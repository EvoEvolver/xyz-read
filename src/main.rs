mod cif;
mod formats;
mod geometry;
mod model;
mod render;
mod selection;
mod xyz;

use std::{
    fs,
    io::{self, Read, Write},
    path::{Path, PathBuf},
};

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use formats::InputFormat;
use geometry::Vec3;
use image::{DynamicImage, ImageFormat, Rgb, RgbImage};
use model::NumberingMode;
use render::RenderSettings;
use serde::{Deserialize, Serialize};
use serde_json::json;
use sha2::{Digest, Sha256};

const SCHEMA_VERSION: u32 = 1;

#[derive(Debug, Parser)]
#[command(
    name = "xyz-read",
    version,
    about = "Inspect, query, and render molecular structures for people and agents",
    long_about = "Read molecular, protein, and crystal structure files; query atoms without modifying the source; and render deterministic PNG images.\n\nxyz-read is deliberately read-only. Use chemistry software such as RDKit to generate or edit structures, then pass MOL/SDF data to xyz-read by file or stdin for inspection and rendering.",
    after_help = "SUPPORTED INPUT:\n  XYZ, PDB, CIF/mmCIF, MOL/SDF, MOL2, and POSCAR/CONTCAR.\n\nAGENT LOOP:\n  xyz-read inspect protein.pdb --json\n  xyz-read query protein.pdb --select 'around:5@ligand'\n  xyz-read views protein.pdb --select 'around:5@ligand' -o candidates\n  xyz-read plan protein.pdb -o view.json --select ligand --atom-numbers --numbering rdkit\n  xyz-read render --spec view.json -o final.png --report final.json\n  xyz-read render crystal.cif -o cell.png --unit-cell\n  Use --frame last to address the final model or trajectory frame.\n\nRDKIT PIPE:\n  python make_conformer.py | xyz-read render - --format mol -o conformer.png --atom-numbers --numbering rdkit\n\nSELECTIONS:\n  all, protein/polymer, ligand, water, ion\n  index:1,3..8              one-based input indices\n  rdkit-index:0,2..7        zero-based RDKit indices\n  element:C,N  chain:A,B  residue:40..60  atom:CA,N\n  within:5@ligand selects exact atoms; around:5@ligand expands complete residues.\n  Use &, |, !, parentheses, and byres(SELECTION).\n\nRun `xyz-read <COMMAND> --help` for command-specific options."
)]
struct Cli {
    /// Emit operation errors as JSON on stderr.
    #[arg(long, global = true)]
    json_errors: bool,

    /// Emit machine JSON without insignificant whitespace.
    #[arg(long, global = true)]
    compact: bool,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Report format, frames, connectivity, protein metadata, and cell data.
    Inspect(InspectArgs),
    /// Return selected atoms plus explicit and resolved bonds as stable JSON.
    Query(QueryArgs),
    /// Render a structure directly, or reproduce a saved view spec.
    Render(RenderArgs),
    /// Validate a view and save its complete, reproducible JSON specification.
    Plan(PlanArgs),
    /// Create a new view spec by changing camera or selection fields.
    Revise(ReviseArgs),
    /// Generate six labeled candidate views, a contact sheet, and a JSON manifest.
    Views(ViewsArgs),
    /// Describe the machine interface, supported formats, and selectors.
    Capabilities,
    /// Print a bundled JSON Schema for an agent-facing document type.
    Schema(SchemaArgs),
    /// Render closer by multiplying an explicitly supplied zoom.
    ZoomIn(ZoomArgs),
    /// Render farther away by dividing an explicitly supplied zoom.
    ZoomOut(ZoomArgs),
    /// Alias for render with explicit rotation options.
    Rotate(RenderArgs),
}

#[derive(Debug, Args)]
struct InspectArgs {
    /// Input file, or - for stdin.
    input: String,
    /// Input format. Content detection works on stdin, but --format is safer in pipelines.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,
    /// Inspect one 1-based frame number, or "last".
    #[arg(long, value_name = "NUMBER|last")]
    frame: Option<String>,
    /// Emit stable machine-readable JSON.
    #[arg(long)]
    json: bool,
}

#[derive(Debug, Args)]
struct QueryArgs {
    /// Input file, or - for stdin.
    input: String,
    /// Input format. Content detection works on stdin, but --format is safer in pipelines.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,
    /// 1-based frame number, or "last".
    #[arg(long, default_value = "1", value_name = "NUMBER|last")]
    frame: String,
    /// Atom selection expression. Query never changes the source structure.
    #[arg(long, default_value = "all")]
    select: String,
    /// Omit per-atom and per-bond arrays; retain counts and residue/element summaries.
    #[arg(long)]
    summary: bool,
}

#[derive(Clone, Debug, Args)]
struct SceneArgs {
    /// 1-based frame number, or "last".
    #[arg(long, default_value = "1", value_name = "NUMBER|last")]
    frame: String,
    /// Show only atoms matching this expression.
    #[arg(long)]
    select: Option<String>,
    /// Add a contrasting outline to visible atoms matching this expression.
    #[arg(long)]
    highlight: Option<String>,
    /// Output width in pixels (128 to 4096).
    #[arg(long, default_value_t = 1000)]
    width: u32,
    /// Output height in pixels (128 to 4096).
    #[arg(long, default_value_t = 750)]
    height: u32,
    /// Absolute zoom multiplier (0.05 to 20).
    #[arg(long, default_value_t = 1.0)]
    zoom: f32,
    /// Base camera direction before applying rotations.
    #[arg(long, value_enum, default_value_t = View::Iso)]
    view: View,
    /// Rotate around the screen-horizontal X axis, in degrees.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    rotate_x: f32,
    /// Rotate around the screen-vertical Y axis, in degrees.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    rotate_y: f32,
    /// Rotate in the image plane around the Z axis, in degrees.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    rotate_z: f32,
    /// Overlay atom numbers. Use --numbering rdkit for RDKit GetIdx() values.
    #[arg(long)]
    atom_numbers: bool,
    /// Number labels as one-based input positions, RDKit zero-based indices, or source serials.
    #[arg(long, value_enum, default_value_t = Numbering::OneBased)]
    numbering: Numbering,
    /// Draw the crystallographic unit-cell boundary when cell data exists.
    #[arg(long)]
    unit_cell: bool,
    /// Draw atoms without explicit or automatically inferred bonds.
    #[arg(long)]
    no_bonds: bool,
}

#[derive(Clone, Debug, Args)]
struct RenderArgs {
    /// Input file, - for stdin, or omit when using --spec.
    input: Option<String>,
    /// Input format for direct rendering.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,
    /// Reproduce a JSON view spec instead of accepting an input path.
    #[arg(
        long,
        value_name = "VIEW.json",
        conflicts_with_all = [
            "input", "format", "frame", "select", "highlight", "width", "height", "zoom", "view",
            "rotate_x", "rotate_y", "rotate_z", "atom_numbers", "numbering", "unit_cell",
            "no_bonds"
        ]
    )]
    spec: Option<PathBuf>,
    /// Output PNG path. Existing files are replaced atomically.
    #[arg(short, long, value_name = "FILE.png")]
    output: PathBuf,
    /// Write a machine-readable render report containing hashes and resolved counts.
    #[arg(long, value_name = "REPORT.json")]
    report: Option<PathBuf>,
    #[command(flatten)]
    scene: SceneArgs,
}

#[derive(Debug, Args)]
struct PlanArgs {
    /// Input structure path. stdin is rejected because a plan must be reproducible.
    input: String,
    /// Input format. Defaults to file-name and content detection.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,
    /// Destination JSON view spec.
    #[arg(short, long, value_name = "VIEW.json")]
    output: PathBuf,
    #[command(flatten)]
    scene: SceneArgs,
}

#[derive(Debug, Args)]
struct ReviseArgs {
    /// Existing JSON view spec.
    spec: PathBuf,
    /// Destination for the revised spec. The structure file is never modified.
    #[arg(short, long, value_name = "VIEW.json")]
    output: PathBuf,
    /// Multiply the saved absolute zoom by this positive factor.
    #[arg(long)]
    zoom_by: Option<f32>,
    /// Add degrees to the saved X rotation.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    rotate_x_by: f32,
    /// Add degrees to the saved Y rotation.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    rotate_y_by: f32,
    /// Add degrees to the saved Z rotation.
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    rotate_z_by: f32,
    /// Replace the saved selection expression.
    #[arg(long)]
    set_selection: Option<String>,
    /// Replace the saved highlight expression.
    #[arg(long, conflicts_with = "clear_highlight")]
    set_highlight: Option<String>,
    /// Remove highlighting from the saved view.
    #[arg(long)]
    clear_highlight: bool,
}

#[derive(Debug, Args)]
struct ViewsArgs {
    /// Input structure file. Candidate manifests require a stable file path.
    input: String,
    /// Input format. Defaults to file-name and content detection.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,
    /// 1-based frame number, or "last".
    #[arg(long, default_value = "1", value_name = "NUMBER|last")]
    frame: String,
    /// Show only atoms matching this expression in every candidate.
    #[arg(long)]
    select: Option<String>,
    /// Outline matching visible atoms in every candidate.
    #[arg(long)]
    highlight: Option<String>,
    /// Directory receiving A.png through F.png, contact-sheet.png, and manifest.json.
    #[arg(short, long, value_name = "DIRECTORY")]
    output: PathBuf,
    /// Width of each candidate image.
    #[arg(long, default_value_t = 480)]
    width: u32,
    /// Height of each candidate image.
    #[arg(long, default_value_t = 360)]
    height: u32,
    /// Overlay atom numbers on every candidate.
    #[arg(long)]
    atom_numbers: bool,
    /// Number labels using one-based, RDKit, or source conventions.
    #[arg(long, value_enum, default_value_t = Numbering::OneBased)]
    numbering: Numbering,
    /// Draw the crystallographic unit-cell boundary.
    #[arg(long)]
    unit_cell: bool,
    /// Draw atoms without explicit or inferred bonds.
    #[arg(long)]
    no_bonds: bool,
}

#[derive(Debug, Args)]
struct ZoomArgs {
    #[command(flatten)]
    render: RenderArgs,
    /// Zoom change factor. Must be greater than 1.
    #[arg(long, default_value_t = 1.5)]
    factor: f32,
}

#[derive(Debug, Args)]
struct SchemaArgs {
    /// Document schema to print.
    #[arg(value_enum)]
    kind: SchemaKind,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SchemaKind {
    View,
    Query,
    RenderReport,
    CandidateManifest,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, ValueEnum)]
#[serde(rename_all = "lowercase")]
enum View {
    /// Three-quarter view that makes depth immediately visible.
    Iso,
    /// Look along the X axis.
    X,
    /// Look along the Y axis.
    Y,
    /// Look along the Z axis with X horizontal and Y vertical.
    Z,
}

impl View {
    fn rotation(self) -> Vec3 {
        match self {
            Self::Iso => Vec3::new(-20.0, 35.0, 0.0),
            Self::X => Vec3::new(0.0, 90.0, 0.0),
            Self::Y => Vec3::new(-90.0, 0.0, 0.0),
            Self::Z => Vec3::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, ValueEnum)]
#[serde(rename_all = "kebab-case")]
enum Numbering {
    /// Stable one-based position in the input frame.
    OneBased,
    /// RDKit-compatible zero-based atom index (Atom.GetIdx()).
    Rdkit,
    /// PDB/mmCIF source serial, falling back to one-based input position.
    Source,
}

impl From<Numbering> for NumberingMode {
    fn from(value: Numbering) -> Self {
        match value {
            Numbering::OneBased => Self::OneBased,
            Numbering::Rdkit => Self::Rdkit,
            Numbering::Source => Self::Source,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct ViewSpec {
    schema_version: u32,
    tool_version: String,
    input: SpecInput,
    frame: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    selection: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    highlight: Option<String>,
    camera: CameraSpec,
    image: ImageSpec,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct SpecInput {
    path: String,
    format: String,
    sha256: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct CameraSpec {
    view: View,
    rotation: [f32; 3],
    zoom: f32,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct ImageSpec {
    width: u32,
    height: u32,
    atom_numbers: bool,
    numbering: Numbering,
    unit_cell: bool,
    bonds: bool,
}

#[derive(Debug, Serialize)]
struct RenderReport {
    schema_version: u32,
    tool_version: &'static str,
    input_sha256: String,
    frame: usize,
    frame_count: usize,
    selection: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    highlight: Option<String>,
    selected_atoms: usize,
    drawn_bonds: usize,
    camera: CameraSpec,
    image: OutputImage,
    warnings: Vec<String>,
}

#[derive(Debug, Serialize)]
struct OutputImage {
    path: String,
    sha256: String,
    width: u32,
    height: u32,
}

#[derive(Debug, Serialize)]
struct CandidateManifest {
    schema_version: u32,
    tool_version: &'static str,
    input: SpecInput,
    frame: usize,
    selection: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    highlight: Option<String>,
    contact_sheet: String,
    candidates: Vec<Candidate>,
}

#[derive(Debug, Serialize)]
struct Candidate {
    label: &'static str,
    file: String,
    spec: String,
    view: View,
    rotation: [f32; 3],
}

struct LoadedInput {
    structure: model::Structure,
    sha256: String,
}

fn main() {
    let cli = Cli::parse();
    let json_errors = cli.json_errors;
    if let Err(error) = run(cli.command, cli.compact) {
        if json_errors {
            eprintln!(
                "{}",
                json!({
                    "schema_version": SCHEMA_VERSION,
                    "error": {"code": "operation_failed", "message": format!("{error:#}")}
                })
            );
        } else {
            eprintln!("xyz-read: {error:#}");
        }
        std::process::exit(1);
    }
}

fn run(command: Command, compact: bool) -> Result<()> {
    match command {
        Command::Inspect(args) => inspect(args, compact),
        Command::Query(args) => query(args, compact),
        Command::Render(args) | Command::Rotate(args) => render_image(args, 1.0),
        Command::Plan(args) => plan(args),
        Command::Revise(args) => revise(args),
        Command::Views(args) => candidate_views(args),
        Command::Capabilities => capabilities(compact),
        Command::Schema(args) => schema(args, compact),
        Command::ZoomIn(args) => {
            validate_zoom_factor(args.factor)?;
            render_image(args.render, args.factor)
        }
        Command::ZoomOut(args) => {
            validate_zoom_factor(args.factor)?;
            render_image(args.render, 1.0 / args.factor)
        }
    }
}

fn load_input(input: &str, format: InputFormat) -> Result<LoadedInput> {
    let source = if input == "-" {
        let mut source = String::new();
        io::stdin()
            .read_to_string(&mut source)
            .context("could not read structure data from stdin")?;
        source
    } else {
        fs::read_to_string(input)
            .with_context(|| format!("could not read structure file {input}"))?
    };
    let hint = (input != "-").then(|| Path::new(input));
    let structure = formats::load_text(&source, hint, format)?;
    Ok(LoadedInput {
        structure,
        sha256: sha256(source.as_bytes()),
    })
}

fn inspect(args: InspectArgs, compact: bool) -> Result<()> {
    let loaded = load_input(&args.input, args.format)?;
    let selected = match args.frame.as_deref() {
        Some(selector) => Some(loaded.structure.frame(selector)?.0),
        None => None,
    };
    let mut inspection = loaded.structure.inspect(selected);
    inspection.input_sha256 = loaded.sha256;
    if args.json {
        print_json(&inspection, compact)?;
        return Ok(());
    }
    println!("format: {}", inspection.format);
    println!("{} frame(s)", inspection.frame_count);
    for frame in inspection.frames {
        let elements = frame
            .elements
            .iter()
            .map(|(element, count)| format!("{element}:{count}"))
            .collect::<Vec<_>>()
            .join(" ");
        println!(
            "frame {}: {} atom(s), {} connectivity [{}]  {}",
            frame.frame, frame.atom_count, frame.connectivity, elements, frame.comment
        );
        if frame.residue_count > 0 || !frame.chains.is_empty() {
            println!(
                "  protein: {} residue(s), chain(s) {}",
                frame.residue_count,
                frame.chains.join(",")
            );
        }
        if let Some(cell) = frame.unit_cell {
            println!(
                "  cell: {:.3} {:.3} {:.3} A; {:.2} {:.2} {:.2} degrees",
                cell.lengths[0],
                cell.lengths[1],
                cell.lengths[2],
                cell.angles[0],
                cell.angles[1],
                cell.angles[2]
            );
        }
    }
    Ok(())
}

fn query(args: QueryArgs, compact: bool) -> Result<()> {
    let loaded = load_input(&args.input, args.format)?;
    let (frame_index, frame) = loaded.structure.frame(&args.frame)?;
    let resolved = render::resolved_bonds(frame);
    let result = selection::query_result(
        &loaded.structure.format,
        frame_index + 1,
        frame,
        &args.select,
        &loaded.sha256,
        &resolved,
        !args.summary,
    )?;
    print_json(&result, compact)
}

fn plan(args: PlanArgs) -> Result<()> {
    if args.input == "-" {
        bail!("plan requires a stable input path; stdin cannot be stored reproducibly");
    }
    let canonical = canonical_input_path(&args.input)?;
    let loaded = load_input(&canonical, args.format)?;
    validate_scene(&args.scene)?;
    let (_, frame) = loaded.structure.frame(&args.scene.frame)?;
    validate_selection(
        frame,
        args.scene.select.as_deref(),
        args.scene.highlight.as_deref(),
    )?;
    let spec = view_spec(&canonical, &loaded, args.scene);
    write_json_atomically(&args.output, &spec)?;
    println!("wrote reproducible view spec to {}", args.output.display());
    Ok(())
}

fn revise(args: ReviseArgs) -> Result<()> {
    let mut spec = read_spec(&args.spec)?;
    if let Some(factor) = args.zoom_by {
        if !factor.is_finite() || factor <= 0.0 {
            bail!("--zoom-by must be finite and greater than zero");
        }
        spec.camera.zoom *= factor;
    }
    for (name, value) in [
        ("rotate-x-by", args.rotate_x_by),
        ("rotate-y-by", args.rotate_y_by),
        ("rotate-z-by", args.rotate_z_by),
    ] {
        if !value.is_finite() {
            bail!("{name} must be finite");
        }
    }
    spec.camera.rotation[0] += args.rotate_x_by;
    spec.camera.rotation[1] += args.rotate_y_by;
    spec.camera.rotation[2] += args.rotate_z_by;
    if let Some(selection) = args.set_selection {
        spec.selection = Some(selection);
    }
    if let Some(highlight) = args.set_highlight {
        spec.highlight = Some(highlight);
    } else if args.clear_highlight {
        spec.highlight = None;
    }
    validate_spec(&spec)?;
    let format = parse_format(&spec.input.format)?;
    let loaded = load_input(&spec.input.path, format)?;
    verify_input_hash(&spec, &loaded)?;
    let (_, frame) = loaded.structure.frame(&spec.frame)?;
    validate_selection(frame, spec.selection.as_deref(), spec.highlight.as_deref())?;
    spec.tool_version = env!("CARGO_PKG_VERSION").to_owned();
    write_json_atomically(&args.output, &spec)?;
    println!("wrote revised view spec to {}", args.output.display());
    Ok(())
}

fn render_image(args: RenderArgs, zoom_modifier: f32) -> Result<()> {
    require_png(&args.output)?;
    let (loaded, scene, stored_camera) = if let Some(path) = &args.spec {
        let spec = read_spec(path)?;
        let format = parse_format(&spec.input.format)?;
        let loaded = load_input(&spec.input.path, format)?;
        verify_input_hash(&spec, &loaded)?;
        let scene = scene_from_spec(&spec);
        (loaded, scene, Some(spec.camera))
    } else {
        let input = args
            .input
            .as_deref()
            .context("render requires an INPUT or --spec VIEW.json")?;
        (load_input(input, args.format)?, args.scene.clone(), None)
    };
    validate_scene(&scene)?;
    let (frame_index, original) = loaded.structure.frame(&scene.frame)?;
    let selected = selected_frame(
        original,
        scene.select.as_deref(),
        scene.highlight.as_deref(),
    )?;
    let rotation =
        scene.view.rotation() + Vec3::new(scene.rotate_x, scene.rotate_y, scene.rotate_z);
    let image = render::render(
        &selected,
        RenderSettings {
            width: scene.width,
            height: scene.height,
            zoom: scene.zoom * zoom_modifier,
            rotation,
            atom_numbers: scene.atom_numbers.then(|| scene.numbering.into()),
            bonds: !scene.no_bonds,
            unit_cell: scene.unit_cell,
        },
    )?;
    let drawn_bonds = if scene.no_bonds {
        0
    } else {
        render::resolved_bonds(&selected).len()
    };
    write_png_atomically(&args.output, image)?;
    if let Some(report_path) = args.report {
        let mut camera = stored_camera.unwrap_or(CameraSpec {
            view: scene.view,
            rotation: [scene.rotate_x, scene.rotate_y, scene.rotate_z],
            zoom: scene.zoom,
        });
        camera.zoom *= zoom_modifier;
        let report = RenderReport {
            schema_version: SCHEMA_VERSION,
            tool_version: env!("CARGO_PKG_VERSION"),
            input_sha256: loaded.sha256,
            frame: frame_index + 1,
            frame_count: loaded.structure.frames.len(),
            selection: scene.select.unwrap_or_else(|| "all".to_owned()),
            highlight: scene.highlight,
            selected_atoms: selected.atoms.len(),
            drawn_bonds,
            camera,
            image: OutputImage {
                path: args.output.display().to_string(),
                sha256: sha256(&fs::read(&args.output)?),
                width: scene.width,
                height: scene.height,
            },
            warnings: Vec::new(),
        };
        write_json_atomically(&report_path, &report)?;
    }
    println!(
        "rendered {} selected atom(s), frame {}/{} to {}",
        selected.atoms.len(),
        frame_index + 1,
        loaded.structure.frames.len(),
        args.output.display()
    );
    Ok(())
}

fn candidate_views(args: ViewsArgs) -> Result<()> {
    if args.input == "-" {
        bail!("views requires a stable input path so its manifest can be reproduced");
    }
    if !(128..=2048).contains(&args.width) || !(128..=2048).contains(&args.height) {
        bail!("candidate width and height must each be between 128 and 2048 pixels");
    }
    let canonical = canonical_input_path(&args.input)?;
    let loaded = load_input(&canonical, args.format)?;
    let (frame_index, original) = loaded.structure.frame(&args.frame)?;
    let selected = selected_frame(original, args.select.as_deref(), args.highlight.as_deref())?;
    fs::create_dir_all(&args.output).with_context(|| {
        format!(
            "could not create candidate directory {}",
            args.output.display()
        )
    })?;
    let definitions = [
        ("A", View::Iso, [0.0, 0.0, 0.0]),
        ("B", View::X, [0.0, 0.0, 0.0]),
        ("C", View::Y, [0.0, 0.0, 0.0]),
        ("D", View::Z, [0.0, 0.0, 0.0]),
        ("E", View::Iso, [0.0, 90.0, 0.0]),
        ("F", View::Iso, [90.0, 0.0, 0.0]),
    ];
    let mut sheet = RgbImage::from_pixel(args.width * 3, args.height * 2, Rgb([247, 249, 252]));
    let mut candidates = Vec::new();
    for (index, (label, view, extra)) in definitions.into_iter().enumerate() {
        let base = view.rotation();
        let rotation = Vec3::new(base.x + extra[0], base.y + extra[1], base.z + extra[2]);
        let mut image = render::render(
            &selected,
            RenderSettings {
                width: args.width,
                height: args.height,
                zoom: 1.0,
                rotation,
                atom_numbers: args.atom_numbers.then(|| args.numbering.into()),
                bonds: !args.no_bonds,
                unit_cell: args.unit_cell,
            },
        )?;
        render::draw_panel_label(&mut image, label);
        let file_name = format!("{label}.png");
        let spec_name = format!("{label}.view.json");
        write_png_atomically(&args.output.join(&file_name), image.clone())?;
        let candidate_scene = SceneArgs {
            frame: args.frame.clone(),
            select: args.select.clone(),
            highlight: args.highlight.clone(),
            width: args.width,
            height: args.height,
            zoom: 1.0,
            view,
            rotate_x: extra[0],
            rotate_y: extra[1],
            rotate_z: extra[2],
            atom_numbers: args.atom_numbers,
            numbering: args.numbering,
            unit_cell: args.unit_cell,
            no_bonds: args.no_bonds,
        };
        write_json_atomically(
            &args.output.join(&spec_name),
            &view_spec(&canonical, &loaded, candidate_scene),
        )?;
        image::imageops::replace(
            &mut sheet,
            &image,
            ((index % 3) as u32 * args.width).into(),
            ((index / 3) as u32 * args.height).into(),
        );
        candidates.push(Candidate {
            label,
            file: file_name,
            spec: spec_name,
            view,
            rotation: extra,
        });
    }
    let contact_name = "contact-sheet.png";
    write_png_atomically(&args.output.join(contact_name), sheet)?;
    let manifest = CandidateManifest {
        schema_version: SCHEMA_VERSION,
        tool_version: env!("CARGO_PKG_VERSION"),
        input: SpecInput {
            path: canonical,
            format: loaded.structure.format.clone(),
            sha256: loaded.sha256,
        },
        frame: frame_index + 1,
        selection: args.select.unwrap_or_else(|| "all".to_owned()),
        highlight: args.highlight,
        contact_sheet: contact_name.to_owned(),
        candidates,
    };
    write_json_atomically(&args.output.join("manifest.json"), &manifest)?;
    println!("generated candidate views in {}", args.output.display());
    Ok(())
}

fn capabilities(compact: bool) -> Result<()> {
    print_json(
        &json!({
            "schema_version": SCHEMA_VERSION,
            "tool_version": env!("CARGO_PKG_VERSION"),
            "read_only": true,
            "structure_editing": false,
            "commands": ["inspect", "query", "render", "plan", "revise", "views", "capabilities", "schema"],
            "formats": ["xyz", "pdb", "cif", "mmcif", "mol", "sdf", "mol2", "poscar"],
            "stdin": true,
            "selectors": {
                "roles": ["all", "protein", "polymer", "ligand", "water", "ion"],
                "properties": ["index", "rdkit-index", "element", "chain", "residue", "atom"],
                "operators": ["&", "|", "!", "()", "within:DISTANCE@SELECTION", "around:DISTANCE@SELECTION", "byres(SELECTION)"]
            },
            "atom_numbering": ["one-based", "rdkit", "source"],
            "render_features": ["selection", "highlight", "unit-cell", "atom-numbers", "candidate-views"],
            "query_modes": ["detailed", "--summary"],
            "view_specs": {"schema_version": SCHEMA_VERSION, "input_hash_verified": true},
            "schemas": {
                "view": "xyz-read schema view",
                "query": "xyz-read schema query",
                "render_report": "xyz-read schema render-report",
                "candidate_manifest": "xyz-read schema candidate-manifest"
            },
            "rdkit": {
                "boundary": "MOL/SDF file or stdin",
                "atom_order_preserved": true,
                "explicit_bond_orders_preserved": true,
                "indexing": "use rdkit-index selectors and --numbering rdkit"
            }
        }),
        compact,
    )
}

fn schema(args: SchemaArgs, compact: bool) -> Result<()> {
    let source = match args.kind {
        SchemaKind::View => include_str!("../schemas/view-v1.schema.json"),
        SchemaKind::Query => include_str!("../schemas/query-v1.schema.json"),
        SchemaKind::RenderReport => include_str!("../schemas/render-report-v1.schema.json"),
        SchemaKind::CandidateManifest => {
            include_str!("../schemas/candidate-manifest-v1.schema.json")
        }
    };
    if compact {
        let value: serde_json::Value = serde_json::from_str(source)?;
        print_json(&value, true)
    } else {
        print!("{source}");
        Ok(())
    }
}

fn print_json(value: &impl Serialize, compact: bool) -> Result<()> {
    if compact {
        println!("{}", serde_json::to_string(value)?);
    } else {
        println!("{}", serde_json::to_string_pretty(value)?);
    }
    Ok(())
}

fn selected_frame(
    frame: &model::Frame,
    expression: Option<&str>,
    highlight: Option<&str>,
) -> Result<model::Frame> {
    let expression = expression.unwrap_or("all");
    let selected = selection::select(frame, expression)?;
    let mut subset = selection::subset(frame, &selected)?;
    if let Some(highlight) = highlight {
        let highlighted = selection::select(frame, highlight)?;
        let mut visible_count = 0;
        for atom in &mut subset.atoms {
            atom.highlighted = highlighted.contains(&atom.input_index);
            visible_count += usize::from(atom.highlighted);
        }
        if visible_count == 0 {
            bail!("highlight expression matched no visible atoms");
        }
    }
    Ok(subset)
}

fn validate_selection(
    frame: &model::Frame,
    expression: Option<&str>,
    highlight: Option<&str>,
) -> Result<()> {
    selected_frame(frame, expression, highlight).map(|_| ())
}

fn validate_scene(scene: &SceneArgs) -> Result<()> {
    if !(128..=4096).contains(&scene.width) || !(128..=4096).contains(&scene.height) {
        bail!("image width and height must each be between 128 and 4096 pixels");
    }
    if !scene.zoom.is_finite() || !(0.05..=20.0).contains(&scene.zoom) {
        bail!("zoom must be between 0.05 and 20");
    }
    for (name, value) in [
        ("rotate-x", scene.rotate_x),
        ("rotate-y", scene.rotate_y),
        ("rotate-z", scene.rotate_z),
    ] {
        if !value.is_finite() {
            bail!("{name} must be finite");
        }
    }
    Ok(())
}

fn view_spec(path: &str, loaded: &LoadedInput, scene: SceneArgs) -> ViewSpec {
    ViewSpec {
        schema_version: SCHEMA_VERSION,
        tool_version: env!("CARGO_PKG_VERSION").to_owned(),
        input: SpecInput {
            path: path.to_owned(),
            format: loaded.structure.format.clone(),
            sha256: loaded.sha256.clone(),
        },
        frame: scene.frame,
        selection: scene.select,
        highlight: scene.highlight,
        camera: CameraSpec {
            view: scene.view,
            rotation: [scene.rotate_x, scene.rotate_y, scene.rotate_z],
            zoom: scene.zoom,
        },
        image: ImageSpec {
            width: scene.width,
            height: scene.height,
            atom_numbers: scene.atom_numbers,
            numbering: scene.numbering,
            unit_cell: scene.unit_cell,
            bonds: !scene.no_bonds,
        },
    }
}

fn scene_from_spec(spec: &ViewSpec) -> SceneArgs {
    SceneArgs {
        frame: spec.frame.clone(),
        select: spec.selection.clone(),
        highlight: spec.highlight.clone(),
        width: spec.image.width,
        height: spec.image.height,
        zoom: spec.camera.zoom,
        view: spec.camera.view,
        rotate_x: spec.camera.rotation[0],
        rotate_y: spec.camera.rotation[1],
        rotate_z: spec.camera.rotation[2],
        atom_numbers: spec.image.atom_numbers,
        numbering: spec.image.numbering,
        unit_cell: spec.image.unit_cell,
        no_bonds: !spec.image.bonds,
    }
}

fn read_spec(path: &Path) -> Result<ViewSpec> {
    let source = fs::read_to_string(path)
        .with_context(|| format!("could not read view spec {}", path.display()))?;
    let spec: ViewSpec = serde_json::from_str(&source)
        .with_context(|| format!("invalid JSON view spec {}", path.display()))?;
    validate_spec(&spec)?;
    Ok(spec)
}

fn validate_spec(spec: &ViewSpec) -> Result<()> {
    if spec.schema_version != SCHEMA_VERSION {
        bail!(
            "unsupported view spec schema {}; this xyz-read supports schema {}",
            spec.schema_version,
            SCHEMA_VERSION
        );
    }
    validate_scene(&scene_from_spec(spec))
}

fn verify_input_hash(spec: &ViewSpec, loaded: &LoadedInput) -> Result<()> {
    if spec.input.sha256 != loaded.sha256 {
        bail!(
            "input hash mismatch for {}; the structure changed after this view was planned",
            spec.input.path
        );
    }
    Ok(())
}

fn parse_format(value: &str) -> Result<InputFormat> {
    InputFormat::from_str(value, true).map_err(|_| anyhow::anyhow!("unknown format {value:?}"))
}

fn canonical_input_path(input: &str) -> Result<String> {
    let path = fs::canonicalize(input)
        .with_context(|| format!("could not resolve structure path {input}"))?;
    path.to_str()
        .map(str::to_owned)
        .context("structure path is not valid UTF-8")
}

fn validate_zoom_factor(factor: f32) -> Result<()> {
    if !factor.is_finite() || factor <= 1.0 {
        bail!("zoom factor must be greater than 1");
    }
    Ok(())
}

fn require_png(path: &Path) -> Result<()> {
    if path.extension().and_then(|value| value.to_str()) != Some("png") {
        bail!("output must use the .png extension");
    }
    Ok(())
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn write_json_atomically(path: &Path, value: &impl Serialize) -> Result<()> {
    let mut bytes = serde_json::to_vec_pretty(value)?;
    bytes.push(b'\n');
    write_bytes_atomically(path, &bytes)
}

fn write_png_atomically(path: &Path, image: RgbImage) -> Result<()> {
    require_png(path)?;
    let parent = output_parent(path);
    fs::create_dir_all(parent)
        .with_context(|| format!("could not create output directory {}", parent.display()))?;
    let mut temporary = tempfile::Builder::new()
        .prefix(".xyz-read-")
        .suffix(".png")
        .tempfile_in(parent)
        .with_context(|| format!("could not create a temporary file in {}", parent.display()))?;
    DynamicImage::ImageRgb8(image)
        .write_to(temporary.as_file_mut(), ImageFormat::Png)
        .with_context(|| format!("could not encode PNG for {}", path.display()))?;
    temporary.as_file_mut().flush()?;
    persist(temporary, path)
}

fn write_bytes_atomically(path: &Path, bytes: &[u8]) -> Result<()> {
    let parent = output_parent(path);
    fs::create_dir_all(parent)
        .with_context(|| format!("could not create output directory {}", parent.display()))?;
    let mut temporary = tempfile::Builder::new()
        .prefix(".xyz-read-")
        .tempfile_in(parent)
        .with_context(|| format!("could not create a temporary file in {}", parent.display()))?;
    temporary.write_all(bytes)?;
    temporary.as_file_mut().flush()?;
    persist(temporary, path)
}

fn persist(temporary: tempfile::NamedTempFile, path: &Path) -> Result<()> {
    temporary
        .persist(path)
        .map_err(|error| error.error)
        .with_context(|| format!("could not replace {}", path.display()))?;
    Ok(())
}

fn output_parent(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}
