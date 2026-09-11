mod cif;
mod formats;
mod geometry;
mod model;
mod render;
mod xyz;

use std::{fs, io::Write, path::PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use formats::InputFormat;
use geometry::Vec3;
use image::{DynamicImage, ImageFormat};
use render::RenderSettings;

#[derive(Debug, Parser)]
#[command(
    name = "xyz-read",
    version,
    about = "Inspect and render molecular, protein, and crystal structures",
    long_about = "Inspect molecular, protein, and crystal structure files and render deterministic PNG images for people and agents.\n\nSupported formats: XYZ, PDB, CIF/mmCIF, MOL/SDF, MOL2, and POSCAR/CONTCAR. Multi-model structures and trajectories are frames numbered from 1. Rendering uses xyz-read's built-in CPU 3D renderer; no GUI or external chemistry software is required.",
    after_help = "QUICK START:\n  xyz-read inspect protein.pdb --json\n  xyz-read render ligand.sdf -o ligand.png --atom-numbers\n  xyz-read render crystal.cif -o cell.png --unit-cell\n  xyz-read render trajectory.xyz -o frame.png --frame 3 --rotate-y 45\n\nCONNECTIVITY:\n  Explicit bonds from PDB CONECT, CIF bond loops, MOL/SDF, and MOL2 are preserved.\n  Missing connectivity is inferred only for formats whose bond tables may be incomplete.\n\nCOMPOSING OPERATIONS:\n  Every image command accepts --format, --frame, --zoom, --rotate-x/y/z, --unit-cell, and --atom-numbers.\n  Use --frame last to render the final model or trajectory frame.\n\nRun `xyz-read <COMMAND> --help` for command-specific options."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Report format, frames, connectivity, protein metadata, and cell data.
    Inspect(InspectArgs),
    /// Render a selected structure frame to PNG.
    Render(RenderArgs),
    /// Render closer by multiplying the current zoom.
    ZoomIn(ZoomArgs),
    /// Render farther away by dividing the current zoom.
    ZoomOut(ZoomArgs),
    /// Render with an explicitly rotated camera view.
    Rotate(RenderArgs),
}

#[derive(Debug, Args)]
struct InspectArgs {
    /// Input structure file.
    input: PathBuf,

    /// Input format. By default, detect it from the name and contents.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,

    /// Inspect one 1-based frame number, or "last".
    #[arg(long, value_name = "NUMBER|last")]
    frame: Option<String>,

    /// Emit stable machine-readable JSON.
    #[arg(long)]
    json: bool,
}

#[derive(Clone, Debug, Args)]
struct RenderArgs {
    /// Input structure file.
    input: PathBuf,

    /// Input format. By default, detect it from the name and contents.
    #[arg(long, value_enum, default_value_t = InputFormat::Auto)]
    format: InputFormat,

    /// Output PNG path. Existing files are replaced atomically.
    #[arg(short, long, value_name = "FILE.png")]
    output: PathBuf,

    /// 1-based frame number, or "last".
    #[arg(long, default_value = "1", value_name = "NUMBER|last")]
    frame: String,

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

    /// Overlay 1-based atom numbers at atom centers.
    #[arg(long)]
    atom_numbers: bool,

    /// Draw the crystallographic unit-cell boundary when cell data exists.
    #[arg(long)]
    unit_cell: bool,

    /// Draw atoms without explicit or automatically inferred bonds.
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

#[derive(Clone, Copy, Debug, ValueEnum)]
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

fn main() {
    if let Err(error) = run(Cli::parse()) {
        eprintln!("xyz-read: {error:#}");
        std::process::exit(1);
    }
}

fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Inspect(args) => inspect(args),
        Command::Render(args) | Command::Rotate(args) => render_image(args, 1.0),
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

fn inspect(args: InspectArgs) -> Result<()> {
    let structure = formats::load(&args.input, args.format)?;
    let selected = match args.frame.as_deref() {
        Some(selector) => Some(structure.frame(selector)?.0),
        None => None,
    };
    let inspection = structure.inspect(selected);

    if args.json {
        println!("{}", serde_json::to_string_pretty(&inspection)?);
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
            "frame {}: {} atom(s), {} bond mode [{}]  {}",
            frame.frame, frame.atom_count, frame.connectivity, elements, frame.comment
        );
        println!(
            "  bounds: x {:.3}..{:.3}, y {:.3}..{:.3}, z {:.3}..{:.3}",
            frame.bounds.min[0],
            frame.bounds.max[0],
            frame.bounds.min[1],
            frame.bounds.max[1],
            frame.bounds.min[2],
            frame.bounds.max[2]
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

fn validate_zoom_factor(factor: f32) -> Result<()> {
    if !factor.is_finite() || factor <= 1.0 {
        bail!("zoom factor must be greater than 1");
    }
    Ok(())
}

fn render_image(args: RenderArgs, zoom_modifier: f32) -> Result<()> {
    if args.output.extension().and_then(|value| value.to_str()) != Some("png") {
        bail!("output must use the .png extension");
    }
    for (name, value) in [
        ("rotate-x", args.rotate_x),
        ("rotate-y", args.rotate_y),
        ("rotate-z", args.rotate_z),
    ] {
        if !value.is_finite() {
            bail!("{name} must be a finite number");
        }
    }

    let structure = formats::load(&args.input, args.format)?;
    let (frame_index, frame) = structure.frame(&args.frame)?;
    let base_rotation = args.view.rotation();
    let image = render::render(
        frame,
        RenderSettings {
            width: args.width,
            height: args.height,
            zoom: args.zoom * zoom_modifier,
            rotation: Vec3::new(
                base_rotation.x + args.rotate_x,
                base_rotation.y + args.rotate_y,
                base_rotation.z + args.rotate_z,
            ),
            atom_numbers: args.atom_numbers,
            bonds: !args.no_bonds,
            unit_cell: args.unit_cell,
        },
    )?;
    write_png_atomically(&args.output, image)?;
    println!(
        "rendered frame {}/{} to {}",
        frame_index + 1,
        structure.frames.len(),
        args.output.display()
    );
    Ok(())
}

fn write_png_atomically(path: &PathBuf, image: image::RgbImage) -> Result<()> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| std::path::Path::new("."));
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
    temporary
        .persist(path)
        .map_err(|error| error.error)
        .with_context(|| format!("could not replace {}", path.display()))?;
    Ok(())
}
