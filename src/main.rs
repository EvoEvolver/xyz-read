mod geometry;
mod render;
mod xyz;

use std::{fs, io::Write, path::PathBuf};

use anyhow::{bail, Context, Result};
use clap::{Args, Parser, Subcommand, ValueEnum};
use geometry::Vec3;
use image::{DynamicImage, ImageFormat};
use render::RenderSettings;
use xyz::Trajectory;

#[derive(Debug, Parser)]
#[command(
    name = "xyz-read",
    version,
    about = "Inspect XYZ trajectories and render molecular images",
    long_about = "Inspect XYZ molecular files and render deterministic PNG images for people and agents.\n\nXYZ trajectories with repeated atom-count/comment/coordinate blocks are supported. Frame numbers are 1-based. Rendering uses xyz-read's built-in CPU 3D renderer; no GUI or external chemistry software is required.",
    after_help = "QUICK START:\n  xyz-read inspect trajectory.xyz --json\n  xyz-read render trajectory.xyz -o frame.png --frame 3 --atom-numbers\n  xyz-read zoom-in molecule.xyz -o close.png --factor 1.5\n  xyz-read rotate molecule.xyz -o turned.png --rotate-y 45 --rotate-x 15\n\nCOMPOSING OPERATIONS:\n  Every image command accepts --frame, --zoom, --rotate-x/y/z, and --atom-numbers.\n  zoom-in multiplies --zoom by --factor; zoom-out divides it by --factor.\n  Use --frame last to render the final trajectory frame.\n\nRun `xyz-read <COMMAND> --help` for command-specific options."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Report frames, atoms, elements, comments, and coordinate bounds.
    Inspect(InspectArgs),
    /// Render a selected XYZ frame to PNG.
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
    /// Input XYZ file. Multi-frame trajectories are supported.
    input: PathBuf,

    /// Inspect one 1-based frame number, or "last".
    #[arg(long, value_name = "NUMBER|last")]
    frame: Option<String>,

    /// Emit stable machine-readable JSON.
    #[arg(long)]
    json: bool,
}

#[derive(Clone, Debug, Args)]
struct RenderArgs {
    /// Input XYZ file. Multi-frame trajectories are supported.
    input: PathBuf,

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

    /// Draw atoms without automatically inferred bonds.
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
    let trajectory = Trajectory::from_path(&args.input)?;
    let selected = match args.frame.as_deref() {
        Some(selector) => Some(trajectory.frame(selector)?.0),
        None => None,
    };
    let inspection = trajectory.inspect(selected);

    if args.json {
        println!("{}", serde_json::to_string_pretty(&inspection)?);
        return Ok(());
    }

    println!("{} frame(s)", inspection.frame_count);
    for frame in inspection.frames {
        let elements = frame
            .elements
            .iter()
            .map(|(element, count)| format!("{element}:{count}"))
            .collect::<Vec<_>>()
            .join(" ");
        println!(
            "frame {}: {} atom(s) [{}]  {}",
            frame.frame, frame.atom_count, elements, frame.comment
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

    let trajectory = Trajectory::from_path(&args.input)?;
    let (frame_index, frame) = trajectory.frame(&args.frame)?;
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
        },
    )?;
    write_png_atomically(&args.output, image)?;
    println!(
        "rendered frame {}/{} to {}",
        frame_index + 1,
        trajectory.frames.len(),
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
