use std::collections::HashMap;

use anyhow::{bail, Result};
use image::{imageops::FilterType, Rgb, RgbImage};

use crate::geometry::Vec3;
use crate::{geometry::rotate_euler, xyz::Frame};

const SUPERSAMPLE: u32 = 2;
const BACKGROUND: [u8; 3] = [247, 249, 252];

#[derive(Clone, Copy, Debug)]
pub struct RenderSettings {
    pub width: u32,
    pub height: u32,
    pub zoom: f32,
    pub rotation: Vec3,
    pub atom_numbers: bool,
    pub bonds: bool,
}

#[derive(Clone, Copy, Debug)]
struct ElementStyle {
    color: [u8; 3],
    radius: f32,
    covalent_radius: f32,
}

#[derive(Clone, Copy, Debug)]
struct ProjectedAtom {
    center: Vec3,
    radius_px: f32,
    color: [u8; 3],
    index: usize,
}

#[derive(Clone, Copy, Debug)]
struct Bond {
    first: usize,
    second: usize,
}

pub fn render(frame: &Frame, settings: RenderSettings) -> Result<RgbImage> {
    if frame.atoms.is_empty() {
        bail!("cannot render an empty XYZ frame");
    }
    if !(128..=4096).contains(&settings.width) || !(128..=4096).contains(&settings.height) {
        bail!("image width and height must each be between 128 and 4096 pixels");
    }
    if !settings.zoom.is_finite() || !(0.05..=20.0).contains(&settings.zoom) {
        bail!("zoom must be between 0.05 and 20");
    }

    let render_width = settings.width * SUPERSAMPLE;
    let render_height = settings.height * SUPERSAMPLE;
    let center = coordinate_center(frame);
    let transformed: Vec<Vec3> = frame
        .atoms
        .iter()
        .map(|atom| {
            rotate_euler(
                atom.position - center,
                settings.rotation.x,
                settings.rotation.y,
                settings.rotation.z,
            )
        })
        .collect();
    let scale = fit_scale(frame, &transformed, render_width, render_height) * settings.zoom;
    let screen_center = Vec3::new(render_width as f32 / 2.0, render_height as f32 / 2.0, 0.0);
    let projected: Vec<ProjectedAtom> = frame
        .atoms
        .iter()
        .zip(&transformed)
        .enumerate()
        .map(|(index, (atom, point))| {
            let style = element_style(&atom.element);
            ProjectedAtom {
                center: Vec3::new(
                    screen_center.x + point.x * scale,
                    screen_center.y - point.y * scale,
                    point.z,
                ),
                radius_px: (style.radius * scale).max(4.0 * SUPERSAMPLE as f32),
                color: style.color,
                index,
            }
        })
        .collect();

    let mut image = RgbImage::from_pixel(render_width, render_height, Rgb(BACKGROUND));
    let mut depth = vec![f32::NEG_INFINITY; (render_width * render_height) as usize];

    if settings.bonds {
        for bond in infer_bonds(frame) {
            draw_bond(
                &mut image,
                &mut depth,
                projected[bond.first],
                projected[bond.second],
                scale,
            );
        }
    }
    for atom in &projected {
        draw_sphere(&mut image, &mut depth, *atom, scale);
    }

    let mut image = image::imageops::resize(
        &image,
        settings.width,
        settings.height,
        FilterType::Lanczos3,
    );
    if settings.atom_numbers {
        let mut label_positions = projected;
        label_positions.sort_by(|left, right| left.center.z.total_cmp(&right.center.z));
        for atom in label_positions {
            draw_number(
                &mut image,
                atom.index + 1,
                atom.center.x / SUPERSAMPLE as f32,
                atom.center.y / SUPERSAMPLE as f32,
                atom.radius_px / SUPERSAMPLE as f32,
            );
        }
    }
    Ok(image)
}

fn coordinate_center(frame: &Frame) -> Vec3 {
    let mut min = frame.atoms[0].position;
    let mut max = min;
    for atom in &frame.atoms[1..] {
        min.x = min.x.min(atom.position.x);
        min.y = min.y.min(atom.position.y);
        min.z = min.z.min(atom.position.z);
        max.x = max.x.max(atom.position.x);
        max.y = max.y.max(atom.position.y);
        max.z = max.z.max(atom.position.z);
    }
    (min + max) / 2.0
}

fn fit_scale(frame: &Frame, points: &[Vec3], width: u32, height: u32) -> f32 {
    let mut min_x = f32::INFINITY;
    let mut max_x = f32::NEG_INFINITY;
    let mut min_y = f32::INFINITY;
    let mut max_y = f32::NEG_INFINITY;
    for (atom, point) in frame.atoms.iter().zip(points) {
        let radius = element_style(&atom.element).radius;
        min_x = min_x.min(point.x - radius);
        max_x = max_x.max(point.x + radius);
        min_y = min_y.min(point.y - radius);
        max_y = max_y.max(point.y + radius);
    }
    let span_x = (max_x - min_x).max(0.5);
    let span_y = (max_y - min_y).max(0.5);
    ((width as f32 * 0.76) / span_x).min((height as f32 * 0.76) / span_y)
}

fn infer_bonds(frame: &Frame) -> Vec<Bond> {
    const CELL_SIZE: f32 = 3.3;
    let mut cells: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    let cell_for = |point: Vec3| {
        (
            (point.x / CELL_SIZE).floor() as i32,
            (point.y / CELL_SIZE).floor() as i32,
            (point.z / CELL_SIZE).floor() as i32,
        )
    };
    for (index, atom) in frame.atoms.iter().enumerate() {
        cells
            .entry(cell_for(atom.position))
            .or_default()
            .push(index);
    }

    let mut bonds = Vec::new();
    for (first, atom) in frame.atoms.iter().enumerate() {
        let cell = cell_for(atom.position);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbor = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                    let Some(indices) = cells.get(&neighbor) else {
                        continue;
                    };
                    for &second in indices {
                        if second <= first {
                            continue;
                        }
                        let other = &frame.atoms[second];
                        let distance = (atom.position - other.position).length();
                        let threshold = (element_style(&atom.element).covalent_radius
                            + element_style(&other.element).covalent_radius)
                            * 1.25;
                        if distance >= 0.1 && distance <= threshold.min(CELL_SIZE) {
                            bonds.push(Bond { first, second });
                        }
                    }
                }
            }
        }
    }
    bonds
}

fn draw_bond(
    image: &mut RgbImage,
    depth: &mut [f32],
    first: ProjectedAtom,
    second: ProjectedAtom,
    scale: f32,
) {
    let midpoint = Vec3::new(
        (first.center.x + second.center.x) / 2.0,
        (first.center.y + second.center.y) / 2.0,
        (first.center.z + second.center.z) / 2.0,
    );
    draw_bond_half(image, depth, first.center, midpoint, first.color, scale);
    draw_bond_half(image, depth, midpoint, second.center, second.color, scale);
}

fn draw_bond_half(
    image: &mut RgbImage,
    depth: &mut [f32],
    start: Vec3,
    end: Vec3,
    color: [u8; 3],
    scale: f32,
) {
    let radius_px = (scale * 0.115).max(2.5 * SUPERSAMPLE as f32);
    let Some((min_x, max_x, min_y, max_y)) = clipped_bounds(
        start.x.min(end.x) - radius_px,
        start.x.max(end.x) + radius_px,
        start.y.min(end.y) - radius_px,
        start.y.max(end.y) + radius_px,
        image.width(),
        image.height(),
    ) else {
        return;
    };
    let delta_x = end.x - start.x;
    let delta_y = end.y - start.y;
    let length_squared = delta_x * delta_x + delta_y * delta_y;

    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let t = if length_squared < 0.001 {
                0.0
            } else {
                (((x as f32 - start.x) * delta_x + (y as f32 - start.y) * delta_y) / length_squared)
                    .clamp(0.0, 1.0)
            };
            let closest_x = start.x + delta_x * t;
            let closest_y = start.y + delta_y * t;
            let pixel_dx = x as f32 + 0.5 - closest_x;
            let pixel_dy = y as f32 + 0.5 - closest_y;
            let distance_squared = pixel_dx * pixel_dx + pixel_dy * pixel_dy;
            if distance_squared > radius_px * radius_px {
                continue;
            }
            let curved_front = (radius_px * radius_px - distance_squared).sqrt() / scale;
            let pixel_depth = start.z + (end.z - start.z) * t + curved_front;
            let offset = (y * image.width() + x) as usize;
            if pixel_depth > depth[offset] {
                depth[offset] = pixel_depth;
                let edge = (1.0 - distance_squared / (radius_px * radius_px)).sqrt();
                *image.get_pixel_mut(x, y) = Rgb(shade(color, 0.48 + 0.42 * edge));
            }
        }
    }
}

fn draw_sphere(image: &mut RgbImage, depth: &mut [f32], atom: ProjectedAtom, scale: f32) {
    let radius = atom.radius_px;
    let Some((min_x, max_x, min_y, max_y)) = clipped_bounds(
        atom.center.x - radius,
        atom.center.x + radius,
        atom.center.y - radius,
        atom.center.y + radius,
        image.width(),
        image.height(),
    ) else {
        return;
    };
    let light = Vec3::new(-0.42, 0.52, 0.74).normalized();

    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let nx = (x as f32 + 0.5 - atom.center.x) / radius;
            let ny = -(y as f32 + 0.5 - atom.center.y) / radius;
            let radial = nx * nx + ny * ny;
            if radial > 1.0 {
                continue;
            }
            let nz = (1.0 - radial).sqrt();
            let pixel_depth = atom.center.z + (nz * radius / scale);
            let offset = (y * image.width() + x) as usize;
            if pixel_depth <= depth[offset] {
                continue;
            }
            depth[offset] = pixel_depth;
            let normal = Vec3::new(nx, ny, nz);
            let diffuse = normal.dot(light).max(0.0);
            let rim = (1.0 - nz).powi(2);
            let highlight = diffuse.powi(18) * 0.34;
            let intensity = 0.29 + diffuse * 0.68 + highlight - rim * 0.08;
            *image.get_pixel_mut(x, y) = Rgb(shade(atom.color, intensity));
        }
    }
}

fn clipped_bounds(
    min_x: f32,
    max_x: f32,
    min_y: f32,
    max_y: f32,
    width: u32,
    height: u32,
) -> Option<(u32, u32, u32, u32)> {
    if max_x < 0.0 || max_y < 0.0 || min_x >= width as f32 || min_y >= height as f32 {
        return None;
    }
    Some((
        min_x.floor().max(0.0) as u32,
        max_x.ceil().min(width as f32 - 1.0) as u32,
        min_y.floor().max(0.0) as u32,
        max_y.ceil().min(height as f32 - 1.0) as u32,
    ))
}

fn shade(color: [u8; 3], intensity: f32) -> [u8; 3] {
    let intensity = intensity.clamp(0.0, 1.18);
    color.map(|channel| ((channel as f32 * intensity).min(255.0)) as u8)
}

fn draw_number(image: &mut RgbImage, number: usize, center_x: f32, center_y: f32, radius: f32) {
    let text = number.to_string();
    let scale = (radius / 8.5).round().clamp(1.0, 4.0) as i32;
    let glyph_width = 5 * scale;
    let gap = scale;
    let total_width = text.len() as i32 * glyph_width + (text.len() as i32 - 1) * gap;
    let start_x = center_x.round() as i32 - total_width / 2;
    let start_y = center_y.round() as i32 - (7 * scale) / 2;

    let mut pixels = Vec::new();
    for (digit_index, character) in text.chars().enumerate() {
        let glyph = digit_glyph(character);
        let origin_x = start_x + digit_index as i32 * (glyph_width + gap);
        for (row, bits) in glyph.iter().enumerate() {
            for column in 0..5 {
                if bits & (1 << (4 - column)) == 0 {
                    continue;
                }
                for sy in 0..scale {
                    for sx in 0..scale {
                        pixels.push((
                            origin_x + column * scale + sx,
                            start_y + row as i32 * scale + sy,
                        ));
                    }
                }
            }
        }
    }

    for &(x, y) in &pixels {
        for oy in -1..=1 {
            for ox in -1..=1 {
                set_pixel_checked(image, x + ox, y + oy, [20, 24, 30]);
            }
        }
    }
    for (x, y) in pixels {
        set_pixel_checked(image, x, y, [255, 255, 255]);
    }
}

fn set_pixel_checked(image: &mut RgbImage, x: i32, y: i32, color: [u8; 3]) {
    if x >= 0 && y >= 0 && x < image.width() as i32 && y < image.height() as i32 {
        *image.get_pixel_mut(x as u32, y as u32) = Rgb(color);
    }
}

fn digit_glyph(character: char) -> [u8; 7] {
    match character {
        '0' => [
            0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110,
        ],
        '1' => [
            0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110,
        ],
        '2' => [
            0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111,
        ],
        '3' => [
            0b11110, 0b00001, 0b00001, 0b01110, 0b00001, 0b00001, 0b11110,
        ],
        '4' => [
            0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010,
        ],
        '5' => [
            0b11111, 0b10000, 0b10000, 0b11110, 0b00001, 0b00001, 0b11110,
        ],
        '6' => [
            0b01110, 0b10000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110,
        ],
        '7' => [
            0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000,
        ],
        '8' => [
            0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110,
        ],
        '9' => [
            0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00001, 0b01110,
        ],
        _ => [0; 7],
    }
}

fn element_style(element: &str) -> ElementStyle {
    let (color, radius, covalent_radius) = match element {
        "H" => ([238, 240, 244], 0.31, 0.31),
        "He" => ([217, 255, 255], 0.31, 0.28),
        "Li" => ([204, 128, 255], 0.48, 1.28),
        "Be" => ([194, 255, 0], 0.44, 0.96),
        "B" => ([255, 181, 181], 0.42, 0.84),
        "C" => ([58, 66, 76], 0.43, 0.76),
        "N" => ([52, 91, 214], 0.41, 0.71),
        "O" => ([222, 54, 64], 0.40, 0.66),
        "F" => ([76, 190, 94], 0.39, 0.57),
        "Ne" => ([179, 227, 245], 0.38, 0.58),
        "Na" => ([171, 92, 242], 0.52, 1.66),
        "Mg" => ([138, 255, 0], 0.50, 1.41),
        "Al" => ([191, 166, 166], 0.48, 1.21),
        "Si" => ([224, 164, 80], 0.47, 1.11),
        "P" => ([237, 117, 42], 0.46, 1.07),
        "S" => ([225, 190, 40], 0.45, 1.05),
        "Cl" => ([53, 181, 72], 0.45, 1.02),
        "Ar" => ([128, 209, 227], 0.44, 1.06),
        "K" => ([143, 64, 212], 0.57, 2.03),
        "Ca" => ([61, 255, 0], 0.55, 1.76),
        "Fe" => ([196, 93, 53], 0.50, 1.32),
        "Co" => ([222, 126, 148], 0.50, 1.26),
        "Ni" => ([72, 165, 88], 0.50, 1.24),
        "Cu" => ([184, 115, 51], 0.50, 1.32),
        "Zn" => ([110, 113, 176], 0.50, 1.22),
        "Br" => ([140, 42, 35], 0.48, 1.20),
        "Ag" => ([156, 160, 171], 0.52, 1.45),
        "I" => ([103, 57, 168], 0.52, 1.39),
        "Au" => ([214, 166, 52], 0.52, 1.36),
        "Hg" => ([143, 145, 166], 0.52, 1.32),
        "Pb" => ([91, 89, 112], 0.54, 1.46),
        _ => ([122, 132, 146], 0.46, 0.77),
    };
    ElementStyle {
        color,
        radius,
        covalent_radius,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::xyz::Trajectory;

    fn water() -> Frame {
        Trajectory::parse("3\nwater\nO 0 0 0\nH .758 .586 0\nH -.758 .586 0\n")
            .unwrap()
            .frames
            .remove(0)
    }

    fn settings() -> RenderSettings {
        RenderSettings {
            width: 320,
            height: 240,
            zoom: 1.0,
            rotation: Vec3::new(-20.0, 35.0, 0.0),
            atom_numbers: false,
            bonds: true,
        }
    }

    #[test]
    fn renders_non_blank_image_at_requested_size() {
        let image = render(&water(), settings()).unwrap();
        assert_eq!(image.dimensions(), (320, 240));
        assert!(image.pixels().any(|pixel| pixel.0 != BACKGROUND));
    }

    #[test]
    fn atom_numbers_change_rendered_pixels() {
        let plain = render(&water(), settings()).unwrap();
        let numbered = render(
            &water(),
            RenderSettings {
                atom_numbers: true,
                ..settings()
            },
        )
        .unwrap();
        assert_ne!(plain.as_raw(), numbered.as_raw());
    }

    #[test]
    fn zoom_changes_the_image() {
        let regular = render(&water(), settings()).unwrap();
        let zoomed = render(
            &water(),
            RenderSettings {
                zoom: 1.5,
                ..settings()
            },
        )
        .unwrap();
        assert_ne!(regular.as_raw(), zoomed.as_raw());
    }

    #[test]
    fn finds_water_bonds() {
        assert_eq!(infer_bonds(&water()).len(), 2);
    }
}
