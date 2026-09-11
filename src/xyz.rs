use anyhow::{bail, Context, Result};

use crate::{
    geometry::Vec3,
    model::{normalize_element, Atom, Frame, Structure},
};

pub fn parse(source: &str) -> Result<Structure> {
    let lines: Vec<&str> = source.lines().collect();
    let mut cursor = 0;
    let mut frames = Vec::new();

    while cursor < lines.len() {
        while cursor < lines.len() && lines[cursor].trim().is_empty() {
            cursor += 1;
        }
        if cursor == lines.len() {
            break;
        }

        let count_line_number = cursor + 1;
        let atom_count: usize = lines[cursor].trim().parse().with_context(|| {
            format!(
                "line {count_line_number}: expected a non-negative atom count, found {:?}",
                lines[cursor]
            )
        })?;
        cursor += 1;

        if cursor >= lines.len() {
            bail!("frame {} is missing its comment line", frames.len() + 1);
        }
        let comment = lines[cursor].to_owned();
        cursor += 1;

        let mut atoms = Vec::with_capacity(atom_count);
        for atom_index in 0..atom_count {
            if cursor >= lines.len() {
                bail!(
                    "frame {} declares {atom_count} atoms but ends after {atom_index}",
                    frames.len() + 1
                );
            }
            let line_number = cursor + 1;
            let fields: Vec<&str> = lines[cursor].split_whitespace().collect();
            if fields.len() < 4 {
                bail!(
                    "line {line_number}: expected ELEMENT X Y Z, found {:?}",
                    lines[cursor]
                );
            }
            let element = normalize_element(fields[0]).with_context(|| {
                format!("line {line_number}: invalid element symbol {:?}", fields[0])
            })?;
            let parse_coordinate = |field: &str, axis: &str| -> Result<f32> {
                let value: f32 = field.parse().with_context(|| {
                    format!("line {line_number}: invalid {axis} coordinate {field:?}")
                })?;
                if !value.is_finite() {
                    bail!("line {line_number}: {axis} coordinate must be finite");
                }
                Ok(value)
            };
            atoms.push(Atom::new(
                element,
                Vec3::new(
                    parse_coordinate(fields[1], "X")?,
                    parse_coordinate(fields[2], "Y")?,
                    parse_coordinate(fields[3], "Z")?,
                ),
            ));
            cursor += 1;
        }

        frames.push(Frame {
            comment,
            atoms,
            bonds: Vec::new(),
            infer_bonds: true,
            cell: None,
        });
    }

    if frames.is_empty() {
        bail!("the file contains no XYZ frames");
    }
    Ok(Structure {
        format: "xyz".to_owned(),
        frames,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const TRAJECTORY: &str = "3\nwater 1\nO 0 0 0\nH 0.758 0.586 0 extra\nH -0.758 0.586 0\n3\nwater 2\nO 0 0 .1\nH .8 .6 0\nH -.8 .6 0\n";

    #[test]
    fn parses_multiple_frames_and_extra_columns() {
        let trajectory = parse(TRAJECTORY).unwrap();
        assert_eq!(trajectory.frames.len(), 2);
        assert_eq!(trajectory.frames[0].atoms.len(), 3);
        assert_eq!(trajectory.frames[0].atoms[0].element, "O");
        assert_eq!(trajectory.frames[1].comment, "water 2");
    }

    #[test]
    fn selects_human_numbered_and_last_frames() {
        let trajectory = parse(TRAJECTORY).unwrap();
        assert_eq!(trajectory.frame("1").unwrap().0, 0);
        assert_eq!(trajectory.frame("last").unwrap().0, 1);
        assert!(trajectory.frame("0").is_err());
        assert!(trajectory.frame("3").is_err());
    }

    #[test]
    fn rejects_truncated_frames() {
        let error = parse("2\ncomment\nH 0 0 0\n").unwrap_err();
        assert!(error.to_string().contains("declares 2 atoms"));
    }

    #[test]
    fn inspection_counts_elements() {
        let trajectory = parse(TRAJECTORY).unwrap();
        let inspection = trajectory.inspect(Some(0));
        assert_eq!(inspection.frame_count, 2);
        assert_eq!(inspection.frames[0].elements["H"], 2);
    }
}
