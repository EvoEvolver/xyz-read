use std::{collections::BTreeMap, fs, path::Path};

use anyhow::{bail, Context, Result};
use serde::Serialize;

use crate::geometry::Vec3;

#[derive(Clone, Debug)]
pub struct Atom {
    pub element: String,
    pub position: Vec3,
}

#[derive(Clone, Debug)]
pub struct Frame {
    pub comment: String,
    pub atoms: Vec<Atom>,
}

#[derive(Clone, Debug)]
pub struct Trajectory {
    pub frames: Vec<Frame>,
}

#[derive(Debug, Serialize)]
pub struct Inspection {
    pub frame_count: usize,
    pub frames: Vec<FrameInspection>,
}

#[derive(Debug, Serialize)]
pub struct FrameInspection {
    pub frame: usize,
    pub comment: String,
    pub atom_count: usize,
    pub elements: BTreeMap<String, usize>,
    pub bounds: Bounds,
}

#[derive(Debug, Serialize)]
pub struct Bounds {
    pub min: [f32; 3],
    pub max: [f32; 3],
}

impl Trajectory {
    pub fn from_path(path: &Path) -> Result<Self> {
        let source = fs::read_to_string(path)
            .with_context(|| format!("could not read XYZ file {}", path.display()))?;
        Self::parse(&source).with_context(|| format!("invalid XYZ file {}", path.display()))
    }

    pub fn parse(source: &str) -> Result<Self> {
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
                atoms.push(Atom {
                    element,
                    position: Vec3::new(
                        parse_coordinate(fields[1], "X")?,
                        parse_coordinate(fields[2], "Y")?,
                        parse_coordinate(fields[3], "Z")?,
                    ),
                });
                cursor += 1;
            }

            frames.push(Frame { comment, atoms });
        }

        if frames.is_empty() {
            bail!("the file contains no XYZ frames");
        }
        Ok(Self { frames })
    }

    pub fn frame(&self, selector: &str) -> Result<(usize, &Frame)> {
        let index = if selector.eq_ignore_ascii_case("last") {
            self.frames.len() - 1
        } else {
            let human_index: usize = selector.parse().with_context(|| {
                format!("frame must be a 1-based number or 'last', found {selector:?}")
            })?;
            if human_index == 0 {
                bail!("frame numbers start at 1");
            }
            human_index - 1
        };
        let frame = self.frames.get(index).with_context(|| {
            format!(
                "frame {} does not exist; this file has {} frame(s)",
                index + 1,
                self.frames.len()
            )
        })?;
        Ok((index, frame))
    }

    pub fn inspect(&self, selected: Option<usize>) -> Inspection {
        let frames = self
            .frames
            .iter()
            .enumerate()
            .filter(|(index, _)| selected.is_none_or(|selected| selected == *index))
            .map(|(index, frame)| frame.inspect(index + 1))
            .collect();
        Inspection {
            frame_count: self.frames.len(),
            frames,
        }
    }
}

impl Frame {
    fn inspect(&self, frame_number: usize) -> FrameInspection {
        let mut elements = BTreeMap::new();
        for atom in &self.atoms {
            *elements.entry(atom.element.clone()).or_insert(0) += 1;
        }

        let (min, max) = if let Some(first) = self.atoms.first() {
            let mut min = first.position;
            let mut max = first.position;
            for atom in &self.atoms[1..] {
                min.x = min.x.min(atom.position.x);
                min.y = min.y.min(atom.position.y);
                min.z = min.z.min(atom.position.z);
                max.x = max.x.max(atom.position.x);
                max.y = max.y.max(atom.position.y);
                max.z = max.z.max(atom.position.z);
            }
            ([min.x, min.y, min.z], [max.x, max.y, max.z])
        } else {
            ([0.0; 3], [0.0; 3])
        };

        FrameInspection {
            frame: frame_number,
            comment: self.comment.clone(),
            atom_count: self.atoms.len(),
            elements,
            bounds: Bounds { min, max },
        }
    }
}

fn normalize_element(value: &str) -> Result<String> {
    if value.is_empty() || value.len() > 3 || !value.bytes().all(|byte| byte.is_ascii_alphabetic())
    {
        bail!("element symbols must contain one to three ASCII letters");
    }
    let mut chars = value.chars();
    let first = chars
        .next()
        .expect("value is not empty")
        .to_ascii_uppercase();
    let rest: String = chars
        .map(|character| character.to_ascii_lowercase())
        .collect();
    Ok(format!("{first}{rest}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    const TRAJECTORY: &str = "3\nwater 1\nO 0 0 0\nH 0.758 0.586 0 extra\nH -0.758 0.586 0\n3\nwater 2\nO 0 0 .1\nH .8 .6 0\nH -.8 .6 0\n";

    #[test]
    fn parses_multiple_frames_and_extra_columns() {
        let trajectory = Trajectory::parse(TRAJECTORY).unwrap();
        assert_eq!(trajectory.frames.len(), 2);
        assert_eq!(trajectory.frames[0].atoms.len(), 3);
        assert_eq!(trajectory.frames[0].atoms[0].element, "O");
        assert_eq!(trajectory.frames[1].comment, "water 2");
    }

    #[test]
    fn selects_human_numbered_and_last_frames() {
        let trajectory = Trajectory::parse(TRAJECTORY).unwrap();
        assert_eq!(trajectory.frame("1").unwrap().0, 0);
        assert_eq!(trajectory.frame("last").unwrap().0, 1);
        assert!(trajectory.frame("0").is_err());
        assert!(trajectory.frame("3").is_err());
    }

    #[test]
    fn rejects_truncated_frames() {
        let error = Trajectory::parse("2\ncomment\nH 0 0 0\n").unwrap_err();
        assert!(error.to_string().contains("declares 2 atoms"));
    }

    #[test]
    fn inspection_counts_elements() {
        let trajectory = Trajectory::parse(TRAJECTORY).unwrap();
        let inspection = trajectory.inspect(Some(0));
        assert_eq!(inspection.frame_count, 2);
        assert_eq!(inspection.frames[0].elements["H"], 2);
    }
}
