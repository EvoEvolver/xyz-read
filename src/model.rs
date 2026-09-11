use std::collections::{BTreeMap, BTreeSet};

use anyhow::{bail, Context, Result};
use serde::Serialize;

use crate::geometry::Vec3;

#[derive(Clone, Debug)]
pub struct Atom {
    pub element: String,
    pub position: Vec3,
    pub label: Option<String>,
    pub residue: Option<String>,
    pub chain: Option<String>,
    pub serial: Option<i32>,
}

impl Atom {
    pub fn new(element: impl Into<String>, position: Vec3) -> Self {
        Self {
            element: element.into(),
            position,
            label: None,
            residue: None,
            chain: None,
            serial: None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bond {
    pub first: usize,
    pub second: usize,
    pub order: u8,
}

impl Bond {
    pub fn new(first: usize, second: usize, order: u8) -> Self {
        Self {
            first: first.min(second),
            second: first.max(second),
            order: order.max(1),
        }
    }
}

#[derive(Clone, Copy, Debug, Serialize)]
pub struct UnitCell {
    pub lengths: [f32; 3],
    pub angles: [f32; 3],
    #[serde(skip)]
    pub vectors: [Vec3; 3],
}

impl UnitCell {
    pub fn from_parameters(lengths: [f32; 3], angles: [f32; 3]) -> Result<Self> {
        if lengths
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
        {
            bail!("unit-cell lengths must be finite and positive");
        }
        if angles
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0 || *value >= 180.0)
        {
            bail!("unit-cell angles must be between 0 and 180 degrees");
        }
        let [a, b, c] = lengths;
        let [alpha, beta, gamma] = angles.map(f32::to_radians);
        let sin_gamma = gamma.sin();
        if sin_gamma.abs() < 1e-6 {
            bail!("unit-cell gamma angle produces a singular lattice");
        }
        let c_x = c * beta.cos();
        let c_y = c * (alpha.cos() - beta.cos() * gamma.cos()) / sin_gamma;
        let c_z_squared = c * c - c_x * c_x - c_y * c_y;
        if c_z_squared < -1e-3 {
            bail!("unit-cell parameters do not form a valid lattice");
        }
        Ok(Self {
            lengths,
            angles,
            vectors: [
                Vec3::new(a, 0.0, 0.0),
                Vec3::new(b * gamma.cos(), b * gamma.sin(), 0.0),
                Vec3::new(c_x, c_y, c_z_squared.max(0.0).sqrt()),
            ],
        })
    }

    pub fn from_vectors(vectors: [Vec3; 3]) -> Result<Self> {
        let lengths = vectors.map(Vec3::length);
        let angle = |first: Vec3, second: Vec3| {
            let cosine = (first.dot(second) / (first.length() * second.length())).clamp(-1.0, 1.0);
            cosine.acos().to_degrees()
        };
        let angles = [
            angle(vectors[1], vectors[2]),
            angle(vectors[0], vectors[2]),
            angle(vectors[0], vectors[1]),
        ];
        let mut cell = Self::from_parameters(lengths, angles)?;
        cell.vectors = vectors;
        Ok(cell)
    }

    pub fn fractional_to_cartesian(&self, fractional: Vec3) -> Vec3 {
        self.vectors[0] * fractional.x
            + self.vectors[1] * fractional.y
            + self.vectors[2] * fractional.z
    }

    pub fn corners(&self) -> [Vec3; 8] {
        let [a, b, c] = self.vectors;
        [Vec3::default(), a, b, c, a + b, a + c, b + c, a + b + c]
    }
}

#[derive(Clone, Debug)]
pub struct Frame {
    pub comment: String,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub infer_bonds: bool,
    pub cell: Option<UnitCell>,
}

#[derive(Clone, Debug)]
pub struct Structure {
    pub format: String,
    pub frames: Vec<Frame>,
}

#[derive(Debug, Serialize)]
pub struct Inspection {
    pub format: String,
    pub frame_count: usize,
    pub frames: Vec<FrameInspection>,
}

#[derive(Debug, Serialize)]
pub struct FrameInspection {
    pub frame: usize,
    pub comment: String,
    pub atom_count: usize,
    pub explicit_bond_count: usize,
    pub connectivity: String,
    pub elements: BTreeMap<String, usize>,
    pub residue_count: usize,
    pub chains: Vec<String>,
    pub bounds: Bounds,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unit_cell: Option<UnitCell>,
}

#[derive(Debug, Serialize)]
pub struct Bounds {
    pub min: [f32; 3],
    pub max: [f32; 3],
}

impl Structure {
    pub fn frame(&self, selector: &str) -> Result<(usize, &Frame)> {
        let index = if selector.eq_ignore_ascii_case("last") {
            self.frames.len().saturating_sub(1)
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
            format: self.format.clone(),
            frame_count: self.frames.len(),
            frames,
        }
    }
}

impl Frame {
    fn inspect(&self, frame_number: usize) -> FrameInspection {
        let mut elements = BTreeMap::new();
        let mut residues = BTreeSet::new();
        let mut chains = BTreeSet::new();
        for atom in &self.atoms {
            *elements.entry(atom.element.clone()).or_insert(0) += 1;
            if let Some(residue) = &atom.residue {
                residues.insert((atom.chain.clone(), residue.clone()));
            }
            if let Some(chain) = &atom.chain {
                chains.insert(chain.clone());
            }
        }

        let (min, max) = if let Some(first) = self.atoms.first() {
            let mut min = first.position;
            let mut max = min;
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
            explicit_bond_count: self.bonds.len(),
            connectivity: match (self.bonds.is_empty(), self.infer_bonds) {
                (true, true) => "inferred",
                (false, true) => "explicit+inferred",
                (false, false) => "explicit",
                (true, false) => "none",
            }
            .to_owned(),
            elements,
            residue_count: residues.len(),
            chains: chains.into_iter().collect(),
            bounds: Bounds { min, max },
            unit_cell: self.cell,
        }
    }
}

pub fn normalize_element(value: &str) -> Result<String> {
    let value = value.trim();
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

pub fn deduplicate_bonds(bonds: &mut Vec<Bond>) {
    for bond in bonds.iter_mut() {
        *bond = Bond::new(bond.first, bond.second, bond.order);
    }
    bonds.sort_by_key(|bond| (bond.first, bond.second));
    bonds.dedup_by(|right, left| {
        if left.first == right.first && left.second == right.second {
            left.order = left.order.max(right.order);
            true
        } else {
            false
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triclinic_cell_converts_fractional_coordinates() {
        let cell = UnitCell::from_parameters([10.0, 11.0, 12.0], [90.0, 90.0, 120.0]).unwrap();
        let point = cell.fractional_to_cartesian(Vec3::new(1.0, 1.0, 1.0));
        assert!((point.x - 4.5).abs() < 0.001);
        assert!((point.y - 9.526).abs() < 0.001);
        assert!((point.z - 12.0).abs() < 0.001);
    }
}
