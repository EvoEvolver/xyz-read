use std::collections::{BTreeMap, BTreeSet, HashMap};

use anyhow::{bail, Context, Result};
use serde::Serialize;

use crate::model::{Atom, AtomRole, Bond, Frame};

#[derive(Debug, Serialize)]
pub struct QueryResult {
    pub schema_version: u32,
    pub tool_version: &'static str,
    pub input_sha256: String,
    pub format: String,
    pub frame: usize,
    pub selection: String,
    pub atom_count: usize,
    pub explicit_bond_count: usize,
    pub resolved_bond_count: usize,
    pub elements: BTreeMap<String, usize>,
    pub roles: BTreeMap<String, usize>,
    pub residue_count: usize,
    pub residues: Vec<ResidueRecord>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub atoms: Option<Vec<AtomRecord>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub explicit_bonds: Option<Vec<BondRecord>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub resolved_bonds: Option<Vec<BondRecord>>,
}

#[derive(Debug, Serialize)]
pub struct ResidueRecord {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain: Option<String>,
    pub residue: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub residue_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub residue_number: Option<i32>,
    pub role: String,
    pub atom_count: usize,
}

#[derive(Debug, Serialize)]
pub struct AtomRecord {
    /// Stable one-based position in the original input frame.
    pub index: usize,
    /// The same atom using RDKit's zero-based indexing convention.
    pub rdkit_index: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_serial: Option<i32>,
    pub element: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub atom_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub residue: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub residue_name: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub residue_number: Option<i32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chain: Option<String>,
    pub role: AtomRole,
    pub position: [f32; 3],
}

#[derive(Debug, Serialize)]
pub struct BondRecord {
    pub first: usize,
    pub second: usize,
    pub first_rdkit_index: usize,
    pub second_rdkit_index: usize,
    pub order: u8,
    pub explicit: bool,
}

pub fn query_result(
    format: &str,
    frame_number: usize,
    frame: &Frame,
    expression: &str,
    input_sha256: &str,
    resolved: &[Bond],
    detailed: bool,
) -> Result<QueryResult> {
    let selected = select(frame, expression)?;
    let mut elements = BTreeMap::new();
    let mut roles = BTreeMap::new();
    let mut residue_groups = BTreeMap::new();
    for &index in &selected {
        let atom = &frame.atoms[index];
        *elements.entry(atom.element.clone()).or_insert(0) += 1;
        *roles.entry(role_name(atom.role).to_owned()).or_insert(0) += 1;
        if let Some(residue) = &atom.residue {
            let key = (
                atom.chain.clone(),
                residue.clone(),
                atom.residue_name.clone(),
                atom.residue_number,
                role_name(atom.role).to_owned(),
            );
            *residue_groups.entry(key).or_insert(0) += 1;
        }
    }
    let residues = residue_groups
        .into_iter()
        .map(
            |((chain, residue, residue_name, residue_number, role), atom_count)| ResidueRecord {
                chain,
                residue,
                residue_name,
                residue_number,
                role,
                atom_count,
            },
        )
        .collect::<Vec<_>>();
    let explicit_bonds = frame
        .bonds
        .iter()
        .filter(|bond| selected.contains(&bond.first) && selected.contains(&bond.second))
        .map(|bond| bond_record(frame, *bond, true))
        .collect::<Vec<_>>();
    let resolved_bonds = resolved
        .iter()
        .filter(|bond| selected.contains(&bond.first) && selected.contains(&bond.second))
        .map(|bond| {
            let explicit = frame
                .bonds
                .iter()
                .any(|source| source.first == bond.first && source.second == bond.second);
            bond_record(frame, *bond, explicit)
        })
        .collect::<Vec<_>>();
    Ok(QueryResult {
        schema_version: 1,
        tool_version: env!("CARGO_PKG_VERSION"),
        input_sha256: input_sha256.to_owned(),
        format: format.to_owned(),
        frame: frame_number,
        selection: expression.to_owned(),
        atom_count: selected.len(),
        explicit_bond_count: explicit_bonds.len(),
        resolved_bond_count: resolved_bonds.len(),
        elements,
        roles,
        residue_count: residues.len(),
        residues,
        atoms: detailed.then(|| {
            selected
                .iter()
                .map(|&index| atom_record(&frame.atoms[index]))
                .collect()
        }),
        explicit_bonds: detailed.then_some(explicit_bonds),
        resolved_bonds: detailed.then_some(resolved_bonds),
    })
}

fn role_name(role: AtomRole) -> &'static str {
    match role {
        AtomRole::Unknown => "unknown",
        AtomRole::Polymer => "polymer",
        AtomRole::Ligand => "ligand",
        AtomRole::Water => "water",
        AtomRole::Ion => "ion",
    }
}

fn atom_record(atom: &Atom) -> AtomRecord {
    AtomRecord {
        index: atom.input_index + 1,
        rdkit_index: atom.input_index,
        source_serial: atom.serial,
        element: atom.element.clone(),
        atom_name: atom.label.clone(),
        residue: atom.residue.clone(),
        residue_name: atom.residue_name.clone(),
        residue_number: atom.residue_number,
        chain: atom.chain.clone(),
        role: atom.role,
        position: [atom.position.x, atom.position.y, atom.position.z],
    }
}

fn bond_record(frame: &Frame, bond: Bond, explicit: bool) -> BondRecord {
    let first = frame.atoms[bond.first].input_index;
    let second = frame.atoms[bond.second].input_index;
    BondRecord {
        first: first + 1,
        second: second + 1,
        first_rdkit_index: first,
        second_rdkit_index: second,
        order: bond.order,
        explicit,
    }
}

pub fn select(frame: &Frame, expression: &str) -> Result<BTreeSet<usize>> {
    let expression = expression.trim();
    if expression.is_empty() {
        bail!("selection expression cannot be empty");
    }
    evaluate(frame, expression)
        .with_context(|| format!("invalid selection expression {expression:?}"))
}

fn evaluate(frame: &Frame, expression: &str) -> Result<BTreeSet<usize>> {
    let expression = strip_outer_parentheses(expression.trim())?;
    if let Some((left, right)) = split_top_level(expression, '|')? {
        let mut result = evaluate(frame, left)?;
        result.extend(evaluate(frame, right)?);
        return Ok(result);
    }
    if let Some((left, right)) = split_top_level(expression, '&')? {
        let left = evaluate(frame, left)?;
        let right = evaluate(frame, right)?;
        return Ok(left.intersection(&right).copied().collect());
    }
    if let Some(rest) = expression.strip_prefix('!') {
        let excluded = evaluate(frame, rest)?;
        return Ok((0..frame.atoms.len())
            .filter(|index| !excluded.contains(index))
            .collect());
    }
    evaluate_term(frame, expression)
}

fn evaluate_term(frame: &Frame, term: &str) -> Result<BTreeSet<usize>> {
    let lower = term.to_ascii_lowercase();
    let role = match lower.as_str() {
        "all" => return Ok((0..frame.atoms.len()).collect()),
        "protein" | "polymer" => Some(AtomRole::Polymer),
        "ligand" => Some(AtomRole::Ligand),
        "water" => Some(AtomRole::Water),
        "ion" => Some(AtomRole::Ion),
        _ => None,
    };
    if let Some(role) = role {
        return Ok(frame
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(index, atom)| (atom.role == role).then_some(index))
            .collect());
    }

    if lower.starts_with("byres(") && term.ends_with(')') {
        let selected = evaluate(frame, &term[6..term.len() - 1])?;
        return Ok(expand_residues(frame, &selected));
    }

    let distance_mode = if lower.starts_with("within:") {
        Some((7, false))
    } else if lower.starts_with("around:") {
        Some((7, true))
    } else {
        None
    };
    if let Some((prefix_length, expand)) = distance_mode {
        let (distance, target) = term[prefix_length..].split_once('@').context(
            "distance selection uses within:DISTANCE@SELECTION or around:DISTANCE@SELECTION",
        )?;
        let distance: f32 = distance
            .parse()
            .context("within distance is not a number")?;
        if !distance.is_finite() || distance <= 0.0 {
            bail!("within distance must be finite and greater than zero");
        }
        let targets = evaluate(frame, target)?;
        if targets.is_empty() {
            return Ok(BTreeSet::new());
        }
        let distance_squared = distance * distance;
        let selected = frame
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(index, atom)| {
                targets
                    .iter()
                    .any(|target| {
                        let delta = atom.position - frame.atoms[*target].position;
                        delta.dot(delta) <= distance_squared
                    })
                    .then_some(index)
            })
            .collect::<BTreeSet<_>>();
        return Ok(if expand {
            expand_residues(frame, &selected)
        } else {
            selected
        });
    }

    let (key, values) = term
        .split_once(':')
        .context("expected a selector such as element:C, chain:A, residue:42, or index:1")?;
    let key = key.trim().to_ascii_lowercase();
    let values = values.trim();
    if values.is_empty() {
        bail!("{key} selector requires a value");
    }
    match key.as_str() {
        "index" => select_indices(frame, values, false),
        "rdkit-index" | "rdkit" => select_indices(frame, values, true),
        "element" => select_text(frame, values, |atom| Some(atom.element.as_str())),
        "chain" => select_text(frame, values, |atom| atom.chain.as_deref()),
        "atom" | "atom-name" => select_text(frame, values, |atom| atom.label.as_deref()),
        "residue" => select_residues(frame, values),
        _ => bail!("unknown selector {key:?}"),
    }
}

fn expand_residues(frame: &Frame, selected: &BTreeSet<usize>) -> BTreeSet<usize> {
    let residue_keys = selected
        .iter()
        .filter_map(|index| {
            let atom = &frame.atoms[*index];
            Some((atom.chain.clone(), atom.residue.clone()?))
        })
        .collect::<BTreeSet<_>>();
    frame
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| {
            if selected.contains(&index)
                || atom.residue.as_ref().is_some_and(|residue| {
                    residue_keys.contains(&(atom.chain.clone(), residue.clone()))
                })
            {
                Some(index)
            } else {
                None
            }
        })
        .collect()
}

fn select_indices(frame: &Frame, values: &str, zero_based: bool) -> Result<BTreeSet<usize>> {
    let requested = parse_integer_set(values)?;
    let mut selected = BTreeSet::new();
    for value in requested {
        let index = if zero_based {
            if value < 0 {
                bail!("RDKit atom indices cannot be negative");
            }
            value as usize
        } else {
            if value <= 0 {
                bail!("atom indices are one-based and must be positive");
            }
            value as usize - 1
        };
        if index >= frame.atoms.len() {
            bail!(
                "atom index {} is outside this frame's {} atoms",
                if zero_based { index } else { index + 1 },
                frame.atoms.len()
            );
        }
        selected.insert(index);
    }
    Ok(selected)
}

fn select_text<'a>(
    frame: &'a Frame,
    values: &str,
    accessor: impl Fn(&'a Atom) -> Option<&'a str>,
) -> Result<BTreeSet<usize>> {
    let requested = values
        .split(',')
        .map(|value| value.trim().to_ascii_lowercase())
        .filter(|value| !value.is_empty())
        .collect::<BTreeSet<_>>();
    if requested.is_empty() {
        bail!("selector value list cannot be empty");
    }
    Ok(frame
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| {
            accessor(atom)
                .is_some_and(|value| requested.contains(&value.to_ascii_lowercase()))
                .then_some(index)
        })
        .collect())
}

fn select_residues(frame: &Frame, values: &str) -> Result<BTreeSet<usize>> {
    let mut numbers = BTreeSet::new();
    let mut names = BTreeSet::new();
    for value in values
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        if value
            .chars()
            .all(|character| character.is_ascii_digit() || matches!(character, '-' | '.'))
        {
            numbers.extend(parse_integer_set(value)?);
        } else {
            names.insert(value.to_ascii_lowercase());
        }
    }
    Ok(frame
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| {
            let number_match = atom
                .residue_number
                .is_some_and(|number| numbers.contains(&number));
            let name_match = atom.residue_name.as_ref().is_some_and(|name| {
                names.contains(&name.to_ascii_lowercase())
                    || atom
                        .residue
                        .as_ref()
                        .is_some_and(|label| names.contains(&label.to_ascii_lowercase()))
            });
            (number_match || name_match).then_some(index)
        })
        .collect())
}

fn parse_integer_set(values: &str) -> Result<BTreeSet<i32>> {
    let mut result = BTreeSet::new();
    for value in values
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        if let Some((start, end)) = value.split_once("..") {
            let start: i32 = start.parse().context("range start is not an integer")?;
            let end: i32 = end.parse().context("range end is not an integer")?;
            if start > end {
                bail!("range start must not be greater than range end");
            }
            if end.saturating_sub(start) > 1_000_000 {
                bail!("selection range is too large");
            }
            result.extend(start..=end);
        } else {
            result.insert(value.parse().context("selector value is not an integer")?);
        }
    }
    if result.is_empty() {
        bail!("selector value list cannot be empty");
    }
    Ok(result)
}

fn strip_outer_parentheses(mut expression: &str) -> Result<&str> {
    loop {
        if !expression.starts_with('(') || !expression.ends_with(')') {
            return Ok(expression);
        }
        let mut depth = 0_i32;
        let mut encloses_all = true;
        for (index, character) in expression.char_indices() {
            match character {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth < 0 {
                        bail!("unmatched closing parenthesis");
                    }
                    if depth == 0 && index + 1 != expression.len() {
                        encloses_all = false;
                        break;
                    }
                }
                _ => {}
            }
        }
        if depth != 0 {
            bail!("unmatched opening parenthesis");
        }
        if !encloses_all {
            return Ok(expression);
        }
        expression = expression[1..expression.len() - 1].trim();
    }
}

fn split_top_level(expression: &str, operator: char) -> Result<Option<(&str, &str)>> {
    let mut depth = 0_i32;
    for (index, character) in expression.char_indices().rev() {
        match character {
            ')' => depth += 1,
            '(' => {
                depth -= 1;
                if depth < 0 {
                    bail!("unmatched opening parenthesis");
                }
            }
            value if value == operator && depth == 0 => {
                let left = expression[..index].trim();
                let right = expression[index + operator.len_utf8()..].trim();
                if left.is_empty() || right.is_empty() {
                    bail!("operator {operator} requires expressions on both sides");
                }
                return Ok(Some((left, right)));
            }
            _ => {}
        }
    }
    if depth != 0 {
        bail!("unmatched closing parenthesis");
    }
    Ok(None)
}

pub fn subset(frame: &Frame, selected: &BTreeSet<usize>) -> Result<Frame> {
    if selected.is_empty() {
        bail!("selection matched no atoms");
    }
    let mut remap = HashMap::new();
    let atoms = selected
        .iter()
        .enumerate()
        .map(|(new_index, old_index)| {
            remap.insert(*old_index, new_index);
            frame.atoms[*old_index].clone()
        })
        .collect();
    let bonds = frame
        .bonds
        .iter()
        .filter_map(|bond| {
            Some(Bond::new(
                *remap.get(&bond.first)?,
                *remap.get(&bond.second)?,
                bond.order,
            ))
        })
        .collect();
    Ok(Frame {
        comment: frame.comment.clone(),
        atoms,
        bonds,
        infer_bonds: frame.infer_bonds,
        cell: frame.cell,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{geometry::Vec3, model::Atom};

    fn frame() -> Frame {
        let mut first = Atom::new("C", Vec3::new(0.0, 0.0, 0.0));
        first.input_index = 0;
        first.chain = Some("A".to_owned());
        first.residue_name = Some("ALA".to_owned());
        first.residue_number = Some(10);
        first.label = Some("CA".to_owned());
        first.role = AtomRole::Polymer;
        let mut second = Atom::new("O", Vec3::new(2.0, 0.0, 0.0));
        second.input_index = 1;
        second.role = AtomRole::Ligand;
        Frame {
            comment: String::new(),
            atoms: vec![first, second],
            bonds: vec![Bond::new(0, 1, 1)],
            infer_bonds: false,
            cell: None,
        }
    }

    #[test]
    fn combines_semantic_and_property_selectors() {
        let selected = select(&frame(), "protein & chain:A & residue:10 & atom:CA").unwrap();
        assert_eq!(selected, BTreeSet::from([0]));
    }

    #[test]
    fn supports_rdkit_indices_and_distance_queries() {
        assert_eq!(
            select(&frame(), "rdkit-index:0").unwrap(),
            BTreeSet::from([0])
        );
        assert_eq!(
            select(&frame(), "within:2.1@ligand").unwrap(),
            BTreeSet::from([0, 1])
        );
    }

    #[test]
    fn around_expands_nearby_atoms_to_complete_residues() {
        let mut frame = frame();
        let mut third = Atom::new("H", Vec3::new(4.0, 0.0, 0.0));
        third.input_index = 2;
        third.chain = Some("A".to_owned());
        third.residue = Some("ALA10".to_owned());
        third.residue_name = Some("ALA".to_owned());
        third.residue_number = Some(10);
        third.role = AtomRole::Polymer;
        frame.atoms[0].residue = Some("ALA10".to_owned());
        frame.atoms.push(third);
        assert_eq!(
            select(&frame, "around:2.1@ligand").unwrap(),
            BTreeSet::from([0, 1, 2])
        );
        assert_eq!(
            select(&frame, "byres(within:2.1@ligand)").unwrap(),
            BTreeSet::from([0, 1, 2])
        );
    }

    #[test]
    fn subset_preserves_input_indices_and_remaps_bonds() {
        let selected = BTreeSet::from([0, 1]);
        let subset = subset(&frame(), &selected).unwrap();
        assert_eq!(subset.atoms[1].input_index, 1);
        assert_eq!(subset.bonds, vec![Bond::new(0, 1, 1)]);
    }
}
