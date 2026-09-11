use std::{collections::HashMap, path::Path};

use anyhow::{bail, Context, Result};
use clap::ValueEnum;

use crate::{
    cif,
    geometry::Vec3,
    model::{
        deduplicate_bonds, normalize_element, Atom, AtomRole, Bond, Frame, Structure, UnitCell,
    },
    xyz,
};

#[derive(Clone, Copy, Debug, Default, ValueEnum)]
pub enum InputFormat {
    #[default]
    Auto,
    Xyz,
    Pdb,
    Cif,
    Mmcif,
    Mol,
    Sdf,
    Mol2,
    Poscar,
}

pub fn load_text(
    source: &str,
    name_hint: Option<&Path>,
    requested: InputFormat,
) -> Result<Structure> {
    let fallback = Path::new("stdin");
    let path = name_hint.unwrap_or(fallback);
    let format = match requested {
        InputFormat::Auto => detect(path, source)?,
        format => format,
    };
    let result = match format {
        InputFormat::Auto => unreachable!(),
        InputFormat::Xyz => xyz::parse(source),
        InputFormat::Pdb => parse_pdb(source),
        InputFormat::Cif | InputFormat::Mmcif => cif::parse(source),
        InputFormat::Mol | InputFormat::Sdf => parse_sdf(source),
        InputFormat::Mol2 => parse_mol2(source),
        InputFormat::Poscar => parse_poscar(source),
    };
    let mut structure = result.with_context(|| {
        format!(
            "invalid {} structure from {}",
            format_name(format),
            path.display()
        )
    })?;
    structure.assign_input_indices();
    Ok(structure)
}

fn detect(path: &Path, source: &str) -> Result<InputFormat> {
    let extension = path
        .extension()
        .and_then(|extension| extension.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase();
    let format = match extension.as_str() {
        "xyz" => Some(InputFormat::Xyz),
        "pdb" | "ent" => Some(InputFormat::Pdb),
        "cif" | "mmcif" | "mcif" => Some(InputFormat::Cif),
        "mol" => Some(InputFormat::Mol),
        "sdf" | "sd" => Some(InputFormat::Sdf),
        "mol2" => Some(InputFormat::Mol2),
        "vasp" | "poscar" => Some(InputFormat::Poscar),
        _ if matches!(file_name.as_str(), "poscar" | "contcar") => Some(InputFormat::Poscar),
        _ => None,
    };
    if let Some(format) = format {
        return Ok(format);
    }

    let first = source
        .lines()
        .find(|line| !line.trim().is_empty())
        .unwrap_or("");
    if source.contains("@<TRIPOS>MOLECULE") {
        Ok(InputFormat::Mol2)
    } else if source.contains("V2000") || source.contains("V3000") {
        Ok(InputFormat::Mol)
    } else if first.starts_with("data_") || source.contains("_atom_site_") {
        Ok(InputFormat::Cif)
    } else if source.lines().any(|line| {
        matches!(
            line.get(..6).unwrap_or("").trim(),
            "ATOM" | "HETATM" | "HEADER" | "MODEL"
        )
    }) {
        Ok(InputFormat::Pdb)
    } else if first.trim().parse::<usize>().is_ok() {
        Ok(InputFormat::Xyz)
    } else if looks_like_poscar(source) {
        Ok(InputFormat::Poscar)
    } else {
        bail!("could not detect the input format; use --format xyz|pdb|cif|mol|sdf|mol2|poscar")
    }
}

fn looks_like_poscar(source: &str) -> bool {
    let lines: Vec<&str> = source
        .lines()
        .filter(|line| !line.trim().is_empty())
        .collect();
    lines.len() >= 8
        && lines[1].trim().parse::<f32>().is_ok()
        && lines[2..5].iter().all(|line| {
            line.split_whitespace()
                .take(3)
                .filter(|value| value.parse::<f32>().is_ok())
                .count()
                == 3
        })
}

fn format_name(format: InputFormat) -> &'static str {
    match format {
        InputFormat::Auto => "auto",
        InputFormat::Xyz => "XYZ",
        InputFormat::Pdb => "PDB",
        InputFormat::Cif => "CIF/mmCIF",
        InputFormat::Mmcif => "mmCIF",
        InputFormat::Mol => "MOL",
        InputFormat::Sdf => "SDF",
        InputFormat::Mol2 => "MOL2",
        InputFormat::Poscar => "POSCAR",
    }
}

fn parse_pdb(source: &str) -> Result<Structure> {
    let mut title = String::new();
    let mut cell = None;
    let mut frames = Vec::new();
    let mut atoms = Vec::new();
    let mut connections = Vec::new();
    let mut model_number = None;

    let finish_frame = |frames: &mut Vec<Frame>,
                        atoms: &mut Vec<Atom>,
                        model: Option<String>,
                        title: &str,
                        cell: Option<UnitCell>| {
        if atoms.is_empty() {
            return;
        }
        frames.push(Frame {
            comment: model
                .map(|number| format!("model {number}"))
                .unwrap_or_else(|| title.trim().to_owned()),
            atoms: std::mem::take(atoms),
            bonds: Vec::new(),
            infer_bonds: true,
            cell,
        });
    };

    for (line_index, line) in source.lines().enumerate() {
        let record = field(line, 0, 6).trim();
        match record {
            "HEADER" | "TITLE" if title.is_empty() || record == "TITLE" => {
                if !title.is_empty() {
                    title.push(' ');
                }
                title.push_str(field(line, 10, 80).trim());
            }
            "CRYST1" => {
                let lengths = [
                    parse_pdb_float(line, 6, 15, line_index + 1, "cell a")?,
                    parse_pdb_float(line, 15, 24, line_index + 1, "cell b")?,
                    parse_pdb_float(line, 24, 33, line_index + 1, "cell c")?,
                ];
                let angles = [
                    parse_pdb_float(line, 33, 40, line_index + 1, "cell alpha")?,
                    parse_pdb_float(line, 40, 47, line_index + 1, "cell beta")?,
                    parse_pdb_float(line, 47, 54, line_index + 1, "cell gamma")?,
                ];
                cell = if lengths.iter().all(|length| (*length - 1.0).abs() < 0.002)
                    && angles.iter().all(|angle| (*angle - 90.0).abs() < 0.01)
                {
                    None
                } else {
                    Some(UnitCell::from_parameters(lengths, angles)?)
                };
            }
            "MODEL" => {
                finish_frame(&mut frames, &mut atoms, model_number.take(), &title, cell);
                model_number = Some(field(line, 10, 20).trim().to_owned());
            }
            "ENDMDL" => finish_frame(&mut frames, &mut atoms, model_number.take(), &title, cell),
            "ATOM" | "HETATM" => {
                let alternate = field(line, 16, 17).trim();
                if !alternate.is_empty() && alternate != "A" && alternate != "1" {
                    continue;
                }
                let atom_name = field(line, 12, 16).trim();
                let element_field = field(line, 76, 78).trim();
                let element = if element_field.is_empty() {
                    pdb_element(atom_name, record == "ATOM")?
                } else {
                    normalize_element(element_field)?
                };
                let serial = field(line, 6, 11).trim().parse::<i32>().ok();
                let mut atom = Atom::new(
                    element,
                    Vec3::new(
                        parse_pdb_float(line, 30, 38, line_index + 1, "X")?,
                        parse_pdb_float(line, 38, 46, line_index + 1, "Y")?,
                        parse_pdb_float(line, 46, 54, line_index + 1, "Z")?,
                    ),
                );
                atom.label = Some(atom_name.to_owned());
                atom.serial = serial;
                let residue_name = field(line, 17, 20).trim();
                let residue_number_text = field(line, 22, 26).trim();
                let insertion_code = field(line, 26, 27).trim();
                if !residue_name.is_empty() || !residue_number_text.is_empty() {
                    atom.residue = Some(format!(
                        "{residue_name}{residue_number_text}{insertion_code}"
                    ));
                }
                atom.residue_name = (!residue_name.is_empty()).then(|| residue_name.to_owned());
                atom.residue_number = residue_number_text.parse().ok();
                atom.role = pdb_atom_role(record, residue_name);
                let chain = field(line, 21, 22).trim();
                if !chain.is_empty() {
                    atom.chain = Some(chain.to_owned());
                }
                atoms.push(atom);
            }
            "CONECT" => {
                let serials: Vec<i32> = line
                    .get(6..)
                    .unwrap_or("")
                    .split_whitespace()
                    .filter_map(|value| value.parse().ok())
                    .collect();
                if let Some((&first, others)) = serials.split_first() {
                    connections.extend(others.iter().map(|&second| (first, second)));
                }
            }
            _ => {}
        }
    }
    finish_frame(&mut frames, &mut atoms, model_number.take(), &title, cell);
    if frames.is_empty() {
        bail!("the file contains no ATOM or HETATM records");
    }

    for frame in &mut frames {
        let indices: HashMap<i32, usize> = frame
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(index, atom)| atom.serial.map(|serial| (serial, index)))
            .collect();
        for &(first, second) in &connections {
            if let (Some(&first), Some(&second)) = (indices.get(&first), indices.get(&second)) {
                frame.bonds.push(Bond::new(first, second, 1));
            }
        }
        deduplicate_bonds(&mut frame.bonds);
    }
    Ok(Structure {
        format: "pdb".to_owned(),
        frames,
    })
}

fn pdb_atom_role(record: &str, residue_name: &str) -> AtomRole {
    if record == "ATOM" {
        return AtomRole::Polymer;
    }
    hetero_atom_role(residue_name)
}

pub(crate) fn hetero_atom_role(residue_name: &str) -> AtomRole {
    let residue = residue_name.to_ascii_uppercase();
    if matches!(residue.as_str(), "HOH" | "WAT" | "DOD" | "H2O") {
        AtomRole::Water
    } else if matches!(
        residue.as_str(),
        "LI" | "NA"
            | "K"
            | "RB"
            | "CS"
            | "MG"
            | "CA"
            | "SR"
            | "BA"
            | "ZN"
            | "FE"
            | "CU"
            | "MN"
            | "CO"
            | "NI"
            | "CD"
            | "HG"
            | "CL"
            | "BR"
            | "IOD"
    ) {
        AtomRole::Ion
    } else {
        AtomRole::Ligand
    }
}

fn field(line: &str, start: usize, end: usize) -> &str {
    line.get(start.min(line.len())..end.min(line.len()))
        .unwrap_or("")
}

fn parse_pdb_float(line: &str, start: usize, end: usize, number: usize, name: &str) -> Result<f32> {
    let value = field(line, start, end).trim();
    value
        .parse()
        .with_context(|| format!("line {number}: invalid {name} value {value:?}"))
}

fn pdb_element(atom_name: &str, polymer_atom: bool) -> Result<String> {
    let letters: String = atom_name
        .chars()
        .filter(|character| character.is_ascii_alphabetic())
        .collect();
    if letters.is_empty() {
        bail!("cannot derive an element from atom name {atom_name:?}");
    }
    let candidate = if polymer_atom || letters.len() == 1 {
        &letters[..1]
    } else {
        &letters[..letters.len().min(2)]
    };
    normalize_element(candidate)
}

fn parse_sdf(source: &str) -> Result<Structure> {
    let mut frames = Vec::new();
    for (record_index, record) in source.split("$$$$").enumerate() {
        let record = record.trim_end_matches(['\r', '\n']);
        // The newline after an SDF delimiter separates records. A leading newline in
        // the first record, however, is a valid empty MOL title as emitted by RDKit.
        let record = if record_index > 0 {
            record
                .strip_prefix("\r\n")
                .or_else(|| record.strip_prefix('\n'))
                .unwrap_or(record)
        } else {
            record
        };
        if record.trim().is_empty() {
            continue;
        }
        frames.push(parse_mol_record(record)?);
    }
    if frames.is_empty() {
        bail!("the file contains no molecule records");
    }
    Ok(Structure {
        format: if source.contains("$$$$") {
            "sdf"
        } else {
            "mol"
        }
        .to_owned(),
        frames,
    })
}

fn parse_mol_record(record: &str) -> Result<Frame> {
    let lines: Vec<&str> = record.lines().collect();
    if lines.len() < 4 {
        bail!("MOL record is shorter than its four-line header");
    }
    if lines[3].contains("V3000") {
        return parse_v3000(&lines);
    }
    let atom_count = field(lines[3], 0, 3).trim().parse::<usize>().or_else(|_| {
        lines[3]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .parse::<usize>()
    })?;
    let bond_count = field(lines[3], 3, 6).trim().parse::<usize>().or_else(|_| {
        lines[3]
            .split_whitespace()
            .nth(1)
            .unwrap_or("")
            .parse::<usize>()
    })?;
    if lines.len() < 4 + atom_count + bond_count {
        bail!("MOL counts line declares more atoms or bonds than the record contains");
    }
    let mut atoms = Vec::with_capacity(atom_count);
    for (offset, line) in lines[4..4 + atom_count].iter().enumerate() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 4 {
            bail!("MOL atom line {} is incomplete", offset + 5);
        }
        let position = Vec3::new(fields[0].parse()?, fields[1].parse()?, fields[2].parse()?);
        let mut atom = Atom::new(normalize_element(fields[3])?, position);
        atom.role = AtomRole::Ligand;
        atoms.push(atom);
    }
    let mut bonds = Vec::with_capacity(bond_count);
    for (offset, line) in lines[4 + atom_count..4 + atom_count + bond_count]
        .iter()
        .enumerate()
    {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            bail!("MOL bond line {} is incomplete", offset + atom_count + 5);
        }
        let first = fields[0]
            .parse::<usize>()?
            .checked_sub(1)
            .context("MOL atom indices start at 1")?;
        let second = fields[1]
            .parse::<usize>()?
            .checked_sub(1)
            .context("MOL atom indices start at 1")?;
        if first >= atoms.len() || second >= atoms.len() {
            bail!("MOL bond references an atom outside the atom block");
        }
        bonds.push(Bond::new(
            first,
            second,
            fields[2].parse::<u8>().unwrap_or(1),
        ));
    }
    deduplicate_bonds(&mut bonds);
    Ok(Frame {
        comment: lines[0].trim().to_owned(),
        atoms,
        bonds,
        infer_bonds: false,
        cell: None,
    })
}

fn parse_v3000(lines: &[&str]) -> Result<Frame> {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut atom_ids = HashMap::new();
    let mut section = "";
    for line in lines {
        let content = line.strip_prefix("M  V30 ").unwrap_or(line).trim();
        match content {
            "BEGIN ATOM" => section = "ATOM",
            "END ATOM" => section = "",
            "BEGIN BOND" => section = "BOND",
            "END BOND" => section = "",
            _ if section == "ATOM" => {
                let fields: Vec<&str> = content.split_whitespace().collect();
                if fields.len() >= 6 {
                    let id: i32 = fields[0].parse()?;
                    atom_ids.insert(id, atoms.len());
                    let mut atom = Atom::new(
                        normalize_element(fields[1])?,
                        Vec3::new(fields[2].parse()?, fields[3].parse()?, fields[4].parse()?),
                    );
                    atom.role = AtomRole::Ligand;
                    atoms.push(atom);
                }
            }
            _ if section == "BOND" => {
                let fields: Vec<&str> = content.split_whitespace().collect();
                if fields.len() >= 4 {
                    let order = fields[1].parse::<u8>().unwrap_or(1);
                    let first = atom_ids
                        .get(&fields[2].parse::<i32>()?)
                        .context("V3000 bond references an unknown atom")?;
                    let second = atom_ids
                        .get(&fields[3].parse::<i32>()?)
                        .context("V3000 bond references an unknown atom")?;
                    bonds.push(Bond::new(*first, *second, order));
                }
            }
            _ => {}
        }
    }
    if atoms.is_empty() {
        bail!("V3000 record contains no atoms");
    }
    deduplicate_bonds(&mut bonds);
    Ok(Frame {
        comment: lines.first().unwrap_or(&"").trim().to_owned(),
        atoms,
        bonds,
        infer_bonds: false,
        cell: None,
    })
}

fn parse_mol2(source: &str) -> Result<Structure> {
    let mut frames = Vec::new();
    for block in source.split("@<TRIPOS>MOLECULE").skip(1) {
        let mut name = String::new();
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut atom_ids = HashMap::new();
        let mut section = "MOLECULE";
        for line in block.lines() {
            let trimmed = line.trim();
            if let Some(next) = trimmed.strip_prefix("@<TRIPOS>") {
                section = next;
                continue;
            }
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            match section {
                "MOLECULE" if name.is_empty() => name = trimmed.to_owned(),
                "ATOM" => {
                    let fields: Vec<&str> = trimmed.split_whitespace().collect();
                    if fields.len() < 6 {
                        bail!("incomplete MOL2 atom line {trimmed:?}");
                    }
                    let id: i32 = fields[0].parse()?;
                    let type_element = fields[5].split('.').next().unwrap_or(fields[5]);
                    let mut atom = Atom::new(
                        normalize_element(type_element)?,
                        Vec3::new(fields[2].parse()?, fields[3].parse()?, fields[4].parse()?),
                    );
                    atom.label = Some(fields[1].to_owned());
                    atom.role = AtomRole::Ligand;
                    if fields.len() >= 8 {
                        atom.residue = Some(fields[7].to_owned());
                        let (name, number) = split_residue_label(fields[7]);
                        atom.residue_name = name;
                        atom.residue_number = number;
                    }
                    atom_ids.insert(id, atoms.len());
                    atoms.push(atom);
                }
                "BOND" => {
                    let fields: Vec<&str> = trimmed.split_whitespace().collect();
                    if fields.len() < 4 {
                        bail!("incomplete MOL2 bond line {trimmed:?}");
                    }
                    let first = atom_ids
                        .get(&fields[1].parse::<i32>()?)
                        .context("MOL2 bond references an unknown atom")?;
                    let second = atom_ids
                        .get(&fields[2].parse::<i32>()?)
                        .context("MOL2 bond references an unknown atom")?;
                    let order = match fields[3] {
                        "2" => 2,
                        "3" => 3,
                        _ => 1,
                    };
                    bonds.push(Bond::new(*first, *second, order));
                }
                _ => {}
            }
        }
        if !atoms.is_empty() {
            deduplicate_bonds(&mut bonds);
            frames.push(Frame {
                comment: name,
                atoms,
                bonds,
                infer_bonds: false,
                cell: None,
            });
        }
    }
    if frames.is_empty() {
        bail!("the file contains no MOL2 molecule blocks");
    }
    Ok(Structure {
        format: "mol2".to_owned(),
        frames,
    })
}

fn split_residue_label(value: &str) -> (Option<String>, Option<i32>) {
    let boundary = value
        .char_indices()
        .find(|(_, character)| character.is_ascii_digit() || *character == '-')
        .map(|(index, _)| index)
        .unwrap_or(value.len());
    let (name, number) = value.split_at(boundary);
    (
        (!name.is_empty()).then(|| name.to_owned()),
        number.parse::<i32>().ok(),
    )
}

fn parse_poscar(source: &str) -> Result<Structure> {
    let lines: Vec<&str> = source
        .lines()
        .filter(|line| !line.trim().is_empty())
        .collect();
    if lines.len() < 8 {
        bail!("POSCAR requires a title, scale, lattice, elements, counts, mode, and positions");
    }
    let title = lines[0].trim().to_owned();
    let requested_scale: f32 = lines[1].trim().parse()?;
    if !requested_scale.is_finite() || requested_scale == 0.0 {
        bail!("POSCAR scale must be finite and non-zero");
    }
    let parse_vector = |line: &str| -> Result<Vec3> {
        let fields: Vec<f32> = line
            .split_whitespace()
            .map(str::parse)
            .collect::<std::result::Result<_, _>>()?;
        if fields.len() < 3 {
            bail!("POSCAR lattice vector requires three values");
        }
        Ok(Vec3::new(fields[0], fields[1], fields[2]))
    };
    let raw_vectors = [
        parse_vector(lines[2])?,
        parse_vector(lines[3])?,
        parse_vector(lines[4])?,
    ];
    let raw_volume = determinant(raw_vectors).abs();
    if raw_volume < 1e-8 {
        bail!("POSCAR lattice vectors form a singular cell");
    }
    let scale = if requested_scale > 0.0 {
        requested_scale
    } else {
        ((-requested_scale) / raw_volume).cbrt()
    };
    let cell = UnitCell::from_vectors(raw_vectors.map(|vector| vector * scale))?;
    let fifth_fields: Vec<&str> = lines[5].split_whitespace().collect();
    let vasp4 = fifth_fields
        .iter()
        .all(|value| value.parse::<usize>().is_ok());
    let (elements, counts_line, mut cursor): (Vec<String>, usize, usize) = if vasp4 {
        let counts: Vec<usize> = fifth_fields
            .iter()
            .map(|value| value.parse())
            .collect::<std::result::Result<_, _>>()?;
        let title_elements: Vec<String> = title
            .split_whitespace()
            .filter_map(|value| normalize_element(value).ok())
            .collect();
        if title_elements.len() != counts.len() {
            bail!(
                "VASP 4 POSCAR has no element line; put {} element symbol(s) in the title",
                counts.len()
            );
        }
        (title_elements, 5, 6)
    } else {
        (
            fifth_fields
                .into_iter()
                .map(normalize_element)
                .collect::<Result<_>>()?,
            6,
            7,
        )
    };
    let counts: Vec<usize> = lines[counts_line]
        .split_whitespace()
        .map(str::parse)
        .collect::<std::result::Result<_, _>>()?;
    if elements.len() != counts.len() {
        bail!("POSCAR element and atom-count lists have different lengths");
    }
    if lines[cursor].trim().to_ascii_lowercase().starts_with('s') {
        cursor += 1;
    }
    let direct = match lines
        .get(cursor)
        .map(|line| line.trim().to_ascii_lowercase())
    {
        Some(mode) if mode.starts_with('d') => true,
        Some(mode) if mode.starts_with('c') || mode.starts_with('k') => false,
        _ => bail!("POSCAR coordinate mode must be Direct or Cartesian"),
    };
    cursor += 1;
    let atom_count: usize = counts.iter().sum();
    if lines.len() < cursor + atom_count {
        bail!("POSCAR declares {atom_count} atoms but has too few coordinate lines");
    }
    let mut atoms = Vec::with_capacity(atom_count);
    for (element, count) in elements.iter().zip(counts) {
        for _ in 0..count {
            let fields: Vec<f32> = lines[cursor]
                .split_whitespace()
                .take(3)
                .map(str::parse)
                .collect::<std::result::Result<_, _>>()?;
            if fields.len() < 3 {
                bail!("POSCAR coordinate line {} is incomplete", cursor + 1);
            }
            let point = Vec3::new(fields[0], fields[1], fields[2]);
            let point = if direct {
                cell.fractional_to_cartesian(point)
            } else {
                point * scale
            };
            atoms.push(Atom::new(element.clone(), point));
            cursor += 1;
        }
    }
    Ok(Structure {
        format: "poscar".to_owned(),
        frames: vec![Frame {
            comment: title,
            atoms,
            bonds: Vec::new(),
            infer_bonds: true,
            cell: Some(cell),
        }],
    })
}

fn determinant(vectors: [Vec3; 3]) -> f32 {
    let [a, b, c] = vectors;
    a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pdb_reads_models_metadata_cell_and_connectivity() {
        let source = "HEADER    TEST PROTEIN\nCRYST1   10.000   11.000   12.000  90.00  90.00 120.00\nMODEL        1\nATOM      1  N   GLY A   1      11.104  13.207   9.101  1.00 20.00           N\nATOM      2  CA  GLY A   1      12.000  13.000   9.500  1.00 20.00           C\nENDMDL\nCONECT    1    2\n";
        let structure = parse_pdb(source).unwrap();
        assert_eq!(structure.frames.len(), 1);
        assert_eq!(structure.frames[0].atoms[1].element, "C");
        assert_eq!(structure.frames[0].bonds.len(), 1);
        assert!(structure.frames[0].cell.is_some());
        assert_eq!(structure.inspect(None).frames[0].chains, ["A"]);
    }

    #[test]
    fn pdb_ignores_dummy_nmr_cell() {
        let source = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1\nATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N\n";
        let structure = parse_pdb(source).unwrap();
        assert!(structure.frames[0].cell.is_none());
    }

    #[test]
    fn pdb_classifies_polymer_water_ion_and_ligand_atoms() {
        let source = "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N\nHETATM    2  O   HOH A 101       2.000   0.000   0.000  1.00 20.00           O\nHETATM    3 NA    NA A 102       4.000   0.000   0.000  1.00 20.00          NA\nHETATM    4  C1  LIG A 103       6.000   0.000   0.000  1.00 20.00           C\n";
        let structure = parse_pdb(source).unwrap();
        let roles = structure.frames[0]
            .atoms
            .iter()
            .map(|atom| atom.role)
            .collect::<Vec<_>>();
        assert_eq!(
            roles,
            vec![
                AtomRole::Polymer,
                AtomRole::Water,
                AtomRole::Ion,
                AtomRole::Ligand
            ]
        );
    }

    #[test]
    fn mol_reads_explicit_bond_orders() {
        let source = "ethene\n  xyz-read\n\n  2  1  0  0  0  0            999 V2000\n    0.0 0.0 0.0 C\n    1.34 0.0 0.0 C\n  1  2  2\nM  END\n";
        let structure = parse_sdf(source).unwrap();
        assert_eq!(structure.frames[0].bonds[0].order, 2);
    }

    #[test]
    fn mol_accepts_the_empty_title_emitted_by_rdkit() {
        let source = "\n     RDKit          3D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3400    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  2  0\nM  END\n";
        let structure = parse_sdf(source).unwrap();
        assert_eq!(structure.frames[0].atoms.len(), 2);
        assert_eq!(structure.frames[0].bonds[0].order, 2);
    }

    #[test]
    fn sdf_reads_multiple_records_after_delimiter_newlines() {
        let mol = "methane\n  xyz-read\n\n  1  0  0  0  0  0            999 V2000\n    0.0 0.0 0.0 C\nM  END\n";
        let source = format!("{mol}$$$$\n{mol}$$$$\n");
        let structure = parse_sdf(&source).unwrap();
        assert_eq!(structure.frames.len(), 2);
    }

    #[test]
    fn mol_v3000_reads_atoms_and_connectivity() {
        let source = "ethene\n  xyz-read\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 C 1.34 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 2 1 2\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        let structure = parse_sdf(source).unwrap();
        assert_eq!(structure.frames[0].atoms.len(), 2);
        assert_eq!(structure.frames[0].bonds[0].order, 2);
    }

    #[test]
    fn mol2_reads_aromatic_connectivity() {
        let source = "@<TRIPOS>MOLECULE\nring\n2 1 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.ar 1 BEN\n2 C2 1.4 0 0 C.ar 1 BEN\n@<TRIPOS>BOND\n1 1 2 ar\n";
        let structure = parse_mol2(source).unwrap();
        assert_eq!(structure.frames[0].atoms[0].element, "C");
        assert_eq!(structure.frames[0].bonds.len(), 1);
        assert!(!structure.frames[0].infer_bonds);
    }

    #[test]
    fn poscar_converts_direct_coordinates() {
        let source = "NaCl\n1.0\n5 0 0\n0 5 0\n0 0 5\nNa Cl\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n";
        let structure = parse_poscar(source).unwrap();
        assert_eq!(
            structure.frames[0].atoms[1].position,
            Vec3::new(2.5, 2.5, 2.5)
        );
        assert!(structure.frames[0].cell.is_some());
    }

    #[test]
    fn poscar_supports_negative_volume_and_vasp4_elements() {
        let source = "Na Cl\n-125.0\n1 0 0\n0 1 0\n0 0 1\n1 1\nDirect\n0 0 0\n0.5 0.5 0.5\n";
        let structure = parse_poscar(source).unwrap();
        assert_eq!(structure.frames[0].atoms[1].element, "Cl");
        let cell = structure.frames[0].cell.unwrap();
        assert!((cell.lengths[0] - 5.0).abs() < 0.001);
    }
}
