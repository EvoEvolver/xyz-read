use std::collections::HashMap;

use anyhow::{bail, Context, Result};

use crate::{
    geometry::Vec3,
    model::{
        deduplicate_bonds, normalize_element, Atom, AtomRole, Bond, Frame, Structure, UnitCell,
    },
};

#[derive(Debug)]
struct CifLoop {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
}

pub fn parse(source: &str) -> Result<Structure> {
    let tokens = tokenize(source)?;
    let mut scalars = HashMap::new();
    let mut loops = Vec::new();
    let mut data_name = String::new();
    let mut cursor = 0;

    while cursor < tokens.len() {
        let token = &tokens[cursor];
        let lower = token.to_ascii_lowercase();
        if let Some(name) = lower.strip_prefix("data_") {
            if data_name.is_empty() {
                data_name = name.to_owned();
            }
            cursor += 1;
        } else if lower == "loop_" {
            cursor += 1;
            let mut headers = Vec::new();
            while cursor < tokens.len() && tokens[cursor].starts_with('_') {
                headers.push(canonical_name(&tokens[cursor]));
                cursor += 1;
            }
            if headers.is_empty() {
                bail!("CIF loop_ has no column names");
            }
            let mut values = Vec::new();
            while cursor < tokens.len()
                && !(is_control_token(&tokens[cursor]) && values.len() % headers.len() == 0)
            {
                values.push(tokens[cursor].clone());
                cursor += 1;
            }
            if values.len() % headers.len() != 0 {
                bail!(
                    "CIF loop {:?} with {} columns has {} values",
                    headers,
                    headers.len(),
                    values.len()
                );
            }
            let rows = values
                .chunks(headers.len())
                .map(<[String]>::to_vec)
                .collect();
            loops.push(CifLoop { headers, rows });
        } else if token.starts_with('_') {
            let value = tokens
                .get(cursor + 1)
                .with_context(|| format!("CIF item {token} has no value"))?;
            scalars.insert(canonical_name(token), value.clone());
            cursor += 2;
        } else {
            cursor += 1;
        }
    }

    let cell = parse_cell(&scalars)?;
    let atom_loop = loops
        .iter()
        .find(|table| {
            has_column(table, "_atom_site_cartn_x") || has_column(table, "_atom_site_fract_x")
        })
        .context("CIF contains no atom_site coordinate loop")?;
    let cartesian = has_column(atom_loop, "_atom_site_cartn_x");
    if !cartesian && cell.is_none() {
        bail!("fractional CIF coordinates require unit-cell parameters");
    }

    let x_name = if cartesian {
        "_atom_site_cartn_x"
    } else {
        "_atom_site_fract_x"
    };
    let y_name = if cartesian {
        "_atom_site_cartn_y"
    } else {
        "_atom_site_fract_y"
    };
    let z_name = if cartesian {
        "_atom_site_cartn_z"
    } else {
        "_atom_site_fract_z"
    };
    let mut grouped: Vec<(String, Vec<Atom>)> = Vec::new();

    for row in &atom_loop.rows {
        let model = value(atom_loop, row, &["_atom_site_pdbx_pdb_model_num"])
            .filter(|value| !missing(value))
            .unwrap_or("1");
        let group_index = grouped
            .iter()
            .position(|(name, _)| name == model)
            .unwrap_or_else(|| {
                grouped.push((model.to_owned(), Vec::new()));
                grouped.len() - 1
            });
        let raw_element = value(
            atom_loop,
            row,
            &[
                "_atom_site_type_symbol",
                "_atom_site_label_atom_id",
                "_atom_site_label",
            ],
        )
        .context("CIF atom row has no element or atom label")?;
        let element = element_from_cif(raw_element)?;
        let raw = Vec3::new(
            parse_number(required(atom_loop, row, x_name)?)?,
            parse_number(required(atom_loop, row, y_name)?)?,
            parse_number(required(atom_loop, row, z_name)?)?,
        );
        let position = if cartesian {
            raw
        } else {
            cell.expect("fractional coordinates checked for cell")
                .fractional_to_cartesian(raw)
        };
        let mut atom = Atom::new(element, position);
        atom.label = value(
            atom_loop,
            row,
            &[
                "_atom_site_label",
                "_atom_site_label_atom_id",
                "_atom_site_auth_atom_id",
            ],
        )
        .filter(|value| !missing(value))
        .map(str::to_owned);
        atom.chain = value(
            atom_loop,
            row,
            &["_atom_site_auth_asym_id", "_atom_site_label_asym_id"],
        )
        .filter(|value| !missing(value))
        .map(str::to_owned);
        let residue_name = value(
            atom_loop,
            row,
            &["_atom_site_auth_comp_id", "_atom_site_label_comp_id"],
        )
        .filter(|value| !missing(value));
        let residue_number = value(
            atom_loop,
            row,
            &["_atom_site_auth_seq_id", "_atom_site_label_seq_id"],
        )
        .filter(|value| !missing(value));
        if residue_name.is_some() || residue_number.is_some() {
            atom.residue = Some(format!(
                "{}{}",
                residue_name.unwrap_or(""),
                residue_number.unwrap_or("")
            ));
        }
        atom.residue_name = residue_name.map(str::to_owned);
        atom.residue_number = residue_number.and_then(|value| value.parse().ok());
        let group = value(atom_loop, row, &["_atom_site_group_pdb"])
            .unwrap_or("")
            .to_ascii_uppercase();
        atom.role = match group.as_str() {
            "ATOM" => AtomRole::Polymer,
            "HETATM" => crate::formats::hetero_atom_role(residue_name.unwrap_or("")),
            _ => AtomRole::Unknown,
        };
        atom.serial =
            value(atom_loop, row, &["_atom_site_id"]).and_then(|value| value.parse().ok());
        grouped[group_index].1.push(atom);
    }

    let mut frames: Vec<Frame> = grouped
        .into_iter()
        .map(|(model, atoms)| Frame {
            comment: if data_name.is_empty() {
                format!("model {model}")
            } else if model == "1" {
                data_name.clone()
            } else {
                format!("{data_name}, model {model}")
            },
            atoms,
            bonds: Vec::new(),
            infer_bonds: true,
            cell,
        })
        .collect();

    if let Some(bond_loop) = loops.iter().find(|table| {
        has_column(table, "_geom_bond_atom_site_label_1")
            && has_column(table, "_geom_bond_atom_site_label_2")
    }) {
        for frame in &mut frames {
            let indices: HashMap<&str, usize> = frame
                .atoms
                .iter()
                .enumerate()
                .filter_map(|(index, atom)| atom.label.as_deref().map(|label| (label, index)))
                .collect();
            for row in &bond_loop.rows {
                let first_label = required(bond_loop, row, "_geom_bond_atom_site_label_1")?;
                let second_label = required(bond_loop, row, "_geom_bond_atom_site_label_2")?;
                if let (Some(&first), Some(&second)) =
                    (indices.get(first_label), indices.get(second_label))
                {
                    let order = value(bond_loop, row, &["_geom_bond_type", "_geom_bond_order"])
                        .map(bond_order)
                        .unwrap_or(1);
                    frame.bonds.push(Bond::new(first, second, order));
                }
            }
            deduplicate_bonds(&mut frame.bonds);
        }
    }

    if frames.is_empty() || frames.iter().all(|frame| frame.atoms.is_empty()) {
        bail!("CIF atom_site loop contains no atoms");
    }
    Ok(Structure {
        format: if atom_loop
            .headers
            .iter()
            .any(|header| header == "_atom_site_group_pdb")
        {
            "mmcif"
        } else {
            "cif"
        }
        .to_owned(),
        frames,
    })
}

fn parse_cell(scalars: &HashMap<String, String>) -> Result<Option<UnitCell>> {
    let names = [
        "_cell_length_a",
        "_cell_length_b",
        "_cell_length_c",
        "_cell_angle_alpha",
        "_cell_angle_beta",
        "_cell_angle_gamma",
    ];
    if names.iter().all(|name| !scalars.contains_key(*name)) {
        return Ok(None);
    }
    let mut values = [0.0; 6];
    for (index, name) in names.iter().enumerate() {
        values[index] = parse_number(
            scalars
                .get(*name)
                .with_context(|| format!("CIF cell is missing {name}"))?,
        )?;
    }
    Ok(Some(UnitCell::from_parameters(
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
    )?))
}

fn value<'a>(table: &CifLoop, row: &'a [String], names: &[&str]) -> Option<&'a str> {
    names.iter().find_map(|name| {
        table
            .headers
            .iter()
            .position(|header| header == name)
            .and_then(|index| row.get(index))
            .map(String::as_str)
    })
}

fn required<'a>(table: &CifLoop, row: &'a [String], name: &str) -> Result<&'a str> {
    value(table, row, &[name])
        .filter(|value| !missing(value))
        .with_context(|| format!("CIF atom row is missing {name}"))
}

fn has_column(table: &CifLoop, name: &str) -> bool {
    table.headers.iter().any(|header| header == name)
}

fn missing(value: &str) -> bool {
    matches!(value, "." | "?")
}

fn parse_number(value: &str) -> Result<f32> {
    let without_uncertainty = value.split('(').next().unwrap_or(value);
    let parsed: f32 = without_uncertainty
        .parse()
        .with_context(|| format!("invalid CIF number {value:?}"))?;
    if !parsed.is_finite() {
        bail!("CIF number must be finite, found {value:?}");
    }
    Ok(parsed)
}

fn element_from_cif(value: &str) -> Result<String> {
    let letters: String = value
        .chars()
        .skip_while(|character| character.is_ascii_digit())
        .take_while(|character| character.is_ascii_alphabetic())
        .take(2)
        .collect();
    normalize_element(&letters).with_context(|| format!("invalid CIF element {value:?}"))
}

fn bond_order(value: &str) -> u8 {
    match value.to_ascii_lowercase().as_str() {
        "2" | "d" | "double" => 2,
        "3" | "t" | "triple" => 3,
        _ => 1,
    }
}

fn is_control_token(value: &str) -> bool {
    let lower = value.to_ascii_lowercase();
    value.starts_with('_')
        || lower == "loop_"
        || lower == "stop_"
        || lower.starts_with("data_")
        || lower.starts_with("save_")
}

fn canonical_name(value: &str) -> String {
    value.to_ascii_lowercase().replace('.', "_")
}

fn tokenize(source: &str) -> Result<Vec<String>> {
    let mut tokens = Vec::new();
    let mut lines = source.lines().peekable();
    while let Some(line) = lines.next() {
        if let Some(initial) = line.strip_prefix(';') {
            let mut value = String::new();
            if !initial.is_empty() {
                value.push_str(initial);
                value.push('\n');
            }
            let mut terminated = false;
            for next in lines.by_ref() {
                if next.starts_with(';') {
                    terminated = true;
                    break;
                }
                value.push_str(next);
                value.push('\n');
            }
            if !terminated {
                bail!("unterminated semicolon-delimited CIF value");
            }
            tokens.push(value.trim_end().to_owned());
            continue;
        }
        let mut characters = line.char_indices().peekable();
        while let Some((start, character)) = characters.next() {
            if character.is_whitespace() {
                continue;
            }
            if character == '#' {
                break;
            }
            if character == '\'' || character == '"' {
                let quote = character;
                let content_start = start + character.len_utf8();
                let mut end = None;
                for (index, next) in characters.by_ref() {
                    if next == quote {
                        end = Some(index);
                        break;
                    }
                }
                let end = end.context("unterminated quoted CIF value")?;
                tokens.push(line[content_start..end].to_owned());
            } else {
                let mut end = line.len();
                while let Some(&(index, next)) = characters.peek() {
                    if next.is_whitespace() || next == '#' {
                        end = index;
                        break;
                    }
                    characters.next();
                }
                tokens.push(line[start..end].to_owned());
                if line.as_bytes().get(end) == Some(&b'#') {
                    break;
                }
            }
        }
    }
    Ok(tokens)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_fractional_crystal_and_explicit_bonds() {
        let source = "data_quartz\n_cell_length_a 5.0\n_cell_length_b 5.0\n_cell_length_c 6.0\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 120\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nSi1 Si 0 0 0\nO1 O 0.5 0.5 0.5\nloop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n_geom_bond_type\nSi1 O1 D\n";
        let structure = parse(source).unwrap();
        assert_eq!(structure.format, "cif");
        assert_eq!(structure.frames[0].atoms.len(), 2);
        assert_eq!(structure.frames[0].bonds.len(), 1);
        assert_eq!(structure.frames[0].bonds[0].order, 2);
        assert!(structure.frames[0].cell.is_some());
    }

    #[test]
    fn parses_mmcif_models_and_protein_metadata() {
        let source = "data_test\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.pdbx_PDB_model_num\nATOM 1 N N GLY A 1 0 0 0 1\nATOM 2 C CA GLY A 1 1.4 0 0 1\nATOM 3 N N GLY A 1 0.1 0 0 2\n";
        let structure = parse(source).unwrap();
        assert_eq!(structure.format, "mmcif");
        assert_eq!(structure.frames.len(), 2);
        assert_eq!(structure.inspect(None).frames[0].residue_count, 1);
    }
}
