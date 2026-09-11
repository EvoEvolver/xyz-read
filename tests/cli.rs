use std::{fs, path::Path};

use assert_cmd::Command;
use predicates::prelude::*;

fn example() -> &'static Path {
    Path::new("examples/benzene-trajectory.xyz")
}

#[test]
fn top_level_help_exposes_agent_workflow() {
    Command::cargo_bin("xyz-read")
        .unwrap()
        .arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("zoom-in"))
        .stdout(predicate::str::contains("zoom-out"))
        .stdout(predicate::str::contains("rotate"))
        .stdout(predicate::str::contains("--atom-numbers"))
        .stdout(predicate::str::contains("--frame last"))
        .stdout(predicate::str::contains("PDB"))
        .stdout(predicate::str::contains("--unit-cell"));
}

#[test]
fn inspect_json_reports_all_frames() {
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args(["inspect", example().to_str().unwrap(), "--json"])
        .assert()
        .success()
        .stdout(predicate::str::contains("\"frame_count\": 2"))
        .stdout(predicate::str::contains("\"C\": 6"));
}

#[test]
fn renders_selected_frame_with_rotation_and_numbers() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("frame.png");
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "render",
            example().to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--frame",
            "last",
            "--width",
            "320",
            "--height",
            "240",
            "--rotate-y",
            "35",
            "--atom-numbers",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("frame 2/2"));

    let image = image::open(&output).unwrap();
    assert_eq!((image.width(), image.height()), (320, 240));
    assert!(fs::metadata(output).unwrap().len() > 1_000);
}

#[test]
fn zoom_and_rotate_commands_produce_distinct_views() {
    let directory = tempfile::tempdir().unwrap();
    let close = directory.path().join("close.png");
    let wide = directory.path().join("wide.png");
    let rotated = directory.path().join("rotated.png");

    for (command, output, extra) in [
        ("zoom-in", &close, vec!["--factor", "1.4"]),
        ("zoom-out", &wide, vec!["--factor", "1.4"]),
        ("rotate", &rotated, vec!["--view", "z", "--rotate-y", "45"]),
    ] {
        let mut arguments = vec![
            command,
            example().to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--width",
            "160",
            "--height",
            "160",
        ];
        arguments.extend(extra);
        Command::cargo_bin("xyz-read")
            .unwrap()
            .args(arguments)
            .assert()
            .success();
    }

    let close = image::open(close).unwrap().into_rgb8();
    let wide = image::open(wide).unwrap().into_rgb8();
    let rotated = image::open(rotated).unwrap().into_rgb8();
    assert_ne!(close.as_raw(), wide.as_raw());
    assert_ne!(close.as_raw(), rotated.as_raw());
    assert_ne!(wide.as_raw(), rotated.as_raw());
}

#[test]
fn zoom_report_contains_the_effective_zoom_once() {
    let directory = tempfile::tempdir().unwrap();
    let image = directory.path().join("close.png");
    let report = directory.path().join("report.json");

    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "zoom-in",
            example().to_str().unwrap(),
            "--output",
            image.to_str().unwrap(),
            "--report",
            report.to_str().unwrap(),
            "--zoom",
            "1.5",
            "--factor",
            "1.4",
        ])
        .assert()
        .success();

    let report: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(report).unwrap()).unwrap();
    assert_eq!(report["camera"]["zoom"], 2.1);
}

#[test]
fn reports_out_of_range_frames() {
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args(["inspect", example().to_str().unwrap(), "--frame", "8"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("this file has 2 frame(s)"));
}

#[test]
fn inspects_protein_metadata_and_crystal_cell() {
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args(["inspect", "examples/peptide.pdb", "--json"])
        .assert()
        .success()
        .stdout(predicate::str::contains("\"format\": \"pdb\""))
        .stdout(predicate::str::contains("\"residue_count\": 2"))
        .stdout(predicate::str::contains("\"A\""));

    Command::cargo_bin("xyz-read")
        .unwrap()
        .args(["inspect", "examples/nacl.cif", "--json"])
        .assert()
        .success()
        .stdout(predicate::str::contains("\"format\": \"cif\""))
        .stdout(predicate::str::contains("\"unit_cell\""));
}

#[test]
fn renders_explicit_connectivity_and_crystal_cell() {
    let directory = tempfile::tempdir().unwrap();
    let ligand = directory.path().join("ligand.png");
    let crystal = directory.path().join("crystal.png");

    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "render",
            "examples/ethene.sdf",
            "-o",
            ligand.to_str().unwrap(),
            "--width",
            "320",
            "--height",
            "240",
        ])
        .assert()
        .success();
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "render",
            "examples/nacl.cif",
            "-o",
            crystal.to_str().unwrap(),
            "--width",
            "320",
            "--height",
            "240",
            "--unit-cell",
        ])
        .assert()
        .success();

    assert!(fs::metadata(ligand).unwrap().len() > 1_000);
    assert!(fs::metadata(crystal).unwrap().len() > 1_000);
}

#[test]
fn format_override_handles_unknown_extensions() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("structure.data");
    fs::copy("examples/ethene.sdf", &input).unwrap();
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "inspect",
            input.to_str().unwrap(),
            "--format",
            "sdf",
            "--json",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("\"explicit_bond_count\": 5"));
}

#[test]
fn detects_extensionless_poscar() {
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args(["inspect", "examples/POSCAR", "--json"])
        .assert()
        .success()
        .stdout(predicate::str::contains("\"format\": \"poscar\""))
        .stdout(predicate::str::contains("\"Si\": 2"));
}

#[test]
fn unit_cell_requires_cell_data() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("invalid.png");
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "render",
            "examples/peptide.pdb",
            "-o",
            output.to_str().unwrap(),
            "--unit-cell",
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("has no unit-cell data"));
    assert!(!output.exists());
}

#[test]
fn query_exposes_rdkit_indices_and_preserves_explicit_bonds() {
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "query",
            "examples/ethene.sdf",
            "--select",
            "rdkit-index:0..1",
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("\"rdkit_index\": 0"))
        .stdout(predicate::str::contains("\"rdkit_index\": 1"))
        .stdout(predicate::str::contains("\"atom_count\": 2"))
        .stdout(predicate::str::contains("\"explicit_bond_count\": 1"))
        .stdout(predicate::str::contains("\"order\": 2"));
}

#[test]
fn query_summary_omits_large_atom_and_bond_arrays() {
    let output = Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "query",
            "examples/peptide.pdb",
            "--select",
            "around:2@residue:1",
            "--summary",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(value.get("atoms").is_none());
    assert!(value.get("explicit_bonds").is_none());
    assert!(value.get("resolved_bonds").is_none());
    assert!(value["atom_count"].as_u64().unwrap() > 0);
    assert!(value["residue_count"].as_u64().unwrap() > 0);
    assert!(value["residues"].as_array().is_some());
}

#[test]
fn stdin_accepts_rdkit_mol_and_sdf_data() {
    let source = fs::read_to_string("examples/ethene.sdf").unwrap();
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args(["inspect", "-", "--format", "sdf", "--json"])
        .write_stdin(source)
        .assert()
        .success()
        .stdout(predicate::str::contains("\"format\": \"sdf\""))
        .stdout(predicate::str::contains("\"schema_version\": 1"));
}

#[test]
fn plan_revise_render_report_is_reproducible() {
    let directory = tempfile::tempdir().unwrap();
    let plan = directory.path().join("view.json");
    let revised = directory.path().join("revised.json");
    let image = directory.path().join("image.png");
    let report = directory.path().join("report.json");

    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "plan",
            "examples/ethene.sdf",
            "-o",
            plan.to_str().unwrap(),
            "--select",
            "rdkit-index:0..1",
            "--highlight",
            "rdkit-index:0",
            "--atom-numbers",
            "--numbering",
            "rdkit",
            "--width",
            "180",
            "--height",
            "160",
        ])
        .assert()
        .success();
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "revise",
            plan.to_str().unwrap(),
            "-o",
            revised.to_str().unwrap(),
            "--zoom-by",
            "1.2",
            "--rotate-y-by",
            "25",
        ])
        .assert()
        .success();
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "render",
            "--spec",
            revised.to_str().unwrap(),
            "-o",
            image.to_str().unwrap(),
            "--report",
            report.to_str().unwrap(),
        ])
        .assert()
        .success()
        .stdout(predicate::str::contains("2 selected atom(s)"));

    let report = fs::read_to_string(report).unwrap();
    assert!(report.contains("\"selected_atoms\": 2"));
    assert!(report.contains("\"drawn_bonds\": 1"));
    assert!(report.contains("\"highlight\": \"rdkit-index:0\""));
    assert!(report.contains("\"zoom\": 1.2"));
    assert!(fs::metadata(image).unwrap().len() > 1_000);
}

#[test]
fn view_spec_rejects_changed_structure_input() {
    use std::io::Write;

    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("molecule.sdf");
    let spec = directory.path().join("view.json");
    let image = directory.path().join("image.png");
    fs::copy("examples/ethene.sdf", &input).unwrap();

    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "plan",
            input.to_str().unwrap(),
            "-o",
            spec.to_str().unwrap(),
        ])
        .assert()
        .success();
    fs::OpenOptions::new()
        .append(true)
        .open(&input)
        .unwrap()
        .write_all(b"\n")
        .unwrap();
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "render",
            "--spec",
            spec.to_str().unwrap(),
            "-o",
            image.to_str().unwrap(),
        ])
        .assert()
        .failure()
        .stderr(predicate::str::contains("input hash mismatch"));
    assert!(!image.exists());
}

#[test]
fn views_produces_labeled_candidates_and_replayable_specs() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("views");
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "views",
            "examples/ethene.sdf",
            "-o",
            output.to_str().unwrap(),
            "--width",
            "128",
            "--height",
            "128",
        ])
        .assert()
        .success();

    for name in [
        "A.png",
        "A.view.json",
        "F.png",
        "F.view.json",
        "contact-sheet.png",
        "manifest.json",
    ] {
        assert!(output.join(name).is_file(), "missing {name}");
    }
    let manifest = fs::read_to_string(output.join("manifest.json")).unwrap();
    assert!(manifest.contains("\"label\": \"A\""));
    assert!(manifest.contains("\"spec\": \"A.view.json\""));
}

#[test]
fn capabilities_declares_read_only_rdkit_boundary() {
    Command::cargo_bin("xyz-read")
        .unwrap()
        .arg("capabilities")
        .assert()
        .success()
        .stdout(predicate::str::contains("\"read_only\": true"))
        .stdout(predicate::str::contains("\"structure_editing\": false"))
        .stdout(predicate::str::contains("\"rdkit-index\""));
}

#[test]
fn bundled_machine_schemas_are_valid_json() {
    for kind in ["view", "query", "render-report", "candidate-manifest"] {
        let output = Command::cargo_bin("xyz-read")
            .unwrap()
            .args(["schema", kind])
            .output()
            .unwrap();
        assert!(output.status.success(), "schema command failed for {kind}");
        let schema: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(
            schema["$schema"],
            "https://json-schema.org/draft/2020-12/schema"
        );
    }
}

#[test]
fn selection_and_highlight_do_not_modify_the_structure() {
    let directory = tempfile::tempdir().unwrap();
    let input = directory.path().join("molecule.sdf");
    let output = directory.path().join("image.png");
    fs::copy("examples/ethene.sdf", &input).unwrap();
    let before = fs::read(&input).unwrap();

    Command::cargo_bin("xyz-read")
        .unwrap()
        .args([
            "render",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "--select",
            "all",
            "--highlight",
            "rdkit-index:0..1",
            "--width",
            "180",
            "--height",
            "160",
        ])
        .assert()
        .success();

    assert_eq!(fs::read(input).unwrap(), before);
    assert!(fs::metadata(output).unwrap().len() > 1_000);
}
