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
        .stdout(predicate::str::contains("--frame last"));
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
fn reports_out_of_range_frames() {
    Command::cargo_bin("xyz-read")
        .unwrap()
        .args(["inspect", example().to_str().unwrap(), "--frame", "8"])
        .assert()
        .failure()
        .stderr(predicate::str::contains("this file has 2 frame(s)"));
}
