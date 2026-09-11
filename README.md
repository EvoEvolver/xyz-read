# xyz-read

`xyz-read` is a small command-line tool for inspecting and rendering molecules, proteins, crystals, and trajectories. It is designed for agent workflows: inspect structured metadata first, render a candidate, adjust the model, zoom, or rotation, and show the result to a human.

It runs headlessly and includes its own CPU 3D renderer. No GUI, browser, or external molecular viewer is needed.

![Numbered benzene rendered by xyz-read](docs/benzene-numbered.png)

PDB protein/peptide structures and crystallographic cells use the same renderer:

![Two-residue PDB peptide](docs/peptide.png)

![NaCl CIF with unit-cell boundary](docs/nacl-cell.png)

## Install

Linux and macOS users can install the latest prebuilt release with:

```sh
curl -fsSL https://raw.githubusercontent.com/EvoEvolver/xyz-read/main/install.sh | sh
```

The installer detects the operating system and CPU, downloads the matching release, verifies its SHA-256 checksum, and installs `xyz-read` to `~/.local/bin`.

To choose another directory or pin a release:

```sh
curl -fsSL https://raw.githubusercontent.com/EvoEvolver/xyz-read/main/install.sh | \
  XYZ_READ_INSTALL_DIR="$HOME/bin" XYZ_READ_VERSION=v0.2.0 sh
```

Prebuilt archives are also available on the [Releases](https://github.com/EvoEvolver/xyz-read/releases) page for:

- Linux x86-64 and ARM64
- macOS Intel and Apple Silicon

## Quick start

Inspect any supported structure file:

```sh
xyz-read inspect protein.pdb --json
xyz-read inspect crystal.cif --json
xyz-read inspect ligand.sdf --json
```

Render a molecule using connectivity stored in the file:

```sh
xyz-read render ligand.sdf -o ligand.png --atom-numbers
```

Render a protein model or trajectory frame:

```sh
xyz-read render ensemble.pdb -o model-7.png --frame 7 --rotate-y 30
```

Render a crystal together with its unit-cell boundary:

```sh
xyz-read render crystal.cif -o crystal.png --unit-cell
```

## Supported formats

| Format | Detection | Preserved data |
| --- | --- | --- |
| XYZ / extended XYZ | `.xyz` | Repeated frames; extra atom columns are accepted |
| PDB | `.pdb`, `.ent` | `MODEL` frames, `CONECT`, `CRYST1`, atom names, residues, chains |
| CIF / mmCIF | `.cif`, `.mmcif`, `.mcif` | Cartesian or fractional coordinates, models, cell parameters, `_geom_bond` connectivity, protein metadata |
| MOL / SDF | `.mol`, `.sdf`, `.sd` | V2000 and V3000 atoms, explicit bonds and bond orders, multiple SDF records as frames |
| MOL2 | `.mol2` | Multiple molecules, atom types, residue names, explicit connectivity |
| POSCAR / CONTCAR | those filenames, `.vasp`, `.poscar` | Direct or Cartesian coordinates and lattice vectors |

Format detection uses both the file name and its contents. Override it when a file has an unusual name:

```sh
xyz-read inspect downloaded.data --format mmcif --json
```

For the CLI, CIF and mmCIF share the `cif` format selector, and MOL/SDF may be selected separately as `mol` or `sdf`.

## Agent workflow

The JSON inspection output is intended to be the first step:

```sh
xyz-read inspect structure-file --json
```

It reports the detected format and total frame count. Each frame includes its 1-based number, comment, atom and explicit-bond counts, connectivity mode, element counts, residues, chains, coordinate bounds, and unit-cell parameters when present. An agent can then select an appropriate model and view without changing the input:

```sh
xyz-read render structure-file -o candidate.png \
  --frame 12 --rotate-y 35
```

After human review, the agent can adjust just the view:

```sh
xyz-read zoom-in structure-file -o closer.png \
  --frame 12 --rotate-y 35 --factor 1.4

xyz-read rotate structure-file -o alternate.png \
  --frame 12 --rotate-x 20 --rotate-y 70
```

Every image-producing command accepts the same format, frame, size, base view, rotation, atom-number, unit-cell, and bond options. This makes revisions explicit and reproducible.

## Connectivity

`xyz-read` distinguishes complete and partial connectivity:

- MOL, SDF, and MOL2 bond tables are treated as complete and rendered exactly as supplied. Double and triple bonds are drawn with parallel strokes.
- PDB `CONECT` and CIF `_geom_bond` records are preserved, then missing local bonds are inferred. This matters because protein PDB files commonly omit ordinary polymer bonds from `CONECT`.
- XYZ, POSCAR, and structures without bond records use spatial hashing and covalent radii to infer bonds efficiently.

Use `--no-bonds` to suppress both explicit and inferred bonds.

## Proteins

PDB `MODEL` blocks and mmCIF `_atom_site.pdbx_PDB_model_num` values become selectable frames. Alternate PDB locations other than blank, `A`, or `1` are skipped to avoid drawing duplicate conformers. Residue and chain summaries are available through `inspect --json`.

Atom numbering is off by default and is usually best left off for a whole protein. It remains useful for a selected small model or ligand:

```sh
xyz-read render peptide.pdb -o peptide.png --atom-numbers
```

## Crystals

CIF fractional coordinates and POSCAR direct coordinates are converted with the full triclinic lattice, including non-orthogonal cell angles. Show the cell boundary explicitly with:

```sh
xyz-read render structure.cif -o structure.png --unit-cell
```

`--unit-cell` participates in automatic framing. It returns an error rather than silently doing nothing when the selected frame has no cell data.

## Frames and trajectories

Standard repeated XYZ blocks are treated as animation or trajectory frames:

```text
3
frame 1
O  0.000  0.000  0.000
H  0.758  0.586  0.000
H -0.758  0.586  0.000
3
frame 2
O  0.000  0.000  0.050
H  0.780  0.600  0.000
H -0.780  0.600  0.000
```

Frame numbers start at 1. Use `--frame NUMBER` or `--frame last`. SDF records, PDB/mmCIF models, and repeated XYZ blocks all use the same frame interface.

## Zoom and rotation

Use an absolute zoom on any image command:

```sh
xyz-read render molecule.xyz -o image.png --zoom 1.8
```

Or use the intent-oriented commands:

```sh
xyz-read zoom-in molecule.xyz -o close.png --factor 1.5
xyz-read zoom-out molecule.xyz -o wide.png --factor 1.5
```

`zoom-in` multiplies `--zoom` by `--factor`; `zoom-out` divides it. The factor must be greater than 1.

Rotation values are degrees and compose across X, Y, and Z:

```sh
xyz-read rotate molecule.xyz -o rotated.png \
  --rotate-x -15 --rotate-y 45 --rotate-z 5
```

The default `--view iso` gives a three-quarter view. `--view x`, `--view y`, and `--view z` provide axis-aligned starting views; rotations are applied on top of that base view.

## Rendering options

```text
--frame NUMBER|last     Select a 1-based trajectory frame
--format FORMAT         auto, xyz, pdb, cif, mol, sdf, mol2, or poscar
--width PIXELS          Output width, 128 to 4096 (default: 1000)
--height PIXELS         Output height, 128 to 4096 (default: 750)
--zoom FACTOR           Absolute zoom, 0.05 to 20 (default: 1)
--view iso|x|y|z        Choose the base camera direction
--rotate-x DEGREES      Rotate around X
--rotate-y DEGREES      Rotate around Y
--rotate-z DEGREES      Rotate in the image plane
--atom-numbers          Overlay stable 1-based atom indices
--unit-cell             Draw and frame the crystallographic cell
--no-bonds              Disable explicit and inferred bonds
```

Run `xyz-read --help` or `xyz-read <command> --help` for the complete command reference.

## Renderer

The built-in renderer is implemented in this repository. It centers and rotates 3D coordinates, resolves explicit and inferred connectivity, auto-fits atoms and optional lattice boundaries, rasterizes shaded atom spheres and bonds with a depth buffer, and downsamples a larger working image for smoother edges. Atom numbers use an embedded regular sans-serif font, so output looks conventional and does not depend on system fonts.

Rendering is deterministic for the same input and arguments. PNG output is written atomically, so an interrupted render does not leave a partial destination file.

## License

MIT
