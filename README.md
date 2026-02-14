# shellqcc

[English](README.md) | [Chinese](README.zh-CN.md)

Shell utilities for high-throughput Quantum ESPRESSO and ORCA workflows.

## Overview

This repository provides Bash scripts to streamline common computational chemistry tasks:

- Generate and convert Quantum ESPRESSO inputs (`qei.sh`)
- Run QE workflows and auto-generate post-processing inputs (`qec.sh`)
- Prepare and batch-run ORCA jobs (`orcai.sh`)
- Collect ORCA Gibbs free energies from frequency outputs (`orcac.sh`)

The project also includes helper tools and examples under `tools/` and `examples/`.

## Repository Layout

- `qei.sh`: QE input generation/conversion/analysis toolkit
- `qec.sh`: QE workflow runner (`pw`, `dos`, `bands`, `w90`, `phonon`, etc.)
- `orcai.sh`: ORCA input preparation and execution helper
- `orcac.sh`: ORCA output parser for Gibbs free energies
- `examples/`: sample molecular and crystal structures, ORCA templates
- `tools/`: helper scripts for environment setup, input generation, and CIF handling

## Requirements

### System

- Bash
- `awk`, `sed`, `grep`, `bc`
- `mpirun` (for QE tasks)

### For `qei.sh`

- `cif2cell`
- ASE CLI (`ase`)
- Access to pseudopotential libraries (paths configurable inside script)

### For `qec.sh`

- Quantum ESPRESSO executables (`pw.x`, `dos.x`, `projwfc.x`, `pp.x`, `bands.x`, `ph.x`, `q2r.x`, `matdyn.x`, etc. as needed by task)
- Optional: `gnuplot`, Wannier90 (`wannier90.x`, `pw2wannier90.x`), `open_grid.x`, `sumpdos.x`, `thermo_pw.x`

### For `orcai.sh`

- ORCA (binary path configured in script)
- ASE CLI (`ase`) for `.mol` to `.xyz` conversion in xtb mode

## Quick Start

```bash
# 1) Generate QE input from CIF
./qei.sh structure.cif qe scf

# 2) Run QE default task (pw)
./qec.sh structure.scf.in

# 3) Run a QE task chain
./qec.sh structure.relax.in pw dos bands

# 4) Prepare ORCA jobs from template
./orcai.sh xtb.tmp

# 5) Extract Gibbs free energies from ORCA freq outputs
./orcac.sh
```

## Script Usage

### `qei.sh`

```bash
# CIF -> input (default: qe scf)
./qei.sh file.cif [qe|vasp|cp2k] [calc_type]

# QE input -> cif/conv
./qei.sh file.in [cif|conv]

# QE output -> analysis
./qei.sh file.out [scf|relax|vc-relax]
```

Notes:
- If a target `.cif` file does not exist, `qei.sh` can launch an interactive CIF generation wizard.
- The script contains interactive prompts for pseudopotential choice and related settings.

### `qec.sh`

```bash
./qec.sh input.in [tasks...]
```

Available tasks:
- `pw` (default), `dos`, `co`, `esp`, `hp`, `bands`, `w90`, `phonon`, `thermo`, `dry`, `all`

Example:

```bash
./qec.sh prefix.relax.in pw dos co esp bands
```

### `orcai.sh`

```bash
./orcai.sh <template_file>
```

Behavior:
- With `xtb.tmp`: converts `.mol` to `.xyz` when needed, builds ORCA inputs, creates job folders, runs ORCA, and copies generated `.xtb.xyz` back.
- With other templates: prepares ORCA input files from existing `.xyz` files.

### `orcac.sh`

```bash
./orcac.sh
```

Scans `*/*.out` files for `Final Gibbs free energy` and writes results to `Gibbs.txt`.

## Examples

- Crystal CIFs: `examples/crystal/`
- Molecular structures: `examples/molecular/`
- ORCA templates: `examples/orcai_tmp/`

## License

This project is released under the MIT License. See `LICENSE`.
