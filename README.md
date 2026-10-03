# GZ-toolkit

Atomistic structure tools, HPC potential testing, and guided EAM potential generation for computational materials science.

`GZ-toolkit` (import name: `gz_toolkit`) provides tools for:

- **Building crystal structures** — seed perfect unit cells (FCC, BCC, HCP, diamond, rocksalt, zinc blende, wurtzite, fluorite, CsCl) and write LAMMPS data files.
- **Creating defects** — voids, ⟨111⟩ and ⟨100⟩ dislocation loops (SIA and vacancy), Frenkel pairs, and interstitial insertion.
- **Merging structures** — carve-and-transplant regions between host/donor configurations with automatic overlap scoring and boundary optimization.
- **Geometric selection** — select atoms by sphere, cylinder, plane, cube, atom type, coordinates, percentage, or atom/molecule ID — all through a unified `region_mask()` dispatcher.
- **Coordinate transforms** — symmetry operations and integer supercell transformations.
- **EAM potential generation** — `eamGen` imports/exports EAM, FS and HE potentials,
  fits selected spline corrections or supported analytic parameters to energy/force/stress references, and hands
  reviewed candidates to potential testing. Includes a Pd–He–Ni starting workflow
  and an agent skill. See [the eamGen manual](examples/eamgen_manual.md).

## Installation

Use Python 3.9 or newer in a dedicated environment. From the project root:

```bash
python -m pip install .
# For development:
python -m pip install -e ".[dev]"
```

Dependencies, including the `myLAMMPS` Git dependency, are installed automatically.
Git and network access are needed for that dependency. LAMMPS/VASP binaries,
plugins, MPI and VASP pseudopotentials must be provided separately on the cluster.
See [INSTALL.md](INSTALL.md) for HPC setup, offline installation considerations,
installation checks and agent-skill installation.

```bash
gz-toolkit-eamgen --help
gz-toolkit-potential-testing --help
```

### Dependencies

| Package | Purpose |
|---------|---------|
| [numpy](https://numpy.org/) | Array operations |
| [pandas](https://pandas.pydata.org/) | Atom table management |
| [scipy](https://scipy.org/) | Neighbor lookups, spline interpolation and bounded optimization |
| [pymatgen](https://pymatgen.org/) | Crystal/lattice operations and VASP reference parsing |
| [myLAMMPS](https://github.com/TaoLiang120/myLAMMPS.git) | LAMMPS data file I/O (`lmpData`, `lmpBox`) |

## Quick Start

### Build a BCC seed crystal

```python
from gz_toolkit import Gen_crystal

g = Gen_crystal(atom_style="atomic")
g.seed_crystal("A2", ["Fe", "Fe"], "bcc_seed.data")
```

### Load and modify a LAMMPS data file

```python
from gz_toolkit import Modlmp_LmpData

data = Modlmp_LmpData.from_file("my_structure.data", "atomic")

# Select atoms in a sphere
mask = data.select_sphere(center=[10, 10, 10], radius=5.0)
print(f"Selected {mask.sum()} atoms")

# Modify atom types in a region
data.mod_atom_type(
    "sphere",
    new_type=2,
    ff_elements=["Fe", "He"],
    atomic_masses=[55.845, 4.003],
    center=[10, 10, 10],
    radius=5.0,
)

data.to_file("modified.data")
```

### Create a void

```python
from gz_toolkit import Modlmp_LmpData, BCCDefect

data = Modlmp_LmpData.from_file("bcc_supercell.data", "atomic")
defect = BCCDefect(data)
n_deleted = defect.void(center=None, radius=8.0, region_type="sphere")
print(f"Deleted {n_deleted} atoms")
data.to_file("with_void.data")
```

### Create dislocation loops

```python
from gz_toolkit import Modlmp_LmpData, BCCDefect

data = Modlmp_LmpData.from_file("bcc_supercell.data", "atomic")
defect = BCCDefect(data)

# Interstitial loop
defect.add_111_loop(
    center=[57.1, 57.1, 57.1],
    habit_plane=(1, 1, 1),
    radius=10.0,
    lattice_const=2.8553,
    ff_elements=["Fe", "Fe"],
    atomic_masses=[55.845, 55.845],
    loop_type="sil",
    n_repeats=1,
)

data.to_file("with_loop.data")
```

## Simulation Workflow

`gz_toolkit` (mainly through `potential_testing/`) drives a LAMMPS potential-testing
run in four stages:

Structures → potential configuration → generated inputs/jobs → reviewed SLURM submission → summaries

1. **Create structure** — `Gen_crystal` (`buildmtx/buildstr.py`) seeds a crystal cell;
   `BCCDefect`/`FCCDefect` (`defect/`) and `structure_ops.py` carve out point defects,
   dislocation loops, alloys, and gas complexes as LAMMPS data files.

2. **Create the input file — manual** — you supply the interatomic potential
   (`pot_lines` and any potential files) in `pot_inputs/<pot>.json`. Everything else —
   `in.reference.lammps`, `in.relax.lammps`, `potential.inc`, etc. — is generated
   automatically by `potential_testing/pipeline.py`.

3. **Create and identify the directory tree** — `prepare_reference_stage()` and
   `build_defects_post_reference()` (`potential_testing/pipeline.py`) lay out
   `<pot>/reference/<element>/`, `<pot>/point_defects/<element>/<case>/`,
   `<pot>/loops/...`, `<pot>/alloy_lc/...`, `<pot>/gas_complexes/...`, etc.
   That tree is later re-*identified* (not just built) by:
   - `DirectoryManager` (`jobs/directory.py`) — generic traversal/filtering of run
     subdirectories by naming convention (`prefix-suffix`, `flat`, or a custom regex).
   - `_list_case_dirs()` (`potential_testing/parallel.py`) — walks the
     `potential_testing`-specific case tree to build job commands.
   - `collect_status()` (`potential_testing/status.py`) — walks the same tree and
     classifies each case as `pending` / `running` / `done` / `failed`.
   - `generate_folder_tree()` (`jobs/tree.py`) — writes a human-readable
     `folder_tree.txt` snapshot of any directory.

4. **Submit** — `run` prepares job files and `submit_all.sh` for every
   `workflow.paral_degree`. Review the generated settings and run the launcher
   explicitly. Degree 1 uses a single pipeline job; degree 2 uses one per
   potential; degree 3 chains stage/group jobs with dependencies. Degree 4
   launches reference and build jobs, then the build/scatter job submits
   individual cases and dependent summaries at runtime. Review its total job
   count and resource settings before launching.

Driven end-to-end via the CLI:

```bash
python -m gz_toolkit.potential_testing.cli init --root . --pot-name my_pot   # scaffold project
python -m gz_toolkit.potential_testing.cli run --pot-inputs pot_inputs       # stage 1 + emit jobs
bash submit_all.sh                                                          # explicitly launch the reviewed jobs (all degrees)
python -m gz_toolkit.potential_testing.cli status --pot-inputs pot_inputs    # check progress anytime
```

## Package Structure

```
gz_toolkit/
├── __init__.py               # public API (lazy imports)
├── core/
│   └── modlmp.py             # Modlmp_LmpData class (selection, add/delete, transform, merge)
├── buildmtx/
│   └── buildstr.py           # Gen_crystal: seed crystal structures
├── defect/
│   ├── bcc_defect.py         # BCCDefect: void, loops, Frenkel pairs, transforms
│   └── fcc_defect.py         # FCCDefect: voids, loops and point defects
├── jobs/                     # SLURM job templates, directory + submission helpers
├── potential_testing/        # potential-testing pipeline (reference → defects → summary)
├── eamGen/                   # EAM/FS/HE fitting, reference adapters and agent skill
├── analyze/                  # formation/binding energies and elastic formulas
└── pot_infobank/             # validated potential metadata

tests/                        # pytest suite
examples/                     # usage scripts, eamGen manual and Pd–He–Ni starter files
```

## eamGen on the HPC

Start with [the full eamGen manual](examples/eamgen_manual.md) and
[the Pd–He–Ni starter files](examples/eamgen_pd_he_ni/README.md).
The initial study retains the supplied Pd–He–H `eam/he` formulation, uses a
Pd–Ni potential for reference data, and adds Ni–He DFT or suitable teacher
references. Cross-potential composition is provisional and requires review.

The workflow is: import/combine a baseline, prepare and label geometries,
fit a fresh revision on the HPC, review its errors with the user, hand it to
potential testing, and compare independent physical results with sourced targets.
The package records inputs/provenance and provides standardized agent commands;
the accompanying skill guides research and human decisions.

Current setfl support is `eam/alloy`, `eam/fs` and `eam/he`. Analytic forms are
limited to the supplied paper's He terms, and composition currently requires an
HE primary. Further generalization remains future work. No validated Pd–He–Ni
potential is included. For this study, run all fitting/reference calculations on
the HPC; examples and installation checks do not start them automatically.

## Running Tests

```bash
python -m pytest tests/ -v
```

The eamGen tests include synthetic optimization. Run them on the HPC for this
study. Installation/help checks in [INSTALL.md](INSTALL.md) do not perform fits.
See [CHANGELOG.md](CHANGELOG.md) for release changes.

## License

MIT
