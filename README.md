# GZ-toolkit

A LAMMPS structure manipulation toolkit for computational materials science.

`GZ-toolkit` (import name: `gz_toolkit`) provides tools for:
- **Building crystal structures** — seed perfect unit cells (FCC, BCC, HCP, diamond, rocksalt, zinc blende, wurtzite, fluorite, CsCl) and write LAMMPS data files.
- **Creating defects** — voids, ⟨111⟩ and ⟨100⟩ dislocation loops (SIA and vacancy), Frenkel pairs, and interstitial insertion.
- **Merging structures** — carve-and-transplant regions between host/donor configurations with automatic overlap scoring and boundary optimization.
- **Geometric selection** — select atoms by sphere, cylinder, plane, cube, atom type, coordinates, percentage, or atom/molecule ID — all through a unified `region_mask()` dispatcher.
- **Coordinate transforms** — symmetry operations and integer supercell transformations.

## Installation

### Prerequisites

`GZ-toolkit` depends on [`myLAMMPS`](https://github.com/TaoLiang120/myLAMMPS.git), a LAMMPS data I/O library. Install it first:

```bash
pip install git+https://github.com/TaoLiang120/myLAMMPS.git
```

### Install GZ-toolkit

```bash
# From the project root
pip install -e .

# With development dependencies (pytest)
pip install -e ".[dev]"
```

### Dependencies

| Package | Purpose |
|---------|---------|
| [numpy](https://numpy.org/) | Array operations |
| [pandas](https://pandas.pydata.org/) | Atom table management |
| [scipy](https://scipy.org/) | Fast neighbor lookups (KD-trees) |
| [pymatgen](https://pymatgen.org/) | Crystal structure generation and lattice operations |
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
    loop_type="interstitial",
    n_layers=1,
)

data.to_file("with_loop.data")
```

## Simulation Workflow

`gz_toolkit` (mainly through `potential_testing/`) drives a LAMMPS potential-testing
run in four stages:

```
1. Create structure          2. Create input file      3. Create + identify        4. Submit
   (buildmtx / defect)          (MANUAL step)              directory tree              ┌─ (a) Python manager: auto-submit
                                                             (jobs / potential_testing)  └─ (b) Job files + .sh: manual submit
```

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

4. **Submit** — controlled by `workflow.paral_degree` in the potential's config:
   - **(a) Python script as manager (automatic, chained submission)** —
     `paral_degree = 3` or `4`. `emit_jobs()` / `scatter_cases()`
     (`potential_testing/parallel.py`) write per-stage or per-case `.job` files
     *and* call `sbatch` themselves, chaining reference → build-defects → cases →
     summary with `--dependency=afterok` so the whole pipeline runs unattended
     once you kick off the first job.
   - **(b) Job files + shell script (manual submission)** — `paral_degree = 1` or
     `2`, or `ScatterSubmitter`/`PatchSubmitter` (`jobs/submission.py`). Every
     `.job` file is written up front alongside one `submit_all.sh`; nothing is
     auto-submitted — `jobs/submission.py` explicitly never calls `sbatch`, so you
     review and run `bash submit_all.sh` yourself.

Driven end-to-end via the CLI:

```bash
python -m gz_toolkit.potential_testing.cli init --root . --pot-name my_pot   # scaffold project
python -m gz_toolkit.potential_testing.cli run --pot-inputs pot_inputs       # stage 1 + emit jobs
bash submit_all.sh                                                          # stage 4(b), if paral_degree in (1,2)
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
│   └── fcc_defect.py         # FCCDefect (placeholder)
├── jobs/                     # SLURM job templates, directory + submission helpers
└── potential_testing/        # potential-testing pipeline (reference → defects → summary)

tests/                        # pytest suite
examples/                     # usage example scripts
```

## Running Tests

```bash
pytest tests/ -v
```

## License

MIT
