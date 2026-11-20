# Step 0 — Working environment

The same workflow runs on very different machines (a Windows laptop, a Linux workstation, an
HPC login node over SSH). Paths to Python, venvs, conda envs and local packages are different on
each. **Never carry paths from one machine to another.**

## 1. Identify the machine

Look at the evidence before asking:

| Evidence | Likely situation |
|----------|------------------|
| Folder connected through the desktop app, paths like `C:\Users\...` | local Windows |
| `/Users/<name>/...` | local macOS |
| A shell where `sbatch`/`squeue` exist, `/home/<user>` or `/scratch/...` | HPC login node (SSH) |
| User says "on the cluster", "on <cluster name>" | SSH to a cluster |

If it is still ambiguous: **ask** — "Is this run on your local machine, or over SSH on a
cluster? If a cluster, which one?"

Quick probes once you have a shell there:

```bash
uname -a 2>/dev/null || ver                 # OS
which sbatch squeue 2>/dev/null             # SLURM present?
module --version 2>&1 | head -1             # Lmod / environment modules?
which python python3 conda 2>/dev/null
```

## 2. Python interpreter and packages

The workflow scripts need only the standard library (+ PyYAML for YAML plans); structure building
needs GZ_toolkit and its deps (usually numpy, maybe ASE/pymatgen); summaries may need `ovito`.

1. Check "Known machines" below. If the machine is listed, use it and confirm with one probe.
2. Otherwise ask the user for:
   - the interpreter: a `.venv` path, a conda env name/prefix, or `module load python/x.y`
   - paths to local packages that are not pip-installed (GZ_toolkit, myLAMMPS, SEAKMC_py, ...)
3. Probe:
   ```bash
   <python> -c "import sys; print(sys.executable, sys.version)"
   <python> -c "import sys; sys.path.insert(0, r'<GZ_toolkit path>'); import <gz package>; print('ok')"
   <python> -c "import yaml, numpy; print('yaml', yaml.__version__)"
   ```
4. Record in `sim_plan.yaml: environment:` (python, python_kind, activate line, modules,
   packages). Offer to add a new machine to the list below so the next chat finds it.

How to call Python on each kind of setup:

| Kind | Call | In a job script |
|------|------|-----------------|
| venv (Linux/mac) | `/path/.venv/bin/python script.py` | `source /path/.venv/bin/activate` |
| venv (Windows) | `C:\path\venv\Scripts\python.exe script.py` | — |
| conda | `conda run -p <prefix> python script.py` or activate first | `source <conda>/etc/profile.d/conda.sh && conda activate <env>` |
| module | `module load python/x.y && python script.py` | same `module load` lines in the header |

Local packages go on `sys.path` at the top of every script that needs them (the structure
template does this from `environment.packages` automatically):

```python
import sys
sys.path.insert(0, r"<path to GZ_toolkit>")
```

### Package re-reading policy (saves effort)

| Package | Policy |
|---------|--------|
| GZ_toolkit | Under active development — list its modules / public functions when starting a structure task, and re-read when an import or attribute fails or the user mentions changes. |
| myLAMMPS, SEAKMC_py(_parallel) | Stable — re-read only on ImportError/AttributeError or if the user says it changed. |

## 3. HPC specifics (only when on a cluster)

Ask once, record in `sim_plan.yaml` (used by the reference job file): partition, account,
cores per node, typical walltime limit, the module lines for the simulation code, and the
scratch filesystem where runs should live. If a code is missing on the cluster, go to
`references/hpc-install.md`.

---

## Known machines

Add an entry whenever the user confirms a new machine. Entries are facts about one machine only.

### zhang — local Windows workstation
- OS: Windows; work is local (no SLURM)
- Python venv: `C:\Users\zhang\Documents\localVenv\pyvenv1`
  - run: `C:\Users\zhang\Documents\localVenv\pyvenv1\Scripts\python.exe script.py`
  - pip: `C:\Users\zhang\Documents\localVenv\pyvenv1\Scripts\pip.exe install <pkg>`
  - activate: `C:\Users\zhang\Documents\localVenv\pyvenv1\Scripts\activate.bat`
- Local packages (add to `sys.path`):
  - GZ-toolkit (import `gz_toolkit`): `C:\Users\zhang\Downloads\pycharm\GZ-toolkit`
  - myLAMMPS: `C:\Users\zhang\Downloads\pycharm\myLAMMPS`
  - SEAKMC_py_parallel: `C:\Users\zhang\Downloads\pycharm\SEAKMC_py_parallel`

<!-- template for a new machine
### <user> — <machine / cluster name>
- OS / access: linux, ssh <user>@<host>; scheduler: slurm
- Python: <venv path | conda env | module>
- activate in jobs: <line>
- Local packages: GZ_toolkit: <path>; ...
- SLURM: partition <p>, account <a>, <n> cores/node, max walltime <t>
- Code modules: lammps: `module load ...`; vasp: `module load ...`
- Scratch: <path>
-->
