# Installing GZ-toolkit on an HPC

Release 0.2.0 includes `eamGen`, potential testing, and the existing structure
tools. The Python package requires Python 3.9 or newer. Use the cluster's
supported Python module and a virtual environment; the Windows environment
configured on Zhang's computer is not portable to Linux.

## Install from your Git checkout

Upload the repository, including `pyproject.toml`, `MANIFEST.in`, `gz_toolkit/`
and `examples/`. On the HPC, load the appropriate Python module, then run:

```bash
cd /path/to/GZ-toolkit
python -m venv "$HOME/venvs/gz-toolkit"
source "$HOME/venvs/gz-toolkit/bin/activate"
python -m pip install --upgrade pip
python -m pip install .
```

Installation resolves NumPy, pandas, SciPy, pymatgen and the `myLAMMPS` Git
dependency automatically. Git and network access are required to fetch that
dependency. For development, replace the final command with
`python -m pip install -e ".[dev]"`. Keep a checkout for the example files.

On clusters without network access, prepare Linux-compatible dependency wheels
on a connected machine using the same Python version/platform. Install the
dependency wheels first, including `myLAMMPS`, then install the toolkit with
`python -m pip install --no-deps /path/to/GZ-toolkit`. The direct Git dependency
cannot be fetched on an offline cluster; copying a Windows virtual environment
will not work.

## Check installation without calculations

```bash
python -m pip check
gz-toolkit-eamgen --help
python -m gz_toolkit.eamGen --help
gz-toolkit-potential-testing --help
python -c 'import gz_toolkit; print(gz_toolkit.__version__)'
```

These checks do not fit a potential or start reference calculations. If a console
command is missing, activate the environment used for installation; the module
form uses that environment's Python directly.

## Configure external calculators

LAMMPS, its required pair styles/plugins, VASP and MPI are separate cluster
installations. This package supplies neither their binaries nor VASP POTCARs.
Load the required modules in your batch jobs and verify that the LAMMPS build
supports `eam/he` for the current Pd–He–Ni study. Provide licensed VASP inputs
with matching POTCAR/POSCAR species order when using DFT references.

Read [the eamGen manual](examples/eamgen_manual.md) before preparing a fit.
The [Pd–He–Ni example](examples/eamgen_pd_he_ni/README.md) describes the starter
files. Edit account, partition, environment and resources in `fit_hpc.job`, and
submit from an environment in which this package is installed:

```bash
sbatch --export=ALL,EAMGEN_PROJECT=/absolute/path/pd_he_ni_fit,EAMGEN_REVISION=r001 examples/eamgen_pd_he_ni/fit_hpc.job
```

Run this only after replacing the project's empty dataset with labelled,
documented reference data and reviewing the fitting options. The example job
performs fitting only; reference calculations and physical testing use separate
jobs. All fitting and reference calculations for this study belong on the HPC.

## Agent skill and release contents

The package ships `gz_toolkit/eamGen/skills/eamgen/SKILL.md` and its full manual.
Copy the entire `eamgen` directory to your host agent's skill location if needed;
Python installation does not install a skill globally. To locate the resources:

```bash
python -c 'from importlib.resources import files; print(files("gz_toolkit.eamGen").joinpath("skills/eamgen"))'
```

Git checkouts and source distributions include `examples/`, this installation
guide and the changelog. Wheels contain the executable package, testing/job
templates, agent skill and its manual; keep the checkout/source distribution
for the example JSON files and batch job. No supplied private potential file,
paper, generated fit or reference dataset is bundled.
