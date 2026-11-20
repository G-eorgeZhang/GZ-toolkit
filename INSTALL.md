# Installing GZ-toolkit and its simulation skill

The package provides atomistic structure tools, HPC potential testing, and a
materials-simulation skill for AI-human collaboration. It requires Python 3.9
or newer. Use a Python environment on the machine where the work will run;
a Windows virtual environment is not portable to Linux.

## Python package

From a checkout containing `pyproject.toml`, `MANIFEST.in`, and `gz_toolkit/`:

```bash
python -m venv "$HOME/venvs/gz-toolkit"
source "$HOME/venvs/gz-toolkit/bin/activate"
python -m pip install .
# Development:
python -m pip install -e ".[dev]"
# YAML plans for the bundled workflow scripts:
python -m pip install ".[workflow]"
```

Dependencies include NumPy, pandas, SciPy, pymatgen and the `myLAMMPS` Git
dependency. Git and network access are required to fetch it. On an offline
cluster, prepare compatible dependency wheels on a connected Linux machine,
install dependencies including `myLAMMPS`, then install the toolkit with
`python -m pip install --no-deps /path/to/GZ-toolkit`.

Check installation without launching calculations:

```bash
python -m pip check
gz-toolkit-potential-testing --help
python -m gz_toolkit.skills --help
python -c 'import gz_toolkit; print(gz_toolkit.__version__)'
```

LAMMPS, VASP, MPI, plugins and licensed VASP POTCAR files must be supplied
separately. Use the site's actual modules, partitions, accounts and resources.
The package does not supply simulation binaries or private potential files.

## Install the skill for Codex

This is optional and separate from Python package installation. Importing the
package also never installs a skill. Agents should first inspect their available
skills or run `python -m gz_toolkit.skills --check` (read-only; exit 1 if missing),
then suggest installation and ask the user if it is missing. Existing explicit
authorization to install does not need another confirmation.

Run on the computer hosting Codex, in the environment containing GZ-toolkit:

```bash
gz-toolkit-install-skill
# Equivalent:
python -m gz_toolkit.skills
```

The destination is `$CODEX_HOME/skills/materials-simulation-skill`, or
`~/.codex/skills/materials-simulation-skill` if CODEX_HOME is unset. The skill
becomes available on the next Codex turn. Invoke `$materials-simulation-skill`
to use it. Python package installation alone does not activate an agent skill.

Pass `--dest /path/to/agent/skills` to choose another host's skills directory.
The installer refuses to overwrite an existing skill unless you pass `--force`.
It copies all guides, scripts and templates without executing them.

Alternatively, copy the complete directory from a checkout:

```text
gz_toolkit/skills/materials-simulation-skill/
```

For a ChatGPT host offering skill uploads, create a zip with that single folder
at its root and use that host's upload interface. Installing a local Codex skill
does not upload it to a separate ChatGPT account or remote workspace.

## Using the package and skill together

For potential testing, use the existing `pot_inputs/` project and CLI. For a
custom campaign, save `sim_plan.yaml` or `sim_plan.json`, copy the skill's
workflow helpers to the campaign's `workflow/`, and prepare inspectable files.
Follow user-selected review checkpoints and reuse recorded decisions.
Preparation does not authorize submission: the agent needs human authorization
before launching the reviewed campaign. Degree 4's launcher includes downstream
runtime submissions, so review those resources and job counts too.

The [integration guide](gz_toolkit/skills/materials-simulation-skill/references/gz-toolkit.md)
documents remaining status, submission-failure, configuration and reproducibility
limits that the skill cannot repair in code.

Git checkouts and source distributions include examples, tests, this guide and
the changelog. Wheels contain the executable package, job/testing templates and
the complete simulation skill. Keep a checkout for example scripts.
