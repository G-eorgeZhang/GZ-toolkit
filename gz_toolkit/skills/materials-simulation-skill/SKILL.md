---
name: materials-simulation-skill
description: "Prepare, troubleshoot, and summarize materials simulation campaigns with GZ-toolkit, LAMMPS, SEAKMC, VASP, SPPARKS, or MOOSE. Use for structure construction, simulation inputs, parameter sweeps, reviewed SLURM jobs, completion checks, result extraction, and HPC setup. Supports AI-human collaboration through saved plans and inspectable intermediate files."
---

# Materials Simulation Workflow

One skill for the whole life of a simulation campaign, adapted from the user's Claude
skill for Codex and bundled with GZ-toolkit. Software guides and environment references
are loaded only when needed. Existing user instructions and configured environments
take precedence; this skill does not replace other installed skills.

```
 0 environment ─▶ 1 dir tree ─▶ 2 structures ─▶ 3 input script ─▶ 4 set up runs
     ─▶ 5 jobs (parallelization level) ─▶ 6 checker ─▶ 7 summarize ─▶ 8 analysis (optional)
                       ▲ user checks ▲             ▲ user checks ▲
```

**Find where the user is.** Someone asking "why does my LAMMPS job say Lost atoms" needs only the
LAMMPS sub-skill, not the whole workflow. Someone saying "I want to study X across T and Y"
starts at step 0. Someone with a finished tree asking for a CSV jumps to step 7. Read
`sim_plan.yaml` in the connected root first if it exists — it records what has been decided.

## Sub-skills and references — read the one you need, when you need it

| Need | File |
|------|------|
| LAMMPS / MD theory, input scripts, potentials, errors, compiling | `software/lammps/GUIDE.md` |
| SEAKMC input.yaml, AV sizing, troubleshooting | `software/seakmc/GUIDE.md` (+ `references/parameters.md`, `examples.md`) |
| VASP / DFT: INCAR, KPOINTS, POTCAR, POSCAR, calc types | `software/vasp/GUIDE.md` (+ `references/calc-types.md`) |
| SPPARKS (lattice KMC / Potts) | `software/spparks/GUIDE.md` |
| MOOSE (FEM / phase field, HIT input files) | `software/moose/GUIDE.md` |
| Any other code (GROMACS, QE, CP2K, ...) | `software/other.md` |
| Step 0 environment, known machines | `references/environment.md` |
| Step 2 structures with GZ_toolkit | `references/structures.md` |
| Step 5 SLURM details, arrays, dependencies, failures | `references/slurm.md` |
| Step 7 extraction recipes, OVITO | `references/summarize.md` |
| Installing/compiling codes on a cluster, module files | `references/hpc-install.md` |
| GZ-toolkit potential-testing CLI, review/submission boundaries, reliability limits | `references/gz-toolkit.md` |

Workflow scripts are in `scripts/` (standard library only; PyYAML if the plan is YAML).
Templates are in `templates/`.

## Ground rules (apply to every step)

1. **One change per step.** Every structure, every input file derives from a checked reference
   by exactly one modification. Never build a final structure from scratch in one go.
2. **Review checkpoints.** Make the directory tree, reference structures, final structures,
   input changes, and job counts inspectable. Honor the user's requested review checkpoints.
   When their existing instructions authorize proceeding, record that authorization and
   continue. Otherwise obtain the needed review before proceeding with dependent work.
   Generate and check all job files before asking for submission approval; a request to
   prepare jobs alone does not authorize `sbatch` or launching simulations.
3. **The user's files win.** Potentials, job files, reference inputs come from `reffiles/`. If
   something is missing, ask for it. Do **not** download potentials or job templates from the
   internet unless the user says so. Documentation and literature lookups for *how to write*
   input are encouraged — cite the doc page.
4. **Nothing is guessed about the machine.** Paths, partitions, accounts, modules differ between
   machines. Get them from `sim_plan.yaml`, `references/environment.md`, or the user.
5. **`sim_plan.yaml` is the memory.** Record every decision there (template:
   `templates/sim_plan.yaml`). Later steps and later chats read it instead of re-asking.
6. Batch missing information into one round. Use the host's available question tool for
   discrete choices, or ask in chat. Reuse choices and authorization already recorded.

## Choosing the workflow

For a GZ-toolkit potential-testing project (`pot_inputs/`), read
`references/gz-toolkit.md` and use the package's existing CLI and directory layout.
For a custom campaign driven by `sim_plan.yaml` or `sim_plan.json`, use the scripts
below. Do not create a second directory tree or independently submit the same cases
with both workflows. Resolve script/template paths relative to this SKILL.md.

---

## Step 0 — Environment

Read `references/environment.md`. In short:

- **Where does the work run?** Local machine or SSH to a cluster? Look at the context (device
  bridge folder, an ssh session, Windows paths). If not sure, **ask**.
- **Python**: which interpreter / `.venv` / conda env runs the workflow scripts on *that*
  machine; which local packages (GZ_toolkit, myLAMMPS, SEAKMC_py ...) must go on `sys.path`.
  Check `references/environment.md` "Known machines" first; if the machine is not listed, ask
  the user for the paths, and offer to add the machine to that list.
- Verify before relying on it: run `<python> -c "import sys; print(sys.version)"` and try the
  imports. Write the result into `sim_plan.yaml: environment:`.
- Copy `scripts/*.py` into `ROOT/workflow/` (the connected folder is ROOT). All workflow
  commands below are run from ROOT with the environment's interpreter.

## Step 1 — Directory tree

Ask the user for the tree: which variables, in which nesting order, their values, the directory
naming, and how many runs (and what a run means: statistical repeat with a different seed, or a
different position/configuration). Typical layout:

```
ROOT/                          level 0  (the folder the user connected)
  sim_plan.yaml  manifest.csv  workflow/
  reffiles/                    reference + unchanging files (potentials, ref structures, ref input, .job)
  T300K/                       level 1  variable 1
    rho1e14/                   level 2  variable 2
      He2V1/                   level 3  variable 3
        run1/ run2/ run3/ ...  level 4  run dirs — simulations happen here
  T400K/ ...
```

Put the answers in `sim_plan.yaml: tree:`, then:

```bash
python workflow/make_tree.py --dry-run     # show the tree to the user  ── GATE
python workflow/make_tree.py               # create dirs + manifest.csv
```

`manifest.csv` (one row per run dir with every variable's value) is what all later scripts
iterate over. If a level's values depend on its parent, use the dict form of `values`
(see the script docstring); for anything stranger, write the manifest yourself.

## Step 2 — Structures (Python + GZ_toolkit)

Read `references/structures.md` before writing any structure code. Essentials:

- Always use **GZ_toolkit** (ask for its path if `sim_plan.yaml`/known machines lacks it), and
  read its actual API before calling it — never invent function names.
- **Stage 0: perfect crystal** into `reffiles/ref0_*`. Verify lattice constant (for the
  potential being used), orientation, box size, atom count. **GATE.**
- **Stage 1..k: one change each** (e.g. add a screw dislocation) into `reffiles/ref1_*` ... Each
  stage reads the previous file. **GATE after each.**
- **Final stage**: the per-run change (e.g. a HeV cluster at position i around the dislocation),
  written **directly into each run dir** under the name in `structures.file`. Parallelize over
  run dirs with `--workers` for big trees.
- Start from `templates/make_structures_template.py` → `ROOT/workflow/make_structures.py`.
  Record each stage in `sim_plan.yaml: structures.chain` and mark `checked: true` when the user
  confirms.
- When all final structures exist: report a sample (atoms, types, box), suggest opening a few in
  OVITO, and **GATE** before step 3.

Worked example (He–V binding around a screw dislocation in BCC Fe): ref0 = perfect BCC Fe with
x=[11-2], y=[111], z=[1-10] → user checks → ref1 = ref0 + ½[111] screw dislocation along y →
user checks → final = ref1 + He_nV_m cluster at site i (one site per run dir) → user checks.

## Step 3 — Reference input script

Open the sub-skill for the code (table above) and write — or take the user's — **one** reference
input in `reffiles/`. Check syntax against the official docs for anything not certain
(LAMMPS: docs.lammps.org, VASP: vasp.at/wiki, SPPARKS: spparks.github.io, MOOSE: mooseframework.inl.gov).

Then turn it into a template: replace each value that varies across the tree with a placeholder
named after the variable — `{{T}}`, `{{cluster}}` — plus built-ins `{{SEED}}` (deterministic
per run), `{{RUN}}`, `{{RUN_INDEX}}`, `{{RELPATH}}`, `{{JOBNAME}}`, and any `constants:` from the
plan. Show the user the diff (reference vs template) — **GATE**. Unfilled placeholders are a
hard error later, so typos surface early.

Alternative when the code supports command-line variables (LAMMPS `-var T 300`, SPPARKS `-var`,
MOOSE `Executioner/end_time=10`): keep the input untouched and put `{{T}}` in the job body
instead. Choose whichever the user finds easier to read.

## Step 4 — Set up every run dir

Everything a run needs comes from `reffiles/`: for LAMMPS that is the structure (already in the
run dir from step 2), the input script, the potential file(s), and the `.job` file. If any is
missing from `reffiles/`, ask the user — do not fetch it online unless told to. Fill
`sim_plan.yaml: reffiles:` (`copy`, `link`, `templates`, `job`, `require`), then:

```bash
python workflow/setup_runs.py --dry-run
python workflow/setup_runs.py --workers 32      # parallel over cores
```

For thousands of run dirs, run the setup itself as a one-node SLURM job so each core handles
a slice of run dirs (`--workers $SLURM_CPUS_ON_NODE`). Already-submitted run dirs are skipped
unless `--force`.

## Step 5 — Jobs and parallelization level

**Parallelization level** = the tree level at which one SLURM job is created.

| Level | One job per | Inside the job | Called |
|-------|-------------|----------------|--------|
| run level (e.g. 4) | run dir | the reference `.job` rendered in the run dir | scatter submission (most nodes) |
| e.g. 3 | variable-3 dir | `cd run1; srun ...; cd ..; cd run2; srun ...` in series | level-3 |
| ... | ... | all runs below it, in series | ... |
| 0 | root | every run in the tree | single job |

`--pack k` runs k run dirs concurrently inside an aggregated job (each srun must request its
share of cores). Walltime of aggregated jobs scales from `job.time_per_run`.

The reference job (`reffiles/<name>.job`) supplies the header (`#SBATCH`, modules) and the per-run
body; mark the split with a `# @RUN` line (see `templates/run.job.example`). Before writing
it, confirm partition, account, nodes/cores, walltime, modules — cluster-specific values are
never guessed (`references/slurm.md`).

If the user asks to run and did **not** give a level, **ask** for it (show the job count each
level would produce: count of distinct dirs at that level in `manifest.csv`). Then:

```bash
python workflow/make_jobs.py --level 3            # writes job files + jobs_lv3.txt
python workflow/submit.py --level 3 --dry-run     # show the user  ── GATE
python workflow/submit.py --level 3               # sbatch from each job's dir; IDs -> jobs_submitted.csv
```

Resubmitting only what failed: `make_jobs.py --level L --only-unfinished`, then submit again.

## Step 6 — Checker

Ask the user what "finished" means — e.g. a file exists, `log.lammps` contains
`Total wall time`, a dump has ≥ N lines, no `ERROR`, ≥ N `KMC_*.dat` files. And whether they
want `.sh` or `.py`. Write them into `sim_plan.yaml: checks:` (types in `scripts/checklib.py`), then:

```bash
python workflow/make_checker.py --lang sh    # -> ROOT/checker.sh   (pure bash)
python workflow/make_checker.py --lang py    # -> ROOT/checker.py   (+ --by LEVEL tallies)
./checker.sh        ./checker.sh -u          # all runs / only unfinished; writes status.csv
```

Suggested defaults per code are in each sub-skill's "Workflow integration" section. If nobody
knows what a finished output looks like, look at one finished test run before writing checks.

## Step 7 — Summarize

When the user asks for results, ask which values from which files, and how to reduce them
(last / min / max / mean). Write `sim_plan.yaml: summary.fields`, then:

```bash
python workflow/summarize.py --only-finished --average   # summary.csv + summary_avg.csv
```

CSV by default. Field types: regex, LAMMPS thermo column, CSV column, VASP energy, line/file
counts, OVITO, or custom Python — see `references/summarize.md` and
`templates/ovito_extract_example.py`. Empty cells mean extraction failed; the warnings say why.

## Step 8 — Analysis (optional, on request)

Plots and fits from `summary.csv` (binding energy vs. position, Arrhenius fits of rates vs. T,
averages ± std over runs). Write the analysis as a script in `ROOT/analysis/` so it can be re-run
when more runs finish.

---

## Talking to the user

- Name the step you're on and what the next gate is.
- Show trees, diffs, job counts, and sample structures — not walls of generated files.
- After a gate, record the decision in `sim_plan.yaml` before moving on.
