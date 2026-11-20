# GZ-toolkit potential testing

Use the installed package (`import gz_toolkit`), its actual API, and the existing
project when the user is testing interatomic potentials. The package prepares
structures and LAMMPS inputs and aggregates physical metrics; it does not fit
potentials. Read the project's README and configuration before changing it.

## Potential names, settings and the pot_infobank

Potential folders are named `<name>_<system>_<year>` (`gao_FeHe_2011`), or
`<name>_<system>_` without a year; short names like `gao2011` are renamed
automatically once their `pot_inputs` JSON exists. The project root's
`settings.json` key `potinfo_savepath` (`false`, `true` = `<root>/pot_infobank/`,
or a path = `<path>/pot_infobank/`) decides whether `summarize-pot` saves each
tested potential as `pot_infobank/<pot_name>/` with `pot_info.json`, the
potential files, `pot.inc` and `README.md`.

When `check` reports missing reference info (paper and DOI; plus GitHub and
Matbench Discovery link for `kind: "universal"`), ask the user for it and add
it under `potential.reference` in that `pot_inputs` JSON before launch. Do not
invent citations.

When the user names a tested potential for a new simulation, find its folder
with `gz_toolkit.pot_infobank.list_potentials` / `pot_dir` (short names work),
read `pot_info.json` for `lc`, `type_map`, masses and crystal structure, copy
the potential files next to the input and `include pot.inc`. Ask the user
where the pot_infobank is if it is not in the project's `settings.json`.

## Preparation, review, and launch

The potential-testing workflow stores one configuration per potential in
`pot_inputs/`. Keep user-supplied potential files in `potentials/<pot_name>/`.
`python -m gz_toolkit.potential_testing.cli run` prepares reference inputs and
job files; it returns a `submit_all.sh` launcher without submitting jobs.

Before preparing a campaign, run `validate` and `check` with the same
`--pot-inputs`/`--run-dir` locations used by the project. These are separate CLI
commands: `run` does not automatically enforce validation. Review the potential
species/type mapping, reference settings, generated inputs, resources, modules,
job counts, and dependencies. Resolve every reported error before launch.

Only execute the launcher when the human has authorized that campaign's
submission. When they already authorized submission and its settings, reuse
that authorization rather than asking again. `scatter-cases` is an execution
command that calls `sbatch`; it is not a preview. For degree 4, approving the
launcher also needs to cover the downstream build job's runtime submissions.
Do not invoke `scatter-cases` separately for cases already launched by it.

Parallelization degrees in this package have fixed meanings:

- 1: one job for all potentials, using the first potential's HPC header.
- 2: one pipeline job per potential.
- 3: reference/build and work-group jobs with dependencies.
- 4: reference/build jobs followed by runtime case submissions, with batching
  when a group exceeds `hpc.max_jobs_per_batch`.

These degrees are different from `make_jobs.py --level`, which selects a custom
campaign's directory-tree depth. Do not interchange them. Alloy-defect cases
are under `point_defects_alloy/<system>/<composition>/repN/<case>/`; each replica
includes `bulk` plus the requested defects. They use `in.relax.lammps`.

## Completion and remaining code limits

### Finite-temperature lattice constants

Pure-element cases now include `thermal_lc/<element>/<temperature>K/` by default:
600, 700, 800, 900, 1000, and 1023 K. Python builds 10 x 10 x 10 supercells,
and `in.thermal.lammps` runs 50,000 NPT steps with a 1 fs timestep, zero external
pressure and isotropic shape-preserving expansion. The last 20,000 steps are
screened for flat length history and temperature consistency. Include these
extra jobs in resource reviews. Set `include_thermal_lc` false to disable;
alloys require `thermal_lc_for_alloys` and configured compositions.

Thermal status is stricter than the legacy classifier below: matching input
fingerprints, complete samples and a clean finish are required. `unconverged`
means the run finished but the averaging window failed its configured screen;
it contributes no accepted lattice constant and is not automatically extended.
Read `thermal_lc.dat` and `thermal_result.json`, inspect the final structure and
phase, and let the human decide whether to lengthen/relaunch. Do not select a
shorter apparently flat tail to make the test pass. The potential tag JSON
keeps thermal values in a separate `thermal_lc` table; `lc` remains the 0 K
lookup. Current potential-testing geometry supports BCC/FCC/HCP only.

### Extended benchmarks and migration

The package now has opt-in benchmark_cases for surfaces, C15 and SFT seeds,
dislocation lines, molecular references, and alloy gas complexes. Read
examples/extended_potential_testing.md in a checkout for schemas, energy
reference rules and limitations. Surfaces can also be enabled using
include_surfaces; loop_radius_factors and loop_types configure BCC size/type
scans. C15 is a guarded BCC seed constructor, and SFT construction is still an
unrelaxed triangular FCC vacancy platelet. Validate the relaxed topology.

SEAKMC case patterns can select extended defects with an explicit active volume.
Its runtime driver requires the completed parent relaxed.data. NEB cases require
reviewed ID/species-matched relaxed endpoints and a replica-aware launch. These
methods are complementary; an unknown SEAKMC destination is not an automatic
NEB endpoint. Check MPI, LAMMPS REPLICA and force-evaluator compatibility.

Python preparation/summary is CPU work, including when run inside a GPU job.
Set backend-specific lammps_args and run_command; no generic GPU suffix is
assumed for MTP/universal potentials. Degrees 3/4 support CPU preparation/summary
profiles via preprocess_partition_key/preprocess_ntasks. GPU SEAKMC needs a
tested template; GPU NEB needs an explicit replica/GPU mapping command.

### Legacy cases

Use `status` as a convenience view, then verify the expected terminal markers,
outputs, and scientific convergence before interpreting or promoting results.
For LAMMPS, require a clean exit (`Total wall time`), absence of errors, expected
RESULT fields, and the required final structure. Record convergence checks
appropriate to the calculation. For custom campaigns, encode these checks in
the saved plan and use `summarize.py --only-finished`.

The current package has limits that instructions cannot repair:

- `classify_case()` treats any RESULT marker as done, even in an unfinished log;
  SEAKMC status also relies on output existence rather than a validated finish.
- Generated jobs now stop on command failures, and failed sbatch calls raise.
  If sbatch is absent, files remain for manual launch. Preserve launcher output
  and job IDs; a summary file alone is not evidence that every case ran.
- Unknown configuration keys are dropped. Cross-check names against the actual
  dataclasses; a clean `validate` result cannot detect a discarded typo.
- Alloy builders use fresh random configurations unless explicitly seeded, and
  workflow generation does not currently expose a base seed. The generic skill
  scripts' deterministic SEED does not change this behavior. Save structures and
  input hashes for reproducibility; matched bulk/defect configurations and a
  pinned dependency revision still require package changes or explicit helpers.

If these limits affect the user's campaign, explain the concrete consequence
and propose a targeted code fix. Do not silently patch all package behaviors,
claim a skill fixes them, or launch a replacement campaign. A missing or failed
result stays missing/failed and is reported as such.
