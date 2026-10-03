# eamGen: guided EAM/FS/HE potential generation

`gz_toolkit.eamGen` provides a working candidate fitter and standardized tools
for agents and humans. It does not promise an accurate potential from a chemical
formula alone. Reference coverage, physical representability and independent
validation determine whether a candidate is useful.

The initial application is Pd–He–Ni based on a supplied Pd–He–H `eam/he` model.
The HE formulation is retained: density depends on both receiver and donor
species and can be negative. A Pd–Ni potential supplies additional reference
data; Ni–He requires suitable DFT/teacher data and explicit physical validation.

For this study, run fitting and reference calculations on the user's HPC only.
The supplied `examples/eamgen_pd_he_ni/fit_hpc.job` is a SLURM template with
explicit project/revision environment variables. Review the cluster's account,
partition, modules, Python environment and resource requests before submission.
It is not submitted automatically. Until HPC access is provided, the delivered
code and inputs are the output; no actual Pd–He–Ni fit is performed locally.

## Installation and working directory

Follow `INSTALL.md` in the repository for HPC installation. From the checkout,
activate the cluster Python environment and run `python -m pip install .`.
Check `gz-toolkit-eamgen --help`; it does not perform fitting. Examples below
assume the checkout is the working directory. Keep input potentials and generated
projects outside the source package, and use fresh output paths/revision names.
Python package installation does not install LAMMPS, VASP or their models.

## Implemented model and fit

The evaluator computes

```text
E = sum_i F_alpha(rho_i) + 1/2 sum_i,j,images phi_alpha,beta(r_ij)
rho_i = sum_j,images f_alpha,beta(r_ij)       (FS / HE)
rho_i = sum_j,images f_beta(r_ij)             (alloy EAM)
```

Cells have lattice vectors as rows. Periodic neighbors include all images inside
the cutoff, including self images; the reference evaluator is deliberately for
small fitting cells. Large/skewed neighbor enumerations are rejected.
Energies are eV, positions Angstrom, forces eV/Angstrom and stress
eV/Angstrom^3, positive in tension. Analytic forces and virial stress are fitted.

Models contain explicit embedding, density and pair tables. Existing potentials
can be imported/exported as `eam/alloy`, `eam/fs`, or `eam/he` setfl files.
No funcfl, CD-EAM or MEAM parser is implied. Setfl pairs store **r times phi**;
the importer preserves the original value at r=0. The internal evaluator
interpolates these products and density/embedding arrays with SciPy cubic
splines. LAMMPS uses its own interpolation, so exported-file parity must be
checked before physical testing. File-format support alone does not establish
that an imported potential behaves identically in this evaluator.

Fitting adds smooth cubic-spline corrections to selected functions. The
baseline core below `radial_min` and radial correction value/slope at the cutoff
are anchored. Embedding corrections preserve F(0), helping retain isolated
species reference energies. All baseline functions outside the selection stay
unchanged. Fixed-density fitting uses bounded linear least squares; optional
density corrections use nonlinear least squares and reject invalid density
domains. Regularization controls coefficient size/nonidentifiability; it is
not proof that a poorly constrained fit is physically correct.

The HE embedding domain may be negative. Choosing density corrections can
change host interactions even if the host pair/embedding tables are frozen.
Directly combining tables from different potentials does not align their gauge.
`compose_model` retains the primary HE grids/functions, imports covered secondary
functions as provisional, and labels missing/domain-incompatible functions as
generic seeds. A missing Ni–He interaction is never represented as fitted data.
Currently composition requires an HE primary; cross-format reference *data*
can always be imported through the reference calculators instead.

## Pd–He–Ni starting point

Run from the repository or an installed package, using your configured Python
environment. `python` below means that environment's interpreter.

```bash
python -m gz_toolkit.eamGen import-potential --file PdHHe.eam.fs.he.t --style eam/he --output pd_he_h.json
python -m gz_toolkit.eamGen import-potential --file PdNi.eam.alloy --style eam/alloy --output pd_ni.json
python -m gz_toolkit.eamGen compose --primary pd_he_h.json --secondary pd_ni.json --elements Pd He Ni --output pd_he_ni_seed.json
python -m gz_toolkit.eamGen init --root pd_he_ni_fit --elements examples/eamgen_pd_he_ni/elements.json --baseline pd_he_ni_seed.json
```

The second potential's filename/style above are examples: choose its actual
format after inspecting the supplied file and paper. Metadata lattice constants
in `elements.json` are geometry starting guesses, not sourced fitting targets.
Prefer the actual baseline's lattice constants and relaxed reference calculations.

The project includes baseline/model JSON, `dataset.json`, `fit_options.json`,
`targets.json`, `project.json`, references and revision folders. Inspect the
composition report's provisional functions before proceeding.
Copy/adapt `examples/eamgen_pd_he_ni/fit_options.json` to fit Ni embedding and
Ni–Ni/Pd–Ni/Ni–He pairs while initially retaining all density functions and
the Pd/He baseline. Add selected `density_keys`, e.g. `Ni<-He` and `He<-Ni`,
only when the reference data supports modifying them.

Generate geometry proposals with:

```bash
python -m gz_toolkit.eamGen sample --spec examples/eamgen_pd_he_ni/sampling.json --output snapshots.json
```

These snapshots have **no** energy/force targets. The example includes EOS,
random Pd–Ni, interstitial He in Ni and a held-out ternary vacancy-gas case.
It is a starting geometry pool, not enough data to validate the intended study.
Use independent arrangements/trajectories; add defects, surfaces, grain boundaries,
He–He configurations and gas pressure data appropriate to the user's study.

## Reference calculators

A LAMMPS preparation spec is a JSON object with:

```json
{
  "backend": "lammps",
  "configuration": {
    "name": "NiHe_tetra", "species": ["Ni", "He"],
    "positions": [[0, 0, 0], [0.88, 0.88, 0.88]],
    "cell": [[10, 0, 0], [0, 10, 0], [0, 0, 10]],
    "split": "train", "group": "NiHe_family",
    "provenance": {"source": "geometry example only; replace with real supercell"}
  },
  "type_map": {"Ni": 1, "He": 2},
  "masses": {"Ni": 58.6934, "He": 4.0026},
  "potential_lines": "pair_style eam/he\npair_coeff * * teacher.eam.he Ni He",
  "potential_files": ["teacher.eam.he"]
}
```

Replace its geometry and teacher with real inputs. Any LAMMPS potential can be
used if it supports **atomic** atom style and **metal** units; the adapter does
not install plugins/models. It rotates triclinic cells to LAMMPS's restricted
frame and converts forces/stress back on import. Periodic snapshots should use
converged cells for the teacher, independently of this fitter's image handling.

A VASP spec replaces backend-specific keys with `"backend": "vasp"`,
`"templates": "path/to/static_inputs"`, and optionally
`"species_order": ["Pd", "Ni", "He"]`. That directory must contain the user's
INCAR, KPOINTS and POTCAR. POTCAR must match the POSCAR order. No pseudopotentials
or VASP binary are supplied. Use static, converged single-point settings for
energy/force snapshots; import rejects unconverged or changed-geometry results.
The importer uses pymatgen `vasprun.xml` and defaults to `e_fr_energy` (free
energy, consistent with forces at finite smearing). Set `energy_kind: "sigma_zero"`
in the VASP preparation spec for `e_0_energy` when reproducing a documented
zero-smearing reference convention; reconcile energy/force consistency and do
not blend the two conventions silently. The importer converts
compression-positive kbar stress to this dataset's tension convention.

```bash
python -m gz_toolkit.eamGen prepare-reference --spec reference_spec.json --directory refs/NiHe01
python -m gz_toolkit.eamGen run-reference --directory refs/NiHe01 --argv-file lammps_command.json
python -m gz_toolkit.eamGen import-reference --directory refs/NiHe01 --output NiHe01.json
python -m gz_toolkit.eamGen collect-references --directories refs/NiHe01 refs/NiHe02 --energy-reference "Documented teacher model, settings and energy zero" --output labelled_dataset.json
```

`lammps_command.json` is an argv list such as `["lmp", "-in", "in.reference"]`.
VASP could use `["mpiexec", "-n", "4", "vasp_std"]`. `run-reference` is synchronous,
uses no shell and captures output/return code/input hashes. Use existing jobs
tools for SLURM preparation; submit explicitly and import after jobs finish.
The runner is not a scheduler monitor. External-run imports mark execution
provenance as unverified unless an execution record exists.

Before blending calculators, reconcile energy references. The dataset accepts
explicit `energy_offsets` (per-element eV subtracted from each target energy).
Offsets never arise from silently fitting away errors. Record magnetic states,
DFT functional/POTCAR/cutoff/k-mesh, teacher version and relative-energy convention.
Stress fitting assumes static references with no thermal kinetic contribution.

## Fit and human-guided iteration

### Analytic equations from the supplied paper

`analytic.py` implements Zhou et al. Eqs. 1/2 (repulsive pair/cutoff), 4
(signed exponential density), 6 (He core polynomial) and 7 (piecewise He
embedding). `examples/eamgen_pd_he_ni/zhou_2021_he_terms.json` transcribes the
FS parameters from Table II, preserving the table's **negative** beta_HeHe.
The original Pd/H functions remain imported tables. This example applies to
the original Pd–H–He model, including its H functions; remove H-related terms
when working on a Pd–He–Ni subset.

To fit a named equation instead of spline corrections, place `analytic_terms.json`
in the fit project's root. Each term declares component, key, form, parameters,
source/hypothesis and `bounds`. Only parameters in `bounds` are free; an empty
mapping freezes them. For a new Ni–He signed density, an explicit term might
use `key: "He<-Ni"`, `form: "zhou_density"`, a negative amplitude starting guess
and an amplitude interval entirely below zero. This is a hypothesis to fit,
not a claim that Ni must use Pd's parameter values.

`fit` detects this file and uses bounded nonlinear parameter fitting, saving
it with the revision. The public API is `fit_analytic_model`. Analytic and spline
parameter families are fitted in separate revisions; a previous candidate can
be the next baseline. Remove the analytic file to select spline mode again.

Use a dataset with `configurations`, `energy_reference`, `energy_offsets` and
`schema_version: 1`. Each configuration carries species, coordinates, cell,
optional energy/forces/stress, weight, group, split and provenance. Duplicate
names or related groups crossing splits are rejected.

After collecting reference results, explicitly populate the fitting project:

```bash
cp examples/eamgen_pd_he_ni/fit_options.json pd_he_ni_fit/fit_options.json
cp labelled_dataset.json pd_he_ni_fit/dataset.json
```

Inspect these files before submitting the HPC fit. Initialization leaves an
empty dataset and a placeholder energy convention; those are deliberately
insufficient for fitting. Edit `project.json` with the intended study, sources
and decisions. `targets.json` contains independent physical comparison targets;
these targets are not automatically converted into optimizer observations.

| Spline option | Meaning |
|---------------|---------|
| `embedding_keys`, `pair_keys` | `null` selects all functions; `[]` freezes the component. Pair keys use alphabetically sorted species, e.g. `Ni-Pd`. |
| `density_keys` | Empty by default; FS/HE keys use `receiver<-donor`, e.g. `He<-Ni`. Changes can alter host behavior. |
| `knots`, `radial_min` | Correction resolution and preserved short-range core. |
| `energy_scale`, `force_scale`, `stress_scale` | Positive residual weighting scales in eV/atom, eV/Angstrom and eV/Angstrom³; smaller values give greater weight. |
| `regularization`, `coefficient_bound` | Positive penalty and symmetric coefficient bounds controlling unconstrained corrections. |
| `max_parameters`, `max_evaluations` | Parameter ceiling and nonlinear optimizer evaluation limit; not a scheduler time limit. |

Configuration `weight` provides additional relative weighting. Analytic mode
uses explicit parameter bounds in its terms file instead of spline selections;
residual weighting and regularization still apply. Start with a small physically
motivated parameter selection rather than freeing every table.

Submit `examples/eamgen_pd_he_ni/fit_hpc.job` with the prepared absolute project
path and a fresh revision name, as described in `INSTALL.md`. The direct `fit`
command below is for an HPC allocation or an equivalent configured batch job.

```bash
python -m gz_toolkit.eamGen fit --root pd_he_ni_fit --revision r001
python -m gz_toolkit.eamGen review --revision pd_he_ni_fit/revisions/r001 --status approved_for_testing --reviewer "Zhang" --reason "Reviewed validation errors; proceed to physical regression tests"
python -m gz_toolkit.eamGen handoff --revision pd_he_ni_fit/revisions/r001 --root testing_r001 --pot-name PdHeNi_r001
```

The review command must represent an actual user's decision. It approves
**testing**, not production use. Each revision saves baseline/dataset/options
snapshots, input hashes, coefficients, model JSON, setfl export and an error
report. Failed optimizer convergence blocks approval for testing. A new revision
name is required for subsequent fits; `--baseline` can point to a prior model.
Change targets, weights or function selections only for an explicit recorded
physical reason. Dataset validation/test points never enter optimization.

The testing bridge creates ordinary `pot_inputs/`, `potentials/` and driver
scripts. He is a gas; Pd/Ni are FCC hosts. It enables elemental elastic tests,
binary alloy lattice tests, single-gas sites and small vacancy-gas complexes.
Interstitial-gas complexes are disabled in the starter config: inspect those
geometries before enabling the existing pipeline's generator. Ternary alloy gas
defects and grain-boundary/bubble dynamics need separate cases. No LAMMPS/SLURM
jobs are submitted by handoff.

Use the existing potential-testing `validate`, `check`, `run`, `status` and
`summarize` commands after reviewing HPC settings. Its gas formation-energy
summary assumes zero gas chemical potential; align the definition before
comparing DFT/experimental values. Extract desired CSV metrics into a JSON mapping
and convert to the target's units before `compare`:

```bash
python -m gz_toolkit.eamGen compare --metrics measured_metrics.json --targets pd_he_ni_fit/targets.json --output physical_comparison_r001.json
python -m gz_toolkit.eamGen evaluate --model pd_he_ni_fit/revisions/r001/model.json --dataset final_holdout.json --include-test --output final_errors.json
```

Each comparison target needs `metric`, `value`, positive `tolerance`, `units`,
`source` (DOI/URL) and `conditions` (temperature, composition, method). It must
come from a verified source, with uncertainty/citation location recorded as
extra fields if available. Missing metrics fail comparison. The comparison
tool does not automatically convert units or correct physical definitions.
Repeated validation-guided iterations can overfit the validation set: preserve
a genuinely independent final test set.

## Command reference and troubleshooting

Use `python -m gz_toolkit.eamGen COMMAND --help` for the exact arguments.

| Commands | Purpose |
|----------|---------|
| `import-potential`, `compose`, `init` | Import a known format, create an HE-primary provisional seed and initialize a project. |
| `sample` | Propose unlabelled geometries with reproducible families/splits. |
| `prepare-reference`, `run-reference` | Write calculator inputs and explicitly execute a synchronous command in an HPC allocation. |
| `import-reference`, `collect-references` | Parse finished calculations and assemble a documented labelled dataset. |
| `fit`, `evaluate` | Fit a fresh revision or evaluate a dataset without changing the model. |
| `review`, `handoff` | Record the user's decision and export an approved candidate for testing. |
| `compare` | Compare measured properties with sourced, unit-aligned targets. |

- **Existing output/revision:** choose a fresh path/name; the CLI generally
  refuses overwriting outputs. Review decisions retain their previous record.
- **Missing targets or placeholder convention:** import actual reference results
  and replace the project's dataset/convention before fitting.
- **Density outside the embedding domain:** inspect seeded densities and sampled
  environments; constrain changes or deliberately redesign the model/domain.
  Do not silently extrapolate a signed HE density.
- **Unconverged VASP or changed input:** rerun a consistent static calculation
  with appropriate convergence settings; retain its provenance.
- **Optimizer failure:** inspect conditioning, bounds, data coverage and selected
  parameters. A failed fit cannot be approved for testing.
- **Poor physical results despite small fit errors:** add relevant reference
  environments, reconsider the functional form or documented weights, and use
  a fresh revision. Keep the final test set independent.
- **LAMMPS mismatch:** verify pair style, species ordering, energy conventions,
  cutoff and table interpolation before attributing the mismatch to fitting.

For each revision retain the input snapshots, `report.json`, `coefficients.json`,
`model.json`, exported `candidate.eam.he` (or the corresponding FS/alloy suffix)
and `review.json`. Archive the testing configuration, calculator versions,
measured properties and comparison report alongside it. Passing numeric tests
does not establish suitability for grain boundaries, gas bubbles or diffusion;
validate the intended phenomena separately.

## Python interface

```python
from gz_toolkit.eamGen import EAMModel, Dataset, FitOptions, fit_model

baseline = EAMModel.load("pd_he_ni_seed.json")
dataset = Dataset.load("labelled_dataset.json")
result = fit_model(baseline, dataset, FitOptions(
    embedding_keys=["Ni"], pair_keys=["He-Ni", "Ni-Pd"], density_keys=[]))
result.model.write_setfl("candidate.eam.he")
print(result.report["metrics"])
```

`python -m gz_toolkit.eamGen --help` lists all standardized commands.
For the host-agent workflow, use `gz_toolkit/eamGen/skills/eamgen/SKILL.md`. This skill and manual
are distributed as package resources under `gz_toolkit/eamGen/skills/eamgen/`.
Copy that entire `eamgen` folder into your host agent's skill directory to install
it. It is shipped with the package, not silently installed globally.

## Sources and practical limits

- [LAMMPS EAM/FS/HE definitions and setfl layout](https://docs.lammps.org/pair_eam.html).
- [NIST Pd–H–He record and Zhou et al., PRB 103, 014108](https://www.ctcms.nist.gov/potentials/entry/2021--Zhou-X-W-Bartelt-N-C-Sills-R-B--Pd-H-He/).
- [VASP reference-training practices](https://vasp.at/wiki/Best_practices_for_machine-learned_force_fields).
- [VASP pressure convention](https://vasp.at/wiki/PSTRESS).

The first release supplies tabulated/spline candidates and the supplied paper's
He-related analytic equations. The original Pd/H functions stay tabulated;
other analytic forms and specialized physical property fitting need explicit
extensions. No self-contained AI agent,
web scraper, universal-model installer, DFT license, GPU provisioning or automatic
"fit until accurate" loop is included. Research and decisions are performed by
the host agent with the accompanying skill and the user.
Synthetic fitting and finite-difference tests verify the implementation; they
do not establish a scientifically accurate Pd–He–Ni potential. Actual baseline
files, papers, reference calculations and validation thresholds are still needed.
