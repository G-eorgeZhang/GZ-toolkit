---
name: eamgen
description: Fit and iteratively validate EAM, Finnis-Sinclair or eam/he candidates with gz_toolkit.eamGen, using sourced references and explicit human guidance. Use for potential generation or modification, not for claiming that an unvalidated candidate is suitable for production.
---

# Guided EAM fitting

Use the standardized `python -m gz_toolkit.eamGen` tools and Python API.
Read [the module manual](references/manual.md) and inspect the project's
`project.json`, `fit_options.json`, dataset provenance, and latest revision report.
Use the configured Python environment. The module does not contain an LLM or
network client: the host agent does literature research and calls these tools.

## Establish the actual study

Record species, composition/temperature/pressure ranges, target phenomena,
reference methods, available potential files, and tolerances in `project.json`.
Respect the user's priorities and unusual physical hypotheses; turn each into
a testable target or proposed model change, explaining the tradeoff.
Do not interpret "generate a potential" as permission to launch unlimited
DFT/ML calculations. Agree on the reference job batch/resources before expensive
execution; prepare inputs and report the command, method and expected scope first.

## Baselines and literature

Search primary papers, their supplements, NIST IPR/OpenKIM model records and
official LAMMPS/VASP documentation. Do not invent numerical experimental targets.
For every datum record URL/DOI, table/figure/page, units, uncertainty, temperature,
composition, reference-energy convention and whether it is fit or validation data.
Put sources in `project.json.sources` and sourced property targets in `targets.json`.
Check licenses and preserve attribution when reusing potential files/functions.

Import with an explicit `--style`; never relabel eam/fs as eam/he or flatten
receiver/donor densities into alloy EAM. Check rho_min/rho_max, radial core,
cutoff continuity, isolated-atom energies and element ordering.
Compare imported baselines against LAMMPS before fitting. The internal evaluator
uses cubic splines and does not promise identical LAMMPS interpolation.

For Pd–He–Ni: preserve the supplied Pd/He HE functions as the starting baseline.
This study's user requires fitting and reference calculations on their HPC,
not on the local workstation. Prepare code/inputs locally; wait for the HPC
connection and execution details before fitting or launching reference jobs.
The Pd–Ni potential is primarily a reference calculator: its absolute energies
and density gauge may differ from the HE baseline. `compose` is only a documented
provisional initialization, never a completed ternary potential. Remove H only
by explicitly selecting Pd, He, Ni. Inspect every provisional function.
Do not assume a universal model includes He or is reliable on dilute noble-gas
defects, gas dimers, or compressed gas; confirm model species/domain support and
benchmark critical Ni–He environments against converged DFT when available.

## Reference dataset

Propose small, converged snapshots using `sample` or user structures. Cover
EOS/strains, alloy compositions/arrangements, displacements, vacancies, gas
tetra/octa sites, gas-vacancy clusters and relevant surfaces/GB environments.
The sampler only supplies FCC/BCC geometry proposals; it does not validate them.
Reserve separate groups/trajectories/compositions for validation and final tests.
Do not split adjacent snapshots of one trajectory across these groups.

Use `prepare-reference` to write LAMMPS atomic/metal run-0 inputs or VASP inputs
from the user's INCAR/KPOINTS/POTCAR. Verify POTCAR species order and DFT settings.
Run only the agreed explicit argv with `run-reference`, or prepare cluster jobs
using the existing jobs tooling and import after completion. Do not pass sbatch
to the synchronous runner. `import-reference`/`collect-references` captures
energies, forces, stress and source hashes; combine only compatible energy zeros.
Different calculator energies need documented per-element offsets and/or
consistent energy differences; never silently blend teacher/DFT energy scales.

## Fit, inspect, ask, revise

For the supplied Zhou model, named analytic forms from Eqs. 1/2/4/6/7 are available
through `analytic_terms.json` and `fit_analytic_model`. Use explicit parameter
bounds, signs and sourced/hypothesis annotations. Analytic and spline corrections
are separate revision modes; inspect which mode is active before fitting.

Start with a small selected set of pair/embedding functions. Keep density functions
frozen initially; opt into specific signed HE density corrections only when data
and the user's model hypothesis justify them. Pair keys are sorted (`Ni-Pd`,
`He-Ni`); density keys are directed (`Ni<-He`, `He<-Ni`).
Explain which functions are frozen, fitted or seeded. Preserve short-range core
and smooth radial corrections. If adding a new equation/function beyond these
splines, implement it explicitly rather than pretending coefficient changes
represent the new physics. Monitor embedding-domain violations and fit rank.

`fit --revision r001` snapshots inputs, coefficients, file hashes and reports.
Check optimizer convergence, train vs validation errors, density ranges, cutoff
behavior and regressions to the primary baseline. Show concrete results and
suggest the next data or parameter change; ask the user to choose when there
is a real physical tradeoff. Save their reasoning in `project.json.decisions`.
Do not record a human review unless the human actually supplied the decision.
Use `review` to record `approved_for_testing` or `rejected`, with reviewer/reason.
Create a new revision for each fit; keep the original baselines and prior outputs.

## Physical validation and stopping

After approval for testing, `handoff` creates a new potential_testing project.
Inspect its species, crystal structures, sizes, gases and HPC settings before
running its normal validation and submission workflow. HE in pure metals is
covered by that suite; ternary alloy gas environments, GBs, diffusion and bubble
behavior need additional explicit reference/production-validation cases.
Gas-complex formation energies in the existing summary assume zero gas chemical
potential; align to sourced DFT/experimental definitions before comparing them.
Interstitial-gas complex generation is disabled in the starter handoff; inspect
those geometries before enabling it.

Use `compare` with numeric metrics already converted to each target's stated
units/conditions. Missing or nonfinite results fail comparison. CSV elastic
constants are in GPa; dataset stress is in eV/Angstrom^3.
Keep the final test set untouched during iterations; use `evaluate --include-test`
once the model/options are selected. Report validation overfitting if the same
holdout has guided many revisions. Do not automatically iterate until thresholds
pass by changing targets/weights, hiding failures, or treating optimizer success
as a scientifically validated potential.

Deliver the candidate, fit report, sourced physical comparisons, intended domain,
remaining failures and next proposed experiment. Label it candidate until the
user accepts the physical evidence for their study.
