# Pd–He–Ni eamGen starter files

Follow [the full manual](../eamgen_manual.md) and
[HPC installation guide](../../INSTALL.md). Run commands from the toolkit
checkout so the manual's `examples/...` paths resolve.

| File | Purpose and required review |
|------|-----------------------------|
| `elements.json` | Species metadata and initial geometry guesses; replace guesses with appropriate reference values. |
| `fit_options.json` | Ni embedding and Ni–Ni/Pd–Ni/Ni–He spline corrections; density functions initially frozen. |
| `sampling.json` | Reproducible geometry proposals with grouped splits; contains no energy/force labels. |
| `zhou_2021_he_terms.json` | Paper's original Pd–H–He analytic parameters; all bounds empty/frozen. Remove H terms and explicitly define bounded Ni terms before a ternary analytic fit. |
| `literature_observations.json` | Paper values with method/definition provenance; not ready-made comparison targets or a fitting dataset. |
| `fit_hpc.job` | Editable SLURM fitting template; requires a prepared project and fresh revision name. |

The Pd–He–H baseline and Pd–Ni teacher files are supplied separately. Composition
produces a provisional seed; it does not establish compatible gauges or fitted
Ni–He interactions. Review all imported/seeded functions, generate reference
labels on the HPC, and keep independent validation/test configurations.

After initialization, copy `fit_options.json` to the project's file of the same
name and replace `dataset.json` with the labelled dataset. Keep the analytic
parameter file outside the project unless intentionally selecting analytic
mode. An `analytic_terms.json` in the project root switches `fit` to that mode.

No command or job is submitted automatically by these files. Reference jobs,
fitting jobs and reviewed physical testing are separate steps. Fitting and
reference calculations for this study must run on the HPC.
