# Changelog

## 0.2.0

- Added `eamGen`: setfl import/export for `eam/alloy`, `eam/fs` and `eam/he`,
  energy/force/stress evaluation, selected spline corrections and the supplied
  Zhou He analytic forms.
- Added provisional composition with an HE primary, geometry proposals,
  LAMMPS/VASP reference adapters, immutable fit revisions, human review records
  and candidate handoff to potential testing.
- Added the `gz-toolkit-eamgen` command, packaged agent workflow, complete
  example manual and Pd–He–Ni starter inputs/SLURM template.
- Added HPC installation instructions and source-distribution documentation.

The initial study preserves the Pd–He `eam/he` formulation. Broader analytic
families and other EAM formats remain future extensions. This release contains
tools and templates, not a fitted or validated Pd–He–Ni potential.
