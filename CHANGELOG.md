# Changelog

## Unreleased

- `BCCDefect.add_111_loop` takes `n_layers=3` (default, perfect 1/2<111> loop)
  or `n_layers=1` (one extra or missing (111) plane, a faulted a/6<111> loop
  placed on the free stacking site). `add_loop` passes `n_layers` through.
- `edge_dislocation` / `screw_dislocation` take `atomic_masses`; when omitted,
  masses are looked up per element. Previously the written Masses section fell
  back to myLAMMPS' 4-element Fe/H/He/Cr table regardless of `ff_elements`.
- Fixed SIL loops with `n_repeats > 1`, which placed inserted atoms exactly on
  top of existing ones; the slab now splits as b/2 and T - b/2.
- Added opt-in extended benchmarks, guarded BCC C15 construction, SFT-seed
  testing, matched-reference line energies, loop size/type scans, alloy-gas
  binding references, molecular references, prescribed NEB paths and broader
  SEAKMC discovery/summary support. Added backend-specific GPU arguments,
  separate CPU preparation/summary headers and fail-fast execution scripts.
- Retained the LAMMPS ELASTIC method and added reviewed imported elastic cells
  with atomic/charge/full styles; no EOS fit replaces the stress/strain tensor.
- Added pure-element finite-temperature lattice-constant tests at 600–1023 K,
  optional alloy tests, Python-built supercells, shape-preserving NPT inputs,
  final-window convergence screening, input fingerprints, diagnostics, status
  and accepted thermal metrics separate from the existing 0 K constants.
- Bundled the materials-simulation skill, software guides, workflow scripts and
  templates, adapted for Codex and AI-human collaboration.
- Added `gz-toolkit-install-skill` and `python -m gz_toolkit.skills` to install the
  skill from a checkout or wheel, with explicit destination/overwrite options.
- Connected alloy-defect bulk/replica cases to all four parallelization degrees
  and status reporting through a shared work-group catalog.
- Removed obsolete potential-fitting exports, command, examples, tests and
  documentation. The toolkit builds structures and tests supplied potentials.
- Included tests and skill resources in source distributions and stopped
  ignoring the test suite in Git.
