# Step 2 — Building structures

Structures are built with Python and **GZ_toolkit**. The output of this step is one structure
file in every run dir, each derived from a chain of user-checked references.

## The reference chain

```
reffiles/ref0_perfect.*      perfect crystal                       ── user checks
reffiles/ref1_<change>.*     ref0 + ONE change (e.g. dislocation)  ── user checks
reffiles/ref2_<change>.*     ref1 + ONE change (optional)          ── user checks
<run dir>/<structures.file>  last ref + ONE per-run change         ── user checks a sample
```

Why: when a final structure is wrong, the chain tells you which single step broke it, and each
reference is reusable (the same ref1 can feed several campaigns).

What counts as "one change": insert one dislocation; add one cluster; apply one strain; swap
species on one sublattice; delete one region. Relaxing a structure with the potential is also a
stage of its own (and its output must be relaxed with **the same potential** that the production
runs use).

## Before writing code

1. Get the GZ_toolkit path (sim_plan / known machines / ask). Add it to `sys.path`.
2. Read the toolkit's API: list its modules and public functions/classes, read docstrings of
   the ones you will use. If the toolkit cannot do something, tell the user and propose a small
   helper — don't silently switch to another library.
3. Get the physics inputs from the user — do not default them: element(s), lattice constant
   (ideally the relaxed value for *their* potential), crystal structure, orientation (x, y, z
   Miller directions — check they are mutually orthogonal and right-handed), box repeats or
   target size, boundary conditions, output format (LAMMPS data `atom_style`, POSCAR, ...).

## Stage 0 — perfect crystal: what to verify and show

- Atom count = n_basis × repeats (BCC 2/cell, FCC 4/cell for cubic cells; for rotated cells,
  count = volume × ρ).
- Box lengths are integer multiples of the orientation's repeat distances (e.g. BCC with
  x=[11-2], y=[111], z=[1-10]: repeats a₀√6, a₀√3/2, a₀√2) — otherwise periodic images don't
  match and the "perfect" crystal has seams.
- Nearest-neighbour distance and coordination equal the ideal values everywhere (BCC: √3/2 a₀,
  8 neighbours; FCC: a₀/√2, 12) — a quick neighbour check or OVITO CNA showing 100 % BCC/FCC.
- Atom types / species order match what the potential's `pair_coeff` (or POTCAR) will expect.

Show the user: these numbers, and suggest opening the file in OVITO. Wait for confirmation.

## Defect stages: typical checks

| Change | What to verify |
|--------|----------------|
| Dislocation (screw/edge) | Burgers vector, line direction, position; DXA in OVITO finds one dislocation of the intended type; boundary handling (fixed layers, periodic arrays, image dislocations) matches what the user wants |
| Vacancy / He–V cluster | Atom count change = expected; position relative to the reference feature (distance from the dislocation core, angle); no overlap (min interatomic distance not below ~0.7 × NN) |
| Interstitial / solute | Site type (octahedral, tetrahedral, dumbbell orientation); min distance |
| Strain | New box = old box × (1+ε); atoms mapped affinely |
| Surface / void | Vacuum thickness, cut plane, no dangling partial cells |

## Final stage — into the run dirs

- The mapping from tree variables to the per-run change lives in code, e.g. `cluster` ←
  level-3 value, `site` ← list indexed by `run_index`. Write the list of sites (coordinates or
  shell/angle labels) to `reffiles/sites.csv` so the user can see and edit it.
- Write directly into `ROOT/<relpath>/<structures.file>`. Use `--workers` to parallelize over
  run dirs; for very large trees submit the generation as a one-node job.
- Report: how many written, a sample of 2–3 (atoms/types/box), and any run dir that failed.
  Ask the user to check a sample (e.g. open run1, a middle one, the last) — **GATE**.

## Template

`templates/make_structures_template.py` — copy to `ROOT/workflow/make_structures.py`. It
handles plan loading, `sys.path`, staging (`--stage 0|1|final`), parallel final stage and a
LAMMPS-data header report. Fill the `TODO(gz)` blocks with real GZ_toolkit calls.
