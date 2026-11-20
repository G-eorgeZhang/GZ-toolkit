# MOOSE — sub-skill of materials-simulation-skill

MOOSE (Multiphysics Object-Oriented Simulation Environment, INL) is a finite-element framework:
phase field, heat conduction, solid mechanics, rate theory, and app-specific physics (e.g.
Marmot, Bison, user apps). Users run an **application executable** built on MOOSE, never
"moose" itself.

**Docs:** https://mooseframework.inl.gov — syntax index:
https://mooseframework.inl.gov/syntax/index.html (every block and object, with parameters).
Which objects exist depends on the **app and its modules** — ask which app the user runs, then
check available syntax with `<app>-opt --dump [search]` or `<app>-opt --show-input`. Do not
use an object without confirming the app has it.

## 1. Input file (HIT format, `.i`)

```
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmax = 100
  ymax = 100
[]

[Variables]
  [c]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = c
  []
  [dt]
    type = TimeDerivative
    variable = c
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = c
    boundary = left
    value = 1
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  dt = 1
  end_time = 100
[]

[Postprocessors]
  [avg_c]
    type = ElementAverageValue
    variable = c
  []
[]

[Outputs]
  exodus = true
  csv = true            # postprocessors -> <input>_out.csv
[]
```

Common top-level blocks: `Mesh`, `Variables`, `AuxVariables`, `Kernels`, `AuxKernels`, `BCs`,
`ICs`, `Materials`, `Functions`, `Postprocessors`, `VectorPostprocessors`, `Preconditioning`,
`Executioner`, `Outputs`, `GlobalParams`. HIT supports `${var}` substitution of top-level
variables defined as `var = value` at the top of the file.

## 2. Running

```bash
<app>-opt -i input.i
mpiexec -n 16 <app>-opt -i input.i
<app>-opt -i input.i Executioner/end_time=200 Mesh/nx=100     # command-line overrides
<app>-opt -i input.i --recover                                  # resume from checkpoint
```

Command-line overrides (`Block/param=value`, nested: `Kernels/diff/variable=c`) make MOOSE a
natural fit for option (b) of workflow step 3: keep one input and put
`Executioner/end_time={{TEND}}` etc. in the job body.

## 3. Workflow integration (steps 3–7)

- **Structure (step 2):** usually a mesh. Either generated in the input (`GeneratedMesh`), or a
  mesh / initial-condition file per run dir (`FileMesh` with `file = ...`, Exodus initial
  conditions, image-based ICs). Ask which.
- **Template candidates:** top-of-file HIT variables are the cleanest place:
  `T = {{T}}` then `${T}` inside blocks.
- **Job body:** `srun <app>-opt -i input.i > moose.out 2>&1`. Add
  `[Outputs] checkpoint = true` for long runs so `--recover` works on resubmission.
- **Checks:** ask the user; reliable choices are that `input_out.csv` exists and its last
  `time` reaches `end_time` (custom Python check, or `regex`/`last_line` on the CSV), and that
  `moose.out` does not contain `*** ERROR ***`. Confirm the exact end-of-run text on a finished
  test run before relying on a `contains` check.
- **Summary:** `csv_column` on `<input>_out.csv` (postprocessors), `which: last`; field data
  from Exodus via a custom Python field (e.g. `netCDF4`/`meshio`).
