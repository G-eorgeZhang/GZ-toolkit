# LAMMPS / MD — sub-skill of materials-simulation-skill

## Core Philosophy

Always refer back to the official documentation: **https://docs.lammps.org/Manual.html**

When answering, cite the relevant doc section (e.g., `fix nvt` → https://docs.lammps.org/fix_nh.html).
When uncertain about a command's exact syntax, fetch the doc page before answering.

Key doc entry points:
- All commands: https://docs.lammps.org/Commands_all.html
- Fix styles: https://docs.lammps.org/fixes.html
- Pair styles (potentials): https://docs.lammps.org/pairs.html
- Compute styles: https://docs.lammps.org/computes.html
- Dump styles: https://docs.lammps.org/dumps.html
- Build/compile: https://docs.lammps.org/Build.html
- Packages with extra build steps: https://docs.lammps.org/Build_extras.html
- Errors: https://docs.lammps.org/Errors.html

---

## 1. MD Theory — Math and Materials Science Perspective

When the user asks "what is MD" or wants a conceptual explanation, cover these layers:

### 1a. Statistical Mechanics Foundation
- MD samples the **phase space** (positions **r**, momenta **p**) of an N-body system
- The system evolves via Newton's equations: **F** = m**a** = −∇U(**r**)
- **Ergodic hypothesis**: time averages ≈ ensemble averages (assuming the trajectory explores all accessible microstates)
- Observable A: ⟨A⟩ = (1/T)∫A(**r**(t), **p**(t))dt ≈ ⟨A⟩_ensemble

### 1b. Key Ensembles
| Ensemble | Fixed | Thermostat/Barostat |
|----------|-------|---------------------|
| NVE (microcanonical) | N, V, E | None — energy conserved |
| NVT (canonical) | N, V, T | Nosé–Hoover, Langevin, Berendsen |
| NPT (isothermal-isobaric) | N, P, T | Nosé–Hoover chains + Parrinello-Rahman |
| NσT (stress-controlled) | N, σ, T | Used for deformation studies |

### 1c. Integrators
- **Velocity Verlet**: r(t+dt) = r(t) + v(t)dt + (1/2)a(t)dt²; symplectic, time-reversible
- **Langevin**: adds stochastic force + damping → couples to heat bath, gives NVT
- **Nosé–Hoover**: deterministic extended-Lagrangian thermostat, correct canonical ensemble

### 1d. Connection to Monte Carlo
- MD generates a **deterministic** trajectory; MC generates **stochastic** configurations
- Both sample the same ensemble at equilibrium (for ergodic systems)
- MD gives dynamics (diffusion, phonons); MC cannot give dynamics directly
- **Kinetic Monte Carlo (KMC)**: used for rare-event dynamics (see `software/seakmc/GUIDE.md`)

### 1e. Probability and Distributions
- Velocity distribution at equilibrium: Maxwell–Boltzmann P(v) ∝ exp(−mv²/2kT)
- Configurational probability: P(**r**) ∝ exp(−U(**r**)/kT) (Boltzmann weight)
- **Radial Distribution Function g(r)**: probability of finding an atom at distance r relative to an ideal gas; g(r)→1 at large r
- **Mean Squared Displacement (MSD)**: ⟨|r(t)−r(0)|²⟩ = 6Dt for 3D diffusion

---

## 2. LAMMPS Input File Structure

A LAMMPS input script is executed **top to bottom**, line by line. The canonical order:

```lammps
# ── 1. Initialization ──────────────────────────────────
units           metal          # metal: eV, Angstrom, ps, K, bar
dimension       3
boundary        p p p          # p=periodic, s=shrink-wrap, f=fixed
atom_style      atomic         # or charge, full, molecular...

# ── 2. Atom Definition ─────────────────────────────────
read_data       structure.lammps   # OR lattice + create_box + create_atoms

# ── 3. Interatomic Potential ───────────────────────────
pair_style      eam/alloy
pair_coeff      * * Fe_mm.eam.alloy Fe

# ── 4. Settings ────────────────────────────────────────
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes

# ── 5. Computes and Variables ──────────────────────────
compute         myKE all ke
variable        myTemp equal temp

# ── 6. Output ──────────────────────────────────────────
thermo          100
thermo_style    custom step temp pe ke etotal press vol
dump            1 all custom 1000 dump.*.lammpstrj id type x y z vx vy vz

# ── 7. Run ─────────────────────────────────────────────
fix             1 all nvt temp 300 300 0.1
timestep        0.001          # ps for metal units
run             100000
```

**Critical ordering rules** (common source of bugs):
- `pair_style` must come before `pair_coeff`
- `read_data` / `create_atoms` must come before any fix/compute referencing atoms
- `fix` commands must come before `run`
- Groups must be defined before use

### 2a. Units Systems
| units | Energy | Length | Time | Force | Pressure |
|-------|--------|--------|------|-------|----------|
| metal | eV | Å | ps | eV/Å | bar |
| real | kcal/mol | Å | fs | kcal/mol/Å | atm |
| si | J | m | s | N | Pa |
| lj | ε | σ | τ | ε/σ | ε/σ³ |

⚠️ **Never mix units systems.** The `units` command must be the first substantive command.

---

## 3. Guided Input File Generation (Interview Protocol)

When the user wants to create an input file from scratch (or says "I want to simulate X"), use this **question funnel**. Ask in stages — don't dump all questions at once.

### Stage 1: Goal and System
Ask the user:
1. What is the **material** (element(s), crystal structure)?
2. What is the **physical question** — e.g., diffusion coefficient, elastic modulus, melting point, radiation damage, crack propagation, interface properties?
3. What **length scale** — single crystal, polycrystal, nanoparticle, surface slab?

### Stage 2: Simulation Details
4. What **potential / force field** do you want to use? (EAM, ReaxFF, Tersoff, MTP, MACE, SNAP, PACE/ACE, CHGNET, etc.)
5. What **ensemble** — NVE, NVT, NPT?
6. What **temperature and pressure**?
7. How long a simulation and what **timestep**?

### Stage 3: Output
8. What **observables** do you need — thermo output, dump trajectory, RDF, MSD, stress, heat flux?
9. Will you run on a **single core, multicore (MPI), or GPU**?

Once you have answers, build the complete input file and explain each section.

---

## 4. Error Diagnosis Protocol

When the user pastes a LAMMPS error message:

### Step 1: Classify the error
- **Initialization errors** (before `run`): syntax, missing files, wrong atom counts
- **Potential errors**: bad `pair_coeff`, missing atom types, cutoff issues
- **Neighbor list errors**: skin distance too small, atoms moved too far
- **Timestep errors**: integration instability (atoms flying apart)
- **Communication errors** (MPI): domain decomposition, ghost atom cutoffs

### Step 2: Common errors and fixes

| Error Pattern | Likely Cause | Fix |
|---------------|--------------|-----|
| `Cannot open file X` | Wrong path or filename | Check `read_data`, `pair_coeff` paths |
| `Lost atoms` | Exploding simulation | Reduce timestep; check potential; verify units |
| `Neighbor list overflow` | Too many neighbors | Increase `one` in `neigh_modify` or reduce density |
| `Invalid atom type` | Mismatch between data file and `pair_coeff` | Count atom types carefully |
| `Expected floating point` | Wrong argument type to a command | Check command syntax at docs |
| `Unrecognized pair style` | Package not compiled in | Recompile with correct package |
| `Fix X requires atom attribute Y` | Wrong `atom_style` | Use `atom_style charge` or `full` |
| `Dangerous builds performed` | Atoms moved > skin/2 in one step | Reduce timestep or increase skin |
| `PPPM grid too large` | Box too big for Ewald | Use `kspace_modify` or `pppm/disp` |
| `Variable X does not exist` | Referencing undefined variable | Define variable before use |

### Step 3: Always check
1. Fetch the relevant LAMMPS error page: https://docs.lammps.org/Errors_details.html
2. Ask for the full error context if only part is shown
3. Ask for the relevant section of the input file

---

## 5. Interatomic Potentials Guide

### 5a. Classical Potentials
| Style | Use Case | Doc |
|-------|----------|-----|
| `eam` | Metals (single element) | https://docs.lammps.org/pair_eam.html |
| `eam/alloy` | Metal alloys | same |
| `eam/fs` | Finnis-Sinclair form | same |
| `tersoff` | Covalent (Si, C, SiC) | https://docs.lammps.org/pair_tersoff.html |
| `sw` | Stillinger-Weber (Si, Ge) | https://docs.lammps.org/pair_sw.html |
| `reax/c` / `reaxff` | Reactive, multi-element | https://docs.lammps.org/pair_reaxff.html |
| `lj/cut` | Simple Lennard-Jones | https://docs.lammps.org/pair_lj.html |
| `morse` | Diatomic / simple bonds | https://docs.lammps.org/pair_morse.html |

### 5b. Machine Learning Potentials (MLPs)
These require special compilation. See Section 7.

| Potential | Package | Notes |
|-----------|---------|-------|
| MTP (moment tensor) | ML-IAP or ML-QUIP | Requires MLIP or QUIP library |
| SNAP / linear SNAP | ML-IAP | Bispectrum descriptors |
| ACE / PACE | ML-PACE | Requires pacemaker library |
| GAP (Gaussian Approx.) | ML-QUIP | Requires QUIP/GAP library |
| NNP (neural network) | ML-HDNNP or ML-IAP | Requires n2p2 or libtorch |
| MACE | ML-IAP (pytorch backend) | Requires libtorch |
| CHGNet / M3GNet | APIP package | Universal potentials |

---

## 6. File Formats

### 6a. LAMMPS Data File (read_data)
```
LAMMPS data file

N atoms
M atom types

xlo xhi
ylo yhi
zlo zhi

Masses
1 55.845   # Fe

Atoms   # atomic style: atom_id mol_id atom_type x y z
1 1 1 0.0 0.0 0.0
...
```
- `atom_style atomic`: id type x y z
- `atom_style charge`: id type charge x y z  
- `atom_style full`: id mol_id type charge x y z

For full format details: https://docs.lammps.org/read_data.html

### 6b. Dump File
```
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1000
ITEM: BOX BOUNDS pp pp pp
0.0 28.7
0.0 28.7
0.0 28.7
ITEM: ATOMS id type x y z
1 1 0.000 0.000 0.000
...
```
Common dump styles: `custom`, `atom`, `xyz`, `cfg`, `dcd`, `h5md`
Full format: https://docs.lammps.org/dump.html

### 6c. log.lammps
Key sections:
- **Thermo header**: column names printed once (Step, Temp, E_pair, E_mol, TotEng, Press...)
- **Thermo data**: one row per `thermo N` steps
- **Performance**: Loop time, CPU usage, MPI stats at end

To extract thermo data with Python:
```python
import pandas as pd

def read_lammps_log(filename):
    data, collecting = [], False
    headers = []
    with open(filename) as f:
        for line in f:
            if line.startswith("Step"):
                headers = line.split()
                collecting = True
                continue
            if collecting:
                if line.startswith("Loop") or line.startswith("WARNING"):
                    collecting = False
                    continue
                try:
                    data.append(list(map(float, line.split())))
                except ValueError:
                    collecting = False
    return pd.DataFrame(data, columns=headers)
```

### 6d. Restart Files
Binary format written by `write_restart` / `restart` fix.
Read with `read_restart` — does **not** restore `fix`, `compute`, `dump` commands; only atomic positions, velocities, topology.

---

## 7. Compiling LAMMPS with Special Packages

**Always use CMake** (traditional make no longer supports packages needing extra build steps as of Sept 2025).
Full reference: https://docs.lammps.org/Build_extras.html

### 7a. Standard Build (no special packages)
```bash
git clone https://github.com/lammps/lammps.git
cd lammps && mkdir build && cd build
cmake ../cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D PKG_MANYBODY=yes \
  -D PKG_EXTRA-FIX=yes
cmake --build . -- -j$(nproc)
```

### 7b. MTP via ML-IAP + MLIP Library
MTP (Moment Tensor Potentials, Shapeev) requires the `MLIP` library:
```bash
# Build MLIP library first (from mlip-3 or mlip-2 source)
# Then:
cmake ../cmake \
  -D PKG_ML-IAP=yes \
  -D PKG_MANYBODY=yes \
  -D MLIAP_ENABLE_PYTHON=yes \       # if Python interface needed
  -D MLIP_LIBRARY=/path/to/mlip/lib
```
For MTP via ML-QUIP: requires QUIP + quippy compiled with MTP support.
Doc: https://docs.lammps.org/Build_extras.html#ml-iap-package

### 7c. ACE / PACE via ML-PACE
```bash
cmake ../cmake \
  -D PKG_ML-PACE=yes
# Downloads and builds pacemaker automatically via FitSNAP/pace-maker
```
Then in input: `pair_style pace` + `pair_coeff * * potential.yace Fe`
Doc: https://docs.lammps.org/Build_extras.html#ml-pace-package

### 7d. GAP / NNP via ML-QUIP
Requires QUIP library (https://github.com/libAtoms/QUIP):
```bash
# Build QUIP first, set QUIP_ROOT
cmake ../cmake \
  -D PKG_ML-QUIP=yes \
  -D QUIP_LIBRARY=/path/to/QUIP/build/linux_x86_64_gfortran/libquip.a
```
In input: `pair_style quip` + `pair_coeff * * gap.xml "Potential xml_label=GAP" 14`
Doc: https://docs.lammps.org/Build_extras.html#ml-quip-package

### 7e. Universal Potentials (CHGNet, M3GNet) via APIP
```bash
cmake ../cmake \
  -D PKG_APIP=yes \
  -D PKG_PYTHON=yes
pip install chgnet m3gnet  # Python-side
```
Doc: https://docs.lammps.org/Build_extras.html#apip-package

### 7f. GPU Acceleration (KOKKOS or GPU package)
```bash
# KOKKOS with CUDA:
cmake ../cmake \
  -D PKG_KOKKOS=yes \
  -D Kokkos_ENABLE_CUDA=yes \
  -D Kokkos_ARCH_VOLTA70=yes   # or AMPERE80, HOPPER90 etc.
  -D CMAKE_CXX_COMPILER=${LAMMPS_DIR}/lib/kokkos/bin/nvcc_wrapper

# GPU package (legacy):
cmake ../cmake \
  -D PKG_GPU=yes \
  -D GPU_API=cuda \
  -D GPU_ARCH=sm_70
```
Run: `lmp -sf gpu -pk gpu 1 < in.script`

### 7g. OpenMP
```bash
cmake ../cmake -D PKG_OPENMP=yes
```
Run: `lmp -sf omp -pk omp 4 < in.script`

---

## 8. Common Simulation Recipes

### 8a. Equilibration Protocol
```lammps
# 1. Minimize first (always)
minimize 1e-6 1e-8 10000 100000

# 2. Ramp temperature NVT
fix 1 all nvt temp 10 300 0.1
run 10000
unfix 1

# 3. NPT to relax box
fix 2 all npt temp 300 300 0.1 iso 0 0 1.0
run 50000
unfix 2

# 4. NVT production
fix 3 all nvt temp 300 300 0.1
run 100000
```

### 8b. Computing RDF
```lammps
compute rdf all rdf 100 1 1   # 100 bins, type 1 - type 1
fix rdf_avg all ave/time 100 1 100 c_rdf[*] file rdf.dat mode vector
```

### 8c. Computing MSD and Diffusion
```lammps
compute msd all msd com yes
fix msd_out all ave/time 10 1 10 c_msd[4] file msd.dat
# D = slope(MSD vs t) / 6  (3D)
```

### 8d. Stress / Elastic Constants
```lammps
compute stress all stress/atom NULL    # per-atom stress tensor
compute pressure all pressure thermo_temp  # global
```

---

## 9. Debugging Checklist

When reviewing a user's input file, always check:

- [ ] `units` declared first and consistent with potential file
- [ ] `atom_style` matches data file format
- [ ] Number of atom types in `create_box` or data file matches `pair_coeff`
- [ ] `pair_style` declared before `pair_coeff`
- [ ] Timestep appropriate for units (metal: 0.001 ps; real: 1 fs typical)
- [ ] Thermostat damping parameter ~ 100× timestep (rule of thumb)
- [ ] `boundary` consistent with geometry (surfaces need `s` or `f`)
- [ ] `neigh_modify delay 0 check yes` set for non-equilibrium runs
- [ ] Dump file naming uses `*` or `%` for multi-file output
- [ ] `run 0` after `read_restart` to re-initialize before production

---

## 10. Reference Architecture

| Topic | Where | When |
|-------|-------|------|
| Potential-specific setup, file formats in depth | fetch the LAMMPS doc page (pair_*, read_data, dump) | User asks about a specific MLP / format detail |
| Parallel / HPC running, job files | `references/slurm.md` + SKILL.md step 5 | Cluster, SLURM, MPI scaling |
| Compiling LAMMPS on a cluster | §7 above + `references/hpc-install.md` | Build questions |
| KMC / rare events | `software/seakmc/GUIDE.md` | Long-timescale dynamics |

Always web_fetch the specific LAMMPS doc page when a user asks about a command you're not 100% certain about. The docs are at https://docs.lammps.org/[command_name].html (e.g., https://docs.lammps.org/fix_nvt.html).

---

## 11. Workflow integration (materials-simulation-skill steps 3–7)

**Files in each run dir:** the structure from step 2 (`read_data structure.data`), the input
template rendered from `reffiles/in.lammps`, the potential file(s) (`copy` or `link` from
reffiles — never downloaded unless the user says so), and the job file.

**Making the reference input a template** — two equivalent styles:

```lammps
# (a) placeholders, filled by setup_runs.py
variable        T     equal {{T}}
variable        seed  equal {{SEED}}
read_data       structure.data
pair_style      eam/alloy
pair_coeff      * * {{POT}} Fe He
velocity        all create ${T} ${seed} mom yes rot yes dist gaussian
```

```bash
# (b) untouched input, values passed on the command line in the job body (below "# @RUN")
srun lmp -in in.lammps -var T {{T}} -var seed {{SEED}} -log log.lammps > lammps.out 2>&1
```

Use (b) when the user wants one identical input file everywhere. Either way, give every run
its own seed for statistics (`{{SEED}}` is deterministic per run dir, so reruns reproduce).

**Make results easy to extract:** end the input with explicit `print` lines, e.g.
```lammps
variable Ef equal pe
print "Efinal = ${Ef}"
```
Energies for binding-energy studies (He–V around a dislocation, etc.) are typically single-point
or minimized energies of several configurations: plan which run dirs/reference files provide
each term (e.g. the cluster near the core, the cluster far from it, the defect-free reference)
before submitting, and record the formula in `sim_plan.yaml` so step 8 can compute it.

**Default checks** (confirm with the user):
```yaml
checks:
  - {type: contains,     file: log.lammps, pattern: "Total wall time"}
  - {type: not_contains, file: log.lammps, pattern: "^ERROR|Lost atoms"}
```

**Default summary fields:** `lammps_thermo` for any thermo column (pick the block with `run:`;
`-1` is the last run/minimize), `regex` for `print` lines.

**Job body:** `srun lmp -in in.lammps -log log.lammps > lammps.out 2>&1`
(executable name differs per cluster: `lmp`, `lmp_mpi`, `lmp_intel_cpu_intelmpi` ... — take it
from the user's reference job).
