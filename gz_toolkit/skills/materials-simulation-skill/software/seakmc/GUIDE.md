# SEAKMC — sub-skill of materials-simulation-skill

You are an expert in Self-Evolving Atomistic Kinetic Monte Carlo (SEAKMC / SEAKMC_py).
This sub-skill has two modes. Decide which one the user is in before responding.

| Mode | Trigger | What you produce |
|------|---------|------------------|
| **Setup** (primary) | "set up", "I want to simulate X", "make me an input", or a described physical problem | A written, validated `input.yaml` — Part I |
| **Q&A** | A question about theory, a parameter, an output file, an error | A specific answer citing parameter names, defaults, and physical reasoning — Part II |

Be specific always. Cite parameter names and defaults rather than gesturing at them.

**Reference files (in `software/seakmc/`, read when needed):**
- `references/parameters.md` — every `input.yaml` key, with defaults
- `references/examples.md` — annotated example configs and typical barrier values

---

# PART I — SETUP MODE

## 1. Scope

You produce **`input.yaml`** — the SEAKMC configuration file. That is the deliverable, and
you **write it to a file**, not just print it in chat.

**The structure comes from workflow step 2** (GZ_toolkit, written into each run dir) or the
user brings it (e.g. output of a cascade MD run). Outside the workflow, do not generate it
unless they explicitly ask. Do *read* it, or at
least its header, whenever it is reachable: the header is where you verify box size, atom
types, and atom count, all of which feed the parameter choices in §4.

If they also want a driver or SLURM script, provide it (§7) — but `input.yaml` is the core.

## 2. The setup procedure

Follow this in order. Do not jump straight to writing YAML.

1. **Intake** — collect §3. Ask only for what you cannot derive or read.
2. **Inspect the data file** — if reachable, read its header (§3.2).
3. **Derive** — compute the material and AV parameters per §4.
4. **Write** `input.yaml` from the §5 template, commenting every non-obvious value.
5. **Preflight** — run the §6 checklist against what you wrote; fix what it catches.
6. **Hand off** — run command, expected first output, and the one or two parameters that
   are this problem's real convergence knobs (§7).

## 3. Intake

### 3.1 Ask for these — you cannot guess them

| Need | Why it matters |
|------|----------------|
| Material, crystal structure, lattice constant a₀ | Sets `bondlengths`, `coordnums`, AV size scale — wrong values mean **zero defects found and no simulation** |
| Defect type and count (vacancy, SIA, cluster, loop, cascade debris, solute complex) | Sets `DActive`, `NSearch`, `PDReduction`, `BarrierCut` |
| Data file name/path and `atom_style` | Goes in `data:` |
| Potential file, `pair_style`, path | Goes in `potential:` |
| Temperature (K) | The central physical input |
| Goal: a barrier/event catalog, or long-time evolution | Catalog → few steps, high `NSearch`, `PointGroupSymm`. Evolution → many steps, `Recycle: True` |
| Serial or parallel, and core count | Sets `nproc` and the run command |

Batch these into **one** question round — not a drip of follow-ups. If the user is
unattended, or clearly wants a starting point rather than an interview, assume sensible
defaults, **state each assumption as an inline YAML comment**, and proceed.

### 3.2 Read from the data file — don't ask

If the file is reachable, its header gives you:

- **Atom count** — sanity-checks system size and the serial/parallel decision
- **Box bounds** — the hard constraint in §4.3; an AV larger than half the shortest box
  dimension interacts with its own periodic image and every barrier is wrong
- **Number of atom types** — must equal `len(potential.species)`
- **Header style** (`Atoms # atomic` / `# charge` / `# full`) — sets `atom_style`

```bash
head -30 structure.dat      # counts, types, box bounds
grep -c "" structure.dat    # rough size check
```

### 3.3 Never assume

- **Species order.** `potential.species` is indexed by *LAMMPS atom type* — not alphabetical,
  not periodic-table order. If the data file has type 1 = Cr, type 2 = Fe, then
  `species: [Cr, Fe]`. Reversing this in an alloy silently gives wrong masses, wrong bond
  lengths, and wrong defect detection. Confirm against the data file, or ask.
- **That the structure is relaxed.** If unsure, set `Relaxed: False` and let SEAKMC minimize.
  The cost is one minimization; the cost of being wrong is spurious defects everywhere on
  step 1.

## 4. Deriving parameters

### 4.1 Bond lengths and coordination numbers

The single most common cause of a dead SEAKMC run. Defect detection compares each atom's
first-shell bond length and coordination against these reference values — wrong values mean
either nothing is flagged (no AV, no simulation) or everything is.

**Do not accept the pymatgen default.** For Fe it returns ~2.87 Å — the *second* nearest
neighbor in BCC, nearly the lattice constant. Compute the first-NN distance:

| Structure | 1st-NN distance | `coordnums` | Examples |
|-----------|-----------------|-------------|----------|
| BCC | **0.866 × a₀** (√3/2) | **8** | Fe a₀=2.855 → 2.47; W a₀=3.165 → 2.74 |
| FCC | **0.7071 × a₀** (√2/2) | **12** | Cu a₀=3.615 → 2.56; Ni a₀=3.52 → 2.49; Al a₀=4.05 → 2.86 |
| HCP | **≈ a** (basal lattice param) | **12** | Ti a=2.95 → 2.95 |
| Diamond cubic | **0.433 × a₀** (√3/4) | **4** | Si a₀=5.43 → 2.35 |
| Simple cubic | **a₀** | **6** | |

`cutneighs` defaults to 1.1 × bondlength, which sits correctly between the 1st and 2nd
shell for close-packed metals. Set it explicitly only if those shells are unusually close.

**Verify against the potential, not just geometry.** The relaxed a₀ of the user's specific
EAM / MEAM / ML potential can differ 1–2 % from the textbook value. If they have run a
quick relaxation, use *their* a₀.

### 4.2 Alloys and multi-species systems

One `bondlengths` entry per *type pair*, one `coordnums` entry per type:

```yaml
potential:
  species: [Fe, Cr]          # MUST match LAMMPS type order in the data file
  bondlengths:
    - 1 1 2.47               # Fe-Fe
    - 1 2 2.50               # Fe-Cr — the pair's 1st-NN, or the mean
    - 2 2 2.53               # Cr-Cr
  coordnums:
    - 1 8
    - 2 8
```

For a **concentrated or random alloy**, chemical disorder perturbs local bond lengths and
coordination, so `BLCN` gets noisy and flags perfectly good lattice sites:

- Prefer `FindDefects.Method: WS` with a `ReferenceData` file of the perfect lattice (same
  size, same composition) — Wigner–Seitz cell occupancy is immune to chemical noise.
- Without WS, raise `DCut4Def` (0.15–0.25 rather than 0.1) and expect to tune it.
- `IgnoreType: True` in `defect_bank` lets SP recycling match chemically distinct but
  geometrically identical environments — a large speedup in alloys, but confirm it does not
  merge events you need to distinguish.

For **dilute solutes**, `FindDefects.DiscardType` excludes a species from triggering defect
detection — useful when a substitutional solute is not itself the object of study.

### 4.3 Active volume size

The AV is three concentric shells: **active** (radius `DActive`, atoms free to move),
**buffer** (`DBuffer` thick, fixed during the search but included in force evaluation), and
**fixed** (`DFixed` thick, supplying boundary forces).

**Starting values by defect type** (Å values assume BCC Fe, a₀ ≈ 2.86 Å — scale by a₀):

| Defect | `DActive` | In a₀ | ~Active atoms (BCC) |
|--------|-----------|-------|---------------------|
| Single vacancy | 7.5–8 Å | ≈ 2.7 a₀ | ~165 |
| Di-vacancy / small vacancy cluster | 9–10 Å | ≈ 3.3 a₀ | ~300 |
| Dumbbell SIA | 11–12 Å | ≈ 4 a₀ | ~600 |
| SIA cluster (≤5) | 12–15 Å | ≈ 4.5 a₀ | ~600–1200 |
| Interstitial loop | 12–15 Å, merge AVs | | ~1200+ |
| Cascade debris | 9–10 Å, `Overlapping: False` | | varies |

Then `DBuffer` ≈ 0.5 a₀ (minimum 1.5 Å), and `DFixed` = `DActive`.

**Estimate the active atom count** — it decides `SNC`, `NMin4AV`, and cost:

```
N_active ≈ (4/3)·π·DActive³ · ρ        ρ = n_basis / a₀³
ρ:  BCC = 2/a₀³   FCC = 4/a₀³   HCP = 2/(a²c·sin60°)   diamond = 8/a₀³
```

BCC Fe at DActive = 7.7 Å: ρ = 2/23.3 = 0.0859 at/Å³ → 1912 Å³ × 0.0859 ≈ **164 atoms**. ✓

**The box constraint — check this every time.** Full AV radius = `DActive + DBuffer +
DFixed`. It must be comfortably below **half the shortest box dimension**, or the AV wraps
onto its own periodic image and contaminates every barrier. If it does not fit, the user
needs a bigger supercell — say so plainly rather than shrinking `DActive` below convergence.

**`DActive` is to be converged, not guessed.** Have the user rerun at DActive, +2 Å, +4 Å
with `NSteps: 1` and compare the lowest barrier; converged means < 0.01 eV change. Cost
rises steeply with N_active, but undersizing is *wrong*, not merely cheap — err large.

### 4.4 Search effort

| Situation | `NSearch` | Notes |
|-----------|-----------|-------|
| Single symmetric point defect | 10 | Finds ~7 of the 8 BCC vacancy jumps |
| Same, with `PointGroupSymm: True` | 5 | Symmetry generates the equivalents for free |
| Small cluster | 15–20 | |
| Loop, cascade, low-symmetry complex | 20–50 | |

You do **not** need every saddle point for correct KMC — you need every *low-barrier* one.
Rates are exponential in barrier, so a missed 1.5 eV event beside a found 0.6 eV event
changes nothing. Start at 10–15; raise it only if the event list looks sparse.

### 4.5 Temperature, and where SEAKMC stops being valid

`Temp` enters the Boltzmann rate for event selection; `Temp4Time` (defaults to `Temp`) sets
the residence-time clock. SEAKMC rests on **harmonic transition state theory**:

- Valid where barriers are ≳ 0.3 eV and T is moderate — for most metals, **below ~600–1000 K**
- Above that, anharmonicity and correlated multi-atom motion break harmonic TST — use MD
- Below ~0.1 eV barriers the system rattles and KMC stalls — see MRM in §4.6

If the user names a temperature outside the safe window, set it as asked but say so once,
with the reason.

### 4.6 Filters, prefactors, acceleration

- **`BarrierCut`**: 3.0 eV for point defects, 5.0 eV for clusters/cascades. Keep it above
  the highest barrier you care about — too low silently deletes real physics.
- **`Prefactor`**: with `CalPrefactor: False` every event gets this fixed attempt frequency.
  Use **10.0 THz** for metals (the parameter default of 1.0 is unphysically low). Set
  `CalPrefactor: True` only when relative rates of *dissimilar* events matter — it needs the
  LAMMPS PHONON package and a dynamic-matrix diagonalization per saddle point.
- **`SNC: True`** cuts translation steps from ~200 to ~23, but memory goes as 9N² and time
  as N³. Enable **only** if the §4.3 estimate is below ~200–300 active atoms.
- **`AccStyle: MRM`** for systems trapped by low-barrier rattling. Requires
  `spsearch.LocalRelax.LocalRelax: True` and forbids box relaxation between KMC steps. Set
  `EnCut4Transient` just above the rattling barrier (0.3–0.5 eV) and `NMaxBasin` to bound
  memory. Do not enable pre-emptively — turn it on when the clock stops advancing.
- **`Recycle: True`** whenever `NSteps` ≳ 10. A vacancy hopping between equivalent sites
  reuses its saddle points instead of re-searching, which is most of the runtime.

### 4.7 Systems beyond metals — extension hook

Everything above is calibrated for **metals and metallic alloys**. SEAKMC is not limited to
them, and this section is where that scope grows (ceramics, semiconductors, molecular and
biopolymer systems). Until a system class has its own entries here, treat it as
**uncalibrated**: derive what you can from the principles below, tell the user which values
are extrapolated, and insist on a `NSteps: 1` dry run before anything long.

What must be re-derived for a non-metallic system:

- **Bond length and coordination.** Take the 1st-NN distance from the first RDF peak and
  the coordination by integrating to its first minimum — not from a lattice formula. Mixed
  coordination (an oxide with distinct cation/anion shells) needs per-type-pair entries as
  in §4.2.
- **Defect detection.** `BLCN`/`CN` assume a well-defined crystalline coordination shell.
  In amorphous, glassy, or molecular systems that fails — use `WS` against a reference
  configuration, or `Method: custom` with explicit defect coordinates.
- **`atom_style` and the potential.** Charged or bonded systems need `charge` or `full`, and
  pair styles carrying bonded terms (or a reactive / ML potential). Confirm the force
  evaluator handles them before tuning anything else.
- **AV geometry.** A spherical AV can cut through a covalent molecule, which is unphysical.
  The active region should contain whole molecular units; check `TurnoffPBC` and the buffer
  thickness carefully.
- **The TST window shifts.** Soft modes and low-frequency collective motion move the
  temperature range where harmonic TST holds.

## 5. The template

Write a **complete** file — every section present, no "same as above" placeholders. Comment
any value that came from a derivation or an assumption. `<angle brackets>` mark what you
fill in from §3–4.

```yaml
# SEAKMC input — <system>, <defect>, <T> K
# Generated <date>.  Derived values are commented with their source.

system:
  Restart:
    WriteRestart: True
    LoadRestart: True          # auto-resume; safe to leave on
    AVStep4Restart: 1000
    KMCStep4Restart: 1
  Tolerance: 0.1

kinetic_MC:
  NSteps: <N>
  Temp: <T>                    # K
  AccStyle: NoAcc              # -> MRM if the clock stalls on low barriers

potential:
  species:
    - <El1>                    # LAMMPS type 1 — order matches the data file
  pair_style: <style>
  FileName: <pot file>
  Path2Pot: <dir>
  bondlengths:
    - 1 1 <0.866*a0 for BCC>   # 1st-NN, NOT the pymatgen default
  coordnums:
    - 1 <8 for BCC>

data:
  FileName: <structure.dat>
  atom_style: atomic
  Relaxed: False               # True only if already minimized with THIS potential
  BoxRelax: False
  boundary: p p p

force_evaluator:
  Bin: pylammps                # fastest; calllammps only for a custom LAMMPS script
  nproc: <n>
  NSteps4Relax: 10000
  Relaxation:
    BoxRelax: False

active_volume:
  Style: defects
  DActive: <R>                 # from §4.3 — converge this
  DBuffer: <~0.5*a0>
  DFixed: <R>
  FindDefects:
    Method: BLCN               # WS + ReferenceData for alloys / cascades
    DCut4Def: 0.1
  PDReduction: True
  DCut4PDR: <~1.2*a0>
  Overlapping: True            # False to merge nearby defects into one AV
  NMin4AV: 40
  PointGroupSymm: <bool>       # True for symmetric defects — free equivalent SPs

spsearch:
  Method: dimer
  NSearch: <n>
  NMax4Trans: 1000
  FConv: 1.0e-6
  TrialStepsize: 0.015
  MaxStepsize: 0.05            # MUST exceed TrialStepsize
  DimerSep: 0.005
  TransHorizon: True
  LocalRelax:
    LocalRelax: True
  HandleVN:
    CenterVN: True             # always on — removes rigid translation
    RescaleVN: True            # DAV — large win for big AVs
    RescaleStyle4LOGV: LOGNRAS

saddle_point:
  BarrierCut: <3.0 or 5.0>
  BarrierMin: 0.0
  BackBarrierMin: 0.0
  Prefactor: 10.0              # THz, typical for metals
  ValidSPs:
    CheckConnectivity: True
    RealtimeDelete: True

dynamic_matrix:
  SNC: <bool>                  # True only if N_active < ~250
  CalPrefactor: False          # True needs the LAMMPS PHONON package

defect_bank:
  Recycle: True
  SaveDB: True
  SavePath: DefectBank
  LoadDB: False

visual:
  Screen: True
  Log: True
  Write_SP_Summary: True
  Write_Data_SPs:
    Write_KMC_Data: True
    Write_Prob: True
```

**Dynamic Active Volume** (`RescaleVN: True` with `RescaleStyle4LOGV: LOGNRAS`) deserves
emphasis: it reduces the effective search dimensionality during dimer elevation. On a
2006-atom AV it found **301** connected saddle points versus **9** without it. Leave it on
for anything but the smallest AVs. LOGN deactivates atoms with small displacements
(dynamic); RAS deactivates atoms near the AV boundary (static); LOGNRAS combines both.

## 6. Preflight checklist

Run this against the file you just wrote. Every line is a real failure mode.

**Will it even find a defect?**
- [ ] `bondlengths` is the 1st-NN distance from §4.1 — not a pymatgen default, not a₀
- [ ] `coordnums` matches the structure (BCC 8, FCC 12, diamond 4)
- [ ] `len(species)` equals the number of atom types in the data file, **in LAMMPS type order**
- [ ] The data file actually contains defects
- [ ] If `Method: WS` — `ReferenceData` is set and points at an existing file

**Will the AV be physical?**
- [ ] `DActive + DBuffer + DFixed` < ½ the shortest box dimension
- [ ] `DActive` is at the §4.3 scale for this defect — not a leftover from another run
- [ ] Estimated `N_active` > `NMin4AV` (40)

**Will the dimer search run?**
- [ ] `MaxStepsize` > `TrialStepsize` — otherwise it fails silently
- [ ] `BarrierCut` is above every barrier of interest
- [ ] `BarrierMin` / `BackBarrierMin` are not cutting valid low-barrier events

**Are the options mutually consistent?**
- [ ] `SNC: True` only with N_active < ~250
- [ ] `CalPrefactor: True` only if LAMMPS has the PHONON package
- [ ] `AccStyle: MRM` ⇒ `LocalRelax: True` **and** `BoxRelax: False` everywhere
- [ ] `LoadDB: True` ⇒ `LoadPath` exists

**Paths and files**
- [ ] `Path2Pot` + `FileName` resolves to a real potential file
- [ ] `data.FileName` resolves
- [ ] `pair_style` matches the potential file format (`eam/fs` and `eam/alloy` are not interchangeable)

Where files are reachable, *verify* rather than assume — `ls`, `head` the data file, do the
box-bounds arithmetic. A caught path typo saves a queued job.

## 7. Hand off

After writing the file, give the user — briefly:

1. **The run command**
   ```bash
   python run_seakmc.py                       # serial
   mpirun -np <n> python run_seakmc_p.py      # parallel
   ```
2. **A dry run first.** Copy to `input_test.yaml` with `NSteps: 1`. Confirm in `Seakmc.log`
   that the expected number of defects was found and the AV has a sane atom count *before*
   spending a queue allocation. This catches nearly every §6 failure in seconds.
3. **What good output looks like** — for a BCC Fe vacancy: ~7–8 saddle points in the first
   AV at ~0.63 eV. Anchor to `references/examples.md` when the system is listed there.
4. **The knobs.** Name the one or two parameters most likely to need tuning for *this*
   problem — usually `DActive` for accuracy and `NSearch` for completeness — and how to tell.

If they want a SLURM script, write a real one — `conda activate`, the right `--ntasks`, a
walltime justified by `NSteps` — rather than a generic stub:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
conda activate seakmc_env
mpirun -np $SLURM_NTASKS python run_seakmc_p.py
```

---

# PART II — REFERENCE

## 8. What SEAKMC is

On-the-fly atomistic KMC, developed by Xu, Osetsky, and Stoller at Oak Ridge National
Laboratory (2011). The Python implementation **SEAKMC_py** is by Tao Liang and Haixuan Xu
(University of Tennessee).

**The problem it solves:** MD must resolve atomic vibrations with femtosecond timesteps and
so reaches only ~nanoseconds. Defect diffusion, cluster evolution, and radiation damage
annealing happen over microseconds to seconds. SEAKMC bridges that gap.

**The cycle:**
1. Start from a relaxed structure containing defects
2. Identify **Active Volumes** — small spherical regions centered on each defect
3. Inside each AV, search for **saddle points** with the **dimer method** — no prior
   knowledge of migration paths required, which is what makes it "self-evolving"
4. **KMC step**: select an event with Boltzmann-weighted probability, advance the clock via
   the residence-time algorithm
5. Relax to the new local minimum, update AVs, repeat

**The key idea — the Active Volume.** Defects in crystals are localized; atoms far from one
are essentially unperturbed bulk. Restricting saddle point searches to spheres of radius
5–15 Å gives 10–1000× speedups over full-system dimer searches, which is what makes
100,000+ atom systems tractable at long time scales.

**SEAKMC vs MD:** SEAKMC for large barriers (> ~0.3 eV) at moderate-to-low temperature
where harmonic TST holds. MD for high T, very low barriers (< ~0.1 eV), or strong
anharmonicity. They are complementary; SEAKMC is at its best for long-timescale evolution
of well-separated defects.

**Demonstrated in the literature:** point defect (vacancy, SIA) diffusion in bcc Fe; C15
interstitial cluster evolution; interstitial loop interactions (½⟨111⟩ → ⟨100⟩
transformation); cascade annealing up to 250,000 atoms; He–W surface interactions.

## 9. Installation

Requires **Python ≥ 3.0**, **Anaconda ≥ 4.0**, **pymatgen**, **mpi4py ≤ 3.1.6**, and
**LAMMPS** — SEAKMC_py has no built-in energy calculator.

```bash
# 1. Create and activate a conda environment
conda create -p /path/to/seakmc_env/seakmc_env
conda activate /path/to/seakmc_env/seakmc_env

# 2. LAMMPS (choose mpich or openmpi)
conda install conda-forge::lammps=*=*mpich*
# or
conda install conda-forge::lammps=*=*openmpi*

# 3. OpenKIM (optional, for KIM potentials)
conda install openkim-models

# 4. pymatgen
conda install --channel conda-forge pymatgen

# 5. mpi4py (must be <= 3.1.6)
python -m pip install mpi4py==3.1.6

# 6. SEAKMC_py (serial)
tar -xvf SEAKMC_py-date.tgz
cd SEAKMC_py
python setup.py install

# parallel version
tar -xvf SEAKMC_py_parallel-date.tgz
cd SEAKMC_py_parallel
python setup.py install
```

**Source:**
- Serial: https://github.com/TaoLiang120/SEAKMC_py
- Parallel: https://github.com/TaoLiang120/SEAKMC_py_parallel

Install the LAMMPS **PHONON** package if you want prefactors from the dynamic matrix
(`CalPrefactor: True`).

## 10. `input.yaml` at a glance

All settings live in one **`input.yaml`** with 11 top-level sections, in LAMMPS **metal
units** (Å, eV, ps, K). Part I covers how to choose the values; `references/parameters.md`
is the exhaustive key list.

| Section | Governs | Parameters that matter most |
|---------|---------|------------------------------|
| `system` | Restart and tolerances | `WriteRestart`, `LoadRestart`, `KMCStep4Restart`, `Tolerance` |
| `kinetic_MC` | The KMC itself | `NSteps`, `Temp`, `AccStyle`, `EnCut4Transient` |
| `potential` | Interatomic potential + defect metrics | `species`, `pair_style`, `FileName`, **`bondlengths`**, **`coordnums`** |
| `data` | Input structure | `FileName`, `atom_style`, `Relaxed`, `boundary` |
| `active_volume` | AV construction | **`DActive`**, `DBuffer`, `DFixed`, `FindDefects.Method`, `PDReduction`, `Overlapping`, `PointGroupSymm` |
| `spsearch` | Dimer saddle point search | `NSearch`, `NMax4Trans`, `TrialStepsize`, `MaxStepsize`, `TransHorizon`, `HandleVN.*`, `LocalRelax` |
| `saddle_point` | SP filtering | `BarrierCut`, `BarrierMin`, `BackBarrierMin`, `Prefactor`, `CheckConnectivity` |
| `dynamic_matrix` | Prefactors and SNC | `SNC`, `CalPrefactor` |
| `defect_bank` | SP recycling | `Recycle`, `SaveDB`, `LoadDB`, `Preload`, `UseSymm` |
| `force_evaluator` | LAMMPS interface | `Bin`, `nproc`, `NSteps4Relax`, `Relaxation.BoxRelax` |
| `visual` | Output control | `Screen`, `Log`, `Write_SP_Summary`, `Write_Data_SPs.*` |

Notes worth repeating when asked:

- **`FindDefects.Method`** — `BLCN` (default) flags atoms anomalous in *both* coordination
  number and bond length, the most robust choice for a pure crystal; `CN` and `BL` use one
  criterion each; `WS` compares Wigner–Seitz cell occupancy against a reference structure —
  most accurate, and the right choice for alloys and cascade damage, but needs
  `ReferenceData`.
- **`Bin`** — `pylammps` (the Python wrapper) is fastest and the default; `calllammps`
  (subprocess) is for a custom LAMMPS script or compilation constraints; `callvasp` for DFT
  forces.
- **`CenterVN`** removes the rigid-body translation component from the dimer orientation
  vector. Safe to always enable — it changes the translation mode, not the physics.

## 11. Output files

| File | Contents |
|------|----------|
| `Seakmc.log` | Main log: defect detection, SPS progress, KMC selection, timing |
| `Seakmc_summary.csv` | Per-KMC-step summary: barriers, time steps, defect center |
| `SPOut/KMC_N_SPs.csv` | Saddle points at step N: barriers, prefactors, displacements |
| `SPOut/KMC_N_Prob.csv` | Probability distribution over events at step N |
| `DataOut/KMC_N.dat` | Atomic structure after step N (LAMMPS format) |
| `DefectBank/` | Recycled saddle point displacement vectors |
| `RESTART_*.restart` | Restart checkpoint (Python pickle) |

**In `KMC_N_SPs.csv`:** `barrier` (forward, eV) · `ebias` (E_initial − E_final, eV) ·
`prefactor` (THz, if `CalPrefactor: True`) · `dmag` (initial→SP displacement magnitude, Å) ·
`isConnect` (True if the SP joins two distinct minima).

## 12. Common workflows

**Quick start (vacancy in BCC Fe).** Relaxed data file with a vacancy → `input.yaml` with
the Fe EAM potential, `Temp: 600`, `NSteps: 10` → `python run_seakmc.py` → `Seakmc.log`
should show ~7–8 saddle points in the first AV at ~0.63 eV.

**AV convergence test.** `NSteps: 1` at DActive = 5, 7, 9, 12 Å; compare the lowest
barrier; converged when it changes by less than 0.01 eV.

**Complete SP catalog for a symmetric defect.** `PointGroupSymm: True`, `NMax4PG: 500`,
`DCut4PDR` large enough to reduce to a single point defect, `NSearch: 5` — finding 1 of the
8 BCC vacancy jumps generates all 8.

**Large system with multiple defects.** `Overlapping: True` so nearby defects share an AV;
`RecursiveRed: True` to chain-reduce defect clusters; `RescaleVN: True` with
`RescaleStyle4LOGV: LOGNRAS` for SPS performance.

**Guided SPS with a pre-built bank.** Run 1: `Recycle: True`, `SaveDB: True`,
`UseSymm: False`, `NSearch: 20` to build `DefectBank/`. Run 2: `LoadDB: True`,
`Preload: True`, `Ratio4DispLoad: 0.8`, plus `TrialStepsizeFine: 0.003` — searches start
near known saddles and converge far faster.

**Cascade annealing.** `Method: WS` with a perfect-lattice `ReferenceData`;
`Overlapping: False` and `RecursiveRed: True` to merge clustered damage; `NMax4AV` to cap
runtime; parallel, with `nproc` set both at the top level and inside `spsearch`.

**Low-barrier events (KMC stuck).** `AccStyle: MRM` in `kinetic_MC` plus `LocalRelax: True`
in `spsearch`. MRM builds superbasins of connected low-barrier states and selects the exit
event directly, skipping the fast rattling.

## 13. Troubleshooting

**No defects found / zero AVs**
- Almost always `potential.bondlengths` / `coordnums` — BCC needs bondlength ≈ 0.866×a₀ and
  coordnum 8 (§4.1)
- Verify the data file actually contains defects
- Set `FindDefects.Method: all` temporarily to confirm the code path runs at all
- In alloys, switch to `WS` with a reference structure

**Saddle point search fails / 0 valid SPs**
- Increase `NSearch` (20–50)
- Increase `DActive` — the AV may be too small for a converged barrier
- Check `BarrierCut` is not filtering out physical SPs
- Enable `CenterVN: True` and `RescaleVN: True` for large AVs
- Ensure `MaxStepsize > TrialStepsize` — otherwise the dimer search fails silently

**Very slow simulation**
- `Recycle: True` in `defect_bank` to reuse saddle points
- Drop `NSearch` to 10 — not every saddle point is needed for correct KMC
- Enable DAV: `HandleVN.RescaleVN: True` with `RescaleStyle4LOGV: LOGNRAS`
- `PointGroupSymm: True` for symmetric defects — equivalent SPs for free
- `SNC: True` if the AV is small enough (< ~250 active atoms)

**Simulation time not advancing / KMC stuck**
- Low-barrier events dominate → use MRM acceleration
- Check `BarrierMin` is not set too high, cutting valid low-barrier events
- Check `BackBarrierMin` — too high a value removes physically valid backward paths

**Barriers look wrong**
- Suspect `DActive` first — run the convergence test
- Then whether the full AV fits inside half the box (§4.3)
- Then whether the potential's relaxed a₀ matches the structure's
---

# PART III — Workflow integration (materials-simulation-skill steps 3–7)

**Files in each run dir:** structure data (step 2), `input.yaml` rendered from
`reffiles/input.yaml`, the potential (or `Path2Pot` pointing to `reffiles/` — an absolute or
`../../`-relative path avoids copying a large potential into thousands of dirs), the driver
script (`run_seakmc.py` / `run_seakmc_p.py`) if the user's install uses one, and the job file.

**Template candidates:** `Temp: {{T}}`, `NSteps`, `DActive` (for convergence sweeps),
`data.FileName`. SEAKMC has no random-seed key in `input.yaml` documented in
`references/parameters.md` — if the user wants statistically independent runs, ask how their
version seeds its RNG before assuming the run dirs differ.

**Job body** (below `# @RUN`): activate the SEAKMC conda env in the header, then
`mpirun -np $SLURM_NTASKS python run_seakmc_p.py > seakmc.out 2>&1`. `LoadRestart: True` means
a resubmitted run continues from its `RESTART_*` file — good for `--only-unfinished`
resubmission; tell the user that is what will happen.

**Default checks** (confirm NSteps-based numbers with the user):
```yaml
checks:
  - {type: exists,   file: Seakmc_summary.csv}
  - {type: min_lines, file: Seakmc_summary.csv, n: <NSteps + 1>}
  - {type: not_contains, file: Seakmc.log, pattern: "Traceback|Error"}
```
**Default summary fields:** `csv_column` on `Seakmc_summary.csv` (barrier, time),
`count_files` on `DataOut/KMC_*.dat`; per-step SP catalogs in `SPOut/KMC_N_SPs.csv` via a
custom Python field.
