# SEAKMC_py Example input.yaml Configurations

---

## Example 1: Single Vacancy in BCC Fe (Quick Start)

Minimal working example for a single vacancy in BCC iron using an EAM/FS potential.

```yaml
system:
  Restart:
    WriteRestart: True
    LoadRestart: True
    KMCStep4Restart: 1

kinetic_MC:
  NSteps: 20
  Temp: 600          # K — good for Fe vacancy diffusion studies

potential:
  species:
    - Fe
  pair_style: eam/fs
  FileName: Fe-P.eam.fs
  Path2Pot: /path/to/potentials/
  bondlengths:
    - 1 1 2.47        # 0.866 * 2.855 Å for BCC Fe
  coordnums:
    - 1 8             # 8 first-nearest neighbors in BCC

force_evaluator:
  Bin: pylammps
  NSteps4Relax: 10000

data:
  FileName: Fe_vacancy.dat
  atom_style: atomic
  Relaxed: False
  boundary: p p p

active_volume:
  Style: defects
  DActive: 7.7        # ~2.7 lattice constants — converged for single vacancy
  DBuffer: 1.5
  DFixed: 7.7
  FindDefects:
    Method: BLCN
    DCut4Def: 0.1
  PDReduction: True
  DCut4PDR: 3.5
  Overlapping: True
  NMin4AV: 40

defect_bank:
  Recycle: True
  SaveDB: True
  SavePath: DefectBank

dynamic_matrix:
  SNC: False
  CalPrefactor: False

spsearch:
  Method: dimer
  NSearch: 10
  NMax4Trans: 1000
  FConv: 1.0e-6
  TrialStepsize: 0.015
  MaxStepsize: 0.05
  DimerSep: 0.005
  TransHorizon: True
  LocalRelax:
    LocalRelax: True
  HandleVN:
    CenterVN: True
    RescaleVN: True
    RescaleStyle4LOGV: LOGNRAS

saddle_point:
  BarrierCut: 3.0
  BarrierMin: 0.0
  BackBarrierMin: 0.0
  Prefactor: 10.0     # Default 10 THz (typical for metals)
  ValidSPs:
    CheckConnectivity: True
    RealtimeDelete: True

visual:
  Screen: True
  Log: True
  Write_SP_Summary: True
  Write_Data_SPs:
    Write_KMC_Data: True
    Write_Data_SPs: True
    Write_Prob: True
```

**Expected output:** SEAKMC_py should find 7–8 saddle points in the first AV (the 8 first-
nearest-neighbor jumps for a BCC vacancy). `Seakmc.log` shows `barrier` ≈ 0.63–0.64 eV
for these events.

---

## Example 2: Large System with Multiple Defects (Cascade Annealing)

For annealing an MD cascade with many defects. Uses parallel version, recycling, and DAV.

```yaml
system:
  Restart:
    WriteRestart: True
    LoadRestart: True
    AVStep4Restart: 500
    KMCStep4Restart: 1
  Tolerance: 0.1

kinetic_MC:
  NSteps: 500
  Temp: 650

potential:
  species:
    - Fe
  pair_style: eam/alloy
  FileName: Fe_Mendelev.eam.alloy
  bondlengths:
    - 1 1 2.47
  coordnums:
    - 1 8

force_evaluator:
  Bin: pylammps
  nproc: 4               # Processors for whole-system relaxation
  NSteps4Relax: 10000

data:
  FileName: cascade_128k.dat
  Relaxed: True          # Cascade output is already relaxed
  boundary: p p p

active_volume:
  Style: defects
  DActive: 9.0           # Larger for interstitial clusters
  DBuffer: 2.0
  DFixed: 9.0
  FindDefects:
    Method: WS           # Wigner-Seitz — best for cascade damage detection
    ReferenceData: Fe_perfect.dat
    DCut4Def: 0.2
  PDReduction: True
  RecursiveRed: True     # Chain-reduce complex defect clusters
  DCut4PDR: 4.0
  Overlapping: False     # Merge nearby defects into one AV
  NMin4AV: 40
  NMax4AV: 5000

defect_bank:
  Recycle: True
  SaveDB: True
  UseSymm: True          # Store one representative per symmetry group

spsearch:
  NSearch: 15
  NMax4Trans: 1500
  TransHorizon: True
  LocalRelax:
    LocalRelax: True
  HandleVN:
    CenterVN: True
    RescaleVN: True
    RescaleStyle4LOGV: LOGNRAS
  force_evaluator:
    nproc: 2             # Processors per AV task

saddle_point:
  BarrierCut: 5.0
  BackBarrierMin: 0.0
  Prefactor: 10.0
  ValidSPs:
    CheckConnectivity: True
    RealtimeDelete: True

visual:
  Screen: True
  Log: True
  Write_Data_SPs:
    Write_KMC_Data: True
```

---

## Example 3: Using Symmetry to Find All Vacancy Saddle Points

For a single symmetric defect where you want a complete SP catalog.

```yaml
# ... (same potential/force_evaluator/data as Example 1) ...

kinetic_MC:
  NSteps: 50
  Temp: 600

active_volume:
  Style: defects
  DActive: 7.7
  DBuffer: 1.5
  DFixed: 7.7
  FindDefects:
    Method: BLCN
  PDReduction: True
  DCut4PDR: 5.0          # Larger cutoff to get single reduced defect
  PointGroupSymm: True   # AUTO-GENERATE symmetry-equivalent SPs
  NMax4PG: 500
  NMin4AV: 40

spsearch:
  NSearch: 5             # Fewer needed — symmetry fills in the rest
  NMax4Trans: 800

defect_bank:
  Recycle: True
  SaveDB: True
  UseSymm: True          # Save only unique SPs

# ... rest same as Example 1 ...
```

With `PointGroupSymm: True`, finding just 1 of the 8 first-NN vacancy jumps will
automatically generate all 8. This is much more efficient than running 10+ searches hoping
to find all variants.

---

## Example 4: Using the Mean Rate Method (MRM) for Low-Barrier Systems

For systems with fast rattling events (small barriers that trap conventional KMC).

```yaml
kinetic_MC:
  NSteps: 100
  Temp: 300
  AccStyle: MRM
  EnCut4Transient: 0.3   # Build superbasins for barriers < 0.3 eV
  NMaxBasin: 50          # Prevent memory overflow
  Handle_no_Backward: OUT

spsearch:
  NSearch: 10
  LocalRelax:
    LocalRelax: True     # REQUIRED for MRM
  HandleVN:
    CenterVN: True
    RescaleVN: True

force_evaluator:
  Relaxation:
    BoxRelax: False      # Box relaxation NOT allowed with MRM

# ... rest same as Example 1 ...
```

---

## Example 5: Guided SPS with Pre-Built Defect Database

Two-step workflow: generate SP database, then use it for subsequent runs.

### Step 1 — Generate database (`input_generate.yaml`)

```yaml
kinetic_MC:
  NSteps: 1
  Temp: 600

active_volume:
  PDReduction: True
  DCut4PDR: 5.0
  PointGroupSymm: True   # Find all symmetry-equivalent SPs

defect_bank:
  Recycle: True
  SaveDB: True
  SavePath: Fe_vacancy_DB
  UseSymm: False         # Save ALL SPs (not just unique)

spsearch:
  NSearch: 20            # High count to ensure complete catalog
```

### Step 2 — Use database (`input_guided.yaml`)

```yaml
defect_bank:
  LoadDB: True
  LoadPath: Fe_vacancy_DB
  Preload: True
  Ratio4DispLoad: 0.8

spsearch:
  NSearch: 10            # First 8 use loaded SPs, last 2 do fresh search
  Preloading:
    Preload: True
    LoadPath: Fe_vacancy_DB
    Ratio4DispLoad: 0.8
  TrialStepsizeFine: 0.003   # Smaller step for guided SPS (closer to SP)
  MaxStepsizeFine: 0.01
```

---

## Typical Barrier Values for Reference

| System | Defect | Barrier (eV) | Notes |
|--------|--------|--------------|-------|
| BCC Fe | Vacancy 1st-NN jump | 0.63–0.64 | Ackland-04 potential |
| BCC Fe | Dumbbell SIA | 0.31 | Ackland-04 potential |
| BCC Fe | C15 di-interstitial → ⟨111⟩ | 0.69 | Ackland-04 |
| BCC Fe | C15 tetra-interstitial → ⟨111⟩ | 0.89 | Ackland-04 |
| BCC W | He₃ interstitial cluster migration | 0.13–0.35 | Juslin-Wirth potential |
| BCC W | He₃ → vacancy complex formation | ~3 eV release | Exothermic |

---

## Common Parameter Combinations

### For accuracy (small defects, few atoms):
```
DActive: 7–8 Å, SNC: True, NSearch: 10, PointGroupSymm: True
```

### For large systems (loops, cascades):
```
DActive: 10–15 Å, RescaleVN: True, RescaleStyle4LOGV: LOGNRAS,
CenterVN: True, NSearch: 15–20, Recycle: True
```

### For high-throughput screening:
```
NSearch: 5, PointGroupSymm: True, Recycle: True, UseSymm: True
```

### For stuck simulations (low-barrier rattling):
```
AccStyle: MRM, EnCut4Transient: 0.3–0.5, LocalRelax: True
```
