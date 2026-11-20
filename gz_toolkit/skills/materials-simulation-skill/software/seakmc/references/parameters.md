# SEAKMC_py Complete Parameter Reference

All parameters from `input.yaml`. Units are LAMMPS metal units (Å, eV, ps, K) unless noted.
Default values shown in brackets.

---

## `system`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TempFiles` | `["tmp0.dat","tmp1.dat","tmp2.dat"]` | Internal temp file names (do not change) |
| `Inteval4ShowProgress` | 10 | Print progress every N tasks |
| `Restart.WriteRestart` | True | Write restart files |
| `Restart.LoadRestart` | True | Auto-load latest restart on start |
| `Restart.LoadFile` | False | Load specific restart file (string or False) |
| `Restart.AVStep4Restart` | 1000 | Write restart every N finished AVs |
| `Restart.KMCStep4Restart` | 1 | Write restart every N KMC steps |
| `Restart.Reset_Simulation_Time` | False | Reset simulation time to 0 on restart |
| `Tolerance` | 0.1 | General distance tolerance (Å) |
| `angle_tolerance` | 5 | Angle tolerance (degrees) |
| `significant_figures` | 6 | Decimal places for output coordinates |
| `float_precision` | 3 | Float precision for general output |
| `VerySmallNumber` | 1.0e-20 | Guard against division by zero |

---

## `kinetic_MC`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NSteps` | 1 | Number of KMC steps to run |
| `Temp` | 800 | Temperature for event selection (K) |
| `Temp4Time` | = Temp | Temperature for KMC time calculation (K) |
| `DispStyle` | FI | Update structure via "SP" or "FI" (final image) |
| `Sorting` | False | Sort event table by rate before KMC |
| `AccStyle` | NoAcc | Acceleration: "NoAcc" or "MRM" |
| `EnCut4Transient` | 0.5 | Max barrier for superbasin in MRM (eV) |
| `NMaxBasin` | NA (unlimited) | Max basins in MRM superbasin |
| `Handle_no_Backward` | OUT | If backward SP missing: "OUT" or "IN" (MRM) |
| `Tol4Disp` | 0.1 | Displacement tolerance for basin matching (Å, MRM) |
| `Tol4Barr` | 0.03 | Barrier tolerance for basin matching (eV, MRM) |

**MRM notes:** Mean Rate Method handles low-barrier rattling. Requires `spsearch.LocalRelax: True`. 
Cannot be used with box relaxation between KMC steps.

---

## `defect_bank`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LoadDB` | False | Load saddle point database from files |
| `Preload` | False | Preload displacements into AV-level search |
| `LoadPath` | DefectBank | Directory for loading DB |
| `Ratio4DispLoad` | 0.8 | Fraction of loaded displacements to use |
| `SortDisps` | False | Load displacements by filename sequence |
| `Recycle` | False | Recycle SPs from previous KMC steps |
| `SaveDB` | DefectBank | Save DB to this directory (or False) |
| `SavePath` | DefectBank | Directory for saving DB |
| `UseSymm` | False | Save only one SP per symmetry group |
| `FileHeader` | DB | Prefix for DB filenames |
| `NMax4DB` | 100 | Max atoms in a DB entry |
| `NMin4DB` | 8 | Min atoms in a DB entry |
| `Tol4Disp` | 0.1 | Tolerance for identifying same structures |
| `Scaling` | 1.0 | Scale factor applied to DB coordinates |
| `IgnoreType` | False | Ignore atom types when matching |

---

## `potential`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `species` | required | List of element symbols (one per LAMMPS type) |
| `pair_style` | required | LAMMPS pair_style (e.g., "eam/fs", "eam/alloy") |
| `FileName` | required | Potential file name |
| `Path2Pot` | False | Absolute path to potential file |
| `pair_coeff` | False | Override pair_coeff string (or False for auto-generate) |
| `masses` | auto | Atomic masses (from pymatgen if not specified) |
| `bondlengths` | pymatgen | Bond lengths per type pair — **must set for BCC metals** |
| `cutneighs` | 1.1×bondlengths | Neighbor cutoff distances |
| `coordnums` | 4 | Preferred coordination numbers — **must set for BCC/FCC** |
| `cutneighmax` | max(cutneighs) | Maximum neighbor cutoff |
| `bondlengths4LAS` | = bondlengths | Bond lengths for local atom strain (LAS) estimation |
| `coordnums4LAS` | = coordnums | Coordination numbers for LAS estimation |
| `OpenKIM.OpenKIM` | False | Use OpenKIM potential |
| `OpenKIM.kim_init` | False | KIM init string (required if OpenKIM=True) |
| `OpenKIM.kim_interaction` | False | KIM interaction string |

**BCC Fe example:**
```yaml
bondlengths:
  - 1 1 2.47   # ≈ 0.866 × 2.855 Å
coordnums:
  - 1 8
```

---

## `force_evaluator`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Bin` | pylammps | Force evaluator: "pylammps", "calllammps", "callvasp" |
| `Path2Bin` | False | Absolute path to binary/script |
| `Style` | = Bin | "pylammps", "lammps", or "vasp" |
| `nproc` | auto | Processors for data relaxation |
| `Nproc4Recal` | auto | Processors for energy recalibration |
| `NSteps4Relax` | 10000 | Max minimization steps (also MD steps) |
| `timestep` | 0.002 | MD timestep (ps) |
| `processors` | False | LAMMPS processors command |
| `Screen` | False | Print LAMMPS output to screen |
| `LogFile` | False | LAMMPS log file (False = no log) |
| `RinputOpt` | False | Custom input file for data relaxation |
| `RinputMD0` | False | Custom input file for energy recalibration |
| `ImportValue4RinputOpt` | False | Import SEAKMC values into RinputOpt |
| `Relaxation.BoxRelax` | False | Relax box shape between KMC steps |
| `Relaxation.InitTemp4Opt` | 0.0 | Initial temperature for NVT during relaxation |
| `Relaxation.TargetTemp4NVT` | 5.0 | Target temperature for NVT |
| `Relaxation.NVTSteps4Opt` | 10000 | NVT steps during relaxation |
| `GPU` | None | KOKKOS GPU arguments for PyLammps |

---

## `data`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `FileName` | required | LAMMPS data file name |
| `atom_style` | atomic | LAMMPS atom style |
| `Relaxed` | False | Is data already relaxed? Skip minimization if True |
| `BoxRelax` | False | Relax box during initial relaxation |
| `RinputOpt` | False | Custom input for initial relaxation |
| `MoleDyn` | False | Run MD on initial structure before KMC |
| `RinputMD` | False | Custom input for initial MD |
| `boundary` | p p p | Boundary conditions (p=periodic, f=fixed, s=shrink) |
| `units` | metal | Only metal units supported |
| `dimension` | 3 | Only 3D supported |

---

## `active_volume`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Style` | defects | AV type: "defects", "all", or "custom" |
| `DActive` | 4.1×cutdefectmax | Radius of active region (Å) |
| `DBuffer` | 0.4×cutdefectmax | Buffer layer thickness (Å) |
| `DFixed` | 4.1×cutdefectmax | Fixed layer thickness (Å) |
| `NMin4AV` | 40 | Min active atoms per AV |
| `NMax4AV` | unlimited | Max active atoms per AV |
| `NMax4Def` | unlimited | Max number of defects |
| `NPredef` | 0 | Number of preset defect types (via molecule-ID) |
| `PredefOnly` | False | Only use preset defects |
| `cutdefectmax` | max(cutneighs) | Max defect search cutoff |
| `FindDefects.Method` | BLCN | Defect detection: CN, BL, BLCN, WS, all, custom |
| `FindDefects.DCut4Def` | 0.1×BL | Deviation threshold for defect detection |
| `FindDefects.DiscardType` | [] | Atom types to exclude from defect search |
| `FindDefects.Defects` | — | Manual defect coordinates (if Method=custom) |
| `FindDefects.ReferenceData` | — | Reference data file (if Method=WS) |
| `PDReduction` | True | Reduce defect groups to point defects |
| `DCut4PDR` | 1.4×cutdefectmax | Distance cutoff for point defect reduction |
| `RecursiveRed` | False | Recursively reduce chained defects |
| `SortD4PDR` | False | Sort defects by distance before reduction |
| `Overlapping` | True | Allow AVs to overlap |
| `DCut4noOverlap` | 2×(DActive+DBuffer) | Non-overlap merge distance |
| `Sorting` | True | Sort atoms in AV by distance from defect center |
| `Sort_by` | Auto → [D,X,Y,Z] | Sorting keys (D=distance, X/Y/Z=rounded coords) |
| `SortingSpacer` | [0.3,0.3,0.3] | Spacing for coordinate rounding in sorting |
| `PointGroupSymm` | False | Generate symmetry-equivalent SPs automatically |
| `NMax4PG` | 1000 | Max atoms for point group symmetry analysis |
| `TurnoffPBC` | [False,False,False] | Turn off PBC per direction for AVs |
| `RT_SetMolID` | False | Real-time set molecule-IDs around defect center |
| `R4RT_SetMolID` | 30.0 | Radius for RT_SetMolID (Å) |

**Convergence test for DActive:**
Typical values by defect type (BCC Fe, a₀=2.855 Å):
- Single vacancy: DActive ≈ 7.7 Å (168 atoms)
- Dumbbell SIA: DActive ≈ 12 Å (610 atoms)
- Small cluster (≤5 SIAs): DActive ≈ 12–15 Å
- Interstitial loop (37 SIAs): DActive ≈ 12 Å + merge all defect AVs

---

## `dynamic_matrix`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `SNC` | False | Scaled Normal Coordinates mode (dramatic speedup for small AVs) |
| `CalPrefactor` | False | Calculate vibrational prefactors (requires PHONON) |
| `OutDynMat` | False | Write dynamic matrix to file |
| `NMax4SNC` | 1000 | Max atoms for SNC (memory scales as 9N²) |
| `displacement` | 1e-6 | Displacement for finite-difference dynamic matrix |
| `Method4Prefactor` | harmonic | Only "harmonic" (Vineyard formula) available |
| `VibCut` | 1e-8 | Vibrational frequency cutoff for prefactor (ignore below this) |

**SNC guidance:** Reduces translation steps from ~200 to ~23 for simple defects. Memory
scales as 9N² and time as N³ — avoid for AVs > ~300 atoms.

---

## `spsearch`

### Main search settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Method` | dimer | Only "dimer" supported |
| `NSearch` | 10 | Search attempts per AV |
| `SearchBuffer` | False | Include buffer atoms in final refinement |
| `NMax4Trans` | 1000 | Max translation steps per attempt |
| `FConv` | 1e-6 | Force convergence criterion (eV/Å) |
| `EnConv` | 1e-5 | Energy convergence criterion (eV) |
| `DimerSep` | 0.005 | Dimer image separation (Å) |
| `DRatio4Relax` | 2.0 | Ratio × DimerSep for pre-relaxation tweak |
| `Tol4Connect` | 0.1 | MAXABSC tolerance for connectivity check (Å) |

### Translation step size

| Parameter | Default | Description |
|-----------|---------|-------------|
| `TrialStepsize` | 0.015 | Trial translation step (Å) |
| `MaxStepsize` | 0.05 | Maximum translation step (Å) — must be > Trial |
| `MinStepsize` | 0.003 | Minimum translation step (Å) |
| `DecayStyle` | fixed | Step decay: "fixed", "staircase", "continuous" |
| `DecayRate` | 0.71 | Base number for decay power function |
| `DecaySteps` | 20 | Steps per decay interval |
| `TrialStepsizeFine` | 0.003 | Trial step for fine-tuning guided SPS |
| `MaxStepsizeFine` | 0.0 | Max step for fine-tuning (0 = same as max) |

### Rotation settings

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NMax4Rot` | 3 | Max dimer rotation steps per translation step |
| `FThres4Rot` | 0.1 | Force threshold to trigger rotation |
| `FMin4Rot` | 0.001 | Rotation force convergence |

### Early termination

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CheckAng` | True | Stop early if angle to known SP is < AngCut |
| `CheckAngSteps` | 50 | Check angle every N steps |
| `AngCut` | 2.0 | Angle threshold for early termination (degrees) |
| `IgnoreSteps` | 4 | Skip angle check for first N steps |
| `TransHorizon` | True | Translate dimer horizontally when trapped |
| `En4TransHorizon` | 0.1 | Energy threshold for TransHorizon (eV) |

### Dynamic Active Volume (DAV) — `HandleVN`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CenterVN` | True | Remove net translation from dimer axis |
| `NSteps4CenterVN` | 1 | Apply CenterVN every N steps |
| `RescaleVN` | True | Enable DAV rescaling |
| `RescaleStyle4LOGV` | LOGNRAS | Rescaling mode: LOGN, RAS, LOGNRAS |
| `XRange4LOGV` | 20 | Range for LOGN sigmoid function |
| `Ratio4Zero4LOGV` | 0.2 | Zero point ratio for LOGN (deactivate atoms below 4.8% of max displacement) |
| `XRange4RAS` | 40 | Range for RAS sigmoid function |
| `Ratio4Zero4RAS` | 0.3 | Zero point ratio for RAS (deactivate outer 30% by sequence) |
| `TakeMin4MixedRescales` | True | Take minimum of LOGN and RAS for LOGNRAS |
| `PowerOnV` | 4.0 | Power applied to displacement before LOGN |
| `Period4MA` | 1 | Moving average period for displacement smoothing |

**LOGNRAS** (default) is the best general-purpose choice: LOGN deactivates atoms with
small displacements (dynamic), RAS deactivates atoms near the AV boundary (static).

### Local relaxation after each SPS

| Parameter | Default | Description |
|-----------|---------|-------------|
| `LocalRelax.LocalRelax` | True | Relax within AV after KMC event |
| `LocalRelax.InitTemp4Opt` | 0.0 | Initial NVT temperature for local relax |
| `LocalRelax.TargetTemp4NVT` | 5.0 | NVT target temperature |
| `LocalRelax.NVTSteps4Opt` | 10000 | NVT steps |

### Preloading (guided SPS)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Preloading.Preload` | False | Load saddle points to guide search |
| `Preloading.LoadPath` | — | Path to preloaded SP files |
| `Preloading.Method` | — | Loading method |
| `Preloading.Ratio4DispLoad` | 0.8 | Fraction of displacement magnitude to use |
| `Preloading.Scaling` | 1.0 | Scale factor for loaded displacements |
| `Preloading.SortDisps` | False | Sort by filename sequence |

---

## `saddle_point`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `BarrierCut` | 5.0 | Discard SPs with forward barrier > this (eV) |
| `BarrierMin` | 0.0 | Discard SPs with forward barrier < this (eV) |
| `EbiasCut` | 5.0 | Discard SPs with |ebias| > this (eV) |
| `EbiasMin` | -5.0 | Discard SPs with ebias < this (eV) |
| `BackBarrierMin` | 0.0 | Min backward barrier (eV) |
| `Prefactor` | 1.0 | Default prefactor if CalPrefactor=False (THz) |
| `DAtomCut` | unlimited | Max number of displaced atoms (nad) |
| `DMagCut` | unlimited | Max displacement magnitude from init→SP (Å) |
| `DtotCut` | unlimited | Max total displacement (Å) |
| `DmaxCut` | unlimited | Max single-atom displacement init→SP (Å) |
| `ValidSPs.CheckConnectivity` | True | Only keep SPs that connect two minima |
| `ValidSPs.RealtimeDelete` | True | Delete invalid SPs immediately (saves memory) |
| `ValidSPs.RealtimeValid` | True | Validate SPs in real time |
| `ValidSPs.NScreenDisp` | 0 | Number of displacement-based screens |
| `ValidSPs.NScreenEng` | 0 | Number of energy-based screens |
| `GroupSP` | False | Group similar SPs as one event |
| `AngCut4GSP` | 30.0 | Angle cutoff for SP grouping |
| `MagCut4GSP` | 0.2 | Magnitude cutoff for SP grouping |
| `EnCut4GSP` | 0.03 | Energy cutoff for SP grouping |

---

## `visual`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `Screen` | True | Print to screen |
| `Log` | True | Write Seakmc.log |
| `Write_SP_Summary` | True | Write SPOut/KMC_N_SPs.csv |
| `Write_Data_SPs.Write_KMC_Data` | True | Write DataOut/KMC_N.dat |
| `Write_Data_SPs.Write_Data_SPs` | True | Write SP structure files |
| `Write_Data_SPs.Write_KMC_Data` | True | Write structure at each KMC step |
| `Write_Data_SPs.Write_Prob` | True | Write probability distribution |
| `DataOutPath` | DataOut | Directory for structure output |
| `SPOutPath` | SPOut | Directory for SP output |
| `Write_AV_SPs.Write_AV_SPs` | False | Write AV-level SP files |
| `Write_AV_SPs.Write_Local_AV` | False | Write AV structure files |
| `OutputStyle` | individual | "individual" or "stack" output format |
| `ShowBuffer` | False | Include buffer atoms in output |
| `ShowFixed` | False | Include fixed atoms in output |
