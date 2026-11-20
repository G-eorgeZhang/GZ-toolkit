# VASP INCAR Templates by Calculation Type

All templates are starting points. Always verify ENCUT against your POTCAR ENMAX values and
converge ENCUT and k-points for your specific system.

---

## 1. Single-Point Energy (SCF only)

```
System = Single-point
ISTART = 0       # Start from scratch
ICHARG = 2       # Superposition of atomic charge densities
ENCUT  = 520     # ≥ 1.3 × max ENMAX in POTCAR
EDIFF  = 1E-6    # SCF convergence criterion (eV)
PREC   = Accurate
ALGO   = Normal
ISMEAR = 0       # Gaussian smearing (insulator/molecule); use 1 for metal
SIGMA  = 0.05
NSW    = 0       # No ionic steps
IBRION = -1
LWAVE  = .FALSE.
LCHARG = .FALSE.
```

---

## 2. Geometry Relaxation — Bulk

```
System = Bulk relaxation
ISTART = 0
ICHARG = 2
ENCUT  = 520
EDIFF  = 1E-6
EDIFFG = -0.02   # Force convergence (eV/Å), negative = force criterion
PREC   = Accurate
ALGO   = Normal
ISMEAR = 1       # Methfessel-Paxton for metals; 0 for insulators
SIGMA  = 0.1
NSW    = 200     # Max ionic steps
IBRION = 2       # Conjugate gradient
ISIF   = 3       # Relax ions + cell shape + volume
LWAVE  = .FALSE.
LCHARG = .FALSE.
```

---

## 3. Geometry Relaxation — Surface Slab

```
System = Surface slab relaxation
ISTART = 0
ICHARG = 2
ENCUT  = 520
EDIFF  = 1E-6
EDIFFG = -0.02
PREC   = Accurate
ALGO   = Normal
ISMEAR = 1       # or 0 for semiconductor surface
SIGMA  = 0.1
NSW    = 200
IBRION = 2
ISIF   = 2       # Relax ions ONLY — do NOT use ISIF=3 for slabs
LWAVE  = .FALSE.
LCHARG = .FALSE.
# Add LDIPOL/IDIPOL for polar slabs:
# IDIPOL = 3
# LDIPOL = .TRUE.
```

> ⚠️ ISIF=3 on a slab will relax the vacuum dimension — almost always wrong.
> Ref: https://vasp.at/wiki/ISIF

---

## 4. Band Structure

Requires two steps:

**Step 1 — SCF to get converged CHGCAR:**
```
ISTART = 0
ICHARG = 2
ENCUT  = 520
EDIFF  = 1E-6
PREC   = Accurate
ALGO   = Normal
ISMEAR = 0
SIGMA  = 0.05
NSW    = 0
IBRION = -1
LWAVE  = .FALSE.
LCHARG = .TRUE.   # Save charge density for step 2
```

**Step 2 — Non-SCF band structure:**
```
ISTART = 1
ICHARG = 11      # Read CHGCAR, keep fixed
ENCUT  = 520
EDIFF  = 1E-6
PREC   = Accurate
ALGO   = Normal
ISMEAR = 0
SIGMA  = 0.05
NSW    = 0
IBRION = -1
LORBIT = 11      # Projected DOS / fat bands
LWAVE  = .FALSE.
LCHARG = .FALSE.
```

KPOINTS for step 2: use explicit k-path (line mode).
Use pymatgen `HighSymmKpath` or VASPKIT option 302 to generate.
Ref: https://vasp.at/wiki/Band-structure_calculation_using_DFT

---

## 5. Density of States (DOS)

```
# After SCF (reuse CHGCAR), run non-SCF:
ISTART = 1
ICHARG = 11
ENCUT  = 520
EDIFF  = 1E-6
PREC   = Accurate
ALGO   = Normal
ISMEAR = -5      # Tetrahedron method — best for DOS; needs dense k-mesh
SIGMA  = 0.05
NSW    = 0
IBRION = -1
LORBIT = 11      # l-decomposed DOS
NEDOS  = 2000    # DOS mesh points
EMIN   = -10     # Energy range (relative to Fermi)
EMAX   = 10
LWAVE  = .FALSE.
LCHARG = .FALSE.
```

> ⚠️ ISMEAR=-5 requires a well-converged k-mesh (≥ 4×4×4 minimum).
> Ref: https://vasp.at/wiki/ISMEAR

---

## 6. DFT+U (LDA+U / GGA+U)

```
# Add to any base INCAR:
LDAU    = .TRUE.
LDAUTYPE = 2     # Dudarev simplified rotational invariant
LDAUL   = 2 -1  # L quantum number per species: 2=d, -1=no U
LDAUU   = 3.5 0  # U values (eV) per species
LDAUJ   = 0.0 0  # J values (eV) — usually 0 for Dudarev
LDAUPRINT = 1    # Print occupancy matrix
LMAXMIX = 4      # Needed for d electrons (6 for f)
```

> ℹ️ LDAUL, LDAUU, LDAUJ entries must match species order in POSCAR/POTCAR.
> Choosing U: consult literature for your system, or use linear response method.
> Ref: https://vasp.at/wiki/LDAU | https://vasp.at/wiki/DFT%2BU:_formalism

---

## 7. Hybrid Functional (HSE06)

```
LHFCALC = .TRUE.
HFSCREEN = 0.2   # HSE06 screening parameter (Å⁻¹)
AEXX    = 0.25   # Fraction of exact exchange (0.25 = HSE06)
ALGO    = All    # Required for hybrids — do NOT use Fast
TIME    = 0.4    # Typical for All algorithm
PRECFOCK = Normal
ENCUT   = 520
EDIFF   = 1E-5
PREC    = Accurate
ISMEAR  = 0
SIGMA   = 0.05
NSW     = 0
```

> ⚠️ Hybrids are 5-10× more expensive than GGA. Pre-converge with GGA first.
> Ref: https://vasp.at/wiki/LHFCALC | https://vasp.at/wiki/HSE06

---

## 8. vdW Correction (DFT-D3)

```
# Grimme D3 with BJ damping (recommended):
IVDW   = 12
# D3 with zero damping:
# IVDW = 11
```

> Add to any GGA INCAR. Compatible with ISMEAR, NSW, IBRION settings.
> Ref: https://vasp.at/wiki/DFT-D3

---

## 9. Spin-Polarised Calculation

```
ISPIN  = 2
MAGMOM = 5 5 0 0  # Initial magnetic moments per atom (match POSCAR count)
# For ferromagnetic Fe bcc (2 atoms): MAGMOM = 2*5.0
# For AFM: alternate signs, e.g., 5 -5
NUPDOWN = -1      # -1 = no constraint on total spin
```

> After SCF, check OUTCAR for "magnetization (x)" to see converged moments.
> Ref: https://vasp.at/wiki/ISPIN | https://vasp.at/wiki/MAGMOM

---

## 10. Noncollinear Magnetism / Spin-Orbit Coupling

```
LSORBIT = .TRUE.
LNONCOLLINEAR = .TRUE.
ISPIN  = 1          # Ignored when LNONCOLLINEAR=.TRUE.
MAGMOM = 0 0 5  0 0 5  # 3 components per atom: mx my mz
SAXIS  = 0 0 1       # Spin quantisation axis
GGA_COMPAT = .FALSE.
LMAXMIX = 4
```

> Ref: https://vasp.at/wiki/LSORBIT | https://vasp.at/wiki/Noncollinear_calculations

---

## 11. NEB (Nudged Elastic Band)

```
# Run in a directory with subdirs 00/ 01/ ... NN/ (images)
ICHAIN = 0    # NEB
IMAGES = 7    # Number of intermediate images
SPRING = -5   # Spring constant (eV/Å²)
IBRION = 3    # Quickmin or LBFGS
POTIM  = 0    # Step size (use 0 with IOPT)
IOPT   = 1    # LBFGS ionic optimizer
EDIFF  = 1E-6
EDIFFG = -0.05
NSW    = 200
LCLIMB = .TRUE.  # CI-NEB for accurate barrier
```

> NEB requires VTST tools (Henkelman group) patched into VASP.
> Ref: https://vasp.at/wiki/NEB | https://theory.cm.utexas.edu/vtsttools/neb.html

---

## 12. Ab Initio Molecular Dynamics (AIMD)

```
IBRION = 0       # MD
MDALGO = 2       # Nosé-Hoover thermostat
NSW    = 5000    # Total MD steps
POTIM  = 2.0     # Timestep in femtoseconds
SMASS  = 0       # Nosé mass parameter (0 = automatic)
TEBEG  = 300     # Starting temperature (K)
TEEND  = 300     # Final temperature (K)
ISYM   = 0       # Disable symmetry for MD
LWAVE  = .FALSE.
LCHARG = .FALSE.
NWRITE = 1       # Write OUTCAR every step
ISMEAR = 0
SIGMA  = 0.1
```

> Ref: https://vasp.at/wiki/IBRION | https://vasp.at/wiki/Molecular_dynamics

---

## 13. Phonons via DFPT (IBRION=8)

```
IBRION = 8       # DFPT linear response
NSW    = 1
EDIFF  = 1E-8    # Very tight for accurate force constants
PREC   = Accurate
ENCUT  = 520
ISMEAR = 0
SIGMA  = 0.01
LEPSILON = .TRUE.  # Also compute dielectric tensor (optional)
```

> Ref: https://vasp.at/wiki/IBRION | https://vasp.at/wiki/Phonons:_Theory

---

## 14. Optical Properties / Dielectric Function

```
ALGO   = Exact   # or CHI for RPA
LOPTICS = .TRUE.
NEDOS  = 2000
NBANDS = <2× default>  # Need many unoccupied bands
CSHIFT = 0.1           # Broadening (eV)
ISMEAR = 0
SIGMA  = 0.01
```

> Ref: https://vasp.at/wiki/LOPTICS | https://vasp.at/wiki/Dielectric_properties

---

## Quick Tag Wiki Index

| Tag | Wiki |
|-----|------|
| ALGO | https://vasp.at/wiki/ALGO |
| EDIFF | https://vasp.at/wiki/EDIFF |
| EDIFFG | https://vasp.at/wiki/EDIFFG |
| ENCUT | https://vasp.at/wiki/ENCUT |
| IBRION | https://vasp.at/wiki/IBRION |
| ICHARG | https://vasp.at/wiki/ICHARG |
| IDIPOL | https://vasp.at/wiki/IDIPOL |
| ISIF | https://vasp.at/wiki/ISIF |
| ISMEAR | https://vasp.at/wiki/ISMEAR |
| ISPIN | https://vasp.at/wiki/ISPIN |
| IVDW | https://vasp.at/wiki/IVDW |
| KPAR | https://vasp.at/wiki/KPAR |
| LDAU | https://vasp.at/wiki/LDAU |
| LHFCALC | https://vasp.at/wiki/LHFCALC |
| LORBIT | https://vasp.at/wiki/LORBIT |
| LSORBIT | https://vasp.at/wiki/LSORBIT |
| LWAVE | https://vasp.at/wiki/LWAVE |
| MAGMOM | https://vasp.at/wiki/MAGMOM |
| NBANDS | https://vasp.at/wiki/NBANDS |
| NELM | https://vasp.at/wiki/NELM |
| NPAR | https://vasp.at/wiki/NPAR |
| NSW | https://vasp.at/wiki/NSW |
| POTIM | https://vasp.at/wiki/POTIM |
| PREC | https://vasp.at/wiki/PREC |
| SIGMA | https://vasp.at/wiki/SIGMA |
