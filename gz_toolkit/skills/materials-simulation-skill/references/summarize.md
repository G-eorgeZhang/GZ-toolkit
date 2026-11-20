# Step 7 — Summarizing results

Output is a CSV (default `summary.csv` in ROOT): one row per run dir with every tree variable,
plus one column per requested value. `--average` also writes `summary_avg.csv` (mean, std, n per
variable combination). Ask the user: which quantity, from which file, which occurrence
(last / first / min / max / mean / all), and whether unfinished runs should be included.

## Field types (`sim_plan.yaml: summary.fields`)

```yaml
summary:
  output: summary.csv
  fields:
    # any text file, first capture group
    - {name: E_bind, type: regex, file: log.lammps, pattern: "^Ebind\\s*=\\s*(\\S+)", which: last}
    # a LAMMPS thermo column; run: which thermo block (-1 = last run/minimize in the log)
    - {name: PE, type: lammps_thermo, file: log.lammps, column: PotEng, which: last, run: -1}
    # a CSV column (SEAKMC summaries, fix ave/time outputs converted to CSV, ...)
    - {name: Ea_min, type: csv_column, file: Seakmc_summary.csv, column: barrier, which: min}
    # VASP: last E0= in OSZICAR (key: "F=" for free energy)
    - {name: E0, type: vasp_energy, file: OSZICAR}
    - {name: n_lines, type: file_lines, file: dump.final}
    - {name: n_kmc, type: count_files, file: "DataOut/KMC_*.dat"}
    # OVITO / any custom Python: func(run_dir: Path, field: dict) -> value
    - {name: n_vac, type: ovito, file: dump.final, script: "reffiles/ovito_extract.py:ws_vacancies",
       reference: reffiles/ref0_perfect.data}
```

Globs are allowed in `file`; the last match (sorted) is used — fine for `dump.*` patterns
where the last is the final frame, but check the sort is right for unpadded numbers
(`dump.100` sorts before `dump.20`). If it isn't, use a custom Python field.

Make the LAMMPS input print what you need in an easy-to-grep form, e.g.
`print "Ebind = ${Eb}"` — this is often cleaner than parsing thermo.

## Per-code quick recipes

| Code | Typical quantities | Where |
|------|--------------------|-------|
| LAMMPS | final PotEng, `print` lines, averaged fix output | `log.lammps`, fix ave/time files, dumps (OVITO) |
| SEAKMC | barriers, KMC time, defect positions | `Seakmc_summary.csv`, `SPOut/KMC_N_SPs.csv`, `DataOut/KMC_N.dat` |
| VASP | E0, magnetization, forces, final volume | `OSZICAR`, `OUTCAR`, `vasprun.xml` (pymatgen/py4vasp in a custom field) |
| SPPARKS | stats columns over time, final grain/cluster stats | `log.spparks`, diag output, dumps |
| MOOSE | postprocessor values | `<input>_out.csv` (csv=true) → `csv_column` |

## OVITO

`pip install ovito` into the workflow environment (check the import works on the target
machine). Typical analyses: Wigner–Seitz (vacancies/interstitials vs. a perfect reference),
CNA / PTM (structure fractions), DXA (dislocation length, Burgers vectors), cluster analysis.
Examples in `templates/ovito_extract_example.py`. OVITO's Python API changes between major
versions — check names at https://docs.ovito.org/python/ for the installed version.

OVITO extraction is slow per file: run `summarize.py --workers N` (or as a one-node job on the
cluster).

## Analysis (step 8, optional)

From `summary.csv`/`summary_avg.csv`: plots (matplotlib), fits (Arrhenius `ln k` vs `1/T` for
rates, binding energy vs distance from the core), tables for papers. Put scripts in
`ROOT/analysis/` and make them read the CSV, so they re-run as more runs finish.
