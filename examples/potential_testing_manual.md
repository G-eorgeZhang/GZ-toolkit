# Potential Testing — Step-by-Step Manual

This is a practical, run-it-yourself guide to `gz_toolkit.potential_testing`:
scaffold a project, fill in one JSON per potential, submit to an HPC cluster,
and collect a comparison table across potentials.

CLI entry point: `gz-toolkit-potential-testing` (see `cli.py`).

---

## Step 1 — Create the project

```bash
cd /path/to/your/hpc/workdir
gz-toolkit-potential-testing init --root . --pot-name my_potential --interactive
```

This scaffolds:

```
<root>/
├── main.py                 # calls run_all() — builds references + job scripts
├── summarize.py             # calls write_grouped_summaries() — writes summary CSV(s), see Step 10
├── check.py                 # user-run pre-flight check — see Step 4
├── pot_inputs/
│   └── my_potential.json    # the config you'll edit in Step 2
└── potentials/
    └── my_potential/        # drop your potential file(s) here
```

Add one more `pot_inputs/<name>.json` per potential you want to compare
side-by-side (they all get summarized into the same `summary.csv`).

**Skipped `init` and hand-wrote `pot_inputs/*.json` instead?** You won't have
`main.py`/`summarize.py`/`check.py` yet — drop them in without touching
anything else you've already set up:

```bash
gz-toolkit-potential-testing write-scripts --root .
```

It's a no-op for any of the three files that already exist.

---

## Step 2 — Fill in `pot_inputs/<name>.json`

The file has three top-level sections: `potential`, `workflow`, `hpc`.
Every field below shows a real example value and what to watch out for.

### `potential` section

| Tag | Example | Notes |
|---|---|---|
| `pot_name` | `"my_potential"` | Must match the JSON filename and the `potentials/<pot_name>/` folder. |
| `pot_lines` | `"pair_style eam/fs\npair_coeff * * Fe.eam.fs Fe"` | **Multi-line tag.** In JSON this is one string with `\n` between LAMMPS lines — it gets written verbatim into `potential.inc` and `#include`d by every generated input script. For MEAM (2 pair_coeff-style lines) or hybrid potentials just add more `\n`-joined lines, e.g. `"pair_style hybrid/overlay eam/fs meam\npair_coeff * * eam/fs Fe.eam.fs Fe\npair_coeff * * meam library.meam Fe FeCr.meam Fe"`. |
| `type_map` | `{"Fe": 1, "Cr": 2, "He": 3}` | Element symbol → LAMMPS atom-type id. Every element used anywhere else in the config (elements, gases, alloys) must appear here. |
| `masses` | `{"1": 55.845, "2": 51.996, "3": 4.0026}` | Keys are the **atom-type ids as strings** (JSON object keys are always strings; loaded back as `int`), values are amu. Must have exactly one entry per `type_map` value. |
| `single_elements` | `["Fe", "Cr"]` | Pure-element reference cells to build (box-relax to find `lc`/`Ecoh`). |
| `gases` | `["He"]` | Interstitial gas species, if any. Leave `[]` if none. |
| `gas_max_n` / `gas_max_m` | `4`, `4` | Max gas atoms / max vacancies-or-interstitials per gas-defect complex (only used if `include_gas_complexes` is on). |
| `gas_complex_pairs` | `[["Fe", "He"]]` | Optional — restrict which (metal, gas) pairs get complexes built. Leave `[]` to default to every `single_elements` × `gases` combination. |
| `gas_with_vacancy` / `gas_with_interstitial` | `true`, `true` | Whether to build `(gas)n(V)m` and/or `(gas)n(I)m` complexes. |
| `two_element_suites` | `[{"A": "Fe", "B": "Cr", "fractions_atpct": [3, 5, 10, 50], "ordering": "random"}]` | Binary alloy compositions to **scan** (a whole curve of B-fractions in an A matrix). `ordering` can be `"random"`, `"B2"`, or `"random\|B2"` (B2 only applies at 50 at.%, falls back to random otherwise). Add one dict per element pair. |
| `multi_element_suites` | `[{"composition": {"Fe": null, "Cr": 3, "Ni": 2}, "ordering": "random"}]` | Fixed-point alloy compositions for **3 or more elements** (also works for 2). `composition` maps element → at.%; use `null` for "whatever's left" — any number of `null` entries split the remainder evenly (e.g. `{"Fe": null, "Cr": null, "Ni": 6}` gives Fe = Cr = 47). Explicit values with no `null` must sum to 100. The **first key is the host lattice species** used to seed/replicate the cell. Only `"random"` ordering is supported (no B2 for 3+ elements). Add one dict per composition point you want tested — unlike `two_element_suites` there's no `fractions_atpct` scan, each entry is one alloy. |
| `crystal_structures` | `{"Fe": {"structure": "bcc"}, "Cr": {"structure": "bcc"}}` | One entry per `single_elements` member. Supported: `"bcc"`, `"fcc"` (hcp not yet wired into defect placement). |
| `lc_initial` | `2.85` | Starting lattice constant (Å) fed into the reference box-relax; gets refined per-element by the LAMMPS run. |
| `size_single` | `20` | N×N×N supercell replication for reference/defect cells. |
| `size_alloy` | `20` | N×N×N replication for alloy composition cells. |
| `size_elastic` | `4` | Smaller replication for elastic-constant cells (Cij converges fast). |
| `size_elastic_alloy` | `6` | Replication for alloy elastic cells (only if `elastic_for_alloys` is on) — a single random realization, treat Cij as an estimate. |
| `potential_files` | `["library.meam", "FeCr.meam"]` | **List, not a single string** — a potential can need more than one file (e.g. MEAM's shared library + per-system parameter file, or a `.mtp` file plus a settings file). Every name listed here must exist under `potentials/<pot_name>/`; each is copied into every case directory. Use `[]` if the potential is self-contained in `pot_lines` (e.g. a built-in pair_style needing no external file). |
| `kind` | `"classical"` | One of `"classical"` (eam/meam/hybrid), `"ml"` (mtp/snap), `"universal"` (chgnet/mace). Affects SEAKMC-compatibility checks. |
| `promote_to_infobank` | `false` | Whether `summarize-pot` (Step 9) also copies this pot's `<pot_name>.json` into `gz_toolkit/pot_infobank/`. `false` (default) never copies. `true` copies to the cluster's `path2gz_toolkit` — requires `hpc.cluster` to be set and that cluster's `gz_toolkit/jobs/cluster_info/<CLUSTER>/meta.json` to have `path2gz_toolkit` filled in (it's per-cluster, not per-pot, since the install location doesn't change across pots on the same machine). A path string copies there directly, no validation. |

### `workflow` section

| Tag | Example | Notes |
|---|---|---|
| `run_find_lc` / `run_reference_supercell` / `run_defect_generation` / `run_energies` / `run_summary` | `true` (all) | Legacy per-stage on/off switches, kept for the single-pot `reference`/`defects` CLI subcommands. Leave `true` for the normal multi-pot flow. |
| `seakmc_enabled` | `false` | Turn on to compute migration barriers via SEAKMC for the defects listed in `seakmc_defects`. |
| `seakmc_temp_K` | `600.0` | KMC temperature. |
| `seakmc_template_path` | `null` | Path to a custom `input.yaml` template; `null` uses the built-in default. |
| `seakmc_defects` | `["1vac", "1int_tetra", "1int_octa"]` | Must be a subset of `defect_catalog` — `validate` will flag it otherwise. |
| `defect_catalog` | `["1vac", "2vac1nn", ..., "3int"]` | Which point defects get built for every element in `single_elements`. Defaults to the full built-in catalog; trim this list to skip specific defects. |
| `include_dislocation_loops` | `false` | BCC-only self-interstitial loops (`SIL111`, `SIL100`). |
| `include_dislocation_lines` | `false` | Not implemented yet — leave `false` (raises `NotImplementedError` if a case tries to build one). |
| `include_alloy_suite` | `false` | Turn on together with a non-empty `two_element_suites` and/or `multi_element_suites`. Builds the lattice-constant scan only (`alloy_lc/`) — no defect energies. |
| `include_alloy_defects` | `false` | Point-defect formation energies for every composition in `two_element_suites`/`multi_element_suites`. Separate opt-in from `include_alloy_suite` — every alloy build is an independently randomized solute arrangement (see the `structure_ops.py` module docstring: seed at lc=1 → apply the defect → scale to the real lc → *then* randomly decide chemistry, so an inserted atom is drawn from the same composition as everyone else), so this runs `alloy_defect_replicas` independent replicas per (composition, case) and averages the resulting `Ef` — meaningfully more expensive than the lc scan above. |
| `alloy_defect_replicas` | `3` | Independent random realizations averaged per (composition, case) under `include_alloy_defects`, to smooth out configurational noise from where the solutes happened to land relative to the defect. |
| `include_gas_complexes` | `false` | Turn on together with non-empty `gases`. |
| `include_elastic` | `true` | Full 6×6 Cij for every pure element (on by default). |
| `elastic_for_alloys` | `false` | Also compute Cij per alloy composition — requires `two_element_suites` and/or `multi_element_suites` non-empty. |
| `elastic_strain` | `1.0e-6` | Finite-deformation magnitude (LAMMPS ELASTIC `up` variable). |
| `dumbbell_sep_factor` | `0.35` | Dumbbell interstitial separation = factor × lc. |
| `loop_radius_factor` | `1.5` | SIL loop radius = factor × lc. |
| `minimize_tol` | `1.0e-18` | etol/ftol for every generated `minimize` command. |
| `paral_degree` | `4` | `1` = one job for everything, `2` = one job per potential, `3` = one job per (potential, work-group) with `--dependency=afterok` chaining, `4` = one job per simulation (scatter — fastest, most queue slots). |

### `hpc` section

| Tag | Example | Notes |
|---|---|---|
| `mode` | `"scatter"` | `"scatter"` (one `.job` per case, pairs with `paral_degree: 4`) or `"batch"` (group cases into fewer `.job` files, capped by `batch_max_per_file`/`max_jobs_per_batch`). |
| `cluster` | `"ISAAC"` | Must match a CSV in `gz_toolkit/jobs/cluster_info/<CLUSTER>/<CLUSTER>.csv` and/or a template in `gz_toolkit/jobs/templates/machines/<CLUSTER>.job`. Use `null` for a generic SLURM script. Each cluster's directory can also hold `<CLUSTER>/meta.json` with `path2gz_toolkit` (see `promote_to_infobank` above). |
| `partition_key` | `"s"` | **Must be quoted** — it's a string key into the cluster CSV, not a bareword. For ISAAC: `"c"` = campus (24:00:00 cap), `"s"` = short (3:00:00 cap), `"g"` = condo (30-00:00:00 cap). Picking a key auto-fills partition/qos/account and caps `walltime` to that partition's max. |
| `walltime` | `"03:00:00"` | `HH:MM:SS` or `D-HH:MM:SS`. Must not exceed the chosen partition's cap (it gets silently capped with a warning if it does). |
| `nodes` | `1` | Number of nodes per job. |
| `ntasks` | `48` | Total MPI ranks. For a full-node job, match the partition's `cores_per_node` from the cluster CSV (e.g. ISAAC `short` = 48 cores/node). |
| `output_name` | `"run.job"` | Base filename for generated job scripts. |
| `batch_max_per_file` | `50` | Only used when `mode: "batch"`. |
| `max_jobs_per_batch` | `200` | Only used when `mode: "batch"` — work groups exceeding this split into `_1`, `_2`, ... files. |
| `run_command` | `"srun -n {{NTASKS}} {{LAMMPS_BIN}} -in {{INPUT}}"` | Placeholders: `{{NTASKS}}` → `hpc.ntasks`, `{{INPUT}}` → the case's input filename (e.g. `in.reference.lammps`), `{{LAMMPS_BIN}}` → `lammps_binary_path` (or plain `lmp_mpi` if unset). **This value can be multiple lines** (`\n`-joined) if you need extra shell commands before the run — see the MTP example below. |
| `device` | `"cpu"` | `"cpu"` or `"gpu"`. `"ml"`/`"universal"` potentials are usually GPU. |
| `gpus_per_node` | `0` | If `device: "gpu"` and this is `> 0`, a `#SBATCH --gres=gpu:<N>` line is auto-inserted. |
| `modules` | `[]` | **List of `module load ...` lines**, one string per line, e.g. `["module load lammps"]`. An **empty list is not "no modules"** on a recognized cluster — it falls back to that cluster's built-in default (e.g. ISAAC's `module unload PE-intel` + `module load anaconda3` + `source $ANACONDA3_SH`), which you usually want anyway for `conda activate` to work. |
| `conda_env` | `"/lustre/isaac24/scratch/qzhang55/myVenvs/mtp_env"` | Full path to a conda environment; emits `conda activate <path>` in the job's env section. `null` to skip. |
| `lammps_binary_path` | `"/lustre/isaac24/scratch/qzhang55/MLIPs/lammps/src/lmp_mpi"` | Full path to a custom-built `lmp_mpi` (e.g. an MTP-enabled build). Leave `null` to use plain `lmp_mpi` from `PATH`/`modules`. |
| `lammps_src_path` | `"/lustre/isaac24/scratch/qzhang55/MLIPs/lammps/src"` | Directory prepended to `LD_LIBRARY_PATH` before the run, for binaries built without RPATH. Often — but not always — the same directory as `lammps_binary_path` (identical here because this build drops `lmp_mpi` directly in `src/` alongside its shared libs; for other builds these two paths can differ). |

---

### Two concrete `hpc` examples

**A — module-based LAMMPS (stock binary on PATH):**
```json
"hpc": {
  "mode": "scatter",
  "cluster": "ISAAC",
  "partition_key": "c",
  "walltime": "24:00:00",
  "nodes": 1,
  "ntasks": 40,
  "run_command": "srun -n {{NTASKS}} {{LAMMPS_BIN}} -in {{INPUT}}",
  "device": "cpu",
  "modules": ["module load lammps"],
  "conda_env": null,
  "lammps_binary_path": null,
  "lammps_src_path": null
}
```
Renders to: `srun -n 40 lmp_mpi -in in.reference.lammps` (plain `lmp_mpi`, since `lammps_binary_path` is unset).

**B — custom-compiled MTP build (own binary + shared libs):**
```json
"hpc": {
  "mode": "scatter",
  "cluster": "ISAAC",
  "partition_key": "s",
  "walltime": "03:00:00",
  "nodes": 1,
  "ntasks": 48,
  "run_command": "srun -n {{NTASKS}} {{LAMMPS_BIN}} -in {{INPUT}}",
  "device": "cpu",
  "modules": [],
  "conda_env": "/lustre/isaac24/scratch/qzhang55/myVenvs/mtp_env",
  "lammps_binary_path": "/lustre/isaac24/scratch/qzhang55/MLIPs/lammps/src/lmp_mpi",
  "lammps_src_path": "/lustre/isaac24/scratch/qzhang55/MLIPs/lammps/src"
}
```
Renders to:
```bash
export LD_LIBRARY_PATH=/lustre/isaac24/scratch/qzhang55/MLIPs/lammps/src:$LD_LIBRARY_PATH
srun -n 48 /lustre/isaac24/scratch/qzhang55/MLIPs/lammps/src/lmp_mpi -in in.reference.lammps
```

---

## Step 3 — Copy the potential file(s) into place

```bash
cp /path/to/library.meam /path/to/FeCr.meam potentials/my_potential/
```
Every name in `potential.potential_files` must exist here — one file per name,
copied verbatim into every case directory that potential builds.

---

## Step 4 — Check the potential files line up

```bash
python check.py
```
This is **user-initiated only** — nothing else in the pipeline calls it, so
run it by hand now, and again any time you touch `pot_lines`,
`potential_files`, or the files under `potentials/<pot_name>/`. It catches
the mistakes that `validate` (Step 5, which only looks at the JSON) can't see
because they involve the filesystem:

* a filename mentioned in `pot_lines` (`pair_style`/`pair_coeff`) that isn't
  listed in `potential_files`, or vice versa — the classic typo that means
  the file you copied is never actually referenced
* a name in `potential_files` with no matching file under
  `potentials/<pot_name>/`
* for MTP potentials, it prints `mlip.ini`'s contents so you can eyeball it,
  and flags it if the `.mtp` filename inside doesn't match what's in
  `potential_files`
* once case directories exist (after Step 6), it also flags any that are
  missing a potential file, or holding a stale copy that no longer matches
  the source under `potentials/<pot_name>/` (e.g. you edited the potential
  file after `main.py` already built the case dirs)

Equivalent CLI form: `gz-toolkit-potential-testing check --pot-inputs pot_inputs --run-dir .`

---

## Step 5 — Validate before touching the scheduler

```bash
gz-toolkit-potential-testing validate --pot-inputs pot_inputs
```
Catches: `type_map`/`masses` mismatches, missing `crystal_structures` entries,
elements referenced but not in `type_map`, bad `paral_degree`, malformed
`walltime`, `seakmc_defects` not present in `defect_catalog`, and more.
Fix everything it reports before proceeding — an empty report means you're
clear to submit.

---

## Step 6 — Build references and generate job scripts

```bash
python main.py
```
This box-relaxes a reference cell per element (writing `relaxed.data` +
`RESULT LC` / `RESULT ECOH` markers), and alongside it also builds and runs
a `reference/<element>/single_atom/` case — one isolated atom of that
element in a 100 Å box — printing `RESULT PE_SINGLE`. This is a sanity check
on the potential itself: an atom with no neighbors in range should read a
well-defined per-atom energy (often ~0, but potential-dependent — some
formalisms don't zero it). It then emits `submit_all.sh` plus every `.job`
file needed for your chosen `paral_degree`.

---

## Step 7 — Submit

```bash
bash submit_all.sh
```
Reference jobs run first; `build-defects` (which creates every point-defect,
elastic, alloy, gas-complex, and SEAKMC case directory) is chained to run
after reference relaxation completes, then every per-case job is submitted
(or grouped, under `mode: "batch"`).

---

## Step 8 — Check progress

```bash
gz-toolkit-potential-testing status --run-dir . -v
```
Classifies every case as `pending` / `running` / `done` / `failed` by
inspecting `log.lammps` (or `Seakmc_summary.csv`/`SPOut/` for SEAKMC cases).

If everything is `pending`/`failed` and never actually ran, re-run
`python check.py` — now that case directories exist it will also catch a
missing or stale potential file copy in them.

---

## Step 9 — `<pot_name>.json` builds itself once each potential is done

`gz_toolkit.pot_infobank` (used later by `buildmtx`/`defect`/`analyze` to
pull a potential's numbers by tag instead of you re-typing `lc` by hand)
needs a standardized per-pot record. That's `<pot_name>/<pot_name>.json`
(e.g. `gao2011/gao2011.json`) — named after the potential, not a fixed
filename, so it's already named correctly once copied into `pot_infobank/`.
It builds itself automatically: `emit_jobs` (Step 6) appends a final
`summarize-pot` job to every potential's job chain, dependent on that
potential's own last stage(s) — it fires the moment that potential is done,
independently of every other potential in `pot_inputs/`. Where it lands per
`paral_degree`:

| `paral_degree` | where the final step is |
|---|---|
| `1` (all pots, one job) | appended inline, end of `run_all.job` |
| `2` (one job per pot) | appended inline, end of `<pot_name>/run.job` |
| `3` (one job per work-group) | `<pot_name>/summarize_pot.job`, chained `--dependency=afterok` on the last group **and** SEAKMC (if enabled) |
| `4` (scatter) | `<pot_name>/summarize_pot.job`, submitted by the `build_defects_scatter` job once it's sbatched every per-case job, depending on all of them |

If that potential's JSON sets `"promote_to_infobank"` (see Step 2), the same
automatic step also copies `<pot_name>.json` into `gz_toolkit/pot_infobank/`.

`<pot_name>.json` is deliberately minimal — just what you need to go from
`pot="gao2011"` to an actual cell: `pot_name`, `type_map`, `masses`, and
`lc` (lattice constant per element, plus any alloy composition built, keyed
by its `<A><B>_<comp>` tag). For example:

```json
{
  "pot_name": "gao2011",
  "type_map": {"Fe": 1, "He": 2},
  "masses": {"1": 55.845, "2": 4.0026},
  "lc": {"Fe": 2.8553, "FeCr_Fe95Cr5": 2.859},
  "updated": "2026-07-10T00:09:03+00:00"
}
```

Everything else `_collect_pot_metrics` can produce — `Ecoh`, elastic
constants (`C11`…`C66` + Voigt `B`/`G`/`E`/`nu`), point-defect formation
energies, SEAKMC migration barriers, gas-defect-complex `Ef`/`Eb`,
dislocation loops — is deliberately left out of this file. None of it is
lost; it's all in `summary.csv` (Step 10) instead (except the off-diagonal
elastic couplings, which aren't persisted anywhere — only available by
calling `_collect_pot_metrics` directly in Python). This file answers one
question only — "what lattice constant do I build this potential at" —
everything about characterizing/comparing potentials belongs in the CSV.

**This step is independent of Step 10 below** — `summarize.py` parses
`log.lammps` files directly, on its own, and doesn't read or require
`<pot_name>.json` at all. So Step 10 works even if this step hasn't run yet
(a potential's own job chain not done, an old project generated before this
feature, `sbatch` missing when `paral_degree: 4` scattered so nothing
auto-chained) — you just won't get a `pot_infobank` entry for that potential
until `<pot_name>.json` exists. To build/refresh it by hand for one potential:

```bash
gz-toolkit-potential-testing summarize-pot --pot-config pot_inputs/my_potential.json --run-dir .
```

---

## Step 10 — Summarize

```bash
python summarize.py
```
Parses every potential's `log.lammps` files directly (independently of Step
9 / `<pot_name>.json`) into one wide CSV — one column per potential, rows for
lattice constants, cohesive energies, the isolated-atom sanity check, defect
formation energies, elastic constants (all 9 main Cij, plus Voigt-averaged
bulk modulus `B`, shear modulus `G`, Young's modulus `E`, and Poisson ratio
`nu` derived from the full 6x6 tensor), alloy lattice constants, gas binding
energies, and SEAKMC migration barriers. You can run this any time, even
partway through — potentials still running just show blank cells for
whatever hasn't finished yet.

Rows are grouped by **system**, not by metric type: every metric for `Fe`
together, then every metric for `Cr`, then any pure gas (`He`), then each
alloy system (`FeCr` from `two_element_suites`, `FeCrNi` from a
`multi_element_suites` entry, etc.), then each gas-defect-complex pair
(`FeHe`) — in the order those systems first appear across `pot_inputs/`. A `# <system>` marker
row precedes each group; empty groups (e.g. a gas that never got its own
metrics because everything attached to a `<metal><gas>` pair) are omitted.
Within a group, the usual `lc` → `Ecoh` → `Esingle` → `elastic` → `alloy_lc`
→ `Ef` → `Emig` → `Eb` ordering still applies:

```
metric,                   pot1,    pot2
# Fe
lc_Fe,                    2.8557,  2.8553
Ecoh_Fe,                  -4.122,  ...
Esingle_Fe,               0.0,     ...
Ef_1vac_Fe,               1.85,    ...
Emig_1vac_Fe,             0.62,    ...
# Cr
lc_Cr,                    2.8841,  ...
...
# FeCr
alloy_lc_FeCr_Fe95Cr5,    2.859,   ...
Ef_1vac_FeCr_Fe95_Cr5,    1.62,    ...   # include_alloy_defects: mean over alloy_defect_replicas
Ef_1vac_FeCr_Fe95_Cr5_std, 0.09,   ...   # replica-to-replica spread
# FeHe
Ef_He_in_Fe_tetra,        ...
Eb_V_He1_V1_Fe,           ...
```

### Splitting the output across several CSVs — `SUMMARY_GROUPS`

`summarize.py` (as scaffolded by `init`/`write-scripts`) has a plain module-level
dict near the top:

```python
SUMMARY_GROUPS: dict[str, list[str] | str] = {}
```

It's just a Python dict literal you edit by hand — each **key** becomes an
output filename (`<key>.csv`), each **value** picks which potentials go in
it: either a list of `pot_name` strings (the exact `pot_name` from that
potential's JSON, not the filename), or the string `"rest"`.

```python
SUMMARY_GROUPS = {}
```
Default. Ignore grouping entirely — one `summary.csv` with every potential in
`pot_inputs/`.

```python
SUMMARY_GROUPS = {
    "fe_pots": ["eam_fs", "mtp_v2", "chgnet"],
}
```
One file, `fe_pots.csv`, containing only those three potentials (any others
in `pot_inputs/` are silently left out — there's no default "everyone else"
unless you add a `"rest"` group too).

```python
SUMMARY_GROUPS = {
    "fe_pots": ["eam_fs", "mtp_v2"],
    "fecr_pots": ["fecr_1", "fecr_2", "fecr_3"],
    "everything_else": "rest",
}
```
Three files: `fe_pots.csv`, `fecr_pots.csv`, and `everything_else.csv` — the
last one automatically gets every potential in `pot_inputs/` that wasn't
named in `fe_pots` or `fecr_pots` above it (order of the other keys doesn't
matter, only which pots got explicitly claimed).

If your `summarize.py` predates this feature (still just calls
`summarize_all()`), upgrade it in place — no copy-pasting needed:

```bash
gz-toolkit-potential-testing write-scripts --root . --force --only summarize.py
```

`--force` overwrites just that one file with the current template (dropping
any hand edits to it); `--only` keeps `main.py`/`check.py` untouched. Drop
`--only` to refresh all three at once.
