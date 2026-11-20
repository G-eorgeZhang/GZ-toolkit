# SPPARKS — sub-skill of materials-simulation-skill

SPPARKS (Sandia) is a parallel **lattice / on-site kinetic Monte Carlo** code: grain growth
(Potts), diffusion, Ising, sintering, additive-manufacturing microstructures, and user apps.
It is *object* KMC on a fixed set of sites — unlike SEAKMC, there is no saddle-point search:
the event list and rates are defined by the chosen application.

**Docs:** https://spparks.github.io/doc/Manual.html — every command has a page
`https://spparks.github.io/doc/<command>.html` (e.g. `app_style.html`, `sweep.html`,
`solve_style.html`, `diag_style.html`). App-specific commands are documented on
`app_<style>.html`. Fetch the page before writing any command you are not sure of; SPPARKS
syntax is app-specific and easy to get subtly wrong.

## 1. How an input script is read

- Read top-to-bottom, one command at a time; each takes effect when read. Order matters
  (a species count cannot be set before the species exists; settings before `run`).
- `#` comments, `&` line continuation, `$x` / `${name}` variable substitution, double quotes
  group words. Commands are case-sensitive and lowercase.
- Typical structure: **(1)** `app_style`, `solve_style` or `sweep` (and lattice/sites);
  **(2)** settings (`seed`, `temperature`, app-specific commands, `stats`, `diag_style`, `dump`);
  **(3)** `run`. Settings + runs can repeat.

## 2. Two ways to evolve the system

| Method | Command | Use |
|--------|---------|-----|
| Rejection KMC / Metropolis sweeps | `sweep random` / `sweep raster` (optionally `mask yes`) | Lattice apps with many sites; time in MC sweeps |
| Rejection-free (BKL / n-fold way) KMC | `solve_style tree` or `linear` (or `group`) | True KMC time from rates; low-rate regimes |

Not every app supports both — check `app_<style>.html`.

## 3. Minimal example (2-d Potts grain growth, the canonical SPPARKS example)

```
seed            {{SEED}}
app_style       potts 100                  # 100 spin values
dimension       2
lattice         sq/8n 1.0
region          box block 0 100 0 100 -0.5 0.5
create_box      box
create_sites    box
set             site range 1 100
sweep           random
diag_style      energy
temperature     {{T}}                      # in units of the app's energy scale
stats           10.0
dump            1 text 10.0 dump.potts
run             100.0
```

`read_sites` instead of `create_sites` when the microstructure comes from step 2 (a site file
written with GZ_toolkit or a converter) — check the expected site-file format on
`read_sites.html` before writing it.

## 4. Running

```bash
spk_mpi < in.potts                                  # stdin by default
mpirun -np 16 spk_mpi -in in.potts -log log.spparks -screen none
spk_mpi -var T 0.5 -var seed 12345 -in in.potts    # -var defines `variable name index value`
```
The executable is named after the build target (`spk_serial`, `spk_mpi`, `spk_<machine>`) —
take it from the user's reference job. Default log: `log.spparks`.

## 5. Workflow integration (steps 3–7)

- **Seeds:** SPPARKS requires `seed`. Use `{{SEED}}` so every run dir is statistically
  independent and reproducible.
- **Template candidates:** `temperature`, app parameters, lattice size, run length. Or keep the
  input fixed and pass `-var` in the job body.
- **Job body:** `srun spk_mpi -in in.spparks -log log.spparks -screen spparks.out`.
- **Checks:** ask the user; good defaults are that `log.spparks` contains `Loop time` (printed
  at the end of a run — confirm on a finished test run) and has no `ERROR`; or that the final
  dump time exists (`count_files` on dumps).
- **Summary:** the `stats` table in `log.spparks` (time, events, diagnostics) — use `regex` or a
  custom Python field; final grain-size statistics from the last dump via custom Python.

## 6. Checklist

- [ ] `seed` set, different per run when statistics are wanted
- [ ] app supports the chosen `sweep` / `solve_style`
- [ ] `stats` and `dump` intervals are in the app's time units and divide the run length sensibly
- [ ] lattice dimensions and `region` match `dimension`
- [ ] number of processors compatible with the lattice decomposition (see `processors.html`)
