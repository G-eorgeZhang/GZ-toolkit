# Other simulation codes (GROMACS, Quantum ESPRESSO, CP2K, ABINIT, NAMD, LAMMPS-like codes ...)

The workflow (steps 0–8) is code-agnostic. For a code without its own guide here:

1. **Find the official documentation** and read the input-format and command-line pages before
   writing anything. Cite the pages you used. For physics settings (cutoffs, k-points,
   functionals, thermostats), check the documentation and, where it matters, recent
   literature for the same material/property — tell the user which values came from where.
2. **Write one reference input** in `reffiles/`, have the user check it (GATE), then template
   it with `{{VAR}}` placeholders or pass values on the command line in the job body.
3. **Identify the per-run file set** (structure, input, parameter/pseudopotential files,
   job) — missing user files are requested from the user, never downloaded unasked
   (pseudopotentials and force fields included).
4. **Run one test** before the campaign: confirm the executable/modules from the user's
   reference job, and look at what a *finished* output looks like — that defines the checks.
5. **Record** in `sim_plan.yaml` (`software:`, `checks:`, `summary:`), and if the code will be
   used again, offer to add a `software/<code>/GUIDE.md` sub-skill following the LAMMPS one's
   "Workflow integration" section.

Quick pointers (verify on the docs):

| Code | Run | Finished marker (verify) | Energies |
|------|-----|--------------------------|----------|
| Quantum ESPRESSO | `srun pw.x -in scf.in > scf.out` | `JOB DONE.` in output | `!    total energy` lines |
| GROMACS | `gmx grompp ...; srun gmx_mpi mdrun -deffnm md` | `Finished mdrun` in md.log | `gmx energy` → .xvg |
| CP2K | `srun cp2k.psmp -i in.inp -o out.log` | `PROGRAM ENDED AT` | `ENERGY| Total FORCE_EVAL` |
| NAMD | `srun namd3 +p$SLURM_NTASKS in.namd > namd.out` | `End of program` | `ENERGY:` lines |
