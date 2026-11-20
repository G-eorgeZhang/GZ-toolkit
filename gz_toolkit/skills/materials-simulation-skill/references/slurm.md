# SLURM reference (workflow step 5)

Expert assistant for writing, submitting, monitoring, and automating HPC simulation jobs on SLURM clusters.

## Scope

- Writing and editing SLURM `.job` (batch) scripts
- Writing bash automation scripts for job submission workflows
- SLURM commands: `sbatch`, `squeue`, `scancel`, `sinfo`, `sacct`, `scontrol`, `srun`
- Job arrays, job dependencies, resubmission loops
- Simulation-specific job patterns: LAMMPS, GROMACS, NAMD (MD); VASP, Quantum ESPRESSO (DFT)
- Advanced Linux usage on HPC (file ops, loops, piping, awk/sed, find, parallel)
- Troubleshooting failed/stuck/OOM jobs

> **SLURM documentation reference:** https://slurm.schedmd.com/documentation.html  
> For any SLURM directive, flag, or behavior that needs verification, fetch the relevant page from docs above. Key pages:
> - Quickstart: https://slurm.schedmd.com/quickstart.html
> - sbatch man page: https://slurm.schedmd.com/sbatch.html
> - Job arrays: https://slurm.schedmd.com/job_array.html
> - Job dependencies: https://slurm.schedmd.com/sbatch.html#OPT_dependency
> - Job state/reason codes: https://slurm.schedmd.com/job_state_codes.html / https://slurm.schedmd.com/job_reason_codes.html
> - GRES (GPU scheduling): https://slurm.schedmd.com/gres.html

---

## Writing .job Files

Always produce well-commented SLURM batch scripts. Use this canonical structure:

```bash
#!/bin/bash
#SBATCH --job-name=<name>
#SBATCH --output=%x_%j.out        # %x=jobname, %j=jobid
#SBATCH --error=%x_%j.err
#SBATCH --partition=<partition>
#SBATCH --nodes=<N>
#SBATCH --ntasks-per-node=<n>
#SBATCH --cpus-per-task=<c>       # for threaded/OpenMP
#SBATCH --mem=<MB>                # or --mem-per-cpu
#SBATCH --time=<HH:MM:SS>
#SBATCH --account=<account>       # if required

# --- environment ---
module purge
module load <software/version>

# --- run ---
cd $SLURM_SUBMIT_DIR
srun <executable> <args>
```

**Rules when writing job files:**
1. Always ask or confirm: partition name, node/core count, walltime, memory, account — if not provided.
2. Use `srun` inside batch scripts for MPI jobs; use direct executable call for serial/threaded-only.
3. Use `$SLURM_SUBMIT_DIR`, `$SLURM_JOB_ID`, `$SLURM_NTASKS` etc. for portability.
4. Always add `--output` and `--error` directives — never let logs go to default locations silently.
5. For GPU jobs, add `#SBATCH --gres=gpu:<type>:<count>`.
6. When walltime is unknown, suggest a conservative overestimate and explain how to check `sacct` afterward.

---

## Simulation-Specific Patterns

### LAMMPS
```bash
srun lmp -in in.lammps -log log.lammps
# For large jobs, add: -sf gpu -pk gpu <N>  (if GPU-enabled build)
```
- Output files: `log.lammps`, dump files, restart files
- Tip: always redirect screen output too: `lmp -in in.lammps > lammps.out 2>&1`

### GROMACS
```bash
# Typical MD pipeline (can split into separate jobs with dependencies)
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
srun gmx_mpi mdrun -v -deffnm md -ntmpi $SLURM_NTASKS
```

### NAMD
```bash
srun namd3 +p $SLURM_NTASKS input.namd > namd.out 2>&1
```

### VASP
```bash
srun vasp_std > vasp.out 2>&1
# Ensure INCAR, POSCAR, POTCAR, KPOINTS are in $SLURM_SUBMIT_DIR
# Check NCORE/NPAR settings match node topology
```
- See also: `software/vasp/GUIDE.md` for VASP input file guidance

### Quantum ESPRESSO
```bash
srun pw.x -in scf.in > scf.out 2>&1
# For k-point parallelism: pw.x -npool <N> -in scf.in
```

---

## Job Arrays

Use job arrays for parameter sweeps or replicate runs:

```bash
#SBATCH --array=1-100             # or 1-100:2, 0-9, 1,3,7
#SBATCH --array=1-100%10          # limit to 10 concurrent

# Inside script: use $SLURM_ARRAY_TASK_ID
cd run_${SLURM_ARRAY_TASK_ID}
srun ./my_sim input_${SLURM_ARRAY_TASK_ID}.dat
```

For parameter mapping from a file:
```bash
PARAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" params.txt)
```

Full docs: https://slurm.schedmd.com/job_array.html

---

## Job Dependencies

Chain jobs so they run in sequence or conditionally:

```bash
JOB1=$(sbatch job1.job | awk '{print $NF}')
JOB2=$(sbatch --dependency=afterok:$JOB1 job2.job | awk '{print $NF}')
sbatch --dependency=afterok:$JOB2 job3.job
```

Key dependency types:
| Type | Meaning |
|------|---------|
| `afterok:<id>` | Start after job succeeds (exit 0) |
| `afternotok:<id>` | Start after job fails |
| `afterany:<id>` | Start after job ends (any status) |
| `after:<id>` | Start after job begins |
| `singleton` | Wait for all jobs with same name+user to finish |

---

## Bash Automation Scripts

When writing submission automation scripts:

```bash
#!/bin/bash
# Submit a sweep of jobs, one per directory

BASE_DIR=$(pwd)
for dir in run_*/; do
    cd "$BASE_DIR/$dir"
    JOB_ID=$(sbatch submit.job | awk '{print $NF}')
    echo "Submitted $dir → Job $JOB_ID"
done
```

**Patterns to know:**
- Capture job ID: `JOB_ID=$(sbatch job.sh | awk '{print $NF}')`
- Wait for completion: use `--dependency` chaining rather than `sleep` polling
- Check status: `squeue -j $JOB_ID -h -o "%T"` → returns state string
- Cancel all user jobs: `scancel -u $USER`
- Resubmit on failure: check `sacct -j $JOB_ID --format=State` then resubmit

---

## Monitoring & Troubleshooting

### Common commands
```bash
squeue -u $USER                        # your jobs
squeue -u $USER --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R"
sinfo -s                               # partition summary
sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS
scontrol show job <jobid>              # full job details
```

### Common failure reasons
| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `OOM` / `OUT_OF_MEMORY` | `--mem` too low | Increase `--mem` or `--mem-per-cpu` |
| `TIMEOUT` | Walltime exceeded | Increase `--time`; use restart files |
| `NODE_FAIL` | Hardware issue | Add `#SBATCH --requeue`; check `sinfo` |
| Job stuck in PD | Wrong partition, resource limit, dependency | Check `squeue -j ID` reason column |
| Exit code ≠ 0 | Application error | Check `.err` file; run interactively with `srun --pty bash` |

For job reason codes: https://slurm.schedmd.com/job_reason_codes.html  
For exit codes: https://slurm.schedmd.com/job_exit_code.html

---

## Useful Linux One-Liners for HPC

```bash
# Find all output files modified in last 24h
find . -name "*.out" -mtime -1

# Watch job queue live
watch -n 30 squeue -u $USER

# Count lines across all output files
wc -l *.out | sort -n

# Check disk usage by directory
du -sh */ | sort -h

# Extract energy from LAMMPS log
grep "^Step" log.lammps | head -1   # header
awk '/^[0-9]/{print $1, $2}' log.lammps > energy.dat

# Monitor a running job's output
tail -f ${SLURM_JOB_ID}.out
```

---

## Interaction Guidelines

- **If partition/account/walltime are unknown**: always ask before writing the job file — these are cluster-specific and wrong values cause immediate failures.
- **If connected to HPC via SSH**: offer to inspect the cluster with `sinfo`, `squeue`, or `sacct` to inform the job script.
- **If the user pastes an existing job file**: review it for correctness, suggest improvements, and flag any missing or risky directives.
- **Always verify non-trivial SLURM behavior** against the official docs at https://slurm.schedmd.com/documentation.html before giving advice.
