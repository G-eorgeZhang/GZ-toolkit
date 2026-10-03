"""Job-script generation for the four parallelization tiers.

The flow looks like::

    python main.py
        ├── prepare_reference (writes reference dirs + LAMMPS inputs)
        ├── emit job scripts via this module
        └── prints `bash submit_all.sh` for the user to run

paral_degree mapping:
    1 — single .job at the project root, runs every potential's full pipeline.
    2 — one .job per potential at <root>/<pot>/run.job; submit_all.sh sbatches each.
    3 — one .job per (potential, work-group) at <root>/<pot>/<group>.job;
        submit_all.sh chains them with --dependency=afterok so reference runs first.
    4 — one .job per simulation. Per-case .job files are emitted ahead of time;
        a build-defects job (depending on reference) creates the defect dirs and
        sbatches every per-case .job. SEAKMC jobs depend on their case relax job.

Notes:
* The "build_defects" step is a Python invocation of
  ``gz_toolkit.potential_testing.cli build-defects --pot-config <json>`` that
  reads relaxed.data from the reference dir and creates every defect / alloy /
  gas case directory.
* GPU support is signalled by ``hpc.device == 'gpu'`` and adds
  ``--gres=gpu:N`` lines to the SLURM header.
"""

from __future__ import annotations

from pathlib import Path
import textwrap

from gz_toolkit.jobs.template_engine import JobTemplate
from gz_toolkit.potential_testing.config import PotentialConfig


# ---------------------------------------------------------------------------
# SLURM header rendering (uses gz_toolkit.jobs.JobTemplate but adds GPU lines)
# ---------------------------------------------------------------------------


def _render_header(
    config: PotentialConfig,
    jobname: str,
    extra_modules: list[str] | None = None,
) -> str:
    hpc = config.hpc
    tpl = JobTemplate(cluster=hpc.cluster) if hpc.cluster else JobTemplate()

    modules = list(hpc.modules)
    if extra_modules:
        modules.extend(extra_modules)
    module_block = "\n".join(modules) if modules else None

    env_lines: list[str] = []
    if hpc.conda_env:
        env_lines.append(f"conda activate {hpc.conda_env}")

    text = tpl.render(
        jobname=jobname,
        partition_key=hpc.partition_key,
        walltime=hpc.walltime,
        nodes=hpc.nodes,
        ntasks=hpc.ntasks,
        module_commands=module_block,
        env_commands=env_lines if env_lines else None,
        run_command="",  # body added separately
    )

    if hpc.device == "gpu" and hpc.gpus_per_node > 0:
        # Insert --gres=gpu:N right after the partition line.
        gres_line = f"#SBATCH --gres=gpu:{hpc.gpus_per_node}"
        if "#SBATCH --partition" in text:
            text = text.replace(
                "#SBATCH --partition",
                f"{gres_line}\n#SBATCH --partition",
                1,
            )
        else:
            # Fall back: prepend after shebang.
            lines = text.splitlines()
            lines.insert(1, gres_line)
            text = "\n".join(lines)
    return text


# ---------------------------------------------------------------------------
# Run-command formatting
# ---------------------------------------------------------------------------


def _run_lammps(config: PotentialConfig, input_name: str) -> str:
    hpc = config.hpc
    cmd = hpc.run_command
    cmd = cmd.replace("{{NTASKS}}", str(hpc.ntasks))
    cmd = cmd.replace("{{INPUT}}", input_name)
    cmd = cmd.replace("{{LAMMPS_BIN}}", hpc.lammps_binary_path or "lmp_mpi")
    if hpc.lammps_src_path:
        cmd = f"export LD_LIBRARY_PATH={hpc.lammps_src_path}:$LD_LIBRARY_PATH\n{cmd}"
    return cmd


def _run_seakmc(config: PotentialConfig) -> str:
    # SEAKMC parallel driver — uses mpirun rather than srun for portability.
    return f"mpirun -np {config.hpc.ntasks} python run_seakmc_p.py"


# ---------------------------------------------------------------------------
# Body builders for each work group
# ---------------------------------------------------------------------------


def _bash_for_dir(case_dir: Path, run_cmd: str, project_root: Path) -> str:
    """Emit ``cd <abs> && <run>``."""
    abs_path = case_dir.resolve()
    return f'cd "{abs_path}"\n{run_cmd}\n'


def _reference_body(config: PotentialConfig, project_root: Path) -> str:
    pot_dir = project_root / config.potential.pot_name
    elements = config.potential.single_elements or [next(iter(config.potential.type_map))]
    lines = ["# --- Stage 1: reference (box/relax per element) ---"]
    for el in elements:
        ref = pot_dir / "reference" / el
        lines.append(_bash_for_dir(ref, _run_lammps(config, "in.reference.lammps"), project_root))
        lines.append(_bash_for_dir(ref / "single_atom", _run_lammps(config, "in.single_atom.lammps"), project_root))
    return "\n".join(lines)


def _build_defects_body(config: PotentialConfig, pot_config_path: Path, project_root: Path) -> str:
    return (
        "# --- Stage 2: build defect / alloy / gas case dirs from relaxed reference ---\n"
        f'cd "{project_root.resolve()}"\n'
        f'python -m gz_toolkit.potential_testing.cli build-defects '
        f'--pot-config "{pot_config_path.resolve()}" --run-dir "{project_root.resolve()}"\n'
    )


def _summarize_pot_body(config: PotentialConfig, pot_config_path: Path, project_root: Path) -> str:
    """Final stage: parse this pot's own logs into <pot_name>.json (and
    promote it, if configured) — runs automatically once every other stage
    for this pot is done, so nobody has to remember to call
    ``summarize-pot`` by hand."""
    return (
        "# --- Final: build this pot's <pot_name>.json (and promote it, if configured) ---\n"
        f'cd "{project_root.resolve()}"\n'
        f'python -m gz_toolkit.potential_testing.cli summarize-pot '
        f'--pot-config "{pot_config_path.resolve()}" --run-dir "{project_root.resolve()}"\n'
    )


def _list_case_dirs(pot_dir: Path, group: str, input_name: str) -> list[tuple[Path, str]]:
    """Return (case_dir, lammps_input_filename) pairs for a work-group."""
    if group == "point_defects":
        if not (pot_dir / "point_defects").is_dir():
            return []
        cases: list[tuple[Path, str]] = []
        for el_dir in sorted((pot_dir / "point_defects").iterdir()):
            if not el_dir.is_dir():
                continue
            for case_dir in sorted(el_dir.iterdir()):
                if case_dir.is_dir():
                    cases.append((case_dir, input_name))
        return cases
    if group == "elastic":
        root = pot_dir / "elastic"
        if not root.is_dir():
            return []
        # one level deep: elastic/<element> and elastic/<alloy_tag>
        return [(c, input_name) for c in sorted(root.iterdir()) if c.is_dir()]
    if group == "loops":
        root = pot_dir / "loops"
        if not root.is_dir():
            return []
        return [(c, input_name) for el in sorted(root.iterdir()) if el.is_dir()
                for c in sorted(el.iterdir()) if c.is_dir()]
    if group == "disloc_lines":
        return []  # not implemented; never executed
    if group == "alloy_lc":
        root = pot_dir / "alloy_lc"
        if not root.is_dir():
            return []
        return [(c, "in.alloy.lammps") for pair in sorted(root.iterdir()) if pair.is_dir()
                for c in sorted(pair.iterdir()) if c.is_dir()]
    if group == "gas_complexes":
        root = pot_dir / "gas_complexes"
        if not root.is_dir():
            return []
        return [(c, input_name) for pair in sorted(root.iterdir()) if pair.is_dir()
                for c in sorted(pair.iterdir()) if c.is_dir()]
    if group == "seakmc":
        cases: list[tuple[Path, str]] = []
        for el_dir in sorted((pot_dir / "point_defects").glob("*")) if (pot_dir / "point_defects").is_dir() else []:
            for case_dir in sorted(el_dir.iterdir()) if el_dir.is_dir() else []:
                meb = case_dir / "meb"
                if meb.is_dir():
                    cases.append((meb, ""))
        return cases
    return []


def _group_body(group: str, cases: list[tuple[Path, str]], config: PotentialConfig) -> str:
    if not cases:
        return f"# --- {group}: no cases ---\n"
    lines = [f"# --- {group}: {len(cases)} case(s) ---"]
    for case_dir, input_name in cases:
        if group == "seakmc":
            cmd = _run_seakmc(config)
        else:
            cmd = _run_lammps(config, input_name)
        lines.append(f'cd "{case_dir.resolve()}"\n{cmd}\n')
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Public entry: emit jobs for one project (all configs)
# ---------------------------------------------------------------------------


def emit_jobs(
    configs: list[PotentialConfig],
    project_root: Path,
    pot_inputs_dir: Path,
) -> Path:
    """Generate every .job file + a top-level submit_all.sh.

    Returns the path of submit_all.sh.
    """
    project_root = project_root.resolve()
    submit_lines: list[str] = ["#!/bin/bash", "set -e", ""]

    # Determine paral_degree from the first config (must be uniform across pots
    # because the master submit script is shared).
    paral = configs[0].workflow.paral_degree if configs else 4
    submit_lines.append(f"# paral_degree = {paral}")
    submit_lines.append("")

    if paral == 1:
        path = _emit_paral1(configs, project_root, pot_inputs_dir)
        submit_lines.append(f'sbatch "{path}"')
    elif paral == 2:
        for cfg in configs:
            jpath = _emit_paral2(cfg, project_root, pot_inputs_dir)
            submit_lines.append(f'sbatch "{jpath}"')
    elif paral == 3:
        for cfg in configs:
            chain = _emit_paral3(cfg, project_root, pot_inputs_dir)
            submit_lines.extend(chain)
            submit_lines.append("")
    elif paral == 4:
        for cfg in configs:
            chain = _emit_paral4(cfg, project_root, pot_inputs_dir)
            submit_lines.extend(chain)
            submit_lines.append("")
    else:
        raise ValueError(f"Unsupported paral_degree: {paral}")

    submit = project_root / "submit_all.sh"
    submit.write_text("\n".join(submit_lines) + "\n", encoding="utf-8")
    submit.chmod(0o755)
    return submit


# --- paral=1: one .job for everything --------------------------------------


def _all_stages_body(cfg: PotentialConfig, pot_config_path: Path, project_root: Path) -> str:
    """Emit a sequential body that runs every stage in one shell.

    Case directories don't exist at python-emission time (they're built by the
    inline ``build_defects`` step), so all post-reference loops use bash globs
    rather than baked-in dir lists.
    """
    pot_dir = project_root / cfg.potential.pot_name
    parts: list[str] = [
        _reference_body(cfg, project_root),
        _build_defects_body(cfg, pot_config_path, project_root),
        _runtime_group_body(cfg, pot_dir, "point_defects", "in.relax.lammps"),
    ]
    if cfg.workflow.include_elastic:
        parts.append(_runtime_group_body(cfg, pot_dir, "elastic", "in.elastic"))
    if cfg.workflow.include_dislocation_loops:
        parts.append(_runtime_group_body(cfg, pot_dir, "loops", "in.relax.lammps"))
    if cfg.workflow.include_alloy_suite:
        parts.append(_runtime_group_body(cfg, pot_dir, "alloy_lc", "in.alloy.lammps"))
    if cfg.workflow.include_gas_complexes:
        parts.append(_runtime_group_body(cfg, pot_dir, "gas_complexes", "in.relax.lammps"))
    if cfg.workflow.seakmc_enabled:
        parts.append(_runtime_group_body(cfg, pot_dir, "seakmc", ""))
    parts.append(_summarize_pot_body(cfg, pot_config_path, project_root))
    return "\n".join(parts)


def _emit_paral1(configs: list[PotentialConfig], project_root: Path, pot_inputs_dir: Path) -> Path:
    """Single master job. NOTE: uses the FIRST config's HPC settings."""
    if not configs:
        raise RuntimeError("No potential configs provided")
    cfg0 = configs[0]
    body_parts: list[str] = []
    for cfg in configs:
        pot_json = pot_inputs_dir / f"{cfg.potential.pot_name}.json"
        body_parts.append(f"# ====== potential: {cfg.potential.pot_name} ======")
        body_parts.append(_all_stages_body(cfg, pot_json, project_root))
    header = _render_header(cfg0, jobname="pot_test_all")
    text = header + "\n" + "\n".join(body_parts) + "\ndate\n"
    out = project_root / "run_all.job"
    out.write_text(text, encoding="utf-8")
    return out


# --- paral=2: one .job per potential ---------------------------------------


def _emit_paral2(cfg: PotentialConfig, project_root: Path, pot_inputs_dir: Path) -> Path:
    pot_dir = project_root / cfg.potential.pot_name
    pot_dir.mkdir(parents=True, exist_ok=True)
    pot_json = pot_inputs_dir / f"{cfg.potential.pot_name}.json"
    body = _all_stages_body(cfg, pot_json, project_root)
    header = _render_header(cfg, jobname=f"pot_{cfg.potential.pot_name}")
    text = header + "\n" + body + "\ndate\n"
    out = pot_dir / "run.job"
    out.write_text(text, encoding="utf-8")
    return out


# --- paral=3: one .job per (potential, work-group) -------------------------


def _emit_paral3(cfg: PotentialConfig, project_root: Path, pot_inputs_dir: Path) -> list[str]:
    """Emit reference + build_defects + per-group jobs with sbatch dependencies."""
    pot_dir = project_root / cfg.potential.pot_name
    pot_dir.mkdir(parents=True, exist_ok=True)
    pot_json = pot_inputs_dir / f"{cfg.potential.pot_name}.json"

    submit: list[str] = [f"# --- {cfg.potential.pot_name}: paral=3 chain ---"]

    # Reference job
    ref_body = _reference_body(cfg, project_root)
    ref_text = _render_header(cfg, jobname=f"{cfg.potential.pot_name}_ref") + "\n" + ref_body + "\ndate\n"
    ref_path = pot_dir / "reference.job"
    ref_path.write_text(ref_text, encoding="utf-8")
    submit.append(f'ref_id=$(sbatch --parsable "{ref_path}")')
    submit.append('echo "  reference: $ref_id"')

    # Build-defects job
    bd_text = _render_header(cfg, jobname=f"{cfg.potential.pot_name}_bd") + "\n" + \
              _build_defects_body(cfg, pot_json, project_root) + "\ndate\n"
    bd_path = pot_dir / "build_defects.job"
    bd_path.write_text(bd_text, encoding="utf-8")
    submit.append(f'bd_id=$(sbatch --parsable --dependency=afterok:$ref_id "{bd_path}")')
    submit.append('echo "  build_defects: $bd_id"')

    # Per-work-group jobs (deferred reading of case lists at runtime — but the
    # cases don't exist yet at python time. We still emit a static body that
    # walks the directory at job-execution time.)
    groups: list[tuple[str, str]] = [("point_defects", "in.relax.lammps")]
    if cfg.workflow.include_elastic:
        groups.append(("elastic", "in.elastic"))
    if cfg.workflow.include_dislocation_loops:
        groups.append(("loops", "in.relax.lammps"))
    if cfg.workflow.include_alloy_suite:
        groups.append(("alloy_lc", "in.alloy.lammps"))
    if cfg.workflow.include_gas_complexes:
        groups.append(("gas_complexes", "in.relax.lammps"))

    last_id_var = "bd_id"
    for group, input_name in groups:
        # Body uses bash globs to discover case dirs at runtime; safer than
        # relying on python having seen them at submit time.
        body = _runtime_group_body(cfg, pot_dir, group, input_name)
        text = _render_header(cfg, jobname=f"{cfg.potential.pot_name}_{group}") + "\n" + body + "\ndate\n"
        path = pot_dir / f"{group}.job"
        path.write_text(text, encoding="utf-8")
        var = f"{group}_id"
        submit.append(f'{var}=$(sbatch --parsable --dependency=afterok:${last_id_var} "{path}")')
        submit.append(f'echo "  {group}: ${var}"')
        last_id_var = var

    dep_ids = [last_id_var]
    if cfg.workflow.seakmc_enabled:
        body = _runtime_group_body(cfg, pot_dir, "seakmc", "")
        text = _render_header(cfg, jobname=f"{cfg.potential.pot_name}_seakmc") + "\n" + body + "\ndate\n"
        path = pot_dir / "seakmc.job"
        path.write_text(text, encoding="utf-8")
        submit.append(f'seakmc_id=$(sbatch --parsable --dependency=afterok:$point_defects_id "{path}")')
        submit.append('echo "  seakmc: $seakmc_id"')
        # seakmc branches off point_defects, in parallel with any later groups —
        # so the final summarize_pot job must wait on both branches, not just
        # whichever "last_id_var" happens to be.
        dep_ids.append("seakmc_id")

    # Final job: build this pot's <pot_name>.json (and promote it, if configured)
    # once every branch above has finished — no manual summarize-pot needed.
    dep_clause = ":".join(f"${v}" for v in dep_ids)
    sm_text = _render_header(cfg, jobname=f"{cfg.potential.pot_name}_summarize") + "\n" + \
              _summarize_pot_body(cfg, pot_json, project_root) + "\ndate\n"
    sm_path = pot_dir / "summarize_pot.job"
    sm_path.write_text(sm_text, encoding="utf-8")
    submit.append(f'summarize_id=$(sbatch --parsable --dependency=afterok:{dep_clause} "{sm_path}")')
    submit.append('echo "  summarize_pot: $summarize_id"')

    return submit


def _runtime_group_body(cfg: PotentialConfig, pot_dir: Path, group: str, input_name: str) -> str:
    """Emit a body that *discovers* case dirs at job-execution time.

    Used for paral=3 where the case dirs only exist after ``build_defects`` runs.
    """
    abs_pot = pot_dir.resolve()
    if group == "seakmc":
        cmd = _run_seakmc(cfg)
        return textwrap.dedent(f"""\
            for d in "{abs_pot}"/point_defects/*/*/meb; do
              [ -d "$d" ] || continue
              cd "$d"
              {cmd}
            done
            """)
    cmd = _run_lammps(cfg, input_name)
    if group == "point_defects":
        glob = f'"{abs_pot}"/point_defects/*/*'
    elif group == "elastic":
        glob = f'"{abs_pot}"/elastic/*'
    elif group == "loops":
        glob = f'"{abs_pot}"/loops/*/*'
    elif group == "alloy_lc":
        glob = f'"{abs_pot}"/alloy_lc/*/*'
    elif group == "gas_complexes":
        glob = f'"{abs_pot}"/gas_complexes/*/*'
    else:
        return f"# unknown group: {group}\n"
    return textwrap.dedent(f"""\
        for d in {glob}; do
          [ -d "$d" ] || continue
          # Skip if no LAMMPS input present (e.g., disloc_lines placeholder).
          [ -f "$d/{input_name}" ] || continue
          cd "$d"
          {cmd}
        done
        """)


# --- paral=4: one .job per simulation --------------------------------------


def _emit_paral4(cfg: PotentialConfig, project_root: Path, pot_inputs_dir: Path) -> list[str]:
    """Reference + build_defects-and-scatter chain.

    The build_defects job ALSO sbatches every per-case .job because the case
    directories don't exist until that point. Per-case .job *templates* are
    emitted by the build_defects step at runtime.
    """
    pot_dir = project_root / cfg.potential.pot_name
    pot_dir.mkdir(parents=True, exist_ok=True)
    pot_json = pot_inputs_dir / f"{cfg.potential.pot_name}.json"
    submit: list[str] = [f"# --- {cfg.potential.pot_name}: paral=4 chain ---"]

    # Reference scatter — one .job per element reference dir.
    elements = cfg.potential.single_elements or [next(iter(cfg.potential.type_map))]
    ref_ids_var = []
    for i, el in enumerate(elements):
        ref_path = pot_dir / "reference" / el / "run.job"
        ref_path.parent.mkdir(parents=True, exist_ok=True)
        body = _bash_for_dir(ref_path.parent, _run_lammps(cfg, "in.reference.lammps"), project_root)
        body += _bash_for_dir(ref_path.parent / "single_atom",
                              _run_lammps(cfg, "in.single_atom.lammps"), project_root)
        text = _render_header(cfg, jobname=f"{cfg.potential.pot_name}_ref_{el}") + "\n" + body + "\ndate\n"
        ref_path.write_text(text, encoding="utf-8")
        var = f"ref{i}_id"
        ref_ids_var.append(var)
        submit.append(f'{var}=$(sbatch --parsable "{ref_path}")')
        submit.append(f'echo "  ref/{el}: ${var}"')
    ref_dep = ",".join(f"${v}" for v in ref_ids_var) if ref_ids_var else ""

    # Build-defects + scatter job (single job that runs python build-defects
    # and then sbatches every per-case job).
    bd_body = _build_defects_body(cfg, pot_json, project_root)
    bd_body += textwrap.dedent(f"""\

        # --- Stage 3: scatter sbatch every case ---
        python -m gz_toolkit.potential_testing.cli scatter-cases \\
            --pot-config "{pot_json.resolve()}" --run-dir "{project_root.resolve()}"
        """)
    text = _render_header(cfg, jobname=f"{cfg.potential.pot_name}_bd_scatter") + "\n" + bd_body + "\ndate\n"
    bd_path = pot_dir / "build_defects_scatter.job"
    bd_path.write_text(text, encoding="utf-8")
    dep_clause = f"--dependency=afterok:{ref_dep}" if ref_dep else ""
    submit.append(f'bd_id=$(sbatch --parsable {dep_clause} "{bd_path}")')
    submit.append('echo "  build_defects_scatter: $bd_id"')

    return submit


# ---------------------------------------------------------------------------
# Per-case scatter (called from inside the bd_scatter job at runtime)
# ---------------------------------------------------------------------------


def scatter_cases(cfg: PotentialConfig, project_root: Path, pot_inputs_dir: Path) -> list[Path]:
    """Write a per-case .job for every directory built post-reference, then sbatch them.

    Honours ``hpc.max_jobs_per_batch`` by splitting any group exceeding that
    limit into ``<group>_1.job``, ``<group>_2.job``, etc. (each .job containing
    a chunk of dirs run sequentially).
    """
    pot_dir = project_root / cfg.potential.pot_name
    written: list[Path] = []
    submitted: list[str] = []

    groups: list[tuple[str, str]] = [("point_defects", "in.relax.lammps")]
    if cfg.workflow.include_elastic:
        groups.append(("elastic", "in.elastic"))
    if cfg.workflow.include_dislocation_loops:
        groups.append(("loops", "in.relax.lammps"))
    if cfg.workflow.include_alloy_suite:
        groups.append(("alloy_lc", "in.alloy.lammps"))
    if cfg.workflow.include_gas_complexes:
        groups.append(("gas_complexes", "in.relax.lammps"))

    case_to_jobid: dict[Path, str] = {}

    for group, input_name in groups:
        cases = _list_case_dirs(pot_dir, group, input_name)
        if not cases:
            continue
        if len(cases) <= cfg.hpc.max_jobs_per_batch:
            # True scatter — one .job per case.
            for case_dir, in_name in cases:
                jpath = case_dir / "run.job"
                cmd = _run_lammps(cfg, in_name)
                body = _bash_for_dir(case_dir, cmd, project_root)
                text = _render_header(cfg, jobname=case_dir.name) + "\n" + body + "\ndate\n"
                jpath.write_text(text, encoding="utf-8")
                written.append(jpath)
                jid = _sbatch(jpath)
                if jid:
                    case_to_jobid[case_dir] = jid
                    submitted.append(jid)
        else:
            # Chunk into groups of max_jobs_per_batch — emit batched .job files.
            chunk = cfg.hpc.max_jobs_per_batch
            for k in range(0, len(cases), chunk):
                slab = cases[k:k + chunk]
                lo, hi = k, k + len(slab) - 1
                lines = [_bash_for_dir(d, _run_lammps(cfg, n), project_root) for d, n in slab]
                body = "\n".join(lines)
                text = _render_header(cfg, jobname=f"{group}_{lo}-{hi}") + "\n" + body + "\ndate\n"
                jpath = pot_dir / f"{group}_{lo}-{hi}.job"
                jpath.write_text(text, encoding="utf-8")
                written.append(jpath)
                jid = _sbatch(jpath)
                if jid:
                    submitted.append(jid)
                    for d, _ in slab:
                        case_to_jobid[d] = jid

    # SEAKMC depends on its specific point-defect job.
    if cfg.workflow.seakmc_enabled:
        for case_dir, _ in _list_case_dirs(pot_dir, "seakmc", ""):
            parent_relax_dir = case_dir.parent  # .../point_defects/<el>/<case>/
            jpath = case_dir / "run.job"
            cmd = _run_seakmc(cfg)
            body = _bash_for_dir(case_dir, cmd, project_root)
            text = _render_header(cfg, jobname=f"seakmc_{case_dir.parent.name}") + "\n" + body + "\ndate\n"
            jpath.write_text(text, encoding="utf-8")
            written.append(jpath)
            dep = case_to_jobid.get(parent_relax_dir)
            jid = _sbatch(jpath, dependency=dep)
            if jid:
                submitted.append(jid)

    # Final job: once every scattered case (+ SEAKMC) job for this pot has
    # finished, build <pot_name>.json (and promote it, if configured)
    # automatically — no manual summarize-pot needed. Always written for a manual fallback;
    # only auto-submitted if we actually collected job ids to depend on
    # (i.e. sbatch was available above).
    pot_json = pot_inputs_dir / f"{cfg.potential.pot_name}.json"
    sm_path = pot_dir / "summarize_pot.job"
    sm_text = _render_header(cfg, jobname=f"{cfg.potential.pot_name}_summarize") + "\n" + \
              _summarize_pot_body(cfg, pot_json, project_root) + "\ndate\n"
    sm_path.write_text(sm_text, encoding="utf-8")
    written.append(sm_path)
    if submitted:
        _sbatch(sm_path, dependency=":".join(submitted))

    return written


def _sbatch(job_path: Path, dependency: str | None = None) -> str | None:
    """Run sbatch and return the job id, or None if sbatch is not available.

    Used at runtime from inside a SLURM job (the ``build_defects_scatter`` job),
    so we expect sbatch on $PATH. If not, we just leave the .job files in place.
    """
    import shutil
    import subprocess
    if not shutil.which("sbatch"):
        return None
    cmd = ["sbatch", "--parsable"]
    if dependency:
        cmd.append(f"--dependency=afterok:{dependency}")
    cmd.append(str(job_path))
    try:
        out = subprocess.check_output(cmd, text=True).strip()
        return out.split(";")[0]  # sbatch --parsable returns "<jobid>" or "<jobid>;<cluster>"
    except subprocess.CalledProcessError:
        return None
