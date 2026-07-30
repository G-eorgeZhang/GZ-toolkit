"""
submission.py - Script generation strategies for HPC jobs.

This module intentionally does NOT submit jobs. It only generates
.job scripts so users can review and submit manually.

Linux-only:
  Designed for Linux HPC workflows and SLURM-like script conventions.
"""

import os

from gz_toolkit.jobs.template_engine import JobTemplate


class ScatterSubmitter:
    """
    Generate one .job file per run directory (scatter style).
    """

    def __init__(self, directory_manager, template=None):
        self.dm = directory_manager
        self.template = template

    def build(
        self,
        run_command,
        selectors=None,
        output_name="slurm.job",
        template=None,
        jobname_func=None,
        dry_run=False,
        **template_kwargs,
    ):
        """
        Generate one script in each selected run directory.

        Parameters
        ----------
        run_command : str
            Command body to run. Supports ``{{DIR_PATH}}`` replacement.
        selectors : selectors type
            Passed to DirectoryManager.get_subdirs().
        output_name : str
            Job filename written in each run directory.
        template : JobTemplate or None
            Overrides ``self.template`` for this call.
        jobname_func : callable or None
            Optional: ``jobname_func(dir_name, dir_path) -> str``.
            Default job name is the run directory name.
        dry_run : bool
            If True, print would-be outputs without writing files.
        **template_kwargs
            Passed to JobTemplate.render().

        Returns
        -------
        list[str]
            Absolute paths of generated (or would-be generated) files.
        """
        tpl = template or self.template or JobTemplate()
        written = []

        for dir_name, dir_path in self.dm.iter_subdirs(selectors):
            jobname = (
                jobname_func(dir_name, dir_path)
                if jobname_func is not None
                else dir_name
            )
            cmd = run_command.replace("{{DIR_PATH}}", dir_path)
            text = tpl.render(
                jobname=jobname,
                run_command=cmd,
                **template_kwargs,
            )

            out_path = os.path.abspath(os.path.join(dir_path, output_name))
            if dry_run:
                print(f"--- [DRY RUN] {out_path} ---")
                print(text)
                print("---")
            else:
                with open(out_path, "w") as f:
                    f.write(text)
                print(f"Wrote scatter script: {out_path}")

            written.append(out_path)

        return written

    def submit(self, *args, **kwargs):
        """
        Backward-compatibility alias.
        This method now only generates scripts and never submits.

        Supported legacy mode:
            submit(sub_file_name, selectors=None, pre_submit_hook=None, dry_run=False)
        In that mode, it validates/optionally edits existing files and
        returns matching paths, but does not call sbatch.
        """
        if args and isinstance(args[0], str) and args[0].endswith(".job"):
            sub_file_name = args[0]
            selectors = kwargs.get("selectors", None)
            pre_submit_hook = kwargs.get("pre_submit_hook", None)
            dry_run = kwargs.get("dry_run", False)

            matched = []
            for dir_name, dir_path in self.dm.iter_subdirs(selectors):
                sub_path = os.path.join(dir_path, sub_file_name)
                if not os.path.isfile(sub_path):
                    print(
                        f"Skipping '{dir_name}': "
                        f"'{sub_file_name}' not found."
                    )
                    continue

                if pre_submit_hook is not None:
                    try:
                        pre_submit_hook(dir_name, dir_path, sub_path)
                    except Exception as e:
                        print(
                            f"pre_submit_hook error for '{dir_name}': {e}"
                        )

                if dry_run:
                    print(f"[DRY RUN] Would use existing script: {sub_path}")
                else:
                    print(f"Validated existing script: {sub_path}")

                matched.append(sub_path)

            print(
                "[ScatterSubmitter] submit() compatibility mode: "
                "no job was submitted."
            )
            return matched

        print(
            "[ScatterSubmitter] submit() no longer submits jobs. "
            "Generating scripts only."
        )
        return self.build(*args, **kwargs)


class PatchSubmitter:
    """
    Generate one or more patch scripts, each containing commands for
    multiple run directories.
    """

    def __init__(self, directory_manager, template=None):
        self.dm = directory_manager
        self.template = template

    def build(
        self,
        run_command,
        selectors=None,
        max_per_file=50,
        output_dir=None,
        output_prefix="patch",
        jobname_base=None,
        template=None,
        dry_run=False,
        **template_kwargs,
    ):
        """
        Generate one or more patch .job files.
        """
        if output_dir is None:
            output_dir = os.getcwd()
        os.makedirs(output_dir, exist_ok=True)

        if jobname_base is None:
            jobname_base = os.path.basename(self.dm.folder_path)

        tpl = template or self.template or JobTemplate()
        subdirs = list(self.dm.iter_subdirs(selectors))
        written = []

        for i in range(0, len(subdirs), max_per_file):
            chunk = subdirs[i: i + max_per_file]
            chunk_lo = i
            chunk_hi = i + len(chunk) - 1
            jobname = f"{jobname_base}_{chunk_lo}-{chunk_hi}"

            body_lines = []
            for _, dir_path in chunk:
                body_lines.append(f"cd {dir_path}")
                body_lines.append(run_command.replace("{{DIR_PATH}}", dir_path))
                body_lines.append("")
            body_text = "\n".join(body_lines)

            text = tpl.render(
                jobname=jobname,
                run_command=body_text,
                **template_kwargs,
            )

            out_name = f"{output_prefix}_{chunk_lo}-{chunk_hi}.job"
            out_path = os.path.abspath(os.path.join(output_dir, out_name))

            if dry_run:
                print(f"--- [DRY RUN] {out_path} ---")
                print(text)
                print("---")
            else:
                with open(out_path, "w") as f:
                    f.write(text)
                print(f"Wrote patch script: {out_path}")

            written.append(out_path)

        return written

    def submit(self, *args, **kwargs):
        """
        Backward-compatibility alias.
        This method now only generates scripts and never submits.
        """
        print(
            "[PatchSubmitter] submit() no longer submits jobs. "
            "Generating scripts only."
        )
        return self.build(*args, **kwargs)
