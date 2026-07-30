"""
Examples for gz_toolkit.jobs (Linux HPC workflow).

This file demonstrates all currently supported jobs APIs:
  - generate_folder_tree
  - DirectoryManager
  - list_available_machines
  - load_cluster_info
  - JobTemplate
  - ScatterSubmitter
  - PatchSubmitter
"""

import re

from gz_toolkit.jobs import (
    DirectoryManager,
    JobTemplate,
    ScatterSubmitter,
    PatchSubmitter,
    generate_folder_tree,
    list_available_machines,
    load_cluster_info,
)


def example_generate_tree():
    # Default: use pwd as base folder.
    generate_folder_tree()

    # Explicit path override.
    generate_folder_tree(
        path="/home/user/my-runs",
        output_file="tree.txt",
        max_depth=5,
        include_files=True,
    )


def example_directory_manager():
    dm = DirectoryManager("/home/user/my-runs", naming="prefix-suffix")
    all_subdirs = dm.get_subdirs()
    one = dm.get_subdirs("12-Fe")
    some = dm.get_subdirs(["12-Fe", "13-He"])
    by_range = dm.get_subdirs((10, 20))
    by_glob = dm.get_subdirs("glob:*-Fe")
    iter_pairs = list(dm.iter_subdirs((10, 12)))
    one_path = dm.subdir_path("12-Fe")
    _ = (all_subdirs, one, some, by_range, by_glob, iter_pairs, one_path)

    # Other naming modes:
    dm_flat = DirectoryManager("/home/user/my-runs-flat", naming="flat")
    dm_regex = DirectoryManager(
        "/home/user/temp-runs",
        naming=re.compile(r"temp_(?P<key>\d+)"),
    )
    _ = (dm_flat, dm_regex)


def example_machine_profiles():
    machines = list_available_machines()
    print("Available machines:", machines)

    isaac = load_cluster_info("ISAAC")
    print("ISAAC keys:", list(isaac.keys()))


def example_job_template():
    # 1) Default template
    tpl_default = JobTemplate()
    script_text = tpl_default.render(
        jobname="myjob",
        run_command="srun -n 40 lmp_mpi -in in.lmp",
    )
    print(script_text[:120])

    # 2) Machine template + machine CSV defaults
    tpl_machine = JobTemplate(machine="ISAAC")
    tpl_machine.write(
        "/home/user/scripts/isaac-c.job",
        partition_key="c",
        jobname="my-isaac-job",
        nodes=1,
        run_command="srun -n 40 lmp_mpi -in in.lmp",
    )

    # 3) Override partition/qos/account/cores directly
    tpl_machine.write(
        "/home/user/scripts/custom.job",
        jobname="custom-job",
        partition="campus",
        qos="campus",
        account="ACF-UTK0035",
        nodes=2,
        cores_per_node=40,
        run_command="srun -n 80 lmp_mpi -in in.lmp",
    )

    # 4) One job file per partition key
    tpl_machine.create_job_set(
        partition_keys=["c", "s", "g"],
        name_prefix="lammps",
        output_dir="/home/user/scripts",
        run_command="srun -n 40 lmp_mpi -in in.lmp",
    )

    # 5) Register and reuse user templates
    JobTemplate.register_template(
        "my-template",
        "/home/user/templates/my_machine.job",
    )
    tpath = JobTemplate.get_registered_template("my-template")
    JobTemplate(template=tpath).write(
        "/home/user/scripts/from-registered.job",
        jobname="reg-job",
        run_command="srun -n 16 python run.py",
    )


def example_scatter_and_patch():
    dm = DirectoryManager("/home/user/my-runs", naming="prefix-suffix")
    tpl = JobTemplate(machine="ISAAC")

    scatter = ScatterSubmitter(dm, template=tpl)
    scatter.build(
        run_command="srun -n 40 lmp_mpi -in in.lmp",
        selectors=(0, 99),
        output_name="slurm.job",
        partition_key="c",
    )

    # Backward-compatible alias (still generation-only).
    scatter.submit(
        run_command="srun -n 40 lmp_mpi -in in.lmp",
        output_name="slurm.job",
        partition_key="c",
    )

    patch = PatchSubmitter(dm, template=tpl)
    patch.build(
        run_command="srun -n 40 lmp_mpi -in in.lmp",
        selectors="glob:*-Fe",
        max_per_file=25,
        output_dir="/home/user/patch-jobs",
        output_prefix="patch",
        partition_key="c",
    )

    # Backward-compatible alias (still generation-only).
    patch.submit(
        run_command="srun -n 40 lmp_mpi -in in.lmp",
        max_per_file=25,
        output_dir="/home/user/patch-jobs",
        partition_key="c",
    )


if __name__ == "__main__":
    # Run only examples you need.
    # Uncomment to execute in your Linux environment.
    # example_generate_tree()
    # example_directory_manager()
    # example_machine_profiles()
    # example_job_template()
    # example_scatter_and_patch()
    pass
