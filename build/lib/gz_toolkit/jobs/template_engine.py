"""
template_engine.py — Create and customise SLURM job files from templates.

A job file has three logical sections:
  1. **Header** — ``#SBATCH`` resource-request directives
  2. **Middle** — module loads and environment setup
  3. **Run**    — conda/venv activation and the actual execution command

Templates use ``{{PLACEHOLDER}}`` markers that get substituted at render time.
Cluster-specific defaults (cores per node, partition, etc.) are loaded from
CSV files shipped under ``gz_toolkit/jobs/cluster_info/``.
"""

import csv
import os


# Path to the built-in templates and cluster_info directories
_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
_TEMPLATES_DIR = os.path.join(_PACKAGE_DIR, "templates")
_MACHINE_TEMPLATES_DIR = os.path.join(_TEMPLATES_DIR, "machines")
_CLUSTER_INFO_DIR = os.path.join(_PACKAGE_DIR, "cluster_info")

# Registry for user-added templates (name -> filepath)
_USER_TEMPLATES = {}


def _walltime_to_secs(wt):
    """Convert a walltime string like '24:00:00' or '30-00:00:00' to seconds."""
    if "-" in wt:
        days_str, rest = wt.split("-", 1)
        days = int(days_str)
    else:
        days, rest = 0, wt
    parts = rest.split(":")
    h, m, s = int(parts[0]), int(parts[1]), int(parts[2])
    return days * 86400 + h * 3600 + m * 60 + s


def _secs_to_walltime(secs):
    """Convert seconds back to HH:MM:SS or D-HH:MM:SS string."""
    days = secs // 86400
    remainder = secs % 86400
    h = remainder // 3600
    m = (remainder % 3600) // 60
    s = remainder % 60
    time_str = f"{h:02d}:{m:02d}:{s:02d}"
    if days > 0:
        return f"{days}-{time_str}"
    return time_str


def load_cluster_info(cluster_name):
    """
    Load cluster partition info from the CSV file.

    Parameters
    ----------
    cluster_name : str
        Name of the cluster (matches a ``.csv`` file in ``cluster_info/``).
        Case-insensitive.  E.g. ``"ISAAC"`` loads ``ISAAC.csv``.

    Returns
    -------
    dict
        Maps partition key (e.g. ``"campus"``, ``"group"``) to a dict of
        ``{partition, qos, account, cores_per_node, max_walltime, key}``.
    """
    fname = f"{cluster_name.upper()}.csv"
    fpath = os.path.join(_CLUSTER_INFO_DIR, fname)
    if not os.path.isfile(fpath):
        raise FileNotFoundError(
            f"Cluster info file '{fpath}' not found.  "
            f"Create a CSV in gz_toolkit/jobs/cluster_info/{fname} with columns: "
            "partition,qos,account,cores_per_node,max_walltime,key"
        )

    presets = {}
    with open(fpath, newline="") as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        for row in reader:
            k = row["key"].strip()
            presets[k] = {
                "partition": row["partition"].strip(),
                "qos": row["qos"].strip(),
                "account": row["account"].strip(),
                "cores_per_node": int(row["cores_per_node"].strip()),
                "max_walltime": row["max_walltime"].strip(),
            }
    return presets


def list_available_machines():
    """
    Return available machine names discovered from:
      - templates/machines/*.job
      - cluster_info/*.csv
    """
    names = set()
    if os.path.isdir(_MACHINE_TEMPLATES_DIR):
        for fn in os.listdir(_MACHINE_TEMPLATES_DIR):
            if fn.lower().endswith(".job"):
                names.add(os.path.splitext(fn)[0].upper())

    if os.path.isdir(_CLUSTER_INFO_DIR):
        for fn in os.listdir(_CLUSTER_INFO_DIR):
            if fn.lower().endswith(".csv"):
                names.add(os.path.splitext(fn)[0].upper())

    return sorted(names)


def _resolve_template_path(template=None, machine=None):
    if template is not None:
        return template

    if machine is not None:
        mpath = os.path.join(_MACHINE_TEMPLATES_DIR, f"{machine.upper()}.job")
        if os.path.isfile(mpath):
            return mpath

    return os.path.join(_TEMPLATES_DIR, "slurm_default.job")


class JobTemplate:
    """
    Template-based SLURM job file generator.

    Parameters
    ----------
    template : str or None
        Path to a custom ``.job`` template file.
        If None, uses the built-in ``slurm_default.job``.
    cluster : str or None
        Name of a cluster whose info CSV is in ``cluster_info/``.
        When provided, :meth:`render` can auto-fill partition, qos,
        account, cores_per_node from the CSV.
    """

    def __init__(self, template=None, cluster=None, machine=None):
        # Load template text
        if machine is None:
            machine = cluster
        tpath = _resolve_template_path(template=template, machine=machine)
        if not os.path.isfile(tpath):
            raise FileNotFoundError(f"Template file '{tpath}' not found.")
        with open(tpath, "r") as f:
            self._template_text = f.read()

        # Load cluster presets (optional)
        self._cluster_presets = None
        self._cluster_name = cluster
        self._machine_name = machine
        profile_name = machine if machine is not None else cluster
        if profile_name is not None:
            try:
                self._cluster_presets = load_cluster_info(profile_name)
            except FileNotFoundError:
                self._cluster_presets = None

            # A named cluster/machine must resolve to SOMETHING — either a
            # machine template or a cluster-info CSV. Silently falling back
            # to the default template would hide typos (e.g. "ISACC").
            has_machine_template = os.path.isfile(
                os.path.join(_MACHINE_TEMPLATES_DIR, f"{str(profile_name).upper()}.job"))
            if self._cluster_presets is None and not has_machine_template:
                raise FileNotFoundError(
                    f"Unknown cluster/machine '{profile_name}': no machine "
                    f"template or cluster_info CSV found. "
                    f"Available: {list_available_machines()}"
                )

    # ------------------------------------------------------------------
    # Render
    # ------------------------------------------------------------------

    def render(
        self,
        jobname="job",
        partition_key=None,
        partition=None,
        qos=None,
        account=None,
        walltime="24:00:00",
        nodes=1,
        cores_per_node=None,
        ntasks=None,
        ntasks_per_node=None,
        module_commands=None,
        env_commands=None,
        run_command="echo 'No run command specified'",
        **extra_placeholders,
    ):
        """
        Substitute placeholders and return the rendered job script string.

        Parameters
        ----------
        jobname : str
            Job name for ``#SBATCH -J``.
        partition_key : str or None
            Key into the cluster CSV (e.g. ``"c"`` for campus, ``"g"`` for
            group).  If the cluster was set in ``__init__``, this auto-fills
            PARTITION, QOS, ACCOUNT, and caps walltime.
        account : str or None
            Override the account from cluster CSV.
        walltime : str
            Walltime string (``"HH:MM:SS"`` or ``"D-HH:MM:SS"``).
        nodes : int
            Number of nodes.
        ntasks : int or None
            Total MPI tasks.  If None, computed as
            ``nodes * cores_per_node`` (from cluster CSV) or ``nodes * 40``.
        ntasks_per_node : int or None
            Tasks per node.  If None, computed as ``ntasks // nodes``.
        module_commands : str or list[str] or None
            Shell commands for the "middle" section (module loads).
            If None, uses a sensible default for ISAAC.
        env_commands : str or list[str] or None
            Shell commands for environment activation
            (e.g. ``"conda activate /path/to/env"``).
            If None, this placeholder is left blank.
        run_command : str
            The actual execution command
            (e.g. ``"srun -n 40 python sim.py"`` or
            ``"srun -n 40 ~/bins/lmp_mpi -in in.lmp"``).
        **extra_placeholders
            Any additional ``{{KEY}}`` → value substitutions.

        Returns
        -------
        str
            The fully rendered job script text.
        """
        # Start with cluster defaults if available
        partition_val = partition or ""
        qos_val = qos or ""
        acct = account or ""
        cores = int(cores_per_node) if cores_per_node is not None else 40

        if self._cluster_presets is not None and partition_key is not None:
            if partition_key not in self._cluster_presets:
                raise ValueError(
                    f"Partition key '{partition_key}' not found in "
                    f"cluster '{self._cluster_name}'.  "
                    f"Available keys: {list(self._cluster_presets.keys())}"
                )
            preset = self._cluster_presets[partition_key]
            partition_val = partition_val or preset["partition"]
            qos_val = qos_val or preset["qos"]
            if not acct:
                acct = preset["account"]
            if cores_per_node is None:
                cores = preset["cores_per_node"]

            # Cap walltime
            max_secs = _walltime_to_secs(preset["max_walltime"])
            if _walltime_to_secs(walltime) > max_secs:
                print(
                    f"⚠ Walltime '{walltime}' exceeds max "
                    f"'{preset['max_walltime']}' — capping."
                )
                walltime = preset["max_walltime"]

        # Compute tasks
        if ntasks is None:
            ntasks = cores * nodes
        if ntasks_per_node is None:
            ntasks_per_node = ntasks // nodes if nodes > 0 else ntasks

        # Module commands (default for ISAAC)
        if module_commands is None:
            module_commands = (
                "module unload PE-intel\n"
                "module load anaconda3\n"
                "source $ANACONDA3_SH"
            )
        elif isinstance(module_commands, list):
            module_commands = "\n".join(module_commands)

        # Environment commands
        if env_commands is None:
            env_commands = ""
        elif isinstance(env_commands, list):
            env_commands = "\n".join(env_commands)

        # Build substitution dict
        subs = {
            "JOBNAME": str(jobname),
            "ACCOUNT": str(acct),
            "WALLTIME": str(walltime),
            "NODES": str(nodes),
            "NTASKS": str(ntasks),
            "NTASKS_PER_NODE": str(ntasks_per_node),
            "PARTITION": str(partition_val),
            "QOS": str(qos_val),
            "MODULE_COMMANDS": str(module_commands),
            "ENV_COMMANDS": str(env_commands),
            "RUN_COMMAND": str(run_command),
        }
        subs.update({k.upper(): str(v) for k, v in extra_placeholders.items()})

        # Substitute
        text = self._template_text
        for key, val in subs.items():
            text = text.replace("{{" + key + "}}", val)

        return text

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------

    def write(self, output_path, **kwargs):
        """
        Render and write the job script to a file.

        Parameters are the same as :meth:`render`.

        Returns
        -------
        str
            The absolute path of the written file.
        """
        text = self.render(**kwargs)
        output_path = os.path.abspath(output_path)
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as f:
            f.write(text)
        print(f"Wrote job file: {output_path}")
        return output_path

    # ------------------------------------------------------------------
    # Create multiple job files at once (the createJobFile.py workflow)
    # ------------------------------------------------------------------

    def create_job_set(
        self,
        partition_keys=None,
        name_prefix="job",
        output_dir=None,
        **common_kwargs,
    ):
        """
        Create one job file per partition key (like ``createJobFile.py``).

        Parameters
        ----------
        partition_keys : list[str] or None
            List of partition keys from the cluster CSV.
            If None and a cluster is loaded, uses all available keys.
        name_prefix : str
            Base name for the output files.  Files are named
            ``{key}-{name_prefix}.job``.
        output_dir : str or None
            Directory to write files into.  Defaults to cwd.
        **common_kwargs
            Passed through to :meth:`render` (walltime, nodes, etc.).

        Returns
        -------
        list[str]
            Paths of written files.
        """
        if output_dir is None:
            output_dir = os.getcwd()

        if partition_keys is None:
            if self._cluster_presets is None:
                raise ValueError(
                    "No cluster loaded and no partition_keys specified."
                )
            partition_keys = list(self._cluster_presets.keys())

        written = []
        for pk in partition_keys:
            out_name = f"{pk}-{name_prefix}.job"
            out_path = os.path.join(output_dir, out_name)
            self.write(out_path, partition_key=pk, **common_kwargs)
            written.append(out_path)

        return written

    # ------------------------------------------------------------------
    # Template registration
    # ------------------------------------------------------------------

    @staticmethod
    def register_template(name, filepath):
        """
        Register a user-defined template for reuse by name.

        Parameters
        ----------
        name : str
            Short name to reference later.
        filepath : str
            Absolute path to the template file.
        """
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"Template file '{filepath}' not found.")
        _USER_TEMPLATES[name] = filepath
        print(f"Registered template '{name}' → {filepath}")

    @staticmethod
    def get_registered_template(name):
        """Return the file path of a previously registered template."""
        if name not in _USER_TEMPLATES:
            raise KeyError(
                f"Template '{name}' not registered.  "
                f"Available: {list(_USER_TEMPLATES.keys())}"
            )
        return _USER_TEMPLATES[name]

    def __repr__(self):
        return (
            f"JobTemplate(cluster='{self._cluster_name}', machine='{self._machine_name}')"
        )
