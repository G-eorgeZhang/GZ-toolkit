# Cluster info

Each cluster gets its own directory: `<CLUSTER>/<CLUSTER>.csv` (partitions) and
an optional `<CLUSTER>/meta.json` (machine-level facts). Both are looked up by
`cluster` name (case-insensitive), e.g. `"ISAAC"` -> `ISAAC/ISAAC.csv` +
`ISAAC/meta.json`.

## `<CLUSTER>.csv` — partitions available on this cluster

| column            | meaning                                            |
|-------------------|----------------------------------------------------|
| `partition`       | SLURM partition name (`#SBATCH --partition=...`)   |
| `qos`             | SLURM QoS                                          |
| `account`         | SLURM account                                      |
| `cores_per_node`  | physical cores per node (used for `--ntasks-per-node`) |
| `max_walltime`    | hard walltime cap (HH:MM:SS or D-HH:MM:SS)         |
| `key`             | short alias (used as `partition_key` in HPC config)|

## `meta.json` — machine-level facts (not partition-specific)

```json
{
  "path2gz_toolkit": "/lustre/isaac/proj/.../gz_toolkit"
}
```

`path2gz_toolkit` is the install location of this `gz_toolkit` checkout on
that cluster. It's defined once per cluster (not per potential) because it
never changes across pots run on the same machine. Used by
`gz_toolkit.pot_infobank.promote_tag_json` to resolve where to copy a
potential's `<pot_name>.json` when `PotentialMetadata.promote_to_infobank = True`
(see `gz_toolkit/pot_infobank/__init__.py`). Leave it `null`/omit the file if
this cluster never promotes.

## Adding a GPU partition

`gz_toolkit.potential_testing` recognises `hpc.device == "gpu"` and `hpc.gpus_per_node`
and injects a `#SBATCH --gres=gpu:<N>` line automatically. You only need to add
a row that points at the GPU partition; for example, on ISAAC:

```csv
gpu_short,gpu,ACF-UTK0035,32,3:00:00,gs
```

Then in your per-potential JSON:

```json
"hpc": {
  "cluster": "ISAAC",
  "partition_key": "gs",
  "device": "gpu",
  "gpus_per_node": 1,
  "modules": ["module load cuda/12.1"],
  "conda_env": "/path/to/mlpot_env",
  "run_command": "srun -n {{NTASKS}} lmp_gpu -in {{INPUT}} -sf gpu -pk gpu 1"
}
```
