Put one machine-specific template per file:
  <MACHINE>.job

Example:
  ISAAC.job
  SUMMIT.job

Use placeholders like:
  {{JOBNAME}} {{ACCOUNT}} {{WALLTIME}} {{NODES}}
  {{NTASKS}} {{NTASKS_PER_NODE}} {{PARTITION}} {{QOS}}
  {{MODULE_COMMANDS}} {{ENV_COMMANDS}} {{RUN_COMMAND}}

Then use:
  JobTemplate(machine="ISAAC")

Machine defaults for partition/qos/account/cores can be stored in:
  gz_toolkit/jobs/cluster_info/<MACHINE>.csv
