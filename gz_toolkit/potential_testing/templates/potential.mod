# NOTE: Adapted from LAMMPS examples/ELASTIC potential.mod.
# The pair style, pair coefficients and per-type masses all come from the
# per-case potential.inc written by gz_toolkit.

include         potential.inc

# Setup neighbor style
neighbor 1.0 bin
neigh_modify once no every 1 delay 0 check yes

# Setup minimization style
min_style       cg
min_modify      dmax ${dmax} line quadratic

# Setup output
thermo          1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
