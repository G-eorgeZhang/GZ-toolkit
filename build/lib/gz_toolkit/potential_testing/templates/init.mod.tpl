# NOTE: Adapted from LAMMPS examples/ELASTIC init.mod.
# The simulation cell is read from structure.data (pre-built by gz_toolkit at
# the reference-relaxed lattice constant) instead of being generated with the
# `lattice` command, so the same template serves pure elements, random alloys,
# and hcp cells.

# Define the finite deformation size. Try several values of this variable
# to verify that results do not depend on it.
variable up equal __UP__

# Define the amount of random jiggle for atoms.
# This prevents atoms from staying on saddle points.
variable atomjiggle equal 1.0e-5

# metal units, elastic constants in GPa
units           metal
variable cfac equal 1.0e-4
variable cunits string GPa

# Define minimization parameters
variable etol equal __ETOL__
variable ftol equal __FTOL__
variable maxiter equal __MAXITER__
variable maxeval equal __MAXEVAL__
variable dmax equal 1.0e-2

boundary        p p p
box tilt large

# Read the pre-built cell. It is orthogonal, so convert it to triclinic
# immediately — the shear deformations in displace.mod (xy/xz/yz change_box)
# require a triclinic box.
read_data       structure.data
change_box all triclinic
