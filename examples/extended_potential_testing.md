# Extended potential benchmarks and migration

These suites are opt-in. The existing LAMMPS ELASTIC stress/strain method is
retained; no energy-volume/EOS fit replaces it. The method is general, but its
numerical accuracy still requires relaxed inputs, force convergence and strain
sensitivity checks. Automatic prototype generation remains BCC/FCC/HCP; reviewed
imported elastic cells support other structures and compounds.

All settings below belong inside the existing per-potential JSON. Examples are
fragments: keep your potential name, type_map, masses, potential files/lines,
reference structures and actual machine settings.

## Surface energies

Set workflow.include_surfaces=true. Defaults generate symmetric, relaxed slabs
for (100), (110), (111), and (112) of each single_elements host. Configure
surface_millers, surface_thickness and surface_vacuum (angstrom). Multiple
symmetry-distinct symmetric terminations are retained. Converge slab thickness,
vacuum and lateral dimensions before treating an energy as a validated value.

The reported value is (E_slab - sum(N_i*mu_i))/(2*A), converted from eV/angstrom^2
to J/m^2. Unary bulk references supply mu; absent references do not become zero.
For an imported stoichiometric ceramic slab, use an explicit surface case and
supply physically justified elemental chemical potentials in eV/atom:

```json
{"workflow": {"benchmark_cases": [
  {"kind": "surface", "name": "SiC_001", "system": "SiC",
   "structure_file": "reffiles/SiC_001_slab.data", "equivalent_surfaces": true,
   "chemical_potentials": {"Si": -5.0, "C": -9.0}}
]}}
```

Those chemical potentials are illustrative placeholders, not reference data.
Choose compatible reservoirs and units for your potential. Slab faces must be
parallel to the first two cell vectors, with two equivalent surfaces separated
by vacuum. An asymmetric/polar slab's two surface energies cannot generally be
separated by simply dividing its excess energy by two. Imported benchmark files
currently use atomic atom style; charged surface benchmarks require extending
that path rather than silently dropping charges.

## C15 and SFT clusters

BCCDefect.add_c15_cluster(nsias, lattice_const=..., ff_elements=...,
atomic_masses=...) wraps myLAMMPS.create_C15_BCC. It uses the Z12 seed coordinates
from the supplied Examples archive. Optional template_0/template_1 accept
reviewed replacement POSCARs. This implementation currently requires a pure,
axis-aligned BCC host and verifies its BCC registry. It checks exact net atom
excess, cell preservation, finite coordinates and minimum periodic separation
before updating the host. Rejection leaves the original host unchanged.

The count is NET interstitials, not all atoms inside the C15 cluster. Successful
construction is a geometric seed, not evidence that the model stabilizes C15.
Relax it, check its structure, and compare formation energies across sizes.
C15 interstitial clusters are established for BCC Fe; this does not justify a
generic HCP constructor. The legacy backend is retained, not claimed to be a
newly validated optimal cluster-growth algorithm.

FCCDefect.stacking_fault_tetrahedron already removes a triangular, single-layer
{111} vacancy platelet. It reports relaxed=false. Relaxation/annealing can form
an SFT, a loop, or another cluster depending on potential, size and temperature;
inspect the relaxed topology before identifying it as an SFT.

```json
{"workflow": {"benchmark_cases": [
  {"kind": "c15", "name": "C15_I8", "element": "Fe", "size": 12, "nsias": 8},
  {"kind": "sft", "name": "SFT_seed", "element": "Cu", "size": 12,
   "edge_length_A": 20, "anneal": {"temperature_K": 300, "steps": 10000}}
]}}
```

Only use the entries compatible with your configured host structure and
potential. SFT element must be FCC, C15 element BCC. Increasing cluster sizes
requires increasing the cell to avoid image interactions. Cluster summaries
report formation energy and energy per net defect; SFT topology is not inferred
from its energy or from a force-converged minimization alone.

## Dislocation lines

The existing BCCDefect/FCCDefect edge and screw constructors are used. A line
benchmark needs a supplied, correctly oriented pristine structure, construction
arguments, and a matched reviewed reference with the SAME final cell geometry.
This makes boundary conditions and surface/strain contributions explicit rather
than guessing a core energy from a bulk energy difference.

```json
{"workflow": {"benchmark_cases": [
  {"kind": "dislocation", "element": "Fe", "name": "edge_100_pair",
   "structure_file": "reffiles/oriented_pristine.data",
   "reference_file": "reffiles/matched_reference.data",
   "line_kind": "edge", "line_axis": 0, "n_lines": 2,
   "boundary": "p p p",
   "construction": {"burgers": [0,1,0], "glide_plane": [0,0,1], "nedges": 2}}
]}}
```

The example assumes the specified directions already coincide with the box
axes. For other orientations provide current_axes/newaxis through construction
as supported by the defect API. Screw cases accept nscrews and handle_pbc.
Single periodic dislocations need a compatible tilted cell or reviewed free
boundaries; a net Burgers vector must not be hidden in a naive periodic cell.

The finite-cell line excess is
(E_def - E_ref - sum((N_def_i-N_ref_i)*mu_i))/(n_lines*L), in eV/angstrom.
Both defect and reference must converge under the intended boundary conditions.
fixed_ids/reference_fixed_ids can hold specified atoms during minimization.
Surface corrections are embodied in the matched reference, not automatically
borrowed from a differently shaped surface calculation. This value contains
finite-size elastic/image contributions and is NOT an isolated core energy.
Changing cell size, boundary conditions, or the line count changes the test.

## Dislocation loops

Keep include_dislocation_loops=true. loop_radius_factors selects a size scan;
loop_types can contain sil and/or vl (SIA/vacancy loops). Both BCC <111> and <100>
families are generated. Empty radius factors preserve the old single-radius
setting and directory names. Summaries add formation energy per net defect.

```json
{"workflow": {"include_dislocation_loops": true,
              "loop_radius_factors": [1.5, 2.0, 3.0], "loop_types": ["sil", "vl"]}}
```

## Gas complexes in alloy matrices

```json
{"workflow": {"benchmark_cases": [
  {"kind": "alloy_gas", "name": "He2V1", "system": "FeCrHe",
   "composition": {"Fe": 90, "Cr": 10}, "gas": "He",
   "n_gas": 2, "m_defect": 1, "defect_kind": "vacancy", "size": 12, "seed": 1729}
]}}
```

The pristine matrix is decorated BEFORE defects are created. Four cases share
one realization: bulk, defect-only, gas-only, complex. Binding energy is
E_complex + E_bulk - E_defect - E_gas (negative means attraction). This avoids
comparing unrelated random alloys or inferring alloy elemental chemical
potentials from a bulk average. Average multiple explicit seeds for disorder
statistics. Gas positions are initial guesses; inspect relaxed configurations
and compare alternative sites. Additional interstitial defect cases use
"defect_kind": "interstitial".

## Molecular references

```json
{"workflow": {"benchmark_cases": [
  {"kind": "molecule", "name": "H2", "elements": ["H","H"], "distance_A": 0.75, "box_A": 40},
  {"kind": "molecule", "name": "He2", "elements": ["He","He"], "distance_A": 2.0, "box_A": 40}
]}}
```

Each dimer gets same-box isolated-atom references. Outputs include relaxed
minimum-image bond distance, binding E_dimer-sum(E_atoms), and its negative
(dissociation energy). A repulsive/unbound dimer is a meaningful failure of a
molecular-binding test, not a reason to assert it has a stable bond. Check
box/cutoff sensitivity, endpoint separation, and whether this potential is
intended to model molecular physics. The sample starting distances are explicit
user-editable choices, not claimed optimized references for your potential.

## SEAKMC discovery and NEB prescribed paths

Both are available; they solve different tasks. SEAKMC discovers local activated
events without specifying their final states. LAMMPS climbing-image NEB tests a
particular path between supplied relaxed endpoints. A failed/incompatible
SEAKMC backend cannot automatically be replaced with a valid NEB path.

For SEAKMC keep seakmc_enabled=true and the original seakmc_defects. To include
other cases, add seakmc_case_patterns relative to the potential directory, e.g.
["gas_complexes/*/*", "benchmarks/Fe/C15_I8"], and set an explicit
seakmc_dactive in angstrom appropriate to the defect size. seakmc_searches
controls NSearch. Review the full input template, defect detection and active
volume settings for extended clusters; the default is still a starting point.

At runtime run_seakmc_p.py requires the parent's cleanly completed relaxed.data
and copies it plus potential files immediately before the search. It never
falls back to an unrelaxed structure. Dependencies cover the selected parent
cases. Summaries retain minimum forward/backward barriers, number of CSV saddle
records and reported search counts when available. Counts across successive
searches are not guaranteed counts of unique physical events.

For NEB:

```json
{"workflow": {"neb_replicas": 8, "neb_force_tolerance": 0.0001,
 "neb_iterations": 5000, "neb_cases": [
  {"name": "vacancy_hop", "system": "Fe", "initial": "reffiles/start.data",
   "final": "reffiles/end.data", "endpoints_relaxed": true}
]}, "hpc": {"ntasks": 8}}
```

Endpoints must already be reviewed/minimized, use the same box, atom IDs and
species, and represent distinct basins. NEB must not be made by independently
renumbering endpoint atoms. The wrapper writes LAMMPS's ID/coordinate endpoint
file, uses fire + fix neb, and launches with -partition NxR. ntasks must be
divisible by the replica count. It requires a LAMMPS build with REPLICA.
Results are accepted only after matching inputs, all image files, a clean
finish, and maximum replica force norm <= the configured force tolerance.
Inspect image continuity, saddle location and replica-count convergence;
linear interpolation is only a starting path. No physical barrier is validated
by unit tests of the input generator/parser.

## CPU preparation, GPU force evaluation, and potential backends

Python geometry, files and summaries remain CPU work. A GPU node also supplies
CPU cores; allocating a GPU does not make NumPy or ordinary Python run on it.
Use a Python environment with the toolkit and its dependencies on the cluster.
The model's LAMMPS pair style or in-process force evaluator decides whether it
can and should use a GPU.

Potential lines/files remain explicit, so classical models, MTP/MLIP plugins
and universal ML pair styles can use the same benchmark geometries. This is a
LAMMPS interface contract, not a claim that every model/build supports every
calculation (isolated atoms, molecules, minimization, stress, species, charge).
Python-only calculators without a working LAMMPS or SEAKMC adapter are not
silently treated as supported LAMMPS models.

hpc.device="gpu" plus gpus_per_node requests GPUs. Set the actual modules,
conda_env, lammps_binary_path, python_executable, run_command and lammps_args
required by YOUR installed backend. For example lammps_args="-k on g 1 -sf kk"
is appropriate only for a compatible Kokkos-enabled pair style/build; a GPU
allocation alone does not enable it. MTP/MLIP implementations may be CPU-only.
Do not automatically add a generic -sf gpu to universal potentials.

Degrees 3/4 allow separate preparation/summary allocations without GPU requests:
set preprocess_partition_key to the site's CPU partition key and
preprocess_ntasks (default 1). Without a CPU partition key, it retains the
configured partition, whose scheduling policy you must check. The same modules
and environment are used. Degrees 1/2 run Python CPU work inside their shared
compute allocation; this works but can leave the GPU idle during preparation.
Generated jobs now stop on command failures instead of hiding them behind date.

GPU SEAKMC needs an explicit, tested seakmc_template_path and, where needed,
seakmc_run_command for its MPI/in-process LAMMPS force evaluator. An executable
GPU LAMMPS working outside SEAKMC does not establish this compatibility.
GPU NEB requires degree 3/4 and an explicit neb_run_command with replica-aware
GPU binding. There is no generic automatic GPU-to-replica mapping; avoid every
replica accidentally selecting GPU 0. All launch commands remain reviewable.

## Imported elastic cells

Keep the LAMMPS ELASTIC method and supply optional workflow.elastic_cases:

```json
{"workflow": {"include_elastic": true, "elastic_cases": [
  {"name": "SiC", "structure_file": "reffiles/relaxed_SiC.data", "atom_style": "atomic"}
]}}
```

The complete supplied cell is copied without geometry/chemistry regeneration.
atomic, charge and full input styles are supported for this imported elastic
path; provide matching potential commands, charges, topology and required
long-range settings. These cells do not become elemental chemical-potential
references for other benchmarks. This extension does not make the existing
point-defect constructors valid for arbitrary ceramic lattices.

## Results and execution

New cases live in benchmarks/<system>/<name>/ or neb/<system>/<name>/ and are
included in all four scheduling tiers, pre-flight file checks and status.
Summaries report explicit status; failed or non-force-converged runs contribute
no accepted energy metric. Default benchmark force tolerance is 1e-6 eV/angstrom
with etol=0; configure benchmark_force_tolerance for the study. Matching input
fingerprints prevent acceptance of stale logs after regenerating a case.

Preparation writes files only. Run the existing launcher only after the campaign
and its resources have been authorized. No tests in this change run LAMMPS,
SEAKMC or NEB calculations or establish physical validation of a potential.

## Method references

- C15 interstitial clusters in BCC Fe: Byggmastar and Granberg,
  [Dynamical stability of radiation-induced C15 clusters in iron](https://arxiv.org/abs/1909.09818).
- [LAMMPS NEB](https://docs.lammps.org/neb.html) and
  [fix neb](https://docs.lammps.org/fix_neb.html): endpoints, replica launch,
  climbing-image optimization and force-norm convergence.
- [LAMMPS GPU package](https://docs.lammps.org/Speed_gpu.html): GPU acceleration
  depends on the force-evaluation implementation and installed build.

C15 Z12 seed coordinate files in gz_toolkit/defect/templates were taken from
the user-supplied Examples.tgz (INPUTs/C15s_paper/Z12_0.POSCAR and Z12_1.POSCAR).
The wrapper calls the installed myLAMMPS implementation; it does not bundle a
copy of that implementation or claim its output has been physically validated.

For imported elastic files, coefficient sections are skipped with read_data
nocoeff; define all force-field coefficients in potential.inc/pot_lines. An
already triclinic supplied cell is preserved without converting it a second
time. The generated elasticity/benchmark scripts use metal units; force
interfaces must return compatible eV/angstrom-based quantities.
