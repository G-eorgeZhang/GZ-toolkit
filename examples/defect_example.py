"""Examples for EVERY defect-creation method in BCCDefect (gz_toolkit.defect.bcc_defect).

Each call lists every parameter on its own line, with a comment explaining it.
Replace file names / coordinates with your own, then call the example
functions you want from __main__ at the bottom.

Common setup used by all examples::

    data = Modlmp_LmpData.from_file("bcc_supercell.data", "atomic")
    defect = BCCDefect(data)
"""

import numpy as np

from gz_toolkit import Modlmp_LmpData, BCCDefect

# --- shared constants (bcc Fe) ---------------------------------------------
A0 = 2.8553                       # lattice constant (Å)
FF_ELEMENTS = ["Fe", "Fe"]        # force-field element per atom type (type 1, type 2)
ATOMIC_MASSES = [55.845, 55.845]  # masses (amu) matching FF_ELEMENTS


########################################################################################################################
# 1. POINT DEFECTS
########################################################################################################################

def ex_vacancy(defect: BCCDefect):
    """Delete atoms sitting at given coordinates."""
    n = defect.vacancy(
        coords=[14.27, 14.27, 14.27],  # one (x,y,z) or a list of (x,y,z) — Cartesian Å positions of atoms to remove
        tolerance=1e-3,                # per-axis match tolerance (Å): atom is deleted if |x-x0|,|y-y0|,|z-z0| all within this
        reinit=True,                   # refresh lmpData bookkeeping after the edit (set False when chaining many edits)
    )
    print(f"vacancy: deleted {n} atoms")


def ex_interstitial(defect: BCCDefect):
    """Insert new atoms at explicit positions."""
    n = defect.interstitial(
        sites=[[1.43, 1.43, 0.0],
               [4.28, 1.43, 0.0]],    # one (x,y,z) or list of (x,y,z) — Cartesian Å positions for the new atoms
        atom_type=2,                  # LAMMPS atom type assigned to every inserted atom
        wrap=False,                   # True: wrap positions back into the periodic box before inserting
        reinit=True,                  # refresh bookkeeping afterwards
    )
    print(f"interstitial: added {n} atoms")
    # interstitial_cluster(...) and impurity_atom(...) take the same parameters —
    # they are aliases of interstitial() for readability.


def ex_random_interstitials(defect: BCCDefect):
    """Add random tetra/octa-site interstitials (alias: random_interstitials)."""
    n = defect.add_interstitial_bcc(
        fraction=0.001,               # number of interstitials = int(fraction * current atom count)
        atom_type=2,                  # atom type of the inserted interstitials
        min_dist=1.2,                 # minimum allowed distance (Å) to ANY existing/new atom (min-image); candidates closer are rejected
        max_trials_per_atom=20,       # random attempts per interstitial before giving up on it
        rng=np.random.default_rng(42),# numpy Generator for reproducibility; None = fresh random seed
        reinit=True,                  # refresh bookkeeping afterwards
    )
    print(f"random interstitials: placed {n}")


def ex_dumbbell(defect: BCCDefect):
    """Replace one lattice atom with a two-atom dumbbell."""
    n = defect.dumbbell(
        center=[14.27, 14.27, 14.27],  # Cartesian Å position of the lattice atom to convert into a dumbbell
        atom_type=1,                   # atom type of the two dumbbell atoms
        direction=(1, 1, 1),           # dumbbell axis (any 3-vector; normalized internally)
        separation=1.0,                # distance (Å) between the two dumbbell atoms (each sits ±separation/2 from center)
        remove_center_atom=True,       # delete the original atom at center first (False: just add 2 atoms on top)
        tolerance=1e-3,                # coordinate-match tolerance (Å) used when removing the center atom
        wrap=True,                     # wrap the two new positions into the periodic box
        reinit=True,                   # refresh bookkeeping afterwards
    )
    print(f"dumbbell: created {n} atoms")
    # Shortcuts with the axis pre-set (all other params identical):
    #   defect.dumbbell_111(center, atom_type=1)   # direction=(1,1,1)
    #   defect.dumbbell_100(center, atom_type=1)   # direction=(1,0,0)
    #   defect.dumbbell_110(center, atom_type=1)   # direction=(1,1,0)


def ex_impurities(defect: BCCDefect):
    """Re-type atoms in a region (substitutional impurities)."""
    n = defect.impurities(
        "sphere",                      # region_type: any region_mask key — "sphere", "cube", "cylinder", "plane", "type", "pct", "atom_id", ...
        new_type=2,                    # atom type the selected atoms are changed to (must exist in Masses)
        ff_elements=["Fe", "He"],      # element symbols per type — asserted into the force field before re-typing
        atomic_masses=[55.845, 4.003], # masses (amu), same length/order as ff_elements
        center=[14.27, 14.27, 14.27],  # region kwarg (sphere): Cartesian center
        radius=5.0,                    # region kwarg (sphere): radius (Å)
        # source_region_type="type",   # optional 2nd region: final selection = region ∩ source region
        # source_kwargs={"types": 1},  # kwargs for the source region (here: only atoms of type 1)
        # reinit=True,                 # refresh bookkeeping afterwards (default True)
    )
    print(f"impurities: re-typed {n} atoms")
    # Random alloying example — re-type 5% of all atoms:
    #   defect.impurities("pct", new_type=2, ff_elements=[...], atomic_masses=[...],
    #                     pct=0.05,                # fraction (0,1] of eligible atoms to pick randomly
    #                     base_region="all")       # region the percentage is drawn from


########################################################################################################################
# 2. LARGER DEFECTS / SINKS
########################################################################################################################

def ex_void(defect: BCCDefect):
    """Carve out a void (sphere / cube / cylinder)."""
    n = defect.void(
        center=None,                   # Cartesian (x,y,z) center; None = geometric center of the system
        radius=8.0,                    # sphere radius (Å); for cube: half side-length; for cylinder: circle radius
        region_type="sphere",          # carve shape: "sphere", "cube"/"box", or "cylinder"
        reinit=True,                   # refresh bookkeeping afterwards
    )
    print(f"void: deleted {n} atoms")
    # Cylinder variant needs two extra kwargs:
    #   defect.void(center=None, radius=8.0, region_type="cylinder",
    #               plane="xy",        # circle plane of the cylinder: "xy", "xz", or "yz"
    #               height=20.0)       # ± half-thickness (Å) about the center along the cylinder axis


def ex_gas(defect: BCCDefect):
    """Insert gas atoms — single interstitial or gas-filled bubble."""
    n = defect.gas(
        mode="interstitial",           # "interstitial"/"single": just insert gas atom(s); "bubble"/"void": carve sphere + 1 gas atom
        atom_type=2,                   # atom type of the gas atom(s) (e.g. He)
        sites=[14.27, 15.70, 14.27],   # gas position(s) for interstitial mode (or pass center=... instead)
    )
    print(f"gas interstitial: {n}")

    report = defect.gas(
        mode="bubble",                 # carve a spherical void, then place one gas atom at its center
        atom_type=2,                   # gas atom type
        center=[28.55, 28.55, 28.55],  # bubble center (Cartesian Å)
        radius=5.0,                    # bubble radius (Å) — all matrix atoms inside are removed
    )
    print(f"gas bubble: {report}")     # dict: {"removed_atoms": ..., "gas_atoms": 1}


def ex_frenkel_pairs(defect: BCCDefect):
    """Create vacancy–interstitial (Frenkel) pairs by displacing random atoms."""
    report = defect.add_FrenkelPair(
        10,                              # n: how many Frenkel pairs to attempt
        min_sep_from_vac=1.5,            # new interstitial must land ≥ this (Å) from its own vacancy site
        min_sep_from_atoms=1.3,          # ... and ≥ this (Å) from every unmoved atom (min-image; KD-tree if box orthogonal)
        min_sep_between_new=1.3,         # ... and ≥ this (Å) from every previously placed interstitial
        offset_range=(10.0, 20.0),       # (lo, hi) Å — displacement distance drawn uniformly in this range
        max_trials=200,                  # random direction/distance attempts per pair before skipping it
        rng=np.random.default_rng(7),    # numpy Generator for reproducibility; None = fresh seed
        eligible_region_type="cube",     # optional: only atoms inside this region may be chosen (e.g. avoid surfaces)
        eligible_region_kwargs={         # kwargs for the eligible region
            "cube_params": {"center": (28.5, 28.5, 28.5), "side_length": 40.0}},
        exclude_region_type=None,        # optional: atoms inside this region are forbidden (e.g. frozen slab)
        exclude_region_kwargs=None,      # kwargs for the exclude region
        atom_mask=None,                  # optional boolean mask (len = n_atoms) ANDed with the region masks
        wrap_new=True,                   # wrap displaced positions back into the periodic box
    )
    print(f"Frenkel: moved {len(report['moved_indices'])}, skipped {len(report['skipped_indices'])}")


########################################################################################################################
# 3. DISLOCATION LOOPS (prismatic, <111> and <100>)
########################################################################################################################

def ex_111_interstitial_loop(defect: BCCDefect):
    """<111> SIL: insert a 3-layer (111) platelet (copy the slab, split ±b/2)."""
    n = defect.add_111_loop(
        loop_type="sil",                   # "sil"/"interstitial" (insert) or "vl"/"vacancy"
        center=[57.106, 57.106, 57.106],   # Cartesian Å center of the loop (REQUIRED)
        habit_plane=(1, 1, 1),             # <111> variant, e.g. (1,-1,1); b = ½<111> -> 3 layers
        radius=10.0,                       # loop radius (Å) in the habit plane
        lattice_const=A0,                  # bcc a (Å) — sets |b| and (111) layer spacing a/(2√3)
        ff_elements=FF_ELEMENTS,           # element symbols per atom type (asserted into force field)
        atomic_masses=ATOMIC_MASSES,       # masses (amu) matching ff_elements
        box_axes=None,                     # matrix coord system (crystal dirs on box x/y/z); None = cubic-aligned
        normal_box=None,                   # OR give the box-frame normal directly (e.g. after a swap_axes)
        shape="round",                     # "round" (disk) or "square"
        n_repeats=1,                       # platelet thickness in Burgers repeats (1 => 3 layers for <111>)
        loop_atom_type=2,                  # type for the NEW atoms; None = inherit each atom's type
        select_types=None,                 # restrict selection to these types/symbols; None = all
        wrap=True,                         # wrap moved/new positions into the periodic box
        reinit=True,                       # refresh bookkeeping afterwards
    )
    print(f"<111> SIL: inserted {n} atoms (3 layers)")


def ex_111_vacancy_loop(defect: BCCDefect):
    """<111> vacancy loop: delete the 3-layer disc + move columns toward the gap."""
    n = defect.add_111_loop(
        loop_type="vl",                    # vacancy branch (separate code path from SIL)
        center=[57.106, 57.106, 57.106],   # Cartesian Å center (REQUIRED)
        habit_plane=(1, 1, 1),             # <111> variant
        radius=10.0,                       # loop radius (Å)
        lattice_const=A0,                  # bcc lattice constant a (Å)
        ff_elements=FF_ELEMENTS,           # element symbols per atom type
        atomic_masses=ATOMIC_MASSES,       # masses (amu)
        n_repeats=1,                       # layers deleted = k*n_repeats (3 for <111>)
        closure_layers=2,                  # intact layers on each side nudged toward the gap
        closure_disp_frac=0.15,            # nudge distance = frac * a per closure layer
        select_types=None,                 # optional type filter
        wrap=True,                         # wrap moved positions into the box
        reinit=True,                       # refresh bookkeeping
    )
    print(f"<111> vacancy loop: {n} atoms deleted")


def ex_100_loop(defect: BCCDefect):
    """<100> loop — same call, only the family differs (2 layers; b = a<100>)."""
    n = defect.add_100_loop(
        loop_type="sil",                   # "sil" or "vl"
        center=[57.106, 57.106, 57.106],   # Cartesian Å center (REQUIRED)
        habit_plane=(1, 0, 0),             # <100> variant: (±1,0,0)/(0,±1,0)/(0,0,±1); 2 layers, spacing a/2
        radius=10.0,                       # loop radius (Å)
        lattice_const=A0,                  # bcc lattice constant a (Å)
        ff_elements=FF_ELEMENTS,           # element symbols per atom type
        atomic_masses=ATOMIC_MASSES,       # masses (amu)
        n_repeats=1,                       # 1 => 2 layers for <100>
        loop_atom_type=2,                  # type for new atoms; None = inherit
        select_types=None,                 # optional type filter
        wrap=True,                         # wrap new/moved positions
        reinit=True,                       # refresh bookkeeping
    )
    print(f"<100> loop: {n}")
    # Oriented cells: pass box_axes (crystal dirs on box x/y/z) so the requested
    # plane is rotated into the box frame — e.g. for a <111> screw/edge cell:
    #   defect.add_111_loop(loop_type="sil", box_axes=[[1,1,1],[1,1,-2],[1,-1,0]], ...)


########################################################################################################################
# 4. DISLOCATION LINES (edge / screw, myLAMMPS backend)
########################################################################################################################

def ex_edge_dislocation(defect: BCCDefect):
    """Edge dislocation: remove a half-plane + shrink the box to close it.

    The data must already be a BUILT SUPERCELL (read seed, scale_data(a),
    make_supercell) — construction stays outside this method.
    """
    report = defect.edge_dislocation(
        ff_elements=["Fe", "Fe"],           # element symbols per atom type (asserted into force field)
        lattice_const=A0,                   # lattice constant a (Å) — auto-computes |b| from `burgers` (½<111> → a√3/2)
        burgers=(1, 1, 1),                  # Burgers vector DIRECTION (½<111> assumed for magnitude)
        glide_plane=(1, 1, 0),              # glide-plane normal: (1,1,0) → <111>{110} edge; (1,1,2) → <111>{112} edge
        current_axes=[[1, 1, 0],            # crystallographic direction currently on x
                      [1, 1, 2],            # ... on y
                      [1, 1, 1]],           # ... on z — used to auto-derive the axis permutation (here: [1,2,0])
        newaxis=None,                       # explicit swap_axes permutation override, e.g. [1,2,0]; skips auto-derivation
        burgerm=None,                       # explicit |b| (Å) override; None = computed from lattice_const + burgers
        nedges=1,                           # number of edge dislocations: 1 (single) or 2 (dipole, opposite halves)
        fy_start=None,                      # fractional y of the removed half-plane; None = 0.5 (or 0.25/0.75 for nedges=2)
        add_vacuum=False,                   # add a vacuum slab after insertion (for free-surface boundary conditions)
        vacuum_direction=2,                 # axis index for the vacuum slab: 0=x, 1=y, 2=z
        lvac=20.0,                          # vacuum thickness (Å) if add_vacuum=True
        restore_axes=False,                 # True: swap axes back to input orientation; False: keep x=line, y=b, z=normal
        reset_ids=True,                     # renumber atom IDs contiguously afterwards
    )
    print(f"edge: {report}")                # dict: burgerm, newaxis used, atom counts before/after, axes_restored


def ex_edge_dislocation_100_cube(defect: BCCDefect):
    """<100>{010} edge dislocation (cube-plane glide, b = a<100>).

    [100](010) and [100](001) are the same dislocation (related by a 90°
    rotation about [100]) — pick either glide_plane. No oriented supercell
    is needed: the plain cubic-axes cell already has line/b/normal on box
    axes, and the permutation is derived automatically.
    """
    report = defect.edge_dislocation(
        ff_elements=["Fe", "Fe"],           # element symbols per atom type (asserted into force field)
        lattice_const=A0,                   # lattice constant a (Å) — auto-computes |b| = a for <100>
        burgers=(1, 0, 0),                  # Burgers vector DIRECTION (full a<100>, not a partial)
        glide_plane=(0, 1, 0),              # glide-plane normal; (0,0,1) is the equivalent variant (line then = [010])
        current_axes=[[1, 0, 0],            # crystallographic direction currently on x
                      [0, 1, 0],            # ... on y
                      [0, 0, 1]],           # ... on z — derives newaxis=[2,0,1]: x=[001] line, y=[100] b, z=[010] normal
        newaxis=None,                       # explicit swap_axes permutation override; None = auto
        burgerm=None,                       # explicit |b| (Å) override; None = a (one a/2-wide (200) plane removed)
        nedges=1,                           # 1 (single) or 2 (dipole)
        fy_start=None,                      # fractional y of the removed half-plane; None = 0.5
        add_vacuum=False,                   # add a vacuum slab after insertion
        vacuum_direction=2,                 # axis index for the vacuum slab
        lvac=20.0,                          # vacuum thickness (Å) if add_vacuum=True
        restore_axes=False,                 # keep x=line, y=b, z=normal afterwards
        reset_ids=True,                     # renumber atom IDs contiguously
    )
    print(f"<100>{{010}} edge: {report}")


def ex_edge_dislocation_100_110(defect: BCCDefect):
    """<100>{110} edge dislocation: b = a[100] on (011), line = [0-11].

    Note: (110) itself does NOT contain [100] (b·n = 1) — the {110} member
    that works with b=[100] is (011) or (0-11). This needs an oriented
    supercell built beforehand, e.g.::

        builder.trans_unit("bcc.POSCAR",
                           [[0, -1, 1], [1, 0, 0], [0, 1, 1]],
                           halve_c=False)   # halve_c is for screw cells only
        # then scale to a and replicate to size
    """
    report = defect.edge_dislocation(
        ff_elements=["Fe", "Fe"],           # element symbols per atom type
        lattice_const=A0,                   # lattice constant a (Å) — auto-computes |b| = a for <100>
        burgers=(1, 0, 0),                  # Burgers vector DIRECTION (full a<100>)
        glide_plane=(0, 1, 1),              # glide-plane normal — the {110} variant containing [100]
        current_axes=[[0, -1, 1],           # direction currently on x (dislocation line)
                      [1, 0, 0],            # ... on y (Burgers vector)
                      [0, 1, 1]],           # ... on z (glide-plane normal) — already ordered, newaxis=[0,1,2]
        newaxis=None,                       # explicit permutation override; None = auto
        burgerm=None,                       # explicit |b| (Å) override; None = a
        nedges=1,                           # 1 (single) or 2 (dipole)
        fy_start=None,                      # fractional y of the removed half-plane; None = 0.5
        add_vacuum=False,                   # add a vacuum slab after insertion
        vacuum_direction=2,                 # axis index for the vacuum slab
        lvac=20.0,                          # vacuum thickness (Å) if add_vacuum=True
        restore_axes=False,                 # keep x=line, y=b, z=normal afterwards
        reset_ids=True,                     # renumber atom IDs contiguously
    )
    print(f"<100>{{110}} edge: {report}")


def ex_screw_dislocation(defect: BCCDefect):
    """Screw dislocation (myLAMMPS backend; OVITO-based variant planned)."""
    report = defect.screw_dislocation(
        ff_elements=["Fe", "Fe"],           # element symbols per atom type
        lattice_const=A0,                   # lattice constant a (Å) — auto-computes |b| = a√3/2 for ½<111>
        burgers=(1, 1, 1),                  # Burgers direction (= line direction for a screw); pinned to z internally
        current_axes=[[1, -1, 0],           # direction currently on x   (bcc recommended frame:
                      [1, 1, -2],           # ... on y                    x=[1,-1,0], y=[1,1,-2], z=[1,1,1])
                      [1, 1, 1]],           # ... on z
        newaxis=None,                       # explicit permutation override; None = derived from current_axes
        burgerm=None,                       # explicit |b| (Å) override
        nscrews=2,                          # number of screw dislocations: 1, 2 (dipole), or 4 (quadrupole) — myLAMMPS-coded
        style="bcc",                        # crystal style: "bcc" or "fcc" (sets layer geometry + glide system)
        handle_pbc="tilt",                  # PBC fix-up for nscrews=1 in an orthogonal box: "tilt" the box or chop atoms
        orientation=True,                   # myLAMMPS internal orientation-handling flag (keep default True)
        add_vacuum=False,                   # add a vacuum slab after insertion
        vacuum_direction=0,                 # axis index for the vacuum slab
        lvac=20.0,                          # vacuum thickness (Å)
        restore_axes=False,                 # swap axes back to input orientation when done
    )
    print(f"screw: {report}")
    # Generic dispatcher: defect.dislocation(kind="edge", **kwargs) / defect.dislocation(kind="screw", **kwargs)


########################################################################################################################
# 5. COORDINATE TRANSFORMATION (helper, not a defect)
########################################################################################################################

def ex_transform(defect: BCCDefect):
    out = defect.transform(
        old_sys=np.eye(3),                  # 3x3 matrix, rows = CURRENT basis vectors
        new_sys=np.array([[1, 1, -2],
                          [1, 1, 1],
                          [1, -1, 0]]),     # 3x3 matrix, rows = TARGET basis vectors
        by="supercell",                     # "supercell": integer replication (returns a NEW lmpData); "transmat": rotate coords in place
        tol=1e-8,                           # tolerance for the integer check of the supercell matrix
    )
    return out


########################################################################################################################
# Legacy scratch — merge example (kept from the original notes)
########################################################################################################################
# data1 = Modlmp_LmpData.from_file("gao_ref111.data", "atomic")  # host
# data2 = Modlmp_LmpData.from_file("111SIA_111.data", "atomic")  # donor
# info = data1.merge_structure_seamless(
#     data2,
#     region_type="cube",
#     center_host=[54.78, 55.225, 58.55],
#     center_donor=[139.881, 97.7032, 80.7604],
#     side_lengths=(48.0, 60.0, 44.0),
#     pad_plus_donor=(0.0, 2.0, 0.0),
#     mapping="relative_cart",
#     donor_shift_cart=(0.0, -1.0, 0.0),
# )
# print(info)


if __name__ == "__main__":
    pot, a = "gao", 2.855312749000520  # Gao potential

    data = Modlmp_LmpData.from_file("gao_ref100.data", "atomic")
    defect = BCCDefect(data)

    # Uncomment the examples you want to run:
    # ex_vacancy(defect)
    # ex_interstitial(defect)
    # ex_random_interstitials(defect)
    # ex_dumbbell(defect)
    # ex_impurities(defect)
    # ex_void(defect)
    # ex_gas(defect)
    # ex_frenkel_pairs(defect)
    # ex_111_interstitial_loop(defect)
    # ex_111_vacancy_loop(defect)
    # ex_100_loop(defect)
    # ex_edge_dislocation(defect)   # needs an oriented supercell, see docstring
    # ex_edge_dislocation_100_cube(defect)  # <100>{010} edge — plain cubic-axes cell works
    # ex_edge_dislocation_100_110(defect)   # <100>{110} edge — needs [0-11]/[100]/[011] oriented cell
    # ex_screw_dislocation(defect)  # needs an oriented supercell
    # ex_transform(defect)

    data.to_file("test.data")
