"""Examples for EVERY public method in Modlmp_LmpData (gz_toolkit.core.modlmp).

Modlmp_LmpData subclasses myLAMMPS' ``lmpData`` and adds:
  * PBC / minimum-image helpers
  * region-based atom SELECTION (sphere, cube, cylinder, plane, type, ...)
  * atom EDITING (add / delete / retype / re-molecule)
  * box SCALING & straining (uniaxial / biaxial / triaxial, stress-driven)
  * crystallographic site identification

Each call lists every parameter on its own line, with a comment explaining it.
Replace file names / coordinates / constants with your own, then call the
example functions you want from ``__main__`` at the bottom.

The ``merge_structure_seamless`` method is intentionally NOT covered here.

Common setup used by all examples::

    data = Modlmp_LmpData.from_file("my_structure.data", "atomic")

Most editing methods return a count or a diagnostics dict, and accept
``reinit=`` to refresh internal bookkeeping (natoms, types, fractional
coords) after the change. Selection methods return a boolean mask aligned
row-for-row with ``data.atoms``.
"""

import numpy as np

from gz_toolkit import Modlmp_LmpData

# --- shared constants (bcc Fe single-type cell; tweak for your system) ------
DATA_FILE = "test.data"            # a LAMMPS data file (atom_style 'atomic')
FF_ELEMENTS = ["Fe", "Fe"]         # force-field element per atom type (type 1, 2)
ATOMIC_MASSES = [55.845, 55.845]   # masses (amu) matching FF_ELEMENTS

# A bcc conventional cell has 2 sites: corner (0,0,0) and body center (.5,.5,.5)
BCC_BASIS_FRAC = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
A0 = 2.8553                        # bcc Fe lattice constant (Å) -> unit-cell length


########################################################################################################################
# 0. LOAD / SAVE / COPY
########################################################################################################################

def ex_load():
    """Read a LAMMPS data file into a Modlmp_LmpData object (inherited ctor)."""
    data = Modlmp_LmpData.from_file(
        DATA_FILE,   # path to the LAMMPS data file
        "atomic",    # atom_style: "atomic" | "charge" | "molecular" | "full" ...
    )
    print("natoms:", data.natoms, "| box lengths:", np.round(data.box.lengths, 4))
    return data


def ex_copy(data: Modlmp_LmpData):
    """Deep-copy the whole object so edits don't touch the original."""
    clone = data.copy()             # independent deep copy
    clone.delete_atoms("all")       # mutate the clone only
    print("original natoms:", data.natoms, "| clone natoms:", clone.natoms)
    return clone


def ex_save(data: Modlmp_LmpData):
    """Write the structure back out to a LAMMPS data file (inherited)."""
    data.to_file(
        "out_structure.data",  # output path
        # to_atom_style=True,  # (default) keep only the columns required by atom_style
    )


########################################################################################################################
# 1. GEOMETRY / PBC HELPERS
########################################################################################################################

def ex_get_center(data: Modlmp_LmpData):
    """Geometric center of the system, in Cartesian or fractional coords."""
    cart = data.get_center(
        coord_type="cart",       # "cart" -> Cartesian Å ; "fraction" -> fractional
        zero_coords=False,       # if True, shift atoms to the box origin first
        thres=(0.1, 0.1, 0.1),   # threshold forwarded to find_center()
    )
    frac = data.get_center(coord_type="fraction")
    print("center (cart):", np.round(cart, 4), "| (frac):", np.round(frac, 4))
    return cart


def ex_wrap_cart(data: Modlmp_LmpData):
    """Wrap arbitrary Cartesian point(s) back into the periodic box (triclinic-safe)."""
    pts = np.array([
        [-1.0, -1.0, -1.0],          # a point just outside the lower corner
        data.box.lengths[0] + 5.0, 0, 0,  # a point past the +x face
    ], dtype=object)
    single = data.wrap_cart([-1.0, -1.0, -1.0])     # one (3,) point -> (3,)
    many = data.wrap_cart(np.array([[-1.0, -1.0, -1.0],
                                    [200.0, 0.0, 0.0]]))  # (N,3) -> (N,3)
    print("wrapped single:", np.round(single, 4))
    print("wrapped many:\n", np.round(many, 4))
    return single


def ex_min_image(data: Modlmp_LmpData):
    """Minimum-image relative vectors and distances under PBC."""
    ref = data.get_center("cart")                       # reference point
    coords = data.atoms[["x", "y", "z"]].to_numpy()[:5]  # first 5 atoms

    rel = data.min_image_rel(
        coords,   # (N,3) or (3,) positions
        ref,      # (N,3) or (3,) reference(s)
    )                                                   # -> folded displacement vectors
    dists = data.min_image_dists(
        coords,   # positions
        ref,      # reference(s)
    )                                                   # -> Euclidean MIC distances
    print("MIC distances of first 5 atoms to center:", np.round(dists, 4))
    return rel, dists


########################################################################################################################
# 2. SELECTION  (each returns a boolean mask aligned with data.atoms)
########################################################################################################################

def ex_select_type(data: Modlmp_LmpData):
    """Select by integer atom type and/or element symbol."""
    mask_int = data.select_type(1)               # type id 1
    mask_list = data.select_type([1, 2])         # several type ids
    mask_sym = data.select_type("Fe")            # by element symbol (needs masses/ff)
    mask_mixed = data.select_type([1, "Fe"])     # mixed ids + symbols
    print("type==1 count:", int(mask_int.sum()))
    return mask_int


def ex_select_coords(data: Modlmp_LmpData):
    """Match atoms sitting at exact coordinate(s), within a tolerance."""
    first_xyz = data.atoms[["x", "y", "z"]].to_numpy()[0]
    mask = data.select_coords(
        first_xyz,          # one (x,y,z) or a list/array of (x,y,z)
        tolerance=1e-3,     # match radius (Å)
    )
    print("matched atoms:", int(mask.sum()))
    return mask


def ex_select_sphere(data: Modlmp_LmpData):
    """Select atoms inside a sphere (Cartesian, no PBC)."""
    mask = data.select_sphere(
        center=data.get_center("cart"),  # (x,y,z) or list of centers
        radius=10.0,                     # radius in Å
    )
    print("atoms in sphere:", int(mask.sum()))
    return mask


def ex_select_cube(data: Modlmp_LmpData):
    """Axis-aligned box selection — three input forms shown."""
    c = data.get_center("cart")

    # form A: center + single side length (cube)
    mask_a = data.select_cube({
        "center": c,            # box center (x,y,z)
        "side_length": 20.0,    # isotropic edge length (Å)
    })

    # form B: center + per-axis lengths, with asymmetric padding
    mask_b = data.select_cube({
        "center": c,                       # box center
        "side_lengths": (20.0, 30.0, 10.0),  # (Lx, Ly, Lz)
        "pad_plus": (2.0, 0.0, 0.0),         # extend +x face by 2 Å
        "pad_minus": (0.0, 0.0, 5.0),        # extend -z face by 5 Å
    })

    # form C: explicit bounds
    mask_c = data.select_cube({
        "x_lo": 0.0, "x_hi": 30.0,
        "y_lo": 0.0, "y_hi": 30.0,
        "z_lo": 0.0, "z_hi": 30.0,
    })
    print("cube A/B/C counts:", int(mask_a.sum()), int(mask_b.sum()), int(mask_c.sum()))
    return mask_a


def ex_select_cylinder(data: Modlmp_LmpData):
    """Select atoms inside a cylinder; circle lies in a chosen plane."""
    c = data.get_center("cart")
    mask = data.select_cylinder(
        plane="xy",        # circle plane: 'xy' | 'xz' | 'yz' (height axis is the 3rd)
        center=c,          # (x,y,z) or list of centers
        radius=8.0,        # cylinder radius (Å)
        height=15.0,       # +/- thickness about center along the height axis (optional)
        # height_lower=0.0,  # OR give explicit bounds along the height axis
        # height_upper=40.0, #   (overrides `height` if provided)
    )
    print("atoms in cylinder:", int(mask.sum()))
    return mask


def ex_select_plane(data: Modlmp_LmpData):
    """Select a thin slab around a plane defined by a normal + point/level."""
    c = data.get_center("cart")
    mask_pt = data.select_plane(
        (0, 0, 1),         # plane normal (h,k,l) / any 3-vector
        point=c,           # a point the plane passes through
        tolerance=2.0,     # slab half-thickness (Å)
    )
    mask_lvl = data.select_plane(
        (1, 1, 0),         # normal
        level=50.0,        # OR constant in n·r = level (instead of point=)
        tolerance=1.5,
    )
    print("slab counts:", int(mask_pt.sum()), int(mask_lvl.sum()))
    return mask_pt


def ex_select_pct(data: Modlmp_LmpData):
    """Randomly select a fraction of atoms (optionally from a sub-region)."""
    base = data.select_sphere(data.get_center("cart"), 15.0)  # eligible pool
    mask = data.select_pct(
        0.10,                              # fraction in (0, 1]; here 10%
        base_mask=base,                    # restrict to this pool (optional)
        rng=np.random.default_rng(0),      # seed for reproducibility (optional)
    )
    print("randomly selected:", int(mask.sum()), "of", int(base.sum()))
    return mask


def ex_region_mask(data: Modlmp_LmpData):
    """Unified dispatcher — pick any region by name; same kwargs as select_*.

    This is what the editing methods (delete/retype) use internally, so the
    region keywords here are exactly what you pass to those methods.
    """
    c = data.get_center("cart")
    m_all = data.region_mask("all")
    m_sphere = data.region_mask("sphere", center=c, radius=10.0)
    m_cube = data.region_mask("cube", cube_params={"center": c, "side_length": 20.0})
    m_type = data.region_mask("type", types=[1])
    m_id = data.region_mask("id", atom_id=[1, 2, 3])     # by atom-ID (DataFrame index)
    # nested "pct of a base region":
    m_pct = data.region_mask(
        "pct",
        pct=0.2,                       # 20%
        base_region="sphere",          # of atoms in this base region
        base_kwargs={"center": c, "radius": 12.0},
        rng=np.random.default_rng(1),
    )
    print("region counts -> all:", int(m_all.sum()),
          "sphere:", int(m_sphere.sum()), "cube:", int(m_cube.sum()),
          "type:", int(m_type.sum()), "id:", int(m_id.sum()), "pct:", int(m_pct.sum()))
    return m_sphere


########################################################################################################################
# 3. EDITING  (mutates data; use data.copy() first if you want to keep the original)
########################################################################################################################

def ex_add_atoms(data: Modlmp_LmpData):
    """Append one or more atoms in Cartesian coordinates."""
    c = data.get_center("cart")
    new_ids = data.add_atoms(
        [c, c + np.array([1.0, 0.0, 0.0])],  # one (3,) or list/array of (3,) positions
        atom_type=1,                          # atom type id (keyword-only)
        molecule_id=None,                     # only used if atom_style has molecule-ID
        charge=None,                          # only used if atom_style has q
        wrap=True,                            # wrap inserted atoms into the box (PBC)
        reinit=True,                          # refresh bookkeeping afterwards
    )
    print("added atom IDs:", new_ids)
    return new_ids


def ex_delete_atoms(data: Modlmp_LmpData):
    """Delete atoms selected by any region_mask region + kwargs."""
    n = data.delete_atoms(
        "sphere",                        # region_type (see ex_region_mask)
        center=data.get_center("cart"),  # region kwargs ...
        radius=5.0,
        reinit=True,                     # refresh bookkeeping afterwards
    )
    print("deleted atoms:", n)
    return n


def ex_mod_atom_type(data: Modlmp_LmpData):
    """Reassign atom type within a region (and register the force field).

    Supports an optional second region: the final selection is
    target_region  ∩  source_region.
    """
    n = data.mod_atom_type(
        "sphere",                          # target region_type
        new_type=2,                        # new atom type id to assign
        ff_elements=FF_ELEMENTS,           # element per type (defines/extends Masses)
        atomic_masses=ATOMIC_MASSES,       # masses matching ff_elements
        center=data.get_center("cart"),    # target region kwargs ...
        radius=8.0,
        source_region_type="type",         # (optional) intersect with this region
        source_kwargs={"types": [1]},      # ... and its kwargs
        reinit=True,
    )
    print("retyped atoms:", n)
    return n


def ex_mod_molecule_id(data: Modlmp_LmpData):
    """Set molecule-ID within a region — requires atom_style='molecular'.

    Load with: Modlmp_LmpData.from_file("...", "molecular")
    """
    if getattr(data, "atom_style", None) != "molecular":
        print("skip: needs atom_style='molecular'")
        return 0
    n = data.mod_molecule_id(
        "cube",                                 # target region_type
        new_value=7,                            # molecule-ID to set
        cube_params={"center": data.get_center("cart"), "side_length": 20.0},
        source_region_type=None,                # optional intersect region
        source_kwargs=None,
        reinit=True,
    )
    print("re-moleculed atoms:", n)
    return n


def ex_wrap_atoms(data: Modlmp_LmpData):
    """Wrap ALL atoms back into the box via fractional coords (triclinic-safe).

    Handy after an affine transform that pushed atoms outside the box.
    """
    data.wrap_atoms(
        style=1,          # fractional normalization style (1 -> [0,1))
        drop_frac=True,   # drop temporary xsn/ysn/zsn columns afterwards
        reinit=False,     # set True to also refresh bookkeeping
    )
    print("wrapped all atoms back into the box")


########################################################################################################################
# 4. BOX SCALING & STRAINING  (scale_box) — 1, 2, or 3 axes
########################################################################################################################

def ex_scale_box_lc(data: Modlmp_LmpData):
    """mode='lc': set a NEW box length on chosen axis/axes."""
    target = data.box.lengths[0] * 1.05            # +5% on x
    info = data.scale_box(
        "x",                  # axis/axes: 'x'/'y'/'z', 0/1/2, or a list e.g. ['x','z']
        target,               # NEW length(s) in Å (scalar broadcasts to all driven axes)
        mode="lc",            # interpret `values` as target lengths
        affine=True,          # carry atoms with the box (preserve fractional coords)
        reinit=True,
    )
    print("lc x ->", np.round(info["new_lengths"], 4), "| vol ratio", round(info["volume_ratio"], 6))
    return info


def ex_scale_box_strain(data: Modlmp_LmpData):
    """mode='strain': apply normal strain(s); optionally conserve volume."""
    # uniaxial tensile strain, free lateral (volume grows)
    info_free = data.copy().scale_box(
        "z",                    # drive z only
        0.02,                   # engineering strain ds -> scale by (1+ds)
        mode="strain",
        conserve_volume=False,  # leave the other axes unchanged
    )

    # SAME strain but volume-conserving: x,y each scale by (1+ds)^(-1/2)
    info_iso = data.copy().scale_box(
        "z",
        0.02,
        mode="strain",
        conserve_volume=True,   # compensate the OTHER axes so volume stays constant
    )

    # biaxial, volume-conserving == Distortion.tetr_dis(ds): x,y up, z down
    info_tetr = data.copy().scale_box(
        ["x", "y"],             # drive two axes
        0.03,                   # same ds on both
        mode="strain",
        conserve_volume=True,   # third axis compensates: z *= 1/(1+ds)^2
    )

    # per-axis strains (orthorhombic): x:+ds, y:-ds, z compensates == orth_dis(ds)
    info_orth = data.copy().scale_box(
        ["x", "y"],
        [0.03, -0.03],          # one value per driven axis
        mode="strain",
        conserve_volume=True,
    )
    print("free vol ratio:", round(info_free["volume_ratio"], 6),
          "| conserve ratios:", round(info_iso["volume_ratio"], 6),
          round(info_tetr["volume_ratio"], 6), round(info_orth["volume_ratio"], 6))
    return info_iso


def ex_scale_box_stress(data: Modlmp_LmpData):
    """mode='stress': apply a stress; strain solved from ε = S·σ (S = C⁻¹).

    Needs elastic constants — cubic C11/C12/C44 OR a full 6×6 matrix.
    Lateral (Poisson) strains come out automatically; conserve_volume is
    ignored in this mode. Stress units must match the elastic constants.
    """
    # cubic Fe-like constants (GPa); apply uniaxial stress σ_xx = 5 GPa
    info_cubic = data.copy().scale_box(
        "x",                                       # axis carrying the applied stress
        5.0,                                       # stress value(s), same units as C (GPa)
        mode="stress",
        elastic_constants={"C11": 243.0, "C12": 145.0, "C44": 116.0},
        # shear_stress=(0.0, 0.0, 0.0),            # optional (σ_yz, σ_xz, σ_xy) -> tilt
    )

    # equivalent call using a full 6×6 stiffness matrix (Voigt: xx,yy,zz,yz,xz,xy)
    C = Modlmp_LmpData._build_stiffness_matrix({"C11": 243.0, "C12": 145.0, "C44": 116.0})
    info_full = data.copy().scale_box(
        "x",
        5.0,
        mode="stress",
        elastic_constants={"C": C},                # full 6×6 (hcp/tetragonal/etc.)
    )

    eps = np.array(info_cubic["strain_tensor"])
    print("strain (diag):", np.round(np.diag(eps), 6),
          "| axial+, lateral- (Poisson):", eps[0, 0] > 0 and eps[1, 1] < 0)
    print("full-C matches cubic dict:",
          np.allclose(info_cubic["strain_tensor"], info_full["strain_tensor"]))
    return info_cubic


########################################################################################################################
# 5. CRYSTALLOGRAPHIC SITE IDENTIFICATION
########################################################################################################################

def ex_identify_sites(data: Modlmp_LmpData):
    """Label each atom by its basis-site ID inside the repeating unit cell.

    The big (orthogonal) box must be an integer number of unit cells along
    each axis. Atoms that don't match any basis site get ID 0 (defect/distorted).

    NOTE: this requires a CUBE-oriented cell whose box lengths are integer
    multiples of the unit-cell lengths. An oriented cell (e.g. a [111]/[110]
    slab whose periodic lengths are irrational multiples of a0) will raise.
    """
    site_ids = data.identify_sites(
        BCC_BASIS_FRAC,              # (Nsite,3) fractional basis in ONE unit cell
        [A0, A0, A0],                # Cartesian unit-cell lengths [lx, ly, lz]
        tol=0.02,                    # tolerance in LOCAL fractional coords (strict=small)
        write_to_molecule_id=False,  # True -> store the ID in the 'molecule-ID' column
        reinit=False,                # True -> refresh bookkeeping after labeling
    )
    uniq, counts = np.unique(site_ids, return_counts=True)
    print("site-ID histogram:", dict(zip(uniq.tolist(), counts.tolist())))
    return site_ids


########################################################################################################################
# MAIN — uncomment the examples you want to run
########################################################################################################################

if __name__ == "__main__":
    data = ex_load()

    # --- geometry / PBC ---
    # ex_get_center(data)
    # ex_wrap_cart(data)
    # ex_min_image(data)

    # --- selection (non-mutating) ---
    # ex_select_type(data)
    # ex_select_coords(data)
    # ex_select_sphere(data)
    # ex_select_cube(data)
    # ex_select_cylinder(data)
    # ex_select_plane(data)
    # ex_select_pct(data)
    # ex_region_mask(data)

    # --- editing (mutating: work on a copy) ---
    # ex_add_atoms(data.copy())
    # ex_delete_atoms(data.copy())
    # ex_mod_atom_type(data.copy())
    # ex_mod_molecule_id(data.copy())   # needs atom_style='molecular'
    # ex_wrap_atoms(data.copy())

    # --- box scaling / straining ---
    # ex_scale_box_lc(data.copy())
    # ex_scale_box_strain(data)
    # ex_scale_box_stress(data)

    # --- site identification ---
    # ex_identify_sites(data)

    # --- save / copy ---
    # ex_copy(data)
    # ex_save(data.copy())
