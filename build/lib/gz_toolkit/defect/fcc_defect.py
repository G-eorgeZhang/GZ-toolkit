"""FCC defect creation tools.

Everything structure-agnostic lives in ``base.LatticeDefect``; this module
supplies the fcc crystallography plus the fcc-only defects.

Method/class tree
-----------------
fcc_defect.py
`-- FCCDefect(LatticeDefect)
    |-- point defects                            [base]
    |   |-- vacancy(coords)
    |   |-- interstitial(sites) / interstitial_cluster(sites)
    |   |-- random_interstitials(...)            octahedral + tetrahedral
    |   |-- dumbbell(center, direction)
    |   |-- dumbbell_100(center)                 <- the fcc SIA ground state
    |   `-- impurities(region_type) / impurity_atom(site)
    |-- larger defects / sinks
    |   |-- void(center, radius)                 [base]
    |   |-- dislocation(kind, ...)               [base] -> edge/screw
    |   |-- edge_dislocation(...)                b = 1/2<110> on {111}
    |   |-- screw_dislocation(...)               style="fcc"
    |   |-- add_loop(family=...)     generic, validates against LOOP_FAMILIES
    |   |-- add_frank_loop(...)      faulted, sessile, b = 1/3<111> on {111}
    |   |-- add_110_loop(...)        perfect, glissile, b = 1/2<110>
    |   |-- stacking_fault(...)      intrinsic SF on {111}
    |   `-- stacking_fault_tetrahedron(...)   1-layer <110>-edged {111}
    |                                         triangle; folds up on relax
    `-- fcc crystallography (this module)
        |-- LOOP_FAMILIES            <111> Frank + <110> perfect
        |-- _default_burgerm(...)    1/2<110>, 1/3<111>, 1/6<112>
        `-- _interstitial_offsets_frac()

Burgers vectors in fcc (a = cubic lattice constant)
---------------------------------------------------
    perfect          1/2<110>   |b| = a/sqrt(2)    glissile, slip on {111}
    Shockley partial 1/6<112>   |b| = a/sqrt(6)    bounds a stacking fault
    Frank partial    1/3<111>   |b| = a/sqrt(3)    sessile, bounds a faulted loop
    stair-rod        1/6<110>   |b| = a/(3*sqrt(2))  SFT edges

A perfect 1/2<110> dislocation dissociates into two Shockley partials
separated by an intrinsic stacking fault. The builders here insert the
*perfect* dislocation; dissociation emerges during relaxation, at a
separation set by the potential's stacking-fault energy.
"""

from __future__ import annotations

import numpy as np

from gz_toolkit.core.modlmp import Modlmp_LmpData      # re-exported for callers
from gz_toolkit.defect.base import LatticeDefect


########################################################################################################################
class FCCDefect(LatticeDefect):
    """FCC defect builder using Modlmp_LmpData selection and edit methods."""

    lattice_name = "fcc"

    #: b = 1/2<110> on {111} is the fcc slip system
    default_burgers = (1, 1, 0)
    default_glide_plane = (1, 1, 1)
    #: myLAMMPS create_screw_dislocation style (recommended frame:
    #: x=[1,1,1], y=[1,1,-2]/2, z=[1,-1,0]/2)
    screw_style = "fcc"

    # ------------------------------------------------------------------
    # fcc crystallography
    # ------------------------------------------------------------------

    @staticmethod
    def _default_burgerm(lattice_const: float, burgers) -> float:
        """|b| from the lattice constant for common fcc Burgers vectors.

        Keyed on the sorted |hkl| signature:
            <110> -> 1/2<110>, |b| = a/sqrt(2)      perfect
            <111> -> 1/3<111>, |b| = a/sqrt(3)      Frank partial
            <112> -> 1/6<112>, |b| = a/sqrt(6)      Shockley partial

        Note <110> is read as the *perfect* dislocation, not the 1/6<110>
        stair-rod -- stair-rods only arise inside an SFT and are never what
        a caller means by ``burgers=(1,1,0)``. Pass burgerm=... for those.
        """
        hkl = sorted(abs(int(round(x))) for x in FCCDefect._parse_hkl(burgers))
        a = float(lattice_const)
        if hkl == [0, 1, 1]:
            return a / np.sqrt(2.0)     # 1/2<110> perfect
        if hkl == [1, 1, 1]:
            return a / np.sqrt(3.0)     # 1/3<111> Frank partial
        if hkl == [1, 1, 2]:
            return a / np.sqrt(6.0)     # 1/6<112> Shockley partial
        raise ValueError(
            f"No default |b| for burgers {burgers}; pass burgerm=... explicitly."
        )

    #: The prismatic-loop families of fcc. Note these are NOT the bcc ones:
    #: fcc grows <110> perfect loops, which bcc does not, and its <111> loop
    #: is the *faulted Frank* loop (b = 1/3<111>, one close-packed layer),
    #: not the bcc 1/2<111> loop.
    #:
    #: "111" -- Frank loop. d_{111} = a/sqrt(3) and |b| = a/sqrt(3), so
    #:          k = 1: a Frank loop is a SINGLE {111} layer. Inserting or
    #:          removing one close-packed layer breaks the ABCABC sequence,
    #:          which is exactly why the loop is faulted and sessile.
    #: "110" -- perfect prismatic loop. The allowed fcc reflection is (220),
    #:          so d = a/(2 sqrt(2)); |b| = a/sqrt(2) gives k = 2. b is a
    #:          full lattice translation, so the loop is unfaulted and
    #:          glissile.
    LOOP_FAMILIES = {
        "111": {
            "allowed_abs": [1, 1, 1],
            "label": "111",
            "default_habit": (1, 1, 1),
            "burgers": lambda a: a / np.sqrt(3.0),         # 1/3<111>
            "k": 1,
            "d_layer": lambda a: a / np.sqrt(3.0),
            "description": "Frank, b = 1/3<111>, faulted sessile, 1 (111) layer",
        },
        "110": {
            "allowed_abs": [0, 1, 1],
            "label": "110",
            "default_habit": (1, 1, 0),
            "burgers": lambda a: a / np.sqrt(2.0),         # 1/2<110>
            "k": 2,
            "d_layer": lambda a: a / (2.0 * np.sqrt(2.0)),
            "description": "perfect, b = 1/2<110>, glissile, 2 (220) layers",
        },
    }

    @staticmethod
    def _parse_111(habit_plane):
        return FCCDefect._parse_family(habit_plane, [1, 1, 1], "habit_plane", label="111")

    @staticmethod
    def _parse_110(habit_plane):
        return FCCDefect._parse_family(habit_plane, [0, 1, 1], "habit_plane", label="110")

    @staticmethod
    def _interstitial_offsets_frac() -> np.ndarray:
        """Octahedral + tetrahedral offsets from an fcc host atom, in units of a.

        Octahedral sites are the 1/2<100> edge/body centres (6 neighbours of
        a lattice atom); tetrahedral sites are the 8 1/4<111> positions.
        """
        octa = np.array([
            [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0], [0.0, -0.5, 0.0],
            [0.0, 0.0, 0.5], [0.0, 0.0, -0.5],
        ])
        tetra = np.array([[sx * 0.25, sy * 0.25, sz * 0.25]
                          for sx in (1, -1) for sy in (1, -1) for sz in (1, -1)])
        return np.vstack([octa, tetra])

    # ------------------------------------------------------------------
    # Dislocation loops
    # ------------------------------------------------------------------

    def add_frank_loop(self, *, loop_type: str = "vl", radius: float, lattice_const: float,
                       ff_elements: list[str], atomic_masses: list[float], center,
                       habit_plane=(1, 1, 1), box_axes=None, normal_box=None,
                       shape: str = "round", n_repeats: int = 1,
                       triangle_rotation_deg: float = 0.0,
                       loop_atom_type: int | None = None, closure_layers: int = 2,
                       closure_disp_frac: float = 0.15, select_types=None,
                       wrap: bool = True, reinit: bool = True):
        """Add a Frank (faulted, sessile) loop on {111}; b = 1/3<111>.

        One {111} layer per repeat, so ``n_repeats`` counts inserted or
        removed close-packed layers directly.

        loop_type="vl" (default) removes a platelet -- the vacancy Frank loop
        that forms under irradiation and is the SFT precursor. "sil" inserts
        one, giving the interstitial Frank loop.

        shape="triangle" gives the triangular habit that real Frank loops
        adopt (edges along <110>); use ``triangle_rotation_deg`` to line the
        edges up with your cell's orientation.
        """
        return self._add_loop(
            family="111", loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit,
            triangle_rotation_deg=triangle_rotation_deg,
        )

    def add_110_loop(self, *, loop_type: str = "sil", radius: float, lattice_const: float,
                     ff_elements: list[str], atomic_masses: list[float], center,
                     habit_plane=(1, 1, 0), box_axes=None, normal_box=None,
                     shape: str = "round", n_repeats: int = 1,
                     loop_atom_type: int | None = None, closure_layers: int = 2,
                     closure_disp_frac: float = 0.15, select_types=None,
                     wrap: bool = True, reinit: bool = True):
        """Add a perfect (unfaulted, glissile) prismatic loop; b = 1/2<110>.

        Two {220} layers per repeat. Because b is a full lattice translation
        the stacking is preserved, so unlike the Frank loop this one carries
        no stacking fault and can glide on its prism.
        """
        return self._add_loop(
            family="110", loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit,
        )

    # ------------------------------------------------------------------
    # Stacking faults
    # ------------------------------------------------------------------

    def stacking_fault(self, *, lattice_const: float, plane_normal=(1, 1, 1),
                       fault_vector=None, plane_point=None, box_axes=None,
                       normal_box=None, wrap: bool = True, reinit: bool = True) -> dict:
        """Create a planar stacking fault by shearing one half of the crystal.

        Every atom on the positive side of the cut plane is displaced by
        ``fault_vector``; atoms below are left alone. The default fault
        vector is a Shockley partial 1/6<112> (|b| = a/sqrt(6)) lying in the
        habit plane, which turns ABCABC into an intrinsic stacking fault.

        Parameters
        ----------
        plane_normal : cubic {111} normal of the fault plane.
        fault_vector : (3,) Cartesian displacement, or None for the default
            Shockley partial. Must lie IN the fault plane -- an out-of-plane
            component would open or close the crystal rather than fault it,
            so it is rejected.
        plane_point : (3,) point on the cut plane (default: box centre).
        box_axes / normal_box : as for the loop builders, to map the cubic
            normal into the box frame of an oriented cell.

        Notes
        -----
        This is a rigid shear of a half-crystal, so it produces ONE fault in
        a non-periodic cell. Under 3-D PBC the displacement is not a lattice
        translation at the box boundary, and a second, unintended fault
        appears where the shear wraps. For a periodic cell either use a
        free surface along the normal, or shear a *slab* by calling this
        twice with opposite vectors.
        """
        data = self.data
        if data.atoms is None or len(data.atoms) == 0:
            return {"n_moved": 0}

        a = float(lattice_const)
        _hkl, n_cubic = self._parse_111(plane_normal)
        n = self._resolve_loop_normal(n_cubic, box_axes, normal_box)

        if plane_point is None:
            p0 = np.asarray(data.get_center("cart"), dtype=float).reshape(3,)
        else:
            p0 = np.asarray(plane_point, dtype=float).reshape(3,)

        if fault_vector is None:
            # a Shockley partial in the fault plane: any in-plane direction of
            # the right magnitude is crystallographically equivalent here.
            u, _v = self._plane_basis_from_normal(n)
            b = (a / np.sqrt(6.0)) * u
        else:
            b = np.asarray(fault_vector, dtype=float).reshape(3,)
            out_of_plane = float(np.dot(b, n))
            if abs(out_of_plane) > 1e-6 * max(1.0, float(np.linalg.norm(b))):
                raise ValueError(
                    f"fault_vector must lie in the fault plane, but its component "
                    f"along the normal is {out_of_plane:.4g}. A normal component "
                    f"opens/closes the crystal instead of faulting it."
                )

        pos = data.atoms[["x", "y", "z"]].to_numpy(dtype=float)
        above = ((pos - p0.reshape(1, 3)) @ n) > 0.0
        if not np.any(above):
            return {"n_moved": 0, "fault_vector": b.tolist(), "normal": n.tolist()}

        idx = data.atoms.index.to_numpy()[above]
        new = pos[above] + b.reshape(1, 3)
        if wrap:
            new = data.wrap_cart(new)
        data.atoms.loc[idx, ["x", "y", "z"]] = new

        if reinit:
            data.initialization(normalization=False, style=1)

        return {
            "n_moved": int(above.sum()),
            "fault_vector": b.tolist(),
            "burgers_magnitude": float(np.linalg.norm(b)),
            "normal": n.tolist(),
            "plane_point": p0.tolist(),
        }

    # ------------------------------------------------------------------
    # Stacking fault tetrahedron
    # ------------------------------------------------------------------

    def stacking_fault_tetrahedron(self, *, edge_length: float, lattice_const: float,
                                   ff_elements: list[str], atomic_masses: list[float],
                                   center, habit_plane=(1, 1, 1), box_axes=None,
                                   normal_box=None, triangle_rotation_deg: float = 0.0,
                                   select_types=None, wrap: bool = True,
                                   reinit: bool = True) -> dict:
        """Build a stacking fault tetrahedron by the Silcox-Hirsch route.

        Removes an **equilateral triangular platelet, exactly one {111}
        atomic layer thick**, with its edges along <110>. That is a
        triangular vacancy Frank loop (b = 1/3<111>).

        **The tetrahedron forms during relaxation, not here.** As written to
        disk this structure is a flat triangular hole. On minimisation /
        anneal the Frank partial bounding each of the three edges
        dissociates into a stair-rod (1/6<110>) plus a Shockley partial
        (1/6<112>); the three Shockleys glide up the inclined {111} planes
        and meet at the apex, closing the platelet into a tetrahedron with
        four faulted {111} faces and six stair-rod edges. Run a LAMMPS
        minimize (and usually a short anneal) before measuring anything.

        Single-layer selection
        ----------------------
        This is the part that has to be right: catching a second {111} plane
        gives a two-layer platelet, which collapses to a faulted loop rather
        than folding into a tetrahedron. The selection uses a half-width of
        0.30 * d_{111} = 0.30 * a/sqrt(3) (0.63 A at a = 3.615), snaps onto
        the actual atomic planes, and takes exactly one layer. The result is
        verified before any atom is deleted -- if the selected slab spans
        more than one interplanar spacing the call raises rather than
        silently producing the wrong defect.

        Parameters
        ----------
        edge_length : float
            Requested triangle edge length L in Angstrom. The platelet is
            quantised: it gains a whole close-packed row at a time, with rows
            spaced (a/sqrt(2))*sqrt(3)/2 apart, so what is removed is the
            largest lattice triangle fitting inside the request. Ask for
            L = 26 A in Cu and you get a 20.45 A triangle of 45 vacancies.
            The report returns ``edge_length_actual`` and ``n_vacancies`` --
            quote those, not the request, and not an area estimate.
        triangle_rotation_deg : float
            Spins the triangle in the habit plane. 0 puts the edges along
            the <110> direction chosen by ``_perpendicular_110``; the three
            <110> in-plane directions are equivalent by symmetry, so this is
            only needed to match a particular cell orientation.

        Returns
        -------
        dict
            n_vacancies, edge_length, layer_spacing, slab_spread,
            plane_tol, normal, and ``relaxed`` (always False -- a reminder
            that the fold-up has not happened yet).
        """
        radius = float(edge_length) / np.sqrt(3.0)

        ctx = self._select_loop_atoms(
            family="111", radius=radius, lattice_const=lattice_const,
            ff_elements=ff_elements, atomic_masses=atomic_masses, center=center,
            habit_plane=habit_plane, shape="triangle", n_repeats=1,
            select_types=select_types, box_axes=box_axes, normal_box=normal_box,
            triangle_rotation_deg=triangle_rotation_deg,
        )
        if ctx is None:
            return {
                "n_vacancies": 0,
                "edge_length_requested": float(edge_length),
                "edge_length_actual": 0.0,
                "layer_spacing": float(lattice_const) / np.sqrt(3.0),
                "slab_spread": 0.0, "plane_tol": None, "normal": None,
                "relaxed": False,
            }

        # Guard the single-layer requirement explicitly: a two-layer platelet
        # relaxes into a faulted loop, not an SFT, and the difference is not
        # obvious by eye afterwards.
        s = ctx["rel_all"][ctx["sel"]] @ ctx["n"]
        spread = float(s.max() - s.min())
        if spread > ctx["d_layer"]:
            raise RuntimeError(
                f"SFT platelet spans {spread:.4f} A along <111>, more than one "
                f"{{111}} interplanar spacing ({ctx['d_layer']:.4f} A). It would "
                f"relax into a faulted loop rather than a tetrahedron. Check that "
                f"the cell really is fcc with lattice_const={lattice_const}, and "
                f"that habit_plane/box_axes describe this cell's orientation."
            )

        # The platelet is quantised by close-packed rows, so the triangle that
        # actually gets removed is the largest lattice triangle fitting inside
        # the requested one -- often noticeably smaller. Report both; for an
        # equilateral triangle the circumradius is the max in-plane distance
        # from the centroid, and edge = R*sqrt(3).
        rel_sel = ctx["rel_all"][ctx["sel"]]
        in_plane = rel_sel - np.outer(rel_sel @ ctx["n"], ctx["n"])
        in_plane = in_plane - in_plane.mean(axis=0)
        edge_actual = float(np.linalg.norm(in_plane, axis=1).max() * np.sqrt(3.0))

        n = self._apply_vacancy_closure(
            ctx, burgers_vector=None,
            # No pre-collapse: the platelet must fold under relaxation. Nudging
            # the closure layers by hand here would bias the resulting geometry.
            closure_layers=0, closure_disp_frac=0.0,
            wrap_move_atoms=wrap, reinit=reinit,
        )

        return {
            "n_vacancies": int(n),
            "edge_length_requested": float(edge_length),
            "edge_length_actual": edge_actual,
            "layer_spacing": float(ctx["d_layer"]),
            "slab_spread": spread,
            "plane_tol": float(ctx["plane_tol"]),
            "normal": np.asarray(ctx["n"], dtype=float).tolist(),
            "relaxed": False,
        }


FccDefect = FCCDefect
