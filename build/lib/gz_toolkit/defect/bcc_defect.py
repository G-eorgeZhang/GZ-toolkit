"""BCC defect creation tools.

Everything structure-agnostic (point defects, voids, dislocation
orientation, the prismatic-loop pipeline, Frenkel pairs, transforms) lives
in ``base.LatticeDefect``. This module supplies only the bcc
crystallography.

Method/class tree
-----------------
bcc_defect.py
`-- BCCDefect(LatticeDefect)
    |-- point defects                            [base]
    |   |-- vacancy(coords)
    |   |-- interstitial(sites) / interstitial_cluster(sites)
    |   |-- random_interstitials(...) / add_interstitial_bcc(...)
    |   |-- dumbbell(center, direction)
    |   |-- dumbbell_111 / dumbbell_100 / dumbbell_110(center)
    |   `-- impurities(region_type) / impurity_atom(site)
    |-- larger defects / sinks
    |   |-- void(center, radius)                 [base]
    |   |-- dislocation(kind, ...)               [base] -> edge/screw
    |   |-- edge_dislocation(...)                [base, bcc |b| + defaults]
    |   |-- screw_dislocation(...)               [base, style="bcc"]
    |   |-- add_loop(family=...)                      [base] validates family
    |   |-- add_111_loop(loop_type="sil"|"vl", ...)   # <111> loop (3 layers)
    |   |-- add_100_loop(loop_type="sil"|"vl", ...)   # <100> loop (2 layers)
    |   |-- precipitate(...) / grain_boundary(...)    # placeholders
    |   `-- gas(mode)
    `-- bcc crystallography (this module)
        |-- LOOP_FAMILIES                 <111> and <100>; no <110> (that is fcc)
        |-- _default_burgerm(...)         1/2<111> -> a*sqrt(3)/2; <100> -> a
        `-- _interstitial_offsets_frac()  tetrahedral + octahedral sites
"""

from __future__ import annotations

import warnings

import numpy as np

from gz_toolkit.core.modlmp import Modlmp_LmpData      # re-exported for callers
from gz_toolkit.defect.base import LatticeDefect


########################################################################################################################
class BCCDefect(LatticeDefect):
    """BCC defect builder using Modlmp_LmpData selection and edit methods."""

    lattice_name = "bcc"

    #: b = 1/2<111> on {110} is the dominant bcc slip system
    default_burgers = (1, 1, 1)
    default_glide_plane = (1, 1, 0)
    #: myLAMMPS create_screw_dislocation style
    screw_style = "bcc"

    # ------------------------------------------------------------------
    # bcc crystallography
    # ------------------------------------------------------------------

    @staticmethod
    def _default_burgerm(lattice_const: float, burgers) -> float:
        """|b| from the lattice constant for common bcc Burgers vectors."""
        hkl = sorted(abs(int(round(x))) for x in BCCDefect._parse_hkl(burgers))
        a = float(lattice_const)
        if hkl == [1, 1, 1]:
            return a * np.sqrt(3.0) / 2.0   # 1/2<111>
        if hkl == [0, 0, 1]:
            return a                         # <100>
        raise ValueError(
            f"No default |b| for burgers {burgers}; pass burgerm=... explicitly."
        )

    #: The two prismatic-loop families observed in bcc metals. <111> loops
    #: (b = 1/2<111>) dominate at low dose; <100> loops (b = a<100>) are the
    #: sessile ones that survive to higher dose in Fe. There is no <110>
    #: family here -- that is an fcc habit.
    LOOP_FAMILIES = {
        "111": {
            "allowed_abs": [1, 1, 1],
            "label": "111",
            "default_habit": (1, 1, 1),
            "burgers": lambda a: a * np.sqrt(3.0) / 2.0,   # 1/2<111>
            "k": 3,
            "d_layer": lambda a: a / (2.0 * np.sqrt(3.0)),
            "description": "b = 1/2<111>, glissile, 3 (111) layers",
        },
        "100": {
            "allowed_abs": [0, 0, 1],
            "label": "100",
            "default_habit": (1, 0, 0),
            "burgers": lambda a: a,                        # a<100>
            "k": 2,
            "d_layer": lambda a: a / 2.0,
            "description": "b = a<100>, sessile, 2 (200) layers",
        },
    }

    @staticmethod
    def _parse_111(habit_plane):
        return BCCDefect._parse_family(habit_plane, [1, 1, 1], "habit_plane", label="111")

    @staticmethod
    def _parse_100(habit_plane):
        return BCCDefect._parse_family(habit_plane, [0, 0, 1], "habit_plane", label="100")

    @staticmethod
    def _interstitial_offsets_frac() -> np.ndarray:
        """Tetrahedral + octahedral offsets from a bcc host atom, in units of a.

        Tetrahedral sites sit at 1/4<210>-type positions (12 per cell),
        octahedral at 1/2<100> face/edge centres (6 per cell).
        """
        base = [(0.25, 0.5, 0.0), (0.25, -0.5, 0.0), (-0.25, 0.5, 0.0), (-0.25, -0.5, 0.0)]
        perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
        tetra = np.unique(
            np.array([[v[p[0]], v[p[1]], v[p[2]]] for p in perms for v in base]), axis=0
        )
        octa = np.array([
            [0.5, 0.0, 0.0], [-0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0], [0.0, -0.5, 0.0],
            [0.0, 0.0, 0.5], [0.0, 0.0, -0.5],
        ])
        return np.vstack([tetra, octa])

    # ------------------------------------------------------------------
    # Legacy alias
    # ------------------------------------------------------------------

    def add_interstitial_bcc(self, **kwargs) -> int:
        """Deprecated alias for :meth:`random_interstitials`."""
        return self.random_interstitials(**kwargs)

    # ------------------------------------------------------------------
    # Prismatic dislocation loops (public entry points)
    #
    # A bcc loop is a coherent platelet of k atomic layers along the habit
    # normal. Both take loop_type = "sil" (self-interstitial) or "vl"
    # (vacancy, a separate code path). For oriented cells, pass box_axes
    # (crystal directions on box x/y/z) so the requested plane is rotated
    # into the box frame, or normal_box to give the box-frame normal.
    # ------------------------------------------------------------------

    def add_111_loop(self, *, loop_type: str = "sil", radius: float, lattice_const: float,
                     ff_elements: list[str], atomic_masses: list[float], center,
                     habit_plane=(1, 1, 1), box_axes=None, normal_box=None,
                     shape: str = "round", n_repeats: int = 1,
                     loop_atom_type: int | None = None, closure_layers: int = 2,
                     closure_disp_frac: float = 0.15, select_types=None,
                     wrap: bool = True, reinit: bool = True):
        """Add a <111> prismatic loop (3 (111) layers; b = 1/2<111>)."""
        return self._add_loop(
            family="111", loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit,
        )

    def add_100_loop(self, *, loop_type: str = "sil", radius: float, lattice_const: float,
                     ff_elements: list[str], atomic_masses: list[float], center,
                     habit_plane=(1, 0, 0), box_axes=None, normal_box=None,
                     shape: str = "round", n_repeats: int = 1,
                     loop_atom_type: int | None = None, closure_layers: int = 2,
                     closure_disp_frac: float = 0.15, select_types=None,
                     wrap: bool = True, reinit: bool = True):
        """Add a <100> prismatic loop (2 (200) layers; b = a<100>)."""
        return self._add_loop(
            family="100", loop_type=loop_type, radius=radius,
            lattice_const=lattice_const, ff_elements=ff_elements,
            atomic_masses=atomic_masses, center=center, habit_plane=habit_plane,
            box_axes=box_axes, normal_box=normal_box, shape=shape,
            n_repeats=n_repeats, loop_atom_type=loop_atom_type,
            closure_layers=closure_layers, closure_disp_frac=closure_disp_frac,
            select_types=select_types, wrap=wrap, reinit=reinit,
        )


BccDefect = BCCDefect
