import copy
import warnings

import numpy as np
import pandas as pd
from pymatgen.core.lattice import Lattice
from mylammps.inputs.data import lmpBox, lmpData, lattice_2_lmpbox
from mylammps.elastic.distortion import Distortion


class Modlmp_LmpData(lmpData):

    # -----------------------
    # Vacuum padding
    # -----------------------

    def add_vacuum(self, lvac=20.0, direction=2, zero_coords=True, thres=[0.1, 0.1, 0.1]):
        """
        Pad the box with vacuum along `direction` WITHOUT moving any atoms.

        Overrides lmpData.add_vacuum (which re-centers atoms by shifting
        them by lvac/2 and stretches the box by only lvac total). Here both
        bounds are pushed out symmetrically by lvac:
            new_lo = old_lo - lvac
            new_hi = old_hi + lvac
        e.g. a 0-20 box along z becomes -20 to 40 for lvac=20. Atom Cartesian
        coordinates (x/y/z) are left exactly as they are; only the box and
        the derived fractional coordinates (xsn/ysn/zsn) are updated.

        `zero_coords`/`thres` are accepted only for call-site compatibility
        with lmpData.create_edge_dislocation/create_screw_dislocation (which
        call self.add_vacuum(..., zero_coords=True, ...)) and are not used.
        """
        bounds = copy.deepcopy(self.box.bounds)
        bounds[direction][0] -= lvac
        bounds[direction][1] += lvac
        self.box = lmpBox(bounds, tilt=self.box.tilt)
        self.coords2fracts(normalization=False)

    # -----------------------
    # PBC helpers (triclinic-safe)
    # -----------------------

    def wrap_cart(self, coords) -> np.ndarray:
        """
        Wrap Cartesian coordinate(s) back into the periodic box.

        coords: (3,) or (N, 3) Cartesian positions.
        Returns the same shape, wrapped via fractional coordinates [0, 1).
        """
        invT = np.asarray(self.box.inv_matrix, dtype=float).T
        matT = np.asarray(self.box.matrix, dtype=float).T
        P = np.asarray(coords, dtype=float)
        single = P.ndim == 1
        f = P.reshape(-1, 3) @ invT
        f -= np.floor(f)
        out = f @ matT
        return out[0] if single else out

    def min_image_rel(self, coords, ref) -> np.ndarray:
        """
        Minimum-image relative vector(s) coords - ref.

        coords: (3,) or (N, 3); ref: (3,) or (N, 3).
        Returns Cartesian displacement(s) folded into [-L/2, L/2).
        """
        invT = np.asarray(self.box.inv_matrix, dtype=float).T
        matT = np.asarray(self.box.matrix, dtype=float).T
        P = np.asarray(coords, dtype=float)
        Q = np.asarray(ref, dtype=float)
        single = P.ndim == 1 and Q.ndim == 1
        df = P.reshape(-1, 3) @ invT - Q.reshape(-1, 3) @ invT
        df -= np.round(df)
        out = df @ matT
        return out[0] if single else out

    def min_image_dists(self, coords, ref) -> np.ndarray:
        """Minimum-image Euclidean distance(s) between coords and ref."""
        rel = self.min_image_rel(coords, ref)
        return np.linalg.norm(rel, axis=-1)

    def add_atoms(self,coords,*,atom_type: int,molecule_id: int | None = None,charge: float | None = None,wrap: bool = False,reinit: bool = True,):
        """
        Append atom(s) to the Atoms section (self.atoms) in the current atom_style.

        coords: (3,) or list/array of (3,) in Cartesian coordinates.
        Returns: list of new atom IDs (DataFrame index values).
        """
        df = self.atoms
        if df is None or len(df.columns) == 0:
            raise ValueError("self.atoms is not initialized.")

        # require x/y/z columns (this matches your select_* and delete workflow)
        if not all(c in df.columns for c in ("x", "y", "z")):
            raise ValueError("Atoms table must contain Cartesian columns: x, y, z.")
        if "type" not in df.columns:
            raise ValueError("Atoms table has no 'type' column.")

        # normalize coords -> list of 3-vectors
        if (isinstance(coords, (list, tuple, np.ndarray))
                and len(coords) == 3
                and np.isscalar(coords[0])):
            targets = [coords]
        else:
            targets = list(coords)

        # start id
        if hasattr(self, "idmax") and self.idmax is not None:
            next_id = int(self.idmax) + 1
        else:
            next_id = int(df.index.to_numpy(dtype=int).max()) + 1 if len(df) else 1

        # Validate and optionally wrap all coordinates
        coords_arr = np.array([np.asarray(p, dtype=float).reshape(3,) for p in targets])
        if not np.isfinite(coords_arr).all():
            bad = coords_arr[~np.isfinite(coords_arr).all(axis=1)]
            raise ValueError(f"Bad coord (non-finite): {bad[0]}")

        if wrap:
            coords_arr = self.wrap_cart(coords_arr)

        # Build all new rows at once (column-wise, no per-row dicts)
        template = self.generate_default_dict()
        new_ids = list(range(next_id, next_id + len(coords_arr)))
        new_df = pd.DataFrame(template, index=new_ids)
        new_df["type"] = int(atom_type)
        new_df["x"] = coords_arr[:, 0]
        new_df["y"] = coords_arr[:, 1]
        new_df["z"] = coords_arr[:, 2]
        if "molecule-ID" in new_df.columns:
            new_df["molecule-ID"] = 0 if molecule_id is None else int(molecule_id)
        if "q" in new_df.columns:
            new_df["q"] = 0.0 if charge is None else float(charge)
        self.atoms = pd.concat([self.atoms, new_df], axis=0)

        # keep idmax consistent
        self.idmax = new_ids[-1] if new_ids else getattr(self, "idmax", 0)

        # optional: refresh internal bookkeeping once
        if reinit:
            self.initialization(normalization=False, style=1)

        return new_ids

    def delete_atoms(self, region_type: str, *, reinit: bool = True, **kwargs):
        """
        Delete atoms selected by region_mask(region_type, **kwargs),
        by directly filtering self.atoms (NOT by atom IDs).
        """
        if self.atoms is None or len(self.atoms) == 0:
            return 0

        mask = self.region_mask(region_type, **kwargs)
        if mask is None or len(mask) != len(self.atoms):
            raise ValueError("region_mask() must return a boolean mask aligned with self.atoms.")

        n_del = int(np.count_nonzero(mask))

        # delete by coords/region: filter rows
        self.atoms = self.atoms.loc[~mask].copy()

        if reinit:
            self.initialization(normalization=False, style=1)

        return n_del

    # ------------------------------------------------------------------
    # Screw dislocation — OVITO-style: atan2 displacement field + box shear.
    #
    # Overrides lmpData.create_screw_dislocation (the old bcc/fcc/hcp
    # glide-type + chop-atoms implementation). Convention: z || Burgers
    # vector == line direction (same frame BCCDefect.screw_dislocation
    # already rotates into via swap_axes before calling this).
    #
    # For each core at (cx, cy) in the local x-y cross-section:
    #     theta = atan2(y - cy, x - cx)
    #     z    += -sign * (burgerm / 2*pi) * theta
    # atan2's branch cut runs along -x from the core, so it always exits
    # the box through the x-periodic boundary. Rather than chopping atoms
    # there, the box is re-expressed with the a-vector carrying a burgerm/2
    # shear along z (so the periodic image across that boundary lines up
    # with the displacement jump); lattice_2_lmpbox + modify_by_symmetry
    # re-canonicalize that sheared matrix back into the restricted
    # LAMMPS triclinic form (a-vector along x) the same way swap_axes does,
    # rotating the atoms along with it so box and atoms stay consistent.
    # A dipole/quadrupole (nscrews=2/4) has canceling signed cores, so the
    # net shear is zero and no box change is needed at all.
    # ------------------------------------------------------------------

    def create_screw_dislocation(self, burgerm, nscrews=1, style="bcc", handle_pbc="tilt",
                                 orientation=True, add_vacuum=False, direction=0, lvac=20.0,
                                 center_offset=(0.0, 0.0)):
        """
        Insert nscrews screw dislocation(s) via an isotropic atan2
        displacement field plus a box shear (see class-level comment above).

        `style`, `handle_pbc`, and `orientation` are accepted only for
        signature compatibility with the call site in BCCDefect.screw_dislocation
        and are not used: this implementation is lattice-agnostic and always
        uses the box-shear (tilt) approach, never a per-style chop-atoms path.

        center_offset : (dx, dy) or list of (dx, dy), Cartesian Angstrom
            Default core placement snaps to the nearest atom to the box
            center (or the dipole/quadrupole fractional positions below).
            center_offset nudges the theta-reference point in the local
            x-y plane AFTER that snap, e.g. center_offset=(0.0, -2.0) moves
            the core 2 A down in y. A single (dx, dy) pair is applied to
            every core; pass a list of nscrews pairs to offset each core
            independently (e.g. for a dipole where the two cores need
            different fine-tuning).
        """
        self.zero_coords()
        buffer = burgerm * 0.03

        if nscrews == 1:
            centers_frac = [[0.5, 0.5, 0.5]]
            signs = [1.0]
        elif nscrews == 2:
            centers_frac = [[0.5, 0.25, 0.5], [0.5, 0.75, 0.5]]
            signs = [1.0, -1.0]
        elif nscrews == 4:
            centers_frac = [[0.25, 0.25, 0.5], [0.25, 0.75, 0.5],
                            [0.75, 0.25, 0.5], [0.75, 0.75, 0.5]]
            signs = [1.0, -1.0, -1.0, 1.0]
        else:
            raise ValueError("Uncoded number of screw dislocations.")

        if (len(center_offset) == 2
                and not isinstance(center_offset[0], (list, tuple, np.ndarray))):
            offsets = [center_offset] * len(centers_frac)
        else:
            offsets = list(center_offset)
            if len(offsets) != len(centers_frac):
                raise ValueError(
                    f"center_offset must be a single (dx, dy) pair or a list of "
                    f"{len(centers_frac)} pairs (one per core), got {len(offsets)}."
                )

        centers = []
        for cf, (ox, oy) in zip(centers_frac, offsets):
            coords = np.dot(np.array(cf), self.box.matrix)
            _, coords, _, _ = self.find_center_atom_coords(
                burgerm + buffer, center=coords, is_cartesian=True, style=1)
            coords = np.array(coords, dtype=float, copy=True)
            coords[0] += float(ox)
            coords[1] += float(oy)
            centers.append(coords)

        x = self.atoms["x"].to_numpy(dtype=float)
        y = self.atoms["y"].to_numpy(dtype=float)
        z = self.atoms["z"].to_numpy(dtype=float)

        dz = np.zeros(len(x))
        for (cx, cy, cz), s in zip(centers, signs):
            theta = np.arctan2(y - cy, x - cx)
            dz += -s * (burgerm / (2.0 * np.pi)) * theta

        self.atoms["z"] = z + dz

        net_sign = float(np.sum(signs))
        if net_sign != 0.0:
            newmatrix = copy.deepcopy(self.box.matrix)
            newmatrix[0, 2] += net_sign * burgerm / 2.0
            newlatt = Lattice(newmatrix)
            self.box, symmop = lattice_2_lmpbox(newlatt)
            self.atoms = lmpData.modify_by_symmetry(self.atoms, symmop)
        else:
            self.coords2fracts(normalization=False)

        if add_vacuum:
            self.add_vacuum(lvac=lvac, direction=direction, zero_coords=True)

########################################################################################################################

    def _require_xyz(self):
        df = self.atoms
        if df is None or len(df) == 0:
            return 0
        if not all(c in df.columns for c in ("x", "y", "z")):
            raise ValueError("Atoms table must contain Cartesian columns: 'x','y','z'.")
        return len(df)

    def _as_list_of_xyz(self, obj):
        """
        Normalize input into a list of (x,y,z) triplets.
        Accepts: single (3,), list of (3,), or (N,3) ndarray.
        """
        if obj is None:
            return []
        if isinstance(obj, np.ndarray) and obj.ndim == 2 and obj.shape[1] == 3:
            return [row for row in obj]
        if (isinstance(obj, (list, tuple, np.ndarray))
                and len(obj) == 3
                and isinstance(obj[0], (int, float, np.floating))):
            return [obj]
        return list(obj)

    def _empty_mask(self):
        n = 0 if self.atoms is None else len(self.atoms)
        return np.zeros(n, dtype=bool)

    def get_center(self,coord_type="cart", zero_coords=False, thres=(0.1, 0.1, 0.1)):
        """
        Return the geometric center of the system.

        Parameters
        ----------
        coord_type : str
            "cart"      → return Cartesian center
            "fraction"  → return fractional center
        zero_coords : bool
            Passed to find_center().
        thres : tuple
            Threshold passed to find_center().

        Returns
        -------
        np.ndarray of shape (3,)
        """
        result = self.find_center(zero_coords=zero_coords,thres=list(thres),)

        if not (isinstance(result, (tuple, list)) and len(result) == 2):
            raise ValueError("find_center() did not return (frac_center, cart_center) as expected.")

        frac_center, cart_center = result
        ct = str(coord_type).lower()

        if ct in ("fraction", "frac", "f"):
            return np.asarray(frac_center, dtype=float).reshape(3,)
        if ct in ("cart", "cartesian", "c"):
            return np.asarray(cart_center, dtype=float).reshape(3,)

        raise ValueError(f"coord_type must be 'fraction' or 'cart', got '{coord_type}'.")

    # -----------------------
    # select methods
    # -----------------------

    def select_coords(self, coords, tolerance: float = 1e-3) -> np.ndarray:
        """
        Match coordinates within tolerance.
        coords can be (x,y,z) or list/array of (x,y,z).
        Returns boolean mask aligned with self.atoms (row order).
        """
        n = self._require_xyz()
        if n == 0:
            return self._empty_mask()

        coords_list = self._as_list_of_xyz(coords)
        if not coords_list:
            return np.zeros(n, dtype=bool)

        xyz = self.atoms[["x", "y", "z"]].to_numpy(dtype=float)
        tol = float(tolerance)

        mask = np.zeros(n, dtype=bool)
        for c in coords_list:
            c = np.asarray(c, dtype=float).reshape(1, 3)
            mask |= np.isclose(xyz, c, atol=tol).all(axis=1)
        return mask

    def select_type(self, types) -> np.ndarray:
        """
        Select by atom type OR by element symbol(s) (via lmpData symbols mapping).

        Examples:
          select_type(1)
          select_type([1,2])
          select_type("Fe")
          select_type(["Fe","He"])
          select_type([1,"He"])   # mixed
        Returns boolean mask aligned with self.atoms.
        """
        df = self.atoms
        if df is None or len(df) == 0:
            return self._empty_mask()
        if "type" not in df.columns:
            raise ValueError("Atoms table has no 'type' column.")
        if types is None:
            return np.zeros(len(df), dtype=bool)

        items = [types] if isinstance(types, (int, np.integer, str)) else list(types)

        type_ids = []
        symbols = []
        for it in items:
            if isinstance(it, (int, np.integer)):
                type_ids.append(int(it))
            else:
                symbols.append(str(it).strip())

        t = df["type"].to_numpy(dtype=int)
        mask = np.zeros(len(df), dtype=bool)

        if type_ids:
            mask |= np.isin(t, np.array(type_ids, dtype=int))

        if symbols:
            # self.symbols is typically set by lmpData.get_data_info()
            if not hasattr(self, "symbols") or not self.symbols:
                self.get_data_info()

            symset = set(symbols)
            sym_type_ids = [i + 1 for i, sym in enumerate(self.symbols) if sym in symset]
            if sym_type_ids:
                mask |= np.isin(t, np.array(sym_type_ids, dtype=int))

        return mask

    def select_sphere(self, center, radius: float) -> np.ndarray:
        """
        Select atoms within a sphere (Cartesian).

        center: (x,y,z) in Cartesian coordinates
        radius: float
        Returns boolean mask aligned with self.atoms.
        """
        df = self.atoms
        if df is None or len(df) == 0:
            return np.zeros(0, dtype=bool)

        if center is None:
            return np.zeros(len(df), dtype=bool)

        if not all(c in df.columns for c in ("x", "y", "z")):
            raise ValueError("Atoms table must contain Cartesian columns: x, y, z")

        # allow single center or list of centers
        is_single = (
            isinstance(center, (list, tuple, np.ndarray))
            and len(center) == 3
            and isinstance(center[0], (int, float, np.floating))
        )
        centers = [center] if is_single else center

        pos = df[["x", "y", "z"]].to_numpy(dtype=float)
        r = float(radius)

        mask = np.zeros(len(df), dtype=bool)
        for c in centers:
            c = np.asarray(c, dtype=float).reshape(1, 3)
            d = np.linalg.norm(pos - c, axis=1)
            mask |= (d <= r)

        return mask

    def select_cylinder(self,plane: str,center,radius: float,height: float = None,height_lower: float = None,height_upper: float = None,):
        """
        Cylinder selection (Cartesian, no PBC/MIC).
        plane: 'xy', 'xz', 'yz' defines circle plane.
        center: (x,y,z) or list of centers.
        radius: radius in the circle plane.
        height: +/- thickness about center along height axis (optional)
        height_lower/height_upper: explicit bounds along height axis (optional, overrides height)
        Returns boolean mask aligned with self.atoms.
        """
        n = self._require_xyz()
        if n == 0:
            return self._empty_mask()

        centers = self._as_list_of_xyz(center)
        if not centers:
            return np.zeros(n, dtype=bool)

        plane = str(plane).lower()
        if plane == "xy":
            c_axes = ("x", "y"); h_axis = "z"
            pick2 = lambda c: (c[0], c[1]); pickh = lambda c: c[2]
        elif plane == "xz":
            c_axes = ("x", "z"); h_axis = "y"
            pick2 = lambda c: (c[0], c[2]); pickh = lambda c: c[1]
        elif plane == "yz":
            c_axes = ("y", "z"); h_axis = "x"
            pick2 = lambda c: (c[1], c[2]); pickh = lambda c: c[0]
        else:
            raise ValueError("plane must be one of: 'xy', 'xz', 'yz'.")

        r = float(radius)
        if r <= 0:
            return np.zeros(n, dtype=bool)
        r2 = r * r

        circ = self.atoms[list(c_axes)].to_numpy(dtype=float)  # (N,2)
        hvals = self.atoms[h_axis].to_numpy(dtype=float)       # (N,)

        mask_total = np.zeros(n, dtype=bool)
        for c in centers:
            c = np.asarray(c, dtype=float).reshape(3,)
            cc = np.array(pick2(c), dtype=float).reshape(1, 2)
            d2 = np.sum((circ - cc) ** 2, axis=1)
            m = d2 <= r2

            # height constraint
            if height_lower is not None or height_upper is not None:
                lo = -np.inf if height_lower is None else float(height_lower)
                hi =  np.inf if height_upper is None else float(height_upper)
                m &= (hvals >= lo) & (hvals <= hi)
            elif height is not None:
                hc = float(pickh(c))
                hh = float(height)
                m &= (np.abs(hvals - hc) <= hh)

            mask_total |= m

        return mask_total

    def select_plane(self,plane,*,point=None,level=None,tolerance: float = 1e-3,):
        """
        Select atoms within +/- tolerance of a plane (Cartesian).
        plane: normal vector (h,k,l) or any 3-vector.
        Provide either:
          - point=(x,y,z) on the plane, OR
          - level: constant in n·r = level
        Returns boolean mask aligned with self.atoms.
        """
        n = self._require_xyz()
        if n == 0:
            return self._empty_mask()

        nvec = np.asarray(plane, dtype=float).reshape(3,)
        norm = np.linalg.norm(nvec)
        if norm == 0:
            raise ValueError("Plane normal cannot be (0,0,0).")
        nvec /= norm

        if point is not None:
            level_val = float(np.dot(nvec, np.asarray(point, dtype=float).reshape(3,)))
        elif level is not None:
            level_val = float(level)
        else:
            raise ValueError("Provide either point=... or level=...")

        pos = self.atoms[["x", "y", "z"]].to_numpy(dtype=float)
        d = pos @ nvec
        return np.abs(d - level_val) <= float(tolerance)

    def select_cube(self, cube_params: dict = None) -> np.ndarray:
        """
        Axis-aligned box selection in Cartesian coordinates.

        center-based forms:
          - {"center": (cx,cy,cz), "side_length": L}
          - {"center": (cx,cy,cz), "side_lengths": (Lx,Ly,Lz)}
          - {"center": (cx,cy,cz), "lx":Lx, "ly":Ly, "lz":Lz}

        Optional asymmetric extension (Å):
          - "pad_plus":  (px,py,pz)  extends only + side  -> x_hi+=px, y_hi+=py, z_hi+=pz
          - "pad_minus": (mx,my,mz)  extends only - side  -> x_lo-=mx, y_lo-=my, z_lo-=mz

        bounds form:
          - {"x_lo":..,"x_hi":..,"y_lo":..,"y_hi":..,"z_lo":..,"z_hi":..}

        Returns boolean mask aligned with self.atoms.
        """


        n = self._require_xyz()
        if n == 0:
            return self._empty_mask()
        if not cube_params:
            return np.zeros(n, dtype=bool)

        pad_plus = np.asarray(cube_params.get("pad_plus", (0.0, 0.0, 0.0)), dtype=float).reshape(3, )
        pad_minus = np.asarray(cube_params.get("pad_minus", (0.0, 0.0, 0.0)), dtype=float).reshape(3, )

        if "center" in cube_params:
            cx, cy, cz = map(float, cube_params["center"])

            if "side_lengths" in cube_params:
                Lx, Ly, Lz = cube_params["side_lengths"]
                Lx, Ly, Lz = float(Lx), float(Ly), float(Lz)
            elif all(k in cube_params for k in ("lx", "ly", "lz")):
                Lx, Ly, Lz = float(cube_params["lx"]), float(cube_params["ly"]), float(cube_params["lz"])
            elif "side_length" in cube_params:
                L = float(cube_params["side_length"])
                Lx = Ly = Lz = L
            else:
                return np.zeros(n, dtype=bool)

            hx, hy, hz = Lx / 2.0, Ly / 2.0, Lz / 2.0
            x_lo, x_hi = cx - hx, cx + hx
            y_lo, y_hi = cy - hy, cy + hy
            z_lo, z_hi = cz - hz, cz + hz

            # apply directional padding (asymmetric)
            x_lo -= pad_minus[0];
            x_hi += pad_plus[0]
            y_lo -= pad_minus[1];
            y_hi += pad_plus[1]
            z_lo -= pad_minus[2];
            z_hi += pad_plus[2]

        elif all(k in cube_params for k in ("x_lo", "x_hi", "y_lo", "y_hi", "z_lo", "z_hi")):
            x_lo, x_hi = float(cube_params["x_lo"]), float(cube_params["x_hi"])
            y_lo, y_hi = float(cube_params["y_lo"]), float(cube_params["y_hi"])
            z_lo, z_hi = float(cube_params["z_lo"]), float(cube_params["z_hi"])

            # (optional) allow pad_* even for bounds form
            x_lo -= pad_minus[0];
            x_hi += pad_plus[0]
            y_lo -= pad_minus[1];
            y_hi += pad_plus[1]
            z_lo -= pad_minus[2];
            z_hi += pad_plus[2]

        else:
            return np.zeros(n, dtype=bool)

        x = self.atoms["x"].to_numpy(dtype=float)
        y = self.atoms["y"].to_numpy(dtype=float)
        z = self.atoms["z"].to_numpy(dtype=float)
        return (x >= x_lo) & (x <= x_hi) & (y >= y_lo) & (y <= y_hi) & (z >= z_lo) & (z <= z_hi)

    def select_pct(self,pct: float,*,base_mask: np.ndarray | None = None,rng: np.random.Generator | None = None,) -> np.ndarray:
        """
        Select a percentage of atoms randomly (optionally from a base_mask).
        Returns boolean mask aligned with self.atoms.
        """

        if self.atoms is None or len(self.atoms) == 0:
            return np.zeros(0, dtype=bool)

        if not (0 < pct <= 1):
            raise ValueError("pct must be in (0,1].")

        n = len(self.atoms)

        if base_mask is None:
            base_mask = np.ones(n, dtype=bool)
        else:
            base_mask = np.asarray(base_mask, dtype=bool)
            if len(base_mask) != n:
                raise ValueError("base_mask length mismatch.")

        eligible_idx = np.where(base_mask)[0]

        if len(eligible_idx) == 0:
            return np.zeros(n, dtype=bool)

        n_select = int(np.floor(pct * len(eligible_idx)))
        n_select = max(1, n_select)

        if rng is None:
            rng = np.random.default_rng()

        chosen = rng.choice(eligible_idx, size=n_select, replace=False)

        mask = np.zeros(n, dtype=bool)
        mask[chosen] = True

        return mask

    def region_mask(self, region_type: str, **kwargs):
        """
        Unified dispatcher: returns boolean mask aligned with self.atoms.
        """
        rt = str(region_type).lower()
        n = 0 if self.atoms is None else len(self.atoms)

        if rt in ("all",):
            return np.ones(n, dtype=bool)

        if rt in ("coordinates", "coords"):
            return self.select_coords(
                coords=kwargs.get("coords"),
                tolerance=kwargs.get("tolerance", 1e-3),
            )

        if rt in ("type", "atom_type"):
            return self.select_type(kwargs.get("types"))

        if rt in ("sphere",):
            return self.select_sphere(
                center=kwargs.get("center"),
                radius=kwargs.get("radius"),
            )

        if rt in ("cylinder",):
            return self.select_cylinder(
                plane=kwargs.get("plane"),
                center=kwargs.get("center"),
                radius=kwargs.get("radius"),
                height=kwargs.get("height"),
                height_lower=kwargs.get("height_lower"),
                height_upper=kwargs.get("height_upper"),
            )

        if rt in ("plane",):
            return self.select_plane(
                plane=kwargs.get("plane"),
                point=kwargs.get("point", None),
                level=kwargs.get("level", None),
                tolerance=kwargs.get("tolerance", 1e-3),
            )

        if rt in ("cube", "box"):
            return self.select_cube(kwargs.get("cube_params", None))

        if rt in ("atom_id", "id"):
            atom_id = kwargs.get("atom_id")
            if atom_id is None:
                return np.zeros(n, dtype=bool)
            ids = [int(atom_id)] if isinstance(atom_id, (int, np.integer)) else [int(i) for i in atom_id]
            return np.isin(self.atoms.index.to_numpy(dtype=int), ids)

        if rt in ("mol_id", "molecule_id", "molecule-id", "molecule"):
            if self.atoms is None or "molecule-ID" not in self.atoms.columns:
                return np.zeros(n, dtype=bool)
            mol_id = kwargs.get("mol_id")
            if mol_id is None:
                return np.zeros(n, dtype=bool)
            mids = [int(mol_id)] if isinstance(mol_id, (int, np.integer)) else [int(i) for i in mol_id]
            return self.atoms["molecule-ID"].astype(int).isin(mids).to_numpy(dtype=bool)

        if rt in ("pct", "percentage"):
            pct = kwargs.get("pct", None)
            if pct is None:
                raise ValueError("pct region requires pct=...")

            # optional base region
            base_region = kwargs.get("base_region", "all")
            base_kwargs = kwargs.get("base_kwargs", {})

            base_mask = self.region_mask(base_region, **base_kwargs)

            return self.select_pct(
                pct=pct,
                base_mask=base_mask,
                rng=kwargs.get("rng", None),
            )

        raise ValueError(f"Unknown region_type '{region_type}'.")

    def _resolve_mask(self,region_type: str,source_region_type: str | None = None,source_kwargs: dict | None = None,**kwargs,) -> np.ndarray:
        """
        Resolve the final selection mask:
            final_mask = target_region ∩ source_region (if source given)
        Returns a boolean mask aligned with self.atoms.
        """
        target_mask = self.region_mask(region_type, **kwargs)
        if target_mask is None or len(target_mask) != len(self.atoms):
            raise ValueError("region_mask() must return a boolean mask aligned with self.atoms.")

        if source_region_type is None:
            return target_mask

        source_kwargs = {} if source_kwargs is None else dict(source_kwargs)
        source_mask = self.region_mask(source_region_type, **source_kwargs)
        if source_mask is None or len(source_mask) != len(self.atoms):
            raise ValueError("source region_mask() invalid.")

        return target_mask & source_mask

########################################################################################################################

    def mod_atom_type(self,region_type: str,*,new_type: int,ff_elements: list[str],atomic_masses: list[float],source_region_type: str | None = None,
    source_kwargs: dict | None = None,reinit: bool = True,**kwargs,):
        """
            Modify atom type for atoms selected by region_mask(region_type, **kwargs).

            Supports two-region logic:
                final_mask = target_region ∩ source_region

            Note:
                pct selection should be handled via region_mask("pct", ...)
            """

        # ------------------------------------------------------------
        # 1) Force field
        # ------------------------------------------------------------
        if ff_elements is None or atomic_masses is None:
            raise ValueError("You must provide ff_elements=... and atomic_masses=...")

        if len(ff_elements) != len(atomic_masses):
            raise ValueError("ff_elements and atomic_masses must have the same length.")

        self.assert_force_field(ff_elements, atomic_masses=atomic_masses)

        # ------------------------------------------------------------
        # 2) Basic checks
        # ------------------------------------------------------------
        if self.atoms is None or len(self.atoms) == 0:
            return 0
        if "type" not in self.atoms.columns:
            raise ValueError("Atoms table has no 'type' column.")

        valid_types = set(int(i) for i in self.masses.index.to_numpy(dtype=int))
        if int(new_type) not in valid_types:
            raise ValueError(
                f"new_type={new_type} not present in Masses after assert_force_field. "
                f"Valid types are: {sorted(valid_types)}"
            )

        # ------------------------------------------------------------
        # 3) Resolve target ∩ source mask and apply
        # ------------------------------------------------------------
        mask = self._resolve_mask(region_type, source_region_type, source_kwargs, **kwargs)

        n_mod = int(np.count_nonzero(mask))
        if n_mod == 0:
            return 0

        self.atoms.loc[mask, "type"] = int(new_type)

        if reinit:
            self.initialization(normalization=False, style=1)

        return n_mod

    def mod_molecule_id(self,region_type: str,*,new_value: int,source_region_type: str | None = None,
    source_kwargs: dict | None = None,reinit: bool = True,**kwargs):
        """
            Set 'molecule-ID' for atoms selected by region_mask(region_type, **kwargs).

            Supports:
                final_mask = target_region ∩ source_region
            """

        if getattr(self, "atom_style", None) != "molecular":
            raise ValueError(
                f"modify_molecule_id_region requires atom_style='molecular', "
                f"but got {getattr(self, 'atom_style', None)!r}."
            )

        if self.atoms is None or len(self.atoms) == 0:
            return 0
        if "molecule-ID" not in self.atoms.columns:
            raise ValueError("Atoms table has no 'molecule-ID' column.")

        mask = self._resolve_mask(region_type, source_region_type, source_kwargs, **kwargs)

        n_mod = int(np.count_nonzero(mask))
        if n_mod == 0:
            return 0

        self.atoms.loc[mask, "molecule-ID"] = int(new_value)

        if reinit:
            self.initialization(normalization=False, style=1)

        return n_mod

########################################################################################################################

    def wrap_atoms(self, *, style: int = 1, drop_frac: bool = True, reinit: bool = False):
        """
        Wrap all atom coordinates back into the current periodic simulation box.

        Triclinic-safe: Cartesian -> fractional (via box.inv_matrix),
        normalize fractional coords into the box, then fractional -> Cartesian.

        Parameters
        ----------
        style : int
            Normalization style forwarded to lmpData.normalize_frac_coords.
            style=1 gives [0,1). Other styles follow your existing implementation.
        drop_frac : bool
            If True, remove temporary fractional columns ('xsn','ysn','zsn') after wrapping.
        reinit : bool
            If True, call self.initialization(...) after wrapping.

        Notes
        -----
        - Does NOT change the box. Only remaps atoms into the existing box with PBC.
        - Intended to be called after transforms like modify_by_symmetry().
        """
        if self.atoms is None or len(self.atoms) == 0:
            return

        if not all(c in self.atoms.columns for c in ("x", "y", "z")):
            raise ValueError("Atoms table must contain Cartesian columns: x, y, z.")

        # This exists in the parent lmpData (data.py):
        # coords2fracts(normalization=True) + fracts2coords()
        self.normalize_coords(style=style)

        if drop_frac:
            dropcols = [c for c in ("xsn", "ysn", "zsn") if c in self.atoms.columns]
            if dropcols:
                self.atoms = self.atoms.drop(columns=dropcols)

        if reinit:
            self.initialization(normalization=False, style=1)


########################################################################################################################
#  scale_box — per-axis box scaling / straining (uniaxial / biaxial / triaxial)
########################################################################################################################

    @staticmethod
    def _build_stiffness_matrix(elastic_constants: dict) -> np.ndarray:
        """
        Build a 6x6 stiffness matrix C (Voigt order [xx, yy, zz, yz, xz, xy]).

        Accepts either:
          - {"C": <(6, 6) array-like>}                          full matrix
          - {"C11": .., "C12": .., "C44": ..}                   cubic
          - {"C11": .., "C12": .., "C13": .., "C33": .., "C44": ..}
                                                                hexagonal (hcp)

        Units are whatever the caller uses; they must match the applied
        stress so that strain (= S . sigma) comes out dimensionless.

        Notes
        -----
        The hexagonal form is *transversely isotropic about z*, so it is only
        correct when the crystal c-axis lies along the box z-axis. That is the
        case for cells built from the ``A3_ORTHO`` prototype (x=[2-1-10],
        y=[01-10], z=[0001]); if you have re-oriented the cell, pass an
        explicit 6x6 ``C`` instead.

        Distinguishing hexagonal from cubic is done on the presence of C13/C33,
        because the cubic key set {C11, C12, C44} is a strict subset of the
        hexagonal one -- an hcp dict would otherwise be silently accepted as
        cubic and give wrong lateral (Poisson) strains.
        """
        ec = dict(elastic_constants)

        if "C" in ec:
            C = np.asarray(ec["C"], dtype=float)
            if C.shape != (6, 6):
                raise ValueError("elastic_constants['C'] must be a 6x6 matrix.")
            return C

        hexagonal = ("C11", "C12", "C13", "C33", "C44")
        if all(k in ec for k in hexagonal):
            c11 = float(ec["C11"])
            c12 = float(ec["C12"])
            c13 = float(ec["C13"])
            c33 = float(ec["C33"])
            c44 = float(ec["C44"])
            # C66 is not independent for a hexagonal crystal, but honour an
            # explicitly supplied value (e.g. a fitted, slightly off-symmetry
            # tensor) rather than silently overwriting it.
            c66 = float(ec.get("C66", 0.5 * (c11 - c12)))
            C = np.zeros((6, 6), dtype=float)
            C[0, 0] = C[1, 1] = c11
            C[2, 2] = c33
            C[0, 1] = C[1, 0] = c12
            C[0, 2] = C[2, 0] = c13
            C[1, 2] = C[2, 1] = c13
            C[3, 3] = C[4, 4] = c44
            C[5, 5] = c66
            return C

        cubic = ("C11", "C12", "C44")
        if all(k in ec for k in cubic):
            c11 = float(ec["C11"])
            c12 = float(ec["C12"])
            c44 = float(ec["C44"])
            C = np.zeros((6, 6), dtype=float)
            # normal block
            for i in range(3):
                for j in range(3):
                    C[i, j] = c11 if i == j else c12
            # shear block
            for i in range(3, 6):
                C[i, i] = c44
            return C

        raise ValueError(
            "elastic_constants must contain either a full 6x6 'C', the cubic "
            "constants 'C11', 'C12', 'C44', or the hexagonal constants "
            "'C11', 'C12', 'C13', 'C33', 'C44'. "
            f"Got: {sorted(ec)}"
        )

    def scale_box(
        self,
        axes,
        values,
        *,
        mode: str = "strain",
        conserve_volume: bool = False,
        elastic_constants: dict | None = None,
        shear_stress=(0.0, 0.0, 0.0),
        affine: bool = True,
        reinit: bool = True,
    ) -> dict:
        """
        Scale / strain the simulation box along one, two, or three axes.

        This is the triclinic-safe, per-axis generalization of
        ``lmpData.scale_data`` (which can only scale x, y, z simultaneously).
        It is intended for applying stress/strain to a piece of material.

        Parameters
        ----------
        axes : str | int | sequence
            Axis or axes to drive. Accepts 'x'/'y'/'z' (case-insensitive),
            integers 0/1/2, or a list/tuple mixing them, e.g. ['x', 'z'].
        values : float | sequence
            One value per driven axis (a scalar is broadcast to all driven
            axes). Meaning depends on ``mode``:
              - 'lc'     : the NEW box length along that axis (Angstrom).
              - 'strain' : the engineering normal strain ds (dimensionless);
                           the axis is scaled by (1 + ds).
              - 'stress' : the applied normal stress along that axis, in the
                           same units as the elastic constants (e.g. GPa).
        mode : {'lc', 'strain', 'stress'}
            How ``values`` is interpreted (see above).
        conserve_volume : bool
            Only used for 'lc' and 'strain'. If True, the axes that are NOT
            driven are scaled so the box volume is preserved, using the same
            construction as ``elastic.distortion.Distortion.tetr_dis`` /
            ``orth_dis`` (where ``ds`` is the strain value): the product of all
            three scale factors is forced to 1, with the compensating strain
            split equally (geometrically) over the remaining axes::

                scale 1 axis by f      -> other two each scale by f ** (-1/2)
                scale 2 axes by f1, f2 -> third axis scales by 1 / (f1 * f2)

            (These reproduce tetr_dis and orth_dis as special cases.)
            Ignored in 'stress' mode and when all three axes are driven.
        elastic_constants : dict, required for mode='stress'
            One of:
              - cubic       ``{"C11": .., "C12": .., "C44": ..}``
              - hexagonal   ``{"C11": .., "C12": .., "C13": .., "C33": ..,
                               "C44": ..}`` -- for hcp; assumes the c-axis is
                            along box z (true for the ``A3_ORTHO`` prototype)
              - explicit    ``{"C": <(6, 6)>}`` full stiffness matrix
            Voigt order [xx, yy, zz, yz, xz, xy]. Units must match the applied
            stress. See ``_build_stiffness_matrix``.
        shear_stress : sequence of 3 floats
            Optional applied shear stresses (sigma_yz, sigma_xz, sigma_xy) for
            mode='stress'. Default (0, 0, 0). Non-zero values produce box tilt.
        affine : bool
            If True (default) atoms are carried with the box (fractional
            coordinates preserved). If False, the box changes but Cartesian
            atom positions are left untouched.
        reinit : bool
            Refresh internal bookkeeping after the transform.

        Returns
        -------
        dict
            Diagnostics: driven axes, per-axis stretch factors, strain tensor,
            deformation matrix, old/new box lengths and volume.

        Notes
        -----
        In 'stress' mode the full strain response is solved from the inverse
        stiffness, eps = S . sigma with S = C^-1 (Voigt). This includes the
        lateral Poisson strains, so the volume generally changes and
        ``conserve_volume`` does not apply.
        """
        if self.atoms is None or len(self.atoms) == 0:
            raise ValueError("No atoms to scale.")
        if not all(c in self.atoms.columns for c in ("x", "y", "z")):
            raise ValueError("Atoms table must contain Cartesian columns: x, y, z.")

        mode = str(mode).lower().strip()
        if mode not in ("lc", "strain", "stress"):
            raise ValueError("mode must be 'lc', 'strain', or 'stress'.")

        # --- normalize axes to a sorted list of unique indices in {0, 1, 2} ---
        axis_map = {"x": 0, "y": 1, "z": 2}
        axes_in = [axes] if isinstance(axes, (str, int, np.integer)) else list(axes)
        driven = set()
        for a in axes_in:
            if isinstance(a, str):
                key = a.strip().lower()
                if key not in axis_map:
                    raise ValueError(f"Invalid axis '{a}'; use 'x'/'y'/'z' or 0/1/2.")
                driven.add(axis_map[key])
            else:
                idx = int(a)
                if idx not in (0, 1, 2):
                    raise ValueError(f"Invalid axis {a}; use 'x'/'y'/'z' or 0/1/2.")
                driven.add(idx)
        driven = sorted(driven)
        if not driven:
            raise ValueError("No axes specified.")

        orgbox = copy.deepcopy(np.asarray(self.box.matrix, dtype=float))
        old_lengths = np.asarray(self.box.lengths, dtype=float).reshape(3,)
        old_volume = float(abs(np.linalg.det(orgbox)))

        # ==============================================================
        # Build the 3x3 deformation matrix D (new_vec_i = D @ old_vec_i)
        # ==============================================================
        if mode in ("lc", "strain"):
            vals = np.atleast_1d(np.asarray(values, dtype=float)).reshape(-1)
            if vals.size == 1:
                vals = np.full(len(driven), float(vals[0]))
            if vals.size != len(driven):
                raise ValueError(
                    f"Got {vals.size} value(s) for {len(driven)} driven axis(es)."
                )

            factors = np.ones(3, dtype=float)
            for ax, v in zip(driven, vals):
                if mode == "lc":
                    if old_lengths[ax] <= 0:
                        raise ValueError(f"Non-positive box length on axis {ax}.")
                    factors[ax] = float(v) / old_lengths[ax]
                else:  # strain
                    factors[ax] = 1.0 + float(v)

            if np.any(factors <= 0):
                raise ValueError("Resulting scale factor(s) must be positive.")

            comp = [a for a in (0, 1, 2) if a not in driven]
            if conserve_volume:
                if not comp:
                    warnings.warn(
                        "conserve_volume ignored: all three axes are driven.",
                        stacklevel=2,
                    )
                else:
                    # ProdAll factors == 1  ->  compensating axes split equally
                    # (geometric). Reproduces tetr_dis / orth_dis.
                    driven_prod = float(np.prod([factors[a] for a in driven]))
                    g = driven_prod ** (-1.0 / len(comp))
                    for a in comp:
                        factors[a] = g

            D = np.diag(factors)
            strain_tensor = D - np.eye(3)

        else:  # stress
            if elastic_constants is None:
                raise ValueError("mode='stress' requires elastic_constants=...")
            if conserve_volume:
                warnings.warn(
                    "conserve_volume is ignored in stress mode; the lateral "
                    "response is the physical Poisson strain from S = C^-1.",
                    stacklevel=2,
                )

            C = self._build_stiffness_matrix(elastic_constants)
            S = np.linalg.inv(C)

            svals = np.atleast_1d(np.asarray(values, dtype=float)).reshape(-1)
            if svals.size == 1:
                svals = np.full(len(driven), float(svals[0]))
            if svals.size != len(driven):
                raise ValueError(
                    f"Got {svals.size} stress value(s) for {len(driven)} driven axis(es)."
                )

            # Applied stress as a Voigt vector [xx, yy, zz, yz, xz, xy]
            sigma = np.zeros(6, dtype=float)
            for ax, v in zip(driven, svals):
                sigma[ax] = float(v)
            sigma[3:6] = np.asarray(shear_stress, dtype=float).reshape(3,)

            # Voigt strain (engineering shear): eps = S . sigma
            eps_v = S @ sigma
            # symmetric strain tensor (engineering shear -> tensor component / 2)
            strain_tensor = np.array([
                [eps_v[0],      eps_v[5] / 2.0, eps_v[4] / 2.0],
                [eps_v[5] / 2.0, eps_v[1],      eps_v[3] / 2.0],
                [eps_v[4] / 2.0, eps_v[3] / 2.0, eps_v[2]],
            ], dtype=float)
            D = np.eye(3) + strain_tensor

        # ==============================================================
        # Apply deformation to the box and (optionally) carry atoms.
        # ==============================================================
        if affine:
            # fractional coords computed from the OLD box (before it changes)
            self.coords2fracts(From_Cart=True, normalization=False)

        newmatrix = Distortion.apply_distortion(D, orgbox)
        self.box, _ = lattice_2_lmpbox(Lattice(newmatrix))

        if affine:
            # remap Cartesian coords through the NEW box (fractional preserved)
            self.fracts2coords()
            dropcols = [c for c in ("xsn", "ysn", "zsn") if c in self.atoms.columns]
            if dropcols:
                self.atoms = self.atoms.drop(columns=dropcols)

        if reinit:
            self.initialization(normalization=False, style=1)

        new_lengths = np.asarray(self.box.lengths, dtype=float).reshape(3,)
        new_volume = float(abs(np.linalg.det(np.asarray(self.box.matrix, dtype=float))))

        return {
            "mode": mode,
            "driven_axes": driven,
            "scale_factors": np.diag(D).tolist(),
            "strain_tensor": strain_tensor.tolist(),
            "deformation_matrix": D.tolist(),
            "old_lengths": old_lengths.tolist(),
            "new_lengths": new_lengths.tolist(),
            "old_volume": old_volume,
            "new_volume": new_volume,
            "volume_ratio": (new_volume / old_volume) if old_volume else None,
        }


########################################################################################################################
#  merge_structure_seamless — deterministic site-identity-based merge
########################################################################################################################

    # ------------------------------------------------------------------
    # Helper: identify site IDs on arbitrary coordinate arrays
    # ------------------------------------------------------------------
    @staticmethod
    def _identify_sites_array(
            xyz: np.ndarray,
            site_basis_frac: np.ndarray,
            unitcell_lengths: np.ndarray,
            *,
            tol: float = 0.15,
    ) -> np.ndarray:
        """
        Assign crystallographic site IDs to an array of Cartesian positions.

        This is a standalone version of identify_sites() that operates on
        raw coordinate arrays rather than self.atoms, so it can be applied
        to arbitrary subsets (boundary atoms, donor chunks, etc.).

        Parameters
        ----------
        xyz : (N, 3) Cartesian positions
        site_basis_frac : (Nsite, 3) fractional basis coordinates in one unit cell
        unitcell_lengths : (3,) Cartesian unit cell lengths [lx, ly, lz]
        tol : float
            Tolerance in LOCAL fractional coordinates for site matching.

        Returns
        -------
        np.ndarray of shape (N,), dtype int
            Site IDs: 1-based for matched sites, 0 for unmatched.
        """
        xyz = np.asarray(xyz, dtype=float)
        if xyz.ndim == 1:
            xyz = xyz.reshape(1, 3)
        basis = np.asarray(site_basis_frac, dtype=float)
        uc = np.asarray(unitcell_lengths, dtype=float).reshape(3,)

        # which small unit cell each atom sits in
        cell_idx = np.floor(xyz / uc.reshape(1, 3)).astype(int)
        local_xyz = xyz - cell_idx * uc.reshape(1, 3)
        frac_local = local_xyz / uc.reshape(1, 3)
        frac_local = frac_local - np.floor(frac_local)  # keep in [0, 1)

        # vectorized distance to all basis sites
        diff = frac_local[:, np.newaxis, :] - basis[np.newaxis, :, :]  # (N, Nbasis, 3)
        diff -= np.round(diff)  # periodic within local unit cell
        d2 = np.sum(diff * diff, axis=2)  # (N, Nbasis)

        best_j = np.argmin(d2, axis=1)
        best_d2 = d2[np.arange(len(xyz)), best_j]

        tol2 = float(tol) ** 2
        return np.where(best_d2 <= tol2, best_j + 1, 0).astype(int)

    # ------------------------------------------------------------------
    # Helper: extract boundary atom indices from a carved region
    # ------------------------------------------------------------------
    @staticmethod
    def _extract_boundary_atoms(
            xyz: np.ndarray,
            center: np.ndarray,
            half_lengths: np.ndarray,
            *,
            mode: str = "corners",
            corner_radius: float | None = None,
            shell_thickness: float = 2.0,
    ) -> np.ndarray:
        """
        Select indices of atoms near the boundary of a carved box.

        Parameters
        ----------
        xyz : (N, 3) atom positions already selected inside the carved region
        center : (3,) center of the carved box
        half_lengths : (3,) half-side-lengths of the carved box (Hx, Hy, Hz)
        mode : str
            Boundary sampling strategy:
            - "corners" : 8 corner zones of the box (default)
            - "shell"   : outer shell (all atoms within shell_thickness of any face)
            Future options: "edges", "faces"
        corner_radius : float or None
            Radius of corner sampling sphere. If None, auto = 0.3 * min(half_lengths).
        shell_thickness : float
            Thickness in Å for shell-based sampling.

        Returns
        -------
        np.ndarray of int
            Indices into xyz of the boundary atoms.
        """
        xyz = np.asarray(xyz, dtype=float)
        center = np.asarray(center, dtype=float).reshape(3,)
        H = np.asarray(half_lengths, dtype=float).reshape(3,)
        mode = str(mode).lower().strip()

        if mode == "corners":
            if corner_radius is None:
                corner_radius = 0.3 * float(np.min(H))
            cr2 = float(corner_radius) ** 2

            # 8 corners of the box
            corners = []
            for sx in (-1, +1):
                for sy in (-1, +1):
                    for sz in (-1, +1):
                        corners.append(center + np.array([sx * H[0], sy * H[1], sz * H[2]]))

            mask = np.zeros(len(xyz), dtype=bool)
            for c in corners:
                d2 = np.sum((xyz - c.reshape(1, 3)) ** 2, axis=1)
                mask |= (d2 <= cr2)
            return np.where(mask)[0]

        elif mode == "shell":
            # atoms within shell_thickness of any face of the box
            rel = np.abs(xyz - center.reshape(1, 3))  # (N, 3)
            # distance from each face = H[i] - |rel[i]|
            face_dist = H.reshape(1, 3) - rel  # positive means inside
            # atom is near boundary if its minimum face distance < shell_thickness
            min_face_dist = np.min(face_dist, axis=1)
            return np.where((min_face_dist >= 0) & (min_face_dist <= shell_thickness))[0]

        else:
            raise ValueError(
                f"Unknown boundary mode '{mode}'. "
                f"Supported: 'corners', 'shell'. "
                f"Future: 'edges', 'faces'."
            )

    # ------------------------------------------------------------------
    # Helper: generate candidate translation vectors from unit cell basis
    # ------------------------------------------------------------------
    @staticmethod
    def _generate_basis_candidates(
            site_basis_frac: np.ndarray,
            unitcell_lengths: np.ndarray,
    ) -> list[np.ndarray]:
        """
        Generate candidate Cartesian translation vectors from the unit cell
        basis positions. Each candidate represents shifting the carve center
        by one basis vector, which cycles through all possible site alignments.

        Parameters
        ----------
        site_basis_frac : (Nsite, 3) fractional basis positions
        unitcell_lengths : (3,) Cartesian unit cell lengths

        Returns
        -------
        list of np.ndarray
            Candidate translation vectors in Cartesian coordinates.
            Always includes [0, 0, 0] as the first candidate.
        """
        basis = np.asarray(site_basis_frac, dtype=float)
        uc = np.asarray(unitcell_lengths, dtype=float).reshape(3,)

        # Each basis position, converted to Cartesian, is a candidate shift
        candidates = [np.zeros(3, dtype=float)]  # identity (no shift)
        for i in range(len(basis)):
            v = basis[i] * uc  # fractional -> Cartesian
            # skip if essentially zero (same as the identity candidate)
            if np.linalg.norm(v) > 1e-6:
                # also skip if duplicate of an existing candidate
                is_dup = False
                for existing in candidates:
                    # check modulo unit cell
                    diff = (v - existing) / uc
                    diff -= np.round(diff)
                    if np.linalg.norm(diff * uc) < 1e-3:
                        is_dup = True
                        break
                if not is_dup:
                    candidates.append(v.copy())
        return candidates

    # ------------------------------------------------------------------
    # Main method: merge_structure_seamless
    # ------------------------------------------------------------------
    def merge_structure_seamless(
            self,
            donor: "Modlmp_LmpData",
            *,
            # ----- Carve region (same shape for both host & donor) -----
            region_type: str = "cube",
            center_host=None,
            center_donor=None,

            # cube/box sizes (shared baseline for both host & donor)
            side_length: float | None = None,
            side_lengths: tuple[float, float, float] | None = None,
            pad_plus=(0.0, 0.0, 0.0),
            pad_minus=(0.0, 0.0, 0.0),

            # Donor-specific padding override (Å).
            # Use when the host defect (e.g. edge dislocation) shifts the
            # boundary inward on one side → donor carve must be larger there
            # to fill the gap.  None = use same pad as host.
            pad_plus_donor: tuple[float, float, float] | None = None,
            pad_minus_donor: tuple[float, float, float] | None = None,

            # sphere
            radius: float | None = None,

            # extra region kwargs (pass-through for future region types)
            host_region_kwargs: dict | None = None,
            donor_region_kwargs: dict | None = None,

            # ----- Site identity matching -----
            site_basis_frac: np.ndarray | None = None,
            unitcell_lengths: np.ndarray | None = None,
            site_tol: float = 0.15,

            # boundary sampling
            boundary_mode: str = "corners",
            corner_radius: float | None = None,
            shell_thickness: float = 2.0,

            # Exclusion zones: boundary atoms within any of these zones are
            # skipped during site-ID comparison.  Use to mask out regions
            # where long-range strain (e.g. edge dislocation) makes site
            # identification unreliable.
            # Format: list of (center_cart, radius) tuples.
            #   e.g. [([x, y, z], 15.0)]  — exclude atoms within 15 Å of (x,y,z)
            boundary_exclude_zones: list[tuple] | None = None,

            # ----- Coordinate mapping -----
            mapping: str = "relative_cart",
            wrap_inserted: bool = True,

            # ----- Human control: manual shift applied AFTER alignment -----
            donor_shift_cart=(0.0, 0.0, 0.0),

            # ----- Human override: skip auto-alignment, use this shift -----
            manual_candidate_shifts_cart: list | None = None,
            skip_auto_alignment: bool = False,

            # ----- Write intermediates -----
            host_deleted_file: str | None = None,
            donor_insert_file: str | None = None,

            # ----- ID / init -----
            reset_ids: bool = True,
            reinit: bool = True,

            # ----- Post-merge overlap removal -----
            # Remove host atoms that end up too close to inserted donor atoms.
            # None = disabled.  Float (e.g. 0.5) = remove host atom if a
            # donor atom is within this distance (Å).
            overlap_tol: float | None = None,
            # Optional spatial filter for overlap detection.
            # Only atoms whose coordinates fall inside this box are checked.
            # Pass a dict with any subset of keys:
            #   {"x": (xlo, xhi), "y": (ylo, yhi), "z": (zlo, zhi)}
            # None (default) = check the entire system.
            overlap_region: dict | None = None,

            # ----- Strictness -----
            require_same_box: bool = False,

            # ----- Logging -----
            verbose: bool = True,
    ) -> dict:
        """
        Merge (transplant) a carved region from `donor` into `self` (host)
        using deterministic site-identity matching at the boundary.

        Workflow
        --------
        1. Carve the same shape from both host and donor.
        2. Generate candidate translation vectors from unit cell basis.
        3. For each candidate, shift the donor carve center, re-carve,
           and compare site IDs of boundary atoms (corners of carved box).
        4. Pick the candidate where boundary site IDs match best.
        5. Transplant the correctly-aligned donor chunk into the host.

        The boundary atoms (corners) are far from the defect center,
        so their site identity is reliable even with defect distortion.

        Parameters
        ----------
        donor : Modlmp_LmpData
            The donor structure to transplant FROM.
        region_type : str
            Shape of the carved region: "cube", "box", "sphere".
            Future: any region_mask-supported type.
        center_host, center_donor : array-like (3,) or None
            Centers of the carved region. None → geometric center.
        side_length / side_lengths : float or (3,)
            Isotropic or anisotropic box sizes for cube/box carving.
        pad_plus, pad_minus : (3,)
            Asymmetric padding (Å) for cube/box carving.
        radius : float
            Radius for sphere carving.
        site_basis_frac : (Nsite, 3) or None
            Fractional basis coordinates in one unit cell.
            Required for site-identity matching.
        unitcell_lengths : (3,) or None
            Cartesian unit cell lengths. Required for site-identity matching.
        site_tol : float
            Tolerance in local fractional coordinates for site matching.
        boundary_mode : str
            How to sample boundary atoms: "corners", "shell".
            Future: "edges", "faces".
        corner_radius : float or None
            Radius for corner sampling. None → auto.
        shell_thickness : float
            Thickness (Å) for shell-based boundary sampling.
        mapping : str
            Coordinate mapping: "relative_cart" or "relative_frac".
        wrap_inserted : bool
            Wrap inserted atoms into host box via PBC.
        donor_shift_cart : (3,)
            Manual Cartesian shift applied to ALL inserted atoms AFTER
            alignment. For human fine-tuning. Default (0,0,0).
        manual_candidate_shifts_cart : list of (3,) or None
            If provided, use these shifts instead of auto-generated
            basis candidates. For human override.
        skip_auto_alignment : bool
            If True, skip site-identity matching entirely and use
            the first manual_candidate_shift (or zero shift).
            For cases where the user knows alignment is correct.
        host_deleted_file, donor_insert_file : str or None
            Write intermediate structures before merge.
        reset_ids, reinit : bool
            Post-merge cleanup.
        require_same_box : bool
            Raise if host and donor boxes differ.
        verbose : bool
            Print progress messages.

        Returns
        -------
        dict with merge statistics and diagnostics.
        """

        def _log(msg: str):
            if verbose:
                print(f"[merge_seamless] {msg}")

        # ==================================================================
        # 1) Validation
        # ==================================================================
        _log("Start")
        if self.atoms is None or len(self.atoms) == 0:
            raise ValueError("Host (self) has no atoms.")
        if donor.atoms is None or len(donor.atoms) == 0:
            raise ValueError("Donor has no atoms.")
        for obj, name in ((self, "Host"), (donor, "Donor")):
            if not all(c in obj.atoms.columns for c in ("x", "y", "z")):
                raise ValueError(f"{name} atoms must contain x, y, z.")

        if require_same_box:
            A = np.asarray(self.box.matrix, dtype=float)
            B = np.asarray(donor.box.matrix, dtype=float)
            if not np.allclose(A, B, rtol=0.0, atol=1e-8):
                raise ValueError("Host and donor box matrices differ.")

        # site identity matching requires basis + unit cell
        do_site_matching = (site_basis_frac is not None
                            and unitcell_lengths is not None
                            and not skip_auto_alignment)
        if do_site_matching:
            site_basis_frac = np.asarray(site_basis_frac, dtype=float)
            if site_basis_frac.ndim != 2 or site_basis_frac.shape[1] != 3:
                raise ValueError("site_basis_frac must have shape (Nsite, 3).")
            unitcell_lengths = np.asarray(unitcell_lengths, dtype=float).reshape(3,)
            if np.any(unitcell_lengths <= 0):
                raise ValueError("unitcell_lengths must be positive.")
            _log(f"Site matching ON: {len(site_basis_frac)} basis sites, "
                 f"unit cell = {unitcell_lengths.tolist()}")
        else:
            if skip_auto_alignment:
                _log("Site matching SKIPPED (skip_auto_alignment=True)")
            else:
                _log("Site matching OFF (no site_basis_frac / unitcell_lengths)")

        # centers
        if center_host is None:
            center_host = self.get_center("cart")
        center_host = np.asarray(center_host, dtype=float).reshape(3,)

        if center_donor is None:
            center_donor = donor.get_center("cart")
        center_donor = np.asarray(center_donor, dtype=float).reshape(3,)

        # mapping
        mapping = str(mapping).lower().strip()
        if mapping in ("relative_cart", "cart", "rel_cart"):
            mapping_norm = "relative_cart"
        elif mapping in ("relative_frac", "frac", "rel_frac"):
            mapping_norm = "relative_frac"
        else:
            raise ValueError(f"mapping must be 'relative_cart' or 'relative_frac', got '{mapping}'.")

        # shifts
        donor_shift_cart = np.asarray(donor_shift_cart, dtype=float).reshape(3,)
        pad_plus = np.asarray(pad_plus, dtype=float).reshape(3,)
        pad_minus = np.asarray(pad_minus, dtype=float).reshape(3,)

        # Donor-specific pads: fall back to shared pads if not specified
        pad_plus_d = np.asarray(pad_plus_donor, dtype=float).reshape(3,) if pad_plus_donor is not None else pad_plus.copy()
        pad_minus_d = np.asarray(pad_minus_donor, dtype=float).reshape(3,) if pad_minus_donor is not None else pad_minus.copy()

        # Parse exclusion zones
        exclude_zones = []
        if boundary_exclude_zones is not None:
            for zone in boundary_exclude_zones:
                zc = np.asarray(zone[0], dtype=float).reshape(3,)
                zr = float(zone[1])
                exclude_zones.append((zc, zr))

        rt = str(region_type).lower().strip()

        _log(f"region_type={rt}, mapping={mapping_norm}")
        _log(f"center_host={center_host.tolist()}")
        _log(f"center_donor={center_donor.tolist()}")
        _log(f"donor_shift_cart={donor_shift_cart.tolist()} (manual, post-alignment)")
        if pad_plus_donor is not None or pad_minus_donor is not None:
            _log(f"Donor-specific pad_plus={pad_plus_d.tolist()}, pad_minus={pad_minus_d.tolist()}")
        if exclude_zones:
            _log(f"Boundary exclude zones: {len(exclude_zones)}")
            for zc, zr in exclude_zones:
                _log(f"  exclude center={zc.tolist()}, radius={zr}")

        # ==================================================================
        # 2) Box matrices for mapping
        # ==================================================================
        host_invT = np.asarray(self.box.inv_matrix, dtype=float).T
        host_matT = np.asarray(self.box.matrix, dtype=float).T
        donor_invT = np.asarray(donor.box.inv_matrix, dtype=float).T
        donor_matT = np.asarray(donor.box.matrix, dtype=float).T

        def _wrap_cart_host(Pcart: np.ndarray) -> np.ndarray:
            f = Pcart @ host_invT
            f -= np.floor(f)
            return f @ host_matT

        def _map_positions(P_cart: np.ndarray, c_host: np.ndarray, c_don: np.ndarray) -> np.ndarray:
            c_host = np.asarray(c_host, float).reshape(3,)
            c_don = np.asarray(c_don, float).reshape(3,)
            if mapping_norm == "relative_cart":
                return P_cart + (c_host - c_don).reshape(1, 3)
            else:
                fP = P_cart @ donor_invT
                fC = c_don.reshape(1, 3) @ donor_invT
                df = fP - fC
                df -= np.round(df)
                dcart = df @ host_matT
                return c_host.reshape(1, 3) + dcart

        # ==================================================================
        # 3) Build region kwargs for carving
        # ==================================================================
        def _build_region_kwargs(
                rt: str, center: np.ndarray, extra_kwargs: dict | None,
                *,
                override_pad_plus: np.ndarray | None = None,
                override_pad_minus: np.ndarray | None = None,
        ) -> dict:
            rk = {} if extra_kwargs is None else dict(extra_kwargs)
            center = np.asarray(center, float).reshape(3,)
            pp = pad_plus if override_pad_plus is None else override_pad_plus
            pm = pad_minus if override_pad_minus is None else override_pad_minus

            if rt == "sphere":
                if radius is None:
                    raise ValueError("Sphere carving requires radius.")
                rk["center"] = center
                rk["radius"] = float(radius)
                return rk

            if rt in ("cube", "box"):
                if "cube_params" not in rk or rk["cube_params"] is None:
                    cp = {"center": center.tolist()}
                    if side_lengths is not None:
                        cp["side_lengths"] = tuple(map(float, side_lengths))
                    elif side_length is not None:
                        cp["side_length"] = float(side_length)
                    else:
                        raise ValueError("Cube/box carving requires side_length or side_lengths.")
                    cp["pad_plus"] = tuple(map(float, pp.tolist()))
                    cp["pad_minus"] = tuple(map(float, pm.tolist()))
                    rk["cube_params"] = cp
                else:
                    # user provided cube_params, just update center
                    cp = dict(rk["cube_params"])
                    cp["center"] = center.tolist()
                    cp.setdefault("pad_plus", tuple(map(float, pp.tolist())))
                    cp.setdefault("pad_minus", tuple(map(float, pm.tolist())))
                    rk["cube_params"] = cp
                return rk

            # generic region types
            rk.setdefault("center", center)
            return rk

        def _shift_cube_center(rk: dict, shift: np.ndarray) -> dict:
            """Return a copy of region_kwargs with cube center shifted."""
            rk = dict(rk)
            if "cube_params" in rk and rk["cube_params"] is not None:
                cp = dict(rk["cube_params"])
                if "center" in cp:
                    c = np.asarray(cp["center"], float).reshape(3,)
                    cp["center"] = (c + shift).tolist()
                rk["cube_params"] = cp
            elif "center" in rk:
                rk["center"] = np.asarray(rk["center"], float).reshape(3,) + shift
            return rk

        # compute half-lengths for boundary extraction
        def _get_half_lengths() -> np.ndarray:
            if side_lengths is not None:
                return np.array(side_lengths, dtype=float) / 2.0 + pad_plus / 2.0 + pad_minus / 2.0
            elif side_length is not None:
                base = np.full(3, float(side_length)) / 2.0
                return base + pad_plus / 2.0 + pad_minus / 2.0
            elif radius is not None:
                return np.full(3, float(radius))
            else:
                return np.full(3, 5.0)  # fallback

        half_lengths = _get_half_lengths()

        # ==================================================================
        # 4) Host carve (one-time)
        # ==================================================================
        _log("Carving host region")
        host_rk = _build_region_kwargs(rt, center_host, host_region_kwargs)
        host_mask = self.region_mask(rt, **host_rk)
        if host_mask is None or len(host_mask) != len(self.atoms):
            raise ValueError("Host region_mask returned invalid mask.")

        n_host_removed = int(np.count_nonzero(host_mask))
        _log(f"Host carved: {n_host_removed} atoms selected for removal")

        # Host boundary atoms (for site comparison)
        host_carved_xyz = self.atoms.loc[host_mask, ["x", "y", "z"]].to_numpy(dtype=float)
        host_boundary_idx = self._extract_boundary_atoms(
            host_carved_xyz, center_host, half_lengths,
            mode=boundary_mode, corner_radius=corner_radius,
            shell_thickness=shell_thickness,
        )
        _log(f"Host boundary atoms: {len(host_boundary_idx)}")

        # Host boundary site IDs
        if do_site_matching and len(host_boundary_idx) > 0:
            host_boundary_xyz = host_carved_xyz[host_boundary_idx]
            host_boundary_sites = self._identify_sites_array(
                host_boundary_xyz, site_basis_frac, unitcell_lengths, tol=site_tol,
            )

            # Apply exclusion zones: mask out boundary atoms near known
            # defects where long-range strain corrupts site identification
            boundary_reliable = np.ones(len(host_boundary_xyz), dtype=bool)
            for zc, zr in exclude_zones:
                d2 = np.sum((host_boundary_xyz - zc.reshape(1, 3)) ** 2, axis=1)
                excluded = d2 <= zr ** 2
                n_excl = int(np.count_nonzero(excluded))
                if n_excl > 0:
                    _log(f"Exclude zone ({zc.tolist()}, r={zr}): "
                         f"masking {n_excl} host boundary atoms")
                boundary_reliable &= ~excluded

            # Zero out site IDs for unreliable atoms so they are
            # skipped during matching (h_site == 0 → skip)
            host_boundary_sites[~boundary_reliable] = 0

            _log(f"Host boundary site IDs: {host_boundary_sites.tolist()}")
            n_host_identified = int(np.count_nonzero(host_boundary_sites > 0))
            n_total_boundary = len(host_boundary_sites)
            n_excluded = int(np.count_nonzero(~boundary_reliable))
            _log(f"Host boundary: {n_host_identified} identified, "
                 f"{n_excluded} excluded, {n_total_boundary} total")

            if n_host_identified == 0:
                _log("WARNING: All host boundary atoms excluded or unidentified! "
                     "Site matching will be ineffective. Consider reducing "
                     "exclude radius or enlarging carve region.")
        else:
            host_boundary_xyz = np.empty((0, 3))
            host_boundary_sites = np.array([], dtype=int)
            boundary_reliable = np.array([], dtype=bool)

        # ==================================================================
        # 5) Generate candidate shifts
        # ==================================================================
        if manual_candidate_shifts_cart is not None:
            candidates = [np.asarray(v, float).reshape(3,) for v in manual_candidate_shifts_cart]
            _log(f"Using {len(candidates)} manual candidate shifts")
        elif do_site_matching:
            candidates = self._generate_basis_candidates(site_basis_frac, unitcell_lengths)
            _log(f"Auto-generated {len(candidates)} candidate shifts from basis")
        else:
            candidates = [np.zeros(3, dtype=float)]
            _log("No site matching → single candidate (zero shift)")

        # ==================================================================
        # 6) Score each candidate: compare boundary site IDs
        # ==================================================================
        best_shift = candidates[0]
        best_score = -1.0
        best_donor_mask = None
        best_donor_center = center_donor
        score_details = []

        for ci, shift in enumerate(candidates):
            cd_try = center_donor + shift

            # carve donor with shifted center (use donor-specific pads)
            donor_rk = _build_region_kwargs(
                rt, cd_try, donor_region_kwargs,
                override_pad_plus=pad_plus_d,
                override_pad_minus=pad_minus_d,
            )
            dm = donor.region_mask(rt, **donor_rk)
            if dm is None or len(dm) != len(donor.atoms):
                raise ValueError("Donor region_mask returned invalid mask.")

            n_donor_sel = int(np.count_nonzero(dm))

            if not do_site_matching or len(host_boundary_idx) == 0:
                # no site matching → use first candidate
                score = 1.0 if n_donor_sel > 0 else 0.0
                score_details.append({
                    "shift": shift.tolist(), "score": score,
                    "donor_selected": n_donor_sel,
                })
                if score > best_score:
                    best_score = score
                    best_shift = shift
                    best_donor_mask = dm
                    best_donor_center = cd_try
                break  # no matching, just use the single candidate

            # -- Site-identity comparison --
            donor_carved_xyz = donor.atoms.loc[dm, ["x", "y", "z"]].to_numpy(dtype=float)

            # Map donor boundary atoms into host frame for comparison
            donor_mapped_xyz = _map_positions(donor_carved_xyz, center_host, cd_try)

            # Extract boundary from the donor chunk (same geometry)
            donor_boundary_idx = self._extract_boundary_atoms(
                donor_mapped_xyz, center_host, half_lengths,
                mode=boundary_mode, corner_radius=corner_radius,
                shell_thickness=shell_thickness,
            )

            if len(donor_boundary_idx) == 0:
                score_details.append({
                    "shift": shift.tolist(), "score": 0.0,
                    "donor_selected": n_donor_sel,
                    "donor_boundary": 0,
                    "reason": "no donor boundary atoms",
                })
                continue

            donor_boundary_mapped = donor_mapped_xyz[donor_boundary_idx]

            # Identify site IDs on the mapped donor boundary atoms
            # (using HOST coordinates, because they're now in host frame)
            donor_boundary_sites = self._identify_sites_array(
                donor_boundary_mapped, site_basis_frac, unitcell_lengths, tol=site_tol,
            )

            # Match: for each host boundary atom, find nearest donor boundary
            # atom and compare site IDs
            n_match = 0
            n_compare = 0
            for hi in range(len(host_boundary_xyz)):
                h_site = host_boundary_sites[hi]
                if h_site == 0:
                    continue  # host atom unidentified, skip

                # find nearest donor boundary atom
                d2 = np.sum((donor_boundary_mapped - host_boundary_xyz[hi]) ** 2, axis=1)
                nearest = int(np.argmin(d2))
                d_site = donor_boundary_sites[nearest]

                n_compare += 1
                if h_site == d_site:
                    n_match += 1

            score = float(n_match) / max(n_compare, 1)

            _log(f"  candidate[{ci}] shift={shift.tolist()}, "
                 f"donor_sel={n_donor_sel}, boundary_match={n_match}/{n_compare}, "
                 f"score={score:.3f}")

            score_details.append({
                "shift": shift.tolist(), "score": score,
                "donor_selected": n_donor_sel,
                "donor_boundary": len(donor_boundary_idx),
                "match": n_match, "compared": n_compare,
            })

            if score > best_score:
                best_score = score
                best_shift = shift
                best_donor_mask = dm
                best_donor_center = cd_try

        _log(f"Best candidate: shift={best_shift.tolist()}, score={best_score:.3f}")

        # If no candidate was selected (all had 0 donor atoms), fallback
        if best_donor_mask is None:
            _log("WARNING: No valid candidate found, using zero-shift fallback")
            best_donor_center = center_donor
            donor_rk = _build_region_kwargs(rt, center_donor, donor_region_kwargs)
            best_donor_mask = donor.region_mask(rt, **donor_rk)

        # ==================================================================
        # 7) Prepare final host_keep and donor_sel
        # ==================================================================
        n_removed = int(np.count_nonzero(host_mask))
        n_selected = int(np.count_nonzero(best_donor_mask))
        _log(f"Final: host_removed={n_removed}, donor_selected={n_selected}")

        # Atom count sanity check
        count_diff = abs(n_selected - n_removed)
        if count_diff > 0:
            pct_diff = 100.0 * count_diff / max(n_removed, 1)
            if pct_diff > 5.0:
                _log(f"WARNING: Atom count mismatch: host_removed={n_removed}, "
                     f"donor_selected={n_selected} (diff={count_diff}, {pct_diff:.1f}%). "
                     f"Consider using pad_plus_donor/pad_minus_donor to extend the "
                     f"donor carve on the defect side.")
            else:
                _log(f"Atom count diff={count_diff} ({pct_diff:.1f}%) — within tolerance")

        host_keep = self.atoms.loc[~host_mask].copy()
        donor_sel = donor.atoms.loc[best_donor_mask].copy()
        if len(donor_sel) == 0:
            raise ValueError("Donor selection is empty; nothing to insert.")

        # Map donor atoms into host frame
        _log("Map donor atoms into host frame")
        P2 = donor_sel[["x", "y", "z"]].to_numpy(dtype=float)
        Pnew = _map_positions(P2, center_host, best_donor_center)

        # Apply manual post-alignment shift
        if np.linalg.norm(donor_shift_cart) > 0.0:
            _log(f"Apply manual donor_shift_cart = {donor_shift_cart.tolist()}")
            Pnew = Pnew + donor_shift_cart.reshape(1, 3)

        if wrap_inserted:
            _log("Wrap inserted atoms into host box")
            Pnew = _wrap_cart_host(Pnew)

        donor_sel.loc[:, ["x", "y", "z"]] = Pnew

        # ==================================================================
        # 8) Write intermediates
        # ==================================================================
        if host_deleted_file is not None:
            _log(f"Write host_deleted_file: {host_deleted_file}")
            tmp = self.copy()
            tmp.atoms = host_keep.copy()
            tmp.to_file(host_deleted_file, to_atom_style=True)

        if donor_insert_file is not None:
            _log(f"Write donor_insert_file: {donor_insert_file}")
            tmp = donor.copy()
            tmp.atoms = donor_sel.copy()
            tmp.box = self.box
            tmp.to_file(donor_insert_file, to_atom_style=True)

        # ==================================================================
        # 9) Remove overlapping atoms in the entire merged system
        # ==================================================================
        n_overlap_removed = 0
        if overlap_tol is not None and float(overlap_tol) > 0:
            _log(f"Overlap removal: tol={overlap_tol} Å (whole-system check)")
            otol = float(overlap_tol)
            otol2 = otol * otol

            # Full arrays for the combined (to-be-merged) system
            donor_xyz = donor_sel[["x", "y", "z"]].to_numpy(dtype=float)
            host_xyz  = host_keep[["x", "y", "z"]].to_numpy(dtype=float)

            # ------------------------------------------------------------------
            # Optional region filter
            # Only atoms that fall inside the user-specified axis-aligned box
            # are considered as *candidates for removal* (host) or as
            # *reference atoms* (donor).  Atoms outside the region are kept as-is.
            # ------------------------------------------------------------------
            if overlap_region is not None:
                def _in_region(xyz_arr: np.ndarray, region: dict) -> np.ndarray:
                    """Return boolean mask: True if atom is inside *region*."""
                    mask = np.ones(len(xyz_arr), dtype=bool)
                    axes = {"x": 0, "y": 1, "z": 2}
                    for ax, col in axes.items():
                        if ax in region:
                            lo, hi = float(region[ax][0]), float(region[ax][1])
                            mask &= (xyz_arr[:, col] >= lo) & (xyz_arr[:, col] <= hi)
                    return mask

                host_in_region  = _in_region(host_xyz,  overlap_region)
                donor_in_region = _in_region(donor_xyz, overlap_region)

                # Summarise the filter for the user
                region_desc = ", ".join(
                    f"{ax}: {overlap_region[ax][0]:.2f}–{overlap_region[ax][1]:.2f}"
                    for ax in ("x", "y", "z") if ax in overlap_region
                )
                _log(f"  overlap_region filter: {region_desc}")
                _log(f"  host atoms in region : {int(np.count_nonzero(host_in_region))} / {len(host_xyz)}")
                _log(f"  donor atoms in region: {int(np.count_nonzero(donor_in_region))} / {len(donor_xyz)}")

                # Only check overlaps within the region
                host_check_xyz  = host_xyz[host_in_region]
                donor_check_xyz = donor_xyz[donor_in_region]
            else:
                host_in_region  = np.ones(len(host_xyz),  dtype=bool)
                donor_in_region = np.ones(len(donor_xyz), dtype=bool)
                host_check_xyz  = host_xyz
                donor_check_xyz = donor_xyz
                _log("  No overlap_region specified → checking entire system")

            # ------------------------------------------------------------------
            # Fast voxel-grid overlap detection
            # Build lookup on *donor* atoms (reference), then sweep host atoms.
            # ------------------------------------------------------------------
            host_too_close = np.zeros(len(host_xyz), dtype=bool)

            if len(donor_check_xyz) > 0 and len(host_check_xyz) > 0:
                # Build voxel map on the donor reference atoms
                vox_d = {}
                vkeys = np.floor(donor_check_xyz / otol).astype(int)
                for i, k in enumerate(vkeys):
                    kk = (int(k[0]), int(k[1]), int(k[2]))
                    vox_d.setdefault(kk, []).append(i)

                # Check only the host atoms that are inside the region
                host_in_indices = np.where(host_in_region)[0]
                hkeys = np.floor(host_check_xyz / otol).astype(int)

                for local_hi, (global_hi, hk) in enumerate(
                        zip(host_in_indices, hkeys)):
                    hkk = (int(hk[0]), int(hk[1]), int(hk[2]))
                    found = False
                    for dx in (-1, 0, 1):
                        for dy in (-1, 0, 1):
                            for dz in (-1, 0, 1):
                                nbr = (hkk[0] + dx, hkk[1] + dy, hkk[2] + dz)
                                didxs = vox_d.get(nbr)
                                if not didxs:
                                    continue
                                d2 = np.sum(
                                    (donor_check_xyz[didxs] - host_check_xyz[local_hi]) ** 2,
                                    axis=1,
                                )
                                if np.any(d2 < otol2):
                                    found = True
                                    break
                            if found:
                                break
                        if found:
                            break
                    if found:
                        host_too_close[global_hi] = True

            n_overlap_removed = int(np.count_nonzero(host_too_close))
            if n_overlap_removed > 0:
                _log(f"Removing {n_overlap_removed} host atoms overlapping with donor atoms")
                host_keep = host_keep.loc[~host_too_close].copy()
            else:
                _log("No overlapping host atoms found")

        # ==================================================================
        # 10) Merge
        # ==================================================================
        _log("Append donor chunk into host")
        self.atoms = pd.concat([host_keep, donor_sel], axis=0)

        if reinit:
            _log("Initialization")
            self.initialization(normalization=False, style=1)

        if reset_ids:
            _log("Reset atom IDs")
            self.reset_atom_ids()

        _log("Done")
        return {
            "host_removed": n_removed,
            "donor_selected": n_selected,
            "final_atoms": int(len(self.atoms)),
            "region_type": rt,
            "mapping": mapping_norm,
            "wrap_inserted": bool(wrap_inserted),
            "manual_donor_shift_cart": donor_shift_cart.tolist(),
            "best_alignment_shift": best_shift.tolist(),
            "best_alignment_score": float(best_score),
            "center_host": center_host.tolist(),
            "center_donor_original": center_donor.tolist(),
            "center_donor_aligned": best_donor_center.tolist(),
            "n_candidates_tested": len(score_details),
            "candidate_scores": score_details,
            "site_matching_used": bool(do_site_matching),
            "boundary_mode": boundary_mode,
            "n_host_boundary_atoms": len(host_boundary_idx),
            "n_boundary_excluded": int(np.count_nonzero(~boundary_reliable)) if len(boundary_reliable) > 0 else 0,
            "pad_plus_donor_used": pad_plus_d.tolist(),
            "pad_minus_donor_used": pad_minus_d.tolist(),
            "n_overlap_removed": n_overlap_removed,
            "wrote_host_deleted_file": host_deleted_file,
            "wrote_donor_insert_file": donor_insert_file,
        }

    def copy(self):
        return copy.deepcopy(self)

    def identify_sites(
            self,
            site_basis_frac,
            unitcell_lengths,
            *,
            tol: float = 0.02,
            write_to_molecule_id: bool = True,
            reinit: bool = False,
    ):
        """
        Identify basis-site IDs by dividing the large orthogonal box into
        repeated small unit cells, then matching each atom by its LOCAL
        fractional coordinate inside the small cell.

        Parameters
        ----------
        site_basis_frac : array-like, shape (Nsite, 3)
            Fractional basis coordinates inside ONE reference unit cell.

        unitcell_lengths : (3,) 
            Cartesian box lengths of the reference unit cell:
            [lx_uc, ly_uc, lz_uc]

        tol : float
            Tolerance in LOCAL fractional coordinates.
            Small value means strict matching.

        write_to_molecule_id : bool
            If True, write site ID into "molecule-ID".

        reinit : bool
            If True, call initialization() after writing labels.

        Returns
        -------
        np.ndarray
            Site IDs aligned with self.atoms.
            0 means unmatched / defect / distorted atom.
        """
        if self.atoms is None or len(self.atoms) == 0:
            return np.array([], dtype=int)

        if not all(c in self.atoms.columns for c in ("x", "y", "z")):
            raise ValueError("Atoms table must contain x, y, z columns.")

        # -------------------------
        # validate basis
        # -------------------------
        basis = np.asarray(site_basis_frac, dtype=float)
        if basis.ndim != 2 or basis.shape[1] != 3:
            raise ValueError("site_basis_frac must have shape (Nsite, 3).")

        # -------------------------
        # validate unit cell lengths
        # -------------------------
        uc = np.asarray(unitcell_lengths, dtype=float).reshape(3, )
        if np.any(uc <= 0):
            raise ValueError("unitcell_lengths must be positive.")

        # -------------------------
        # large box lengths
        # -------------------------
        bigL = np.asarray(self.box.lengths, dtype=float)

        # this method assumes orthogonal replication logic
        if not self.box.is_orthogonal:
            raise ValueError(
                "identify_sites(unitcell_lengths=...) currently assumes an orthogonal box."
            )

        # number of small cells along each direction
        ncell = bigL / uc
        ncell_round = np.rint(ncell)

        if not np.allclose(ncell, ncell_round, atol=1e-6):
            raise ValueError(
                f"Large box lengths {bigL.tolist()} are not integer multiples of "
                f"unitcell_lengths {uc.tolist()}."
            )

        ncell = ncell_round.astype(int)

        # -------------------------
        # atom positions
        # -------------------------
        xyz = self.atoms[["x", "y", "z"]].to_numpy(dtype=float)

        # which small cell the atom is in
        cell_index = np.floor(xyz / uc.reshape(1, 3)).astype(int)

        # clamp edge cases for atoms exactly on upper boundary
        for i in range(3):
            cell_index[:, i] = np.clip(cell_index[:, i], 0, ncell[i] - 1)

        # local Cartesian coord inside that small cell
        local_xyz = xyz - cell_index * uc.reshape(1, 3)

        # local fractional coord inside one small cell
        frac_local = local_xyz / uc.reshape(1, 3)

        # numerical cleanup: keep in [0,1)
        frac_local = frac_local - np.floor(frac_local)

        # -------------------------
        # assign site IDs
        # -------------------------
        tol2 = float(tol) ** 2
        n_atoms = len(frac_local)
        n_basis = len(basis)

        # Vectorized: compute all atom-basis distances at once
        # diff[i, j, :] = frac_local[i] - basis[j] (broadcast)
        diff = frac_local[:, np.newaxis, :] - basis[np.newaxis, :, :]  # (N, Nbasis, 3)
        diff -= np.round(diff)  # periodic compare inside local unit cell
        d2 = np.sum(diff * diff, axis=2)  # (N, Nbasis)

        best_j = np.argmin(d2, axis=1)  # (N,)
        best_d2 = d2[np.arange(n_atoms), best_j]  # (N,)

        site_ids = np.where(best_d2 <= tol2, best_j + 1, 0).astype(int)

        # -------------------------
        # write labels
        # -------------------------
        if write_to_molecule_id:
            if "molecule-ID" not in self.atoms.columns:
                self.atom_style = "molecular"
                self.atoms["molecule-ID"] = 0
            self.atoms["molecule-ID"] = site_ids

        if reinit:
            self.initialization(normalization=False, style=1)

        return site_ids



















