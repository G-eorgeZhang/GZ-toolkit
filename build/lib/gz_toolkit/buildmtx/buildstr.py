from pymatgen.core import Lattice, Structure
from pymatgen.io.lammps.data import LammpsData
from mylammps.inputs.data import lmpData
from typing import Literal
import math
import numpy as np
import os
import copy

class Gen_crystal:

    def __init__(self,atom_style: Literal["atomic", "charge", "full"] = "atomic",):
        self.atom_style = atom_style

        self.PROTOTYPES = {

            # ---------- A TYPES ----------
            "A1": {  # FCC
                "lattice": lambda: Lattice.cubic(1.0),
                "coords": [
                    [0, 0, 0],
                    [0, 0.5, 0.5],
                    [0.5, 0, 0.5],
                    [0.5, 0.5, 0],
                ],
            },

            "A2": {  # BCC
                "lattice": lambda: Lattice.cubic(1.0),
                "coords": [
                    [0, 0, 0],
                    [0.5, 0.5, 0.5],
                ],
            },

            "A3": {  # HCP
                "lattice": lambda: Lattice.hexagonal(
                    1.0,
                    math.sqrt(8 / 3),
                ),
                "coords": [
                    [0, 0, 0],
                    [2 / 3, 1 / 3, 0.5],
                ],
            },

            "A4": {  # Diamond
                "lattice": lambda: Lattice.cubic(1.0),
                "coords": [
                    [0, 0, 0],
                    [0.25, 0.25, 0.25],
                    [0.5, 0.5, 0],
                    [0.75, 0.75, 0.25],
                    [0.5, 0, 0.5],
                    [0.75, 0.25, 0.75],
                    [0, 0.5, 0.5],
                    [0.25, 0.75, 0.75],
                ],
            },

            # ---------- B TYPES ----------
            "B1": {  # Rocksalt
                "lattice": lambda: Lattice.cubic(1.0),
                "coords": [
                    [0, 0, 0],
                    [0.5, 0.5, 0],
                    [0.5, 0, 0.5],
                    [0, 0.5, 0.5],
                    [0.5, 0, 0],
                    [0, 0.5, 0],
                    [0, 0, 0.5],
                    [0.5, 0.5, 0.5],
                ],
            },

            "B2": {  # CsCl
                "lattice": lambda: Lattice.cubic(1.0),
                "coords": [
                    [0, 0, 0],
                    [0.5, 0.5, 0.5],
                ],
            },

            "B3": {  # Zinc blende
                "lattice": lambda: Lattice.cubic(1.0),
                "coords": [
                    [0, 0, 0],
                    [0.25, 0.25, 0.25],
                    [0.5, 0.5, 0],
                    [0.75, 0.75, 0.25],
                    [0.5, 0, 0.5],
                    [0.75, 0.25, 0.75],
                    [0, 0.5, 0.5],
                    [0.25, 0.75, 0.75],
                ],
            },

            "B4": {  # Wurtzite
                "lattice": lambda: Lattice.hexagonal(
                    1.0,
                    math.sqrt(8 / 3),
                ),
                "coords": [
                    [0, 0, 0],
                    [0, 0, 0.5],
                    [1 / 3, 2 / 3, 0.25],
                    [2 / 3, 1 / 3, 0.75],
                ],
            },

            # ---------- C TYPES ----------
            "C1": {  # Fluorite (CaF2)
                "lattice": lambda: Lattice.cubic(1.0),
                "coords": [
                    [0, 0, 0],
                    [0.5, 0.5, 0],
                    [0.5, 0, 0.5],
                    [0, 0.5, 0.5],
                    [0.25, 0.25, 0.25],
                    [0.25, 0.25, 0.75],
                    [0.25, 0.75, 0.25],
                    [0.75, 0.25, 0.25],
                    [0.75, 0.75, 0.75],
                    [0.75, 0.75, 0.25],
                    [0.75, 0.25, 0.75],
                    [0.25, 0.75, 0.75],
                ],
            },
        }

    # =====================================================
    # Main entry
    # =====================================================

    def seed_crystal(
        self,
        structure_type: str,
        elements: list[str],
        filename: str,
    ):

        structure_type = structure_type.upper()

        if structure_type not in self.PROTOTYPES:
            raise ValueError(
                f"Unknown structure type: {structure_type}"
            )

        proto = self.PROTOTYPES[structure_type]
        coords = proto["coords"]

        if len(elements) != len(coords):
            raise ValueError(
                f"{structure_type} requires "
                f"{len(coords)} elements, "
                f"but got {len(elements)}."
            )

        struct = Structure(
            proto["lattice"](),
            elements,
            coords,
        )


        LammpsData.from_structure(
            struct,
            atom_style=self.atom_style,
        ).write_file(filename)


    def trans_unit(
        self,
        seed_file: str,
        directions: list[list[int]],
        filename: str | None = None,
        *,
        ff_elements: list[str] = ("Fe",),
        atom_style: str = "atomic",
        is_sort: bool = False,
        halve_c: bool = True,
    ):
        """
        Transform a seed unit cell to a new crystallographic orientation.

        Parameters
        ----------
        seed_file : str
            Path to a POSCAR file (e.g. "bcc.POSCAR").
        directions : list of 3 lists of 3 ints
            Supercell transform matrix.
        filename : str or None
            Output file name. If None, derived from seed_file as
            "<seed_stem>_transformed.data".
        ff_elements : list of str
            Force-field element symbols passed to lmpData.from_POSCAR.
        atom_style : str
            Atom style passed to lmpData.from_POSCAR.
        is_sort : bool
            Sorting flag forwarded to lmpData.from_POSCAR.
        halve_c : bool
            If True (historical default), halve the c-axis length of the
            transformed box (used for screw-dislocation cells).

        Returns
        -------
        str
            Path to the written file.
        """
        thisdata = lmpData.from_POSCAR(
            seed_file, atom_style, ff_elements=list(ff_elements), is_sort=is_sort
        )

        outdata = thisdata.make_supercell(directions)
        if halve_c:
            newmatrix = copy.deepcopy(outdata.box.matrix)
            newmatrix[2][2] /= 2.0
            outdata.modify_lmpbox(newmatrix, style=2)
        outdata.reset_atom_ids()
        if filename is None:
            stem = os.path.splitext(os.path.basename(seed_file))[0]
            filename = f"{stem}_transformed.data"
        outdata.to_file(filename)

        return filename
    # =====================================================
    # Replicate — fast numpy supercell builder
    # =====================================================

    def replicate(
        self,
        seed_file: str,
        nx: int,
        ny: int,
        nz: int,
        lattice_parameter: float,
        filename: str,
    ):
        """
        Read a seed .data file, tile it nx × ny × nz times using
        numpy vector operations, and write the result.

        Parameters
        ----------
        seed_file : str
            Path to a LAMMPS data file (the seed unit cell).
        nx, ny, nz : int
            Repetitions along each lattice vector.
        lattice_parameter : float
            Real lattice parameter (Å).  Seed is stored with a=1.
        filename : str
            Output file.  Detected by name:
              *.data          → LAMMPS data
              POSCAR* / *.vasp → VASP POSCAR
        """
        # ----- read seed -----
        seed = LammpsData.from_file(
            seed_file, atom_style=self.atom_style
        ).structure

        elements = [str(sp) for sp in seed.species]
        lattice_matrix = np.array(seed.lattice.matrix)    # (3,3)
        frac_basis = np.array(seed.frac_coords)            # (n_basis, 3)

        # ----- translation grid -----
        gi, gj, gk = np.meshgrid(
            np.arange(nx, dtype=np.float64),
            np.arange(ny, dtype=np.float64),
            np.arange(nz, dtype=np.float64),
            indexing="ij",
        )
        translations = np.column_stack(
            [gi.ravel(), gj.ravel(), gk.ravel()]
        )
        n_cells = len(translations)

        # ----- broadcast: basis + translations → all atoms -----
        all_frac = translations[:, None, :] + frac_basis[None, :, :]
        all_frac = all_frac.reshape(-1, 3)

        # ----- frac → Cartesian -----
        cart = all_frac @ (lattice_matrix * lattice_parameter)

        # ----- species list -----
        all_species = elements * n_cells

        # ----- build pymatgen Structure -----
        super_lattice = Lattice(
            lattice_matrix * lattice_parameter
            * np.array([[nx], [ny], [nz]])
        )
        struct = Structure(
            super_lattice,
            all_species,
            cart,
            coords_are_cartesian=True,
        )

        # ----- write output -----
        fname_lower = filename.lower()
        if fname_lower.endswith(".vasp") or "poscar" in fname_lower:
            struct.to(filename=filename, fmt="poscar")
        else:
            LammpsData.from_structure(
                struct, atom_style=self.atom_style,
            ).write_file(filename)



            
