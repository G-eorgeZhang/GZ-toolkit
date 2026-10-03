"""Formation/binding/cohesive energy and elastic-modulus formulas.

Every function here is pure: energies and atom counts in, a derived number
out. None of them touch a filesystem or know about the potential-testing
directory layout — that extraction lives in
``gz_toolkit.potential_testing.summary``.
"""

from __future__ import annotations


def formation_energy(e_cell: float, n_ref_species_in_cell: float, e_ref: float, n_ref: int) -> float:
    """Defect/complex formation energy relative to a bulk reference cell.

        E_f = E_cell - (n_ref_species_in_cell / n_ref) * E_ref

    ``e_cell``: total potential energy of the cell containing the
    defect/complex.
    ``n_ref_species_in_cell``: how many atoms of the *reference* species are
    in that cell — e.g. the metal atom count for a gas-defect complex, not
    the cell's total atom count if other species (gas) are mixed in.
    ``e_ref`` / ``n_ref``: total energy / atom count of the pure reference
    cell (same species as ``n_ref_species_in_cell`` counts).
    """
    return e_cell - (n_ref_species_in_cell / n_ref) * e_ref


def binding_energy(ef_complex: float, ef_reduced: float, ef_single: float) -> float:
    """Binding energy of one constituent to a complex.

        E_b = E_f(complex) - E_f(complex minus one constituent) - E_f(single constituent)

    ``ef_reduced`` is the formation energy of the same complex with one
    fewer of that constituent; ``ef_single`` is the formation energy of that
    constituent on its own (e.g. one vacancy, or one gas atom in bulk).
    """
    return ef_complex - ef_reduced - ef_single


def cohesive_energy(e_per_atom: float, e_single_atom: float) -> float:
    """True cohesive energy: bulk energy/atom minus the isolated-atom energy.

    LAMMPS's own "ECOH" print is really just energy-per-atom of the bulk
    cell — this is what turns that into the actual cohesive energy, using
    the isolated single-atom energy from the same potential.
    """
    return e_per_atom - e_single_atom


# The 9 diagonal + coupled normal/shear Cij every cubic-or-lower symmetry
# elastic-constant run reports; anything else (C14, C15, ...) is an
# off-diagonal coupling not used by voigt_moduli.
MAIN_CIJ = ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66")


def voigt_moduli(cij: dict[str, float]) -> dict[str, float]:
    """Voigt-average bulk/shear/Young's modulus + Poisson ratio from the 9 main Cij.

    Returns ``{}`` if ``cij`` doesn't have all 9 of ``MAIN_CIJ`` (e.g. a
    lower-symmetry structure where some components weren't computed).
    """
    if not all(k in cij for k in MAIN_CIJ):
        return {}
    diag = cij["C11"] + cij["C22"] + cij["C33"]
    off = cij["C12"] + cij["C13"] + cij["C23"]
    shear = cij["C44"] + cij["C55"] + cij["C66"]
    B = (diag + 2.0 * off) / 9.0
    G = (diag - off + 3.0 * shear) / 15.0
    out = {"B": B, "G": G}
    denom = 3.0 * B + G
    if denom != 0:
        out["E"] = 9.0 * B * G / denom
        out["nu"] = (3.0 * B - 2.0 * G) / (2.0 * denom)
    return out
