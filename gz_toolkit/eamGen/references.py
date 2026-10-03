"""Reference jobs with explicit commands, captured outputs and input hashes.

LAMMPS runs arbitrary user-supplied pair styles (including universal ML).
VASP preparation copies user-provided INCAR/KPOINTS/POTCAR; no pseudopotentials
are downloaded or invented. VASP results are imported only after convergence.
"""
from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np

from .data import Configuration


def _write_json(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def prepare_lammps(c: Configuration, directory, type_map: dict[str, int], masses: dict[str, float],
                   potential_lines: str, potential_files=()) -> Path:
    """Write a run-0 energy/force/stress job without running/submitting it.

    Input units and atom_style are metal/atomic. Potential lines must use
    exactly the supplied type map; charge/reaction ML styles need a custom
    external reference adapter instead of this atomic adapter.
    """
    c.validate()
    if set(type_map) != set(masses) or sorted(type_map.values()) != list(range(1, len(type_map) + 1)):
        raise ValueError("type_map must be contiguous from 1 with matching masses.")
    if not set(c.species).issubset(type_map) or not potential_lines.strip():
        raise ValueError("Missing species mapping or potential lines.")
    if any(not np.isfinite(m) or m <= 0 for m in masses.values()):
        raise ValueError("Masses must be finite and positive.")
    cell = np.asarray(c.cell)
    if np.linalg.det(cell) <= 0:
        raise ValueError("LAMMPS adapter requires a right-handed cell.")
    q, r = np.linalg.qr(cell.T)
    signs = np.sign(np.diag(r)); q = q @ np.diag(signs); restricted = (np.diag(signs) @ r).T
    coords = np.asarray(c.positions) @ q
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=False)
    lines = ["eamGen reference", "", f"{len(c.species)} atoms", f"{len(type_map)} atom types", "",
             f"0 {restricted[0,0]:.16g} xlo xhi", f"0 {restricted[1,1]:.16g} ylo yhi", f"0 {restricted[2,2]:.16g} zlo zhi",
             f"{restricted[1,0]:.16g} {restricted[2,0]:.16g} {restricted[2,1]:.16g} xy xz yz", "", "Masses", ""]
    for s, t in sorted(type_map.items(), key=lambda item: item[1]):
        lines.append(f"{t} {masses[s]:.16g}")
    lines += ["", "Atoms # atomic", ""]
    for atom_id, (s, xyz) in enumerate(zip(c.species, coords), 1):
        lines.append(f"{atom_id} {type_map[s]} " + " ".join(f"{v:.16g}" for v in xyz))
    (root / "structure.data").write_text("\n".join(lines) + "\n", encoding="utf-8")
    boundary = " ".join("p" if p else "f" for p in c.pbc)
    script = f"""units metal
atom_style atomic
boundary {boundary}
read_data structure.data
{potential_lines}
thermo_modify norm no
thermo_style custom step pe pxx pyy pzz pxy pxz pyz
thermo_modify format float %.16g
run 0
write_dump all custom forces.dump id fx fy fz modify sort id format float %.16g
variable E equal pe
print "${{E}}" file energy.txt screen no
variable Sxx equal -pxx/1602176.6208
variable Syy equal -pyy/1602176.6208
variable Szz equal -pzz/1602176.6208
variable Sxy equal -pxy/1602176.6208
variable Sxz equal -pxz/1602176.6208
variable Syz equal -pyz/1602176.6208
print "${{Sxx}} ${{Syy}} ${{Szz}} ${{Sxy}} ${{Sxz}} ${{Syz}}" file stress.txt screen no
"""
    (root / "in.reference").write_text(script, encoding="utf-8")
    for file in potential_files:
        source = Path(file)
        if (root / source.name).exists():
            raise ValueError(f"Duplicate/reserved potential filename: {source.name}")
        shutil.copy2(source, root / source.name)
    _write_json(root / "reference.json", {"backend": "lammps", "configuration": asdict(c),
                                          "rotation": q.tolist(), "type_map": type_map})
    return root


def prepare_vasp(c: Configuration, directory, templates, species_order=None, energy_kind="free_energy") -> Path:
    """Prepare a single-point VASP job; POTCAR order must match POSCAR order.

    Verify POTCAR identities against this order yourself before execution.
    For fitting forces, use unrelaxed snapshots and static INCAR settings.
    """
    c.validate()
    if energy_kind not in {"free_energy", "sigma_zero"}:
        raise ValueError("VASP energy_kind must be free_energy or sigma_zero.")
    if not all(c.pbc):
        raise ValueError("VASP snapshots must be periodic (use explicit vacuum for isolated systems).")
    if np.linalg.det(c.cell) <= 0:
        raise ValueError("VASP requires a right-handed cell.")
    order = list(dict.fromkeys(c.species)) if species_order is None else species_order
    if len(set(order)) != len(order) or set(order) != set(c.species):
        raise ValueError("species_order must contain each present species once.")
    template_root = Path(templates)
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        if not (template_root / name).is_file():
            raise FileNotFoundError(template_root / name)
    root = Path(directory); root.mkdir(parents=True, exist_ok=False)
    permutation = [i for s in order for i, element in enumerate(c.species) if element == s]
    lines = [c.name, "1.0"] + [" ".join(f"{v:.16g}" for v in row) for row in c.cell]
    lines += [" ".join(order), " ".join(str(c.species.count(s)) for s in order), "Cartesian"]
    lines += [" ".join(f"{v:.16g}" for v in c.positions[i]) for i in permutation]
    (root / "POSCAR").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for name in ("INCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(template_root / name, root / name)
    _write_json(root / "reference.json", {"backend": "vasp", "configuration": asdict(c),
                                          "permutation": permutation, "species_order": order, "energy_kind": energy_kind})
    return root


def run_reference(directory, command: list[str], timeout=3600) -> Path:
    """Execute exactly the supplied argv, with shell=False and no auto submission.

    Use only for synchronous execution (e.g. ['mpiexec','-n','4','vasp_std']);
    scheduler submissions need a separate job and import after completion.
    No executable/library is installed by this tool.
    """
    if not command or not all(isinstance(x, str) and x for x in command) or timeout <= 0:
        raise ValueError("Supply a nonempty argv and positive timeout.")
    root = Path(directory).resolve()
    meta = json.loads((root / "reference.json").read_text(encoding="utf-8"))
    if (root / "execution.json").exists():
        raise FileExistsError("Reference already executed; prepare a new directory for a new attempt.")
    names = ["INCAR", "KPOINTS", "POTCAR", "POSCAR"] if meta["backend"] == "vasp" else [p.name for p in root.iterdir() if p.is_file()]
    record = {"command": command, "input_sha256": {name: _hash(root / name) for name in names}, "status": "running"}
    _write_json(root / "execution.json", record)
    try:
        with (root / "stdout.txt").open("w", encoding="utf-8") as stdout, (root / "stderr.txt").open("w", encoding="utf-8") as stderr:
            result = subprocess.run(command, cwd=root, shell=False, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
        record.update(returncode=result.returncode, status="completed" if result.returncode == 0 else "failed")
        _write_json(root / "execution.json", record)
        if result.returncode != 0:
            raise RuntimeError(f"Reference command failed ({result.returncode}); inspect {root / 'stderr.txt'}")
    except (OSError, subprocess.TimeoutExpired) as exc:
        record.update(status="failed", error=str(exc)); _write_json(root / "execution.json", record)
        raise
    return root


def import_reference(directory) -> Configuration:
    root = Path(directory).resolve()
    meta = json.loads((root / "reference.json").read_text(encoding="utf-8"))
    c = Configuration(**meta["configuration"])
    provenance = {**c.provenance, "backend": meta["backend"], "directory": str(root)}
    execution = root / "execution.json"
    if execution.exists():
        run = json.loads(execution.read_text(encoding="utf-8"))
        if run["status"] != "completed":
            raise ValueError("Reference execution did not complete successfully.")
        if any(not (root / name).is_file() or _hash(root / name) != digest for name, digest in run["input_sha256"].items()):
            raise ValueError("Reference inputs changed after execution; cannot attribute these outputs.")
        provenance["execution"] = run
    else:
        provenance["execution"] = "external run; confirm command and input provenance"
    if meta["backend"] == "lammps":
        c.energy = float((root / "energy.txt").read_text().strip())
        dump = (root / "forces.dump").read_text().splitlines()
        header = next(i for i, line in enumerate(dump) if line.startswith("ITEM: ATOMS"))
        columns = dump[header].split()[2:]
        rows = np.asarray([[float(v) for v in line.split()] for line in dump[header + 1:] if line.strip()])
        ids = rows[:, columns.index("id")].astype(int)
        if sorted(ids.tolist()) != list(range(1, len(c.species) + 1)):
            raise ValueError("Reference force dump lost/reordered atom identities.")
        local = rows[np.argsort(ids)][:, [columns.index(k) for k in ("fx", "fy", "fz")]]
        q = np.asarray(meta["rotation"])
        c.forces = (local @ q.T).tolist()
        values = np.loadtxt(root / "stress.txt").reshape(6)
        xx, yy, zz, xy, xz, yz = values
        local_stress = np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        c.stress = (q @ local_stress @ q.T).tolist() if all(c.pbc) else None
        outputs = ["energy.txt", "forces.dump", "stress.txt"]
    elif meta["backend"] == "vasp":
        from pymatgen.io.vasp.outputs import Vasprun
        v = Vasprun(root / "vasprun.xml", parse_dos=False, parse_eigen=False, parse_potcar_file=False)
        if not v.converged_electronic or not v.converged_ionic:
            raise ValueError("VASP output is not converged.")
        step = v.ionic_steps[-1]
        positions = np.asarray(v.final_structure.cart_coords)
        expected = np.asarray(c.positions)[meta["permutation"]]
        if [str(s) for s in v.final_structure.species] != [c.species[i] for i in meta["permutation"]]:
            raise ValueError("VASP species do not match the prepared snapshot.")
        if not np.allclose(v.final_structure.lattice.matrix, c.cell, atol=1e-5) or not np.allclose(positions, expected, atol=1e-5):
            raise ValueError("VASP geometry changed; use a static run for this snapshot.")
        energy_kind = meta.get("energy_kind", "free_energy")
        energy_key = "e_fr_energy" if energy_kind == "free_energy" else "e_0_energy"
        c.energy = float(step[energy_key])
        forces = np.zeros((len(c.species), 3)); forces[meta["permutation"]] = step["forces"]
        c.forces = forces.tolist()
        # VASP stress: compression-positive kbar -> tension-positive eV/A^3.
        c.stress = (-np.asarray(step["stress"]) / 1602.1766208).tolist()
        provenance["energy_kind"] = f"VASP {energy_key} ({energy_kind})"
        outputs = ["vasprun.xml"]
    else:
        raise ValueError("Unknown reference backend.")
    provenance["output_sha256"] = {name: _hash(root / name) for name in outputs}
    c.provenance = provenance
    c.validate()
    return c
