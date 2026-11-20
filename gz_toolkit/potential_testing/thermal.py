"""NPT thermal-expansion cases and explicit acceptance of a stationary final window."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import PotentialConfig, resolve_composition

INPUT_NAME = "in.thermal.lammps"


def validate_thermal_options(config: PotentialConfig) -> list[str]:
    wf = config.workflow
    if not wf.include_thermal_lc:
        return (["thermal_lc_for_alloys requires include_thermal_lc."] if wf.thermal_lc_for_alloys else [])
    errors = []
    integers = {
        "size_thermal_lc": config.potential.size_thermal_lc,
        **{name: getattr(wf, name) for name in (
            "thermal_lc_steps", "thermal_lc_average_steps", "thermal_lc_sample_every", "thermal_lc_block_steps")},
    }
    for name, value in integers.items():
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            errors.append(f"{name} must be a positive integer.")
    for name in ("thermal_lc_timestep_ps", "thermal_lc_tdamp_steps", "thermal_lc_pdamp_steps",
                 "thermal_lc_max_relative_drift", "thermal_lc_max_block_spread", "thermal_lc_temperature_tolerance"):
        value = getattr(wf, name)
        if not isinstance(value, (int, float)) or not np.isfinite(value) or value <= 0:
            errors.append(f"{name} must be finite and positive.")
    if not isinstance(wf.thermal_lc_pressure_bar, (int, float)) or not np.isfinite(wf.thermal_lc_pressure_bar):
        errors.append("thermal_lc_pressure_bar must be finite (metal units: bar).")
    temps = wf.thermal_lc_temperatures
    if not isinstance(temps, list) or not temps or any(
        isinstance(t, bool) or not isinstance(t, (int, float)) or not np.isfinite(t) or t <= 0 for t in temps
    ):
        errors.append("thermal_lc_temperatures must be a nonempty list of positive finite temperatures in K.")
    elif len(set(temps)) != len(temps):
        errors.append("thermal_lc_temperatures must not contain duplicates.")
    if wf.thermal_lc_barostat not in ("iso", "aniso"):
        errors.append("thermal_lc_barostat must be iso (preserve shape) or aniso.")
    if isinstance(wf.thermal_lc_seed, bool) or not isinstance(wf.thermal_lc_seed, int) or not 1 <= wf.thermal_lc_seed < 900000000:
        errors.append("thermal_lc_seed must be an integer from 1 to 899999999.")
    if not any("positive integer" in error for error in errors):
        total, window, sample, block = (wf.thermal_lc_steps, wf.thermal_lc_average_steps,
                                         wf.thermal_lc_sample_every, wf.thermal_lc_block_steps)
        if window >= total:
            errors.append("thermal_lc_average_steps must be smaller than thermal_lc_steps.")
        if any(value % sample for value in (total, window, block)):
            errors.append("Thermal run, averaging window, and block lengths must be multiples of thermal_lc_sample_every.")
        if window % block or window // block < 4 or block // sample < 2:
            errors.append("The thermal averaging window needs at least four equal blocks, each with at least two samples.")
    if wf.thermal_lc_for_alloys and not (config.potential.two_element_suites or config.potential.multi_element_suites):
        errors.append("thermal_lc_for_alloys needs configured alloy compositions.")
    return errors


def render_thermal_input(config: PotentialConfig, temperature: float, seed: int) -> str:
    wf = config.workflow
    tdamp = wf.thermal_lc_tdamp_steps * wf.thermal_lc_timestep_ps
    pdamp = wf.thermal_lc_pdamp_steps * wf.thermal_lc_timestep_ps
    return f'''# NPT thermal lattice constant; structures are built by Python.
units metal
dimension 3
boundary p p p
atom_style atomic
read_data structure.data
include potential.inc
neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
thermo {wf.thermal_lc_sample_every}
thermo_style custom step temp press pe vol lx ly lz
thermo_modify lost error flush yes
min_style cg
minimize {wf.minimize_tol:g} {wf.minimize_tol:g} 100000 100000
reset_timestep 0
timestep {wf.thermal_lc_timestep_ps:g}
velocity all create {temperature:g} {seed} mom yes rot no dist gaussian
fix thermal_npt all npt temp {temperature:g} {temperature:g} {tdamp:g} {wf.thermal_lc_barostat} {wf.thermal_lc_pressure_bar:g} {wf.thermal_lc_pressure_bar:g} {pdamp:g}
variable box_x equal lx
variable box_y equal ly
variable box_z equal lz
variable measured_temp equal temp
variable measured_press equal press
fix thermal_history all ave/time {wf.thermal_lc_sample_every} 1 {wf.thermal_lc_sample_every} v_box_x v_box_y v_box_z v_measured_temp v_measured_press file thermal_lc.dat ave one format "%20.16g"
run {wf.thermal_lc_steps}
write_data thermal_final.data
print "THERMAL_LC_COMPLETE {wf.thermal_lc_steps}"
'''


def build_thermal_cases(config: PotentialConfig, pot_dir: Path, project_root: Path,
                        per_element_lc: dict[str, float]) -> list[dict[str, Any]]:
    from .pipeline import _copy_potential_files, _write_potential_include
    from .structure_ops import build_reference_structure, build_alloy_structure, build_multi_alloy_structure
    from mylammps.inputs.data import lmpData

    if not config.workflow.include_thermal_lc:
        return []
    errors = validate_thermal_options(config)
    if errors:
        raise ValueError("Invalid thermal-lc configuration: " + "; ".join(errors))
    pot, wf = config.potential, config.workflow
    n = pot.size_thermal_lc
    specs = []
    for element in pot.single_elements or [next(iter(pot.type_map))]:
        specs.append((element, {element: 100.0}, pot.crystal_structures.get(element, {}).get("structure", "bcc"), "random"))
    if wf.thermal_lc_for_alloys:
        for suite in pot.two_element_suites:
            a, b = suite["A"], suite["B"]
            for fraction in suite.get("fractions_atpct", []):
                composition = {a: 100 - float(fraction), b: float(fraction)}
                tag = f"{a}{b}_" + "_".join(f"{el}{pct:g}" for el, pct in composition.items())
                specs.append((tag, composition, suite.get("structure") or pot.crystal_structures.get(a, {}).get("structure", "bcc"), suite.get("ordering", "random")))
        for suite in pot.multi_element_suites:
            composition = resolve_composition(suite["composition"])
            host = next(iter(composition))
            tag = "".join(composition) + "_" + "_".join(f"{el}{pct:g}" for el, pct in composition.items())
            specs.append((tag, composition, suite.get("structure") or pot.crystal_structures.get(host, {}).get("structure", "bcc"), "random"))
    manifest = []
    for tag, composition, structure, ordering in specs:
        host = next(iter(composition))
        if host not in per_element_lc:
            raise ValueError(f"Thermal case {tag} needs a 0 K reference for {host}.")
        for temperature in wf.thermal_lc_temperatures:
            case = pot_dir / "thermal_lc" / tag / f"{temperature:g}K"
            case.mkdir(parents=True, exist_ok=True)
            digest = int(hashlib.sha256(f"{pot.pot_name}/{tag}/{temperature:g}/{wf.thermal_lc_seed}".encode()).hexdigest()[:12], 16)
            seed = digest % 899999999 + 1
            # Use the same alloy arrangement at each temperature; velocity seeds differ.
            alloy_seed = int(hashlib.sha256(f"{tag}/{wf.thermal_lc_seed}".encode()).hexdigest()[:8], 16)
            lc = per_element_lc[host]
            if len(composition) == 1:
                build_reference_structure(case, config, host, structure, lc_value=lc, n_replicate=n)
            elif len(composition) == 2:
                a, b = composition
                build_alloy_structure(case, config, a, b, composition[b], structure, ordering=ordering,
                                      relaxed_lc_A=lc, n_replicate=n, seed=alloy_seed)
            else:
                build_multi_alloy_structure(case, config, composition, structure, ordering="random",
                                            relaxed_lc_host=lc, n_replicate=n, seed=alloy_seed)
            _write_potential_include(case, pot.pot_lines, pot.masses)
            _copy_potential_files(project_root, pot.pot_name, case, pot.potential_files)
            text = render_thermal_input(config, float(temperature), seed)
            file_hashes = {name: hashlib.sha256((case / name).read_bytes()).hexdigest()
                           for name in ("structure.data", "potential.inc", *pot.potential_files)
                           if (case / name).is_file()}
            case_id = hashlib.sha256((text + json.dumps(file_hashes, sort_keys=True)).encode()).hexdigest()
            text = f'print "THERMAL_LC_CASE {case_id}"\n' + text
            (case / INPUT_NAME).write_text(text, encoding="utf-8")
            file_hashes[INPUT_NAME] = hashlib.sha256((case / INPUT_NAME).read_bytes()).hexdigest()
            data = lmpData.from_file(str(case / "structure.data"), "atomic")
            initial_lengths = np.asarray(data.box.lengths, dtype=float).reshape(3)
            metadata = {
                "system": "".join(composition), "tag": tag, "composition": composition,
                "structure": structure, "temperature_K": temperature, "replications": n,
                "initial_lc_A": lc, "initial_box_A": initial_lengths.tolist(), "velocity_seed": seed,
                "case_id": case_id, "file_hashes": file_hashes,
                "options": {k: v for k, v in asdict(wf).items() if k.startswith("thermal_lc_")},
            }
            (case / "thermal_case.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
            manifest.append({"system": tag, "temperature_K": temperature, "dir": str(case.relative_to(project_root))})
    return manifest


def analyze_thermal_case(case_dir: Path) -> dict[str, Any]:
    """Read raw samples; never accept a truncated, drifting, or shorter selected window."""
    case = Path(case_dir)
    result: dict[str, Any] = {"state": "pending", "converged": False, "reasons": []}
    try:
        meta = json.loads((case / "thermal_case.json").read_text(encoding="utf-8"))
        result.update({key: meta[key] for key in ("system", "tag", "temperature_K")})
        opts = meta["options"]
        log = case / "log.lammps"
        if not log.is_file():
            return result
        text = log.read_text(encoding="utf-8", errors="replace")
        for name, digest in meta.get("file_hashes", {}).items():
            if hashlib.sha256((case / name).read_bytes()).hexdigest() != digest:
                raise ValueError(f"Input {name} changed since this thermal case was generated.")
        if f"THERMAL_LC_CASE {meta['case_id']}" not in text:
            raise ValueError("The log does not belong to the current thermal case inputs.")
        if "ERROR" in text:
            return {**result, "state": "failed", "reasons": ["LAMMPS reported an error."]}
        if "Total wall time" not in text:
            return {**result, "state": "running"}
        total = opts["thermal_lc_steps"]
        if f"THERMAL_LC_COMPLETE {total}" not in text or not (case / "thermal_final.data").is_file():
            raise ValueError("Missing completion marker or thermal_final.data.")
        history = np.loadtxt(case / "thermal_lc.dat", comments="#", ndmin=2)
        if history.shape[1] != 6 or not np.isfinite(history).all():
            raise ValueError("Thermal history must contain six finite columns: step, lx, ly, lz, temperature, pressure.")
        sample, window, block = (opts["thermal_lc_sample_every"], opts["thermal_lc_average_steps"], opts["thermal_lc_block_steps"])
        expected = np.arange(sample, total + 1, sample)
        if not np.array_equal(history[:, 0], expected):
            raise ValueError("Thermal history is incomplete, duplicated, or has unexpected timesteps.")
        if np.any(history[:, 1:4] <= 0):
            raise ValueError("Thermal history contains nonpositive cell lengths.")
        tail = history[history[:, 0] > total - window]
        lengths = tail[:, 1:4]
        means = lengths.mean(axis=0)
        blocks = lengths.reshape(window // block, block // sample, 3).mean(axis=1)
        times = tail[:, 0] - tail[:, 0].mean()
        slopes = (times[:, None] * lengths).sum(axis=0) / (times * times).sum()
        drift = np.abs(slopes) * window / means
        spread = np.ptp(blocks, axis=0) / means
        reasons = []
        if np.any(drift > opts["thermal_lc_max_relative_drift"]):
            reasons.append("The final averaging window still has excessive length drift.")
        if np.any(spread > opts["thermal_lc_max_block_spread"]):
            reasons.append("Block means in the final averaging window are not flat enough.")
        measured_temp = float(tail[:, 4].mean())
        if abs(measured_temp / meta["temperature_K"] - 1) > opts["thermal_lc_temperature_tolerance"]:
            reasons.append("The final-window mean temperature is outside the target tolerance.")
        scales = lengths / np.asarray(meta["initial_box_A"])
        if opts["thermal_lc_barostat"] == "iso" and not np.allclose(scales, scales[:, [0]], rtol=1e-6, atol=0):
            reasons.append("The cell shape changed despite requesting isotropic expansion.")
        result.update({
            "state": "unconverged" if reasons else "done", "converged": not reasons, "reasons": reasons,
            "window_start_exclusive": total - window, "window_end": total, "samples": len(tail),
            "mean_temperature_K": measured_temp, "mean_pressure_bar": float(tail[:, 5].mean()),
            "relative_drift_xyz": drift.tolist(), "relative_block_spread_xyz": spread.tolist(),
            "candidate_box_A": means.tolist(), "block_means_box_A": blocks.tolist(),
            "box_std_A": lengths.std(axis=0, ddof=1).tolist(),
        })
        if not reasons:
            result["box_A"] = means.tolist()
            result["unit_cell_edges_A"] = (means / meta["replications"]).tolist()
            # Shape-preserving thermal expansion scales the prototype's a0;
            # for non-cubic cells, raw repeat edges are not three copies of a0.
            result["lattice_axes_A"] = (means / np.asarray(meta["initial_box_A"]) * meta["initial_lc_A"]).tolist()
            if opts["thermal_lc_barostat"] == "iso":
                result["lc_A"] = float(np.mean(result["lattice_axes_A"]))
        return result
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return {**result, "state": "failed", "reasons": [str(exc)]}
