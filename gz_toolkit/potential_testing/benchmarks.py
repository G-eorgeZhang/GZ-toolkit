"""Opt-in extended benchmarks with explicit geometry, references and normalization.

All geometries are prepared in Python. No external calculators run here.
Dislocation energies are finite-cell excesses, not isolated core energies.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import numpy as np
from pymatgen.core import Element, Lattice, Structure
from pymatgen.core.surface import SlabGenerator

from gz_toolkit.core.modlmp import Modlmp_LmpData
from gz_toolkit.defect import BCCDefect, FCCDefect
from .config import PotentialConfig, resolve_composition


def _safe_name(value):
    value = str(value)
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", value) or value in (".", ".."):
        raise ValueError(f"Invalid benchmark name: {value!r}")
    return value


def _counts(data, cfg):
    return {el: int((data.atoms['type'] == tid).sum()) for el, tid in cfg.potential.type_map.items()}


def _from_structure(structure, cfg):
    from .structure_ops import _assert_full_force_field
    data = Modlmp_LmpData.from_pmg_structure(structure, 'atomic')
    original = sorted(Element(el) for el in structure.symbol_set)
    mapping = {i + 1: cfg.potential.type_map[el.symbol] for i, el in enumerate(original)}
    data.atoms['type'] = data.atoms['type'].map(mapping).astype(int)
    _assert_full_force_field(data, cfg)
    return data


def _pmg(data, cfg):
    symbols = {v: k for k, v in cfg.potential.type_map.items()}
    lo = np.asarray(data.box.bounds)[:, 0]
    return Structure(Lattice(data.box.matrix), [symbols[int(t)] for t in data.atoms['type']],
                     data.atoms[['x', 'y', 'z']].to_numpy() - lo, coords_are_cartesian=True)


def _input(cfg, *, boundary='p p p', fixed_ids=(), anneal=None):
    from .pipeline import render_min_input
    text = render_min_input(box_relax=False, minimize_tol=cfg.workflow.minimize_tol)
    old_tol = f'{cfg.workflow.minimize_tol:g}'
    text = text.replace(f'minimize        {old_tol} {old_tol} 100000 100000',
                        f'minimize 0.0 {cfg.workflow.benchmark_force_tolerance:g} 100000 100000')
    text = text.replace('boundary        p p p', f'boundary        {boundary}')
    extra = ''
    if fixed_ids:
        extra += 'group fixed id ' + ' '.join(str(i) for i in fixed_ids) + '\nfix hold fixed setforce 0 0 0\n'
        extra += 'group mobile subtract all fixed\n'
    if anneal:
        t, steps = float(anneal['temperature_K']), int(anneal['steps'])
        if t <= 0 or steps < 1:
            raise ValueError('Anneal temperature and steps must be positive.')
        group = 'mobile' if fixed_ids else 'all'
        extra += f'velocity {group} create {t:g} 1729 mom yes rot no\ntimestep 0.001\nfix anneal {group} nvt temp {t:g} {t:g} 0.1\nrun {steps}\nunfix anneal\n'
    text = text.replace('min_style       cg', extra + 'min_style       cg')
    return text


def _save(case, data, cfg, root, metadata, *, boundary='p p p', fixed_ids=(), anneal=None):
    from .pipeline import _write_potential_include, _copy_potential_files
    from .structure_ops import _write_structure
    case.mkdir(parents=True, exist_ok=True)
    if cfg.workflow.benchmark_force_tolerance <= 0:
        raise ValueError('benchmark_force_tolerance must be positive.')
    if any(int(i) not in data.atoms.index for i in fixed_ids):
        raise ValueError('fixed_ids must refer to atoms present in the constructed cell.')
    _write_structure(data, case / 'structure.data')
    _write_potential_include(case, cfg.potential.pot_lines, cfg.potential.masses)
    _copy_potential_files(root, cfg.potential.pot_name, case, cfg.potential.potential_files)
    metadata['counts'] = _counts(data, cfg)
    hashes = {name: hashlib.sha256((case / name).read_bytes()).hexdigest()
              for name in ('structure.data', 'potential.inc', *cfg.potential.potential_files)
              if (case / name).is_file()}
    text = _input(cfg, boundary=boundary, fixed_ids=fixed_ids, anneal=anneal)
    case_id = hashlib.sha256((text + json.dumps(hashes, sort_keys=True)).encode()).hexdigest()
    text = f'print "BENCHMARK_CASE {case_id}"\n' + text + '\nprint "BENCHMARK_COMPLETE"\n'
    (case / 'in.benchmark.lammps').write_text(text, encoding='utf-8')
    hashes['in.benchmark.lammps'] = hashlib.sha256((case / 'in.benchmark.lammps').read_bytes()).hexdigest()
    metadata.update(case_id=case_id, file_hashes=hashes)
    (case / 'benchmark.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    return {'kind': metadata['kind'], 'system': metadata['system'], 'dir': str(case.relative_to(root))}


def _host(case, cfg, element, lc, size, structure=None):
    from .structure_ops import _build_unit_cell, _scale_cell_to_lc
    case.mkdir(parents=True, exist_ok=True)
    structure = structure or cfg.potential.crystal_structures.get(element, {}).get('structure', 'bcc')
    data = _build_unit_cell(case, cfg, structure, element, size)
    _scale_cell_to_lc(data, lc)
    return data, structure


def build_benchmarks(cfg: PotentialConfig, pot_dir: Path, root: Path, lattice: dict) -> list[dict]:
    """Build configured surfaces, C15/SFT, dislocations, molecules and alloy gas complexes."""
    from .structure_ops import build_gas_complex_structure, _decorate_alloy_types
    manifest = []
    if cfg.workflow.include_surfaces:
        if cfg.workflow.surface_thickness <= 0 or cfg.workflow.surface_vacuum <= 0:
            raise ValueError('Surface thickness and vacuum must be positive.')
        for el in cfg.potential.single_elements:
            base = pot_dir / 'benchmarks' / el
            seed, _ = _host(base / '_surface_seed', cfg, el, lattice[el], 1)
            bulk = _pmg(seed, cfg)
            for miller in cfg.workflow.surface_millers:
                slabs = SlabGenerator(bulk, tuple(miller), cfg.workflow.surface_thickness,
                    cfg.workflow.surface_vacuum, center_slab=True, primitive=False).get_slabs(symmetrize=True)
                if not slabs:
                    raise ValueError(f'No symmetric slab found for {el} {miller}.')
                for termination, slab in enumerate(slabs):
                    tag = 'surface_' + '_'.join(map(str, miller)) + f'_t{termination}'
                    area = float(np.linalg.norm(np.cross(slab.lattice.matrix[0], slab.lattice.matrix[1])))
                    manifest.append(_save(base / tag, _from_structure(slab, cfg), cfg, root,
                        {'kind': 'surface', 'system': el, 'miller': list(miller), 'area_A2': area, 'n_surfaces': 2}))
    used = set()
    for index, spec in enumerate(cfg.workflow.benchmark_cases):
        kind = spec['kind']
        el = spec.get('element')
        system = _safe_name(spec.get('system') or el or ''.join(spec.get('elements', [])))
        name = _safe_name(spec.get('name', f'{kind}_{index}'))
        case = pot_dir / 'benchmarks' / system / name
        if case in used or (case / 'benchmark.json').is_file() and name.startswith('surface_'):
            raise ValueError(f'Duplicate benchmark case: {case}')
        used.add(case)
        meta = {'kind': kind, 'system': system, 'name': name}
        if kind == 'molecule':
            elements = list(spec['elements'])
            if len(elements) != 2:
                raise ValueError('A molecular reference requires two elements.')
            distance, box = float(spec.get('distance_A', 1.0)), float(spec.get('box_A', 40))
            if distance <= 0 or box <= 2 * distance:
                raise ValueError('Require positive bond distance and box_A > 2 * distance_A.')
            data = _from_structure(Structure(Lattice.cubic(box), elements,
                [[box / 2 - distance / 2, box / 2, box / 2], [box / 2 + distance / 2, box / 2, box / 2]], coords_are_cartesian=True), cfg)
            meta['elements'] = elements
            # Separate isolated atoms use exactly the same box and potential.
            for element in set(elements):
                isolated = _from_structure(Structure(Lattice.cubic(box), [element], [[box / 2] * 3], coords_are_cartesian=True), cfg)
                atom_dir = case.parent / f'{name}_atom_{element}'
                manifest.append(_save(atom_dir, isolated, cfg, root, {'kind': 'isolated_atom', 'system': system, 'element': element}))
            meta['atom_reference_dirs'] = {e: f'{name}_atom_{e}' for e in set(elements)}
        elif kind == 'alloy_gas':
            composition = resolve_composition(spec['composition'])
            host = next(iter(composition))
            gas = spec['gas']
            if gas in composition:
                raise ValueError('The gas must not be part of the matrix composition.')
            case.mkdir(parents=True, exist_ok=True)
            structure = spec.get('structure') or cfg.potential.crystal_structures[host]['structure']
            seed = int(spec.get('seed', 1729))
            n_gas, n_def = int(spec['n_gas']), int(spec['m_defect'])
            if n_gas < 1 or n_def < 0:
                raise ValueError('alloy_gas needs n_gas >= 1 and m_defect >= 0.')
            defect_kind = spec.get('defect_kind', 'vacancy')
            # Build the same chemical host first, then decorate only metal atoms.
            for suffix, ng, nd in (('', n_gas, n_def), ('_bulk', 0, 0), ('_defect', 0, n_def), ('_gas', n_gas, 0)):
                out = case if not suffix else case.parent / (name + suffix)
                out.mkdir(parents=True, exist_ok=True)
                build_gas_complex_structure(out, host, gas, ng, nd, defect_kind, cfg, structure,
                    relaxed_lc=lattice[host], n_replicate=int(spec.get('size', cfg.potential.size_alloy)), composition_atpct=composition, seed=seed)
                candidate = Modlmp_LmpData.from_file(str(out / 'structure.data'), 'atomic')
                extra = {**meta, 'counts': _counts(candidate, cfg), 'composition': composition}
                if suffix:
                    extra['kind'] = 'alloy_gas_reference'
                else:
                    extra.update(reference_dirs={key: name + val for key, val in [('bulk', '_bulk'), ('defect', '_defect'), ('gas', '_gas')]}, gas=gas,
                                 chemical_potentials=spec.get('chemical_potentials', {}))
                manifest.append(_save(out, candidate, cfg, root, extra))
            continue
        elif kind == 'surface':
            if spec.get('equivalent_surfaces') is not True:
                raise ValueError('Imported surface energies require equivalent_surfaces=true; asymmetric terminations cannot be split equally.')
            data = Modlmp_LmpData.from_file(str(root / spec['structure_file']), 'atomic')
            area = float(np.linalg.norm(np.cross(data.box.matrix[0], data.box.matrix[1])))
            meta.update(area_A2=area, n_surfaces=2, chemical_potentials=spec.get('chemical_potentials', {}))
        elif kind in ('c15', 'sft', 'dislocation'):
            size = int(spec.get('size', cfg.potential.size_single))
            if el not in lattice:
                raise ValueError(f'{kind} needs an elemental reference for {el}.')
            data, structure = _host(case, cfg, el, lattice[el], size)
            before = len(data.atoms)
            ff = [e for e, _ in sorted(cfg.potential.type_map.items(), key=lambda kv: kv[1])]
            masses = [cfg.potential.masses[cfg.potential.type_map[e]] for e in ff]
            if kind == 'c15':
                if structure != 'bcc':
                    raise ValueError('C15 seeds are currently supported in BCC hosts only.')
                meta['construction'] = BCCDefect(data).add_c15_cluster(int(spec['nsias']), lattice_const=lattice[el], ff_elements=ff, atomic_masses=masses)
            elif kind == 'sft':
                if structure != 'fcc':
                    raise ValueError('SFT vacancy platelets require FCC.')
                meta['construction'] = FCCDefect(data).stacking_fault_tetrahedron(edge_length=float(spec['edge_length_A']), center=spec.get('center', data.get_center('cart')), lattice_const=lattice[el], ff_elements=ff, atomic_masses=masses)
                meta['structural_validation_required'] = 'Confirm the relaxed structure formed an SFT, not a loop or void.'
            else:
                # Use explicitly oriented or imported cells; preserve a matched
                # pristine reference BEFORE inserting a line defect.
                if 'structure_file' not in spec:
                    raise ValueError('Dislocation benchmark requires an oriented structure_file and current_axes.')
                data = Modlmp_LmpData.from_file(str(root / spec['structure_file']), 'atomic')
                if data.atoms['type'].nunique() != 1 or structure not in ('bcc', 'fcc'):
                    raise ValueError('Generated line defects currently need a pure BCC/FCC host.')
                before = len(data.atoms)
                builder = (BCCDefect if structure == 'bcc' else FCCDefect)(data)
                kwargs = dict(spec.get('construction', {}))
                report = builder.dislocation(spec.get('line_kind', 'screw'), lattice_const=lattice[el], ff_elements=ff, **kwargs)
                # An arbitrary surface-containing reference cannot be inferred
                # from the cut geometry. Require a reviewed matched reference.
                if 'reference_file' not in spec or 'line_axis' not in spec or 'n_lines' not in spec:
                    raise ValueError('Dislocations need reference_file, line_axis and n_lines for matched excess energies.')
                reference = Modlmp_LmpData.from_file(str(root / spec['reference_file']), 'atomic')
                if not np.allclose(reference.box.matrix, data.box.matrix):
                    raise ValueError('The dislocation and its matched pristine reference must have the same box.')
                reference_dir = case.parent / (name + '_reference')
                boundary = spec.get('boundary', 'p p p')
                manifest.append(_save(reference_dir, reference, cfg, root, {'kind': 'line_reference', 'system': system}, boundary=boundary, fixed_ids=spec.get('reference_fixed_ids', [])))
                axis, count = int(spec['line_axis']), int(spec['n_lines'])
                if axis not in (0, 1, 2) or count < 1:
                    raise ValueError('line_axis must be 0/1/2 and n_lines positive.')
                meta.update(reference_dir=reference_dir.name, reference_counts=_counts(reference, cfg),
                            length_A=float(np.linalg.norm(data.box.matrix[axis])), n_lines=count, construction=report)
            meta['net_defects'] = len(data.atoms) - before
            if kind in ('c15', 'sft') and meta['net_defects'] == 0:
                raise ValueError('The requested cluster selected no atoms; increase its size or the cell.')
        else:
            raise ValueError(f'Unsupported benchmark kind: {kind}')
        manifest.append(_save(case, data, cfg, root, meta, boundary=spec.get('boundary', 'p p p'),
                              fixed_ids=spec.get('fixed_ids', []), anneal=spec.get('anneal')))
    return manifest


def read_result(case):
    """Only complete matching runs become physical metrics; no failure-to-zero conversion."""
    from .summary import _extract_results
    case = Path(case)
    meta = json.loads((case / 'benchmark.json').read_text())
    log = case / 'log.lammps'
    if not log.exists():
        return meta, None, 'pending'
    text = log.read_text(errors='replace')
    if 'ERROR' in text:
        return meta, None, 'failed'
    if 'Total wall time' not in text:
        return meta, None, 'running'
    if f"BENCHMARK_CASE {meta['case_id']}" not in text or 'BENCHMARK_COMPLETE' not in text or not (case / 'relaxed.data').exists():
        return meta, None, 'failed'
    for name, digest in meta['file_hashes'].items():
        if not (case / name).is_file() or hashlib.sha256((case / name).read_bytes()).hexdigest() != digest:
            return meta, None, 'failed'
    results = _extract_results(log)
    if not all(key in results and np.isfinite(results[key]) for key in ('pe', 'natoms')):
        return meta, None, 'failed'
    if int(results['natoms']) != sum(meta['counts'].values()):
        return meta, None, 'failed'
    force_stops = re.findall(r'Stopping criterion\s*=\s*(.*)', text)
    if not force_stops or 'force tolerance' not in force_stops[-1].lower():
        return meta, None, 'unconverged'
    return meta, results, 'done'


def collect_benchmark_metrics(cfg, pot_dir, reference_energies):
    values, systems = {}, {}
    def put(key, value, system):
        values[key], systems[key] = value, system
    for case in sorted((pot_dir / 'benchmarks').glob('*/*')):
        if not (case / 'benchmark.json').exists():
            continue
        meta, results, state = read_result(case)
        tag, system, kind = case.name, meta['system'], meta['kind']
        put(f'benchmark_status_{system}_{tag}', state, system)
        if not results or kind in ('isolated_atom', 'alloy_gas_reference', 'line_reference'):
            continue
        e = results['pe']
        mu = dict(reference_energies)
        mu.update(meta.get('chemical_potentials', {}))
        if kind == 'molecule':
            atom_refs = {el: read_result(case.parent / path)[1] for el, path in meta['atom_reference_dirs'].items()}
            if all(atom_refs.values()):
                binding = e - sum(atom_refs[el]['pe'] for el in meta['elements'])
                put(f'molecule_binding_{system}_{tag}', binding, system)
                put(f'molecule_dissociation_{system}_{tag}', -binding, system)
            relaxed = Modlmp_LmpData.from_file(str(case / 'relaxed.data'), 'atomic')
            xyz = relaxed.atoms[['x', 'y', 'z']].to_numpy()
            if len(xyz) == 2:
                delta = (xyz[1] - xyz[0]) @ np.linalg.inv(relaxed.box.matrix)
                delta -= np.rint(delta)
                put(f'molecule_bond_length_{system}_{tag}', float(np.linalg.norm(delta @ relaxed.box.matrix)), system)
        elif kind == 'alloy_gas':
            refs = {key: read_result(case.parent / path)[1] for key, path in meta['reference_dirs'].items()}
            if all(refs.values()):
                put(f'alloy_gas_binding_{system}_{tag}', e + refs['bulk']['pe'] - refs['defect']['pe'] - refs['gas']['pe'], system)
        elif kind == 'dislocation':
            _, reference, _ = read_result(case.parent / meta['reference_dir'])
            if reference and all(el in mu for el, n in meta['counts'].items() if n != meta['reference_counts'].get(el, 0)):
                excess = e - reference['pe'] - sum((n - meta['reference_counts'].get(el, 0)) * mu.get(el, 0) for el, n in meta['counts'].items())
                put(f'line_excess_per_length_{system}_{tag}', excess / meta['length_A'] / meta['n_lines'], system)
        elif all(el in mu for el, count in meta['counts'].items() if count):
            excess = e - sum(count * mu.get(el, 0) for el, count in meta['counts'].items())
            if kind == 'surface':
                put(f'surface_energy_J_m2_{system}_{tag}', excess / (meta['n_surfaces'] * meta['area_A2']) * 16.02176634, system)
            else:
                put(f'cluster_formation_{system}_{tag}', excess, system)
                put(f'cluster_formation_per_defect_{system}_{tag}', excess / abs(meta['net_defects']), system)
    return values, systems
