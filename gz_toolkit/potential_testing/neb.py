"""LAMMPS climbing-image NEB for reviewed, atom-ID-matched relaxed endpoints."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
from gz_toolkit.core.modlmp import Modlmp_LmpData


def build_neb_cases(cfg, pot_dir, root):
    from .pipeline import _write_potential_include, _copy_potential_files
    from .benchmarks import _safe_name
    wf = cfg.workflow
    if not wf.neb_cases:
        return []
    if wf.neb_replicas < 3 or cfg.hpc.ntasks % wf.neb_replicas:
        raise ValueError('NEB requires >=3 replicas and ntasks divisible by neb_replicas.')
    if cfg.hpc.device == 'gpu' and cfg.workflow.paral_degree in (1, 2):
        raise ValueError('GPU NEB needs degree 3/4 and an explicitly reviewed replica/GPU launch command.')
    if cfg.hpc.device == 'gpu':
        from .config import HPCOptions
        if cfg.hpc.neb_run_command == HPCOptions().neb_run_command:
            raise ValueError('Set neb_run_command with replica-aware GPU binding; no default GPU sharing is assumed.')
    if wf.neb_force_tolerance <= 0 or wf.neb_iterations < 100 or wf.neb_iterations % 100:
        raise ValueError('NEB needs a positive force tolerance and iterations that are a positive multiple of 100.')
    manifest, seen = [], set()
    for spec in wf.neb_cases:
        name, system = _safe_name(spec['name']), _safe_name(spec['system'])
        case = pot_dir / 'neb' / system / name
        if case in seen:
            raise ValueError(f'Duplicate NEB case {case}')
        seen.add(case)
        if spec.get('endpoints_relaxed') is not True:
            raise ValueError('Review and relax both NEB endpoints, then explicitly set endpoints_relaxed=true.')
        first = Modlmp_LmpData.from_file(str(root / spec['initial']), 'atomic')
        last = Modlmp_LmpData.from_file(str(root / spec['final']), 'atomic')
        a, b = first.atoms.sort_index(), last.atoms.sort_index()
        if not a.index.equals(b.index) or not np.array_equal(a['type'], b['type']):
            raise ValueError('NEB endpoints must preserve every atom ID and atom type.')
        if not np.allclose(first.box.matrix, last.box.matrix) or not np.allclose(first.box.bounds, last.box.bounds):
            raise ValueError('NEB endpoints must have identical cell geometry and origin.')
        if not np.isfinite(a[['x','y','z']].to_numpy()).all() or not np.isfinite(b[['x','y','z']].to_numpy()).all():
            raise ValueError('NEB endpoints contain nonfinite coordinates.')
        if np.allclose(a[['x','y','z']], b[['x','y','z']]):
            raise ValueError('NEB endpoints are identical; supply a transition between distinct states.')
        case.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / spec['initial'], case / 'initial.data')
        with (case / 'coords.final').open('w', encoding='utf-8') as handle:
            handle.write(f'{len(b)}\n')
            for atom_id, row in b.iterrows():
                handle.write(f'{int(atom_id)} {row.x:.16g} {row.y:.16g} {row.z:.16g}\n')
        _write_potential_include(case, cfg.potential.pot_lines, cfg.potential.masses)
        _copy_potential_files(root, cfg.potential.pot_name, case, cfg.potential.potential_files)
        text = f'''units metal
atom_style atomic
atom_modify map array sort 0 0.0
boundary p p p
read_data initial.data
include potential.inc
timestep 0.001
min_style fire
fix path all neb 1.0
thermo 100
neb 0.0 {wf.neb_force_tolerance:g} {wf.neb_iterations} {wf.neb_iterations} 100 final coords.final
variable replica equal part
write_data image.${{replica}}.data
print "NEB_COMPLETE"
'''
        hashes = {name: hashlib.sha256((case / name).read_bytes()).hexdigest()
                  for name in ('initial.data', 'coords.final', 'potential.inc', *cfg.potential.potential_files) if (case / name).is_file()}
        case_id = hashlib.sha256((text + json.dumps(hashes, sort_keys=True)).encode()).hexdigest()
        text = f'print "NEB_CASE {case_id}"\n' + text
        (case / 'in.neb.lammps').write_text(text, encoding='utf-8')
        hashes['in.neb.lammps'] = hashlib.sha256((case / 'in.neb.lammps').read_bytes()).hexdigest()
        (case / 'neb.json').write_text(json.dumps({'system': system, 'name': name, 'replicas': wf.neb_replicas,
            'force_tolerance': wf.neb_force_tolerance, 'case_id': case_id, 'file_hashes': hashes}, indent=2))
        manifest.append({'system': system, 'dir': str(case.relative_to(root))})
    return manifest


def analyze_neb_case(case):
    case = Path(case)
    meta = json.loads((case / 'neb.json').read_text())
    result = {'state': 'pending', 'system': meta['system'], 'name': meta['name']}
    if not (case / 'log.lammps').exists():
        return result
    text = (case / 'log.lammps').read_text(errors='replace')
    if 'ERROR' in text:
        return {**result, 'state': 'failed'}
    if 'Total wall time' not in text:
        return {**result, 'state': 'running'}
    if 'NEB_COMPLETE' not in text or f"NEB_CASE {meta['case_id']}" not in text:
        return {**result, 'state': 'failed'}
    for name, digest in meta['file_hashes'].items():
        if not (case / name).is_file() or hashlib.sha256((case / name).read_bytes()).hexdigest() != digest:
            return {**result, 'state': 'failed'}
    if not all((case / f'image.{i}.data').exists() for i in range(meta['replicas'])):
        return {**result, 'state': 'failed'}
    rows = []
    for line in text.splitlines():
        tokens = line.split()
        if len(tokens) != 9 + 2 * meta['replicas']:
            continue
        try:
            row = np.asarray([float(x) for x in tokens])
        except ValueError:
            continue
        if np.isfinite(row).all():
            rows.append(row)
    if not rows:
        return {**result, 'state': 'failed'}
    final = rows[-1]
    # ftol is the maximum replica 3N force norm, not max force component.
    if final[1] > meta['force_tolerance']:
        return {**result, 'state': 'unconverged', 'max_replica_force': float(final[1])}
    energies = final[10::2]
    if len(energies) != meta['replicas']:
        return {**result, 'state': 'failed'}
    return {**result, 'state': 'done', 'forward_barrier': float(max(energies) - energies[0]),
            'backward_barrier': float(max(energies) - energies[-1]), 'energies': energies.tolist(),
            'reaction_coordinates': final[9::2].tolist(), 'max_replica_force': float(final[1])}
