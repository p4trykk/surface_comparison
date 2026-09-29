#!/usr/bin/env python3
"""
how to use: 
python align_nylon.py --reference nylon_aligned_for_reference.pdb --inputs input_nylon_to_align.pdb --output-dir aligned_nylon_output_dir

outputs:
- aligned PDB files in specified output directory
- JSON report files in specified output directory
"""


import argparse
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
import numpy as np

DEFAULT_NAMES = ('N', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O')

@dataclass
class Structure:
    path: Path
    lines: list
    atom_lines: list
    xyz: np.ndarray
    residues: dict

def read_pdb(path):
    path = Path(path)
    lines = path.read_text().splitlines(keepends=True)
    atoms, xyz, residues = [], [], {}
    if sum(line[:6].strip() == 'MODEL' for line in lines) > 1:
        raise ValueError(f'{path}: multiple MODELs are unsupported; select one first.')
    for line_number, line in enumerate(lines):
        record = line[:6].strip()
        if record in {'CRYST1', 'SCALE1', 'SCALE2', 'SCALE3', 'ORIGX1', 'ORIGX2',
                      'ORIGX3', 'MTRIX1', 'MTRIX2', 'MTRIX3'}:
            raise ValueError(f'{path}: {record} is frame-dependent; prepare a finite-block PDB without these records first.')
        if record not in {'ATOM', 'HETATM'}:
            continue
        if len(line.rstrip('\r\n')) < 54:
            raise ValueError(f'{path}:{line_number+1}: truncated coordinate record.')
        if line[16] != ' ':
            raise ValueError(f'{path}:{line_number+1}: select one alternate conformer first.')
        try:
            point = [float(line[a:a+8]) for a in (30, 38, 46)]
            resid = int(line[22:26])
        except ValueError as exc:
            raise ValueError(f'{path}:{line_number+1}: invalid XYZ or residue number.') from exc
        key = (line[21], resid, line[26])
        name = line[12:16].strip()
        residue = residues.setdefault(key, {'resname': line[17:20].strip(), 'atoms': {}})
        if residue['resname'] != line[17:20].strip() or name in residue['atoms']:
            raise ValueError(f'{path}: ambiguous residue/duplicate atom {key}, {name}; use unique chain and residue IDs.')
        residue['atoms'][name] = len(xyz)
        atoms.append(line_number)
        xyz.append(point)
    array = np.asarray(xyz, dtype=float)
    if len(array) < 3 or not np.isfinite(array).all():
        raise ValueError(f'{path}: insufficient atoms or nonfinite coordinates.')
    return Structure(path, lines, atoms, array, residues)


def parse_residue(text):
    match = re.fullmatch(r'([^:]):(-?\d+)([A-Za-z]?)', text)
    chain, number, insertion = match.groups()
    return (' ' if chain == '_' else chain, int(number), insertion or ' ')


def residue_label(key):
    return f"{key[0] if key[0] != ' ' else '_'}:{key[1]}{key[2].strip()}"


def anchor_points(structure, key, names, resname):
    residue = structure.residues.get(key)
    if residue is None:
        raise ValueError(f'{structure.path}: anchor {residue_label(key)} not found; specify --anchor explicitly.')
    if residue['resname'] != resname:
        raise ValueError(f'{structure.path}: anchor {residue_label(key)} is {residue["resname"]}, expected {resname}.')
    missing = set(names) - set(residue['atoms'])
    if missing:
        raise ValueError(f'{structure.path}: anchor {residue_label(key)} lacks atoms {sorted(missing)}.')
    return structure.xyz[[residue['atoms'][name] for name in names]]


def fit_orientation(mobile, reference, pairs, names, resname, rmsd_limit):
    source, target = [], []
    for mkey, rkey in pairs:
        a = anchor_points(mobile, mkey, names, resname)
        b = anchor_points(reference, rkey, names, resname)
        source.append(a - a.mean(axis=0))
        target.append(b - b.mean(axis=0))
    a, b = np.vstack(source), np.vstack(target)
    for label, arr in [('mobile', a), ('reference', b)]:
        sv = np.linalg.svd(arr, compute_uv=False)
        if sv[0] <= 1e-12 or sv[1] / sv[0] < 1e-3:
            raise ValueError(f'{label}: anchor geometry is almost collinear; choose a wider chemical motif.')
    u, singular, vt = np.linalg.svd(a.T @ b)
    correction = np.eye(3)
    correction[2, 2] = 1. if np.linalg.det(u @ vt) > 0 else -1.
    rotation = u @ correction @ vt  
    if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-10) or not np.isclose(np.linalg.det(rotation), 1.):
        raise ValueError('A proper rigid rotation could not be determined.')
    errors = []
    for pair, aa, bb in zip(pairs, source, target):
        rmsd = float(np.sqrt(np.mean(np.sum((aa @ rotation - bb)**2, axis=1))))
        errors.append({'mobile': residue_label(pair[0]), 'reference': residue_label(pair[1]),
                       'centred_motif_RMSD_A': rmsd})
        if rmsd > rmsd_limit:
            raise ValueError(f'Anchor {errors[-1]["mobile"]}: centred motif RMSD {rmsd:.4f} Å exceeds {rmsd_limit:g} Å. '
                             'Check anchor equivalence/conformation; do not increase tolerance without inspection.')
    return rotation, errors, singular


def transform_lines(structure, rotation, shift):
    out = structure.lines.copy()
    transformed = structure.xyz @ rotation + shift
    for line_index, point in zip(structure.atom_lines, transformed):
        fields = [f'{v:8.3f}' for v in point]
        if any(len(field) != 8 for field in fields):
            raise ValueError('Transformed coordinates exceed fixed-width PDB limits.')
        line = out[line_index]
        out[line_index] = line[:30] + ''.join(fields) + line[54:]
    for i, line in enumerate(out):
        if line[:6].strip() == 'ANISOU':
            try:
                q = [int(line[j:j+7]) for j in (28, 35, 42, 49, 56, 63)]
            except ValueError as exc:
                raise ValueError('Malformed ANISOU record.') from exc
            tensor = np.array([[q[0], q[3], q[4]], [q[3], q[1], q[5]], [q[4], q[5], q[2]]], float)
            rotated = rotation.T @ tensor @ rotation
            values = [rotated[0, 0], rotated[1, 1], rotated[2, 2], rotated[0, 1], rotated[0, 2], rotated[1, 2]]
            fields = [f'{int(np.rint(v)):7d}' for v in values]
            if any(len(f) != 7 for f in fields):
                raise ValueError('Rotated ANISOU exceeds PDB field limits.')
            out[i] = line[:28] + ''.join(fields) + line[70:]
    return out, transformed


def align_one(mobile, reference, pairs, names, resname='NYL', center='reference', rmsd_limit=.15):
    rotation, errors, singular = fit_orientation(mobile, reference, pairs, names, resname, rmsd_limit)
    source_center = mobile.xyz.mean(axis=0)
    reference_center = reference.xyz.mean(axis=0)
    target_center = reference_center if center == 'reference' else np.zeros(3)
    shift = target_center - source_center @ rotation
    lines, aligned = transform_lines(mobile, rotation, shift)
    bounds = np.array([aligned.min(axis=0), aligned.max(axis=0)])
    ray_anchor = np.array([bounds[0, 0], target_center[1], target_center[2]])
    report = {
        'mobile': str(mobile.path), 'reference': str(reference.path),
        'mobile_sha256': hashlib.sha256(mobile.path.read_bytes()).hexdigest(),
        'reference_sha256': hashlib.sha256(reference.path.read_bytes()).hexdigest(),
        'atom_count': len(mobile.xyz), 'reference_atom_count': len(reference.xyz),
        'method': 'Proper Kabsch rotation of chemically matched, individually centred repeat units',
        'anchor_atom_names': list(names), 'anchors': errors,
        'transform_equation': 'xyz_aligned = xyz_input @ rotation_row + translation_A',
        'rotation_row': rotation.tolist(), 'translation_A': shift.tolist(),
        'rotation_determinant': float(np.linalg.det(rotation)),
        'source_geometric_centroid_A': source_center.tolist(),
        'target_geometric_centroid_A': target_center.tolist(),
        'center_mode': center, 'centering_atom_selection': 'all ATOM and HETATM records',
        'aligned_bounds_A': bounds.tolist(), 'aligned_extent_A': np.ptp(aligned, axis=0).tolist(),
        'audited_SURFMAP_central_direction_aligned': [-1., 0., 0.],
        'equivalent_direction_in_input_frame': (np.array([-1., 0., 0.]) @ rotation.T).tolist(),
        'negative_X_atomic_support_anchor_A': ray_anchor.tolist(),
        }
    return lines, report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--reference', type=Path, required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--inputs', type=Path, nargs='+')
    group.add_argument('--input-dir', type=Path)
    parser.add_argument('--pattern', default='*.pdb', help='Nonrecursive pattern for --input-dir.')
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--anchor', action='append', help='MOBILE=REFERENCE residue pair; repeat for several anchors. Default A:1=A:1.')
    parser.add_argument('--atoms', default=','.join(DEFAULT_NAMES), help='Comma-separated atom names common to every anchor.')
    parser.add_argument('--resname', default='NYL')
    parser.add_argument('--center', choices=['reference', 'origin'], default='reference')
    parser.add_argument('--max-anchor-rmsd', type=float, default=.15)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args(argv)
    if not np.isfinite(args.max_anchor_rmsd) or args.max_anchor_rmsd <= 0:
        parser.error('--max-anchor-rmsd must be positive and finite.')
    names = tuple(n.strip() for n in args.atoms.split(',') if n.strip())
    if len(names) < 3 or len(set(names)) != len(names):
        parser.error('--atoms must contain at least three distinct names.')
    try:
        pairs = []
        for spec in args.anchor or ['A:1=A:1']:
            if spec.count('=') != 1:
                raise ValueError('Each --anchor must have MOBILE=REFERENCE form.')
            a, b = spec.split('=')
            pairs.append((parse_residue(a), parse_residue(b)))
        if len({a for a, b in pairs}) != len(pairs) or len({b for a, b in pairs}) != len(pairs):
            raise ValueError('Anchor mapping must be one-to-one without repeated residues.')
        reference = read_pdb(args.reference)
        for _, refkey in pairs:
            anchor_points(reference, refkey, names, args.resname)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    files = sorted(args.input_dir.glob(args.pattern)) if args.input_dir else args.inputs
    files = list(dict.fromkeys(p.resolve() for p in files if p.is_file()))
    if args.inputs and len(files) != len(set(p.resolve() for p in args.inputs)):
        parser.error('One or more explicit input paths are missing or not regular files.')
    if not files:
        parser.error('No input PDB files were found.')
    if len({p.stem for p in files}) != len(files):
        parser.error('Input filenames have duplicate stems; use separate output folders.')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    protected = set(files) | {args.reference.resolve()}
    print('Reference geometric centroid (Å):', reference.xyz.mean(axis=0))
    print('Anchors:', ', '.join(f'{residue_label(a)}={residue_label(b)}' for a, b in pairs))
    print('Use equivalent repeat-unit positions/phases; local RMSD does not resolve crystal symmetry.')
    failures = 0
    for path in files:
        output = args.output_dir / f'{path.stem}_aligned.pdb'
        audit = args.output_dir / f'{path.stem}_alignment.json'
        try:
            if output.resolve() in protected or audit.resolve() in protected:
                raise ValueError('Output would overwrite an input/reference; choose a different output folder.')
            if not args.overwrite and (output.exists() or audit.exists()):
                raise FileExistsError(f'{output.name} or its report already exists; use --overwrite deliberately.')
            mobile = read_pdb(path)
            lines, report = align_one(mobile, reference, pairs, names, args.resname, args.center, args.max_anchor_rmsd)
            # Validate formatting/reconstruction before writing either deliverable.
            coords = np.array([[float(lines[k][j:j+8]) for j in (30, 38, 46)] for k in mobile.atom_lines])
            expected = mobile.xyz @ np.array(report['rotation_row']) + report['translation_A']
            if np.max(np.abs(coords-expected)) > .000501:
                raise ValueError('PDB coordinate rounding check failed.')
            report['written_centroid_A'] = coords.mean(axis=0).tolist()
            output.write_text(''.join(lines))
            audit.write_text(json.dumps(report, indent=2) + '\n')
            maximum = max(row['centred_motif_RMSD_A'] for row in report['anchors'])
            print(f'OK {path.name}: {len(coords)} atoms; max anchor RMSD {maximum:.6f} Å -> {output}')
        except (OSError, ValueError, np.linalg.LinAlgError) as exc:
            failures += 1
            print(f'FAILED {path.name}: {exc}', file=sys.stderr)
    return 1 if failures else 0

if __name__ == '__main__':
    raise SystemExit(main())
