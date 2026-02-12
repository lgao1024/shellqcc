#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

import argparse
import json
import math
import re
import sys

BOHR_TO_ANGSTROM = 0.529177210903

_ELEMENTS = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
    'Es', 'Fm', 'Md', 'No', 'Lr',
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',
    'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'
]
_ATOMIC_NUMBERS = {symbol: idx + 1 for idx, symbol in enumerate(_ELEMENTS)}


def getMetaData():
    metaData = {}
    metaData['outputFormat'] = 'cjson'
    metaData['operations'] = ['read']
    metaData['identifier'] = 'CP2K Input Format'
    metaData['name'] = 'CP2K INP'
    metaData['description'] = 'Read CP2K input and convert CELL/COORD blocks to cjson.'
    metaData['fileExtensions'] = ['inp']
    metaData['mimeTypes'] = ['chemical/x-cp2k-input']
    return metaData


def _clean_line(raw):
    text = raw.split('!', 1)[0].split('#', 1)[0]
    return text.strip()


def _clean_symbol(token):
    m = re.match(r'^([A-Za-z]{1,3})', token)
    if not m:
        return token
    sym = m.group(1)
    return sym[0].upper() + sym[1:].lower()


def _parse_bracket_unit(line):
    m = re.search(r'\[\s*([^\]]+)\s*\]', line)
    if not m:
        return ''
    return m.group(1).strip().lower()


def _unit_scale(unit):
    u = unit.lower()
    if 'bohr' in u or u == 'au':
        return BOHR_TO_ANGSTROM
    return 1.0


def _extract_floats(text):
    nums = re.findall(r'[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?', text)
    out = []
    for n in nums:
        out.append(float(n.replace('d', 'e').replace('D', 'E')))
    return out


def _lattice_vectors_from_abc(a, b, c, alpha_deg, beta_deg, gamma_deg):
    alpha = math.radians(alpha_deg)
    beta = math.radians(beta_deg)
    gamma = math.radians(gamma_deg)

    va = [a, 0.0, 0.0]
    vb = [b * math.cos(gamma), b * math.sin(gamma), 0.0]

    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / max(math.sin(gamma), 1.0e-14)
    cz2 = c * c - cx * cx - cy * cy
    cz = math.sqrt(max(cz2, 0.0))
    vc = [cx, cy, cz]

    return [va, vb, vc]


def _parse_cp2k(text):
    lines = text.splitlines()
    in_cell = False
    in_coord = False

    cell_a = None
    cell_b = None
    cell_c = None
    abc = None
    abg = None

    coord_scaled = False
    coord_unit = 'angstrom'
    atoms_raw = []

    for raw in lines:
        line = _clean_line(raw)
        if not line:
            continue

        upper = line.upper()

        if upper.startswith('&CELL'):
            in_cell = True
            continue
        if upper.startswith('&COORD'):
            in_coord = True
            unit = _parse_bracket_unit(line)
            if unit:
                coord_unit = unit
            continue

        if upper.startswith('&END'):
            if in_cell:
                in_cell = False
            elif in_coord:
                in_coord = False
            continue

        if in_cell:
            key = upper.split()[0]
            unit = _parse_bracket_unit(line)
            scale = _unit_scale(unit)
            vals = _extract_floats(line)

            if key == 'A' and len(vals) >= 3:
                cell_a = [vals[0] * scale, vals[1] * scale, vals[2] * scale]
            elif key == 'B' and len(vals) >= 3:
                cell_b = [vals[0] * scale, vals[1] * scale, vals[2] * scale]
            elif key == 'C' and len(vals) >= 3:
                cell_c = [vals[0] * scale, vals[1] * scale, vals[2] * scale]
            elif key == 'ABC' and len(vals) >= 3:
                abc = [vals[0] * scale, vals[1] * scale, vals[2] * scale]
            elif key == 'ALPHA_BETA_GAMMA' and len(vals) >= 3:
                abg = [vals[0], vals[1], vals[2]]
            continue

        if in_coord:
            parts = line.split()
            if not parts:
                continue

            key = parts[0].upper()
            if key == 'UNIT' and len(parts) >= 2:
                coord_unit = parts[1].lower()
                continue
            if key == 'SCALED' and len(parts) >= 2:
                v = parts[1].upper()
                coord_scaled = v in ('T', 'TRUE', '.TRUE.', '1', 'YES', 'Y')
                continue

            symbol = _clean_symbol(parts[0])
            vals = _extract_floats(' '.join(parts[1:]))
            if len(vals) < 3:
                continue
            atoms_raw.append((symbol, [vals[0], vals[1], vals[2]]))

    cell = None
    if cell_a is not None and cell_b is not None and cell_c is not None:
        cell = [cell_a, cell_b, cell_c]
    elif abc is not None:
        if abg is None:
            abg = [90.0, 90.0, 90.0]
        cell = _lattice_vectors_from_abc(abc[0], abc[1], abc[2], abg[0], abg[1], abg[2])

    return cell, atoms_raw, coord_scaled, coord_unit


def _frac_to_cart(cell, frac):
    a = cell[0]
    b = cell[1]
    c = cell[2]
    return [
        frac[0] * a[0] + frac[1] * b[0] + frac[2] * c[0],
        frac[0] * a[1] + frac[1] * b[1] + frac[2] * c[1],
        frac[0] * a[2] + frac[1] * b[2] + frac[2] * c[2],
    ]


def _mat_inv3(m):
    a, b, c = m[0]
    d, e, f = m[1]
    g, h, i = m[2]
    det = (
        a * (e * i - f * h)
        - b * (d * i - f * g)
        + c * (d * h - e * g)
    )
    if abs(det) < 1.0e-14:
        return None

    return [
        [(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
        [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
        [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det],
    ]


def _mat_vec_mul(m, v):
    return [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]


def _norm(v):
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def _angle_deg(v1, v2):
    n1 = _norm(v1)
    n2 = _norm(v2)
    if n1 < 1.0e-14 or n2 < 1.0e-14:
        return 90.0
    c = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2)
    c = max(-1.0, min(1.0, c))
    return math.degrees(math.acos(c))


def _convert_atoms(cell, atoms_raw, coord_scaled, coord_unit):
    atoms_cart = []
    atoms_frac = []

    if coord_scaled and cell is not None:
        for symbol, frac in atoms_raw:
            atoms_frac.append((symbol, frac))
            atoms_cart.append((symbol, _frac_to_cart(cell, frac)))
        return atoms_cart, atoms_frac

    scale = _unit_scale(coord_unit)
    for symbol, coords in atoms_raw:
        atoms_cart.append((symbol, [coords[0] * scale, coords[1] * scale, coords[2] * scale]))

    if cell is None:
        return atoms_cart, atoms_frac

    inv = _mat_inv3(cell)
    if inv is None:
        return atoms_cart, atoms_frac

    for symbol, cart in atoms_cart:
        atoms_frac.append((symbol, _mat_vec_mul(inv, cart)))
    return atoms_cart, atoms_frac


def _build_cjson(cell, atoms_cart, atoms_frac):
    numbers = []
    coords3d = []
    coords3d_frac = []

    for symbol, cart in atoms_cart:
        numbers.append(_ATOMIC_NUMBERS.get(_clean_symbol(symbol), 0))
        coords3d.extend([cart[0], cart[1], cart[2]])

    for _, frac in atoms_frac:
        coords3d_frac.extend([frac[0], frac[1], frac[2]])

    cjson = {
        'chemicalJson': 1,
        'name': 'Generated from CP2K input',
        'atoms': {
            'elements': {'number': numbers},
            'coords': {'3d': coords3d}
        }
    }

    if coords3d_frac:
        cjson['atoms']['coords']['3dFractional'] = coords3d_frac

    if cell is not None:
        cjson['unitCell'] = {
            'spaceGroup': 'P1',
            'a': _norm(cell[0]),
            'b': _norm(cell[1]),
            'c': _norm(cell[2]),
            'alpha': _angle_deg(cell[1], cell[2]),
            'beta': _angle_deg(cell[0], cell[2]),
            'gamma': _angle_deg(cell[0], cell[1]),
            'cellVectors': [
                cell[0][0], cell[0][1], cell[0][2],
                cell[1][0], cell[1][1], cell[1][2],
                cell[2][0], cell[2][1], cell[2][2],
            ]
        }

    return json.dumps(cjson)


def read():
    text = sys.stdin.read()
    cell, atoms_raw, coord_scaled, coord_unit = _parse_cp2k(text)
    atoms_cart, atoms_frac = _convert_atoms(cell, atoms_raw, coord_scaled, coord_unit)
    return _build_cjson(cell, atoms_cart, atoms_frac)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('CP2K INP file format script.')
    parser.add_argument('--metadata', action='store_true')
    parser.add_argument('--read', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    if args['metadata']:
        print(json.dumps(getMetaData()))
    elif args['display_name']:
        print(getMetaData()['name'])
    elif args['read']:
        print(read())
