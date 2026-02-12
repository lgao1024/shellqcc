#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

import argparse
import json
import math
import re
import sys

BOHR_TO_ANGSTROM = 0.529177210903


def getMetaData():
    metaData = {}
    metaData['outputFormat'] = 'cjson'
    metaData['operations'] = ['read']
    metaData['identifier'] = 'QE Input Format'
    metaData['name'] = 'QE Input'
    metaData['description'] = "Read QE input and convert ATOMIC_POSITIONS/CELL_PARAMETERS to cjson."
    metaData['fileExtensions'] = ['in', 'pwi']
    metaData['mimeTypes'] = ['chemical/x-quantum-espresso-input']
    return metaData


def _parse_unit(header, keyword):
    line = header.strip().lower()
    m = re.search(r'\{\s*([^}]+)\s*\}', line)
    if m:
        return m.group(1).strip()

    tail = line[len(keyword):].strip()
    if tail:
        return tail.split()[0].strip()
    return ''


def _extract_system_value(text, key):
    m = re.search(r'(?im)^\s*' + re.escape(key) + r'\s*=\s*([^,\s!/]+)', text)
    if not m:
        return None
    token = m.group(1).strip()
    token = token.replace('d', 'e').replace('D', 'E')
    try:
        return float(token)
    except ValueError:
        return None


def _is_block_header(line):
    up = line.strip().upper()
    if not up:
        return True
    return (
        up.startswith('&')
        or up.startswith('/')
        or up.startswith('K_POINTS')
        or up.startswith('ATOMIC_SPECIES')
        or up.startswith('ATOMIC_POSITIONS')
        or up.startswith('CELL_PARAMETERS')
        or up.startswith('CONSTRAINTS')
        or up.startswith('OCCUPATIONS')
        or up.startswith('HUBBARD')
    )


def _parse_three_floats(parts):
    vals = []
    for i in range(3):
        token = parts[i].replace('d', 'e').replace('D', 'E')
        vals.append(float(token))
    return vals


def _clean_symbol(token):
    m = re.match(r'^([A-Za-z]{1,3})', token)
    if not m:
        return token
    sym = m.group(1)
    return sym[0].upper() + sym[1:].lower()


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


def _parse_qe(text):
    lines = text.splitlines()
    celldm1 = _extract_system_value(text, 'celldm(1)')
    alat = _extract_system_value(text, 'A')
    if alat is None and celldm1 is not None:
        alat = celldm1 * BOHR_TO_ANGSTROM

    cell = None
    cell_scale = 1.0
    atoms = []
    atom_pos_mode = 'crystal'

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        up = line.upper()

        if up.startswith('CELL_PARAMETERS'):
            unit = _parse_unit(line, 'cell_parameters').lower()
            if unit in ('angstrom',):
                cell_scale = 1.0
            elif unit in ('bohr', 'au'):
                cell_scale = BOHR_TO_ANGSTROM
            else:
                # default QE behavior is alat when not explicit
                cell_scale = alat if alat is not None else 1.0

            vectors = []
            for j in range(1, 4):
                if i + j >= len(lines):
                    break
                parts = lines[i + j].split()
                if len(parts) < 3:
                    break
                try:
                    v = _parse_three_floats(parts)
                except ValueError:
                    break
                vectors.append(v)

            if len(vectors) == 3:
                cell = [
                    [vectors[0][0] * cell_scale, vectors[0][1] * cell_scale, vectors[0][2] * cell_scale],
                    [vectors[1][0] * cell_scale, vectors[1][1] * cell_scale, vectors[1][2] * cell_scale],
                    [vectors[2][0] * cell_scale, vectors[2][1] * cell_scale, vectors[2][2] * cell_scale],
                ]
            i += 4
            continue

        if up.startswith('ATOMIC_POSITIONS'):
            unit = _parse_unit(line, 'atomic_positions').lower()
            if unit in ('crystal',):
                atom_pos_mode = 'crystal'
            elif unit in ('angstrom',):
                atom_pos_mode = 'angstrom'
            elif unit in ('bohr', 'au'):
                atom_pos_mode = 'bohr'
            else:
                atom_pos_mode = 'alat'

            i += 1
            while i < len(lines):
                raw = lines[i]
                s = raw.strip()
                if not s:
                    break
                if _is_block_header(raw):
                    break
                parts = s.split()
                if len(parts) >= 4:
                    try:
                        xyz = _parse_three_floats(parts[1:4])
                    except ValueError:
                        i += 1
                        continue
                    atoms.append((parts[0], xyz))
                i += 1
            continue

        i += 1

    return cell, atoms, atom_pos_mode, alat


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

    inv = [
        [(e * i - f * h) / det, (c * h - b * i) / det, (b * f - c * e) / det],
        [(f * g - d * i) / det, (a * i - c * g) / det, (c * d - a * f) / det],
        [(d * h - e * g) / det, (b * g - a * h) / det, (a * e - b * d) / det],
    ]
    return inv


def _mat_vec_mul(m, v):
    return [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]


def _to_cartesian_and_fractional(cell, atoms, atom_pos_mode, alat):
    if atom_pos_mode == 'crystal':
        frac_atoms = [(label, coords) for label, coords in atoms]
        if cell is None:
            return [], frac_atoms
        cart_atoms = []
        for label, frac in frac_atoms:
            cart_atoms.append((label, _frac_to_cart(cell, frac)))
        return cart_atoms, frac_atoms

    scale = 1.0
    if atom_pos_mode == 'angstrom':
        scale = 1.0
    elif atom_pos_mode == 'bohr':
        scale = BOHR_TO_ANGSTROM
    elif atom_pos_mode == 'alat':
        scale = alat if alat is not None else 1.0

    cart_atoms = []
    for label, coords in atoms:
        cart_atoms.append((label, [coords[0] * scale, coords[1] * scale, coords[2] * scale]))

    if cell is None:
        return cart_atoms, []

    inv = _mat_inv3(cell)
    if inv is None:
        return cart_atoms, []

    frac_atoms = []
    for label, cart in cart_atoms:
        frac_atoms.append((label, _mat_vec_mul(inv, cart)))
    return cart_atoms, frac_atoms


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


def _build_cjson(cell, atoms_cart, atoms_frac):
    atomic_numbers = []
    coords3d = []
    coords3d_frac = []

    for label, cart in atoms_cart:
        symbol = _clean_symbol(label)
        atomic_numbers.append(_ATOMIC_NUMBERS.get(symbol, 0))
        coords3d.extend([cart[0], cart[1], cart[2]])

    if atoms_frac:
        for _, frac in atoms_frac:
            coords3d_frac.extend([frac[0], frac[1], frac[2]])

    cjson = {
        'chemicalJson': 1,
        'name': 'Generated from QE input',
        'atoms': {
            'elements': {
                'number': atomic_numbers
            },
            'coords': {
                '3d': coords3d
            }
        }
    }

    if coords3d_frac:
        cjson['atoms']['coords']['3dFractional'] = coords3d_frac

    if cell is not None:
        a = _norm(cell[0])
        b = _norm(cell[1])
        c = _norm(cell[2])
        alpha = _angle_deg(cell[1], cell[2])
        beta = _angle_deg(cell[0], cell[2])
        gamma = _angle_deg(cell[0], cell[1])
        cjson['unitCell'] = {
            'a': a,
            'b': b,
            'c': c,
            'alpha': alpha,
            'beta': beta,
            'gamma': gamma,
            'cellVectors': [
                cell[0][0], cell[0][1], cell[0][2],
                cell[1][0], cell[1][1], cell[1][2],
                cell[2][0], cell[2][1], cell[2][2]
            ]
        }

    return json.dumps(cjson)


def read():
    text = sys.stdin.read()
    cell, atoms, atom_pos_mode, alat = _parse_qe(text)
    atoms_cart, atoms_frac = _to_cartesian_and_fractional(cell, atoms, atom_pos_mode, alat)
    return _build_cjson(cell, atoms_cart, atoms_frac)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('QE input file format script.')
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
