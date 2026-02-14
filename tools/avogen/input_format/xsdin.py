#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET

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
    meta_data = {}
    meta_data['outputFormat'] = 'cjson'
    meta_data['operations'] = ['read']
    meta_data['identifier'] = 'Materials Studio XSD Format'
    meta_data['name'] = 'Materials Studio XSD'
    meta_data['description'] = 'Read Materials Studio XSD and preserve cell/symmetry as cjson.'
    meta_data['fileExtensions'] = ['xsd']
    meta_data['mimeTypes'] = ['chemical/x-materials-studio-xsd']
    return meta_data


def _clean_symbol(token):
    if not token:
        return ''
    token = token.strip()
    if not token:
        return ''
    return token[0].upper() + token[1:].lower()


def _parse_csv_floats(text):
    if not text:
        return []
    vals = []
    for part in text.split(','):
        s = part.strip()
        if not s:
            continue
        vals.append(float(s.replace('D', 'E').replace('d', 'e')))
    return vals


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


def _frac_to_cart(cell, frac):
    a = cell[0]
    b = cell[1]
    c = cell[2]
    return [
        frac[0] * a[0] + frac[1] * b[0] + frac[2] * c[0],
        frac[0] * a[1] + frac[1] * b[1] + frac[2] * c[1],
        frac[0] * a[2] + frac[1] * b[2] + frac[2] * c[2],
    ]


def _parse_operators(raw):
    if not raw:
        return []
    ops = []
    chunks = [chunk.strip() for chunk in raw.split(':') if chunk.strip()]
    for chunk in chunks:
        vals = _parse_csv_floats(chunk)
        if len(vals) == 12:
            ops.append(vals)
    return ops


def _format_screw_axis(token):
    if len(token) == 2 and token[0] in ('2', '3', '4', '6') and token[1].isdigit():
        return token[0] + '_' + token[1]
    return token


def _normalize_space_group(group_name, long_name):
    # Prefer compact GroupName (e.g. "P42/MNM"), normalize into
    # a Hermann-Mauguin-like text used by Avogadro (e.g. "P 4_2/m n m").
    g = (group_name or '').strip()
    if g:
        g = g.replace(' ', '')
        if len(g) >= 2 and g[0].isalpha():
            lattice = g[0].upper()
            tail = g[1:]

            if '/' in tail:
                first, rest = tail.split('/', 1)
                first = _format_screw_axis(first)
                rest_tokens = [ch.lower() for ch in rest if ch.isalpha() or ch.isdigit()]
                if rest_tokens:
                    return '{} {}/{}'.format(lattice, first, rest_tokens[0]) + (
                        ' ' + ' '.join(rest_tokens[1:]) if len(rest_tokens) > 1 else ''
                    )
                return '{} {}'.format(lattice, first)

            # No slash present: just split lattice and body.
            return '{} {}'.format(lattice, _format_screw_axis(tail))

    # Fallback: keep long name if present.
    l = (long_name or '').strip()
    if l:
        return l
    return 'P1'


def _first_or_none(nodes):
    if not nodes:
        return None
    return nodes[0]


def _resolve_atom_symbol(atom_node, atom_by_id):
    comp = atom_node.get('Components', '')
    if comp:
        return _clean_symbol(comp)

    image_of = atom_node.get('ImageOf', '')
    if image_of and image_of in atom_by_id:
        ref = atom_by_id[image_of]
        return _clean_symbol(ref.get('Components', ''))

    name = atom_node.get('Name', '')
    letters = ''.join(ch for ch in name if ch.isalpha())
    return _clean_symbol(letters)


def _parse_xsd(text):
    root = ET.fromstring(text)

    space_group = _first_or_none(root.findall('.//SpaceGroup'))
    symm_system = _first_or_none(root.findall('.//SymmetrySystem'))

    vectors = None
    symmetry = {}
    if space_group is not None:
        av = _parse_csv_floats(space_group.get('AVector', ''))
        bv = _parse_csv_floats(space_group.get('BVector', ''))
        cv = _parse_csv_floats(space_group.get('CVector', ''))
        if len(av) == 3 and len(bv) == 3 and len(cv) == 3:
            vectors = [av, bv, cv]

        symmetry = {
            'name': space_group.get('Name', ''),
            'groupName': space_group.get('GroupName', ''),
            'longName': space_group.get('LongName', ''),
            'qualifier': space_group.get('Qualifier', ''),
            'itNumber': space_group.get('ITNumber', ''),
            'schoenfliesName': space_group.get('SchoenfliesName', ''),
            'crystalSystem': space_group.get('System', ''),
            'crystalClass': space_group.get('Class', ''),
            'lattice': space_group.get('Lattice', ''),
            'centering': space_group.get('Centering', ''),
            'operatorsRaw': space_group.get('Operators', ''),
            'operators': _parse_operators(space_group.get('Operators', '')),
        }

    atoms_parent = symm_system if symm_system is not None else root
    atom_nodes = atoms_parent.findall('.//Atom3d')

    atom_by_id = {}
    for atom in atom_nodes:
        atom_id = atom.get('ID', '')
        if atom_id:
            atom_by_id[atom_id] = atom

    mapping_by_id = {}
    for mapping_tag in ('IdentityMapping', 'ImageMapping'):
        for mapping in atoms_parent.findall('.//' + mapping_tag):
            mid = mapping.get('ID', '')
            if mid:
                mapping_by_id[mid] = mapping

    def mapping_translation(atom_node):
        mapping_id = atom_node.get('Mapping', '')
        if not mapping_id:
            return None
        mapping = mapping_by_id.get(mapping_id)
        if mapping is None:
            return None

        # Prefer the affine element translation; fallback to constraint translation.
        m12 = _parse_csv_floats(mapping.get('Element', ''))
        if len(m12) == 12:
            return [m12[3], m12[7], m12[11]]
        c12 = _parse_csv_floats(mapping.get('Constraint', ''))
        if len(c12) == 12:
            return [c12[3], c12[7], c12[11]]
        return None

    records = []
    for atom in atom_nodes:
        # Keep only asymmetric-unit atoms and rely on space-group metadata.
        if atom.get('ImageOf', ''):
            continue
        xyz = _parse_csv_floats(atom.get('XYZ', ''))
        if len(xyz) != 3:
            xyz = mapping_translation(atom) or []
        records.append({
            'id': atom.get('ID', ''),
            'name': atom.get('Name', ''),
            'imageOf': '',
            'symbol': _resolve_atom_symbol(atom, atom_by_id),
            'xyz': xyz if len(xyz) == 3 else None,
        })

    atoms = []
    for rec in records:
        frac = rec['xyz']
        if frac is None or not rec['symbol']:
            continue
        atoms.append({
            'symbol': rec['symbol'],
            'fractional': frac,
            'id': rec['id'],
            'name': rec['name'],
            'imageOf': rec['imageOf'],
        })

    atom_ids = set([atom['id'] for atom in atoms if atom['id']])
    id_to_index = {}
    for idx, atom in enumerate(atoms):
        if atom['id']:
            id_to_index[atom['id']] = idx

    bonds = []
    seen = set()
    bond_nodes = atoms_parent.findall('.//Bond')
    for bond in bond_nodes:
        connects = bond.get('Connects', '')
        if not connects:
            continue
        parts = [p.strip() for p in connects.split(',') if p.strip()]
        if len(parts) != 2:
            continue
        a_id, b_id = parts[0], parts[1]
        if a_id not in atom_ids or b_id not in atom_ids:
            continue
        ia = id_to_index[a_id]
        ib = id_to_index[b_id]
        if ia == ib:
            continue
        key = (ia, ib) if ia < ib else (ib, ia)
        if key in seen:
            continue
        seen.add(key)

        order_attr = bond.get('Order', '').strip()
        order = 1
        if order_attr:
            try:
                order = int(float(order_attr))
            except ValueError:
                order = 1
        bonds.append((key[0], key[1], order))

    return root, vectors, symmetry, atoms, bonds


def _build_cjson(root, vectors, symmetry, atoms, bonds):
    atomic_numbers = []
    coords_3d = []
    coords_3d_frac = []

    for atom in atoms:
        symbol = atom['symbol']
        frac = atom['fractional']
        cart = _frac_to_cart(vectors, frac) if vectors is not None else frac

        atomic_numbers.append(_ATOMIC_NUMBERS.get(symbol, 0))
        coords_3d.extend([cart[0], cart[1], cart[2]])
        coords_3d_frac.extend([frac[0], frac[1], frac[2]])

    cjson = {
        'chemicalJson': 1,
        'name': 'Generated from XSD',
        'atoms': {
            'elements': {'number': atomic_numbers},
            'coords': {'3d': coords_3d}
        },
        'properties': {
            'sourceFormat': 'xsd',
            'xsdVersion': root.get('Version', ''),
            'xsdWrittenBy': root.get('WrittenBy', ''),
        },
    }

    # In asymmetric-unit mode with space-group expansion, cjson bond indices
    # cannot reference symmetry-generated atoms in neighboring images.
    # Emitting only asymmetric bonds would show incomplete connectivity
    # (e.g., only one Ti-O in TiO2). Let Avogadro perceive bonds instead.
    emit_bonds = not (vectors is not None and (symmetry.get('groupName') or symmetry.get('name')))
    if emit_bonds and bonds:
        conn_index = []
        orders = []
        for ia, ib, order in bonds:
            conn_index.extend([ia, ib])
            orders.append(order)
        cjson['bonds'] = {
            'connections': {'index': conn_index},
            'order': orders,
        }

    if coords_3d_frac:
        cjson['atoms']['coords']['3dFractional'] = coords_3d_frac

    if vectors is not None:
        cjson['unitCell'] = {
            'spaceGroup': _normalize_space_group(
                symmetry.get('groupName', ''),
                symmetry.get('longName', ''),
            ),
            'a': _norm(vectors[0]),
            'b': _norm(vectors[1]),
            'c': _norm(vectors[2]),
            'alpha': _angle_deg(vectors[1], vectors[2]),
            'beta': _angle_deg(vectors[0], vectors[2]),
            'gamma': _angle_deg(vectors[0], vectors[1]),
            'cellVectors': [
                vectors[0][0], vectors[0][1], vectors[0][2],
                vectors[1][0], vectors[1][1], vectors[1][2],
                vectors[2][0], vectors[2][1], vectors[2][2],
            ],
            'itNumber': symmetry.get('itNumber', ''),
            'schoenfliesName': symmetry.get('schoenfliesName', ''),
            'lattice': symmetry.get('lattice', ''),
            'centering': symmetry.get('centering', ''),
            'crystalSystem': symmetry.get('crystalSystem', ''),
            'crystalClass': symmetry.get('crystalClass', ''),
            'symmetryOperators': symmetry.get('operators', []),
            'symmetryOperatorsRaw': symmetry.get('operatorsRaw', ''),
        }

    cjson['properties']['xsdSymmetry'] = {
        'name': symmetry.get('name', ''),
        'groupName': symmetry.get('groupName', ''),
        'longName': symmetry.get('longName', ''),
        'qualifier': symmetry.get('qualifier', ''),
    }

    atom_meta = []
    for atom in atoms:
        atom_meta.append({
            'id': atom['id'],
            'name': atom['name'],
            'imageOf': atom['imageOf'],
            'symbol': atom['symbol'],
        })
    cjson['properties']['xsdAtoms'] = atom_meta

    return json.dumps(cjson)


def read():
    text = sys.stdin.read()
    root, vectors, symmetry, atoms, bonds = _parse_xsd(text)
    return _build_cjson(root, vectors, symmetry, atoms, bonds)


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Materials Studio XSD file format script.')
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
