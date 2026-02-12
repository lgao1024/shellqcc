"""
/*****************************************************************************

  This source file is not yet a part of the Avogadro project.

Copyright 2015 Tomislav Subic <tomislav.subic@gmail.com>
Copyright 2026 Lei Gao <lgao1024@outlook.com>

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

******************************************************************************/
"""

import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET


TARGET_NAME = "CP2K2025"
DEBUG = False
FRACTIONAL_TOL = 1.0e-6

# Number of valence electrons for selected chemical elements.
VALENCE_ELECTRONS = {
  "H": 1,
  "He": 2,
  "Li": 3,
  "Be": 4,
  "B": 3,
  "C": 4,
  "N": 5,
  "O": 6,
  "F": 7,
  "Ne": 8,
  "Na": 9,
  "Mg": 10,
  "Al": 3,
  "Si": 4,
  "P": 5,
  "S": 6,
  "Cl": 7,
  "Ar": 8,
  "K": 9,
  "Ca": 10,
  "Sc": 11,
  "Ti": 12,
  "V": 13,
  "Cr": 14,
  "Mn": 15,
  "Fe": 16,
  "Co": 17,
  "Ni": 18,
  "Cu": 11,
  "Zn": 12,
  "Ga": 13,
  "Ge": 4,
  "As": 5,
  "Se": 6,
  "Br": 7,
  "Kr": 8,
  "Sr": 10,
  "Y": 11,
  "Zr": 12,
  "Mo": 14,
  "Ru": 16,
  "Rh": 17,
  "Pd": 18,
  "Ag": 11,
  "In": 13,
  "Sb": 5,
  "Te": 6,
  "I": 7,
  "Ba": 10,
  "W": 14,
  "Au": 11,
  "Bi": 15,
}


def getOptions():
  user_options = {
    "tabName": "Basic",
    "Title": {
      "type": "string",
      "default": "",
      "toolTip": "Title of the input file",
    },
    "Filename Base": {
      "type": "string",
      "default": "job",
    },
    "Charge": {
      "type": "integer",
      "default": 0,
      "minimum": -9,
      "maximum": 9,
      "toolTip": "Total charge of the system",
    },
    "Multiplicity": {
      "type": "integer",
      "default": 1,
      "minimum": 1,
      "maximum": 6,
      "toolTip": "Total spin multiplicity of the system",
    },
    "Grid number": {
      "type": "integer",
      "default": 4,
      "minimum": 1,
    },
    "Grid cutoff": {
      "type": "integer",
      "default": 200,
      "minimum": 1,
      "maximum": 999,
    },
    "Grid rel cutoff": {
      "type": "integer",
      "default": 40,
      "minimum": 1,
      "maximum": 999,
    },
    "Run Type": {
      "type": "stringList",
      "default": 0,
      "values": [
        "Energy",
        "Energy and forces",
        "Geometry Optimization",
        "Molecular dynamics",
      ],
      "toolTip": "Type of calculation to perform",
    },
    "Method": {
      "type": "stringList",
      "default": 0,
      "values": [
        "Electronic structure methods (DFT)",
        "Molecular Mechanics",
        "Hybrid quantum classical (Not yet supported)",
      ],
    },
    "Basis Set": {
      "type": "stringList",
      "default": 2,
      "values": [
        "minix",
        "SZV-MOLOPT-SR-GTH",
        "DZVP-MOLOPT-SR-GTH",
        "TZVP-MOLOPT-SR-GTH",
        "TZV2P-MOLOPT-SR-GTH",
      ],
    },
    "Functional": {
      "type": "stringList",
      "default": 0,
      "values": [
        "PBE",
        "BLYP",
        "BP",
        "HCTH120",
        "PADE",
      ],
    },
  }

  return {
    "userOptions": user_options,
    "inputMoleculeFormat": "cml",
  }


def parse_cml(cml):
  return ET.fromstring(cml)


def parse_atoms(cml, vectors=None):
  root = parse_cml(cml)
  atom_array = root.find(".//{*}atomArray")
  atoms = []
  if atom_array is None:
    return atoms

  for atom in atom_array.findall("{*}atom"):
    element = atom.get("elementType")
    x = atom.get("x3")
    y = atom.get("y3")
    z = atom.get("z3")

    if element is None:
      continue

    if x is not None and y is not None and z is not None:
      atoms.append((element, float(x), float(y), float(z)))
      continue

    fx = atom.get("xFract")
    fy = atom.get("yFract")
    fz = atom.get("zFract")
    if fx is not None and fy is not None and fz is not None and vectors is not None:
      cart = _fractional_to_cartesian([float(fx), float(fy), float(fz)], vectors)
      atoms.append((element, cart[0], cart[1], cart[2]))

  return atoms


def generateElements(cml, unique=0):
  root = parse_cml(cml)
  atom_array = root.find(".//{*}atomArray")
  elements = []
  if atom_array is None:
    return set() if unique == 1 else elements

  for atom in atom_array.findall("{*}atom"):
    element = atom.get("elementType")
    if element:
      elements.append(element)

  return set(elements) if unique == 1 else elements


def _read_crystal_scalars(cml):
  root = parse_cml(cml)
  crystal = root.find(".//{*}crystal")
  if crystal is None:
    return None

  scalar_values = {
    "a": None,
    "b": None,
    "c": None,
    "alpha": None,
    "beta": None,
    "gamma": None,
  }

  for scalar in crystal.findall("{*}scalar"):
    title = scalar.get("title")
    if title not in scalar_values or scalar.text is None:
      continue

    value = float(scalar.text)
    if title in ("alpha", "beta", "gamma"):
      value = math.radians(value)
    scalar_values[title] = value

  if any(value is None for value in scalar_values.values()):
    return None

  return scalar_values


def _lattice_vectors_from_scalars(scalars):
  a = scalars["a"]
  b = scalars["b"]
  c = scalars["c"]
  alpha = scalars["alpha"]
  beta = scalars["beta"]
  gamma = scalars["gamma"]

  vector_a = [a, 0.0, 0.0]

  bx = b * math.cos(gamma)
  by = b * math.sin(gamma)
  vector_b = [bx, by, 0.0]

  cx = c * math.cos(beta)
  cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
  cz = math.sqrt(c**2 - cx**2 - cy**2)
  vector_c = [cx, cy, cz]

  return vector_a, vector_b, vector_c


def generateCellVectors(cml):
  scalars = _read_crystal_scalars(cml)
  if scalars is None:
    return None
  return _lattice_vectors_from_scalars(scalars)


def generateMolecularCell(cml, vacuum_padding=10.0, min_size=12.0):
  root = parse_cml(cml)
  atom_array = root.find(".//{*}atomArray")
  fallback = ([20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0])

  if atom_array is None:
    return fallback

  coords = []
  for atom in atom_array.findall("{*}atom"):
    x = atom.get("x3")
    y = atom.get("y3")
    z = atom.get("z3")
    if x is None or y is None or z is None:
      continue
    coords.append((float(x), float(y), float(z)))

  if not coords:
    return fallback

  xs, ys, zs = zip(*coords)
  lx = max(max(xs) - min(xs) + vacuum_padding, min_size)
  ly = max(max(ys) - min(ys) + vacuum_padding, min_size)
  lz = max(max(zs) - min(zs) + vacuum_padding, min_size)

  return ([lx, 0.0, 0.0], [0.0, ly, 0.0], [0.0, 0.0, lz])


def format_vector(vector):
  return " ".join(f"{value:12.7f}" for value in vector)


def format_atom_line(atom):
  element, x, y, z = atom
  return f"      {element:<2s} {x:16.10f} {y:16.10f} {z:16.10f}"


def _cross(v1, v2):
  return [
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0],
  ]


def _dot(v1, v2):
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


def _cartesian_to_fractional(cart, vectors, tol):
  vector_a, vector_b, vector_c = vectors
  b_cross_c = _cross(vector_b, vector_c)
  c_cross_a = _cross(vector_c, vector_a)
  a_cross_b = _cross(vector_a, vector_b)
  volume = _dot(vector_a, b_cross_c)
  if abs(volume) < tol:
    raise ValueError("Invalid cell vectors: near-zero cell volume.")

  fx = _dot(cart, b_cross_c) / volume
  fy = _dot(cart, c_cross_a) / volume
  fz = _dot(cart, a_cross_b) / volume
  return [fx, fy, fz]


def _fractional_to_cartesian(frac, vectors):
  vector_a, vector_b, vector_c = vectors
  fx, fy, fz = frac
  return [
    fx * vector_a[0] + fy * vector_b[0] + fz * vector_c[0],
    fx * vector_a[1] + fy * vector_b[1] + fz * vector_c[1],
    fx * vector_a[2] + fy * vector_b[2] + fz * vector_c[2],
  ]


def _wrap_fractional(value, tol):
  wrapped = value - math.floor(value)
  if wrapped < tol or wrapped > 1.0 - tol:
    return 0.0
  return wrapped


def deduplicate_periodic_atoms(atoms, vectors, tol=FRACTIONAL_TOL):
  unique_atoms = []
  seen = set()

  for element, x, y, z in atoms:
    frac = _cartesian_to_fractional([x, y, z], vectors, tol)
    wrapped = [_wrap_fractional(value, tol) for value in frac]
    key = (
      element,
      int(round(wrapped[0] / tol)),
      int(round(wrapped[1] / tol)),
      int(round(wrapped[2] / tol)),
    )
    if key in seen:
      continue
    seen.add(key)
    cart = _fractional_to_cartesian(wrapped, vectors)
    unique_atoms.append((element, cart[0], cart[1], cart[2]))

  return unique_atoms


def _append_global(lines, title, calculate):
  run_type_map = {
    "Energy": "ENERGY",
    "Energy and forces": "ENERGY_FORCE",
    "Molecular dynamics": "MOLECULAR_DYNAMICS",
    "Geometry Optimization": "GEO_OPT",
  }

  try:
    run_type = run_type_map[calculate]
  except KeyError as exc:
    raise ValueError(f"Invalid calculation type: {calculate}") from exc

  lines.extend([
    "&GLOBAL",
    f"  PROJECT {title} ",
    f"  RUN_TYPE {run_type}",
    "  PRINT_LEVEL LOW",
    "  ELPA_KERNEL GENERIC_SIMPLE",
    "&END GLOBAL",
    "",
  ])


def _append_dft_block(lines, opts, functional, is_periodic_system):
  lines.extend([
    "  &DFT",
    "    BASIS_SET_FILE_NAME  BASIS_MOLOPT",
    "    POTENTIAL_FILE_NAME  POTENTIAL",
    f"    CHARGE {opts['Charge']}",
    f"    MULTIPLICITY {opts['Multiplicity']}",
  ])

  if not is_periodic_system:
    lines.extend([
      "    &POISSON",
      "      PERIODIC NONE",
      "      PSOLVER WAVELET",
      "    &END POISSON",
    ])

  lines.extend([
    "    &QS",
    "      EPS_DEFAULT 1.0E-10",
    "    &END QS",
    "    &MGRID",
    f"      NGRIDS {opts['Grid number']}",
    f"      CUTOFF {opts['Grid cutoff']}",
    f"      REL_CUTOFF {opts['Grid rel cutoff']}",
    "    &END MGRID",
    "    &XC",
    f"      &XC_FUNCTIONAL {functional}",
    "      &END XC_FUNCTIONAL",
    "    &END XC",
    "    &SCF",
    "      EPS_SCF 5.0E-06",
    "      MAX_SCF 200",
    "      SCF_GUESS ATOMIC",
    "      &DIAGONALIZATION",
    "        ALGORITHM STANDARD",
    "      &END DIAGONALIZATION",
    "      &MIXING",
    "        METHOD BROYDEN_MIXING",
    "        ALPHA 0.4",
    "        NBROYDEN 8",
    "      &END MIXING",
    "    &END SCF",
    "  &END DFT",
    "",
    "  &PRINT",
    "    &FORCES ON",
    "    &END FORCES",
    "  &END PRINT",
    "",
  ])


def _append_subsys_dft_block(lines, atoms, basis_set, functional, vectors, is_periodic_system):
  vector_a, vector_b, vector_c = vectors
  lines.append("  &SUBSYS")

  for element in sorted({atom[0] for atom in atoms}):
    if element not in VALENCE_ELECTRONS:
      raise ValueError(f"No valence electron mapping for element: {element}")

    lines.extend([
      f"    &KIND {element}",
      f"      ELEMENT   {element}",
      f"      BASIS_SET   {basis_set}",
      f"      POTENTIAL   GTH-{functional}-q{VALENCE_ELECTRONS[element]}",
      "    &END KIND",
    ])

  lines.extend([
    "    &CELL",
    "      PERIODIC XYZ" if is_periodic_system else "      PERIODIC NONE",
    f"      A     {format_vector(vector_a)}",
    f"      B     {format_vector(vector_b)}",
    f"      C     {format_vector(vector_c)}",
    "    &END CELL ",
    "    &COORD",
  ])
  lines.extend([format_atom_line(atom) for atom in atoms])
  lines.append("    &END COORD")

  if not is_periodic_system:
    lines.extend([
      "    &TOPOLOGY",
      "      &CENTER_COORDINATES",
      "      &END CENTER_COORDINATES",
      "    &END TOPOLOGY",
    ])

  lines.append("  &END SUBSYS")


def _append_mm_block(lines):
  lines.extend([
    "  &SUBSYS",
    "    &CELL",
    "    A     10.00000000    0.000000000    0.000000000",
    "    B     0.000000000    10.00000000    0.000000000",
    "    C     0.000000000    0.000000000    10.00000000",
    "    &END CELL ",
    "    &COORD",
    "$$coords:          S      x    y    z$$",
    "    &END COORD",
    "",
    "    &TOPOLOGY",
    "      CHARGE_BETA",
    "      CONNECTIVITY AMBER",
    "      CONN_FILE_NAME ! Add file name that contains connectivity data",
    "    &END TOPOLOGY ",
    "    &PRINT",
    "      &TOPOLOGY_INFO",
    "        AMBER_INFO",
    "      &END",
    "    &END ",
    "  &END SUBSYS",
    "  &MM",
    "    &FORCEFIELD",
    "      ! Add file name that contains force field parameters",
    "      PARMTYPE AMBER",
    "      &SPLINE",
    "        EMAX_SPLINE 10000",
    "      &END SPLINE",
    "    &END FORCEFIELD",
    "    &POISSON",
    "      &EWALD",
    "        EWALD_TYPE SPME",
    "        ALPHA .36",
    "        GMAX 128",
    "      &END EWALD",
    "    &END POISSON",
    "    &PRINT",
    "      &FF_INFO",
    "      $END",
    "      &FF_PARAMETER_FILE",
    "      &END",
    "    &END PRINT",
    "  &END MM",
  ])


def generateInputFile(cml, opts):
  title = opts["Title"]
  calculate = opts["Run Type"]
  method = opts["Method"]
  basis_set = opts["Basis Set"]
  functional = opts["Functional"]

  method_map = {
    "Electronic structure methods (DFT)": "QUICKSTEP",
    "Hybrid quantum classical (Not yet supported)": "QMMM",
    "Molecular Mechanics": "FIST",
  }
  try:
    cp2k_method = method_map[method]
  except KeyError as exc:
    raise ValueError(f"Invalid method type: {method}") from exc

  cell_vectors = generateCellVectors(str(cml))
  is_periodic_system = cell_vectors is not None
  vectors = cell_vectors if is_periodic_system else generateMolecularCell(str(cml))
  atoms = parse_atoms(str(cml), vectors=vectors)
  if is_periodic_system:
    atoms = deduplicate_periodic_atoms(atoms, vectors)

  lines = []
  _append_global(lines, title, calculate)

  lines.extend([
    "&FORCE_EVAL",
    f"  METHOD {cp2k_method}",
  ])

  if method == "Electronic structure methods (DFT)":
    _append_dft_block(lines, opts, functional, is_periodic_system)
    _append_subsys_dft_block(lines, atoms, basis_set, functional, vectors, is_periodic_system)
  elif method == "Molecular Mechanics":
    _append_mm_block(lines)

  lines.append("&END FORCE_EVAL")
  return "\n".join(lines)


def generateInput():
  stdin_str = sys.stdin.read()
  payload = json.loads(stdin_str)

  inp = generateInputFile(payload["cml"], payload["options"])
  base_name = payload["options"]["Filename Base"]

  files = [{"filename": f"{base_name}.inp", "contents": inp}]
  if DEBUG:
    files.append({"filename": "debug_info", "contents": stdin_str})

  return {
    "files": files,
    "mainFile": f"{base_name}.inp",
  }


if __name__ == "__main__":
  parser = argparse.ArgumentParser(f"Generate a {TARGET_NAME} input file.")
  parser.add_argument("--debug", action="store_true")
  parser.add_argument("--print-options", action="store_true")
  parser.add_argument("--generate-input", action="store_true")
  parser.add_argument("--display-name", action="store_true")
  parser.add_argument("--lang", nargs="?", default="en")
  args = vars(parser.parse_args())

  DEBUG = args["debug"]

  if args["display_name"]:
    print(TARGET_NAME)
  if args["print_options"]:
    print(json.dumps(getOptions()))
  elif args["generate_input"]:
    print(json.dumps(generateInput()))
