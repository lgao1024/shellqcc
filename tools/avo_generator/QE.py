#!/usr/bin/env python3
"""
Quantum ESPRESSO input generator for Avogadro2.
Implements the core QE input-generation logic from qei.sh (cif2qe branch)
in a pure-Python plugin style.
"""

import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET


TARGET_NAME = "QuantumESPRESSO"
DEBUG = False
FRACTIONAL_TOL = 1.0e-6

VALENCE_ELECTRONS = {
  "H": 1, "He": 2,
  "Li": 3, "Be": 4, "B": 3, "C": 4, "N": 5, "O": 6, "F": 7, "Ne": 8,
  "Na": 9, "Mg": 10, "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Ar": 8,
  "K": 9, "Ca": 10,
  "Sc": 11, "Ti": 12, "V": 13, "Cr": 14, "Mn": 15, "Fe": 16, "Co": 17,
  "Ni": 18, "Cu": 11, "Zn": 12,
  "Ga": 13, "Ge": 4, "As": 5, "Se": 6, "Br": 7, "Kr": 8,
  "Y": 11, "Zr": 12, "Nb": 13, "Mo": 14,
  "Ru": 16, "Rh": 17, "Pd": 18, "Ag": 11, "Cd": 12,
  "In": 13, "Sn": 4, "Sb": 5, "Te": 6, "I": 7,
  "Ba": 10, "W": 14, "Au": 11, "Hg": 12, "Bi": 15,
}

ATOMIC_MASSES = {
  "H": 1.00794, "He": 4.00260,
  "Li": 6.941, "Be": 9.01218, "B": 10.811, "C": 12.0107, "N": 14.0067,
  "O": 15.999, "F": 18.9984, "Ne": 20.1797,
  "Na": 22.9898, "Mg": 24.305, "Al": 26.9815, "Si": 28.0855, "P": 30.9738,
  "S": 32.065, "Cl": 35.453, "Ar": 39.948,
  "K": 39.0983, "Ca": 40.078,
  "Sc": 44.9559, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961, "Mn": 54.9380,
  "Fe": 55.845, "Co": 58.9332, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38,
  "Ga": 69.723, "Ge": 72.64, "As": 74.9216, "Se": 78.96, "Br": 79.904,
  "Kr": 83.798,
  "Y": 88.9058, "Zr": 91.224, "Nb": 92.9064, "Mo": 95.96,
  "Ru": 101.07, "Rh": 102.9055, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.411,
  "In": 114.818, "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.9045,
  "Ba": 137.327, "W": 183.84, "Au": 196.9666, "Hg": 200.59, "Bi": 208.9804,
}


def getOptions():
  job = {
    "tabName": "Job",
    "Title": {"type": "string", "default": "qe_job", "toolTip": "Used as QE prefix"},
    "Filename Base": {"type": "string", "default": "qe"},
    "Calculation Type": {
      "type": "stringList",
      "default": 0,
      "values": ["SCF", "Relax", "VC-Relax"],
    },
    "Charge": {"type": "integer", "default": 0, "minimum": -9, "maximum": 9},
    "Total Magnetization": {
      "type": "integer",
      "default": 0,
      "minimum": -20,
      "maximum": 20,
      "toolTip": "0 means non-spin-polarized by default",
    },
    "Input DFT": {
      "type": "string",
      "default": "PBEsol",
      "toolTip": "Leave PBEsol to omit input_dft; use hsesol for sla+pw+hse+psc",
    },
  }

  pseudo = {
    "tabName": "Pseudo",
    "Pseudo Directory": {"type": "string", "default": "./pseudo"},
    "Pseudo Suffix": {
      "type": "string",
      "default": ".UPF",
      "toolTip": "default pseudo filename pattern: <Elem><suffix>",
    },
    "Pseudo Type": {
      "type": "stringList",
      "default": 0,
      "values": ["US", "PAW", "NC"],
      "toolTip": "Used for default ecutrho multiplier",
    },
    "Ecutwfc (Ry)": {"type": "integer", "default": 50, "minimum": 1, "maximum": 9999},
    "Use DFT-D3": {"type": "boolean", "default": True},
    "DFT-D3 Version": {"type": "integer", "default": 4, "minimum": 2, "maximum": 5},
  }

  kpoints = {
    "tabName": "K-Points",
    "K-Point Mode": {
      "type": "stringList",
      "default": 0,
      "values": ["Auto by Resolution", "Manual Grid", "Gamma"],
    },
    "K-Resolution (1/A)": {"type": "number", "default": 0.35, "minimum": 0.01, "maximum": 10.0},
    "K-Point Grid": {"type": "string", "default": "4 4 4"},
    "K-Point Shift": {"type": "string", "default": "0 0 0"},
    "nqx1": {"type": "integer", "default": 1, "minimum": 1, "maximum": 48},
    "nqx2": {"type": "integer", "default": 1, "minimum": 1, "maximum": 48},
    "nqx3": {"type": "integer", "default": 2, "minimum": 1, "maximum": 48},
  }

  scf = {
    "tabName": "SCF",
    "Set nbnd": {"type": "boolean", "default": True},
    "nbnd Factor": {"type": "number", "default": 1.25, "minimum": 1.0, "maximum": 3.0},
    "nbnd Extra": {"type": "integer", "default": 4, "minimum": 0, "maximum": 200},
    "Mixing Beta": {"type": "number", "default": 0.7, "minimum": 0.0, "maximum": 1.0},
    "SCF Conv Thr": {"type": "string", "default": "1.0d-6"},
    "Electron Max Steps": {"type": "integer", "default": 200, "minimum": 1, "maximum": 5000},
  }

  relax = {
    "tabName": "Relax",
    "Relax Max Steps": {"type": "integer", "default": 200, "minimum": 1, "maximum": 5000},
    "forc_conv_thr": {"type": "string", "default": "1.0d-3"},
    "etot_conv_thr": {"type": "string", "default": "1.0d-4"},
  }

  runtime = {
    "tabName": "Runtime",
    "Outdir": {"type": "string", "default": "./tmp"},
  }

  return {
    "userOptions": [job, pseudo, kpoints, scf, relax, runtime],
    "inputMoleculeFormat": "cml",
    "highlightStyles": _get_highlight_styles(),
  }


def _get_highlight_styles():
  rules = []

  rules.append({
    "patterns": [{"regexp": "\\b[+-]?[.0-9]+(?:[eEdD][+-]?[.0-9]+)?\\b"}],
    "format": {"preset": "literal"},
  })

  rules.append({
    "patterns": [{"regexp": "&[A-Z_]+"}],
    "format": {"preset": "title"},
  })

  rules.append({
    "patterns": [{"regexp": "^\\s*[A-Z_]{2,}\\b"}],
    "format": {"preset": "title"},
  })

  rules.append({
    "patterns": [{"regexp": "^\\s*[a-z_][a-z0-9_]*\\s*="}],
    "format": {"preset": "keyword"},
  })

  rules.append({
    "patterns": [{"regexp": "^\\s*[a-z_][a-z0-9_]*\\s*=\\s*([^#\\n]+)"}],
    "format": {"preset": "property"},
  })

  rules.append({
    "patterns": [{"regexp": "^\\s*[A-Z_]{2,}(\\s+[^#\\n]+)"}],
    "format": {"preset": "property"},
  })

  rules.append({
    "patterns": [{"regexp": "#[^\\n]*"}, {"regexp": "![^\\n]*"}],
    "format": {"preset": "comment"},
  })

  return [{"style": "qe-default", "rules": rules}]


def parse_cml(cml):
  return ET.fromstring(cml)


def _cross(v1, v2):
  return [
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0],
  ]


def _dot(v1, v2):
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


def _vector_norm(v):
  return math.sqrt(_dot(v, v))


def _cartesian_to_fractional(cart, vectors, tol=FRACTIONAL_TOL):
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


def _wrap_fractional(value, tol=FRACTIONAL_TOL):
  wrapped = value - math.floor(value)
  if wrapped < tol or wrapped > 1.0 - tol:
    return 0.0
  return wrapped


def parse_cell(cml):
  root = parse_cml(cml)
  crystal = root.find(".//{*}crystal")
  if crystal is None:
    return None

  values = {"a": None, "b": None, "c": None, "alpha": None, "beta": None, "gamma": None}
  for scalar in crystal.findall("{*}scalar"):
    title = scalar.get("title")
    if title not in values or scalar.text is None:
      continue
    value = float(scalar.text)
    if title in ("alpha", "beta", "gamma"):
      value = math.radians(value)
    values[title] = value

  if any(values[k] is None for k in values):
    return None

  a = values["a"]
  b = values["b"]
  c = values["c"]
  alpha = values["alpha"]
  beta = values["beta"]
  gamma = values["gamma"]

  v1 = [a, 0.0, 0.0]
  v2 = [b * math.cos(gamma), b * math.sin(gamma), 0.0]
  cx = c * math.cos(beta)
  cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
  cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
  v3 = [cx, cy, cz]
  return v1, v2, v3


def parse_atoms_with_fractional(cml, vectors=None):
  root = parse_cml(cml)
  atom_array = root.find(".//{*}atomArray")
  atoms = []
  if atom_array is None:
    return atoms

  for atom in atom_array.findall("{*}atom"):
    element = atom.get("elementType")
    if element is None:
      continue

    x = atom.get("x3")
    y = atom.get("y3")
    z = atom.get("z3")
    if x is not None and y is not None and z is not None:
      cart = [float(x), float(y), float(z)]
      frac = _cartesian_to_fractional(cart, vectors) if vectors is not None else None
      atoms.append({"element": element, "cart": cart, "frac": frac})
      continue

    fx = atom.get("xFract")
    fy = atom.get("yFract")
    fz = atom.get("zFract")
    if fx is not None and fy is not None and fz is not None:
      frac = [float(fx), float(fy), float(fz)]
      cart = _fractional_to_cartesian(frac, vectors) if vectors is not None else None
      atoms.append({"element": element, "cart": cart, "frac": frac})

  return atoms


def deduplicate_periodic_atoms(atoms, vectors, tol=FRACTIONAL_TOL):
  unique_atoms = []
  seen = set()
  for atom in atoms:
    element = atom["element"]
    frac = atom["frac"]
    if frac is None:
      frac = _cartesian_to_fractional(atom["cart"], vectors, tol)

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
    unique_atoms.append({"element": element, "cart": cart, "frac": wrapped})

  return unique_atoms


def _parse_triplet(text, fallback):
  try:
    values = [int(x) for x in str(text).split()]
  except ValueError:
    values = []
  if len(values) != 3:
    return fallback
  return values


def _k_grid_from_resolution(vectors, resolution):
  lengths = [_vector_norm(v) for v in vectors]
  grid = []
  for length in lengths:
    if length <= 0:
      grid.append(1)
      continue
    k = int(math.ceil((2.0 * math.pi) / (length * resolution)))
    grid.append(max(1, k))
  return grid


def _estimate_nbnd(atoms, charge, factor, extra):
  nelec = -float(charge)
  for atom in atoms:
    element = atom["element"]
    if element not in VALENCE_ELECTRONS:
      raise ValueError(f"No valence data for {element}")
    nelec += VALENCE_ELECTRONS[element]

  return int((nelec / 2.0) * float(factor) + int(extra))


def _calc_to_qe(calc_type):
  mapping = {
    "SCF": "scf",
    "Relax": "relax",
    "VC-Relax": "vc-relax",
  }
  try:
    return mapping[calc_type]
  except KeyError as exc:
    raise ValueError(f"Invalid calculation type: {calc_type}") from exc


def _ecutrho_from_type(ecutwfc, pseudo_type):
  p = pseudo_type.upper()
  if p == "US":
    return int(ecutwfc * 10)
  if p == "PAW":
    return int(ecutwfc * 6)
  return int(ecutwfc * 4)


def _pseudo_filename(element, suffix):
  return f"{element}{suffix}"


def _format_species_line(element, suffix):
  mass = ATOMIC_MASSES.get(element, 1.0)
  return f"{element:>3s} {mass:9.5f}  {_pseudo_filename(element, suffix)}"


def _opt(opts, key, default):
  return opts.get(key, default)


def build_qe_input(cml, opts):
  calculation = _calc_to_qe(_opt(opts, "Calculation Type", "SCF"))
  cell = parse_cell(cml)
  is_periodic = cell is not None

  if not is_periodic:
    cell = ([20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0])

  atoms = parse_atoms_with_fractional(cml, vectors=cell)
  if is_periodic:
    atoms = deduplicate_periodic_atoms(atoms, cell)

  atoms = _sorted_atoms_by_element(atoms, cell)

  if not atoms:
    raise ValueError("No atoms parsed from CML (x3/y3/z3 or xFract/yFract/zFract).")

  species = sorted({atom["element"] for atom in atoms})
  nat = len(atoms)
  ntyp = len(species)

  ecutwfc = int(_opt(opts, "Ecutwfc (Ry)", 50))
  ecutrho = _ecutrho_from_type(ecutwfc, _opt(opts, "Pseudo Type", "NC"))

  lines = [
    "&CONTROL",
    f"  calculation = '{calculation}'",
    f"  prefix = '{_opt(opts, 'Title', 'qe_job')}'",
  ]

  if calculation == "scf":
    lines.append("  verbosity = 'high'")

  lines.extend([
    f"  outdir = '{_opt(opts, 'Outdir', './tmp')}'",
    f"  pseudo_dir = '{_opt(opts, 'Pseudo Directory', './pseudo')}'",
    "  !restart_mode = 'restart'",
  ])

  if calculation in ("relax", "vc-relax"):
    lines.extend([
      f"  nstep = {_opt(opts, 'Relax Max Steps', 200)}",
      f"  forc_conv_thr = {_opt(opts, 'forc_conv_thr', '1.0d-3')}",
      f"  etot_conv_thr = {_opt(opts, 'etot_conv_thr', '1.0d-4')}",
    ])

  lines.extend([
    "/",
    "&SYSTEM",
    "  ibrav = 0",
    f"  A = {cell[0][0]:.5f}",
    f"  nat = {nat}",
    f"  ntyp = {ntyp}",
  ])

  if _opt(opts, "Set nbnd", True):
    nbnd = _estimate_nbnd(
      atoms,
      int(_opt(opts, "Charge", 0)),
      float(_opt(opts, "nbnd Factor", 1.25)),
      int(_opt(opts, "nbnd Extra", 4)),
    )
    lines.append(f"  nbnd = {nbnd}")

  lines.extend([
    f"  ecutwfc = {ecutwfc}",
    f"  ecutrho = {ecutrho}",
  ])

  charge = int(_opt(opts, "Charge", 0))
  if charge != 0:
    lines.append(f"  tot_charge = {charge}")

  total_magnetization = int(_opt(opts, "Total Magnetization", 0))
  if total_magnetization != 0:
    lines.extend([
      "  nspin = 2",
      f"  tot_magnetization = {total_magnetization}",
      "  !starting_magnetization(1) = 0.8",
      "  occupations = 'smearing'",
      "  smearing = 'mv'",
      "  degauss = 0.001",
    ])

  input_dft = str(_opt(opts, "Input DFT", "PBEsol")).strip()
  if input_dft and input_dft.lower() != "pbesol":
    if input_dft.lower() == "hsesol":
      lines.append("  input_dft = 'sla+pw+hse+psc'")
    else:
      lines.append(f"  input_dft = '{input_dft}'")

  if _opt(opts, "Use DFT-D3", True):
    lines.extend([
      "  vdw_corr = 'dft-d3'",
      f"  dftd3_version = {_opt(opts, 'DFT-D3 Version', 4)}",
    ])

  lines.extend([
    f"  nqx1 = {_opt(opts, 'nqx1', 1)}",
    f"  nqx2 = {_opt(opts, 'nqx2', 1)}",
    f"  nqx3 = {_opt(opts, 'nqx3', 2)}",
    "/",
    "&ELECTRONS",
    f"  mixing_beta = {_opt(opts, 'Mixing Beta', 0.7)}",
    f"  conv_thr = {_opt(opts, 'SCF Conv Thr', '1.0d-6')}",
    f"  electron_maxstep = {_opt(opts, 'Electron Max Steps', 200)}",
    "  !diagonalization = 'cg'",
    "  !startingwfc = 'random'",
  ])

  if calculation in ("relax", "vc-relax"):
    lines.append("  scf_must_converge = .false.")

  lines.append("/")

  if calculation in ("relax", "vc-relax"):
    lines.extend([
      "&IONS",
      "  ion_dynamics = 'bfgs'",
      "/",
    ])

  if calculation == "vc-relax":
    lines.extend([
      "&CELL",
      "  cell_dynamics = 'bfgs'",
      "  press_conv_thr = 0.5d0",
      "  cell_dofree = 'all'",
      "/",
    ])

  lines.append("CELL_PARAMETERS {alat}")
  alat = cell[0][0]
  for vector in cell:
    lines.append(
      "  {:.15f} {:.15f} {:.15f}".format(vector[0] / alat, vector[1] / alat, vector[2] / alat)
    )

  lines.append("ATOMIC_SPECIES")
  suffix = str(_opt(opts, "Pseudo Suffix", ".UPF"))
  for element in species:
    lines.append("  " + _format_species_line(element, suffix))

  lines.append("ATOMIC_POSITIONS {crystal}")
  for atom in atoms:
    frac = atom["frac"]
    if frac is None:
      frac = _cartesian_to_fractional(atom["cart"], cell)
    lines.append(
      "{:<2s} {: .15f} {: .15f} {: .15f}".format(atom["element"], frac[0], frac[1], frac[2])
    )

  k_mode = _opt(opts, "K-Point Mode", "Auto by Resolution")
  if k_mode == "Gamma":
    lines.append("K_POINTS gamma")
  else:
    lines.append("K_POINTS automatic")
    if k_mode == "Auto by Resolution":
      resolution = float(_opt(opts, "K-Resolution (1/A)", 0.35))
      kgrid = _k_grid_from_resolution(cell, resolution)
    else:
      kgrid = _parse_triplet(_opt(opts, "K-Point Grid", "4 4 4"), [4, 4, 4])
    shift = _parse_triplet(_opt(opts, "K-Point Shift", "0 0 0"), [0, 0, 0])
    lines.append(f"{kgrid[0]} {kgrid[1]} {kgrid[2]}  {shift[0]} {shift[1]} {shift[2]}")

  return "\n".join(lines) + "\n"


_METALS = {
  "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
  "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
  "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
  "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
  "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
  "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
  "Rg", "Cn", "Nh", "Fl", "Mc", "Lv",
}


def _sorted_atoms_by_element(atoms, vectors):
  def atom_key(atom):
    frac = atom["frac"]
    if frac is None:
      frac = _cartesian_to_fractional(atom["cart"], vectors)
    element = atom["element"]
    metal_group = 0 if element in _METALS else 1
    return (metal_group, element, frac[0], frac[1], frac[2])

  return sorted(atoms, key=atom_key)


def generateInput():
  payload = json.loads(sys.stdin.read())
  cml = payload["cml"]
  opts = payload["options"]

  text = build_qe_input(cml, opts)
  base = opts["Filename Base"]

  files = [{"filename": f"{base}.in", "contents": text, "highlightStyles": ["qe-default"]}]
  if DEBUG:
    files.append({"filename": "debug_info", "contents": json.dumps(payload)})

  return {
    "files": files,
    "mainFile": f"{base}.in",
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
