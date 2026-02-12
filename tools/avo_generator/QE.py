#!/usr/bin/env python3
"""
Quantum ESPRESSO input generator for Avogadro2.
Implements the core QE input-generation logic from qei.sh (cif2qe branch)
in a pure-Python plugin style.

High-level flow:
1. Parse CML cell + atoms (cartesian/fractional both supported).
2. Normalize atom list (dedupe for periodic systems, stable sorting).
3. Map UI options to QE blocks (&CONTROL/&SYSTEM/&ELECTRONS/...).
4. Emit optional extensions (DFT-D, Hubbard U, k-points).
"""

import argparse
import json
import math
import sys
import xml.etree.ElementTree as ET


TARGET_NAME = "QE"
DEBUG = False
FRACTIONAL_TOL = 1.0e-6
AUTO_K_RESOLUTION = 0.35

# Minimal valence electron reference used for nbnd estimation.
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

# Atomic masses used in ATOMIC_SPECIES output (fallback: 1.0 if unknown).
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

HUBBARD_U_EXPLICIT = {
  # Keep qei.sh presets for common 3d metals.
  "Ti": ("3d", 4.0),
  "Fe": ("3d", 5.0),
  "Cu": ("3d", 4.0),
  "Ni": ("3d", 6.0),
  "Co": ("3d", 5.0),
  "V": ("3d", 3.0),
  "Cr": ("3d", 3.0),
  "Mn": ("3d", 4.0),
  "Zn": ("3d", 8.0),
  "Sc": ("3d", 2.0),
}


def _build_hubbard_u_presets():
  """Build complete Hubbard U presets, keeping explicit values as priority."""
  # Start from explicit (highest-priority) element presets.
  presets = dict(HUBBARD_U_EXPLICIT)

  def add(elements, orbital, u_value):
    # Only fill elements not explicitly specified above.
    for element in elements:
      presets.setdefault(element, (orbital, u_value))

  # d-block metals
  add(["Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd"], "4d", 3.0)
  add(["Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"], "5d", 2.5)
  add(["Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"], "6d", 2.5)

  # f-block metals
  add(["La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"], "4f", 6.0)
  add(["Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"], "5f", 4.0)

  # s-block metals
  add(["Li", "Be"], "2s", 1.5)
  add(["Na", "Mg"], "3s", 1.5)
  add(["K", "Ca"], "4s", 1.5)
  add(["Rb", "Sr"], "5s", 1.5)
  add(["Cs", "Ba"], "6s", 1.5)
  add(["Fr", "Ra"], "7s", 1.5)

  # p-block metals
  add(["Al"], "3p", 2.0)
  add(["Ga"], "4p", 2.0)
  add(["In", "Sn"], "5p", 2.0)
  add(["Tl", "Pb", "Bi", "Po"], "6p", 2.0)
  add(["Nh", "Fl", "Mc", "Lv"], "7p", 2.0)

  return presets


HUBBARD_U_PRESETS = _build_hubbard_u_presets()


def getOptions():
  """Return Avogadro UI schema for this generator."""
  # Basic tab: core chemistry and run setup.
  job = {
    "tabName": "Basic",
    "Title": {"type": "string", "default": "", "hide": True, "toolTip": "Used as QE prefix"},
    "Filename Base": {"type": "string", "default": "qe"},
    "Calculation Type": {
      "type": "stringList",
      "default": 0,
      "values": ["scf", "relax", "vc-relax"],
    },
    "Charge": {"type": "integer", "default": 0, "minimum": -9, "maximum": 9},
    "Multiplicity": {
      "type": "integer",
      "default": 1,
      "minimum": 1,
      "maximum": 41,
      "toolTip": "Converted to tot_magnetization = multiplicity - 1",
    },
    "Theory": {
      "type": "string",
      "default": "pbesol",
      "toolTip": "Used for input_dft and pseudo filename; hsesol -> sla+pw+hse+psc",
    },
    "Basis": {
      "type": "stringList",
      "default": 2,
      "values": ["US", "PAW", "NC"],
      "toolTip": "Used for pseudo filename and default ecutrho multiplier",
    },
    "Ecutwfc (Ry)": {"type": "integer", "default": 50, "minimum": 1, "maximum": 9999},
    "DFT-D": {
      "type": "stringList",
      "default": 0,
      "values": ["None", "D3BJ", "D3ZERO", "D2"],
    },
  }

  # Advanced tab: optional expert features.
  advanced = {
    "tabName": "Advanced",
    "Use Hubbard U": {"type": "boolean", "default": False},
  }

  return {
    "userOptions": [job, advanced],
    "inputMoleculeFormat": "cml",
    "highlightStyles": _get_highlight_styles(),
  }


def _get_highlight_styles():
  """Define syntax highlight rules for generated QE input text."""
  rules = []

  # Numeric literals.
  rules.append({
    "patterns": [{"regexp": "\\b[+-]?[.0-9]+(?:[eEdD][+-]?[.0-9]+)?\\b"}],
    "format": {"preset": "literal"},
  })

  # QE section headers (&CONTROL, &SYSTEM, ...).
  rules.append({
    "patterns": [{"regexp": "&[A-Z_]+"}],
    "format": {"preset": "title"},
  })

  # Standalone uppercase block labels (e.g., K_POINTS, ATOMIC_SPECIES).
  rules.append({
    "patterns": [{"regexp": "^\\s*[A-Z_]{2,}\\b"}],
    "format": {"preset": "title"},
  })

  # Parameter keys on assignment lines.
  rules.append({
    "patterns": [{"regexp": "^\\s*[a-z_][a-z0-9_]*\\s*="}],
    "format": {"preset": "keyword"},
  })

  # Assigned values (right-hand side).
  rules.append({
    "patterns": [{"regexp": "^\\s*[a-z_][a-z0-9_]*\\s*=\\s*([^#\\n]+)"}],
    "format": {"preset": "property"},
  })

  # Trailing properties in uppercase data blocks.
  rules.append({
    "patterns": [{"regexp": "^\\s*[A-Z_]{2,}(\\s+[^#\\n]+)"}],
    "format": {"preset": "property"},
  })

  # Both '#' and '!' comments.
  rules.append({
    "patterns": [{"regexp": "#[^\\n]*"}, {"regexp": "![^\\n]*"}],
    "format": {"preset": "comment"},
  })

  return [{"style": "qe-default", "rules": rules}]


def parse_cml(cml):
  """Parse CML XML text into an ElementTree root element."""
  return ET.fromstring(cml)


def _cross(v1, v2):
  """3D vector cross product."""
  return [
    v1[1] * v2[2] - v1[2] * v2[1],
    v1[2] * v2[0] - v1[0] * v2[2],
    v1[0] * v2[1] - v1[1] * v2[0],
  ]


def _dot(v1, v2):
  """3D vector dot product."""
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]


def _vector_norm(v):
  """Euclidean norm of a 3D vector."""
  return math.sqrt(_dot(v, v))


def _cartesian_to_fractional(cart, vectors, tol=FRACTIONAL_TOL):
  """Convert Cartesian coordinates to fractional coordinates for a general cell."""
  # Build reciprocal-like basis via cross products.
  vector_a, vector_b, vector_c = vectors
  b_cross_c = _cross(vector_b, vector_c)
  c_cross_a = _cross(vector_c, vector_a)
  a_cross_b = _cross(vector_a, vector_b)
  # Cell volume is the triple product.
  volume = _dot(vector_a, b_cross_c)
  if abs(volume) < tol:
    raise ValueError("Invalid cell vectors: near-zero cell volume.")

  # Project Cartesian coordinate onto reciprocal basis.
  fx = _dot(cart, b_cross_c) / volume
  fy = _dot(cart, c_cross_a) / volume
  fz = _dot(cart, a_cross_b) / volume
  return [fx, fy, fz]


def _fractional_to_cartesian(frac, vectors):
  """Convert fractional coordinates to Cartesian coordinates."""
  vector_a, vector_b, vector_c = vectors
  fx, fy, fz = frac
  return [
    fx * vector_a[0] + fy * vector_b[0] + fz * vector_c[0],
    fx * vector_a[1] + fy * vector_b[1] + fz * vector_c[1],
    fx * vector_a[2] + fy * vector_b[2] + fz * vector_c[2],
  ]


def _wrap_fractional(value, tol=FRACTIONAL_TOL):
  """Wrap fractional value to [0,1), snapping boundary noise to 0."""
  wrapped = value - math.floor(value)
  if wrapped < tol or wrapped > 1.0 - tol:
    return 0.0
  return wrapped


def parse_cell(cml):
  """Parse crystal parameters from CML and return 3x3 Cartesian cell vectors."""
  root = parse_cml(cml)
  # Crystal block is optional in CML; missing means non-periodic molecule.
  crystal = root.find(".//{*}crystal")
  if crystal is None:
    return None

  # Collect a, b, c and alpha, beta, gamma.
  values = {"a": None, "b": None, "c": None, "alpha": None, "beta": None, "gamma": None}
  for scalar in crystal.findall("{*}scalar"):
    title = scalar.get("title")
    if title not in values or scalar.text is None:
      continue
    value = float(scalar.text)
    # QE math uses radians for angles.
    if title in ("alpha", "beta", "gamma"):
      value = math.radians(value)
    values[title] = value

  # Any missing scalar means we cannot build a valid periodic cell.
  if any(values[k] is None for k in values):
    return None

  # Convert lengths/angles to Cartesian lattice vectors.
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
  """Parse atoms from CML, preserving both cart and frac forms when possible."""
  root = parse_cml(cml)
  atom_array = root.find(".//{*}atomArray")
  atoms = []
  # Empty atom array yields empty atom list.
  if atom_array is None:
    return atoms

  for atom in atom_array.findall("{*}atom"):
    # Skip malformed atoms without element symbol.
    element = atom.get("elementType")
    if element is None:
      continue

    # Prefer Cartesian coordinates when available.
    x = atom.get("x3")
    y = atom.get("y3")
    z = atom.get("z3")
    if x is not None and y is not None and z is not None:
      cart = [float(x), float(y), float(z)]
      # Derive fractional coords if lattice vectors are known.
      frac = _cartesian_to_fractional(cart, vectors) if vectors is not None else None
      atoms.append({"element": element, "cart": cart, "frac": frac})
      continue

    # Fallback to fractional coordinates if Cartesian was absent.
    fx = atom.get("xFract")
    fy = atom.get("yFract")
    fz = atom.get("zFract")
    if fx is not None and fy is not None and fz is not None:
      frac = [float(fx), float(fy), float(fz)]
      # Derive Cartesian coords if lattice vectors are known.
      cart = _fractional_to_cartesian(frac, vectors) if vectors is not None else None
      atoms.append({"element": element, "cart": cart, "frac": frac})

  return atoms


def deduplicate_periodic_atoms(atoms, vectors, tol=FRACTIONAL_TOL):
  """Remove periodic duplicates by wrapped fractional coordinates + element key."""
  unique_atoms = []
  seen = set()
  for atom in atoms:
    # Ensure fractional coordinates exist for periodic keying.
    element = atom["element"]
    frac = atom["frac"]
    if frac is None:
      frac = _cartesian_to_fractional(atom["cart"], vectors, tol)

    # Wrap to unit cell and quantize by tolerance to avoid float jitter.
    wrapped = [_wrap_fractional(value, tol) for value in frac]
    key = (
      element,
      int(round(wrapped[0] / tol)),
      int(round(wrapped[1] / tol)),
      int(round(wrapped[2] / tol)),
    )
    # Keep first occurrence only.
    if key in seen:
      continue
    seen.add(key)
    # Store wrapped position in both fractional and Cartesian forms.
    cart = _fractional_to_cartesian(wrapped, vectors)
    unique_atoms.append({"element": element, "cart": cart, "frac": wrapped})

  return unique_atoms


def _k_grid_from_resolution(vectors, resolution):
  """Compute automatic Monkhorst-Pack grid from target k-resolution (1/Angstrom)."""
  # Use real-space cell lengths to estimate reciprocal sampling counts.
  lengths = [_vector_norm(v) for v in vectors]
  grid = []
  for length in lengths:
    # Degenerate axis fallback.
    if length <= 0:
      grid.append(1)
      continue
    # k ~= 2pi / (L * resolution), rounded up.
    k = int(math.ceil((2.0 * math.pi) / (length * resolution)))
    grid.append(max(1, k))
  return grid


def _estimate_nbnd(atoms, charge, factor, extra):
  """Estimate QE nbnd from valence electrons, charge and empirical margin."""
  # Start from total electrons corrected by charge sign convention.
  nelec = -float(charge)
  for atom in atoms:
    element = atom["element"]
    if element not in VALENCE_ELECTRONS:
      raise ValueError(f"No valence data for {element}")
    nelec += VALENCE_ELECTRONS[element]

  # Convert electrons to occupied bands and add safety margin.
  return int((nelec / 2.0) * float(factor) + int(extra))


def _calc_to_qe(calc_type):
  """Normalize UI calc type labels to QE calculation keywords."""
  key = str(calc_type).strip()
  # Keep backward compatibility with old capitalized labels.
  mapping = {
    "scf": "scf",
    "relax": "relax",
    "vc-relax": "vc-relax",
    "SCF": "scf",
    "Relax": "relax",
    "VC-Relax": "vc-relax",
  }
  try:
    return mapping[key]
  except KeyError as exc:
    raise ValueError(f"Invalid calculation type: {calc_type}") from exc


def _ecutrho_from_type(ecutwfc, pseudo_type):
  """Return recommended ecutrho multiplier based on pseudopotential family."""
  p = pseudo_type.upper()
  # US/PAW typically require higher augmentation cutoffs.
  if p == "US":
    return int(ecutwfc * 10)
  if p == "PAW":
    return int(ecutwfc * 6)
  return int(ecutwfc * 4)


def _pseudo_filename(element, functional, pseudo_type):
  """Build pseudo filename as '<Elem>-<functional>-<Basis>.UPF'."""
  func = (str(functional).strip() or "pbesol").lower()
  ptype = str(pseudo_type).strip().upper() or "NC"
  return f"{element}-{func}-{ptype}.UPF"


def _format_species_line(element, functional, pseudo_type):
  """Format one ATOMIC_SPECIES row with aligned element and mass columns."""
  mass = ATOMIC_MASSES.get(element, 1.0)
  return f"{element:<2s} {mass:9.5f}  {_pseudo_filename(element, functional, pseudo_type)}"


def _species_in_atom_order(atoms):
  """Return unique element list in first-appearance order from atom sequence."""
  species = []
  seen = set()
  for atom in atoms:
    element = atom["element"]
    # Preserve order of first appearance.
    if element in seen:
      continue
    seen.add(element)
    species.append(element)
  return species


def _hubbard_u_lines(species):
  """Generate HUBBARD block lines for species that have preset values."""
  # Emit only for species present in this structure.
  hubbard_species = [element for element in species if element in HUBBARD_U_PRESETS]
  if not hubbard_species:
    return []
  # Use ortho-atomic syntax as in qei.sh workflow.
  lines = ["HUBBARD {ortho-atomic}"]
  for element in hubbard_species:
    orbital, u_value = HUBBARD_U_PRESETS[element]
    lines.append(f"  U  {element}-{orbital}  {u_value:.1f}")
  return lines


def _opt(opts, key, default):
  """Safe option fetch with default fallback."""
  return opts.get(key, default)


def build_qe_input(cml, opts):
  """Build full Quantum ESPRESSO input text from CML + UI options."""
  # 1) Parse calculation mode and structural model.
  calculation = _calc_to_qe(_opt(opts, "Calculation Type", "scf"))
  cell = parse_cell(cml)
  is_periodic = cell is not None

  # For molecules without crystal info, use a large cubic box.
  if not is_periodic:
    cell = ([20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 20.0])

  # Parse atoms and remove periodic duplicates only for periodic systems.
  atoms = parse_atoms_with_fractional(cml, vectors=cell)
  if is_periodic:
    atoms = deduplicate_periodic_atoms(atoms, cell)

  # Sort atoms for stable output and derived sections.
  atoms = _sorted_atoms_by_element(atoms, cell)
  kgrid = _k_grid_from_resolution(cell, AUTO_K_RESOLUTION)
  shift = [0, 0, 0]

  if not atoms:
    raise ValueError("No atoms parsed from CML (x3/y3/z3 or xFract/yFract/zFract).")

  # Pull composition metrics used in SYSTEM block.
  species = _species_in_atom_order(atoms)
  nat = len(atoms)
  ntyp = len(species)

  # Basic cutoff setup from basis family and selected ecutwfc.
  ecutwfc = int(_opt(opts, "Ecutwfc (Ry)", 50))
  ecutrho = _ecutrho_from_type(ecutwfc, _opt(opts, "Basis", "NC"))

  # CONTROL block.
  lines = [
    "&CONTROL",
    f"  calculation = '{calculation}'",
    f"  prefix = '{(str(_opt(opts, 'Filename Base', 'qe')).strip() or 'qe')}'",
  ]

  # For SCF runs, keep verbose output by default.
  if calculation == "scf":
    lines.append("  verbosity = 'high'")

  # Fixed runtime paths and restart hint (kept as comment).
  lines.extend([
    "  outdir = './tmp'",
    "  !pseudo_dir = './pseudo'",
    "  !restart_mode = 'restart'",
  ])

  # Relax-specific controls.
  if calculation in ("relax", "vc-relax"):
    lines.extend([
      "  nstep = 200",
      "  forc_conv_thr = 1.0d-3",
      "  etot_conv_thr = 1.0d-4",
    ])

  # SYSTEM block header with cell/atom counts.
  lines.extend([
    "/",
    "&SYSTEM",
    "  ibrav = 0",
    f"  A = {cell[0][0]:.5f}",
    f"  nat = {nat}",
    f"  ntyp = {ntyp}",
  ])

  if calculation not in ("relax", "vc-relax"):
    # Keep nbnd auto only for SCF-like runs.
    nbnd = _estimate_nbnd(
      atoms,
      int(_opt(opts, "Charge", 0)),
      1.25,
      4,
    )
    lines.append(f"  nbnd = {nbnd}")

  # Core kinetic cutoffs.
  lines.extend([
    f"  ecutwfc = {ecutwfc}",
    f"  ecutrho = {ecutrho}",
  ])

  # Optional net charge.
  charge = int(_opt(opts, "Charge", 0))
  if charge != 0:
    lines.append(f"  tot_charge = {charge}")

  # Spin-polarized setup when multiplicity implies non-zero 2S.
  multiplicity = max(1, int(_opt(opts, "Multiplicity", 1)))
  # QE expects 2S (= multiplicity - 1) as tot_magnetization.
  total_magnetization = multiplicity - 1
  if total_magnetization != 0:
    lines.extend([
      "  nspin = 2",
      f"  tot_magnetization = {total_magnetization}",
      "  !starting_magnetization(1) = 0.8",
      "  occupations = 'smearing'",
      "  smearing = 'mv'",
      "  degauss = 0.001",
    ])

  # XC functional mapping.
  functional = (str(_opt(opts, "Theory", "pbesol")).strip() or "pbesol").lower()
  # pbesol is treated as implicit default (omit input_dft).
  if functional != "pbesol":
    if functional == "hsesol":
      lines.append("  input_dft = 'sla+pw+hse+psc'")
    else:
      lines.append(f"  input_dft = '{functional}'")

  # DFT-D dispersion mapping.
  dftd = str(_opt(opts, "DFT-D", "None")).strip().upper()
  if dftd == "D3BJ":
    lines.extend([
      "  vdw_corr = 'dft-d3'",
      "  dftd3_version = 4",
    ])
  elif dftd == "D3ZERO":
    lines.extend([
      "  vdw_corr = 'dft-d3'",
      "  dftd3_version = 3",
    ])
  elif dftd == "D2":
    lines.append("  vdw_corr = 'dft-d'")

  # ELECTRONS block (fixed defaults).
  lines.extend([
    "/",
    "&ELECTRONS",
    "  mixing_beta = 0.7",
    "  conv_thr = 1.0d-6",
    "  electron_maxstep = 200",
    "  !diagonalization = 'cg'",
    "  !startingwfc = 'random'",
  ])

  # Relax runs can proceed even if intermediate SCF is not fully converged.
  if calculation in ("relax", "vc-relax"):
    lines.append("  scf_must_converge = .false.")

  lines.append("/")

  # IONS block for geometry optimization.
  if calculation in ("relax", "vc-relax"):
    lines.extend([
      "&IONS",
      "  ion_dynamics = 'bfgs'",
      "/",
    ])

  # CELL block only for variable-cell relaxation.
  if calculation == "vc-relax":
    lines.extend([
      "&CELL",
      "  cell_dynamics = 'bfgs'",
      "  press_conv_thr = 0.5d0",
      "  cell_dofree = 'all'",
      "/",
    ])

  # Cell vectors in alat units.
  lines.append("CELL_PARAMETERS {alat}")
  alat = cell[0][0]
  for vector in cell:
    lines.append(
      "  {:.15f} {:.15f} {:.15f}".format(vector[0] / alat, vector[1] / alat, vector[2] / alat)
    )

  # Species table with generated pseudo filenames.
  lines.append("ATOMIC_SPECIES")
  pseudo_type = str(_opt(opts, "Basis", "NC")).strip().upper() or "NC"
  for element in species:
    lines.append("  " + _format_species_line(element, functional, pseudo_type))

  # Atomic positions in crystal coordinates.
  lines.append("ATOMIC_POSITIONS {crystal}")
  for atom in atoms:
    frac = atom["frac"]
    if frac is None:
      frac = _cartesian_to_fractional(atom["cart"], cell)
    lines.append(
      "  {:<2s} {: .15f} {: .15f} {: .15f}".format(atom["element"], frac[0], frac[1], frac[2])
    )

  if kgrid == [1, 1, 1] and shift == [0, 0, 0]:
    # Keep gamma-only shortcut for the minimal automatic mesh.
    lines.append("K_POINTS gamma")
  else:
    lines.append("K_POINTS automatic")
    lines.append(f"  {kgrid[0]} {kgrid[1]} {kgrid[2]}  {shift[0]} {shift[1]} {shift[2]}")

  if _opt(opts, "Use Hubbard U", False):
    # Append HUBBARD block after k-point section when enabled.
    lines.append("")
    lines.extend(_hubbard_u_lines(species))

  return "\n".join(lines) + "\n"


# Metal list used by sorting (metals first) and U-preset coverage checks.
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
  """Sort atoms deterministically: metals first, then element and position."""
  # Include fractional coordinates in key for stable intra-element ordering.
  def atom_key(atom):
    frac = atom["frac"]
    if frac is None:
      frac = _cartesian_to_fractional(atom["cart"], vectors)
    element = atom["element"]
    metal_group = 0 if element in _METALS else 1
    return (metal_group, element, frac[0], frac[1], frac[2])

  return sorted(atoms, key=atom_key)


def generateInput():
  """Avogadro plugin entry: consume JSON payload and return generated files."""
  # Plugin payload contract: {"cml": ..., "options": ...}
  payload = json.loads(sys.stdin.read())
  cml = payload["cml"]
  opts = payload["options"]

  # Default output naming convention: <base>.<calc>.in
  text = build_qe_input(cml, opts)
  base = str(opts["Filename Base"]).strip() or "qe"
  calculation = _calc_to_qe(_opt(opts, "Calculation Type", "scf"))
  filename = f"{base}.{calculation}.in"

  # Return file list + main file marker as required by Avogadro.
  files = [{"filename": filename, "contents": text, "highlightStyles": ["qe-default"]}]
  if DEBUG:
    files.append({"filename": "debug_info", "contents": json.dumps(payload)})

  return {
    "files": files,
    "mainFile": filename,
  }


if __name__ == "__main__":
  """CLI shim used by Avogadro plugin protocol and local debugging."""
  # Expose plugin-compatible command switches.
  parser = argparse.ArgumentParser(f"Generate a {TARGET_NAME} input file.")
  parser.add_argument("--debug", action="store_true")
  parser.add_argument("--print-options", action="store_true")
  parser.add_argument("--generate-input", action="store_true")
  parser.add_argument("--display-name", action="store_true")
  parser.add_argument("--lang", nargs="?", default="en")
  args = vars(parser.parse_args())

  DEBUG = args["debug"]

  # Entrypoint dispatch by selected command flag.
  if args["display_name"]:
    print(TARGET_NAME)
  if args["print_options"]:
    print(json.dumps(getOptions()))
  elif args["generate_input"]:
    print(json.dumps(generateInput()))
