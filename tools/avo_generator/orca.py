"""
/*****************************************************************************
  This source file is part of the Avogadro project. Version 260212

  This source code is released under the New BSD License, (the "License").
******************************************************************************/
"""

import argparse
import json
import sys


TARGET_NAME = "ORCA"
DEBUG = False


def getOptions():
  """Return Avogadro option schema for ORCA input generation."""
  # Basic tab: core electronic structure setup.
  user_options = {
    "tabName": "Basic",
    "Title": {
      "type": "string",
      "default": "",
      "toolTip": "Title of the input file",
      "hide": False,
    },
    "Processor Cores": {
      "type": "integer",
      "default": 24,
      "minimum": 1,
    },
    "Memory": {
      "type": "integer",
      "default": 240,
      "minimum": 1,
      "maximum": 999,
    },
    "Calculation Type": {
      "type": "stringList",
      "default": 1,
      "values": [
        "Single Point",
        "Geometry Optimization",
        "Frequencies",
        "Transition State",
        "Transition State (NEB method)",
        "Intrinsic Reaction Coordinate",
        "Dynamics",
      ],
      "toolTip": "Type of calculation to perform",
    },
    "Theory": {
      "type": "stringList",
      "default": 1,
      "values": [
        "XTB",
        "r2SCAN-3C",
        "PBEh-3C",
        "PBE",
        "revPBE",
        "PBE0",
        "TPSS",
        "TPSSh",
        "TPSS0",
        "M06L",
        "M062X",
        "BLYP",
        "B3LYP",
        "B3LYP/G",
        "CAM-B3LYP",
        "B97-3C",
        "wB97X-D4",
        "wB97M-V",
        "wB97X-2",
        "PWPB95",
        "HF",
        "MP2",
        "CCSD",
        "CCSD(T)",
      ],
    },
    "RI Approximation": {
      "type": "stringList",
      "default": 0,
      "values": ["None", "NORI", "RIJK", "RIJONX", "RIJCOSX"],
    },
    "Dispersion Correction": {
      "type": "stringList",
      "default": 0,
      "values": ["None", "D3ZERO", "D3BJ", "D4"],
      "toolTip": "Any added dispersion corrections",
      "hide": False,
    },
    "Basis": {
      "type": "stringList",
      "default": 2,
      "values": [
        "def2-SVP",
        "def2-TZVP",
        "def2-TZVPP",
        "def2-TZVPPD",
        "def2-QZVPP",
        "def2-QZVPPD",
        "ma-def2-SVP",
        "ma-def2-TZVP",
        "ma-def2-QZVPP",
        "cc-pVDZ",
        "cc-pVTZ",
        "cc-pVQZ",
        "aug-cc-pVDZ",
        "aug-cc-pVTZ",
        "aug-cc-pVQZ",
        "6-31G(d)",
        "6-311G(d,p)",
      ],
      "toolTip": "Gaussian basis set",
    },
    "Solvent": {
      "type": "stringList",
      "default": 0,
      "values": [
        "None (gas)",
        "-",
        "Water",
        "Acetonitrile",
        "Acetone",
        "Ethanol",
        "Methanol",
        "CCl4",
        "CH2Cl2",
        "Chloroform",
        "DMSO",
        "DMF",
        "Hexane",
        "Toluene",
        "Pyridine",
        "THF",
        "Toluene",
      ],
      "toolTip": "Solvent",
    },
    "Solvation Model": {
      "type": "stringList",
      "default": 0,
      "values": ["CPCM", "SMD"],
      "toolTip": "Solvation Method",
    },
    "Filename Base": {
      "type": "string",
      "default": "orca",
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
    "Keywords": {
      "type": "string",
      "default": "SlowConv DefGrid3 TightSCF",
    },
    "Excitation States": {
      "type": "boolean",
      "default": False,
    },
  }

  # Dynamics tab: AIMD-specific controls.
  aimd_options = {
    "tabName": "Dynamics",
    "AIMD TimeStep": {
      "type": "string",
      "default": "0.5_fs",
    },
    "AIMD Initvel": {
      "type": "integer",
      "default": 300,
      "maximum": 9999,
      "minimum": 0,
    },
    "AIMD Thermostat Type": {
      "type": "stringList",
      "default": 1,
      "values": ["Berendsen", "CSVR", "NHC", "None"],
    },
    "AIMD Thermostat Temp": {
      "type": "integer",
      "default": 300,
      "maximum": 1000,
      "minimum": 0,
    },
    "AIMD Thermostat Time": {
      "type": "string",
      "default": "30_fs",
    },
    "AIMD Restart": {
      "type": "stringList",
      "default": 0,
      "values": ["Not Restart", "Restart", "Restart IfExists"],
    },
    "AIMD Output Trajectory": {
      "type": "stringList",
      "default": 0,
      "values": ["Position", "Position, Velocity and Force", "All", "None"],
    },
    "AIMD RunTime": {
      "type": "integer",
      "default": 2000,
      "maximum": 999999999,
      "minimum": 0,
    },
  }

  return {
    "userOptions": [user_options, aimd_options],
    "inputMoleculeFormat": "cjson",
    "highlightStyles": _get_highlight_styles(),
  }


def _resolve_calc_mode(calculation_type):
  """Map UI calculation label to ORCA route keyword and geometry-control flag."""
  mapping = {
    "Single Point": ("SP", False),
    "Geometry Optimization": ("Opt", True),
    "Frequencies": ("Opt Freq", True),
    "Dynamics": ("MD", True),
    "Transition State": ("OptTS", True),
    "Transition State (NEB method)": ("NEB-TS Freq", True),
    "Intrinsic Reaction Coordinate": ("IRC", True),
  }
  try:
    return mapping[calculation_type]
  except KeyError as exc:
    raise ValueError(f"Unhandled calculation type: {calculation_type}") from exc


def _resolve_solvation(solvent, solvation_model):
  """Resolve solvent/solvation model settings into ORCA route fragments."""
  # Explicit gas-phase selection disables solvent keyword.
  if "None" in solvent:
    return ""
  if solvation_model == "CPCM":
    return f"CPCM({solvent})"
  if solvation_model == "SMD":
    return "CPCM"
  return ""


def _resolve_basis_with_aux(ri, basis):
  """Resolve RI keyword and required auxiliary basis for the selected main basis."""
  # RIJ-like approximations use Coulomb-fitting auxiliary basis sets.
  rij_basis = {
    "6-31G(d)": "AutoAux",
    "6-311G(d,p)": "AutoAux",
    "cc-pVDZ": "Def2/J",
    "cc-pVTZ": "Def2/J",
    "cc-pVQZ": "Def2/J",
    "aug-cc-pVDZ": "AutoAux",
    "aug-cc-pVTZ": "AutoAux",
    "aug-cc-pVQZ": "AutoAux",
    "def2-SVP": "Def2/J",
    "def2-TZVP": "Def2/J",
    "def2-TZVPP": "Def2/J",
    "def2-TZVPPD": "AutoAux",
    "def2-QZVPP": "Def2/J",
    "def2-QZVPPD": "AutoAux",
    "ma-def2-SVP": "AutoAux",
    "ma-def2-TZVP": "AutoAux",
    "ma-def2-QZVPP": "AutoAux",
  }

  # RIJK-like approximations require JK auxiliary basis sets.
  rijk_basis = {
    "6-31G(d)": "AutoAux",
    "6-311G(d,p)": "AutoAux",
    "cc-pVDZ": "cc-pVDZ/JK",
    "cc-pVTZ": "cc-pVTZ/JK",
    "cc-pVQZ": "cc-pVQZ/JK",
    "aug-cc-pVDZ": "aug-cc-pVDZ/JK",
    "aug-cc-pVTZ": "aug-cc-pVTZ/JK",
    "aug-cc-pVQZ": "aug-cc-pVQZ/JK",
    "def2-SVP": "Def2/JK",
    "def2-TZVP": "Def2/JK",
    "def2-TZVPP": "Def2/JK",
    "def2-TZVPPD": "aug-cc-pVTZ/JK",
    "def2-QZVPP": "Def2/JK",
    "def2-QZVPPD": "aug-cc-pVQZ/JK",
    "ma-def2-SVP": "aug-cc-pVDZ/JK",
    "ma-def2-TZVP": "aug-cc-pVTZ/JK",
    "ma-def2-QZVPP": "aug-cc-pVQZ/JK",
  }

  # No RI approximation: keep original basis only.
  if ri == "None":
    return "", basis
  # Explicit NORI still adds a route token but no aux basis suffix.
  if ri == "NORI":
    return " NORI", basis

  # Choose auxiliary basis family by RI method class.
  if ri in ("RIJONX", "RIJCOSX"):
    aux_basis = rij_basis[basis]
  else:
    aux_basis = rijk_basis[basis]

  return f" {ri}", f"{basis} {aux_basis}"


def _build_orca_code(opts, calc_str):
  """Build the ORCA route line from theory/basis/RI/dispersion/solvation options."""
  theory = opts["Theory"]
  basis = opts["Basis"]
  addopts = opts["Keywords"]
  dispersion = opts["Dispersion Correction"]
  ri = opts["RI Approximation"]

  # Resolve RI and possibly append auxiliary basis info.
  ri_str, basis = _resolve_basis_with_aux(ri, basis)
  dispersion_str = "" if dispersion == "None" else f" {dispersion}"

  # Composite methods already bundle basis/dispersion choices.
  if "-3c" in theory or "-3C" in theory or "XTB" in theory:
    code = f"{calc_str} {theory}{ri_str}"
  else:
    code = f"{calc_str} {theory}{dispersion_str}{ri_str} {basis}"

  # Solvation keyword is appended when a solvent is selected.
  solvation = _resolve_solvation(opts["Solvent"], opts["Solvation Model"])
  if solvation:
    code = f"{code} {solvation}"
  # Keep user-entered free-form options at the end.
  if addopts:
    code = f"{code} {addopts}"

  return code


def _append_md_block(lines, opts):
  """Append ORCA %md block from AIMD options."""
  lines.extend([
    "%md",
    f"   Timestep {opts['AIMD TimeStep']}",
  ])

  # Restart mode changes how initial velocity is handled.
  restart_mode = str(opts["AIMD Restart"])
  init_vel = opts["AIMD Initvel"]
  if restart_mode != "Not Restart":
    lines.append(f"   Initvel {init_vel}_K No_Overwrite")
  else:
    lines.append(f"   Initvel {init_vel}_K")

  # Thermostat control line.
  lines.append(
    "   Thermostat "
    + f"{opts['AIMD Thermostat Type']} "
    + f"{opts['AIMD Thermostat Temp']}_K Timecon {opts['AIMD Thermostat Time']}"
  )

  # Optional trajectory dumps.
  traj = str(opts["AIMD Output Trajectory"])
  if traj != "None":
    lines.append('   Dump Position Stride 1 Filename "position.xyz"')
    if traj != "Position":
      lines.append('   Dump Force Stride 1 Filename "force.xyz"')
      lines.append('   Dump Velocity Stride 1 Filename "velocity.xyz"')
      if traj != "Position, Velocity and Force":
        lines.append('   Dump GBW Stride 10 Filename "wfc"')
        lines.append('   Dump EnGrad Stride 10 Filename "EnGrad"')

  # Restart directives are emitted at the end of the block.
  if restart_mode != "Not Restart":
    lines.append(f"   {restart_mode}")

  lines.extend([
    f"   Run {opts['AIMD RunTime']} CenterCOM",
    "end",
    "",
  ])


def _append_geom_constraints(lines, cjson):
  """Append %geom constraint entries from cjson constraints/frozen atom data."""
  has_constraints = "constraints" in cjson
  has_frozen = "atoms" in cjson and "frozen" in cjson["atoms"]
  if not has_constraints and not has_frozen:
    return

  lines.append("   Constraints")

  # Distance/angle/dihedral constraints.
  if has_constraints:
    for constraint in cjson["constraints"]:
      if len(constraint) == 3:
        value, atom1, atom2 = constraint
        lines.append(f"      {{ B {atom1} {atom2} {value:.6f} C }}")
      elif len(constraint) == 4:
        value, atom1, atom2, atom3 = constraint
        lines.append(f"      {{ A {atom1} {atom2} {atom3} {value:.6f} C }}")
      elif len(constraint) == 5:
        value, atom1, atom2, atom3, atom4 = constraint
        lines.append(f"      {{ D {atom1} {atom2} {atom3} {atom4} {value:.6f} C }}")

  # Frozen atom modes:
  # - length == atom_count: freeze whole atoms
  # - length == 3*atom_count: freeze x/y/z selectively
  if has_frozen:
    frozen = cjson["atoms"]["frozen"]
    atom_count = len(cjson["atoms"]["elements"]["number"])

    if len(frozen) == atom_count:
      for i, value in enumerate(frozen):
        if value == 1:
          lines.append(f"      {{ C {i} C }}")
    elif len(frozen) == 3 * atom_count:
      for i in range(0, len(frozen), 3):
        if frozen[i] == 0:
          lines.append(f"      {{ X {i} C }}")
        if frozen[i + 1] == 0:
          lines.append(f"      {{ Y {i} C }}")
        if frozen[i + 2] == 0:
          lines.append(f"      {{ Z {i} C }}")

  lines.append("   end")


def _append_geom_block(lines, calculate, cjson):
  """Append ORCA %geom control block including optional constraints."""
  lines.append("%geom")
  _append_geom_constraints(lines, cjson)
  # Transition-state optimization benefits from periodic Hessian refresh.
  if calculate == "Transition State":
    lines.append("   ReCalc_Hess 100")
  lines.extend([
    "   MaxIter 1024",
    "end",
    "",
  ])


def generateInputFile(opts, cjson):
  """Generate ORCA input text from options + cjson molecular data."""
  title = opts["Title"]
  calculate = opts["Calculation Type"]
  charge = opts["Charge"]
  multiplicity = opts["Multiplicity"]
  # Convert total memory (GB) into ORCA's per-core MB setting.
  n_cores = int(opts["Processor Cores"])
  memory = int((opts["Memory"] * 1024) / n_cores)

  # Resolve route keyword and whether geometry-control block is needed.
  calc_str, geomcontrol = _resolve_calc_mode(calculate)
  code = _build_orca_code(opts, calc_str)

  # Header and core runtime controls.
  lines = [
    f"# {title}",
    f"%maxcore {memory}",
  ]

  # Enable parallel block only when more than one core is requested.
  if n_cores > 1:
    lines.append(f"%pal nprocs {n_cores} end")

  # Geometry jobs use MiniPrint, SP jobs keep explicit orbital/basis prints.
  if geomcontrol:
    lines.append(f"! {code} MiniPrint")
  else:
    lines.extend([
      f"! {code} NoAutoStart PrintMOs PrintBasis",
      "# ! MORead",
      f"# %Moinp \"{opts['Filename Base']}_restart.gbw\"",
    ])

  lines.append("")

  # SMD requires an explicit %cpcm block with smd=true.
  if "None" not in opts["Solvent"] and opts["Solvation Model"] == "SMD":
    lines.extend([
      "%cpcm",
      "   smd true",
      f"   SMDSolvent \"{opts['Solvent']}\"",
      "end",
      "",
    ])

  # Optional TDDFT excited-state setup.
  if opts["Excitation States"]:
    lines.extend([
      "%tddft",
      "  nroots 30",
      "  triplets false",
      "  TDA false",
      "end",
      "",
    ])

  # Optional MD control block.
  if calc_str == "MD":
    _append_md_block(lines, opts)

  # Optional geometry control block.
  if geomcontrol:
    _append_geom_block(lines, calculate, cjson)

  # Optional NEB control.
  if calculate == "Transition State (NEB method)":
    lines.extend([
      "%neb",
      f"  product \"{opts['Filename Base']}_FS.xyz\"",
      "end",
      "",
    ])

  # Optional IRC control.
  if calculate == "Intrinsic Reaction Coordinate":
    lines.extend([
      "%irc",
      "  MaxIter 1024",
      "  Direction both",
      "  InitHess read",
      f"  Hess_Filename \"{opts['Filename Base']}_TS.hess\"",
      "end",
      "",
    ])

  # Coordinate body in xyz format expected by ORCA.
  lines.extend([
    f"* xyz {charge} {multiplicity}",
    "$$coords:___Sxyz$$",
    "*",
    "",
    "",
  ])

  return "\n".join(lines)


def generateInput():
  """Avogadro plugin entry: parse payload and return ORCA input files."""
  stdin_str = sys.stdin.read()
  payload = json.loads(stdin_str)

  inp = generateInputFile(payload["options"], payload["cjson"])
  base_name = payload["options"]["Filename Base"]

  # Output file follows <Filename Base>.inp convention.
  files = [{"filename": f"{base_name}.inp", "contents": inp, "highlightStyles": ["orca-default"]}]
  if DEBUG:
    files.append({"filename": "debug_info", "contents": stdin_str})

  return {
    "files": files,
    "mainFile": f"{base_name}.inp",
  }


def _get_highlight_styles():
  """Define highlight rules for ORCA input preview."""
  rules = []

  rules.append({
    "patterns": [
      {"regexp": "^!.*$"},
      {"regexp": "\\bnprocs\\b"},
      {"regexp": "\\bMaxIter\\b"},
      {"regexp": "\\bxyz\\b"},
      {"regexp": "\\btddft\\b"},
      {"regexp": "\\bmd\\b"},
      {"regexp": "\\bconstraints\\b"},
    ],
    "format": {"preset": "keyword"},
  })

  rules.append({
    "patterns": [
      {"regexp": "^%maxcore\\s+([0-9]+)\\s*$"},
      {"regexp": "^%pal\\s+nprocs\\s+([0-9]+)\\s+end\\s*$"},
      {"regexp": "^\\s*MaxIter\\s+([0-9]+)\\s*$"},
    ],
    "format": {"preset": "property"},
  })

  rules.append({
    "patterns": [{"regexp": "\\b[+-]?[.0-9]+(?:[eEdD][+-]?[.0-9]+)?\\b"}],
    "format": {"preset": "literal"},
  })

  rules.append({
    "patterns": [{"regexp": "^#.*$"}],
    "format": {"preset": "comment"},
  })

  # Keep this rule last so `%block` headers (e.g. `%geom`, `%pal`) and `end`
  # are not overridden by generic keyword rules.
  rules.append({
    "patterns": [
      {"regexp": "^%\\w+"},
      {"regexp": "^\\s*end\\b"},
      {"regexp": "\\bend\\s*$"},
    ],
    "format": {"preset": "title"},
  })

  return [{"style": "orca-default", "rules": rules}]


if __name__ == "__main__":
  """CLI shim for plugin protocol and local testing."""
  parser = argparse.ArgumentParser(f"Generate a {TARGET_NAME} input file.")
  parser.add_argument("--debug", action="store_true")
  parser.add_argument("--print-options", action="store_true")
  parser.add_argument("--generate-input", action="store_true")
  parser.add_argument("--display-name", action="store_true")
  parser.add_argument("--lang", nargs="?", default="en")
  args = vars(parser.parse_args())

  DEBUG = args["debug"]

  # Dispatch selected command.
  if args["display_name"]:
    print(TARGET_NAME)
  if args["print_options"]:
    print(json.dumps(getOptions()))
  elif args["generate_input"]:
    print(json.dumps(generateInput()))
