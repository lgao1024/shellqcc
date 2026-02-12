"""
/*****************************************************************************
  This source file is part of the Avogadro project. Version 260212

  This source code is released under the New BSD License, (the "License").

******************************************************************************/
"""

import argparse
import json
import sys


TARGET_NAME = "Gaussian"
DEBUG = False

# Option type schema for the primary tab.
OPTION_TYPES = {
  "Title": "string",
  "Calculation Type": "stringList",
  "Theory": "stringList",
  "Basis": "stringList",
  "Processor Cores": "integer",
  "Memory": "integer",
  "Multiplicity": "integer",
  "Charge": "integer",
}

# UI calculation labels mapped to Gaussian route fragments.
CALC_TYPE = {
  "Single Point": " SP",
  "Equilibrium Geometry": " Opt Freq",
  "Transition State": " Opt(ts,noeigen,calcfc) Freq",
}

THEORIES = [
  "PBEPBE",
  "PBE1PBE",
  "TPSSTPSS",
  "TPSSh",
  "BLYP",
  "B3LYP",
  "CAM-B3LYP",
  "M06L",
  "M062X",
  "WB97XD",
  "PM7",
  "HF",
  "MP2",
  "CCSD",
]

BASES = [
  "6-31G(d)",
  "6-311G(d,p)",
  "6-31+G(d)",
  "6-311+G(d,p)",
  "LANL2DZ",
  "SDD",
  "Def2SVP",
  "Def2TZVP",
  "Def2TZVPP",
  "Def2QZVPP",
  "cc-pVDZ",
  "cc-pVTZ",
  "cc-pVQZ",
  "aug-cc-pVDZ",
  "aug-cc-pVTZ",
  "aug-cc-pVQZ",
]


def _build_standard_options():
  """Build the primary Gaussian option tab."""
  # Start from the shared option type map and then assign defaults/ranges.
  user_options = {option: {"type": kind} for option, kind in OPTION_TYPES.items()}
  user_options["tabName"] = "Standard"

  # Hidden title is still kept in the input as a comment header line.
  user_options["Title"].update({"default": "", "hide": False})

  # Calculation route base.
  user_options["Calculation Type"].update({"default": 1, "values": list(CALC_TYPE)})

  # Electronic structure choices.
  user_options["Theory"].update({"default": 1, "values": THEORIES})
  user_options["Basis"].update({"default": 6, "values": BASES})

  # Runtime and file naming controls.
  user_options["Filename Base"] = {"type": "string", "default": "gaussian"}
  user_options["Processor Cores"].update({"default": 24, "minimum": 1})
  user_options["Memory"].update({"default": 240, "minimum": 1, "maximum": 9999})

  # Spin/charge controls.
  user_options["Multiplicity"].update({"default": 1, "minimum": 1, "maximum": 5})
  user_options["Charge"].update({"default": 0, "minimum": -9, "maximum": 9})

  # Free-form route suffix entered by users.
  user_options["Keywords"] = {"type": "string", "default": "scf=(xqc, maxcycle=1024)"}

  return user_options


def _build_alternate_options():
  """Build the secondary tab with override and environment options."""
  return {
    "tabName": "Alternate",
    "Alternate Theory": {
      "type": "string",
      "default": "",
    },
    "Alternate Basis Set": {
      "type": "string",
      "default": "",
    },
    "Dispersion Correction": {
      "type": "stringList",
      "values": ["GD3BJ", "GD3", "None"],
      "default": 0,
    },
    "Write Checkpoint File": {
      "type": "boolean",
      "default": True,
    },
    "Output Format": {
      "type": "stringList",
      "values": ["Standard", "Molden", "Molekel"],
      "default": 0,
    },
    "Solvation": {
      "type": "stringList",
      "values": [
        "None (gas)",
        "Water",
        "Acetonitrile",
        "Acetone",
        "Ethanol",
        "Methanol",
        "CCl4",
        "CH2Cl2",
        "DMSO",
        "Toluene",
        "THF",
        "Toluene",
        "DiethylEther",
        "DiChloroEthane",
        "Benzene",
        "ChloroBenzene",
        "CH3NO2",
        "Heptane",
        "CycloHexane",
        "Aniline",
        "n-Octanol",
      ],
      "default": 0,
    },
    "Solvation Type": {
      "type": "stringList",
      "values": ["IEFPCM", "SMD"],
      "default": 0,
    },
  }


def getOptions():
  """Return Avogadro option schema for Gaussian input generation."""
  standard_options = _build_standard_options()
  alternate_options = _build_alternate_options()

  return {
    "userOptions": [standard_options, alternate_options],
    "highlightStyles": _get_highlight_styles(),
  }


def _build_solvation(solvation_type, solvent):
  """Build the Gaussian SCRF route fragment from solvent settings."""
  # "None (gas)" means explicit gas-phase calculation.
  if "None" in solvent:
    return ""

  if solvation_type == "IEFPCM":
    return f" scrf=(IEFPCM, solvent={solvent})"

  if solvation_type == "SMD":
    return f" scrf=(SMD, solvent={solvent})"

  return ""


def _build_theory_basis(theory, basis, warnings):
  """Build theory/basis route head and collect compatibility warnings."""
  # Semi-empirical methods do not take an explicit basis keyword in Gaussian.
  if theory in ("AM1", "PM3"):
    warnings.append("Ignoring basis set for semi-empirical calculation.")
    return f"#p {theory}"

  basis_clean = str(basis).replace(" ", "")
  return f"#p {theory}/{basis_clean}"


def _build_output_format_suffix(output_format):
  """Map output format selection to extra Gaussian route keywords."""
  if output_format == "Standard":
    return ""

  if output_format == "Molden":
    return " gfprint pop=full"

  if output_format == "Molekel":
    return " gfoldprint pop=full"

  raise ValueError(f"Invalid output format: {output_format}")


def _resolve_theory_and_basis(opts):
  """Resolve alternate overrides for theory and basis when provided."""
  # Alternate values are free-form overrides and take precedence.
  theory = opts["Alternate Theory"] or opts["Theory"]
  basis = opts["Alternate Basis Set"] or opts["Basis"]
  return theory, basis


def _build_route_line(opts, warnings):
  """Compose a complete Gaussian route line from all route-affecting options."""
  calculate = opts["Calculation Type"]
  theory, basis = _resolve_theory_and_basis(opts)

  route = _build_theory_basis(theory, basis, warnings)

  # Optional empirical dispersion term.
  dispersion = opts["Dispersion Correction"]
  if dispersion in ("GD3BJ", "GD3"):
    route += f" em={dispersion}"

  # Optional solvent model term.
  route += _build_solvation(opts["Solvation Type"], opts["Solvation"])

  # Mandatory run type fragment.
  try:
    route += CALC_TYPE[calculate]
  except KeyError as exc:
    raise ValueError(f"Invalid calculation type: {calculate}") from exc

  # Optional output format modifiers.
  route += _build_output_format_suffix(opts["Output Format"])

  # User-provided additional keywords are appended last.
  keywords = str(opts["Keywords"]).strip()
  if keywords:
    route += f" {keywords}"

  return route


def generateInputFile(opts, warnings=None):
  """Generate Gaussian input text from options and return it as a string."""
  warnings = [] if warnings is None else warnings

  title = opts["Title"]
  multiplicity = opts["Multiplicity"]
  charge = opts["Charge"]
  n_cores = opts["Processor Cores"]
  memory_gb = opts["Memory"]

  # Header directives controlling resources and checkpoint IO.
  lines = [
    f"%NProcShared={n_cores}",
    f"%mem={memory_gb}GB",
  ]

  if opts["Write Checkpoint File"]:
    lines.append(f"%Chk={opts['Filename Base']}.chk")

  # Main route line.
  route = _build_route_line(opts, warnings)

  # Geometry/title/body section.
  lines.extend([
    route,
    "",
    f" {title}",
    "",
    f"{charge} {multiplicity}",
    "$$coords:Sxyz$$",
    "",
  ])

  # Gaussian silently fails without a trailing newline.
  return "\n".join(lines)


def generateInput():
  """Avogadro plugin entry: parse payload and return generated Gaussian files."""
  stdin_str = sys.stdin.read()
  payload = json.loads(stdin_str)

  warnings = []
  inp = generateInputFile(payload["options"], warnings=warnings)
  base_name = payload["options"]["Filename Base"]

  result = {
    "mainFile": f"{base_name}.gjf",
    "files": [{"filename": f"{base_name}.gjf", "contents": inp, "highlightStyles": ["gaussian-default"]}],
  }

  if DEBUG:
    result["files"].append({"filename": "debug_info", "contents": stdin_str})

  if warnings:
    result["warnings"] = warnings

  return result


def _get_highlight_styles():
  """Define highlight rules for generated Gaussian input text."""
  rules = []

  # Link0 directives.
  rules.append({
    "patterns": [{"regexp": "^%[A-Za-z_][A-Za-z0-9_]*"}],
    "format": {"preset": "title"},
  })

  # Link0 assigned values.
  rules.append({
    "patterns": [{"regexp": "^%[A-Za-z_][A-Za-z0-9_]*=([^\\n]+)$"}],
    "format": {"preset": "property"},
  })

  # Route line.
  rules.append({
    "patterns": [{"regexp": "^#.*$"}],
    "format": {"preset": "keyword"},
  })

  # Pipe-style inline comments/diagnostics.
  rules.append({
    "patterns": [{"regexp": "^\\s+[^\\n]*\\|[^\\n]*$"}],
    "format": {"preset": "comment"},
  })

  # Numeric literals.
  rules.append({
    "patterns": [{"regexp": "\\b[+-]?[.0-9]+(?:[eEdD][+-]?[.0-9]+)?\\b"}],
    "format": {"preset": "literal"},
  })

  return [{"style": "gaussian-default", "rules": rules}]


if __name__ == "__main__":
  """CLI shim for Avogadro plugin protocol and local debugging."""
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
