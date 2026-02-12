"""
/*****************************************************************************
  This source file is part of the Avogadro project. Version 251027

  This source code is released under the New BSD License, (the "License").

******************************************************************************/
"""

import argparse
import json
import sys


TARGET_NAME = "Gaussian"
DEBUG = False
WARNINGS = []

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


def getOptions():
  user_options = {option: {"type": kind} for option, kind in OPTION_TYPES.items()}
  user_options["tabName"] = "Standard"
  user_options["Title"]["default"] = ""

  user_options["Calculation Type"]["default"] = 1
  user_options["Calculation Type"]["values"] = list(CALC_TYPE)

  user_options["Theory"]["default"] = 1
  user_options["Theory"]["values"] = THEORIES

  user_options["Basis"]["default"] = 6
  user_options["Basis"]["values"] = BASES

  user_options["Filename Base"] = {
    "type": "string",
    "default": "job",
  }

  user_options["Processor Cores"]["default"] = 24
  user_options["Processor Cores"]["minimum"] = 1

  user_options["Memory"]["default"] = 240
  user_options["Memory"]["minimum"] = 1
  user_options["Memory"]["maximum"] = 9999

  user_options["Multiplicity"]["default"] = 1
  user_options["Multiplicity"]["minimum"] = 1
  user_options["Multiplicity"]["maximum"] = 5

  user_options["Charge"]["default"] = 0
  user_options["Charge"]["minimum"] = -9
  user_options["Charge"]["maximum"] = 9

  user_options["Keywords"] = {
    "type": "string",
    "default": "scf=(xqc, maxcycle=1024)",
  }

  alt_options = {
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

  return {"userOptions": [user_options, alt_options], "highlightStyles": _get_highlight_styles()}


def _build_solvation(solvation_type, solvent):
  if "None" in solvent:
    return ""
  if solvation_type == "IEFPCM":
    return f" scrf=(IEFPCM, solvent={solvent})"
  if solvation_type == "SMD":
    return f" scrf=(SMD, solvent={solvent})"
  return ""


def _build_theory_basis(theory, basis):
  if theory in ("AM1", "PM3"):
    WARNINGS.append("Ignoring basis set for semi-empirical calculation.")
    return f"#p {theory}"
  return f"#p {theory}/{basis.replace(' ', '')}"


def generateInputFile(opts):
  title = opts["Title"]
  calculate = opts["Calculation Type"]
  theory = opts["Alternate Theory"] or opts["Theory"]
  basis = opts["Alternate Basis Set"] or opts["Basis"]
  multiplicity = opts["Multiplicity"]
  charge = opts["Charge"]
  output_format = opts["Output Format"]
  checkpoint = opts["Write Checkpoint File"]
  n_cores = opts["Processor Cores"]
  memory_gb = opts["Memory"]
  keywords = opts["Keywords"]
  dispersion = opts["Dispersion Correction"]
  solvation_type = opts["Solvation Type"]
  solvent = opts["Solvation"]

  lines = [
    f"%NProcShared={n_cores}",
    f"%mem={memory_gb}GB",
  ]

  if checkpoint:
    lines.append(f"%Chk={opts['Filename Base']}.chk")

  route = _build_theory_basis(theory, basis)
  if dispersion in ("GD3BJ", "GD3"):
    route += f" em={dispersion}"

  route += _build_solvation(solvation_type, solvent)

  try:
    route += CALC_TYPE[calculate]
  except KeyError as exc:
    raise ValueError(f"Invalid calculation type: {calculate}") from exc

  if output_format == "Molden":
    route += " gfprint pop=full"
  elif output_format == "Molekel":
    route += " gfoldprint pop=full"
  elif output_format != "Standard":
    raise ValueError(f"Invalid output format: {output_format}")

  route += f" {keywords}"

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
  stdin_str = sys.stdin.read()
  payload = json.loads(stdin_str)

  WARNINGS.clear()
  inp = generateInputFile(payload["options"])
  base_name = payload["options"]["Filename Base"]

  result = {
    "mainFile": f"{base_name}.gjf",
    "files": [{"filename": f"{base_name}.gjf", "contents": inp, "highlightStyles": ["gaussian-default"]}],
  }

  if DEBUG:
    result["files"].append({"filename": "debug_info", "contents": stdin_str})

  if WARNINGS:
    result["warnings"] = WARNINGS

  return result


def _get_highlight_styles():
  rules = []

  rules.append({
    "patterns": [{"regexp": "^%[A-Za-z_][A-Za-z0-9_]*"}],
    "format": {"preset": "title"},
  })

  rules.append({
    "patterns": [{"regexp": "^%[A-Za-z_][A-Za-z0-9_]*=([^\\n]+)$"}],
    "format": {"preset": "property"},
  })

  rules.append({
    "patterns": [{"regexp": "^#.*$"}],
    "format": {"preset": "keyword"},
  })

  rules.append({
    "patterns": [{"regexp": "^\\s+[^\\n]*\\|[^\\n]*$"}],
    "format": {"preset": "comment"},
  })

  rules.append({
    "patterns": [{"regexp": "\\b[+-]?[.0-9]+(?:[eEdD][+-]?[.0-9]+)?\\b"}],
    "format": {"preset": "literal"},
  })

  return [{"style": "gaussian-default", "rules": rules}]


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
