"""
/******************************************************************************
  This source file is part of the Avogadro project.
  
  This source code is released under the New BSD License, (the "License").
  
  version 251023
******************************************************************************/
"""

import argparse
import json
import sys

# Some globals:
targetName = "ORCA"
debug = False


def getOptions():
    userOptions = {"tabName": "Basic"}

    userOptions["Filename Base"] = {
        "type": "string",
        "default": "job",
    }

    userOptions["Processor Cores"] = {
        "type": "integer",
        "default": 24,
        "minimum": 1,
    }

    userOptions["Calculation Type"] = {
        "type": "stringList",
        "default": 1,
        "values": [
            "Single Point",
            "Geometry Optimization",
            "Transition State",
            "Transition State Scan",
            "Transition State (NEB method)",
            "Intrinsic Reaction Coordinate",
            "Excitation States", 
            "Molecular Dynamics",
        ],
        "toolTip": "Type of calculation to perform",
    }

    userOptions["Theory"] = {
        "type": "stringList",
        "default": 3,
        "values": [
            "XTB",
            "HF-3c",
            "B97-3c",
            'r2SCAN-3c',
            'PBEh-3c',
            "PBE",
            "revPBE",
            "BLYP",
            "M06L",
            'TPSS',
            "PBE0",
            "M062X",
            "B3LYP",
            'B3LYP/G',
            'CAM-B3LYP',
            "wB97X-D4",
            'wB97M-V',
            'TPSSh',
            'TPSS0',
            'PWPB95',
            'wB97X-2',
            "HF",
            "MP2",
            "CCSD",
            "CCSD(T)",
        ],
    }

    userOptions["Basis"] = {
        "type": "stringList",
        "default": 0,
        "values": [
            'None',
            "def2-SVP",
            "def2-TZVP",
            "def2-TZVPP",
            "def2-QZVPP",
            "cc-pVDZ",
            "cc-pVTZ",
            "cc-pVQZ",
            "6-31G(d)",
            '6-311G(d,p)',
            "def2-TZVPPD",
            "def2-QZVPPD",
            "ma-def2-SVP",
            "ma-def2-TZVP",
            "ma-def2-QZVP",
            "aug-cc-pVDZ",
            "aug-cc-pVTZ",
            "aug-cc-pVQZ",
        ],
        "toolTip": "Gaussian basis set",
    }

    userOptions["Charge"] = {
        "type": "integer",
        "default": 0,
        "minimum": -9,
        "maximum": 9,
        "toolTip": "Total charge of the system",
    }

    userOptions["Multiplicity"] = {
        "type": "integer",
        "default": 1,
        "minimum": 1,
        "maximum": 6,
        "toolTip": "Total spin multiplicity of the system",
    }

    userOptions["Keywords"] = {
        "type": "string",
        "default": "noautostart miniprint nopop slowconv defgrid3 tightSCF",
    }

    userOptions["Memory"] = {
        "type": "integer",
        "default": 240,
        "minimum": 1,
        "maximum": 999,
    }
    
    userOptions["Solvation"] = {
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
        "toolTip": "Solvent Model",
    }

    advOptions = {"tabName": "Advanced"}

    advOptions["Constrain Atoms"] = {
        "type": "string",
        "default": "",
    }

    advOptions["Scan"] = {
        "type": "string",
        "default": "",
    }

    advOptions["Hybrid Hess"] = {
        "type": "string",
        "default": "",
    }

    advOptions["Dispersion Correction"] = {
        "type": "stringList",
        "default": 0,
        "values": ["None", "D3ZERO", "D3BJ", "D4"],
        "toolTip": "Any added dispersion corrections",
        "hide": False,
    }

    advOptions["RI Approximation"] = {
        "type": "stringList",
        "default": 0,
        "values": ["Default", "NORI", "RIJ", "RIJK", "RIJONX", "RIJCOSX"],
    }

    advOptions["Solvation Type"] = {
        "type": "stringList",
        "default": 0,
        "values": ["CPCM", "SMD"],
        "toolTip": "Solvent model",
    }

    aimdOptions = {"tabName": "Dynamics"}
    aimdOptions["AIMD TimeStep"] = {
        "type": "string",
        "default": "0.5_fs",
    }
    aimdOptions["AIMD Initvel"] = {
        "type": "integer",
        "default": 300,
        "maximum": 9999,
        "minimum": 0,
    }
    aimdOptions["AIMD Thermostat Type"] = {
        "type": "stringList",
        "default": 1,
        "values": ["Berendsen", "CSVR", "NHC", "None"],
    }
    aimdOptions["AIMD Thermostat Temp"] = {
        "type": "integer",
        "default": 300,
        "maximum": 1000,
        "minimum": 0,
    }
    aimdOptions["AIMD Thermostat Time"] = {
        "type": "string",
        "default": "30_fs",
    }
    aimdOptions["AIMD Restart"] = {
        "type": "stringList",
        "default": 0,
        "values": ["Not Restart", "Restart", "Restart IfExists"],
    }
    aimdOptions["AIMD Output Trajectory"] = {
        "type": "stringList",
        "default": 0,
        "values": ["Position", "Position, Velocity and Force", "All", "None"],
    }
    aimdOptions["AIMD RunTime"] = {
        "type": "integer",
        "default": 2000,
        "maximum": 999999999,
        "minimum": 0,
    }

    opts = {"userOptions": [userOptions, advOptions, aimdOptions]}

    return opts


def generateInputFile(opts):
    calculate = opts["Calculation Type"]
    theory = opts["Theory"]
    basis = opts["Basis"]
    charge = opts["Charge"]
    multiplicity = opts["Multiplicity"]
    nCores = int(opts["Processor Cores"])
    memory = int((opts["Memory"] * 1024) / nCores)
    addopts = opts["Keywords"]
    solvtype = opts["Solvation Type"]
    solvent = opts["Solvation"]
    disp = opts["Dispersion Correction"]
    ri = opts["RI Approximation"]
    rijbasis = {
        "None": "None",
        "def2-SVP": "Def2/J",
        "def2-TZVP": "Def2/J",
        "def2-TZVPP": "Def2/J",
        "def2-QZVPP": "Def2/J",
        "cc-pVDZ": "Def2/J",
        "cc-pVTZ": "Def2/J",
        "cc-pVQZ": "Def2/J",
        "6-31G(d)": "AutoAux",
        '6-311G(d,p)': 'AutoAux',
        "def2-TZVPPD": "AutoAux",
        "def2-QZVPPD": "AutoAux",
        "ma-def2-SVP": "AutoAux",
        "ma-def2-TZVP": "AutoAux",
        "ma-def2-QZVP": "AutoAux",
        "aug-cc-pVDZ": "AutoAux",
        "aug-cc-pVTZ": "AutoAux",
        "aug-cc-pVQZ": "AutoAux",
    }
    rijkbasis = {
        "None": "None",
        "def2-SVP": "Def2/JK",
        "def2-TZVP": "Def2/JK",
        "def2-TZVPP": "Def2/JK",
        "def2-QZVPP": "Def2/JK",
        "cc-pVDZ": "cc-pVDZ/JK",
        "cc-pVTZ": "cc-pVTZ/JK",
        "cc-pVQZ": "cc-pVQZ/JK",
        "6-31G(d)": "AutoAux",
        '6-311G(d,p)': 'AutoAux',
        "def2-TZVPPD": "aug-cc-pVTZ/JK",
        "def2-QZVPPD": "aug-cc-pVQZ/JK",
        "ma-def2-SVP": "aug-cc-pVDZ/JK",
        "ma-def2-TZVP": "aug-cc-pVTZ/JK",
        "ma-def2-QZVP": "aug-cc-pVQZ/JK",
        "aug-cc-pVDZ": "aug-cc-pVDZ/JK",
        "aug-cc-pVTZ": "aug-cc-pVTZ/JK",
        "aug-cc-pVQZ": "aug-cc-pVQZ/JK",
    }

    # Convert to code-specific strings
    calcStr = ""
    geomcontrol = True
    scancontrol = False
    if calculate == "Single Point":
        calcStr = "SP"
        geomcontrol = False
    elif calculate == "Geometry Optimization":
        calcStr = "Opt Freq"
        scancontrol = True
    elif calculate == "Transition State":
        calcStr = "OptTS Freq"
    elif calculate == "Transition State Scan":
        calcStr = "ScanTS Freq"
        scancontrol = True
    elif calculate == "Transition State (NEB method)":
        calcStr = "NEB-TS Freq"
    elif calculate == "Intrinsic Reaction Coordinate":
        calcStr = "IRC"
    elif calculate == "Excitation States":
        calcStr = "SP"
        geomcontrol = False
    elif calculate == "Molecular Dynamics":
        calcStr = "MD"
        geomcontrol = False
    else:
        raise Exception("Unhandled calculation type: %s" % calculate)

    if basis == 'None':
        basis = ""
        if ri == "NORI":
            basis = " " + ri 
    else:
        if ri in ["RIJ", "RIJONX", "RIJCOSX"]:
            basis = " " + ri + " " + basis + " " + rijbasis[basis]
        elif ri == "RIJK":
            basis = " " + ri + " " + basis + " " + rijkbasis[basis]
        elif ri == "NORI":
            basis = " " + ri + " " + basis
        else:
            basis = " " + basis

    if disp == "None":
        theory = theory + basis
    else:
        theory = theory + " " + disp + basis
    
    solvation = ""
    if not "None" in opts["Solvation"] and solvtype == "CPCM":
        solvation = "CPCM(" + solvent + ")"
    elif not "None" in opts["Solvation"] and solvtype == "SMD":
        solvation = "CPCM"

    code = f"{calcStr} {theory}"
    if not solvation == "":
        code = f"{code} {solvation}"
    if not addopts == "":
        code = f"{code} {addopts}"




    output = ""
    output += "%maxcore " + str(memory) + "\n"
    output += "%pal nprocs "+ str(nCores) + " end\n"
    output += f"! {code}\n"
    if calculate == "Single Point":
        output += "# ! MORead\n"
        output += "# %Moinp \"" + str(opts["Filename Base"]) + "_restart.gbw\"\n"
    output += "\n"

    if not "None" in opts["Solvation"] and solvtype == "SMD":
        output += "%cpcm\n"
        output += "   smd true\n"
        output += '   SMDSolvent "' + solvent + '"\n'
        output += "end\n\n"


    if geomcontrol == True :
        output += "%geom\n"

        if str(opts["Constrain Atoms"]) != "":
            output += "  Constraints\n"
            for Atoms in str(opts["Constrain Atoms"]).split(";"):
                output += "    { " + str(Atoms) + " C}\n"
            output += "  end\n"
        
        if str(opts["Scan"]) != "" and scancontrol == True:
            output += "  Scan\n"
            output += "    " + str(opts["Scan"]) + "\n"
            output += "  end\n"

        if str(opts["Hybrid Hess"]) != "":
            output += "  Hybrid_Hess\n"
            output += "    { " + str(opts["Hybrid Hess"]) + " }\n"
            output += "  end\n"

        if calculate == "Transition State":
            output += "  ReCalc_Hess 50\n"
        
        if calculate == "Transition State Scan":
            output += "  FullScan true\n"

        output += "  MaxIter 1024\n"
        output += "end\n\n"

    if calculate == "Transition State (NEB method)":
        output += "%neb\n"
        output += "  product \"" + str(opts["Filename Base"]) + "_FS.xyz\"\n"
        output += "end\n\n"

    if calculate == "Intrinsic Reaction Coordinate":
        output += "%irc\n"
        output += "  MaxIter 1024\n"
        output += "  Direction both\n"
        output += "  InitHess read\n"
        output += "  Hess_Filename \"" + str(opts["Filename Base"]) + "_TS.hess\"\n"
        output += "end\n\n"
        
    if calculate == "Excitation States":
        output += "%tddft\n"
        output += "  nroots 30\n"
        output += "  triplets false\n"
        output += "  TDA false\n"
        output += "end\n\n"

    if calcStr == "MD":
        output += "%md\n"
        output += "   Timestep " + opts["AIMD TimeStep"] + "\n"

        if str(opts["AIMD Restart"]) != "Not Restart":
            output += "   Initvel " + str(opts["AIMD Initvel"]) + "_K No_Overwrite\n"
        else:
            output += "   Initvel " + str(opts["AIMD Initvel"]) + "_K\n"

        output += (
            "   Thermostat "
            + str(opts["AIMD Thermostat Type"]) + " "
            + str(opts["AIMD Thermostat Temp"])
            + "_K Timecon "
            + opts["AIMD Thermostat Time"]
            + "\n"
        )

        if str(opts["AIMD Output Trajectory"]) != "None":
            output += '   Dump Position Stride 1 Filename "position.xyz"\n'
            if str(opts["AIMD Output Trajectory"]) != "Position":
                output += '   Dump Force Stride 1 Filename "force.xyz"\n'
                output += '   Dump Velocity Stride 1 Filename "velocity.xyz"\n'
                if str(opts["AIMD Output Trajectory"]) != "Position, Velocity and Force":
                    output += '   Dump GBW Stride 10 Filename "wfc"\n'
                    output += '   Dump EnGrad Stride 10 Filename "EnGrad"\n'

        if str(opts["AIMD Restart"]) != "Not Restart":
            output += "   " + str(opts["AIMD Restart"]) + "\n"

        output += "   run " + str(opts["AIMD RunTime"]) + " CenterCOM\n"
        output += "end\n\n"

    output += f"* xyz {charge} {multiplicity}\n"
    output += "$$coords:___Sxyz$$\n"
    output += "*\n\n\n"

    return output


def generateInput():
    # Read options from stdin
    stdinStr = sys.stdin.read()

    # Parse the JSON strings
    opts = json.loads(stdinStr)

    # Generate the input file
    inp = generateInputFile(opts["options"])

    # Basename for input files:
    baseName = opts["options"]["Filename Base"]

    # Prepare the result
    result = {}
    # Input file text -- will appear in the same order in the GUI as they are
    # listed in the array:
    files = []
    files.append({"filename": "%s.inp" % baseName, "contents": inp})
    if debug:
        files.append({"filename": "debug_info", "contents": stdinStr})
    result["files"] = files
    # Specify the main input file. This will be used by MoleQueue to determine
    # the value of the $$inputFileName$$ and $$inputFileBaseName$$ keywords.
    result["mainFile"] = "%s.inp" % baseName
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Generate a %s input file." % targetName)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--print-options", action="store_true")
    parser.add_argument("--generate-input", action="store_true")
    parser.add_argument("--display-name", action="store_true")
    parser.add_argument("--lang", nargs="?", default="en")
    args = vars(parser.parse_args())

    debug = args["debug"]

    if args["display_name"]:
        print(targetName)
    if args["print_options"]:
        print(json.dumps(getOptions()))
    elif args["generate_input"]:
        print(json.dumps(generateInput()))
