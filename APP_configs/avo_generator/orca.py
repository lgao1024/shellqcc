"""
/******************************************************************************
  This source file is part of the Avogadro project. Version 251027

  This source code is released under the New BSD License, (the "License").
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

    userOptions["Title"] = {
        "type": "string",
        "default": "",
        "toolTip": "Title of the input file",
        "hide": False,
    }

    userOptions["Processor Cores"] = {
        "type": "integer",
        "default": 24,
        "minimum": 1,
    }

    userOptions["Memory"] = {
        "type": "integer",
        "default": 240,
        "minimum": 1,
        "maximum": 999,
    }

    userOptions["Calculation Type"] = {
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
    }

    userOptions["Theory"] = {
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
    }

    userOptions["RI Approximation"] = {
        "type": "stringList",
        "default": 0,
        "values": ["None", "NORI", "RIJK", "RIJONX", "RIJCOSX"],
    }

    userOptions["Dispersion Correction"] = {
        "type": "stringList",
        "default": 0,
        "values": ["None", "D3ZERO", "D3BJ", "D4"],
        "toolTip": "Any added dispersion corrections",
        "hide": False,
    }

    userOptions["Basis"] = {
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
    }

    userOptions["Solvent"] = {
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
    }

    userOptions["Solvation Model"] = {
        "type": "stringList",
        "default": 0,
        "values": ["CPCM", "SMD"],
        "toolTip": "Solvation Method",
    }

    userOptions["Filename Base"] = {
        "type": "string",
        "default": "job",
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
        "default": "SlowConv DefGrid3 TightSCF",
    }
    
    userOptions["Excitation States"] = {
        "type": "boolean",
        "default": False,
    }
    
    #userOptions["AutoAux"] = {
    #    "type": "boolean",
    #    "default": False,
    #    "toolTip": "Automatically select auxiliary basis set",
    #}

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

    opts = {"userOptions": [userOptions, aimdOptions]}
    opts['inputMoleculeFormat'] = 'cjson'

    return opts


def generateInputFile(opts,cjson):
    # Extract options:
    title = opts["Title"]
    calculate = opts["Calculation Type"]
    theory = opts["Theory"]
    basis = opts["Basis"]
    charge = opts["Charge"]
    multiplicity = opts["Multiplicity"]
    nCores = int(opts["Processor Cores"])
    memory = int((opts["Memory"] * 1024) / nCores)
    addopts = opts["Keywords"]
    solvtype = opts["Solvation Model"]
    solvent = opts["Solvent"]
    excitation = opts["Excitation States"]
    # autoaux = opts["AutoAux"]
    disp = opts["Dispersion Correction"]
    ri = opts["RI Approximation"]
    auxbasis = "None"

    rijbasis = {
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

    rijkbasis = {
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

    # Convert to code-specific strings
    calcStr = ""
    geomcontrol = True
    if calculate == "Single Point":
        calcStr = "SP"
        geomcontrol = False
    elif calculate == "Geometry Optimization":
        calcStr = "Opt"
    elif calculate == "Frequencies":
        calcStr = "Opt Freq"
    elif calculate == "Dynamics":
        calcStr = "MD"
    elif calculate == "Transition State":
        calcStr = "OptTS"
    elif calculate == "Transition State (NEB method)":
        calcStr = "NEB-TS Freq"
    elif calculate == "Intrinsic Reaction Coordinate":
        calcStr = "IRC"
    else:
        raise Exception("Unhandled calculation type: %s" % calculate)



    solvation = ""
    if not "None" in opts["Solvent"] and solvtype == "CPCM":
        solvation = "CPCM(" + solvent + ")"
    elif not "None" in opts["Solvent"] and solvtype == "SMD":
        solvation = "CPCM"

    if disp == "None":
        disp = ""
    else:
        disp = " " + disp

    if ri in ["None"]:
        autoaux = False
        ri = ""
    # see https://discuss.avogadro.cc/t/orca-input-generator-does-not-print-nori-when-selected/5489
    elif ri in ["NORI"]:
        autoaux = False
        ri = " " + ri
    else:
        if ri in ["RIJONX", "RIJCOSX"]:
            auxbasis = rijbasis[basis]
        else:
            auxbasis = rijkbasis[basis]
        ri = " " + ri

    #if autoaux == True:
    #    auxbasis = "AutoAux"

    if auxbasis != "None":
        basis = basis + " " + auxbasis

    if "-3c" in theory or "-3C" in theory or "XTB" in theory:
        # -3c composite methods have everything together
        code = f"{calcStr} {theory}{ri}"
    else:
        theory = theory + disp + ri
        # put the pieces together
        code = f"{calcStr} {theory} {basis}"
    if not solvation == "":
        code = f"{code} {solvation}"
    if not addopts == "":
        code = f"{code} {addopts}"



    output = ""

    output += "# " + title + "\n"
    output += "%maxcore " + str(memory) + "\n"
    if nCores > 1:
        output += "%pal nprocs "+ str(nCores) + " end\n"
    
    if geomcontrol == True:
        output += f"! {code}" + " MiniPrint\n"
    else:
        output += f"! {code}" + " NoAutoStart PrintMOs PrintBasis\n"
        output += "# ! MORead\n"
        output += "# %Moinp \"" + str(opts["Filename Base"]) + "_restart.gbw\"\n"
    output += "\n"

    if not "None" in opts["Solvent"] and solvtype == "SMD":
        output += "%cpcm\n"
        output += "   smd true\n"
        output += '   SMDSolvent "' + solvent + '"\n'
        output += "end\n\n"

    if excitation == True:
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

        output += "   Run " + str(opts["AIMD RunTime"]) + " CenterCOM\n"
        output += "end\n\n"



    if geomcontrol == True:
        # check for constraints and frozen atoms in cjson
        output += "%geom\n"
        
        if "constraints" in cjson or "frozen" in cjson["atoms"]:
            output += "   Constraints\n"
            
            # look for bond, angle, torsion constraints
            if "constraints" in cjson:
                # loop through the output
                # e.g. "{ B N1 N2 value C }"
                for constraint in cjson["constraints"]:
                    if len(constraint) == 3:
                        # distance
                        value, atom1, atom2 = constraint
                        output += "      " + ""f"{{ B {atom1} {atom2} {value:.6f} C }}\n"
                    if len(constraint) == 4:
                        # angle
                        value, atom1, atom2, atom3 = constraint
                        output += "      " + f"{{ A {atom1} {atom2} {atom3} {value:.6f} C }}\n"
                    if len(constraint) == 5:
                        # torsion / dihedral
                        value, atom1, atom2, atom3, atom4 = constraint
                        output += "      " + f"{{ D {atom1} {atom2} {atom3} {atom4} {value:.6f} C }}\n"

            # look for frozen atoms
            if "frozen" in cjson["atoms"]:
                # two possibilities - same number of atoms
                # or .. 3*number of atoms
                frozen = cjson["atoms"]["frozen"]
                atomCount = len(cjson["atoms"]["elements"]["number"])
                if len(frozen) == atomCount:
                    # look for 1 or 0
                    for i in range(len(frozen)):
                        if frozen[i] == 1:
                            output += "      " + f"{{ C {i} C }}\n"
                elif len(frozen) == 3 * atomCount:
                    # look for 1 or 0 - x, y, z for each atom
                    for i in range(0, len(frozen), 3):
                        if frozen[i] == 0:
                            output += "      " + f"{{ X {i} C }}\n"
                        if frozen[i + 1] == 0:
                            output += "      " + f"{{ Y {i} C }}\n"
                        if frozen[i + 2] == 0:
                            output += "      " + f"{{ Z {i} C }}\n"

            output += "   end\n"
        if calculate == "Transition State":
            output += "   ReCalc_Hess 100\n"
        output += "   MaxIter 1024\n"
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
    inp = generateInputFile(opts["options"], opts["cjson"])

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
