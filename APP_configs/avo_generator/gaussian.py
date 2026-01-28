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
targetName = "Gaussian"
debug = False
warnings = []

option_types={
    "Title":"string",
    "Calculation Type":"stringList", 
    "Theory":"stringList",
    "Basis":"stringList", 
    "Processor Cores":"integer",
    "Memory":"integer", 
    "Multiplicity":"integer", 
    "Charge":"integer",
    }

calc_type={
    "Single Point": " SP",
    "Equilibrium Geometry": " Opt Freq",
    "Transition State": " Opt(ts,noeigen,calcfc) Freq"
    }

theories=['PBEPBE', 'PBE1PBE', 'TPSSTPSS', 'TPSSh', 'BLYP', 'B3LYP', 'CAM-B3LYP', 'M06L','M062X', 'WB97XD', 'PM7', 'HF', 'MP2', 'CCSD']

bases=['6-31G(d)', '6-311G(d,p)', '6-31+G(d)', '6-311+G(d,p)', 'LANL2DZ', 'SDD', 
       'Def2SVP', 'Def2TZVP', 'Def2TZVPP', 'Def2QZVPP', 
       'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ']

def getOptions():
    userOptions = {option:{"type":t} for option,t in option_types.items()}
    
    userOptions["tabName"]="Standard"
    userOptions['Title']['default'] = ''

    userOptions['Calculation Type']['default'] = 1
    userOptions['Calculation Type']['values'] = [c for c in calc_type]

    userOptions['Theory']['default'] = 1
    userOptions['Theory']['values'] = theories

    userOptions['Basis']['default'] = 6
    userOptions['Basis']['values'] = bases

    userOptions['Filename Base'] = {}
    userOptions['Filename Base']['type'] = 'string'
    userOptions['Filename Base']['default'] = 'job'

    userOptions['Processor Cores']['default'] = 24
    userOptions['Processor Cores']['minimum'] = 1

    userOptions['Memory']['default'] = 240
    userOptions['Memory']['minimum'] = 1
    userOptions['Memory']['maximum'] = 9999

    userOptions['Multiplicity']['default'] = 1
    userOptions['Multiplicity']['minimum'] = 1
    userOptions['Multiplicity']['maximum'] = 5

    userOptions['Charge']['default'] = 0
    userOptions['Charge']['minimum'] = -9
    userOptions['Charge']['maximum'] = 9

    userOptions['Keywords'] = {}
    userOptions['Keywords']['type'] = 'string'
    userOptions['Keywords']['default'] = 'scf=(xqc, maxcycle=1024)'

    altOptions={f"Alternate {k}":{} for k in ("Theory", "Basis Set")}
    altOptions["tabName"]="Alternate"
    altOptions['Alternate Theory']['type'] = 'string'
    altOptions['Alternate Theory']['default'] = ''

    altOptions['Alternate Basis Set']['type'] = 'string'
    altOptions['Alternate Basis Set']['default'] = ''
    
    altOptions['Dispersion Correction'] = {}
    altOptions['Dispersion Correction']['type'] = 'stringList'
    altOptions['Dispersion Correction']['values'] = ['GD3BJ', 'GD3', 'None']
    altOptions['Dispersion Correction']['default'] = 0

    altOptions['Write Checkpoint File'] = {}
    altOptions['Write Checkpoint File']['type'] = "boolean"
    altOptions['Write Checkpoint File']['default'] = True

    altOptions['Output Format'] = {}
    altOptions['Output Format']['type'] = 'stringList'
    altOptions['Output Format']['values'] = ['Standard', 'Molden', 'Molekel']
    altOptions['Output Format']['default'] = 0
    
    altOptions['Solvation'] = {}
    altOptions['Solvation']['type'] = 'stringList'
    altOptions['Solvation']['values'] = [
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
            "n-Octanol"
            ]
    altOptions['Solvation']['default'] = 0
    
    altOptions['Solvation Type'] = {}
    altOptions['Solvation Type']['type'] = 'stringList'
    altOptions['Solvation Type']['values'] = ['IEFPCM', 'SMD']
    altOptions['Solvation Type']['default'] = 0
    
    # userOptions=[userOptions,altOptions]
    opts = {'userOptions': [userOptions, altOptions]}

    return opts

def generateInputFile(opts):
    # Extract options:
    title = opts['Title']
    calculate = opts['Calculation Type']
    theory = opts['Theory']
    if opts["Alternate Theory"]:
        theory = opts["Alternate Theory"]
    basis = opts['Alternate Basis Set'] or opts['Basis']
    multiplicity = opts['Multiplicity']
    charge = opts['Charge']
    outputFormat = opts['Output Format']
    checkpoint = opts['Write Checkpoint File']
    nCores = opts['Processor Cores']
    keywords = opts['Keywords']
    disp = opts['Dispersion Correction']
    solvtype = opts["Solvation Type"]
    solvent = opts["Solvation"]
    
    
    output = ''

    # Number of cores
    output += f"%NProcShared={nCores}\n"
    output += f"%mem={opts['Memory']}GB\n"

    # Checkpoint
    if checkpoint:
        baseName=opts["Filename Base"]
        output += f'%Chk={baseName}.chk\n'

    # Theory/Basis
    if theory in ('AM1','PM3'):
        output += f'#p {theory}'
        warnings.append('Ignoring basis set for semi-empirical calculation.')
    else:
        output += f"#p {theory}/{basis.replace(' ', '')}" 
    
    if disp in ('GD3BJ','GD3'):
        output += f' em={disp}'
    
    solvation = ""
    if not "None" in opts["Solvation"] and solvtype == "IEFPCM":
        solvation = " scrf=(IEFPCM, solvent=" + solvent + ")"
    elif not "None" in opts["Solvation"] and solvtype == "SMD":
        solvation = " scrf=(SMD, solvent=" + solvent + ")"
    output += f'{solvation}'
    
    # Calculation type
    try:
        output += calc_type[calculate]
    except KeyError:
        raise Exception(f'Invalid calculation type: {calculate}')

    # Output format
    if outputFormat == 'Standard':
        pass
    elif outputFormat == 'Molden':
        output += ' gfprint pop=full'
    elif outputFormat == 'Molekel':
        output += ' gfoldprint pop=full'
    else:
        raise Exception(f'Invalid output format: {outputFormat}')
    
    # Additional keywords
    output += f" {keywords}"

    # Title
    output += f'\n\n {title}\n\n'
    
    # Charge/Multiplicity
    output += f"{charge} {multiplicity}\n"

    # Coordinates
    output += '$$coords:Sxyz$$\n'

    # The gaussian code is irritatingly fickle -- it *will* silently crash if
    # this extra, otherwise unnecessary newline is not present at the end of the
    # file.
    output += '\n'

    return output


def generateInput():
    # Read options from stdin
    stdinStr = sys.stdin.read()

    # Parse the JSON strings
    opts = json.loads(stdinStr)

    # Generate the input file
    inp = generateInputFile(opts['options'])

    # Basename for input files:
    baseName = opts['options']['Filename Base']

    # Prepare the result
    # Specify the main input file. This will be used by MoleQueue to determine
    # the value of the $$inputFileName$$ and $$inputFileBaseName$$ keywords.
    result = {"mainFile": f'{baseName}.gjf'}
    # Input file text -- will appear in the same order in the GUI as they are
    # listed in the array:
    files = [{'filename': f'{baseName}.gjf', 'contents': inp}]
    if debug:
        files.append({'filename': 'debug_info', 'contents': stdinStr})
    result['files'] = files

    if len(warnings) > 0:
        result['warnings'] = warnings

    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Generate a %s input file." % targetName)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--print-options', action='store_true')
    parser.add_argument('--generate-input', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    debug = args['debug']

    if args['display_name']:
        print(targetName)
    if args['print_options']:
        print(json.dumps(getOptions()))
    elif args['generate_input']:
        print(json.dumps(generateInput()))
