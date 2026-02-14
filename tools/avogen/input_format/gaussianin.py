#  This source file is part of the Avogadro project.
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").

import argparse
import json
import re
import sys


def getMetaData():
    metaData = {}
    metaData['inputFormat'] = 'gaussian gjf'
    metaData['outputFormat'] = 'xyz'
    metaData['operations'] = ['read']
    metaData['identifier'] = 'Gaussian GJF Format'
    metaData['name'] = 'Gaussian GJF'
    metaData['description'] = "Read coordinates from Gaussian input files by parsing " +\
                              "the geometry block after charge/multiplicity."
    metaData['fileExtensions'] = ['gjf', 'com']
    metaData['mimeTypes'] = ['chemical/x-gaussian-input']
    return metaData


def read():
    lines = sys.stdin.readlines()
    atoms = []
    in_geom_block = False
    charge_mult_re = re.compile(r'^\s*[-+]?\d+\s+[-+]?\d+\s*$')

    for raw in lines:
        line = raw.strip()

        if not in_geom_block:
            if charge_mult_re.match(line):
                in_geom_block = True
            continue

        if not line:
            break

        words = line.split()
        if len(words) < 4:
            continue

        element = words[0]
        try:
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
        except ValueError:
            continue

        atoms.append((element, x, y, z))

    result = str(len(atoms)) + '\n'
    result += 'Generated from Gaussian GJF\n'
    for element, x, y, z in atoms:
        result += '%-3s %12.6f %12.6f %12.6f\n' % (element, x, y, z)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Gaussian GJF file format script.')
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
