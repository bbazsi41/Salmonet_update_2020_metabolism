#!/bin/python

import argparse
import sys
import os
from time import strftime


def parse_args(args):
    
    help_text = \
    """
    Salmonet metabolic layer formatter
    """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input-file",
        help="<Path to input file [mandatory]>",
        type=str,
        dest="input_file",
        action="store",
        required=True)
    parser.add_argument("-o", "--output-file",
        help="<path to an output file> [mandatory]",
        type=str,
        dest="output_file",
        action="store",
        required=True)

    results = parser.parse_args(args)

    return results.input_file, results.output_file


def main():
    input_file, output_file = parse_args(sys.argv[1:])
    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Starting merge')
    
    kizar = [
        'P0A263',
        'P0A263',
        'P30705',
        'Q8ZPL4',
    ]

    NET = {}

    files = {
        input_file : 'STM',
        'BioModels_EEI.csv': 'BioModels',
    }

    for fn in files:
        with open(fn) as f:
            for line in f:
                cell = line.strip().split(";")
                if cell[0] in kizar or cell[1] in kizar:
                    ms = cell[2].split(",")
                    if not len(ms) > 2 or not "co-catalysis" in ms:
                        continue
                fe = "locustag:%s\tlocustag:%s" % (cell[0], cell[1])
                re = "locustag:%s\tlocustag:%s" % (cell[1], cell[0])
                if fe in NET:
                    NET[fe].append(files[fn])
                elif re in NET:
                    NET[re].append(files[fn])
                else:
                    NET[fe] = []
                    NET[fe].append(files[fn])

    with open("%s_out" % output_file, "w") as f:
        for e in NET:
            f.write("%s\t%s\n" % (e, ",".join(NET[e])))

if __name__ == '__main__':
    main()