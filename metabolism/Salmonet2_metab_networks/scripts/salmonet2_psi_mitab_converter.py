import argparse
import sys
import os
from time import strftime


def parse_args(args):
    help_text = \
        """
        === SalmoNet 2.0 PSI-Mitab Converter script ===
        
        
        **Description:**
        
        This script takes a TSV or CSV file, which contains protein-protein
        interactions and converts it to PSI-MITAB 2.7 file format.
        
        The script will delete duplicates!
        
        It does not matter whether the input file has header or not!
        
        IMPORTANT: In this script the numbering of the columns starts with 0!!!
        
        
        **Parameters:**
        
        -i, --input-file <path>          : path to an existing csv/tsv file [mandatory]
        
        -hea, --header <boolen>          : if this parameter is given, then the input file has a header, default: NULL [optional]
        
        -t, --tax-id <int>               : taxonomy identifier of interactor A and interactor B [mandatory]
        
        -nc, --needed-columns <list>     : comma separated list of column numbers, which columns does the user want
                                           to integrate to the output mitab file [mandatory]
        
        -mc, --mitab-columns <list>      : comma separated list of column numbers, in which column of the mitab file
                                           does the user want to keep the informations [mandatory]
        
        -o, --output-file <path>         : path to an output mitab file [mandatory]
        
        
        **Exit codes**
        
        Exit code 1: The specified input file does not exists!
        Exit code 2: The number of the needed_columns parameter and the number of the mitab_columns parameter must be equal!
        """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input-file",
                        help="<path to an existing csv/tsv file> [mandatory]",
                        type=str,
                        dest="input_file",
                        action="store",
                        required=True)

    parser.add_argument("-hea", "--header",
                        help="<if this parameter is given, then the input file has a header, default: NULL> [optional]",
                        dest="header",
                        action="store_true",
                        default=None)

    parser.add_argument("-t", "--tax-id",
                        help="<taxonomy identifier of interactor A and interactor B> [mandatory]",
                        type=int,
                        dest="tax_id",
                        action="store",
                        required=True)

    parser.add_argument("-nc", "--needed-columns",
                        help="<comma separated list of column numbers, which columns does the user want to integrate to"
                             "the output mitab file> [mandatory]",
                        type=str,
                        dest="needed_columns",
                        action="store",
                        required=True)

    parser.add_argument("-mc", "--mitab-columns",
                        help="<comma separated list of column numbers, in which column of the mitab file does the user"
                             "want to keep the informations> [mandatory]",
                        type=str,
                        dest="mitab_columns",
                        action="store",
                        required=True)

    parser.add_argument("-o", "--output-file",
                        help="<path to an output mitab file> [mandatory]",
                        type=str,
                        dest="output_file",
                        action="store",
                        required=True)

    results = parser.parse_args(args)

    return results.input_file, results.header, results.tax_id, results.needed_columns, results.mitab_columns,\
           results.output_file


def check_params(input_file, needed_columns_numbers, mitab_columns_numbers):

    if not os.path.isfile(input_file):
        sys.stderr.write(f'ERROR MESSAGE [{strftime("%H:%M:%S")}]: The specified input file does not exists: {input_file}')
        sys.exit(1)

    if len(needed_columns_numbers) != len(mitab_columns_numbers):
        sys.stderr.write(f'ERROR MESSAGE [{strftime("%H:%M:%S")}]: The number of the needed_columns parameter and '
                         f'the number of the mitab_columns parameter must be equal:'
                         f'{len(needed_columns_numbers)} != {len(mitab_columns_numbers)}')
        sys.exit(2)


def write_to_output(line, out, needed_columns, mitab_columns, tax_id):

    into_mitab = {}

    for r in range(0, 42):
        into_mitab[r] = "-"
        if r == 9 or r == 10:
            into_mitab[r] = f"taxid:{tax_id}"

    for x in range(0, len(needed_columns)):
        into_mitab[mitab_columns[x]] = line[needed_columns[x]]

    for keys, value in into_mitab.items():
        if keys == 41:
            out.write(value)
        else:
            out.write(value + '\t')
    out.write('\n')


def main():

    input_file, header, tax_id, needed_columns, mitab_columns, output_file = parse_args(sys.argv[1:])

    needed_columns_number = needed_columns.split(",")
    needed_columns_numbers = []
    for m in needed_columns_number:
        needed_columns_numbers.append(int(m))
    mitab_columns_number = mitab_columns.split(",")
    mitab_columns_numbers = []
    for n in mitab_columns_number:
        mitab_columns_numbers.append(int(n))

    check_params(input_file, needed_columns_numbers, mitab_columns_numbers)
    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Parameters are fine, starting...')

    abspath_output_file = os.path.abspath(output_file)

    with open(input_file, 'r') as i, open(output_file, 'w') as out:

        if header:
            i.readline()

        avoid_duplicates = set()

        print(f'MESSAGE [{strftime("%H:%M:%S")}]: Write results to the output mitab file: {abspath_output_file}')
        for line in i:
            line = line.strip().split('\t')
            interactions = (f'{line[0]}', f'{line[1]}')
            if interactions not in avoid_duplicates:
                avoid_duplicates.add(interactions)
                write_to_output(line, out, needed_columns_numbers, mitab_columns_numbers, tax_id)

    print(f'MESSAGE [{strftime("%H:%M:%S")}]: SalmoNet2 PSI-Mitab Converter script finished successfully!')


if __name__ == '__main__':
    main()
