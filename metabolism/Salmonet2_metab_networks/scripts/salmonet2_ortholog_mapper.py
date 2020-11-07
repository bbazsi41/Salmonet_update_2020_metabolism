import argparse
import sys
import os
from time import strftime


strain_tax_ids = {
    "ECOLI": 83333,
    "SALA4": 454166,
    "SALAR": 41514,
    "SALBC": 218493,
    "SALCH": 321314,
    "SALDC": 439851,
    "SALEP": 550537,
    "SALG2": 550538,
    "SALHS": 454169,
    "SALNS": 423368,
    "SALPK": 554290,
    "SALPA": 295319,
    "SALPB": 1016998,
    "SALPC": 476213,
    "SALPS": 1081093,
    "SALSV": 439843,
    "SALTI": 90370,
    "SALT1": 588858,
    "SALT4": 909946,
    "SALTD": 568708,
    "SALTY": 99287,
    "SALTS": 216597
}


def parse_args(args):
    help_text = \
        """
        === SalmoNet 2.0 Ortholog Mapper script ===
        
        
        **Description:**

        This script takes a PSI-Mitab 2.7 file, which contains protein-protein
        interactions and mapping the IDs in the first two columns to appropriate
        orholog OMA identifiers.
        
        The script need a mapping file (-m) which contains the OMA IDs of different
        identifiers (Uniprot, Ensembl, etc..). One example file of this is in the
        github repo.
        
        The script also need an ortholog file (-of) which contains the OMA groups of
        the different identifiers. One example file of this is in the github repo.
        
        Taxonomy identifiers, what the user gives, apply to both interaction partners!
        
        
        **Parameters:**
        
        -i, --input-file <path>          : path to an existing mitab file [mandatory]
        
        -m, --mapping-file <path>        : path to an existing ortholog mapping file [optional]
        
        -of, --ortholog-file <path>      : path to an existing ortholog group file [optional]
        
        -t, --tax-id <int>               : taxonomy identifier of the species what the user wants to investigate [mandatory]
        
        -to, --tax-ids-to-output <list>  : comma separated list of taxonomy identifiers of the species what the user wants to put into the output file [mandatory]
        
        -o, --output-file <path>         : path to an output mitab file [mandatory]
        
        
        **Exit codes**
        
        Exit code 1: The specified input file does not exists!
        Exit code 2: The specified mapping file does not exists!
        Exit code 3: The specified ortholog file does not exists!
        """

    parser = argparse.ArgumentParser(description=help_text)

    parser.add_argument("-i", "--input-file",
                        help="<path to an existing mitab file> [mandatory]",
                        type=str,
                        dest="input_file",
                        action="store",
                        required=True)

    parser.add_argument("-m", "--mapping-file",
                        help="<path to an existing ortholog mapping file> [optional]",
                        type=str,
                        dest="mapping_file",
                        action="store",
                        required=False)

    parser.add_argument("-of", "--ortholog-file",
                        help="<path to an existing ortholog group file> [optional]",
                        type=str,
                        dest="ortholog_file",
                        action="store",
                        required=False)

    parser.add_argument("-t", "--tax-id",
                        help="<taxonomy identifier of the species what the user wants to investigate> [mandatory]",
                        type=int,
                        dest="tax_id",
                        action="store",
                        required=True)

    parser.add_argument("-to", "--tax-ids-to-output",
                        help="<comma separated list of taxonomy identifiers of the species what the user wants to put into the output file> [mandatory]",
                        type=str,
                        dest="tax_ids_to_output",
                        action="store",
                        required=False)

    parser.add_argument("-o", "--output-file",
                        help="<path to an output mitab file> [mandatory]",
                        type=str,
                        dest="output_file",
                        action="store",
                        required=True)

    results = parser.parse_args(args)

    return results.input_file, results.mapping_file, results.ortholog_file, results.tax_id, results.tax_ids_to_output, \
           results.output_file


def check_params(input_file, mapping_file, ortholog_file):

    if not os.path.isfile(input_file):
        sys.stderr.write(f'ERROR MESSAGE [{strftime("%H:%M:%S")}]: The specified input file does not exists: {input_file}')
        sys.exit(1)

    if mapping_file:
        if not os.path.isfile(mapping_file):
            sys.stderr.write(f'ERROR MESSAGE [{strftime("%H:%M:%S")}]: The specified mapping file does not exists: {mapping_file}')
            sys.exit(2)

    if ortholog_file:
        if not os.path.isfile(ortholog_file):
            sys.stderr.write(f'ERROR MESSAGE [{strftime("%H:%M:%S")}]: The specified ortholog file does not exists: {ortholog_file}')
            sys.exit(3)


def collect_oma_ids(mapping_file):

    oma_ids = {}

    with open(mapping_file, 'r') as mapping:

        for line in mapping:
            line = line.strip().split('\t')

            if len(line) > 2:
                ids = line[2].strip().split("|")
                other_ids = []
                for id in ids:
                    new_id = id.replace(" ", "")
                    other_ids.append(new_id)
                oma_id = ids[0].strip()
                oma_ids[oma_id] = other_ids[1:-1]

    mapping.close()
    return oma_ids


def collect_orthologs(ortholog_file, tax_id):

    orthologs = {}

    with open(ortholog_file, 'r') as ort:

        for line in ort:

            if "#" in line:
                continue

            if tax_id == 83333 and "ECOLI" not in line:
                continue

            line = line.strip().split('\t')

            oma_group_id = line[0]
            if oma_group_id not in orthologs:
                orthologs[oma_group_id] = []

            for x in range(1, len(line)):
                oma_id = line[x].split("|")[0].split(":")[1].strip()
                orthologs[oma_group_id].append(oma_id)

    ort.close()
    return orthologs


def mapping_oma_ids(oma_ids, interactor_a_ids, interactor_b_ids):

    interactor_a = {}
    interactor_b = {}

    for keys, values in oma_ids.items():

        for ids_a in interactor_a_ids:
            if ids_a in values:
                if keys not in interactor_a:
                    interactor_a[keys] = None

        for ids_b in interactor_b_ids:
            if ids_b in values:
                if keys not in interactor_b:
                    interactor_b[keys] = None

    return interactor_a, interactor_b


def collect_orthologs_for_interactors(interactor, orthologs):

    interactor_orthologs = {}

    for key, values in orthologs.items():

        for oma_id_a in interactor:
            if oma_id_a in values:
                for value in values:
                    if value not in interactor_orthologs:
                        interactor_orthologs[value] = None

    return interactor_orthologs


def write_to_output(line, interactor_a, interactor_b, orthologs, tax_id, out, strain_tax_ids, tax_ids):

    a_orthologs = collect_orthologs_for_interactors(interactor_a, orthologs)
    b_orthologs = collect_orthologs_for_interactors(interactor_b, orthologs)

    for id_a in a_orthologs:
        for id_b in b_orthologs:

            line[0] = f'OMA:{id_a}'
            line[1] = f'OMA:{id_b}'

            if id_a[0:5] != id_b[0:5]:
                continue
            line[9] = f'taxid:{strain_tax_ids[id_a[0:5]]}'
            line[10] = f'taxid:{strain_tax_ids[id_b[0:5]]}'

            if len(tax_ids) != 0:
                if str(strain_tax_ids[id_a[0:5]]) not in tax_ids and str(strain_tax_ids[id_b[0:5]]) not in tax_ids:
                    continue

            out.write("\t".join(line) + '\n')


def collect_other_ids(column, interactor_ids):

    alt_ids = []
    if column != '-':
        alt_ids = column.strip().split("|")

    for alt_id in alt_ids:
        id = alt_id.split(":")[1]
        if id not in interactor_ids:
            interactor_ids.append(id)


def main():

    input_file, mapping_file, ortholog_file, tax_id, tax_ids_to_output, output_file = parse_args(sys.argv[1:])
    tax_ids = []
    if tax_ids_to_output:
        tax_ids = tax_ids_to_output.split(",")

    check_params(input_file, mapping_file, ortholog_file)
    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Parameters are fine, starting...')

    if not mapping_file:
        mapping_file = "../files/Map-SeqNum-ID.txt"
    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Collect OMA IDs from mapping file: {mapping_file}')
    oma_ids = collect_oma_ids(mapping_file)

    if not ortholog_file:
        ortholog_file = "../files/OrthologousGroups.txt"
    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Collect orthologs from ortholog file: {ortholog_file}')
    orthologs = collect_orthologs(ortholog_file, tax_id)

    print(f'MESSAGE [{strftime("%H:%M:%S")}]: Writing results to output file: {output_file}')
    with open(input_file, 'r') as i, open(output_file, 'w') as out:

        for line in i:

            if line.startswith("#"):
                continue

            line = line.strip().split('\t')

            if f"taxid:{tax_id}" not in line[9] and f"taxid:{tax_id}" not in line[10]:
                continue

            interactor_a_ids = []
            interactor_b_ids = []

            interactor_a_ids.append(line[0].split(":")[1])
            interactor_b_ids.append(line[1].split(":")[1])

            collect_other_ids(line[2], interactor_a_ids)
            collect_other_ids(line[3], interactor_b_ids)
            collect_other_ids(line[4], interactor_a_ids)
            collect_other_ids(line[5], interactor_b_ids)

            if line[2] == '-':
                line[2] = f'{line[0]}'
            else:
                line[2] = f'{line[0]}|{line[2]}'

            if line[3] == '-':
                line[3] = f'{line[1]}'
            else:
                line[3] = f'{line[1]}|{line[3]}'

            interactor_a, interactor_b = mapping_oma_ids(oma_ids, interactor_a_ids, interactor_b_ids)

            write_to_output(line, interactor_a, interactor_b, orthologs, tax_id, out, strain_tax_ids, tax_ids)

    print(f'MESSAGE [{strftime("%H:%M:%S")}]: SalmoNet 2.0 Ortholog Mapper script finished successfully!')


if __name__ == '__main__':
    main()
