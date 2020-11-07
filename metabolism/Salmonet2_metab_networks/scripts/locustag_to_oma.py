import os
orthologs = {}

with open('OrthologousGroups.txt', 'r') as infile:
    for line in infile:

        if '#' in line:
            continue

        line = line.strip().split('\t')

        oma_group_id = line[0]

        if oma_group_id not in orthologs:
            orthologs[oma_group_id] = []

        for x in range(1, len(line)):
            lt_id = line[x].split(" | ")
            omaID = lt_id[0].split(':')
            # if 'ECOLI' in lt_id[0] or 'SALTS' in lt_id[0]:
            orthologs[oma_group_id].append(lt_id[2])

def invert_dict(d): 
    inverse = dict() 
    for key in d: 
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse: 
                # If not create a new list
                inverse[item] = [key] 
            else: 
                inverse[item].append(key) 
    return inverse

locustag2oma = invert_dict(orthologs)

for filename in os.listdir('../results/'):
    with open(f'../results/{filename}', 'r') as f:
        with open(f'../results/oma_{filename}', 'w') as outfile:
            for lines in f:
                cells = lines.strip().split('\t')
                n1 = cells[0].split(':')
                n2 = cells[1].split(':')
                try:
                    outfile.write(f"{''.join(locustag2oma[n1[1]])}\t{''.join(locustag2oma[n2[1]])}\t{''.join(locustag2oma[n1[1]])}-{''.join(locustag2oma[n2[1]])}\n")
                except IndexError:
                    continue

