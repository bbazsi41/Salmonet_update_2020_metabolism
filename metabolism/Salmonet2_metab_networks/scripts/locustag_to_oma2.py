import os
for filename in os.listdir('../results/'):
    strain_ints = []
    if filename.startswith('oma_'):
        with open(f'../results/{filename}', 'r') as f:
            for lines in f:
                cells = lines.strip().split('\t')
                ids = cells[0]+'\t'+cells[1]
                strain_ints.append(ids)


        all_interactions = []
        with open('all_met_interactions.tsv', 'r') as f:
            for lines in f:
                cells = lines.strip().split('\t')
                ids2 = cells[0]+'\t'+cells[1]
                all_interactions.append(ids2)

        intmatrix = []
        for i in all_interactions:
            if i in strain_ints:
                #print(i+'\t'+'1')
                intmatrix.append(i+'\t'+'1')
            elif i not in strain_ints:
                #print(i)
                intmatrix.append(i+'\t'+'0')

        with open(f'../results/metnet_{filename}', 'w') as outfile:
            outfile.write(f"node1\tnode2\tPresence\tSource\n")
            for j in intmatrix:
                outfile.write(f"{j}\tPMID:31455646\n")

import pandas as pd
import os


df = pd.read_csv('../results/metnet_oma_1016998', sep = '\t')

for filename in os.listdir('../results/'):
    if filename.startswith('metnet_'):
        strain = filename.split('_')
        df2 = pd.read_csv(f'../results/{filename}', sep = '\t')
        df[f'{strain[2]}'] = df2['Presence']

df.to_csv('metnet_matrix.tsv', sep='\t', index=False)