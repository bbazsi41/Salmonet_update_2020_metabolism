import pandas as pd
import os


df = pd.read_csv('../results/metnet_oma_1016998', sep = '\t')

for filename in os.listdir('../results/'):
    if filename.startswith('metnet_'):
        strain = filename.split('_')
        df2 = pd.read_csv(f'../results/{filename}', sep = '\t')
        df[f'{strain[2]}'] = df2['Presence']

df.to_csv('metnet_matrix.tsv', sep='\t', index=False)