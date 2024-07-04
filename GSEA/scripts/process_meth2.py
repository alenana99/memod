#!/usr/bin/env python

import sys
import pandas as pd

def find_gene(position, pos_ID_df):
    for _, row in pos_ID_df.iterrows():
        if row['pos_in'] <= position <= row['pos_end']:
            return row['ID_gene']
    return None

def main(pos_ID_path, pmsm_path, output_path):
    pos_ID_df = pd.read_csv(pos_ID_path, sep='\t', header=None)
    pos_ID_df.columns = ['pos_in', 'pos_end', 'ID_gene']
    
    pmsm_df = pd.read_csv(pmsm_path, sep='\t')
    pmsm_df['ID_gene'] = pmsm_df['Position'].apply(lambda pos: find_gene(pos, pos_ID_df))
    pmsm_df['ID_gene'] = pmsm_df['ID_gene'].str.replace('_gene', '')
    
    pmsm_count_df = pmsm_df.groupby(['ID_gene'])['ID_gene'].count().reset_index(name='n_meth')
    pmsm_count_df = pmsm_count_df.sort_values(by='n_meth', ascending=False)
    pmsm_count_df['ID_gene'] = pmsm_count_df['ID_gene'].str.replace('ID=', '')
    
    pmsm_count_df.to_csv(output_path, index=False)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: ./process_methylations.py <pos_ID_path> <pmsm_path> <output_path>")
        sys.exit(1)
    
    pos_ID_path = sys.argv[1]
    pmsm_path = sys.argv[2]
    output_path = sys.argv[3]
    
    main(pos_ID_path, pmsm_path, output_path)
