#!/usr/bin/env python
# count_meth.py

#import library
import pandas as pd
import argparse
import os

def find_gene(position, pos_ID_df):
    for _, row in pos_ID_df.iterrows():
        if row['pos_in'] <= position <= row['pos_end']:
            return row['ID_gene']
    return None

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process methylation data and gene positions.')
    parser.add_argument('-p', '--pos_id_file', required=True, help='File with the start and end positions of each gene')
    parser.add_argument('-m', '--methylation_file', required=True, help='Methylation data file')
    parser.add_argument('-o', '--output_file', required=True, help='Output file for the results')
    
    args = parser.parse_args()

    # Load the file with the start and end positions of each gene
    pos_ID_df = pd.read_csv(args.pos_id_file, header=None, sep='\t')
    pos_ID_df.columns = ['pos_in', 'pos_end', 'ID_gene']

    # Load methylation data
    pmsm_df = pd.read_csv(args.methylation_file, sep='\t')

    # Find gene for each methylated position
    pmsm_df['ID_gene'] = pmsm_df['Position'].apply(lambda pos: find_gene(pos, pos_ID_df))
    pmsm_df['ID_gene'] = pmsm_df['ID_gene'].str.replace('_gene','')

    # Count the number of methylations per gene
    pmsm_count_df = pmsm_df.groupby(['ID_gene'])['ID_gene'].count().reset_index(name='n_meth')

    # Sort by number of methylations
    pmsm_count_df = pmsm_count_df.sort_values(by='n_meth', ascending=False)
    pmsm_count_df['ID_gene'] = pmsm_count_df['ID_gene'].str.replace('ID=', '')

    # Save the results
    pmsm_count_df.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    main()
