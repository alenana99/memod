#process meth
#!/usr/bin/env python

## - Libraries - ##
import os
import sys
import argparse
import pandas as pd
from pathlib import Path

## - Functions - ##
def parse_gff(gff_path):
    pos_ID_data = []
    with open(gff_path, 'r') as gff_file:
        for line in gff_file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if parts[2] == 'gene':
                pos_in = int(parts[3])
                pos_end = int(parts[4])
                attributes = parts[8]
                id_gene = attributes.split(';')[0].replace('ID=', '')
                pos_ID_data.append((pos_in, pos_end, id_gene))
    pos_ID_df = pd.DataFrame(pos_ID_data, columns=['pos_in', 'pos_end', 'ID_gene'])
    return pos_ID_df

def parse_tsv(tsv_path):
    L_df = pd.read_csv(tsv_path, sep='\t', usecols=[3, 4, 5, 18])
    L_df.columns = ['Position', 'Modification', 'Strand', 'Motif']
    return L_df

def find_gene(position, pos_ID_df):
    for _, row in pos_ID_df.iterrows():
        if row['pos_in'] <= position <= row['pos_end']:
            return row['ID_gene']
    return None

def process_methylations(gff_path, tsv_path, output_path):
    pos_ID_df = parse_gff(gff_path)
    L_df = parse_tsv(tsv_path)
    
    L_df['ID_gene'] = L_df['Position'].apply(lambda pos: find_gene(pos, pos_ID_df))
    L_df['ID_gene'] = L_df['ID_gene'].str.replace('_gene', '')
    
    L_count_df = L_df.groupby(['ID_gene'])['ID_gene'].count().reset_index(name='n_meth')
    L_count_df = L_count_df.sort_values(by='n_meth', ascending=False)
    L_count_df['ID_gene'] = L_count_df['ID_gene'].str.replace('ID=', '')
    
    L_count_df.to_csv(output_path, index=False)

def main():
    # Argument parser setup
    parser = argparse.ArgumentParser(
        prog="Methylation Processor",
        description="Process methylation data",
        epilog="Version 1.0"
    )
    
    parser.add_argument('gff_path', help='Path to the GFF3 annotation file')
    parser.add_argument('tsv_path', help='Path to the TSV file from MicrobeMod')
    parser.add_argument('output_path', help='Path to the output file')
    
    args = parser.parse_args()
    
    # Process methylations
    process_methylations(args.gff_path, args.tsv_path, args.output_path)

if __name__ == '__main__':
    main()
