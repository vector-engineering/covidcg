#!/usr/bin/env python3
# coding: utf-8

"""Assign influenza serotype

Authors: 
    - David Favela - Vector Engineering Team (dfavela@broadinstitute.org)
    - Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import pandas as pd
import pysam

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--in-bam', type=str, nargs='+', required=True, help='Input BAM file')
    parser.add_argument('--out-csv', type=str, required=True, help='Output CSV file')

    args = parser.parse_args()

    all_serotypes = []
    
    for bamfile in args.in_bam:
        bamfile = pysam.AlignmentFile(bamfile, "r", check_sq=False)  # pylint: disable=no-member
        
        for read in bamfile.fetch(until_eof=True):
            # print(read.query_name, read.reference_name)
            template = read.reference_name.split('/')[0]
            serotype = template.split('_')[0]
            genus = template.split('_')[1]

            accession_id = read.query_name.split('|')[0]

            all_serotypes.append((accession_id, genus, serotype))

        bamfile.close()

    serotype_df = pd.DataFrame.from_records(
        all_serotypes, columns=["Accession ID", "genus", "serotype"]
    )
    serotype_df.to_csv(args.out_csv, index=False)

if __name__ == '__main__':
    main()