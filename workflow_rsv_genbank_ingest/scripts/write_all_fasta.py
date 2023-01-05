#!/usr/bin/env python3
# coding: utf-8

import argparse
import csv
import gzip
import pandas as pd
import sys

csv.field_size_limit(sys.maxsize)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--feed', type=str, required=True, help='Feed CSV file')
    parser.add_argument('--metadata-in', type=str, required=True, help='Path to metadata csv file')
    parser.add_argument('--out', type=str, required=True, help='Output fasta file suffix')

    args = parser.parse_args()

    # Load metadata
    metadata = pd.read_csv(args.metadata_in, index_col='Accession ID')

    with open(args.feed, "r", newline="") as fp_in, gzip.open(args.out, "at") as fp_out:

        feed_reader = csv.DictReader(fp_in, delimiter=",", quotechar='"')
        for row in feed_reader:

            # If this entry isn't present in the cleaned metadata, then skip
            accession_id = row['genbank_accession']
            if accession_id not in metadata.index:
                continue

            fp_out.write('>{}\n{}\n'.format(accession_id, row['sequence']))


if __name__ == '__main__':
    main()