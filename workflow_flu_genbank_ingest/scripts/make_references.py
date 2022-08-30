#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
from collections import OrderedDict


def read_fasta_file(lines):
    """Read a reference FASTA file. This function is built for many entries in
    the FASTA file but we should be loading in one reference at a time for 
    simplicity's sake. If we need to extract from multiple references then we
    should just run the program multiple times instead.

    Parameters
    ----------
    lines: list of str
        - Output from file.readlines()
    
    Returns
    -------
    entries: dict
        key: Name of the FASTA entry
        value: DNA sequence of the FASTA entry
    """

    entries = OrderedDict()

    # Read sequences
    cur_entry = ""
    cur_seq = ""
    for i, line in enumerate(lines):
        # Strip whitespace
        line = line.strip()

        # If not the name of an entry, add this line to the current sequence
        # (some FASTA files will have multiple lines per sequence)
        if ">" not in line:

            # Skip if empty
            if not line:
                pass
            else:
                cur_seq = cur_seq + line

                # Force to uppercase
                _cur_seq = list(cur_seq.upper())

                # Throw an error for non canonical bases
                # https://www.bioinformatics.org/sms/iupac.html
                # for char in _cur_seq:
                #     if char not in 'ATCGURYSWKMBDHVN':
                #         error_msg = 'Non canonical base: \"{}\" in sequence {} on line {}.'.format(char, line, i)
                #         raise Exception(error_msg)

                # IUPAC also defines gaps as '-' or '.',
                # but the reference shouldn't have gaps.
                # Maybe I can add this in later...

                # Replace the sequence with the edited one
                cur_seq = "".join(_cur_seq)

        # Start of another entry = end of the previous entry
        if ">" in line or i == (len(lines) - 1):
            # Avoid capturing the first one and pushing an empty sequence
            if cur_entry:
                entries[cur_entry] = cur_seq

            # Clear the entry and sequence
            cur_entry = line[1:]
            # Ignore anything past the first whitespace
            if cur_entry:
                cur_entry = cur_entry.split()[0]
            cur_seq = ""

    return entries


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("ref", type=str, help="Path to reference JSON file")

    args = parser.parse_args()

    df = pd.read_json(args.ref)

    # Fix strain naming
    df.loc[df['strain'] == 'A/Hong Kong/1073/99 (H9N2)', 'strain'] = 'A/Hong Kong/1073/99'
    df.loc[df['strain'] == 'A/Goose/Guangdong/1/96(H5N1)', 'strain'] = 'A/goose/Guangdong/1/1996'
    df.loc[df['strain'] == 'B/Lee/40', 'strain'] = 'B/Lee/1940'
    df.loc[df['strain'] == 'A/Puerto Rico/8/1934(H1N1)', 'strain'] = 'A/Puerto Rico/8/1934'

    # Fill in missing serotypes
    df.loc[df['strain'] == 'A/goose/Guangdong/1/1996', 'serotype'] = 'H5N1'
    df.loc[df['strain'] == 'A/Puerto Rico/8/1934', 'serotype'] = 'H1N1'
    df.loc[df['strain'] == 'A/Hong Kong/1073/99', 'serotype'] = 'H9N2'
    df.loc[df['strain'] == 'A/California/07/2009', 'serotype'] = 'H1N1pdm'

    df.loc[df['strain'] == 'C/Ann Arbor/1/50', 'serotype'] = 'C'
    df.loc[df['strain'] == 'B/Lee/1940', 'serotype'] = 'B-lee'
    df.loc[df['genbank_accession'] == 'NC_006306', 'serotype'] = 'C'
    df.loc[df['genbank_accession'] == 'NC_006306', 'strain'] = 'C/Ann Arbor/1/50'
    
    # Fix segments
    df.loc[df['segments'] == 'RNA 1', 'segments'] = '1'

    # Assign segments
    df.loc[df['genbank_accession'] == 'NC_007359', 'segments'] = '3' # PA
    df.loc[df['genbank_accession'] == 'NC_004908', 'segments'] = '4' # HA
    df.loc[df['genbank_accession'] == 'NC_004909', 'segments'] = '6' # NA
    df.loc[df['genbank_accession'] == 'NC_004910', 'segments'] = '1' # PB2
    df.loc[df['genbank_accession'] == 'NC_004911', 'segments'] = '2' # PB1
    df.loc[df['genbank_accession'] == 'NC_004912', 'segments'] = '3' # PA
    df.loc[df['genbank_accession'] == 'NC_007357', 'segments'] = '1' # PB2
    df.loc[df['genbank_accession'] == 'NC_007360', 'segments'] = '5' # NP
    df.loc[df['genbank_accession'] == 'NC_007361', 'segments'] = '6' # NA

    df['segments'] = df['segments'].astype(int)

    df = df.sort_values('strain')

    # print(df.loc[df['serotype'].isna()])

    # Load B references (yamagata, victoria)
    with open('b-vic.fa', 'r') as fp:
        vic = read_fasta_file(fp.readlines())

    with open('b-yam.fa', 'r') as fp:
        yam = read_fasta_file(fp.readlines())
    
    for i in range(8):
        with open('references/{}.fa'.format(i+1), 'w') as fp:
            dff = df.loc[df['segments'] == (i + 1)]
            for _, row in dff.iterrows():
                fp.write('>{}_{}\n'.format(row['serotype'], row['strain']))
                fp.write(row['sequence'] + '\n')

            fp.write('>B-yam_B/Yamagata/16/88\n')
            fp.write(list(yam.items())[i][1] + '\n')

            # Vic missing segment 8
            if i < 7:
                fp.write('>B-vic_B/Victoria/2/87\n')
                fp.write(list(vic.items())[i][1] + '\n')

    


if __name__ == '__main__':
    main()