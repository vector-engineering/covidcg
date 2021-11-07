# coding: utf-8

"""Load SNP data, create SNP signatures

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import io
import numpy as np
import pandas as pd

from pathlib import Path


def extract_ids(fasta_file):
    """Extract Accession IDs (entry names) from a fasta file

    Parameters
    ----------
    fasta_file: str

    Returns
    -------
    out: list of tuples, (Accession ID, date)

    """

    out = []

    # Read sequences
    cur_entry = ""
    cur_seq = ""

    # Get the date from the fasta file name, as a string
    file_date = Path(fasta_file).name.replace(".fa.gz", "")

    with gzip.open(fasta_file, "rt") as fp:
        lines = fp.readlines()
        for i, line in enumerate(lines):
            # Strip whitespace
            line = line.strip()

            # Ignore empty lines that aren't the last line
            if not line and i < (len(lines) - 1):
                continue

            # If not the name of an entry, add this line to the current sequence
            # (some FASTA files will have multiple lines per sequence)
            if line[0] != ">":
                cur_seq = cur_seq + line

            # Start of another entry = end of the previous entry
            if line[0] == ">" or i == (len(lines) - 1):
                # Avoid capturing the first one and pushing an empty sequence
                if cur_entry:
                    out.append((cur_entry, file_date,))

                # Clear the entry and sequence
                cur_entry = line[1:]
                # Ignore anything past the first whitespace
                if cur_entry:
                    cur_entry = cur_entry.split()[0]
                cur_seq = ""

    # print("Read {} entries for file {}".format(len(out), fasta_file))

    return out


def process_snps(
    processed_fasta_files,
    snp_files,
    # SNPs must occur at least this many times to pass filters
    count_threshold=3,
    mode="dna",  # dna, gene_aa, protein_aa
):

    manifest = []
    for fasta_file in sorted(Path(processed_fasta_files).glob("*.fa.gz")):
        manifest.extend(extract_ids(fasta_file))
    manifest = pd.DataFrame.from_records(manifest, columns=["Accession ID", "date"])
    pruned_manifest = manifest.drop_duplicates(["Accession ID"], keep="last")

    # Dump all SNP chunks into a text buffer
    snp_df_io = io.StringIO()
    for i, chunk in enumerate(snp_files):
        file_date = Path(chunk).name.replace("_" + mode + "_snp.csv", "")
        with open(chunk, "r") as fp_in:
            # Write dates, so we can remove duplicate sequences
            # and default to the mutations of the latest sequence, by date
            for j, line in enumerate(fp_in):
                # Write the header of the first file
                if i == 0 and j == 0:
                    snp_df_io.write(line.strip() + ",date\n")
                # Or write any line that's not the header
                # (to avoid writing the header more than once)
                elif j > 0:
                    snp_df_io.write(line.strip() + "," + file_date + "\n")

    # Read the buffer into a dataframe, then discard the buffer
    snp_df_io.seek(0)
    snp_df = pd.read_csv(snp_df_io)
    snp_df_io.close()

    # --------------------------
    # Remove duplicate sequences
    # --------------------------

    snp_df = pruned_manifest.set_index(["Accession ID", "date"]).join(
        snp_df.set_index(["Accession ID", "date"]), how="left"
    )
    snp_df.reset_index(inplace=True)

    # Set aside sequences with no mutations
    # We'll add them back in later
    no_snp_seqs = snp_df.loc[snp_df["pos"].isna(), "Accession ID"].values

    # Remove sequences with no mutations before counting mutations
    snp_df = snp_df.loc[~snp_df["pos"].isna()]
    snp_df.loc[:, "pos"] = snp_df["pos"].astype(int)
    # Replace NaNs in the 'ref' and 'alt' column with '-'
    snp_df.fillna("-", inplace=True)

    groupby_cols = ["pos", "ref", "alt"]
    if mode == "gene_aa":
        groupby_cols.insert(0, "gene")
    elif mode == "protein_aa":
        groupby_cols.insert(0, "protein")

    # Generate SNP strings for initial dataframes
    # e.g., for DNA, its pos|ref|alt
    # snp_df['snp_str'] = np.nan
    snp_df = snp_df.assign(
        snp_str=(
            snp_df[groupby_cols[0]]
            .astype(str)
            .str.cat([snp_df[col].astype(str) for col in groupby_cols[1:]], sep="|")
        )
    )

    # Collapse by Accession ID and count occurrences
    # Filter out low global occurrence SNPs
    valid_snps = snp_df.groupby("snp_str")["Accession ID"].count()
    valid_snps = valid_snps[valid_snps > count_threshold]

    # Filter SNPs by valid SNPs
    snp_df = snp_df.join(valid_snps.rename("count"), how="right", on="snp_str")

    # Map SNPs to integer IDs
    snp_map = pd.Series(snp_df["snp_str"].unique())
    # Flip index and values
    snp_map = pd.Series(snp_map.index.values, index=snp_map)
    snp_df["snp_id"] = snp_df["snp_str"].map(snp_map)

    snp_group_df = snp_df.groupby("Accession ID")["snp_id"].agg(list).reset_index()
    # Add back the sequences with no mutations
    snp_group_df = pd.concat(
        [
            snp_group_df,
            pd.DataFrame({"Accession ID": no_snp_seqs}).assign(
                snp_id=[[]] * len(no_snp_seqs)
            ),
        ],
        axis=0,
        ignore_index=True,
    ).set_index("Accession ID")

    return snp_group_df, snp_map
