# coding: utf-8

"""Load mutation data, create mutation signatures

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import gzip
import io
import pandas as pd
import numpy as np

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


def process_mutations(
    processed_fasta_files,
    mutation_files,
    references,
    # Mutations must occur at least this many times to pass filters
    count_threshold=3,
    mode="dna",  # dna, gene_aa, protein_aa
):
    """Process mutation data
    
    Parameters
    ----------
    processed_fasta_files: string
        - Path to processed FASTA files folder
    mutation_files: list of strings
        - Paths to mutation CSV files
    references: list of strings
        - Names of reference FASTA files
    count_threshold: int
        - Mutations must occur at least this many times to pass filters
    mode: string
        - dna, gene_aa, protein_aa
    
    Returns
    -------
    out: tuple of pandas.DataFrames
        - Mutation group dataframe
        - Mutation map dataframe
    """

    manifest = []
    for fasta_file in sorted(Path(processed_fasta_files).glob("*.fa.gz")):
        manifest.extend(extract_ids(fasta_file))
    manifest = pd.DataFrame.from_records(manifest, columns=["Accession ID", "date"])
    pruned_manifest = manifest.drop_duplicates(["Accession ID"], keep="last")

    # Add references to manifest
    pruned_manifest[
        "reference"
    ] = np.nan  # Have to first assign a singular value before setting to array
    # pruned_manifest['reference'] = pruned_manifest['reference'].apply(lambda _: ['NC_038235.1', 'NC_001781.1', 'KX858757.1', 'KX858756.1'])
    pruned_manifest["reference"] = pruned_manifest["reference"].apply(
        lambda _: references
    )
    pruned_manifest = pruned_manifest.explode("reference")

    # Dump all mutation chunks into a text buffer
    mutation_df_io = io.StringIO()
    for i, chunk in enumerate(mutation_files):
        file_date = Path(chunk).name.replace("_" + mode + "_mutation.csv", "")
        with open(chunk, "r") as fp_in:
            # Write dates, so we can remove duplicate sequences
            # and default to the mutations of the latest sequence, by date
            for j, line in enumerate(fp_in):
                # Write the header of the first file
                if i == 0 and j == 0:
                    mutation_df_io.write(line.strip() + ",date\n")
                # Or write any line that's not the header
                # (to avoid writing the header more than once)
                elif j > 0:
                    mutation_df_io.write(line.strip() + "," + file_date + "\n")

    # Read the buffer into a dataframe, then discard the buffer
    mutation_df_io.seek(0)
    mutation_df = pd.read_csv(mutation_df_io)
    mutation_df_io.close()

    # --------------------------
    # Remove duplicate sequences
    # --------------------------

    # Also has the effect of adding rows for sequences without mutations
    # (pos, ref, alt, etc filled with NaNs)
    mutation_df = pruned_manifest.set_index(["Accession ID", "reference", "date"]).join(
        mutation_df.set_index(["Accession ID", "reference", "date"]), how="left"
    )
    mutation_df.reset_index(inplace=True)

    # Remove sequences with no mutations before counting mutations
    mutation_df = mutation_df.loc[~mutation_df["pos"].isna()]
    mutation_df.loc[:, "pos"] = mutation_df["pos"].astype(int)
    # Replace NaNs in the 'ref' and 'alt' column with '-'
    mutation_df.fillna("-", inplace=True)

    groupby_cols = ["pos", "ref", "alt"]
    if mode == "gene_aa":
        groupby_cols.insert(0, "gene")
    elif mode == "protein_aa":
        groupby_cols.insert(0, "protein")

    # Generate mutation strings for initial dataframes
    # e.g., for DNA, its pos|ref|alt
    # mutation_df['mutation_str'] = np.nan
    mutation_df = mutation_df.assign(
        mutation_str=(
            mutation_df[groupby_cols[0]]
            .astype(str)
            .str.cat(
                [mutation_df[col].astype(str) for col in groupby_cols[1:]], sep="|"
            )
        )
    )

    # Collapse by Accession ID and count occurrences
    # Filter out low global occurrence mutations
    valid_mutations = mutation_df.groupby("mutation_str")["Accession ID"].count()
    valid_mutations = valid_mutations[valid_mutations > count_threshold]

    # Filter mutations by valid mutations
    mutation_df = mutation_df.join(
        valid_mutations.rename("count"), how="right", on="mutation_str"
    )

    # Map mutations to integer IDs
    mutation_map = pd.Series(mutation_df["mutation_str"].unique())
    # Flip index and values
    mutation_map = pd.Series(mutation_map.index.values, index=mutation_map)
    mutation_df["mutation_id"] = mutation_df["mutation_str"].map(mutation_map)

    mutation_group_df = mutation_df.groupby(
        ["Accession ID", "reference"], as_index=False
    )["mutation_id"].agg(list)

    mutation_group_df = (
        pruned_manifest.drop(columns=["date"])
        .merge(
            mutation_group_df,
            left_on=["Accession ID", "reference"],
            right_on=["Accession ID", "reference"],
            how="left",
        )
        .sort_values(["Accession ID", "reference"])
    )

    # Fill NaNs with empty arrays
    mutation_group_df.loc[
        mutation_group_df["mutation_id"].isna(), "mutation_id"
    ] = pd.Series(
        [[]] * mutation_group_df["mutation_id"].isna().sum(),
        index=mutation_group_df.index[mutation_group_df["mutation_id"].isna()],
    )

    return mutation_group_df, mutation_map
