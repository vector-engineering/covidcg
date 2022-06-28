# coding: utf-8

"""Load mutation data, create mutation signatures

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import io
import pandas as pd

from pathlib import Path


def process_mutations(
    manifest,
    mutation_files,
    # Mutations must occur at least this many times to pass filters
    count_threshold=3,
    mode="dna",  # dna, gene_aa, protein_aa
):
    """Process mutation data
    
    Parameters
    ----------
    manifest: pandas.DataFrame
        - Sequence manifest (all sequence-reference pairs)
    mutation_files: list of strings
        - Paths to mutation CSV files
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

    # Dump all mutation chunks into a text buffer
    mutation_df_io = io.StringIO()
    for i, chunk in enumerate(mutation_files):
        file_name = Path(chunk).name.replace("_" + mode + "_mutation.csv", "")
        with open(chunk, "r") as fp_in:
            # Write dates, so we can remove duplicate sequences
            # and default to the mutations of the latest sequence, by date
            for j, line in enumerate(fp_in):
                # Write the header of the first file
                if i == 0 and j == 0:
                    mutation_df_io.write(line.strip() + ",file_name\n")
                # Or write any line that's not the header
                # (to avoid writing the header more than once)
                elif j > 0:
                    mutation_df_io.write(line.strip() + "," + file_name + "\n")

    # Read the buffer into a dataframe, then discard the buffer
    mutation_df_io.seek(0)
    mutation_df = pd.read_csv(mutation_df_io)
    mutation_df_io.close()

    # --------------------------
    # Remove duplicate sequences
    # --------------------------

    # Also has the effect of adding rows for sequences without mutations
    # (pos, ref, alt, etc filled with NaNs)
    mutation_df = manifest.merge(
        mutation_df,
        how="left",
        left_on=["Accession ID", "reference", "file_name"],
        right_on=["Accession ID", "reference", "file_name"],
    )

    # Remove sequences with no mutations before counting mutations
    mutation_df = mutation_df.loc[~mutation_df["pos"].isna()]
    mutation_df.loc[:, "pos"] = mutation_df["pos"].astype(int)
    # Replace NaNs in the 'ref' and 'alt' column with '-'
    mutation_df.fillna("-", inplace=True)

    groupby_cols = ["pos", "ref", "alt"]
    if mode == "dna":
        groupby_cols.insert(0, "segment")
    elif mode == "gene_aa":
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
        manifest.drop(columns=["file_name", "date", "subtype", "segment"])
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
