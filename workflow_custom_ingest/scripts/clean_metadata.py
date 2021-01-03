# coding: utf-8

"""Clean metadata from Genbank

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import datetime
import pandas as pd


def clean_metadata(metadata_in, lineages_in, metadata_out):
    """Clean metadata from GenBank

    Required columns:
        "Accession ID": index
        # For later processing this must be unique
        # within the first 8 characters.

        # date information
        "collection_date"
        "submission_date"

        # location information
        "region"
        "country"
        "division"
        "location"

    Parameters
    ----------
    metadata_in: str
    lineages_in: str
    metadata_out: str

    Returns
    -------
    None
    """

    df = pd.read_csv(metadata_in, index_col='Accession ID')

    # Fields:
    #   Accession ID,
    #   strain,
    #   region,
    #   country,
    #   division,
    #   location,
    #   collection_date,
    #   submission_date,
    #   host,

    # Remove sequences without region, collection date, or submission date
    remove_rows = (
        (df["region"].isna())
        | (df["submission_date"].isna())
        | (df["collection_date"].isna())
    )
    df = df.loc[~remove_rows]

    # Fill in missing values with empty strings, for remaining
    # metadata columns
    fill_in_cols = [
        "strain",
        "host",
    ]
    for col in fill_in_cols:
        df.loc[:, col] = df[col].fillna("Unknown")

    # Load lineages and join to dataframe
    lineages_df = pd.read_csv(lineages_in)
    lineages_df = lineages_df.rename(columns={"taxon": "Accession ID"}).set_index(
        "Accession ID"
    )
    df = df.join(lineages_df[["lineage"]], how="left")
    # Fill in missing values
    df.loc[:, "lineage"] = df["lineage"].fillna("None")

    df.to_csv(metadata_out)
