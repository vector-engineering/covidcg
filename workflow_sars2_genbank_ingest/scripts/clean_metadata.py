# coding: utf-8

"""Clean metadata from Genbank

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import datetime
import pandas as pd


def main():
    """Clean metadata from GenBank

    Required columns:
        "genbank_accession": index

        # date information
        "collection_date"
        "submission_date"

        # location information
        "region"
        "country"
        "division"
        "location"
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--metadata-in",
        type=str,
        required=True,
        help="Path to input metadata CSV file",
    )

    parser.add_argument(
        "--lineages",
        type=str,
        required=True,
        help="Path to lineages CSV file",
    )
    parser.add_argument(
        "--quality",
        type=str,
        required=True,
        help="Path to quality CSV file",
    )

    parser.add_argument(
        "-o",
        "--metadata-out",
        type=str,
        required=True,
        help="Path to output metadata CSV file",
    )

    args = parser.parse_args()

    df = pd.read_csv(args.metadata_in)

    # Fields:
    #   genbank_accession,
    #   database,
    #   strain,
    #   region,
    #   location,
    #   collected,
    #   submitted,
    #   length,
    #   host,
    #   isolation_source,
    #   biosample_accession,
    #   title,
    #   authors,
    #   publications

    # Remove some unnecessary columns
    df = df.drop(columns=["title", "length"])

    # Rename columns, set index
    df = df.rename(
        columns={
            "genbank_accession": "Accession ID",
            "strain": "virus_name",
            "submitted": "submission_date",
            "collected": "collection_date",
        }
    )
    df = df.set_index("Accession ID")

    # Remove sequences without region, collection date, or submission date
    remove_rows = (
        (df["region"].isna())
        | (df["submission_date"].isna())
        | (df["collection_date"].isna())
    )
    df = df.loc[~remove_rows]

    # Remove "Z" from the end of the submission date string, and convert from
    # ISO datetime to ISO date
    def datetime_to_date(x):
        return datetime.datetime.fromisoformat(x[:-1]).strftime("%Y-%m-%d")

    df.loc[:, "submission_date"] = df["submission_date"].apply(datetime_to_date)

    # Parse location data
    def parse_genbank_location(s):
        """Convert a Genbank location string into a tuple of (country, division, location)
        Fill missing data with "-1" for groupbys later
        """

        # No additional location data
        if not s or type(s) is not str:
            return (-1, -1, -1)

        # The country is always the first part
        country_chunks = s.split(":")
        country = country_chunks[0].strip()

        # Only the country is defined
        if len(country_chunks) == 1:
            return (country, -1, -1)

        # If the second part exists, the division and location might
        # be comma-delimited, i.e., "division, location"

        division_chunks = country_chunks[1].split(",")
        division = division_chunks[0].strip()

        # Only the division is defined
        if len(division_chunks) == 1:
            return (country, division, -1)

        # All are defined
        return (country, division, division_chunks[1].strip())

    loc_tuples = df["location"].apply(parse_genbank_location)
    # Drop the original "location" column, and join on new columns
    df = df.drop(columns=["location"]).join(
        pd.DataFrame(
            [[a, b, c] for a, b, c in loc_tuples.values],
            columns=["country", "division", "location"],
            index=loc_tuples.index,
        )
    )

    # Fill in missing values with empty strings, for remaining
    # metadata columns
    fill_in_cols = [
        "database",
        "virus_name",
        "host",
        "isolation_source",
        "biosample_accession",
        "authors",
        "publications",
    ]
    for col in fill_in_cols:
        df.loc[:, col] = df[col].fillna("Unknown")

    # Load lineages and join to dataframe
    lineages_df = pd.read_csv(args.lineages)
    lineages_df = lineages_df.rename(columns={"taxon": "Accession ID"}).set_index(
        "Accession ID"
    )
    df = df.join(lineages_df[["lineage"]], how="left")
    # Fill in missing values
    df.loc[:, "lineage"] = df["lineage"].fillna("None")

    # Isolate ID = same as Accession ID
    df["isolate_id"] = df.index.values
    # Segment = 1
    df["segment"] = 1

    # Load quality and join to dataframe
    quality_df = pd.read_csv(args.quality, index_col="Accession ID")
    df = df.join(quality_df, how="left")
    df["length"] = df["length"].fillna(0).astype(int)
    # Calculate percent ambiguous, drop the num_ambiguous column
    df["num_ambiguous"] = ((df["num_ambiguous"] / df["length"]) * 100).fillna(0)
    df.rename(columns={"num_ambiguous": "percent_ambiguous"}, inplace=True)

    # Filter out entries without any sequence
    df = df.loc[df["length"] > 0]

    df.to_csv(args.metadata_out)


if __name__ == "__main__":
    main()
