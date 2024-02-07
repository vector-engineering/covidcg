# coding: utf-8

"""Clean metadata from Genbank

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import pandas as pd


def main():
    """Clean metadata from GenBank

    Required columns:
        "Accession ID": index

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

    parser.add_argument("--metadata-in", type=str, required=True, help="Metadata in")
    parser.add_argument("--metadata-out", type=str, required=True, help="Metadata out")

    args = parser.parse_args()

    df = pd.read_csv(
        args.metadata_in,
        usecols=[
            "genbank_accession",
            "database",
            "strain",
            "region",
            "location",
            "collected",
            "submitted",
            "host",
            "isolation_source",
            "biosample_accession",
            "authors",
            "publications",
            "protein_names",
            # 'sequence', 'length', 'title'
        ],
    )

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
    #   publications,
    #   serotype,
    #   protein_names

    # Rename columns, set index
    df = df.rename(
        columns={
            "genbank_accession": "Accession ID",
            "strain": "virus_name",
            "submitted": "submission_date",
            "collected": "collection_date",
        }
    )

    # Drop duplicate accession IDs
    df.drop_duplicates("Accession ID", keep="first", inplace=True)

    df = df.set_index("Accession ID")

    # Remove sequences without region, collection date, or submission date
    # Remove sequences without F and G proteins
    remove_rows = (
        (df["region"].isna())
        | (df["submission_date"].isna())
        | (df["collection_date"].isna())
    )
    df = df.loc[~remove_rows]

    # Remove "Z" from the end of the submission date string, and convert from
    # ISO datetime to ISO date
    df.loc[:, "submission_date"] = pd.to_datetime(
        df["submission_date"], errors="coerce"
    ).apply(lambda x: x.strftime("%Y-%m-%d"))
    # Remove rows with invalid submission dates
    df = df.loc[~(df["submission_date"].isna())]

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

    # Isolate ID = same as Accession ID
    df["isolate_id"] = df.index.values
    # Segment = 1
    df["segment"] = 1

    df.to_csv(args.metadata_out)


if __name__ == "__main__":
    main()
