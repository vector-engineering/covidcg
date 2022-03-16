#!/usr/bin/env python3
# coding: utf-8

"""Clean metadata

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import numpy as np
import pandas as pd


def clean_df(df):
    df.drop(
        columns=[
            # Unsupported segments
            "HE Segment_Id",
            "P3 Segment_Id",
            # Either irrelevant or way too little data (or just empty)
            "Isolate_Id",
            # 'Lineage',
            "Note",
            "Antigen_Character",
            "Animal_Vaccin_Product",
            "Zip_Code",
            "Pathogen_Test_Info",
            "Animal_Health_Status",
            "Domestic_Status",
            "PMID",
            # Include these later if we want to link to GenBank?
            "PB2 INSDC_Upload",
            "PB1 INSDC_Upload",
            "PA INSDC_Upload",
            "HA INSDC_Upload",
            "NP INSDC_Upload",
            "NA INSDC_Upload",
            "MP INSDC_Upload",
            "NS INSDC_Upload",
            "HE INSDC_Upload",
            "P3 INSDC_Upload",
            # Include these later?
            "Submitting_Sample_Id",
            "Originating_Sample_Id",
            "Adamantanes_Resistance_geno",
            "Oseltamivir_Resistance_geno",
            "Zanamivir_Resistance_geno",
            "Peramivir_Resistance_geno",
            "Other_Resistance_geno",
            "Adamantanes_Resistance_pheno",
            "Oseltamivir_Resistance_pheno",
            "Zanamivir_Resistance_pheno",
            "Peramivir_Resistance_pheno",
            "Other_Resistance_pheno",
            "Host_Age",
            "Host_Age_Unit",
            "Outbreak",
            "Is_Vaccinated",
            "Human_Specimen_Source",
            "Animal_Specimen_Source",
            "Update_Date",
        ],
        inplace=True,
    )

    # Rename columns
    df.rename(
        columns={
            "Isolate_Name": "virus_name",
            "Subtype": "serotype",
            "Passage_History": "passage",
            "Host": "host",
            "Isolate_Submitter": "isolate_submitter",
            "Submitting_Lab": "submitting_lab",
            "Authors": "authors",
            "Publication": "publication",
            "Originating_Lab": "originating_lab",
            "Collection_Date": "collection_date",
            "Submission_Date": "submission_date",
            # 'Update_Date': 'update_date',
            "Host_Gender": "gender",
            "Patient_Status": "patient_status",
        },
        inplace=True,
    )

    segment_cols = [
        "PB2 Segment_Id",
        "PB1 Segment_Id",
        "PA Segment_Id",
        "HA Segment_Id",
        "NP Segment_Id",
        "NA Segment_Id",
        "MP Segment_Id",
        "NS Segment_Id",
    ]

    df.insert(
        0,
        "accession_ids",
        (
            df[segment_cols]
            .fillna("")
            .apply(lambda x: x.str.split("|"))
            .applymap(lambda x: [x[0]] if len(x) > 1 else [])
            .apply(lambda x: sum(x, []), axis=1)
        ),
    )

    df.insert(
        1,
        "segments",
        (
            # Mask missing segments
            (
                df[segment_cols]
                # .applymap(lambda x: 1 if not np.isnan(x) else np.nan)
                .isna()
            )
            .astype(int)
            .applymap(lambda x: 0 if x == 1 else 1)
            # Segments are already in order from 1-8, so assign number based on column position
            .multiply(np.array([1, 2, 3, 4, 5, 6, 7, 8]))
            # Join into one array per isolate
            .apply(lambda x: [_x for _x in list(x) if _x > 0], axis=1)
        ),
    )
    df.drop(columns=segment_cols, inplace=True)

    # Clean up subtype
    # Replace B serotypes with the lineage, if exists
    b_serotype = df["serotype"] == "B"
    df.loc[b_serotype, "serotype"] = df.loc[b_serotype, "Lineage"].map(
        {"Victoria": "B-vic", "Yamagata": "B-yam"}
    )
    # Clean A serotypes
    df.loc[~b_serotype, "serotype"] = (
        df.loc[~b_serotype, "serotype"].str.split("/").apply(lambda x: x[1]).str.strip()
    )

    # Expand location metadata
    # print("Processing location data...", end="", flush=True)
    # Location data is stored in one column, "region / country / division / location"
    location_df = (
        # Unset index for faster filtering
        df["Location"]
        .fillna('')
        .str.split("/", expand=True)
        .iloc[:, :4]  # Only take 4 columns
        # Rename columns
        .rename(columns={0: "region", 1: "country", 2: "division", 3: "location"})
        .applymap(lambda x: x.strip() if x else x)
        # Placeholder for missing values, so that it will still
        # be caught by groupby() later on
        .fillna(-1)
    )

    # With sparse data, sometimes none of the sequences fed in will have a
    # "division" or "location" entry.
    # Make them manually now, if they don't already exist
    if "division" not in location_df.columns:
        location_df["division"] = -1
    if "location" not in location_df.columns:
        location_df["location"] = -1

    # Clean location data
    # location_df = clean_location_data(location_df, location_corretions)

    # Join values to master dataframe
    # Do an inner join... but none of the rows should've been filtered out
    # so this doesn't really matter
    df = df.join(location_df, how="inner")
    # Drop original column
    df.drop(columns=["Location"], inplace=True)

    # Fill in "Unknown" in sparse cols
    df.loc[:, "isolate_submitter"] = (
        df["isolate_submitter"].fillna("Unknown").astype(str).str.strip()
    )
    df.loc[:, "submitting_lab"] = (
        df["submitting_lab"].fillna("Unknown").astype(str).str.strip()
    )
    df.loc[:, "authors"] = df["authors"].fillna("Unknown").astype(str).str.strip()
    df.loc[:, "publication"] = (
        df["publication"].fillna("Unknown").astype(str).str.strip()
    )
    df.loc[:, "originating_lab"] = (
        df["originating_lab"].fillna("Unknown").astype(str).str.strip()
    )
    df.loc[:, "gender"] = df["gender"].fillna("Unknown").astype(str).str.strip()
    df.loc[:, "patient_status"] = (
        df["patient_status"].fillna("Unknown").astype(str).str.strip()
    )

    # Convert dates to datetime
    df.loc[:, "collection_date"] = pd.to_datetime(df["collection_date"], yearfirst=True)
    df.loc[:, "submission_date"] = pd.to_datetime(df["submission_date"], yearfirst=True)
    # Backfill submission date with collection date
    df.loc[:, 'submission_date'] = df['submission_date'].combine_first(df['collection_date'])

    # Enforce column order for easier concatenation later
    df = df[
        [
            "accession_ids",
            "segments",
            "virus_name",
            "serotype",
            "passage",
            "host",
            "isolate_submitter",
            "submitting_lab",
            "authors",
            "publication",
            "originating_lab",
            "collection_date",
            "submission_date",
            "gender",
            "patient_status",
            "region",
            "country",
            "division",
            "location",
        ]
    ]

    # Collapse by isolate
    df = df.groupby("virus_name").agg(
        accession_ids=("accession_ids", sum),
        segments=("segments", sum),
        serotype=("serotype", "first"),
        passage=("passage", "first"),
        host=("host", "first"),
        isolate_submitter=("isolate_submitter", "first"),
        submitting_lab=("submitting_lab", "first"),
        authors=("authors", "first"),
        publication=("publication", "first"),
        originating_lab=("originating_lab", "first"),
        collection_date=("collection_date", "first"),
        submission_date=("submission_date", "first"),
        gender=("gender", "first"),
        patient_status=("patient_status", "first"),
        region=("region", "first"),
        country=("country", "first"),
        division=("division", "first"),
        location=("location", "first"),
    )

    return df


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--metadata-in",
        type=str,
        nargs="+",
        required=True,
        help="Path to input metadata CSV files",
    )

    parser.add_argument(
        "--metadata-out",
        type=str,
        required=True,
        help="Path to output metadata CSV file",
    )
    parser.add_argument(
        "--metadata-virus-out",
        type=str,
        required=True,
        help="Path to output metadata virus CSV file",
    )

    args = parser.parse_args()

    dfs = []
    for f in args.metadata_in:
        print("Processing {}...".format(f))
        df = pd.read_excel(f, sheet_name=0)
        dfs.append(clean_df(df))
    dfs = pd.concat(dfs, axis=0)

    dfs.to_csv(args.metadata_virus_out)

    # Expand by Accession ID
    df = dfs.explode(["accession_ids", "segments"]).rename(
        columns={"accession_ids": "Accession ID", "segments": "segment",}
    )
    df.to_csv(args.metadata_out)


if __name__ == "__main__":
    main()
