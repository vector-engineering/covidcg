#!/usr/bin/env python3
# coding: utf-8

"""Clean metadata

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import numpy as np
import pandas as pd

"""All columns:

Isolate_Id
PB2 Segment_Id
PB1 Segment_Id
PA Segment_Id
HA Segment_Id
NP Segment_Id
NA Segment_Id
MP Segment_Id
NS Segment_Id
HE Segment_Id
P3 Segment_Id
Isolate_Name
Subtype
Lineage -- sometimes this column doesn't exist
Clade -- sometimes this column doesn't exist
Passage_History
Location
Host
Isolate_Submitter
Submitting_Lab
Submitting_Sample_Id
Authors
Publication
Originating_Lab
Originating_Sample_Id
Collection_Date
Note
Update_Date
Submission_Date
Antigen_Character
Animal_Vaccin_Product
Adamantanes_Resistance_geno
Oseltamivir_Resistance_geno
Zanamivir_Resistance_geno
Peramivir_Resistance_geno
Other_Resistance_geno
Adamantanes_Resistance_pheno
Oseltamivir_Resistance_pheno
Zanamivir_Resistance_pheno
Peramivir_Resistance_pheno
Other_Resistance_pheno
Host_Age
Host_Age_Unit
Host_Gender
Patient_Status
Zip_Code
Outbreak
Pathogen_Test_Info
Is_Vaccinated
Human_Specimen_Source
Animal_Specimen_Source
Animal_Health_Status
Domestic_Status
PMID
PB2 INSDC_Upload
PB1 INSDC_Upload
PA INSDC_Upload
HA INSDC_Upload
NP INSDC_Upload
NA INSDC_Upload
MP INSDC_Upload
NS INSDC_Upload
HE INSDC_Upload
P3 INSDC_Upload
"""


def clean_df(df):
    df.drop(
        columns=[
            # Unsupported segments
            "HE Segment_Id",
            "P3 Segment_Id",
            # Either irrelevant or way too little data (or just empty)
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
            # "Adamantanes_Resistance_geno",
            # "Oseltamivir_Resistance_geno",
            # "Zanamivir_Resistance_geno",
            # "Peramivir_Resistance_geno",
            # "Other_Resistance_geno",
            # "Adamantanes_Resistance_pheno",
            # "Oseltamivir_Resistance_pheno",
            # "Zanamivir_Resistance_pheno",
            # "Peramivir_Resistance_pheno",
            # "Other_Resistance_pheno",
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

    # Add columns if they don't exist yet
    if "Clade" not in df.columns:
        df.insert(3, "Clade", np.nan)

    if "Lineage" not in df.columns:
        df.insert(3, "Lineage", np.nan)

    # Rename columns
    df.rename(
        columns={
            "Isolate_Id": "isolate_id",
            "Isolate_Name": "virus_name",
            "Subtype": "serotype",
            "Clade": "clade",
            "Lineage": "lineage",
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
            "Adamantanes_Resistance_geno": "adamantanes_resistance_geno",
            "Oseltamivir_Resistance_geno": "oseltamivir_resistance_geno",
            "Zanamivir_Resistance_geno": "zanamivir_resistance_geno",
            "Peramivir_Resistance_geno": "peramivir_resistance_geno",
            "Other_Resistance_geno": "other_resistance_geno",
            "Adamantanes_Resistance_pheno": "adamantanes_resistance_pheno",
            "Oseltamivir_Resistance_pheno": "oseltamivir_resistance_pheno",
            "Zanamivir_Resistance_pheno": "zanamivir_resistance_pheno",
            "Peramivir_Resistance_pheno": "peramivir_resistance_pheno",
            "Other_Resistance_pheno": "other_resistance_pheno",
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
            .apply(lambda x: [str(_x) for _x in list(x) if _x > 0], axis=1)
        ),
    )
    df.drop(columns=segment_cols, inplace=True)

    # Clean up subtype
    # Replace B serotypes with the lineage, if exists
    b_serotype = df["serotype"] == "B"
    df.loc[b_serotype, "serotype"] = df.loc[b_serotype, "lineage"].map(
        {"Victoria": "B-vic", "Yamagata": "B-yam"}
    )
    # Clean A serotypes
    df.loc[~b_serotype, "serotype"] = (
        df.loc[~b_serotype, "serotype"].str.split("/").apply(lambda x: x[1]).str.strip()
    )
    # Generalize H5, H7, H9, H10 serotypes
    # But first save the original serotype
    df["original_serotype"] = df["serotype"]
    df.loc[df["serotype"].str.startswith("H5"), "serotype"] = "H5NX"
    df.loc[df["serotype"].str.startswith("H7"), "serotype"] = "H7NX"
    # df['serotype'].str.replace(r'^H7N?[1-9]?$', 'H7NX', regex=True)
    df.loc[df["serotype"].str.startswith("H9"), "serotype"] = "H9NX"
    df.loc[df["serotype"].str.startswith("H10"), "serotype"] = "H10NX"

    # Extract N subtype
    df["n_subtype"] = "NA"
    df.loc[~b_serotype, "n_subtype"] = (
        df.loc[~b_serotype, "original_serotype"]
        .str.extract(r".*N(\d+)$", expand=False)
        .fillna("Unknown")
    )
    df.loc[b_serotype, "n_subtype"] = "NA"

    # Remove rows without segments
    df = df.loc[df["segments"].apply(len) > 0, :]

    # Expand location metadata
    # print("Processing location data...", end="", flush=True)
    # Location data is stored in one column, "region / country / division / location"
    location_df = (
        # Unset index for faster filtering
        df["Location"]
        .fillna("")
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
    # df.drop(columns=["Location"], inplace=True)
    # Rename original column
    df.rename(columns={"Location": "location_original"}, inplace=True)

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
    df.loc[:, "collection_date"] = pd.to_datetime(
        df["collection_date"], yearfirst=True, errors="coerce"
    )
    df.loc[:, "submission_date"] = pd.to_datetime(
        df["submission_date"], yearfirst=True, errors="coerce"
    )
    # Backfill submission date with collection date
    df.loc[:, "submission_date"] = df["submission_date"].combine_first(
        df["collection_date"]
    )

    # Remove rows without collection dates
    df = df.loc[~pd.isna(df["collection_date"]), :]

    # Infer "Original", "Cell", or "Egg" passage from the "passage" field

    passage_clean = df["passage"].fillna("").str.lower().str.strip()

    # Get rid of nuisance terms
    replace_map = {
        r"passage\s?(details)?:\s?": "",
        r".*mdck.*": "cell",
        r".*293t.*": "cell",
        r".*siat[123]?.*": "cell",
        r".*qmc[12].*": "cell",
        r".*spf[123].*": "cell",
        r".*rhmk[123]?.*": "cell",
        # r"p[0-9]{1}": "cell",
        # r"s[0-9]{1}": "cell",
        # r"C[0-9]{1}\+C[1-9]{1}": "cell",
        r"[0-9]+ passages?.+cells?.+": "cell",
        r"cell passage [0-9]+": "cell",
        r".*ax-4.*": "cell",
        # r"hCK": "cell",
        r".*caco-2.*": "cell",
        r".*md[123].*": "cell",
        r".*pmk1.*": "cell",
        r".*spfe1.*": "cell",
        r"^[csmr][0-9x]+[\+\/].*": "cell",
        r"^[csmr][0-9x]+[csm]{1}[0-9x]{1}.*": "cell",
        r"^[xm][0-9x]?\/[cs].*": "cell",
        r".*egg.*": "egg",
        r"^e[0-9x]{1}.*": "egg",
        r"^p[0-9]{1}\,e[0-9]+.*": "egg",
        r"^ece.*": "egg",
        r".*original.*": "original",
        r".*clinical.*": "original",
        r".*direct.*": "original",
        r".*autopsy.*": "original",
        r".*swab.*": "original",
        r".*organ.*": "original",
        r".*tissue.*": "original",
    }
    for k, v in replace_map.items():
        passage_clean = passage_clean.str.replace(k, v, regex=True)

    passage_map = {
        "original": [
            "original",
            "origina",
            "orginal",
            "org",
            "orginal_sample",
            "p0",
            "p0-ori",
            "blank",
            "cs",
            "cs-ori",
            "cs_ori",
            "or_ir",
            "or-ir",
            "ori",
            "sample",
            "pooled lungs and oropharyngeal-tracheal swab",
            "human",
            "human, unpassaged",
            "or",
            "initial",
            "primary specimen",
            "no passage",
            "swab",
            "nasal swab",
            "first",
            "orignal",
            "rna",
            "tissues",
            "isolated directly from host; no passage",
            "ferret",
            "op&np",
        ],
        "cell": [
            "c",
            "c0",
            "cell",
            "с1",
            "c1",
            "c1-ori",
            "c2",
            "c3",
            "c4",
            "c5",
            "cx",
            "c1 in allantoidal liquid",
            "s1",
            "s2",
            "s3",
            "s4",
            "sx",
            "c4",
            "c2hck2/c1",
            "p1",
            "p-1",
            "p2",
            "p3",
            "pi",
            # "x",
            "x1",
            "x2",
            "x3",
            "r0",
            "cxs1",
            "i-cs",
            "m1",
            "m2",
            "m3",
            "m4",
            "c1s1",
            "paepc1",
            "paepc2",
            "rmk 1st passage",
            "gmkc1",
        ],
        "egg": [
            "egg",
            "e0",
            "e1",
            "e2",
            "e3",
            "e4",
            "p1,e1",
            "p2,e1",
            "p2,e2",
            "p2,e3",
            "c1e1",
            "ce1",
        ],
    }
    passage_map_reverse = {}
    for k, v in passage_map.items():
        for vv in v:
            passage_map_reverse[vv] = k

    df["passage_category"] = passage_clean.map(passage_map_reverse)

    num_unmapped = (df["passage_category"].isna() & ~df["passage"].isna()).sum()
    print(f"Unmapped passage values: {num_unmapped} / {len(df)}")
    print(
        df["passage"][
            (df["passage_category"].isna() & ~df["passage"].isna())
        ].value_counts()
    )
    # print(
    #     passage_clean[
    #         (df["passage_category"].isna() & ~df["passage"].isna())
    #     ].value_counts()
    # )

    df["passage_category"].fillna("unknown", inplace=True)

    # Enforce column order for easier concatenation later
    df = df[
        [
            "accession_ids",
            "segments",
            "isolate_id",
            "virus_name",
            "serotype",
            "n_subtype",
            "lineage",
            "clade",
            "passage",
            "passage_category",
            "host",
            "isolate_submitter",
            "submitting_lab",
            "authors",
            "publication",
            "originating_lab",
            "collection_date",
            "submission_date",
            "adamantanes_resistance_geno",
            "oseltamivir_resistance_geno",
            "zanamivir_resistance_geno",
            "peramivir_resistance_geno",
            "other_resistance_geno",
            "adamantanes_resistance_pheno",
            "oseltamivir_resistance_pheno",
            "zanamivir_resistance_pheno",
            "peramivir_resistance_pheno",
            "other_resistance_pheno",
            "gender",
            "patient_status",
            "region",
            "country",
            "division",
            "location",
            "location_original",
        ]
    ]

    # Collapse by isolate

    # For these columns, create a list of values for each isolate
    agg_dict = {
        "accession_ids": ("accession_ids", sum),
        "segments": ("segments", sum),
    }
    # For all other columns, take the first value
    for col in df.columns[3:]:
        agg_dict[col] = (col, "first")

    df = df.groupby("isolate_id").agg(**agg_dict)

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
        help="Path to output metadata virus JSON file",
    )

    args = parser.parse_args()

    dfs = []
    for f in args.metadata_in:
        print("Processing {}...".format(f))
        df = pd.read_excel(f, sheet_name=0)
        dfs.append(clean_df(df))
    dfs = pd.concat(dfs, axis=0)

    dfs.to_json(args.metadata_virus_out, orient="records")

    # Expand by Accession ID
    df = dfs.explode(["accession_ids", "segments"]).rename(
        columns={
            "accession_ids": "Accession ID",
            "segments": "segment",
        }
    )
    df.to_csv(args.metadata_out)


if __name__ == "__main__":
    main()
