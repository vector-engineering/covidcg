# coding: utf-8

"""Get representative sequences from each lineage

Code is derived from CoVizu:
  * https://filogeneti.ca/covizu/
  * https://github.com/PoonLab/covizu

Please find the attached license at LICENSE_COVIZU

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd
import numpy as np
import json

from datetime import date


def date2float(isodate):
    """ Convert ISO date string to float (years) """
    year, month, day = map(int, isodate.split("-"))
    dt = date(year, month, day)
    origin = date(dt.year, 1, 1)
    td = (dt - origin).days
    return dt.year + td / 365.25


def get_representative_seqs(
    isolate_data_path,
    metadata_map_path,
    ref_seq_path,
    fasta_out_path,
    datefile_out_path,
    table_out_path,
):

    # Load DNA mutation ID map
    with open(metadata_map_path, "r") as fp:
        metadata_map = json.loads(fp.read())
    dna_mutation = metadata_map["dna_mutation"]
    id_to_dna_mutation = {v: k for k, v in dna_mutation.items()}

    # Load data and select a representative sequence for each lineage
    # by choosing the last seen sequence from each lineage
    df = pd.read_json(isolate_data_path)
    df["collection_date"] = pd.to_datetime(df["collection_date"])

    # Join region string
    # Map integer IDs to strings
    df.loc[:, "region"] = df["region"].map(
        {int(k): v for k, v in metadata_map["region"].items()}
    )

    # https://github.com/cov-lineages/pango-designation
    # See: DOI 10.1038/s41564-020-0770-5, 10.1093/ve/veab064
    pango_lineages = pd.read_csv(
        "https://github.com/cov-lineages/pango-designation/raw/master/lineages.csv"
    )

    reps = (
        df.loc[
            (df["collection_date"] >= pd.to_datetime("2019-12-15"))
            & (df["collection_date"] <= pd.to_datetime(date.today()))
            & (
                df["virus_name"]
                .str.replace("hCoV-19/", "")
                .isin(pango_lineages["taxon"])
            )
        ]
        .drop_duplicates("lineage", keep="first")[
            ["isolate_id", "collection_date", "lineage", "dna_mutation"]
        ]
        .set_index("isolate_id")
        # Join first and last collection dates
        .join(
            (
                df.loc[
                    (df["collection_date"] > pd.to_datetime("2019-12-15"))
                    & (df["collection_date"] < pd.to_datetime(date.today()))
                ]
                .groupby("lineage")
                .agg(
                    date_min=("collection_date", np.min),
                    date_max=("collection_date", np.max),
                )
            ),
            on="lineage",
        )
        # Join region with the most counts of this lineage
        .join(
            (
                df.groupby(["lineage", "region"], as_index=False)[["isolate_id"]]
                .agg("count")
                .groupby("lineage")
                .apply(lambda x: x["region"].values[np.argmax(x["isolate_id"])])
                .rename("region_most_common")
            ),
            on="lineage",
        )
    )

    # Load reference sequence
    with open(ref_seq_path, "r") as fp:
        ref = json.loads(fp.read())
    ref_seq = ref["WIV04"]["segments"]["1"]["sequence"]

    # For each representative sequence, use its DNA mutations to reconstruct its
    # "MSA" genome from the reference genome
    # Ignore insertions, since this would require a true MSA
    reps["sequence"] = ref_seq
    for accession_id, row in reps.iterrows():
        dna_mutations = row["dna_mutation"]
        seq = list(ref_seq)

        for mut in dna_mutations:
            mut = id_to_dna_mutation[mut].split("|")
            pos = int(mut[1])
            ref = mut[2]
            alt = mut[3]

            # Skip insertions
            if ref == "-" or len(alt) > len(ref):
                continue

            # Deletions
            if alt == "-":
                seq[pos - 1 : (pos - 1) + len(ref)] = ["-"] * len(ref)
            # Mutations
            else:
                seq[pos - 1 : (pos - 1) + len(ref)] = list(alt)

        reps.loc[accession_id, "sequence"] = "".join(seq)

    # Write "MSA" sequences as a FASTA file
    with open(fasta_out_path, "w") as fp:
        for accession_id, row in reps.iterrows():
            fp.write(">{}\n{}\n".format(accession_id, row["sequence"]))

    # Write files for treetime
    # Extract dates from sequence headers
    reps["collection_date"] = reps["collection_date"].astype(str)
    with open(datefile_out_path, "w") as fp:
        fp.write("name,date\n")
        for accession_id, row in reps.iterrows():
            # Date format should be fine: https://github.com/neherlab/treetime#metadata-and-date-format
            fp.write("{},{}\n".format(accession_id, date2float(row["collection_date"])))
        fp.close()

    # Write dataframe for future
    reps.drop(columns=["dna_mutation"]).to_csv(table_out_path)
