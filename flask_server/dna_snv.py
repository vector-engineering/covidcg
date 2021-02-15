# coding: utf-8

import pandas as pd

from .color import get_snv_colors


def format_nt_snv(snv_str):
    # Print as REF POS ALT
    # i.e., 23403|A|G -> A23403G
    chunks = snv_str.split("|")
    return "{}{}{}".format(chunks[1], chunks[0], chunks[2])


def process_dna_snvs(dna_snv_dict):
    dna_snv = (
        pd.DataFrame.from_dict(dna_snv_dict, orient="index", columns=["id"])
        .reset_index()
        .rename(columns={"index": "snp_str"})
    )
    dna_snv = dna_snv.join(
        dna_snv["snp_str"]
        .str.split("|", expand=True)
        .rename(columns={0: "pos", 1: "ref", 2: "alt"})
    ).set_index("id")
    dna_snv.loc[:, "pos"] = dna_snv["pos"].astype(int)
    dna_snv["color"] = get_snv_colors(dna_snv.index.values)
    dna_snv["snv_name"] = dna_snv["snp_str"].apply(format_nt_snv)

    return dna_snv
