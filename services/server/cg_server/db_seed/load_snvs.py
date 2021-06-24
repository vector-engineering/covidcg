# coding: utf-8

"""Add additional fields to SNVs
Used only by database seeder

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd

# from https://personal.sron.nl/~pault/#sec:qualitative
# 'vibrant', then 'muted'
snv_colors = [
    "#0077bb",
    "#33bbee",
    "#009988",
    "#ee7733",
    "#cc3311",
    "#ee3377",
    "#332288",
    "#88ccee",
    "#44aa99",
    "#117733",
    "#999933",
    "#ddcc77",
    "#cc6677",
    "#882255",
    "#aa4499",
]


def get_snv_colors(snvs):
    colors = []
    for i, snv in enumerate(snvs):
        colors.append(snv_colors[i % len(snv_colors)])
    return colors


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


def format_aa_snv(snv_str):
    # Print as GENE/PROTEIN · REF POS ALT
    # i.e., S|614|D|G -> S · D614G
    chunks = snv_str.split("|")
    return "{}:{}{}{}".format(chunks[0], chunks[2], chunks[1], chunks[3])


def process_aa_snvs(aa_snv_dict, name, defs):
    """
    Parameters
    ----------
    aa_snv_dict: dict
    name: str
    defs: pandas.DataFrame
        - genes or proteins dataframe from genes_and_proteins.py
    """
    aa_snv = (
        pd.DataFrame.from_dict(aa_snv_dict, orient="index", columns=["id"])
        .reset_index()
        .rename(columns={"index": "snp_str"})
    )
    aa_snv = aa_snv.join(
        aa_snv["snp_str"]
        .str.split("|", expand=True)
        .rename(columns={0: name, 1: "pos", 2: "ref", 3: "alt"})
    ).set_index("id")
    aa_snv.loc[:, "pos"] = aa_snv["pos"].astype(int)
    aa_snv["color"] = get_snv_colors(aa_snv.index.values)
    aa_snv["snv_name"] = aa_snv["snp_str"].apply(format_aa_snv)

    def get_nt_pos(snv):
        segment_ind = [
            i
            for i, segment in enumerate(defs.at[snv[name], "aa_ranges"])
            if snv["pos"] >= segment[0] and snv["pos"] <= segment[1]
        ][0]
        return defs.at[snv[name], "segments"][segment_ind][0] + (
            (snv["pos"] - defs.at[snv[name], "aa_ranges"][segment_ind][0]) * 3
        )

    aa_snv["nt_pos"] = aa_snv.apply(get_nt_pos, axis=1)

    return aa_snv
