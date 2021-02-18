# coding: utf-8

import pandas as pd

from .color import get_snv_colors


def format_aa_snv(snv_str):
    # Print as GENE/PROTEIN · REF POS ALT
    # i.e., S|614|D|G -> S · D614G
    chunks = snv_str.split("|")
    return "{}:{}{}{}".format(chunks[0], chunks[2], chunks[1], chunks[3])


# TODO: maybe this should be moved into workflow_main?
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
            for i, segment in enumerate(defs.at[snv[name], "aa_segments"])
            if snv["pos"] >= segment[0] and snv["pos"] <= segment[1]
        ][0]
        return defs.at[snv[name], "segments"][segment_ind][0] + (
            (snv["pos"] - defs.at[snv[name], "aa_segments"][segment_ind][0]) * 3
        )

    aa_snv["nt_pos"] = aa_snv.apply(get_nt_pos, axis=1)

    return aa_snv
