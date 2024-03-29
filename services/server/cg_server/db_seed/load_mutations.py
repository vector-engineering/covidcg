# coding: utf-8

"""Add additional fields to mutations
Used only by database seeder

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd

# from https://personal.sron.nl/~pault/#sec:qualitative
# 'vibrant', then 'muted'
mutation_colors = [
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


def get_mutation_colors(mutations):
    colors = []
    for i in range(len(mutations)):
        colors.append(mutation_colors[i % len(mutation_colors)])
    return colors


def format_nt_mutation(dna_mutation_str):
    # Print as REF POS ALT
    # i.e., 23403|A|G -> A23403G
    chunks = dna_mutation_str.split("|")
    return "{}{}{}".format(chunks[2], chunks[1], chunks[3])


def process_dna_mutations(dna_mutation_dict):
    dna_mutation = (
        pd.DataFrame.from_dict(dna_mutation_dict, orient="index", columns=["id"])
        .reset_index()
        .rename(columns={"index": "mutation_str"})
    )
    dna_mutation = dna_mutation.join(
        dna_mutation["mutation_str"]
        .str.split("|", expand=True)
        .rename(columns={0: "segment", 1: "pos", 2: "ref", 3: "alt"})
    ).set_index("id")
    dna_mutation.loc[:, "pos"] = dna_mutation["pos"].astype(int)
    dna_mutation["color"] = get_mutation_colors(dna_mutation.index.values)
    dna_mutation["mutation_name"] = dna_mutation["mutation_str"].apply(
        format_nt_mutation
    )

    return dna_mutation


def format_aa_mutation(aa_mutation_str):
    # Print as GENE/PROTEIN · REF POS ALT
    # i.e., S|614|D|G -> S · D614G
    chunks = aa_mutation_str.split("|")
    return "{}:{}{}{}".format(chunks[0], chunks[2], chunks[1], chunks[3])


def process_aa_mutations(aa_mutation_dict):
    """
    Parameters
    ----------
    aa_mutation_dict: dict
    """
    aa_mutation = (
        pd.DataFrame.from_dict(aa_mutation_dict, orient="index", columns=["id"])
        .reset_index()
        .rename(columns={"index": "mutation_str"})
    )
    aa_mutation = aa_mutation.join(
        aa_mutation["mutation_str"]
        .str.split("|", expand=True)
        .rename(columns={0: "feature", 1: "pos", 2: "ref", 3: "alt"})
    ).set_index("id")
    aa_mutation.loc[:, "pos"] = aa_mutation["pos"].astype(int)
    aa_mutation["color"] = get_mutation_colors(aa_mutation.index.values)
    aa_mutation["mutation_name"] = aa_mutation["mutation_str"].apply(format_aa_mutation)

    return aa_mutation
