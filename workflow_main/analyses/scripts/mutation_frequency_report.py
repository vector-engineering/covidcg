# coding: utf-8

"""Add mutation metadata to group mutation frequencies file
Also spit out more readable CSVs for easier ingestion

Author: Albert Chen (albert_chen@g.harvard.edu)
"""

import argparse
import json

import pandas as pd


def dna_mutation_to_name(mut_str):
    split = mut_str.split("|")
    # REF POS ALT
    return split[2] + split[1] + split[3]


def aa_mutation_to_name(mut_str):
    split = mut_str.split("|")
    # REF POS ALT
    return split[2] + split[1] + split[3]


def main():
    """Command-line entry point"""

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--group-mutation-frequencies",
        type=str,
        required=True,
        help="Path to group mutation frequencies JSON file",
    )
    parser.add_argument(
        "--metadata-map", type=str, required=True, help="Metadata map JSON file"
    )
    parser.add_argument(
        "--output-json", type=str, required=True, help="Path to output JSON file"
    )
    parser.add_argument(
        "--output-dna-csv", type=str, required=True, help="Path to output DNA CSV file"
    )
    parser.add_argument(
        "--output-gene-aa-csv",
        type=str,
        required=True,
        help="Path to output gene AA CSV file",
    )
    parser.add_argument(
        "--output-protein-aa-csv",
        type=str,
        required=True,
        help="Path to output protein AA CSV file",
    )
    args = parser.parse_args()

    # Open mutation frequency data
    with open(args.group_mutation_frequencies, "r") as fp:
        group_mutation_frequencies = json.loads(fp.read())

    # Load mutation ID maps
    with open(args.metadata_map, "r") as fp:
        metadata_map = json.loads(fp.read())

    id_to_dna_mutation = {v: k for k, v in metadata_map["dna_mutation"].items()}
    id_to_gene_aa_mutation = {v: k for k, v in metadata_map["gene_aa_mutation"].items()}
    id_to_protein_aa_mutation = {
        v: k for k, v in metadata_map["protein_aa_mutation"].items()
    }

    # -----------------------------
    #     MUTATION DEFINITIONS
    # -----------------------------

    dna_mutation_def = (
        pd.Series(id_to_dna_mutation)
        .rename("mutation_str")
        .rename_axis("mutation_id")
        .to_frame()
    )
    dna_mutation_def.insert(
        0, "segment", dna_mutation_def["mutation_str"].apply(lambda x: x.split("|")[0])
    )
    dna_mutation_def.insert(
        1, "pos", dna_mutation_def["mutation_str"].apply(lambda x: int(x.split("|")[1]))
    )
    dna_mutation_def.insert(
        2, "ref", dna_mutation_def["mutation_str"].apply(lambda x: x.split("|")[2])
    )
    dna_mutation_def.insert(
        3, "alt", dna_mutation_def["mutation_str"].apply(lambda x: x.split("|")[3])
    )
    dna_mutation_def.insert(
        0,
        "mutation_name",
        dna_mutation_def["mutation_str"].apply(dna_mutation_to_name),
    )

    gene_aa_mutation_def = (
        pd.Series(id_to_gene_aa_mutation)
        .rename("mutation_str")
        .rename_axis("mutation_id")
        .to_frame()
    )
    gene_aa_mutation_def.insert(
        1, "gene", gene_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[0])
    )
    gene_aa_mutation_def.insert(
        2,
        "pos",
        gene_aa_mutation_def["mutation_str"].apply(lambda x: int(x.split("|")[1])),
    )
    gene_aa_mutation_def.insert(
        3, "ref", gene_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[2])
    )
    gene_aa_mutation_def.insert(
        4, "alt", gene_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[3])
    )
    gene_aa_mutation_def.insert(
        0,
        "mutation_name",
        gene_aa_mutation_def["mutation_str"].apply(aa_mutation_to_name),
    )

    protein_aa_mutation_def = (
        pd.Series(id_to_protein_aa_mutation)
        .rename("mutation_str")
        .rename_axis("mutation_id")
        .to_frame()
    )
    protein_aa_mutation_def.insert(
        1,
        "protein",
        protein_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[0]),
    )
    protein_aa_mutation_def.insert(
        2,
        "pos",
        protein_aa_mutation_def["mutation_str"].apply(lambda x: int(x.split("|")[1])),
    )
    protein_aa_mutation_def.insert(
        3,
        "ref",
        protein_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[2]),
    )
    protein_aa_mutation_def.insert(
        4,
        "alt",
        protein_aa_mutation_def["mutation_str"].apply(lambda x: x.split("|")[3]),
    )
    protein_aa_mutation_def.insert(
        0,
        "mutation_name",
        protein_aa_mutation_def["mutation_str"].apply(aa_mutation_to_name),
    )

    mutation_defs = {
        "dna": dna_mutation_def.to_dict(orient="index"),
        "gene_aa": gene_aa_mutation_def.to_dict(orient="index"),
        "protein_aa": protein_aa_mutation_def.to_dict(orient="index"),
    }

    out_dfs = {"dna": [], "gene_aa": [], "protein_aa": []}

    for group in group_mutation_frequencies.keys():
        for reference_name in group_mutation_frequencies[group].keys():
            for mutation_type in group_mutation_frequencies[group][
                reference_name
            ].keys():
                mutations = group_mutation_frequencies[group][reference_name][
                    mutation_type
                ]
                defs = [
                    mutation_defs[mutation_type][mutation["mutation_id"]]
                    for mutation in mutations
                ]
                group_mutation_frequencies[group][reference_name][mutation_type] = [
                    {**a, **b} for a, b in zip(mutations, defs)
                ]
                mutation_df = pd.DataFrame.from_dict(
                    group_mutation_frequencies[group][reference_name][mutation_type]
                )
                mutation_df.insert(1, "reference", reference_name)
                out_dfs[mutation_type].append(mutation_df)

    for mutation_type in ["dna", "gene_aa", "protein_aa"]:
        out_dfs[mutation_type] = pd.concat(
            out_dfs[mutation_type], axis=0, ignore_index=True
        )

    with open(args.output_json, "w") as fp:
        fp.write(json.dumps(group_mutation_frequencies))

    out_dfs["dna"].to_csv(args.output_dna_csv, index=False)
    out_dfs["gene_aa"].to_csv(args.output_gene_aa_csv, index=False)
    out_dfs["protein_aa"].to_csv(args.output_protein_aa_csv, index=False)


if __name__ == "__main__":
    main()
