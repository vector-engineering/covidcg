# coding: utf-8

import gzip
import json
import pandas as pd
import urllib.request

from flask_server.config import config
from flask_server.genes_and_proteins import load_genes_or_proteins
from flask_server.dna_snv import process_dna_snvs
from flask_server.aa_snv import process_aa_snvs

from flask_server.color import get_categorical_colormap

genes = load_genes_or_proteins("static_data/genes.json")
proteins = load_genes_or_proteins("static_data/proteins.json")


def load_data(url):

    with urllib.request.urlopen(url) as f:
        print("Download complete")
        f = gzip.decompress(f.read())
        print("Decompression complete")
        f = json.loads(f)
        print("Loaded into memory")

    case_data = (
        pd.DataFrame.from_dict(f["case_data"], orient="columns")
        .set_index("Accession ID")
        .rename(
            columns={
                "dna_snp_str": "dna_snp",
                "gene_aa_snp_str": "gene_aa_snp",
                "protein_aa_snp_str": "protein_aa_snp",
            }
        )
    )
    case_data["collection_date"] = pd.to_datetime(case_data["collection_date"])
    # print(case_data)

    # Load metadata map
    metadata_map = f["metadata_map"]
    for meta_col in config["metadata_cols"].keys():
        metadata_map[meta_col] = {int(k): v for k, v in metadata_map[meta_col].items()}
    # print(list(metadata_map.keys()))

    # Load location map
    location_map = pd.DataFrame(f["location_map"])
    # print(location_map)

    # Load global group counts
    global_group_counts = f["global_group_counts"]
    # Convert SNV keys from strings to integers
    global_group_counts["dna_snp"] = {
        int(k): v for k, v in global_group_counts["dna_snp"].items()
    }
    global_group_counts["gene_aa_snp"] = {
        int(k): v for k, v in global_group_counts["gene_aa_snp"].items()
    }
    global_group_counts["protein_aa_snp"] = {
        int(k): v for k, v in global_group_counts["protein_aa_snp"].items()
    }
    # print(global_group_counts)

    # Group consensus SNVs
    group_consensus_snvs = f["group_consensus_snps"]
    # print(group_consensus_snps)

    # Build colormaps
    sequence_groups = list(group_consensus_snvs.keys())
    group_colormaps = dict()
    for group in sequence_groups:
        group_colormaps[group] = get_categorical_colormap(
            list(group_consensus_snvs[group].keys())
        )

    dna_snp = process_dna_snvs(metadata_map["dna_snp"])
    # print(dna_snp)

    gene_aa_snp = process_aa_snvs(metadata_map["gene_aa_snp"], "gene", genes)
    protein_aa_snp = process_aa_snvs(
        metadata_map["protein_aa_snp"], "protein", proteins
    )
    # print(gene_aa_snp)
    # print(protein_aa_snp)

    return {
        "case_data": case_data,
        "metadata_map": metadata_map,
        "location_map": location_map,
        "global_group_counts": global_group_counts,
        "group_consensus_snvs": group_consensus_snvs,
        "group_colormaps": group_colormaps,
        "dna_snp": dna_snp,
        "gene_aa_snp": gene_aa_snp,
        "protein_aa_snp": protein_aa_snp,
        "data_date": f["data_date"],
        "num_sequences": len(case_data),
        "country_score": f["country_score"],
        "geo_select_tree": f["geo_select_tree"],
    }

