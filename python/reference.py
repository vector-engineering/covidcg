# coding: utf-8

"""Load and parse the Wuhan-Hu-1 reference genome

Author: Albert Chen (Deverman Lab - Broad Institute)
"""

import numpy as np
import pandas as pd

from pathlib import Path

from fasta import read_fasta_file
from util import translate, data_dir, static_data_dir

# Load the reference sequence
ref_fasta_path = static_data_dir / "reference.fasta"
with ref_fasta_path.open("r") as fp:
    lines = fp.readlines()
    ref = read_fasta_file(lines)
    ref_seq = list(ref.values())[0]

# Load genes
genes_path = static_data_dir / "genes.csv"
genes_df = pd.read_csv(genes_path)

# Write to JSON
genes_df.to_json(static_data_dir / "genes.json", orient="records")

# Dict of gene: (start, end) nucleotide positions
# Positions are 1-indexed and inclusive, and ranges
# are inclusive [start, end]
genes = {}
for i, gene in genes_df.iterrows():
    # Skip non-protein-coding
    if gene["protein_coding"] == 0:
        continue

    start = gene["start"]
    end = gene["end"]

    genes[gene["gene"]] = (start, end)

# Reference translated genes
# From Wuhan-Hu-1, NCBI: NC_045512.2
# With stop codons added onto the ends
gene_aa = {}
for gene, rnge in genes.items():
    # [start, end], so end = end + 1
    # and because ranges are 1-indexed,
    # start = start - 1 and end = end - 1
    # so the range is [start - 1, end)
    gene_aa[gene] = translate(ref_seq[(rnge[0] - 1) : rnge[1]])

# Load proteins
proteins_path = static_data_dir / "proteins.csv"
proteins_df = pd.read_csv(proteins_path, comment="#")

# Write to JSON
proteins_df.to_json(static_data_dir / "proteins.json", orient="records")

