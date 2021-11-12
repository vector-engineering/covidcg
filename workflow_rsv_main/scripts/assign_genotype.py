# coding: utf-8

"""Assign RSV genotype (A or B)

Author: David Favela - Vector Engineering Team (dfavela@broadinstitute.org)
"""
import pandas as pd
import pysam


def assign_genotype(sam_file):
    samfile = pysam.AlignmentFile(sam_file, "r")  # pylint: disable=no-member

    all_genotypes = []
    for read in samfile.fetch(until_eof=True):
        # Skip if unmapped
        if read.is_unmapped:
            continue

        if read.reference_name == "NC_038235.1":
            all_genotypes.append([read.query_name, "A"])
        elif read.reference_name == "NC_001781.1":
            all_genotypes.append([read.query_name, "B"])

    samfile.close()

    genotype_df = pd.DataFrame.from_records(
        all_genotypes, columns=["Accession ID", "genotype"]
    )

    return genotype_df
