# coding: utf-8

"""Assign RSV genotype (A or B)

Author: David Favela - Vector Engineering Team (dfavela@broadinstitute.org)
"""
import pandas as pd
import pysam


def assign_genotype(sam_file):
    samfile = pysam.AlignmentFile(sam_file, "r", check_sq=False)  # pylint: disable=no-member
    all_genotypes = []
    for read in samfile.fetch(until_eof=True):
        nameArr = read.query_name.split('_')
        if nameArr[0] == 'RSVA':
            print('Assigned A')
            all_genotypes.append([read.reference_name, "A", nameArr[1]])
        elif nameArr[0] == 'RSVB':
            print('Assigned B')
            all_genotypes.append([read.reference_name, "B", nameArr[1]])
        else:
            print("Reference not recognized for query: ", read.query_name, read.reference_name)

    samfile.close()

    genotype_df = pd.DataFrame.from_records(
        all_genotypes, columns=["Accession ID", "subtype", "genotype"]
    )

    return genotype_df
