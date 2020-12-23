# coding: utf-8

"""Combine all lineage assignments into one file

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import io
import pandas as pd

from pathlib import Path


def combine_lineages(lineages, lineage_out):
    # Dump all SNP chunks into a text buffer
    df_io = io.StringIO()
    for i, chunk in enumerate(lineages):

        # For some reason, snakemake likes to pass folders in
        # just skip these
        if not Path(chunk).is_file():
            continue

        with open(chunk, "r") as fp_in:
            for j, line in enumerate(fp_in):
                # Write the header of the first file
                # Or write any line that's not the header
                # (to avoid writing the header more than once)
                if (i == 0 and j == 0) or j > 0:
                    df_io.write(line)

    # Read the buffer into a dataframe, then discard the buffer
    df_io.seek(0)
    df = pd.read_csv(df_io)
    df_io.close()

    df.to_csv(lineage_out, index=False)
