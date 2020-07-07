#!/usr/bin/env python3
# coding: utf-8

"""Process acknowledgement files

Author: Albert Chen (Deverman Lab - Broad Institute)
"""

import numpy as np
import pandas as pd

from pathlib import Path

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (
    project_root_path / "data"
).resolve()  # Resolve any symlinks --> absolute path


def process_ack():
    """COLLECT ACKNOWLEDGEMENTS	
    """

    ack_files = sorted((data_dir / "acknowledgements").glob("*.xls"))
    print(
        "Collecting {} acknowledgement files...".format(len(ack_files)),
        end="",
        flush=True,
    )
    ack_df = pd.DataFrame()
    for f in ack_files:
        _df = pd.read_excel(f, skiprows=[0, 1, 3])
        ack_df = pd.concat([ack_df, _df], ignore_index=True)

    print("done")

    print("Factorizing acknowledgements into IDs...", end="", flush=True)
    code, uniques = pd.factorize(
        list(
            zip(ack_df["Originating lab"], ack_df["Submitting lab"], ack_df["Authors"])
        )
    )

    # Create map of ID -> acknowledgement entry
    unique_ack_df = pd.DataFrame.from_records(
        uniques, columns=["Originating lab", "Submitting lab", "Authors"]
    )

    # Append code to acknowledgement dataframe, and take
    # subset of columns
    ack_df = pd.concat([ack_df, pd.Series(code, name="ack_id")], axis=1)[
        ["Accession ID", "ack_id"]
    ]
    ack_df = ack_df.set_index("Accession ID")
    # Cast ack_id to integer
    ack_df["ack_id"] = ack_df["ack_id"].astype(int)

    print("done")
    # print(ack_df)

    print("Saving acknowledgement files...", end="", flush=True)

    # ack_df.to_csv(data_dir / 'taxon_acknowledgements.csv', header=['ack_id'], index_label='gisaid_id')
    # ack_df.to_json(data_dir / 'taxon_acknowledgements.json')

    # Drop MultiIndex to columns, make index the real index again
    unique_ack_df.to_csv(
        data_dir / "acknowledgement_map.csv", index=True, index_label="ack_id"
    )
    unique_ack_df.to_json(data_dir / "acknowledgement_map.json", orient="index")

    print("done", flush=True)

    return ack_df


if __name__ == "__main__":
    process_ack()
