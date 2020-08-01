#!/usr/bin/env python3
# coding: utf-8

"""Process acknowledgement files

Author: Albert Chen (Deverman Lab - Broad Institute)
"""

import numpy as np
import pandas as pd


def process_ack(ack_files):
    """Process acknowledgement files
    """

    print(
        "Loading {} acknowledgement files...".format(len(ack_files)),
        end="",
        flush=True,
    )
    ack_df = pd.DataFrame()
    for f in ack_files:
        df = pd.read_excel(f, skiprows=[0, 1, 3])
        ack_df = pd.concat([ack_df, df], ignore_index=True)

    # Strip strings
    ack_df["Originating lab"] = ack_df["Originating lab"].str.strip()
    ack_df["Submitting lab"] = ack_df["Submitting lab"].str.strip()
    ack_df["Authors"] = ack_df["Authors"].str.strip()

    print("done")

    print("Factorizing acknowledgements into IDs...", end="", flush=True)
    code, uniques = pd.factorize(
        list(
            zip(ack_df["Originating lab"], ack_df["Submitting lab"], ack_df["Authors"])
        )
    )

    # Create map of ID -> acknowledgement entry
    ack_map = pd.DataFrame.from_records(
        uniques, columns=["Originating lab", "Submitting lab", "Authors"]
    )

    # Create an abbreviated authors column - with (First author) Last name et al.

    # Default to the most common format - First-Name Last-Name, ...
    # Exclude commas/semicolons since we want to capture all characters, including weird utf-8 ones
    # Also have to specify both the english and chinese comma characters
    ack_map["authors_abbrev"] = (
        ack_map["Authors"]
        .str.extract(r"^(?:[^,，;\s]+\s)*([^,，;]+)[,，;]")
        .squeeze()  # DataFrame of one column -> Series
        .str.strip()
    )
    # If the format is Last-Name Initials, ...
    # i.e., FirstName L-N, or FirstName L.N, or FirstName L.N., or FirstName LN
    # Then extract the last name
    last_name_initial_inds = (
        ack_map["Authors"]
        .str.contains(r"^[^,，;]+[\s,，;][A-Z]{1}(?:[-\.]{0,2}[A-Z])?\.?[,，;]")
        .fillna(False)
    )
    ack_map.loc[last_name_initial_inds, "authors_abbrev"] = (
        ack_map.loc[last_name_initial_inds, "Authors"]
        .str.extract(
            # This is a mess. I can't believe it kind of works.
            # TODO: Split this up into more specific cases
            r"^([\w’\'\-]+\s?[\w’\'\-]*)[,;\s](?:[^,，;]+)?[A-Z]{1}[-\.]{0,2}[A-Z]?\.?[,，;]"
        )
        .squeeze()
        .str.strip()
    )

    # If the format is already in "LastName et al", then keep the Authors
    et_al_inds = ack_map["Authors"].str.contains(r"etl? al\.?$").fillna(False)
    ack_map.loc[et_al_inds, "authors_abbrev"] = (
        ack_map.loc[et_al_inds, "Authors"]
        .str.extract(r"^(.+)etl? al\.$")
        .squeeze()
        .str.strip()
    )

    # Capitalize name and add "et al."
    ack_map["authors_abbrev"] = ack_map["authors_abbrev"].str.capitalize() + " et al."

    # If the format is one name: "LastName, F."
    one_name_first_initial_inds = (
        ack_map["Authors"]
        .str.contains(r"^([\w]+),?(?:[^,，;]+)?[A-Za-z]{1}[-\.]{0,2}[A-Za-z]?\.?$")
        .fillna(False)
    )
    # ack_map.loc[one_name_first_initial_inds, "authors_abbrev"] = (
    #     ack_map.loc[one_name_first_initial_inds, "Authors"]
    #     .str.extract(r"^([\w]+),?(?:[^,，;]+)?[A-Za-z]{1}[-\.]{0,2}[A-Za-z]?\.?$")
    #     .squeeze()
    #     .str.strip()
    # )
    ack_map.loc[one_name_first_initial_inds, "authors_abbrev"] = ack_map.loc[
        one_name_first_initial_inds, "Authors"
    ]

    # If the format is one name: "FirstName LastName"
    # I guess I can't know whether the first string is the first name or not
    # without making assumptions based on the names themselves
    # But I'll assume english naming convention here (extract the last part as last name)
    one_name_first_last_inds = (
        ack_map["Authors"].str.contains(r"^(?:[^,，;]{2,})\s([^,，;]{2,})$").fillna(False)
    ) & (~one_name_first_initial_inds)
    # ack_map.loc[one_name_first_last_inds, "authors_abbrev"] = (
    #     ack_map.loc[one_name_first_last_inds, "Authors"]
    #     .str.extract(r"^(?:[^\.,，;]{2,})\s([^\.,，;]{2,})$")
    #     .squeeze()
    #     .str.strip()
    # )
    ack_map.loc[one_name_first_last_inds, "authors_abbrev"] = ack_map.loc[
        one_name_first_last_inds, "Authors"
    ]

    # Default to the original Author's column if the abbreviation is empty
    ack_map["authors_abbrev"] = ack_map["authors_abbrev"].combine_first(
        ack_map["Authors"]
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

    return ack_df, ack_map


# if __name__ == "__main__":
#     process_ack()