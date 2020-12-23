# coding: utf-8

"""Abbreviate author strings

2020-12-21: DEPRECATED

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import numpy as np
import pandas as pd


def abbreviate_authors(authors):
    """Process acknowledgement files
    """

    # Create an abbreviated authors column - with (First author) Last name et al.

    # Default to the most common format - First-Name Last-Name, ...
    # Exclude commas/semicolons since we want to capture all characters, including weird utf-8 ones
    # Also have to specify both the english and chinese comma characters
    authors_abbrev = (
        authors.str.extract(r"^(?:[^,，;\s]+\s)*([^,，;]+)[,，;]")
        .squeeze()  # DataFrame of one column -> Series
        .str.strip()
    )
    # If the format is Last-Name Initials, ...
    # i.e., FirstName L-N, or FirstName L.N, or FirstName L.N., or FirstName LN
    # Then extract the last name
    last_name_initial_inds = authors.str.contains(
        r"^[^,，;]+[\s,，;][A-Z]{1}(?:[-\.]{0,2}[A-Z])?\.?[,，;]"
    ).fillna(False)
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
    et_al_inds = authors.str.contains(r"etl? al\.?$").fillna(False)
    ack_map.loc[et_al_inds, "authors_abbrev"] = (
        ack_map.loc[et_al_inds, "Authors"]
        .str.extract(r"^(.+)etl? al\.$")
        .squeeze()
        .str.strip()
    )

    # Capitalize name and add "et al."
    authors_abbrev = authors_abbrev.str.capitalize() + " et al."

    # If the format is one name: "LastName, F."
    one_name_first_initial_inds = authors.str.contains(
        r"^([\w]+),?(?:[^,，;]+)?[A-Za-z]{1}[-\.]{0,2}[A-Za-z]?\.?$"
    ).fillna(False)
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
        authors.str.contains(r"^(?:[^,，;]{2,})\s([^,，;]{2,})$").fillna(False)
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
    authors_abbrev = authors_abbrev.combine_first(authors)

    return authors_abbrev
