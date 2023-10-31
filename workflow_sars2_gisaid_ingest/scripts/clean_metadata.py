#!/usr/bin/env python3
# coding: utf-8

"""Clean metadata
Merger of old scripts for cleaning patient metadata, sequencing metadata, and 
acknowledgements separately. Now they're all merged into this big file

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import numpy as np
import pandas as pd
import re

from process_location_metadata import process_location_metadata


def clean_name_metadata(df):
    df["virus_name"] = df["covv_virus_name"].astype(str).str.strip()
    return df


def clean_host_metadata(df):
    """Clean host metadata"""
    # print("Cleaning host metadata", end="", flush=True)
    df["host"] = df["covv_host"].astype(str).str.strip()
    # In the future - collapse common + taxonomy names?

    # print("done")
    return df


def clean_gender_metadata(df):
    """Clean patient gender metadata"""

    # print("Cleaning patient gender metadata...", end="", flush=True)

    # Make a copy, strip whitespace
    df["gender"] = df["covv_gender"].fillna("Unknown").astype(str).str.strip()

    # replace_map = [
    #     (r"^female", "Female", False),
    #     (r"^male", "Male", False),
    #     (r"^f$", "Female", False),
    #     (r"^m$", "Male", False),
    #     (r"unknown", "Unknown", False),
    # ]

    # for pair in replace_map:
    #     df["gender"] = df["gender"].str.replace(
    #         pair[0], pair[1], -1, pair[2] if len(pair) > 2 else True
    #     )

    # gender_key_map = {
    #     "Female": ["Woman", "Femal", "Famle", "Famale", "Fenale", "Femele"],
    #     "Male": ["M"],
    #     "Unknown": [
    #         "unknwon",
    #         "Unkown",
    #         "U",
    #         "nan",
    #         "Not reported to protect privacy",
    #         "unknow",
    #         "?",
    #     ],
    # }

    # gender_map = {}
    # for k, v in gender_key_map.items():
    #     # Add self
    #     gender_map[k] = k
    #     for _v in v:
    #         gender_map[_v] = k

    # df["gender"] = (
    #     df["gender"].map(gender_map).combine_first(df["gender"]).fillna("Unknown")
    # )

    # print("done")
    # print('"Gender" values: {}'.format(", ".join(df["gender"].unique())))
    # df["gender"] = df["gender"].fillna("Unknown")

    return df


def clean_age_metadata(df):
    """Clean patient age metadata

    For each age value we want to define an age range that the value
    corresponds to. This is necessary since the metadata is provided in
    various different specificities.

    i.e., exact age (72.343), year (42), or age range (20-30)

    Define a range [start, end), for each age
    This ranges can then be filtered over in a way that includes as much data
    as possible

    """

    # print("Cleaning patient age metadata...", end="", flush=True)

    # Do some basic cleanup before we start
    df["age_clean"] = df["covv_patient_age"].astype(str)
    df["age_clean"] = df["age_clean"].fillna("Unknown")
    df["age_clean"] = df["age_clean"].str.strip()

    df["age_start"] = np.nan
    df["age_end"] = np.nan

    # for i, v in df["age_clean"].iteritems():

    #     # Skip if Unknown
    #     if (
    #         v
    #         in [
    #             "Unknown",
    #             "unknown",
    #             "unkown",
    #             "unknwon",
    #             "unavailable",
    #             "uknown",
    #             "no data",
    #             "Male",
    #             "uknown",
    #         ]
    #         or not v
    #     ):
    #         # df.loc[i, 'Patient age clean'] = 'Unknown'
    #         continue

    #     # Parse "Adult" and "Over 18" to be 18 -> 100
    #     elif v in ["Over 18", "over 18", "> 18", "Adult"]:
    #         df.loc[i, "age_start"] = 18.0
    #         df.loc[i, "age_end"] = 100.0

    #     # Parse clean integers. e.g., "42"
    #     elif re.match(r"^[0-9]+$", v):
    #         df.loc[i, "age_start"] = float(int(v))
    #         df.loc[i, "age_end"] = float(int(v) + 1)

    #     # Parse fractions. e.g., "72.323"
    #     elif re.match(r"^[0-9]*\.[0-9]+$", v):
    #         df.loc[i, "age_start"] = float(v)
    #         df.loc[i, "age_end"] = float(v)

    #     # Parse N months (less than 1 years old). e.g., "7 months"
    #     elif m := re.match(r"^([1]*[0-9])+\smonth(s)?$", v, re.IGNORECASE):
    #         month = int(m.groups()[0])
    #         # Convert months to fraction of a year, then round to an integer
    #         df.loc[i, "age_start"] = month / 12.0
    #         df.loc[i, "age_end"] = (month + 1) / 12.0

    #     # Same, but days. e.g., "17 days"
    #     elif m := re.match(r"^([1]*[0-9])+\sday(s)?$", v, re.IGNORECASE):
    #         day = int(m.groups()[0])
    #         # Convert days to fraction of a year, then round to an integer
    #         df.loc[i, "age_start"] = day / 365.0
    #         df.loc[i, "age_end"] = (day + 1) / 365.0

    #     # Same, but weeks. e.g., "6 weeks"
    #     elif m := re.match(r"^([1]*[0-9])+\sweek(s)?$", v, re.IGNORECASE):
    #         week = int(m.groups()[0])
    #         # Convert weeks to fraction of a year, then round to an integer
    #         df.loc[i, "age_start"] = week / 52.0
    #         df.loc[i, "age_end"] = (week + 1) / 52.0

    #     # Remove '-year old', 'years' or 'age' at end. e.g., "44-year old"
    #     # None of the years/age entries are fractions, so treat as a valid integer
    #     # This might easily break in the future
    #     elif m := re.match(r"^([0-9]+)\s?(-year\sold|years|age)$", v, re.IGNORECASE):
    #         year = int(m.groups()[0])
    #         df.loc[i, "age_start"] = float(year)
    #         df.loc[i, "age_end"] = float(year + 1)

    #     # Match year-month, e.g., "6 years 2 months"
    #     elif m := re.match(
    #         r"^([0-9]+)(?:,\s|\syears\s)([1]?[0-9]+)\smonths$", v, re.IGNORECASE
    #     ):
    #         # Extract years, months
    #         years = int(m.groups()[0])
    #         months = int(m.groups()[1])
    #         # Round to nearest year
    #         df.loc[i, "age_start"] = years + (months / 12.0)
    #         df.loc[i, "age_end"] = years + ((months + 1) / 12.0)

    #     # Match "(year), (x) month". e.g., "62, 1 month"
    #     elif m := re.match(r"^([0-9]+),\s([0-9]+)\smonth$", v, re.IGNORECASE):
    #         # Extract years, months
    #         years = int(m.groups()[0])
    #         months = int(m.groups()[1])
    #         # Round to nearest year
    #         df.loc[i, "age_start"] = years + (months / 12.0)
    #         df.loc[i, "age_end"] = years + ((months + 1) / 12.0)

    #     # Extract decade ranges. e.g., "30s"
    #     elif m := re.match(r"^([0-9]+)\'?s$", v):
    #         decade = int(m.groups()[0])
    #         df.loc[i, "age_start"] = float(decade)
    #         df.loc[i, "age_end"] = float(decade + 10)

    #     # Extract year ranges, e.g., "10-20"
    #     elif m := re.match(r"^([0-9]+)\s?[-|–]\s?([0-9]+)$", v):
    #         start = int(m.groups()[0])
    #         end = int(m.groups()[1])
    #         df.loc[i, "age_start"] = float(start)
    #         df.loc[i, "age_end"] = float(
    #             end + 1
    #         )  # Assume that the provided range is [start, end]

    #     # Extract year ranges with text separator, e.g., "10 to 20"
    #     elif m := re.match(r"^([0-9]+)\s?to\s?([0-9]+)$", v, re.IGNORECASE):
    #         start = int(m.groups()[0])
    #         end = int(m.groups()[1])
    #         df.loc[i, "age_start"] = float(start)
    #         df.loc[i, "age_end"] = float(
    #             end + 1
    #         )  # Assume that the provided range is [start, end]

    #     # Extract inequalities, e.g., ">60" or "over 60"
    #     elif m := re.match(r"^(?:>|over)\s?([0-9]+)$", v, re.IGNORECASE):
    #         start = int(m.groups()[0])

    #         # Assume that the range is just one year
    #         df.loc[i, "age_start"] = float(start)
    #         df.loc[i, "age_end"] = float(start + 1)

    #     # Extract <18. This could be anywhere from 0 - 18
    #     elif re.match(r"^<\s?18$", v):
    #         df.loc[i, "age_start"] = 0.0
    #         df.loc[i, "age_end"] = 18.0

    #     # Extract "(year)+" or "(year)x". i.e., "30+" or "30x"
    #     # Assume that the age_end is just the next year
    #     elif m := re.match(r"^([0-9]+)[\+|x]$", v):
    #         start = int(m.groups()[0])
    #         # Assume that the range is just one year
    #         df.loc[i, "age_start"] = float(start)
    #         df.loc[i, "age_end"] = float(start + 1)

    #     # Extract "(year) plus"
    #     # Assume that the age_end is just the next year
    #     elif m := re.match(r"^([0-9]+)\s?plus$", v, re.IGNORECASE):
    #         start = int(m.groups()[0])
    #         # Assume that the range is just one year
    #         df.loc[i, "age_start"] = float(start)
    #         df.loc[i, "age_end"] = float(start + 1)

    #     # Extract "(x) Year"
    #     # Assume that the age_end is just the next year
    #     elif m := re.match(r"^([0-9]+)\syears?$", v, re.IGNORECASE):
    #         start = int(m.groups()[0])
    #         # Assume that the range is just one year
    #         df.loc[i, "age_start"] = float(start)
    #         df.loc[i, "age_end"] = float(start + 1)

    # print("done")
    # print(
    #     'Could not parse the following "Patient Age" values: {}'.format(
    #         ", ".join(
    #             [
    #                 '"{}"'.format(x)
    #                 for x in df["age_clean"][pd.isnull(df["age_start"])]
    #                 .astype(str)
    #                 .str.strip()
    #                 .unique()
    #                 .astype(str)
    #             ]
    #         )
    #     )
    # )

    return df


def clean_patient_status_metadata(df):
    # print("Cleaning patient status metadata...", end="", flush=True)

    # Strip whitespace
    df["patient_status"] = (
        df["covv_patient_status"].fillna("Unknown").astype(str).str.strip()
    )

    # replace_map = [
    #     (r"hospitalized", "Hospitalized", False),
    #     (r"^fever", "Fever", False),
    #     (r"^live", "Live", False),
    #     (r"outpatient", "Outpatient", False),
    #     (r"released", "Released", False),
    #     (r"asymptomatic", "Asymptomatic", False),
    #     (r"unknown", "Unknown", False),
    # ]

    # for pair in replace_map:
    #     df["patient_status"] = df["patient_status"].str.replace(
    #         pair[0], pair[1], -1, pair[2] if len(pair) > 2 else True
    #     )

    # status_key_map = {
    #     "Unknown": [
    #         "unknwon",
    #         "Unkown",
    #         "unkown",
    #         "unknow",
    #         "Oro-pharyngeal swab",
    #         "Z039",
    #         "-",
    #         "Unknow",
    #         "uknown",
    #         "\ufeffUnknown",
    #         "uncknown",
    #     ],
    #     "Hospitalized": [
    #         "Hospitaized",
    #         "Hospitalized patient",
    #         "Hopsitalized",
    #         "Hospitalised",
    #         "Still Hospitalized",
    #     ],
    #     "Symptomatic": ["symptomatic"],
    #     "Fever": ["fever"],
    #     "Intensive Care Unit": [
    #         "Intensive Care Unit",
    #         "Hospitalized in ICU",
    #         "ICU",
    #         "Hospitalized (Intensive care unit)",
    #     ],
    #     "Mild Case": ["Mild case", "Mild"],
    #     "Deceased": ["Death", "deceased",],
    # }

    # status_map = {}
    # for k, v in status_key_map.items():
    #     # Add self
    #     status_map[k] = k
    #     for _v in v:
    #         status_map[_v] = k

    # df["patient_status"] = (
    #     df["patient_status"]
    #     .map(status_map)
    #     .combine_first(df["patient_status"])
    #     .fillna("Unknown")
    # )

    # print("done")
    # print(
    #     '"Patient Status" values: {}'.format(
    #         ", ".join(['"{}"'.format(x) for x in df["patient_status"].unique()])
    #     )
    # )

    return df


def clean_passage_metadata(df):
    """Clean cell passage metadata"""

    print("Cleaning cell passage metadata...", end="", flush=True)

    # Basic cleaning
    df["passage"] = df["covv_passage"].astype(str).str.strip()

    passage_key_map = {
        "Original": [
            "Orginal",
            "Orignal",
            "Origional",
            "Original (first sample)",
            "original",
            "Orignial",
            "Oroginal",
            "Orginial",
            "Original, from nasopharyngeal aspirate",
            "Nasopharyngeal",
            "ORIGINAL",
            "origin",
            "Original,",
            "Original isolate",
            "Origina",
            "ORiginal",
            "Origin",
            "Direct from swab",
            "Original Isolate",
            "Original (Clinical Sample)",
            "Original (Clinical sample)",
            "Original(Clinical Sample)",
            "Clinical Specimen",
            "Clinical specimen",
            "Clinicial specimen",
            "Originzl",
            "Originaal",
            "Orig9inal",
            "Origginal",
            "Origianl",
            "Nasopharyngeal swab",
            "Orifinal",
            "Oiginal",
            "Oringinal",
            "Clinical nasopharyngeal swab specimen",
            "Originla",
            "Oriiginal",
            "Orriginal",
            "Origina;",
            "Originale",
            "origianl",
            "Original, clinical sample",
            "Originall",
            "origina",
            "Human",
            "direct sequencing",
            "e.g. Original",
            "Original isolate isolate",
            "Original swab",
            "Originalo",
            "joriginal",
            "originale",
            "Original sample in syrian hamster_P2",
            "Original sample in syrian hamster_P1",
            "Originjal",
            "New variant",
            "Priginal",
            "Criblage sauvage",
            "Original.",
            "Oriignal",
            "originall",
            "isolate",
            "Oiriginal",
            "Nasal swab",
        ]
        # "Vero": ["Vero cells"],
        # "Vero P1": [
        #     "Vero p1",
        #     "Vero CCL81 isolate P1",
        #     "Vero cell P1",
        #     "Vero C1",
        #     "Vero 1",
        #     "Vero1",
        # ],
        # "Vero P2": ["Vero cell P2", "P2, Vero", "Vero2"],
        # "Vero P3": ["Vero cell P3"],
        # "Vero P4": ["Vero cell P4"],
        # "LLC-MK2": [],
        # "P1": ["Virus Isolate, Passage 1", "C1", "Viral culture P1", "V1"],
        # "P2": ["Passage 2"],
        # "P3": [],
        # "P4": [],
        # "P5": [],
        # "Vero E6": [
        #     "VeroE6",
        #     "Vero-E6",
        #     "VeroE6 cells",
        #     "Initial isolation on VeroE6 cells",
        # ],
        # "Vero E6 P1": [
        #     "VeroE6/P1",
        #     "Vero E6, 1st passage",
        #     "P1/Vero E6",
        #     "Vero E6 #1",
        #     "Pi/Vero E6",
        #     "Vero E6, 1 passage",
        #     "Vero E6 P1, mouse-P14",
        #     "P1 Vero E6",
        #     "Passage 1 / Vero-E6",
        #     "P1 / Vero-E6",
        #     "P1/ Vero-E6",
        #     "VeroE6 P1",
        # ],
        # "Vero E6 P2": [
        #     "VeroE6/P2",
        #     "VERO E6 / P2",
        #     "Vero E6 cells, P2",
        #     "Vero E6 / P2",
        #     "Vero E6/P2",
        #     "VeroE6 P2",
        #     "Vero C2",
        # ],
        # "Vero E6 P3": ["VeroE6/P3", "Vero E6-P3", "VeroE6, passage 3"],
        # "Vero E6 P4": ["Passage 4 in Vero E6 cells"],
        # "Vero/hSLAM P1": [],
        # "Vero E6/TMPRSS2": ["VeroE6/TMPRSS2", "Vero E6 / TMPRSS2"],
        # "Caco-2": ["Caco2"],
        # "Caco-2 P1": [],
        # "Caco-2 P2": ["Passage 2, Caco2"],
        # "Culture (unknown)": ["Virus culture", "Viral culture"],
    }

    passage_map = {}
    for k, v in passage_key_map.items():
        # Add self
        passage_map[k] = k
        for _v in v:
            passage_map[_v] = k

    print("done")
    print(
        'Unmapped "Passage" values: {}'.format(
            ", ".join(
                [
                    '"{}"'.format(x)
                    for x in df["covv_passage"][
                        pd.isnull(df["passage"].map(passage_map))
                    ]
                    .astype(str)
                    .str.strip()
                    .unique()
                    .astype(str)
                ]
            )
        )
    )

    df["passage"] = (
        (df["passage"].map(passage_map)).combine_first(df["passage"]).fillna("Unknown")
    )

    # Only keep rows that are Original (throw out unknown as well)
    original_passage = df["passage"].isin(
        ["Original", "Original (Clinical Sample)", "Clinical Specimen"]
    )
    print("Excluding {} isolates passaged in culture".format((~original_passage).sum()))
    df = df.loc[original_passage, :]

    return df


def clean_specimen_metadata(df):
    # print("Cleaning specimen metadata...", end="", flush=True)

    # Basic cleanup
    df.loc[:, "specimen"] = (
        df["covv_specimen"].fillna("Unknown").astype(str).str.strip()
    )

    # specimen_key_map = {
    #     "Alveolar lavage fluid": [],
    #     "Anal swab": ["Anus", "Anal swabs"],
    #     "Aspirate": [],
    #     "Autopsy, Bronch and lung": [],
    #     "Blood": ["blood sample", "Blood sample"],
    #     "Breathing air using VIVAs Air sampler": [],
    #     "Bronchial scraper": [],
    #     "Bronchoalveolar lavage fluid": [
    #         "bronchoalveolar lavage fluid",
    #         "Bronchoalveolar Lavage NCIT:C51913",
    #         "Bronchoalveolar lavage",
    #         "Bronchoalveolar fluid",
    #         "Broncho-alveolar lavage",
    #         "Bronchoalveolar lavage fluid (BALF)",
    #         "Bronchoalveolar-lavage fluid",
    #         "Broncho-alveolar fluid",
    #         "Broncoalveolar",
    #         "BAL",
    #     ],
    #     "Bronchial secretion": [],
    #     "Buccal swab": ["buccal swab"],
    #     "Deep throat saliva": [],
    #     "Dry swab": [],
    #     "EDTA-Plasma": [],
    #     "Endotracheal aspirate": [
    #         "endotracheal aspirates",
    #         "Endotracheal aspirate (ETA)",
    #         "endotracheal aspirate (ETA)",
    #         "Endotracheales Aspirat",
    #     ],
    #     "Enviromental swab": ["Door handle", "Air", "Air sample", "Home environment"],
    #     "Fecal swab": [],
    #     "Intestine tissue": [],
    #     "Lung tissue": ["Lung", "lung tissue", "lung biopsy"],
    #     "Lung & intenstine tissue": ["Mixed lung and intestine tissue"],
    #     "Lung secretion": ["Lung secretion"],
    #     "Lung swab": ["lung swab"],
    #     "Mid-turbinate nasal swab": ["Mid-Turbinate nasal swab"],
    #     "Muscle tissue": ["muscle tissue"],
    #     "Nasal swab": [
    #         "Nose",
    #         "Nose swab",
    #         "nasal swab",
    #         "Nasal Swab",
    #         "Mid-nasal swab",
    #         "nose swab",
    #         "Nasal swab specimen",
    #     ],
    #     "Nasal swab & blood": ["Nasal Swab, Blood"],
    #     "Nasal swab & serum": ["Nasal Swab, Serum", "Nasal Swab, Serum,"],
    #     "Nasal swab, serum, blood": [
    #         "Nasal swab, Serum, Blood",
    #         "Nasal Swab, Blood, Serum",
    #         "Nasal Swab, Serum, Blood",
    #         "Nasal Swab, Serum, Blood,",
    #     ],
    #     "Nasal swab, serum, blood, plasma": [
    #         "Nasal Swab, Blood, Serum, Plasma",
    #         "Nasal Swab, Blood, Serum, Plasma",
    #         "Nasal Swab, Serum, Blood, Plasma",
    #     ],
    #     "Nasal swab, serum, blood, sputum": ["Nasal Swab, Serum, Blood, Sputum"],
    #     "Nasal swab, serum, blood, sputum, urine": [
    #         "Nasal Swab, Serum, Blood, Sputum, Urine"
    #     ],
    #     "Nasal & throat swab": [
    #         "Nose/throat swab",
    #         "Nost and Throad swab",
    #         "Nose and Throat swab",
    #         "Nose-Throat",
    #         "Throat swab Nasal swab",
    #         "Mixture of Throat swab Nasal swab",
    #         "Mixture of Throat swab and Nasal swab",
    #         "Mixture of throat and nasal swab",
    #         "Nasal swab and Throat swab",
    #         "Nasal and throat swab",
    #         "Throat and Nasal Swab",
    #         "Combined nasal and throat swab",
    #     ],
    #     "Nasal & oro-pharyngeal swab": [
    #         "Nasal, oropharyngeal swab",
    #         "Oro-pharyngeal swab, Nasal swab",
    #     ],
    #     "Nasal swab, oral swab, tracheal wash": [
    #         "Oral swab; Nasal swab; Tracheal wash"
    #     ],
    #     "Nasal washing": ["Nasal Washing"],
    #     "Nasopharyngeal swab": [
    #         "nasopharyngeal swab",
    #         "Nasopharyngeal swab",
    #         "Naso-pharyngeal swab",
    #         "naso-pharyngeal swab",
    #         "Nasopharyngeal swap",
    #         "Nasopharyngeal Swab",
    #         "Nasopharynx swab",
    #         "Nasopharingeal swab",
    #         "NP Swab",
    #         "Nasopharnygeal swab",
    #         "NP swab",
    #         "nasopharyngial swab",
    #         "NP",
    #         "Naso-pharingeal swab",
    #         "Naso-pharygeal swab",
    #         "Nasopharynx",
    #         "nasopharingeal swab",
    #         "Rhino-pharyngeal swab",
    #         "Upper Resp NPS",
    #         "Nasal-pharyngeal swab",
    #         "Nasopharyngeal",
    #         "Nasopharygeal Swab",
    #         "Nasopharyngeal swa",
    #     ],
    #     "Nasopharyngeal VTM": [],
    #     "Nasopharyngeal & oro-pharyngeal swab": [
    #         "Nasopharyngeal swab and oropharyngeal swab",
    #         "Nasopharyngeal swab, oropharyngeal swab",
    #         "Combined nasopharyngeal and oropharyngeal swab",
    #         "Nasopharyngeal/oropharyngeal swab",
    #         "Naso-Oropharyngeal swab",
    #         "Swab nasopharyngeal oropharyngeal",
    #         "Naso and oro-pharyngeal swabs",
    #         "Oro and naso-pharyngeal swabs",
    #         "Oro and naso-pharyngeal swab",
    #         "Nasopharyngeal / oropharyngeal swab",
    #         "Nasopharyngeal and Oropharyngeal swab",
    #         "Nasopharyngeal / Oro-pharyngeal swab",
    #         "nasopharyngeal / oropharyngeal swab",
    #         "Nasopharyngeal / Oropharyngeal swab",
    #         "Oro-nasopharyngeal",
    #         "nasopharyngeal and oropharyngeal swab",
    #         "np/op",
    #         "Nasopharyngeal and oropharyngeal swab",
    #         "Nasopharyngeal and Oro-pharyngeal swab",
    #         "Naso-pharyngeal swab/Oro-pharyngeal swab",
    #         "Combined naso- and oropharyngeal swab",
    #         "combined nasopharyngeal and oropharyngeal swab",
    #         "Nasopharyngeal swab/Oropharyngeal swab",
    #         "Naso and/or oro-pharyngeal swab",
    #         "Naso and/or oropharyngeal swab",
    #         "Oro/naso-pharyngeal swab",
    #         "Oro/naso-pharyngeal swabs",
    #         "Oro-naso-pharyngeal swab",
    #         "Oro-naso-pharyngeal swabs",
    #         "NP/OP swab",
    #         "nasopharyngeal/oropharyngeal swab",
    #     ],
    #     "Nasopharyngeal & pharyngeal swab": [
    #         "Naso and pharyngeal swab",
    #         "Pharyngeal and Nasopharyngeal",
    #         "Pharyngeal and Nasopharyngeal swab",
    #     ],
    #     "Nasopharyngeal & throat swab": [
    #         "Nasopharyngeal swab and Throat swab",
    #         "Nasopharyngeal/throat swab",
    #         "Nasopharyngeal swab and Throat swab",
    #         "Nasopharyngeal swab & Throat swab",
    #         "Nasopharyngeal swab/throat swab",
    #         "NPS & TS",
    #         "NPS+TS",
    #     ],
    #     "Nasopharyngeal & tracheal swab": ["Nasopharyngeal & Tracheal swab"],
    #     "Nasopharyngeal aspirate": [
    #         "Naso-pharyngeal aspirate",
    #         "Nasopharyngeal Aspirate",
    #     ],
    #     "Nasopharyngeal aspirate & throat swab": [
    #         "Nasopharyngeal aspirate & Throat swab",
    #         "Nasopharyngeal aspirate/throat swab",
    #         "NPA & TS",
    #     ],
    #     "Nasopharyngeal aspirate & tracheal swab": [
    #         "Nasopharyngeal aspirate & Tracheal swab",
    #         "Nasopharygeal aspirate & Tracheal swab",
    #     ],
    #     "Nasopharyngeal exudate": ["Naso-pharyngeal exudate", "Nasopharingeal exudate"],
    #     "Nasopharyngeal washings": ["Nasopharyngeal (Throat) washings"],
    #     "Oral swab": ["mouth swab", "oral swab"],
    #     "Oro-pharyngeal swab": [
    #         "oropharyngeal swab",
    #         "Oro-Pharyngeal swab",
    #         "Orao-pharungeal swab",
    #         "Oropharyngeal swab",
    #         "Oro-pharyngael swab",
    #         "Oro-pharingeal swab",
    #         "Oro-pharyngel swab",
    #         "Oro-pharyngeal",
    #         "Oro-pharyngeal swab,",
    #         "ora-pharyngeal swab",
    #         "Ora-pharyngeal swab",
    #         "OF swab",
    #         "Orol-pharyngeal swab",
    #         "Throat (Oropharyngeal) swab",
    #         "Oro-pharangeal swab",
    #         "Oro-Pharyngeal Swab",
    #         "Posterior oropharyngeal swab",
    #         "Oro-pharyngeal Swab",
    #         "OPS",
    #         "Oropharynx",
    #     ],
    #     "Oro-pharyngeal & conjunctival swab": ["Oropharyngeal/conjunctival swab"],
    #     "Oronasopharynx": ["Oronasopahrynx", "oronasopharynx", "Oronasopharynx",],
    #     "Pharyngeal swab": ["pharyngeal swab", "pharynx swab"],
    #     "Pharyngeal exudate": ["Pharingeal exudate"],
    #     "Plasma": [],
    #     "Pulmonary swab": ["pulmonary swab"],
    #     "Respiratory swab": ["respiratory swab"],
    #     "Respiratory secretion": [],
    #     "Saliva": [],
    #     "Saliva & nasal swab": ["Saliva and Nasal Swab", "Saliva and Nasal swab"],
    #     "Saliva & nasopharyngeal swab": [
    #         "Saliva and Nasopharyngeal swab",
    #         "Saliva and Nasopharyngeal Swab",
    #     ],
    #     "Scale": [],
    #     "Serum": [],
    #     "Sputum": ["sputum", "Sputum/PBS"],
    #     "Stool": ["Faeces", "Feces"],
    #     "Swab (unspecified)": ["swab", "Swab"],
    #     "Throat saliva": [],
    #     "Throat swab": [
    #         "throat swab",
    #         "Throat Swab",
    #         "Throat",
    #         "Throat swabs",
    #         "Throat swab; passaged twice in Vero E6 cells",
    #     ],
    #     "Throat swab, serum, blood, sputum, urine": [
    #         "Throat Swab, Serum, Blood, Sputum, Urine"
    #     ],
    #     "Throat washing": ["throat washing"],
    #     "Tracheal aspirate": ["Tracheal aspirate sample"],
    #     "Tracheal secretion": ["tracheal secretion"],
    #     "Tracheal swab": ["Thracheal Swab", "tracheal"],
    #     "Upper respiratory tract": [
    #         "Upper respiratory specimen",
    #         "upper respiratory specimen",
    #         "Upper respiratory tract swab",
    #         "upper respiratory tract swab",
    #     ],
    #     "Upper respiratory secretion": ["upper respiratory secretion"],
    #     "Urine": [],
    #     "Urine & blood": ["Blood, urine"],
    #     "Vomit fluid": [],
    #     "Wastewater": [
    #         "Untreated wastewater",
    #         "Wastewater; sewage sample",
    #         "Raw sewage",
    #     ],
    # }

    # specimen_map = {}
    # for k, v in specimen_key_map.items():
    #     # Add self
    #     specimen_map[k] = k
    #     for _v in v:
    #         specimen_map[_v] = k

    # df.loc[:, "specimen"] = df["specimen"].map(specimen_map)

    # df.loc[:, "specimen"] = (
    #     df["specimen"].map(specimen_map).combine_first(df["specimen"]).fillna("Unknown")
    # )

    # print("done")
    # print(
    #     '"Specimen" values: {}'.format(
    #         ", ".join(['"{}"'.format(x) for x in df["specimen"].unique()])
    #     )
    # )

    return df


def clean_date_metadata(df):
    """Clean the collection and submission date metadata"""

    df.loc[:, "collection_date"] = df["covv_collection_date"].astype(str).str.strip()
    df.loc[:, "submission_date"] = df["covv_subm_date"].astype(str).str.strip()

    # Filter out really unspecific collection dates
    # If the date is 4 characters or less (a year, like "2019", or "2020"), then remove it
    df = df.loc[df["collection_date"].str.len() > 4, :]

    # df["collection_date"] = df["collection_date"].fillna(
    #     "Unknown"
    # )

    # Fix weird year formats
    df.loc[:, "collection_date"] = df["collection_date"].str.replace(r"^0020", "2020")
    df.loc[:, "collection_date"] = df["collection_date"].str.replace(r"^0021", "2021")

    # Convert dates to datetime
    df.loc[:, "collection_date"] = pd.to_datetime(df["collection_date"], yearfirst=True)
    df.loc[:, "submission_date"] = pd.to_datetime(df["submission_date"], yearfirst=True)

    return df


def clean_lineage_metadata(df):
    df["lineage"] = df["covv_lineage"].astype(str).str.strip()

    # Filter out "None" lineages
    #remove_seqs = (df["lineage"] == "None") | (df["lineage"] == "nan")
    #df = df.loc[~remove_seqs, :]
    #print("Removed {} sequences without a lineage assignment".format(remove_seqs.sum()))

    return df


def clean_clade_metadata(df):
    df["clade"] = df["covv_clade"].astype(str).str.strip()
    return df


def clean_seq_tech_metadata(df):
    """Clean "Sequencing Technology" values"""

    # print("Cleaning sequencing technology metadata...", end="", flush=True)

    # Basic cleaning
    df["sequencing_tech"] = (
        df["covv_seq_technology"].fillna("Unknown").astype(str).str.strip()
    )

    # replace_map = [
    #     (r"illumina", "Illumina", False),
    #     (r"Ilumina", "Illumina"),
    #     (r"^llumina", "Illumina"),
    #     (r"Illumina\'s", "Illumina"),
    #     (r"Illumina_", "Illumina "),
    #     (r"Illumina technology", "Illumina"),
    #     (r"Immumina", "Illumina"),
    #     (r"iSeq([0-9]+)", lambda m: "iSeq {}".format(m.groups()[0])),
    #     (r"hiseq", "HiSeq", False),
    #     (r"miseq", "MiSeq", False),
    #     (r"MiSeq\.", "MiSeq"),
    #     (r"^MiSeq", "Illumina MiSeq"),
    #     (r"miniseq", "MiniSeq", False),
    #     (r"nextseq", "NextSeq", False),
    #     (r"Next\sSeq", "NextSeq"),
    #     (r"NextSeq([0-9]+)", lambda m: "NextSeq {}".format(m.groups()[0])),
    #     (r"NextSeq 5[01]{1}[0-9]{1}", "NextSeq 500"),
    #     (r"^NextSeq", "Illumina NextSeq"),
    #     (r"novaseq", "NovaSeq", False),
    #     (r"Noveseq", "NovaSeq"),
    #     (r"NovaSeq([0-9+])", lambda m: "NovaSeq {}".format(m.groups()[0])),
    #     (
    #         r"^(NovaSeq|MiSeq|NextSeq|iSeq|HiSeq|MiniSeq)",
    #         lambda m: "Illumina {}".format(m.groups()[0]),
    #     ),
    #     (r"OXFORD_NANOPORE", "Nanopore"),
    #     (r"nanopore", "Nanopore", False),
    #     (r"Nanpore", "Nanopore"),
    #     (r"Nanopre", "Nanopore"),
    #     (r"Nanopore technology", "Nanopore"),
    #     (r"oxford nanopore technology", "Nanopore", False),
    #     (r"Oxoford", "Oxford"),
    #     (r"Oxford Nanopore", "Nanopore"),
    #     (r"minion", "MinION", False),
    #     (r"MinION[,\.]$", "MinION"),
    #     (r"GRIDION", "GridION"),
    #     (r"^(MinION|GridION)", lambda m: "Nanopore {}".format(m.groups()[0])),
    #     (r"IonTorrent", "Ion Torrent", False),
    #     (r"IonTorren", "Ion Torrent"),
    #     (r"ion torrent", "Ion Torrent", False),
    #     (r"Ion Torren", "Ion Torrent"),
    #     (r"Ion Torrentt", "Ion Torrent"),
    #     (r"artic", "ARTIC", False),
    #     (r"ARTIC\sprotocol", "ARTIC sequencing protocol"),
    #     (r"ARTIC\sProtocol", "ARTIC sequencing protocol"),
    #     (r"\-\-\sARTIC", "- ARTIC"),
    #     (r"MGISEQ-?([0-9]+)", lambda m: "MGISEQ {}".format(m.groups()[0])),
    #     (r"Sanger dideoxy sequencing", "Sanger"),
    #     (r"Sanger Sequencing Method", "Sanger"),
    #     (r", assembled sequences", ""),
    #     (r";", ","),
    #     (r"\sand\s", ", "),
    #     (r"\s\&\s", ", "),
    #     (r"+", ", "),
    #     (r",\s+", ", "),
    #     (r"_", " "),
    # ]

    # for pair in replace_map:
    #     df["sequencing_tech"] = df["sequencing_tech"].str.replace(
    #         pair[0], pair[1], -1, pair[2] if len(pair) > 2 else True
    #     )

    # seq_key_map = {
    #     "DNBSEQ-T7": [],
    #     "DNBSEQ-G400": ["MGI Tech. DNBSEQ-G400"],
    #     "DNBSEQ-G400RS": ["DNBSEQ-G400RS, MGI Tech Co., Ltd"],
    #     "Illumina, Nanopore": [],
    #     "Nanopore GridION": ["Nanopore - GridION"],
    #     "Nanopore MinION - ARTIC sequencing protocol": [
    #         "Nanopore MinION ARTIC sequencing protocol"
    #     ],
    #     "Nanopore MinION ARTIC Network V1": ["Nanopore MinION ARTIC v1 primers"],
    #     "Nanopore, Sanger": [],
    #     "Sanger": ["Sanger dideoxy sequencing"],
    #     "Unspecified NGS": ["NGS"],
    #     "Unknown": ["unknown", "Oro-pharyngeal swab", "unkown", "contacted submitter"],
    # }

    # seq_map = {}
    # for k, v in seq_key_map.items():
    #     # Add self
    #     seq_map[k] = k
    #     for _v in v:
    #         seq_map[_v] = k

    # df["sequencing_tech"] = (
    #     (df["sequencing_tech"].map(seq_map))
    #     .combine_first(df["sequencing_tech"])
    #     .fillna("Unknown")
    # )

    # print("done")

    return df


def clean_assembly_metadata(df):
    """Clean "Assembly Method" column"""

    print("Cleaning assembly method metadata...", end="", flush=True)

    df["assembly_method"] = (
        df["covv_assembly_method"].fillna("Unknown").astype(str).str.strip()
    )

    # replace_map = [
    #     # Aliases
    #     (r"artic", "ARTIC", False),
    #     (r"ArticNetwork", "ARTIC Network"),
    #     (r"bcftools", "BCFtools", False),
    #     (r"bowtie", "bowtie", False),
    #     (r"bwa", "BWA", False),
    #     (r"bwa[-|\s]mem", "BWA-MEM", False),
    #     (r"(Burrows-Wheeler Aligner method)", "", False),
    #     (r"Burrows-Wheeler Aligner Tool (BWA)", "BWA", False),
    #     (r"custom", "Custom", False),
    #     (r"dragen", "DRAGEN", False),
    #     (r"geneious prime", "Geneious Prime", False),
    #     (r"megahit", "MEGAHIT", False),
    #     (r"minimap2", "minimap2", False),
    #     (r"mpileup", "mpileup", False),
    #     (r"samtools", "samtools", False),
    #     (r"PAdes", "SPAdes"),
    #     (r"SSPAdes", "SPAdes"),
    #     (r"spades", "SPAdes", False),
    #     (r"vcftools", "vcftools", False),
    #     (r"Workbench ([0-9\.]+)", lambda m: "Workbench v{}".format(m.groups()[0])),
    #     (r"v. ([0-9\.]+)", lambda m: "v{}".format(m.groups()[0])),
    #     ("\ufeff", ""),
    #     (r" assembly method", ""),
    #     # Separators
    #     (r";", ","),
    #     (r"\sand\s", ", "),
    #     (r"\s\&\s", ", "),
    #     (r"+", ", "),
    #     (r",\s+", ", "),
    #     (r"--", "-"),
    #     (r"\s\/\s", " - "),
    #     # These are Excel errors
    #     (r"Sequencher 5.4.[789]", "Sequencher 5.4.6"),
    #     (r"Sequencher 5.4.[123]{1}[0-9]{1}", "Sequencher 5.4.6"),
    #     (r"Sequencher 5.4.4[12]", "Sequencher 5.4.40"),
    #     # More Excel errors
    #     (r"mapped to NC_045512v[3-9]", "mapped to NC_045512v2"),
    #     (r"mapped to NC_045512v[0-9]{1}[0-9]{1}", "mapped to NC_045512v2"),
    #     (r"IRMA v1.0.[1-9]", "IRMA v1.0.0"),
    #     (r"IRMA v1.0.[0-9]{1}[0-9]{1}", "IRMA v1.0.0"),
    #     (r"IME-BJ0[2-9]", "IME-BJ01"),
    #     (r"Geneious Prime 2020.1.[3-9]{1}", "Geneious Prime 2020.1.2"),
    #     (r"Geneious Prime 2020.1.[0-9]{1}[0-9]{1}", "Geneious Prime 2020.1.2"),
    #     (r"Geneious Prime 2020.0.[7-9]", "Geneious Prime 2020.0.6"),
    #     (r"bowtie[3-9]", "bowtie2"),
    #     (r"bowtie[0-9]{1}[0-9]{1}", "bowtie2"),
    #     (r"minimap[3-9]", "minimap2"),
    #     (r"minimap[0-9]{1}[0-9]{1}", "minimap2"),
    # ]

    # for pair in replace_map:
    #     df["assembly_method"] = df["assembly_method"].str.replace(
    #         pair[0], pair[1], -1, pair[2] if len(pair) > 2 else True
    #     )

    # method_key_map = {"Unknown": ["unknown"]}

    # method_map = {}
    # for k, v in method_key_map.items():
    #     # Add self
    #     method_map[k] = k
    #     for _v in v:
    #         method_map[_v] = k

    # df["assembly_method"] = (
    #     (df["assembly_method"].map(method_map))
    #     .combine_first(df["assembly_method"])
    #     .fillna("Unknown")
    # )

    # print("done")

    return df


def clean_comment_type_metadata(df):
    df["comment_type"] = df["covv_comment_type"].astype(str).str.strip()
    df["comment_type"] = df["comment_type"].fillna("None")
    return df


def clean_author_metadata(df):
    df["authors"] = df["covv_authors"].astype(str).str.strip()
    return df


def clean_originating_lab_metadata(df):
    df["originating_lab"] = df["covv_orig_lab"].astype(str).str.strip()
    return df


def clean_submitting_lab_metadata(df):
    df["submitting_lab"] = df["covv_subm_lab"].astype(str).str.strip()
    return df


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--metadata-in",
        type=str,
        required=True,
        help="Path to input metadata CSV file",
    )

    parser.add_argument(
        "-l",
        "--location-corrections",
        type=str,
        required=True,
        help="Path to location corrections CSV file",
    )

    parser.add_argument(
        "--lineages",
        type=str,
        required=True,
        help="Path to lineages CSV file",
    )
    parser.add_argument(
        "--quality",
        type=str,
        required=True,
        help="Path to quality CSV file",
    )

    parser.add_argument(
        "-o",
        "--metadata-out",
        type=str,
        required=True,
        help="Path to output metadata CSV file",
    )

    args = parser.parse_args()

    # Load metadata
    df = pd.read_csv(args.metadata_in, low_memory=False)
    df = df.rename(columns={"covv_accession_id": "Accession ID"}).set_index(
        "Accession ID"
    )

    df = clean_name_metadata(df)
    df = clean_host_metadata(df)
    df = clean_gender_metadata(df)
    df = clean_age_metadata(df)
    df = clean_patient_status_metadata(df)
    df = clean_passage_metadata(df)
    df = clean_specimen_metadata(df)
    df = clean_date_metadata(df)
    df = clean_lineage_metadata(df)
    df = clean_clade_metadata(df)
    df = clean_seq_tech_metadata(df)
    df = clean_assembly_metadata(df)
    df = clean_comment_type_metadata(df)
    df = clean_author_metadata(df)
    df = clean_originating_lab_metadata(df)
    df = clean_submitting_lab_metadata(df)

    df = process_location_metadata(df, args.location_corrections)

    # Take subset of columns
    df = df[
        [
            "virus_name",
            "collection_date",
            "submission_date",
            "host",
            "gender",
            "age_start",
            "age_end",
            "patient_status",
            "passage",
            "specimen",
            "lineage",
            "clade",
            "sequencing_tech",
            "assembly_method",
            "comment_type",
            "authors",
            "originating_lab",
            "submitting_lab",
            # Location data
            "region",
            "country",
            "division",
            "location",
        ]
    ]

    # Isolate ID = same as Accession ID
    df["isolate_id"] = df.index.values
    # Segment = 1
    df["segment"] = 1

    # Load lineages and join to dataframe
    lineages_df = pd.read_csv(args.lineages)
    lineages_df = lineages_df.rename(columns={"taxon": "Accession ID"}).set_index(
        "Accession ID"
    )
    df = df.rename(columns={"lineage": "gisaid_lineage"}).join(
        lineages_df[
            [
                "lineage",
                "conflict",
                "ambiguity_score",
                "scorpio_call",
                "scorpio_support",
                "scorpio_conflict",
                "scorpio_notes",
                # "version",
                # "pangolin_version",
                # "scorpio_version",
                # "constellation_version",
                "is_designated",
                "qc_status",
                "qc_notes",
                "note",
            ]
        ].rename(
            columns={
                "is_designated": "pangolin_is_designated",
                "qc_status": "pangolin_qc_status",
                "qc_notes": "pangolin_qc_notes",
                "note": "pangolin_note",
            }
        ),
        how="left",
    )
    # Fill in missing values with GISAID lineages
    df.loc[:, "lineage"] = df["lineage"].combine_first(df["gisaid_lineage"])

    # Load quality and join to dataframe
    quality_df = pd.read_csv(args.quality, index_col="Accession ID")
    df = df.join(quality_df, how="left")
    # Calculate percent ambiguous, drop the num_ambiguous column
    df["num_ambiguous"] = (
        ((df["num_ambiguous"] / df["length"]) * 100).fillna(0).astype(int)
    )
    df.rename(columns={"num_ambiguous": "percent_ambiguous"}, inplace=True)

    df.to_csv(args.metadata_out)


if __name__ == "__main__":
    main()
