# coding: utf-8

"""Clean patient metadata

Author: Albert Chen (Deverman Lab - Broad Institute)
"""

import numpy as np
import pandas as pd
import re


def clean_gender_metadata(patient_meta_df):
    """Clean patient gender metadata
    """

    print("Cleaning patient gender metadata...", end="", flush=True)

    # Make a copy, strip whitespace
    patient_meta_df["gender"] = patient_meta_df["Gender"].str.strip()

    replace_map = [
        (r"^female", "Female", False),
        (r"^male", "Male", False),
        (r"^f$", "Female", False),
        (r"^m$", "Male", False),
        (r"unknown", "Unknown", False),
    ]

    for pair in replace_map:
        patient_meta_df["gender"] = patient_meta_df["gender"].str.replace(
            pair[0], pair[1], -1, pair[2] if len(pair) > 2 else True
        )

    gender_key_map = {
        "Female": ["Woman", "Femal"],
        "Male": ["M"],
        "Unknown": ["unknwon", "Unkown", "U"],
    }

    gender_map = {}
    for k, v in gender_key_map.items():
        # Add self
        gender_map[k] = k
        for _v in v:
            gender_map[_v] = k

    patient_meta_df["gender"] = (
        patient_meta_df["gender"]
        .map(gender_map)
        .combine_first(patient_meta_df["gender"])
        .fillna("Unknown")
    )

    print("done")
    print('"Gender" values: {}'.format(", ".join(patient_meta_df["gender"].unique())))

    patient_meta_df["gender"] = patient_meta_df["gender"].fillna("Unknown")

    return patient_meta_df


def clean_age_metadata(patient_meta_df):
    """Clean patient age metadata

    For each age value we want to define an age range that the value
    corresponds to. This is necessary since the metadata is provided in
    various different specificities. 

    i.e., exact age (72.343), year (42), or age range (20-30) 

    Define a range [start, end), for each age
    This ranges can then be filtered over in a way that includes as much data
    as possible

    """

    print("Cleaning patient age metadata...", end="", flush=True)

    # Do some basic cleanup before we start
    patient_meta_df["age_clean"] = patient_meta_df["Patient age"]
    patient_meta_df["age_clean"] = patient_meta_df["age_clean"].fillna("Unknown")
    patient_meta_df["age_clean"] = patient_meta_df["age_clean"].str.strip()

    patient_meta_df["age_start"] = np.nan
    patient_meta_df["age_end"] = np.nan

    for i, v in patient_meta_df["age_clean"].iteritems():

        # Skip if Unknown
        if v == "Unknown":
            continue

        # Merge "Unknowns"
        # "Adult" is too vague
        # "Over 18" is too vague
        elif v in [
            "unknown",
            "unkown",
            "unknwon",
            "unavailable",
            "uknown",
            "no data",
            "Male",
            "Adult",
            "Over 18",
            "over 18",
        ]:
            # patient_meta_df.loc[i, 'Patient age clean'] = 'Unknown'
            continue

        # Parse clean integers. e.g., "42"
        elif re.match(r"^[0-9]+$", v):
            patient_meta_df.loc[i, "age_start"] = float(int(v))
            patient_meta_df.loc[i, "age_end"] = float(int(v) + 1)

        # Parse fractions. e.g., "72.323"
        elif re.match(r"^[0-9]*\.[0-9]+$", v):
            patient_meta_df.loc[i, "age_start"] = float(v)
            patient_meta_df.loc[i, "age_end"] = float(v)

        # Parse N months (less than 1 years old). e.g., "7 months"
        elif re.match(r"^([1]*[0-9])+\smonth(s)?$", v):
            # Re-run regex to extract months
            m = re.match(r"^([1]*[0-9])+\smonth(s)?$", v)
            month = int(m.groups()[0])
            # Convert months to fraction of a year, then round to an integer
            patient_meta_df.loc[i, "age_start"] = month / 12.0
            patient_meta_df.loc[i, "age_end"] = (month + 1) / 12.0

        # Same, but days. e.g., "17 days"
        elif re.match(r"^([1]*[0-9])+\sday(s)?$", v):
            # Re-run regex to extract days
            m = re.match(r"^([1]*[0-9])+\sday(s)?$", v)
            day = int(m.groups()[0])
            # Convert days to fraction of a year, then round to an integer
            patient_meta_df.loc[i, "age_start"] = day / 365.0
            patient_meta_df.loc[i, "age_end"] = (day + 1) / 365.0

        # Same, but weeks. e.g., "6 weeks"
        elif re.match(r"^([1]*[0-9])+\sweek(s)?$", v):
            # Re-run regex to extract weeks
            m = re.match(r"^([1]*[0-9])+\sweek(s)?$", v)
            week = int(m.groups()[0])
            # Convert weeks to fraction of a year, then round to an integer
            patient_meta_df.loc[i, "age_start"] = week / 52.0
            patient_meta_df.loc[i, "age_end"] = (week + 1) / 52.0

        # Remove '-year old', 'years' or 'age' at end. e.g., "44-year old"
        # None of the years/age entries are fractions, so treat as a valid integer
        # This might easily break in the future
        elif re.match(r"^([0-9]+)\s?(-year\sold|years|age)$", v):
            # Re-run regex to extract age
            m = re.match(r"^([0-9]+)\s?(-year\sold|years|age)$", v)
            year = int(m.groups()[0])
            patient_meta_df.loc[i, "age_start"] = float(year)
            patient_meta_df.loc[i, "age_end"] = float(year + 1)

        # Match year-month, e.g., "6 years 2 months"
        elif re.match(r"^([0-9]+)(?:,\s|\syears\s)([1]?[0-9]+)\smonths$", v):
            # Re-run to extract years and months
            m = re.match(r"^([0-9]+)(?:,\s|\syears\s)([1]?[0-9]+)\smonths$", v)
            # Extract years, months
            years = int(m.groups()[0])
            months = int(m.groups()[1])
            # Round to nearest year
            patient_meta_df.loc[i, "age_start"] = years + (months / 12.0)
            patient_meta_df.loc[i, "age_end"] = years + ((months + 1) / 12.0)

        # Extract decade ranges. e.g., "30s"
        elif re.match(r"^([0-9]+)\'?s$", v):
            # Re-run to extract decade
            m = re.match(r"^([0-9]+)\'?s$", v)
            decade = int(m.groups()[0])
            patient_meta_df.loc[i, "age_start"] = float(decade)
            patient_meta_df.loc[i, "age_end"] = float(decade + 10)

        # Extract year ranges, e.g., "10-20"
        elif re.match(r"^([0-9]+)\s?-\s?([0-9])+$", v):
            # Re-run to extract range
            m = re.match(r"^([0-9]+)\s?-\s?([0-9])+$", v)
            start = int(m.groups()[0])
            end = int(m.groups()[1])

            patient_meta_df.loc[i, "age_start"] = float(start)
            patient_meta_df.loc[i, "age_end"] = float(
                end + 1
            )  # Assume that the provided range is [start, end]

        # Extract inequalities, e.g., ">60"
        elif re.match(r"^>([0-9]+)$", v):
            # Re-run to extract year
            m = re.match(r"^>([0-9]+)$", v)
            start = int(m.groups()[0])

            # Assume that the range is just one year
            patient_meta_df.loc[i, "age_start"] = float(start)
            patient_meta_df.loc[i, "age_end"] = float(start + 1)

    print("done")
    print(
        'Could not parse the following "Patient Age" values: {}'.format(
            ", ".join(
                patient_meta_df["age_clean"][pd.isnull(patient_meta_df["age_start"])]
                .str.strip()
                .unique()
                .astype(str)
            )
        )
    )

    return patient_meta_df


def clean_patient_status_metadata(patient_meta_df):

    print("Cleaning patient status metadata...", end="", flush=True)

    # Strip whitespace
    patient_meta_df["patient_status"] = patient_meta_df["Patient status"].str.strip()

    replace_map = [
        (r"hospitalized", "Hospitalized", False),
        (r"^fever", "Fever", False),
        (r"^live", "Live", False),
        (r"outpatient", "Outpatient", False),
        (r"released", "Released", False),
        (r"asymptomatic", "Asymptomatic", False),
        (r"unknown", "Unknown", False),
    ]

    for pair in replace_map:
        patient_meta_df["patient_status"] = patient_meta_df[
            "patient_status"
        ].str.replace(pair[0], pair[1], -1, pair[2] if len(pair) > 2 else True)

    status_key_map = {
        "Unknown": [
            "unknwon",
            "Unkown",
            "unkown",
            "unknow",
            "Oro-pharyngeal swab",
            "Z039",
            "-",
            "Unknow",
            "uknown",
            "\ufeffUnknown",
        ],
        "Alive": ["Live", "Live, physical examination"],
        "Hospitalized": [
            "Hospitalized; Stable",
            "Hospitalized or to be Hospitalized",
            "Hospitalized, Live",
            "In-hospital",
            "Hospitaized",
            "inpatient",
            "Hospitalized patient",
        ],
        "Recovered": ["Cured"],
        "Discharged": [
            "Released",
            "Recovered and Released",
            "Initially Hospitalized, but now improved and discharged",
            "Outpatient",
            "Hospitalized/Released",
            "Released, Live",
            "Moderate/Outpatient",
            "Discharged after recovery",
        ],
        "Pneumonia": [
            "Pneumonia (chest X-ray), not critical",
            "Pneumonia (chest X-ray)",
        ],
        "Asymptomatic": [
            "Asymptomatic/Released",
            "Mild/Contact exposure/Asymptomatic",
            "Healthy",
        ],
        "Symptomatic": [
            "Mild symptoms (fever, cardiovascular disorders)",
            "Mild symptoms inpatient for observation",
            "Mild clinical signs without hospitalization",
            "Mild",
            "fever",
            "Live, acute respiratory infection",
        ],
        "Intensive Care Unit": [
            "Intensive Care Unit",
            "Hospitalized in ICU",
            "ICU; Serious",
            "Severe/ICU",
        ],
        "Mild Case": ["Mild case"],
        "Deceased": ["Death", "Hospitalized/Deceased"],
        "LTC": ["EHPAD", "EHPAD_IRA"],
        "Quarantine": ["Stable in quarantine", "Quarantined", "Isolation"],
        "Home": ["Mild, at home.", "Live, mild symptoms, at home"],
    }

    status_map = {}
    for k, v in status_key_map.items():
        # Add self
        status_map[k] = k
        for _v in v:
            status_map[_v] = k

    patient_meta_df["patient_status"] = (
        patient_meta_df["patient_status"]
        .map(status_map)
        .combine_first(patient_meta_df["patient_status"])
        .fillna("Unknown")
    )

    print("done")
    print(
        '"Patient Status" values: {}'.format(
            ", ".join(patient_meta_df["patient_status"].unique())
        )
    )

    return patient_meta_df


def clean_passage_metadata(patient_meta_df):
    """Clean cell passage metadata
    """

    print("Cleaning cell passage metadata...", end="", flush=True)

    # Basic cleaning
    patient_meta_df["passage"] = patient_meta_df["Passage"].str.strip()

    passage_key_map = {
        "Original": [
            "Orginal",
            "Orignal",
            "Origional",
            "Original (first sample)",
            "original",
            "Clinical specimen",
            "Clinicial specimen",
            "Original (Clinical sample)",
            "Orignial",
            "Original (Clinical Sample)",
            "Liver",
            "Original (Brain tissue)",
            "Oroginal",
            "Original(Clinical Sample)",
            "Orginial",
            "Original, from nasopharyngeal aspirate",
        ],
        "Vero": [],
        "Vero P1": [
            "Vero p1",
            "Vero CCL81 isolate P1",
            "Vero cell P1",
            "Vero C1",
            "Vero 1",
            "Vero1",
        ],
        "Vero P2": ["Vero cell P2", "P2, Vero"],
        "Vero P3": ["Vero cell P3"],
        "Vero P4": ["Vero cell P4"],
        "LLC-MK2": [],
        "P1": ["Virus Isolate, Passage 1", "C1"],
        "P2": [],
        "P3": [],
        "Vero E6": [
            "VeroE6",
            "Vero-E6",
            "VeroE6 cells",
            "Initial isolation on VeroE6 cells",
        ],
        "Vero E6 P1": [
            "VeroE6/P1",
            "Vero E6, 1st passage",
            "P1/Vero E6",
            "Vero E6 #1",
            "Pi/Vero E6",
            "Vero E6, 1 passage",
            "Vero E6 P1, mouse-P14",
        ],
        "Vero E6 P2": ["VeroE6/P2", "VERO E6 / P2", "Vero E6 cells, P2"],
        "Vero E6 P3": ["VeroE6/P3", "Vero E6-P3"],
        "Vero E6 P4": ["Passage 4 in Vero E6 cells"],
        "Vero/hSLAM P1": [],
        "Vero E6/TMPRSS2": ["VeroE6/TMPRSS2"],
        "Caco-2": ["C2"],
        "Caco-2 P1": [],
        "Culture (unknown)": ["Virus culture"],
    }

    passage_map = {}
    for k, v in passage_key_map.items():
        # Add self
        passage_map[k] = k
        for _v in v:
            passage_map[_v] = k

    patient_meta_df["passage"] = patient_meta_df["passage"].map(passage_map)

    print("done")
    print(
        'Setting "Passage" values to "Unknown": {}'.format(
            ", ".join(
                patient_meta_df["Passage"][pd.isnull(patient_meta_df["passage"])]
                .str.strip()
                .unique()
                .astype(str)
            )
        )
    )

    patient_meta_df["passage"] = patient_meta_df["passage"].fillna("Unknown")

    return patient_meta_df


def clean_specimen_metadata(patient_meta_df):

    print("Cleaning specimen metadata...", end="", flush=True)

    # Basic cleanup
    patient_meta_df["specimen"] = patient_meta_df["Specimen"].str.strip()

    specimen_key_map = {
        "Alveolar lavage fluid": [],
        "Anal swab": ["Anus", "Anal swabs"],
        "Aspirate": [],
        "Autopsy, Bronch and lung": [],
        "Blood": ["blood sample", "Blood sample"],
        "Breathing air using VIVAs Air sampler": [],
        "Bronchial scraper": [],
        "Bronchoalveolar lavage fluid": [
            "BALF sample of 2019 Wuhan peumonia patient 01",
            "BALF sample of 2019 Wuhan peumonia patient 02",
            "bronchoalveolar lavage fluid",
            "Bronchoalveolar Lavage NCIT:C51913",
            "Bronchoalveolar lavage",
            "Bronchoalveolar fluid",
            "Broncho-alveolar lavage",
            "Bronchoalveolar lavage fluid (BALF)",
            "Bronchoalveolar-lavage fluid",
            "Broncho-alveolar fluid",
            "Broncoalveolar",
        ],
        "Buccal swab": ["buccal swab"],
        "Deep throat saliva": [],
        "Dry swab": [],
        "EDTA-Plasma": [],
        "Endotracheal aspirate": [
            "endotracheal aspirates",
            "Endotracheal aspirate (ETA)",
            "endotracheal aspirate (ETA)",
        ],
        "Enviromental swab": ["Door handle", "Air"],
        "Fecal swab": [],
        "Intestine tissue": [],
        "Lung tissue": ["Lung", "lung tissue", "lung biopsy"],
        "Lung & intenstine tissue": ["Mixed lung and intestine tissue"],
        "Lung secretion": ["Lung secretion"],
        "Lung swab": ["lung swab"],
        "Mid-turbinate nasal swab": ["Mid-Turbinate nasal swab"],
        "Nasal swab": [
            "Nose",
            "Nose swab",
            "nasal swab",
            "Nasal Swab",
            "Mid-nasal swab",
            "nose swab",
            "Nasal swab specimen",
        ],
        "Nasal swab & blood": ["Nasal Swab, Blood"],
        "Nasal swab & serum": ["Nasal Swab, Serum", "Nasal Swab, Serum,"],
        "Nasal swab, serum, blood": [
            "Nasal swab, Serum, Blood",
            "Nasal Swab, Blood, Serum",
            "Nasal Swab, Serum, Blood",
            "Nasal Swab, Serum, Blood,",
        ],
        "Nasal swab, serum, blood, plasma": [
            "Nasal Swab, Blood, Serum, Plasma",
            "Nasal Swab, Blood, Serum, Plasma",
            "Nasal Swab, Serum, Blood, Plasma",
        ],
        "Nasal swab, serum, blood, sputum": ["Nasal Swab, Serum, Blood, Sputum"],
        "Nasal swab, serum, blood, sputum, urine": [
            "Nasal Swab, Serum, Blood, Sputum, Urine"
        ],
        "Nasal & throat swab": [
            "Nose/throat swab",
            "Nost and Throad swab",
            "Nose and Throat swab",
            "Nose-Throat",
            "Throat swab Nasal swab",
            "Mixture of Throat swab Nasal swab",
            "Mixture of Throat swab and Nasal swab",
            "Mixture of throat and nasal swab",
            "Nasal swab and Throat swab",
            "Nasal and throat swab",
            "Throat and Nasal Swab",
        ],
        "Nasal & oro-pharyngeal swab": [
            "Nasal, oropharyngeal swab",
            "Oro-pharyngeal swab, Nasal swab",
            "oronasopharynx",
        ],
        "Nasal swab, oral swab, tracheal wash": [
            "Oral swab; Nasal swab; Tracheal wash"
        ],
        "Nasopharyngeal swab": [
            "nasopharyngeal swab",
            "Nasopharyngeal swab",
            "Naso-pharyngeal swab",
            "naso-pharyngeal swab",
            "Nasopharyngeal swap",
            "Nasopharyngeal Swab",
            "Nasopharynx swab",
            "Nasopharingeal swab",
            "NP Swab",
            "Nasopharnygeal swab",
            "NP swab",
            "nasopharyngial swab",
            "NP",
            "Naso-pharingeal swab",
            "Naso-pharygeal swab",
            "Nasopharynx",
            "nasopharingeal swab",
            "Rhino-pharyngeal swab",
            "Upper Resp NPS",
        ],
        "Nasopharyngeal VTM": [],
        "Nasopharyngeal & oro-pharyngeal swab": [
            "Nasopharyngeal swab and oropharyngeal swab",
            "Nasopharyngeal swab, oropharyngeal swab",
            "Combined nasopharyngeal and oropharyngeal swab",
            "Nasopharyngeal/oropharyngeal swab",
            "Naso-Oropharyngeal swab",
            "Swab nasopharyngeal oropharyngeal",
            "Naso and oro-pharyngeal swabs",
            "Oro and naso-pharyngeal swabs",
            "Oro and naso-pharyngeal swab",
            "Nasopharyngeal / oropharyngeal swab",
            "Nasopharyngeal and Oropharyngeal swab",
            "Nasopharyngeal / Oro-pharyngeal swab",
            "nasopharyngeal / oropharyngeal swab",
            "Nasopharyngeal / Oropharyngeal swab",
            "Oro-nasopharyngeal",
            "nasopharyngeal and oropharyngeal swab",
            "np/op",
            "Nasopharyngeal and oropharyngeal swab",
            "Nasopharyngeal and Oro-pharyngeal swab",
            "Naso-pharyngeal swab/Oro-pharyngeal swab",
            "Combined naso- and oropharyngeal swab",
            "combined nasopharyngeal and oropharyngeal swab",
            "Nasopharyngeal swab/Oropharyngeal swab",
            "Naso and/or oro-pharyngeal swab",
            "Naso and/or oropharyngeal swab",
            "Oronasopharynx",
        ],
        "Nasopharyngeal & pharyngeal swab": [
            "Naso and pharyngeal swab",
            "Pharyngeal and Nasopharyngeal",
            "Pharyngeal and Nasopharyngeal swab",
        ],
        "Nasopharyngeal & throat swab": [
            "Nasopharyngeal swab and Throat swab",
            "Nasopharyngeal/throat swab",
            "Nasopharyngeal swab and Throat swab",
            "Nasopharyngeal swab & Throat swab",
            "Nasopharyngeal swab/throat swab",
            "NPS & TS",
            "NPS+TS",
        ],
        "Nasopharyngeal & tracheal swab": ["Nasopharyngeal & Tracheal swab"],
        "Nasopharyngeal aspirate": [
            "Naso-pharyngeal aspirate",
            "Nasopharyngeal Aspirate",
        ],
        "Nasopharyngeal aspirate & throat swab": [
            "Nasopharyngeal aspirate & Throat swab",
            "Nasopharyngeal aspirate/throat swab",
            "NPA & TS",
        ],
        "Nasopharyngeal aspirate & tracheal swab": [
            "Nasopharyngeal aspirate & Tracheal swab",
            "Nasopharygeal aspirate & Tracheal swab",
        ],
        "Nasopharyngeal exudate": ["Naso-pharyngeal exudate", "Nasopharingeal exudate"],
        "Nasopharyngeal washings": ["Nasopharyngeal (Throat) washings"],
        "Oral swab": ["mouth swab", "oral swab"],
        "Oro-pharyngeal swab": [
            "oropharyngeal swab",
            "Oro-Pharyngeal swab",
            "Orao-pharungeal swab",
            "Oropharyngeal swab",
            "Oro-pharyngael swab",
            "Oro-pharingeal swab",
            "Oro-pharyngel swab",
            "Oro-pharyngeal",
            "Oro-pharyngeal swab,",
            "ora-pharyngeal swab",
            "Ora-pharyngeal swab",
            "OF swab",
            "Orol-pharyngeal swab",
            "Throat (Oropharyngeal) swab",
            "Oro-pharangeal swab",
            "Oro-Pharyngeal Swab",
            "Posterior oropharyngeal swab",
            "Oro-pharyngeal Swab",
            "OPS",
        ],
        "Oro-pharyngeal & conjunctival swab": ["Oropharyngeal/conjunctival swab"],
        "Pharyngeal swab": ["pharyngeal swab", "pharynx swab"],
        "Pharyngeal exudate": ["Pharingeal exudate"],
        "Pulmonary swab": ["pulmonary swab"],
        "Respiratory swab": ["respiratory swab"],
        "Respiratory secretion": [],
        "Saliva": [],
        "Saliva & nasal swab": ["Saliva and Nasal Swab", "Saliva and Nasal swab"],
        "Saliva & nasopharyngeal swab": [
            "Saliva and Nasopharyngeal swab",
            "Saliva and Nasopharyngeal Swab",
        ],
        "Scale": [],
        "Serum": [],
        "Sputum": ["sputum", "Sputum/PBS"],
        "Stool": ["Faeces", "Feces"],
        "Swab (unspecified)": ["swab", "Swab"],
        "Throat saliva": [],
        "Throat swab": ["throat swab", "Throat Swab", "Throat", "Throat swabs"],
        "Throat swab, serum, blood, sputum, urine": [
            "Throat Swab, Serum, Blood, Sputum, Urine"
        ],
        "Throat washing": ["throat washing"],
        "Tracheal aspirate": ["Tracheal aspirate sample"],
        "Tracheal secretion": ["tracheal secretion"],
        "Tracheal swab": ["Thracheal Swab", "tracheal"],
        "Upper respiratory tract": [
            "Upper respiratory specimen",
            "upper respiratory specimen",
            "Upper respiratory tract swab",
            "upper respiratory tract swab",
        ],
        "Upper respiratory secretion": ["upper respiratory secretion"],
        "Urine": [],
        "Urine & blood": ["Blood, urine"],
        "Wastewater": [
            "Untreated wastewater",
            "Wastewater; sewage sample",
            "Raw sewage",
        ],
    }

    specimen_map = {}
    for k, v in specimen_key_map.items():
        # Add self
        specimen_map[k] = k
        for _v in v:
            specimen_map[_v] = k

    patient_meta_df["specimen"] = patient_meta_df["specimen"].map(specimen_map)

    print("done")
    print(
        'Setting "Specimen" values to "Unknown": {}'.format(
            ", ".join(
                patient_meta_df["Specimen"][pd.isnull(patient_meta_df["specimen"])]
                .str.strip()
                .unique()
                .astype(str)
            )
        )
    )

    patient_meta_df["specimen"] = patient_meta_df["specimen"].fillna("Unknown")

    return patient_meta_df


def clean_collection_date_metadata(patient_meta_df):

    patient_meta_df["collection_date"] = patient_meta_df["Collection date"].str.strip()

    # patient_meta_df["collection_date"] = patient_meta_df["collection_date"].fillna(
    #     "Unknown"
    # )

    # Convert Collection date to datetime
    patient_meta_df["collection_date"] = pd.to_datetime(
        patient_meta_df["collection_date"], yearfirst=True
    )

    return patient_meta_df


def clean_lineage_metadata(patient_meta_df):
    patient_meta_df["lineage"] = patient_meta_df["Lineage"].str.strip()
    return patient_meta_df


def clean_clade_metadata(patient_meta_df):
    patient_meta_df["clade"] = patient_meta_df["Clade"].str.strip()
    return patient_meta_df


def clean_patient_metadata(patient_meta_df):
    patient_meta_df = clean_gender_metadata(patient_meta_df)
    patient_meta_df = clean_age_metadata(patient_meta_df)
    patient_meta_df = clean_patient_status_metadata(patient_meta_df)
    patient_meta_df = clean_passage_metadata(patient_meta_df)
    patient_meta_df = clean_specimen_metadata(patient_meta_df)
    patient_meta_df = clean_collection_date_metadata(patient_meta_df)
    patient_meta_df = clean_lineage_metadata(patient_meta_df)
    patient_meta_df = clean_clade_metadata(patient_meta_df)

    # Take subset of columns
    patient_meta_df = patient_meta_df[
        [
            "collection_date",
            "Location",
            "gender",
            "age_start",
            "age_end",
            "patient_status",
            "passage",
            "specimen",
            "lineage",
            "clade",
        ]
    ]

    return patient_meta_df
