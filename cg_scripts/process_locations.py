#!/usr/bin/env python3
# coding: utf-8

"""Clean location metadata
Build hierarchical location tree for the UI

Author: Albert Chen (Deverman Lab - Broad Institute)
"""

import json
import numpy as np
import pandas as pd

from functools import reduce
from cg_scripts.util import human_format


def process_location_metadata(case_df):

    # Unset index for faster filtering
    case_df = case_df.reset_index()

    print("Processing location data...", end="", flush=True)
    # Location data is stored in one column, "region / country / division / location"
    location_df = (
        case_df["Location"]
        .str.split("/", expand=True)
        .iloc[:, :4]  # Only take 4 columns
        # Rename columns
        .rename(columns={0: "region", 1: "country", 2: "division", 3: "location"})
        .applymap(lambda x: x.strip() if x else x)
        # Placeholder for missing values, so that it will still
        # be caught by groupby() later on
        .fillna(-1)
    )

    # Clean location data
    location_df = clean_location_data(location_df)

    # Create complete location column from the separate parts
    # This time with no padding
    location_df["loc_str"] = location_df["region"].str.cat(
        [
            location_df["country"].astype(str),
            location_df["division"].astype(str),
            location_df["location"].astype(str),
        ],
        sep="/",
    )

    location_map_df = (
        location_df["loc_str"]
        .drop_duplicates()
        .sort_values(ignore_index=True)
        .reset_index()  # Produce an index column
    )
    # Location data is stored in one column, "region/country/division/location"
    location_map_df = pd.concat(
        [
            location_map_df,
            (
                location_map_df["loc_str"]
                .str.split("/", expand=True)
                .iloc[:, :4]  # Only take 4 columns
                # Rename columns
                .rename(
                    columns={0: "region", 1: "country", 2: "division", 3: "location"}
                )
            ),
        ],
        axis=1,
    )

    # Map location IDs back to taxon_locations dataframe
    location_df["location_id"] = location_df["loc_str"].map(
        pd.Series(
            location_map_df["index"].values, index=location_map_df["loc_str"].values,
        )
    )

    # Take subset of columns, re-index
    location_df = location_df[
        ["location_id", "region", "country", "division", "location"]
    ].reset_index(drop=True)
    print("done")

    # Re-apply Accession ID index
    location_df["Accession ID"] = case_df["Accession ID"]
    location_df = location_df.set_index("Accession ID")

    return location_df, location_map_df


def clean_location_data(location_df):
    """Fix typos, unify nomenclature in location data
    """

    # Make rules, which are mappings of input region/country/division/location ->
    # output region/country/division/location
    # Each rule is a tuple, where the first element is the input rule
    # (which entries to match)
    # and the second element is the output rule
    # (what fields to change)

    rules = [
        # DRC
        # ---
        # Abbreviate name to save space
        ({"country": "Democratic Republic of the Congo"}, {"country": "DRC"}),
        # Fix typos
        (
            {"country": "DRC", "division": "Kongo Central"},
            {"division": "Kongo-Central"},
        ),
        # MOROCCO
        # -------
        # Move to Africa
        ({"country": "Morocco"}, {"region": "Africa"}),
        # Move Cadiz to Andalusia, Spain
        (
            {"country": "Morocco", "location": "Cadiz_h"},
            {
                "region": "Europe",
                "country": "Spain",
                "division": "Andalusia",
                "location": "Cadiz",
            },
        ),
        # Fix typos
        ({"country": "Marocco"}, {"country": "Morocco"}),
        # SENEGAL
        # -------
        # Fix typos
        (
            {"country": "Senegal", "division": ["St Louis", "St Louis (RTL)"]},
            {"division": "St-Louis"},
        ),
        ({"country": "Senegal", "division": ["Dakarhuman"]}, {"division": "Dakar"}),
        # SOUTH AFRICA
        # ------------
        # Unabbreviate province names
        ({"country": "South Africa", "division": "EC"}, {"division": "Eastern Cape"}),
        ({"country": "South Africa", "division": "KZN"}, {"division": "KwaZulu-Natal"}),
        ({"country": "South Africa", "division": "GP"}, {"division": "Gauteng"}),
        ({"country": "South Africa", "division": "LP"}, {"division": "Limpopo"}),
        ({"country": "South Africa", "division": "MP"}, {"division": "Mpumalanga"}),
        ({"country": "South Africa", "division": "NW"}, {"division": "North West"}),
        ({"country": "South Africa", "division": "FS"}, {"division": "Free State"}),
        (
            {"country": "South Africa", "division": ["WC", "Western Cape Province"]},
            {"division": "Western Cape"},
        ),
        # Remove Unknown division
        ({"country": "South Africa", "division": "Unknown"}, {"division": -1}),
        # ASIA
        # ----
        # Fix typos
        ({"region": "Asian"}, {"region": "Asia"}),
        # BANGLADESH
        # ----------
        ({"country": "Bangladesh", "division": "Ragshahi"}, {"division": "Rajshahi"}),
        # CHINA
        # -----
        # Move NanChang into Jiangxi province
        (
            {"country": "China", "division": "NanChang"},
            {"division": "Jiangxi", "location": "Nanchang"},
        ),
        # Move Guangzhou into Guangdong province
        (
            {"country": "China", "division": "Guangzhou"},
            {"division": "Guangdong", "location": "Guangzhou"},
        ),
        # Move Meizhou into Guangdong province
        (
            {"country": "China", "division": "Meizhou"},
            {"division": "Guangdong", "location": "Meizhou"},
        ),
        # Move Hangzhou into Zhejiang province
        (
            {"country": "China", "division": "Hangzhou"},
            {"division": "Zhejiang", "location": "Hangzhou"},
        ),
        # Move Wuhan into Hubei province
        (
            {"country": "China", "division": "Wuhan"},
            {"division": "Hubei", "location": "Wuhan"},
        ),
        # Move Keqiao into Zhejiang
        # Move Shangyu into Zhejiang
        # Move Yuecheng into Zhejiang
        (
            {"country": "China", "division": "Keqiao"},
            {"division": "Zhejiang", "location": "Keqiao"},
        ),
        (
            {"country": "China", "division": "Shangyu"},
            {"division": "Zhejiang", "location": "Shangyu"},
        ),
        (
            {"country": "China", "division": "Yuecheng"},
            {"division": "Zhejiang", "location": "Yuecheng"},
        ),
        # Move Harbin into Heilongjiang
        (
            {"country": "China", "division": "Harbin"},
            {"division": "Heilongjiang", "location": "Harbin"},
        ),
        # Move Shenzhen into Guangdong
        (
            {"country": "China", "division": "Shenzhen"},
            {"division": "Guangdong", "location": "Shenzhen"},
        ),
        # Fix Typos
        ({"country": "China", "division": "Chongqinq"}, {"division": "Chongqing"}),
        # Move HuaShang into Guangdong (province of submitting lab),
        # can't find a location for this. The sequences are super short
        # so these won't get displayed anyways
        ({"country": "China", "division": "HuaShang"}, {"division": "Guangdong"}),
        # INDIA
        # -----
        # Fix typos
        ({"country": "India", "division": "West_Bengal"}, {"division": "West Bengal"}),
        ({"country": "India", "division": "Telengana"}, {"division": "Telangana"}),
        ({"country": "India", "division": "Rajsathan"}, {"division": "Rajasthan"}),
        # Jammu --> Jammu and Kashmir
        (
            {"country": "India", "division": "Jammu"},
            {"division": "Jammu and Kashmir", "location": "Jammu"},
        ),
        # Kargil --> Ladakh/Kargil
        (
            {"country": "India", "division": "Kargil"},
            {"division": "Ladakh", "location": "Kargil"},
        ),
        # Mumbai --> Maharashtra/Mumbai
        (
            {"country": "India", "division": "Mumbai"},
            {"division": "Maharashtra", "location": "Mumbai"},
        ),
        # MP -> Madhya Pradesh
        ({"country": "India", "division": "MP"}, {"division": "Madhya Pradesh"}),
        # ISRAEL
        # ------
        # Fix typos
        ({"country": "ISRAEL"}, {"country": "Israel"}),
        # South Coast District --> South District
        (
            {
                "country": "Israel",
                "division": ["South Coast District", "Southern District"],
            },
            {"division": "South District"},
        ),
        # Move cities into districts
        (
            {"country": "Israel", "division": ["ELAD"]},
            {"division": "Center District", "location": "El'ad"},
        ),
        (
            {"country": "Israel", "division": ["GANI TIKVA"]},
            {"division": "Center District", "location": "Ganei Tikva"},
        ),
        (
            {"country": "Israel", "division": ["GEDERA"]},
            {"division": "Center District", "location": "Gedera"},
        ),
        (
            {"country": "Israel", "division": ["KFAR-HABAD"]},
            {"division": "Center District", "location": "Kfar Chabad"},
        ),
        (
            {"country": "Israel", "division": ["NETANYA"]},
            {"division": "Center District", "location": "Netanya"},
        ),
        (
            {"country": "Israel", "division": ["PARDESIA"]},
            {"division": "Center District", "location": "Pardesia"},
        ),
        (
            {"country": "Israel", "division": ["PETAH-TIKVA"]},
            {"division": "Center District", "location": "Petah Tikva"},
        ),
        (
            {"country": "Israel", "division": ["RAANANA"]},
            {"division": "Center District", "location": "Ra'anana"},
        ),
        (
            {"country": "Israel", "division": ["RAMLA"]},
            {"division": "Center District", "location": "Ramla"},
        ),
        (
            {"country": "Israel", "division": ["SHOAM"]},
            {"division": "Center District", "location": "Shoham"},
        ),
        (
            {"country": "Israel", "division": ["YAVNE"]},
            {"division": "Center District", "location": "Yavne"},
        ),
        (
            {"country": "Israel", "division": ["YEHUD"]},
            {"division": "Center District", "location": "Yehud"},
        ),
        (
            {"country": "Israel", "division": ["Kefar Sava", "Kfar Saba"]},
            {"division": "Center District", "location": "Kfar Saba"},
        ),
        (
            {"country": "Israel", "division": ["HAIFA"]},
            {"division": "Haifa District", "location": "Haifa"},
        ),
        (
            {"country": "Israel", "division": ["BEIT SHEMESH"]},
            {"division": "Jerusalem District", "location": "Beit Shemesh"},
        ),
        (
            {"country": "Israel", "division": ["JERUSALEM"]},
            {"division": "Jerusalem District", "location": "Jerusalem"},
        ),
        (
            {"country": "Israel", "division": ["ELKANA"]},
            {"division": "Judea and Samaria Area", "location": "Elkana"},
        ),
        (
            {"country": "Israel", "division": ["HUKUK"]},
            {"division": "North District", "location": "Hukok"},
        ),
        (
            {"country": "Israel", "division": ["MAALOT-TARSHICHA"]},
            {"division": "North District", "location": "Ma'alot-Tarshiha"},
        ),
        (
            {"country": "Israel", "division": ["TIRAT ZVI"]},
            {"division": "North District", "location": "Tirat Zvi"},
        ),
        (
            {"country": "Israel", "division": ["TVERIA"]},
            {"division": "North District", "location": "Tiberias"},
        ),
        (
            {"country": "Israel", "division": ["ZEFAT"]},
            {"division": "North District", "location": "Safed"},
        ),
        (
            {"country": "Israel", "division": ["ASHDOD"]},
            {"division": "South District", "location": "Ashdod"},
        ),
        (
            {"country": "Israel", "division": ["KIRYAT-GAT"]},
            {"division": "South District", "location": "Kiryat Gat"},
        ),
        (
            {"country": "Israel", "division": ["NETIVOT"]},
            {"division": "South District", "location": "Netivot"},
        ),
        (
            {"country": "Israel", "division": ["BAT-YAM"]},
            {"division": "Tel Aviv District", "location": "Bat Yam"},
        ),
        (
            {"country": "Israel", "division": ["BENI-BRAK"]},
            {"division": "Tel Aviv District", "location": "Bnei Brak"},
        ),
        (
            {"country": "Israel", "division": ["HOLON"]},
            {"division": "Tel Aviv District", "location": "Holon"},
        ),
        (
            {"country": "Israel", "division": ["KIRYAT ONO"]},
            {"division": "Tel Aviv District", "location": "Kiryat Ono"},
        ),
        (
            {"country": "Israel", "division": ["RAMAT-GAN", "Ramat Gan"]},
            {"division": "Tel Aviv District", "location": "Ramat Gan"},
        ),
        (
            {"country": "Israel", "division": ["TEL-AVIV"]},
            {"division": "Tel Aviv District", "location": "Tel Aviv"},
        ),
        # JAPAN
        # -----
        # Remove unknown division
        ({"country": "Japan", "division": "unknown"}, {"division": -1}),
        # OMAN
        # ----
        # Fix typos
        ({"country": "Oman", "division": ["Muscta", "Musca"]}, {"division": "Muscat"}),
        # PAKISTAN
        # --------
        # Unabbreviate province names
        (
            {"country": "Pakistan", "division": "KPK"},
            {"division": "Khyber Pakhtunkhwa"},
        ),
        # SOUTH KOREA
        # -----------
        # Korea --> South Korea
        # I assume North Korea is not submitting genomes...
        ({"country": "Korea"}, {"country": "South Korea"}),
        # TAIWAN
        # ------
        # Fix typos
        (
            {"country": "Taiwan", "division": "New Taipei city"},
            {"division": "New Taipei City"},
        ),
        # THAILAND
        # --------
        # Fix typos
        (
            {"country": "Thailand", "division": "Phatum thani"},
            {"division": "Pathum Thani"},
        ),
        # VIETNAM
        # -------
        # Fix typos
        ({"country": "Vietnam", "division": "Quangning"}, {"division": "Quangninh"}),
        ({"country": "Vietnam", "division": "Thanh Hoa"}, {"division": "Thanhhoa"}),
        # BELGIUM
        # -------
        # Move Belgian cities into provinces
        # Since it is getting out of hand. the list is too big
        # And merge towns into parent municipalities
        # Antwerp
        (
            {
                "country": "Belgium",
                "division": [
                    "Antwerp",
                    "Antwerpen",
                    "Berchem",
                    "Berendrech-Zandvliet-Lillo",
                    "Borgerhout",
                    "Deurne",
                    "Ekeren",
                    "Hoboken",
                    "Merksem",
                    "Wilrijk",
                ],
            },
            {"division": "Antwerp", "location": "Antwerp"},
        ),
        (
            {"country": "Belgium", "division": "Balen"},
            {"division": "Antwerp", "location": "Balen"},
        ),
        (
            {"country": "Belgium", "division": ["Beerse", "Vlimmeren"]},
            {"division": "Antwerp", "location": "Beerse"},
        ),
        (
            {"country": "Belgium", "division": "Bonheiden"},
            {"division": "Antwerp", "location": "Bonheiden"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Bornem", "Hingene", "Wintam", "Mariekerke", "Weert"],
            },
            {"division": "Antwerp", "location": "Bornem"},
        ),
        (
            {"country": "Belgium", "division": "Brasschaat"},
            {"division": "Antwerp", "location": "Brasschaat"},
        ),
        (
            {"country": "Belgium", "division": ["Grobbendonk", "Bouwel"]},
            {"division": "Antwerp", "location": "Grobbendonk"},
        ),
        (
            {"country": "Belgium", "division": ["Hemisem", "Hemiksem"]},
            {"division": "Antwerp", "location": "Hemiksem"},
        ),
        (
            {"country": "Belgium", "division": ["Herselt", "Ramsel"]},
            {"division": "Antwerp", "location": "Herselt"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Laakdal",
                    "Laakdaal",
                    "Veerle",
                    "Eindhout",
                    "Vorst",
                    "Varendonk",
                    "Vorst-Meerlaar",
                    "Klein-Vorst",
                ],
            },
            {"division": "Antwerp", "location": "Laakdal"},
        ),
        (
            {"country": "Belgium", "division": ["Lille", "Gierle"]},
            {"division": "Antwerp", "location": "Lille"},
        ),
        (
            {"country": "Belgium", "division": "Kalmthout"},
            {"division": "Antwerp", "location": "Kalmthout"},
        ),
        (
            {"country": "Belgium", "division": "Kasterlee"},
            {"division": "Antwerp", "location": "Kasterlee"},
        ),
        (
            {"country": "Belgium", "division": ["Malle", "Oostmalle", "Westmalle"]},
            {"division": "Antwerp", "location": "Malle"},
        ),
        (
            {"country": "Belgium", "division": ["Nijlen", "Bevel", "Kessel"]},
            {"division": "Antwerp", "location": "Nijlen"},
        ),
        (
            {"country": "Belgium", "division": "Mol"},
            {"division": "Antwerp", "location": "Mol"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Olen",
                    "Olen-Centrum",
                    "Saint-Martins",
                    "Boekel",
                    "Onze-Lieve-Vrouw-Olen",
                    "Achter-Olen",
                    "St.-Jozef-Olen",
                    "Olen-Fabriek",
                ],
            },
            {"division": "Antwerp", "location": "Olen"},
        ),
        (
            {"country": "Belgium", "division": "Oud-Turnhout"},
            {"division": "Antwerp", "location": "Oud-Turnhout"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Puurs", "Breendonk", "Liezele", "Kalfort", "Ruisbroek"],
            },
            {"division": "Antwerp", "location": "Puurs"},
        ),
        (
            {"country": "Belgium", "division": ["Ravels", "Poppel", "Weelde"]},
            {"division": "Antwerp", "location": "Ravels"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Ranst", "Broechem", "Emblem", "Oelegem"],
            },
            {"division": "Antwerp", "location": "Ranst"},
        ),
        (
            {"country": "Belgium", "division": ["Retie", "Rethy", "Schoonbroek"]},
            {"division": "Antwerp", "location": "Retie"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Rijkevorsel",
                    "Achtel",
                    "Sint-Jozef-Rijkevorsel",
                    "Gammel",
                ],
            },
            {"division": "Antwerp", "location": "Rijkevorsel"},
        ),
        (
            {"country": "Belgium", "division": "Schoten"},
            {"division": "Antwerp", "location": "Schoten"},
        ),
        (
            {"country": "Belgium", "division": ["Sint-Amands", "Lippelo", "Oppuurs"]},
            {"division": "Antwerp", "location": "Sint-Amands"},
        ),
        (
            {"country": "Belgium", "division": "Turnhout"},
            {"division": "Antwerp", "location": "Turnhout"},
        ),
        (
            {"country": "Belgium", "division": "Vorselaar"},
            {"division": "Antwerp", "location": "Vorselaar"},
        ),
        (
            {"country": "Belgium", "division": "Oevel"},
            {"division": "Antwerp", "location": "Westerlo"},
        ),
        (
            {"country": "Belgium", "division": "Wijnegem"},
            {"division": "Antwerp", "location": "Wijnegem"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Willebroek",
                    "Blaasveld",
                    "Heindonk",
                    "Tisselt",
                    "Klein Willebroek",
                ],
            },
            {"division": "Antwerp", "location": "Willebroek"},
        ),
        (
            {"country": "Belgium", "division": ["Zoersel", "Halle", "St. Antonius"]},
            {"division": "Antwerp", "location": "Zoersel"},
        ),
        # East Flanders
        (
            {"country": "Belgium", "division": ["Aalter", "Lotenhulle"]},
            {"division": "East Flanders", "location": "Aalter"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Beveren",
                    "Doel",
                    "Haasdonk",
                    "Kallo",
                    "Kieldrecht",
                    "Melsele",
                    "Verrebroek",
                    "Vrasene",
                    "Beveren-Waas",
                ],
            },
            {"division": "East Flanders", "location": "Beveren"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Buggenhout", "Briel", "Opdorp", "Opstal"],
            },
            {"division": "East Flanders", "location": "Buggenhout"},
        ),
        (
            {"country": "Belgium", "division": "Deinze"},
            {"division": "East Flanders", "location": "Deinze"},
        ),
        (
            {"country": "Belgium", "division": "Dendermonde"},
            {"division": "East Flanders", "location": "Dendermonde"},
        ),
        (
            {"country": "Belgium", "division": ["Asper", "Gavere"]},
            {"division": "East Flanders", "location": "Gavere"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Gent",
                    "Ghent",
                    "Afsnee",
                    "Desteldonk",
                    "Drongen",
                    "Gentbrugge",
                    "Ledeberg",
                    "Mariakerke",
                    "Mendonk",
                    "Oostakker",
                    "Sint-Amandsberg",
                    "Sint-Denijs-Westrem",
                    "Sint-Kruis-Winkel",
                    "Wondelgem",
                    "Zwijnaarde",
                ],
            },
            {"division": "East Flanders", "location": "Ghent"},
        ),
        (
            {"country": "Belgium", "division": "Geraardsbergen"},
            {"division": "East Flanders", "location": "Geraardsbergen"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Hamme", "Kastel-Moerzeke", "St-Anna", "Zogge"],
            },
            {"division": "East Flanders", "location": "Hamme"},
        ),
        (
            {"country": "Belgium", "division": ["Kruibeke", "Bazel", "Rupelmonde"]},
            {"division": "East Flanders", "location": "Kruibeke"},
        ),
        (
            {"country": "Belgium", "division": ["Laarne", "Kalken"]},
            {"division": "East Flanders", "location": "Laarne"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Lievegem", "Waarschoot", "Lovendegem", "Zomergem"],
            },
            {"division": "East Flanders", "location": "Lievegem"},
        ),
        (
            {"country": "Belgium", "division": "Madegem"},
            {"division": "East Flanders", "location": "Madegem"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Merelbeke",
                    "Bottelare",
                    "Lemberge",
                    "Melsen",
                    "Munte",
                    "Schelderode",
                ],
            },
            {"division": "East Flanders", "location": "Merelbeke"},
        ),
        (
            {"country": "Belgium", "division": "Nevele"},
            {"division": "East Flanders", "location": "Nevele"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Oudenaarde",
                    "Bevere",
                    "Edelare",
                    "Eine",
                    "Ename",
                    "Herune",
                    "Leupegem",
                    "Mater",
                    "Melden",
                    "Mullem",
                    "Nederename",
                    "Volkegem",
                    "Welden",
                    "Ooike",
                ],
            },
            {"division": "East Flanders", "location": "Oudenaarde"},
        ),
        (
            {"country": "Belgium", "division": "Ronse"},
            {"division": "East Flanders", "location": "Ronse"},
        ),
        (
            {"country": "Belgium", "division": "Sint-Gillis-Waas"},
            {"division": "East Flanders", "location": "Sint-Gillis-Waas"},
        ),
        (
            {"country": "Belgium", "division": "Sint-Niklaas"},
            {"division": "East Flanders", "location": "Sint-Niklaas"},
        ),
        (
            {"country": "Belgium", "division": "Stekene"},
            {"division": "East Flanders", "location": "Stekene"},
        ),
        (
            {"country": "Belgium", "division": "Temse"},
            {"division": "East Flanders", "location": "Temse"},
        ),
        (
            {"country": "Belgium", "division": "Waasmunster"},
            {"division": "East Flanders", "location": "Waasmunster"},
        ),
        # Flemish Brabant
        (
            {"country": "Belgium", "division": "Aarschot"},
            {"division": "Flemish Brabant", "location": "Aarschot"},
        ),
        (
            {"country": "Belgium", "division": "Asse"},
            {"division": "Flemish Brabant", "location": "Asse"},
        ),
        (
            {"country": "Belgium", "division": "Beersel"},
            {"division": "Flemish Brabant", "location": "Beersel"},
        ),
        (
            {"country": "Belgium", "division": "Bierbeek"},
            {"division": "Flemish Brabant", "location": "Bierbeek"},
        ),
        (
            {"country": "Belgium", "division": "Boutersem"},
            {"division": "Flemish Brabant", "location": "Boutersem"},
        ),
        (
            {"country": "Belgium", "division": ["Diest", "Kaggevinne", "Molenstede"]},
            {"division": "Flemish Brabant", "location": "Diest"},
        ),
        (
            {"country": "Belgium", "division": "Dilbeek"},
            {"division": "Flemish Brabant", "location": "Dilbeek"},
        ),
        (
            {"country": "Belgium", "division": "Geetbets"},
            {"division": "Flemish Brabant", "location": "Geetbets"},
        ),
        (
            {"country": "Belgium", "division": "Grimbergen"},
            {"division": "Flemish Brabant", "location": "Grimbergen"},
        ),
        (
            {"country": "Belgium", "division": ["Haacht", "Wespelaar"]},
            {"division": "Flemish Brabant", "location": "Haacht"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Winksele", "Veltem-beisem", "Veltem-Beisem"],
            },
            {"division": "Flemish Brabant", "location": "Herent"},
        ),
        (
            {"country": "Belgium", "division": "Hoegaarden"},
            {"division": "Flemish Brabant", "location": "Hoegaarden"},
        ),
        (
            {"country": "Belgium", "division": "Holsbeek"},
            {"division": "Flemish Brabant", "location": "Holsbeek"},
        ),
        (
            {"country": "Belgium", "division": "Huldenberg"},
            {"division": "Flemish Brabant", "location": "Huldenberg"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Kampenhout", "Berg", "Buken", "Nederokkerzeel"],
            },
            {"division": "Flemish Brabant", "location": "Kampenhout"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Kortenaken",
                    "Hoeleden",
                    "Kersbeek-Miskom",
                    "Kersbeek-miskom",
                    "Ransberg",
                    "Waanrode",
                ],
            },
            {"division": "Flemish Brabant", "location": "Kortenaken"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Kortenberg", "Erps-Kwerps", "Everberg", "Meerbeek"],
            },
            {"division": "Flemish Brabant", "location": "Kortenberg"},
        ),
        (
            {"country": "Belgium", "division": "Kraainem"},
            {"division": "Flemish Brabant", "location": "Kraainem"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Landen",
                    "Attenhoven",
                    "Eliksem",
                    "Ezemaal",
                    "Laar",
                    "Neerlanden",
                    "Neerwinden",
                    "Overwinden",
                    "Rumsdorp",
                    "Waasmont",
                    "Walsbets",
                    "Walshoutem",
                    "Wange",
                    "Wezeren",
                ],
            },
            {"division": "Flemish Brabant", "location": "Landen"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Leuven", "Heverlee", "Kessel-Lo", "Ladeuze", "Wilsele"],
            },
            {"division": "Flemish Brabant", "location": "Leuven"},
        ),
        (
            {"country": "Belgium", "division": "Linter"},
            {"division": "Flemish Brabant", "location": "Linter"},
        ),
        (
            {"country": "Belgium", "division": "Lubbeek"},
            {"division": "Flemish Brabant", "location": "Lubbeek"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Oud-Heverlee",
                    "Blanden",
                    "Haasrode",
                    "Sint-Joris-Weert",
                    "Vaalbeek",
                ],
            },
            {"division": "Flemish Brabant", "location": "Oud-Heverlee"},
        ),
        (
            {"country": "Belgium", "division": ["Overijse"]},
            {"division": "Flemish Brabant", "location": "Overijse"},
        ),
        (
            {"country": "Belgium", "division": "Machelen"},
            {"division": "Flemish Brabant", "location": "Machelen"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Scherpenheuvel-Zichem",
                    "Scherpenheuvel-Zchem",
                    "Averbode",
                    "Messelbroek",
                    "Okselaar",
                    "Scherpenheuvel",
                    "Schoonderbuken",
                    "Keiberg",
                    "Testelt",
                    "Zichem",
                ],
            },
            {"division": "Flemish Brabant", "location": "Scherpenheuvel-Zichem"},
        ),
        (
            {"country": "Belgium", "division": "Sint-Genesius-Rode"},
            {"division": "Flemish Brabant", "location": "Sint-Genesius-Rode"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Steenokkerzeel",
                    "Humelgem",
                    "Wambeek",
                    "Perk",
                    "Huinhoven",
                    "Boekt",
                    "Melsbroek",
                ],
            },
            {"division": "Flemish Brabant", "location": "Steenokkerzeel"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Tervuren", "Duisberg", "Vossem", "Moorsel"],
            },
            {"division": "Flemish Brabant", "location": "Tervuren"},
        ),
        (
            {"country": "Belgium", "division": "Tielt-Winge"},
            {"division": "Flemish Brabant", "location": "Tielt-Winge"},
        ),
        (
            {"country": "Belgium", "division": "Tienen"},
            {"division": "Flemish Brabant", "location": "Tienen"},
        ),
        (
            {"country": "Belgium", "division": "Vilvoorde"},
            {"division": "Flemish Brabant", "location": "Vilvoorde"},
        ),
        (
            {"country": "Belgium", "division": "Wezembeek-Oppem"},
            {"division": "Flemish Brabant", "location": "Wezembeek-Oppem"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Zaventem",
                    "Sterrebeek",
                    "Nossegem",
                    "Sint-Stevens-Woluwe",
                ],
            },
            {"division": "Flemish Brabant", "location": "Zaventem"},
        ),
        (
            {"country": "Belgium", "division": "Zoutleeuw"},
            {"division": "Flemish Brabant", "location": "Zoutleeuw"},
        ),
        # Limburg
        (
            {"country": "Belgium", "division": "Alken"},
            {"division": "Limburg", "location": "Alken"},
        ),
        (
            {"country": "Belgium", "division": ["Beringen", "Koersel", "Paal"]},
            {"division": "Limburg", "location": "Beringen"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Bilzen",
                    "Beverst",
                    "Eigenbilzen",
                    "Grote-Spouwen",
                    "Hees",
                    "Hoelbeek",
                    "Kleine-Spouwen",
                    "Martenslinde",
                    "Mopertinen",
                    "Munsterbilzen",
                    "Rijkhoven",
                    "Rosmeer",
                    "Waltwilder",
                ],
            },
            {"division": "Limburg", "location": "Bilzen"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Bocholt", "Kaulillecity", "Reppel", "Lozen", "Kaulile"],
            },
            {"division": "Limburg", "location": "Bocholt"},
        ),
        (
            {"country": "Belgium", "division": "Bree"},
            {"division": "Limburg", "location": "Bree"},
        ),
        (
            {"country": "Belgium", "division": "Genk"},
            {"division": "Limburg", "location": "Genk"},
        ),
        (
            {"country": "Belgium", "division": ["Halen", "Loksbergen", "Zelem"]},
            {"division": "Limburg", "location": "Halen"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Hasselt",
                    "Sint-Lambrechts-Herk",
                    "Wimmertingen",
                    "Kermt",
                    "Spalbeek",
                    "Kuringen",
                    "Stokrooie",
                    "Stevoort",
                    "Runkst",
                    "Kiewit",
                    "Godsheide",
                    "Rapertingen",
                ],
            },
            {"division": "Limburg", "location": "Hasselt"},
        ),
        (
            {"country": "Belgium", "division": "Houthalen-Helchteren"},
            {"division": "Limburg", "location": "Houthalen-Helchteren"},
        ),
        (
            {"country": "Belgium", "division": "Lanaken"},
            {"division": "Limburg", "location": "Lanaken"},
        ),
        (
            {"country": "Belgium", "division": ["Leopoldsburg", "Heppen", "Beverloo"]},
            {"division": "Limburg", "location": "Leopoldsburg"},
        ),
        (
            {"country": "Belgium", "division": "Lommel"},
            {"division": "Limburg", "location": "Lommel"},
        ),
        (
            {"country": "Belgium", "division": "Lummen"},
            {"division": "Limburg", "location": "Lummen"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Maaseik",
                    "Neeroeteren",
                    "Opoeteren",
                    "Aldeneik",
                    "Heppeneert",
                    "Wurfeld",
                    "'t Ven",
                    "Gremeslo",
                    "Berg",
                    "Schootsheide",
                    "Voorshoven",
                    "Waterloos",
                    "De Riet",
                    "Dorne",
                ],
            },
            {"division": "Limburg", "location": "Maaseik"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Maasmechelen",
                    "Mechelen-aan-de-Maas",
                    "Vucht",
                    "Leut",
                    "Meeswijk",
                    "Uikhoven",
                    "Eisden",
                    "Opgrimbie",
                    "Boorsem",
                    "Kotem",
                ],
            },
            {"division": "Limburg", "location": "Maasmechelen"},
        ),
        (
            {"country": "Belgium", "division": ["Neerpelt", "Overpelt"]},
            {"division": "Limburg", "location": "Pelt"},
        ),
        (
            {"country": "Belgium", "division": "Nieuwerkerken"},
            {"division": "Limburg", "location": "Nieuwerkerken"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Peer",
                    "Grote-Brogel",
                    "Kleine-Brogel",
                    "Wauberg",
                    "Erpekom",
                    "Wijchmaal",
                ],
            },
            {"division": "Limburg", "location": "Peer"},
        ),
        (
            {"country": "Belgium", "division": "Riemst"},
            {"division": "Limburg", "location": "Riemst"},
        ),
        (
            {"country": "Belgium", "division": "Sint-Truiden"},
            {"division": "Limburg", "location": "Sint-Truiden"},
        ),
        (
            {"country": "Belgium", "division": "Tessenderlo"},
            {"division": "Limburg", "location": "Tessenderlo"},
        ),
        (
            {"country": "Belgium", "division": ["Heusden-Zolder", "Zolder"]},
            {"division": "Limburg", "location": "Heusden-Zolder"},
        ),
        (
            {"country": "Belgium", "division": "Zonhoven"},
            {"division": "Limburg", "location": "Zonhoven"},
        ),
        # West Flanders
        (
            {"country": "Belgium", "division": "Avelgem"},
            {"division": "West Flanders", "location": "Avelgem"},
        ),
        (
            {"country": "Belgium", "division": "Anzegem"},
            {"division": "West Flanders", "location": "Anzegem"},
        ),
        (
            {"country": "Belgium", "division": ["Blankenberge", "Uitkerke"]},
            {"division": "West Flanders", "location": "Blankenberge"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Bruges",
                    "Brugge",
                    "Assebroek",
                    "Koolkerke",
                    "Sint-Andries",
                    "Sint-Michiels",
                    "Sint-Kruis",
                    "Dudzele",
                    "Lissewege",
                    "Sint-Jozef",
                    "Sint-Pieters",
                ],
            },
            {"division": "West Flanders", "location": "Bruges"},
        ),
        (
            {"country": "Belgium", "division": "Deerlijk"},
            {"division": "West Flanders", "location": "Deerlijk"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Diksmuide",
                    "Beerst",
                    "Esen",
                    "Kaaskerke",
                    "Keiem",
                    "Lampernisse",
                    "Leke",
                    "Nieuwkapelle",
                    "Oostkerke",
                    "Oudakapelle",
                    "Pervijze",
                    "Sint-Jacobs-Kapelle",
                    "Stuivekenskerke",
                    "Vladslo",
                    "Woumen",
                ],
            },
            {"division": "West Flanders", "location": "Diksmuide"},
        ),
        (
            {"country": "Belgium", "division": ["Harelbeke", "Bavikhove", "Hulste"]},
            {"division": "West Flanders", "location": "Harelbeke"},
        ),
        (
            {"country": "Belgium", "division": "Knokke-Heist"},
            {"division": "West Flanders", "location": "Knokke-Heist"},
        ),
        (
            {"country": "Belgium", "division": "Kuurne"},
            {"division": "West Flanders", "location": "Kuurne"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Lo-Reninge",
                    "Lo",
                    "Noordschote",
                    "Pollinkhove",
                    "Reninge",
                ],
            },
            {"division": "West Flanders", "location": "Lo-Reninge"},
        ),
        (
            {"country": "Belgium", "division": "Oostrozebeke"},
            {"division": "West Flanders", "location": "Oostrozebeke"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Ostend", "Oostende", "Raversijde", "Stene", "Zandvoorde"],
            },
            {"division": "West Flanders", "location": "Ostend"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Oudenburg", "Ettelgem", "Roksem", "Westkerke"],
            },
            {"division": "West Flanders", "location": "Oudenburg"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Poperinge",
                    "Krombeke",
                    "Proven",
                    "Reningelst",
                    "Roesbrugge-Haringe",
                    "Roesbrugge",
                    "Haringe",
                    "Watou",
                    "Sint-Jan-Ter-Biezen",
                    "Abele",
                ],
            },
            {"division": "West Flanders", "location": "Poperinge"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Torhout",
                    "Wijnendale",
                    "Sint-Henricus",
                    "De Driekoningen",
                ],
            },
            {"division": "West Flanders", "location": "Torhout"},
        ),
        (
            {"country": "Belgium", "division": ["Wervik", "Geluwe"]},
            {"division": "West Flanders", "location": "Wervik"},
        ),
        (
            {"country": "Belgium", "division": "Wevelgem"},
            {"division": "West Flanders", "location": "Wevelgem"},
        ),
        (
            {"country": "Belgium", "division": ["Wingene", "Zwevezele"]},
            {"division": "West Flanders", "location": "Wingene"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Ypres",
                    "Ieper",
                    "Boezinge",
                    "Brielen",
                    "Dikkebus",
                    "Elverdinge",
                    "Hollebekke",
                    "Sint-Jan-Vlamertinge",
                    "Vlamertinge",
                    "Voormezele",
                    "Zillebeke",
                    "Zuidschote",
                ],
            },
            {"division": "West Flanders", "location": "Ypres"},
        ),
        # Hainaut
        (
            {"country": "Belgium", "division": ["Aiseau", "Presles"]},
            {"division": "Hainaut", "location": "Aiseau-Presles"},
        ),
        (
            {"country": "Belgium", "division": ["Ath", "Arbre"]},
            {"division": "Hainaut", "location": "Ath"},
        ),
        (
            {"country": "Belgium", "division": "Bassilly"},
            {"division": "Hainaut", "location": "Silly"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Quevaucamps",
                    "Wadelincourt",
                    "Basecles",
                    "Ramegnies",
                    "Thumaide",
                    "Aubechies",
                    "Ellignies-Sainte-Anne",
                ],
            },
            {"division": "Hainaut", "location": "Belœil"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Bernissart",
                    "Blaton",
                    "Harchies",
                    "Pommeroeul",
                    "Poomerœul",
                ],
            },
            {"division": "Hainaut", "location": "Bernissart"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Binche",
                    "Bray",
                    "Buvrinnes",
                    "Epinois",
                    "Leval-Trahegnies",
                    "Peronnes-lez-Binche",
                    "Ressaix",
                    "Waudrez",
                ],
            },
            {"division": "Hainaut", "location": "Binche"},
        ),
        (
            {"country": "Belgium", "division": "Boussu"},
            {"division": "Hainaut", "location": "Boussu"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Brugelette",
                    "Cambron-Casteau",
                    "Attre",
                    "Mevergnies-lez-Lens",
                    "Gages",
                ],
            },
            {"division": "Hainaut", "location": "Brugelette"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Braine-Le-Comte",
                    "Hennuyères",
                    "Henripont",
                    "Petit-Roeulx-lez-Braine",
                    "Ronquières",
                    "Steenkerque",
                    "Steenkerke",
                ],
            },
            {"division": "Hainaut", "location": "Braine-Le-Comte"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Brunehaut",
                    "Bléharies",
                    "Guignies",
                    "Hollain",
                    "Jollain-Merlin",
                    "Wez-Velvain",
                    "Lesdain",
                    "Laplaigne",
                    "Rongy",
                    "Howardries",
                ],
            },
            {"division": "Hainaut", "location": "Brunehaut"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Péruwelz",
                    "Bon-Secours",
                    "Roucourt",
                    "Braffe",
                    "Bury",
                    "Baugnies",
                    "Wasmes-Audemez-Briffoeil",
                    "Brasménil",
                    "Wiers",
                ],
            },
            {"division": "Hainaut", "location": "Péruwelz"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Charleroi",
                    "Dampremy",
                    "Lodenlinsart",
                    "Gilly",
                    "Montignies-Sur-Sambre",
                    "Montignies-sur-sambre",
                    "Couillet",
                    "Marcinelle",
                    "Mont-sur-Marchienne",
                    "Marchienne-au-Pont",
                    "Monceau-sur-Sambre",
                    "Monceau-sur-sambre",
                    "Goutroux",
                    "Roux",
                    "Jumet",
                    "Gosselies",
                    "Ransart",
                ],
            },
            {"division": "Hainaut", "location": "Charleroi"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Chatelet", "Châtelet", "Bouffioulx", "Chatelineau"],
            },
            {"division": "Hainaut", "location": "Châtelet"},
        ),
        (
            {"country": "Belgium", "division": ["Colfontaine", "Confontaine"]},
            {"division": "Hainaut", "location": "Colfontaine"},
        ),
        (
            {"country": "Belgium", "division": "Cuesmes"},
            {"division": "Hainaut", "location": "Cuesmes"},
        ),
        (
            {"country": "Belgium", "division": "Dour"},
            {"division": "Hainaut", "location": "Dour"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Ecaussinnes",
                    "Écaussinees",
                    "Écaussinees-d'Enghien",
                    "Marche-lez-Écaussinnes",
                    "Écaussinees-Lalaing",
                ],
            },
            {"division": "Hainaut", "location": "Écaussinees"},
        ),
        (
            {"country": "Belgium", "division": "Edingen"},
            {"division": "Hainaut", "location": "Edingen"},
        ),
        (
            {"country": "Belgium", "division": "Ellezelles"},
            {"division": "Hainaut", "location": "Ellezelles"},
        ),
        (
            {"country": "Belgium", "division": ["Flobecq", "Vloesberg"]},
            {"division": "Hainaut", "location": "Flobecq"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Fontaine-l Eveque", "Forchies-la-marche"],
            },
            {"division": "Hainaut", "location": "Fontaine-l'Évêque"},
        ),
        (
            {"country": "Belgium", "division": "Frameries"},
            {"division": "Hainaut", "location": "Frameries"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Frasnes-lez-Anvaing",
                    "Arc",
                    "Ainières",
                    "Wattripont",
                    "Anvaing",
                    "Buissenal",
                    "Moustier",
                    "Oeudeghien",
                    "Ellignies",
                    "Dergneau",
                    "Saint-Sauveur",
                    "Hacquegnies",
                    "Herquegnies",
                ],
            },
            {"division": "Hainaut", "location": "Frasnes-lez-Anvaing"},
        ),
        (
            {"country": "Belgium", "division": "Haulchin"},
            {"division": "Hainaut", "location": "Haulchin"},
        ),
        (
            {"country": "Belgium", "division": ["Hensies", "Thulin"]},
            {"division": "Hainaut", "location": "Hensies"},
        ),
        (
            {"country": "Belgium", "division": "Honnelles"},
            {"division": "Hainaut", "location": "Honnelles"},
        ),
        (
            {"country": "Belgium", "division": "Jurbise"},
            {"division": "Hainaut", "location": "Jurbise"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Comines-Warneton",
                    "Komen",
                    "Comines",
                    "Comines-ten-Brielen",
                    "Houthem",
                    "Warneton",
                    "Waasten",
                    "Bas-Warneton",
                    "Neer-Waasten",
                    "Ploegsteert",
                    "Le Bizet",
                ],
            },
            {"division": "Hainaut", "location": "Comines-Warneton"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Lens",
                    "Bauffe",
                    "Cambron-Saint-Vincent",
                    "Lombise",
                    "Montignies-lez-Lens",
                ],
            },
            {"division": "Hainaut", "location": "Lens"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "La Louvriere",
                    "La Louviere",
                    "La louvière",
                    "La Louvière",
                    "Haine-Saint-Paul",
                    "Haine-Saint-Pierre",
                    "Saint-Vaast",
                    "Trivieres",
                    "Boussoit",
                    "Houdeng-Aimeries",
                    "Houdeng-Gœgnies",
                    "Maurage",
                    "Strepy-Bracquegnies",
                ],
            },
            {"division": "Hainaut", "location": "La Louvière"},
        ),
        (
            {"country": "Belgium", "division": "Le roeulx"},
            {"division": "Hainaut", "location": "Le Rœulx"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Les Bons Villers",
                    "Frasnes-lez-Gosselies",
                    "Mellet",
                    "Reves",
                    "Villers-Pierwin",
                    "Wayaux",
                ],
            },
            {"division": "Hainaut", "location": "Les Bons Villers"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Lessines",
                    "Ghoy",
                    "Bois-de-Lessines",
                    "Lessenbos",
                    "Ogy",
                    "Oseke",
                    "Papignies",
                    "Papegem",
                    "Deux-Acren",
                    "Acren",
                    "Wannebecq",
                    "Wannabeek",
                    "Ollignies",
                    "Woelingen",
                ],
            },
            {"division": "Hainaut", "location": "Lessines"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Merbes-le-Chateau",
                    "Merbes-le-Château",
                    "Fontaine-Valmont",
                    "Labuissiere",
                    "Merbes-Sainte-Marie",
                ],
            },
            {"division": "Hainaut", "location": "Merbes-le-Château"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Mons",
                    "Mon",
                    "Courcelles",
                    "Havre",
                    "Jemappes",
                    "Trazegnies",
                    "Flénu",
                    "Ghlin",
                    "Shape",
                ],
            },
            {"division": "Hainaut", "location": "Mons"},
        ),
        (
            {"country": "Belgium", "division": "Montigny-le-Tilleul"},
            {"division": "Hainaut", "location": "Montigny-le-Tilleul"},
        ),
        (
            {"country": "Belgium", "division": "Morlanwelz"},
            {"division": "Hainaut", "location": "Morlanwelz"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Mouscron",
                    "Dottignies",
                    "Dottenijs",
                    "Luingne",
                    "Lowingen",
                    "Herseaux",
                    "Herzeeuw",
                ],
            },
            {"division": "Hainaut", "location": "Mouscron"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Pont-à-Celles",
                    "Buzet",
                    "Liberchies",
                    "Luttre",
                    "Obaix",
                    "Thimeon",
                    "Viesville",
                ],
            },
            {"division": "Hainaut", "location": "Pont-à-Celles"},
        ),
        (
            {"country": "Belgium", "division": "Quaregnon"},
            {"division": "Hainaut", "location": "Quaregnon"},
        ),
        (
            {"country": "Belgium", "division": ["Quevy", "Givry", "Quévy"]},
            {"division": "Hainaut", "location": "Quevy"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Saint-Ghislain",
                    "Tertre",
                    "Villerot",
                    "Baudour",
                    "Neufmaison",
                ],
            },
            {"division": "Hainaut", "location": "Saint-Ghislain"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Soignies",
                    "Casteau",
                    "Chaussee-Notre-Dame-Louvignies",
                    "Horrues",
                    "Neufvilles",
                    "Naast",
                    "Thieusies",
                    "Saisinne",
                    "Zinnik",
                ],
            },
            {"division": "Hainaut", "location": "Soignies"},
        ),
        (
            {"country": "Belgium", "division": "Seneffe"},
            {"division": "Hainaut", "location": "Seneffe"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Tournai",
                    "Doornik",
                    "Ere",
                    "Saint-Maur",
                    "Orcq",
                    "Esplechin",
                    "Froyennes",
                    "Froidmont",
                    "Willemeau",
                    "Ramegnies-Chin",
                    "Templeuve",
                    "Chercq",
                    "Blandain",
                    "Hertain",
                    "Lamain",
                    "Marquain",
                    "Gaurain-Ramecroix",
                    "Havinnes",
                    "Beclers",
                    "Thimougnies",
                    "Barry",
                    "Maulde",
                    "Vaulx",
                    "Vezon",
                    "Kain",
                    "Melles",
                    "Quartes",
                    "Rumillies",
                    "Mont-Saint-Aubert",
                    "Mourcourt",
                    "Warchin",
                ],
            },
            {"division": "Hainaut", "location": "Tournai"},
        ),
        # Liege
        (
            {"country": "Belgium", "division": "Ans"},
            {"division": "Liege", "location": "Ans"},
        ),
        (
            {"country": "Belgium", "division": ["Baelen", "Membach"]},
            {"division": "Liege", "location": "Baelen"},
        ),
        (
            {"country": "Belgium", "division": "Berloz"},
            {"division": "Liege", "location": "Berloz"},
        ),
        (
            {"country": "Belgium", "division": "Bierset"},
            {"division": "Liege", "location": "Grâce-Hollogne"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Bullingen",
                    "Honsfeld",
                    "Hunningen",
                    "Murringen",
                    "Rocherath",
                    "Krinkelt",
                    "Wirtzfeld",
                    "Manderfeld",
                    "Afst",
                    "Allmuthen",
                    "Andlermuhle",
                    "Berterath",
                    "Buchholz",
                    "Eimerscheid",
                    "Hasenvenn",
                    "Hergersberg",
                    "Holzheim",
                    "Hullscheid",
                    "Igelmonder Hof",
                    "Igelmondermuhle",
                    "Kehr",
                    "Krewinkel",
                    "Lanzerath",
                    "Losheimergraben",
                    "Medendorf",
                    "Merlscheid",
                    "Weckerath",
                ],
            },
            {"division": "Liege", "location": "Bullingen"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Burdinne",
                    "Hanneche",
                    "Lamontzee",
                    "Marneffe",
                    "Oteppe",
                    "Vissoul",
                ],
            },
            {"division": "Liege", "location": "Burdinne"},
        ),
        (
            {"country": "Belgium", "division": "Couthuin"},
            {"division": "Liege", "location": "Heron"},
        ),
        (
            {"country": "Belgium", "division": "Dalhem"},
            {"division": "Liege", "location": "Dalhem"},
        ),
        (
            {"country": "Belgium", "division": "Eupen"},
            {"division": "Liege", "location": "Eupen"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Hannut",
                    "Abolens",
                    "Avernas-le-Bauduin",
                    "Avin",
                    "Bertrée",
                    "Blehen",
                    "Cras-Avernas",
                    "Crehen",
                    "Grand-Hallet",
                    "Les-Saint-Remy",
                    "Merdorp",
                    "Moxhe",
                    "Petit-Hallet",
                    "Poucet",
                    "Thisnes",
                    "Trognée",
                    "Truielingen",
                    "Villers-le-Peuplier",
                    "Wansin",
                ],
            },
            {"division": "Liege", "location": "Hannut"},
        ),
        (
            {"country": "Belgium", "division": "Henri-chapelle"},
            {"division": "Liege", "location": "Welkenraedt"},
        ),
        (
            {"country": "Belgium", "division": "Herstal"},
            {"division": "Liege", "location": "Herstal"},
        ),
        (
            {"country": "Belgium", "division": ["Herve", "Chaineux"]},
            {"division": "Liege", "location": "Herve"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Juprelle",
                    "Fexhe-Slins",
                    "Lantin",
                    "Paifve",
                    "Slins",
                    "Villers-Saint-Simeon",
                    "Voroux-lez-Liers",
                    "Voroux-les-Liers",
                    "Wihogne",
                    "Nudorp",
                ],
            },
            {"division": "Liege", "location": "Juprelle"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Kelmis",
                    "La Calamine",
                    "Hergenrath",
                    "Neu-Moresnet",
                    "Vieille Montagne",
                ],
            },
            {"division": "Liege", "location": "Kelmis"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Liege",
                    "Liège",
                    "Angleur",
                    "Bressoux",
                    "Chenee",
                    "Glain",
                    "Grivegnee",
                    "Jupille-sur-Meuse",
                    "Rocourt",
                    "Wandre",
                ],
            },
            {"division": "Liege", "location": "Liège"},
        ),
        (
            {"country": "Belgium", "division": "Limbourg"},
            {"division": "Liege", "location": "Limbourg"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Remicourt", "Momalle", "Hodeige", "Lamine", "Pousset"],
            },
            {"division": "Liege", "location": "Remicourt"},
        ),
        (
            {"country": "Belgium", "division": "Saive"},
            {"division": "Liege", "location": "Blegny"},
        ),
        (
            {"country": "Belgium", "division": "Seraing"},
            {"division": "Liege", "location": "Seraing"},
        ),
        (
            {"country": "Belgium", "division": "Soumagne"},
            {"division": "Liege", "location": "Soumagne"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Sprimont",
                    "Dolembreux",
                    "Gomze-Andoumont",
                    "Rouvreux",
                    "Louveigne",
                    "Banneux",
                    "Damre",
                    "Florze",
                    "Lince",
                    "Ogne",
                    "Presseux",
                    "Rivage",
                ],
            },
            {"division": "Liege", "location": "Sprimont"},
        ),
        (
            {"country": "Belgium", "division": "Theux"},
            {"division": "Liege", "location": "Theux"},
        ),
        (
            {"country": "Belgium", "division": ["Antheit", "Wanze"]},
            {"division": "Liege", "location": "Wanze"},
        ),
        (
            {"country": "Belgium", "division": "Yernee-fraineux"},
            {"division": "Liege", "location": "Nandrin"},
        ),
        # Luxembourg
        # ...
        # Namur
        (
            {"country": "Belgium", "division": "Andenne"},
            {"division": "Namur", "location": "Andenne"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Ciney",
                    "Achene",
                    "Braibant",
                    "Chevetogne",
                    "Conneux",
                    "Leignon",
                    "Pessoux",
                    "Serinchamps",
                    "Sovet",
                    "Chapois",
                ],
            },
            {"division": "Namur", "location": "Ciney"},
        ),
        (
            {"country": "Belgium", "division": "Couvin"},
            {"division": "Namur", "location": "Couvin"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Aisemont",
                    "Le Roux",
                    "Sart-Eustache",
                    "Sart-Saint-Laurent",
                    "Vitrival",
                ],
            },
            {"division": "Namur", "location": "Fosses-la-Ville"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Gembloux",
                    "Beuzet",
                    "Bossiere",
                    "Bothey",
                    "Corroy-le-Chateau",
                    "Ernage",
                    "Grand-Leez",
                    "Grand-Manil",
                    "Isnes",
                    "Lonzee",
                    "Mazy",
                    "Sauveniere",
                ],
            },
            {"division": "Namur", "location": "Gembloux"},
        ),
        (
            {"country": "Belgium", "division": "Havelange"},
            {"division": "Namur", "location": "Havelange"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Namur",
                    "Beez",
                    "Belgrade",
                    "Saint-Servais",
                    "Saint-Marc",
                    "Bouge",
                    "Champion",
                    "Daussoulx",
                    "Flawinne",
                    "Malonne",
                    "Suarlee",
                    "Temploux",
                    "Vedrin",
                    "Boninne",
                    "Cognelee",
                    "Gelbressee",
                    "Marche-les-Dames",
                    "Dave",
                    "Jambes",
                    "Naninne",
                    "Wepion",
                    "Wierde",
                    "Erpent",
                    "Lives-sur-Meuse",
                    "Loyers",
                ],
            },
            {"division": "Namur", "location": "Namur"},
        ),
        # Walloon Brabant
        (
            {
                "country": "Belgium",
                "division": [
                    "Braine-l alleud",
                    "Braine-l'alleud",
                    "Braine-L'Alleud",
                    "Ophain-Bois-Seigneur-Isaac",
                    "Lillois-Witterzée",
                    "Lillois-witterzée",
                ],
            },
            {"division": "Walloon Brabant", "location": "Braine-l'Alleud"},
        ),
        (
            {"country": "Belgium", "division": "Braine-le-Chateau"},
            {"division": "Walloon Brabant", "location": "Braine-le-Château"},
        ),
        (
            {"country": "Belgium", "division": "Genappe"},
            {"division": "Walloon Brabant", "location": "Genappe"},
        ),
        (
            {"country": "Belgium", "division": "Lasne"},
            {"division": "Walloon Brabant", "location": "Lasne"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Rebecq",
                    "Bierghes",
                    "Bierk",
                    "Rebecq-Rognon",
                    "Roosbeek",
                    "Quenast",
                    "Kenast",
                ],
            },
            {"division": "Walloon Brabant", "location": "Rebecq"},
        ),
        (
            {
                "country": "Belgium",
                "division": ["Rixensart", "Rosieres", "Rozieren", "Genval"],
            },
            {"division": "Walloon Brabant", "location": "Rixensart"},
        ),
        (
            {"country": "Belgium", "division": "Tubize"},
            {"division": "Walloon Brabant", "location": "Tubize"},
        ),
        (
            {"country": "Belgium", "division": "Waterloo"},
            {"division": "Walloon Brabant", "location": "Waterloo"},
        ),
        # Brussels-Capital Region
        (
            {"country": "Belgium", "division": "Anderlecht"},
            {"division": "Brussels-Capital Region", "location": "Anderlecht"},
        ),
        (
            {"country": "Belgium", "division": ["Brussel", "Brussels", "Laken"]},
            {"division": "Brussels-Capital Region", "location": "Brussels"},
        ),
        (
            {"country": "Belgium", "division": ["Elsene", "Ixelles"]},
            {"division": "Brussels-Capital Region", "location": "Ixelles"},
        ),
        (
            {"country": "Belgium", "division": "Etterbeek"},
            {"division": "Brussels-Capital Region", "location": "Etterbeek"},
        ),
        (
            {"country": "Belgium", "division": "Vorst"},
            {"division": "Brussels-Capital Region", "location": "Forest"},
        ),
        (
            {"country": "Belgium", "division": "Ganshoren"},
            {"division": "Brussels-Capital Region", "location": "Ganshoren"},
        ),
        (
            {"country": "Belgium", "division": "Jette"},
            {"division": "Brussels-Capital Region", "location": "Jette"},
        ),
        (
            {"country": "Belgium", "division": "Oudergem"},
            {"division": "Brussels-Capital Region", "location": "Auderghem"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Molenbeek-Saint-Jean",
                    "Sint-Jans-Molenbeek",
                    "Molenbeek",
                ],
            },
            {"division": "Brussels-Capital Region", "location": "Molenbeek-Saint-Jean"},
        ),
        (
            {"country": "Belgium", "division": "Sint-Agatha-Berchem"},
            {"division": "Brussels-Capital Region", "location": "Sint-Agatha-Berchem"},
        ),
        (
            {"country": "Belgium", "division": "Ukkel"},
            {"division": "Brussels-Capital Region", "location": "Uccle"},
        ),
        (
            {"country": "Belgium", "division": "Watermaal-Bosvoorde"},
            {"division": "Brussels-Capital Region", "location": "Watermael-Boitsfort"},
        ),
        (
            {
                "country": "Belgium",
                "division": [
                    "Sint-Pieter-Woluwe",
                    "Sint-Pieters-Woluwe",
                    "Sint-pieters-woluwe",
                ],
            },
            {"division": "Brussels-Capital Region", "location": "Woluwe-Saint-Pierre"},
        ),
        # Can't find a division for this location
        ({"country": "Belgium", "division": "Mepel"}, {"division": -1}),
        # BOSNIA AND HERZEGOVINA
        # ----------------------
        # Fix typos
        ({"country": "Bosnia and Hercegovina"}, {"country": "Bosnia and Herzegovina"}),
        # CROATIA
        # -------
        # Fix typos
        (
            {"country": "Croatia", "division": "Varazdin"},
            {"division": "Varaždin County"},
        ),
        # CZECHIA
        # -------
        # Rename from "Czech Republic" --> "Czechia"
        ({"country": "Czech Republic"}, {"country": "Czechia"}),
        # Move "Prague" --> "Central Czhechia and Prague"
        (
            {"country": "Czechia", "division": "Prague"},
            {"division": "Central Czechia and Prague"},
        ),
        # Fix typos
        (
            {"country": "Czechia", "division": "Vysocina Region"},
            {"division": "Vysocina"},
        ),
        # DENMARK
        # -------
        # Remove Unknown region
        ({"country": "Denmark", "division": "Unknown"}, {"division": -1}),
        # FRANCE
        # ------
        # Move Saint-Saulve from Belgium to France
        (
            {"country": "Belgium", "division": "Saint-Saulve"},
            {"country": "France", "division": "Hauts-de-France"},
        ),
        # Move Roucourt from belgium to France
        (
            {"country": "Belgium", "division": "Roucourt"},
            {"country": "France", "division": "Hauts-de-France"},
        ),
        # Move Flines-les-Mortagne (Flines-les-Morta) from Belgium into France
        (
            {"country": "Belgium", "division": ["Flines-Les-Morta"]},
            {
                "country": "France",
                "division": "Hauts-de-France",
                "location": "Flines-lès-Mortagne",
            },
        ),
        # Fix typos
        (
            {
                "country": "France",
                "division": [
                    "Bourgogne Franche comte",
                    "Bourgogne-France-Comté",
                    "Bourgogne-Franche Comte",
                    "Bourgogne-France-Comte",
                    "Bourgogne",
                ],
            },
            {"division": "Bourgogne-Franche-Comté"},
        ),
        # Unabbreviate province names, unify province names
        (
            {"country": "France", "division": "ARA"},
            {"division": "Auvergne-Rhône-Alpes"},
        ),
        (
            {"country": "France", "division": "Centre - Val de Loire"},
            {"division": "Centre-Val de Loire"},
        ),
        (
            {"country": "France", "division": ["Grand-Est", "Grand-est"]},
            {"division": "Grand Est"},
        ),
        (
            {"country": "France", "division": ["Hauts De France", "Hauts de France"]},
            {"division": "Hauts-de-France"},
        ),
        (
            {
                "country": "France",
                "division": ["IDF", "Ile De France", "Ile de France", "Ile-de-France"],
            },
            {"division": "Île-de-France"},
        ),
        # Move Toulouse to Occitanie
        (
            {"country": "France", "division": "Toulouse"},
            {"division": "Occitanie", "location": "Toulouse"},
        ),
        # Fix typos
        (
            {"country": "France", "location": ["Chateau-Thierry", "Chateau-thierry"]},
            {"location": "Château-Thierry"},
        ),
        (
            {"country": "France", "location": ["Crépy-en -Valois", "Crépy en Valois"]},
            {"location": "Crépy-en-Valois"},
        ),
        (
            {"country": "France", "location": ["Asnieres sur Seine"]},
            {"location": "Asnières-sur-Seine"},
        ),
        (
            {"country": "France", "location": ["Moussy le Neuf"]},
            {"location": "Moussy Le Neuf"},
        ),
        (
            {"country": "France", "location": ["Gondrecourt-le-chateau"]},
            {"location": "Gondrecourt-le-Chateau"},
        ),
        ({"country": "France", "location": ["Compiegne"]}, {"location": "Compiègne"}),
        # Move Paris into Ile-de-France
        (
            {"country": "France", "division": "Paris"},
            {"division": "Île-de-France", "location": "Paris"},
        ),
        # GEORGIA
        # -------
        # Move to Asia
        ({"region": "Europe", "country": "Georgia"}, {"region": "Asia"}),
        # GERMANY
        # -------
        # Move Munich division into Bavaria division
        (
            {"country": "Germany", "division": "Munich"},
            {"division": "Bavaria", "location": "Munich"},
        ),
        # Duesseldorf -> North Rhine Westphalia
        (
            {"country": "Germany", "division": "Duesseldorf"},
            {"division": "North Rhine-Westphalia", "location": "Duesseldorf"},
        ),
        # Frankfurt -> Hesse
        (
            {"country": "Germany", "division": "Frankfurt"},
            {"division": "Hesse", "location": "Frankfurt"},
        ),
        # Rostock -> Mecklenburg-Vorpommern
        (
            {"country": "Germany", "division": "Rostock"},
            {"division": "Mecklenburg-Vorpommern", "location": "Rostock"},
        ),
        # Fix typos
        (
            {"country": "Germany", "division": "Baden-Wuerttemberg"},
            {"division": "Baden-Württemberg"},
        ),
        (
            {"country": "Germany", "division": "North Rhine Westphalia"},
            {"division": "North Rhine-Westphalia"},
        ),
        # ITALY
        # -----
        # Fix typos
        (
            {"country": "Italy", "division": ["Lombarida", "Lombardia"]},
            {"division": "Lombardy"},
        ),
        # Castel di Sangro -> Abruzzo
        (
            {"country": "Italy", "division": "Castel di Sangro"},
            {"division": "Abruzzo", "location": "Castel di Sangro"},
        ),
        # Teramo -> Abruzzo
        (
            {"country": "Italy", "division": "Teramo"},
            {"division": "Abruzzo", "location": "Teramo"},
        ),
        # Cagliari -> Sardinia
        (
            {"country": "Italy", "division": "Cagliari"},
            {"division": "Sardinia", "location": "Cagliari"},
        ),
        # Milan -> Lombardy
        (
            {"country": "Italy", "division": "Milan"},
            {"division": "Lombardy", "location": "Milan"},
        ),
        # Rename Trento
        # NVM, can't have that separator in the name...
        # ({'country': 'Italy', 'division': 'PA Trento'}, {'division': 'Trentino-Alto Adige/Südtirol'}),
        # Palermo -> Silicy
        (
            {"country": "Italy", "division": "Palermo"},
            {"division": "Sicily", "location": "Palermo"},
        ),
        # Rome -> Lazio
        (
            {"country": "Italy", "division": "Rome"},
            {"division": "Lazio", "location": "Rome"},
        ),
        # NETHERLANDS
        # -----------
        # Move municipalities into provinces
        # Drenthe
        (
            {"country": "Netherlands", "division": ["Coevorden", "Dalen"]},
            {"division": "Drenthe", "location": "Coevorden"},
        ),
        ({"country": "Netherlands", "division": "Drente"}, {"division": "Drenthe"}),
        # Flevoland
        (
            {"country": "Netherlands", "division": "Zeewolde"},
            {"division": "Flevoland", "location": "Zeewolde"},
        ),
        # North Brabant
        (
            {"country": "Netherlands", "division": "Andel"},
            {"division": "North Brabant", "location": "Altena"},
        ),
        (
            {"country": "Netherlands", "division": "Berlicum"},
            {"division": "North Brabant", "location": "Sint-Michielsgestel"},
        ),
        (
            {"country": "Netherlands", "division": "Eindhoven"},
            {"division": "North Brabant", "location": "Eindhoven"},
        ),
        (
            {"country": "Netherlands", "division": "Helmond"},
            {"division": "North Brabant", "location": "Helmond"},
        ),
        (
            {"country": "Netherlands", "division": "Loon op zand"},
            {"division": "North Brabant", "location": "Loop op Zand"},
        ),
        (
            {"country": "Netherlands", "division": "Milheeze"},
            {"division": "North Brabant", "location": "Milheeze"},
        ),
        (
            {"country": "Netherlands", "division": "Oisterwijk"},
            {"division": "North Brabant", "location": "Oisterwijk"},
        ),
        (
            {"country": "Netherlands", "division": "Oss"},
            {"division": "North Brabant", "location": "Oss"},
        ),
        (
            {"country": "Netherlands", "division": "Tilburg"},
            {"division": "North Brabant", "location": "Tilburg"},
        ),
        (
            {"country": "Netherlands", "division": "Noord Brabant"},
            {"division": "North Brabant"},
        ),
        # North Holland
        (
            {"country": "Netherlands", "division": "Blaricum"},
            {"division": "North Holland", "location": "Blaricum"},
        ),
        (
            {"country": "Netherlands", "division": "Diemen"},
            {"division": "North Holland", "location": "Diemen"},
        ),
        (
            {"country": "Netherlands", "division": "Haarlem"},
            {"division": "North Holland", "location": "Haarlem"},
        ),
        (
            {"country": "Netherlands", "division": "Naarden"},
            {"division": "North Holland", "location": "Gooise Meren"},
        ),
        (
            {"country": "Netherlands", "division": "Noord Holland"},
            {"division": "North Holland"},
        ),
        # South Holland
        (
            {"country": "Netherlands", "division": "Delft"},
            {"division": "South Holland", "location": "Delft"},
        ),
        (
            {"country": "Netherlands", "division": "Hardinxveld Giessendam"},
            {"division": "South Holland", "location": "Hardinxveld-Giessendam"},
        ),
        (
            {"country": "Netherlands", "division": "Leiden"},
            {"division": "South Holland", "location": "Leiden"},
        ),
        (
            {"country": "Netherlands", "division": "Nieuwendijk"},
            {"division": "South Holland", "location": "Nieuwendijk"},
        ),
        (
            {"country": "Netherlands", "division": "Nootdorp"},
            {"division": "South Holland", "location": "Nootdorp"},
        ),
        (
            {"country": "Netherlands", "division": "Rotterdam"},
            {"division": "South Holland", "location": "Rotterdam"},
        ),
        (
            {"country": "Netherlands", "division": "Zuid Holland"},
            {"division": "South Holland"},
        ),
        # Utrecht
        (
            {"country": "Netherlands", "division": "Houten"},
            {"division": "Utrecht", "location": "Houten"},
        ),
        # POLAND
        # ------
        # Fix typos
        (
            {"country": "Poland", "division": "Dolnoslakie"},
            {"division": "Dolnoslaskie"},
        ),
        # Don't use anglicized names here
        (
            {"country": "Poland", "division": ["Pomorze", "Pomerania"]},
            {"division": "Pomorskie"},
        ),
        ({"country": "Poland", "division": "Malopolska"}, {"division": "Malopolskie"}),
        (
            {"country": "Poland", "division": "Wielkopolska"},
            {"division": "Wielkopolskie"},
        ),
        # Zielonogorskie -> Lubusz
        ({"country": "Poland", "division": "Zielonogorskie"}, {"division": "Lubusz"}),
        # Fix typos
        ({"country": "Poland", "location": "Krakow"}, {"location": "Kraków"}),
        # ROMANIA
        # -------
        # Fix weird encoding issue
        ({"country": ["ÄéRomania", "\u200eRomania"]}, {"country": "Romania"}),
        # Anglicize
        ({"country": "Romania", "division": "Bucuresti"}, {"division": "Bucharest"}),
        # RUSSIA
        # ------
        # Fix typos
        ({"country": "Russia", "division": "Moscow"}, {"division": "Moscow Region"}),
        (
            {
                "country": "Russia",
                "division": ["Saint-Petersburg", "St. Petersburg", "St.Petersburg"],
            },
            {"division": "Saint Petersburg"},
        ),
        # SPAIN
        # -----
        # Move Andalusia from Sweden to Spain
        ({"country": "Sweden", "division": "Andalusia"}, {"country": "Spain"}),
        # Fix typos, unify provinces/locations
        (
            {"country": "Spain", "division": ["BasqueCountry", "Basque_Country"]},
            {"division": "Basque Country"},
        ),
        (
            {"country": "Spain", "division": "Castilla La Mancha"},
            {"division": "Castilla-La Mancha"},
        ),
        (
            {"country": "Spain", "division": "Castilla y Leon"},
            {"division": "Castilla y León"},
        ),
        ({"country": "Spain", "division": "Catalunya"}, {"division": "Catalonia"}),
        (
            {"country": "Spain", "division": "Comunitat_Valenciana"},
            {"division": "Comunitat Valenciana"},
        ),
        (
            {"country": "Spain", "location": "Bonrepos_i_Mirambell"},
            {"location": "Bonrepos i Mirambell"},
        ),
        (
            {"country": "Spain", "location": "Canet_d'En_Berenguer"},
            {"location": "Canet d'En Berenguer"},
        ),
        ({"country": "Spain", "location": "El_Puig"}, {"location": "El Puig"}),
        (
            {"country": "Spain", "location": "Grau_de_Sagunt"},
            {"location": "Grau de Sagunt"},
        ),
        (
            {"country": "Spain", "location": "Palma_de_Gandia"},
            {"location": "Palma de Gandia"},
        ),
        (
            {"country": "Spain", "location": "Tavernes_de_la_Valldigna"},
            {"location": "Tavernes de la Valldigna"},
        ),
        ({"country": "Spain", "location": "Valencia_h"}, {"location": "Valencia"}),
        ({"country": "Spain", "location": "Malaga_h"}, {"location": "Malaga"}),
        (
            {"country": "Spain", "division": ["LaRioja", "La_Rioja"]},
            {"division": "La Rioja"},
        ),
        # Barcelona -> Catalonia
        ({"country": "Spain", "division": "Barcelona"}, {"division": "Catalonia"}),
        # Fix more typos
        (
            {"country": "Spain", "location": "Alhaurin_de_la_Torre"},
            {"location": "Alhaurin de la Torre"},
        ),
        (
            {"country": "Spain", "location": "Jerez_de_la_Frontera"},
            {"location": "Jerez de la Frontera"},
        ),
        ({"country": "Spain", "location": "La_Linea"}, {"location": "La Linea"}),
        (
            {"country": "Spain", "location": "Mairena_Aljarafe"},
            {"location": "Mairena Aljarafe"},
        ),
        (
            {"country": "Spain", "location": "Puerto_Santa_Maria"},
            {"location": "Puerto Santa Maria"},
        ),
        (
            {"country": "Spain", "location": "Rincon_de_la_Victoria"},
            {"location": "Rincon de la Victoria"},
        ),
        (
            {"country": "Spain", "location": "San_Fernando"},
            {"location": "San Fernando"},
        ),
        ({"country": "Spain", "location": "San_Roque"}, {"location": "San Roque"}),
        ({"country": "Spain", "location": "Vitoria_h"}, {"location": "Vitoria"}),
        (
            {"country": "Spain", "location": "Donostia-San_Sebastian"},
            {"location": "Donostia-San Sebastian"},
        ),
        (
            {"country": "Spain", "location": "Simat_de_la_Valldigna"},
            {"location": "Simat de la Valldigna"},
        ),
        (
            {"country": "Spain", "location": "Torres_de_Elorz"},
            {"location": "Torres de Elorz"},
        ),
        (
            {"country": "Spain", "location": "Velez_Malaga"},
            {"location": "Velez Malaga"},
        ),
        # More typos
        (
            {"country": "Spain", "division": "Balear_Islands"},
            {"division": "Balearic Islands"},
        ),
        (
            {"country": "Spain", "division": "Canary_Islands"},
            {"division": "Canary Islands"},
        ),
        (
            {"country": "Spain", "location": "Rio_San_Pedro"},
            {"location": "Rio San Pedro"},
        ),
        ({"country": "Spain", "location": "Las_Palmas"}, {"location": "Las Palmas"}),
        ({"country": "Spain", "location": "La_Estrada"}, {"location": "La Estrada"}),
        ({"country": "Spain", "location": "Santa_Comba"}, {"location": "Santa Comba"}),
        (
            {"country": "Spain", "location": "Santiago_de_Compostela"},
            {"location": "Santiago de Compostela"},
        ),
        ({"country": "Spain", "location": "LeganA(c)s"}, {"location": "Leganés"}),
        ({"country": "Spain", "location": "Santa_Fe"}, {"location": "Santa Fe"}),
        ({"country": "Spain", "location": "Campo_Real"}, {"location": "Campo Real"}),
        # More typos
        (
            {"country": "Spain", "location": "Palma_de_Mallorca_h"},
            {"location": "Palma de Mallorca"},
        ),
        (
            {"country": "Spain", "location": "Santa_Maria_del_Camino"},
            {"location": "Santa Maria del Camino"},
        ),
        ({"country": "Spain", "location": "Son_Servera"}, {"location": "Son Servera"}),
        ({"country": "Spain", "location": "Logrono_h"}, {"location": "Logrono"}),
        # More typos
        ({"country": "Spain", "location": "Barcelona_h"}, {"location": "Barcelona"}),
        (
            {"country": "Spain", "location": "Alcala_de_Henares"},
            {"location": "Alcala de Henares"},
        ),
        (
            {"country": "Spain", "location": "Arganda_del_Rey"},
            {"location": "Arganda del Rey"},
        ),
        # Move Donostia-San Sebastian into Basque Country
        (
            {
                "country": "Spain",
                "division": ["Donostia-San Sebatian", "Donostia-San Sebastian"],
            },
            {
                "country": "Spain",
                "division": "Basque Country",
                "location": "Donostia-San Sebastián",
            },
        ),
        # SWEDEN
        # ------
        # Fix typos
        (
            {"country": "Sweden", "division": "Vasterbotten"},
            {"division": "Västerbotten"},
        ),
        (
            {"country": "Sweden", "division": "Gavleborgs lan"},
            {"division": "Gavleborg"},
        ),
        ({"country": "Sweden", "division": "Orebro lan"}, {"division": "Orebro"}),
        (
            {"country": "Sweden", "division": "Jamtland Harjedalen"},
            {"division": "Jamtland"},
        ),
        # SWITZERLAND
        # -----------
        # Fix typos
        (
            {
                "country": "Switzerland",
                "division": ["Basel", "Basel Stadt", "Basel Land"],
            },
            {"division": "Basel-Stadt"},
        ),
        (
            {"country": "Switzerland", "division": ["Genève", "Geneve"]},
            {"division": "Geneva"},
        ),
        ({"country": "Switzerland", "division": "Luzern"}, {"division": "Lucerne"}),
        ({"country": "Switzerland", "division": "Argovie"}, {"division": "Aargau"}),
        ({"country": "Switzerland", "division": "Zurich"}, {"division": "Zürich"}),
        (
            {"country": "Switzerland", "division": "Graubunden"},
            {"division": "Graubünden"},
        ),
        (
            {"country": "Switzerland", "division": "Sankt Gallen"},
            {"division": "St Gallen"},
        ),
        # UNITED KINGDOM
        # --------------
        # Fix typos
        (
            {"country": "United Kingdom", "location": "Northamtonshire"},
            {"location": "Northamptonshire"},
        ),
        # Move England into UK
        ({"country": "England"}, {"country": "United Kingdom", "division": "England"}),
        # CANADA
        # ------
        # Unabbreviate province names
        ({"country": "Canada", "division": "MB"}, {"division": "Manitoba"}),
        ({"country": "Canada", "division": "NB"}, {"division": "New Brunswick"}),
        (
            {"country": "Canada", "division": "NL"},
            {"division": "Newfoundland and Labrador"},
        ),
        ({"country": "Canada", "division": "NS"}, {"division": "Nova Scotia"}),
        ({"country": "Canada", "division": "SK"}, {"division": "Saskatchewan"}),
        # Toronto -> Ontario
        (
            {"country": "Canada", "division": "Toronto"},
            {"division": "Ontario", "location": "Toronto"},
        ),
        # MEXICO
        # ------
        # Unabbreviate province names
        ({"country": "Mexico", "division": "CDMX"}, {"division": "Mexico City"}),
        ({"country": "Mexico", "division": "BC"}, {"division": "Baja California"}),
        ({"country": "Mexico", "division": "BCS"}, {"division": "Baja California Sur"}),
        # North America
        # -------------
        # Who misspelled north??
        ({"region": ["Noth America", "North america"]}, {"region": "North America"}),
        # USA
        # ---
        # Merge with "United States"
        ({"country": "United States"}, {"country": "USA"}),
        # Merge with USA region
        ({"region": "USA"}, {"region": "North America", "country": "USA"}),
        # Washington DC
        # -------------
        # Unify with "District of Columbia"
        (
            {"country": "USA", "division": "District of Columbia"},
            {"division": "Washington DC"},
        ),
        # Arizona
        # -------
        # Fix typos
        # Move Phoenix into Maricopa County
        (
            {
                "country": "USA",
                "division": "Arizona",
                "location": ["Maricopa county", "Phoenix"],
            },
            {"location": "Maricopa County"},
        ),
        # California
        # ----------
        # Fix typos
        ({"country": "USA", "division": "Califonia"}, {"division": "California"}),
        # Unify county names
        (
            {
                "country": "USA",
                "division": "California",
                "location": ["Grand Princess", "Grand Princess cruise ship"],
            },
            {"location": "Grand Princess Cruise Ship"},
        ),
        (
            {"country": "USA", "division": "California", "location": "San Diego"},
            {"location": "San Diego County"},
        ),
        (
            {"country": "USA", "division": "California", "location": "San Francisco"},
            {"location": "San Francisco County"},
        ),
        # Davis -> Yolo County
        (
            {"country": "USA", "division": "California", "location": "Davis"},
            {"location": "Yolo County"},
        ),
        # LA -> LA County
        (
            {
                "country": "USA",
                "division": "California",
                "location": ["Los Angeles", "Los Angeles county"],
            },
            {"location": "Los Angeles County"},
        ),
        # Los Angeles division -> California
        (
            {"country": "USA", "division": "Los Angeles"},
            {"division": "California", "location": "Los Angeles County"},
        ),
        (
            {"country": "USA", "division": "California", "location": "Alameda"},
            {"location": "Alameda County"},
        ),
        (
            {"country": "USA", "division": "California", "location": "Solano"},
            {"location": "Solano County"},
        ),
        # Move San Diego division into California
        (
            {"country": "USA", "division": "San Diego"},
            {"division": "California", "location": "San Diego County"},
        ),
        # Imperial -> Imperial County
        (
            {"country": "USA", "division": "California", "location": "Imperial"},
            {"location": "Imperial County"},
        ),
        # Fix typos
        (
            {"country": "USA", "division": "California", "location": ["Marin county"]},
            {"location": "Marin County"},
        ),
        (
            {
                "country": "USA",
                "division": "California",
                "location": ["San Diego county"],
            },
            {"location": "San Diego County"},
        ),
        (
            {
                "country": "USA",
                "division": "California",
                "location": ["Santa Clara county"],
            },
            {"location": "Santa Clara County"},
        ),
        # Colorado
        # --------
        # Move Colorado Springs from Wisconsin to Colorado
        (
            {"country": "USA", "division": "Wisconsin", "location": "Coloardo Springs"},
            {"division": "Colorado", "location": "Colorado Springs"},
        ),
        # Colorado Springs -> El Paso County
        (
            {"country": "USA", "division": "Colorado", "location": "Colorado Springs"},
            {"location": "El Paso County"},
        ),
        # Fix typos
        (
            {"country": "USA", "division": "Colorado", "location": ["Denver county"]},
            {"location": "Denver County"},
        ),
        # Connecticut
        # -----------
        # CT --> Connecticut
        ({"country": "USA", "division": "CT"}, {"division": "Connecticut"}),
        # Move towns to counties
        (
            {
                "country": "USA",
                "division": "Connecticut",
                "location": [
                    "Fairfield county",
                    "Fairfield",
                    "Greenwich",
                    "NEWTOWN",
                    "Newtown",
                    "Bridgeport",
                    "New Fairfield",
                    "Port Chester",
                    "Riverside",
                    "Shelton",
                    "Stratford",
                    "Trumbull",
                ],
            },
            {"location": "Fairfield County"},
        ),
        (
            {
                "country": "USA",
                "division": "Connecticut",
                "location": ["Bristol", "Glastonbury", "Manchester"],
            },
            {"location": "Hartford County"},
        ),
        # ({'country': 'USA', 'division': 'Connecticut', 'location': []}, {'location': 'Litchfield County'}),
        # ({'country': 'USA', 'division': 'Connecticut', 'location': []}, {'location': 'Middlesex County'}),
        (
            {
                "country": "USA",
                "division": "Connecticut",
                "location": [
                    "Beacon Falls",
                    "BRANFORD",
                    "Branford",
                    "North-Branford",
                    "Bethany",
                    "Cheshire",
                    "East-Haven",
                    "East Haven",
                    "HAMDEN",
                    "Hamden",
                    "MADISON",
                    "Madison",
                    "MERIDEN",
                    "Meriden",
                    "Milford",
                    "NORTHFORD",
                    "Northford",
                    "New-Haven",
                    "New Haven",
                    "North-Haven",
                    "North Haven",
                    "Prospect",
                    "SEYMOUR",
                    "Seymour",
                    "Seymor",
                    "WALLINGFORD",
                    "Wallingford",
                    "Waterbury",
                    "WOODBRIDGE",
                    "Woodbridge",
                    "West-Haven",
                    "West Haven",
                ],
            },
            {"location": "New Haven County"},
        ),
        # ({'country': 'USA', 'division': 'Connecticut', 'location': []}, {'location': 'New London County'}),
        (
            {"country": "USA", "division": "Connecticut", "location": ["Mount-Vernon"]},
            {"location": "Tolland County"},
        ),
        # ({'country': 'USA', 'division': 'Connecticut', 'location': []}, {'location': 'Windham County'}),
        # Florida
        # -------
        # Fix typos
        (
            {
                "country": "USA",
                "division": "Florida",
                "location": ["Miami", "Miami-dade county"],
            },
            {"location": "Miami-Dade County"},
        ),
        (
            {
                "country": "USA",
                "division": "Florida",
                "location": ["Palm Beach county"],
            },
            {"location": "Palm Beach County"},
        ),
        # Illinois
        # --------
        # Fix typos
        # Move Chicago -> Cook County
        (
            {
                "country": "USA",
                "division": "Illinois",
                "location": ["Chicago", "Cook county"],
            },
            {"location": "Cook County"},
        ),
        (
            {"country": "USA", "division": "Illinois", "location": ["Dupage county"],},
            {"location": "Dupage County"},
        ),
        (
            {"country": "USA", "division": "Illinois", "location": ["Kane county"],},
            {"location": "Kane County"},
        ),
        (
            {"country": "USA", "division": "Illinois", "location": ["Lake county"],},
            {"location": "Lake County"},
        ),
        (
            {
                "country": "USA",
                "division": "Illinois",
                "location": ["Mchenry county", "McHenry county"],
            },
            {"location": "McHenry County"},
        ),
        (
            {
                "country": "USA",
                "division": "Illinois",
                "location": ["Winnebago county"],
            },
            {"location": "Winnebago County"},
        ),
        # Iowa
        # ----
        # Fix typos
        (
            {"country": "USA", "division": "Iowa", "location": ["Jackson county"],},
            {"location": "Jackson County"},
        ),
        # Louisiana
        # ---------
        # Assume LA is Louisiana
        ({"country": "USA", "division": "LA"}, {"division": "Louisiana"}),
        # New Jersey
        # ----------
        # Fix typos
        (
            {"country": "USA", "division": "New Jersey", "location": "Hudson"},
            {"location": "Hudson County"},
        ),
        # Massachusetts
        # -------------
        # Fix typos
        (
            {
                "country": "USA",
                "division": "Massachusetts",
                "location": ["Middlesex county"],
            },
            {"location": "Middlesex County"},
        ),
        (
            {
                "country": "USA",
                "division": "Massachusetts",
                "location": ["Norfolk county"],
            },
            {"location": "Norfolk County"},
        ),
        # Minnesota
        # ---------
        # Fix typos
        (
            {
                "country": "USA",
                "division": "Minnesota",
                "location": ["Stearns county"],
            },
            {"location": "Stearns County"},
        ),
        # New Jersey
        # ----------
        # Fix typos
        (
            {
                "country": "USA",
                "division": "New Jersey",
                "location": ["Bergen county"],
            },
            {"location": "Bergen County"},
        ),
        (
            {"country": "USA", "division": "New Jersey", "location": ["Essex county"],},
            {"location": "Essex County"},
        ),
        # New York
        # --------
        # Assume NY is NY State
        ({"country": "USA", "division": "NY"}, {"division": "New York"}),
        # Move NYC division into NY State division
        (
            {"country": "USA", "division": "New York City"},
            {"division": "New York", "location": "New York City"},
        ),
        # Remove redundant "New York"
        (
            {"country": "USA", "division": "New York", "location": "New York"},
            {"location": -1},
        ),
        # Merge NYC boroughs into one NYC
        (
            {
                "country": "USA",
                "division": "New York",
                "location": [
                    "Bronx",
                    "Brooklyn",
                    "Manhattan",
                    "Queens",
                    "Staten Island",
                ],
            },
            {"location": "New York City"},
        ),
        # Unify county names
        (
            {
                "country": "USA",
                "division": "New York",
                "location": ["Nassau", "Nassau county"],
            },
            {"location": "Nassau County"},
        ),
        (
            {"country": "USA", "division": "New York", "location": "Rockland"},
            {"location": "Rockland County"},
        ),
        (
            {
                "country": "USA",
                "division": "New York",
                "location": ["Suffolk", "Suffolk county"],
            },
            {"location": "Suffolk County"},
        ),
        (
            {"country": "USA", "division": "New York", "location": "Westchester"},
            {"location": "Westchester County"},
        ),
        # Remove empty location
        ({"country": "USA", "division": "New York", "location": ""}, {"location": -1}),
        # Move Nassau County into NY
        (
            {"country": "USA", "division": "Nassau County"},
            {"division": "New York", "location": "Nassau County"},
        ),
        # Move Yorktown Heights from CT into NY
        (
            {"country": "USA", "location": "Yorktown Heights"},
            {"division": "New York", "location": "Westchester County"},
        ),
        # Move towns into counties
        (
            {"country": "USA", "division": "New York", "location": "New Rochelle"},
            {"location": "Westchester County"},
        ),
        # Fix typos
        (
            {
                "country": "USA",
                "division": "New York",
                "location": ["Onondaga county"],
            },
            {"location": "Onondaga County"},
        ),
        (
            {
                "country": "USA",
                "division": "New York",
                "location": ["Rockland county"],
            },
            {"location": "Rockland County"},
        ),
        (
            {
                "country": "USA",
                "division": "New York",
                "location": ["Westchester county"],
            },
            {"location": "Westchester County"},
        ),
        # Pennsylvania
        # ------------
        # Fix typos
        (
            {
                "country": "USA",
                "division": "Pennsylvania",
                "location": ["Montgomery county"],
            },
            {"location": "Montgomery County"},
        ),
        # Pittsburgh -> Alleghany County
        (
            {"country": "USA", "division": "Pennsylvania", "location": ["Pittsburgh"],},
            {"location": "Alleghany County"},
        ),
        # South Carolina
        # --------------
        # Fix typos
        (
            {"country": "USA", "division": "South Carolina", "location": ["Richland"],},
            {"location": "Richland County"},
        ),
        (
            {
                "country": "USA",
                "division": "South Carolina",
                "location": ["Fairfield"],
            },
            {"location": "Fairfield County"},
        ),
        (
            {"country": "USA", "division": "South Carolina", "location": ["Berkeley"],},
            {"location": "Berkeley County"},
        ),
        # Vermont
        # -------
        # VT -> Vermont
        ({"country": "USA", "division": "VT"}, {"division": "Vermont"}),
        # Washington
        # ----------
        # Move towns into counties
        (
            {
                "country": "USA",
                "division": "Washington",
                "location": ["Kirkland", "Seattle"],
            },
            {"location": "King County"},
        ),
        (
            {"country": "USA", "division": "Washington", "location": "Tacoma"},
            {"location": "Pierce County"},
        ),
        # Remove Unknown County
        (
            {"country": "USA", "division": "Washington", "location": "Unknown County"},
            {"location": -1},
        ),
        # Wisconsin
        # ---------
        # Unify county names
        (
            {"country": "USA", "division": "Wisconsin", "location": "Campbellsp"},
            {"location": "Campbellsport"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Jackson"},
            {"location": "Jackson County"},
        ),
        # Move towns into counties
        # If a town straddles a county line, move into the one its in more, population-wise
        (
            {"country": "USA", "division": "Wisconsin", "location": "Bayside"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Belleville"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Beloit"},
            {"location": "Rock County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Blanchardville"},
            {"location": "Lafayette County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Brookfield"},
            {"location": "Waukesha County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Brooklyn"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Brown Deer"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Cambridge"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Campbellsport"},
            {"location": "Fond du Lac County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Chippewa Falls"},
            {"location": "Chippewa County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Columbus"},
            {"location": "Columbia County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Cottage Grove"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Cross Plains"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Cudahy"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "DeForest"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Elm Grove"},
            {"location": "Waukesha County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Fitchburg"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Fitchburg"},
            {"location": "Dane County"},
        ),
        # Assume this Franklin is in milwaukee
        (
            {"country": "USA", "division": "Wisconsin", "location": "Franklin"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Glendale"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Grafton"},
            {"location": "Ozaukee County"},
        ),
        # Just dump this all in Milwaukee
        (
            {
                "country": "USA",
                "division": "Wisconsin",
                "location": "Greater Milwaukee Area",
            },
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Greenfield"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Hartford"},
            {"location": "Washington County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Hillpoint"},
            {"location": "Sauk County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Holy Cross"},
            {"location": "Ozaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Janesville"},
            {"location": "Rock County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Madison"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Marshall"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Mequon"},
            {"location": "Ozaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Middleton"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Milwaukee"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Monona"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Mount Horeb"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Muscoda"},
            {"location": "Grant County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "New Berlin"},
            {"location": "Waukesha County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Oak Creek"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Oregon"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Pewaukee"},
            {"location": "Waukesha County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Port Washi"},
            {"location": "Ozaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Racine"},
            {"location": "Racine County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Richland Center"},
            {"location": "Richland County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "River Hills"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Saukville"},
            {"location": "Ozaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Slinger"},
            {"location": "Washington County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "South Milwaukee"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Spring Green"},
            {"location": "Sauk County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Stoughton"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Sun Prarie"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Thiensvill"},
            {"location": "Ozaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Tomah"},
            {"location": "Monroe County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Verona"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Waunakee"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Wauwatosa"},
            {"location": "Milwaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Whitefish"},
            {"location": "Milwaukee County"},
        ),
        # Fix typos
        (
            {"country": "USA", "division": "Wisconsin", "location": "Ozaukee county"},
            {"location": "Ozaukee County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Rock county"},
            {"location": "Rock County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Lafayette county"},
            {"location": "Lafayette County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Columbia county"},
            {"location": "Columbia County"},
        ),
        (
            {
                "country": "USA",
                "division": "Wisconsin",
                "location": "Fond du Lac county",
            },
            {"location": "Fond du Lac County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Jefferson county"},
            {"location": "Jefferson County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Milwaukee county"},
            {"location": "Milwaukee County"},
        ),
        (
            {
                "country": "USA",
                "division": "Wisconsin",
                "location": "Washington county",
            },
            {"location": "Washington County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Waukesha county"},
            {"location": "Waukesha County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Brown county"},
            {"location": "Brown County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Dane county"},
            {"location": "Dane County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Dodge county"},
            {"location": "Dodge County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Door county"},
            {"location": "Door County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Jackson county"},
            {"location": "Jackson County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Manitowoc county"},
            {"location": "Manitowoc County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Outagamie county"},
            {"location": "Outagamie County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Sauk county"},
            {"location": "Sauk County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Sheboygan county"},
            {"location": "Sheboygan County"},
        ),
        (
            {
                "country": "USA",
                "division": "Wisconsin",
                "location": "Trempealeau county",
            },
            {"location": "Trempealeau County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Vernon county"},
            {"location": "Vernon County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Waupaca county"},
            {"location": "Waupaca County"},
        ),
        (
            {"country": "USA", "division": "Wisconsin", "location": "Winnebago county"},
            {"location": "Winnebago County"},
        ),
        # Australia
        # ---------
        # Fix typos, unabbreviate province names
        ({"country": "Australia", "division": "NSW"}, {"division": "New South Wales"}),
        (
            {"country": "Australia", "division": "Northern territory"},
            {"division": "Northern Territory"},
        ),
        # New Zealand
        # -----------
        # Fix typos
        (
            {"country": "New Zealand", "division": "Combined Wellington"},
            {"division": "Wellington"},
        ),
        (
            {"country": "New Zealand", "division": "Counties Manukau"},
            {"division": "Auckland"},
        ),
        # Argentina
        # ---------
        # Fix Typos
        (
            {"country": "Argentina", "division": "Ciudad Autonoma de Buenos Aires"},
            {"division": "Buenos Aires"},
        ),
        # Brazil
        # ------
        # Fix typos
        (
            {"country": "Brazil", "division": "Minas gerais"},
            {"division": "Minas Gerais"},
        ),
        ({"country": "Brazil", "division": "São Paulo"}, {"division": "Sao Paulo"}),
        # Remove redundant Sao Paolo location
        ({"country": "Brazil", "location": "Sao Paulo"}, {"location": -1}),
        # Fix more typos
        (
            {"country": "Brazil", "location": "DUQUE DE CAXIAS"},
            {"location": "Duque de Caxias"},
        ),
        ({"country": "Brazil", "location": "NITEROI"}, {"location": "Niteroi"}),
        ({"country": "Brazil", "location": "NOVA IGUACU"}, {"location": "Nova Iguacu"}),
        ({"country": "Brazil", "location": "PETROPOLIS"}, {"location": "Petropolis"}),
        ({"country": "Brazil", "location": "QUEIMADOS"}, {"location": "Queimados"}),
        (
            {"country": "Brazil", "location": "RIO DE JANEIRO"},
            {"location": "Rio de Janeiro"},
        ),
        # More typos
        ({"country": "Brazil", "division": "Amazonas"}, {"division": "Amazonas State"}),
        # Fix typos
        ({"country": "Brazil", "division": "Ceara"}, {"division": "Ceará"}),
        ({"country": "Brazil", "division": "Goiais"}, {"division": "Goiás"}),
        ({"country": "Brazil", "division": "Maranhao"}, {"division": "Maranhão"}),
        ({"country": "Brazil", "division": "Parana"}, {"division": "Paraná"}),
        # Colombia
        # --------
        # Fix typos
        (
            {"country": "Colombia", "division": "Norte de santander"},
            {"division": "Norte de Santander"},
        ),
        ({"country": "Colombia", "location": "Cúcuta"}, {"location": "Cucuta"}),
        ({"country": "Colombia", "location": "El cerrito"}, {"location": "El Cerrito"}),
        # Clean up, move cities to departments
        (
            {"country": "Colombia", "division": "Armenia"},
            {"division": "Quindio", "location": "Armenia"},
        ),
        (
            {"country": "Colombia", "division": "Barrancabermeja"},
            {"division": "Santander", "location": "Barrancabermeja"},
        ),
        (
            {"country": "Colombia", "division": "Barranquilla"},
            {"division": "Atlantico", "location": "Barranquilla"},
        ),
        (
            {"country": "Colombia", "division": "Bello"},
            {"division": "Antioquia", "location": "Bello"},
        ),
        (
            {"country": "Colombia", "division": "Bucaramanga"},
            {"division": "Santander", "location": "Bucaramanga"},
        ),
        (
            {"country": "Colombia", "division": "Cali"},
            {"division": "Valle del Cauca", "location": "Cali"},
        ),
        (
            {"country": "Colombia", "division": "Cartagena"},
            {"division": "Bolivar", "location": "Cartagena"},
        ),
        (
            {"country": "Colombia", "division": "Cartago"},
            {"division": "Valle del Cauca", "location": "Cartago"},
        ),
        (
            {"country": "Colombia", "division": "Cienaga"},
            {"division": "Magdalena", "location": "Cienaga"},
        ),
        (
            {"country": "Colombia", "division": "Cucuta"},
            {"division": "Norte de Santander", "location": "Cucuta"},
        ),
        (
            {"country": "Colombia", "division": "Ibague"},
            {"division": "Tolima", "location": "Ibague"},
        ),
        (
            {"country": "Colombia", "division": "Leticia"},
            {"division": "Amazonas", "location": "Leticia"},
        ),
        (
            {"country": "Colombia", "division": "Manizales"},
            {"division": "Caldas", "location": "Manizales"},
        ),
        (
            {"country": "Colombia", "division": "Medellin"},
            {"division": "Antioquia", "location": "Medellin"},
        ),
        (
            {"country": "Colombia", "division": "Palmira"},
            {"division": "Valle del Cauca", "location": "Palmira"},
        ),
        (
            {"country": "Colombia", "division": "Pereira"},
            {"division": "Risaralda", "location": "Pereira"},
        ),
        (
            {"country": "Colombia", "division": "Popayan"},
            {"division": "Cauca", "location": "Popayan"},
        ),
        (
            {"country": "Colombia", "division": "Santa Marta"},
            {"division": "Magdalena", "location": "Santa Marta"},
        ),
        (
            {"country": "Colombia", "division": "Tumaco"},
            {"division": "Narino", "location": "Tumaco"},
        ),
        (
            {"country": "Colombia", "division": "Villavicencio"},
            {"division": "Meta", "location": "Villavicencio"},
        ),
        # Central America
        # ---------------
        # Move to North America
        ({"region": "Central America"}, {"region": "North America"}),
    ]

    for rule in rules:
        # print(rule)
        input_rule = rule[0]
        output_rule = rule[1]

        # Get matching entries for the input rule
        # by creating a logical mask
        # Start out with matching everything
        loc_mask = pd.Series(np.repeat(True, len(location_df)))
        for key in input_rule.keys():
            vals = input_rule[key]
            # Make it a list if it's just a single value
            if type(vals) is not list:
                vals = [vals]

            # Turn each value into a logical mask
            vals = [location_df[key] == v for v in vals]
            # Combine logical masks with logical ORs, and merge into the master mask with AND
            loc_mask = loc_mask & reduce(lambda x, y: (x | y), vals)

        # Set the output rules on the matching entries from loc_mask
        for out_key in output_rule.keys():
            location_df.loc[loc_mask, out_key] = output_rule[out_key]

    # Done
    return location_df


def build_select_tree(
    location_df, location_map_df, emoji_map_file="static_data/country_to_emoji.xls"
):
    """Build tree for ReactDropdownTreeSelect

    data
    Type: Object or Array

    Data for rendering the tree select items. The object requires the following structure:

    {
    label,          // required: Checkbox label
    value,          // required: Checkbox value
    children,       // optional: Array of child objects
    checked,        // optional: Initial state of checkbox. if true, checkbox is selected and corresponding pill is rendered.
    disabled,       // optional: Selectable state of checkbox. if true, the checkbox is disabled and the node is not selectable.
    expanded,       // optional: If true, the node is expanded (children of children nodes are not expanded by default unless children nodes also have expanded: true).
    className,      // optional: Additional css class for the node. This is helpful to style the nodes your way
    tagClassName,   // optional: Css class for the corresponding tag. Use this to add custom style the pill corresponding to the node.
    actions,        // optional: An array of extra action on the node (such as displaying an info icon or any custom icons/elements)
    dataset,        // optional: Allows data-* attributes to be set on the node and tag elements
    isDefaultValue, // optional: Indicate if a node is a default value. When true, the dropdown will automatically select the node(s) when there is no other selected node. Can be used on more than one node.
    ...             // optional: Any extra properties that you'd like to receive during `onChange` event
    }
    The action object requires the following structure:

    {
    className, // required: CSS class for the node. e.g. `fa fa-info`
    title,     // optional: HTML tooltip text
    text,      // optional: Any text to be displayed. This is helpful to pass ligatures if you're using ligature fonts
    ...        // optional: Any extra properties that you'd like to receive during `onChange` event
    }
    An array renders a tree with multiple root level items whereas an object renders a tree with a single root element (e.g. a Select All root node).


    Example:
    const data = {
    label: 'search me',
    value: 'searchme',
    children: [
        {
        label: 'search me too',
        value: 'searchmetoo',
        children: [
            {
            label: 'No one can get me',
            value: 'anonymous',
            },
        ],
        },
    ],
    }
    """

    # Set unspecified locations to None so that they don't get
    # caught up in the groupby
    location_df.loc[location_df["region"] == "-1", "region"] = None
    location_df.loc[location_df["country"] == "-1", "country"] = None
    location_df.loc[location_df["division"] == "-1", "division"] = None
    location_df.loc[location_df["location"] == "-1", "location"] = None

    # Count sequences per grouping level
    region_counts = dict(location_df.groupby("region")["location_id"].count())
    country_counts = dict(
        location_df.groupby(["region", "country"])["location_id"].count()
    )
    division_counts = dict(
        location_df.groupby(["region", "country", "division"])["location_id"].count()
    )
    location_counts = dict(
        location_df.groupby(["region", "country", "division", "location"])[
            "location_id"
        ].count()
    )

    # Load country -> emoji map
    emoji_map = pd.read_excel(emoji_map_file, skiprows=1)
    # Expand country aliases, remove whitespace from each alias
    emoji_map["aliases"] = (
        emoji_map["aliases"].str.split(",").apply(lambda x: [y.strip() for y in x])
    )

    # Root node
    select_tree = {"label": "All", "value": "All", "children": []}

    for i, loc in location_map_df.iterrows():
        # Add region node
        if loc["region"] == "-1":
            continue

        region_node = [
            c for c in select_tree["children"] if c["value"] == loc["region"]
        ]
        if region_node:
            region_node = region_node[0]
        else:
            region_node = {
                "label": loc["region"],
                "value": loc["region"],
                "level": "region",
                "location_id": loc["index"],
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(region_counts[loc["region"]]) + " sequences",
                        "text": human_format(region_counts[loc["region"]]),
                    }
                ],
                "children": [],
            }
            select_tree["children"].append(region_node)

        # Add country --> region
        if loc["country"] == "-1":
            continue

        country_node = [
            c for c in region_node["children"] if c["value"] == loc["country"]
        ]
        if country_node:
            country_node = country_node[0]
        else:

            # Look for an emoji for this country
            country_emoji = ""
            emoji_entry = emoji_map.loc[
                emoji_map["aliases"].apply(lambda x: loc["country"] in x), :
            ]
            # Fill the country emoji, if it exists
            if len(emoji_entry) == 1:
                country_emoji = emoji_entry.iat[0, 1] + " "

            country_node = {
                "label": country_emoji + loc["country"],
                "value": loc["country"],
                "region": loc["region"],
                "level": "country",
                "location_id": loc["index"],
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(country_counts[(loc["region"], loc["country"])])
                        + " sequences",
                        "text": human_format(
                            country_counts[(loc["region"], loc["country"])]
                        ),
                    }
                ],
                "children": [],
            }
            region_node["children"].append(country_node)

        # Add division --> country
        if loc["division"] == "-1":
            continue

        division_node = [
            c for c in country_node["children"] if c["value"] == loc["division"]
        ]
        if division_node:
            division_node = division_node[0]
        else:
            division_node = {
                "label": loc["division"],
                "value": loc["division"],
                "region": loc["region"],
                "country": loc["country"],
                "level": "division",
                "location_id": loc["index"],
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(
                            division_counts[
                                (loc["region"], loc["country"], loc["division"])
                            ]
                        )
                        + " sequences",
                        "text": human_format(
                            division_counts[
                                (loc["region"], loc["country"], loc["division"])
                            ]
                        ),
                    }
                ],
                "children": [],
            }
            country_node["children"].append(division_node)

        # Add location --> division
        if loc["location"] == "-1":
            continue

        location_node = [
            c for c in division_node["children"] if c["value"] == loc["location"]
        ]
        if location_node:
            location_node = location_node[0]
        else:
            location_node = {
                "label": loc["location"],
                "value": loc["location"],
                "region": loc["region"],
                "country": loc["country"],
                "division": loc["division"],
                "level": "location",
                "location_id": loc["index"],
                "actions": [
                    {
                        "className": "fa fa-info",
                        "title": str(
                            location_counts[
                                (
                                    loc["region"],
                                    loc["country"],
                                    loc["division"],
                                    loc["location"],
                                )
                            ]
                        )
                        + " sequences",
                        "text": human_format(
                            location_counts[
                                (
                                    loc["region"],
                                    loc["country"],
                                    loc["division"],
                                    loc["location"],
                                )
                            ]
                        ),
                    }
                ],
                "children": [],
            }
            division_node["children"].append(location_node)

    # print(loc_tree.nodes)
    # print(loc_tree.in_edges('New York'))
    return select_tree


# if __name__ == "__main__":
#     patient_meta_df = load_patient_metadata()
#     location_df, location_map_df = process_location_metadata(patient_meta_df)
