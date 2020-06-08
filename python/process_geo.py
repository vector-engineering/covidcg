#!/usr/bin/env python3
# coding: utf-8

'''Process clade information
Author: Albert Chen (Deverman Lab, Broad Institute)
'''

import json
import numpy as np
import networkx as nx
import pandas as pd
import re

from functools import reduce
from pathlib import Path

project_root_path = Path(__file__).resolve().parent.parent
data_dir = (project_root_path / 'data').resolve() # Resolve any symlinks --> absolute path

def load_geo_data():
    '''Load location data
    '''

    patient_meta_files = sorted((data_dir / 'patient_meta').glob('*.tsv'))
    print('Collecting {} patient metadata files...'.format(len(patient_meta_files)), end='', flush=True)
    patient_meta_df = pd.DataFrame()
    for f in patient_meta_files:
        _df = pd.read_csv(f, sep='\t', skiprows=2)
        patient_meta_df = pd.concat([patient_meta_df, _df], ignore_index=True)

    # Save dataframe
    patient_meta_df.to_csv(data_dir / 'patient_meta.csv', index=False)
    print('done', flush=True)

    # Location data is stored in one column, "region / country / division / location"
    location_df = (
        patient_meta_df['Location'].str.split('/', expand=True)
        .iloc[:, :4] # Only take 4 columns
        # Rename columns
        .rename(columns={0: 'region', 1: 'country', 2: 'division', 3: 'location'})
        .applymap(lambda x: x.strip() if x else x)
        # Placeholder for missing values, so that it will still 
        # be caught by groupby() later on
        .fillna(-1)
    )
    # Re-add metadata columns
    location_df['name'] = patient_meta_df['Virus name']
    location_df['gisaid_id'] = patient_meta_df['Accession ID']
    location_df['sample_date'] = patient_meta_df['Collection date']

    # Convert sample_date to datetime
    location_df['sample_date'] = pd.to_datetime(location_df['sample_date'], yearfirst=True)

    # Clean location data
    location_df = clean_location_data(location_df)

    # Create complete location column from the separate parts
    # This time with no padding
    location_df['loc_str'] = location_df['region'].str.cat(
        [location_df['country'].astype(str), location_df['division'].astype(str), location_df['location'].astype(str)],
        sep='/'
    )

    unique_location_df = (
        location_df['loc_str']
        .drop_duplicates()
        .sort_values(ignore_index=True)
        .reset_index() # Produce an index column
    )
    # Location data is stored in one column, "region/country/division/location"
    unique_location_df = pd.concat([unique_location_df, (
        unique_location_df['loc_str'].str.split('/', expand=True)
        .iloc[:, :4] # Only take 4 columns
        # Rename columns
        .rename(columns={0: 'region', 1: 'country', 2: 'division', 3: 'location'})
    )], axis=1)

    # Map location IDs back to taxon_locations dataframe
    location_df['location_id'] = location_df['loc_str'].map(pd.Series(
        unique_location_df['index'].values, index=unique_location_df['loc_str'].values
    ))

    # Take subset of columns, re-index
    location_df = (
        location_df
        [['name', 'gisaid_id', 'sample_date', 'location_id', 'region', 'country', 'division', 'location']]
        .sort_values('location_id')
        .reset_index(drop=True)
    )

    # print(location_df)

    print('Saving taxon locations')
    location_df.to_csv(data_dir / 'taxon_locations.csv', index=False)
    # location_df.to_json(data_dir / 'taxon_locations.json', orient='records')

    print('Saving unique locations')
    unique_location_df.drop(columns=['loc_str']).to_csv(data_dir / 'location_map.csv', index=False)
    unique_location_df.drop(columns=['loc_str']).to_json(data_dir / 'location_map.json', orient='records')

    return location_df, unique_location_df


def clean_location_data(location_df):
    '''Fix typos, unify nomenclature in location data
    '''

    # Make rules, which are mappings of input region/country/division/location -> 
    # output region/country/division/location
    # Each rule is a tuple, where the first element is the input rule
    # (which entries to match)
    # and the second element is the output rule
    # (what fields to change)

    rules = [
        # SOUTH AFRICA
        # ------------

        # Unabbreviate province names
        ({'country': 'South Africa', 'division': 'EC'}, {'division': 'Eastern Cape'}),
        ({'country': 'South Africa', 'division': 'KZN'}, {'division': 'KwaZulu-Natal'}),
        ({'country': 'South Africa', 'division': 'GP'}, {'division': 'Gauteng'}),
        ({'country': 'South Africa', 'division': 'LP'}, {'division': 'Limpopo'}),
        ({'country': 'South Africa', 'division': 'MP'}, {'division': 'Mpumalanga'}),
        # Remove Unknown division
        ({'country': 'South Africa', 'division': 'Unknown'}, {'division': -1}),


        # CHINA
        # -----

        # Move NanChang into Jiangxi province
        (
            {'country': 'China', 'division': 'NanChang'}, 
            {'division': 'Jiangxi', 'location': 'Nanchang'}
        ),
        # Move Guangzhou into Guangdong province
        (
            {'country': 'China', 'division': 'Guangzhou'}, 
            {'division': 'Guangdong', 'location': 'Guangzhou'}
        ),
        # Move Meizhou into Guangdong province
        (
            {'country': 'China', 'division': 'Meizhou'}, 
            {'division': 'Guangdong', 'location': 'Meizhou'}
        ),
        # Move Hangzhou into Zhejiang province
        (
            {'country': 'China', 'division': 'Hangzhou'}, 
            {'division': 'Zhejiang', 'location': 'Hangzhou'}
        ),
        # Move Wuhan into Hubei province
        (
            {'country': 'China', 'division': 'Wuhan'}, 
            {'division': 'Hubei', 'location': 'Wuhan'}
        ),

        # INDIA
        # -----

        # Fix typos
        ({'country': 'India', 'division': 'West_Bengal'}, {'division': 'West Bengal'}),
        ({'country': 'India', 'division': 'Telengana'}, {'division': 'Telangana'}),
        ({'country': 'India', 'division': 'Rajsathan'}, {'division': 'Rajasthan'}),
        # Jammu --> Jammu and Kashmir
        (
            {'country': 'India', 'division': 'Jammu'}, 
            {'division': 'Jammu and Kashmir', 'location': 'Jammu'}
        ),
        # Kargil --> Ladakh/Kargil
        (
            {'country': 'India', 'division': 'Kargil'}, 
            {'division': 'Ladakh', 'location': 'Kargil'}
        ),
        # Mumbai --> Maharashtra/Mumbai
        (
            {'country': 'India', 'division': 'Mumbai'}, 
            {'division': 'Maharashtra', 'location': 'Mumbai'}
        ),

        # ISRAEL
        # ------
        # South Coast District --> South District
        ({'country': 'Israel', 'division': 'South Coast District'}, {'division': 'South District'}),

        # JAPAN
        # -----
        # Remove unknown division
        ({'country': 'Japan', 'division': 'unknown'}, {'division': -1}),

        # OMAN
        # ----
        # Fix typos
        ({'country': 'Oman', 'division': ['Muscta', 'Musca']}, {'division': 'Muscat'}),

        # PAKISTAN
        # --------
        # Unabbreviate province names
        ({'country': 'Pakistan', 'division': 'KPK'}, {'division': 'Khyber Pakhtunkhwa'}),

        # SOUTH KOREA
        # -----------
        # Korea --> South Korea
        # I assume North Korea is not submitting genomes...
        ({'country': 'Korea'}, {'country': 'South Korea'}),

        # TAIWAN
        # ------
        # Fix typos
        ({'country': 'Taiwan', 'division': 'New Taipei city'}, {'division': 'New Taipei City'}),

        # THAILAND
        # --------
        # Fix typos
        ({'country': 'Thailand', 'division': 'Phatum thani'}, {'division': 'Pathum Thani'}),

        # VIETNAM
        # -------
        # Fix typos
        ({'country': 'Vietnam', 'division': 'Quangning'}, {'division': 'Quangninh'}),


        # BELGIUM
        # -------
        # Move Belgian cities into provinces
        # Since it is getting out of hand. the list is too big
        # And merge towns into parent municipalities

        # Antwerp
        ({'country': 'Belgium', 'division': 'Antwerp'}, {'division': 'Antwerp', 'location': 'Antwerp'}),
        ({'country': 'Belgium', 'division': 'Bonheiden'}, {'division': 'Antwerp', 'location': 'Bonheiden'}),
        ({'country': 'Belgium', 'division': ['Gierle', 'Lille']}, {'division': 'Antwerp', 'location': 'Lille'}),
        ({'country': 'Belgium', 'division': 'Hemisem'}, {'division': 'Antwerp', 'location': 'Hemiksem'}),
        ({'country': 'Belgium', 'division': 'Herselt'}, {'division': 'Antwerp', 'location': 'Herselt'}),
        ({'country': 'Belgium', 'division': 'Kalmthout'}, {'division': 'Antwerp', 'location': 'Kalmthout'}),
        ({'country': 'Belgium', 'division': 'Kasterlee'}, {'division': 'Antwerp', 'location': 'Kasterlee'}),
        ({'country': 'Belgium', 'division': 'Oevel'}, {'division': 'Antwerp', 'location': 'Westerlo'}),
        ({'country': 'Belgium', 'division': 'Poppel'}, {'division': 'Antwerp', 'location': 'Ravels'}),
        ({'country': 'Belgium', 'division': 'Schoten'}, {'division': 'Antwerp', 'location': 'Schoten'}),
        ({'country': 'Belgium', 'division': 'Schoten'}, {'division': 'Antwerp', 'location': 'Schoten'}),

        # East Flanders
        ({'country': 'Belgium', 'division': 'Deinze'}, {'division': 'East Flanders', 'location': 'Deinze'}),
        ({'country': 'Belgium', 'division': 'Dendermonde'}, {'division': 'East Flanders', 'location': 'Dendermonde'}),
        ({'country': 'Belgium', 'division': ['Gent', 'Ghent']}, {'division': 'East Flanders', 'location': 'Ghent'}),
        ({'country': 'Belgium', 'division': 'Geraardsbergen'}, {'division': 'East Flanders', 'location': 'Geraardsbergen'}),
        ({'country': 'Belgium', 'division': 'Merelbeke'}, {'division': 'East Flanders', 'location': 'Merelbeke'}),
        ({'country': 'Belgium', 'division': 'Nevele'}, {'division': 'East Flanders', 'location': 'Nevele'}),
        ({'country': 'Belgium', 'division': 'Ronse'}, {'division': 'East Flanders', 'location': 'Ronse'}),
        ({'country': 'Belgium', 'division': 'Sint-Gillis-Waas'}, {'division': 'East Flanders', 'location': 'Sint-Gillis-Waas'}),
        ({'country': 'Belgium', 'division': 'Sint-Niklaas'}, {'division': 'East Flanders', 'location': 'Sint-Niklaas'}),
        ({'country': 'Belgium', 'division': 'Stekene'}, {'division': 'East Flanders', 'location': 'Stekene'}),
        ({'country': 'Belgium', 'division': 'Waasmunster'}, {'division': 'East Flanders', 'location': 'Waasmunster'}),

        # Flemish Brabant
        ({'country': 'Belgium', 'division': 'Asse'}, {'division': 'Flemish Brabant', 'location': 'Asse'}),
        ({'country': 'Belgium', 'division': 'Beersel'}, {'division': 'Flemish Brabant', 'location': 'Beersel'}),
        ({'country': 'Belgium', 'division': 'Bierbeek'}, {'division': 'Flemish Brabant', 'location': 'Bierbeek'}),
        ({'country': 'Belgium', 'division': 'Boutersem'}, {'division': 'Flemish Brabant', 'location': 'Boutersem'}),
        ({'country': 'Belgium', 'division': 'Dilbeek'}, {'division': 'Flemish Brabant', 'location': 'Dilbeek'}),
        ({'country': 'Belgium', 'division': 'Grimbergen'}, {'division': 'Flemish Brabant', 'location': 'Grimbergen'}),
        ({'country': 'Belgium', 'division': 'Hoegaarden'}, {'division': 'Flemish Brabant', 'location': 'Hoegaarden'}),
        ({'country': 'Belgium', 'division': 'Holsbeek'}, {'division': 'Flemish Brabant', 'location': 'Holsbeek'}),
        ({'country': 'Belgium', 'division': 'Huldenberg'}, {'division': 'Flemish Brabant', 'location': 'Huldenberg'}),
        ({'country': 'Belgium', 'division': 'Kraainem'}, {'division': 'Flemish Brabant', 'location': 'Kraainem'}),
        ({'country': 'Belgium', 'division': ['Leuven', 'Heverlee', 'Kessel-Lo', 'Ladeuze']}, {'division': 'Flemish Brabant', 'location': 'Leuven'}),
        ({'country': 'Belgium', 'division': 'Linter'}, {'division': 'Flemish Brabant', 'location': 'Linter'}),
        ({'country': 'Belgium', 'division': 'Lubbeek'}, {'division': 'Flemish Brabant', 'location': 'Lubbeek'}),
        ({'country': 'Belgium', 'division': 'Sint-Genesius-Rode'}, {'division': 'Flemish Brabant', 'location': 'Sint-Genesius-Rode'}),
        ({'country': 'Belgium', 'division': 'Tervuren'}, {'division': 'Flemish Brabant', 'location': 'Tervuren'}),
        ({'country': 'Belgium', 'division': 'Tielt-Winge'}, {'division': 'Flemish Brabant', 'location': 'Tielt-Winge'}),
        ({'country': 'Belgium', 'division': 'Tienen'}, {'division': 'Flemish Brabant', 'location': 'Tienen'}),
        ({'country': 'Belgium', 'division': 'Vilvoorde'}, {'division': 'Flemish Brabant', 'location': 'Vilvoorde'}),
        ({'country': 'Belgium', 'division': 'Winksele'}, {'division': 'Flemish Brabant', 'location': 'Herent'}),
        ({'country': 'Belgium', 'division': 'Zoutleeuw'}, {'division': 'Flemish Brabant', 'location': 'Zoutleeuw'}),

        # Limburg
        ({'country': 'Belgium', 'division': 'Alken'}, {'division': 'Limburg', 'location': 'Alken'}),
        ({'country': 'Belgium', 'division': ['Beringen', 'Koersel']}, {'division': 'Limburg', 'location': 'Beringen'}),
        ({'country': 'Belgium', 'division': 'Bree'}, {'division': 'Limburg', 'location': 'Bree'}),
        ({'country': 'Belgium', 'division': 'Hasselt'}, {'division': 'Limburg', 'location': 'Hasselt'}),
        ({'country': 'Belgium', 'division': 'Houthalen-Helchteren'}, {'division': 'Limburg', 'location': 'Houthalen-Helchteren'}),
        ({'country': 'Belgium', 'division': 'Kaulile'}, {'division': 'Limburg', 'location': 'Kaulille'}),
        ({'country': 'Belgium', 'division': 'Lanaken'}, {'division': 'Limburg', 'location': 'Lanaken'}),
        ({'country': 'Belgium', 'division': 'Lummen'}, {'division': 'Limburg', 'location': 'Lummen'}),
        ({'country': 'Belgium', 'division': ['Neerpelt', 'Overpelt']}, {'division': 'Limburg', 'location': 'Pelt'}),
        ({'country': 'Belgium', 'division': 'Nieuwerkerken'}, {'division': 'Limburg', 'location': 'Nieuwerkerken'}),
        ({'country': 'Belgium', 'division': 'Riemst'}, {'division': 'Limburg', 'location': 'Riemst'}),
        ({'country': 'Belgium', 'division': 'Sint-Truiden'}, {'division': 'Limburg', 'location': 'Sint-Truiden'}),
        ({'country': 'Belgium', 'division': 'Zolder'}, {'division': 'Limburg', 'location': 'Heusden-Zolder'}),
        ({'country': 'Belgium', 'division': 'Zonhoven'}, {'division': 'Limburg', 'location': 'Zonhoven'}),

        # West Flanders
        ({'country': 'Belgium', 'division': 'Geluwe'}, {'division': 'West Flanders', 'location': 'Wervik'}),
        ({'country': 'Belgium', 'division': 'Knokke-Heist'}, {'division': 'West Flanders', 'location': 'Knokke-Heist'}),
        ({'country': 'Belgium', 'division': 'Oostrozebeke'}, {'division': 'West Flanders', 'location': 'Oostrozebeke'}),
        ({'country': 'Belgium', 'division': 'Wevelgem'}, {'division': 'West Flanders', 'location': 'Wevelgem'}),

        # Hainaut
        ({'country': 'Belgium', 'division': ['Ath', 'Arbre']}, {'division': 'Hainaut', 'location': 'Ath'}),
        ({'country': 'Belgium', 'division': 'Bassilly'}, {'division': 'Hainaut', 'location': 'Silly'}),
        ({'country': 'Belgium', 'division': ['Blandain', 'Doornik', 'Kain', 'Tournai']}, {'division': 'Hainaut', 'location': 'Tournai'}),
        ({'country': 'Belgium', 'division': 'Boussu'}, {'division': 'Hainaut', 'location': 'Boussu'}),
        ({'country': 'Belgium', 'division': ['Charleroi', 'Montignies-Sur-Sambre', 'Couillet', 'Jumet']}, {'division': 'Hainaut', 'location': 'Charleroi'}),
        ({'country': 'Belgium', 'division': ['Colfontaine', 'Confontaine']}, {'division': 'Hainaut', 'location': 'Colfontaine'}),
        ({'country': 'Belgium', 'division': 'Cuesmes'}, {'division': 'Hainaut', 'location': 'Cuesmes'}),
        ({'country': 'Belgium', 'division': 'Ellezelles'}, {'division': 'Hainaut', 'location': 'Ellezelles'}),
        ({'country': 'Belgium', 'division': ['Fontaine-l Eveque', 'Forchies-la-marche']}, {'division': 'Hainaut', 'location': 'Fontaine-l\'Évêque'}),
        ({'country': 'Belgium', 'division': 'Frameries'}, {'division': 'Hainaut', 'location': 'Frameries'}),
        ({'country': 'Belgium', 'division': 'Haulchin'}, {'division': 'Hainaut', 'location': 'Haulchin'}),
        ({'country': 'Belgium', 'division': 'Honnelles'}, {'division': 'Hainaut', 'location': 'Honnelles'}),
        ({'country': 'Belgium', 'division': 'Komen'}, {'division': 'Hainaut', 'location': 'Comines-Warneton'}),
        ({'country': 'Belgium', 'division': ['La Louvriere', 'La Louviere']}, {'division': 'Hainaut', 'location': 'La Louvière'}),
        ({'country': 'Belgium', 'division': ['Mons', 'Courcelles', 'Havre', 'Jemappes', 'Trazegnies']}, {'division': 'Hainaut', 'location': 'Mons'}),
        ({'country': 'Belgium', 'division': 'Montigny-le-Tilleul'}, {'division': 'Hainaut', 'location': 'Montigny-le-Tilleul'}),
        ({'country': 'Belgium', 'division': 'Mouscron'}, {'division': 'Hainaut', 'location': 'Mouscron'}),
        ({'country': 'Belgium', 'division': 'Quaregnon'}, {'division': 'Hainaut', 'location': 'Quaregnon'}),
        ({'country': 'Belgium', 'division': 'Quevaucamps'}, {'division': 'Hainaut', 'location': 'Belœil'}),
        ({'country': 'Belgium', 'division': 'Quevy'}, {'division': 'Hainaut', 'location': 'Quevy'}),
        ({'country': 'Belgium', 'division': ['Saint-Ghislain', 'Tertre', 'Villerot']}, {'division': 'Hainaut', 'location': 'Saint-Ghislain'}),
        ({'country': 'Belgium', 'division': ['Saisinne', 'Soignies', 'Zinnik']}, {'division': 'Hainaut', 'location': 'Soignies'}),

        # Liege
        ({'country': 'Belgium', 'division': 'Berloz'}, {'division': 'Liege', 'location': 'Berloz'}),
        ({'country': 'Belgium', 'division': 'Bierset'}, {'division': 'Liege', 'location': 'Grâce-Hollogne'}),
        ({'country': 'Belgium', 'division': 'Couthuin'}, {'division': 'Liege', 'location': 'Heron'}),
        ({'country': 'Belgium', 'division': 'Dalhem'}, {'division': 'Liege', 'location': 'Dalhem'}),
        ({'country': 'Belgium', 'division': 'Eupen'}, {'division': 'Liege', 'location': 'Eupen'}),
        ({'country': 'Belgium', 'division': 'Henri-chapelle'}, {'division': 'Liege', 'location': 'Welkenraedt'}),
        ({'country': 'Belgium', 'division': 'Herstal'}, {'division': 'Liege', 'location': 'Herstal'}),
        ({'country': 'Belgium', 'division': ['Herve', 'Chaineux']}, {'division': 'Liege', 'location': 'Herve'}),
        ({'country': 'Belgium', 'division': ['Liege', 'Liège']}, {'division': 'Liege', 'location': 'Liège'}),
        ({'country': 'Belgium', 'division': 'Limbourg'}, {'division': 'Liege', 'location': 'Limbourg'}),
        ({'country': 'Belgium', 'division': 'Saive'}, {'division': 'Liege', 'location': 'Blegny'}),
        ({'country': 'Belgium', 'division': 'Seraing'}, {'division': 'Liege', 'location': 'Seraing'}),
        ({'country': 'Belgium', 'division': 'Soumagne'}, {'division': 'Liege', 'location': 'Soumagne'}),
        ({'country': 'Belgium', 'division': 'Theux'}, {'division': 'Liege', 'location': 'Theux'}),
        ({'country': 'Belgium', 'division': 'Wandre'}, {'division': 'Liege', 'location': 'Wandre'}),
        ({'country': 'Belgium', 'division': 'Wanze'}, {'division': 'Liege', 'location': 'Wanze'}),
        ({'country': 'Belgium', 'division': 'Yernee-fraineux'}, {'division': 'Liege', 'location': 'Nandrin'}),

        # Luxembourg
        # ...

        # Namur
        ({'country': 'Belgium', 'division': 'Andenne'}, {'division': 'Namur', 'location': 'Andenne'}),
        ({'country': 'Belgium', 'division': 'Couvin'}, {'division': 'Namur', 'location': 'Couvin'}),
        ({'country': 'Belgium', 'division': 'Havelange'}, {'division': 'Namur', 'location': 'Havelange'}),
        ({'country': 'Belgium', 'division': 'Temploux'}, {'division': 'Namur', 'location': 'Namur'}),

        # Walloon Brabant
        ({'country': 'Belgium', 'division': 'Braine-l alleud'}, {'division': 'Walloon Brabant', 'location': 'Braine-l\'Alleud'}),
        ({'country': 'Belgium', 'division': 'Braine-le-Chateau'}, {'division': 'Walloon Brabant', 'location': 'Braine-le-Château'}),
        ({'country': 'Belgium', 'division': 'Genappe'}, {'division': 'Walloon Brabant', 'location': 'Genappe'}),
        ({'country': 'Belgium', 'division': 'Lasne'}, {'division': 'Walloon Brabant', 'location': 'Lasne'}),
        ({'country': 'Belgium', 'division': 'Rixensart'}, {'division': 'Walloon Brabant', 'location': 'Rixensart'}),
        ({'country': 'Belgium', 'division': 'Tubize'}, {'division': 'Walloon Brabant', 'location': 'Tubize'}),
        ({'country': 'Belgium', 'division': 'Waterloo'}, {'division': 'Walloon Brabant', 'location': 'Waterloo'}),

        # Brussels-Capital Region
        ({'country': 'Belgium', 'division': 'Anderlecht'}, {'division': 'Brussels-Capital Region', 'location': 'Anderlecht'}),
        ({'country': 'Belgium', 'division': ['Brussel', 'Brussels']}, {'division': 'Brussels-Capital Region', 'location': 'Brussels'}),
        ({'country': 'Belgium', 'division': 'Elsene'}, {'division': 'Brussels-Capital Region', 'location': 'Ixelles'}),
        ({'country': 'Belgium', 'division': 'Etterbeek'}, {'division': 'Brussels-Capital Region', 'location': 'Etterbeek'}),
        ({'country': 'Belgium', 'division': 'Ganshoren'}, {'division': 'Brussels-Capital Region', 'location': 'Ganshoren'}),
        ({'country': 'Belgium', 'division': 'Oudergem'}, {'division': 'Brussels-Capital Region', 'location': 'Auderghem'}),
        ({'country': 'Belgium', 'division': 'Sint-Agatha-Berchem'}, {'division': 'Brussels-Capital Region', 'location': 'Sint-Agatha-Berchem'}),
        ({'country': 'Belgium', 'division': ['Sint-Pieter-Woluwe', 'Sint-Pieters-Woluwe']}, {'division': 'Brussels-Capital Region', 'location': 'Woluwe-Saint-Pierre'}),
        ({'country': 'Belgium', 'division': 'Ukkel'}, {'division': 'Brussels-Capital Region', 'location': 'Uccle'}),
        ({'country': 'Belgium', 'division': 'Vorst'}, {'division': 'Brussels-Capital Region', 'location': 'Forest'}),
        ({'country': 'Belgium', 'division': 'Watermaal-Bosvoorde'}, {'division': 'Brussels-Capital Region', 'location': 'Watermael-Boitsfort'}),

        # CROATIA
        # -------
        # Fix typos
        ({'country': 'Croatia', 'division': 'Varazdin'}, {'division': 'Varaždin County'}),

        # CZECHIA
        # -------
        # Rename from "Czech Republic" --> "Czechia"
        ({'country': 'Czech Republic'}, {'country': 'Czechia'}),
        # Move "Prague" --> "Central Czhechia and Prague"
        ({'country': 'Czechia', 'division': 'Prague'}, {'division': 'Central Czechia and Prague'}),

        # DENMARK
        # -------
        # Remove Unknown region
        ({'country': 'Denmark', 'division': 'Unknown'}, {'division': -1}),

        # FRANCE
        # ------
        
        # Move Saint-Saulve from Belgium to France
        ({'country': 'Belgium', 'division': 'Saint-Saulve'}, {'country': 'France', 'division': 'Hauts-de-France'}),
        # Fix typos
        (
            {'country': 'France', 'division': ['Bourgogne Franche comte', 'Bourgogne-France-Comté', 'Bourgogne-Franche Comte', 'Bourgogne']}, 
            {'division': 'Bourgogne-Franche-Comté'}
        ),
        # Unabbreviate province names, unify province names
        ({'country': 'France', 'division': 'ARA'}, {'division': 'Auvergne-Rhône-Alpes'}),
        ({'country': 'France', 'division': 'Centre - Val de Loire'}, {'division': 'Centre-Val de Loire'}),
        ({'country': 'France', 'division': ['Grand-Est', 'Grand-est']}, {'division': 'Grand Est'}),
        ({'country': 'France', 'division': ['Hauts De France', 'Hauts de France']}, {'division': 'Hauts-de-France'}),
        ({'country': 'France', 'division': ['IDF', 'Ile De France', 'Ile de France', 'Ile-de-France']}, {'division': 'Île-de-France'}),
        # Move Toulouse to Occitanie
        ({'country': 'France', 'division': 'Toulouse'}, {'division': 'Occitanie', 'location': 'Toulouse'}),


        # GEORGIA
        # -------
        # Move to Asia
        ({'region': 'Europe', 'country': 'Georgia'}, {'region': 'Asia'}),

        # GERMANY
        # -------
        
        # Move Munich division into Bavaria division
        ({'country': 'Germany', 'division': 'Munich'}, {'division': 'Bavaria', 'location': 'Munich'}),
        # Duesseldorf -> North Rhine Westphalia
        ({'country': 'Germany', 'division': 'Duesseldorf'}, {'division': 'North Rhine Westphalia', 'location': 'Duesseldorf'}),
        # Frankfurt -> Hesse
        ({'country': 'Germany', 'division': 'Frankfurt'}, {'division': 'Hesse', 'location': 'Frankfurt'}),
        # Rostock -> Mecklenburg-Vorpommern
        ({'country': 'Germany', 'division': 'Rostock'}, {'division': 'Mecklenburg-Vorpommern', 'location': 'Rostock'}),

        # ITALY
        # -----
        
        # Fix typos
        ({'country': 'Italy', 'division': ['Lombarida', 'Lombarida']}, {'division': 'Lombardy'}),
        # Castel di Sangro -> Abruzzo
        ({'country': 'Italy', 'division': 'Castel di Sangro'}, {'division': 'Abruzzo', 'location': 'Castel di Sangro'}),
        # Teramo -> Abruzzo
        ({'country': 'Italy', 'division': 'Teramo'}, {'division': 'Abruzzo', 'location': 'Teramo'}),
        # Cagliari -> Sardinia
        ({'country': 'Italy', 'division': 'Cagliari'}, {'division': 'Sardinia', 'location': 'Cagliari'}),
        # Milan -> Lombardy
        ({'country': 'Italy', 'division': 'Milan'}, {'division': 'Lombardy', 'location': 'Milan'}),
        # Rename Trento
        # NVM, can't have that separator in the name...
        # ({'country': 'Italy', 'division': 'PA Trento'}, {'division': 'Trentino-Alto Adige/Südtirol'}),
        # Palermo -> Silicy
        ({'country': 'Italy', 'division': 'Palermo'}, {'division': 'Sicily', 'location': 'Palermo'}),
        # Rome -> Lazio
        ({'country': 'Italy', 'division': 'Rome'}, {'division': 'Lazio', 'location': 'Rome'}),

        # MOROCCO
        # -------
        # Move to Africa
        # Move Cadiz to Andalusia, Spain

        ({'country': 'Morocco'}, {'region': 'Africa'}),
        ({'country': 'Morocco', 'location': 'Cadiz_h'}, {'region': 'Europe', 'country': 'Spain', 'division': 'Andalusia', 'location': 'Cadiz'}),

        # NETHERLANDS
        # -----------
        # Move municipalities into provinces

        # Drenthe
        ({'country': 'Netherlands', 'division': ['Coevorden', 'Dalen']}, {'division': 'Drenthe', 'location': 'Coevorden'}),
        ({'country': 'Netherlands', 'division': 'Drente'}, {'division': 'Drenthe'}),

        # Flevoland
        ({'country': 'Netherlands', 'division': 'Zeewolde'}, {'division': 'Flevoland', 'location': 'Zeewolde'}),

        # North Brabant
        ({'country': 'Netherlands', 'division': 'Andel'}, {'division': 'North Brabant', 'location': 'Altena'}),
        ({'country': 'Netherlands', 'division': 'Berlicum'}, {'division': 'North Brabant', 'location': 'Sint-Michielsgestel'}),
        ({'country': 'Netherlands', 'division': 'Eindhoven'}, {'division': 'North Brabant', 'location': 'Eindhoven'}),
        ({'country': 'Netherlands', 'division': 'Helmond'}, {'division': 'North Brabant', 'location': 'Helmond'}),
        ({'country': 'Netherlands', 'division': 'Loon op zand'}, {'division': 'North Brabant', 'location': 'Loop op Zand'}),
        ({'country': 'Netherlands', 'division': 'Milheeze'}, {'division': 'North Brabant', 'location': 'Milheeze'}),
        ({'country': 'Netherlands', 'division': 'Oisterwijk'}, {'division': 'North Brabant', 'location': 'Oisterwijk'}),
        ({'country': 'Netherlands', 'division': 'Oss'}, {'division': 'North Brabant', 'location': 'Oss'}),
        ({'country': 'Netherlands', 'division': 'Tilburg'}, {'division': 'North Brabant', 'location': 'Tilburg'}),
        ({'country': 'Netherlands', 'division': 'Noord Brabant'}, {'division': 'North Brabant'}),

        # North Holland
        ({'country': 'Netherlands', 'division': 'Blaricum'}, {'division': 'North Holland', 'location': 'Blaricum'}),
        ({'country': 'Netherlands', 'division': 'Diemen'}, {'division': 'North Holland', 'location': 'Diemen'}),
        ({'country': 'Netherlands', 'division': 'Haarlem'}, {'division': 'North Holland', 'location': 'Haarlem'}),
        ({'country': 'Netherlands', 'division': 'Naarden'}, {'division': 'North Holland', 'location': 'Gooise Meren'}),
        ({'country': 'Netherlands', 'division': 'Noord Holland'}, {'division': 'North Holland'}),

        # South Holland
        ({'country': 'Netherlands', 'division': 'Delft'}, {'division': 'South Holland', 'location': 'Delft'}),
        ({'country': 'Netherlands', 'division': 'Hardinxveld Giessendam'}, {'division': 'South Holland', 'location': 'Hardinxveld-Giessendam'}),
        ({'country': 'Netherlands', 'division': 'Leiden'}, {'division': 'South Holland', 'location': 'Leiden'}),
        ({'country': 'Netherlands', 'division': 'Nieuwendijk'}, {'division': 'South Holland', 'location': 'Nieuwendijk'}),
        ({'country': 'Netherlands', 'division': 'Nootdorp'}, {'division': 'South Holland', 'location': 'Nootdorp'}),
        ({'country': 'Netherlands', 'division': 'Rotterdam'}, {'division': 'South Holland', 'location': 'Rotterdam'}),
        ({'country': 'Netherlands', 'division': 'Zuid Holland'}, {'division': 'South Holland'}),

        # Utrecht
        ({'country': 'Netherlands', 'division': 'Houten'}, {'division': 'Utrecht', 'location': 'Houten'}),

        # POLAND
        # ------
        # Fix typos
        # Don't use anglicized names here
        ({'country': 'Poland', 'division': 'Dolnoslakie'}, {'division': 'Dolnoslaskie'}),
        ({'country': 'Poland', 'division': ['Pomorze', 'Pomerania']}, {'division': 'Pomorskie'}),
        ({'country': 'Poland', 'division': 'Malopolska'}, {'division': 'Malopolskie'}),
        ({'country': 'Poland', 'division': 'Wielkopolska'}, {'division': 'Wielkopolskie'}),
        # Zielonogorskie -> Lubusz
        ({'country': 'Poland', 'division': 'Zielonogorskie'}, {'division': 'Lubusz'}),

        # RUSSIA
        # ------
        # Fix typos

        ({'country': 'Russia', 'division': 'Moscow'}, {'division': 'Moscow Region'}),
        ({'country': 'Russia', 'division': ['Saint-Petersburg', 'St. Petersburg', 'St.Petersburg']}, {'division': 'Saint Petersburg'}),

        # SPAIN
        # -----

        # Move Andalusia from Sweden to Spain
        ({'country': 'Sweden', 'division': 'Andalusia'}, {'country': 'Spain'}),

        # Fix typos, unify provinces/locations
        ({'country': 'Spain', 'division': ['BasqueCountry', 'Basque_Country']}, {'division': 'Basque Country'}),
        ({'country': 'Spain', 'division': 'Castilla La Mancha'}, {'division': 'Castilla-La Mancha'}),
        ({'country': 'Spain', 'division': 'Castilla y Leon'}, {'division': 'Castilla y León'}),
        ({'country': 'Spain', 'division': 'Catalunya'}, {'division': 'Catalonia'}),
        ({'country': 'Spain', 'division': 'Comunitat_Valenciana'}, {'division': 'Comunitat Valenciana'}),
        ({'country': 'Spain', 'location': 'Bonrepos_i_Mirambell'}, {'location': 'Bonrepos i Mirambell'}),
        ({'country': 'Spain', 'location': 'Canet_d\'En_Berenguer'}, {'location': 'Canet d\'En Berenguer'}),
        ({'country': 'Spain', 'location': 'El_Puig'}, {'location': 'El Puig'}),
        ({'country': 'Spain', 'location': 'Grau_de_Sagunt'}, {'location': 'Grau de Sagunt'}),
        ({'country': 'Spain', 'location': 'Palma_de_Gandia'}, {'location': 'Palma de Gandia'}),
        ({'country': 'Spain', 'location': 'Tavernes_de_la_Valldigna'}, {'location': 'Tavernes de la Valldigna'}),
        ({'country': 'Spain', 'location': 'Valencia_h'}, {'location': 'Valencia'}),
        ({'country': 'Spain', 'location': 'Malaga_h'}, {'location': 'Malaga'}),
        ({'country': 'Spain', 'division': ['LaRioja', 'La_Rioja']}, {'division': 'La Rioja'}),
        # Barcelona -> Catalonia
        ({'country': 'Spain', 'division': 'Barcelona'}, {'division': 'Catalonia'}),
        # Fix more typos
        ({'country': 'Spain', 'location': 'Alhaurin_de_la_Torre'}, {'location': 'Alhaurin de la Torre'}),
        ({'country': 'Spain', 'location': 'Jerez_de_la_Frontera'}, {'location': 'Jerez de la Frontera'}),
        ({'country': 'Spain', 'location': 'La_Linea'}, {'location': 'La Linea'}),
        ({'country': 'Spain', 'location': 'Mairena_Aljarafe'}, {'location': 'Mairena Aljarafe'}),
        ({'country': 'Spain', 'location': 'Puerto_Santa_Maria'}, {'location': 'Puerto Santa Maria'}),
        ({'country': 'Spain', 'location': 'Rincon_de_la_Victoria'}, {'location': 'Rincon de la Victoria'}),
        ({'country': 'Spain', 'location': 'San_Fernando'}, {'location': 'San Fernando'}),
        ({'country': 'Spain', 'location': 'San_Roque'}, {'location': 'San Roque'}),
        ({'country': 'Spain', 'location': 'Vitoria_h'}, {'location': 'Vitoria'}),
        ({'country': 'Spain', 'location': 'Donostia-San_Sebastian'}, {'location': 'Donostia-San Sebastian'}),
        ({'country': 'Spain', 'location': 'Simat_de_la_Valldigna'}, {'location': 'Simat de la Valldigna'}),
        ({'country': 'Spain', 'location': 'Torres_de_Elorz'}, {'location': 'Torres de Elorz'}),
        

        # SWEDEN
        # ------

        # Fix typos
        ({'country': 'Sweden', 'division': 'Vasterbotten'}, {'division': 'Västerbotten'}),
        ({'country': 'Sweden', 'division': 'Gavleborgs lan'}, {'division': 'Gavleborg'}),
        ({'country': 'Sweden', 'division': 'Orebro lan'}, {'division': 'Orebro'}),
        

        # SWITZERLAND
        # -----------
        # Fix typos

        ({'country': 'Switzerland', 'division': ['Basel', 'Basel Stadt', 'Basel Land']}, {'division': 'Basel-Stadt'}),
        ({'country': 'Switzerland', 'division': ['Genève', 'Geneve']}, {'division': 'Geneva'}),
        ({'country': 'Switzerland', 'division': 'Luzern'}, {'division': 'Lucerne'}),
        ({'country': 'Switzerland', 'division': 'Argovie'}, {'division': 'Aargau'}),

        # UNITED KINGDOM
        # --------------
        # Fix typos

        ({'country': 'United Kingdom', 'location': 'Northamtonshire'}, {'location': 'Northamptonshire'}),

        # CANADA
        # ------
        # Unabbreviate province names

        ({'country': 'Canada', 'division': 'MB'}, {'division': 'Manitoba'}),
        ({'country': 'Canada', 'division': 'NB'}, {'division': 'New Brunswick'}),
        ({'country': 'Canada', 'division': 'NL'}, {'division': 'Newfoundland and Labrador'}),
        ({'country': 'Canada', 'division': 'NS'}, {'division': 'Nova Scotia'}),
        ({'country': 'Canada', 'division': 'SK'}, {'division': 'Saskatchewan'}),

        # MEXICO
        # ------
        # Unabbreviate province names

        ({'country': 'Mexico', 'division': 'CDMX'}, {'division': 'Mexico City'}),

        # North America
        # -------------
        # Who misspelled north??

        ({'region': ['Noth America', 'North america']}, {'region': 'North America'}),

        # USA
        # ---

        # Washington DC
        # -------------
        # Unify with "District of Columbia"

        ({'country': 'USA', 'division': 'District of Columbia'}, {'division': 'Washington DC'}),

        # California
        # ----------
        # Unify county names

        ({'country': 'USA', 'division': 'California', 'location': ['Grand Princess', 'Grand Princess cruise ship']}, {'location': 'Grand Princess Cruise Ship'}),
        ({'country': 'USA', 'division': 'California', 'location': 'San Diego'}, {'location': 'San Diego County'}),
        ({'country': 'USA', 'division': 'California', 'location': 'San Francisco'}, {'location': 'San Francisco County'}),
        # Davis -> Yolo County
        ({'country': 'USA', 'division': 'California', 'location': 'Davis'}, {'location': 'Yolo County'}),

        # Colorado
        # --------
        # Move Colorado Springs from Wisconsin to Colorado
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Coloardo Springs'}, {'division': 'Colorado', 'location': 'Colorado Springs'}),

        # Louisiana
        # ---------
        # Assume LA is Louisiana

        ({'country': 'USA', 'division': 'LA'}, {'division': 'Louisiana'}),

        # New Jersey
        # ----------
        # Fix typos
        ({'country': 'USA', 'division': 'New Jersey', 'location': 'Hudson'}, {'location': 'Hudson County'}),


        # New York
        # --------
        # Assume NY is NY State
        # Unify county names
        # Move NYC division into NY State division
        # Merge NYC boroughs into one NYC

        ({'country': 'USA', 'division': 'NY'}, {'division': 'New York'}),
        
        ({'country': 'USA', 'division': 'New York City'}, {'division': 'New York', 'location': 'New York City'}),
        # Remove redundant "New York"
        ({'country': 'USA', 'division': 'New York', 'location': 'New York'}, {'location': -1}),
        ({'country': 'USA', 'division': 'New York', 'location': ['Bronx', 'Brooklyn', 'Manhattan', 'Queens', 'Staten Island']}, {'location': 'New York City'}),
        ({'country': 'USA', 'division': 'New York', 'location': ['Nassau', 'Nassau county']}, {'location': 'Nassau County'}),
        ({'country': 'USA', 'division': 'New York', 'location': 'Rockland'}, {'location': 'Rockland County'}),
        ({'country': 'USA', 'division': 'New York', 'location': ['Suffolk', 'Suffolk county']}, {'location': 'Suffolk County'}),
        ({'country': 'USA', 'division': 'New York', 'location': 'Westchester'}, {'location': 'Westchester County'}),
        # Remove empty location
        ({'country': 'USA', 'division': 'New York', 'location': ''}, {'location': -1}),


        # Washington
        # ----------
        # Move towns into counties
        ({'country': 'USA', 'division': 'Washington', 'location': ['Kirkland', 'Seattle']}, {'location': 'King County'}),
        ({'country': 'USA', 'division': 'Washington', 'location': 'Tacoma'}, {'location': 'Pierce County'}),
        # Remove Unknown County
        ({'country': 'USA', 'division': 'Washington', 'location': 'Unknown County'}, {'location': -1}),

        # Wisconsin
        # ---------
        # Unify county names
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Campbellsp'}, {'location': 'Campbellsport'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Jackson'}, {'location': 'Jackson County'}),
        # Move towns into counties
        # If a town straddles a county line, move into the one its in more, population-wise
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Bayside'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Belleville'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Beloit'}, {'location': 'Rock County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Blanchardville'}, {'location': 'Lafayette County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Brookfield'}, {'location': 'Waukesha County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Brooklyn'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Brown Deer'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Cambridge'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Campbellsport'}, {'location': 'Fond du Lac County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Chippewa Falls'}, {'location': 'Chippewa County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Columbus'}, {'location': 'Columbia County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Cottage Grove'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Cross Plains'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Cudahy'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'DeForest'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Elm Grove'}, {'location': 'Waukesha County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Fitchburg'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Fitchburg'}, {'location': 'Dane County'}),
        # Assume this Franklin is in milwaukee
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Franklin'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Glendale'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Grafton'}, {'location': 'Ozaukee County'}),
        # Just dump this all in Milwaukee
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Greater Milwaukee Area'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Greenfield'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Hartford'}, {'location': 'Washington County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Hillpoint'}, {'location': 'Sauk County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Holy Cross'}, {'location': 'Ozaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Janesville'}, {'location': 'Rock County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Madison'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Marshall'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Mequon'}, {'location': 'Ozaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Middleton'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Milwaukee'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Monona'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Mount Horeb'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Muscoda'}, {'location': 'Grant County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'New Berlin'}, {'location': 'Waukesha County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Oak Creek'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Oregon'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Pewaukee'}, {'location': 'Waukesha County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Port Washi'}, {'location': 'Ozaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Racine'}, {'location': 'Racine County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Richland Center'}, {'location': 'Richland County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'River Hills'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Saukville'}, {'location': 'Ozaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Slinger'}, {'location': 'Washington County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'South Milwaukee'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Spring Green'}, {'location': 'Sauk County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Stoughton'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Sun Prarie'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Thiensvill'}, {'location': 'Ozaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Tomah'}, {'location': 'Monroe County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Verona'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Waunakee'}, {'location': 'Dane County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Wauwatosa'}, {'location': 'Milwaukee County'}),
        ({'country': 'USA', 'division': 'Wisconsin', 'location': 'Whitefish'}, {'location': 'Milwaukee County'}),
        

        # Australia
        # ---------
        # Fix typos, unabbreviate province names

        ({'country': 'Australia', 'division': 'NSW'}, {'division': 'New South Wales'}),
        ({'country': 'Australia', 'division': 'Northern territory'}, {'division': 'Northern Territory'}),

        # New Zealand
        # -----------
        # Fix typos
        ({'country': 'New Zealand', 'division': 'Combined Wellington'}, {'division': 'Wellington'}),
        ({'country': 'New Zealand', 'division': 'Counties Manukau'}, {'division': 'Auckland'}),

        # Argentina
        # ---------
        # Fix Typos
        ({'country': 'Argentina', 'division': 'Ciudad Autonoma de Buenos Aires'}, {'division': 'Buenos Aires'}),

        # Brazil
        # ------
        # Fix typos
        # Remove redundant Sao Paolo location

        ({'country': 'Brazil', 'division': 'Minas gerais'}, {'division': 'Minas Gerais'}),
        ({'country': 'Brazil', 'division': 'São Paulo'}, {'division': 'Sao Paulo'}),
        ({'country': 'Brazil', 'location': 'Sao Paulo'}, {'location': -1}),

        # Colombia
        # --------
        # Fix typos
        ({'country': 'Colombia', 'division': 'Norte de santander'}, {'division': 'Norte de Santander'}),
        ({'country': 'Colombia', 'location': 'Cúcuta'}, {'location': 'Cucuta'}),
        ({'country': 'Colombia', 'location': 'El cerrito'}, {'location': 'El Cerrito'}),
        # Clean up, move cities to departments
        ({'country': 'Colombia', 'division': 'Armenia'}, {'division': 'Quindio', 'location': 'Armenia'}),
        ({'country': 'Colombia', 'division': 'Barrancabermeja'}, {'division': 'Santander', 'location': 'Barrancabermeja'}),
        ({'country': 'Colombia', 'division': 'Barranquilla'}, {'division': 'Atlantico', 'location': 'Barranquilla'}),
        ({'country': 'Colombia', 'division': 'Bello'}, {'division': 'Antioquia', 'location': 'Bello'}),
        ({'country': 'Colombia', 'division': 'Bucaramanga'}, {'division': 'Santander', 'location': 'Bucaramanga'}),
        ({'country': 'Colombia', 'division': 'Cali'}, {'division': 'Valle del Cauca', 'location': 'Cali'}),
        ({'country': 'Colombia', 'division': 'Cartagena'}, {'division': 'Bolivar', 'location': 'Cartagena'}),
        ({'country': 'Colombia', 'division': 'Cartago'}, {'division': 'Valle del Cauca', 'location': 'Cartago'}),
        ({'country': 'Colombia', 'division': 'Cienaga'}, {'division': 'Magdalena', 'location': 
        'Cienaga'}),
        ({'country': 'Colombia', 'division': 'Cucuta'}, {'division': 'Norte de Santander', 'location': 'Cucuta'}),
        ({'country': 'Colombia', 'division': 'Ibague'}, {'division': 'Tolima', 'location': 'Ibague'}),
        ({'country': 'Colombia', 'division': 'Leticia'}, {'division': 'Amazonas', 'location': 'Leticia'}),
        ({'country': 'Colombia', 'division': 'Manizales'}, {'division': 'Caldas', 'location': 'Manizales'}),
        ({'country': 'Colombia', 'division': 'Medellin'}, {'division': 'Antioquia', 'location': 'Medellin'}),
        ({'country': 'Colombia', 'division': 'Palmira'}, {'division': 'Valle del Cauca', 'location': 'Palmira'}),
        ({'country': 'Colombia', 'division': 'Pereira'}, {'division': 'Risaralda', 'location': 'Pereira'}),
        ({'country': 'Colombia', 'division': 'Popayan'}, {'division': 'Cauca', 'location': 'Popayan'}),
        ({'country': 'Colombia', 'division': 'Santa Marta'}, {'division': 'Magdalena', 'location': 'Santa Marta'}),
        ({'country': 'Colombia', 'division': 'Tumaco'}, {'division': 'Narino', 'location': 'Tumaco'}),
        ({'country': 'Colombia', 'division': 'Villavicencio'}, {'division': 'Meta', 'location': 'Villavicencio'}),

        # Central America
        # ---------------
        # Move to North America

        ({'region': 'Central America'}, {'region': 'North America'})
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
            loc_mask = (loc_mask & reduce(lambda x, y: (x | y), vals)) 
            
        # Set the output rules on the matching entries from loc_mask
        for out_key in output_rule.keys():
            location_df.loc[loc_mask, out_key] = output_rule[out_key]
    

    # Done
    return location_df


# Thanks to user "rtaft" from https://stackoverflow.com/questions/579310/formatting-long-numbers-as-strings-in-python
def human_format(num):
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '({}{})'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])


def build_select_tree(location_df, unique_location_df):
    '''Build tree for ReactDropdownTreeSelect

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
    '''

    # Set unspecified locations to None so that they don't get
    # caught up in the groupby
    location_df.loc[location_df['region'] == '-1', 'region'] = None
    location_df.loc[location_df['country'] == '-1', 'country'] = None
    location_df.loc[location_df['division'] == '-1', 'division'] = None
    location_df.loc[location_df['location'] == '-1', 'location'] = None

    # Count sequences per grouping level
    region_counts = dict(location_df.groupby('region')['gisaid_id'].count())
    country_counts = dict(location_df.groupby(['region', 'country'])['gisaid_id'].count())
    division_counts = dict(location_df.groupby(['region', 'country', 'division'])['gisaid_id'].count())
    location_counts = dict(location_df.groupby(['region', 'country', 'division', 'location'])['gisaid_id'].count())

    # Root node
    select_tree = {
        'label': 'All',
        'value': 'All',
        'children': []
    }

    for i, loc in unique_location_df.iterrows():
        # Add region node
        if loc['region'] == '-1':
            continue
        
        region_node = [c for c in select_tree['children'] if c['value'] == loc['region']]
        if region_node:
            region_node = region_node[0]
        else:
            region_node = {
                'label': loc['region'],
                'value': loc['region'],
                'level': 'region',
                'actions': [{
                    'className': 'fa fa-info',
                    'title': str(region_counts[loc['region']]) + ' sequences',
                    'text': human_format(region_counts[loc['region']])
                }],
                'children': []
            }
            select_tree['children'].append(region_node)

        # Add country --> region
        if loc['country'] == '-1':
            continue

        country_node = [c for c in region_node['children'] if c['value'] == loc['country']]
        if country_node:
            country_node = country_node[0]
        else:
            country_node = {
                'label': loc['country'],
                'value': loc['country'],
                'region': loc['region'],
                'level': 'country',
                'actions': [{
                    'className': 'fa fa-info',
                    'title': str(country_counts[(loc['region'], loc['country'])]) + ' sequences',
                    'text': human_format(country_counts[(loc['region'], loc['country'])])
                }],
                'children': []
            }
            region_node['children'].append(country_node)
        
        # Add division --> country
        if loc['division'] == '-1':
            continue

        division_node = [c for c in country_node['children'] if c['value'] == loc['division']]
        if division_node:
            division_node = division_node[0]
        else:
            division_node = {
                'label': loc['division'],
                'value': loc['division'],
                'region': loc['region'],
                'country': loc['country'],
                'level': 'division',
                'actions': [{
                    'className': 'fa fa-info',
                    'title': str(division_counts[(loc['region'], loc['country'], loc['division'])]) + ' sequences',
                    'text': human_format(division_counts[(loc['region'], loc['country'], loc['division'])])
                }],
                'children': []
            }
            country_node['children'].append(division_node)

        # Add location --> division
        if loc['location'] == '-1':
            continue
        
        location_node = [c for c in division_node['children'] if c['value'] == loc['location']]
        if location_node:
            location_node = location_node[0]
        else:
            location_node = {
                'label': loc['location'],
                'value': loc['location'],
                'region': loc['region'],
                'country': loc['country'],
                'division': loc['division'],
                'level': 'location',
                'actions': [{
                    'className': 'fa fa-info',
                    'title': str(location_counts[(loc['region'], loc['country'], loc['division'], loc['location'])]) + ' sequences',
                    'text': human_format(location_counts[(loc['region'], loc['country'], loc['division'], loc['location'])])
                }],
                'children': []
            }
            division_node['children'].append(location_node)


    # Save tree as json file
    print('Saving geo select tree')
    select_tree_path = data_dir / 'geo_select_tree.json'
    with select_tree_path.open('w') as fp:
        fp.write(json.dumps(select_tree))

    # print(loc_tree.nodes)
    # print(loc_tree.in_edges('New York'))
    return select_tree


def main():
    location_df, unique_location_df = load_geo_data()
    select_tree = build_select_tree(location_df, unique_location_df)
    # print(select_tree)


if __name__ == '__main__':
    main()

