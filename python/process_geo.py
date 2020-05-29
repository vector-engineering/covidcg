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
        [['name', 'gisaid_id', 'sample_date', 'location_id']]
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

    # SOUTH AFRICA
    # ------------
    # Unabbreviate province names

    location_df.loc[
        (location_df['country'] == 'South Africa') &
        (location_df['division'] == 'EC'),
        'division'
    ] = 'Eastern Cape'

    location_df.loc[
        (location_df['country'] == 'South Africa') &
        (location_df['division'] == 'KZN'),
        'division'
    ] = 'KwaZulu-Natal'

    location_df.loc[
        (location_df['country'] == 'South Africa') &
        (location_df['division'] == 'GP'),
        'division'
    ] = 'Gauteng'

    location_df.loc[
        (location_df['country'] == 'South Africa') &
        (location_df['division'] == 'LP'),
        'division'
    ] = 'Limpopo'

    location_df.loc[
        (location_df['country'] == 'South Africa') &
        (location_df['division'] == 'MP'),
        'division'
    ] = 'Mpumalanga'

    # CHINA
    # -----
    # Move NanChang into Jiangxi province
    # Move Guangzhou into Guangdong province
    # Move Meizhou into Guangdong province
    # Move Hangzhou into Zhejiang province
    # Move Wuhan into Hubei province

    location_df.loc[
        (location_df['country'] == 'China') &
        (location_df['division'] == 'NanChang'),
        ['division', 'location']
    ] = ['Jiangxi', 'Nanchang']

    location_df.loc[
        (location_df['country'] == 'China') &
        (location_df['division'] == 'Guangzhou'),
        ['division', 'location']
    ] = ['Guangdong', 'Guangzhou']

    location_df.loc[
        (location_df['country'] == 'China') &
        (location_df['division'] == 'Meizhou'),
        ['division', 'location']
    ] = ['Guangdong', 'Meizhou']

    location_df.loc[
        (location_df['country'] == 'China') &
        (location_df['division'] == 'Hangzhou'),
        ['division', 'location']
    ] = ['Zhejiang', 'Hangzhou']

    location_df.loc[
        (location_df['country'] == 'China') &
        (location_df['division'] == 'Wuhan'),
        ['division', 'location']
    ] = ['Hubei', 'Wuhan']

    # INDIA
    # -----
    # Fix typos
    # Jammu --> Jammu and Kashmir
    # Kargil --> Ladakh/Kargil
    # Mumbai --> Maharashtra/Mumbai

    location_df.loc[
        (location_df['country'] == 'India') &
        (location_df['division'] == 'West_Bengal'),
        'division'
    ] = 'West Bengal'

    location_df.loc[
        (location_df['country'] == 'India') &
        (location_df['division'] == 'Telengana'),
        'division'
    ] = 'Telangana'

    location_df.loc[
        (location_df['country'] == 'India') &
        (location_df['division'] == 'Jammu'),
        ['division', 'location']
    ] = ['Jammu and Kashmir', 'Jammu']

    location_df.loc[
        (location_df['country'] == 'India') &
        (location_df['division'] == 'Kargil'),
        ['division', 'location']
    ] = ['Ladakh', 'Kargil']

    location_df.loc[
        (location_df['country'] == 'India') &
        (location_df['division'] == 'Mumbai'),
        ['division', 'location']
    ] = ['Maharashtra', 'Mumbai']

    # ISRAEL
    # ------
    # South Coast District --> South District

    location_df.loc[
        (location_df['country'] == 'Israel') &
        (location_df['division'] == 'South Coast District'),
        'division'
    ] = 'South District'

    # PAKISTAN
    # --------
    # Unabbreviate province names

    location_df.loc[
        (location_df['country'] == 'Pakistan') &
        (location_df['division'] == 'KPK'),
        'division'
    ] = 'Khyber Pakhtunkhwa'

    # SOUTH KOREA
    # -----------
    # Korea --> South Korea
    # I assume North Korea is not submitting genomes...

    location_df.loc[
        (location_df['country'] == 'Korea'),
        'country'
    ] = 'South Korea'

    # TAIWAN
    # ------
    # Fix typos

    location_df.loc[
        (location_df['division'] == 'New Taipei city'),
        'division'
    ] = 'New Taipei City'

    # BELGIUM
    # -------
    # Move Belgian cities into provinces
    # Since it is getting out of hand. the list is too big
    # And merge towns into parent municipalities

    # Antwerp

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Antwerp'),
        ['division', 'location']
    ] = ['Antwerp', 'Antwerp']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Bonheiden'),
        ['division', 'location']
    ] = ['Antwerp', 'Bonheiden']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Gierle') |
            (location_df['division'] == 'Lille')
        ),
        ['division', 'location']
    ] = ['Antwerp', 'Lille']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Hemisem'),
        ['division', 'location']
    ] = ['Antwerp', 'Hemiksem']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Herselt'),
        ['division', 'location']
    ] = ['Antwerp', 'Herselt']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Kalmthout'),
        ['division', 'location']
    ] = ['Antwerp', 'Kalmthout']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Kasterlee'),
        ['division', 'location']
    ] = ['Antwerp', 'Kasterlee']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Oevel')
        ),
        ['division', 'location']
    ] = ['Antwerp', 'Westerlo']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Poppel')
        ),
        ['division', 'location']
    ] = ['Antwerp', 'Ravels']

    location_df.loc[
        (location_df['country'] == 'Belgium') & 
        (location_df['division'] == 'Schoten'),
        ['division', 'location']
    ] = ['Antwerp', 'Schoten']

    # East Flanders

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Deinze'),
        ['division', 'location']
    ] = ['East Flanders', 'Deinze']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Dendermonde'),
        ['division', 'location']
    ] = ['East Flanders', 'Dendermonde']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Gent') |
            (location_df['division'] == 'Ghent')
        ),
        ['division', 'location']
    ] = ['East Flanders', 'Ghent']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Geraardsbergen'),
        ['division', 'location']
    ] = ['East Flanders', 'Geraardsbergen']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Merelbeke'),
        ['division', 'location']
    ] = ['East Flanders', 'Merelbeke']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Nevele'),
        ['division', 'location']
    ] = ['East Flanders', 'Nevele']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Ronse'),
        ['division', 'location']
    ] = ['East Flanders', 'Ronse']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Sint-Gillis-Waas'),
        ['division', 'location']
    ] = ['East Flanders', 'Sint-Gillis-Waas']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Sint-Niklaas'),
        ['division', 'location']
    ] = ['East Flanders', 'Sint-Niklaas']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Stekene'),
        ['division', 'location']
    ] = ['East Flanders', 'Stekene']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Waasmunster'),
        ['division', 'location']
    ] = ['East Flanders', 'Waasmunster']

    # Flemish Brabant

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Asse'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Asse']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Beersel'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Beersel']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Bierbeek'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Bierbeek']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Boutersem'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Boutersem']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Dilbeek'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Dilbeek']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Grimbergen'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Grimbergen']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Hoegaarden'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Hoegaarden']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Holsbeek'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Holsbeek']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Huldenberg'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Huldenberg']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Kraainem'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Kraainem']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Leuven') |
            (location_df['division'] == 'Heverlee') |
            (location_df['division'] == 'Kessel-Lo') |
            (location_df['division'] == 'Ladeuze')
        ),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Leuven']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Linter'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Linter']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Lubbeek'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Lubbeek']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Sint-Genesius-Rode'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Sint-Genesius-Rode']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Tervuren'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Tervuren']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Tielt-Winge'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Tielt-Winge']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Tienen'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Tienen']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Vilvoorde'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Vilvoorde']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Winksele'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Herent']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Zoutleeuw'),
        ['division', 'location']
    ] = ['Flemish Brabant', 'Zoutleeuw']

    # Limburg

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Alken'),
        ['division', 'location']
    ] = ['Limburg', 'Alken']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Beringen') |
            (location_df['division'] == 'Koersel')
        ),
        ['division', 'location']
    ] = ['Limburg', 'Beringen']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Bree'),
        ['division', 'location']
    ] = ['Limburg', 'Bree']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Hasselt'),
        ['division', 'location']
    ] = ['Limburg', 'Hasselt']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Houthalen-Helchteren'),
        ['division', 'location']
    ] = ['Limburg', 'Houthalen-Helchteren']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Kaulile'),
        ['division', 'location']
    ] = ['Limburg', 'Kaulille']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Lanaken'),
        ['division', 'location']
    ] = ['Limburg', 'Lanaken']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Lummen'),
        ['division', 'location']
    ] = ['Limburg', 'Lummen']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Neerpelt') | 
            (location_df['division'] == 'Overpelt')
        ),
        ['division', 'location']
    ] = ['Limburg', 'Pelt']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Nieuwerkerken'),
        ['division', 'location']
    ] = ['Limburg', 'Nieuwerkerken']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Riemst'),
        ['division', 'location']
    ] = ['Limburg', 'Riemst']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Sint-Truiden'),
        ['division', 'location']
    ] = ['Limburg', 'Sint-Truiden']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Zolder'),
        ['division', 'location']
    ] = ['Limburg', 'Heusden-Zolder']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Zonhoven'),
        ['division', 'location']
    ] = ['Limburg', 'Zonhoven']

    # West Flanders

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Geluwe'),
        ['division', 'location']
    ] = ['West Flanders', 'Wervik']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Knokke-Heist'),
        ['division', 'location']
    ] = ['West Flanders', 'Knokke-Heist']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Oostrozebeke'),
        ['division', 'location']
    ] = ['West Flanders', 'Oostrozebeke']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Wevelgem'),
        ['division', 'location']
    ] = ['West Flanders', 'Wevelgem']

    # Hainaut


    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Ath') |
            (location_df['division'] == 'Arbre')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'Ath']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Bassilly'),
        ['division', 'location']
    ] = ['Hainaut', 'Silly']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Blandain') |
            (location_df['division'] == 'Doornik') |
            (location_df['division'] == 'Kain')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'Tournai']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Boussu'),
        ['division', 'location']
    ] = ['Hainaut', 'Boussu']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Charleroi') |
            (location_df['division'] == 'Montignies-Sur-Sambre') |
            (location_df['division'] == 'Couillet') |
            (location_df['division'] == 'Jumet')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'Charleroi']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Colfontaine') |
            (location_df['division'] == 'Confontaine')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'Colfontaine']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Cuesmes'),
        ['division', 'location']
    ] = ['Hainaut', 'Cuesmes']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Ellezelles'),
        ['division', 'location']
    ] = ['Hainaut', 'Ellezelles']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Fontaine-l Eveque') | 
            (location_df['division'] == 'Forchies-la-marche')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'Fontaine-l\'Évêque']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Frameries'),
        ['division', 'location']
    ] = ['Hainaut', 'Frameries']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Haulchin'),
        ['division', 'location']
    ] = ['Hainaut', 'Haulchin']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Honnelles'),
        ['division', 'location']
    ] = ['Hainaut', 'Honnelles']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Komen'),
        ['division', 'location']
    ] = ['Hainaut', 'Comines-Warneton']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'La Louvriere') | 
            (location_df['division'] == 'La Louviere')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'La Louvière']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Mons') |
            (location_df['division'] == 'Courcelles') | 
            (location_df['division'] == 'Havre') |
            (location_df['division'] == 'Jemappes') | 
            (location_df['division'] == 'Trazegnies')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'Mons']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Montigny-le-Tilleul'),
        ['division', 'location']
    ] = ['Hainaut', 'Montigny-le-Tilleul']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Mouscron'),
        ['division', 'location']
    ] = ['Hainaut', 'Mouscron']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Quaregnon'),
        ['division', 'location']
    ] = ['Hainaut', 'Quaregnon']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Quevaucamps'),
        ['division', 'location']
    ] = ['Hainaut', 'Belœil']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Quevy'),
        ['division', 'location']
    ] = ['Hainaut', 'Quevy']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Saint-Ghislain') |
            (location_df['division'] == 'Tertre') |
            (location_df['division'] == 'Villerot')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'Saint-Ghislain']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Saisinne') | 
            (location_df['division'] == 'Soignies') |
            (location_df['division'] == 'Zinnik')
        ),
        ['division', 'location']
    ] = ['Hainaut', 'Soignies']
    

    # Liege

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Berloz'),
        ['division', 'location']
    ] = ['Liege', 'Berloz']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Bierset'),
        ['division', 'location']
    ] = ['Liege', 'Grâce-Hollogne']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Couthuin'),
        ['division', 'location']
    ] = ['Liege', 'Heron']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Dalhem'),
        ['division', 'location']
    ] = ['Liege', 'Dalhem']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Eupen'),
        ['division', 'location']
    ] = ['Liege', 'Eupen']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Henri-chapelle'),
        ['division', 'location']
    ] = ['Liege', 'Welkenraedt']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Herstal'),
        ['division', 'location']
    ] = ['Liege', 'Herstal']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Herve') |
            (location_df['division'] == 'Chaineux')
        ),
        ['division', 'location']
    ] = ['Liege', 'Herve']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Liege') |
            (location_df['division'] == 'Liège')
        ),
        ['division', 'location']
    ] = ['Liege', 'Liège']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Limbourg'),
        ['division', 'location']
    ] = ['Liege', 'Limbourg']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Saive'),
        ['division', 'location']
    ] = ['Liege', 'Blegny']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Seraing'),
        ['division', 'location']
    ] = ['Liege', 'Seraing']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Soumagne'),
        ['division', 'location']
    ] = ['Liege', 'Soumagne']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Theux'),
        ['division', 'location']
    ] = ['Liege', 'Theux']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Wandre'),
        ['division', 'location']
    ] = ['Liege', 'Wandre']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Wanze'),
        ['division', 'location']
    ] = ['Liege', 'Wanze']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Yernee-fraineux'),
        ['division', 'location']
    ] = ['Liege', 'Nandrin']

    # Luxembourg
    # Namur

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Andenne'),
        ['division', 'location']
    ] = ['Namur', 'Andenne']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Couvin'),
        ['division', 'location']
    ] = ['Namur', 'Couvin']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Havelange'),
        ['division', 'location']
    ] = ['Namur', 'Havelange']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Temploux'),
        ['division', 'location']
    ] = ['Namur', 'Namur']

    # Walloon Brabant

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Braine-l alleud'),
        ['division', 'location']
    ] = ['Walloon Brabant', 'Braine-l\'Alleud']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Braine-le-Chateau'),
        ['division', 'location']
    ] = ['Walloon Brabant', 'Braine-le-Château']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Genappe'),
        ['division', 'location']
    ] = ['Walloon Brabant', 'Genappe']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Lasne'),
        ['division', 'location']
    ] = ['Walloon Brabant', 'Lasne']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Rixensart'),
        ['division', 'location']
    ] = ['Walloon Brabant', 'Rixensart']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Tubize'),
        ['division', 'location']
    ] = ['Walloon Brabant', 'Tubize']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Waterloo'),
        ['division', 'location']
    ] = ['Walloon Brabant', 'Waterloo']

    # Brussels-Capital Region

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Anderlecht'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Anderlecht']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Brussel') |
            (location_df['division'] == 'Brussels')
        ),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Brussels']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Elsene'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Ixelles']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Etterbeek'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Etterbeek']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Ganshoren'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Ganshoren']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Oudergem'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Auderghem']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Sint-Agatha-Berchem'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Sint-Agatha-Berchem']

    location_df.loc[
        (location_df['country'] == 'Belgium') & (
            (location_df['division'] == 'Sint-Pieter-Woluwe') |
            (location_df['division'] == 'Sint-Pieters-Woluwe')
        ),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Woluwe-Saint-Pierre']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Ukkel'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Uccle']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Vorst'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Forest']

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Watermaal-Bosvoorde'),
        ['division', 'location']
    ] = ['Brussels-Capital Region', 'Watermael-Boitsfort']
    
    # CZECHIA
    # -------
    # Rename from "Czech Republic" --> "Czechia"
    # Move "Prague" --> "Central Czhechia and Prague"
    
    location_df.loc[
        (location_df['country'] == 'Czech Republic'),
        'country'
    ] = 'Czechia'

    location_df.loc[
        (location_df['country'] == 'Czech Republic') &
        (location_df['division'] == 'Prague'),
        'division'
    ] = 'Central Czechia and Prague'


    # FRANCE
    # ------
    # Fix typos, unify province names
    # Unabbreviate province names
    # Move Saint-Saulve from Belgium to France

    location_df.loc[
        (location_df['country'] == 'Belgium') &
        (location_df['division'] == 'Saint-Saulve'),
        ['country', 'division']
    ] = ['France', 'Hauts-de-France']


    location_df.loc[
        (location_df['division'] == 'Bourgogne Franche comte') |
        (location_df['division'] == 'Bourgogne-France-Comté') |
        (location_df['division'] == 'Bourgogne-Franche Comte') |
        (location_df['division'] == 'Bourgogne'),
        'division'
    ] = 'Bourgogne-Franche-Comté'

    location_df.loc[
        (location_df['division'] == 'ARA'),
        'division'
    ] = 'Auvergne-Rhône-Alpes'

    location_df.loc[
        (location_df['division'] == 'Centre - Val de Loire'),
        'division'
    ] = 'Centre-Val de Loire'

    location_df.loc[
        (location_df['division'] == 'Grand-Est') |
        (location_df['division'] == 'Grand-est'),
        'division'
    ] = 'Grand Est'

    location_df.loc[
        (location_df['division'] == 'Hauts De France') |
        (location_df['division'] == 'Hauts de France'),
        'division'
    ] = 'Hauts-de-France'

    location_df.loc[
        (location_df['division'] == 'IDF') |
        (location_df['division'] == 'Ile De France') | 
        (location_df['division'] == 'Ile de France') | 
        (location_df['division'] == 'Ile-de-France'),
        'division'
    ] = 'Île-de-France'

    # GEORGIA
    # -------
    # Move to Asia

    location_df.loc[
        (location_df['country'] == 'Georgia') &
        (location_df['region'] == 'Europe'),
        'region'
    ] = 'Asia'

    # GERMANY
    # -------
    # Move Munich division into Bavaria division
    # Duesseldorf --> North Rhine Westphalia

    location_df.loc[
        (location_df['country'] == 'Germany') &
        (location_df['division'] == 'Munich'),
        ['division', 'location']
    ] = ['Bavaria', 'Munich'] 

    location_df.loc[
        (location_df['country'] == 'Germany') &
        (location_df['division'] == 'Duesseldorf'),
        ['division', 'location']
    ] = ['North Rhine Westphalia', 'Duesseldorf'] 

    # ITALY
    # -----
    # Fix typos
    # Castel di Sangro --> Abruzzo
    # Teramo --> Abruzzo

    location_df.loc[
        (location_df['country'] == 'Italy') &
        (location_df['division'] == 'Lombardia'),
        'division'
    ] = 'Lombardy'

    location_df.loc[
        (location_df['country'] == 'Italy') &
        (location_df['division'] == 'Castel di Sangro'),
        ['division', 'location']
    ] = ['Abruzzo', 'Castel di Sangro']

    location_df.loc[
        (location_df['country'] == 'Italy') &
        (location_df['division'] == 'Teramo'),
        ['division', 'location']
    ] = ['Abruzzo', 'Teramo']

    # MOROCCO
    # -------
    # Move to Africa
    # Fix typos

    location_df.loc[
        (location_df['country'] == 'Morocco'),
        'region'
    ] = 'Africa'

    location_df.loc[
        (location_df['country'] == 'Morocco') &
        (location_df['location'] == 'Cadiz_h'),
        'location'
    ] = 'Cadiz'

    # NETHERLANDS
    # -----------
    # Move municipalities into provinces

    # Drenthe
    location_df.loc[
        (location_df['country'] == 'Netherlands') & (
            (location_df['division'] == 'Coevorden') |
            (location_df['division'] == 'Dalen')
        ),
        ['division', 'location']
    ] = ['Drenthe', 'Coevorden']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Drente'),
        'division'
    ] = 'Drenthe'

    # Flevoland
    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Zeewolde'),
        ['division', 'location']
    ] = ['Flevoland', 'Zeewolde']

    # North Brabant
    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Andel'),
        ['division', 'location']
    ] = ['North Brabant', 'Altena']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Berlicum'),
        ['division', 'location']
    ] = ['North Brabant', 'Sint-Michielsgestel']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Eindhoven'),
        ['division', 'location']
    ] = ['North Brabant', 'Eindhoven']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Helmond'),
        ['division', 'location']
    ] = ['North Brabant', 'Helmond']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Loon op zand'),
        ['division', 'location']
    ] = ['North Brabant', 'Loop op Zand']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Milheeze'),
        ['division', 'location']
    ] = ['North Brabant', 'Milheeze']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Oisterwijk'),
        ['division', 'location']
    ] = ['North Brabant', 'Oisterwijk']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Oss'),
        ['division', 'location']
    ] = ['North Brabant', 'Oss']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Tilburg'),
        ['division', 'location']
    ] = ['North Brabant', 'Tilburg']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Noord Brabant'),
        'division'
    ] = 'North Brabant'

    # North Holland
    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Blaricum'),
        ['division', 'location']
    ] = ['North Holland', 'Blaricum']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Diemen'),
        ['division', 'location']
    ] = ['North Holland', 'Diemen']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Haarlem'),
        ['division', 'location']
    ] = ['North Holland', 'Haarlem']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Naarden'),
        ['division', 'location']
    ] = ['North Holland', 'Gooise Meren']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Noord Holland'),
        'division'
    ] = 'North Holland'

    # South Holland
    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Delft'),
        ['division', 'location']
    ] = ['South Holland', 'Delft']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Hardinxveld Giessendam'),
        ['division', 'location']
    ] = ['South Holland', 'Hardinxveld-Giessendam']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Leiden'),
        ['division', 'location']
    ] = ['South Holland', 'Leiden']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Nieuwendijk'),
        ['division', 'location']
    ] = ['South Holland', 'Nieuwendijk']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Nootdorp'),
        ['division', 'location']
    ] = ['South Holland', 'Nootdorp']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Rotterdam'),
        ['division', 'location']
    ] = ['South Holland', 'Rotterdam']

    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Zuid Holland'),
        'division'
    ] = 'South Holland'

    # Utrecht
    location_df.loc[
        (location_df['country'] == 'Netherlands') &
        (location_df['division'] == 'Houten'),
        ['division', 'location']
    ] = ['Utrecht', 'Houten']

    # POLAND
    # ------
    # Fix typos

    location_df.loc[
        (location_df['country'] == 'Poland') &
        (location_df['division'] == 'Dolnoslakie'),
        'division'
    ] = 'Dolnoslaskie'

    location_df.loc[
        (location_df['country'] == 'Poland') & (
            (location_df['division'] == 'Pomorskie') |
            (location_df['division'] == 'Pomorze')
        ),
        'division'
    ] = 'Pomerania'

    location_df.loc[
        (location_df['country'] == 'Poland') &
        (location_df['division'] == 'Malopolskie'),
        'division'
    ] = 'Malopolska'

    # RUSSIA
    # ------
    # Fix typos

    location_df.loc[
        (location_df['division'] == 'Moscow'),
        'division'
    ] = 'Moscow Region'

    location_df.loc[
        (location_df['country'] == 'Russia') & (
            (location_df['division'] == 'Saint-Petersburg') |
            (location_df['division'] == 'St. Petersburg') |
            (location_df['division'] == 'St.Petersburg')
        ),
        'division'
    ] = 'Saint Petersburg'

    # SPAIN
    # -----
    # Fix typos, unify provinces/locations

    location_df.loc[
        (location_df['division'] == 'BasqueCountry') | 
        (location_df['division'] == 'Basque_Country'),
        'division'
    ] = 'Basque Country'

    location_df.loc[
        (location_df['division'] == 'Castilla La Mancha'),
        'division'
    ] = 'Castilla-La Mancha'

    location_df.loc[
        (location_df['division'] == 'Castilla y Leon'),
        'division'
    ] = 'Castilla y León'

    location_df.loc[
        (location_df['division'] == 'Catalunya'),
        'division'
    ] = 'Catalonia'

    location_df.loc[
        (location_df['division'] == 'Catalunya'),
        'division'
    ] = 'Catalonia'

    location_df.loc[
        (location_df['division'] == 'Comunitat_Valenciana'),
        'division'
    ] = 'Comunitat Valenciana'

    location_df.loc[
        (location_df['division'] == 'Comunitat Valenciana') &
        (location_df['location'] == 'Bonrepos_i_Mirambell'),
        'location'
    ] = 'Bonrepos i Mirambell'

    location_df.loc[
        (location_df['division'] == 'Comunitat Valenciana') &
        (location_df['location'] == 'Canet_d\'En_Berenguer'),
        'location'
    ] = 'Canet d\'En Berenguer'

    location_df.loc[
        (location_df['division'] == 'Comunitat Valenciana') &
        (location_df['location'] == 'El_Puig'),
        'location'
    ] = 'El Puig'

    location_df.loc[
        (location_df['division'] == 'Comunitat Valenciana') &
        (location_df['location'] == 'Grau_de_Sagunt'),
        'location'
    ] = 'Grau de Sagunt'

    location_df.loc[
        (location_df['division'] == 'Comunitat Valenciana') &
        (location_df['location'] == 'Palma_de_Gandia'),
        'location'
    ] = 'Palma de Gandia'

    location_df.loc[
        (location_df['division'] == 'Comunitat Valenciana') &
        (location_df['location'] == 'Tavernes_de_la_Valldigna'),
        'location'
    ] = 'Tavernes de la Valldigna'

    location_df.loc[
        (location_df['division'] == 'Comunitat Valenciana') &
        (location_df['location'] == 'Valencia_h'),
        'location'
    ] = 'Valencia'

    location_df.loc[
        (location_df['country'] == 'Spain') & (
            (location_df['division'] == 'LaRioja') |
            (location_df['division'] == 'La_Rioja')
        ),
        'division'
    ] = 'La Rioja'

    # SWEDEN
    # ------
    # Fix typos

    location_df.loc[
        (location_df['division'] == 'Vasterbotten'),
        'division'
    ] = 'Västerbotten'

    # SWITZERLAND
    # -----------
    # Fix typos

    location_df.loc[
        (location_df['country'] == 'Switzerland') & (
            (location_df['division'] == 'Basel') |
            (location_df['division'] == 'Basel Stadt') |
            (location_df['division'] == 'Basel Land')
        ),
        'division'
    ] = 'Basel-Stadt'

    location_df.loc[
        (location_df['country'] == 'Switzerland') & (
            (location_df['division'] == 'Genève') |
            (location_df['division'] == 'Geneve')
        ),
        'division'
    ] = 'Geneva'

    location_df.loc[
        (location_df['country'] == 'Switzerland') & (
            (location_df['division'] == 'Luzern')
        ),
        'division'
    ] = 'Lucerne'


    # UNITED KINGDOM
    # --------------
    # Fix typos

    location_df.loc[
        (location_df['country'] == 'United Kingdom') &
        (location_df['location'] == 'Northamtonshire'),
        'location'
    ] = 'Northamptonshire'

    # CANADA
    # ------
    # Unabbreviate province names

    location_df.loc[
        (location_df['country'] == 'Canada') &
        (location_df['division'] == 'MB'),
        'division'
    ] = 'Manitoba'

    location_df.loc[
        (location_df['country'] == 'Canada') &
        (location_df['division'] == 'NB'),
        'division'
    ] = 'New Brunswick'

    location_df.loc[
        (location_df['country'] == 'Canada') &
        (location_df['division'] == 'NL'),
        'division'
    ] = 'Newfoundland and Labrador'

    location_df.loc[
        (location_df['country'] == 'Canada') &
        (location_df['division'] == 'NS'),
        'division'
    ] = 'Nova Scotia'

    location_df.loc[
        (location_df['country'] == 'Canada') &
        (location_df['division'] == 'SK'),
        'division'
    ] = 'Saskatchewan'

    # MEXICO
    # ------
    # Unabbreviate province names

    location_df.loc[
        (location_df['country'] == 'Mexico') &
        (location_df['division'] == 'CDMX'),
        'division'
    ] = 'Mexico City'

    # North America
    # -------------
    # Who misspelled north??

    location_df.loc[
        location_df['region'] == 'Noth America',
        'region'
    ] = 'North America'

    # USA
    # ---

    # Washington DC
    # -------------
    # Unify with "District of Columbia"

    location_df.loc[
        (location_df['division'] == 'District of Columbia'),
        'division'
    ] = 'Washington DC'

    # California
    # ----------
    # Unify county names

    location_df.loc[
        (location_df['division'] == 'California') & (
            (location_df['location'] == 'Grand Princess') |
            (location_df['location'] == 'Grand Princess cruise ship')
        ),
        'location'
    ] = 'Grand Princess Cruise Ship'

    location_df.loc[
        (location_df['division'] == 'California') & (
            (location_df['location'] == 'San Diego')
        ),
        'location'
    ] = 'San Diego County'

    location_df.loc[
        (location_df['division'] == 'California') & (
            (location_df['location'] == 'San Francisco')
        ),
        'location'
    ] = 'San Francisco County'

    # Louisiana
    # ---------
    # Assume LA is Louisiana

    location_df.loc[
        (location_df['division'] == 'LA'),
        'division'
    ] = 'Louisiana'

    # New York
    # --------
    # Assume NY is NY State
    # Unify county names
    # Move NYC division into NY State division
    # Merge NYC boroughs into one NYC

    location_df.loc[
        (location_df['division'] == 'NY'),
        'division'
    ] = 'New York'

    location_df.loc[
        (location_df['division'] == 'New York City'),
        ['division', 'location']
    ] = ['New York', 'New York City']

    # Remove redundant NY
    location_df.loc[
        (location_df['division'] == 'New York') & (
            (location_df['location'] == 'New York')
        ),
        'location'
    ] = -1

    location_df.loc[
        (location_df['division'] == 'New York') & (
            (location_df['location'] == 'Bronx') |
            (location_df['location'] == 'Brooklyn') |
            (location_df['location'] == 'Manhattan') |
            (location_df['location'] == 'Queens') |
            (location_df['location'] == 'Staten Island')
        ),
        'location'
    ] = 'New York City'

    location_df.loc[
        (location_df['division'] == 'New York') & (
            (location_df['location'] == 'Nassau') |
            (location_df['location'] == 'Nassau county')
        ),
        'location'
    ] = 'Nassau County'

    location_df.loc[
        (location_df['division'] == 'New York') & (
            (location_df['location'] == 'Rockland')
        ),
        'location'
    ] = 'Rockland County'

    location_df.loc[
        (location_df['division'] == 'New York') & (
            (location_df['location'] == 'Suffolk') |
            (location_df['location'] == 'Suffolk county')
        ),
        'location'
    ] = 'Suffolk County'

    location_df.loc[
        (location_df['division'] == 'New York') & (
            (location_df['location'] == 'Westchester')
        ),
        'location'
    ] = 'Westchester County'

    # Washington
    # ----------
    # Move towns into counties

    location_df.loc[
        (location_df['country'] == 'USA') &
        (location_df['division'] == 'Washington') & (
            (location_df['location'] == 'Kirkland') &
            (location_df['location'] == 'Seattle')
        ),
        'location'
    ] = 'King County'

    location_df.loc[
        (location_df['country'] == 'USA') &
        (location_df['division'] == 'Washington') & (
            (location_df['location'] == 'Tacoma')
        ),
        'location'
    ] = 'Pierce County'

    # Wisconsin
    # ---------
    # Unify county names

    location_df.loc[
        (location_df['division'] == 'Wisconsin') & (
            (location_df['location'] == 'Campbellsp')
        ),
        'location'
    ] = 'Campbellsport'

    location_df.loc[
        (location_df['division'] == 'Wisconsin') & (
            (location_df['location'] == 'Jackson')
        ),
        'location'
    ] = 'Jackson County'

    # Australia
    # ---------
    # Fix typos, unabbreviate province names

    location_df.loc[
        (location_df['country'] == 'Australia') &
        (location_df['division'] == 'NSW'),
        'division'
    ] = 'New South Wales'

    location_df.loc[
        (location_df['country'] == 'Australia') &
        (location_df['division'] == 'Northern territory'),
        'division'
    ] = 'Northern Territory'

    # Brazil
    # ------
    # Fix typos
    # Remove redundant Sao Paolo location

    location_df.loc[
        (location_df['country'] == 'Brazil') &
        (location_df['division'] == 'Minas gerais'),
        'division'
    ] = 'Minas Gerais'

    location_df.loc[
        (location_df['country'] == 'Brazil') &
        (location_df['division'] == 'São Paulo'),
        'division'
    ] = 'Sao Paulo'

    # Remove Sao Paulo location
    location_df.loc[
        (location_df['country'] == 'Brazil') &
        (location_df['location'] == 'Sao Paulo'),
        'location'
    ] = -1

    # Central America
    # ---------------
    # Move to North America

    location_df.loc[
        (location_df['country'] == 'Costa Rica'),
        'region'
    ] = 'North America'

    location_df.loc[
        (location_df['country'] == 'Panama'),
        'region'
    ] = 'North America'

    # Done
    return location_df


def build_select_tree(unique_location_df):
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
    select_tree = build_select_tree(unique_location_df)
    # print(select_tree)


if __name__ == '__main__':
    main()

