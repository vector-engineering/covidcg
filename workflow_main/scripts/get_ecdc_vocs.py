#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import requests
from bs4 import BeautifulSoup


def get_ecdc_vocs():
    url = 'https://www.ecdc.europa.eu/en/covid-19/variants-concern'
    level_order = ['VOC', 'VOI', 'Other']
    variant_list = []

    vocPage = requests.get(url)
    soup = BeautifulSoup(vocPage.content, 'html.parser')

    # Find all tables (VOC, VOI, and Other)
    variantTables = soup.find_all('tbody')
    level_ind = 0

    for table in variantTables:
        if level_ind >= len(level_order):
            # Break if there is a table we are not interested in
            break
        level = level_order[level_ind]
        for row in table.find_all('tr'):
            name = list(row.find_all('td')[1].stripped_strings)[0]
            variant = {'name': name, 'level': level}
            variant_list.append(variant)
        level_ind += 1

    df = pd.DataFrame(variant_list)
    df.set_index('name')

    return df


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Path to output file")

    args = parser.parse_args()

    df = get_ecdc_vocs()
    df.to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
