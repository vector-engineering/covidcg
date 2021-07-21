#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import requests
from bs4 import BeautifulSoup


def get_who_vocs():
    url = 'https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/'
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
            label = ''
            names = []
            if level == 'VOC' or level == 'VOI':
                # Find all of the pango lineages in each row of the table
                names = list(row.find_all('td')[1].stripped_strings)
                # Determine who label
                label = list(row.find_all('td')[0].stripped_strings)[0]
            else:
                names = list(row.find_all('td')[0].stripped_strings)

            for name in names:
                variant = {}
                if label:
                    variant = {
                               'name': name,
                               'level': level,
                               'who_label': label
                               }
                else:
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

    df = get_who_vocs()
    df.to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
