#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
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
            # Find all of the pango lineages in each row of the table
            names = list(row.find_all('td')[0].stripped_strings)

            for name in names:
                # Strip any non alphanumeric char that is not a .
                name = ''.join(e for e in name if e.isalnum() or e == '.')
                # If the name is valid, add the variant
                if name:
                    variant = {'name': name, 'level': level}
                    variant_list.append(variant)
        level_ind += 1

    return variant_list


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Path to output file")

    args = parser.parse_args()

    variant_list = get_who_vocs()

    with open(args.output, 'w') as fp:
        fp.write(json.dumps(variant_list, indent=2))


if __name__ == "__main__":
    main()
