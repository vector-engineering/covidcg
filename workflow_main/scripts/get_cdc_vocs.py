#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import requests
from bs4 import BeautifulSoup


def get_cdc_vocs():
    url = ('https://www.cdc.gov/coronavirus/2019-ncov/'
           'variants/variant-info.html')
    level_order = ['VOI', 'VOC']
    variant_list = []

    vocPage = requests.get(url)
    soup = BeautifulSoup(vocPage.content, 'html.parser')

    # Find all tables (VOI, VOC)
    variantTables = soup.find_all('div', class_='step-table m-3')
    level_ind = 0

    for table in variantTables:
        if level_ind >= len(level_order):
            # Break if there is a table we are not interested in
            break
        level = level_order[level_ind]
        rowgroup = table.div
        rows = rowgroup.find_all('div', class_='col-md-12 p-md-3 pt-0 pb-3')
        for row in rows:
            variants = list(row.stripped_strings)
            variant = {'name': variants[0],
                       'level': level}
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

    df = get_cdc_vocs()
    df.to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
