#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import requests
import re
from bs4 import BeautifulSoup


def clean_text(s):
    s = s.replace('\u2013', '-')
    s = s.replace('\u00a0', '')
    s = s.replace('\\', '')
    s = s.replace('\n', '')
    return s


def update_vocs():
    url = 'https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html'
    level_order = ['VOI', 'VOC', 'VOHC']
    columns = ['spike_subs',
               'nextstrain',
               'who_label',
               'first_detection',
               'ref_isolate',
               'attributes']

    vocPage = requests.get(url)

    soup = BeautifulSoup(vocPage.content, 'html.parser')

    variantTables = soup.find_all('div', class_='table-responsive')

    # level_ind keeps trach of which table is being processed
    level_ind = 0
    dict_list = []

    # Iterate through tables, constructing a dict for each variant
    for table in variantTables:
        # Tables are ordered: First is VOI, then VOC, then VOHC
        for row in table.tbody.find_all('tr'):
            dict = {'level': level_order[level_ind]}
            dict['name'] = row.th.text

            # columnInd tracks which column of the table we are in
            columnInd = 0
            for td in row.find_all('td'):
                # Check for columns that require more processing
                text = ''
                if columnInd == 0:
                    # Turn the list of spike substitutions into an array
                    text = td.text.replace('Spike: ', '')
                    text = text.split(', ')
                elif columnInd == 4:
                    # Remove "external icon" from string
                    text = td.text.replace('external icon', '')
                elif columnInd == 5:
                    # Turn list of attributes into an array
                    arr = []
                    for li in td.find_all('li'):
                        text = re.sub(r"[\d-]", '', li.text)
                        text = text.strip(' ,')
                        text = clean_text(text)
                        arr.append(text)
                    text = arr
                else:
                    text = td.text

                if isinstance(text, str):
                    dict[columns[columnInd]] = clean_text(text)
                else:
                    dict[columns[columnInd]] = text

                columnInd += 1

            dict_list.append(dict)
        level_ind += 1

    df = pd.DataFrame(dict_list)
    df.set_index('name')

    return df


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o", "--output", required=True, type=str, help="Path to output file"
    )

    args = parser.parse_args()

    df = update_vocs()
    df.to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
