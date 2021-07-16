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
    # Start with CDC website and construct an initial dict of variants
    url_list = ['https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html',
                'https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/']
    org_list = ["CDC", "WHO"]
    columns = ['spike_subs',
               'nextstrain',
               'who_label',
               'first_detection',
               'ref_isolate',
               'attributes']
   variant_list = []

   for idx, url in enumerate(url_list):
        vocPage = requests.get(url)

        soup = BeautifulSoup(vocPage.content, 'html.parser')

        if org_list[idx] == "CDC":
            variantTables = soup.find_all('div', class_='table-responsive')
            level_order = ['VOI', 'VOC', 'VOHC']
            # level_ind keeps trach of which table is being processed
            level_ind = 0

            # Iterate through tables, constructing a dict for each variant
            for table in variantTables:
                # Tables are ordered: First is VOI, then VOC, then VOHC
                for row in table.tbody.find_all('tr'):
                    variant = {'level': {"CDC": level_order[level_ind]}}
                    variant['name'] = row.th.text

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
                            variant[columns[columnInd]] = clean_text(text)
                        else:
                            variant[columns[columnInd]] = text

                        columnInd += 1

                    variant_list.append(variant)
                level_ind += 1
        elif org_list[idx] == "WHO":
            variantTables = soup.find_all('tbody')
            level_order = ['VOC', 'VOI', 'AfFM']
            # level_ind keeps trach of which table is being processed
            level_ind = 0
            for table in variantTables:
                if level_order[level_ind] == "AfFM":
                    break
                # Tables are ordered: First is VOC, then VOI, then Alerts for Further Monitoring
                for row in table.find_all('tr'):
                    name = row.find_all('td')[1].split(' ')[0]
                    in_list = False
                    for variant in variant_list:
                        # Check if WHO variant is in the variant_list
                        if variant['name'] == name:
                            # If so, add to the level key of the variant dict
                            variant['level']['WHO'] = level_order[level_ind]
                            in_list = True
                            break
                    if not in_list:
                        # Add to variant_list
                level_ind++


    df = pd.DataFrame(variant_list)
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
