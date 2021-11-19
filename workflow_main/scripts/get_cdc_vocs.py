#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
import requests
from bs4 import BeautifulSoup


def get_cdc_vocs():
    url = ('https://www.cdc.gov/coronavirus/2019-ncov/'
           'variants/variant-info.html')
    variant_list = []

    vocPage = requests.get(url)
    soup = BeautifulSoup(vocPage.content, 'html.parser')
    lineageRowArr = None

    # Get VBM table
    table = soup.find('table', class_='table table-bordered nein-scroll')
    level = 'Other'

    vbmRows = table.find_all('tr')

    for i, row in enumerate(vbmRows):
        if i > 0:
            vbmVariantStr = row.find_all('td')[1].text
            vbmVariantArr = vbmVariantStr.split('and')
            mainVBMVariant = vbmVariantArr[0]
            vbmNameArr = clean_text(mainVBMVariant)
            for vbmName in vbmNameArr:
                variant = {'name': vbmName, 'level': level}
                variant_list.append(variant)

            if len(vbmVariantArr) > 1:
                # Logic for CDC claiming whole lineages as VBMs
                # i.e. "B.1.617.2 and descendant lineages"
                if "lineages" in vbmVariantArr[1].strip():
                    if lineageRowArr is None:
                        # Get all covid lineages
                        covLinURL = "https://cov-lineages.org/lineage_list.html"
                        covLinPage = requests.get(covLinURL)
                        covLinSoup = BeautifulSoup(covLinPage.content, 'html.parser')
                        lineageTable = covLinSoup.find('table', class_='table')
                        lineageRowArr = lineageTable.find_all('tr')

                    target = vbmVariantArr[1].strip().split(' ')[0]
                    if target == "descendent":
                        target = vbmNameArr[0]

                    variant_list = get_all_lineages(target, variant_list,
                                                    level, lineageRowArr)

    # Get VOC table
    table = soup.find('div', class_='mb-0 pt-3')
    level = 'VOC'

    vocRows = table.find_all('p')
    vocArr = list(vocRows[1].stripped_strings)
    name = vocArr[1].split(' ')[0]
    variant = {"name": name, "level": level}
    variant_list.append(variant)

    # The CDC classifies all AY lineages as VOCs but we cannot display 120
    # buttons in the columns we have so we're just writing a note instead

    # variant_list = get_all_lineages(vocArr[1].split(' ')[-1], variant_list,
    #                                level, lineageRowArr)

    return variant_list


def clean_text(str):
    # Remove trailing whitespace and check for commas
    str = str.replace('\xa0', "")
    arr = str.strip().split(',')

    if len(arr) == 1:
        arr = arr[0].split('\n')

    return arr


def get_all_lineages(target, variant_list, level, lineageRowArr):
    for i, lineageRow in enumerate(lineageRowArr):
        if i == 0:
            continue
        currLineage = lineageRow.find_all('td')[0].text
        if (len(currLineage) > len(target) and
           currLineage[0:len(target) + 1] == target + "."):
            variant = {"name": currLineage, "level": level}
            variant_list.append(variant)

    return variant_list


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Path to output file")
    args = parser.parse_args()

    variant_list = get_cdc_vocs()

    with open(args.output, 'w') as fp:
        fp.write(json.dumps(variant_list, indent=2))


if __name__ == "__main__":
    main()
