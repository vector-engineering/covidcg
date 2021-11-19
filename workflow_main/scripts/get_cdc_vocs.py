#!/usr/bin/env python3
# coding: utf-8

import re
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

    # Get VBM table
    vbmTable = soup.find('table', class_='table table-bordered nein-scroll')
    level = 'Other'

    vbmRows = vbmTable.find_all('tr')

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
                if "lineages" in vbmVariantArr[1].strip():
                    # Get all covid lineages
                    covLinURL = "https://cov-lineages.org/lineage_list.html"
                    covLinPage = requests.get(covLinURL)
                    covLinSoup = BeautifulSoup(covLinPage.content, 'html.parser')
                    lineageTable = covLinSoup.find('table', class_='table')
                    lineageRowArr = lineageTable.find_all('tr')

                    # Get target lineage
                    # i.e. the Q if vbmVariantArr[1] = 'Q lineages'
                    target = vbmVariantArr[1].strip().split(' ')[0]

                    if target == "descendent":
                        target = vbmNameArr[0]

                    for i, lineageRow in enumerate(lineageRowArr):
                        if i == 0:
                            continue
                        currLineage = lineageRow.find_all('td')[0].text
                        if (len(currLineage) > len(target) and currLineage[0:len(target) + 1] == target + "."):
                            variant = {"name": currLineage, "level": level}
                            variant_list.append(variant)

    # Get VOC table

    return variant_list


def clean_text(str):
    # Remove trailing whitespace and check for commas
    str = str.replace('\xa0', "")
    arr = str.strip().split(',')

    if len(arr) == 1:
        arr = arr[0].split('\n')

    return arr


def main():
    print(get_cdc_vocs())
#    parser = argparse.ArgumentParser()
#
#    parser.add_argument("-o", "--output", required=True, type=str,
#                        help="Path to output file")
#    args = parser.parse_args()
#
#    variant_list = get_cdc_vocs()
#
#    with open(args.output, 'w') as fp:
#        fp.write(json.dumps(variant_list, indent=2))


if __name__ == "__main__":
    main()
