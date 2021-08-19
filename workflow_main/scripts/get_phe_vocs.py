#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
import requests
from bs4 import BeautifulSoup


def get_phe_vocs():
    url = ("https://www.gov.uk/government/publications/"
           "covid-19-variants-genomically-confirmed-case-numbers")

    homePage = requests.get(url)
    homeSoup = BeautifulSoup(homePage.content, 'html.parser')

    # Find link to newest report
    variant_div = homeSoup.find('div', class_='attachment-details')
    href = variant_div.find('a', href=True)['href']
    variant_url = "https://www.gov.uk" + href

    # Scrape newest variant report
    variantPage = requests.get(variant_url)
    soup = BeautifulSoup(variantPage.content, 'html.parser')
    variant_list = []

    variant_table = soup.find('tbody')
    for row in variant_table.find_all('tr'):
        level = []
        # Variant levels for the PHE are VOC or VUI
        # col 1 of the variant row has the variant's level in the form:
        # Level-DDMM-YY where the date is the day the level was assigned
        level = list(row.find_all('td')[1].stripped_strings)
        # Col 1 sometimes has N/A. In these cases, the level is in col 0
        if level[0] == 'N/A':
            level = list(row.find_all('td')[0].stripped_strings)
        level = level[0]

        # Rename VUI to VOI for ease of processing
        if level == 'VUI':
            level = 'VOI'

        # Lineage is in col 2. Some have multiple lineages separated by spaces
        lineages = list(row.find_all('td')[2].stripped_strings)
        # Remove any accidental spaces
        lineages = lineages[0].replace('. ', '.')
        lineage_list = lineages.split(' ')
        for index, lineage in enumerate(lineage_list):
            if lineage.find('.') != -1:
                variant = {'name': lineage, 'level': level}
                variant_list.append(variant)

    return variant_list


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Path to output file")

    args = parser.parse_args()

    variant_list = get_phe_vocs()
    with open(args.output, 'w') as fp:
        fp.write(json.dumps(variant_list, indent=2))


if __name__ == "__main__":
    main()
