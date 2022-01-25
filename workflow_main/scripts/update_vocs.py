#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
import re
from pathlib import Path


def clean_text(str):
    str = str.replace(" ", "")
    str = str.split('+')[0]
    str = str.split('/')
    return str


def update_vocs(voc_files):
    variant_dict = {}
    for filepath in voc_files:
        org = Path(filepath).stem
        with open(filepath) as voc_list:
            temp_dict = json.load(voc_list)
            for variant in temp_dict:
                name = clean_text(variant['name'])
                if isinstance(name, list):
                    for n in name:
                        variant_dict = addVariantToDict(variant_dict,
                                                        variant,
                                                        n, org)
                else:
                    variant_dict = addVariantToDict(variant_dict, variant,
                                                    name, org)

    return variant_dict


def addVariantToDict(variant_dict, variant, name, org):
    # If the name does not start in a letter, have any number of
    # '.'s followed by numbers, and end in a number
    # This is a malformed lineage name, do nothing
    if not re.match("^[A-Z]+(.[0-9])+$", name):
        return variant_dict

    if name not in variant_dict:
        # Initialize variant into dict
        variant_dict[name] = {}
        variant_dict[name][org] = variant['level']
    elif org in variant_dict[name]:
        # If variant already in dict and org already in variant
        # Keep the highest level a variant is given per org
        if compareLevels(variant['level'], variant_dict[name][org]):
            variant_dict[name][org] = variant['level']
    else:
        # If variant in dict but org not in variant
        variant_dict[name][org] = variant['level']

    return variant_dict


def compareLevels(newLvl, currLvl):
    if currLvl == "VOC":
        # Never overwrite VOC assignment
        return False

    if newLvl == "VOC" or (newLvl == "VOI" and currLvl == "Other"):
        return True
    else:
        return False


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, nargs='+',
                        help="List of paths to voc jsons")
    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Path to output file")

    args = parser.parse_args()

    variant_dict = update_vocs(args.input)

    with open(args.output, 'w') as fp:
        fp.write(json.dumps(variant_dict, indent=2))


if __name__ == "__main__":
    main()
