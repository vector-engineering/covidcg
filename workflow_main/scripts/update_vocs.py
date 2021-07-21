#!/usr/bin/env python3
# coding: utf-8

import argparse
import json
from pathlib import Path


def update_vocs(voc_files):
    variant_dict = {}
    for filepath in voc_files:
        org = Path(filepath).stem
        with open(filepath) as voc_list:
            temp_dict = json.load(voc_list)
            for variant in temp_dict:
                name = variant['name']
                if name not in variant_dict:
                    variant_dict[name] = {}
                variant_dict[name][org] = variant['level']
                if org == 'who' and variant['who_label']:
                    variant_dict[name]['who_label'] = variant['who_label']

    return variant_dict


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
