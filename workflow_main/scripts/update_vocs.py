#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd
import json
from os import listdir, path


def update_vocs(voc_dir):
    variant_dict = {}
    filenames = listdir(voc_dir)
    for filename in filenames:
        org = filename.split('.')[0]
        filepath = path.join(voc_dir, filename)
        with open(filepath) as voc_list:
            temp_dict = json.load(voc_list)
            for variant in temp_dict:
                name = variant['name']
                if name in variant_dict:
                    return
                else:
                    return

    return


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str,
                        help="Path to dir with seperate voc lists")
    parser.add_argument("-o", "--output", required=True, type=str,
                        help="Path to output file")

    args = parser.parse_args()

    df = update_vocs(args.input)
    df.to_json(args.output, orient="records")


if __name__ == "__main__":
    main()
