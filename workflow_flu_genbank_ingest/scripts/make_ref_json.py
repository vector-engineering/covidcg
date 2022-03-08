#!/usr/bin/env python3
# coding: utf-8

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("data", type=str, help="Path to whole data file")
    parser.add_argument("out", type=str, help="Path to output file")

    args = parser.parse_args()

    df = pd.read_csv(args.data)

    dff = df.loc[df['database'] == 'RefSeq']
    dff.to_json(args.out)


if __name__ == '__main__':
    main()