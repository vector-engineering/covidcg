# coding: utf-8

"""Create standalone global sequencing vega spec, with inlined data

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import json


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--in-spec", type=str, required=True, help="Path to input vega spec",
    )
    parser.add_argument(
        "--country-score", type=str, required=True, help="Path to country score data",
    )
    parser.add_argument(
        "--out-spec", type=str, required=True, help="Path to output vega spec",
    )

    args = parser.parse_args()

    standalone_map_spec(args.in_spec, args.country_score, args.out_spec)


def standalone_map_spec(in_spec, data, out_spec):
    # Load vega spec JSON
    with open(in_spec, "r") as fp:
        spec = json.loads(fp.read())

    # Load data JSON
    with open(data, "r") as fp:
        data = json.loads(fp.read())

    # Inject data array into spec object
    spec["data"][0]["values"] = data

    # Save spec with injected data
    with open(out_spec, "w") as fp:
        fp.write(json.dumps(spec))


if __name__ == "__main__":
    main()
