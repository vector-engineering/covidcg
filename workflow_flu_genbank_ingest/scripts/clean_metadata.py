#!/usr/bin/env python3
# coding: utf-8

"""Clean metadata from GenBank

Required columns:
    "Accession ID": index

    # date information
    "collection_date"
    "submission_date"

    # location information
    "region"
    "country"
    "division"
    "location"

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import argparse
import datetime
import json
import pandas as pd


def protein_to_segment(x):
    rename_map = {
        "pB1": "PB1",
        "Pb1": "PB1",
        "\(PB1\)": "PB1",
        "Pb2": "PB2",
        "\(PB2\)": "PB2",
        "Polymerase": "polymerase",
        "poylmerase": "polymerase",
        "plymerase": "polymerase",
        "polymerae": "polymerase",
        "polemerase": "polymerase",
        "polymelase": "polymerase",
        "polymerse": "polymerase",
        "\(P3\)": "P3",
        "\(PA\)": "PA",
        "Hemagglutinin": "hemagglutinin",
        "Haemagglutinin": "hemagglutinin",
        "hemaggluinin": "hemagglutinin",
        "hemaglutinin": "hemagglutinin",
        "haemagglutinin": "hemagglutinin",
        "hemagglurinin": "hemagglutinin",
        "heamagglutinin": "hemagglutinin",
        "hemagglutnin": "hemagglutinin",
        "hemagglutinine": "hemagglutinin",
        "hemaggulutinin": "hemagglutinin",
        "haemaglutinin": "hemagglutinin",
        "hemagglutunin": "hemagglutinin",
        "hemmaglutinin": "hemagglutinin",
        "hemegglutinin": "hemagglutinin",
        "hemaggultinin": "hemagglutinin",
        "prehemagglutinin": "hemagglutinin",
        "haemagluttinin": "hemagglutinin",
        "haemagglutin": "hemagglutinin",
        r"hemagglutin$": "hemagglutinin",
        "hemagluttinin": "hemagglutinin",
        "hemagglutinn": "hemagglutinin",
        "hemagglutintin": "hemagglutinin",
        "hemmagulutinin": "hemagglutinin",
        "hamagglutinin": "hemagglutinin",
        "hemaggluitnin": "hemagglutinin",
        "haemmagglutinin": "hemagglutinin",
        "hemmagglutinin": "hemagglutinin",
        "HAY": "HA",
        "neuraminidase-like": "neuraminidase",
        "neuramindase": "neuraminidase",
        "Neuraminidase": "neuraminidase",
        "neuramidase": "neuraminidase",
        "nueraminidase": "neuraminidase",
        "neuraminadase": "neuraminidase",
        "neuramidinase": "neuraminidase",
        "neuroaminidase": "neuraminidase",
        r"neuraminidas$": "neuraminidase",
        "nueuraminidase": "neuraminidase",
        "neruaminidase": "neuraminidase",
        "nucelar": "nuclear",
        "Nucleoprotein": "nucleoprotein",
        "nuleoprotein": "nucleoprotein",
        "necleoprtotein": "nucleoprotein",
        "Nucleocapsid": "nucleocapsid",
        "nucelocapsid": "nucleocapsid",
        "Nonstructrual": "nonstructural",
        "non-structrual": "nonstructural",
        "nonstractual": "nonstructural",
        "Non structural": "nonstructural",
        "Non-structural": "nonstructural",
        "Nonstructural": "nonstructural",
        "nonstrucral": "nonstructural",
        "Structural": "structural",
        "structual": "structural",
        "Protein": "protein",
        "prtoein": "protein",
        "proptein": "protein",
        "protien": "protein",
        "glytoprotein": "glycoprotein",
        "NS-2": "NS2",
        "NS-1": "NS1",
        "ns2": "NS2",
        "ns1": "NS1",
        "\(NS1\)": "NS1",
        "\(NS2\)": "NS2",
        "martix": "matrix",
        "matirx": "matrix",
        "matriX": "matrix",
        "matix": "matrix",
        "mstrix": "matrix",
        "Matrix": "matrix",
        "\(M1\)": "M1",
        "\(M2\)": "M2",
        "truncated ": "",
        "putative ": "",
        "pre-": "",
        "\sprecursor": "",
        "\sprepropeptide": "",
        "unspliced product of ": "",
        "\s\(partial\)": "",
        "incomplete ": "",
        "baseic": "basic",
        " \(high yield phenotype\)": "",
        " \(HA1 domain\)": "",
        " precusor HA1 chain": "",
        " \(partial[\s0-9A-Za-z,-]+\)$": "",
    }

    protein_to_segment_map = {
        r"polymerase(\scomplex)?(\ssubunit)? PB2": "1",
        r"(RNA\s)?(polymerase|protein)(\spolymerase|\sprotein)? basic(\s(protein|subunit|complex))?(\s(protein|subunit|complex))?\s?(PB)?2": "1",
        r"(cap-binding protein\s)?PB2(\spolymerase)?(\ssubunit|\sprotein)?": "1",
        r"polymerase basic 2 protein": "1",
        r"polymerase protein(\s2)? PB2": "1",
        r"polymerase B2": "1",
        r"basic(\spolymerase)?(\ssubunit|\sprotein)? 2": "1",
        r"PB protein2": "1",
        r"PB2 polymerase(\ssubunit|\sprotein|\s2)": "1",
        r"polymerase 3 \(pb2\)": "1",
        r"polymerase(\scomplex)?(\ssubunit)? PB1": "2",
        r"(RNA\s)?(polymerase|protein)(\spolymerase|\sprotein)? basic(\s(protein|subunit|complex))?(\s(protein|subunit|complex))?\s?(PB)?1([-\s]F2|\sPB1)?": "2",
        r"(polymerase\s)?PB1(-F2)?(\sprotein)?": "2",
        r"PB1(-F2)? polymerase(\ssubunit|\sprotein|\s1)?": "2",
        r"(2-)?F2": "2",
        r"polymerase basic 1 protein": "2",
        r"(polymerase|basic) protein 1(\sPB1)?": "2",
        r"RNA-directed RNA polymerase subunit P1": "2",
        r"polymerase protein PB1": "2",
        r"polymerase B1": "2",
        r"PB1(, PB1-F2|-N40| F2)?(\sprotein)?": "2",
        r"basic(\spolymerase)?(\ssubunit|\sprotein)? 1": "2",
        r"N40 protein": "2",
        r"PB protein1": "2",
        r"polymerase(\ssubunit|\sprotein|\scomplex subunit|\s3)? PA(-X)?(\sprotein)?": "3",
        r"PA(-X)? protein": "3",
        r"PA(, PA-X|-X)?": "3",
        r"PA polymerase(\ssubunit|\sprotein)?(\s1)?": "3",
        r"(RNA\s)?polymerase(\sprotein)? A": "3",
        r"polymerase (acid|A|acidic)(\sprotein|\ssubunit)?(\s2)?(\s[A-Z0-9]{3})?": "3",
        r"acidic (polymerase|protein)(\ssubunit)?(\s[12])?": "3",
        r"RNA-directed RNA polymerase subunit P2": "3",
        r"polymerase(\ssubunit)? P3": "3",
        r"polymerase protein 3 P3": "3",
        r"P3": "3",
        r"(H[35]\s)?hemagglutinin(\sprotein|\sprecursor|[-\s]esterase|\sglycoprotein|\sgene|\sHA[12]?|\sHA-?1? (chain|subunit|region)|\sH[0-9](\ssubtype)?|\s5|\sA)?": "4",
        r"4 HA1 domain": "4",
        r"4 glycoprotein": "4",
        r"hemagglutinin(\ssubunit)? HA1(\sdomain|\sregion)?": "4",
        r"hemagglutinin 1(\ssubunit|\schain)?": "4",
        r"hemagglutinin (subunit|subtype) H?[15]": "4",
        r"pre-hemagglutinin": "4",
        r"HA(\sprotein|1|0)?(\shemagglutinin)?": "4",
        r"hemagglutinin protein subunit 1": "4",
        r"H3HA1 surface glycoprotein": "4",
        r"HA subunit of hemagglutinin": "4",
        r"hemagglutinin-esterase-fusion(\sglycoprotein|\sprotein)?": "4",
        r"hemagglutinins HA1 and HA2": "4",
        r"receptor binding and fusion protein": "4",
        r"nucleocapsid(\sprotein)?": "5",
        r"nucleoprotein(\sNP)?": "5",
        r"NP(\sprotein)?": "5",
        r"nucleocapsid protein \(NP\)": "5",
        r"nucleoprotein [A-Z0-9\/]{3,7}": "5",
        r"(N2\s)?neuraminidase(\ssubtype)?(\sglycoprotein|\ssurface protein|\sprotein|\sNA|\sN?[12]|\scell surface glycoprotein|\sNAN1)?(\sN1|\s2)?": "6",
        r"NB (glyco)?protein": "6",
        r"NB": "6",
        r"CM2(\sprotein)?": "6",
        r"NA glycoprotein": "6",
        r"NA(\sprotein)?": "6",
        r"P42(\sprotein)?": "6",
        r"(membrane\s)?matrix\s?protein(\sM?[12])?": "7",
        r"BM2(\sprotein)?": "7",
        r"membrane ion channel(\;?\sM?2|\sprotein)?(\sM[12])?": "7",
        r"(matrix|membrane)(\s?M?[12])?(\sprotein)?(\s?[12M])?": "7",
        r"(7\s)?M[124]*(\s(matrix|membrane))?(\sprotein|\sion channel|\sgene)?": "7",
        r"membrane(\smatrix)? protein M[12]": "7",
        r"(trans)?membrane protein": "7",
        r"matrix\s?protein(\s[12])?\s?M?[12]": "7",
        r"MP": "7",
        r"nuclear export protein(\s2|\s1|\sNEP)?": "8",
        r"(NS1\s|NS2\s)?non[-\s]?structural(\s[12])?(\sprotein)?(\s?(NS)?[123])?(\sNS[12]?)?": "8",
        r"(8\s)?NS[12]?(\sprotein)?": "8",
        r"NEP(\/NS2)?": "8",
        r"NS2/NEP": "8",
        r"nuclear export protein NS2": "8",
        r"NSP?": "8",
        r"nuclear export protein/nonstructural protein 2": "8",
    }

    blacklist = [
        r"(RNA\s)?polymerase(\sprotein)?(\s[123])?",
        r"polymerase subunit(\s[12])?",
        r"polymerase [12] protein",
        r"polymerase basic protein",
        r"polyprotein",
        r"nuclear protein",
        r"glycoprotein",
        r"RNA polymerase subunit",
        r"MA",  # matrix?
        r"hypothetical NB-NA hybrid protein",
        r"ORF",
        r"DI-[123] protein",
        r"L protein",
        r"hemagglutinin\/neuraminidase",
        r"integral membrane glycoprotein",
        r"RBD-NEP fusion protein",
    ]

    for k, v in rename_map.items():
        x = x.str.replace(k, v, regex=False)

    for r in blacklist:
        match = x.str.match("^" + r + "$")
        match[match.isna()] = False
        x[match] = None

    for k, v in protein_to_segment_map.items():
        x = x.str.replace("^" + k + "$", v, regex=True)

    # Convert to integers
    is_num = x.str.match("^[0-9]$")
    is_num[is_num.isna()] = False
    x[is_num] = x[is_num].astype(int)

    return x


def collapse_segment_list(segments):
    """Choose consensus segment from list of derived segments"""
    # Remove non-integers from this list
    segments = [x for x in segments if type(x) is int]

    # If empty, return None
    if len(segments) == 0:
        return None
    # If just one, choose that
    elif len(segments) == 1:
        return segments[0]
    else:
        # If segments differ, then return None
        # otherwise, return the consensus
        if sum(segments, 0) / len(segments) != segments[0]:
            return None
        else:
            return segments[0]


def parse_genbank_location(s):
    """Convert a Genbank location string into a tuple of (division, location)
    Fill missing data with "-1" for groupbys later
    """

    # No additional location data
    if not s or type(s) is not str:
        return (-1, -1)

    # The country is always the first part
    country_chunks = s.split(":")
    # country = country_chunks[0].strip()

    # Only the country is defined
    if len(country_chunks) == 1:
        return (-1, -1)

    # If the second part exists, the division and location might
    # be comma-delimited, i.e., "division, location"

    division_chunks = country_chunks[1].split(",")
    division = division_chunks[0].strip()

    # Only the division is defined
    if len(division_chunks) == 1:
        return (division, -1)

    # All are defined
    return (division, division_chunks[1].strip())


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--metadata-in", type=str, required=True, help="Metadata in")
    parser.add_argument(
        "--quality",
        type=str,
        required=True,
        help="Path to sequence quality input CSV file",
    )
    parser.add_argument("--metadata-out", type=str, required=True, help="Metadata out")
    parser.add_argument(
        "--metadata-virus-out", type=str, required=True, help="Aggregate metadata out"
    )

    args = parser.parse_args()

    df = pd.read_csv(
        args.metadata_in,
        usecols=[
            "genbank_accession",
            "set_id",
            "database",
            "genus",
            "serotype",
            "strain",
            "length",
            "segments",
            "genome_coverage",
            "region",
            "country",
            "location",
            "collected",
            "submitted",
            "updated",
            "host",
            "isolation_source",
            "biosample_accession",
            "authors",
            "publications",
        ],
    )

    # Fields:
    #   genbank_accession: string,
    #   set_id: string,
    #   database: string,
    #   genus: string (Alphainfluenzavirus, Betainfluenzavirus)
    #   serotype: string
    #   strain: string,
    #   length: int
    #   is_segmented: bool (dropped)
    #   complete: bool (dropped)
    #   segments: string
    #   genome_coverage: JSON
    #   region: string,
    #   country: string,
    #   location: string,
    #   collected: string (ISO?),
    #   submitted: ISO datetime,
    #   updated: ISO datetime,
    #   host: string,
    #   isolation_source: string,
    #   biosample_accession: string,
    #   title: string (dropped),
    #   authors: string,
    #   publications: string
    #   sequence: string (dropped)

    """
    {"genbank_accession":"NC_002016","database":"RefSeq","genus":"Alphainfluenzavirus","serotype":"H1N1","strain":"A\/Puerto Rico\/8\/1934","length":1027,"is_segmented":true,"complete":"complete","segments":"7","genome_coverage":"[{\"id\": \"NC_002016.1\", \"segment\": \"7\", \"proteins\": [{\"id\": \"NP_040979.2\", \"name\": \"matrix protein 2\", \"location\": \"join(26..51,740..1007)\"}, {\"id\": \"NP_040978.1\", \"name\": \"matrix protein 1\", \"location\": \"26..784\"}]}]","region":"North America","country":"USA","location":"Puerto Rico","collected":"1934","submitted":"1982-06-09T00:00:00Z","updated":"2018-08-13T00:00:00Z","host":null,"isolation_source":null,"biosample_accession":null,"title":"Influenza A virus (A\/Puerto Rico\/8\/1934(H1N1)) segment 7, complete sequence","authors":"Winter,G., Fields,S.","publications":"6927841","sequence":"AGCGAAAGCAGGTAGATATTGAAAGATGAGTCTTCTAACCGAGGTCGAAACGTACGTTCTCTCTATCATCCCGTCAGGCCCCCTCAAAGCCGAGATCGCACAGAGACTTGAAGATGTCTTTGCAGGGAAGAACACCGATCTTGAGGTTCTCATGGAATGGCTAAAGACAAGACCAATCCTGTCACCTCTGACTAAGGGGATTTTAGGATTTGTGTTCACGCTCACCGTGCCCAGTGAGCGAGGACTGCAGCGTAGACGCTTTGTCCAAAATGCCCTTAATGGGAACGGGGATCCAAATAACATGGACAAAGCAGTTAAACTGTATAGGAAGCTCAAGAGGGAGATAACATTCCATGGGGCCAAAGAAATCTCACTCAGTTATTCTGCTGGTGCACTTGCCAGTTGTATGGGCCTCATATACAACAGGATGGGGGCTGTGACCACTGAAGTGGCATTTGGCCTGGTATGTGCAACCTGTGAACAGATTGCTGACTCCCAGCATCGGTCTCATAGGCAAATGGTGACAACAACCAACCCACTAATCAGACATGAGAACAGAATGGTTTTAGCCAGCACTACAGCTAAGGCTATGGAGCAAATGGCTGGATCGAGTGAGCAAGCAGCAGAGGCCATGGAGGTTGCTAGTCAGGCTAGGCAAATGGTGCAAGCGATGAGAACCATTGGGACTCATCCTAGCTCCAGTGCTGGTCTGAAAAATGATCTTCTTGAAAATTTGCAGGCCTATCAGAAACGAATGGGGGTGCAGATGCAACGGTTCAAGTGATCCTCTCGCTATTGCCGCAAATATCATTGGGATCTTGCACTTGATATTGTGGATTCTTGATCGTCTTTTTTTCAAATGCATTTACCGTCGCTTTAAATACGGACTGAAAGGAGGGCCTTCTACGGAAGGAGTGCCAAAGTCTATGAGGGAAGAATATCGAAAGGAACAGCAGAGTGCTGTGGATGCTGACGATGGTCATTTTGTCAGCATAGAGCTGGAGTAAAAAACTACCTTGTTTCTACT"}
    """

    # Rename columns, set index
    df.rename(
        columns={
            "genbank_accession": "Accession ID",
            "set_id": "isolate_id",
            "strain": "virus_name",
            "submitted": "submission_date",
            "collected": "collection_date",
            "updated": "updated_date",
        },
        inplace=True,
    )
    df.set_index("Accession ID", inplace=True)

    # Backfill set_id with strain
    df["isolate_id"].fillna(df["virus_name"], inplace=True)
    # Backfill more with Accession ID
    df["isolate_id"].fillna(pd.Series(df.index.values, index=df.index), inplace=True)

    # Remove sequences without isolate_id, virus_name, region, collection date, or submission date
    # TODO: fill in missing isolate_ids with accession IDs?
    # TODO: investigate whether missing isolate IDs are biased towards
    #       certain segments. i.e., if they're all 4 (HA), then they
    #       can probably be treated as standalone isolates
    # print(f"Missing isolate_id: {df['isolate_id'].isna().sum()}")
    # print(f"Missing virus_name: {df['virus_name'].isna().sum()}")
    print(f"Missing region: {df['region'].isna().sum()}")
    print(f"Missing submission_date: {df['submission_date'].isna().sum()}")
    print(f"Missing collection_date: {df['collection_date'].isna().sum()}")

    remove_rows = (
        (df["isolate_id"].isna())  # Shouldn't happen
        # | (df["virus_name"].isna())
        | (df["region"].isna())
        | (df["submission_date"].isna())
        | (df["collection_date"].isna())
    )

    print(f"Removing {remove_rows.sum()} rows")

    df.drop(df.index[remove_rows], inplace=True)

    # Redundant, but only keep alpha+beta
    df.drop(
        df.index[~df["genus"].isin(["Alphainfluenzavirus", "Betainfluenzavirus"])],
        inplace=True,
    )

    # Rename serotypes
    # TODO: preserve original H and N types here?
    """
    All serotypes from metadata:
    [
    'H1', 'H1N1', 'H1N2', 'HIN2', 'H1N2v', 'H1N3', 'H1N4', 'H1N5', 'H1N6', 'H1N7', 'H1N8', 'H1N9', 'H01N2',
    'H2', 'H2N1', 'H2N2', 'H2N3', 'H2N4', 'H2N5', 'H2N6', 'H2N7', 'H2N8', 'H2N9',
    'H3', 'H3N-', 'H3N1', 'H3N2', 'H3N2v', 'H3N3', 'H3N4', 'H3N5', 'H3N6', 'H3N7', 'H3N8', 'H3N9', 'H3Nx',
    'H4', 'H4N1', 'H4N2', 'H4N3', 'H4N4', 'H4N4,8', 'H4N5', 'H4N6', 'H4N7', 'H4N8', 'H4N9', 'H4Nx',
    'H5', 'H5N1', 'h5n1', 'H5N2', 'H5N2?', 'H5N3', 'H5N4', 'H5N5', 'H5N6', 'H5N7', 'H5N8', 'H5N9', 'H5N?', 'H5Nx',
    'H6', 'H6N', 'H6N1', 'H6N1,4', 'H6N2', 'H6N3', 'H6N4', 'H6N5', 'H6N6', 'H6N7', 'H6N8', 'H6N9', 'H6N?', 'H6Nx',
    'H7', 'H7N1', 'H7N2', 'H7N3', 'H7N4', 'H7N5', 'H7N6', 'H7N7', 'H7N8', 'H7N9', 'H7Nx',
    'H8', 'H8N1', 'H8N2', 'H8N3', 'H8N4', 'H8N5', 'H8N6', 'H8N7', 'H8N8', 'H8Nx',
    'H9', 'H91', 'H9N1', 'H9N1,4', 'H9N1,9', 'H9N2', 'H9N3', 'H9N4', 'H9N5', 'H9N6', 'H9N7', 'H9N8', 'H9N9',
    'H10', 'H10N1', 'H10N2', 'H10N3', 'H10N4', 'H10N5', 'H10N6', 'H10N7', 'H10N8', 'H10N9',
    'H11', 'H11N1', 'H11N2', 'H11N3', 'H11N4', 'H11N5', 'H11N6', 'H11N7', 'H11N8', 'H11N9', 'H11N9/N2', 'H11N?',
    'H12', 'H12N1', 'H12N2', 'H12N3', 'H12N4', 'H12N5', 'H12N6', 'H12N7', 'H12N8', 'H12N9',
    'H13', 'H13N1', 'H13N2', 'H13N3', 'H13N4', 'H13N6', 'H13N8', 'H13N9', 'H13N?',
    'H14N2', 'H14N3', 'H14N4', 'H14N5', 'H14N6', 'H14N7', 'H14N8',
    'H15N2', 'H15N4', 'H15N5', 'H15N6', 'H15N7', 'H15N8', 'H15N9',
    'H16', 'H16N3', 'H16N8', 'H16N9',
    'H17N10',
    'H18N11',
    'H?N1', 'H?N2', 'H?N4', 'H?N6',
    'HxN1', 'HxN3', 'HxNx',
    'Mixed', 'mixed',
    'N1', 'n1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9',
    'B',
    'Victoria',
    'Yamagata', 'Yamagata-like',
    'unidentified', 'unknown'
    ]

    """
    serotype_rename_map = {
        "H1N1pdm09": "H1N1",
        "H3N2v": "H3N2",
        "H5": "H5NX",
        "H5N1": "H5NX",
        "h5n1": "H5NX",
        "H5N2": "H5NX",
        "H5N2?": "H5NX",
        "H5N3": "H5NX",
        "H5N4": "H5NX",
        "H5N5": "H5NX",
        "H5N6": "H5NX",
        "H5N7": "H5NX",
        "H5N8": "H5NX",
        "H5N9": "H5NX",
        "H5N?": "H5NX",
        "H5Nx": "H5NX",
        "H7": "H7NX",
        "H7N1": "H7NX",
        "H7N2": "H7NX",
        "H7N3": "H7NX",
        "H7N4": "H7NX",
        "H7N5": "H7NX",
        "H7N6": "H7NX",
        "H7N7": "H7NX",
        "H7N8": "H7NX",
        "H7N9": "H7NX",
        "H7Nx": "H7NX",
        "H9": "H9NX",
        "H91": "H9NX",
        "H9N1": "H9NX",
        "H9N1,4": "H9NX",
        "H9N1,9": "H9NX",
        "H9N2": "H9NX",
        "H9N3": "H9NX",
        "H9N4": "H9NX",
        "H9N5": "H9NX",
        "H9N6": "H9NX",
        "H9N7": "H9NX",
        "H9N8": "H9NX",
        "H9N9": "H9NX",
        "Yamagata-like": "B-yam",
        "Yamagata": "B-yam",
        "Victoria": "B-vic",
    }
    df["serotype"] = df["serotype"].replace(serotype_rename_map)

    # Filter on serotypes
    # Only do this for alpha, leave beta alone
    valid_serotypes = ["H1N1", "H3N2", "H5NX", "H7NX", "H9NX", "B-yam", "B-vic"]
    valid = (df["serotype"].isin(valid_serotypes)) | (df["serotype"].isna())
    df.drop(df.index[~valid], inplace=True)

    # Segment extraction
    # Clean segments
    df["genome_coverage"] = df["genome_coverage"].apply(json.loads)
    df["proteins"] = df["genome_coverage"].apply(
        lambda x: (
            [p["name"] for p in x[0]["proteins"]]
            if x[0]["proteins"] is not None
            else []
        )
    )
    segment_df = df[["segments", "proteins"]].explode("proteins")
    segment_df["protein_segment"] = protein_to_segment(segment_df["proteins"])
    segment_df = (
        segment_df.reset_index()
        .groupby("Accession ID")
        .agg(segments=("segments", "first"), protein_segment=("protein_segment", list))
    )
    segment_df["protein_segment"] = segment_df["protein_segment"].apply(
        collapse_segment_list
    )
    segment_map = {
        "NS": 8,
        "M": 7,
        "HA": 4,
        "NP": 8,
        "PB2": 1,
        "PB1": 2,
        "RNA 4": 4,
        "RNA4": 4,
        "PA": 3,
        "segment 7": 7,
        "MA": 7,
        "segment 4": 4,
        "RNA 1": 1,
        "segment 6": 6,
    }
    segment_df["segments"] = segment_df["segments"].apply(
        lambda x: segment_map[x] if x in segment_map else x
    )
    segment_df["segment_complete"] = segment_df["segments"].combine_first(
        segment_df["protein_segment"]
    )
    segment_df.drop(
        segment_df.index[segment_df["segment_complete"].isna()], inplace=True
    )
    segment_df["segment_complete"] = segment_df["segment_complete"].astype(int)

    df = (
        df.drop(columns=["segments"])
        .merge(
            segment_df.drop(columns=["segments", "protein_segment"]),
            how="inner",
            left_index=True,
            right_index=True,
            copy=False,
        )
        .rename(columns={"segment_complete": "segment"})
        .drop(columns=["proteins", "genome_coverage"])
    )

    # Remove "Z" from the end of the submission date string, and convert from
    # ISO datetime to ISO date
    def datetime_to_date(x):
        return datetime.datetime.fromisoformat(x[:-1]).strftime("%Y-%m-%d")

    df.loc[:, "submission_date"] = df["submission_date"].apply(datetime_to_date)
    df.loc[:, "updated_date"] = df["updated_date"].apply(datetime_to_date)

    # Parse location data
    loc_tuples = df["location"].apply(parse_genbank_location)
    # Drop the original "location" column, and join on new columns
    df = df.drop(columns=["location"]).join(
        pd.DataFrame(
            [[a, b] for a, b in loc_tuples.values],
            columns=["division", "location"],
            index=loc_tuples.index,
        )
    )

    # Fill in missing values with empty strings, for remaining
    # metadata columns
    fill_in_cols = [
        "database",
        "host",
        "isolation_source",
        "biosample_accession",
        "authors",
        "publications",
    ]
    for col in fill_in_cols:
        df.loc[:, col] = df[col].fillna("Unknown")

    # Fill in whitespace in the virus name with hyphens
    # (For future minimap2 alignment, query names with spaces
    # are not correctly parsed and are cutoff when written to BAM files)
    df.loc[:, "virus_name"] = df["virus_name"].str.replace(r"\s", "-")

    # Load quality and join to dataframe
    quality_df = pd.read_csv(args.quality, index_col="Accession ID")
    df = df.join(quality_df, how="left")
    df["length"] = df["length"].fillna(0).astype(int)
    # Calculate percent ambiguous, drop the num_ambiguous column
    df["num_ambiguous"] = ((df["num_ambiguous"] / df["length"]) * 100).fillna(0)
    df.rename(columns={"num_ambiguous": "percent_ambiguous"}, inplace=True)

    # Filter out entries without any sequence
    df = df.loc[df["length"] > 0]

    df.to_csv(args.metadata_out)

    # Write summary dataframe
    virus_df = (
        df.reset_index()
        .groupby("virus_name")
        .agg(
            accession_ids=("Accession ID", list),
            segments=("segment", list),
            genus=("genus", "first"),
            serotype=("serotype", "first"),
            region=("region", "first"),
            country=("country", "first"),
            division=("division", "first"),
            location=("location", "first"),
            collection_date=("collection_date", "first"),
            submission_date=("submission_date", "first"),
            updated_date=("updated_date", "first"),
            host=("host", "first"),
            isolation_source=("isolation_source", "first"),
            biosample_accession=("biosample_accession", "first"),
            authors=("authors", "first"),
            publications=("publications", "first"),
        )
        .assign(
            n_segments=lambda x: x["segments"].apply(len),
            has_ha=lambda x: x["segments"].apply(lambda y: 4 in y),
            has_na=lambda x: x["segments"].apply(lambda y: 6 in y),
        )
    )

    # Convert lists to JSON
    virus_df.loc[:, "accession_ids"] = virus_df["accession_ids"].apply(json.dumps)
    virus_df.loc[:, "segments"] = virus_df["segments"].apply(json.dumps)

    virus_df.to_csv(args.metadata_virus_out)


if __name__ == "__main__":
    main()
