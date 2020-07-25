# coding: utf-8

"""Clean sequencing technology metadata

Author: Albert Chen (Deverman Lab - Broad Institute)
"""

import numpy as np
import pandas as pd
import re


def clean_seq_tech_metadata(seq_meta_df):
    """Clean "Sequencing Technology" values
    """

    print("Cleaning sequencing technology metadata...", end="", flush=True)

    # Basic cleaning
    seq_meta_df["sequencing_tech"] = seq_meta_df["Sequencing technology"].str.strip()

    replace_map = [
        (r"illumina", "Illumina", False),
        (r"Ilumina", "Illumina"),
        (r"^llumina", "Illumina"),
        (r"Illumina\'s", "Illumina"),
        (r"Illumina_", "Illumina "),
        (r"Illumina technology", "Illumina"),
        (r"Immumina", "Illumina"),
        (r"iSeq([0-9]+)", lambda m: "iSeq {}".format(m.groups()[0])),
        (r"hiseq", "HiSeq", False),
        (r"miseq", "MiSeq", False),
        (r"MiSeq\.", "MiSeq"),
        (r"^MiSeq", "Illumina MiSeq"),
        (r"miniseq", "MiniSeq", False),
        (r"nextseq", "NextSeq", False),
        (r"Next\sSeq", "NextSeq"),
        (r"NextSeq([0-9]+)", lambda m: "NextSeq {}".format(m.groups()[0])),
        (r"NextSeq 5[01]{1}[0-9]{1}", "NextSeq 500"),
        (r"^NextSeq", "Illumina NextSeq"),
        (r"novaseq", "NovaSeq", False),
        (r"Noveseq", "NovaSeq"),
        (r"NovaSeq([0-9+])", lambda m: "NovaSeq {}".format(m.groups()[0])),
        (
            r"^(NovaSeq|MiSeq|NextSeq|iSeq|HiSeq|MiniSeq)",
            lambda m: "Illumina {}".format(m.groups()[0]),
        ),
        (r"OXFORD_NANOPORE", "Nanopore"),
        (r"nanopore", "Nanopore", False),
        (r"Nanpore", "Nanopore"),
        (r"Nanopre", "Nanopore"),
        (r"Nanopore technology", "Nanopore"),
        (r"oxford nanopore technology", "Nanopore", False),
        (r"Oxoford", "Oxford"),
        (r"Oxford Nanopore", "Nanopore"),
        (r"minion", "MinION", False),
        (r"MinION[,\.]$", "MinION"),
        (r"GRIDION", "GridION"),
        (r"^(MinION|GridION)", lambda m: "Nanopore {}".format(m.groups()[0])),
        (r"IonTorrent", "Ion Torrent", False),
        (r"IonTorren", "Ion Torrent"),
        (r"ion torrent", "Ion Torrent", False),
        (r"Ion Torren", "Ion Torrent"),
        (r"Ion Torrentt", "Ion Torrent"),
        (r"artic", "ARTIC", False),
        (r"ARTIC\sprotocol", "ARTIC sequencing protocol"),
        (r"ARTIC\sProtocol", "ARTIC sequencing protocol"),
        (r"\-\-\sARTIC", "- ARTIC"),
        (r"MGISEQ-?([0-9]+)", lambda m: "MGISEQ {}".format(m.groups()[0])),
        (r"Sanger dideoxy sequencing", "Sanger"),
        (r"Sanger Sequencing Method", "Sanger"),
        (r", assembled sequences", ""),
        (r";", ","),
        (r"\sand\s", ", "),
        (r"\s\&\s", ", "),
        (r"+", ", "),
        (r",\s+", ", "),
        (r"_", " "),
    ]

    for pair in replace_map:
        seq_meta_df["sequencing_tech"] = seq_meta_df["sequencing_tech"].str.replace(
            pair[0], pair[1], -1, pair[2] if len(pair) > 2 else True
        )

    seq_key_map = {
        "DNBSEQ-T7": [],
        "DNBSEQ-G400": ["MGI Tech. DNBSEQ-G400"],
        "DNBSEQ-G400RS": ["DNBSEQ-G400RS, MGI Tech Co., Ltd"],
        "Illumina, Nanopore": [],
        "Nanopore GridION": ["Nanopore - GridION"],
        "Nanopore MinION - ARTIC sequencing protocol": [
            "Nanopore MinION ARTIC sequencing protocol"
        ],
        "Nanopore MinION ARTIC Network V1": ["Nanopore MinION ARTIC v1 primers"],
        "Nanopore, Sanger": [],
        "Sanger": ["Sanger dideoxy sequencing"],
        "Unspecified NGS": ["NGS"],
        "Unknown": ["unknown", "Oro-pharyngeal swab", "unkown", "contacted submitter"],
    }

    seq_map = {}
    for k, v in seq_key_map.items():
        # Add self
        seq_map[k] = k
        for _v in v:
            seq_map[_v] = k

    seq_meta_df["sequencing_tech"] = (
        (seq_meta_df["sequencing_tech"].map(seq_map))
        .combine_first(seq_meta_df["sequencing_tech"])
        .fillna("Unknown")
    )

    print("done")

    return seq_meta_df


def clean_assembly_metadata(seq_meta_df):
    """Clean "Assembly Method" column
    """

    print("Cleaning assembly method metadata...", end="", flush=True)

    seq_meta_df["assembly_method"] = seq_meta_df["Assembly method"].str.strip()

    replace_map = [
        # Aliases
        (r"artic", "ARTIC", False),
        (r"ArticNetwork", "ARTIC Network"),
        (r"bcftools", "BCFtools", False),
        (r"bowtie", "bowtie", False),
        (r"bwa", "BWA", False),
        (r"bwa[-|\s]mem", "BWA-MEM", False),
        (r"(Burrows-Wheeler Aligner method)", "", False),
        (r"Burrows-Wheeler Aligner Tool (BWA)", "BWA", False),
        (r"custom", "Custom", False),
        (r"dragen", "DRAGEN", False),
        (r"geneious prime", "Geneious Prime", False),
        (r"megahit", "MEGAHIT", False),
        (r"minimap2", "minimap2", False),
        (r"mpileup", "mpileup", False),
        (r"samtools", "samtools", False),
        (r"PAdes", "SPAdes"),
        (r"SSPAdes", "SPAdes"),
        (r"spades", "SPAdes", False),
        (r"vcftools", "vcftools", False),
        (r"Workbench ([0-9\.]+)", lambda m: "Workbench v{}".format(m.groups()[0])),
        (r"v. ([0-9\.]+)", lambda m: "v{}".format(m.groups()[0])),
        ("\ufeff", ""),
        (r" assembly method", ""),
        # Separators
        (r";", ","),
        (r"\sand\s", ", "),
        (r"\s\&\s", ", "),
        (r"+", ", "),
        (r",\s+", ", "),
        (r"--", "-"),
        (r"\s\/\s", " - "),
        # These are Excel errors
        (r"Sequencher 5.4.[789]", "Sequencher 5.4.6"),
        (r"Sequencher 5.4.[123]{1}[0-9]{1}", "Sequencher 5.4.6"),
        (r"Sequencher 5.4.4[12]", "Sequencher 5.4.40"),
        # More Excel errors
        (r"mapped to NC_045512v[3-9]", "mapped to NC_045512v2"),
        (r"mapped to NC_045512v[0-9]{1}[0-9]{1}", "mapped to NC_045512v2"),
        (r"IRMA v1.0.[1-9]", "IRMA v1.0.0"),
        (r"IRMA v1.0.[0-9]{1}[0-9]{1}", "IRMA v1.0.0"),
        (r"IME-BJ0[2-9]", "IME-BJ01"),
        (r"Geneious Prime 2020.1.[3-9]{1}", "Geneious Prime 2020.1.2"),
        (r"Geneious Prime 2020.1.[0-9]{1}[0-9]{1}", "Geneious Prime 2020.1.2"),
        (r"Geneious Prime 2020.0.[7-9]", "Geneious Prime 2020.0.6"),
        (r"bowtie[3-9]", "bowtie2"),
        (r"bowtie[0-9]{1}[0-9]{1}", "bowtie2"),
        (r"minimap[3-9]", "minimap2"),
        (r"minimap[0-9]{1}[0-9]{1}", "minimap2"),
    ]

    for pair in replace_map:
        seq_meta_df["assembly_method"] = seq_meta_df["assembly_method"].str.replace(
            pair[0], pair[1], -1, pair[2] if len(pair) > 2 else True
        )

    method_key_map = {"Unknown": ["unknown"]}

    method_map = {}
    for k, v in method_key_map.items():
        # Add self
        method_map[k] = k
        for _v in v:
            method_map[_v] = k

    seq_meta_df["assembly_method"] = (
        (seq_meta_df["assembly_method"].map(method_map))
        .combine_first(seq_meta_df["assembly_method"])
        .fillna("Unknown")
    )

    print("done")

    return seq_meta_df


def clean_comment_type_metadata(seq_meta_df):
    seq_meta_df["comment_type"] = seq_meta_df["Comment type"].str.strip()

    seq_meta_df["comment_type"] = seq_meta_df["comment_type"].fillna("Unknown")

    return seq_meta_df


def clean_seq_metadata(seq_meta_df):

    seq_meta_df = clean_seq_tech_metadata(seq_meta_df)
    seq_meta_df = clean_assembly_metadata(seq_meta_df)
    seq_meta_df = clean_comment_type_metadata(seq_meta_df)

    # Comment type is already clean - no need to do anything else

    # Take subset of columns
    seq_meta_df = seq_meta_df[["sequencing_tech", "assembly_method", "comment_type"]]

    return seq_meta_df

