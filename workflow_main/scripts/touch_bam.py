#!/usr/bin/env python3
# coding: utf-8

"""Migrate BAM files over after a reprocess step,
but only if the underlying fasta entries did not change
"""

import gzip
import shutil
from pathlib import Path

data_folder = Path("../data/")

old_fasta = data_folder / "fasta_processed_backup"
new_fasta = data_folder / "fasta_processed"

old_seqs = {}
for f in old_fasta.iterdir():
    old_seqs[f.name] = []
    with gzip.open(str(f), "rt") as fp:
        for line in fp:
            line = line.strip()
            if line and line[0] == ">":
                old_seqs[f.name].append(line[1:])


new_seqs = {}
for f in new_fasta.iterdir():
    new_seqs[f.name] = []
    with gzip.open(str(f), "rt") as fp:
        for line in fp:
            line = line.strip()
            if line and line[0] == ">":
                new_seqs[f.name].append(line[1:])


same = []
for k in new_seqs.keys():
    if k not in old_seqs:
        continue

    old = set(old_seqs[k])
    new = set(new_seqs[k])

    if old == new:
        same.append(k)

print(len(same))

old_bam_path = Path("../data/bam_old")
new_bam_path = Path("../data/bam")
for f in same:
    f = f.replace(".fa.gz", "") + ".bam"
    if not (old_bam_path / f).exists():
        continue
    shutil.copy2(str(old_bam_path / f), str(new_bam_path / f))
    (new_bam_path / f).touch()  # Update mtime for snakemake

