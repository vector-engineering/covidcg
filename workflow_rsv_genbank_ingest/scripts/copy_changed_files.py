# coding: utf-8

"""
Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import re
import shutil

from pathlib import Path


def copy_changed_files(files, data_folder):
    """Using the `get_changed_chunks` function, only copy fasta files which have changed
    from the purgatory `fasta_temp` folder to the `fasta_raw` folder. By copying over the files,
    it will flag to snakemake that they (and only they - not the others) will need to be
    reprocessed and realigned.

    Parameters
    ----------
    files: list of str
    data_folder: str

    Returns
    -------
    None
    """
    # Make the fasta_raw folder if it doesn't already exist yet
    # snakemake won't do this automatically since no fasta_raw files
    # are explicitly defined in any output
    (Path(data_folder) / "fasta_raw").mkdir(exist_ok=True)

    # For each changed chunk (as defined by `get_changed_chunks`),
    # copy over the fasta file from the `fasta_temp` folder to the `fasta_raw` folder
    for chunk in files:
        chunk = re.sub(r"\.fa\.gz$", "", Path(chunk).name)
        fasta_temp_path = Path(data_folder) / "fasta_temp" / (chunk + ".fa.gz")
        fasta_raw_path = Path(data_folder) / "fasta_raw" / (chunk + ".fa.gz")

        shutil.copyfile(fasta_temp_path, fasta_raw_path)
        shutil.copystat(fasta_temp_path, fasta_raw_path)
