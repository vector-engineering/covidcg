#!/bin/bash

DATA_FOLDER="../data_genbank2"

# For development: clean up files for re-ingest, without re-downloading feed

rm -rf ${DATA_FOLDER}/fasta_raw
rm -rf ${DATA_FOLDER}/fasta_temp
rm -rf ${DATA_FOLDER}/lineages
rm -f ${DATA_FOLDER}/status/merge_sequences*
rm -f ${DATA_FOLDER}/metadata.csv