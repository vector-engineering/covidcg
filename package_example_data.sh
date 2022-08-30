#!/bin/bash

set -e

rm -f example_data_genbank/sars2.tar.gz
rm -f example_data_genbank/rsv.tar.gz
rm -f example_data_genbank/flu.tar.gz

cd example_data_genbank
tar -czf sars2.tar.gz sars2/fasta_raw sars2/metadata.csv
tar -czf rsv.tar.gz rsv/fasta_raw rsv/metadata.csv
tar -czf flu.tar.gz flu/fasta_raw flu/metadata.csv