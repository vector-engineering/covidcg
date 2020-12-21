#!/bin/bash

rclone copy data "mydrive:covid_data" -P --transfers=20 --checkers=20 --exclude "/data/feed.json" --exclude "/data/fasta_temp/**"

# Push package to google cloud storage
gsutil -h "Cache-Control: no-cache" cp data/data_package.json gs://ve-public/data_package.json
gsutil -h "Cache-Control: no-cache" cp data/data_package.json.gz gs://ve-public/data_package.json.gz
gsutil -h "Cache-Control: no-cache" cp data/map_combined_standalone.vg.json gs://ve-public/map_combined_standalone.vg.json

