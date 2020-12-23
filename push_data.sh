#!/bin/bash

rclone copy data "mydrive:covid_data" --transfers=10 -P

# Push package to google cloud storage
gsutil -h "Cache-Control: no-cache" cp data/data_package.json.gz gs://ve-public/data_package.json.gz
gsutil -h "Cache-Control: no-cache" cp data/map_combined_standalone.vg.json gs://ve-public/map_combined_standalone.vg.json

# Push genbank data
rclone copy data_genbank "mydrive:covid_data/genbank" --transfers 10 -P

# Push genbank packages
gsutil -h "Cache-Control: no-cache" cp data_genbank/data_package.json.gz gs://ve-public/genbank/data_package.json.gz
gsutil -h "Cache-Control: no-cache" cp data_genbank/map_combined_standalone.vg.json gs://ve-public/genbank/map_combined_standalone.vg.json