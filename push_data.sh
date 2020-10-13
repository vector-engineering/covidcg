#!/bin/bash

rclone copy data "shared:VE computational/covid_data" -P --transfers=20 --checkers=20

# Push package to google cloud storage
gsutil -h "Cache-Control: no-cache" cp data/data_package.json gs://ve-public/data_package.json
gsutil -h "Cache-Control: no-cache" cp data/data_package.json.gz gs://ve-public/data_package.json.gz
