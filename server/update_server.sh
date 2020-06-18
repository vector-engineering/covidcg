#!/usr/bin/bash

gcloud compute scp --recurse ~/covid_ui/dist chena@covid-server:/var/www/covid
