#!/usr/bin/bash

gcloud compute scp --recurse ~/covid_ui/dist chena@covid-server:/home/chena/covid_ui/dist
