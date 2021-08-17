#!/bin/bash

export FLASK_APP=cg_server/main.py
export FLASK_ENV=development
# export CONFIGFILE=/path/to/configfile.yml # Set config file here

flask run --host 0.0.0.0 --port=5000