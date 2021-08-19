#!/bin/bash

export FLASK_APP=cg_server/main.py
export FLASK_ENV=development
# export CONFIGFILE=/path/to/configfile.yml # Set config file here

# export DATA_PATH=... # Optional - defaults to data path in config file
# export STATIC_DATA_PATH=... # Optional - defaults to static data path in config file

# POSTGRES CONFIG
export POSTGRES_USER=postgres
export POSTGRES_PASSWORD=cg
export POSTGRES_DB=cg_dev
export POSTGRES_HOST=127.0.0.1
export POSTGRES_PORT=5432
export POSTGRES_MAX_CONN=20

flask run --host 0.0.0.0 --port=5001