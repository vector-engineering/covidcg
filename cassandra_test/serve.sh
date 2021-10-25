#!/bin/bash

export FLASK_APP=cg_server/main.py
export FLASK_ENV=development
export CONFIGFILE=/Users/chena/covidcg/config/config_genbank_example.yaml
# export CONFIGFILE=/path/to/configfile.yml # Set config file here

# export DATA_PATH=... # Optional - defaults to data path in config file
# export STATIC_DATA_PATH=... # Optional - defaults to static data path in config file

# Cassandra config
export CASSANDRA_USER=cassandra
export CASSANDRA_PASSWORD=cassandra
export CASSANDRA_KEYSPACE=cg_gisaid
export CASSANDRA_HOST=127.0.0.1
export CASSANDRA_PORT=9042

flask run --host 0.0.0.0 --port=5001