#!/bin/bash

export FLASK_APP=flask_server/main.py
export FLASK_ENV=development

flask run --port=5000
