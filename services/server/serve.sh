#!/bin/bash

export FLASK_APP=cg_server/main.py
export FLASK_ENV=development

flask run --host 0.0.0.0 --port=5000