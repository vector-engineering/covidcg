#!/bin/bash

FLASK_APP=cg_server/main.py
FLASK_ENV=development

flask run --host 0.0.0.0 --port=5000