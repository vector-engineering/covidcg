#!/usr/bin/bash

# Instructions from https://stackoverflow.com/questions/16420092/how-to-make-python-script-run-as-service

systemctl daemon-reload && systemctl enable covid && systemctl start covid --no-block

# To view logs:
systemctl status covid

