# coding: utf-8

"""Configuration settings import

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import os

from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


config = {}

# Load app configuration

# Print warning if configfile is not defined
# Only an issue in non-docker deployments
if "CONFIGFILE" not in os.environ:
    print(
        'CONFIGFILE environment variable not defined. Defaulting config file path to "/opt/config.yaml"'
    )

config_file_path = os.getenv("CONFIGFILE", "/opt/config.yaml")
with open(config_file_path, "r") as fp:
    config = load(fp.read(), Loader=Loader)

# print(config)
