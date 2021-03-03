# coding: utf-8

import os

from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


config = {}

# Load app configuration
config_file_path = os.getenv("CONFIGFILE", "/opt/config.yaml")
with open(config_file_path, "r") as fp:
    config = load(fp.read(), Loader=Loader)

# print(config)
