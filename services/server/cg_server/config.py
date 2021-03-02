# coding: utf-8

import os

from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


config = {}

# Load app configuration
with open("/opt/config.yaml", "r") as fp:
    config = load(fp.read(), Loader=Loader)

# print(config)
