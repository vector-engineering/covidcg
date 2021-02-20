# coding: utf-8

import os

from yaml import load, dump

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

if os.environ.get("CONFIGFILE", None) is None:
    print("NO CONFIG FILE FOUND")

config = {}

# Load app configuration
with open(os.environ.get("CONFIGFILE"), "r") as fp:
    config = load(fp.read(), Loader=Loader)

# print(config)
