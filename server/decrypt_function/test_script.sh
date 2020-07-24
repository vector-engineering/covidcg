#!/bin/bash

functions-framework --target=decrypt_accession_id --debug

# TEST REQUEST:

# Type: POST
# Headers: Content-Type: application/json
# Body:
# {
# 	"accession_ids": [
#       "62da3a1bdaec4495d689f7211287d70e2801cd5e899ea0afd84c0fdb20b92a8c", # 
#       "bb1c5afb8a003a053f5668bfec64710e32ef7559e0d8d929da5856be8c9014ff",
# 	    "asdf"
#   ]
# }
# Should return:
# {"accession_ids": ["EPI_ISL_402120", "EPI_ISL_402123", null]}