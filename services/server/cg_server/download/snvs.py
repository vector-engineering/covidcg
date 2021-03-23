# coding: utf-8

"""Download selected SNVs in long form

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import pandas as pd

from flask import make_response

from cg_server.query.selection import query_sequences


def download_snvs(conn, req):
    res_df, res_snv = query_sequences(conn, req)
    return make_response(
        res_snv.drop(
            columns=[
                "sequence_id",
                "snp_id",
                "location_id",
                "collection_date",
                "snp_str",
                "color",
            ]
        ).to_csv(index=False),
        200,
        {"Content-Type": "text/csv"},
    )
