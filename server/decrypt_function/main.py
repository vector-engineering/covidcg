import csv
import io
import json
from google.cloud import storage

client = storage.Client()
DATA_BUCKET = "covid-cg"
bucket = client.get_bucket(DATA_BUCKET)

blob = bucket.get_blob("accession_hashmap.csv")
csvdata = blob.download_as_string().decode("utf-8")

accession_id_map = {}
csvreader = csv.reader(csvdata.splitlines(), delimiter=",", quotechar='"')
for row in csvreader:
    # Map hash --> real Accession ID
    accession_id_map[row[1]] = row[0]


def decrypt_accession_id(request):
    """HTTP Cloud Function.
    Args:
        request (flask.Request): The request object.
        <http://flask.pocoo.org/docs/1.0/api/#flask.Request>
    Returns:
        The response text, or any set of values that can be turned into a
        Response object using `make_response`
        <http://flask.pocoo.org/docs/1.0/api/#flask.Flask.make_response>.
    """

    # For more information about CORS and CORS preflight requests, see
    # https://developer.mozilla.org/en-US/docs/Glossary/Preflight_request
    # for more information.

    # Set CORS headers for the preflight request
    if request.method == "OPTIONS":
        # Allows GET requests from any origin with the Content-Type
        # header and caches preflight response for an 3600s
        headers = {
            "Access-Control-Allow-Origin": "*",
            "Access-Control-Allow-Methods": "GET",
            "Access-Control-Allow-Headers": "Content-Type",
            "Access-Control-Max-Age": "3600",
        }

        return ("", 204, headers)

    # Set CORS headers for the main request
    headers = {"Access-Control-Allow-Origin": "*"}

    request_json = request.get_json(silent=True)
    request_args = request.args

    # print("json:", request_json)
    # print("args:", request_args)

    if not request_json:
        headers["ContentType"] = "text/plain"
        return "Error: No JSON data provided", 400, headers

    if "accession_ids" not in request_json:
        headers["ContentType"] = "text/plain"
        return ("Error: No Accession IDs provided in JSON data", 400, headers)

    hashed_ids = request_json["accession_ids"]
    accession_ids = []
    for hashed_id in hashed_ids:
        if hashed_id not in accession_id_map:
            accession_ids.append(None)
        else:
            accession_ids.append(accession_id_map[hashed_id])

    headers["ContentType"] = "application/json"
    return (json.dumps({"accession_ids": accession_ids}), 200, headers)

