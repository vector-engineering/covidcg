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
    request_json = request.get_json(silent=True)
    request_args = request.args

    # print("json:", request_json)
    # print("args:", request_args)

    if not request_json:
        return "Error: No JSON data provided", 400, {"ContentType": "text/plain"}

    if "accession_ids" not in request_json:
        return (
            "Error: No Accession IDs provided in JSON data",
            400,
            {"ContentType": "text/plain"},
        )

    hashed_ids = request_json["accession_ids"]
    accession_ids = []
    for hashed_id in hashed_ids:
        if hashed_id not in accession_id_map:
            accession_ids.append(None)
        else:
            accession_ids.append(accession_id_map[hashed_id])

    return (
        json.dumps({"accession_ids": accession_ids}),
        200,
        {"ContentType": "application/json"},
    )

