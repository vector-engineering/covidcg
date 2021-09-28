from google.cloud import bigquery
from google.cloud import bigtable

# Uploads a csv file from uri using a load job.
def uploadFile(client, uri):
    job_config = bigquery.LoadJobConfig(
        schema=[
            # TODO: Determine schema.
            bigquery.SchemaField("name", "STRING"),
            bigquery.SchemaField("post_abbr", "STRING"),
        ],
        skip_leading_rows=1,
        # The source format defaults to CSV, so the line below is optional.
        source_format=bigquery.SourceFormat.CSV,
    )

    project_id = "test-project-1-2-3"
    instance_id = "covid-cg"
    table_id = "gisaid"

    load_job = client.load_table_from_uri(
        uri, table_id, job_config=job_config
    )  # Make an API request.

    load_job.result()  # Waits for the job to complete.

    destination_table = client.get_table(table_id)  # Make an API request.
    print("Loaded {} rows.".format(destination_table.num_rows))


# TODO: Figure out whether we want to overwrite or delete and write new data.