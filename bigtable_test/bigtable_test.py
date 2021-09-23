# pip install google-cloud-bigtable

import argparse

from google.cloud import bigtable


def main():

    project_id="test-project-1-2-3"
    instance_id="covid-cg"
    table_id="gisaid"

    # Create a Cloud Bigtable client.
    client = bigtable.Client(project=project_id)

    # Connect to an existing Cloud Bigtable instance.
    instance = client.instance(instance_id)

    # Open an existing table.
    table = instance.table(table_id)

    row_key = 'r1'
    row = table.read_row(row_key.encode('utf-8'))

    column_family_id = 'metadata'
    column_id = 'c1'.encode('utf-8')
    value = row.cells[column_family_id][column_id][0].value.decode('utf-8')

    print('Row key: {}\nData: {}'.format(row_key, value))


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(
    #     description=__doc__,
    #     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('project_id', help='Your Cloud Platform project ID.')
    # parser.add_argument(
    #     'instance_id', help='ID of the Cloud Bigtable instance to connect to.')
    # parser.add_argument(
    #     '--table',
    #     help='Existing table used in the quickstart.',
    #     default='my-table')

    # args = parser.parse_args()
    # main(args.project_id, args.instance_id, args.table)
    main()