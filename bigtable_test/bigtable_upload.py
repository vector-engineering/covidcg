
import datetime
import time
import math

from google.cloud import bigtable
from google.cloud.bigtable import column_family
from google.cloud.bigtable import row_filters
# [END bigtable_hw_imports]
import pandas as pd

import multiprocessing as mp
from functools import partial


family_map = {'metadata': [
    'collection_date', 'submission_date', 'Accession ID', 'location', 'location_id']}

def upload_df_to_bigtable(ind, df):
    project_id = 'test-project-1-2-3'
    instance_id = 'covid-cg'
    table_id = 'test'

    client = bigtable.Client(project=project_id, admin=True)
    instance = client.instance(instance_id)
    table = instance.table(table_id)


    batcher = table.mutations_batcher(flush_count=100000, max_row_bytes=50242880)
    rows = []
    column_family_id = 'metadata'
    start = time.time()
    for i, dfrow in df.iterrows():
        row_key = '{}'.format(dfrow['Accession ID']).encode('utf-8')

        row = table.direct_row(row_key)
        row.set_cell(column_family_id,'dummy',1,timestamp=datetime.datetime.utcnow())

        for mut in dfrow['dna_snp_str']:
            column = mut.encode('utf-8')
            row.set_cell(column_family_id,
                        column,
                        1,
                        timestamp=datetime.datetime.utcnow())
        
        if i > 0 and i % 10000 == 0:
            print(i)
        
        rows.append(row)

    batcher.mutate_rows(rows)
    batcher.flush()
    print('{} time: {}'.format(ind, time.process_time() - start))

    return 1


# def main(project_id, instance_id, table_id):
def main():

    start = time.time()

    project_id = 'test-project-1-2-3'
    instance_id = 'covid-cg'
    table_id = 'test'
    
    df = pd.read_json('test_mutations.csv')
    # TODO: Test raw JSON .loads

    """
    
    # [START bigtable_hw_connect]
    # The client must be created with admin=True because it will create a
    # table.
    client = bigtable.Client(project=project_id, admin=True)
    instance = client.instance(instance_id)
    # [END bigtable_hw_connect]
    # [START bigtable_hw_create_table]
    print('Creating the {} table.'.format(table_id))
    table = instance.table(table_id)
    # batcher = table.mutations_batcher(flush_count=100000, max_row_bytes=50242880)

    print('Creating column family cf1 with Max Version GC rule...')
    # Create a column family with GC policy : most recent N versions
    # Define the GC policy to retain only the most recent 2 versions
    max_versions_rule = column_family.MaxVersionsGCRule(1)
    column_family_id = 'metadata'
    column_families = {column_family_id: max_versions_rule}
    if not table.exists():
        table.create(column_families=column_families)
    else:
        print("Table {} already exists.".format(table_id))
    # [END bigtable_hw_create_table]
    """

    # [START bigtable_hw_write_rows]
    print('Writing some data to the table.')
    #greetings = ['Hello World!', 'Hello Cloud Bigtable!', 'Hello Python!']
    rows = []
    # column = 'greeting'.encode()
    # for i, value in enumerate(greetings):
    start = time.process_time()

    ## MULTIPROCESSING
    n_processes = 3
    chunk_size = 5000
    n_chunks = math.ceil(len(df) / chunk_size)

    ## BREAK UP DATAFRAME INTO CHUNKS
    df_chunks = []
    for i in range(n_chunks):
        start = i * chunk_size
        end = start + chunk_size
        if end > (len(df) - 1):
            end = len(df) - 1
        df_chunks.append(df.iloc[start:end])

    with mp.get_context("spawn").Pool(processes=n_processes) as pool:

        def process_callback(ret):
            pass

        def handle_error(e):
            raise e

        for i, chunk in enumerate(df_chunks):
            pool.apply_async(
                partial(upload_df_to_bigtable, i, chunk),
                callback=process_callback,
                error_callback=handle_error,
            )

        pool.close()
        pool.join()
    

    """
    for i, dfrow in df.iterrows():

        # if i > 10000:
        #     break

        # Note: This example uses sequential numeric IDs for simplicity,
        # but this can result in poor performance in a production
        # application.  Since rows are stored in sorted order by key,
        # sequential keys can result in poor distribution of operations
        # across nodes.
        #
        # For more information about how to design a Bigtable schema for
        # the best performance, see the documentation:
        #
        #     https://cloud.google.com/bigtable/docs/schema-design
        row_key = '{}'.format(dfrow['Accession ID']).encode('utf-8')

        row = table.direct_row(row_key)
        row.set_cell(column_family_id,'dummy',1,timestamp=datetime.datetime.utcnow())

        for mut in dfrow['dna_snp_str']:
            column = mut.encode('utf-8')
            row.set_cell(column_family_id,
                        column,
                        1,
                        timestamp=datetime.datetime.utcnow())
        
        if i > 0 and i % 10000 == 0:
            print(i)
        
        rows.append(row)
    batcher.mutate_rows(rows)
    batcher.flush()
    print('time: {}'.format(time.process_time() - start))
    # [END bigtable_hw_write_rows]
    """

    project_id = 'test-project-1-2-3'
    instance_id = 'covid-cg'
    table_id = 'test'

    client = bigtable.Client(project=project_id, admin=True)
    instance = client.instance(instance_id)
    table = instance.table(table_id)

    
    # [START bigtable_hw_create_filter]
    # Create a filter to only retrieve the most recent version of the cell
    # for each column accross entire row.
    row_filter = row_filters.CellsColumnLimitFilter(1)
    # [END bigtable_hw_create_filter]

    # [START bigtable_hw_get_with_filter]
    print('Getting a single greeting by row key.')
    key = 'MN996528'.encode()

    row = table.read_row(key, row_filter)
    print(row)
    # cell = row.cells[column_family_id][column][0]
    # print(cell.value.decode('utf-8'))
    # [END bigtable_hw_get_with_filter]

    # [START bigtable_hw_scan_with_filter]
    print('Scanning for all greetings:')
    partial_rows = table.read_rows(limit=100, filter_=row_filter)

    for row in partial_rows:
        print(row)
        #cell = row.cells[column_family_id][column][0]
        #print(cell.value.decode('utf-8'))
    # [END bigtable_hw_scan_with_filter]
    

    # [START bigtable_hw_delete_table]
    # print('Deleting the {} table.'.format(table_id))
    # table.delete()
    # [END bigtable_hw_delete_table]

    print('total runtime: {:.3} ms'.format(time.time() - start))


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(
    #     description=__doc__,
    #     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # parser.add_argument('project_id', help='Your Cloud Platform project ID.')
    # parser.add_argument(
    #     'instance_id', help='ID of the Cloud Bigtable instance to connect to.')
    # parser.add_argument(
    #     '--table',
    #     help='Table to create and destroy.',
    #     default='Hello-Bigtable')

    # args = parser.parse_args()
    # main(args.project_id, args.instance_id, args.table)
    main()
