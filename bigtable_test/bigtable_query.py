from google.cloud import bigtable
from google.cloud.bigtable import column_family
from google.cloud.bigtable import row_filters

import time
from collections import defaultdict


def main():

    project_id = 'test-project-1-2-3'
    instance_id = 'covid-cg'
    table_id = 'test'

    client = bigtable.Client(project=project_id, read_only=True)
    instance = client.instance(instance_id)
    table = instance.table(table_id)

    start = time.time()

    
    # [START bigtable_hw_create_filter]
    # Create a filter to only retrieve the most recent version of the cell
    # for each column accross entire row.
    row_filter = row_filters.CellsColumnLimitFilter(1)
    # [END bigtable_hw_create_filter]

    # print('Getting a single greeting by row key.')
    # key = 'MN996528'.encode()
    # row = table.read_row(key, row_filter)
    # print(row.cells)

    # [START bigtable_hw_scan_with_filter]
    print('Scanning for all greetings:')
    n_rows = 10000
    partial_rows = table.read_rows(limit=n_rows, filter_=row_filter)
    print('fetched {} rows in {:.3f} ms'.format(n_rows, (time.time() - start)*1000))

    dna_mut_accumulator = defaultdict(lambda: 0)
    for row in partial_rows:
        for k, v in row.cells['metadata'].items():
            if k.decode('utf-8') == 'dummy':
                continue
            
            dna_mut_accumulator[k.decode('utf-8')] += 1

    #print(dna_mut_accumulator)
        #cell = row.cells[column_family_id][column][0]
        #print(cell.value.decode('utf-8'))
    # [END bigtable_hw_scan_with_filter]
    

    # [START bigtable_hw_delete_table]
    # print('Deleting the {} table.'.format(table_id))
    # table.delete()
    # [END bigtable_hw_delete_table]

    print('accumulating {} rows | total runtime: {:.3f} ms'.format(n_rows, (time.time() - start)*1000))


if __name__ == '__main__':
    main()