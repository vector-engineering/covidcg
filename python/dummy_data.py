#!/usr/bin/env python3
# coding: utf-8

'''Construct dummy testing data
Author: Albert Chen (Deverman Lab, Broad Institute)
'''

import pandas as pd
import networkx as nx
import numpy as np
import random
import tqdm

from datetime import date, timedelta
from scipy.stats import uniform, expon, norm

from process_clades import load_clades, build_clade_tree, build_clade_sequences
from process_geo import load_geo_data, build_geo_graph

start_date = date.fromisoformat('2020-02-01')
end_date = date.fromisoformat('2020-05-10')

def main():
    clades_df = load_clades()
    clade_tree = build_clade_tree(clades_df)
    clades = [n for n in list(clade_tree.nodes) if n != 'root']

    clade_data_df = pd.read_csv('processed_data/clade_data.csv')

    clade_seqs, clade_positions = build_clade_sequences(clades_df, clade_tree)

    location_df = load_geo_data()
    loc_tree = build_geo_graph(location_df)

    # Output dataframe
    case_df = pd.DataFrame()

    # For each location:
    for i in tqdm.tqdm(range(len(location_df))):
    #for i in range(1):
        # Choose how many clades we want
        # Uniform(3, 12)
        n_clades = random.randint(3, 12)

        # Choose a random subset of the clades
        loc_clades = random.sample(clades, n_clades)
        # print(loc_clades)

        # Choose a random start date, offset from the start_date
        # by N days. Uniform(0, 30)
        day_offset = random.randint(0, 30)
        loc_start_date = start_date + timedelta(days=day_offset)
        num_days = int((end_date - loc_start_date).days)
        loc_date_range = [loc_start_date + timedelta(days=i) for i in range(num_days)]

        # Choose growth rates for each of the clades
        # Uniform(1, 1.1)
        growth_rates = uniform.rvs(loc=0, scale=5, size=n_clades)
        # Exp(0, 0.01)
        # growth_rates = 1 + expon.rvs(loc=0, scale=0.1, size=n_clades)
        # print(growth_rates)

        # Each clade starts out with 1
        initial_num = [1] * n_clades
        
        # From the start until the end date, increase the number of
        # each clade by their respective growth rate
        # Rows = clades, Cols = days
        
        case_mat = np.zeros((n_clades, num_days))
        case_mat[:, 0] = initial_num
        for j in range(1, num_days):
            # Linear growth
            case_mat[:, j] = case_mat[:, j - 1] + growth_rates
            # Add some noise
            case_mat[:, j] += norm.rvs(loc=0, scale=1)
            # Make sure it's above 0
            case_mat[case_mat[:, j] < 0, j] = 0

        # To integers
        case_mat = case_mat.astype(int)
        # print(case_mat)

        # Create dataframe for this location
        loc_df = pd.DataFrame({
            'location_id': location_df.at[i, 'index'],
            'clade_id': np.repeat(loc_clades, num_days),
            'date': [d.isoformat() for d in loc_date_range] * n_clades,
            'cases': case_mat.flatten() # Flatten row-wise
        })
        # Map clade names to ids
        loc_df['clade_id'] = loc_df['clade_id'].map(pd.Series(
            clade_data_df['index'].values,
            index=clade_data_df['clade'].values
        ))

        # Randomly drop 50% of the rows
        #keep_rows = random.sample(range(1, len(loc_df)), int(len(loc_df) * 0.5))
        #loc_df = loc_df.loc[keep_rows, :]

        # Append to master dataframe
        case_df = case_df.append(loc_df, ignore_index=True)

    print(case_df)
    case_df.to_csv('processed_data/simulated_case_data.csv', index=False)
    case_df.to_json('processed_data/simulated_case_data.json', orient='records')

    '''
    # Create dataframe of positions for the entropy graph
    entropy_df = case_df.copy()
    entropy_df['position'] = case_df['clade'].map(clade_positions)
    # Unnest (explode) the positions column
    entropy_df = entropy_df.explode(column='position')
    print(entropy_df)
    entropy_df.to_csv('processed_data/simulated_entropy_data.csv', index=False)
    entropy_df.to_json('processed_data/simulated_entropy_data.json', orient='records')
    '''
        

if __name__ == '__main__':
    main()
