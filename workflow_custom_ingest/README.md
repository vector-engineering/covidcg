# Custom ingestion for COVID-CG

This processes custom data that has been inserted into the `data_custom/` directory by a user. No chunking on sequence data is performed, it is assumed that the data will be organized by the user of this workflow outside of covidcg. Lineages are assigned to each sequence with [pangolin](https://github.com/cov-lineages/pangolin).

## Running

### Example dataset

An example dataset is included with this code in the folder [test-data/data_custom](test-data/data_custom). This data was derived from a portion of data downloaded from the `workflow_genbank_ingest` version of covidcg (so all sequences and metadata are available on Genbank).

To make use of this data please do the following:

1. Copy test data to the `data_custom/` directory in the covidcg root:

   ```bash
   > cd workflow_custom_ingest
   > cp -r test-data/data_custom ../
   ```

2. Review the configuration settings in [config/config_custom.yaml](../config/config_custom.yaml). In particular, for this test you may want to change `mutation_count_threshold` to `0` (so that covidcg does not remove any low-frequency mutations since the test dataset is very small).

3. Run the workflow in `workflow_custom_ingest` to assign lineages with pangolin and clean up the metadata.

   ```bash
   > snakemake --cores 1 --use-conda
   ```

4. Run main workflow to call mutations/combine all metadata together (see the [Main covidcg instructions](../README.md#main-analysis) for more details).

   ```bash
   > cd ../workflow_main
   > CONFIGFILE=config/config_custom.yaml snakemake --cores 1 --use-conda
   ```

5. Modify `docker-compose.yaml` to mount the produced data folder (as defined in `config_custom.yaml`) to the server container. You'll then be able to start COVID-19 CG and seed your database (see [Main covidcg instructions](../README.md#installation)):

   ```bash
   > docker-compose build
   > docker-compose up -d
   > curl http://localhost:5000/seed
   ```

### Your own data

To run on your own data, please follow the template in the [test-data/data_custom](test-data/data_custom) directory. In particular, copy your gzipped FASTA sequences to `data_custom/fasta_raw` (you can have as many files as you like) and the metadata for all sequences in the `data_custom/metadata_raw.csv` file.

Once you've copied your own data and metadata, you can run with:

```bash
snakemake --cores 1 --use-conda
```

Make sure to include `--use-conda` so that `pangolin` gets automatically installed from `envs/pangolin.yaml`

## Configuration

All configuration options, and their descriptions, are available in the `config/config_custom.yaml` file.

Metadata columns (`metadata_cols`) and sequence groupings (`group_cols`) specific to this pipeline are defined in the `config/config_custom.yaml` file.

## Metadata Requirements

The following fields are **required** by the downstream `workflow_main`:

- `Accession ID`
- `submission_date`
- `collection_date`
- `region`
- `country`
- `division`
- `location`

Sequences without `submission_date`, `collection_date`, or `region` are filtered out.

For more granular location metadata (`country`, `division`, `location`), missing or undetermined values are replaced by the integer -1. This is to facilitate easier groupby-aggregate operations downstream.

Please see [test-data/data_custom/metadata_raw.csv](test-data/data_custom/metadata_raw.csv) for an example.

## Acknowledgements

The code for this method of ingestion of data is derived from the `workflow_genbank_ingest` code.

The `pangolin` lineage assignment tool is hosted on GitHub: [https://github.com/cov-lineages/pangolin](https://github.com/cov-lineages/pangolin).
