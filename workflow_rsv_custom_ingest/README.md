# Custom ingestion for RSVG

This processes custom data that has been inserted into the `data_custom/` directory by a user. No chunking on sequence data is performed, it is assumed that the data will be organized by the user of this workflow outside of rsvg.

## Running

### Example dataset

An example dataset is included with this code in the folder [test-data/data_custom](test-data/data_custom). This data was derived from a portion of data downloaded from the `workflow_rsv_genbank_ingest` version of rsvg (so all sequences and metadata are available on Genbank).

To make use of this data please do the following:

1. Copy test data to the `data_custom/rsv` directory in the repository root:

   ```bash
   > cd workflow_rsv_custom_ingest
   > cp -r test-data/data_custom ../
   ```

2. Review the configuration settings in [config/config_rsv_custom.yaml](../config/config_rsv_custom.yaml). In particular, for this test you may want to change `mutation_count_threshold` to `0` (so that rsvg does not remove any low-frequency mutations since the test dataset is very small).

3. Run the workflow in `workflow_custom_ingest` to clean up the metadata.

   ```bash
   > snakemake --cores 1
   ```

4. Run main workflow to call mutations/combine all metadata together (see the [Main covidcg instructions](../README.md#main-analysis) for more details).

   ```bash
   > cd ../workflow_rsv_main
   > CONFIGFILE=config/config_rsv_custom.yaml snakemake --cores 1
   ```

5. Modify `docker-compose.yaml` to mount the produced data folder (as defined in `config_rsv_custom.yaml`) to the server container. You'll then be able to start COVID-19 CG and seed your database (see [Main covidcg instructions](../README.md#installation)):

   ```bash
   > docker-compose build
   > docker-compose up -d
   > curl http://localhost:5000/seed
   ```

### Your own data

To run on your own data, please follow the template in the [test-data/data_custom](test-data/data_custom) directory. In particular, copy your gzipped FASTA sequences to `data_custom/fasta_raw` (you can have as many files as you like) and the metadata for all sequences in the `data_custom/metadata_raw.csv` file.

Once you've copied your own data and metadata, you can run with:

```bash
snakemake --cores 1
```

## Configuration

All configuration options, and their descriptions, are available in the `config/config_rsv_custom.yaml` file.

Metadata columns (`metadata_cols`) and sequence groupings (`group_cols`) specific to this pipeline are defined in the `config/config_rsv_custom.yaml` file.

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
