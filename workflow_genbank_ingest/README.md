# GenBank ingestion for COVID-CG

This workflow downloads data from GenBank, chunks data into files by submission date, and cleans the metadata as provided by GenBank. Additionally, lineages are assigned to each sequence with [pangolin](https://github.com/cov-lineages/pangolin).

## Running

Please include `--use-conda` in the snakemake call, to automatically install `pangolin` and its dependent packages/libraries. For example:

```
snakemake --cores 8 --config data_folder=../data_genbank --use-conda
```

The environment file for `pangolin` is placed in `envs/pangolin.yaml`

## Configuration

All configuration options, and their descriptions, are available in `config.yaml` in the project root

Metadata columns (`metadata_cols`) and sequence groupings (`group_cols`) specific to this pipeline are defined in the `config.yaml` file under the `genbank` key.

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

## Acknowledgements

The GenBank download code is derived from the ncov-ingest tool ([https://github.com/nextstrain/ncov-ingest](https://github.com/nextstrain/ncov-ingest)) from the [Nextstrain](https://nextstrain.org/) team. The license for this code can be found in the `LICENSE_NEXTSTRAIN` file.

The `pangolin` lineage assignment tool is hosted on GitHub: [https://github.com/cov-lineages/pangolin](https://github.com/cov-lineages/pangolin).