# Genbank (NCBI Virus) ingestion for Flu PathMut

## Running

```
snakemake --cores 8
```

## Known Issues

If running into an error during the `copy_changed_files` step, try raising your open files limit with:

```
ulimit -n 100000
```

## Configuration

All configuration options, and their descriptions, are available in the `config/config_flu_genbank.yaml` file.

Metadata columns (`metadata_cols`) and sequence groupings (`group_cols`) specific to this pipeline are defined in the `config/config_flu_genbank.yaml` file.

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
