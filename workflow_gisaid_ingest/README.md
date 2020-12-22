# GISAID ingestion for COVID-CG

**NOTE: While this code is open-source, this workflow is intended for internal use only. It utilizes a GISAID data source that is not intended for use by the general public. Please see the `workflow_genbank_ingest` to receive data from GenBank, or use either workflows as a template to ingest custom/in-house data.**

This workflow downloads data from a GISAID endpoint, and chunks data by submission date, in order to avoid unnecessarily repeating the computationally-heavy step of aligning genomes to the reference.

Metadata from GISAID is cleaned and processed into the common format required by `workflow_main`.

## Running

```
snakemake --cores 8 --config data_folder=../data_gisaid
```

Lineage/clade assignments are provided by GISAID, so therefore we don't run `pangolin` again to determine sequence lineages. Please see `workflow_genbank_ingest` for an example of how to integrate `pangolin` into the data ingestion

## Configuration

All configuration options, and their descriptions, are available in `config.yaml`

Metadata columns (`metadata_cols`) and sequence groupings (`group_cols`) specific to this pipeline are defined in the `config.yaml` file inside the `workflow_main` folder, under the "GISAID" heading.

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

We are extremely grateful to the [GISAID Initiative](https://www.gisaid.org/) and all its data contributors, i.e. the Authors from the Originating laboratories responsible for obtaining the speciments and the Submitting laboratories where genetic sequence data were generated and shared via the GISAID Initiative, on which this research is based.

Elbe, S., and Buckland-Merrett, G. (2017) Data, disease and diplomacy: GISAID’s innovative contribution to global health. _Global Challenges_, 1:33-46. DOI:[10.1002/gch2.1018](https://doi.org/10.1002/gch2.1018) PMCID: [31565258](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6607375/)