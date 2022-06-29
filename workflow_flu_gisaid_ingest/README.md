# GISAID ingestion for Flu PathMut

**NOTE: While this code is open-source, this workflow is intended for internal use only. It utilizes a GISAID data source that is not intended for use by the general public. Please see the `workflow_flu_genbank_ingest` to receive data from GenBank, or use either workflows as a template to ingest custom/in-house data.**

## Running

```
snakemake --cores 4
```

## Known Issues

If running into an error during the `copy_changed_files` step, try raising your open files limit with:

```
ulimit -n 100000
```
