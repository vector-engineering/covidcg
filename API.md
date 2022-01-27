# COVID CG API

- [COVID CG API](#covid-cg-api)
- [Data enabling COVID CG](#data-enabling-covid-cg)
- [Metadata mappings](#metadata-mappings)
- [Aggregate data](#aggregate-data)
  - [Mutation mode](#mutation-mode)
  - [Group mode (lineage mode)](#group-mode-lineage-mode)
- [Group mutation frequencies](#group-mutation-frequencies)
- [Dynamic group mutation frequencies](#dynamic-group-mutation-frequencies)

## Data enabling COVID CG

We are extremely grateful to the [GISAID Initiative](https://www.gisaid.org/) and all its data contributors, i.e. the Authors from the Originating laboratories responsible for obtaining the speciments and the Submitting laboratories where genetic sequence data were generated and shared via the GISAID Initiative, on which this research is based.

Elbe, S., and Buckland-Merrett, G. (2017) Data, disease and diplomacy: GISAID’s innovative contribution to global health. _Global Challenges_, 1:33-46. DOI:[10.1002/gch2.1018](https://doi.org/10.1002/gch2.1018) PMCID: [31565258](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6607375/)

Note: When using results from these analyses in your manuscript, ensure that you acknowledge the contributors of data, i.e. _We gratefully acknowledge all the Authors from the Originating laboratories responsible for obtaining the speciments and the Submitting laboratories where genetic sequence data were generated and shared via the GISAID Initiative, on which this research is based_.

## Metadata mappings

Much of the input to the API and output from the API is encoded into integers for faster query times. For example, mutations from the aggregate data API are returned as integer IDs instead of human-readable mutations. To map these back from IDs, first get the metadata map with:

```
curl https://covidcg.org/init
```

This map contains integer mappings that the server is currently using, for mutations and all other metadata.

**Update your metadata map every time before running a real query. The metadata mappings change every day**

## Aggregate data

### Mutation mode

```
curl --header "Content-Type: application/json" --request POST --data '{
  "group_key": "mutation",
  "dna_or_aa": "AA",
  "coordinate_mode": "gene",
  "coordinate_ranges": [[21563, 25384]],
  "selected_gene": "S",
  "region": [0, 1, 2],
  "country": [10, 11, 12],
  "division": [3, 4, 5],
  "location": [19, 20, 21],
  "selected_group_fields": { "lineage": ["AY.4", "BA.1"] },
  "selected_metadata_fields": { "host": [0, 1] },
  "start_date": "2021-01-01",
  "end_date": "2021-06-01",
  "subm_start_date": "2021-03-01",
  "subm_end_date": "2021-05-01"
}' https://covidcg.org/data
```

#### Parameters

- `group_key`: string - required
  - "mutation" to aggregate over sequence mutations
  - "lineage" to aggregate over sequence lineage (see section below)
- `dna_or_aa`: string - required
  - "DNA": Mutations are on the nucleotide level
  - "AA": Mutations are on the amino acid level
- `coordinate_mode`: string - required only for AA mode
  - "gene": Amino acid mutations are relative to canonical genes (may be missing ORFs within some genes). i.e., will return mutations relative to Orf1a coding frame instead of relative to nsp1's or nsp2's
  - "protein": Amino acid mutations are relative to all ORFs
- `coordinate_ranges`: array of array of ints - required
  - Only extracts mutations within this set of ranges. Each range is an array of two integers, `[a, b]`, where a and b are both nucleotide positions relative to the WIV04 reference sequence MN996528.1.
- `selected_gene`: string - required only for gene coordinate mode
  - Specify the name of the gene to extract mutations from. A list of gene names is in `static_data/genes.json`
- `selected_protein`: string - required only for protein coordinate mode
  - Specify the name of the protein to extract mutations from. A list of protein names is in `static_data/proteins.json`
- `region`, `country`, `division`, `location`: array of ints - at least one is required
  - Geographical IDs from which to filter sequences from. _These are integers representing IDs of geographical locations, not strings!_
  - To acquire mappings of geographical IDs -> names, see the previous [Metadata mappings](#metadata-mappings) section.
- `selected_group_fields`: object
  - Has the form: `{ group: [group names...], ... }`
  - e.g., to filter for specific PANGO lineages corresponding to Omicron: `{ "lineage": ["BA.1", "BA.2", "BA.3"] }`
- `selected_metadata_fields`: object
  - Has the form: `{ metadata_field: [metadata_value_IDs...], ... }`
  - Available metadata fields are defined in the relevant config YAML file in `config/`. These can also be obtained from the [metadata map](#metadata-mappings)
  - Metadata values are integers, _not strings_, and are defined by the [metadata map](#metadata-mappings).
- `start_date`, `end_date`: string, required
  - Start and end date of sample collection. Dates are strings in ISO format (YYYY-MM-DD)
- `subm_start_date`, `subm_end_date`: string
  - Start and end date of sample submission. Dates are strings in ISO format (YYYY-MM-DD)

#### Returns

Returns a JSON object, with the format:

```
[
  {
    "location": string
      - Name of the specified location. e.g., "North America"
    "collection_date": integer
      - Collection date, in javascript time (milliseconds since Unix epoch)
    "group_id": array of integers or null
      - Represents a group of co-occurring mutation IDs. These mutation IDs can be mapped with the metadata map
      - null value represents no mutations within the specified genomic region
    "counts": integer
      - Occurrences of this group of mutations, at the location, at the given collection date
  },
  ...
]
```

Important to note:

- Mutations are only reported within the genomic coordinates specified in the request parameters. For example, if requesting AA mutations within the spike gene, mutations in the N gene will not be reported.
- Mutations can be double-counted if request contains overlapping locations. For example, if requesting mutations in "North America" and "United States", the result will contain entries for both these locations separately.
- If you desire to count individual mutations, then unfold the list of mutations in `group_id` and then re-aggregate

### Group mode (lineage mode)

```
curl --header "Content-Type: application/json" --request POST --data '{
  "group_key": "lineage",
  "region": [0, 1, 2],
  "country": [10, 11, 12],
  "division": [3, 4, 5],
  "location": [19, 20, 21],
  "selected_group_fields": { "lineage": ["AY.4", "BA.1"] },
  "selected_metadata_fields": { "host": [0, 1] },
  "start_date": "2021-01-01",
  "end_date": "2021-06-01",
  "subm_start_date": "2021-03-01",
  "subm_end_date": "2021-05-01"
}' https://covidcg.org/data
```

Same as mutation mode, except `group_key` is set to `lineage` (PANGO lineage designation) or another grouping as defined by the `group_cols` setting in the config YAML file. The GISAID site has an additional `clade` phylogenetic grouping. Any fields relating to mutatins or genomic coordinates can be omitted in group mode.

## Group mutation frequencies

This function queries all mutations associated with a group (i.e., PANGO lineage) in all of our data. For more targeted queries, see [Dynamic group mutation frequencies](#dynamic-group-mutation-frequencies)

```
curl --header "Content-Type: application/json" --request POST --data '{
  "group": "lineage",
  "mutation_type": "gene_aa",
  "consensus_threshold": 0.9
}' https://covidcg.org/group_mutation_frequencies
```

#### Parameters

`group` determines the field with which sequences are grouped. Typically this is set to `lineage`.

Valid choices for `group` are:

- `lineage`: PANGO lineage designation
- `clade`: GISAID clade designation (GISAID site only, not available on Genbank site)

`mutation_type` determines the format of mutations to return.

Valid choices for `mutation_type` are:

- `dna`: Mutations on the nucleotide level.
- `gene_aa`: Mutations on the AA level, assigned to genes (not all ORFs, i.e., will return mutations relative to Orf1a coding frame instead of relative to nsp1's or nsp2's)
- `protein_aa`: Mutations on the AA level, assigned to proteins (all ORFs)

`consensus_threshold` determines the cutoff for reporting mutations. i.e., mutations that occur less than this frequency will be excluded from the results. Set this to `0.0` to include all mutations associated with the particular grouping.

#### Returns

Returns a JSON object, with the format:

```
[
  {
    "name": string
        - Name of the grouping, i.e., for `group` of `lineage`, this will be the lineage name,
    "count": integer
        - Number of occurrences of this mutation within the group
    "fraction": float
        - Fraction of occurrences of this mutation within the group
    "gene": string
        - Name of gene - only for `mutation_type` of `gene_aa`
    "protein": string
        - Name of protein/ORF - only for `mutation_type` of `protein_aa`
    "mutation_id": integer
        - Mutation ID in our database - these IDs change frequently
    "pos": integer
        - Position of the mutation. With `mutation_type` of `dna`, this is the position in nucleotides relative to the WIV04 reference genome MN996528.1
    "ref": string
        The reference base (for `dna` mode) or amino acid (for `gene_aa` or `protein_aa` mode). A `_` character designates a stop codon, and a `-` character represents a gap, i.e., an insertion when ref = `-`
    "alt": string
        - The alternate base (for `dna` mode) or amino acid (for `gene_aa` or `protein_aa` mode). A `_` character designates a stop codon, and a `-` character represents a gap, i.e., a deletion when alt = `-`
    "mutation_name": string
        - Human readable name of the mutation, in the form ref:pos:alt if in `dna` mode, or gene/protein:ref:pos:alt if in `gene_aa` or `protein_aa` mode.
  },
  ...,
]
```

## Dynamic group mutation frequencies

```
curl --header "Content-Type: application/json" --request POST --data '{
  "group_key": "lineage",
  "mutation_type": "gene_aa",
  "consensus_threshold": 0.1,
  "region": [0, 1, 2, 3, 4, 5],
  "start_date": "2021-12-01",
  "end_date": "2022-01-01",
  "selected_group_fields": { "lineage": ["BA.1"] }
}' https://covidcg.org/group_mutation_frequencies_dynamic
```

#### Parameters

The `group`, `mutation_type`, and `consensus_threshold` parameters are the same as described in the above [Group mutation frequencies](#group-mutation-frequencies) function.

Additional sequence filtering parameters can be provided. All options described in the [Aggregate data](#aggregate-data) are supported, except for `group_key`, `dna_or_aa`, `coordinate_mode`, `coordinate_ranges`, `selected_gene`, and `selected_protein`.

#### Returns

Return data is the same as the above [Group mutation frequencies](#group-mutation-frequencies) function.
