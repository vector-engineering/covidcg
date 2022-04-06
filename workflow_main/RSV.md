This README will briefly explain the differences between the RSV pipeline and
the CovidCG pipeline

# Differences From CovidCG

## preprocess_sequences

The RSV site only has one condition to filter sequences: No sequence can have
more than 5% of its total length comprised of ambiguous NTs

## assign/combine_genotypes

The assign_genotypes rule aligns each sequence to the references listed in
/static_data/rsv/genotype_references.json. These reference sequences were
provided by Astrazeneca. Each sequence that aligns to any genotype is assigned
that genotype. This means that in the current implementation, a sequence can
have multiple genotypes. So far, each sequence that was assigned multiple
genotypes have been assigned genotypes within the same subtype but some
additional analysis may be required so that each sequence is only assigned its
"best" genotype match. The assign_genotypes rule produces .bam files which are
passed into the combine_genotypes rule.

The combine_genotypes rule creates a DataFrame indexed by Accession ID, with
each sequences' assigned genotype and subtype as columns. Subtype is determined
by the genotype assigned in the assign_genotype step

## align_sequences

No difference from covidcg pipeline

## extract_dna_mutations

As there are multiple reference sequences to work with in RSV, the mutation
extraction steps are different between RSV and CovidCG. The reference sequence
passed into ReadExtractor is determined by the reference_name field of the read.
The ReadExtractor itself has also changed slightly to include the subtype of the
read in its produced DataFrame. The DataFrame is also filtered to only include
sequences that have been assigned an A or B subtype

## process_genes_and_proteins

The process_genes_and_proteins rule has been updated to account for having
multiple reference sequences. The processed gene and protein file is an object
with the reference sequence name as keys. The values of the object are the usual
JSON files. The unprocessed gene/protein files have the same form

## extract_aa_mutations

Like extract_dna_mutations, extract_aa_mutations has been updated to work with
multiple reference sequences. The extract_aa_mutations now loops through the
reference names, assigning gene_or_protein_df and the subtype within the loop
as appropriate. In line 84 of the extract_aa_mutations function, we filter the
dna_mutation_df by the subtype being considered

## combine_all_data

The genotypes/subtypes are joined to case_data

## NO CHANGES IN THE REST OF THE PIPELINE
