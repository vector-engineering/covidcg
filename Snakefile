from cg_scripts.get_aa_snps import get_aa_snps
from cg_scripts.get_dna_snps import get_dna_snps
from cg_scripts.preprocess_sequences import preprocess_sequences

SAMPLES, = glob_wildcards("example_data/fasta_raw/{sample}.fasta")

rule all:
    input:
        expand("example_data/gene_aa_snp/{sample}_gene_aa_snp.csv", sample=SAMPLES)

rule preprocess_sequences:
    input:
        fasta = "example_data/fasta_raw/{sample}.fasta"
    params:
        nextstrain_exclude = "static_data/nextstrain_exclude_20200520.txt"
    output:
        fasta = "example_data/fasta_processed/{sample}.fasta"
    run:
        preprocess_sequences(input.fasta, params.nextstrain_exclude, output.fasta)

rule bt2build:
    input: "static_data/reference.fasta"
    params:
        basename="example_data/reference_index/reference"
    output:
        output1="example_data/reference_index/reference.1.bt2",
        output2="example_data/reference_index/reference.2.bt2",
        output3="example_data/reference_index/reference.3.bt2",
        output4="example_data/reference_index/reference.4.bt2",
        outputrev1="example_data/reference_index/reference.rev.1.bt2",
        outputrev2="example_data/reference_index/reference.rev.2.bt2"
    shell:
        """
        bowtie2-build {input} {params.basename}
        """
        
rule align_sequences:
    input:
        fasta = "example_data/fasta_processed/{sample}.fasta",
        bt2_1 = "example_data/reference_index/reference.1.bt2",
        bt2_2 = "example_data/reference_index/reference.2.bt2",
        bt2_3 = "example_data/reference_index/reference.3.bt2",
        bt2_4 = "example_data/reference_index/reference.4.bt2",
        bt2_rev1="example_data/reference_index/reference.rev1.bt2",
        bt2_rev2="example_data/reference_index/reference.rev2.bt2"
    params:
        index_name = "example_data/reference_index/reference"
    output:
        sam = "example_data/sam/{sample}.sam"
    shell:
        """
        bowtie2 --end-to-end --very-fast --xeq --reorder --sam-no-qname-trunc -x {params.index_name} -f -U {input.fasta} -S {output.sam} --threads {workflow.cores}
        """

rule get_dna_snps:
    input:
        sam = "example_data/sam/{sample}.sam"
    output:
        dna_snp = "example_data/dna_snp/{sample}_dna_snp.csv"
    run:
        get_dna_snps(input.sam, output.dna_snp)

rule get_aa_snps:
    input:
        dna_snp = "example_data/dna_snp/{sample}_dna_snp.csv"
    params:
        genes_file = "static_data/genes.csv",
        proteins_file = "static_data/proteins.csv"
    output:
        gene_aa_snp = "example_data/gene_aa_snp/{sample}_gene_aa_snp.csv",
        protein_aa_snp = "example_data/protein_aa_snp/{sample}_protein_aa_snp.csv"
    run:
        get_aa_snps(
            input.dna_snp, 
            params.genes_file, 
            params.proteins_file,
            output.gene_aa_snp,
            output.protein_aa_snp
        )

rule combine_snps:
    input:
        dna_snp = expand(
            "example_data/dna_snp/{sample}_dna_snp.csv", 
            sample=SAMPLES
        ),
        gene_aa_snp = expand(
            "example_data/gene_aa_snp/{sample}_gene_aa_snp.csv", 
            sample=SAMPLES
        ),
        protein_aa_snp = expand(
            "example_data/protein_aa_snp/{sample}_protein_aa_snp.csv", 
            sample=SAMPLES
        )
    output:
        dna_snp = "example_data/dna_snp.csv",
        gene_aa_snp = "example_data/gene_aa_snp.csv",
        protein_aa_snp = "example_data/protein_aa_snp.csv"
    shell:
        """
        # https://apple.stackexchange.com/questions/80611/merging-multiple-csv-files-without-merging-the-header
        awk '(NR == 1) || (FNR > 1)' {input.dna_snp} > {output.dna_snp}
        awk '(NR == 1) || (FNR > 1)' {input.gene_aa_snp} > {output.gene_aa_snp}
        awk '(NR == 1) || (FNR > 1)' {input.protein_aa_snp} > {output.protein_aa_snp}
        """



# rule all:
#     input:
#         "concat.txt"

# rule concat:
#     input:
#         expand("test_files/{sample}.txt", sample=SAMPLES)
#     output:
#         "concat.txt"
#     shell:
#         "cat {input} > {output}"