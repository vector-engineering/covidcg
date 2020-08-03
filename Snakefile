import json
import pandas as pd

from cg_scripts.fasta import read_fasta_file
from cg_scripts.get_aa_snps import get_aa_snps
from cg_scripts.get_dna_snps import get_dna_snps
from cg_scripts.process_ack import process_ack
from cg_scripts.process_artic_primers import process_artic_primers
from cg_scripts.process_lineages import get_consensus_snps
from cg_scripts.process_locations import process_location_metadata, build_select_tree
from cg_scripts.process_patient_metadata import process_patient_metadata
from cg_scripts.process_seq_metadata import process_seq_metadata
from cg_scripts.process_snps import process_snps
from cg_scripts.preprocess_sequences import preprocess_sequences
from cg_scripts.util import hash_accession_id

# For updating Accession ID hashmaps on GCS
# from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

data_folder = "data"
static_data_folder = "static_data"

SAMPLES, = glob_wildcards(data_folder + "/fasta_raw/{sample}.fasta")

rule all:
    input:
        expand(
            data_folder + "/dna_snp/{sample}_dna_snp.csv", 
            sample=SAMPLES
        ),
        expand(
            data_folder + "/gene_aa_snp/{sample}_gene_aa_snp.csv", 
            sample=SAMPLES
        ),
        expand(
            data_folder + "/protein_aa_snp/{sample}_protein_aa_snp.csv", 
            sample=SAMPLES
        ),
        data_folder + "/case_data.json",
        # Generate reference-related data
        static_data_folder + "/reference.json", 
        static_data_folder + "/genes.json",
        static_data_folder + "/proteins.json",
        static_data_folder + "/primers.json",
        # Run consensus SNPs
        data_folder + "/lineage_snp.json", 
        data_folder + "/clade_snp.json",
        # Get global group counts
        data_folder + "/global_group_counts.json"

rule preprocess_sequences:
    input:
        fasta = data_folder + "/fasta_raw/{sample}.fasta",
        nextstrain_exclude = static_data_folder + "/nextstrain_exclude_20200520.txt"
    output:
        fasta = data_folder + "/fasta_processed/{sample}.fasta"
    run:
        preprocess_sequences(input.fasta, input.nextstrain_exclude, output.fasta)

rule bt2build:
    input: static_data_folder + "/reference.fasta"
    params:
        basename=data_folder + "/reference_index/reference"
    output:
        output1=data_folder + "/reference_index/reference.1.bt2",
        output2=data_folder + "/reference_index/reference.2.bt2",
        output3=data_folder + "/reference_index/reference.3.bt2",
        output4=data_folder + "/reference_index/reference.4.bt2",
        outputrev1=data_folder + "/reference_index/reference.rev.1.bt2",
        outputrev2=data_folder + "/reference_index/reference.rev.2.bt2"
    shell:
        """
        bowtie2-build {input} {params.basename}
        """
        
rule align_sequences:
    input:
        fasta = data_folder + "/fasta_processed/{sample}.fasta",
        bt2_1 = data_folder + "/reference_index/reference.1.bt2",
        bt2_2 = data_folder + "/reference_index/reference.2.bt2",
        bt2_3 = data_folder + "/reference_index/reference.3.bt2",
        bt2_4 = data_folder + "/reference_index/reference.4.bt2",
        bt2_rev1=data_folder + "/reference_index/reference.rev.1.bt2",
        bt2_rev2=data_folder + "/reference_index/reference.rev.2.bt2"
    params:
        index_name = data_folder + "/reference_index/reference"
    output:
        sam = data_folder + "/sam/{sample}.sam"
    # bowtie2 is really memory intensive, so make sure it doesn't crash by only letting
    # one instance run at a time, and only give it 75% of the cores
    threads: workflow.cores * 0.75 
    shell:
        """
        bowtie2 --end-to-end --very-fast --xeq --reorder --sam-no-qname-trunc -x {params.index_name} -f -U {input.fasta} -S {output.sam} --threads {threads}
        """

rule get_dna_snps:
    input:
        reference = static_data_folder + "/reference.fasta",
        sam = data_folder + "/sam/{sample}.sam"
    output:
        dna_snp = data_folder + "/dna_snp/{sample}_dna_snp.csv"
    run:
        dna_snp_df = get_dna_snps(input.sam, input.reference)
        dna_snp_df.to_csv(output.dna_snp, index=False)


rule get_aa_snps:
    input:
        dna_snp = data_folder + "/dna_snp/{sample}_dna_snp.csv",
        reference = static_data_folder + "/reference.fasta",
        genes_file = static_data_folder + "/genes.csv",
        proteins_file = static_data_folder + "/proteins.csv"
    output:
        gene_aa_snp = data_folder + "/gene_aa_snp/{sample}_gene_aa_snp.csv",
        protein_aa_snp = data_folder + "/protein_aa_snp/{sample}_protein_aa_snp.csv"
    run:
        gene_aa_snp_df = get_aa_snps(
            input.dna_snp, 
            input.genes_file, 
            input.reference, 
            mode="gene"
        )
        protein_aa_snp_df = get_aa_snps(
            input.dna_snp, 
            input.proteins_file, 
            input.reference, 
            mode="protein"
        )

        gene_aa_snp_df.to_csv(output.gene_aa_snp, index=False)
        protein_aa_snp_df.to_csv(output.protein_aa_snp, index=False)


rule combine_snps:
    input:
        dna_snp = expand(
            data_folder + "/dna_snp/{sample}_dna_snp.csv", 
            sample=SAMPLES
        ),
        gene_aa_snp = expand(
            data_folder + "/gene_aa_snp/{sample}_gene_aa_snp.csv", 
            sample=SAMPLES
        ),
        protein_aa_snp = expand(
            data_folder + "/protein_aa_snp/{sample}_protein_aa_snp.csv", 
            sample=SAMPLES
        )
    output:
        dna_snp = data_folder + "/dna_snp.csv",
        gene_aa_snp = data_folder + "/gene_aa_snp.csv",
        protein_aa_snp = data_folder + "/protein_aa_snp.csv"
    shell:
        """
        # https://apple.stackexchange.com/questions/80611/merging-multiple-csv-files-without-merging-the-header
        awk '(NR == 1) || (FNR > 1)' {input.dna_snp} > {output.dna_snp}
        awk '(NR == 1) || (FNR > 1)' {input.gene_aa_snp} > {output.gene_aa_snp}
        awk '(NR == 1) || (FNR > 1)' {input.protein_aa_snp} > {output.protein_aa_snp}
        """

rule process_snps:
    input:
        dna_snp = data_folder + "/dna_snp.csv",
        gene_aa_snp = data_folder + "/gene_aa_snp.csv",
        protein_aa_snp = data_folder + "/protein_aa_snp.csv"
    params:
        count_threshold = 3
    output:
        dna_snp_group = data_folder + "/dna_snp_group.csv",
        gene_aa_snp_group = data_folder + "/gene_aa_snp_group.csv",
        protein_aa_snp_group = data_folder + "/protein_aa_snp_group.csv",

        dna_snp_map = data_folder + "/dna_snp_map.json",
        gene_aa_snp_map = data_folder + "/gene_aa_snp_map.json",
        protein_aa_snp_map = data_folder + "/protein_aa_snp_map.json"
    run:
        dna_snp_group_df, dna_snp_map = process_snps(
            input.dna_snp, 
            mode="dna", 
            count_threshold=params.count_threshold
        )
        gene_aa_snp_group_df, gene_aa_snp_map = process_snps(
            input.gene_aa_snp, 
            mode="gene_aa", 
            count_threshold=params.count_threshold
        )
        protein_aa_snp_group_df, protein_aa_snp_map = process_snps(
            input.protein_aa_snp, 
            mode="protein_aa", 
            count_threshold=params.count_threshold
        )

        # Save files
        dna_snp_group_df.to_csv(output.dna_snp_group, index=False)
        gene_aa_snp_group_df.to_csv(output.gene_aa_snp_group, index=False)
        protein_aa_snp_group_df.to_csv(output.protein_aa_snp_group, index=False)

        # Save maps
        # snp_map.to_csv(data_dir / "snp_map.csv", index_label="snp", header=["id"])
        dna_snp_map.to_json(output.dna_snp_map, orient="index")
        gene_aa_snp_map.to_json(output.gene_aa_snp_map, orient="index")
        protein_aa_snp_map.to_json(output.protein_aa_snp_map, orient="index")

PATIENT_META_FILES, = glob_wildcards("data/patient_meta/{patient_meta_file}.tsv")

rule process_patient_metadata:
    input:
        patient_meta = expand(
            "data/patient_meta/{patient_meta_file}.tsv", 
            patient_meta_file=PATIENT_META_FILES
        )
    output:
        patient_meta = data_folder + "/patient_meta.csv"
    run:
        patient_meta_df = process_patient_metadata(input.patient_meta)
        patient_meta_df.to_csv(output.patient_meta)


SEQ_META_FILES, = glob_wildcards("data/seq_meta/{seq_meta_file}.tsv")

rule process_seq_metadata:
    input:
        seq_meta = expand(
            "data/seq_meta/{seq_meta_file}.tsv", 
            seq_meta_file=SEQ_META_FILES
        )
    output:
        seq_meta = data_folder + "/seq_meta.csv"
    run:
        seq_meta_df = process_seq_metadata(input.seq_meta)
        seq_meta_df.to_csv(output.seq_meta)


ACK_FILES, = glob_wildcards("data/acknowledgements/{ack_file}.xls")

rule process_acknowledgements:
    input:
        ack = expand(
            "data/acknowledgements/{ack_file}.xls", 
            ack_file=ACK_FILES
        )
    output:
        ack_meta = data_folder + "/ack_meta.csv",
        ack_map = data_folder + "/ack_map.json"
    run:
        ack_df, ack_map = process_ack(input.ack)
        ack_df.to_csv(output.ack_meta)
        ack_map.to_json(output.ack_map, orient="index")
        
# Main rule for generating the data files for the browser
# Mostly just a bunch of joins
rule generate_ui_data:
    input:
        patient_meta = data_folder + "/patient_meta.csv",
        seq_meta = data_folder + "/seq_meta.csv",
        ack_meta = data_folder + "/ack_meta.csv",
        dna_snp_group = data_folder + "/dna_snp_group.csv",
        gene_aa_snp_group = data_folder + "/gene_aa_snp_group.csv",
        protein_aa_snp_group = data_folder + "/protein_aa_snp_group.csv",
        emoji_map_file = static_data_folder + "/country_to_emoji.xls"
    output:
        accession_hashmap = data_folder + "/accession_hashmap.csv",
        metadata_map = data_folder + "/metadata_map.json",
        location_map = data_folder + "/location_map.json",
        geo_select_tree = data_folder + "/geo_select_tree.json",
        case_data = data_folder + "/case_data.json",
        # CSV is just for excel/debugging
        case_data_csv = data_folder + "/case_data.csv"
    run:
        patient_meta_df = pd.read_csv(input.patient_meta, index_col="Accession ID")
        seq_meta_df = pd.read_csv(input.seq_meta, index_col="Accession ID")
        ack_meta_df = pd.read_csv(input.ack_meta, index_col="Accession ID")

        dna_snp_group_df = pd.read_csv(input.dna_snp_group, index_col="Accession ID")
        gene_aa_snp_group_df = pd.read_csv(input.gene_aa_snp_group, index_col="Accession ID")
        protein_aa_snp_group_df = pd.read_csv(input.protein_aa_snp_group, index_col="Accession ID")

        # Join patient and sequencing metadata on Accession ID
        df = patient_meta_df.join(
            seq_meta_df, on="Accession ID", how="left", sort=True
        )
        # Filter out "None" lineages
        df = df.loc[df["lineage"] != "None", :]

        # Join acknowledgement IDs onto main metadata dataframe
        df = df.join(ack_meta_df, on="Accession ID", how="left", sort=False)
        # Replace missing acknowledgement IDs with -1, then cast to integer
        df["ack_id"] = df["ack_id"].fillna(-1).astype(int)

        # Join SNPs to main dataframe
        # inner join to exclude filtered out sequences
        df = df.join(
            dna_snp_group_df[["snp_str"]],
            on="Accession ID",
            how="inner",
            sort=False,
        ).rename(columns={"snp_str": "dna_snp_str"})
        df = df.join(
            gene_aa_snp_group_df[["snp_str"]],
            on="Accession ID",
            how="inner",
            sort=False,
        ).rename(columns={"snp_str": "gene_aa_snp_str"})
        df = df.join(
            protein_aa_snp_group_df[["snp_str"]],
            on="Accession ID",
            how="inner",
            sort=False,
        ).rename(columns={"snp_str": "protein_aa_snp_str"})

        # Semicolon-delimited string to array of SNPs
        df[["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]] = (
            df[["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]]
            .astype(str)
            .applymap(
                lambda x: [int(_x) for _x in x.split(";")]
            )
        )

        # Process location metadata
        location_df, location_map_df = process_location_metadata(df)

        location_map_df.drop(columns=["loc_str"]).to_json(
            output.location_map, orient="records"
        )
        # Save tree as json file
        geo_select_tree = build_select_tree(
            location_df, 
            location_map_df, 
            emoji_map_file=input.emoji_map_file
        )
        with open(output.geo_select_tree, "w") as fp:
            fp.write(json.dumps(geo_select_tree))

        # Join location IDs onto main metadata dataframe, then drop original Location column
        df = df.join(
            location_df[["location_id"]], on="Accession ID", how="inner", sort=False
        ).drop(columns=["Location"])

        # Hash Accession IDs
        df["hashed_id"] = df.index.to_series().apply(hash_accession_id)
        # Create map of hash -> Accession ID
        accession_hash_df = df[["hashed_id"]]
        accession_hash_df.to_csv(output.accession_hashmap, index_label="Accession ID")

        # Delete old accession ID column, reassign to hashed ID
        df = (
            df.reset_index()
            .drop(columns=["Accession ID"])
            .rename(columns={"hashed_id": "Accession ID"})
            .set_index("Accession ID")
        )

        # Factorize some more metadata columns
        map_cols = [
            "gender",
            "patient_status",
            "passage",
            "specimen",
            "sequencing_tech",
            "assembly_method",
            "comment_type",
        ]
        metadata_maps = {}

        for i, col in enumerate(map_cols):
            factor = pd.factorize(df[col])

            id_col = col + "_id"
            df[id_col] = factor[0]

            metadata_maps[col] = pd.Series(factor[1]).to_dict()

        # Drop the original metadata columns
        df = df.drop(columns=map_cols)

        # Write the metadata map to a JSON file
        with open(output.metadata_map, "w") as fp:
            fp.write(json.dumps(metadata_maps))

        # Write final dataframe
        df.to_csv(output.case_data_csv, index_label="Accession ID")
        df.reset_index().to_json(output.case_data, orient="records")


rule write_reference_files:
    input:
        genes = static_data_folder + "/genes.csv",
        proteins = static_data_folder + "/proteins.csv",
        reference = static_data_folder + "/reference.fasta",
        primers = static_data_folder + "/primers.csv"
    output:
        # Write data to JSON for the JS/UI to handle
        genes = static_data_folder + "/genes.json",
        proteins = static_data_folder + "/proteins.json",
        reference = static_data_folder + "/reference.json",
        primers = static_data_folder + "/primers.json"
    run:
        # Write the reference fasta file to json
        # Load the reference sequence
        with open(input.reference, "r") as fp:
            lines = fp.readlines()
            ref = read_fasta_file(lines)
            ref_seq = list(ref.values())[0]

        with open(output.reference, "w") as fp:
            fp.write(json.dumps({"ref_seq": ref_seq}))
        
        # Load genes, write to JSON
        genes_df = pd.read_csv(input.genes, comment="#")
        genes_df.to_json(output.genes, orient="records")

        # Load proteins, write to JSON
        proteins_df = pd.read_csv(input.proteins, comment="#")
        proteins_df.to_json(output.proteins, orient="records")

        # Load primers, write to JSON
        primers_df = pd.read_csv(input.primers, comment="#")
        # Only take a subset of the data to kee file sizes down
        primers_df[["Institution", "Name", "Sequence", "Reverse", "Start", "End"]].to_json(
            output.primers, orient="records"
        )


rule get_consensus_snps:
    input:
        case_data = data_folder + "/case_data.csv"
    params:
        consensus_fraction = 0.9
    output:
        lineage_snp = data_folder + "/lineage_snp.json",
        clade_snp = data_folder + "/clade_snp.json"
    run:
        case_df = pd.read_csv(input.case_data, index_col="Accession ID")

        # Serialized list back to list 
        cols = ["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]
        for col in cols:
            case_df[col] = (
                case_df[col]
                .str.strip("[]")
                .str.split(",")
                .apply(lambda x: [int(_x) for _x in x])
            )

        lineage_snp_df = get_consensus_snps(
            case_df, "lineage",
            consensus_fraction=params.consensus_fraction
        )
        lineage_snp_df.to_json(output.lineage_snp, orient="records")

        clade_snp_df = get_consensus_snps(
            case_df, "clade",
            consensus_fraction=params.consensus_fraction
        )
        clade_snp_df.to_json(output.clade_snp, orient="records")


rule get_global_group_counts:
    input:
        case_data = data_folder + "/case_data.csv"
    output:
        global_group_counts = data_folder + "/global_group_counts.json"
    run:
        from itertools import chain
        from collections import Counter

        case_df = pd.read_csv(input.case_data, index_col="Accession ID")

        # Serialized list back to list 
        cols = ["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str"]
        for col in cols:
            case_df[col] = (
                case_df[col]
                .str.strip("[]")
                .str.split(",")
                .apply(lambda x: [int(_x) for _x in x])
            )

        global_group_counts = {}
        # Count lineages and clades
        global_group_counts["lineage"] = case_df["lineage"].value_counts().to_dict()
        global_group_counts["clade"] = case_df["clade"].value_counts().to_dict()

        # Count global SNV frequencies
        # Collapse list of lists into one list, then count individual
        # occurrences, then cast to a regular dict
        global_group_counts["dna_snp"] = dict(Counter(
            list(chain.from_iterable(case_df["dna_snp_str"]))
        ))
        global_group_counts["gene_aa_snp"] = dict(Counter(
            list(chain.from_iterable(case_df["gene_aa_snp_str"]))
        ))
        global_group_counts["protein_aa_snp"] = dict(Counter(
            list(chain.from_iterable(case_df["protein_aa_snp_str"]))
        ))

        with open(output.global_group_counts, 'w') as fp:
            fp.write(json.dumps(global_group_counts))



rule process_artic_primers:
    input:
        reference_file = static_data_folder + "/reference.fasta"
    params:
        artic_files = [
            "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V1/nCoV-2019.tsv",
            "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V2/nCoV-2019.tsv",
            "https://raw.githubusercontent.com/artic-network/artic-ncov2019/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv",
        ],
        artic_prefixes = ["V1", "V2", "V3"],
        artic_dates = ["2020-01-24", "2020-03-13", "2020-03-20"]
    output:
        artic_primers = static_data_folder + "/artic_primers.csv"
    run:
        artic_df = process_artic_primers(
            params.artic_files, 
            params.artic_prefixes, 
            params.artic_dates,
            input.reference_file
        )
        artic_df.to_csv(output.artic_primers, index=False)


# # This is only for site maintainers
# if "upload_hashmap" in config and config["upload_hashmap"]:
#     GS = GSRemoteProvider()
#     rule update_gcs_hashmap:
#         input:
#             data_folder + "/accession_hashmap.csv"
#         output:
#             GS.remote("covid-cg/accession_hashmap.csv")
#         shell:
#             """
#             cp {input[0]} {output[0]}
#             """
