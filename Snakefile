import datetime
import json
import pandas as pd
import os
import numpy as np
import shutil

from collections import defaultdict
from pathlib import Path

from cg_scripts.clean_metadata import clean_metadata
from cg_scripts.fasta import read_fasta_file
from cg_scripts.get_aa_snps import get_aa_snps
from cg_scripts.get_dna_snps import get_dna_snps
from cg_scripts.process_ack import process_ack
from cg_scripts.process_artic_primers import process_artic_primers
from cg_scripts.process_lineages import get_consensus_snps
from cg_scripts.process_locations import process_location_metadata, build_select_tree
from cg_scripts.process_snps import process_snps
from cg_scripts.preprocess_sequences import preprocess_sequences
from cg_scripts.util import hash_accession_id

# For updating Accession ID hashmaps on GCS
# from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

data_folder = "data"
static_data_folder = "static_data"

# Get today's date in ISO format (YYYY-MM-DD)
today_str = datetime.date.today().isoformat()

rule all:
    input:
        # Download latest data feed
        os.path.join(data_folder, "status", "download_" + today_str + ".done"),
        os.path.join(data_folder, "status", "merge_sequences_" + today_str + ".done"),
        # Process SNVs
        os.path.join(data_folder, "dna_snp.csv"),
        os.path.join(data_folder, "gene_aa_snp.csv"),
        os.path.join(data_folder, "protein_aa_snp.csv"),
        # Process case data
        os.path.join(data_folder, "case_data.json"),
        # Generate reference-related data
        static_data_folder + "/reference.json", 
        static_data_folder + "/genes.json",
        static_data_folder + "/proteins.json",
        static_data_folder + "/primers.json",
        # Calculate global sequencing stats?
        country_seq_stats = data_folder + "/country_score.json",
        # Packaged data
        data_package = data_folder + "/data_package.json",
        data_package_gz = data_folder + "/data_package.json.gz",
        # Standalone map spec
        standalone_map_spec = os.path.join(data_folder, "map_combined_standalone.vg.json")

# Throw an error if the data feed credentials don't exist
if not Path("credentials/data_feed_credentials").exists():
    error_msg = "Credentials file in \"credentials/data_feed_credentials\" (contents of username:password) does not exist. GISAID data feed credentials are required to process the data for this app. If you are working with in-house datasets, please contact the authors at covidcg@broadinstitute.org"
    raise Exception(error_msg)

# Get the login:password string from the "data_feed_credentials" file in the "credentials/" folder
with open("credentials/data_feed_credentials", "r") as fp:
    cred_str = fp.read().strip()

with open("credentials/data_feed_url", "r") as fp:
    feed_url = fp.read().strip()


rule download_data_feed:
    """Download the data feed JSON object from the GISAID database, using our data feed credentials. The resulting file will need to be decompressed by `decompress_data_feed`
    """
    output:
        feed = temp(os.path.join(data_folder, "feed.json")),
        status = touch(os.path.join(data_folder, "status", "download_" + today_str + ".done"))
    params:
        cred_str = cred_str,
        feed_url = feed_url
    threads: workflow.cores
    shell:
        "curl -L https://{params.cred_str}@{params.feed_url} | unxz --threads={threads} -c - > {output.feed}"


checkpoint rewrite_data_feed:
    """Split up the data feed's individual JSON objects into metadata and fasta files. Chunk the fasta files so that every day we only reprocess the subset of fasta files that have changed. The smaller the chunk size, the more efficient the updates, but the more files on the filesystem.
    On a 48-core workstation with 128 GB RAM, aligning 200 sequences takes about 10 minutes, and this is more acceptable than having to align 1000 sequences, which takes ~1 hour. We end up with thousands of files, but the filesystem seems to be handling it well.
    """
    input:
        data_feed = os.path.join(data_folder, "feed.json")
    output:
        fasta = directory(os.path.join(data_folder, "fasta_temp")),
        metadata = os.path.join(data_folder, "metadata.csv")
    params:
        chunk_size = 10000
    run:

        output_path = Path(output.fasta)

        # Make the output directory, if it hasn't been made yet
        # Snakemake won't make the directory itself, since it's a special
        # directory output
        if not output_path.exists():
            output_path.mkdir(exist_ok=True)
        else:
            # Erase all files in the output directory
            for fasta_file in output_path.iterdir():
                if fasta_file.is_file():
                    fasta_file.unlink()

        # Keep track of which chunk we're on and how far we're along the
        # current chunk
        chunk_i = 0
        cur_i = 0

        # Get fields for each isolate
        fields = []
        with open(input.data_feed, "r") as fp_in:
            isolate = json.loads(fp_in.readline().strip())
            for i, key in enumerate(isolate.keys()):
                # Skip the special sequence column
                if key == "sequence":
                    continue
                fields.append(key)

        # Store metadata entries as a list of dictionaries, for now
        # We'll wrap it in a pandas DataFrame later for easier serialization
        metadata_df = []
        
        with open(input.data_feed, "r") as fp_in:


            # Open up the initial fasta file for the first chunk
            # fasta_out = open(output.fasta + "/" + str(chunk_i) + ".fa", "w")
            fasta_by_subm_date = defaultdict(list)

            # I'm pretty sure the output file is write-buffered,
            # so we can just spam the write() function
            line_counter = 0
            for line in fp_in:

                # Flush results if chunk is full, or if we hit the end of the feed
                if not line or cur_i == params.chunk_size:
                    for date, seqs in fasta_by_subm_date.items():
                        # Open the output fasta file for this date chunk
                        fasta_out_path = output_path / (date + ".fa")
                        with fasta_out_path.open("a") as fp_out:
                            for seq in seqs:
                                fp_out.write(">" + seq[0] + "\n" + seq[1] + "\n")
                    # Reset chunk counter
                    chunk_i += 1
                    cur_i = 0
                    # Reset sequence dictionary
                    fasta_by_subm_date = defaultdict(list)

                try:
                    isolate = json.loads(line.strip())
                except json.JSONDecodeError as err:
                    print("ERROR PARSING LINE", line_counter)
                    print(line)
                    continue

                # Add to metadata list
                metadata_df.append({k:isolate[k] for k in fields})

                # Store sequence in dictionary
                fasta_by_subm_date[isolate['covv_subm_date']].append((
                    isolate["covv_accession_id"], isolate["sequence"]
                ))

                # Iterate the intra-chunk counter
                cur_i += 1

                line_counter += 1

        # Cast the list of dictionaries (list of metadata entries) into a pandas
        # DataFrame, and then serialize it to disk
        # Do this step since pandas can handle some special serialization options
        # that I didn't want to implement manually (such as wrapping certain strings 
        # in double quotes)
        metadata_df = pd.DataFrame(metadata_df, columns=fields)
        metadata_df.to_csv(output.metadata, index=False)


def get_changed_chunks(wildcards):
    """Helper function for detecting which chunks have changed in terms of their contents 
    (measured in equality by bytes of disk space occupied). Only re-process and re-align chunks which have changed. This will save us a ton of computational time, as now that there are 200K+
    isolates on GISAID, aligning them would take 1 week for the whole batch.
    """
    
    # Get all chunks from the fasta_temp directory
    checkpoint_output = checkpoints.rewrite_data_feed.get(**wildcards).output[0]
    chunks, = glob_wildcards(os.path.join(checkpoint_output, "{i}.fa"))

    # Keep track of which chunks have changed
    changed_chunks = []

    for chunk in chunks:
        fasta_temp_path = Path(data_folder) / "fasta_temp" / (chunk + ".fa")
        fasta_raw_path = Path(data_folder) / "fasta_raw" / (chunk + ".fa")

        # The chunk has changed if:
        # 1) The current chunk does not exist yet
        # 2) The new chunk and the current chunk are different sizes (in bytes)
        if (
                not fasta_raw_path.exists() or 
                not fasta_raw_path.is_file() or 
                fasta_temp_path.stat().st_size != fasta_raw_path.stat().st_size
            ):
            
            changed_chunks.append(chunk)

    # Return a list of fasta_temp files that have changed, so that they can be copied
    # over to fasta_raw by the below `copy_changed_files` rule
    return expand(os.path.join(data_folder, "fasta_temp", "{i}.fa"), i=changed_chunks)


checkpoint copy_changed_files:
    """Using the `get_changed_chunks` function, only copy fasta files which have changed
    from the purgatory `fasta_temp` folder to the `fasta_raw` folder. By copying over the files,
    it will flag to snakemake that they (and only they - not the others) will need to be
    reprocessed and realigned.
    """
    input:
        get_changed_chunks
    output:
        # Instead of explicitly defining the fasta_raw outputs
        # (and risking touching fasta files that haven't actually changed)
        # Have the output be a flag instead, that the "all" rule checks for
        # to make sure that we actually run this rule
        touch(os.path.join(data_folder, "status", "merge_sequences_" + today_str + ".done"))
    run:
        # Make the fasta_raw folder if it doesn't already exist yet
        # snakemake won't do this automatically since no fasta_raw files
        # are explicitly defined in any output
        (Path(data_folder) / "fasta_raw").mkdir(exist_ok=True)

        # For each changed chunk (as defined by `get_changed_chunks`),
        # copy over the fasta file from the `fasta_temp` folder to the `fasta_raw` folder
        for chunk in input:
            chunk = Path(chunk).stem
            fasta_temp_path = Path(data_folder) / "fasta_temp" / (chunk + ".fa")
            fasta_raw_path = Path(data_folder) / "fasta_raw" / (chunk + ".fa")

            shutil.copyfile(fasta_temp_path, fasta_raw_path)
            shutil.copystat(fasta_temp_path, fasta_raw_path)


rule preprocess_sequences:
    """Filter out sequences (adapted from van Dorp et al, 2020)
    1. Filter against nextstrain exclusion list
    2. Remove animal/environmental isolates (bat, pangolin, mink, tiger, cat, canine, env)
    3. Can't be less than 29700NT
	4. Can't have more than 5% ambiguous NT
    """
    input:
        fasta = os.path.join(data_folder, "fasta_raw", "{chunk}.fa"),
        nextstrain_exclude = os.path.join(static_data_folder, "nextstrain_exclude_20200520.txt")
    output:
        fasta = os.path.join(data_folder, "fasta_processed", "{chunk}.fa")
    run:
        preprocess_sequences(input.fasta, input.nextstrain_exclude, output.fasta)


rule bt2build:
    """Build the bowtie2 index for the reference sequence
    This should only be run once - it shouldn't ever change
    """
    input: os.path.join(static_data_folder, "reference.fasta")
    params:
        basename = os.path.join(data_folder, "reference_index", "reference")
    output:
        output1 = os.path.join(data_folder, "reference_index", "reference.1.bt2"),
        output2 = os.path.join(data_folder, "reference_index", "reference.2.bt2"),
        output3 = os.path.join(data_folder, "reference_index", "reference.3.bt2"),
        output4 = os.path.join(data_folder, "reference_index", "reference.4.bt2"),
        outputrev1 = os.path.join(data_folder, "reference_index", "reference.rev.1.bt2"),
        outputrev2 = os.path.join(data_folder, "reference_index", "reference.rev.2.bt2")
    shell:
        """
        bowtie2-build {input} {params.basename}
        """
        
rule align_sequences:
    """Align each sequence to the reference using bowtie2
    The alignment is costly but this allows us to get SNVs for each sequence
    by using a widely used tool, without having to write any bespoke code ourselves
    """
    input:
        fasta = os.path.join(data_folder, "fasta_processed", "{sample}.fa"),
        bt2_1 = os.path.join(data_folder, "reference_index", "reference.1.bt2"),
        bt2_2 = os.path.join(data_folder, "reference_index", "reference.2.bt2"),
        bt2_3 = os.path.join(data_folder, "reference_index", "reference.3.bt2"),
        bt2_4 = os.path.join(data_folder, "reference_index", "reference.4.bt2"),
        bt2_rev1 = os.path.join(data_folder, "reference_index", "reference.rev.1.bt2"),
        bt2_rev2 = os.path.join(data_folder, "reference_index", "reference.rev.2.bt2")
    params:
        index_name = os.path.join(data_folder, "reference_index", "reference")
    threads: workflow.cores
    output:
        sam = os.path.join(data_folder, "sam", "{sample}.sam")
    # bowtie2 is really memory intensive (10GB per thread), so make sure it 
    # doesn't crash by allocating a set number of cores, where ncores = RAM / 10GB
    shell:
        """
        bowtie2 --end-to-end --very-fast --xeq --reorder --sam-no-qname-trunc -x {params.index_name} -f -U {input.fasta} -S {output.sam} --threads {threads}
        """

rule get_dna_snps:
    """Find SNVs on the NT level for each sequence
    """
    input:
        reference = os.path.join(static_data_folder, "reference.fasta"),
        sam = os.path.join(data_folder, "sam", "{sample}.sam")
    output:
        dna_snp = os.path.join(data_folder, "dna_snp", "{sample}_dna_snp.csv")
    run:
        dna_snp_df = get_dna_snps(input.sam, input.reference)
        dna_snp_df.to_csv(output.dna_snp, index=False)


rule get_aa_snps:
    """Using the NT SNVs, translate genes/proteins and find SNVs on the AA level
    using both gene and protein definitions
    """
    input:
        dna_snp = os.path.join(data_folder, "dna_snp", "{sample}_dna_snp.csv"),
        reference = os.path.join(static_data_folder, "reference.fasta"),
        genes_file = os.path.join(static_data_folder, "genes.json"),
        proteins_file = os.path.join(static_data_folder, "proteins.json")
    output:
        gene_aa_snp = os.path.join(data_folder, "gene_aa_snp", "{sample}_gene_aa_snp.csv"),
        protein_aa_snp = os.path.join(data_folder, "protein_aa_snp", "{sample}_protein_aa_snp.csv")
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


# These functions will run after the `copy_changed_files` rule/checkpoint has completed
# This will return a list of all chunks (not just those that have changed), so that the
# `combine_snps` rule below can aggregate all SNV data. This is required since we don't
# know how many chunks to expect after each download, so there's no way to define the list
# of chunks during the initial compilation of the DAG. We'll only know the final structure
# of the DAG after running the `copy_changed_files` rule.

def get_all_dna_snp_chunks(wildcards):
    # Only do this to trigger the DAG recalculation
    checkpoint_output = checkpoints.copy_changed_files.get().output[0]
    return expand(
        os.path.join(data_folder, "dna_snp", "{chunk}_dna_snp.csv"),
        chunk=glob_wildcards(os.path.join(data_folder, "fasta_raw", "{i}.fa")).i
    )

def get_all_gene_aa_snp_chunks(wildcards):
    # Only do this to trigger the DAG recalculation
    checkpoint_output = checkpoints.copy_changed_files.get().output[0]
    return expand(
        os.path.join(data_folder, "gene_aa_snp", "{chunk}_gene_aa_snp.csv"),
        chunk=glob_wildcards(os.path.join(data_folder, "fasta_raw", "{i}.fa")).i
    )

def get_all_protein_aa_snp_chunks(wildcards):
    # Only do this to trigger the DAG recalculation
    checkpoint_output = checkpoints.copy_changed_files.get().output[0]
    return expand(
        os.path.join(data_folder, "protein_aa_snp", "{chunk}_protein_aa_snp.csv"),
        chunk=glob_wildcards(os.path.join(data_folder, "fasta_raw", "{i}.fa")).i
    )

rule combine_snps:
    """Merge the NT and AA SNVs from each chunk back into one file
    """
    input:
        dna_snp = get_all_dna_snp_chunks,
        gene_aa_snp = get_all_gene_aa_snp_chunks,
        protein_aa_snp = get_all_protein_aa_snp_chunks
    output:
        dna_snp = os.path.join(data_folder, "dna_snp.csv"),
        gene_aa_snp = os.path.join(data_folder, "gene_aa_snp.csv"),
        protein_aa_snp = os.path.join(data_folder, "protein_aa_snp.csv")
    run:
        # This was originally done with this shell command:
        """
        # https://apple.stackexchange.com/questions/80611/merging-multiple-csv-files-without-merging-the-header
        awk '(NR == 1) || (FNR > 1)' {input.dna_snp} > {output.dna_snp}
        awk '(NR == 1) || (FNR > 1)' {input.gene_aa_snp} > {output.gene_aa_snp}
        awk '(NR == 1) || (FNR > 1)' {input.protein_aa_snp} > {output.protein_aa_snp}
        """
        # But this did not handle empty files well. So instead just writing this part in
        # python. It'll be a bit slower but whatever

        snp_types = ["dna_snp", "gene_aa_snp", "protein_aa_snp"]

        for snp_type in snp_types:
            with open(output[snp_type], "w") as fp_out:
                chunks = input[snp_type]
                for i, chunk in enumerate(chunks):
                    with open(chunk, "r") as fp_in:
                        for j, line in enumerate(fp_in):
                            # Write the header of the first file
                            # Or write any line that's not the header
                            # (to avoid writing the header more than once)
                            if (i == 0 and j == 0) or j > 0:
                                fp_out.write(line)

        
rule process_snps:
    """Filter out low-occurrence SNVs and assign each
    SNV an integer ID, to be mapped back on the front-end
    This saves a lot of data transfer space within the data package
    """
    input:
        dna_snp = os.path.join(data_folder, "dna_snp.csv"),
        gene_aa_snp = os.path.join(data_folder, "gene_aa_snp.csv"),
        protein_aa_snp = os.path.join(data_folder, "protein_aa_snp.csv")
    params:
        count_threshold = 3
    output:
        dna_snp_group = os.path.join(data_folder, "dna_snp_group.csv"),
        gene_aa_snp_group = os.path.join(data_folder, "gene_aa_snp_group.csv"),
        protein_aa_snp_group = os.path.join(data_folder, "protein_aa_snp_group.csv"),

        dna_snp_map = os.path.join(data_folder, "dna_snp_map.json"),
        gene_aa_snp_map = os.path.join(data_folder, "gene_aa_snp_map.json"),
        protein_aa_snp_map = os.path.join(data_folder, "protein_aa_snp_map.json")
    run:
        snp_types = ["dna_snp", "gene_aa_snp", "protein_aa_snp"]
        snp_modes = ["dna", "gene_aa", "protein_aa"]
        for snp_type, snp_mode in zip(snp_types, snp_modes):
            snp_group_df, snp_map = process_snps(
                input[snp_type], 
                mode=snp_mode, 
                count_threshold=params.count_threshold
            )

            # Save files
            snp_group_df.to_csv(output[snp_type + "_group"], index=False)
            snp_map.to_json(output[snp_type + "_map"], orient="index")


rule generate_ui_data:
    """Main rule for generating the data files for the browser
    Mostly just a bunch of joins
    """
    input:
        metadata = os.path.join(data_folder, "metadata.csv"),
        dna_snp_group = os.path.join(data_folder, "dna_snp_group.csv"),
        gene_aa_snp_group = os.path.join(data_folder, "gene_aa_snp_group.csv"),
        protein_aa_snp_group = os.path.join(data_folder, "protein_aa_snp_group.csv"),
        location_corrections = os.path.join(static_data_folder, "location_corrections.csv"),
        emoji_map_file = os.path.join(static_data_folder, "country_to_emoji.xls")
    output:
        accession_hashmap = os.path.join(data_folder, "accession_hashmap.csv"),
        metadata_map = os.path.join(data_folder, "metadata_map.json"),
        location_map = os.path.join(data_folder, "location_map.json"),
        ack_map = os.path.join(data_folder, "ack_map.json"),
        geo_select_tree = os.path.join(data_folder, "geo_select_tree.json"),
        case_data = os.path.join(data_folder, "case_data.json"),
        # CSV is just for excel/debugging
        case_data_csv = os.path.join(data_folder, "case_data.csv")
    run:
        # Load metadata
        df = pd.read_csv(input.metadata)
        df = (
            df
            .rename(columns={'covv_accession_id': 'Accession ID'})
            .set_index("Accession ID")
        )

        # Clean metadata
        df = clean_metadata(df)

        # Filter out "None" lineages
        df = df.loc[df["lineage"] != "None", :]
        df = df.loc[df["lineage"] != "nan", :]
        # Exclude sequences without a lineage/clade assignment
        df = df.loc[~pd.isnull(df["lineage"]), :]
        df = df.loc[~pd.isnull(df["clade"]), :]

        dna_snp_group_df = pd.read_csv(input.dna_snp_group, index_col="Accession ID")
        gene_aa_snp_group_df = pd.read_csv(input.gene_aa_snp_group, index_col="Accession ID")
        protein_aa_snp_group_df = pd.read_csv(input.protein_aa_snp_group, index_col="Accession ID")

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
        location_df, location_map_df = process_location_metadata(df, input.location_corrections)

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
        ).drop(columns=["covv_location"])

        # Process acknowledgement metadata
        df, ack_map = process_ack(df)
        # Save acknowledgements map
        ack_map.to_json(output.ack_map, orient="records")

        # Hash Accession IDs. Only take the first 8 chars, that's good enough
        df["hashed_id"] = np.random.rand(len(df))
        df["hashed_id"] = df["hashed_id"].astype(str).apply(hash_accession_id).str.slice(stop=8)
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
    """Write some of the reference sequence data as JSON
    files that can be easily loaded by the front-end
    """
    input:
        reference = static_data_folder + "/reference.fasta",
        primers = static_data_folder + "/primers.csv"
    output:
        # Write data to JSON for the JS/UI to handle
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

        # Load primers, write to JSON
        primers_df = pd.read_csv(input.primers, comment="#")
        # Only take a subset of the data to kee file sizes down
        primers_df[["Institution", "Name", "Sequence", "Reverse", "Start", "End"]].to_json(
            output.primers, orient="records"
        )


rule get_consensus_snps:
    """For each lineage and clade, get the lineage/clade-defining SNVs,
    on both the NT and AA level
    Lineage/clade-defining SNVs are defined as SNVs which occur in
    >= [consensus_fraction] of sequences within that lineage/clade.
    [consensus_fraction] is a parameter which can be adjusted here
    """
    input:
        case_data = os.path.join(data_folder, "case_data.csv")
    params:
        consensus_fraction = 0.9
    output:
        lineage_snp = os.path.join(data_folder, "lineage_snp.json"),
        clade_snp = os.path.join(data_folder, "clade_snp.json")
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
    """Get the number of sequences in each group
    Doing this in the pipeline just saves some work for the browser later
    """
    input:
        case_data = os.path.join(data_folder, "case_data.csv")
    output:
        global_group_counts = os.path.join(data_folder, "global_group_counts.json")
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

        with open(output.global_group_counts, "w") as fp:
            fp.write(json.dumps(global_group_counts))


rule process_artic_primers:
    """Write ARTIC primer data from GitHub into our primer.csv format
    This should only be run once - and not part of the pipeline generally
    """
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


rule calc_global_sequencing_efforts:
    """Merge sequence data with case counts and ISO geographical data,
    to produce the "Global Sequencing Effort" plot in the web app
    """
    input:
        case_data = os.path.join(data_folder, "case_data.csv"),
        location_map = os.path.join(data_folder, "location_map.json")
    output:
        country_seq_stats = data_folder + "/country_score.json"
    run:
        # Load case counts by country
        case_count_df = pd.read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")

        # Upgrade some province/states to country/regions
        upgrade_provinces = [
            "Hong Kong", "Macau", 
            "Faroe Islands", "Greenland", 
            "French Guiana", "French Polynesia", "Guadeloupe", "Martinique",
            "Mayotte", "New Caledonia", "Reunion", "Saint Barthelemy",
            "Saint Pierre and Miquelon", "St Martin", "Aruba",
            "Bonaire, Sint Eustatius and Saba", "Curacao", "Sint Maarten",
            "Anguilla", "Bermuda", "British Virgin Islands", "Cayman Islands",
            "Falkland Islands (Malvinas)",
            "Gibraltar", "Isle of Man", "Channel Islands",
            "Montserrat", "Turks and Caicos Islands",
            "American Samoa",
            "Guam", "Northern Mariana Islands", "Virgin Islands",
            "Puerto Rico"
        ]
        upgrade_province_inds = case_count_df["Province/State"].isin(upgrade_provinces)
        case_count_df.loc[upgrade_province_inds, "Country/Region"] = (
            case_count_df.loc[upgrade_province_inds, "Province/State"]
        )

        # Group by country/region
        case_count_df = (
            case_count_df
            .drop(columns=["Lat", "Long"])
            .groupby("Country/Region")
            .agg(np.sum)
            .reset_index()
        )
        # Unpivot table
        case_count_df = pd.melt(
            case_count_df, 
            id_vars=["Country/Region"], 
            var_name="date", 
            value_name="cumulative_cases"
        )
        # Convert date strings to datetime objects
        case_count_df["date"] = pd.to_datetime(case_count_df["date"])
        case_count_df["month"] = case_count_df["date"].dt.to_period("M")

        case_df = pd.read_csv(input.case_data, index_col="Accession ID")
        location_map = pd.read_json(input.location_map)
        location_map = location_map.set_index("index")
        case_df = case_df.join(location_map, on="location_id", how="left")

        case_df["collection_date"] = pd.to_datetime(case_df["collection_date"], errors="coerce")
        case_df["submission_date"] = pd.to_datetime(
            case_df["submission_date"], errors="coerce"
        )

        # Remove failed date parsing
        case_df = case_df.loc[
            (~pd.isnull(case_df["collection_date"])) & 
            (~pd.isnull(case_df["submission_date"]))
        ]

        # Only take dates from 2019-12-15
        case_df = case_df.loc[
            case_df["collection_date"] > pd.to_datetime("2019-12-15")
        ]

        # Calculate time deltas
        case_df["turnaround_days"] = (
            case_df["submission_date"] - case_df["collection_date"]
        ).dt.days
        # Extract month
        case_df["year_month"] = case_df["collection_date"].dt.to_period("M")

        # Remove invalid submission dates (negative turnaround times)
        case_df = case_df.loc[case_df["turnaround_days"] >= 0]

        # Upgrade provinces to countries
        upgrade_inds = case_df["division"].isin(upgrade_provinces)
        case_df.loc[upgrade_inds, "country"] = case_df.loc[upgrade_inds, "division"]

        # Load UID ISO FIPS lookup table
        iso_lookup_df = pd.read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/UID_ISO_FIPS_LookUp_Table.csv")
        # Upgrade provinces to country/regions
        upgrade_inds = iso_lookup_df["Province_State"].isin(upgrade_provinces)
        iso_lookup_df.loc[upgrade_inds, "Country_Region"] = iso_lookup_df.loc[upgrade_inds, "Province_State"]

        # Only take countries, then set as the index
        iso_lookup_df = (
            iso_lookup_df
            .loc[
                (upgrade_inds & pd.isnull(iso_lookup_df["Admin2"])) | 
                (pd.isnull(iso_lookup_df["Province_State"]))
            ]
            .set_index("Country_Region")
            .rename({
                "US": "USA",
                "Congo (Kinshasa)": "Democratic Republic of the Congo",
                "Congo (Brazzaville)": "Republic of the Congo",
                "Korea, South": "South Korea",
                "Taiwan*": "Taiwan",
                "Czechia": "Czech Republic",
                "Burma": "Myanmar"
            })
        )

        # Combine everything together
        country_df = (
            case_df
            # .loc[
            #     (nextmeta_df["date"] > pd.to_datetime("2020-01-01")) &
            #     (nextmeta_df["date"] < pd.to_datetime("2020-07-01"))
            # ]
            .reset_index()
            .groupby("country").agg(
                median_turnaround_days=pd.NamedAgg(column="turnaround_days", aggfunc=np.median),
                min_turnaround_days=pd.NamedAgg(column="turnaround_days", aggfunc=np.min),
                max_turnaround_days=pd.NamedAgg(column="turnaround_days", aggfunc=np.max),
                num_sequences=pd.NamedAgg(column="Accession ID", aggfunc="count")
            )
            .join(
                case_count_df
                .groupby("Country/Region")
                ["cumulative_cases"]
                .agg(np.max)
                .rename({
                    "US": "USA",
                    "Congo (Kinshasa)": "Democratic Republic of the Congo",
                    "Congo (Brazzaville)": "Republic of the Congo",
                    "Korea, South": "South Korea",
                    "Taiwan*": "Taiwan",
                    "Czechia": "Czech Republic",
                    "Burma": "Myanmar"
                })
            ).join(
                iso_lookup_df, 
                how="right"
            )
            .reset_index()
            .rename(columns={
                "index": "country",
                "cumulative_cases": "cases"
            })
        )

        # Fill some column"s missing values with 0
        country_df["num_sequences"] = country_df["num_sequences"].fillna(0)
        country_df["sequences_per_case"] = (
            country_df["num_sequences"] / country_df["cases"]
        ).fillna(0)

        # Only take some columns
        country_df = country_df.loc[:, [
            "UID", "Country_Region",
            "median_turnaround_days","min_turnaround_days","max_turnaround_days",
            "num_sequences", "cases", "sequences_per_case"
        ]]

        # Write to disk
        # First write JSON to string
        country_df_str = country_df.to_json(orient='records')
        # Manually add some missing records
        country_df_str = country_df_str[:-1] + (
            ',{"UID":260,"Country_Region":"Fr. S. Antarctic Lands","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}' + 
            ',{"UID":795,"Country_Region":"Turkmenistan","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}' + 
            ',{"UID":10,"Country_Region":"Antarctica","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}' + 
            ',{"UID":408,"Country_Region":"North Korea","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}' + 
            ',{"UID":90,"Country_Region":"Solomon Islands","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}' + 
            ',{"UID":548,"Country_Region":"Vanuatu","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}' + 
            # GISAID really wants French Guiana separate from france, so in my custom geojson I made French Guiana ID: -98
            ',{"UID":-98,"Country_Region":"French Guiana","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}' + 
            # Northern Cyprus
            ',{"UID":-99,"Country_Region":"Northern Cyprus","median_turnaround_days":null,"min_turnaround_days":null,"max_turnaround_days":null,"num_sequences":null,"cases":null,"sequences_per_case":null}' + 
            "]"
        )
        with open(output.country_seq_stats, "w") as fp:
            fp.write(country_df_str)


rule assemble_data_package:
    """Assemble the complete data package, that will be downloaded
    by the app upon initial load
    """
    input:
        case_data = os.path.join(data_folder, "case_data.json"),
        ack_map = os.path.join(data_folder, "ack_map.json"),
        clade_snp = os.path.join(data_folder, "clade_snp.json"),
        country_score = os.path.join(data_folder, "country_score.json"),
        dna_snp_map = os.path.join(data_folder, "dna_snp_map.json"),
        gene_aa_snp_map = os.path.join(data_folder, "gene_aa_snp_map.json"),
        geo_select_tree = os.path.join(data_folder, "geo_select_tree.json"),
        global_group_counts = os.path.join(data_folder, "global_group_counts.json"),
        lineage_snp = os.path.join(data_folder, "lineage_snp.json"),
        location_map = os.path.join(data_folder, "location_map.json"),
        metadata_map = os.path.join(data_folder, "metadata_map.json"),
        protein_aa_snp_map = os.path.join(data_folder, "protein_aa_snp_map.json")
    output:
        data_package = os.path.join(data_folder, "data_package.json")
    run:
        data_package = {
            "data_date": datetime.date.today().isoformat()
        }

        with open(input.case_data, "r") as fp:
            data_package["case_data"] = json.loads(fp.read())
        with open(input.ack_map, "r") as fp:
            data_package["ack_map"] = json.loads(fp.read())
        with open(input.clade_snp, "r") as fp:
            data_package["clade_snp"] = json.loads(fp.read())
        with open(input.country_score, "r") as fp:
            data_package["country_score"] = json.loads(fp.read())
        with open(input.dna_snp_map, "r") as fp:
            data_package["dna_snp_map"] = json.loads(fp.read())
        with open(input.gene_aa_snp_map, "r") as fp:
            data_package["gene_aa_snp_map"] = json.loads(fp.read())
        with open(input.geo_select_tree, "r") as fp:
            data_package["geo_select_tree"] = json.loads(fp.read())
        with open(input.global_group_counts, "r") as fp:
            data_package["global_group_counts"] = json.loads(fp.read())
        with open(input.lineage_snp, "r") as fp:
            data_package["lineage_snp"] = json.loads(fp.read())
        with open(input.location_map, "r") as fp:
            data_package["location_map"] = json.loads(fp.read())
        with open(input.metadata_map, "r") as fp:
            data_package["metadata_map"] = json.loads(fp.read())
        with open(input.protein_aa_snp_map, "r") as fp:
            data_package["protein_aa_snp_map"] = json.loads(fp.read())

        with open(output.data_package, "w") as fp:
            fp.write(json.dumps(data_package))

rule compress_data_package:
    """Compress the above data package for quicker transport
    """
    input:
        os.path.join(data_folder, "data_package.json")
    output:
        os.path.join(data_folder, "data_package.json.gz")
    shell:
        """
        gzip -9 -k {input} -c > {output}
        """


rule create_standalone_map_spec:
    input:
        data = os.path.join(data_folder, "country_score.json"),
        spec = os.path.join("src", "vega_specs", "map_combined.vg.json")
    output:
        standalone_spec = os.path.join(data_folder, "map_combined_standalone.vg.json")
    run:
        # Load vega spec JSON
        with open(input.spec, "r") as fp:
            spec = json.loads(fp.read())
        
        # Load data JSON
        with open(input.data, "r") as fp:
            data = json.loads(fp.read())
        
        # Inject data array into spec object
        spec["data"][0]["values"] = data

        # Save spec with injected data
        with open(output.standalone_spec, "w") as fp:
            fp.write(json.dumps(spec))


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
