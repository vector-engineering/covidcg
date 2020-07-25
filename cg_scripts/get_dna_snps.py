import pysam
from pathlib import Path

from cg_scripts.read_extractor_lite import ReadExtractor


def get_dna_snps(sam_file, dna_snp_file):
    samfile = pysam.AlignmentFile(sam_file, "r")  # pylint: disable=no-member

    all_dna_snps = []
    for read in samfile.fetch(until_eof=True):
        # Skip if unmapped
        if read.is_unmapped:
            continue

        read_extractor = ReadExtractor(read)

        # print(read.query_name)
        dna_snps = read_extractor.process_all()
        all_dna_snps.extend(dna_snps)

    samfile.close()

    # Write to disk
    # Open output files for writing
    fp_dna = Path(dna_snp_file).open("w")
    # Write headers
    fp_dna.write("taxon,pos,ref,alt\n")

    # Write DNA SNPs/indels
    for snp in all_dna_snps:
        fp_dna.write(
            "{},{},{},{}\n".format(snp[0], str(snp[1]), str(snp[2]), str(snp[3]))
        )

    fp_dna.close()


def load_dna_snp_file(dna_snp_file):
    dna_snp_df = pd.read_csv(dna_snp_file)
    # Extract the GISAID ID
    dna_snp_df["gisaid_id"] = dna_snp_df["taxon"].str.split("|", expand=True)[1]

    # Fill NaN values
    dna_snp_df["ref"].fillna("", inplace=True)
    dna_snp_df["alt"].fillna("", inplace=True)

    # Drop duplicate entries
    # dna_snp_df.drop_duplicates(["taxon", "pos"], inplace=True)

    dna_snp_df = dna_snp_df.reset_index(drop=True)
    print("Loaded {} DNA SNPs".format(len(dna_snp_df)), flush=True)
    # dna_snp_df.to_csv('dna_snp.csv', index=False)

    return dna_snp_df
