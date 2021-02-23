#!/usr/bin/env python3
# coding: utf-8

import json
import os
import pandas as pd
import sqlite3
import sys
import yaml

from pathlib import Path
from yaml import load, dump, Loader, Dumper

sys.path.append(Path(__file__).parent)
from load_snvs import process_dna_snvs, process_aa_snvs
from gene_protein_defs import load_genes_or_proteins

DB_NAME = "covidcg.db"

if os.environ.get("CONFIGFILE", None) is None:
    print("NO CONFIG FILE FOUND")

with open(os.environ.get("CONFIGFILE"), "r") as fp:
    config = load(fp.read(), Loader=Loader)

data_path = Path(config["data_folder"])
static_data_path = Path(config["static_data_folder"])

genes = load_genes_or_proteins(str(static_data_path / "genes.json"))
proteins = load_genes_or_proteins(str(static_data_path / "proteins.json"))


def main():

    # Drop existing database if it exists
    db_path = Path("../") / DB_NAME
    if db_path.exists() and db_path.is_file():
        db_path.unlink()

    # Create new database
    conn = sqlite3.connect(str(db_path))
    c = conn.cursor()

    # Load metadata map
    with open(data_path / "metadata_map.json", "r") as fp:
        metadata_map = json.loads(fp.read())

    # Make a table for each metadata field
    for field in config["metadata_cols"].keys():
        table_name = "metadata_{}".format(field)
        c.execute(
            """
            CREATE TABLE "{table_name}" (
                id INTEGER PRIMARY KEY,
                value TEXT NOT NULL
            );
        """.format(
                table_name=table_name
            )
        )
        c.executemany(
            """
            INSERT INTO {table_name} VALUES (?, ?)
        """.format(
                table_name=table_name
            ),
            list(metadata_map[field].items()),
        )
        c.execute(
            """
            CREATE INDEX "ix_{table_name}_id" ON "{table_name}"("id");
            """.format(
                table_name=table_name
            )
        )

    # DNA SNVs
    dna_snp = process_dna_snvs(metadata_map["dna_snp"])
    dna_snp.to_sql("dna_snp", conn, index=True, index_label="id")
    c.execute('CREATE INDEX "ix_dna_snp_pos" ON "dna_snp"("pos");')

    gene_aa_snp = process_aa_snvs(metadata_map["gene_aa_snp"], "gene", genes)
    gene_aa_snp.to_sql("gene_aa_snp", conn, index=True, index_label="id")
    c.execute('CREATE INDEX "ix_gene_aa_snp_pos" ON "gene_aa_snp"("pos");')
    c.execute('CREATE INDEX "ix_gene_aa_snp_nt_pos" ON "gene_aa_snp"("nt_pos");')
    c.execute('CREATE INDEX "ix_gene_aa_snp_gene" ON "gene_aa_snp"("gene");')

    protein_aa_snp = process_aa_snvs(
        metadata_map["protein_aa_snp"], "protein", proteins
    )
    protein_aa_snp.to_sql("protein_aa_snp", conn, index=True, index_label="id")
    c.execute('CREATE INDEX "ix_protein_aa_snp_pos" ON "protein_aa_snp"("pos");')
    c.execute('CREATE INDEX "ix_protein_aa_snp_nt_pos" ON "protein_aa_snp"("nt_pos");')
    c.execute(
        'CREATE INDEX "ix_protein_aa_snp_protein" ON "protein_aa_snp"("protein");'
    )

    # Locations
    location_map = pd.read_json(data_path / "location_map.json")
    location_map.to_sql("location", conn, index=True, index_label="id")

    # Consensus SNVs
    with (data_path / "group_consensus_snps.json").open("r") as fp:
        group_consensus_snps = json.loads(fp.read())

    snp_fields = ["dna", "gene_aa", "protein_aa"]
    for grouping in group_consensus_snps.keys():
        for snp_field in snp_fields:
            # Create tables
            c.execute(
                """
                CREATE TABLE "{grouping}_consensus_{snp_field}_snp" (
                    name TEXT NOT NULL,
                    snp_id INTEGER NOT NULL
                );
                """.format(
                    grouping=grouping, snp_field=snp_field
                )
            )

            # Collect tuples of (group, snp_id)
            group_snps = []
            for group in group_consensus_snps[grouping].keys():
                for snp_id in group_consensus_snps[grouping][group][
                    snp_field + "_snp_ids"
                ]:
                    group_snps.append((group, snp_id))

            # Inserts
            c.executemany(
                'INSERT INTO "{grouping}_consensus_{snp_field}_snp" VALUES (?, ?);'.format(
                    grouping=grouping, snp_field=snp_field
                ),
                group_snps,
            )

    # Sequence metadata
    case_data = pd.read_json(data_path / "case_data.json")
    # case_data = case_data.set_index("Accession ID")
    case_data["collection_date"] = pd.to_datetime(case_data["collection_date"])
    case_data["submission_date"] = pd.to_datetime(case_data["submission_date"])

    case_data.drop(
        columns=["dna_snp_str", "gene_aa_snp_str", "protein_aa_snp_str",]
    ).to_sql("sequence", conn, index=True, index_label="id")
    # Create indices
    c.execute(
        'CREATE INDEX "ix_sequence_collection_date" ON "sequence"("collection_date");'
    )
    c.execute(
        'CREATE INDEX "ix_sequence_submission_date" ON "sequence"("submission_date");'
    )
    c.execute('CREATE INDEX "ix_sequence_location_id" ON "sequence"("location_id");')
    for field in config["metadata_cols"].keys():
        c.execute(
            'CREATE INDEX "ix_sequence_{field}" ON "sequence"("{field}");'.format(
                field=field
            )
        )

    # Sequence SNV data
    snp_fields = ["dna", "gene_aa", "protein_aa"]
    for snp_field in snp_fields:
        snp_col = snp_field + "_snp_str"
        c.execute(
            """
            CREATE TABLE "sequence_{snp_field}_snp" (
                sequence_id INTEGER NOT NULL,
                snp_id INTEGER NOT NULL
            );
            """.format(
                snp_field=snp_field
            )
        )
        c.executemany(
            'INSERT INTO "sequence_{snp_field}_snp" VALUES (?, ?);'.format(
                snp_field=snp_field
            ),
            case_data[[snp_col]].explode(snp_col).to_records().tolist(),
        )
        c.execute(
            'CREATE INDEX "ix_sequence_{snp_field}_sequence_id" ON "sequence_{snp_field}_snp"("sequence_id");'.format(
                snp_field=snp_field
            )
        )
        c.execute(
            'CREATE INDEX "ix_sequence_{snp_field}_snp_id" ON "sequence_{snp_field}_snp"("snp_id");'.format(
                snp_field=snp_field
            )
        )

    conn.commit()
    conn.close()


if __name__ == "__main__":
    main()
