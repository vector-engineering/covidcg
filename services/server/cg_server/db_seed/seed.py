# coding: utf-8

"""Database seeding

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import datetime
import io
import json
import os
import pandas as pd

from pathlib import Path
from psycopg2 import sql
from psycopg2.extras import Json

from cg_server.color import get_categorical_colormap
from cg_server.config import config
from .load_snvs import process_dna_snvs, process_aa_snvs

# root/services/server/cg_server/db_seed/seed.py
project_root = Path(__file__).parent.parent.parent.parent.parent
data_path = Path(os.getenv("DATA_PATH", project_root / config["example_data_folder"]))
static_data_path = Path(
    os.getenv("STATIC_DATA_PATH", project_root / config["static_data_folder"])
)

genes = pd.read_json(str(static_data_path / "genes_processed.json"))
genes = genes.set_index("name")
proteins = pd.read_json(str(static_data_path / "proteins_processed.json"))
proteins = proteins.set_index("name")


def df_to_sql(cur, df, table, index_label=None):
    buffer = io.StringIO()
    if index_label:
        df.to_csv(buffer, index_label=index_label, header=False)
    else:
        df.to_csv(buffer, index=False, header=False)
    buffer.seek(0)
    # cur.copy_from(buffer, table, sep=",")
    cur.copy_expert(
        """
        COPY "{table}" FROM STDIN WITH (FORMAT CSV)
        """.format(
            table=table
        ),
        buffer,
    )


def seed_database(conn, schema="public"):
    with conn.cursor() as cur:

        cur.execute(sql.SQL("SET search_path TO {};").format(sql.Identifier(schema)))

        print("Writing metadata maps...", end="", flush=True)

        # Load metadata map
        with open(data_path / "metadata_map.json", "r") as fp:
            metadata_map = json.loads(fp.read())

        # Make a table for each metadata field
        for field in config["metadata_cols"].keys():
            table_name = "metadata_{}".format(field)
            cur.execute(
                'DROP TABLE IF EXISTS "{table_name}";'.format(table_name=table_name)
            )
            cur.execute(
                """
                CREATE TABLE "{table_name}" (
                    id INTEGER PRIMARY KEY,
                    value TEXT NOT NULL
                );
            """.format(
                    table_name=table_name
                )
            )
            metadata_df = pd.DataFrame.from_records(
                list(metadata_map[field].items()), columns=["id", "value"]
            ).set_index("id")
            df_to_sql(cur, metadata_df, table_name, index_label="id")
            cur.execute(
                """
                CREATE INDEX "ix_{table_name}_id" ON "{table_name}"("id");
                """.format(
                    table_name=table_name
                )
            )

        print("done")

        print("Writing SNV maps...", end="", flush=True)

        # DNA SNVs
        dna_snp = process_dna_snvs(metadata_map["dna_snp"])
        cur.execute('DROP TABLE IF EXISTS "dna_snp";')
        cur.execute(
            """
            CREATE TABLE "dna_snp" (
                id          INTEGER  PRIMARY KEY,
                snp_str     TEXT     NOT NULL,
                pos         INTEGER  NOT NULL,
                ref         TEXT     NOT NULL,
                alt         TEXT     NOT NULL,
                color       TEXT     NOT NULL,
                snv_name    TEXT     NOT NULL
            );
            """
        )
        df_to_sql(cur, dna_snp, "dna_snp", index_label="id")
        cur.execute('CREATE INDEX "ix_dna_snp_pos" ON "dna_snp"("pos");')

        # AA SNVs
        gene_aa_snp = process_aa_snvs(metadata_map["gene_aa_snp"], "gene", genes)
        cur.execute('DROP TABLE IF EXISTS "gene_aa_snp";')
        cur.execute(
            """
            CREATE TABLE "gene_aa_snp" (
                id          INTEGER  PRIMARY KEY,
                snp_str     TEXT     NOT NULL,
                gene        TEXT     NOT NULL,
                pos         INTEGER  NOT NULL,
                ref         TEXT     NOT NULL,
                alt         TEXT     NOT NULL,
                color       TEXT     NOT NULL,
                snv_name    TEXT     NOT NULL,
                nt_pos      INTEGER  NOT NULL
            );
            """
        )
        df_to_sql(cur, gene_aa_snp, "gene_aa_snp", index_label="id")
        cur.execute('CREATE INDEX "ix_gene_aa_snp_pos" ON "gene_aa_snp"("pos");')
        cur.execute('CREATE INDEX "ix_gene_aa_snp_nt_pos" ON "gene_aa_snp"("nt_pos");')
        cur.execute('CREATE INDEX "ix_gene_aa_snp_gene" ON "gene_aa_snp"("gene");')

        protein_aa_snp = process_aa_snvs(
            metadata_map["protein_aa_snp"], "protein", proteins
        )
        cur.execute('DROP TABLE IF EXISTS "protein_aa_snp";')
        cur.execute(
            """
            CREATE TABLE "protein_aa_snp" (
                id          INTEGER  PRIMARY KEY,
                snp_str     TEXT     NOT NULL,
                protein     TEXT     NOT NULL,
                pos         INTEGER  NOT NULL,
                ref         TEXT     NOT NULL,
                alt         TEXT     NOT NULL,
                color       TEXT     NOT NULL,
                snv_name    TEXT     NOT NULL,
                nt_pos      INTEGER  NOT NULL
            );
            """
        )
        df_to_sql(cur, protein_aa_snp, "protein_aa_snp", index_label="id")
        cur.execute('CREATE INDEX "ix_protein_aa_snp_pos" ON "protein_aa_snp"("pos");')
        cur.execute(
            'CREATE INDEX "ix_protein_aa_snp_nt_pos" ON "protein_aa_snp"("nt_pos");'
        )
        cur.execute(
            'CREATE INDEX "ix_protein_aa_snp_protein" ON "protein_aa_snp"("protein");'
        )

        print("done")

        print("Writing location map...", end="", flush=True)

        # Locations
        location_map = pd.read_json(data_path / "location_map.json")
        cur.execute('DROP TABLE IF EXISTS "location";')
        cur.execute(
            """
            CREATE TABLE "location" (
                id        INTEGER  PRIMARY KEY,
                region    TEXT     NOT NULL,
                country   TEXT     NOT NULL,
                division  TEXT     NOT NULL,
                location  TEXT     NOT NULL
            );
            """
        )
        df_to_sql(cur, location_map, "location", index_label="id")

        print("done")

        print("Writing groups and consensus SNVs...", end="", flush=True)

        # Consensus SNVs
        with (data_path / "group_consensus_snps.json").open("r") as fp:
            group_consensus_snps = json.loads(fp.read())

        snp_fields = ["dna", "gene_aa", "protein_aa"]
        for grouping in group_consensus_snps.keys():
            for snp_field in snp_fields:
                table_name = "{grouping}_consensus_{snp_field}_snp".format(
                    grouping=grouping, snp_field=snp_field
                )
                # Create tables
                cur.execute(
                    'DROP TABLE IF EXISTS "{table_name}";'.format(table_name=table_name)
                )
                cur.execute(
                    """
                    CREATE TABLE "{table_name}" (
                        name    TEXT     NOT NULL,
                        snp_id  INTEGER  NOT NULL
                    );
                    """.format(
                        table_name=table_name
                    )
                )

                # Collect tuples of (group, snp_id)
                group_snps = []
                for group in group_consensus_snps[grouping].keys():
                    for snp_id in group_consensus_snps[grouping][group][
                        snp_field + "_snp_ids"
                    ]:
                        group_snps.append((group, snp_id))

                group_snp_df = pd.DataFrame.from_records(
                    group_snps, columns=["name", "snp_id"]
                ).set_index("name")
                df_to_sql(cur, group_snp_df, table_name, index_label="name")

                cur.execute(
                    """
                    CREATE INDEX "idx_{table_name}_name" ON "{table_name}"("name");
                    """.format(
                        table_name=table_name
                    )
                )

        print("done")

        # SNV frequencies
        print("Writing group SNV frequencies...", end="", flush=True)
        with (data_path / "group_snv_frequencies.json").open("r") as fp:
            group_snv_frequencies = json.loads(fp.read())

        snp_fields = ["dna", "gene_aa", "protein_aa"]
        for grouping in group_snv_frequencies.keys():
            for snp_field in snp_fields:
                table_name = "{grouping}_frequency_{snp_field}_snp".format(
                    grouping=grouping, snp_field=snp_field
                )
                # Create tables
                cur.execute(
                    'DROP TABLE IF EXISTS "{table_name}";'.format(table_name=table_name)
                )
                cur.execute(
                    """
                    CREATE TABLE "{table_name}" (
                        name      TEXT     NOT NULL,
                        snv_id    INTEGER  NOT NULL,
                        count     INTEGER  NOT NULL,
                        fraction  REAL     NOT NULL
                    );
                    """.format(
                        table_name=table_name
                    )
                )
                group_snv_frequency_df = (
                    pd.DataFrame.from_records(
                        group_snv_frequencies[grouping][snp_field]
                    )
                    .rename(columns={"group": "name"})
                    .set_index("name")
                )
                df_to_sql(cur, group_snv_frequency_df, table_name, index_label="name")

                cur.execute(
                    """
                    CREATE INDEX "idx_{table_name}_name" ON "{table_name}"("name");
                    """.format(
                        table_name=table_name
                    )
                )

        # Grouping tables

        # Build colormaps
        for grouping in group_consensus_snps.keys():
            group_df = pd.DataFrame(
                {"name": list(group_consensus_snps[grouping].keys())}
            )
            group_df["color"] = group_df["name"].map(
                get_categorical_colormap(group_df["name"])
            )

            cur.execute('DROP TABLE IF EXISTS "{grouping}";'.format(grouping=grouping))
            cur.execute(
                """
                CREATE TABLE "{grouping}" (
                    id     INTEGER  PRIMARY KEY,
                    name   TEXT     NOT NULL,
                    color  TEXT     NOT NULL
                );
                """.format(
                    grouping=grouping
                )
            )
            df_to_sql(cur, group_df, grouping, index_label="id")
            cur.execute(
                """
                CREATE INDEX "idx_{grouping}_name" on "{grouping}"("name")
                """.format(
                    grouping=grouping
                )
            )

        print("done")

        print("Writing sequence metadata...", end="", flush=True)

        # Sequence metadata
        case_data = pd.read_json(data_path / "case_data.json")
        # case_data = case_data.set_index("Accession ID")
        case_data["collection_date"] = pd.to_datetime(case_data["collection_date"])
        case_data["submission_date"] = pd.to_datetime(case_data["submission_date"])
        # print(case_data.columns)

        # Make a column for each metadata field
        metadata_cols = []
        metadata_col_defs = ""
        for field in config["metadata_cols"].keys():
            metadata_col_defs += "{field} INTEGER NOT NULL,\n".format(field=field)
            metadata_cols.append(field)

        # Make a column for each grouping
        grouping_cols = []
        grouping_col_defs = ""
        for grouping in group_consensus_snps.keys():
            grouping_col_defs += "{grouping} TEXT NOT NULL,\n".format(grouping=grouping)
            grouping_cols.append(grouping)

        cur.execute('DROP TABLE IF EXISTS "metadata";')
        cur.execute(
            """
            CREATE TABLE "metadata" (
                id               INTEGER    PRIMARY KEY,
                "Accession ID"   TEXT       NOT NULL,
                collection_date  TIMESTAMP  NOT NULL,
                submission_date  TIMESTAMP  NOT NULL,
                {metadata_col_defs}
                {grouping_col_defs}
                location_id      INTEGER    NOT NULL
            );
            """.format(
                metadata_col_defs=metadata_col_defs, grouping_col_defs=grouping_col_defs
            )
        )
        df_to_sql(
            cur,
            case_data[
                ["Accession ID", "collection_date", "submission_date"]
                + metadata_cols
                + grouping_cols
                + ["location_id"]
            ],
            "metadata",
            index_label="id",
        )
        # Create indices
        cur.execute(
            'CREATE INDEX "ix_metadata_collection_date" ON "metadata"("collection_date");'
        )
        cur.execute(
            'CREATE INDEX "ix_metadata_submission_date" ON "metadata"("submission_date");'
        )
        cur.execute(
            'CREATE INDEX "ix_metadata_location_id" ON "metadata"("location_id");'
        )
        for field in config["metadata_cols"].keys():
            cur.execute(
                'CREATE INDEX "ix_metadata_{field}" ON "metadata"("{field}");'.format(
                    field=field
                )
            )

        print("done")

        print("Writing sequence SNVs...", end="", flush=True)

        # Sequence SNV data
        snp_fields = ["dna", "gene_aa", "protein_aa"]
        for snp_field in snp_fields:
            snp_col = snp_field + "_snp_str"
            table_name = "sequence_{snp_field}_snp".format(snp_field=snp_field)
            cur.execute(
                'DROP TABLE IF EXISTS "{table_name}";'.format(table_name=table_name)
            )
            cur.execute(
                """
                CREATE TABLE "{table_name}" (
                    sequence_id  INTEGER  NOT NULL,
                    snp_id       INTEGER  NOT NULL
                );
                """.format(
                    table_name=table_name
                )
            )
            snv_df = case_data[[snp_col]].explode(snp_col)
            snv_df = snv_df.loc[~snv_df[snp_col].isna()]
            df_to_sql(
                cur, snv_df, table_name, index_label="sequence_id",
            )
            cur.execute(
                'CREATE INDEX "ix_{table_name}_sequence_id" ON "{table_name}"("sequence_id");'.format(
                    table_name=table_name
                )
            )
            cur.execute(
                'CREATE INDEX "ix_{table_name}_snp_id" ON "{table_name}"("snp_id");'.format(
                    table_name=table_name
                )
            )

        print("done")

        print("Writing JSONs...", end="", flush=True)

        # Stats table
        cur.execute('DROP TABLE IF EXISTS "stats";')
        cur.execute(
            """
            CREATE TABLE "stats" (
                value JSON NOT NULL
            );
            """
        )

        stats = {
            "num_sequences": len(case_data),
            "data_date": datetime.date.today().isoformat(),
        }

        cur.execute(
            """
            INSERT INTO "stats" (value) VALUES (%s)
            """,
            [Json(stats)],
        )

        # Country score
        # Just dump this as a big JSON
        cur.execute('DROP TABLE IF EXISTS "country_score";')
        cur.execute(
            """
            CREATE TABLE "country_score" (
                value JSON NOT NULL
            );
            """
        )
        with (data_path / "country_score.json").open("r") as fp:
            country_score = json.loads(fp.read())

        cur.execute(
            """
            INSERT INTO "country_score" (value) VALUES (%s)
            """,
            [Json(country_score)],
        )

        # Geo select tree
        # Just dump this as a JSON string in a table
        cur.execute('DROP TABLE IF EXISTS "geo_select_tree";')
        cur.execute(
            """
            CREATE TABLE "geo_select_tree" (
                value JSON NOT NULL
            )
            """
        )
        with (data_path / "geo_select_tree.json").open("r") as fp:
            geo_select_tree = json.loads(fp.read())

        cur.execute(
            """
            INSERT INTO "geo_select_tree" (value) VALUES (%s)
            """,
            [Json(geo_select_tree)],
        )

        print("done")
