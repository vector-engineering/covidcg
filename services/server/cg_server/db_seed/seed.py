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
from .load_mutations import process_dna_mutations, process_aa_mutations

# root/services/server/cg_server/db_seed/seed.py
project_root = Path(__file__).parent.parent.parent.parent.parent
virus = Path(config["virus"])
data_path = Path(os.getenv("DATA_PATH", project_root / config["example_data_folder"]))
static_data_path = Path(
    os.getenv("STATIC_DATA_PATH", project_root / config["static_data_folder"])
)


genes = pd.read_json(str(static_data_path / "genes_processed.json")).set_index("name")
proteins = pd.read_json(str(static_data_path / "proteins_processed.json")).set_index(
    "name"
)

loc_levels = [
    "region",
    "country",
    "division",
    "location",
]


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
        cur.execute("CREATE EXTENSION IF NOT EXISTS intarray;")

        cur.execute("DROP EXTENSION IF EXISTS intarray;")
        cur.execute("CREATE EXTENSION intarray;")

        print("Writing metadata maps...", end="", flush=True)

        # Load metadata map
        with open(data_path / "metadata_map.json", "r") as fp:
            metadata_map = json.loads(fp.read())

        # Make a table for each metadata field
        for field in list(config["metadata_cols"].keys()) + loc_levels:
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
                CREATE INDEX "ix_{table_name}_def_id" ON "{table_name}"("id");
                """.format(
                    table_name=table_name
                )
            )

        print("done")

        print("Writing mutation maps...", end="", flush=True)

        # DNA mutations
        dna_mutation = process_dna_mutations(metadata_map["dna_mutation"])
        cur.execute('DROP TABLE IF EXISTS "dna_mutation";')
        cur.execute(
            """
            CREATE TABLE "dna_mutation" (
                id             INTEGER  PRIMARY KEY,
                mutation_str   TEXT     NOT NULL,
                pos            INTEGER  NOT NULL,
                ref            TEXT     NOT NULL,
                alt            TEXT     NOT NULL,
                color          TEXT     NOT NULL,
                mutation_name  TEXT     NOT NULL
            );
            """
        )
        df_to_sql(cur, dna_mutation, "dna_mutation", index_label="id")
        cur.execute('CREATE INDEX "ix_dna_mutation_pos" ON "dna_mutation"("pos");')

        # AA mutations
        gene_aa_mutation = process_aa_mutations(
            metadata_map["gene_aa_mutation"], "gene", genes
        )
        cur.execute('DROP TABLE IF EXISTS "gene_aa_mutation";')
        cur.execute(
            """
            CREATE TABLE "gene_aa_mutation" (
                id             INTEGER  PRIMARY KEY,
                mutation_str   TEXT     NOT NULL,
                gene           TEXT     NOT NULL,
                pos            INTEGER  NOT NULL,
                ref            TEXT     NOT NULL,
                alt            TEXT     NOT NULL,
                color          TEXT     NOT NULL,
                mutation_name  TEXT     NOT NULL,
                nt_pos         INTEGER  NOT NULL
            );
            """
        )
        df_to_sql(cur, gene_aa_mutation, "gene_aa_mutation", index_label="id")
        cur.execute(
            'CREATE INDEX "ix_gene_aa_mutation_pos" ON "gene_aa_mutation"("pos");'
        )
        cur.execute(
            'CREATE INDEX "ix_gene_aa_mutation_nt_pos" ON "gene_aa_mutation"("nt_pos");'
        )
        cur.execute(
            'CREATE INDEX "ix_gene_aa_mutation_gene" ON "gene_aa_mutation"("gene");'
        )

        protein_aa_mutation = process_aa_mutations(
            metadata_map["protein_aa_mutation"], "protein", proteins
        )
        cur.execute('DROP TABLE IF EXISTS "protein_aa_mutation";')
        cur.execute(
            """
            CREATE TABLE "protein_aa_mutation" (
                id             INTEGER  PRIMARY KEY,
                mutation_str   TEXT     NOT NULL,
                protein        TEXT     NOT NULL,
                pos            INTEGER  NOT NULL,
                ref            TEXT     NOT NULL,
                alt            TEXT     NOT NULL,
                color          TEXT     NOT NULL,
                mutation_name  TEXT     NOT NULL,
                nt_pos         INTEGER  NOT NULL
            );
            """
        )
        df_to_sql(cur, protein_aa_mutation, "protein_aa_mutation", index_label="id")
        cur.execute(
            'CREATE INDEX "ix_protein_aa_mutation_pos" ON "protein_aa_mutation"("pos");'
        )
        cur.execute(
            'CREATE INDEX "ix_protein_aa_mutation_nt_pos" ON "protein_aa_mutation"("nt_pos");'
        )
        cur.execute(
            'CREATE INDEX "ix_protein_aa_mutation_protein" ON "protein_aa_mutation"("protein");'
        )

        print("done")

        # mutation frequencies
        print("Writing group mutation frequencies...", end="", flush=True)
        with (data_path / "group_mutation_frequencies.json").open("r") as fp:
            group_mutation_frequencies = json.loads(fp.read())

        mutation_fields = ["dna", "gene_aa", "protein_aa"]
        for grouping in group_mutation_frequencies.keys():
            for mutation_field in mutation_fields:
                table_name = "{grouping}_frequency_{mutation_field}_mutation".format(
                    grouping=grouping, mutation_field=mutation_field
                )
                # Create tables
                cur.execute(
                    'DROP TABLE IF EXISTS "{table_name}";'.format(table_name=table_name)
                )
                cur.execute(
                    """
                    CREATE TABLE "{table_name}" (
                        name         TEXT     NOT NULL,
                        mutation_id  INTEGER  NOT NULL,
                        count        INTEGER  NOT NULL,
                        fraction     REAL     NOT NULL
                    );
                    """.format(
                        table_name=table_name
                    )
                )
                group_mutation_frequency_df = (
                    pd.DataFrame.from_records(
                        group_mutation_frequencies[grouping][mutation_field]
                    )
                    .rename(columns={"group": "name"})
                    .set_index("name")
                )
                df_to_sql(
                    cur, group_mutation_frequency_df, table_name, index_label="name"
                )

                cur.execute(
                    """
                    CREATE INDEX "idx_{table_name}_name" ON "{table_name}"("name");
                    """.format(
                        table_name=table_name
                    )
                )
        print("done")

        # Grouping tables
        print("Writing groups...", end="", flush=True)
        with (data_path / "global_group_counts.json").open("r") as fp:
            global_group_counts = json.loads(fp.read())

        # Build colormaps
        for grouping in config["group_cols"].keys():
            group_df = pd.DataFrame(
                {"name": list(global_group_counts[grouping].keys())}
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

        # Partition settings
        delta = 31
        if config["mutation_partition_break"] == "Y":
            delta *= 12

        min_date = case_data["collection_date"].min()
        max_date = case_data["collection_date"].max() + pd.Timedelta(delta, unit="D")

        partition_dates = [
            d.isoformat()
            for d in pd.period_range(
                start=min_date, end=max_date, freq=config["mutation_partition_break"],
            )
            .to_timestamp()
            .date
        ]

        # Make a column for each metadata field
        metadata_cols = []
        metadata_col_defs = []
        for field in list(config["metadata_cols"].keys()) + loc_levels:
            metadata_col_defs.append("{field} INTEGER NOT NULL".format(field=field))
            metadata_cols.append(field)
        metadata_col_defs = ",\n".join(metadata_col_defs)

        # Make a column for each grouping
        grouping_cols = []
        grouping_col_defs = []
        for grouping in config["group_cols"].keys():
            grouping_col_defs.append(
                "{grouping} TEXT NOT NULL".format(grouping=grouping)
            )
            grouping_cols.append(grouping)
        grouping_col_defs = ",\n".join(grouping_col_defs)

        cur.execute('DROP TABLE IF EXISTS "metadata";')
        cur.execute(
            """
            CREATE TABLE "metadata" (
                sequence_id               INTEGER    NOT NULL,
                "Accession ID"   TEXT       NOT NULL,
                collection_date  TIMESTAMP  NOT NULL,
                submission_date  TIMESTAMP  NOT NULL,
                {metadata_col_defs},
                {grouping_col_defs}
            )
            PARTITION BY RANGE(collection_date);
            """.format(
                metadata_col_defs=metadata_col_defs, grouping_col_defs=grouping_col_defs
            )
        )
        for i in range(len(partition_dates) - 1):
            start = partition_dates[i]
            end = partition_dates[i + 1]
            cur.execute(
                """
                CREATE TABLE "{partition_table_name}" PARTITION OF "metadata"
                FOR VALUES FROM ('{start}') TO ('{end}')
                WITH (parallel_workers = 4);
                """.format(
                    partition_table_name="metadata_{}".format(i), start=start, end=end,
                )
            )

        df_to_sql(
            cur,
            case_data[
                ["Accession ID", "collection_date", "submission_date"]
                + metadata_cols
                + grouping_cols
            ],
            "metadata",
            index_label="sequence_id",
        )
        # Create indices
        cur.execute(
            'CREATE INDEX "ix_metadata_sequence_id" ON "metadata"("sequence_id");'
        )
        cur.execute(
            'CREATE INDEX "ix_metadata_collection_date" ON "metadata"("collection_date");'
        )
        cur.execute(
            'CREATE INDEX "ix_metadata_submission_date" ON "metadata"("submission_date");'
        )
        for field in list(config["metadata_cols"].keys()) + loc_levels:
            cur.execute(
                'CREATE INDEX "ix_metadata_{field}" ON "metadata"("{field}");'.format(
                    field=field
                )
            )
        for grouping in config["group_cols"].keys():
            cur.execute(
                'CREATE INDEX "ix_metadata_{grouping}" ON "metadata"("{grouping}");'.format(
                    grouping=grouping
                )
            )

        print("done")

        print("Writing sequence mutations...", end="", flush=True)

        # Sequence mutation data
        mutation_fields = ["dna", "gene_aa", "protein_aa"]
        for mutation_field in mutation_fields:
            mutation_col = mutation_field + "_mutation_str"
            table_name = "sequence_{mutation_field}_mutation".format(
                mutation_field=mutation_field
            )
            cur.execute(
                'DROP TABLE IF EXISTS "{table_name}";'.format(table_name=table_name)
            )
            cur.execute(
                """
                CREATE TABLE "{table_name}" (
                    sequence_id      INTEGER    NOT NULL,
                    collection_date  TIMESTAMP  NOT NULL,
                    submission_date  TIMESTAMP  NOT NULL,
                    {metadata_col_defs},
                    {grouping_col_defs},
                    mutations        INTEGER[]  NOT NULL
                )
                PARTITION BY RANGE(collection_date);
                """.format(
                    table_name=table_name,
                    metadata_col_defs=metadata_col_defs,
                    grouping_col_defs=grouping_col_defs,
                )
            )

            for i in range(len(partition_dates) - 1):
                start = partition_dates[i]
                end = partition_dates[i + 1]
                cur.execute(
                    """
                    CREATE TABLE "{partition_table_name}" PARTITION OF "{table_name}"
                    FOR VALUES FROM ('{start}') TO ('{end}')
                    WITH (parallel_workers = 4);
                    """.format(
                        partition_table_name="{}_{}".format(table_name, i),
                        table_name=table_name,
                        start=start,
                        end=end,
                    )
                )

            mutation_df = case_data[
                ["collection_date", "submission_date"]
                + metadata_cols
                + grouping_cols
                + [mutation_col]
            ]
            mutation_df = mutation_df.loc[~mutation_df[mutation_col].isna()]
            # Serialize list of integers
            mutation_df.loc[:, mutation_col] = mutation_df[mutation_col].apply(
                lambda x: "{" + ",".join([str(_x) for _x in sorted(x)]) + "}"
            )

            df_to_sql(
                cur, mutation_df, table_name, index_label="sequence_id",
            )

            # Create indices
            cur.execute(
                'CREATE INDEX "ix_{table_name}_sequence_id" ON "{table_name}"("sequence_id");'.format(
                    table_name=table_name
                )
            )
            cur.execute(
                'CREATE INDEX "ix_{table_name}_collection_date" ON "{table_name}"("collection_date");'.format(
                    table_name=table_name
                )
            )
            cur.execute(
                'CREATE INDEX "ix_{table_name}_submission_date" ON "{table_name}"("submission_date");'.format(
                    table_name=table_name
                )
            )
            for field in list(config["metadata_cols"].keys()) + loc_levels:
                cur.execute(
                    'CREATE INDEX "ix_{table_name}_{field}" ON "{table_name}"("{field}");'.format(
                        table_name=table_name, field=field
                    )
                )
            for grouping in config["group_cols"].keys():
                cur.execute(
                    'CREATE INDEX "ix_{table_name}_{grouping}" ON "{table_name}"("{grouping}");'.format(
                        table_name=table_name, grouping=grouping
                    )
                )

        print("done")

        print("Writing JSONs...", end="", flush=True)

        # Stats table
        cur.execute('DROP TABLE IF EXISTS "jsons";')
        cur.execute(
            """
            CREATE TABLE "jsons" (
                key TEXT PRIMARY KEY,
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
            INSERT INTO "jsons" (key, value) VALUES (%s, %s);
            """,
            ["stats", Json(stats)],
        )
        if config["virus"] == "sars2":
            with (data_path / "country_score.json").open("r") as fp:
                country_score = json.loads(fp.read())
            cur.execute(
                """
                INSERT INTO "jsons" (key, value) VALUES (%s, %s);
                """,
                ["country_score", Json(country_score)],
            )
        with (data_path / "geo_select_tree.json").open("r") as fp:
            geo_select_tree = json.loads(fp.read())
        cur.execute(
            """
            INSERT INTO "jsons" (key, value) VALUES (%s, %s);
            """,
            ["geo_select_tree", Json(geo_select_tree)],
        )
        with (data_path / "surveillance" / "group_counts2.json").open("r") as fp:
            surv_group_counts = json.loads(fp.read())
        cur.execute(
            """
            INSERT INTO "jsons" (key, value) VALUES (%s, %s);
            """,
            ["surv_group_counts", Json(surv_group_counts)],
        )
        with (data_path / "surveillance" / "group_regression2.json").open("r") as fp:
            surv_group_regression = json.loads(fp.read())
        cur.execute(
            """
            INSERT INTO "jsons" (key, value) VALUES (%s, %s);
            """,
            ["surv_group_regression", Json(surv_group_regression)],
        )

        if config["virus"] == "sars2":
            with (data_path / "vocs" / "vocs.json").open("r") as fp:
                vocs = json.loads(fp.read())
            cur.execute(
                """
                INSERT INTO "jsons" (key, value) VALUES (%s, %s);
                """,
                ["vocs", Json(vocs)],
            )

        # Metadata map
        table_queries = []
        for mutation_type in ["dna", "gene_aa", "protein_aa"]:
            table_queries.append(
                sql.SQL(
                    """
                    SELECT {mutation_type_name} as "field", "mutation_str" as "id", "id"::text as "value"
                    FROM {mutation_type}
                    """
                ).format(
                    mutation_type_name=sql.Literal(mutation_type + "_mutation"),
                    mutation_type=sql.Identifier(mutation_type + "_mutation"),
                )
            )
        table_queries = sql.SQL(" UNION ALL ").join(table_queries)
        cur.execute(
            sql.SQL(
                """
                SELECT "field", json_object_agg("id", "value") as "map"
                FROM (
                    {table_queries}
                ) a
                GROUP BY "field";
                """
            ).format(table_queries=table_queries)
        )
        metadata_map = dict(cur.fetchall())
        cur.execute(
            """
            INSERT INTO "jsons" (key, value) VALUES (%s, %s);
            """,
            ["metadata_map", Json(metadata_map)],
        )

        # Groups
        group_queries = []
        for group in config["group_cols"].keys():
            group_queries.append(
                sql.SQL(
                    """
                SELECT
                    "name",
                    "color",
                    {group_name} as "group"
                FROM {group_table}
                """
                ).format(
                    group_name=sql.Literal(group), group_table=sql.Identifier(group)
                )
            )

        cur.execute(sql.SQL(" UNION ALL ").join(group_queries))
        groups = [
            {"name": rec[0], "color": rec[1], "group": rec[2]} for rec in cur.fetchall()
        ]
        cur.execute(
            """
            INSERT INTO "jsons" (key, value) VALUES (%s, %s);
            """,
            ["groups", Json(groups)],
        )

        print("done")
