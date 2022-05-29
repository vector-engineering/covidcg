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
data_path = Path(os.getenv("DATA_PATH", project_root / config["example_data_folder"]))
static_data_path = Path(
    os.getenv("STATIC_DATA_PATH", project_root / config["static_data_folder"])
)

with open(str(static_data_path / "genes_processed.json")) as f:
    genes = json.load(f)
with open(str(static_data_path / "proteins_processed.json")) as f:
    proteins = json.load(f)
with open(str(static_data_path / "reference.json")) as f:
    references = json.load(f)

loc_levels = [
    "region",
    "country",
    "division",
    "location",
]


def df_to_sql(cur, df, table, index_label=None):
    n = 500  # chunk row size
    list_df = [df[i : i + n] for i in range(0, df.shape[0], n)]

    for chunk in list_df:
        buffer = io.StringIO()
        if index_label:
            chunk.to_csv(buffer, index_label=index_label, header=False)
        else:
            chunk.to_csv(buffer, index=False, header=False)

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

        cur.execute("DROP EXTENSION IF EXISTS intarray;")
        cur.execute("CREATE EXTENSION IF NOT EXISTS intarray;")

        print("Writing metadata maps...", end="", flush=True)

        # Load metadata map
        with open(data_path / "metadata_map.json", "r") as fp:
            metadata_map = json.loads(fp.read())

        # Make a table for each metadata field
        for field in list(config["metadata_cols"].keys()) + loc_levels:
            table_name = "metadata_{}".format(field)
            cur.execute(
                sql.SQL(
                    """
                DROP TABLE IF EXISTS {table_name};
                CREATE TABLE {table_name} (
                    id INTEGER PRIMARY KEY,
                    value TEXT NOT NULL
                );
            """
                ).format(table_name=sql.Identifier(table_name))
            )
            metadata_df = pd.DataFrame.from_records(
                list(metadata_map[field].items()), columns=["id", "value"]
            ).set_index("id")
            df_to_sql(cur, metadata_df, table_name, index_label="id")
            cur.execute(
                sql.SQL(
                    """
                CREATE INDEX {index_name} ON {table_name}("id");
                """
                ).format(
                    index_name=sql.Identifier(f"ix_{table_name}_def_id"),
                    table_name=sql.Identifier(table_name),
                )
            )

        print("done")

        print("Writing mutation maps...", end="", flush=True)

        # DNA mutations
        dna_mutation = process_dna_mutations(metadata_map["dna_mutation"])
        cur.execute(
            """
            DROP TABLE IF EXISTS "dna_mutation";
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
            metadata_map["gene_aa_mutation"], "gene"
        )
        cur.execute(
            """
            DROP TABLE IF EXISTS "gene_aa_mutation";
            CREATE TABLE "gene_aa_mutation" (
                id             INTEGER  PRIMARY KEY,
                mutation_str   TEXT     NOT NULL,
                gene           TEXT     NOT NULL,
                pos            INTEGER  NOT NULL,
                ref            TEXT     NOT NULL,
                alt            TEXT     NOT NULL,
                color          TEXT     NOT NULL,
                mutation_name  TEXT     NOT NULL
            );
            """
        )
        df_to_sql(cur, gene_aa_mutation, "gene_aa_mutation", index_label="id")
        cur.execute(
            """
            CREATE INDEX "ix_gene_aa_mutation_pos" ON "gene_aa_mutation"("pos");
            CREATE INDEX "ix_gene_aa_mutation_gene" ON "gene_aa_mutation"("gene");
            """
        )

        protein_aa_mutation = process_aa_mutations(
            metadata_map["protein_aa_mutation"], "protein"
        )
        cur.execute(
            """
            DROP TABLE IF EXISTS "protein_aa_mutation";
            CREATE TABLE "protein_aa_mutation" (
                id             INTEGER  PRIMARY KEY,
                mutation_str   TEXT     NOT NULL,
                protein        TEXT     NOT NULL,
                pos            INTEGER  NOT NULL,
                ref            TEXT     NOT NULL,
                alt            TEXT     NOT NULL,
                color          TEXT     NOT NULL,
                mutation_name  TEXT     NOT NULL
            );
            """
        )
        df_to_sql(cur, protein_aa_mutation, "protein_aa_mutation", index_label="id")
        cur.execute(
            """
            CREATE INDEX "ix_protein_aa_mutation_pos" ON "protein_aa_mutation"("pos");
            CREATE INDEX "ix_protein_aa_mutation_protein" ON "protein_aa_mutation"("protein");
            """
        )

        print("done")

        # mutation frequencies
        print("Writing group mutation frequencies...", end="", flush=True)
        with (data_path / "group_mutation_frequencies.json").open("r") as fp:
            group_mutation_frequencies = json.loads(fp.read())

        mutation_fields = ["dna", "gene_aa", "protein_aa"]
        for grouping in group_mutation_frequencies.keys():

            # Get references
            reference_names = sorted(group_mutation_frequencies[grouping].keys())

            for mutation_field in mutation_fields:
                table_name = f"{grouping}_frequency_{mutation_field}_mutation"
                # Create tables
                cur.execute(
                    sql.SQL(
                        """
                    DROP TABLE IF EXISTS {table_name};
                    CREATE TABLE {table_name} (
                        name         TEXT     NOT NULL,
                        reference    TEXT     NOT NULL,
                        mutation_id  INTEGER  NOT NULL,
                        count        INTEGER  NOT NULL,
                        fraction     REAL     NOT NULL
                    );
                    """
                    ).format(table_name=sql.Identifier(table_name))
                )

                group_mutation_frequency_df = []
                for reference_name in reference_names:
                    group_mutation_frequency_df.append(
                        pd.DataFrame.from_records(
                            group_mutation_frequencies[grouping][reference_name][
                                mutation_field
                            ]
                        )
                        .rename(columns={"group": "name"})
                        .assign(reference=lambda _: reference_name)
                    )
                group_mutation_frequency_df = pd.concat(
                    group_mutation_frequency_df, axis=0
                )[["name", "reference", "mutation_id", "count", "fraction"]]

                df_to_sql(cur, group_mutation_frequency_df, table_name)

                for field in ["name", "reference", "mutation_id", "fraction"]:
                    cur.execute(
                        sql.SQL(
                            """
                        CREATE INDEX {index_name} ON {table_name}({field});
                        """
                        ).format(
                            index_name=sql.Identifier(f"idx_{table_name}_{field}"),
                            table_name=sql.Identifier(table_name),
                            field=sql.Identifier(field),
                        )
                    )

        print("done")

        # Grouping tables
        print("Writing groups...", end="", flush=True)
        with (data_path / "global_group_counts.json").open("r") as fp:
            global_group_counts = json.loads(fp.read())

        # Build colormaps
        for grouping in config["group_cols"].keys():

            # Collect unique group names
            group_names = []
            for reference in global_group_counts.keys():
                group_names.extend(
                    list(global_group_counts[reference][grouping].keys())
                )
            group_names = sorted(set(group_names))

            group_df = pd.DataFrame({"name": group_names})
            group_df["color"] = group_df["name"].map(
                get_categorical_colormap(group_df["name"])
            )

            cur.execute(
                sql.SQL(
                    """
                DROP TABLE IF EXISTS {grouping};
                CREATE TABLE {grouping} (
                    id     INTEGER  PRIMARY KEY,
                    name   TEXT     NOT NULL,
                    color  TEXT     NOT NULL
                );
                """
                ).format(grouping=sql.Identifier(grouping))
            )
            df_to_sql(cur, group_df, grouping, index_label="id")
            cur.execute(
                sql.SQL(
                    """
                CREATE INDEX {index_name} on {grouping}("name")
                """
                ).format(
                    index_name=sql.Identifier(f"idx_{grouping}_name"),
                    grouping=sql.Identifier(grouping),
                )
            )
        print("done")

        print("Writing sequence metadata...", end="", flush=True)
        # Make a column for each metadata field
        metadata_cols = []
        metadata_col_defs = []
        for field in list(config["metadata_cols"].keys()) + loc_levels:
            metadata_col_defs.append(sql.SQL(f"{field} INTEGER NOT NULL"))
            metadata_cols.append(field)
        metadata_col_defs = sql.SQL(",\n").join(metadata_col_defs)

        # Make a column for each grouping
        grouping_cols = []
        grouping_col_defs = []
        for grouping in config["group_cols"].keys():
            grouping_col_defs.append(sql.SQL(f"{grouping} TEXT NOT NULL"))
            grouping_cols.append(grouping)
        grouping_col_defs = sql.SQL(",\n").join(grouping_col_defs)

        cur.execute(
            sql.SQL(
                """
            DROP TABLE IF EXISTS "metadata";
            CREATE TABLE "metadata" (
                sequence_id      INTEGER    NOT NULL,
                "Accession ID"   TEXT       NOT NULL,
                collection_date  TIMESTAMP  NOT NULL,
                submission_date  TIMESTAMP  NOT NULL,
                {metadata_col_defs},
                {grouping_col_defs}
            )
            PARTITION BY RANGE(collection_date);
            """
            ).format(
                metadata_col_defs=metadata_col_defs, grouping_col_defs=grouping_col_defs
            )
        )

        # Only save what we want from the case_data json
        mutation_fields = ["dna", "gene_aa", "protein_aa"]
        mutation_cols = []
        for mutation_field in mutation_fields:
            mutation_cols.append(mutation_field + "_mutation_str")
        case_data_columns = (
            ["sequence_id", "Accession ID", "collection_date", "submission_date"]
            + metadata_cols
            + grouping_cols
            + ["reference",]
            + mutation_cols
        )

        case_data = pd.read_json(data_path / "case_data.json")[case_data_columns]
        # case_data = case_data.set_index("Accession ID")
        case_data["collection_date"] = pd.to_datetime(case_data["collection_date"])
        case_data["submission_date"] = pd.to_datetime(case_data["submission_date"])
        # print(case_data.columns)

        # Partition settings
        min_date = case_data["collection_date"].min()
        # Round latest sequence to the nearest partition break
        if config["mutation_partition_break"] == "M":
            max_date = (
                case_data["collection_date"].max().normalize() + pd.offsets.MonthBegin()
            )
        elif config["mutation_partition_break"] == "Y":
            max_date = (
                case_data["collection_date"].max().normalize() + pd.offsets.YearBegin()
            )

        partition_dates = [
            d.isoformat()
            for d in pd.period_range(
                start=min_date.to_period(config["mutation_partition_break"]),
                end=max_date.to_period(config["mutation_partition_break"]),
                freq=config["mutation_partition_break"],
            )
            .to_timestamp()
            .date
        ]

        for i in range(len(partition_dates) - 1):
            date_start = partition_dates[i]
            date_end = partition_dates[i + 1]
            cur.execute(
                sql.SQL(
                    """
                CREATE TABLE {partition_table_name} PARTITION OF "metadata"
                FOR VALUES FROM ({date_start}) TO ({date_end})
                WITH (parallel_workers = 4);
                """
                ).format(
                    partition_table_name=sql.Identifier(f"metadata_{i}"),
                    date_start=sql.Literal(date_start),
                    date_end=sql.Literal(date_end),
                )
            )

        df_to_sql(
            cur,
            # Remove duplicates associated with multiple references
            case_data[
                ["sequence_id", "Accession ID", "collection_date", "submission_date"]
                + metadata_cols
                + grouping_cols
            ].drop_duplicates("Accession ID", keep="first"),
            "metadata",
        )
        # Create indices
        for field in (
            ["sequence_id", "collection_date", "submission_date"]
            + list(config["metadata_cols"].keys())
            + loc_levels
            + list(config["group_cols"].keys())
        ):
            cur.execute(
                sql.SQL(
                    """
                CREATE INDEX {index_name} ON "metadata"({field});
                """
                ).format(
                    index_name=sql.Identifier(f"idx_metadata_{field}"),
                    field=sql.Identifier(field),
                )
            )
        print("done")

        print("Writing sequence mutations...", end="", flush=True)

        # Sequence mutation data
        mutation_fields = ["dna", "gene_aa", "protein_aa"]
        for mutation_field in mutation_fields:
            mutation_col = mutation_field + "_mutation_str"
            table_name = f"sequence_{mutation_field}_mutation"
            cur.execute(
                sql.SQL(
                    """
                DROP TABLE IF EXISTS {table_name};
                CREATE TABLE {table_name} (
                    sequence_id      INTEGER    NOT NULL,
                    collection_date  TIMESTAMP  NOT NULL,
                    submission_date  TIMESTAMP  NOT NULL,
                    {metadata_col_defs},
                    {grouping_col_defs},
                    reference        TEXT       NOT NULL,
                    mutations        INTEGER[]  NOT NULL
                )
                PARTITION BY LIST(reference);
                """
                ).format(
                    table_name=sql.Identifier(table_name),
                    metadata_col_defs=metadata_col_defs,
                    grouping_col_defs=grouping_col_defs,
                )
            )

            reference_names = list(references.keys())

            # Create table partitions
            for reference_name in reference_names:
                # Clean up the reference name as a SQL ident - no dots
                reference_name_sql = reference_name.replace(".", "_")

                reference_partition_name = f"{table_name}_{reference_name_sql}"

                # Create reference partition
                cur.execute(
                    sql.SQL(
                        """
                    CREATE TABLE {reference_partition_name} PARTITION OF {table_name}
                    FOR VALUES IN ({reference_name})
                    PARTITION BY RANGE(collection_date);
                    """
                    ).format(
                        reference_partition_name=sql.Identifier(
                            reference_partition_name
                        ),
                        reference_name=sql.Literal(reference_name),
                        table_name=sql.Identifier(table_name),
                    )
                )

                for i in range(len(partition_dates) - 1):
                    date_start = partition_dates[i]
                    date_end = partition_dates[i + 1]
                    cur.execute(
                        sql.SQL(
                            """
                        CREATE TABLE {date_partition_name} PARTITION OF {reference_partition_name}
                        FOR VALUES FROM ({date_start}) TO ({date_end})
                        WITH (parallel_workers = 4);
                        """
                        ).format(
                            date_partition_name=sql.Identifier(
                                f"{table_name}_{reference_name_sql}_{i}"
                            ),
                            reference_partition_name=sql.Identifier(
                                reference_partition_name
                            ),
                            date_start=sql.Literal(date_start),
                            date_end=sql.Literal(date_end),
                        )
                    )

            mutation_df = case_data[
                ["sequence_id", "collection_date", "submission_date"]
                + metadata_cols
                + grouping_cols
                + ["reference", mutation_col]
            ]
            mutation_df = mutation_df.loc[~mutation_df[mutation_col].isna()]
            # Serialize list of integers
            mutation_df.loc[:, mutation_col] = mutation_df[mutation_col].apply(
                lambda x: "{" + ",".join([str(_x) for _x in sorted(x)]) + "}"
            )

            df_to_sql(cur, mutation_df, table_name)

            # Create indices
            for field in (
                ["sequence_id", "collection_date", "submission_date", "reference"]
                + list(config["metadata_cols"].keys())
                + loc_levels
                + list(config["group_cols"].keys())
            ):

                cur.execute(
                    sql.SQL(
                        "CREATE INDEX {index_name} ON {table_name}({field});"
                    ).format(
                        index_name=sql.Identifier(f"ix_{table_name}_{field}"),
                        table_name=sql.Identifier(table_name),
                        field=sql.Identifier(field),
                    )
                )

        print("done")

        print("Writing coverage tables...", end="", flush=True)

        mutation_fields = ["dna", "gene_aa", "protein_aa"]
        for mutation_field in mutation_fields:
            table_name = f"{mutation_field}_coverage"

            coverage_cols = ["sequence_id", "reference", "start", "end"]
            gene_or_protein_col = sql.SQL("")
            if mutation_field == "gene_aa":
                gene_or_protein_col = sql.SQL("gene TEXT NOT NULL,")
                coverage_cols.insert(2, "gene")
            elif mutation_field == "protein_aa":
                gene_or_protein_col = sql.SQL("protein TEXT NOT NULL,")
                coverage_cols.insert(2, "protein")

            cur.execute(
                sql.SQL(
                    """
                DROP TABLE IF EXISTS {table_name};
                CREATE TABLE {table_name} (
                    sequence_id  INTEGER  NOT NULL,
                    reference    TEXT     NOT NULL,
                    {gene_or_protein_col}
                    range_start  INTEGER  NOT NULL,
                    range_end    INTEGER  NOT NULL
                )
                PARTITION BY LIST(reference);
                """
                ).format(
                    table_name=sql.Identifier(table_name),
                    gene_or_protein_col=gene_or_protein_col,
                )
            )

            reference_names = list(references.keys())
            # Create table partitions
            for reference_name in reference_names:
                # Clean up the reference name as a SQL ident - no dots
                reference_name_sql = reference_name.replace(".", "_")

                # For gene/protein, additionally partition by gene/protein name
                additional_partition = ""
                if mutation_field == "gene_aa":
                    additional_partition = "PARTITION BY LIST(gene)"
                elif mutation_field == "protein_aa":
                    additional_partition = "PARTITION BY LIST(protein)"

                reference_partition_name = f"{table_name}_{reference_name_sql}"

                cur.execute(
                    sql.SQL(
                        """
                    CREATE TABLE {reference_partition_name} PARTITION OF {table_name}
                    FOR VALUES IN ({reference_name})
                    {additional_partition}
                    ;
                    """
                    ).format(
                        reference_partition_name=sql.Identifier(
                            reference_partition_name
                        ),
                        table_name=sql.Identifier(table_name),
                        reference_name=sql.Literal(reference_name),
                        additional_partition=sql.SQL(additional_partition),
                    )
                )

                # Create gene/protein name partitions
                if mutation_field in ["gene_aa", "protein_aa"]:
                    if mutation_field == "gene_aa":
                        features = genes[reference_name]
                    elif mutation_field == "protein_aa":
                        features = proteins[reference_name]

                    for feature in features:
                        cur.execute(
                            sql.SQL(
                                """
                            CREATE TABLE {gene_protein_partition_name} PARTITION OF {reference_partition_name}
                            FOR VALUES IN ({gene_protein})
                            WITH (parallel_workers = 4);
                            """
                            ).format(
                                gene_protein_partition_name=sql.Identifier(
                                    f'{reference_partition_name}_{feature["name"]}'
                                ),
                                reference_partition_name=sql.Identifier(
                                    reference_partition_name
                                ),
                                gene_protein=sql.Literal(feature["name"]),
                            )
                        )

            coverage_df = pd.read_csv(data_path / f"coverage_{mutation_field}.csv")
            df_to_sql(
                cur, coverage_df[coverage_cols], table_name,
            )

            # Create indices
            index_fields = ["sequence_id", "reference", "range_start", "range_end"]
            if mutation_field == "gene_aa":
                index_fields.append("gene")
            elif mutation_field == "protein_aa":
                index_fields.append("protein")

            for field in index_fields:
                cur.execute(
                    sql.SQL(
                        "CREATE INDEX {index_name} ON {table_name}({field});"
                    ).format(
                        index_name=sql.Identifier(f"ix_{table_name}_{field}"),
                        table_name=sql.Identifier(table_name),
                        field=sql.Identifier(field),
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
            "num_sequences": len(
                case_data.drop_duplicates("Accession ID", keep="first")
            ),
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
        cur.close()
