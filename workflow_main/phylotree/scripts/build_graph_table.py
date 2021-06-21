# coding: utf-8

"""Build tabular data from treetime output

Code is derived from CoVizu:
  * https://filogeneti.ca/covizu/
  * https://github.com/PoonLab/covizu

Please find the attached license at LICENSE_COVIZU

Author: Albert Chen - Vector Engineering Team (chena@broadinstitute.org)
"""

import networkx as nx
import pandas as pd

from Bio import Phylo


def build_graph_table(
    case_data_path, nexus_file_path, representative_table_path, graph_table_out_path
):
    df = pd.read_json(case_data_path)

    # Load tree
    phy = Phylo.read(nexus_file_path, format="nexus")
    net = Phylo.to_networkx(phy)

    reps = pd.read_csv(representative_table_path, index_col="Accession ID")

    nodes = pd.DataFrame(
        [(node.name, node.branch_length) for node in net.nodes],
        columns=["name", "branch_length"],
    )
    nodes = nodes.join(reps[["lineage"]], on="name")
    nodes["leaf"] = ~nodes["lineage"].isna()

    edges = pd.DataFrame(
        [(edge[0].name, edge[1].name) for edge in list(net.edges)],
        columns=["parent", "child"],
    )

    nodes = nodes.set_index("name").join(edges.set_index("child"))
    nodes["cumulative_length"] = nodes["branch_length"]

    # Recursively add length to children
    def add_length(node_name):
        children = nodes.loc[nodes["parent"] == node_name]

        for child_name, child in children.iterrows():
            nodes.loc[child_name, "cumulative_length"] += nodes.loc[
                node_name, "cumulative_length"
            ]
            add_length(child_name)

    add_length("NODE_0000001")

    # Recursively add counts to children
    counts_per_lineage = df["lineage"].value_counts()
    nodes["num_seqs"] = nodes["lineage"].map(counts_per_lineage).fillna(0).astype(int)

    def add_seqs_to_parent(node_name):
        num_seqs = nodes.at[node_name, "num_seqs"]
        children = nodes.loc[nodes["parent"] == node_name]

        for child_name, child in children.iterrows():
            num_seqs += add_seqs_to_parent(child_name)

        nodes.at[node_name, "num_seqs"] = num_seqs
        return num_seqs

    add_seqs_to_parent("NODE_0000001")

    nodes.reset_index().to_json(graph_table_out_path, orient="records")
