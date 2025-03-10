# type: ignore
# flake8: noqa
#
#
#
#
#
#
#
#| label: Read Data 
# Set up 
from pprint import pprint
import pandas as pd
from pathlib import Path
import numpy as np

# Get the current directory
# current_dir = pathlib.Path(__file__).parent.resolve()
current_directory = Path.cwd()

print(current_directory)

# Get folder names
dir_dat = current_directory /"P30 pilot"/ "0_data" / "clean_data"
dir_res = current_directory /"P30 pilot"/ "2_results"


# Read in proteomics metadata
# gene_data = pd.read_csv(dir_res / 'Combined_proteomics_scRNAseq_sig_results_073024.csv')
gene_data = pd.read_csv(dir_res / 'Combined_proteomics_scRNAseq_sig_fdr_results_090924.csv')

# Get names of genes in a list
gene_list = gene_data['EntrezGeneSymbol'].tolist()

# Combine Genes into str
genes = "[" + ", ".join(["'" + symbol + "'" for symbol in gene_list]) + "]"

overlapping_gene_list = ["TMEM87B", "THOC1", "SPRED1", "MFN1"]
genes_both_human_invitro = gene_data[gene_data['sig_overall_simplified'] == 'Both']['EntrezGeneSymbol'].tolist()
genes_human_only = gene_data[gene_data['sig_overall_simplified'] == 'TODAY']['EntrezGeneSymbol'].tolist()
genes_in_vitro_only = gene_data[gene_data['sig_overall_simplified'] == 'scRNAseq']['EntrezGeneSymbol'].tolist()

#print(split_genes)
#print(split_genes)
# type(gene_data)
#print(type(prot_metadata.EntrezGeneSymbol))
#
#
#
#
#
#
#| label: ComptoxAI query
# Connect to neo4j database
import comptox_ai
from comptox_ai.db.graph_db import GraphDB
db = GraphDB()
#db = GraphDB(username="neo4j_user", password="12345", hostname="localhost:7687")

# Create cypher query
# Start and end strings
start_string = "MATCH (c:Chemical {xrefDTXSID: 'DTXSID8031863'}) MATCH (d:Disease {commonName: 'Diabetic Nephropathy'}) MATCH (node1:Gene) WHERE node1.geneSymbol IN "

end_string = "WITH collect(id(node1))+collect(c)+collect(d) as nodes CALL apoc.algo.cover(nodes) YIELD rel RETURN  startNode(rel), rel, endNode(rel);"

# Combine the start and end strings with the unique_gene_symbols_list
query_string = start_string + "[" + ", ".join(["'" + symbol + "'" for symbol in gene_list]) + "] " + end_string

# Run Cypher Query
data = db.run_cypher(query_string)
# print(query_string)
#
#
#
#
#
#
#
#
#| label: Create network diagram
import networkx as nx
import matplotlib.pyplot as plt
# Create a new graph
G = nx.DiGraph()

# List mediator proteins

# Function to compute the combined 'type' attribute
def compute_type(node):
    # if node.get('commonName') == 'PFNA':
    if node.get('xrefDTXSID') == 'DTXSID8031863':
        return 'PFAS'
    if node.get('commonName') == 'Diabetic Nephropathy':
        return 'Disease'
    if node.get('geneSymbol') in genes_both_human_invitro:
        return 'Gene: Human and in-vitro'
    if node.get('geneSymbol') in genes_human_only:
        return 'Gene: Human only'
    if node.get('geneSymbol') in genes_in_vitro_only:
        return 'Gene: In-vitro only'
    else:
        return 'Gene: from ComptoxAI'

    
# Add nodes and edges with combined 'type' attribute
for entry in data:
    start_node = entry['startNode(rel)']
    end_node = entry['endNode(rel)']
    rel_type = entry['rel'][1]  # Relationship type is the second item in the tuple
    
    # Set combined 'type' attribute
    start_node['type'] = compute_type(start_node)
    end_node['type'] = compute_type(end_node)

    # Node identifiers
    start_node_id = start_node.get('geneSymbol')  or start_node.get('commonName')
    end_node_id = end_node.get('geneSymbol') or end_node.get('xrefUmlsCUI') 

    # Add nodes with combined 'type' attribute
    G.add_node(start_node_id, **start_node)
    G.add_node(end_node_id, **end_node)

    # Add edge
    G.add_edge(start_node_id, end_node_id, relationship=rel_type, weight = 1)

del data 
#
#
#
#
#
from collections import Counter

# Color mapping
color_map = {
    'PFAS': 'magenta',
    'Disease': 'red',
    'Gene: Human and in-vitro': 'green',
    'Gene: Human only': 'blue', 
    'Gene: In-vitro only': 'grey'
}

# Compute node colors based on 'type' attribute
node_colors = [color_map[G.nodes[node]['type']] for node in G]

print(Counter(node_colors))
print(len(G))
#print(len(node_colors))

# Draw the graph
#plt.figure(figsize=(12, 8))  # Set the figure size
#pos = nx.spring_layout(G)  # Layout for the nodes
#nx.draw_networkx(G, pos, with_labels=True, node_color=node_colors, node_size=700, edge_color='k', linewidths=1, font_size=10, arrows=True)
#nx.draw_networkx_edge_labels(G, pos, edge_labels=nx.get_edge_attributes(G, 'relationship'))

#plt.show()
#
#
#
#
#
#| label: Check self loops
# Check graph after adding metadata
isolated_nodes = list(nx.isolates(G))
print(nx.number_of_isolates(G))
print(len(G.nodes))

# Examine self loop edges:
self_loop_edges = list(nx.selfloop_edges(G))
# Extract nodes that have self-loops
self_loop_nodes = set(node for node, _ in self_loop_edges)
# Identify nodes that only have self-loops and no other connections
isolated_self_loop_nodes = [
    node for node in self_loop_nodes
    if G.degree(node) == 1  # Total degree (in-degree + out-degree) is 1, indicating only a self-loop
]
print(f"Nodes only connected to themselves: {isolated_self_loop_nodes}")
#
#
#
#
#| label: identify disconnected graphs
# Generate weakly connected components
weakly_connected_components = nx.weakly_connected_components(G)

# Get the size of each component (number of nodes)
component_sizes = [len(component) for component in weakly_connected_components]
# print(Counter(component_sizes))

# Identify genes not in the large componenet -----------

# Find all weakly connected components
weakly_connected_components = list(nx.weakly_connected_components(G))

# Find the largest weakly connected component (by number of nodes)
largest_component = max(weakly_connected_components, key=len)

# Find all nodes in the graph
all_nodes = set(G.nodes())

# Find nodes not in the largest component
nodes_not_in_largest = all_nodes - largest_component

# Prepare the data for the table
data = []
for node in nodes_not_in_largest:
    node_type = G.nodes[node].get('type', 'Unknown')  # Use 'Unknown' if no type is provided
    data.append((node, node_type))

# Create a DataFrame
df = pd.DataFrame(data, columns=['Node', 'Type'])

# Sort the DataFrame by the 'Type' column
df_sorted = df.sort_values(by='Type')

# Display the sorted DataFrame
#print(df_sorted)

# Remove these nodes from the graph
G_trim = G.copy()
G_trim.remove_nodes_from(nodes_not_in_largest)

print(len(G.nodes))
print(len(G_trim.nodes))
print(len(nodes_not_in_largest))
print(len(G_trim.nodes) == len(G.nodes)-len(nodes_not_in_largest))
print(nodes_not_in_largest)
#
#
#
#
#
#
#| label: compute leidenalg communities 
import igraph as ig
import leidenalg as la
import cairocffi as cairo
import matplotlib.pyplot as plt


source_node = "DTXSID8031863"
target_node = "C0011881"

# convert to igraph
h = ig.Graph.from_networkx(G_trim)
# Create a new attribute 'name' from '_nx_name'
h.vs['name'] = h.vs['_nx_name']

# Create a dictionary to map gene symbols to sig_overall values
gene_to_sig = dict(zip(gene_data['EntrezGeneSymbol'], gene_data['sig_overall']))

# Initialize the cluster assignments for nodes based on 'sig_overall'
initial_clusters = []
for vertex in h.vs:
    gene_symbol = vertex['name']
    if gene_symbol in gene_to_sig:
        initial_clusters.append(gene_to_sig[gene_symbol])
    elif gene_symbol == source_node:
        initial_clusters.append(source_node)
    elif gene_symbol == target_node:
        initial_clusters.append(target_node)    
    else:
        initial_clusters.append("error")  # Assign -1 if no match is found

# Convert the list of strings to a numeric list using pandas' factorize function
import pandas as pd
initial_clusters_num, _ = pd.factorize(pd.Series(initial_clusters))
initial_clusters = initial_clusters_num.tolist()

# Set the initial cluster assignment as a vertex attribute
h.vs['initial_cluster'] = initial_clusters
Counter(initial_clusters)

# Identify cluster partitions
partition = la.find_partition(graph=h, partition_type=la.ModularityVertexPartition, initial_membership=initial_clusters, n_iterations=-1, seed = 3787)
print(len(partition))
# Optimize partitions
optimiser = la.Optimiser()
diff = optimiser.optimise_partition(partition, n_iterations=1000)
print(len(partition))
print(diff)

#| label: Get leidenalg community assignments
leidenalg_community = partition.membership
# Create a concatenated variable that combines groups 3 and above
leidenalg_community_4 = [i if i < 3 else 3 for i in leidenalg_community]
leidenalg_community_8 = [i if i < 7 else 7 for i in leidenalg_community]

# Get info about source and target nodes 
source_vertex = h.vs.find(name=source_node)
target_vertex = h.vs.find(name=target_node)
source_cluster = partition.membership[source_vertex.index]
target_cluster = partition.membership[target_vertex.index]
source_cluster_size = Counter(leidenalg_community)[source_cluster]
target_cluster_size = Counter(leidenalg_community)[target_cluster]
print(f'PFNA Cluster: {source_cluster}, cluster size: {source_cluster_size}')
print(f'DKD Cluster: {target_cluster}, cluster size:  {target_cluster_size}')
Counter(leidenalg_community_4)
Counter(leidenalg_community)
# Create a concatenated variable that combines groups 3 and above

#
#
#
#
#| label: map with G_trim  
# Create a mapping from igraph vertices to NetworkX nodes
# Assuming nodes in G_trim are labeled from 0 to n-1
mapping = {v.index: node for v, node in zip(h.vs, G_trim.nodes)}

# Add leidenalg_community information to NetworkX graph
for idx, community in enumerate(leidenalg_community):
    nx_node = mapping[idx]
    G_trim.nodes[nx_node]['leidenalg_community'] = community

# Add leidenalg_community_red information to NetworkX graph
for idx, community in enumerate(leidenalg_community_4):
    nx_node = mapping[idx]
    G_trim.nodes[nx_node]['leidenalg_community_4'] = community

# Add leidenalg_community_red information to NetworkX graph
for idx, community in enumerate(leidenalg_community_8):
    nx_node = mapping[idx]
    G_trim.nodes[nx_node]['leidenalg_community_8'] = community

# Add "z" to G to control depth on Gephi figure
for idx, community in enumerate(leidenalg_community_4):
    nx_node = mapping[idx]
    if community == source_cluster: 
        G_trim.nodes[nx_node]['z'] = 2
        G_trim.nodes[nx_node]['clust_type1'] = "Source_Cluster"
    elif community == target_cluster: 
        G_trim.nodes[nx_node]['z'] = 1
        G_trim.nodes[nx_node]['clust_type1'] = "Target_Cluster"
    else:
        G_trim.nodes[nx_node]['z'] = 0
        G_trim.nodes[nx_node]['clust_type1'] = "Other"
#
#
#
#

# Create a subgraph with nodes that have 'z' == 1
# filtered_nodes = [node for node, data in G_trim.nodes(data=True) if data.get('z') == 1]
# G_subset = G_trim.subgraph(filtered_nodes)

# # Verify by printing only a few nodes
# for node, data in G_subset.nodes(data=True):
#     print(f"Node: {node}, z: {data.get('z')}, clust_type: {data.get('clust_type')}, clust: {data.get('leidenalg_community_min')}")
#
#
#
#
# Find all simple paths from source_node to target_node
all_paths = list(nx.shortest_path(G_trim, source=source_node, target=target_node))
print(all_paths)
#
#
#
#
#| label: Create DF of cluster results
# Create a DataFrame linking node names with the community identities
node_ids = list(G_trim.nodes)
df_embeddings = pd.DataFrame(node_ids, columns=['gene'])

# Merge with clusters 
df_embeddings['leidenalg_community']     = leidenalg_community
df_embeddings['leidenalg_community_min'] = leidenalg_community_min


# Before merge: record the number of rows in each dataframe
rows_before = {'df_embeddings': len(df_embeddings), 'gene_data': len(gene_data)}

# Merging with gene_data metadata, by 'gene' and 'EntrezGeneSymbol' columns
df_merged = pd.merge(df_embeddings, gene_data, how='left', left_on='gene', right_on='EntrezGeneSymbol')

# After merge: record the number of rows in the merged dataframe
rows_after = len(df_merged)

# Print the merge details
print(f"Number of rows in df_embeddings before merge: {rows_before['df_embeddings']}")
print(f"Number of rows in gene_data before merge: {rows_before['gene_data']}")
print(f"Number of rows in df_merged after merge: {rows_after}")
print(f"Number of rows matched: {rows_after}")
print(f"Number of rows lost from df_embeddings: {rows_before['df_embeddings'] - rows_after}")
print(f"Number of rows lost from original gene_data: {rows_before['gene_data'] - rows_after}")
# print(df_merged.columns)

df_merged.to_csv(dir_res / 'network_graph_analysis.csv', index=False)
#
#
#
#
#
#| label:  Cross tabs
cross_table = pd.crosstab(df_merged['sig_overall'], df_merged['leidenalg_community_min'])
# Convert counts to percentages
# cross_table_percent = cross_table.div(cross_table.sum(axis=1), axis=0) * 100
cross_table
#
#
#
#
#
#| label:  save graph

# Save file
# file_path = dir_res / 'ComptoxAI' / 'PFAS_prot_in_vitro_sig_fdr_trimmed_073024.graphml'
file_path = dir_res / 'ComptoxAI' / 'PFAS_prot_in_vitro_sig_100924.graphml'

# Write the graph to a GraphML file
nx.write_graphml(G_trim, file_path)
print(f"Network saved to {file_path}")
#
#
#
