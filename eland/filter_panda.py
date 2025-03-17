import os
import pandas as pd
import networkx as nx
import numpy as np
import scipy.sparse as sp
import communiy.community_louvain as community_louvain

def load_networks(prior_file, panda_file, delimiter='\t'):
    """
    Loads prior and PANDA networks from files as directed graphs.
    
    Parameters:
    prior_file (str): Path to the prior network edge list file.
    panda_file (str): Path to the PANDA network edge list file.
    delimiter (str): Delimiter used in the edge list files.
    
    Returns:
    tuple: A tuple containing the prior network graph and the PANDA network graph.
    """
    prior = nx.read_weighted_edgelist(prior_file, delimiter=delimiter, create_using=nx.DiGraph)
    panda = nx.read_weighted_edgelist(panda_file, delimiter=delimiter, create_using=nx.DiGraph)
    return prior, panda

def intersect_networks(prior, panda, prior_only=False):
    """
    Intersects the edges of prior and PANDA networks, creating a new directed graph with common edges.
    
    Parameters:
    prior (networkx.DiGraph): Directed graph for the prior network.
    panda (networkx.DiGraph): Directed graph for the PANDA network.
    
    Returns:
    scipy.sparse.csr_matrix: The adjacency matrix of the intersected network.
    list: The list of nodes in the intersected network.
    """
    common_edges = set(prior.edges()).intersection(panda.edges())
    weighted_prior = nx.DiGraph()
    
    for u, v in common_edges:
        if panda.has_edge(u, v):
            weighted_prior.add_edge(u, v, weight=panda[u][v]['weight'])
    
    if prior_only:
        return weighted_prior
    else:
        node_list = list(weighted_prior.nodes())
        adjacency_matrix = nx.adjacency_matrix(weighted_prior, nodelist=node_list)
        return adjacency_matrix, node_list

def filter_edges_by_max_values(adjacency_matrix):
    """
    Filters edges by retaining only those with weights greater than the minimum of the maximum weights per row or column.
    
    Parameters:
    adjacency_matrix (scipy.sparse.csr_matrix): The adjacency matrix of the network.
    
    Returns:
    scipy.sparse.csr_matrix: The filtered adjacency matrix.
    """
    max_per_tf = adjacency_matrix.max(axis=1).toarray().flatten()
    max_per_gene = adjacency_matrix.max(axis=0).toarray().flatten()
    min_edge = min(max_per_gene.min(), max_per_tf.min())
    selected_edges_mask = adjacency_matrix > min_edge
    selected_edges = adjacency_matrix.multiply(selected_edges_mask)
    return selected_edges

def convert_sparse_matrix_to_edgelist(selected_edges, node_list):
    """
    Converts a sparse matrix back to an edgelist and creates a DataFrame for better visualization or saving.
    
    Parameters:
    selected_edges (scipy.sparse.csr_matrix): The sparse matrix containing the selected edges.
    node_list (list): The list of node names corresponding to the indices in the sparse matrix.
    
    Returns:
    pd.DataFrame: A DataFrame representing the edgelist with columns 'source', 'target', and 'weight'.
    """
    selected_edges_coo = selected_edges.tocoo()
    edges_list = [
        (node_list[selected_edges_coo.row[i]], node_list[selected_edges_coo.col[i]], selected_edges_coo.data[i])
        for i in range(len(selected_edges_coo.data))
    ]
    edges_df = pd.DataFrame(edges_list, columns=['source', 'target', 'weight'])
    return edges_df

def calculate_modularity(edge_list):
    """
    Calculates the modularity of a network given its edge list using the Louvain method.
    
    Parameters:
    edge_list (pd.DataFrame): DataFrame representing the edgelist with columns 'source', 'target', and 'weight'.
    
    Returns:
    float: The modularity of the network.
    """
    # Ensure the edge list has at least three columns for source, target, and weight
    if edge_list.shape[1] < 3:
        raise ValueError("Edge list must have at least three columns for source, target, and weight.")
    
    # Create a graph from the edge list using the first two columns for source and target
    G = nx.from_pandas_edgelist(edge_list, source=edge_list.columns[0], target=edge_list.columns[1], edge_attr=edge_list.columns[2])
    
    # Detect communities using the Louvain method
    partition = community_louvain.best_partition(G)
    
    # Calculate the modularity of the network
    modularity_value = community_louvain.modularity(partition, G)
    
    return modularity_value

def generate_random_network(G):
    """
    Generates a random bipartite network with the same degree distribution as the given graph.
    
    Parameters:
    G (networkx.Graph): The bipartite graph.
    
    Returns:
    networkx.Graph: A random bipartite graph with the same degree distribution.
    """
    # Get the degree sequences for the two sets of nodes
    top_nodes = [n for n in G.nodes() if n.endswith('_TF')]
    bottom_nodes = [n for n in G.nodes() if n.endswith('_tar')]
    top_degrees = [G.degree(n) for n in top_nodes]
    bottom_degrees = [G.degree(n) for n in bottom_nodes]
    
    # Generate a random bipartite graph with the same degree distribution
    random_bipartite = nx.bipartite.configuration_model(top_degrees, bottom_degrees)
    random_bipartite = nx.Graph(random_bipartite)  # Convert to simple graph
    random_bipartite.remove_edges_from(nx.selfloop_edges(random_bipartite))  # Remove self-loops
    
    return random_bipartite


def filter_panda(prior_file, panda_file, output_file = None, delimiter='\t', compare_with_random = False, 
                prior_only=False):
    """
    Filters the PANDA network by processing prior and PANDA networks, finding their intersection,
    filtering the edges, and converting the filtered adjacency matrix to an edgelist DataFrame.
    
    Parameters:
    prior_file (str): Path to the prior network edge list file.
    panda_file (str): Path to the PANDA network edge list file.
    output_file (str): (Optional) Path to the output file for the filtered edgelist. If not provided, the DataFrame is returned.
    delimiter (str): Delimiter used in the edge list files.
    compare_with_random (bool): Whether to compare the modularity of the filtered network with a random network with preserved degree distribution.
    
    Returns:
    pd.DataFrame: A DataFrame representing the filtered edgelist with columns 'source', 'target', and 'weight'.
    """
    prior, panda = load_networks(prior_file, panda_file, delimiter)
    
    if prior_only:
        filtered_edges = intersect_networks(prior, panda, prior_only=True)
        edges_df = nx.to_pandas_edgelist(filtered_edges)
    else:
        adjacency_matrix, node_list = intersect_networks(prior, panda)
        filtered_edges = filter_edges_by_max_values(adjacency_matrix)
        edges_df = convert_sparse_matrix_to_edgelist(filtered_edges, node_list)
    
    if compare_with_random:
        # Append suffixes to node names to ensure bipartite structure
        edges_df['source'] = edges_df['source'].astype(str) + '_TF'
        edges_df['target'] = edges_df['target'].astype(str) + '_tar'
    
        # Create a bipartite graph from the edge list DataFrame
        G = nx.from_pandas_edgelist(edges_df, source='source', target='target')
        
        # Calculate modularity of the filtered PANDA network
        filtered_modularity = calculate_modularity(G)
        log_message = f"Modularity of the filtered PANDA network: {filtered_modularity}\n"
    
        # Generate a random bipartite network with the same degree distribution
        random_bipartite = generate_random_network(G)
    
        # Calculate modularity of the random network
        random_modularity = calculate_modularity(random_bipartite)
        log_message += f"Modularity of the random network: {random_modularity}\n"
        
        # Remove the extension from the output file name and add ".log"
        log_file = os.path.splitext(output_file)[0] + ".log"
        
        with open(log_file, 'a') as f:
                f.write(log_message)
        
        # Remove suffixes from node names
        edges_df['source'] = edges_df['source'].str.replace('_TF$', '', regex=True)
        edges_df['target'] = edges_df['target'].str.replace('_tar$', '', regex=True)
    
    # save file
    if output_file is not None:
        edges_df.to_csv(output_file, sep=delimiter, index=False, header=False)
    else:
        return edges_df