import pandas as pd
import networkx as nx
import numpy as np
import scipy.sparse as sp

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

def intersect_networks(prior, panda):
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

def filter_panda(prior_file, panda_file, delimiter='\t'):
    """
    Filters the PANDA network by processing prior and PANDA networks, finding their intersection,
    filtering the edges, and converting the filtered adjacency matrix to an edgelist DataFrame.
    
    Parameters:
    prior_file (str): Path to the prior network edge list file.
    panda_file (str): Path to the PANDA network edge list file.
    delimiter (str): Delimiter used in the edge list files.
    
    Returns:
    pd.DataFrame: A DataFrame representing the filtered edgelist with columns 'source', 'target', and 'weight'.
    """
    prior, panda = load_networks(prior_file, panda_file, delimiter)
    adjacency_matrix, node_list = intersect_networks(prior, panda)
    filtered_edges = filter_edges_by_max_values(adjacency_matrix)
    edges_df = convert_sparse_matrix_to_edgelist(filtered_edges, node_list)
    return edges_df