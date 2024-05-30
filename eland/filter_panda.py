def filter_edges(prior_file, panda_file, delimiter="\t"):
    """
    Filters edges from the panda network based on the intersection with the prior network,
    and retains only edges with weights greater than the minimum of the maximum weights
    per row or column.

    Parameters:
    prior_file (str): Path to the prior network edgelist file.
    panda_file (str): Path to the panda network edgelist file.
    delimiter (str): Delimiter used in the edgelist files.

    Returns:
    pd.DataFrame: DataFrame containing the filtered edgelist.
    """
    # Load prior and panda networks as directed graphs
    prior = nx.read_weighted_edgelist(prior_file, delimiter=delimiter, create_using=nx.DiGraph)
    panda = nx.read_weighted_edgelist(panda_file, delimiter=delimiter, create_using=nx.DiGraph)

    # Get the intersection of edges between prior and panda networks
    common_edges = set(prior.edges()).intersection(panda.edges())

    # Create a new directed graph to store the subset of the panda network
    weighted_prior = nx.DiGraph()

    # Add only the common edges to the panda subset with weights
    for u, v in common_edges:
        if panda.has_edge(u, v):
            weighted_prior.add_edge(u, v, weight=panda[u][v]['weight'])

    # Get the node list to map indices back to node names
    node_list = list(weighted_prior.nodes())
    node_index_map = {node: i for i, node in enumerate(node_list)}

    # Transform the weighted prior to adjacency matrix
    wt_prior_adj = nx.adjacency_matrix(weighted_prior, nodelist=node_list)

    # Find the maximum value for each row (TF)
    max_per_tf = wt_prior_adj.max(axis=1).toarray().flatten()

    # Find the maximum value for each column (gene)
    max_per_gene = wt_prior_adj.max(axis=0).toarray().flatten()

    # Determine the minimum edge selected by above
    min_edge = min(max_per_gene.min(), max_per_tf.min())

    # Select all edges higher than the minimum
    selected_edges_mask = wt_prior_adj > min_edge

    # Apply the mask to get the selected edges
    selected_edges = wt_prior_adj.multiply(selected_edges_mask)

    # Convert the sparse matrix back to an edgelist
    selected_edges_coo = selected_edges.tocoo()
    edges_list = [(node_list[selected_edges_coo.row[i]], node_list[selected_edges_coo.col[i]], selected_edges_coo.data[i])
                  for i in range(len(selected_edges_coo.data))]

    # Create a DataFrame for better visualization or saving to file
    edges_df = pd.DataFrame(edges_list, columns=['source', 'target', 'weight'])

    return edges_df