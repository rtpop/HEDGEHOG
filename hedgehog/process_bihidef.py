import pandas as pd
from datetime import datetime

def gmt_from_bihidef(df, output_file):
    """
    Create a GMT file from a Pandas DataFrame.

    Parameters:
    - df: DataFrame where:
        - The first column (index 0) is the gene set name.
        - The second column (index 1) is the size.
        - The third column (index 2) is a list or string of gene identifiers.
    - output_file: File path to save the GMT file.
    """
    with open(output_file, 'w') as file:
        for _, row in df.iterrows():
            gene_set_name = row[0]
            # Size (row[1]) can be used if needed for checks
            genes_column = row[2]
            
            # Check if genes are in a list or a space-separated string
            if isinstance(genes_column, list):
                genes = genes_column
            elif isinstance(genes_column, str):
                genes = genes_column.split(' ')
            else:
                raise ValueError("Genes must be a list or a space-separated string")

            gene_set_description = str(row[1]) + " genes"
            line_content = [gene_set_name, gene_set_description] + genes
            file.write("\t".join(line_content) + '\n')
            
def select_communities(filename, min_genes, max_genes, log_file=None):
    """
    Select communities with a number of genes within the specified range.

    Parameters:
    - filename: Path to a BIHIDEF .nodes file.
    - min_genes: Minimum number of genes in the community.
    - max_genes: Maximum number of genes in the community.
    - log_file: Path to a log file to save the number of selected communities.If None, no log file is saved.

    Returns:
    - DataFrame with the selected communities.
    """
    # Read the .nodes file into a DataFrame
    communities = pd.read_csv(filename, sep='\t', header=None)
    
    # drop the last column (I don't know what it is) \o/
    communities = communities.drop(columns=[3])
    
    # Assign column names to the DataFrame
    communities.columns = ['Community', 'Size', 'Genes']
    
    # Split the 'Genes' column into a list of genes
    communities['Genes'] = communities['Genes'].str.split(' ')
    
    # Filter communities based on the size constraints
    selected_communities = communities[(communities['Size'] >= min_genes) & (communities['Size'] <= max_genes)]
    
    if log_file:
        with open(log_file, 'w') as log:
            current_time = datetime.now().strftime("%d-%m-%Y %H:%M:%S")
            log.write(f"Time and date: {current_time}\n")
            log.write(f"Maximum number of genes per community: {max_genes}\n")
            log.write(f"Minimum number of genes per community: {min_genes}\n")
            log.write(f"Selected communities: {len(selected_communities)}\n")
            log.write(f"Total communities: {len(communities)}\n")
    
    # Return the filtered DataFrame
    return selected_communities