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