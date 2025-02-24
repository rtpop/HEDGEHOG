# import libraries
import pandas as pd

def process_edge_list(input_file, output_file, sep = ' '):
    """
    Processes an edge list file as output by netZooPy PANDA by removing the third column
    (the prior edge) and saves the result to a new file in a tab-delimited format suitable
    for loading as a weighted edgelist in NetworkX.

    Parameters:
    input_file (str): The path to the input edge list file.
    output_file (str): The path to the output file where the processed edge list will be saved.

    Returns:
    None
    """
    # Read the file into a pandas DataFrame
    df = pd.read_csv(input_file, sep = sep, header=None)
    
    # Ensure there are sufficient columns in the DataFrame
    if df.shape[1] < 4:
        raise ValueError("The input file does not have at least four columns.")
    
    # Remove the third column
    df = df.drop(columns=[2])
    
    # Write the modified DataFrame to a new file with tab separation
    df.to_csv(output_file, sep='\t', index=False, header=False)