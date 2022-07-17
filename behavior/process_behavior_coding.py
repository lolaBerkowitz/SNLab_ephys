import pandas as pd
import numpy as np
import glob
import os


def load_data(data_dir, pattern = 'ol_scoring.csv',var_names = ['file','frame','action']):
    """
    Load behavior_coding.py output csv from a directory
    """
    # find output using pattern 
    files = glob.glob(data_dir + os.path.sep +  "**\*" + pattern, recursive=True)
    if len(files) == 0: 
        raise Exception("No files found matching pattern: " + pattern)

    # load data
    df = pd.DataFrame()
    for file in files:
        print("Loading file: " + file)
        df = pd.concat([df ,pd.read_csv(file,usecols = var_names)], ignore_index=True)
    return df

