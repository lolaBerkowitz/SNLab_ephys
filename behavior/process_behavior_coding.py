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

def add_metadata(df, df_meta):
    """
    Add metadata to dataframe
    """
    vars = df_meta.keys()

    for video in df_meta.vidname.unique():
        temp_df = df_meta[df_meta.vidname == video].copy()

        idx = df["file"].str.contains(temp_df.vidname.values[0])
        for var in vars:
            df.loc[idx, var] = temp_df[var].values[0]

    return df