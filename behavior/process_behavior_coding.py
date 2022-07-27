import pandas as pd
import numpy as np
import glob
import os
from logging import raiseExceptions
import nelpy as nel
from cmath import nan

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
        # print("Loading file: " + file)
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
        if idx.sum() == 0:
            print("Video not found: " + video )
            continue # skip if video not found
        for var in vars:
            df.loc[idx, var] = temp_df[var].values[0]

    return df

def check_start_stop(df):
    '''
    Returns index of rows where a stop does not follow a start
    '''
    df_check = df.copy()
    check_index = []
    # get index of rows that contain start 
    start_index = df_check.action.str.contains('start')
    start_index = start_index[start_index == True].index
    for i in start_index:
        if 'start' not in df.loc[0, 'action']:
            raiseExceptions('Start frame not found')
        if 'stop' not in df.loc[i+1, 'action']:
            check_index.append( i+1 )

    return check_index

def frames_to_seconds(df, fs = 60):
    '''
    Convert frames to seconds
    input: 
        df: dataframe with columns 'start', 'stop' and columns indicating trial start end frame
        fs: video frame rate

    output:
        df_new: converts all columns with frames to seconds
    '''
    df_new = pd.DataFrame(columns = df.columns)
    for file in df.file.unique():
        temp_df = df[df.file == file].copy()
        start_ts = temp_df.trial_start_1.values[0]
        temp_df.loc[temp_df.file == file, 'start'] = (temp_df.start.values - start_ts) / fs
        temp_df.loc[temp_df.file == file, 'stop'] = (temp_df.stop.values - start_ts) / fs
        temp_df[temp_df.keys()[temp_df.keys().str.contains('trial')]] = (temp_df[temp_df.keys()[temp_df.keys().str.contains('trial')]].values - start_ts) / fs
        
        df_new = pd.concat([df_new, temp_df], ignore_index=True)

    df_new.reset_index(drop = True,inplace=True)
    return df_new

def compute_dr(object_explore, moved_object_explore):
    """
    Compute the discrimination ratio DR from behavioral coding data

    inputs:
    object_explore: time in seconds when the still object was explored
    moved_object_explore: time in seconds when the moved object was explored

    outputs: 
    DR: discrimination ratio - the ratio of the time spent on the moved object
    """
    # compute DR
    if (moved_object_explore + object_explore) == 0:
        DR = nan
    else:
        DR = (moved_object_explore - object_explore) / (moved_object_explore + object_explore)
    return DR

def compute_explore_time(df, epoch_name = 'all_baseline_5min_test'):
    """
    Compute the time spent exploring each object within a given epoch

    inputs:
    df: dataframe with object_exploration (seconds), and indicator variables for when object exploration occurred. 
    epoch_name: name of indicator variable to use for object exploration

    outputs:
    df_out: dataframe with object_exploration (seconds), and indicator variables for when object exploration occurred.
   
    """
    # setup output dataframe    
    df_out = pd.DataFrame()
    d = []

    # compute object exploration for each object across all files within a given epoch
    for file in df.file.unique():
        temp_df = df[df.file == file].copy()
        
        # If there is no data for a given file, skip it
        if ~np.any(temp_df[epoch_name] == True):
            obj1 = 0
            obj2 = 0
            dr = nan
             # append data
            d.append({
                'file': file,
                'subid': temp_df.subid.iloc[0],
                'session_date': temp_df.session_date.iloc[0],
                'age': temp_df.age.iloc[0],
                'exposure': temp_df.exposure.iloc[0],
                'genotype': temp_df.iloc[0].genotype,
                'vidname': temp_df.vidname.unique()[0],
                'condition': temp_df.condition.unique()[0],
                'treatment': temp_df.treatment.unique()[0],
                'paradigm': temp_df.paradigm.unique()[0],
                'objects': temp_df.objects.unique()[0],
                'object1_explore': obj1,
                'object2_explore': obj2,
                'total_explore': obj1+obj2,
                'DR': dr       
                })
            continue
            
        temp_df = temp_df[temp_df[epoch_name] == True]
        obj1 = np.sum(temp_df[temp_df['object'] == '1'].stop - temp_df[temp_df['object'] == '1'].start)
        obj2 = np.sum(temp_df[temp_df['object'] == '2'].stop - temp_df[temp_df['object'] == '2'].start)

        if temp_df['moved_object'].str.contains('1').any():
            dr = compute_dr(obj2,obj1)
        else: 
            dr = compute_dr(obj1,obj2)
        print(temp_df.vidname.unique()[0])
        # append data
        d.append({
            'file': file,
            'subid': temp_df.subid.iloc[0],
            'session_date': temp_df.session_date.iloc[0],
            'age': temp_df.age.iloc[0],
            'exposure': temp_df.exposure.iloc[0],
            'genotype': temp_df.iloc[0].genotype,
            'vidname': temp_df.vidname.unique()[0],
            'condition': temp_df.condition.unique()[0],
            'treatment': temp_df.treatment.unique()[0],
            'paradigm': temp_df.paradigm.unique()[0],
            'objects': temp_df.objects.unique()[0],
            'object1_explore': obj1,
            'object2_explore': obj2,
            'total_explore': obj1+obj2,
            'DR': dr       
            })

    df_out = pd.concat([df_out, pd.DataFrame(d)], ignore_index=True)
    df_out.reset_index(drop = True,inplace=True)
    
    return df_out
