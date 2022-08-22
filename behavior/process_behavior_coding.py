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
        new_df = pd.read_csv(file,usecols = var_names)
        new_df["basepath"] = os.path.dirname(file)
        new_df["basename"] = os.path.basename(os.path.dirname(file))
        df = pd.concat([df ,new_df], ignore_index=True)
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
                'basepath': temp_df.basepath.iloc[0],
                'basename': temp_df.basename.iloc[0],
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
            'basepath': temp_df.basepath.iloc[0],
            'basename': temp_df.basename.iloc[0],
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

def reorganize_df(df):
    '''
    Reorganize dataframe to have columns for each action, and checks start and stop are balanced times for each action
    '''

    # make two new columns, one for object and one for action made
    df[['object','action']] = df.action.str.split('_', expand=True)
    df['object'] = df.object.str.replace('object','')

    # create dataframe where each row is a visit to object 1 or object 2
    df_visit = pd.DataFrame(columns = ['file','start','stop','object'])
    for file in df.file.unique():
        temp_df = df[df.file == file].copy()
        object_1_idx = temp_df.object == '1'
        object_2_idx = temp_df.object == '2'

        if np.sum(object_1_idx) % 2 == 1:
            raise Exception('Object 1 has an odd number of starts and stops')
        if np.sum(object_2_idx) % 2 == 1:
            raise Exception('Object 2 has an odd number of starts and stops')
        
        d = []
        for idx, i in enumerate(range(0,len(temp_df),2)):
            d.append(
                {'file': temp_df.iloc[i].file,
                'basepath':temp_df.iloc[i].basepath,
                'basename':temp_df.iloc[i].basename,
                'start': temp_df.iloc[i].frame,
                'stop': temp_df.iloc[i+1].frame,
                'object': temp_df.iloc[i].object
                }
            )
        
        df_visit = pd.concat([df_visit, pd.DataFrame(d)], ignore_index=True)

    df_visit.reset_index(drop = True,inplace=True)
    
    return df_visit

def add_indicators(df):
    """
    Adds logical variables for parsing dataframe
    first_5min = first 5 minutes of visit
    first_3min = last 5 minutes of visit    
    """

    df['object_exploration_time'] = df.stop - df.start
    df['first_5min'] = df.stop.apply(lambda x: x < 300)
    df['first_3min'] = df.stop.apply(lambda x: x < 180)
    df['whole_session'] = True
    df['all_baseline_3min_test'] = True
    df.loc[df['condition'] == 'test','all_baseline_3min_test'] = df.loc[df['condition'] == 'test','first_3min']

    df['all_baseline_5min_test'] = True
    df.loc[df['condition'] == 'test','all_baseline_5min_test'] = df.loc[df['condition'] == 'test','first_5min']

    df['trial_n'] = np.zeros(len(df))
    df.loc[df['stop'] < df['trial_stop_1'],'trial_n'] = 1
    df.loc[(df['stop'] > df['trial_stop_1']) & (df['stop'] < df['trial_stop_2']),'trial_n'] = 2
    df.loc[(df['stop'] > df['trial_stop_2']) & (df['stop'] < df['trial_stop_3']),'trial_n'] = 3
    df["in_trial"] = df.trial_n.apply(lambda x: x > 0)
    return df