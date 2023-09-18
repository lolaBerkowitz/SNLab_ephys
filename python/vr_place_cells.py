from ripple_heterogeneity.utils import (
    functions,
    loading,
)
from skimage import measure

from ripple_heterogeneity.place_cells import maps
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import nelpy as nel
import nelpy.plotting as npl



def load_data(basepath):
    # load position
    position_df = loading.load_animal_behavior(basepath)
    # restrict to non-nan times
    position_df_no_nan = position_df.query("not x.isnull() & not y.isnull()")
    position_df_no_nan['x'] = position_df_no_nan['x'] + np.abs(position_df_no_nan['x'].min())
    # create pos object
    pos = nel.AnalogSignalArray(
        data=position_df_no_nan["x"].values.T,
        timestamps=position_df_no_nan.timestamps.values,
    )
    # load in spike data 
    st,cm = loading.load_spikes(basepath,putativeCellType="Pyr|Int",brainRegion="CA1|Dentate")

    # load epochs
    epoch_df = loading.load_epoch(basepath)

    # put behavior epochs into nel.EpochArray
    beh_epochs = nel.EpochArray(np.array([epoch_df.startTime,epoch_df.stopTime]).T)

    return pos,st,cm,beh_epochs,epoch_df

def main(basepath):

    # load data
    pos,st,cm,beh_epochs,epoch_df = load_data(basepath)

    # find behavior_idx from epoch_df
    behavior_idx = int(epoch_df[epoch_df['name'].str.contains('vr')].index[0])

    # get outbound and inbound epochs
    (outbound_epochs, inbound_epochs) = functions.get_linear_track_lap_epochs(
        pos.abscissa_vals, pos.data[0], newLapThreshold=20)
    
    spatial_maps = maps.SpatialMap(
        pos[beh_epochs[behavior_idx]],
        st[beh_epochs[behavior_idx]],
        dir_epoch=outbound_epochs,
        dim=1,
        place_field_min_size=5,
        tuning_curve_sigma=3,
        place_field_sigma=.1,
        smooth_mode="wrap",
        speed_thres=10, 
        )
    spatial_maps.find_fields()
    
    # determine if fields are place fields
    pvals = spatial_maps.shuffle_spatial_information()

    # get IC for lag
    ic = spatial_maps.spatial_information()
    # get seletivity
    seletivity = spatial_maps.spatial_selectivity()

    # save to dataframe
    temp_df = pd.DataFrame(columns = ["IC","pval","seletivity","UID","basepath"])
    temp_df['IC'] = ic
    temp_df['pvals'] = pvals
    temp_df['seletivity'] = seletivity
    temp_df['UID'] = cm.UID
    temp_df['basepath'] = basepath
    temp_df['brainRegion'] = cm.brainRegion

    return temp_df
