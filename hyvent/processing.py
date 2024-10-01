#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 12:52:37 2024

@author: jmette@uni-bremen.de
"""

def corr_mapr_depth(data,lat):
    """

    Calculates a corrected depth based on the 200 first pressure values larger then 5 dB which are assumed as measured before deploying. The new, corrected depth is calculated with Gibbs Seawater Toolbox (gsw.conversions.z_from_p()) with the pressure corrected by the mean of the pre-deployment values and written to the dataframe as a new column.


    Parameters
    ----------
    data : pandas dataframe
        Data of one station of one MAPR as a pandas dataframe.
    lat : interger or float
        Latitude where the measurement was taken to calculate the depth from pressure.

    Returns
    -------
    depth_corr : pandas series
        Column containining the corrected depth.

    """

    import gsw

    p_deck = data[:200]     #take first 200 measurements which are negative
    p_deck = p_deck[p_deck['Press(dB)']<5]      #exclude values which are already in water measured
    if len(p_deck) < 50:        #warn if only a few measurements are pre-deplyoment
        print('Warning: Less then 50 values used for deck reference pressure!')
    p_deck_mean = p_deck['Press(dB)'].mean()
    #print('Mean deck pressure: '+str(p_deck_mean)+' dB')
    data['Depth_corr(m)'] = -data.apply(lambda x: gsw.conversions.z_from_p(x['Press(dB)']-p_deck_mean,lat), axis=1)         #calculates the corrected depth

    depth_corr = data['Depth_corr(m)']    #select only the corrected depth column

    return depth_corr

def calc_mean_profile(data, var, station_list):
    """
    This function calculates a mean profile of the data of mutliple stations for a given variable. Values a binned to 1m depth intervals between the minimum and maximum depth present in the data.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data of one or multiple CTD stations with variables as columns. One column has to be the station designation.
    var : string
        Variable to calculate the mean of. Must be a column key in data.
    station_list : list of strings
        List of station designations, over which the mean should be calculated.

    Returns
    -------
    mean_profile : pandas dataframe
        Dataframe containing two columns which are the binned depth and the associated mean of the variable.

    """
    import numpy as np
    import pandas as pd

    min_depth = data['DEPTH'].min()
    max_depth = data['DEPTH'].max()
    depth_vec = np.arange(min_depth,max_depth,1)        #create depth vector with step 1m, where variabel is binned to
    all_profiles = pd.DataFrame(depth_vec,columns=['DEPTH'])
    for station in station_list:
        profile = data[data['Station']==station][['DEPTH',var]]
        new_key = var+station
        profile.rename(columns={var:new_key},inplace=True)
        profile.sort_values(by='DEPTH',ascending=True,inplace=True)
        profile.dropna(inplace=True)
        all_profiles = pd.merge_asof(all_profiles,profile,on='DEPTH',direction='nearest')       # for every station merge variable by nearest depth value -> binning
    del all_profiles['DEPTH']
    mean_name = str(var)+'_mean'
    all_profiles[mean_name] = all_profiles.mean(axis=1)    #calculate mean variable value
    all_profiles['DEPTH'] = depth_vec
    mean_profile = all_profiles[['DEPTH',mean_name]]

    return mean_profile