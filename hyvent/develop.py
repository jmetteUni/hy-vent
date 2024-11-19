#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:59:34 2024

@author: jmette@uni-bremen.de
"""
import os
working_dir = '/home/jonathan/Dokumente/Masterarbeit/python/hy-vent/hyvent/'
os.chdir(working_dir)

from lat_lon_parser import parse

from hyvent.io import (read_from_csv,
                       read_gebco)
from hyvent.plotting_maps import *
from hyvent.plotting import *
from hyvent.misc import keys_to_data, add_castno
from hyvent.processing import *

import pandas as pd
import gsw
import matplotlib.pyplot as plt

#%% set paths
cnv_path = '/home/jonathan/Dokumente/SHK Maren/PS137/CTD_data_processing/csv_files_noprocessing/'
btl_path = '/home/jonathan/Dokumente/SHK Maren/PS137/CTD_data_processing/csv_files_noprocessing/'
mapr_path = '/home/jonathan/Dokumente/SHK Maren/PS137/MAPR/csv_files/'
gebcopath = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84.nc'     # path for bathymetry files
gebcopath_hillshade = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84_hillshade.nc'
fig_path = '/home/jonathan/Dokumente/Masterarbeit/Thesis_Latex_Project/figures/4-results/obs/'

lonW = -6.45
lonE = -6.175
latS = 82.875
latN = 82.92
vent_loc = (parse("6° 15.32'W"),parse("82° 53.83'N"),3900)
aurora_stations = ['022_01', '026_01', '028_01', '033_01', '036_01', '041_01', '054_01', '055_01']

#%% import data
btl_data = read_from_csv(btl_path, 'btl')
profile_data = read_from_csv(cnv_path, 'cnv')
mapr_data = read_from_csv(mapr_path, 'mapr')
bathy = read_gebco(gebcopath,lonE,lonW,latS,latN)

profile_data = keys_to_data(profile_data, 'cnv')
btl_data = keys_to_data(btl_data, 'btl')
mapr_data = keys_to_data(mapr_data, 'mapr')

profile_data = pd.concat(profile_data.values(),ignore_index=True)
btl_data = pd.concat(btl_data.values(),ignore_index=True)
mapr_data = pd.concat(mapr_data.values(),ignore_index=True)

#%% processing

#try to remove the density spike in 028-01
#import numpy as np
#profile_data[profile_data['Station']=='028_01'] = profile_data[profile_data['Station']=='028_01'].mask((profile_data[profile_data['Station']=='028_01']['DEPTH']>3180) & (profile_data[profile_data['Station']=='028_01']['DEPTH']<3220) & (profile_data[profile_data['Station']=='028_01']['density']>41.98396) & (profile_data[profile_data['Station']=='028_01']['density'] <41.9846),other=np.nan)

profile_data = derive_CTD(profile_data)
mapr_data = mapr_data.rename(columns={'Press(dB)':'PRES','Temp(°C)':'TEMP','Depth_corr(m)':'DEPTH'})
mapr_data = mapr_data[['datetime','PRES','TEMP','DEPTH','Neph(volts)','Dship_lon','Dship_lat','CTD_lon','CTD_lat','dORP','Cruise','Station','SN','Operation']]
mapr_data = derive_mapr(mapr_data, profile_data, aurora_stations)

#%% unify profile and mapr

mapr_formerge = mapr_data[mapr_data['Operation']=='CTD']
profile_mapr = pd.concat((profile_data,mapr_formerge),axis=0)
profile_mapr['datetime'] = pd.to_datetime(profile_mapr['datetime'])

#%% sort stations
profile_background = profile_data[profile_data['Station']=='018_01']
btl_background = btl_data[btl_data['Station']=='018_01']
mapr_background = mapr_data[mapr_data['Station']=='018_01']

profile_explo = profile_data[profile_data['Station']=='049_01']
mapr_explo = mapr_data[mapr_data['Station']=='049_01']
btl_explo = btl_data[mapr_data['Station']=='049_01']

profile_data = profile_data[profile_data['Station'].isin(aurora_stations)]
btl_data = btl_data[btl_data['Station'].isin(aurora_stations)]
mapr_data = mapr_data[mapr_data['Station'].isin(aurora_stations)]
profile_mapr = profile_mapr[profile_mapr['Station'].isin(aurora_stations)]


#%% histogramms & 2D plots & sections
# plot_hist_he(btl_data, lower_depth=3700, upper_depth=3000, bins=30)
# plot_hist_dorp(profile_data, lower_depth=3700, upper_depth=3000, bins=100, ranges=(-0.2,0.2))

# plot2D_all_stations(profile_data, 'dORP', vent_loc, 2000, 5000,bathy)

# plot2D_all_stations_btl(btl_data, 'delta3He', vent_loc, 2000, 5000,bathy)
# plot2D_station(profile_data, 'PS137_036_01', 'dORP', vent_loc, 2000, 5000,bathy)

# plot_section(profile_data, 'PS137_028_01', 'CTD_lat', 'DEPTH', 'dORP', 2400, 20)
# plot_section(profile_mapr[profile_mapr['Station']=='028_01'], '', 'CTD_lat', 'density', 'dORP', 2400, 20)

#%% separate casts and group data for 2D plotting
min_dep = 2000
max_dep = 4600

control_plot =False

var = 'potemperature'
data_list = [d for _, d in profile_data.groupby(['Station','SN'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, profile_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=control_plot)
profile_delta_potemperature = pd.concat(data_list)


var = 'potemperature'
data_list = [d for _, d in mapr_data.groupby(['Station','SN'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, mapr_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=control_plot)
mapr_delta_potemperature = pd.concat(data_list)
mapr_delta_potemperature = mapr_delta_potemperature[(mapr_delta_potemperature['SN']=='74') | (mapr_delta_potemperature['SN']=='72')]

var = 'Neph_smoo(volts)'
data_list = [d for _, d in mapr_data.groupby(['Station','SN'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_turb_delta(station, var, (2000,3000), (3500,4500))
mapr_delta_turb = pd.concat(data_list)


#%%

data = mapr_delta_turb
window_size = 700
var = 'potemperature'
#data = data[(data['Cast']<13) | (data['Cast']>20)]      #remove casts during constant depth in station 028-01


def plot_2Dpercast(data, var, min_dep, max_dep, window_size=1000, bathy='None', vent_loc='None'):
    """
    This function plots the data of one or mutliple CTD stations on a 2D Map, with one point per up- or downcast and the value being the maximum value of the variable during this cast. It is plotted on the mean position calculated from the position data of the cast. The data can be limited by a depth range and optionaly bathymetry and a vent location can be plooted.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data of the CTD station or multiple stations with variable as columns.
    var : string
        Variable to plot. Must be a key for a column in data and is only plotted for values with coordinates in the columns "CTD_lon"/"CTD_lat" or "Dship_lon"/"Dship_lat".
    min_dep : int
        Minimum depth below which datapoints should be plotted.
    max_dep :  int
        Maximum depth above which datapoints should be plotted.
    window_size : int
        Size of the rolling window, to find local extrema within. The optimal size is dependent on the number of measurements per up or down cast. The default is 1000.
    vent_loc : tuple, optional
        Longitude and latitude of a point of interest, plotted as a red star. Default is 'None"'
    bathy : tuple, optional
        Tuple of longitude, latitude and elevation of bathymetry data. The default is 'None'.

    Returns
    -------
    None.

    """
    from hyvent.misc import get_var
    from hyvent.misc import add_castno
    import numpy as np

    label, color, cmap = get_var(var)

    data = add_castno(data, window_size)
    data = data[(data['DEPTH']<min_dep) & data['DEPTH']<max_dep]

    fig, ax = plt.subplots(figsize=(8,6))

    #plot bathymetry
    if bathy != 'None':
        contourlines = ax.contour(bathy[0],bathy[1],-bathy[2],levels=40, colors='black',linestyles='solid',linewidths=0.5,alpha=0.3)
        ax.clabel(contourlines, inline=True, fontsize=6, fmt='%d', colors = 'black')

    #plot vent
    if vent_loc != 'None':
        ax.scatter(vent_loc[0],vent_loc[1],color='red',marker='*',s=100,label='Aurora Vent Site')

    # initialize an empty DataFrame to collect values
    data_binned = pd.DataFrame(columns=['Station','lon', 'lat', 'var_max'])

    # loop through the list of dataframes
    cast_list = [d for _, d in data.groupby(['Cast','Station'])]
    for i, cast in enumerate(cast_list):
        # get the values from the current dataframe
        station = cast['Station'].iloc[0]
        lon = cast['CTD_lon'].mean()
        lat = cast['CTD_lat'].mean()
        if np.isnan(lat) == True:
            lat = cast['Dship_lat'].mean()
        if np.isnan(lon) == True:
            lon = cast['Dship_lon'].mean()
        var_max = cast['Delta_'+var].max()

        # add the values to the DataFrame
        data_binned.loc[i] = [station, lon, lat, var_max]

    var_plot = plt.scatter(data_binned['lon'], data_binned['lat'], c=data_binned['var_max'], cmap = get_var(var)[2])

    fig.colorbar(var_plot,label=label)
    #change colorbar ticks settings at some point?

    plt.show()

plot_2Dpercast(mapr_delta_turb, 'Neph_smoo(volts)', 2000, 4600, window_size=700, bathy=bathy, vent_loc=vent_loc)


#%%

