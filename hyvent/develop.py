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

profile_data = add_castno(profile_data)
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
plot_hist_he(btl_data, lower_depth=3700, upper_depth=3000, bins=30)
plot_hist_dorp(profile_data, lower_depth=3700, upper_depth=3000, bins=100, ranges=(-0.2,0.2))

plot2D_all_stations(profile_data, 'dORP', vent_loc, 2000, 5000,bathy)

plot2D_all_stations_btl(btl_data, 'delta3He', vent_loc, 2000, 5000,bathy)
plot2D_station(profile_data, 'PS137_036_01', 'dORP', vent_loc, 2000, 5000,bathy)

plot_section(profile_data, 'PS137_028_01', 'CTD_lat', 'DEPTH', 'dORP', 2400, 20)
plot_section(profile_mapr[profile_mapr['Station']=='028_01'], '', 'CTD_lat', 'density', 'dORP', 2400, 20)

#%% separate casts and group data for 2D plotting
min_dep = 2000
max_dep = 4600

control_plot =False

var = 'potemperature'
data_list = [d for _, d in profile_data.groupby(['Station'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, profile_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=control_plot)
profile_delta_potemperature = pd.concat(data_list)

#%% enumerates all casts

profile_delta_potemperature = profile_delta_potemperature[(profile_delta_potemperature['Cast']<13) | (profile_delta_potemperature['Cast']>20)]      #remove casts during constant depth in station 028-01

#%% plot binned cast data

data = profile_delta_potemperature

from hyvent.misc import get_var
import numpy as np

label, color, cmap = get_var(var)

fig, ax = plt.subplots(figsize=(8,6))

#plot bathymetry
if bathy != 'None':
    contourlines = ax.contour(bathy[0],bathy[1],-bathy[2],levels=40, colors='black',linestyles='solid',linewidths=0.5,alpha=0.3)
    ax.clabel(contourlines, inline=True, fontsize=6, fmt='%d', colors = 'black')

#plot vent
if vent_loc != 'None':
    ax.scatter(vent_loc[0],vent_loc[1],color='red',marker='*',s=100,label='Aurora Vent Site')

# initialize an empty DataFrame to collect values
data_binned = pd.DataFrame(columns=['lon', 'lat', 'var_max'])

# loop through the list of dataframes
cast_list = [d for _, d in data.groupby(['Cast'])]
for i, cast in enumerate(cast_list):
    # get the values from the current dataframe
    lon = cast['CTD_lon'].mean()
    lat = cast['CTD_lat'].mean()
    if np.isnan(lat) == True:
        lat = cast['Ship_lat'].mean()
    if np.isnan(lon) == True:
        lon = cast['Ship_lon'].mean()
    var_max = cast['Delta_'+var].max()

    # add the values to the DataFrame
    data_binned.loc[i] = [lon, lat, var_max]

# cast_list = [d for _, d in data.groupby(['Cast'])]
# for cast in cast_list:

#         lat = cast['CTD_lat'].mean()
#         lon = cast['CTD_lon'].mean()
#         if np.isnan(lat) == True:
#             lat = cast['Ship_lat'].mean()
#         if np.isnan(lon) == True:
#             lon = cast['Ship_lon'].mean()

#         value = cast['Delta_'+var].max()
#         var_plot = plt.scatter(lon, lat, c=value)

#         print(lon,lat,value)

        #ax.text(lon.iloc[0]+0.007,lat.iloc[0]+0.0005,cast['cast'].iloc[0][1:-3],fontsize='small')

fig.colorbar(var_plot,label=label)






