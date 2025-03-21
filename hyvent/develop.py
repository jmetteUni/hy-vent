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

#remove level part from 018 and 028 (no reindexing, maybe needed?) and stores it separetely -> Todo also MAPR
sta018 = profile_data['PS137_018_01']
sta028 = profile_data['PS137_028_01']
sta018_level = sta018.loc[(sta018['datetime'] > '2023-07-01 02:20:00') & (sta018['datetime'] < '2023-07-01 04:03:00')]
sta028_level = sta028.loc[(sta028['datetime'] > '2023-07-06 00:39:00') & (sta028['datetime'] < '2023-07-06 03:11:00')]
sta018 = pd.concat([sta018[sta018['datetime'] <= '2023-07-01 02:20:00'],sta018[sta018['datetime'] >= '2023-07-01 04:03:00']])
sta028 = pd.concat([sta028[sta028['datetime'] <= '2023-07-06 00:39:00'],sta028[sta028['datetime'] >= '2023-07-06 03:11:00']])

profile_data['PS137_018_01'] = sta018
profile_data['PS137_028_01'] = sta028
profile_data['PS137_018_le'] = sta018_level
profile_data['PS137_028_le'] = sta028_level

profile_data = keys_to_data(profile_data, 'cnv')
btl_data = keys_to_data(btl_data, 'btl')
mapr_data = keys_to_data(mapr_data, 'mapr')

profile_data = pd.concat(profile_data.values(),ignore_index=True)
btl_data = pd.concat(btl_data.values(),ignore_index=True)
mapr_data = pd.concat(mapr_data.values(),ignore_index=True)

#%% processing

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

sta018_level = profile_data[profile_data['Station']=='018_le']
sta028_level = profile_data[profile_data['Station']=='028_le']

profile_explo = profile_data[profile_data['Station']=='049_01']
mapr_explo = mapr_data[mapr_data['Station']=='049_01']
btl_explo = btl_data[mapr_data['Station']=='049_01']

profile_data = profile_data[profile_data['Station'].isin(aurora_stations)]
btl_data = btl_data[btl_data['Station'].isin(aurora_stations)]
mapr_data = mapr_data[mapr_data['Station'].isin(aurora_stations)]
profile_mapr = profile_mapr[profile_mapr['Station'].isin(aurora_stations)]

#%% calc deltas for profile
min_dep = 2000
max_dep = 4600

control_plot =False

var = 'potemperature'
data_list = [d for _, d in profile_data.groupby(['Station'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, profile_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=control_plot)
profile_delta_potemperature = pd.concat(data_list)

var = 'Sigma3'
data_list = [d for _, d in profile_data.groupby(['Station'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, profile_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=control_plot)
profile_delta_sigma3 = pd.concat(data_list)

var = 'PSAL'
data_list = [d for _, d in profile_data.groupby(['Station'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, profile_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=control_plot)
profile_delta_psal = pd.concat(data_list)

var = 'delta3He'
data_list = [d for _, d in btl_data.groupby(['Station'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_helium_delta(station, btl_background, var, 990, max_dep,control_plot=control_plot)
btl_delta_he = pd.concat(data_list)

#cal deltas for mapr
var = 'potemperature'
data_list = [d for _, d in mapr_data.groupby(['Station','SN'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, mapr_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=control_plot)
mapr_delta_potemperature = pd.concat(data_list)
mapr_delta_potemperature = mapr_delta_potemperature[(mapr_delta_potemperature['SN']=='74') | (mapr_delta_potemperature['SN']=='72')]

var = 'Sigma3'
data_list = [d for _, d in mapr_data.groupby(['Station','SN'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, mapr_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=control_plot)
mapr_delta_sigma3 = pd.concat(data_list)
mapr_delta_sigma3 = mapr_delta_sigma3[(mapr_delta_sigma3['SN']=='74') | (mapr_delta_sigma3['SN']=='72')]


# var = 'Neph_outl(volts)'
# data_list = [d for _, d in mapr_data.groupby(['Station','SN'])]
# for i, station in enumerate(data_list):
#     data_list[i] = calc_delta_by_bgfit(station, mapr_background, var, min_dep, max_dep, fit='poly', param=(10),control_plot=control_plot)
# mapr_delta_turb = pd.concat(data_list)

#Turb deltas by steation mean for a smaller range
var = 'Neph_smoo(volts)'
data_list = [d for _, d in mapr_data.groupby(['Station','SN'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_turb_delta(station, var, (2000,3000), (3500,4500))
mapr_delta_turb = pd.concat(data_list)


#%% histogramms & 2D plots & sections
# plot_hist_he(btl_data, lower_depth=3700, upper_depth=3000, bins=30)
# plot_hist_dorp(profile_data, lower_depth=3700, upper_depth=3000, bins=100, ranges=(-0.2,0.2))

# plot2D_all_stations(profile_data, 'dORP', vent_loc, 2000, 5000,bathy)

# plot2D_all_stations_btl(btl_data, 'delta3He', vent_loc, 2000, 5000,bathy)
# plot2D_station(profile_data, 'PS137_036_01', 'dORP', vent_loc, 2000, 5000,bathy)

# plot_section(profile_data, 'PS137_028_01', 'CTD_lat', 'DEPTH', 'dORP', 2400, 20)
# plot_section(profile_mapr[profile_mapr['Station']=='028_01'], '', 'CTD_lat', 'density', 'dORP', 2400, 20)


#%% crosssection

data = mapr_delta_turb
var = 'Delta_Neph_smoo(volts)'
min_dep = 2500
x = 'CTD_lat'

data = data[data['Station']=='041_01']
data = add_castno(data)
data = data[data['DEPTH']>min_dep]

# data_list = [d for _, d in data.groupby(['Cast'])]

# for i, cast in enumerate(data_list):
#     cast['Lat_mean'] = cast['CTD_lat'].mean()
#     data_list[i] = cast

# data = pd.concat(data_list)
data = data.dropna(subset=[x,'DEPTH',var])

fig, ax = plt.subplots()
fig.set_tight_layout(True)


contourf = ax.tricontourf(data[x], data['DEPTH'], data[var], levels=50)

ax.invert_yaxis()
plt.colorbar(contourf)  # r'$\Delta \theta$'
