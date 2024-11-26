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

#%% calc deltas for profile
min_dep = 2000
max_dep = 4600

control_plot =False

var = 'potemperature'
data_list = [d for _, d in profile_data.groupby(['Station'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_by_bgfit(station, profile_background, var, min_dep, max_dep, fit='poly',param=(10),control_plot=True)
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

data = data[data['Station']=='028_01']
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


#%%  temp anomaly baker 2004

data = profile_data[profile_data['Station']=='022_01']
#data = data.sort_values(by='DEPTH')
# maybe take only first downcast
lim_above = (1000,2200)
lim_below = (3200,5000)
fit_order = 1

import numpy as np

data = sep_casts(data)[0]
data = data.sort_values(by='DEPTH')

data = data[(data['DEPTH']>lim_above[0]) & (data['DEPTH']<lim_below[1])]

import matplotlib.pyplot as plt
plt.scatter(data['potemperature'],data['Sigma3'])
#plt.scatter(data['potemperature'],data['DEPTH'])
plt.gca().invert_yaxis()

plt.colorbar()
plt.scatter(data_fit['Sigma3'],data_fit['potemperature'])
plt.scatter(below['Sigma3'],below['potemperature'])

above = data[(data['DEPTH']>lim_above[0]) & (data['DEPTH']<lim_above[1])]
below = data[(data['DEPTH']>lim_below[0]) & (data['DEPTH']<lim_below[1])]

data_fit = pd.concat([above,below])
data_fit = data_fit.dropna(subset=['potemperature','Sigma3'])

lin_theta = np.linspace(data_fit['potemperature'].min(),data_fit['potemperature'].max(),data_fit.index.max()-data_fit.index.min())


coef, res, rank, singular_val, rcond = np.polyfit(data_fit['potemperature'],data_fit['Sigma3'],fit_order,full=True)

plt.plot(data_fit['Sigma3'],data_fit['potemperature'])
fit_func = np.poly1d(coef)
fit = fit_func(lin_theta)
plt.plot(fit,lin_theta)

data['Bg_potemperature'] = (coef[0]*data['Sigma3']**1+coef[1]*data['Sigma3']**0)

plt.plot(data['Bg_potemperature'],data['DEPTH'])
plt.plot(data['potemperature'],data['DEPTH'])
plt.plot(data['potemperature']-data['Bg_potemperature'],data['DEPTH'])

#%%

min_dep = 2000
max_dep = 4600

control_plot =False

var = 'potemperature'
data_list = [d for _, d in profile_data.groupby(['Station'])]
for i, station in enumerate(data_list):
    data_list[i] = calc_delta_stafit(station, var, (1000,2500), (3500,4600), fit='uni',param=(5,0.005),control_plot=True)
profile_delta_potemperature = pd.concat(data_list)

