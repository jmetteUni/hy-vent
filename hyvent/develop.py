#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:59:34 2024

@author: jonathan
"""
import os
working_dir = '/home/jonathan/Dokumente/Masterarbeit/python/hy-vent/'
os.chdir(working_dir)

from lat_lon_parser import parse

from hyvent.io import (read_from_csv,
                       read_gebco)
from hyvent.plotting_maps import *
from hyvent.plotting import *
from hyvent.misc import keys_to_data
from hyvent.processing import calc_mean_profile, derive_mapr, derive_CTD

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
vent_loc = (parse("6° 15.32'W"),parse("82° 53.83'N"))
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

profile_data = derive_CTD(profile_data)
mapr_data = mapr_data.rename(columns={'Press(dB)':'PRES','Temp(°C)':'TEMP','Depth_corr(m)':'DEPTH'})
mapr_data = mapr_data[['datetime','PRES','TEMP','DEPTH','Neph_outl(volts)','Dship_lon','Dship_lat','CTD_lon','CTD_lat','dORP','Cruise','Station','SN','Operation']]
mapr_data = derive_mapr(mapr_data, profile_data, aurora_stations)

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

#%% histogramms & 2D plots & sections
plot_hist_he(btl_data, lower_depth=3700, upper_depth=3000, bins=30)
plot_hist_dorp(profile_data, lower_depth=3700, upper_depth=3000, bins=100, ranges=(-0.2,0.2))

plot2D_all_stations(profile_data, 'dORP', vent_loc, 2000, 5000,bathy)

plot2D_all_stations_btl(btl_data, 'delta3He', vent_loc, 2000, 5000,bathy)
plot2D_station(profile_data, 'PS137_036_01', 'dORP', vent_loc, 2000, 5000,bathy)

plot_section(profile_data, 'PS137_028_01', 'CTD_lat', 'DEPTH', 'dORP', 2400, 20)
plot_section(profile_mapr[profile_mapr['Station']=='028_01'], '', 'CTD_lat', 'density', 'dORP', 2400, 20)

# %% substract background

def get_bg_fit(bg, var, min_dep, max_dep, fit_order=10):
    """
    This funtions fits a polynomial to a variable of a profile (usually to create a smooth background profile), with depth as x values for the fit.

    Parameters
    ----------
    bg : pandas dataframe
        Dataset which should be used as background.
    var : string
        Variable, where the deviation should be calculated.
    min_dep : int
        Minimum depth for the fit's depth range.
    max_dep : int
        Maximum depth for the fit's depth range.
    fit_order : TYPE, optional
        Order of the polynomial fit. The default is 10.

    Returns
    -------
    df_fit : pandas dataframe
        Datframe containing the depth and the fit result.

    """
    import numpy as np
    import pandas as pd
    from hyvent.processing import sep_casts

    #get first downcast of background and subset by depth range
    bg = sep_casts(bg, window_size=5000)
    if len(bg) > 2:
        print('Warning: Your background station was separated into more than two up or down casts. Check the casts returned by sep_casts and the used window size.')
    bg = bg[0]
    bg = bg[(bg['DEPTH']>=min_dep) & (bg['DEPTH']<=max_dep)]

    coef = np.polyfit(bg['DEPTH'],bg[var],fit_order)
    fit = np.poly1d(coef)
    fit = fit(bg['DEPTH'])
    df = {'DEPTH':bg['DEPTH'],'Fit_'+var:fit}
    df_fit = pd.DataFrame(df)

    return df_fit


def calc_delta_by_fit(data, bg, var, min_dep, max_dep): #casts: data, order: Order of PolyFit

    bg_fit = get_bg_fit(bg, var, min_dep, max_dep)
    pd.merge_asof(data[pd.notna(data['DEPTH'])], bg_fit, on='DEPTH', direction='nearest')
    data['Delta_'+var] = data[var] - data['Fit_'+var]

    return data

delta_fit036 = calc_delta_by_fit(sta036, profile_background, 'potemperature', 2000, 5000)


#%% test bg fit

plt.plot(bg['potemperature'],bg['DEPTH'])
plt.plot(Bg_fit['Bg_fit_potemperature'],Bg_fit['DEPTH'])


#%% plot for anomaly

sta036 = profile_data[profile_data['Station']=='036_01']
#delta_dens036 = subtract_bg_by_press(sta036, profile_background, 'potemperature', 1800, 5000)
delta_presslayers036 = calcTempAnomalyPress(sta036, profile_background, 100)

#%%
plt.figure()
plt.plot(sta036['potemperature'],sta036['DEPTH'],label='theta 036')
plt.plot(profile_background['potemperature'],profile_background['DEPTH'],label='theta 018')

plt.gca().invert_yaxis()
plt.legend()

#%%
plt.figure()
plt.plot(delta036['TempAnomaly'],delta036['DEPTH'],label='delta theta')
plt.plot(delta_presslayers036['TempAnomaly'],delta_presslayers036['DEPTH'],label='delta theta press')
#plt.plot(delta_dens036['Delta_potemperature'],delta_dens036['DEPTH'],label='delta theta dens')


plt.gca().invert_yaxis()
plt.legend()






