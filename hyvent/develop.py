#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 10:59:34 2024

@author: jonathan
"""
from lat_lon_parser import parse

from hyvent.io import (read_from_csv,
                       read_gebco)

from hyvent.plotting_maps import *
from hyvent.plotting import *
from hyvent.misc import *

#%%
cnv_path = '/home/jonathan/Dokumente/SHK Maren/PS137/CTD_data_processing/csv_files_noprocessing/'
btl_path = '/home/jonathan/Dokumente/SHK Maren/PS137/CTD_data_processing/csv_files_noprocessing/'
mapr_path = '/home/jonathan/Dokumente/SHK Maren/PS137/MAPR/csv_files/'
gebcopath = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84.nc'     # path for bathymetry files
gebcopath_hillshade = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84_hillshade.nc'

vent_loc = (parse("6° 15.32'W"),parse("82° 53.83'N"),3888)

#%%
btl_data = read_from_csv(btl_path, 'btl')
profile_data = read_from_csv(cnv_path, 'cnv')
mapr_data = read_from_csv(mapr_path, 'mapr')
lats, lons, elev = read_gebco(gebcopath)
bathy = (lons,lats,elev)

profile_data = keys_to_data(profile_data, 'cnv')
btl_data = keys_to_data(btl_data, 'btl')
mapr_data = keys_to_data(mapr_data, 'mapr')
#%%
plot_largemap_lambert(vent_loc,'Aurora Vent Site')
plot_largemap_stereo(vent_loc,'Aurora Vent Site')

#%% histogramms
plot_hist_he(btl_data, lower_depth=3700, upper_depth=3000, bins=30)
plot_hist_dorp(profile_data, lower_depth=3700, upper_depth=3000, bins=100, ranges=(-0.2,0.2))

#%%

plot2D_station(profile_data, 'PS137_036_01', 'dORP', vent_loc, 2000, 5000,bathy)
#%%
plot2D_all_stations(profile_data, 'dORP', vent_loc, 2000, 5000,bathy)

plot2D_all_stations_btl(btl_data, 'delta3He', vent_loc, 2000, 5000,bathy)
#%% merge profile and mapr data

import pandas as pd

all_profile = pd.concat(profile_data.values(),ignore_index=True)
all_mapr = pd.concat(mapr_data.values(),ignore_index=True)
all_mapr = all_mapr.rename(columns={'Press(dB)':'PRES','Temp(°C)':'TEMP','Depth(m)':'DEPTH'})
mapr_for_merge = all_mapr[['datetime','PRES','TEMP','DEPTH','Neph(volts)','Dship_lon','Dship_lat','CTD_lon','CTD_lat','dORP','Cruise','Station','SN','Operation']]
mapr_for_merge = mapr_for_merge[mapr_for_merge['Operation']=='CTD']
profile_mapr = pd.concat((all_profile,mapr_for_merge),axis=0)
profile_mapr.sort_values(by='datetime',ascending=True,inplace=True)

#profile_mapr = profile_mapr.groupby('Station')[['dORP']].apply(lambda g: g.values.tolist()).to_dict()

#%%
import matplotlib.pyplot as plt

plot_section(profile_data, 'PS137_028_01', 'CTD_lat', 'DEPTH', 'dORP', 2400, 20)
plot_section(profile_mapr[profile_mapr['Station']=='028_01'], '', 'CTD_lat', 'density', 'dORP', 2400, 20)


#%%

sta036 = profile_mapr[profile_mapr['Station']=='036_01']
ctd036 = sta036[sta036['SN']=='SBE9']
sn74_036 = sta036[sta036['SN']=='74']
sn72_036 = sta036[sta036['SN']=='72']
sn73_036 = sta036[sta036['SN']=='73']
#%%
plt.plot(ctd036['CTD_lat'],ctd036['DEPTH'])
plt.plot(sn74_036['CTD_lat'],sn74_036['DEPTH'])
plt.plot(sn73_036['CTD_lat'],sn73_036['DEPTH'])
plt.plot(sn72_036['CTD_lat'],sn72_036['DEPTH'])

#%%
from hyvent.misc import get_var
get_var('dORPs')


