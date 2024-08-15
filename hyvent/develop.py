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

#%%
cnv_path = '/home/jonathan/Dokumente/SHK Maren/PS137/CTD_data_processing/csv_files_noprocessing/'
btl_path = '/home/jonathan/Dokumente/SHK Maren/PS137/CTD_data_processing/csv_files_noprocessing/'
gebcopath = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84.nc'     # path for bathymetry files
gebcopath_hillshade = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84_hillshade.nc'

vent_loc = (parse("6° 15.32'W"),parse("82° 53.83'N"),3888)

#%%
btl_data = read_from_csv(btl_path, 'btl')
profile_data = read_from_csv(cnv_path, 'cnv')
lats, lons, elev = read_gebco(gebcopath)
bathy = (lons,lats,elev)
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

plot2D_all_stations_he(btl_data, 'delta3He', vent_loc, 2000, 5000,bathy)
#%%
