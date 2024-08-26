#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 14:57:38 2024

@author: jmette@uni.bremen.de
"""

from lat_lon_parser import parse

from hyvent.io import (read_from_csv,
                       read_gebco)

from hyvent.plotting_maps import *
from hyvent.plotting import *

#%% set paths
cnv_path = '/home/jonathan/Dokumente/SHK Maren/PS137/CTD_data_processing/csv_files_noprocessing/'
btl_path = '/home/jonathan/Dokumente/SHK Maren/PS137/CTD_data_processing/csv_files_noprocessing/'
gebcopath = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84.nc'     # path for bathymetry files
gebcopath_hillshade = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84_hillshade.nc'

vent_loc = (parse("6° 15.32'W"),parse("82° 53.83'N"))

#%% import data
btl_data = read_from_csv(btl_path, 'btl')
profile_data = read_from_csv(cnv_path, 'cnv')
lats, lons, elev = read_gebco(gebcopath)
bathy = (lons,lats,elev)
#%% Aurora overview maps
plot_map(profile_data, btl_data, lats, lons, elev, tracer_type='None',path_save='/home/jonathan/Dokumente/SHK Maren/PS137/fuer_Chris_Gunter/Paper_Chris/AuroraMap.png')
plot_map(profile_data, btl_data, lats, lons, elev, tracer_type='dna_poc',path_save='/home/jonathan/Dokumente/SHK Maren/PS137/fuer_Chris_Gunter/Paper_Chris/AuroraMap_dnapoc.png')
plot_map(profile_data, btl_data, lats, lons, elev, tracer_type='he',path_save='/home/jonathan/Dokumente/SHK Maren/PS137/fuer_Chris_Gunter/Paper_Chris/AuroraMap_he.png')

#%% He samples

plot2D_all_stations_btl(btl_data, 'delta3He', vent_loc, 0, 6000,bathy)
plot2D_all_stations_btl(btl_data, 'delta3He', vent_loc, 2800, 6000,bathy)

