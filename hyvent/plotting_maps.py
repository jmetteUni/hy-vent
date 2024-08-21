#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 13:27:30 2024

@author: jmette@uni-bremen.de
"""

def plot_map(profile_data, btl_data, lats, lons, elev, tracer_type='None', path_save='None'):
    """
    Plots CTD tow-yo tracks on to a 2D map together with bathymety and optional samples of tracer. Tracks are a combination of acoustic tracker position of the CTD if available (solid line), otherwise ship position (dashed line).

    Some sections contain code, which makes this function specific for PS137 data for now.

    Parameters
    ----------
    profile_data : dictionary
        Dictionary of CTD profile data where the keys are unique station identifiers and
        values are data in pandas dataframes.
    btl_data : dictionary
        Dictionary of SeaBird bottle data where the keys are unique station identifiers and
        values are data in pandas dataframes.
    lats : numpy ndarrays
        Latitude coordinates of bathymetry.
    lons : numpy ndarrays
        Longitude coordinates of bathymetry.
    elevs : numpy ndarrays
        Elevation of bathymetry.
    tracer_type : string, optional
        Type of sample to plot on the map. Can be "dna" for DNA samples, "poc" for poc samples, "dna_poc" for both and "he" for helium samples The default is 'None'.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """

    import cartopy.crs as ccrs
    from lat_lon_parser import parse
    import numpy as np
    import matplotlib.pyplot as plt

    #use Arial from custom path
    import matplotlib.font_manager

    font_dir = ['/home/jonathan/Downloads/arial_font/']
    for font in matplotlib.font_manager.findSystemFonts(font_dir):
        matplotlib.font_manager.fontManager.addfont(font)

    plt.rcParams['font.family'] = 'Arial'#'sans-serif'
    #plt.rcParams.update({'font.size': 10})

    #%% subset profile data for plots

    profile_plot = dict()

    for key in profile_data:
        profile_plot[key] = profile_data[key]
        #profile_plot[key] = profile_plot[key][profile_data[key]['DEPTH']>200]
        #profile_plot[key].reset_index(inplace=True)

    del profile_plot['PS137_018_01']        #background station, south Lena Trough
    del profile_plot['PS137_025_01']        #not deep, not interested in
    del profile_plot['PS137_040_01']        #aborted at 1600m because of ice drift
    del profile_plot['PS137_049_01']        #ATWAICE calibration cast
    del profile_plot['PS137_052_01']        #calibration cast, only down to 1000m
    del profile_plot['PS137_053_01']        #aborted at 170m due to ice situation
    del profile_plot['PS137_058_01']        #Lucky B
    del profile_plot['PS137_061_01']        #Lucky B

    slice028 = profile_plot['PS137_028_01'].iloc[6925:6975, :]
    slice028.loc[:, 'CTD_lon'] = np.nan
    profile_plot['PS137_028_01'].update(slice028)

    start_lon = dict()
    start_lat = dict()
    end_lon = dict()
    end_lat = dict()

    start_lon['PS137_022_01'] = parse('6W,15.10W')   #no towyo
    start_lat['PS137_022_01'] = parse('82N,53.89')
    end_lon['PS137_022_01'] = parse('6W,12.955W')
    end_lat['PS137_022_01'] = parse('82N,53.703')
    #start_lon['PS137_025_01'] = parse('6W,15.136W')  #no towo
    #start_lat['PS137_025_01'] = parse('82N,54.17')

    start_lon['PS137_026_01'] = parse('6W,13.097W')  #no towo
    start_lat['PS137_026_01'] = parse('82N,53.649')
    end_lon['PS137_026_01'] = parse('6W,11.152W')
    end_lat['PS137_026_01'] = parse('82N,53.392N')

    start_lon['PS137_028_01'] = parse('6W,24.96W')   #towyo with posidonia
    start_lat['PS137_028_01'] = parse('82N,54.59N')
    end_lon['PS137_028_01'] = parse('6W,21.784W')
    end_lat['PS137_028_01'] = parse('82N,54.08N')

    start_lon['PS137_033_01'] = parse('6W,24.036W')   #towyo with posidonia
    start_lat['PS137_033_01'] = parse('82N,54.514N')
    end_lon['PS137_033_01'] = parse('6W,16.63W')
    end_lat['PS137_033_01'] = parse('82N,52.69N')

    start_lon['PS137_036_01'] = parse('6W,19.02W')    #towyo without posidonia
    start_lat['PS137_036_01'] = parse('82N,54.83')
    end_lon['PS137_036_01'] = parse('6W,15.8W')
    end_lat['PS137_036_01'] = parse('82N,53.18N')

    start_lon['PS137_041_01'] = parse('6W,14.80W')   #towyo with posidonia
    start_lat['PS137_041_01'] = parse('82N,54.03N')
    end_lon['PS137_041_01'] = parse('6W,14.25W')
    end_lat['PS137_041_01'] = parse('82N,53.59N')

    start_lon['PS137_054_01'] = parse('6W,18.36W')    #no towyo
    start_lat['PS137_054_01'] = parse('82N,53.84')
    end_lon['PS137_054_01'] = parse('6W,11.935W')
    end_lat['PS137_054_01'] = parse('82N,53.149N')

    start_lon['PS137_055_01'] = parse('6W,16.121W')    #no towyo
    start_lat['PS137_055_01'] = parse('82N,54.062')
    end_lon['PS137_055_01'] = parse('6W,13.519W')
    end_lat['PS137_055_01'] = parse('82N,53.728N')

    #%%  DNA and POC positions
    poc = {'PS137_026_01' :13,
           'PS137_028_01' :10,
           'PS137_036_01' :7,
           'PS137_041_01' :9,
           'PS137_054_01' :5}

    poc_lon = list()
    poc_lat = list()

    for key in poc:
        if key == 'PS137_022_01' or key == 'PS137_026_01':
            poc_lon.append(btl_data[key]['Dship_lon'][btl_data[key]['Bottle']==poc[key]])
            poc_lat.append(btl_data[key]['Dship_lat'][btl_data[key]['Bottle']==poc[key]])
        else:
            poc_lon.append(btl_data[key]['CTD_lon'][btl_data[key]['Bottle']==poc[key]])
            poc_lat.append(btl_data[key]['CTD_lat'][btl_data[key]['Bottle']==poc[key]])
        #two samples without position

    dna = {'PS137_018_01' :[13],
           'PS137_022_01' :[12],
           'PS137_026_01' :[12],
           'PS137_028_01' :[6,8,10,20,23],
           'PS137_033_01' :[10],
           'PS137_036_01' :[7],
           'PS137_041_01' :[9,10,14,16,20],
           'PS137_054_01' :[5,20,22]}

    dna_lon = list()
    dna_lat = list()

    for key in dna:
        for i in dna[key]:
            if key == 'PS137_022_01' or key == 'PS137_026_01':
                dna_lon.append(btl_data[key]['Dship_lon'][btl_data[key]['Bottle']==i].to_numpy()[0])
                dna_lat.append(btl_data[key]['Dship_lat'][btl_data[key]['Bottle']==i].to_numpy()[0])
            else:
                dna_lon.append(btl_data[key]['CTD_lon'][btl_data[key]['Bottle']==i].to_numpy()[0])
                dna_lat.append(btl_data[key]['CTD_lat'][btl_data[key]['Bottle']==i].to_numpy()[0])

    fig = plt.figure()
    #ax = plt.axes(projection=ccrs.Mercator())
    colors = plt.rcParams["axes.prop_cycle"]()
    fig.set_tight_layout(True)

    contours = plt.contourf(lons,lats,elev,levels=40,alpha=0.7)
    contourlines = plt.contour(lons,lats,elev,levels = 40,colors='black',linestyles='solid',linewidths=0.5,alpha=0.3)
    plt.colorbar(contours,label='Depth in m')

    plt.scatter(parse("6° 15.32'W"),parse("82° 53.83'N"),color='red',marker='*',s=100,label='Aurora Vent Site')

    station=['22','26','28','33','36','41','54','55']
    i=0
    for key in start_lat:
        plt.scatter(start_lon[key],start_lat[key],marker='v',color='black')
        plt.text(start_lon[key]+0.007,start_lat[key]+0.00005,station[i])
        plt.scatter(end_lon[key],end_lat[key],marker='^',color='black')
        if key == 'PS137_022_01' or  key == 'PS137_026_01':
            plt.plot(profile_plot[key]['Dship_lon'],profile_plot[key]['Dship_lat'],linewidth='0.8',linestyle='--',color='black')  #plot estimated track
        else:
            plt.plot([start_lon[key],profile_data[key]['CTD_lon'].iloc[profile_data[key]['CTD_lon'].first_valid_index()]],
                     [start_lat[key],profile_data[key]['CTD_lat'].iloc[profile_data[key]['CTD_lat'].first_valid_index()]],linewidth='0.8',linestyle='--',color='black') #estimated track between start & 1st Posi
            plt.plot([end_lon[key],profile_data[key]['CTD_lon'].iloc[profile_data[key]['CTD_lon'].last_valid_index()]],
                     [end_lat[key],profile_data[key]['CTD_lat'].iloc[profile_data[key]['CTD_lat'].last_valid_index()]],linewidth='0.8',linestyle='--',color='black') #estimated track between last Posi & end
            plt.plot(profile_plot[key]['CTD_lon'],profile_plot[key]['CTD_lat'],color='black',linewidth='0.8')       #plot posi track
        i=i+1

    match tracer_type:      #mark different tracer types on the map
        case 'dna':
            plt.scatter(dna_lon,dna_lat,marker='x',color='blue',label='DNA sample')
        case 'poc':
            plt.scatter(poc_lon,poc_lat,marker='x',color='red',label='POC sample')
        case 'dna_poc':
            plt.scatter(dna_lon,dna_lat,marker='x',color='blue',label='DNA sample')
            plt.scatter(poc_lon,poc_lat,marker='x',color='red',label='POC sample')
        case 'he':
            btl_data_he = dict()
            for key in btl_data:
                he_station = btl_data[key].dropna(subset=['delta3He'])      # subsets rows with helium samples
                if he_station.size != 0:
                    btl_data_he[key] = he_station       # takes all stations with helium samples in new dict

            del btl_data_he['PS137_018_01']        #background station, south Lena Trough
            del btl_data_he['PS137_049_01']        #atwaice cal
            del btl_data_he['PS137_058_01']        #Lucky B
            del btl_data_he['PS137_061_01']        #Lucky B

            for key in btl_data_he:
                try:
                    plt.scatter(btl_data_he[key]['CTD_lon'],btl_data_he[key]['CTD_lat'], marker='x', color='orange')
                except:
                    print(key+'CTD posi not available, using Dship posi')
                    plt.scatter(btl_data_he[key]['Dship_lon'],btl_data_he[key]['Dship_lat'], marker='x', color='orange')


    plt.ylabel('Latitude')
    plt.xlabel('Longitude')
    plt.legend(loc='lower left')
    #plt.title('Aurora Vent Site')

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot_largemap_lambert(vent_loc,vent_label):

    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import matplotlib.path as mpath
    import xarray as xa

    bathy_path = '/home/jonathan/Dokumente/model/inputs/Gridbuilder/Bathymetry/GEBCO_17_Apr_2024_b985334cc502/gebco_2023_n88.4619_s77.6514_w-24.1699_e13.1836.nc'
    bathy = xa.open_dataset(bathy_path)

    ax = plt.axes(projection=ccrs.LambertConformal(0))#central_longitude=vent_loc[1],central_latitude=vent_loc[0]))
    #ax.set_extent([0,360,60,90], crs=ccrs.PlateCarree())     #W,E, S, N
    ax.coastlines()

    # Lon and Lat Boundaries
    xlim = [-40, 40]
    ylim = [60, 90]
    lower_space = 20 # this needs to be manually increased if the lower arched is cut off by changing lon and lat lims

    rect = mpath.Path([[xlim[0], ylim[0]],
                       [xlim[1], ylim[0]],
                       [xlim[1], ylim[1]],
                       [xlim[0], ylim[1]],
                       [xlim[0], ylim[0]],
                       ]).interpolated(20)

    proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax) - ax.transData
    rect_in_target = proj_to_data.transform_path(rect)

    ax.set_boundary(rect_in_target)
    ax.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

    #plot vent location
    data_projection = ccrs.PlateCarree()
    ax.plot(vent_loc[0], vent_loc[1], marker='*', color='tab:red', transform=data_projection)
    ax.contourf(lons, lats, data, transform=data_projection)

    plt.show()

def plot_largemap_stereo(vent_loc,vent_label):

    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import matplotlib.path as mpath

    ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=vent_loc[1]))
    #ax.set_extent([0,360,60,90], crs=ccrs.NorthPolarStereo())     #W,E, S, N
    ax.set_ylim(bottom=60)
    ax.coastlines()

    #plt.scatter(vent_loc[0],vent_loc[1],color='red',marker='*',s=100,label=vent_label)

    plt.show()


