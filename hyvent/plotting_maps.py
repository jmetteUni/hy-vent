#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 13:27:30 2024

@author: jmette@uni-bremen.de
"""

def plot_map(profile_data, btl_data, bathy, tracer_type='None', path_save='None'):
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

    from lat_lon_parser import parse
    import numpy as np
    import matplotlib.pyplot as plt

    #use Arial from custom path
    # import matplotlib.font_manager

    # font_dir = ['/home/jonathan/Downloads/arial_font/']
    # for font in matplotlib.font_manager.findSystemFonts(font_dir):
    #     matplotlib.font_manager.fontManager.addfont(font)

    # plt.rcParams['font.family'] = 'Arial'#'sans-serif'
    #plt.rcParams.update({'font.size': 10})

    # subset profile data for plots

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

    #  DNA and POC positions
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

    lons, lats, elev = bathy
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
    plt.tight_layout()
    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot_track(profile_data, btl_data, vent_loc='None', bathy='None'):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    #plot bathymetry
    if bathy != 'None':
        contourlines = ax.contour(bathy[0],bathy[1],-bathy[2],levels=23 ,colors='black',linestyles='solid',linewidths=0.5,alpha=0.3)
        ax.clabel(contourlines, inline=True, fontsize=6, fmt='%d', colors = 'black')       #use this line to plot with lables instead colorbar
        contours = plt.contourf(bathy[0],bathy[1],bathy[2],levels=40,alpha=0.7)         #use this two lines to plot with colorbar instead labels
        plt.colorbar(contours,label='Depth in m')

    #plot vent
    if vent_loc != 'None':
        ax.scatter(vent_loc[0],vent_loc[1],color='red',marker='*',s=100,label='Aurora Vent Site')

    #plot tracks
    data_list = [d for _, d in profile_data.groupby(['Station','SN'])]
    for profile in data_list:
        lat = profile['CTD_lat'].fillna(profile['Dship_lat'])       #fill gaps in acoustic position with Dship positions
        lon = profile['CTD_lon'].fillna(profile['Dship_lon'])
        #lat = profile['Dship_lat']
        #lon = profile['Dship_lon']
        ax.scatter(lon,lat,color='black',s=0.001)

        #plot labels & marker
        ax.scatter(lon.iloc[0],lat.iloc[0],marker='v',color='black')
        ax.text(lon.iloc[0]+0.007,lat.iloc[0]+0.00005,profile['Station'].iloc[0][1:-3],fontsize='small')
        ax.scatter(lon.iloc[-1],lat.iloc[-1],marker='^',color='black')

    plt.tight_layout()
    plt.show()

def plot_var_in_2D(data,var,min_dep,max_dep,nth_point,vent_loc='None',bathy='None',path_save='None'):
    """
    This function plots a variable in x-y-coordinates with colorcoding. It can be filtered by a depth range and only every n-th point or line segments can be plotted for better visibility. Bathymetry lines can be plotted additonally.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data of the CTD station or multiple stations with variable as columns as columns. One column has to be the station designation.
    var : string
        Variable to plot. Must be a key for a column in data and is only plotted for values with coordinates in the columns "CTD_lon"/"CTD_lat" or "Dship_lon"/"Dship_lat".
    min_dep : int
        Minimum depth below which datapoints should be plotted.
    max_dep :  int
        Maximum depth above which datapoints should be plotted.
    nth_point : int
        Every n-th point is plotted, for example nth_point=1 for every measurement. nth_point=0 activates the plotting of line segments.
    vent_loc : tuple, optional
        Longitude and latitude of a point of interest, plotted as a red star. Default is 'None"'
    bathy : tuple, optional
        Tuple of longitude, latitude and elevation of bathymetry data. The default is 'None'.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    import warnings
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import LineCollection
    from hyvent.misc import get_var

    def colored_line(x, y, c, ax, **lc_kwargs):
        """
        Plot a line with a color specified along the line by a third value.

        It does this by creating a collection of line segments. Each line segment is
        made up of two straight lines each connecting the current (x, y) point to the
        midpoints of the lines connecting the current point with its two neighbors.
        This creates a smooth line with no gaps between the line segments.
        Function by https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html#sphx-glr-gallery-lines-bars-and-markers-multicolored-line-py

        Parameters
        ----------
        x, y : array-like
            The horizontal and vertical coordinates of the data points.
        c : array-like
            The color values, which should be the same size as x and y.
        ax : Axes
            Axis object on which to plot the colored line.
        **lc_kwargs
            Any additional arguments to pass to matplotlib.collections.LineCollection
            constructor. This should not include the array keyword argument because
            that is set to the color argument. If provided, it will be overridden.

        Returns
        -------
        matplotlib.collections.LineCollection
            The generated line collection representing the colored line.
        """
        if "array" in lc_kwargs:
            warnings.warn('The provided "array" keyword argument will be overridden')

        # Default the capstyle to butt so that the line segments smoothly line up
        default_kwargs = {"capstyle": "butt"}
        default_kwargs.update(lc_kwargs)

        # Compute the midpoints of the line segments. Include the first and last points
        # twice so we don't need any special syntax later to handle them.
        x = np.asarray(x)
        y = np.asarray(y)
        x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
        y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

        # Determine the start, middle, and end coordinate pair of each line segment.
        # Use the reshape to add an extra dimension so each pair of points is in its
        # own list. Then concatenate them to create:
        # [
        #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
        #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
        #   ...
        # ]
        coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
        coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
        coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
        segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

        lc = LineCollection(segments, **default_kwargs)
        lc.set_array(c)  # set the colors of each segment

        return ax.add_collection(lc)

    #get properties for variable
    label, color, cmap = get_var(var)

    fig, ax = plt.subplots(figsize=(8,6))

    #plot bathymetry
    if bathy != 'None':
        contourlines = ax.contour(bathy[0],bathy[1],-bathy[2],levels=40, colors='black',linestyles='solid',linewidths=0.5,alpha=0.3)
        ax.clabel(contourlines, inline=True, fontsize=6, fmt='%d', colors = 'black')

    #plot vent
    if vent_loc != 'None':
        ax.scatter(vent_loc[0],vent_loc[1],color='red',marker='*',s=100,label='Aurora Vent Site')

    if var == 'dORP':       #ommit positve dORP values, as only negativ ones are of interest
        data = data[data[var]<=0]
    if var == 'Delta_potemperature':    #ommit negative delta theta values, as only positiv ones are of interest
        data = data[data[var]>=0]
    if var == 'delta3He':       #for delta3He, data is bottle data, therefore rename depth column
        data = data.rename({'DepSM_mean':'DEPTH'},axis=1)

    data_nb = data[(data['DEPTH']>min_dep) & (data['DEPTH']<max_dep)]       #subset by expected plume depth

    #plot data
    if 'delta3He' in var:
        lat = data_nb['CTD_lat'].fillna(data_nb['Dship_lat'])       #fill gaps in acoustic position with Dship positions
        lon = data_nb['CTD_lon'].fillna(data_nb['Dship_lon'])
        var_plot = ax.scatter(lon,lat,c=data_nb[var],cmap=cmap, edgecolors='black',linewidth=0.1)
        fig.colorbar(var_plot,label=label)

    else:
        #plot station labels for continous measurements
        station_list = [d for _, d in data_nb.groupby(['Station'])]
        for station in station_list:
            lat = station['CTD_lat'].fillna(station['Dship_lat'])       #fill gaps in acoustic position with Dship positions
            lon = station['CTD_lon'].fillna(station['Dship_lon'])
            ax.text(lon.iloc[0]+0.007,lat.iloc[0]+0.0005,station['Station'].iloc[0][1:-3],fontsize='small')

        #plot data, as colored line segments (nth_point=0) or as scatter points(nth_point>0)
        if nth_point == 0:      #plot line segments
            data_list = [d for _, d in data_nb.groupby(['Station','SN'])]
            for profile in data_list:
                lat = profile['CTD_lat'].fillna(profile['Dship_lat'])       #fill gaps in acoustic position with Dship positions
                lon = profile['CTD_lon'].fillna(profile['Dship_lon'])
                lines = colored_line(lon,lat,profile[var],ax,linewidth=2,cmap=cmap)
            fig.colorbar(lines)
        else:       #plot single dots
            data_nb = data_nb.iloc[::nth_point,:]        #plot only every n-th point for better visibility
            data_list = [d for _, d in data_nb.groupby(['Station','SN'])]
            for profile in data_list:
                lat = profile['CTD_lat'].fillna(profile['Dship_lat'])       #fill gaps in acoustic position with Dship positions
                lon = profile['CTD_lon'].fillna(profile['Dship_lon'])
                var_plot = ax.scatter(lon,lat,c=profile[var],s=15,cmap=cmap, edgecolors='black',linewidth=0)
            fig.colorbar(var_plot,label=label)
            #ax.plot(data_nb['Dship_lon'],data_nb['Dship_lat'])

    #set plot limits
    # lat = profile['CTD_lat'].fillna(profile['Dship_lat'])       #fill gaps in acoustic position with Dship positions
    # lon = profile['CTD_lon'].fillna(profile['Dship_lon'])
    # lon_limits = [lon.min()-(lon.max()-lon.min())/15,lon.max()+(lon.max()-lon.min())/15]
    # lat_limits = [lat.min()-(lat.max()-lat.min())/15,lat.max()+(lat.max()-lat.min())/15]
    # ax.set_xlim(lon_limits)
    # ax.set_ylim(lat_limits)
    #plt.gca().set_aspect(5)
    plt.xlabel('Longitude')
    plt.ylabel('Longitude')
    plt.tight_layout()

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

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
    import matplotlib.pyplot as plt
    import pandas as pd
    from hyvent.misc import get_var
    from hyvent.misc import add_castno
    import numpy as np

    label, color, cmap = get_var(var)

    data = add_castno(data, window_size)
    data = data[(data['DEPTH']>min_dep) & (data['DEPTH']<max_dep)]

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
        if var == 'dORP':
            var_max = cast[var].min()
        else:
            var_max = cast[var].max()

        # add the values to the DataFrame
        data_binned.loc[i] = [station, lon, lat, var_max]

    if var == 'Delta_potemperature':
        var_min = 0
        var_max = 0.0138
    elif var == 'dORP':
        var_min = -0.135
        var_max = 0
    else:
        var_min = data_binned['var_max'].min()
        var_max = data_binned['var_max'].max()

    var_plot = plt.scatter(data_binned['lon'], data_binned['lat'], c=data_binned['var_max'], cmap = get_var(var)[2],edgecolors='black',linewidth=0.3, vmin=var_min, vmax=var_max)

    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig.colorbar(var_plot,label=label)
    #change colorbar ticks settings at some point?
    plt.tight_layout()

    plt.show()


#%% lambert projections tests
# def plot_largemap_lambert(vent_loc,vent_label):

#     import cartopy.crs as ccrs
#     import matplotlib.pyplot as plt
#     import matplotlib.path as mpath
#     import xarray as xa

#     bathy_path = '/home/jonathan/Dokumente/model/inputs/Gridbuilder/Bathymetry/GEBCO_17_Apr_2024_b985334cc502/gebco_2023_n88.4619_s77.6514_w-24.1699_e13.1836.nc'
#     bathy = xa.open_dataset(bathy_path)

#     ax = plt.axes(projection=ccrs.LambertConformal(0))#central_longitude=vent_loc[1],central_latitude=vent_loc[0]))
#     #ax.set_extent([0,360,60,90], crs=ccrs.PlateCarree())     #W,E, S, N
#     ax.coastlines()

#     # Lon and Lat Boundaries
#     xlim = [-40, 40]
#     ylim = [60, 90]
#     lower_space = 20 # this needs to be manually increased if the lower arched is cut off by changing lon and lat lims

#     rect = mpath.Path([[xlim[0], ylim[0]],
#                        [xlim[1], ylim[0]],
#                        [xlim[1], ylim[1]],
#                        [xlim[0], ylim[1]],
#                        [xlim[0], ylim[0]],
#                        ]).interpolated(20)

#     proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax) - ax.transData
#     rect_in_target = proj_to_data.transform_path(rect)

#     ax.set_boundary(rect_in_target)
#     ax.set_extent([xlim[0], xlim[1], ylim[0] - lower_space, ylim[1]])

#     #plot vent location
#     data_projection = ccrs.PlateCarree()
#     ax.plot(vent_loc[0], vent_loc[1], marker='*', color='tab:red', transform=data_projection)
#     ax.contourf(lons, lats, data, transform=data_projection)

#     plt.show()

# def plot_largemap_stereo(vent_loc,vent_label):

#     import cartopy.crs as ccrs
#     import matplotlib.pyplot as plt
#     import matplotlib.path as mpath

#     ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=vent_loc[1]))
#     #ax.set_extent([0,360,60,90], crs=ccrs.NorthPolarStereo())     #W,E, S, N
#     ax.set_ylim(bottom=60)
#     ax.coastlines()

#     #plt.scatter(vent_loc[0],vent_loc[1],color='red',marker='*',s=100,label=vent_label)

#     plt.show()


