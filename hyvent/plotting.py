#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:25:58 2024

@author: jonathan
"""

def plot_hist_he(data,depth_max,depth_min,bins,path_save='None'):
    """
    Plots histogram of helium isotope data.

    Parameters
    ----------
    data : dictionary or pandas dataframe
        dictionary of pandas dataframes with station, bottle number, Ne/He, delta3He and delta22Ne as columns or one pandas dataframe.
    depth_max : int
        Maximum depth of the samples to plot.
    depth_min : int
        Minimum depth of the samples to plot..
    bins : int
        Number of bins in the histogram.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    import pandas as pd
    import matplotlib.pyplot as plt

    label_sel = 'Samples between '+str(depth_min)+' - '+str(depth_max)+'m'
    if type(data) == dict:
        he_data = pd.concat(data.values(),ignore_index=True)
    else:
        he_data = data

    plt.figure()
    plt.hist(he_data['delta3He'],bins=bins,label='All samples')
    plt.hist(he_data['delta3He'][(he_data['DepSM_mean']>depth_min) & (he_data['DepSM_mean']<depth_max)],bins=bins,label=label_sel,alpha=0.7)

    plt.xlabel('Delta 3He in %')
    plt.ylabel('Number of samples')
    plt.legend()

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot_hist_dorp(data,depth_max,depth_min,bins,ranges,path_save='None'):
    """
    Plots histogram of dORP of CTD profile data.

    Parameters
    ----------
    data : dictionary or pandas dataframe
        dictionary of pandas dataframes or one pandas dataframe with CTD profile data.
    depth_max : int
        Maximum depth of the samples to plot.
    depth_min : int
        Minimum depth of the samples to plot..
    bins : int
        Number of bins in the histogram.
    ranges : tuple of int
        Upper and lower limit of bins in the histogram.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    import pandas as pd
    import matplotlib.pyplot as plt

    label_sel = 'Samples between '+str(depth_min)+' - '+str(depth_max)+'m'
    if type(data) == dict:
        dorp_data = pd.concat(data.values(),ignore_index=True)
    else:
        dorp_data = data

    plt.figure()
    plt.hist(dorp_data['dORP'],bins=bins,range=ranges,label='All samples',log=True)
    plt.hist(dorp_data['dORP'][(dorp_data['DEPTH']>depth_min) & (dorp_data['DEPTH']<depth_max)],bins=bins,range=ranges,label=label_sel,log=True)

    plt.xlabel('dORP')
    plt.ylabel('Number of samples')
    plt.legend()

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot_section(profile_data,key,xvar,yvar,zvar,depth,levels,path_save='None'):
    """
    Plots interpolated section of a selected variable from a dataframe/station of profile data along an x-variable. The dataframe is selected by a key in the dictionary.

    Parameters
    ----------
    profile_data : dictionary
        Dictionary of pandas dataframes or one pandas dataframe with CTD profile data.
    key : string
        Key of the dataframe/station to plot.
    xvar : string
        Key of dataframe column to use as x-coordinate.
    yvar : string
        Key of dataframe column to use as y-coordinate.
    zvar : string
        Key of dataframe column to use as z-coordinate.
    depth : int
        Upper depth limit of the plot.
    levels : int
        Number of individual contour lines.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """

    from matplotlib import colors
    import matplotlib.pyplot as plt

    label_title = 'Interpolated crossection for '+key

    if type(profile_data) == dict:
        station = profile_data[key]
    else:
        station = profile_data

    station = station[station['DEPTH']>depth]
    if station.empty:
        print('Station '+key+' has no values lower then '+str(depth))
        return
    station = station[[xvar,yvar,zvar]]
    station = station.dropna(subset=[xvar,yvar]) #drop rows with NAN in x or y coordinates
    if station.empty:
        print('Station '+key+' has no values in one of the variables used')
        return

    plt.figure()
    plt.title(label_title)

    contourf = plt.tricontourf(station[xvar],station[yvar],station[zvar],levels=levels,cmap='RdYlBu',norm=colors.CenteredNorm())
    track = plt.scatter(station[xvar],station[yvar],s=0.05,color='black',alpha=0.3)

    plt.gca().invert_yaxis()
    plt.colorbar(contourf,label=zvar)
    plt.xlabel(xvar)
    plt.ylabel(yvar)

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot3Ddistr(profile_data,key,var,bathy,vent_loc,path_save='None'):
    """
    Plots positions of a variable from one dataframe of a dictionary of profile data in 3D space with the data as color coding and the bathymetry as a 3D surface.

    Parameters
    ----------
    profile_data : dictionary
        Dictionary of pandas dataframes or one pandas dataframe with CTD profile data.
    key : string
        Key of the dataframe/station to plot.
    var : string
        Key of dataframe column to plot.
    bathy : tuple of three int
        Bathymetry data with longitude and latitude coordinates and elevation.
    vent_loc : tuple of three float
        Longitude, latitude and depth of a position of interest to mark in the plot.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    import matplotlib.pyplot as plt
    from matplotlib import colors

    if type(profile_data) == dict:
        station = profile_data[key]
    else:
        station = profile_data

    label_title = '3D distribution for '+key

    station = profile_data[key]
    station = station[station.DEPTH>2500]
    lons = bathy[0]
    lats = bathy[1]
    elev = bathy[2]*(-1)

    n=5        #plot only nth point
    station = station.iloc[::n,:]

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection='3d')

    vent = ax.scatter(vent_loc[0],vent_loc[1], vent_loc[2], marker='*', color='r', s=100)
    orp1 = ax.scatter(station['CTD_lon'], station['CTD_lat'],station['DEPTH'], c=station[var],marker='.',alpha=1,cmap='RdYlBu',norm=colors.CenteredNorm(),edgecolor='black',linewidth=0.2)
    surf = ax.plot_surface(lons,lats,elev,alpha=0.5,linewidth=0.5,color='white',antialiased=False)

    ax.invert_zaxis()
    fig.colorbar(orp1, label=var,fraction=0.02)

    plt.title(label_title)
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')
    ax.set_zlabel('Depth in m')
    ax.xaxis.labelpad=10
    ax.yaxis.labelpad=10
    ax.zaxis.labelpad=10
    ax.view_init(8, -173)

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot2D_station(profile_data,key,var,vent_loc,depth_min,depth_max,bathy=False,path_save='None'):
    """
    Plots positions of a variable from one dataframe/station of a dictionary of profile data on a 2D map. Samples to plot can be selected by depth. Plot can be combined with bathymetry. If no CTD positions are found, ship positions are used instead.

    Parameters
    ----------
    profile_data : dictionary
        Dictionary of pandas dataframes with CTD profile data.
    key : string
        Key of the dataframe/station to plot.
    var : string
        Key of dataframe column to plot.
    vent_loc : tuple of three float
        Longitude, latitude and depth of a position of interest to mark in the plot.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.
    depth_min : int
        Minimum depth of samples to plot.
    depth_max : int
        Maximum dpeth of samples to plot.
    bathy : tuple of three int, optional
        Bathymetry data with longitude and latitude coordinates and elevation.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    import matplotlib.pyplot as plt

    label_title = '2D Distribution for '+key+'\n'+'in '+str(depth_min)+' - '+str(depth_max)+' m'

    if bathy != False:
        lons = bathy[0]
        lats = bathy[1]
        elev = bathy[2]
    station = profile_data[key]
    station = station[(station['DEPTH'] > depth_min) & (station['DEPTH'] < depth_max)]

    lon = station['CTD_lon']
    lat = station['CTD_lat']
    if (lon.dropna().empty) | (lat.dropna().empty):
        print('No CTD position data, falling back to ship position')
        lon = station['Dship_lon']
        lat = station['Dship_lat']

    fig, ax = plt.subplots()
    fig.set_tight_layout(True)

    if bathy != False:
        contours = plt.contourf(lons,lats,elev,levels=40,alpha=0.7)
        contourlines = plt.contour(lons,lats,elev,levels = 40,colors='black',linestyles='solid',linewidths=0.5,alpha=0.3)

    vent = ax.scatter(vent_loc[0],vent_loc[1], marker='*', color='r', s=100)
    var_plot = ax.scatter(lon, lat,c=station[var],s=5)


    #plt.colorbar(contours,label='Depth in m')

    ax.set_title(label_title)
    fig.colorbar(var_plot, label=var)
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot2D_all_stations(profile_data,var,vent_loc,depth_min,depth_max,bathy=False,path_save='None'):
    """
    Plots positions of a variable from all dataframes/stations of a dictionary of profile data on a 2D map. Samples to plot can be selected by depth. Plot can be combined with bathymetry. If no CTD positions are found, ship positions are used instead.

    Parameters
    ----------
    profile_data : dictionary
        Dictionary of pandas dataframes with CTD profile data.
    key : string
        Key of the dataframe/station to plot.
    var : string
        Key of dataframe column to plot.
    vent_loc : tuple of three float
        Longitude, latitude and depth of a position of interest to mark in the plot.
    depth_min : int
        Minimum depth of samples to plot.
    depth_max : int
        Maximum dpeth of samples to plot.
    bathy : tuple of three int, optional
        Bathymetry data with longitude and latitude coordinates and elevation.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    import matplotlib.pyplot as plt

    label_title = '2D Distribution in '+str(depth_min)+' - '+str(depth_max)+' m'

    if bathy != False:
        lons = bathy[0]
        lats = bathy[1]
        elev = bathy[2]

    fig, ax = plt.subplots()
    fig.set_tight_layout(True)

    if bathy != False:
        contours = plt.contourf(lons,lats,elev,levels=40,alpha=0.7)
        contourlines = plt.contour(lons,lats,elev,levels = 40,colors='black',linestyles='solid',linewidths=0.5,alpha=0.3)

    for key in profile_data:

        station = profile_data[key]
        station = station[(station['DEPTH'] > depth_min) & (station['DEPTH'] < depth_max)]

        lon = station['CTD_lon']
        lat = station['CTD_lat']
        if (lon.dropna().empty) | (lat.dropna().empty):
            print('No CTD position data, falling back to ship position')
            lon = station['Dship_lon']
            lat = station['Dship_lat']

        var_plot = ax.scatter(lon, lat,c=station[var],s=5)

    vent = ax.scatter(vent_loc[0],vent_loc[1], marker='*', color='r', s=100)

    ax.set_title(label_title)
    fig.colorbar(var_plot, label=var)
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot2D_all_stations_btl(btl_data,var,vent_loc,depth_min,depth_max,bathy=False,path_save='None'):
    """
    Plots positions of a variable of bottle data from all dataframes/stations of a dictionary on a 2D map. Samples to plot can be selected by depth. Plot can be combined with bathymetry. If no CTD positions are found, ship positions are used instead.

    Parameters
    ----------
    btl_data : dictionary
        Dictionary of pandas dataframes with CTD bottle data.
    key : string
        Key of the dataframe/station to plot.
    var : string
        Key of dataframe column to plot.
    vent_loc : tuple of three float
        Longitude, latitude and depth of a position of interest to mark in the plot.
    depth_min : int
        Minimum depth of samples to plot.
    depth_max : int
        Maximum dpeth of samples to plot.
    bathy : tuple of three int, optional
        Bathymetry data with longitude and latitude coordinates and elevation.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    import matplotlib.pyplot as plt

    label_title = '2D Distribution in '+str(depth_min)+' - '+str(depth_max)+' m'

    if bathy != False:
        lons = bathy[0]
        lats = bathy[1]
        elev = bathy[2]

    fig, ax = plt.subplots()
    fig.set_tight_layout(True)

    if bathy != False:
        #contours = plt.contourf(lons,lats,elev,levels=40,alpha=0.3)
        contourlines = plt.contour(lons,lats,elev,levels = 40,colors='black',linestyles='solid',linewidths=0.5,alpha=0.1)

    for key in btl_data:

        station = btl_data[key]
        station = station[(station['DepSM_mean'] > depth_min) & (station['DepSM_mean'] < depth_max)]

        lon = station['CTD_lon']
        lat = station['CTD_lat']
        if (lon.dropna().empty) | (lat.dropna().empty):
            print('No CTD position data, falling back to ship position')
            lon = station['Dship_lon']
            lat = station['Dship_lat']

        var_plot = ax.scatter(lon, lat,c=station[var],s=5)

    vent = ax.scatter(vent_loc[0],vent_loc[1], marker='*', color='r', s=100)


    ax.set_title(label_title)
    fig.colorbar(var_plot, label=var)
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()








