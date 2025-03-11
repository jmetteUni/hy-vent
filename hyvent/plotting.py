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
    if isinstance(data,dict):
        he_data = pd.concat(data.values(),ignore_index=True)
    else:
        he_data = data

    plt.figure()
    plt.hist(he_data['delta3He'],bins=bins,label='All samples')
    plt.hist(he_data['delta3He'][(he_data['DepSM_mean']>depth_min) & (he_data['DepSM_mean']<depth_max)],bins=bins,label=label_sel,alpha=0.7)

    plt.xlabel('Delta 3He in %')
    plt.ylabel('Number of samples')
    plt.legend()
    plt.tight_layout()
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
    if isinstance(data,dict):
        dorp_data = pd.concat(data.values(),ignore_index=True)
    else:
        dorp_data = data

    plt.figure()
    plt.hist(dorp_data['dORP'],bins=bins,range=ranges,label='All samples',log=True)
    plt.hist(dorp_data['dORP'][(dorp_data['DEPTH']>depth_min) & (dorp_data['DEPTH']<depth_max)],bins=bins,range=ranges,label=label_sel,log=True)

    plt.xlabel('dORP')
    plt.ylabel('Number of samples')
    plt.legend()
    plt.tight_layout()
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

    if isinstance(profile_data,dict):
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
    plt.tight_layout()
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

    if isinstance(profile_data,dict):
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
    plt.tight_layout()
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
    plt.tight_layout()
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
        contourlines = plt.contour(lons,lats,elev,levels = 40,colors='black',linestyles='solid',linewidths=0.5,alpha=0.5)

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
    plt.tight_layout()
    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def depth_plot(data,xvar,yvar,depth_min,background='None',path_save='None'):
    """
    This function plots a variable of one or multiple CTD stations against depth. Additionally a special station called background is plotted separetely.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data of the CTD station or multiple stations with variable as columns. One column has to be the station designation.
    xvar : string
        Variable to plot as x-variable. Must be a column key in data and background. If xvar is "delta3He" it is plotted as a scatterplot, else as a lineplot.
    yvar : string
        Variable to plot as y-variable. Must be a column key in data and background. Should be a type of depth or pressure coordinate and is plotted inverted.
    depth_min : int
        Minimal cutoff depth to plot.
    background : pandas dataframe, optional
        Dataframe with the data of the station which should be plotted separately. The default is 'None'.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    from hyvent.misc import get_var
    import matplotlib.pyplot as plt
    import pandas as pd

    xlabel, xcolor, cmap = get_var(xvar)

    data = data[data[yvar]>depth_min]
    if isinstance(background,pd.DataFrame):
        background = background[background[yvar]>depth_min]
        bg_list = [d for _, d in background.groupby(['SN'])]

    data_list = [d for _, d in data.groupby(['Station','SN'])]

    plt.figure(figsize=(6,6))
    plt.tight_layout()
    plt.rcParams['axes.autolimit_mode'] = 'data'
    #plt.rcParams['axes.autolimit_mode'] = 'round_numbers'


    if (xvar == 'delta3He') | (xvar == 'Delta_delta3He'):
        for station in data_list:
            plt.scatter(station[xvar],station[yvar],color=get_var(xvar)[1])
        if isinstance(background,pd.DataFrame):
            for device in bg_list:
                plt.scatter(device[xvar],device[yvar],color='black',linewidth=1)
    else:
        for station in data_list:
            station = station.sort_index()
            plt.plot(station[xvar],station[yvar],color=xcolor,linewidth=1)
            #plt.plot(station[xvar],station[yvar],linewidth=1,label=station['SN'].iloc[0])  #plot every station, SN in different color

            #plt.title(station['Station'].iloc[0]+', '+station['SN'].iloc[0])   #plotting station and SN for debugging
        if isinstance(background,pd.DataFrame):
            for device in bg_list:
                plt.plot(device[xvar],device[yvar],color='black',linewidth=1)
                #plt.plot(device[xvar],device[yvar],linewidth=1,linestyle='--',label=device['SN'].iloc[0])  #plot every station, SN in different color


    if (xvar == 'Delta_potemperature') | (xvar == 'Delta_Sigma3') | (xvar == 'Delta_delta3He') | (xvar == 'Delta_Neph_outl(volts)') | (xvar == 'Delta_Neph_smoo(volts)') | (xvar == 'Delta_PSAL'):
        plt.axvline(0, color = 'black',alpha=0.3)
    plt.gca().invert_yaxis()
    plt.ylabel(get_var(yvar)[0])
    plt.xlabel(get_var(xvar)[0])
    #plt.legend()

    #plt.locator_params(axis='x', nbins=8)
    plt.tight_layout()
    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def time_plot(data,station,depth_min,path_save='None'):
    """
    This function plots four quantities (dpeht, potential temperature, dORP, Sigma3) of one CTD station against time.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data of the CTD station or multiple stations with variable as columns as columns. One column has to be the station designation.
    station : string
        Indentifier of the station, which should be plotted. Needs to be a column value of all rows of this station in data.
    depth_min : int
        Minimal cutoff depth to plot.
    path_save : string, optional
        Path to save the plot as a png with dpi=300. The default is 'None'.

    Returns
    -------
    None.

    """
    from hyvent.misc import get_var
    import pandas as pd
    import matplotlib.pyplot as plt

    data = data[data['DEPTH']>depth_min]

    station = data[data['Station']==station]
    station['datetime'] = pd.to_datetime(station['datetime'])

    fig, ax = plt.subplots()
    fig.set_tight_layout(True)
    fig.set_figwidth(12)
    fig.set_figheight(5)

    var = 'DEPTH'
    color = get_var(var)[1]
    label = get_var(var)[0]
    depth, = ax.plot(station['datetime'], station[var], color = color, label = label)
    ax.set_ylabel(label)
    #ax.yaxis.label.set_color(color)

    twin1 = ax.twinx()
    var = 'potemperature'
    color = get_var(var)[1]
    label = get_var(var)[0]
    temp, = twin1.plot(station['datetime'], station[var], color = color, label = label)
    twin1.set_ylabel(label)
    #twin1.yaxis.label.set_color(color)

    twin2 = ax.twinx()
    twin2.spines.right.set_position(("axes", 1.12))
    var = 'dORP'
    color = get_var(var)[1]
    label = get_var(var)[0]
    orp, = twin2.plot(station['datetime'], station[var], color = color, label = label)
    twin2.set_ylabel(label)
    #twin2.yaxis.label.set_color(color)

    twin3 = ax.twinx()
    twin3.spines.right.set_position(("axes", 1.25))
    var = 'Sigma3'
    color = get_var(var)[1]
    label = get_var(var)[0]
    dens, = twin3.plot(station['datetime'], station[var],color = color, label = label)
    twin3.set_ylabel(label)
    #twin3.yaxis.label.set_color(color)

    ax.invert_yaxis()

    ax.set_xlabel('Time')
    fig.legend(bbox_to_anchor=(1,0.3), bbox_transform=ax.transAxes)
    plt.tight_layout()
    if path_save != 'None':
        plt.savefig(path_save, dpi=300)
    plt.show()

def plot_ts_zoom(data, c_var, min_dep, max_dep, p_ref, lon, lat):
    """
    This function plots a zoomed version of the data of a CTD cast as a T-S diagram with potential temperature and practical salinity. The samples can be filtered by a depth range. For the isopycnal calculation a reference density, longitude and latitude is used.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with CTD data containing the potential temperature and practical salinity.
    c_var : string
        Variable which should be used as third dimension with colorbar.
    min_dep : int
        Minimal depth of the samples used.
    max_dep : int
        Maximal depth of the samples used.
    p_ref : int
        Reference pressure used for calculating the isopycnals of the potential density anomaly in dbar.
    lon : float
        Longitude of the samples.
    lat : float
        Latitude of the samples.

    Returns
    -------
    None.

    """
    import matplotlib.pyplot as plt
    import numpy as np
    import gsw
    from hyvent.misc import get_var, add_castno
    import pandas as pd

    # calculate the density values for the grid
    def calculate_density(pt_grid, psal_grid, p, lon, lat):
        SA = gsw.conversions.SA_from_SP(psal_grid, p, lon, lat)
        CT = gsw.conversions.CT_from_pt(SA, pt_grid)
        return gsw.density.sigma3(SA, CT)

    plt.figure(figsize=(8,6))

    data = data[(data['DEPTH']>min_dep) & (data['DEPTH']<max_dep)]

    # create a grid of temperatures and salinities
    pt_grid = np.linspace(data['potemperature'].min()-0.005, data['potemperature'].max()+0.005, 100)
    psal_grid = np.linspace(data['PSAL'].min()-0.0005, data['PSAL'].max()+0.0005, 100)
    pt_grid, psal_grid = np.meshgrid(pt_grid, psal_grid)

    density = calculate_density(pt_grid, psal_grid, p_ref, lon, lat)

    #plt.figure(figsize=(8, 6))
    # add constant density lines
    contours = plt.contour(psal_grid, pt_grid, density, levels=np.arange(np.min(density), np.max(density), (np.max(density)-np.min(density))/20), colors='black')       #for all 10, for zoom 20
    plt.clabel(contours, inline=True, fontsize=10)

    #iterate through stations
    if c_var == 'Station':
        station_list = [d for _, d in data.groupby(['Station'])]

        #iterate colors
        ncolors = len(station_list)
        colors = plt.cm.jet(np.linspace(0,1,ncolors))# Initialize holder for trajectories

        for i, station in enumerate(station_list):
                # create the T-S Diagram
                plt.scatter(station['PSAL'], station['potemperature'],s=20, label=station['Station'].iloc[0], color=colors[i])

        plt.legend(markerscale=2, ncol=2, loc='lower left')  #for markerscale=2, else =5

        #format x axis without exponetial
        from matplotlib.ticker import ScalarFormatter
        ax = plt.gca()
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)

    #iterate through casts
    elif c_var == 'Cast':
        station_list = [d for _, d in data.groupby(['Station'])]
        for i, station in enumerate(station_list):
            station_list[i] = add_castno(station_list[i])
        all_casts = pd.concat(station_list)
        all_casts = [d for _, d in all_casts.groupby(['Station','Cast'])]
        print(len(all_casts))

        #iterate colors
        ncolors = len(all_casts)
        colors = plt.cm.jet(np.linspace(0,1,ncolors))# Initialize holder for trajectories

        for i, cast in enumerate(all_casts):
                # create the T-S Diagram
                plt.scatter(cast['PSAL'], cast['potemperature'],s=5, label=cast['Station'].iloc[0], color=colors[i])

        plt.legend(markerscale=5, ncol=2,loc='lower left')

        #format x axis without exponetial
        from matplotlib.ticker import ScalarFormatter
        ax = plt.gca()
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)

    else:
        # create the T-S Diagram
        plt.scatter(data['PSAL'], data['potemperature'], c=data[c_var], cmap='viridis_r',s=20)

        cbar = plt.colorbar()
        cbar.set_label(get_var(c_var)[0])

        #format x axis without exponetial
        from matplotlib.ticker import ScalarFormatter
        ax = plt.gca()
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)

    #plt.colorbar(label='Potential Density Anomaly (kg/m³)')
    plt.ylabel(get_var('potemperature')[0])
    plt.xlabel(get_var('PSAL')[0])

    plt.tight_layout()
    plt.show()

def plot_ts(data, c_var, min_dep, max_dep, p_ref, lon, lat):
    """
    This function plots data of a CTD cast as a T-S diagram with potential temperature and practical salinity. The samples can be filtered by a depth range. For the isopycnal calculation a reference density, longitude and latitude is used.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with CTD data containing the potential temperature and practical salinity.
    c_var : string
        Variable which should be used as third dimension with colorbar.
    min_dep : int
        Minimal depth of the samples used.
    max_dep : int
        Maximal depth of the samples used.
    p_ref : int
        Reference pressure used for calculating the isopycnals of the potential density anomaly in dbar.
    lon : float
        Longitude of the samples.
    lat : float
        Latitude of the samples.

    Returns
    -------
    None.

    """
    import matplotlib.pyplot as plt
    import numpy as np
    import gsw
    from hyvent.misc import get_var, add_castno
    import pandas as pd

    # calculate the density values for the grid
    def calculate_density(pt_grid, psal_grid, p, lon, lat):
        SA = gsw.conversions.SA_from_SP(psal_grid, p, lon, lat)
        CT = gsw.conversions.CT_from_pt(SA, pt_grid)
        return gsw.density.sigma3(SA, CT)

    plt.figure(figsize=(8,6))

    data = data[(data['DEPTH']>min_dep) & (data['DEPTH']<max_dep)]

    # create a grid of temperatures and salinities
    pt_grid = np.linspace(data['potemperature'].min()-0.005, data['potemperature'].max()+0.005, 100)
    psal_grid = np.linspace(data['PSAL'].min()-0.0005, data['PSAL'].max()+0.0005, 100)
    pt_grid, psal_grid = np.meshgrid(pt_grid, psal_grid)

    density = calculate_density(pt_grid, psal_grid, p_ref, lon, lat)

    #plt.figure(figsize=(8, 6))
    # add constant density lines
    contours = plt.contour(psal_grid, pt_grid, density, levels=np.arange(np.min(density), np.max(density), (np.max(density)-np.min(density))/10), colors='black')
    plt.clabel(contours, inline=True, fontsize=10)

    #iterate through stations
    if c_var == 'Station':
        station_list = [d for _, d in data.groupby(['Station'])]

        #iterate colors
        ncolors = len(station_list)
        colors = plt.cm.jet(np.linspace(0,1,ncolors))# Initialize holder for trajectories

        for i, station in enumerate(station_list):
                # create the T-S Diagram
                plt.scatter(station['PSAL'], station['potemperature'],s=5, label=station['Station'].iloc[0], color=colors[i])

        plt.legend(markerscale=5, ncol=2, loc='lower left')

        #format x axis without exponetial
        from matplotlib.ticker import ScalarFormatter
        ax = plt.gca()
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)

    #iterate through casts
    elif c_var == 'Cast':
        station_list = [d for _, d in data.groupby(['Station'])]
        for i, station in enumerate(station_list):
            station_list[i] = add_castno(station_list[i])
        all_casts = pd.concat(station_list)
        all_casts = [d for _, d in all_casts.groupby(['Station','Cast'])]
        print(len(all_casts))

        #iterate colors
        ncolors = len(all_casts)
        colors = plt.cm.jet(np.linspace(0,1,ncolors))# Initialize holder for trajectories

        for i, cast in enumerate(all_casts):
                # create the T-S Diagram
                plt.scatter(cast['PSAL'], cast['potemperature'],s=5, label=cast['Station'].iloc[0], color=colors[i])

        plt.legend(markerscale=5, ncol=2,loc='lower left')

        #format x axis without exponetial
        from matplotlib.ticker import ScalarFormatter
        ax = plt.gca()
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)

    else:
        # create the T-S Diagram
        plt.scatter(data['PSAL'], data['potemperature'], c=data[c_var], cmap='viridis_r',s=5)

        cbar = plt.colorbar()
        cbar.set_label(get_var(c_var)[0])

        #format x axis without exponetial
        from matplotlib.ticker import ScalarFormatter
        ax = plt.gca()
        formatter = ScalarFormatter(useOffset=False)
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)

    #plt.colorbar(label='Potential Density Anomaly (kg/m³)')
    plt.ylabel(get_var('potemperature')[0])
    plt.xlabel(get_var('PSAL')[0])

    plt.tight_layout()
    plt.show()

def plot3Ddistr(data, var, depth_min, bathy=None ,vent_loc=None):
    """
    This functions plots a variable of a dataset of one or a list of stations in 3D space together with the bathymetry.

    Parameters
    ----------
    data : pandas dataframe or list of pandas dataframes
        Dataframes containing the data with 3D coordinates.
    var : string
        Key of dataframe column to plot.
    depth_min : int
        Minimal cutoff depth to plot.
    bathy : tuple of three int, OPTIONAL
        Bathymetry data with longitude and latitude coordinates and elevation. Default is None.
    vent_loc : tuple of three float, OPTIONAL
        Longitude, latitude and depth of a position of interest to mark in the plot. Default is None.

    Returns
    -------
    None.

    """
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    from hyvent.misc import get_var

    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(projection='3d')

    #prepare and plot data for one dataframe
    if isinstance(data, pd.DataFrame):
        #subset depth
        dataset =data[data.DEPTH>depth_min]
        #mask negative anomalies with nan, because not interesting
        if var =='dORP':
            dataset = dataset.mask(dataset[var]>0, np.nan)
        else:
            dataset = dataset.mask(dataset[var]<0, np.nan)

        #plot only nth point
        n=1
        dataset = dataset.iloc[::n,:]

        data = ax.scatter(dataset['CTD_lon'], dataset['CTD_lat'],dataset['DEPTH'], c=dataset[var],marker='.',alpha=1,cmap=get_var(var)[2],linewidth=0.2)

    #prepare and plot data for a list of dataframes
    if isinstance(data, list):
        for station in data:
            #subset depth
            dataset =station[station.DEPTH>depth_min]
            #mask negative anomalies with nan, because not interesting
            if var =='dORP':
                dataset = dataset.mask(dataset[var]>0, np.nan)
            else:
                dataset = dataset.mask(dataset[var]<0, np.nan)
            #plot only nth point
            n=1
            dataset = dataset.iloc[::n,:]

            data = ax.scatter(dataset['CTD_lon'], dataset['CTD_lat'],dataset['DEPTH'], c=dataset[var],marker='.',alpha=1,cmap=get_var(var)[2],linewidth=0.2)

    #plot bathymetry
    if bathy != None:
        surf = ax.plot_surface(bathy[0], bathy[1], -bathy[2],alpha=0.3,linewidth=0.5,color='white',antialiased=False)
    #plot vent location
    if vent_loc != None:
        vent = ax.scatter(vent_loc[0], vent_loc[1], vent_loc[2], marker='*', color='r', s=300)

    ax.invert_zaxis()
    fig.colorbar(data, label=get_var(var)[0],fraction=0.02)

    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')
    ax.set_zlabel('Depth in m')

    ax.xaxis.labelpad=10
    ax.yaxis.labelpad=10
    ax.zaxis.labelpad=10

    ax.view_init(8, -173)

    plt.show()

def plot_contourf(data, var, xvar, depth_min, vent_loc=None):
    """
    This funtion plots interpolated crossections along an x coordinate and depth.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the CTD data in two dimensions.
    var : string
        Key of dataframe column to plot as the interpolated data.
    xvar : string
        Key of dataframe column to use as an x coordinate. Must be latitude or longitude.
    depth_min : int
        Minimal cutoff depth to plot.
    vent_loc : tuple of three float, OPTIONAL
        Longitude, latitude and depth of a position of interest to mark in the plot. Default is None..

    Returns
    -------
    None.

    """

    import matplotlib.pyplot as plt
    from hyvent.misc import get_var

    #subset data by depth and remove nan in var
    data = data[data['DEPTH']>depth_min]
    data = data.dropna(subset=[var,xvar])

    #set vmin and vmax for colobar to ignore not interesting anomaly values
    if var =='dORP':
        vmin = data[var].min()
        vmax = 0
    else:
        vmin = 0
        vmax = data[var].max()

    fig,ax = plt.subplots()

    #plot data
    contourf = ax.tricontourf(data[xvar],data['DEPTH'],data[var],levels=50,cmap=get_var(var)[2],vmin=vmin,vmax=vmax)
    #plot sample locations
    ax.scatter(data[xvar],data['DEPTH'],color='black',s=0.2,alpha=0.2)

    ax.set_ylabel(get_var('DEPTH')[0])
    if xvar == 'CTD_lat':
        ax.set_xlabel('Latitude')
    if xvar == 'CTD_lon':
        ax.set_xlable('Longitude')
    fig.colorbar(contourf, label=get_var(var)[0])

    #plot vent location
    if vent_loc != None:
        if xvar == 'CTD_lon':
            vent = ax.scatter(vent_loc[0], vent_loc[2], marker='*', color='r', s=300)
        if xvar == 'CTD_lat':
            vent = ax.scatter(vent_loc[1], vent_loc[2], marker='*', color='r', s=300)

    ax.invert_yaxis()

    fig.set_tight_layout(True)
    plt.show()


