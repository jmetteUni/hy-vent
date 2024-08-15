#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:25:58 2024

@author: jonathan
"""

def plot_hist_he(data,depth_max,depth_min,bins,path_save='None'):
    import pandas as pd
    import matplotlib.pyplot as plt

    label_sel = 'Samples between '+str(depth_min)+' - '+str(depth_max)+'m'
    he_data = pd.concat(data.values(),ignore_index=True)

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
    import pandas as pd
    import matplotlib.pyplot as plt

    label_sel = 'Samples between '+str(depth_min)+' - '+str(depth_max)+'m'
    dorp_data = pd.concat(data.values(),ignore_index=True)

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
    from matplotlib import colors
    import matplotlib.pyplot as plt

    label_title = 'Interpolated crossection for '+key

    station = profile_data[key]
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
    import matplotlib.pyplot as plt
    from matplotlib import colors

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

def plot2D_all_stations_he(btl_data,var,vent_loc,depth_min,depth_max,bathy=False,path_save='None'):
    import matplotlib.pyplot as plt

    label_title = '2D Distribution in '+str(depth_min)+' - '+str(depth_max)+' m'

    if bathy != False:
        lons = bathy[0]
        lats = bathy[1]
        elev = bathy[2]

    fig, ax = plt.subplots()
    fig.set_tight_layout(True)

    if bathy != False:
        contours = plt.contourf(lons,lats,elev,levels=40,alpha=0.3)
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








