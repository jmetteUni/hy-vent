#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 13:30:29 2024

@author: jmette@uni-bremen.de
"""

def keys_to_data(data,datatype,project_name=None):
    """

    Adds meta information which is stored in the key of the "data" dictionary
    as column values for every row. The use of this function is recommended if
    you want to use all station as one dataframe and not organized in a
    dictionary

    Parameters
    ----------
    data : dictionary
        Dictionary of pandas dataframes. The keys are unique station
        identifiers similar to "PS137_018_01" with "PS137" as the cruise name,
        "018" as the 3-digit station number and "01" as the two digit cast
        number. The values are pandas dataframes with the data organised as
        column variables.
        For NOAA PMEL MAPR data the format is "PS137_NUI57-1_72" where "NUI"
        describes the type of operation and "72" is the device's serial number.
    datatype : string
        String which describes the type of data. Can be "cnv" for full CTD
        profiles, "btl" for SeaBird bottle file data or "mapr" for NOAA PMEL
        MAPR data.
    project_name : string, optional
        Designates the project name. If this is not 'None', it is assumed that the project name is not included in the keys of the data sets. Therefore it is added in the beginning of each key. The default is 'None'.

    Returns
    -------
    data : dictionary
        Updated dictionary with additional columns in the pandas dataframes.

    """

    if datatype == 'cnv':
        for key in data:
            #check and eventually add project name
            if project_name != None:
                new_key = project_name+'_'+key
            else:
                new_key = key
            new_key = new_key.replace('-','_')
            new_key = str.split(new_key,'_')
            cruise = new_key[0]
            station = new_key[1]+'_'+new_key[2]

            data[key]['Cruise'] = cruise
            data[key]['Station'] = station
            data[key]['SN'] = 'SBE9'
            data[key]['Operation'] = 'profile'

    if datatype == 'btl':
        for key in data:
            #check and eventually add project name
            if project_name != None:
                new_key = project_name+'_'+key
            else:
                new_key = key
            new_key = new_key.replace('-','_')
            new_key = str.split(new_key,'_')
            cruise = new_key[0]
            station = new_key[1]+'_'+new_key[2]

            data[key]['Cruise'] = cruise
            data[key]['Station'] = station
            data[key]['SN'] = 'SBE32'
            data[key]['Operation'] = 'water sample'

    if datatype == 'mapr':
        for key in data:
            #check and eventually add project name
            if project_name != None:
                new_key = project_name+'_'+key
            else:
                new_key = key
            new_key = new_key.replace('-','_')
            new_key = str.split(new_key,'_')
            cruise = new_key[0]
            sn = new_key[2]
            operation = new_key[1][:-4]
            station = new_key[1][-4:]
            station = str.split(station,'-')
            station = '0'+station[0]+'_0'+station[1]

            data[key]['Cruise'] = cruise
            data[key]['Station'] = station
            data[key]['SN'] = sn
            data[key]['Operation'] = operation

    return data

def add_castno(data, window_size=1000):
    """
    This funcion adds to each measurement a cast identification number, increasing for all up and down cast which are found in the data per station, ignoring station. It identifies casts by local extrema using the function hyvent.processing.sep_casts(). Check, if the window size in sep_casts is fitting.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe containing data with several up and downcasts. Must contain the column "Station".
    window_size : int, optional
        Size of the window use by the sep_casts() function. Default is 1000.

    Returns
    -------
    data_cast : pandas dataframe
        Dataframe similar as the input, with an additonal column containing the cast number.
    """
    from hyvent.processing import sep_casts
    import pandas as pd

    if isinstance(data, pd.DataFrame):
        casts = []
        data_list = [d for _, d in data.groupby(['Station','SN'])]
        for data in data_list:
            cast_no = 1
            casts_list = sep_casts(data, window_size)
            # import matplotlib.pyplot as plt   #used for controlling number of casts
            # plt.figure()
            for cast in casts_list:
                # plt.plot(cast['DEPTH'])
                cast['Cast'] = cast_no
                cast_no = cast_no +1
            casts.append(pd.concat(casts_list))
        data_cast = pd.concat(casts)

        return data_cast

def get_var(var):
    """
    Returns porperties such as plot color and axis label for a given variable.

    Parameters
    ----------
    var : string
        Variable, for which properies should be returned. Must be defined in the function.

    Returns
    -------
    label : string
        Axis label for the variable.
    color : string
        Color to be used for plotting the variable.
    """
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np

    def truncate_colormap(cmap, minval=0, maxval=0.8, n=100):
        cmap = plt.get_cmap(cmap)
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

#temperature
    if var == 'potemperature':
        color = 'red'
        label = '$\Theta$ in $^{\circ}$C'
        cmap = 'Reds'
    elif var == 'CT':
        color = 'red'
        label = 'Conservative Temperature in $^{\circ}$C'
        cmap = 'Reds'
    elif var == 'TEMP':
        color = 'red'
        label = 'Temperature in $^{\circ}$C'
        cmap = 'Reds'
    elif var == 'Delta_potemperature':
        color = 'red'
        label = '$\Delta$$\Theta$ in $^{\circ}$C'
        cmap = 'Reds'
#turbidity
    elif var == 'Neph(volts)' or var == 'Neph_outl(volts)' or var == 'Neph_smoo(volts)':
        color = 'blue'
        label = 'Turbidity in NTU'
        cmap = 'Blues'
    elif var == 'Delta_Neph(volts)' or var == 'Delta_Neph_outl(volts)' or var == 'Delta_Neph_smoo(volts)':
        color = 'blue'
        label = '$\Delta$ Turbidity in NTU'
        cmap = 'Blues'
#ORP
    elif var == 'dORP':
        color = 'tab:purple'
        label = '$dE_h/dt$ in mV s$^{-1}$'
        cmap = 'Purples_r'
        cmap = truncate_colormap(cmap)
    elif var == 'upoly0':
        color = 'tab:purple'
        label = '$E$ in mV'
        cmap = 'Purples_r'
#salinity
    elif var == 'PSAL' or var == 'PSAL_mean':
        color = 'cyan'
        label = 'Practical Salinity in 1'
        cmap = 'None'
    elif var == 'SA':
        color = 'cyan'
        label = 'Absolute Salinity in ?'
        cmap = 'None'
    elif var == 'Delta_PSAL':
        color = 'cyan'
        label = 'Practical Salinity in 1'
        cmap = 'None'
#density
    elif var == 'Sigma0':
        color = 'tab:olive'
        label = '$\sigma_0$ in kg m$^{-3}$'
        cmap = 'YlGn_r'
    elif var == 'Sigma3':
        color = 'tab:olive'
        label = r'$\sigma_{3000}$ in kg m$^{-3}$'
        cmap = 'YlGn_r'
    elif var == 'Delta_Sigma3':
        color = 'tab:olive'
        label = r'$\Delta$$\sigma_{3000}$ in kg m$^{-3}$'
        cmap = 'YlGn_r'
    elif var == 'Rho' or var == 'density':
        color = 'tab:olive'
        label = '$\rho$ in kg m$^{-3}$'
        cmap = 'YlGn_r'
    elif var == 'delta_podensity':
        color = 'tab:olive'
        label = 'Potential Density Anomaly in kg m$^{-3}$'
        cmap = 'YlGn_r'
#depth
    elif var == 'DEPTH' or var == 'Depth_corr(m)' or var == 'DepSM_mean':
        color = 'black'
        label = 'Depth in m'
        cmap = 'None'
#pressure
    elif var == 'PRES':
        color = 'black'
        label = 'Pressure in dbar'
        cmap = 'None'
#helium
    elif var == 'delta3He':
        color = 'orange'
        label = '$\delta^3$He in %'
        cmap = 'Oranges'
        cmap = truncate_colormap(cmap, minval=0.2, maxval=1)
    elif var == 'Delta_delta3He':
        color = 'orange'
        label = '$\Delta$ $\delta^3$He in %'
        cmap = 'Oranges'
        cmap = truncate_colormap(cmap, minval=0.2, maxval=1)
    elif var == 'Dist_vent':
        color = 'black'
        label = 'Distance to vent in m'
        cmap = 'None'

    else:
        print('Properties for variable '+var+' not defined yet.')

    return label,color,cmap

def dist(lon1, lat1, lon2, lat2, dep1='None', dep2='None'):
    """
    This function calculates the geodesic distance between two geographical positions. If depth values is also given, it calculates the 3D geometric distance instead. The second method is only recommended for shorter distances.

    Parameters
    ----------
    lon1 : float
        Longitude of first position.
    lat1 : float
        Latitude of first position.
    lon2 : float
        Longitude of second position.
    lat2 : float
        Latitude of second position.
    dep1 : float, optional
        Depth (positive) of first position. Default is None.
    dep2 : float, optional
        Depth (positive) of second position. Default is None.

    Returns
    -------
    d : float
        The calculated distance.

    """
    #from geopy.distance import geodesic
    from geographiclib.geodesic import Geodesic
    import pandas as pd
    import numpy as np

    if pd.isna(lat1) == True or pd.isna(lon1) == True or pd.isna(lat2) == True or pd.isna(lon2) == True or pd.isna(dep1) == True or pd.isna(dep2):
        d = np.nan
    else:
        #d = geodesic((lat1, lon1), (lat2, lon2)).m
        if dep1 !='None':
            d = np.sqrt((lat2+lat1)**2+(lon2+lon1)**2+(dep2-dep1)**2)
        else:
            d = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)['s12']
    return(d)

def bearing(lon1, lat1, lon2, lat2):
    """
    This function calculates the bearing between two geographical positions, from the view of the first point.

    Parameters
    ----------
    lon1 : float
        Longitude of first position.
    lat1 : float
        Latitude of first position.
    lon2 : float
        Longitude of second position.
    lat2 : float
        Latitude of second position.

    Returns
    -------
    bearing : float
        The calculated bearing.
    """
    from geographiclib.geodesic import Geodesic
    import pandas as pd

    import numpy as np

    if pd.isna(lat1) == True or pd.isna(lon1) == True or pd.isna(lat2) == True or pd.isna(lon2) == True:
        bearing = np.nan
    else:
        bearing = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)['azi2']
    return(bearing)

def calc_delta_pos(data):
    """
    This function calculates the differences in distance and bearing between lon/lat position in one row and the next row in a dataframe.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with consectuive position data with longitude ("CTD_lon") and latitude ("CTD_lat") as one column each.

    Returns
    -------
    data : pandas dataframe
        Dataframe containing the longitude, latitude, distance and bearing between position in this row and the next row.
    """

    data = data[['CTD_lon','CTD_lat']]
    data['Lat_next'] = data['CTD_lat'].shift()
    data['Lon_next'] = data['CTD_lon'].shift()
    data['dist'] = data.apply(lambda x: dist(x['CTD_lon'], x['CTD_lat'], x['Lon_next'], x['Lat_next']), axis=1)
    data['bearing'] = data.apply(lambda x: bearing(x['CTD_lon'], x['CTD_lat'], x['Lon_next'], x['Lat_next']), axis=1)

    return data[['CTD_lon','CTD_lat','dist','bearing']]

def calc_vel_from_posi(data):
    """
    This functions calculates the velocity between each data point of a moving CTD (or another device) based on lon/lat position data.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe containing the longitude ("CTD_lon") and the latitude ("CTD_lat") of the CTD. Must also contain the column "datetime".

    Returns
    -------
    data : pandas dataframe
        Dataframe similar to the input but with u and v components and the complex velocity each as additional column.
    """
    from hyvent.misc import calc_delta_pos
    import pandas as pd
    import numpy as np

    #make sure, "datetime" is a datetime object
    data['datetime'] = pd.to_datetime(data['datetime'])

    #calculate the time difference in seconds
    data['timediff'] = (data['datetime'] - data['datetime'].shift()).dt.total_seconds()
    #calculate the difference in distance between measurements
    delta = calc_delta_pos(data)

    #calculate complex velocities between measurements
    data['u_mov'] = delta['dist'] / data['timediff'] *np.cos(np.radians(delta['bearing']))
    data['v_mov'] = delta['dist'] / data['timediff'] *np.sin(np.radians(delta['bearing']))
    data['cv_mov'] = data['u_mov'] + 1j* data['v_mov']

    return data


