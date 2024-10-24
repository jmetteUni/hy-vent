#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 13:30:29 2024

@author: jmette@uni-bremen.de
"""

def keys_to_data(data,datatype):
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

    Returns
    -------
    data : dictionary
        Updated dictionary with additional columns in the pandas dataframes.

    """

    if datatype == 'cnv':
        for key in data:
            new_key = str.split(key,'_')
            cruise = new_key[0]
            station = new_key[1]+'_'+new_key[2]

            data[key]['Cruise'] = cruise
            data[key]['Station'] = station
            data[key]['SN'] = 'SBE9'
            data[key]['Operation'] = 'profile'

    if datatype == 'btl':
        for key in data:
            new_key = str.split(key,'_')
            cruise = new_key[0]
            station = new_key[1]+'_'+new_key[2]

            data[key]['Cruise'] = cruise
            data[key]['Station'] = station
            data[key]['SN'] = 'SBE32'
            data[key]['Operation'] = 'water sample'

    if datatype == 'mapr':
        for key in data:
            new_key = str.split(key,'_')
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
        color = 'green'
        label = '$dE/dt$ in mV/sec'
        cmap = 'Greens_r'
        cmap = truncate_colormap(cmap)
    elif var == 'upoly0':
        color = 'green'
        label = '$E$ in mV'
        cmap = 'Greens_r'
#salinity
    elif var == 'PSAL' or var == 'PSAL_mean':
        color = 'cyan'
        label = 'Practical Salinity in 1'
        cmap = 'None'
    elif var == 'SA':
        color = 'cyan'
        label = 'Absolute Salinity in ?'
        cmap = 'None'
#density
    elif var == 'Sigma0':
        color = 'purple'
        label = '$\sigma_0$ in kg/m$^3$'
        cmap = 'Purples_r'
    elif var == 'Sigma3':
        color = 'purple'
        label = r'$\sigma_{3000}$ in kg/m$^3$'
        cmap = 'Purples_r'
    elif var == 'Delta_Sigma3':
        color = 'purple'
        label = r'$\Delta$$\sigma_{3000}$ in kg/m$^3$'
        cmap = 'Purples_r'
    elif var == 'Rho':
        color = 'purple'
        label = '$\rho$ in kg/m$^3$'
        cmap = 'Purples_r'
    elif var == 'delta_podensity':
        color = 'purple'
        label = 'Potential Density Anomaly in kg/m$^3$'
        cmap = 'Purples_r'
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
    else:
        print('Properties for variable '+var+' not defined yet.')

    return label,color,cmap


