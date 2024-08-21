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
