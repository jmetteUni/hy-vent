#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 13:30:29 2024

@author: jmette@uni-bremen.de
"""

def keys_to_data(data,datatype):

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
