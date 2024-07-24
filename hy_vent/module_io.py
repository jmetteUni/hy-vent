#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 17:04:38 2024

@author: jmett@uni-bremen.de
"""

def read_cnv(path):
    import pandas as pd
    from seabird.cnv import fCNV
    import os
    import datetime as dt

    file_list = []                  #get all cnv filenames
    for file in os.listdir(path):
        if file.endswith(".cnv"):
           file_list.append(file)
           file_list.sort()

    cnv_data = dict()           #read cnv file into pandas dataframe with fCNV package
    for file in file_list:
        key = file.rstrip('.cnv')
        cnv_data[key] = fCNV(path+file).as_DataFrame()
        cnv_data[key]['basedate'] = dt.datetime(2000,1,1)      #create datetime object column
        cnv_data[key]['datetime'] = cnv_data[key]['basedate']+pd.to_timedelta(cnv_data[key]['timeQ'],unit='seconds')
        del cnv_data[key]['basedate']
        print('Read '+file+' sucessfully')

    return cnv_data

def write_to_csv(cnv_data,out_path,datatype):    #writes all stations to individual csv sheets
    import pandas as pd
    for key in cnv_data:
        station = pd.DataFrame(cnv_data[key])

        station.to_csv(out_path+key+'_'+datatype+'.csv',index=False)
        print('Wrote '+key+' to csv sucessfully')

def read_from_csv(csv_path,datatype):        #reads in stations in individual csv sheets into a dictionary
    import pandas as pd
    import os

    file_list = []                  #get all cnv filenames
    for file in os.listdir(csv_path):
        if datatype == 'cnv':
            if file.endswith("cnv.csv"):
               file_list.append(file)
               file_list.sort()
        if datatype == 'btl':
            if file.endswith("btl.csv"):
               file_list.append(file)
               file_list.sort()

    cnv_data = dict()
    for file in file_list:
        if file.endswith("cnv.csv"):
            cnv_data[file.rstrip('_cnv.csv')] = pd.read_csv(csv_path+file)       #,encoding='iso-8859-1' sometimes needed
        if file.endswith("btl.csv"):
            cnv_data[file.rstrip('_btl.csv')] = pd.read_csv(csv_path+file)       #,encoding='iso-8859-1' sometimes needed
        print('Read '+file+' sucessfully')

    return cnv_data

def read_dshippos(dship_path):
    import pandas as pd

    dship = pd.read_fwf(dship_path,skiprows=3,names=['Date','Time','Dship_lat','Dship_lon'],encoding='iso-8859-1')
    dship['datetime'] = pd.to_datetime(dship['Date']+' '+dship['Time'])
    del dship['Date']
    del dship['Time']

    return dship

def read_btl(btl_path):
    import pandas as pd
    import os

    btl_data = dict()

    file_list = []                  #get all filenames of the positiondata csv sheets
    for file in os.listdir(btl_path):
        if file.endswith(".btl"):
           file_list.append(file)
           file_list.sort()

    for file in file_list:
        f = open(btl_path+file,'r',encoding='latin-1')
        first_data_line = 0
        while True:
            line = f.readline()
            if line[0] != '#' and line[0] != '*':
                f.close()
                break
            first_data_line = first_data_line+1
        df = pd.read_fwf(btl_path+file, skiprows=first_data_line,)
        df = df[1:]
        n_bottles = int(len(df)/4)
        columns = list(df.columns)
        columns.remove('Date')
        columns.remove('Bottle')
        columns.remove('Latitude')
        columns.remove('Longitude')
        new_columns = []
        for column in columns:
            new_columns.append(column+'_mean')
            new_columns.append(column+'_sd')
            new_columns.append(column+'_min')
            new_columns.append(column+'_max')

        df_output = pd.DataFrame(columns=['Bottle','Date','Time','Latitude','Longitude']+new_columns,index=range(1,n_bottles+1))

        for i in range(n_bottles):
            index = i*4
            date = df['Date'][index+1]
            time = df['Date'][index+2]
            lat = df['Latitude'][index+1]
            long = df['Longitude'][index+1]
            bottle = df['Bottle'][index+1]
            for column in columns:
                df_output[column+'_mean'][i+1] = df[column][index+1]
                df_output[column+'_sd'][i+1] = df[column][index+2]
                df_output[column+'_min'][i+1] = df[column][index+3]
                df_output[column+'_max'][i+1] = df[column][index+4]
                df_output['Date'][i+1] = date
                df_output['Time'][i+1] = time
                df_output['Latitude'][i+1] = lat
                df_output['Longitude'][i+1] = long
                df_output['Bottle'][i+1] = bottle

        df_output= df_output.iloc[: , :-4]
        df_output['datetime'] = pd.to_datetime(df_output['Date']+' '+df_output['Time'],dayfirst=True)      #creates datetime object column
        df_output['Bottle'] = df_output['Bottle'].astype(int)
        file = file.rstrip('.btl')
        btl_data[file] = df_output
        print('Read '+file+' successfully!')

    return btl_data