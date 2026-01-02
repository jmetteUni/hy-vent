#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 18:10:45 2024

@author: jmette@uni-bremen.de
"""


def read_cnv(path,suffix=None):
    """
    Reads in data from SeaBird CTD processing software which ends with ".cnv"
    to a dictionary with the filename as the key and a pandas dataframe per
    station as the value. It exspects one file per station in the given
    directory.This function uses the seabird package by castelao: https://github.com/castelao/seabird

    Parameters
    ----------
    path : string
        Path to the directory where the data is stored in one station per file.

    Returns
    -------
    cnv_data : dictionary
        Dictionary where the keys are unique station identifiers and
        values are data in pandas dataframes.
    suffix : string, optional
        Suffix at the end of the filename, to filter the files which are read in. Default is None

    """
    import pandas as pd
    from seabird.cnv import fCNV
    import os
    import datetime as dt

    file_list = []  # get all cnv filenames
    for file in os.listdir(path):
        if suffix != None:
            if file.endswith(suffix+".cnv"):
                file_list.append(file)
                file_list.sort()
        else:
            if file.endswith(".cnv"):
                file_list.append(file)
                file_list.sort()

    cnv_data = dict()  # read cnv file into pandas dataframe with fCNV package
    for file in file_list:
        key = file.rstrip('.cnv')
        if suffix in key:
            key = key.replace(suffix, '')
        try:
            cnv_data[key] = fCNV(os.path.join(path,file)).as_DataFrame()
            cnv_data[key]['basedate'] = dt.datetime(
                2000, 1, 1)  # create datetime object column
            if cnv_data[key]['timeQ'].isna().any() == True:
                cnv_data[key]['timeQ'] = cnv_data[key]['timeQ'].interpolate(
                    axis=0)  # interpolates nan values
                if cnv_data[key]['timeQ'].duplicated().any() == True:
                    # average rows with the same timestamp, NaNs already filled
                    cnv_data[key] = cnv_data[key].groupby(
                        'timeQ').mean().reset_index()

            cnv_data[key]['datetime'] = cnv_data[key]['basedate'] + \
                pd.to_timedelta(cnv_data[key]['timeQ'], unit='seconds')
            del cnv_data[key]['basedate']
            print('Read '+file+' sucessfully')
        except:
            print('Error reading '+file+', skipping...')

    # infer datetypes
    for key in cnv_data:
        cnv_data[key] = cnv_data[key].infer_objects(copy=False)

    return cnv_data


# writes all stations to individual csv sheets
def write_to_csv(data_in, out_path, datatype):
    """
    Writes out data from a dictionary with the filename as the key and a pandas
    dataframe per station as the value to csv files, where the key is the
    filename and the dataframe is the csv files content. If the input data is one dataframe with several stations, the filename is created from column data.


    Parameters
    ----------
    cnv_data : dictionary or pandas dataframe
        Dictionary where the keys are unique station identifiers and
        values are data in pandas dataframes or one dataframe

    out_path : string
        Path to the directory to write to.

    datatype : string
        String which describes the type of data. Can be "cnv" for full CTD
        profiles, "btl" for SeaBird bottle file data or "mapr" for NOAA PMEL
        MAPR data.

    Returns
    -------
    None

    """

    import pandas as pd
    import os

    # if data is in one dataframe
    if isinstance(data_in, pd.DataFrame):
        # differentiate between datatype to create the the corresponding filenames
        if datatype == 'mapr':
            data_list = [d for _, d in data_in.groupby(['Station', 'SN'])]
            for data in data_list:
                filename = data['Cruise'].iloc[0]+'_' + \
                    data['Station'].iloc[0]+'_' + \
                    data['SN'].iloc[0]+'_'+datatype+'.csv'
                write_path = os.path.join(out_path, filename)
                data.to_csv(write_path, index=False)
                print('Wrote '+filename+' sucessfully')
        else:
            data_list = [d for _, d in data_in.groupby(['Station'])]
            for data in data_list:
                filename = data['Cruise'].iloc[0]+'_' + \
                    data['Station'].iloc[0]+'_'+datatype+'.csv'
                write_path = os.path.join(out_path, filename)
                data.to_csv(write_path, index=False)
                print('Wrote '+filename+' sucessfully')

    # if data is a dict of dataframes
    elif isinstance(data_in, dict):
        for key in data_in:
            station = pd.DataFrame(data_in[key])
            filename = key+'_'+datatype+'.csv'
            write_path = os.path.join(out_path, filename)
            station.to_csv(write_path, index=False)
            print('Wrote '+filename+' sucessfully')


# reads in stations in individual csv sheets into a dictionary
def read_from_csv(csv_path, datatype):
    """
    Reads data from csv files to a dictionary with the filename as the key and
    a pandas dataframe per station as the value. It exspects one file per
    station in the given directory.


    Parameters
    ----------
    csv_path : string
        Path to the directory where the csv files are stored.

    datatype : string
        String which describes the type of data. Can be "cnv" for full CTD
        profiles, "btl" for SeaBird bottle file data or "mapr" for NOAA PMEL
        MAPR data.

    Returns
    -------
    cnv_data : dictionary
        Dictionary where the keys are unique station identifiers and
        values are data in pandas dataframes.

    """

    import pandas as pd
    import os

    file_list = []  # get all cnv filenames
    for file in os.listdir(csv_path):
        if datatype == 'cnv':
            if file.endswith("cnv.csv"):
                file_list.append(file)
                file_list.sort()
        if datatype == 'btl':
            if file.endswith("btl.csv"):
                file_list.append(file)
                file_list.sort()
        if datatype == 'mapr':
            if file.endswith("mapr.csv"):
                file_list.append(file)
                file_list.sort()

    cnv_data = dict()
    for file in file_list:
        if file.endswith("cnv.csv"):
            in_path = os.path.join(csv_path,file)
            cnv_data[file.rstrip('_cnv.csv')] = pd.read_csv(
                in_path)  # ,encoding='iso-8859-1' sometimes needed
        if file.endswith("btl.csv"):
            in_path = os.path.join(csv_path,file)
            cnv_data[file.rstrip('_btl.csv')] = pd.read_csv(
                in_path)  # ,encoding='iso-8859-1' sometimes needed
        if file.endswith("mapr.csv"):
            in_path = os.path.join(csv_path,file)
            cnv_data[file.rstrip('_mapr.csv')] = pd.read_csv(
                in_path)  # ,encoding='iso-8859-1' sometimes needed
        print('Read '+file+' sucessfully')

    for key in cnv_data:
        cnv_data[key]['datetime'] = pd.to_datetime(
            cnv_data[key]['datetime'], format='mixed')

    return cnv_data


def read_dshippos(dship_path):
    """
    Reads in tab separated file with data from "D-SHIP" portal. Expects 3 rows
    of headers and 4 columns with date, time, latitude and longitude.

    Parameters
    ----------
    dship_path : string
        Path to the file.

    Returns
    -------
    dship : pandas dataframe
        Pandas Dataframe with date, time, latitude and longitude as columns.

    """

    import pandas as pd

    dship = pd.read_fwf(dship_path, skiprows=3, names=[
                        'Date', 'Time', 'Dship_lat', 'Dship_lon'], encoding='iso-8859-1')
    dship['datetime'] = pd.to_datetime(dship['Date']+' '+dship['Time'])
    del dship['Date']
    del dship['Time']

    return dship


# position data as a collection of csv sheets, posidonia, qgis
def posidata_PS137(posi_path):
    """
    Reads in comma separated files with data from the Posidonia acoustic position tracker on RV Polarstern from the cruise PS137. Megrges all files in the "path" directory. Corrects data specific errors.

    Parameters
    ----------
    posi_path : string
        Path to the file.

    Returns
    -------
    posi : pandas dataframe
        Pandas Dataframe with datetime, Poalrstern latitude, Polarstern longitude, CTD latitude, CTD longitude and CTD depth as columns.
    """

    import pandas as pd
    import os
    import datetime as dt
    # from scipy import stats

    file_list = []  # get all filenames of the positiondata csv sheets
    for file in os.listdir(posi_path):
        if file.endswith(".csv"):
            file_list.append(file)
            file_list.sort()

    posi = pd.DataFrame()  # initializes dataframe

    for file in file_list:  # reads in position data
        posi_sheets = pd.read_csv(os.path.join(posi_path,file), sep='\t', float_precision='high')
        posi_sheets['datetime'] = pd.to_datetime(
            # creates datetime object column
            posi_sheets['Date']+' '+posi_sheets['Time'], dayfirst=True)
        posi = pd.concat([posi, posi_sheets], axis=0)  # merge all posi sheets

    # splits position data in 5 sections and advances time in two of them by 2/1 hours based on depth correlation with CTD Data
    part1 = posi[posi['datetime'] < dt.datetime.strptime(
        '2023-07-05 18:00:00', '%Y-%m-%d %H:%M:%S')]
    part2 = posi[(posi['datetime'] >= dt.datetime.strptime('2023-07-05 18:00:00', '%Y-%m-%d %H:%M:%S'))
                 & (posi['datetime'] < dt.datetime.strptime('2023-07-06 03:00:00', '%Y-%m-%d %H:%M:%S'))]
    part3 = posi[(posi['datetime'] >= dt.datetime.strptime('2023-07-06 03:00:00', '%Y-%m-%d %H:%M:%S'))
                 & (posi['datetime'] < dt.datetime.strptime('2023-07-08 13:00:00', '%Y-%m-%d %H:%M:%S'))]
    part4 = posi[(posi['datetime'] >= dt.datetime.strptime('2023-07-08 13:00:00', '%Y-%m-%d %H:%M:%S'))
                 & (posi['datetime'] < dt.datetime.strptime('2023-07-10 14:00:00', '%Y-%m-%d %H:%M:%S'))]
    part5 = posi[posi['datetime'] >= dt.datetime.strptime(
        '2023-07-10 15:00:00', '%Y-%m-%d %H:%M:%S')]

    part2['datetime'] = part2['datetime'] + pd.Timedelta(2, 'h')
    part4['datetime'] = part4['datetime'] + pd.Timedelta(1, 'h')

    posi = pd.concat([part1, part2, part3, part4, part5], axis=0)

    # posi = posi[(abs(stats.zscore(posi['CTD Lon'])) < 2)]
    # posi = posi[(abs(stats.zscore(posi['CTD Lat'])) < 2)]
    # posi = posi[(abs(stats.zscore(posi['CTD Depth'])) < 2)]

    posi = posi[['datetime', 'Polarstern Lat', 'Polarstern Lon',
                 'CTD Lat', 'CTD Lon', 'CTD Depth']]  # subset relevant columns
    posi.rename({'datetime': 'datetime', 'Polarstern Lat': 'Ship_lat', 'Polarstern Lon': 'Ship_lon',
                'CTD Lat': 'CTD_lat', 'CTD Lon': 'CTD_lon', 'CTD Depth': 'CTD_depth'}, axis=1, inplace=True)
    posi = posi.sort_values(by='datetime', axis=0)

    return (posi)


def posidata_SO301(posi_path):  # position data as one csv sheet, Ranger output
    """
    Reads in comma separated files with data from the Sonardyne Ranger acoustic position tracker on RV Sonne from the cruise SO301 exported from the "D-Ship" portal. Megrges all files in the "path" directory.

    Parameters
    ----------
    posi_path : string
        Path to the file.

    Returns
    -------
    posi : pandas dataframe
        Pandas Dataframe with datetime, instrument latitude,instrument longitude and intrument depth as columns.
    """

    import pandas as pd
    from lat_lon_parser import parse

    posi = pd.read_csv(posi_path, sep=';', float_precision='high',
                       # reads posi data
                       encoding='iso-8859-1', skiprows=3, names=['datetime', 'Instr_depth', 'Instr_lat', 'Instr_lon'])
    posi['datetime'] = pd.to_datetime(
        posi['datetime'])  # creates datetime object
    posi['Instr_lat'] = posi['Instr_lat'].astype(str).apply(lambda x: parse(x))
    posi['Instr_lon'] = posi['Instr_lon'].astype(str).apply(lambda x: parse(x))

    return (posi)


def posidata_M210(posi_path):
    """
    This function reads in position data of the ship and from four underwater acoustic position transponders for the cruise M210. The data was exported from the DSHIP system. The position data is grouped by transponder and NaN values as well as not needed data is removed.

    Parameters
    ----------
    posi_path : string
        Path to the data file.

    Returns
    -------
    transp_list : dict
        Dictionary containing key-value pairs for four different transponders as pandas dataframes and one key-value pair with unit information.
    """
    import numpy as np
    import pandas as pd
    from lat_lon_parser import parse

    dship = pd.read_csv(posi_path, sep=';',
                        encoding='iso-8859-1')  # reads posi data

    # remove all rows with all nan
    dship_nonan = dship
    dship_units = dship.iloc[1]
    dship_nonan = dship_nonan.replace(' ', np.nan)
    dship_nonan = dship_nonan.replace('#', np.nan)
    dship_nonan = dship_nonan[2:]
    dship_nonan = dship_nonan.dropna(axis=1, how='all')

    # sort transponders into dict
    transp_list = dict()
    transp_list['transp0'] = dship_nonan[['date time', 'SYS.DISP.ActName', 'SYS.STR.DPT', 'SYS.STR.PosLat', 'SYS.STR.PosLon',
                                          'POSI.PTSAG.0.Depth_BUC',
                                          'POSI.PTSAG.0.position_latitude',
                                          'POSI.PTSAG.0.NS',
                                          'POSI.PTSAG.0.position_longitude',
                                          'POSI.PTSAG.0.EW',
                                          'POSI.PTSAG.0.raw_time',
                                          'POSI.PTSAG.0.transponder_No']]

    transp_list['transp1'] = dship_nonan[['date time', 'SYS.DISP.ActName', 'SYS.STR.DPT', 'SYS.STR.PosLat', 'SYS.STR.PosLon',
                                          'POSI.PTSAG.1.Depth_BUC',
                                          'POSI.PTSAG.1.position_latitude',
                                          'POSI.PTSAG.1.NS',
                                          'POSI.PTSAG.1.position_longitude',
                                          'POSI.PTSAG.1.EW',
                                          'POSI.PTSAG.1.raw_time',
                                          'POSI.PTSAG.1.transponder_No']]

    transp_list['transp2'] = dship_nonan[['date time', 'SYS.DISP.ActName', 'SYS.STR.DPT', 'SYS.STR.PosLat', 'SYS.STR.PosLon',
                                          'POSI.PTSAG.2.Depth_BUC',
                                          'POSI.PTSAG.2.position_latitude',
                                          'POSI.PTSAG.2.NS',
                                          'POSI.PTSAG.2.position_longitude',
                                          'POSI.PTSAG.2.EW',
                                          'POSI.PTSAG.2.transponder_No']]

    transp_list['transp3'] = dship_nonan[['date time', 'SYS.DISP.ActName', 'SYS.STR.DPT', 'SYS.STR.PosLat', 'SYS.STR.PosLon',
                                          'POSI.PTSAG.3.Depth_BUC',
                                          'POSI.PTSAG.3.position_latitude',
                                          'POSI.PTSAG.3.NS',
                                          'POSI.PTSAG.3.position_longitude',
                                          'POSI.PTSAG.3.EW',
                                          'POSI.PTSAG.3.transponder_No']]

    transp_list['transp4'] = dship_nonan[['date time', 'SYS.DISP.ActName', 'SYS.STR.DPT', 'SYS.STR.PosLat', 'SYS.STR.PosLon',
                                          'POSI.PTSAG.4.Depth_BUC',
                                          'POSI.PTSAG.4.position_latitude',
                                          'POSI.PTSAG.4.NS',
                                          'POSI.PTSAG.4.position_longitude',
                                          'POSI.PTSAG.4.EW',
                                          'POSI.PTSAG.4.transponder_No']]

    # rename and convert columns
    for key in transp_list:
        transp = transp_list[key]
        transp.columns = transp.columns.str.replace('POSI.PTSAG.', '')
        no = key[-1]
        transp.columns = transp.columns.str.replace(no, '')
        transp.columns = transp.columns.str.replace('.', '')
        transp = transp.rename(columns={'date time': 'Datetime', 'SYSDISPActName': 'Station', 'SYSSTRDPT': 'Water_depth', 'SYSSTRPosLat': 'Ship_lat',
                               'SYSSTRPosLon': 'Ship_lon', 'Depth_BUC': 'Depth', 'position_latitude': 'Instr_lat', 'position_longitude': 'Instr_lon', 'transponder_No': 'Transp_no'})

        # drop nan rows
        transp = transp.dropna(
            subset=['Depth', 'Instr_lat', 'Instr_lon'], how='all')
        # convert datatypes
        transp['Depth'] = transp['Depth'].astype(float)
        transp['Datetime'] = pd.to_datetime(transp['Datetime'])
        transp['Water_depth'] = transp['Water_depth'].astype(float)
        try:
            transp['Ship_lat'] = transp['Ship_lat'].astype(
                str).apply(lambda x: parse(x))
        except:
            transp['Ship_lat'] = np.nan
        try:
            transp['Ship_lon'] = transp['Ship_lon'].astype(
                str).apply(lambda x: parse(x))
        except:
            transp['Ship_lon'] = np.nan
        transp['Instr_lat'] = transp['Instr_lat'].astype(
            str).apply(lambda x: parse(x))
        transp['Instr_lon'] = transp['Instr_lon'].astype(
            str).apply(lambda x: parse(x))
        transp['Transp_no'] = transp['Transp_no'].astype(int)

        transp_list[key] = transp
    transp_list['units'] = dship_units

    return transp_list


def read_btl(btl_path):
    """
    [Depreciated, new read_btl using a different seabird package]
    Reads in bottle data from SeaBird CTD processing software which ends with ".btl"
    to a dictionary with the filename as the key and a pandas dataframe per
    station as the value. It exspects one file per station in the given
    directory. This function uses the seabird package by castelao: https://github.com/castelao/seabird

    Parameters
    ----------
    btl_path : string
        Path to the directory where the data is stored ideally in one station per file.

    Returns
    -------
    btl_data : dictionary
        Dictionary where the keys are unique station identifiers and
        values are data in pandas dataframes.

    Reference
    ---------
    Code provided by N. Froehberg.
    """

    import pandas as pd
    import os

    btl_data = dict()

    file_list = []  # get all filenames of the btl files sheets
    for file in os.listdir(btl_path):
        if file.endswith(".btl"):
            file_list.append(file)
            file_list.sort()

    for file in file_list:
        f = open(os.path.join(btl_path,file), 'r', encoding='latin-1')
        first_data_line = 0
        while True:
            line = f.readline()
            if line[0] != '#' and line[0] != '*':
                f.close()
                break
            first_data_line = first_data_line+1
        df = pd.read_fwf(os.path.join(btl_path,file), skiprows=first_data_line,)
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

        df_output = pd.DataFrame(columns=[
                                 'Bottle', 'Date', 'Time', 'Latitude', 'Longitude']+new_columns, index=range(1, n_bottles+1))

        for i in range(n_bottles):
            index = i*4
            date = df.loc[index+1,'Date']
            time = df.loc[index+2,'Date']
            lat = df.loc[index+1,'Latitude']
            long = df.loc[index+1,'Longitude']
            bottle = df.loc[index+1,'Bottle']
            for column in columns:
                df_output.loc[i+1,column+'_mean'] = df.loc[index+1,column]
                df_output.loc[i+1,column+'_sd'] = df.loc[index+2,column]
                df_output.loc[i+1,column+'_min'] = df.loc[index+3,column]
                df_output.loc[i+1,column+'_max'] = df.loc[index+4,column]
                df_output.loc[i+1,'Date'] = date
                df_output.loc[i+1,'Time'] = time
                df_output.loc[i+1,'Latitude'] = lat
                df_output.loc[i+1,'Longitude'] = long
                df_output.loc[i+1,'Bottle'] = bottle

        df_output = df_output.iloc[:, :-4]
        df_output['datetime'] = pd.to_datetime(
            # creates datetime object column
            df_output['Date']+' '+df_output['Time'], dayfirst=True)
        df_output['Bottle'] = df_output['Bottle'].astype(int)
        file = file.replace('.btl','')
        btl_data[file] = df_output
        print('Read '+file+' successfully!')

    # infer datetypes
    for key in btl_data:
        btl_data[key] = btl_data[key].infer_objects(copy=False)

    return btl_data

def read_btl_iow(btl_path, stats_as_columns=True):
    """
    Reads in bottle data from SeaBird CTD processing software which ends with ".btl"
    to a dictionary with the filename as the key and a pandas dataframe per
    station as the value. It exspects one file per station in the given
    directory. This function uses the seabirdfilehandler module: https://git.io-warnemuende.de/CTD-Software/SeabirdFileHandler .

    Parameters
    ----------
    btl_path : string
        Path to the directory where the data is stored ideally in one station per file.
    stats_as_columns : Boolean, Optional
        Parameter which controls if statistical values for each bottle closure are displayed as additional columns for every parameter or as additional rows with a statistics column (as in the btl file). Default is True.

    Returns
    -------
    btl_data : dictionary
        Dictionary where the keys are unique station identifiers and
        values are data in pandas dataframes.

    """
    from seabirdfilehandler import BottleFile
    import os
    import pandas as pd

    btl_data = dict()

    # get all filenames of the btl files sheets
    file_list = []
    for file in os.listdir(btl_path):
        if file.endswith(".btl"):
            file_list.append(file)
            file_list.sort()

    #iterate over filenames and read data into dict
    for file in file_list:
        #key for dict
        key = file.replace('.btl','')
        #load data
        data = BottleFile(os.path.join(btl_path,file)).create_dataframe()
        #remove empty column from two-row header
        del data['Position']

        #append statistics data as columns
        if stats_as_columns==True:
            if 'Statistic' in data.columns:
                # subset columns which have not meaningful statistics
                if 'Longitude' in data.columns:
                    data_static = data[['Bottle','Date','Time','Longitude','Latitude']]
                    data = data.drop(['Bottle','Date','Time','Longitude','Latitude'],axis=1)
                else:
                    data_static = data[['Bottle','Date','Time']]
                    data = data.drop(['Bottle','Date','Time'],axis=1)
                #subset every 4th column to account for less rows
                data_static = data_static.iloc[::4, :]
                data_static = data_static.reset_index(drop=True)
                #split by statistic
                data_mean = data[data['Statistic']=='avg']
                data_sdev = data[data['Statistic']=='sdev']
                data_min = data[data['Statistic']=='min']
                data_max = data[data['Statistic']=='max']
                for df in [data_mean,data_sdev,data_min,data_max]:
                    #append statistic to column names
                    if df['Statistic'].unique() =='avg':
                        stat = '_mean'
                    else:
                        stat = '_'+df['Statistic'].unique()
                    del df['Statistic']
                    df.columns = df.columns+stat
                #reset index for concat
                data_mean = data_mean.reset_index(drop=True)
                data_sdev = data_sdev.reset_index(drop=True)
                data_min  = data_min.reset_index(drop=True)
                data_max  = data_max.reset_index(drop=True)
                #concat back to one df
                data = pd.concat([data_static,data_mean,data_sdev,data_min,data_max],axis=1)

        #create datetime object column
        data['datetime'] = pd.to_datetime(data['Date']+' '+data['Time'], dayfirst=True)
        #write data into dictionary
        btl_data[key] = data
        #infer data types
        btl_data[key] = btl_data[key].infer_objects(copy=False)
        print('Read '+file+' successfully!')


    return btl_data

def read_mapr(mapr_path, raw=False, error_sheets=[]):
    """
    Reads in MAPR data from NOAA PMEL MAPR sensor excels sheets
    to a dictionary with the filename as the key and a pandas dataframe per
    station as the value. It exspects one file per sensor and station in the given
    directory.

    Parameters
    ----------
    mapr_path : string
        Path to the directory where the data is stored in one sensor and station per file.
    raw : bool, optional
        Flag which controls if metadata is stripped and the format is reduced to a table containing only relevant data. If raw = True, the original format and metadata is kept. This is intended to produce readible raw data files for archiving. The defaul is False.
    error_sheets : list, optional
        If the excel document contains single sheets with errornous data, the corresponding sheet names can be given to the function as string in this list and are skipped. This can be expecially helpfull as pandas seems to read sometimes one additional errornous sheet called 'MAPR_files module'. Default is an empty list.

    Returns
    -------
    mapr_data : dictionary
        Dictionary where the keys are unique sensor and station identifiers and
        values are data in pandas dataframes.
    """

    import pandas as pd

    mapr_all = pd.read_excel(mapr_path, sheet_name=None, header=None)
    sheets = list(mapr_all.keys())

    mapr = dict()
    column_names = dict()

    for sheet in sheets:
        print(sheet)
        # exclude sheets wich are falsely read by pandas
        if sheet in error_sheets:
            print('Skipping sheet '+sheet)
        else:
            try:
                if raw == False:
                    # select only data
                    mapr[sheet] = mapr_all[sheet].iloc[3:, 1:14]
                    # mapr[sheet] = mapr[sheet].iloc[:,[0,1,2,3,4,5,10,11,12]] # only quantities with physical units
                    # keys/quantities & units from first row,
                    column_names[sheet] = mapr[sheet].loc[3, :]
                    column_names[sheet] = column_names[sheet].str.split(
                        '\n', expand=True)
                    column_names[sheet] = column_names[sheet][0] + \
                        column_names[sheet][1]
                    mapr[sheet].columns = column_names[sheet]
                    mapr[sheet].rename(
                        {'date/time(dd/mm/yyyy hh:mm:ss)': 'datetime'}, axis=1, inplace=True)
                    # mapr[sheet]['datetime'] = pd.to_datetime(mapr[sheet]['datetime'])
                    mapr[sheet] = mapr[sheet].drop(3, axis=0)
                    mapr[sheet]['datetime'] = pd.to_datetime(
                        mapr[sheet]['datetime'])
                else:
                    mapr[sheet] = mapr_all[sheet]
                    #remove empty columns
                    mapr[sheet] = mapr[sheet].iloc[:,:14]
                print('Read '+sheet + ' successfully')
                # to check for proper datetimes
                # if (mapr[sheet]['datetime'].is_monotonic_increasing and mapr[sheet]['datetime'].is_unique)==False:
                #     print('Monotonic: '+mapr[sheet]['datetime'].is_monotonic_increasing)
                #     print('Unique: '+mapr[sheet]['datetime'].is_unique())

            except:
                print('Error while reading sheet '+sheet+', skipping')

    # mapr['PS137_OFOBS59-1_74'] = mapr['PS137_OFOBS59-1_74'][mapr['PS137_OFOBS59-1_74']['date/time']<pd.to_datetime('2023-07-26 04:50:00')]
    # what for? -> check times

    # infer datetypes
    for key in mapr:
        mapr[key] = mapr[key].infer_objects(copy=False)

    return mapr


def read_he(he_path):
    """
    Reads in helium isotope data form excel sheets from the Bremen Helium Isotope Lab.

    Parameters
    ----------
    he_path : string
        Path to the file.

    Returns
    -------
    he_data_full : pandas dataframe
        Pandas dataframe with all quantities as columns.
    """

    import pandas as pd
    import numpy as np

    try:
        he_data_full = pd.read_csv(he_path)
    except:
        he_data_full = pd.read_csv(he_path, encoding="ISO-8859-1")

    # del he_data_full['Unnamed: 54']
    # del he_data_full['Unnamed: 55']

    he_data_full.replace('kommt noch', np.nan, inplace=True)
    # he_data_full['%fname'].iloc[65] = np.nan
    # he_data['datetime'] = pd.to_datetime(he_data['jultime'],unit='D',origin='julian')
    # he_data_full['station'] = 'PS137_0' + he_data_full['station'].astype(str)+'_01'
    he_data_full['station'] = he_data_full['station'].astype(str)
    for index, row in he_data_full.iterrows():
        if len(row['station']) == 1:
            he_data_full.loc[index,'station'] = '00'+he_data_full.loc[index,'station']
        if len(row['station']) == 2:
            he_data_full.loc[index,'station'] = '0'+he_data_full.loc[index,'station']
    he_data_full['station'] = he_data_full['station']+'_01'

    he_data_full.rename({'bottle': 'Bottle'}, axis=1, inplace=True)

    return he_data_full


def read_gebco(gebcopath, lonE, lonW, latS, latN):  # hillshade_path
    """
    Reads in bathymetry from NetCDF file.

    Parameters
    ----------
    gebcopath : string
        Path to the file.

    Returns
    -------
    lats,lons,elev : tuple of three numpy ndarrays
        lats: latitude coordinates
        lons: longitude coordinates
        elevs: elevation
    """

    import xarray as xa
    import numpy as np

    # gebcopath = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84.nc'     # path for bathymetry files
    # gebcopath_hillshade = '/home/jonathan/Dokumente/SHK Maren/PS137/PS137_AuroraVent_25m_bilinear_WGS84_hillshade.nc'

    gebco = xa.open_dataset(gebcopath)

    # lonW = -6.45
    # lonE = -6.175
    # latS = 82.875
    # latN = 82.92
    # res = '10m'

    gebco = gebco.sel(lat=slice(latS, latN), lon=slice(lonW, lonE))

    elev = gebco.variables['Band1'][:]
    lats = gebco.variables['lat'][:]
    lons = gebco.variables['lon'][:]
    lons, lats = np.meshgrid(lons, lats)

    return (lons, lats, elev)
