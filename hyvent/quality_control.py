#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:03:01 2024

@author: jmette@uni-bremen.de
"""

def qc_lat_lon_IQR(profile, vars_to_qc, threshold=1.5, boxplots=False):        #qc on no of lat lon in a profile with the IQR method #threshold example =1.5
    """
    This function removes outliers based on the interquartile ranges method on selected lat lon coordinates in a dataframe. Optional, box plots are produced for controlling the outlier detection. If the selected variabels are acoustic tracker positions names "CTD_lat" or "CTD_lon", the corresponing coordinate for an outlier is also removed for consistency.

    Parameters
    ----------
    profile : pandas dataframe
        Data wich should be quality controlled.
    vars_to_qc : list of strings
        List of variables were outliers should be removed. If one of the is called "CTD_lat" or "CTD_lon" the corresponding coordinate to the outlier is also removed.
    threshold : int or float, optional
        Threshold of the outlier removal, which defines the upper and lower bounds. Default is 1.5.
    boxplots : boolian, optional
        Boolian which controlles if boxplots are produced or not. The default is False.

    Returns
    -------
    profile: pandas dataframe
        Returns the data with outliers as NaN

    """

    import matplotlib.pyplot as plt

    df = profile[vars_to_qc]       #select variables which should be qc

    df_qc = df.copy(deep=True)
    # Calculate 25% percentile and 75% percentile
    for var in vars_to_qc:
        Q1 = df[var].quantile(0.25)
        Q3 = df[var].quantile(0.75)
        # Calculate Interquartile Range (IQR)
        IQR = Q3 - Q1
        # Define lower and upper bounds
        lower_bound = Q1 - threshold * IQR
        upper_bound = Q3 + threshold * IQR

        df_qc[var] = df[var][(df[var] >= lower_bound) & (df[var] <= upper_bound)]       #selects only the data which meets the conditions

    if boxplots == True:        #boxplot of qc variables
        for var in vars_to_qc:
            plt.figure()
            df[var] = df[var].astype(float)   #convert to float for boxplot function
            df.boxplot(column=var) #,label='w/o qc')
            plt.title('w/o qc')
            plt.show()

            plt.figure()
            df_qc[var] = df_qc[var].astype(float)
            df_qc.boxplot(column=var) #,label='with qc, '+str(threshold))
            plt.title('with qc, '+str(threshold))
            plt.show()

    if ('CTD_lat' in vars_to_qc) or ('CTD_lon' in vars_to_qc):      #if var_to_qc contains posi lon lat data, drop all rows with some nan, to only get valid lon-lat combinations
        df_qc.dropna(inplace=True)

    for var in vars_to_qc:         #writes qc vars back into profile
        profile[var] = df_qc[var]

    return profile          #returns the profile with the original vars overwritten by the qc vars



def qc_IQR(profile, var_to_qc, threshold=1.5, boxplots=False):        #qc on no of variables in a profile with the IQR method, originally designed for posidonia latlon data #threshold example =1.5
    """
    This function removes outliers based on the interquartile ranges method on selected variables in a dataframe. Optional, box plots are produced for controlling the outlier detection.

    Parameters
    ----------
    profile : pandas dataframe
        Data wich should be quality conttrolled.
    vars_to_qc : list of strings
        List of variables were outliers should be removed. If one of the is called "CTD_lat" or "CTD_lon" the corresponding coordinate to the outlier is also removed.
    threshold : int or float, optional
        Threshold of the outlier removal, which defines the upper and lower bounds. Default is 1.5.
    boxplots : boolian, optional
        Boolian which controlles if boxplots are produced or not. The default is False.

    Returns
    -------
    profile: pandas series
        Returns the data with outliers as NaN

    """

    import matplotlib.pyplot as plt

    df = profile

    df_qc = df.copy(deep=True)
    # Calculate 25% percentile and 75% percentile
    Q1 = df[var_to_qc].quantile(0.25)
    Q3 = df[var_to_qc].quantile(0.75)
    # Calculate Interquartile Range (IQR)
    IQR = Q3 - Q1
    # Define lower and upper bounds
    lower_bound = Q1 - threshold * IQR
    upper_bound = Q3 + threshold * IQR

    df_qc[var_to_qc] = df[var_to_qc][(df[var_to_qc] >= lower_bound) & (df[var_to_qc] <= upper_bound)]       #selects only the data which meets the conditions

    if boxplots == True:        #boxplot of qc variables
        plt.figure()
        df[var_to_qc] = df[var_to_qc].astype(float)   #convert to float for boxplot function
        df.boxplot(column=var_to_qc) #,label='w/o qc')
        plt.title('w/o qc')
        plt.show()

        plt.figure()
        df_qc[var_to_qc] = df_qc[var_to_qc].astype(float)
        df_qc.boxplot(column=var_to_qc) #,label='with qc, '+str(threshold))
        plt.title('with qc, '+str(threshold))
        plt.show()

    df_qc = df_qc[var_to_qc]   #select only quantites which were quality controlled

    return df_qc          #returns the profile with the original vars overwritten by the qc vars

def cut_prepost_deploy(data, pres_limit=10, window_limit=10, plot=False):
    """
    This functions cuts the beginning and end of CTD or MAPR stations of a dataset based on the first and the last occurence of consecutive values below a pressure limit. Plots to control the result cna be printed optionally.

    Parameters
    ----------
    data : pandas dataframe or dictionary
        Data which should be cutted. Can be one station as a dataframe or multiple dataframes as values in adicitonary. The pressure column is exspected to be called "Press(dB".
    pres_limit : int, optional
        Consecutive values which mark the start and end of a deployment have to be below this limit. The default is 10dbar.
    window_limit : int, optional
        The number of consecutive values which mark the start and end of a deployment. The default is 10.
    plot : boolean, optional
        Controls if a plot with the original pressure data and the cutted pressure data is printed. The default is False.

    Returns
    -------
    data : pandas dataframe or dictionary
        Dataset with the cutted data. Datatype will be the same as the input.

    """
    import pandas as pd
    import matplotlib.pyplot as plt

    if isinstance(data,pd.DataFrame):
        station = data.copy(deep=True)
        # define the condition
        condition = station['Press(dB)'] >= pres_limit

        # create a rolling window of 10 values
        rolling_condition = condition.rolling(window=window_limit, min_periods=window_limit).mean()

        # get the indices of the first row of each group of 10 consecutive values that meet the condition
        indices = rolling_condition[rolling_condition == 1].index
        station_part = station.iloc[indices[0]:indices[-1]+window_limit]

        if plot == True:
            plt.figure()
            plt.plot(station['Press(dB)'])
            plt.plot(station_part['Press(dB)'])
            plt.ylabel('Pressure in dbar')
            plt.show()

        data = station_part

    elif isinstance(data,dict):
        for key in data:
            station = data[key].copy(deep=True)
            # define the condition
            condition = station['Press(dB)'] >= pres_limit

            # create a rolling window of 10 values
            rolling_condition = condition.rolling(window=window_limit, min_periods=window_limit).mean()

            # get the indices of the first row of each group of 10 consecutive values that meet the condition
            indices = rolling_condition[rolling_condition == 1].index
            station_part = station.iloc[indices[0]:indices[-1]+10]

            if plot == True:
                plt.figure()
                plt.plot(station['Press(dB)'])
                plt.plot(station_part['Press(dB)'])
                plt.ylabel('Pressure in dbar')
                plt.title(key)
                plt.show()

            data[key] = station_part

    else:
        print('Could not cut dataset, must be pd.DataFrame or dict.')

    return data