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

        df_qc[var] = df[var][(df[var] >= lower_bound) & (df[var] <= upper_bound)]       #selects only the data which meets the conditions.

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

def cut_prepost_deploy(data, pres_limit=10, window_limit=10, control_plot=False):
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
    control_plot : boolean, optional
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
        station = station.reset_index(drop=True)

        #infer the pressure variable from the data ('Press(dB)' for MAPR, 'PRES' for CTD)
        if 'Press(dB)' in data.keys():
            pres_var = 'Press(dB)'
        if 'PRES' in data.keys():
            pres_var = 'PRES'

        # Define condition
        condition = station[pres_var] >= pres_limit

        # Rolling mean to find sequences of 10 consecutive True values
        rolling_condition = condition.rolling(window=window_limit, min_periods=window_limit).mean()

        # Indices where all 10 consecutive values meet the condition
        valid_indices = rolling_condition[rolling_condition == 1].index

        if len(valid_indices) > 0:
            # Find the first and last valid index block
            start_index = valid_indices[0] - (window_limit - 1)  # include the first full window
            if start_index < 0:
                start_index = 0

            end_index = valid_indices[-1] + 1  # include up to the last valid window
            if end_index > len(station):
                end_index = len(station)

            # Trim the dataset
            station_part = station.iloc[start_index:end_index].copy()
        else:
            # No segment meets the requirement — return empty DataFrame
            print('No values removed.')

        if control_plot == True:
            plt.figure()
            plt.plot(station['datetime'],station[pres_var])
            plt.plot(station_part['datetime'],station_part[pres_var])
            plt.ylabel('Pressure in dbar')
            plt.title(station['Station'].iloc[0])
            plt.show()

        data = station_part

    elif isinstance(data,dict):
        for key in data:
            station = data[key].copy(deep=True)

            #infer the pressure variable from the data ('Press(dB)' for MAPR, 'PRES' for CTD)
            if 'Press(dB)' in station.keys():
                pres_var = 'Press(dB)'
            if 'PRES' in station.keys():
                pres_var = 'PRES'

            # Define condition
            condition = station[pres_var] >= pres_limit

            # Rolling mean to find sequences of 10 consecutive True values
            rolling_condition = condition.rolling(window=window_limit, min_periods=window_limit).mean()

            # Indices where all 10 consecutive values meet the condition
            valid_indices = rolling_condition[rolling_condition == 1].index

            if len(valid_indices) > 0:
                # Find the first and last valid index block
                start_index = valid_indices[0] - (window_limit - 1)  # include the first full window
                if start_index < 0:
                    start_index = 0

                end_index = valid_indices[-1] + 1  # include up to the last valid window
                if end_index > len(station):
                    end_index = len(station)

                # Trim the dataset
                station_part = station.iloc[start_index:end_index].copy()
            else:
                # No segment meets the requirement — return empty DataFrame
                print('No values removed.')

            if control_plot == True:
                plt.figure()
                plt.plot(station['datetime'],station[pres_var])
                plt.plot(station_part['datetime'],station_part[pres_var])
                plt.ylabel('Pressure in dbar')
                plt.title(key)
                plt.show()

            data[key] = station_part

    else:
        print('Could not cut dataset, must be pd.DataFrame or dict.')

    return data

def despike_pressure(series, window, threshold, control_plot=False):
    """
    Despike a pandas pressure series using a rolling median + MAD (Median Absolute Deviation) method. Spikes are replaced with np.nan and then interpolated.

    Parameters
    ----------
    series : pd.Series
        Input time series data.
    window : int
        Sliding window size (should be odd).
    threshold : float
        Number of MADs from the rolling median to consider a spike.
    control_plot : boolean, optional
        Flag, to plot the original and the despiked series. Default is False.

    Returns
    -------
    pd.Series
        Despiked series with np.nan replacing spikes.
    """
    import numpy as np

    s = series.copy()

    roll_median = s.rolling(window, center=True, min_periods=1).median()
    roll_mad = (s - roll_median).abs().rolling(window, center=True, min_periods=1).median()
    roll_mad[roll_mad == 0] = np.nan  # prevent divide-by-zero

    deviation = (s - roll_median).abs() / roll_mad
    spikes = deviation > threshold

    if control_plot == True:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(s.index,s.values,label='original data')
        plt.plot(s[spikes].index,s[spikes].values,'x',color='r', label='removed spikes')
        plt.xlabel('Index')
        plt.ylabel('Variable')
        plt.legend()

    s[spikes] = np.nan
    #interpolate nans
    s = s.interpolate()

    return s

def check_timestamps(data, f_s = 1.0):
    """
    Checks for consistent timestamps if they are strictly monotonic increasing. Single errors, which are likely due to bitflips are corrected, assuming a regular frequency in the timestamps. It corrects errors of the following nature:

    12:01:24
    12:01:25
    12:01:27
    12:01:28
    12:01:29

    The function operates on subsets by station, so a 'Station' column is required. The function prints the result of the checks to the console

    Parameters
    ----------
    data : pandas dataframe
        Data containing a 'datetime' and a 'Station' column.

    f_s : float, optional
        Frequency of the timestamps in Hz. The default is 1 Hz.

    Returns
    -------
    data : pandas dataframe
        Corrected data.
    """
    import pandas as pd

    data_list = [d for _, d in data.groupby(['Station'])]
    for station in data_list:
        # Compute time difference
        station["time_diff"] = station["datetime"].diff()

        # Define conditions to catch single errors which are likely due to bitflips:
        diff_before = station["time_diff"] == pd.Timedelta(seconds=2*1/f_s)
        diff_after = station["time_diff"].shift(-1) == pd.Timedelta(seconds=0)

        # Rows where both conditions hold
        mask = diff_before & diff_after

        # Print results
        n_errors = len(mask[mask==True])
        print('Station '+station['Station'].iloc[0]+':')
        print(str(n_errors)+' single timestamp errors found and corrected.')

        # Subtract one period from those datetime values
        station.loc[mask, "datetime"] -= pd.Timedelta(seconds=1/f_s)

        # check if strictly monotonic increasing
        is_strict = (station['datetime'].diff().dropna() > pd.Timedelta(0)).all()
        # Print results
        if is_strict == True:
            print('Datetime is strictly monotonic increasing.')
        if is_strict == False:
            print('Dateim is NOT strictly monotonic increasing, still datetime errors present!')

    # update original data with new values from station df
    data.update(station)

    return data
