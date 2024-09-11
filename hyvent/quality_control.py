#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:03:01 2024

@author: jmette@uni-bremen.de
"""

def qc_lat_lon_IQR(profile, vars_to_qc, threshold, boxplots=False):        #qc on no of lat lon in a profile with the IQR method #threshold example =1.5
    """
    This function removes outliers based on the interquartile ranges method on selected lat lon coordinates in a dataframe. Optional, box plots are produced for controlling the outlier detection. If the selected variabels are acoustic tracker positions names "CTD_lat" or "CTD_lon", the corresponing coordinate for an outlier is also removed for consistency.

    Parameters
    ----------
    profile : pandas dataframe
        Data wich should be quality conttrolled.
    vars_to_qc : list of strings
        List of variables were outliers should be removed. If one of the is called "CTD_lat" or "CTD_lon" the corresponding coordinate to the outlier is also removed.
    threshold : int or float
        Threshold of the outlier removal, which defines the upper and lower bounds.
    boxplots : boolian, optional
        Boolian which controlles if boxplots are produced or not. The default is False.

    Returns
    -------
    profile: pandas dataframe
        Returns the data with outliers as NaN

    """

    import matplotlib.pyplot as plt

    df = profile[vars_to_qc]       #select variables which should be qc

    df_qc = df.copy()
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

    if ('CTD_lat' in vars_to_qc) or ('CTD_lon' in vars_to_qc):      #if var_to_qc contains posi lon lat data, drop all rows with some nan, to only get valif lon-lat combinations
        df_qc.dropna(inplace=True)

    for var in vars_to_qc:         #writes qc vars back into profile
        profile[var] = df_qc[var]

    return profile          #returns the profile with the original vars overwritten by the qc vars



def qc_IQR(profile, var_to_qc, threshold, boxplots=False):        #qc on no of variables in a profile with the IQR method, originally designed for posidonia latlon data #threshold example =1.5
    """
    This function removes outliers based on the interquartile ranges method on selected variables in a dataframe. Optional, box plots are produced for controlling the outlier detection.

    Parameters
    ----------
    profile : pandas dataframe
        Data wich should be quality conttrolled.
    vars_to_qc : list of strings
        List of variables were outliers should be removed. If one of the is called "CTD_lat" or "CTD_lon" the corresponding coordinate to the outlier is also removed.
    threshold : int or float
        Threshold of the outlier removal, which defines the upper and lower bounds.
    boxplots : boolian, optional
        Boolian which controlles if boxplots are produced or not. The default is False.

    Returns
    -------
    profile: pandas series
        Returns the data with outliers as NaN

    """

    import matplotlib.pyplot as plt

    df = profile

    df_qc = df.copy()
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

