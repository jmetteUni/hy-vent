#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:03:01 2024

@author: jonathan
"""

def qc_IQR(profile, vars_to_qc, threshold, boxplots=False):        #qc on no of variables in a profile with the IQR method, originally designed for posidonia latlon data

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
            df.boxplot(column=var) #,label='w/o qc')
            plt.title('w/o qc')
            plt.show()

            plt.figure()
            df_qc.boxplot(column=var) #,label='with qc, '+str(threshold))
            plt.title('with qc, '+str(threshold))
            plt.show()

    if ('CTD_lat' in vars_to_qc) or ('CTD_lon' in vars_to_qc):      #if var_to_qc contains posi lon lat data, drop all rows with some nan, to only get valif lon-lat combinations
        df_qc.dropna(inplace=True)

    for var in vars_to_qc:         #writes qc vars back into profile
        profile[var] = df_qc[var]

    return profile          #returns the profile with the original vars overwritten by the qc vars
