#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 12:52:37 2024

@author: jmette@uni-bremen.de
"""

def corr_mapr_depth(data,lat):
    """

    Calculates a corrected depth based on the 200 first pressure values larger then 5 dbar which are assumed as measured before deploying. The new, corrected depth is calculated with Gibbs Seawater Toolbox (gsw.conversions.z_from_p()) with the pressure corrected by the mean of the pre-deployment values and written to the dataframe as a new column.


    Parameters
    ----------
    data : pandas dataframe
        Data of one station of one MAPR as a pandas dataframe.
    lat : interger or float
        Latitude where the measurement was taken to calculate the depth from pressure.

    Returns
    -------
    depth_corr : pandas series
        Column containining the corrected depth.

    """

    import gsw

    p_deck = data[:200]     #take first 200 measurements which are negative
    p_deck = p_deck[p_deck['Press(dB)']<5]      #exclude values which are already in water measured
    if len(p_deck) < 50:        #warn if only a few measurements are pre-deplyoment
        print('Warning: Less then 50 values used for deck reference pressure!')
    p_deck_mean = p_deck['Press(dB)'].mean()
    #print('Mean deck pressure: '+str(p_deck_mean)+' dB')
    data['Depth_corr(m)'] = -data.apply(lambda x: gsw.conversions.z_from_p(x['Press(dB)']-p_deck_mean,lat), axis=1)         #calculates the corrected depth

    depth_corr = data['Depth_corr(m)']    #select only the corrected depth column

    return depth_corr

def calc_mean_profile(data, var, station_list):
    """
    This function calculates a mean profile of the data of mutliple stations for a given variable. Values a binned to 1m depth intervals between the minimum and maximum depth present in the data.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data of one or multiple CTD stations with variables as columns. One column has to be the station designation.
    var : string
        Variable to calculate the mean of. Must be a column key in data.
    station_list : list of strings
        List of station designations, over which the mean should be calculated.

    Returns
    -------
    mean_profile : pandas dataframe
        Dataframe containing two columns which are the binned depth and the associated mean of the variable.

    """
    import numpy as np
    import pandas as pd

    min_depth = data['DEPTH'].min()
    max_depth = data['DEPTH'].max()
    depth_vec = np.arange(min_depth,max_depth,1)        #create depth vector with step 1m, where variable is binned to
    all_profiles = pd.DataFrame(depth_vec,columns=['DEPTH'])
    for station in station_list:
        profile = data[data['Station']==station][['DEPTH',var]]
        new_key = var+station
        profile.rename(columns={var:new_key},inplace=True)
        profile.sort_values(by='DEPTH',ascending=True,inplace=True)
        profile.dropna(inplace=True)
        all_profiles = pd.merge_asof(all_profiles,profile,on='DEPTH',direction='nearest')       # for every station merge variable by nearest depth value -> binning
    del all_profiles['DEPTH']
    mean_name = str(var)+'_mean'
    all_profiles[mean_name] = all_profiles.mean(axis=1)    #calculate mean variable value
    all_profiles['DEPTH'] = depth_vec
    mean_profile = all_profiles[['DEPTH',mean_name]]

    return mean_profile

def derive_mapr(data, data_for_mean, station_list):
    """
    This function calculates derived physical properties with the Gibbs Seawater Toolbox for a dataset (usually MAPR) without salinity measurements. The salinity values are calculated as a mean over the stations given from a second dataset (usually CTD stations).

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data without salinity values of one or multiple CTD stations or MAPR operations with variables as columns.
    data_for_mean : pandas dataframe
        Dataframe with the data with salinity values of one or multiple CTD stations with variables as columns. One column has to be the station designation. The salinity variable is assumed to be called "PSAL" and be in Practical Salinity units.
    station_list : list of strings
        List of station designations, over which the mean for salinity should be calculated.

    Returns
    -------
    data_der : pandas dataframe
        Dataframe with the original data and the derived properties.

    """
    from hyvent.processing import calc_mean_profile
    import pandas as pd
    import gsw

    var = 'PSAL'

    mean_profile = calc_mean_profile(data_for_mean, var, station_list)       #calculates a mean profile of PSAL from the list of stations given

    data.sort_values(by='DEPTH',ascending=True,inplace=True)
    data_der = pd.merge_asof(data, mean_profile, on='DEPTH', direction='nearest',tolerance=10)     #merges mean PSAL profile with data based on nearest depth value with maximum tolerance of 10m

    data_der['SA'] = gsw.conversions.SA_from_SP(data_der['PSAL_mean'], data_der['PRES'], data_der['Dship_lon'], data_der['Dship_lat'])     #calculates derived parameters
    data_der['potemperature'] = gsw.conversions.pt0_from_t(data_der['SA'],data_der['TEMP'],data_der['PRES'])
    data_der['CT'] = gsw.conversions.CT_from_pt(data_der['SA'], data_der['potemperature'])
    data_der['Rho'] = gsw.density.rho(data_der['SA'],data_der['CT'],data_der['PRES'])
    data_der['Sigma3'] = gsw.density.sigma3(data_der['SA'],data_der['CT'])

    data_der.sort_values(by='datetime',ascending=True,inplace=True)
    del data_der['PSAL_mean']

    return data_der

def derive_CTD(data):
    """
    This function calculates the derived properties absolute salinity, conservative temperature, density and density anomaly of CTD data with the Gibbs Seawater Toolbox and appends them as columns.

    Parameters
    ----------
    data : pandas dataframe or dictionary
        Dataset as one dataframe or a dictonary with dataframes as values. Must contain practical salinity ("PSAL"), pressure ("PRES"), longitude ("LONGITUDE"), latitude ("LATITUDE") and potential temperature ("potemperature")

    Returns
    -------
    data : pandas dataframe or dictionary
        Data with the derived properties as additonal columns.

    """
    import gsw
    import pandas as pd

    if isinstance(data,pd.DataFrame):

        data['SA'] = gsw.conversions.SA_from_SP(data['PSAL'], data['PRES'], data['LONGITUDE'], data['LATITUDE'])
        data['CT'] = gsw.conversions.CT_from_pt(data['SA'], data['potemperature'])
        data['Rho'] = gsw.density.rho(data['SA'],data['CT'],data['PRES'])
        data['Sigma3'] = gsw.density.sigma3(data['SA'],data['CT'])

    elif isinstance(data,dict):

        for key in data:
            data[key]['SA'] = gsw.conversions.SA_from_SP(data[key]['PSAL'], data[key]['PRES'], data[key]['LONGITUDE'], data[key]['LATITUDE'])
            data[key]['CT'] = gsw.conversions.CT_from_pt(data[key]['SA'], data[key]['potemperature'])
            data[key]['Rho'] = gsw.density.rho(data[key]['SA'],data[key]['CT'],data[key]['PRES'])
            data[key]['Sigma3'] = gsw.density.sigma3(data[key]['SA'],data[key]['CT'])

    else:
        print('Could not derive properties, must be pd.DataFrame or dict.')

    return data

def sep_casts(data, window_size=1000):
    """
    This function separates a station into up and down casts by local extrema.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data of one station.
    window_size : TYPE, optional
        Size of the rolling window, to find local extrema within. The optimal size is dependent on the number of measurements per up or down cast. The default is 1000.

    Returns
    -------
    casts : list of pandas dataframes
        A list were each entry is a dataframe with the data of the consecutive up and down casts.
    """
    import pandas as pd

    local_min_vals = data.loc[data['DEPTH'] == data['DEPTH'].rolling(window_size, center=True).min()]
    local_max_vals = data.loc[data['DEPTH'] == data['DEPTH'].rolling(window_size, center=True).max()]
    #concat all local min, max and the start and end rows
    extrema = pd.concat([local_max_vals,local_min_vals])
    extrema = pd.concat([extrema,data[data.index==data.index.min()]])
    extrema = pd.concat([extrema,data[data.index==data.index.max()]])
    extrema = extrema.sort_index()

    #separate data by the extrema indices and store it in list
    casts = list()
    indexes = extrema.index.to_list()
    for i in  range(len(indexes)-1):
        start = indexes[i]
        end = indexes[i+1]

        cast = data[(data.index >= start) & (data.index < end)]
        casts.append(cast)

    return casts

def sel_bg_cast(data):
    """
    This function subsets all data until the deepest measurement which is expected to be the first downcast. This can be used as a background or standard profile.

    Parameters
    ----------
    data : pandas dataframe
        Dataframe with the data of one station. Must contain a column "DEPTH" with depth values.

    Returns
    -------
    data : pandas dataframe
        Dataframe with the data until the deepest value.

    """
    data = data.reset_index(drop=True)
    i_deepest = data['DEPTH'].idxmax()
    data = data.iloc[:i_deepest]

    return data



def subtract_bg_by_press(data, bg, var, min_dep, max_dep, control_plot=False):
    """
    This functions calculates the deviation of a variable from a number of profiles from a background by comparing equal density layers and then substracting the variable from the background. For the background measurement, all measurements until the deepest one (downcast) are used. Values outside a depth range are replaces with NaNs. Optionally a control plot comparing the data and the background is shown.

    Parameters
    ----------
    data : pandas dataframe
        Dataset of one or several profiles.
    bg : pandas dataframe
        Dataset which should be used as background.
    var : string
        Variable, where the deviation should be calculated.
    min_dep : int
        Minimum depth, above which the calculated deviation is replaced with NaN
    max_dep : int
        Maximum depth, below which the calculated deviation is replaced with NaN
    control_plot : boolean, optional
        Boolian which controlles if control plots are produced or not. The default is False.

    Returns
    -------
    data : pandas dataframe
        Dataset like the input but with the calculated deviation as an additional column.

    """

    import matplotlib.pyplot as plt
    import pandas as pd

    #selects first downcast as background
    bg = sel_bg_cast(bg)

    #prepares background for merge
    bg = bg[['Rho',var]]
    bg = bg.dropna(subset=['Rho'])
    bg = bg.sort_values(by='Rho')

    #prepares data for merge
    data = data.reset_index(drop=True)
    data['datetime'] = pd.to_datetime(data['datetime'])
    merge = data[['datetime','DEPTH','Rho','Station','SN',var]]
    merge = merge.dropna(subset=['Rho'])
    merge = merge.sort_values(by='Rho')

    #merges data and background on density "Rho" with a tolerance and calculates delta of the variable
    merge = pd.merge_asof(merge, bg, on='Rho',tolerance=0.01,direction='nearest')
    merge['Delta_'+var] = merge[var+'_x'] - merge[var+'_y']
    if merge['Delta_'+var].isna().any() == True:
        print('Delta_'+var+' was not calculated for all values, because there no equal background density value was found inside the tolerance limit.')
    merge = merge.sort_values(by='datetime')

    for column in ['Delta_'+var,var+'_x',var+'_y']:
        merge[column] = merge[column].where((merge['DEPTH']>=min_dep) & (merge['DEPTH']<=max_dep))

    if control_plot == True:
        merge_list = [d for _, d in merge.groupby(['Station','SN'])]
        for station in merge_list:
            plt.figure()
            plt.plot(station['datetime'],station[var+'_x'])
            plt.plot(station['datetime'],station[var+'_y'])
            #plt.plot(station['datetime'],station['Delta_'+var)         #can be used to plot Delta
            plt.title(station['Station'].iloc[0]+', '+station['SN'].iloc[0])
            plt.show()

    merge = merge[['Delta_'+var]]
    data = data.join(merge)

    return data

# calculates Temperature anomaly based on background and unique density values as layers
def calcTempAnomaly(dataset, background):
    """
    [Depreciated, just for dicumentation] This function calculates a temperature anomaly between a dataset and a background dataset based on nearest common density values. The background curve is extracted from the background dataset by taking all values until the deepest pressure value.

    Parameters
    ----------
    dataset : pandas dataframe
        Dataframe with the data of one station.
    background : pandas dataframe
        Dataset which should be used as background.

    Returns
    -------
    data : pandas dataframe
        Dataset like the input but with the background temperature and the calculated temperature anomaly as an additional column.

    """
    import pandas as pd

    bg = background[background['PRES'] > 1800]
    # select downcast to the highest pressure value
    bg = bg[bg.index <= int(bg[['PRES']].idxmax())].reset_index()
    # get unique density and corellated temperature
    density_layers = bg.groupby('density')['potemperature'].mean().reset_index()
    density_layers.rename(
        {'density': 'density_l', 'potemperature': 'Tempbg'}, axis=1, inplace=True)

    dataset = dataset[dataset['PRES'] > 1800]
    dataset = pd.merge_asof(dataset.sort_values(by='density'), density_layers, left_on='density', right_on='density_l', direction='nearest').sort_values(
        by='datetime')  # merge dataset and density layers from backgroundcast on nearest density values
    # calculate Temperature anomaly by substracting temperature - background temperature
    dataset['TempAnomaly'] = dataset['potemperature'] - dataset['Tempbg']

    return(dataset)


# calculates Temperature anomaly based on background and No of pressure layers
def calcTempAnomalyPress(dataset, background, No_layers):
    """
    [Depreciated, just for documentation] This function calculates a temperature anomaly between a dataset and a background dataset based for common vertical pressure layer. The background curve is extracted from the background dataset by taking all values until the deepest pressure value. The number of pressure layer (e.g. resolution) must be specified.

    Parameters
    ----------
    dataset : pandas dataframe
        Dataframe with the data of one station.
    background : pandas dataframe
        Dataset which should be used as background.
    No_layers : int
        Number of pressure layers for which one temperature anomaly is calculated.

    Returns
    -------
     data : pandas dataframe
         Dataset like the input but with the background temperature and the calculated temperature anomaly as an additional column.

    """
    import numpy as np
    import pandas as pd

    bg = background[background['PRES'] > 1800]
    # select downcast to the highest pressure value
    bg = bg[bg.index <= int(bg[['PRES']].idxmax())].reset_index()

    press_layers = pd.DataFrame()
    press_layers['PRES_l'] = np.linspace(
        bg['PRES'].iloc[0], bg['PRES'].iloc[-1], No_layers)  # linear spaced pressure layers
    press_layers['PRES_mean'] = np.nan
    press_layers['Tempbg'] = np.nan
    # for every layer group bg by Press values between top of layer and next layer pressure
    for index, row in press_layers.iterrows():
        if index != len(press_layers)-1:  # accounts for out-of-bounds index
            layer = bg[((bg['PRES'] >= press_layers['PRES_l'].iloc[index]) & (
                bg['PRES'] < press_layers['PRES_l'].iloc[index+1]))]
        else:
            # bottom layer equals last row in background,
            layer = bg[(bg['PRES'] >= press_layers['PRES_l'].iloc[index])]

        # mean Pressure per layer for checking algorithm
        press_layers['PRES_mean'][index] = layer['PRES'].mean()
        # mean Temperature for this layer
        press_layers['Tempbg'][index] = layer['potemperature'].mean()

    dataset = dataset[dataset['PRES'] > 1800]
    dataset = pd.merge_asof(dataset.sort_values(by='PRES'), press_layers, left_on='PRES', right_on='PRES_l', direction='nearest').sort_values(
        by='datetime')  # merge dataset and density layers from backgroundcast on nearest density values
    # calculate Temperature anomaly by substracting temperature - background temperature
    dataset['TempAnomaly'] = dataset['potemperature'] - dataset['Tempbg']

    return(dataset)