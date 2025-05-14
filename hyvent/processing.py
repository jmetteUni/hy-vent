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
        profile = profile.dropna()
        all_profiles = pd.merge_asof(all_profiles,profile,on='DEPTH',direction='nearest')       # for every station merge variable by nearest depth value -> binning
    del all_profiles['DEPTH']
    mean_name = str(var)+'_mean'
    all_profiles[mean_name] = all_profiles.mean(axis=1)    #calculate mean variable value
    all_profiles['DEPTH'] = depth_vec
    mean_profile = all_profiles[['DEPTH',mean_name]]

    return mean_profile

def derive_mapr(data, data_for_mean, station_list):
    """
    This function calculates derived physical properties with the Gibbs Seawater Toolbox for a dataset (usually MAPR) without salinity measurements. The salinity values are calculated as a mean over the stations given from a second dataset (usually CTD stations). The function also removes outliers ("Neph_outl(volts)") and additionally smoothes the turbdity data ("Neph_smoo(volts)").

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
    from hyvent.quality_control import qc_IQR
    import pandas as pd
    import gsw
    import numpy as np
    from scipy.signal import savgol_filter

    var = 'PSAL'

    mean_profile = calc_mean_profile(data_for_mean, var, station_list)       #calculates a mean profile of PSAL from the list of stations given

    data.sort_values(by='DEPTH',ascending=True,inplace=True)
    data_der = pd.merge_asof(data, mean_profile, on='DEPTH', direction='nearest',tolerance=10)     #merges mean PSAL profile with data based on nearest depth value with maximum tolerance of 10m

    data_der['SA'] = gsw.conversions.SA_from_SP(data_der['PSAL_mean'], data_der['PRES'], data_der['Dship_lon'], data_der['Dship_lat'])     #calculates derived parameters
    data_der['potemperature'] = gsw.conversions.pt0_from_t(data_der['SA'],data_der['TEMP'],data_der['PRES'])
    data_der['CT'] = gsw.conversions.CT_from_pt(data_der['SA'], data_der['potemperature'])
    data_der['Rho'] = gsw.density.rho(data_der['SA'],data_der['CT'],data_der['PRES'])
    data_der['Sigma3'] = gsw.density.sigma3(data_der['SA'],data_der['CT'])

    #remove outliers and smooth turbidity below 100m
    data_list = [d for _, d in data_der.groupby(['Station','SN'])]
    for data in data_list:
        data_part = data.copy(deep=True)
        data_part['Neph(volts)'] = data_part['Neph(volts)'].mask(data_part['DEPTH']<100,np.nan)      #mask all vlaues in Neph above 100m
        data_part['Neph_outl(volts)'] = qc_IQR(data_part, 'Neph(volts)', 5)     #remove outliers
        data_part['Neph_outl(volts)'] = data_part['Neph_outl(volts)'].interpolate()
        data_part['Neph_smoo(volts)'] = savgol_filter(data_part['Neph_outl(volts)'], 70, 3)       #smooth
        data['Neph_outl(volts)'] = data_part['Neph_outl(volts)']     #write processed columns back to data_partframe
        data['Neph_smoo(volts)'] = data_part['Neph_smoo(volts)']

    data_der = pd.concat(data_list)

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
    import numpy as np

    if isinstance(data,pd.DataFrame):

        data['SA'] = gsw.conversions.SA_from_SP(data['PSAL'], data['PRES'], data['LONGITUDE'], data['LATITUDE'])
        data['CT'] = gsw.conversions.CT_from_pt(data['SA'], data['potemperature'])
        data['Rho'] = gsw.density.rho(data['SA'],data['CT'],data['PRES'])
        data['Sigma3'] = gsw.density.sigma3(data['SA'],data['CT'])
        data['Sigma0'] = gsw.density.sigma0(data['SA'],data['CT'])
        #data['Nsquared'] = np.append(gsw.stability.Nsquared(data['SA'],data['CT'],data['PRES'])[0], np.nan)

    elif isinstance(data,dict):

        for key in data:
            data[key]['SA'] = gsw.conversions.SA_from_SP(data[key]['PSAL'], data[key]['PRES'], data[key]['LONGITUDE'], data[key]['LATITUDE'])
            data[key]['CT'] = gsw.conversions.CT_from_pt(data[key]['SA'], data[key]['potemperature'])
            data[key]['Rho'] = gsw.density.rho(data[key]['SA'],data[key]['CT'],data[key]['PRES'])
            data[key]['Sigma3'] = gsw.density.sigma3(data[key]['SA'],data[key]['CT'])
            #data[key]['Nsquared'] = np.append(gsw.stability.Nsquared(data['SA'],data['CT'],data['PRES'])[0], np.nan)


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

    #for mapr: sort by datetime and reindex.
    #would be also good for ctd data, needs testing
    # if 'Operation' in data.columns:
    #     if 'CTD' in data['Operation']:          # this is a lazy and bad check for MAPR data, needs improvement and testing
    #         if 'datetime' in data.columns:
    #             if data['datetime'].isna().any() == True:
    #                 nan_rows = len(data[data['datetime'].isna()==True])
    #                 print('Warning: Dataset contains '+str(nan_rows)+' rows without datetime. These are removed.')
    #                 data = data[data['datetime'].isna()==False]

    #             data = data.sort_values(by='datetime')
    #         data = data.reset_index(drop = True)
    if 'datetime' in data.columns:
        if data['datetime'].isna().any() == True:
            print('Warning: NaN values in datetime. Remove this prior for correct cast separation!')

    data = data.sort_values(by='datetime')
    data = data.reset_index(drop = True)

    try:
        local_min_vals = data.loc[data['DEPTH'] == data['DEPTH'].rolling(window_size, center=True).min()]
        local_max_vals = data.loc[data['DEPTH'] == data['DEPTH'].rolling(window_size, center=True).max()]
    except:
        try:
            print('Did not find a depth coordinate "DEPTH", trying with "Depth(m))" ')
            local_min_vals = data.loc[data['Depth(m)'] == data['Depth(m)'].rolling(window_size, center=True).min()]
            local_max_vals = data.loc[data['Depth(m)'] == data['Depth(m)'].rolling(window_size, center=True).max()]
        except:
            print('Did not find a depth coordinate "DEPTH", trying with "DepSM_mean" ')
            local_min_vals = data.loc[data['DepSM_mean'] == data['DepSM_mean'].rolling(window_size, center=True).min()]
            local_max_vals = data.loc[data['DepSM_mean'] == data['DepSM_mean'].rolling(window_size, center=True).max()]


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

        if i == len(indexes)-1:
            cast = data[(data.index >= start) & (data.index <= end)]    #include the last row in the last cast
        else:
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

def calc_delta_by_dens(data, bg, var, min_dep, max_dep, control_plot=False):
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
    density_layers = density_layers.rename(
        {'density': 'density_l', 'potemperature': 'Tempbg'}, axis=1)

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

def get_bg_polyfit(bg, dep_vec, var, min_dep, max_dep, fit_order=10):
    """
    This funtions fits a polynomial (np.polyfit) to a variable of a profile (usually to create a smooth background profile), with depth of the main dataset as x values for the fit.

    Parameters
    ----------
    bg : pandas dataframe
        Dataset which should be used as background.
    dep_vec: pandas datafame
        Dataframe with one column ("DEPTH") containing the depth coordinates from the main dataset.
    var : string
        Variable, where the deviation should be calculated.
    min_dep : int or float
        Minimum depth for the fit's depth range.
    max_dep : int or float
    fit_order : int or float, optional
        Order of the polynomial fit. The default is 10.

    Returns
    -------
    df_fit : pandas dataframe
        Datframe containing the depth and the fit result.

    """
    import numpy as np
    import pandas as pd

    #calculate the polynomial fit
    coef, res, rank, singular_val, rcond = np.polyfit(bg['DEPTH'],bg[var],fit_order,full=True)
    print('Residual Polyfit: '+str(res))
    fit_func = np.poly1d(coef)
    fit = fit_func(dep_vec)

    #construct dataframe with depth vector and fit
    df = {'DEPTH':dep_vec['DEPTH'],'Bgfit_'+var:np.ravel(fit)}
    df_fit = pd.DataFrame(df)
    df_fit = df_fit.sort_values(by='DEPTH',ascending=True).dropna(subset=['DEPTH'])

    return df_fit

def get_bg_unifit(bg, dep_vec, var, min_dep, max_dep, k=None, s=0.005):
    """
    This funtions fits a univariate spline (scipy.interpolate.UnivariateSpline) to a variable of a profile (usually to create a smooth background profile), with depth of the main dataset as x values for the fit.

    Parameters
    ----------
    bg : pandas dataframe
        Dataset which should be used as background.
    dep_vec: pandas datafame
        Dataframe with one column ("DEPTH") containing the depth coordinates from the main dataset.
    var : string
        Variable, where the deviation should be calculated.
    min_dep : int or float
        Minimum depth for the fit's depth range.
    max_dep : int or float
        Maximum depth for the fit's depth range.
    k : int or float, optional
        Degree of the smoothing spline. The default is None.
    s : int or float, optional
        Smoothing paramter used to construct the spline function. The default is 0.005.

    Returns
    -------
    df_fit : pandas dataframe
        Datframe containing the depth and the fit result.

    """
    import pandas as pd
    import numpy as np
    from scipy.interpolate import UnivariateSpline

    #calculate the univariate spline fit
    fit_func = UnivariateSpline(bg['DEPTH'], bg[var],k=k,s=s)
    res = fit_func.get_residual()
    print('Residual UnivariateSpline: '+str(res))
    fit = fit_func(dep_vec)

    df = {'DEPTH':dep_vec['DEPTH'],'Bgfit_'+var:np.ravel(fit)}
    df_fit = pd.DataFrame(df)
    df_fit = df_fit.sort_values(by='DEPTH',ascending=True).dropna(subset=['DEPTH'])

    return df_fit

def calc_delta_by_bgfit(data, bg, var, min_dep, max_dep, fit, param, control_plot=False):
    """
    This functions calculates the deviation of a variable from a polynomial fit of the variable from a different (background) station in a specified depth range. For this, first a polynomial fit of the background variable is calculated in the depth range and then merged with the main dataset based on similar depth values. Then this fit is substracted from the variable's values in the main dataset.

    Parameters
    ----------
    data : pandas dataframe or dictionary
        Dataset containing one station. Must have the variables depth ("DEPTH") and datetime ("datetime").
    bg : pandas dataframe
        Dataset which should be used as background.
    var : string
        Variable, where the deviation should be calculated.
    min_dep : int
        Minimum depth for the fit's depth range. It is advisable to set this to the region of interest to produce a good fit.
    max_dep : int
        Maximum depth for the fit's depth range. It is advisable to set this to the region of interest to produce a good fit.
    fit : string
        Keyword for the fit method used. Can be "poly" for polynomial fit (np.polyfit) or "uni" for univariate spline fit (scipy.interpolate.UnivariateSpline).
    param : int or tuple
        If the poly fit is used: Order of the polynomial fit.
        If the univariate Spline is used: Tuple of k (degree of the smoothing spline) and s (smoothing factor)
    control_plot : boolean, optional
        This option controls if a control plot is shown, containing the dataset variable, the background variable and the calculated fit in the specified depth range. The default is False.

    Returns
    -------
    data : pandas dataframe
        Dataframe similar to the original dataframe with the background fit and the deviation of the variable as additional columns. Values outside the depth range in these columns are NaN.
    """
    import pandas as pd

    #check if there is the same serial number in background data, if not then warn
    if len(bg['SN'].unique()) > 1:
        sn = data['SN'].iloc[0]
        if sn in bg['SN']:
            bg = bg[bg['SN'] == sn]
        else:
            sn = bg['SN'].iloc[0]
            bg = bg[bg['SN'] == sn]
            print('Warning: No profile with the same serial number found in the background. Falling back to any other serial number: '+sn)

    #for continous data, the fit based on the background is calculated, merged with the dataset and subtracted
    data['datetime'] = pd.to_datetime(data['datetime'])
    dep_vec = data[['DEPTH']]
    dep_vec = dep_vec[(dep_vec['DEPTH']>=min_dep) & (dep_vec['DEPTH']<=max_dep)]
    #dep_vec = sep_casts(dep_vec)[0]
    dep_vec = dep_vec.sort_values(by='DEPTH')

    #get first downcast of background
    bg = sep_casts(bg, window_size=5000)
    if len(bg) > 2:
        print('Warning: Your background station was separated into more than two up or down casts. Check the casts returned by sep_casts and the used window size.')
    bg = bg[0]
    bg = bg[(bg['DEPTH']>=min_dep) & (bg['DEPTH']<=max_dep)]
    bg = bg.sort_values(by='DEPTH')
    bg = bg[bg[var].notna()]

    if fit == 'poly':
        bg_fit = get_bg_polyfit(bg, dep_vec, var, min_dep, max_dep, param)
    if fit == 'uni':
        bg_fit = get_bg_unifit(bg, dep_vec, var, min_dep, max_dep, param[0],param[1])

    #data maybe needs droppping DEPTH nan rows?
    data = data.sort_values(by='DEPTH')
    data = data.reset_index().merge(bg_fit.drop_duplicates(subset='DEPTH'), how='left', on='DEPTH').set_index('index')
    data = data.sort_index()
    data['Delta_'+var] = data[var] - data['Bgfit_'+var]

    #control plot, showing the variable in the dataset, the variable in the background and the calculated fit
    if control_plot == True:
        import matplotlib.pyplot as plt
        from hyvent.misc import get_var

        plt.figure()

        plt.plot(data[var],data['DEPTH'], label=get_var(var)[0])
        plt.plot(bg[var],bg['DEPTH'], label='Var from background')
        plt.plot(data['Bgfit_'+var],data['DEPTH'], label='Fit')

        plt.gca().invert_yaxis()
        plt.xlabel(get_var(var)[0])
        plt.ylabel('Depth in m')
        plt.ylim((max_dep,min_dep))
        #get x limits
        data_cut = pd.concat([data,bg])
        data_cut = data_cut[(data_cut['DEPTH']>=min_dep) & (data_cut['DEPTH']<=max_dep)]
        plt.xlim((data_cut[var].min()-abs(data_cut[var].min()/50),data_cut[var].max()+abs(data_cut[var].max()/50)))
        plt.title(fit+' s/n='+str(param)+'; Station = '+str(data['Station'].iloc[0])+', '+str(data['SN'].iloc[0]))
        plt.legend()
        plt.show()

    return data

def calc_delta_stafit(data, var, lim_above, lim_below, fit, param, control_plot=False):
    """
    [Not working yet] This functions calculates the deviation of a variable from a polynomial fit of the variable from a different (background) station in a specified depth range. For this, first a polynomial fit of the background variable is calculated in the depth range and then merged with the main dataset based on similar depth values. Then this fit is substracted from the variable's values in the main dataset. If the variable is "delta3He", that means only discrete samples, just a mean is calculated in the depth range and substracted from the dataset values.

    Parameters
    ----------
    data : pandas dataframe or dictionary
        Dataset containing one station. Must have the variables depth ("DEPTH") and datetime ("datetime").
    var : string
        Variable, where the deviation should be calculated.
    min_dep : int
        Minimum depth for the fit's depth range. It is advisable to set this to the region of interest to produce a good fit.
    max_dep : int
        Maximum depth for the fit's depth range. It is advisable to set this to the region of interest to produce a good fit.
    fit : string
        Keyword for the fit method used. Can be "poly" for polynomial fit (np.polyfit) or "uni" for univariate spline fit (scipy.interpolate.UnivariateSpline).
    param : int or tuple
        If the poly fit is used: Order of the polynomial fit.
        If the univariate Spline is used: Tuple of k (degree of the smoothing spline) and s (smoothing factor)
    control_plot : boolean, optional
        This option controls if a control plot is shown, containing the dataset variable, the background variable and the calculated fit in the specified depth range. The default is False.

    Returns
    -------
    data : pandas dataframe
        Dataframe similar to the original dataframe with the background fit and the deviation of the variable as additional columns. Values outside the depth range in these columns are NaN.
    """
    import pandas as pd

    #for continous data, the fit based on the background is calculated, merged with the dataset and subtracted
    data['datetime'] = pd.to_datetime(data['datetime'])
    dep_vec = data[['DEPTH']]
    dep_vec = dep_vec[(dep_vec['DEPTH']>=lim_above[1]) & (dep_vec['DEPTH']<=lim_below[0])]
    #dep_vec = sep_casts(dep_vec)[0]
    dep_vec = dep_vec.sort_values(by='DEPTH')

    above = data[(data['DEPTH']>=lim_above[0]) & (data['DEPTH']<=lim_above[1])]
    below = data[(data['DEPTH']>=lim_below[0]) & (data['DEPTH']<=lim_below[1])]
    data_fit = pd.concat([above,below])
    data_fit = data_fit.sort_values(by='DEPTH')
    data_fit = data_fit[data_fit[var].notna()]

    if fit == 'poly':
        bg_fit = get_bg_polyfit(data_fit, dep_vec, var, lim_above[1], lim_below[0], param)
    if fit == 'uni':

        bg_fit = get_bg_unifit(data_fit, dep_vec, var, lim_above[1], lim_below[0], param[0],param[1])

    #data maybe needs droppping DEPTH nan rows?
    data = data.sort_values(by='DEPTH')
    data = data.reset_index().merge(bg_fit.drop_duplicates(subset='DEPTH'), how='left', on='DEPTH').set_index('index')
    data = data.sort_index()
    data['Delta_'+var] = data[var] - data['Bgfit_'+var]

    #control plot, showing the variable in the dataset, the variable in the background and the calculated fit
    if control_plot == True:
        import matplotlib.pyplot as plt
        from hyvent.misc import get_var

        plt.figure()

        plt.plot(data[var],data['DEPTH'], label=get_var(var)[0])
        plt.scatter(data_fit[var],data_fit['DEPTH'], label='Fit data',s=0.1)
        plt.plot(data['Bgfit_'+var],data['DEPTH'], label='Fit')

        plt.gca().invert_yaxis()
        plt.xlabel(get_var(var)[0])
        plt.ylabel('Depth in m')
        plt.ylim((lim_below[1]),lim_above[0])
        #get x limits
        data_cut = data
        data_cut = data_cut[(data_cut['DEPTH']>=lim_above[1]) & (data_cut['DEPTH']<=lim_below[0])]
        plt.xlim((data_cut[var].min()-abs(data_cut[var].min()/50),data_cut[var].max()+abs(data_cut[var].max()/50)))
        plt.title(fit+' s/n='+str(param)+'; Station = '+str(data['Station'].iloc[0])+', '+str(data['SN'].iloc[0]))
        plt.legend()
        plt.show()

    return data

def calc_helium_delta(data, bg, var, min_dep, max_dep, control_plot=False):
    """
    This functions calculates the deviation of a delta3He from the mean of a different (background) station in a specified depth range. For this, the mean of the background variable is calculated and substracted from the variable's values in the main dataset.

    Parameters
    ----------
    data : pandas dataframe or dictionary
        Dataset containing one station. Must have the variable depth ("DEPTH").
    bg : pandas dataframe
        Dataset which should be used as background.
    var : string
        Variable, where the deviation should be calculated.
    min_dep : int
        Minimum depth for the depth range where the mean is calculated. It is advisable to exclude surface measurements.
    max_dep : int
        Maximum depth for the depth range where the mean is calulated.
    control_plot : boolean, optional
        This option controls if a control plot is shown, containing the dataset variable, the background variable and the calculated fit in the specified depth range. The default is False.

    Returns
    -------
    data : pandas dataframe
        Dataframe similar to the original dataframe with the bakcground's mean and the deviation as additional columns. Values outside the depth range in these columns are NaN.
    """
    import pandas as pd

    #for delta3He, data is bottle data, therefore take only mean and subtract it
    bg = bg.rename({'DepSM_mean':'DEPTH'},axis=1)
    bg = bg[(bg['DEPTH']>=min_dep) & (bg['DEPTH']<=max_dep)]
    bg_mean = bg['delta3He'].mean()
    data = data.rename({'DepSM_mean':'DEPTH'},axis=1)
    data['Bgmean_'+var] = bg_mean
    data['Delta_'+var] = data[var] - data['Bgmean_'+var]

    #control plot, showing the variable in the dataset and the variable in the background
    if control_plot == True:
        import matplotlib.pyplot as plt
        from hyvent.misc import get_var

        plt.figure()
        if var == 'delta3He':
            plt.scatter(data[var],data['DEPTH'], label=get_var(var)[0])
            plt.scatter(bg[var],bg['DEPTH'], label='Var from background')
            plt.vlines(bg_mean,min_dep,max_dep, label='Mean from background')

        plt.gca().invert_yaxis()
        plt.xlabel(get_var(var)[0])
        plt.ylabel('Depth in m')
        plt.ylim((max_dep,min_dep))
        #get x limits
        data_cut = pd.concat([data,bg])
        data_cut = data_cut[(data_cut['DEPTH']>=min_dep) & (data_cut['DEPTH']<=max_dep)]
        plt.xlim((data_cut[var].min()-abs(data_cut[var].min()/50),data_cut[var].max()+abs(data_cut[var].max()/50)))
        plt.legend()
        plt.show()

    return data

def calc_turb_delta(data, var, upper_layer, lower_layer):
    """
    This functions calculates the deviation of turbidity. This is done by excluding a plume layer and defining a layer above and below the plume, where the mean is calculated and then subtracted.

    Parameters
    ----------
    data : pandas dataframe or dictionary
        Dataset containing one station. Must have the variable depth ("DEPTH").
    var : string
        Variable, where the deviation should be calculated.
    upper_layer : tuple of float
        Upper and lower depth of the layer above the plume, which is is used as the background.
    lower_layer : tuple of float
        Upper and lower depth of the layer below the plume, which is is used as the background.

    Returns
    -------
    data : pandas dataframe
        Dataframe similar to the original dataframe with the the deviation as an additional column. Values outside the depth range in these columns are NaN.
    """
    import numpy as np

    above = data[(data['DEPTH']<upper_layer[1]) & (data['DEPTH']>upper_layer[0])]
    below = data[(data['DEPTH']<lower_layer[1]) & (data['DEPTH']<lower_layer[0])]
    mean = np.mean((above[var].mean(),below[var].mean()))

    data['Delta_'+var] = data[var] - mean

    return data

def zmax_theta(data, min_dep, vent_depth, uncert, control_plot=False):
    """
    This function finds the maximum rise height of a plume based on the potential temperature anomaly. That is the vent's depth minus the minimal depth below a limit ("min_dep") where the anomaly is higher then 0 + a threshold ("uncert"). Optionally a control plot is printed, which shows the profile and and zmax.

    Parameters
    ----------
    data : pandas dataframe
        Dataset of one or several profiles. Must contain the potential temperature anomaly ("Delta_potemperature") and depth ("DEPTH").
    min_dep : int
        Minimum depth below which data is exmained for zmax.
    vent_depth : int
        Depth of the hydrothermal vent or the seafloor.
    thres: float
        Threshold of the potential temperature which is added to zero. This defines the threshold a value has to pass to be considered for zmax.
    control_plot : boolean, optional
        Boolian which controlles if control plots are produced or not. The default is False.

    Returns
    -------
    zmax : float
        The maximum rise height of the plume above the vent depth.
    """
    import matplotlib.pyplot as plt
    from hyvent.misc import get_var

    var = 'Delta_potemperature'

    #Find the minimum depth where the anomaly is larger then 0+the uncertainty and below a minimum depth
    zmax = data['DEPTH'][(data[var]>0+uncert) & (data['DEPTH'] > min_dep)].min()
    #print(vent_loc[2]-zmax)

    #control plot with profile and zmax
    if control_plot == True:
        plt.figure()

        data_list = [d for _, d in data.groupby(['Station','SN'])]
        for station in data_list:
            plt.plot(station[var],station['DEPTH'],color=get_var(var)[1])
        plt.axhline(zmax)
        plt.axvline(0, color = 'black',alpha=0.3)
        plt.gca().invert_yaxis()
        plt.show()

    #calculate the rise height as height above the vent
    zmax = vent_depth - zmax

    return zmax

def calc_flux(zmax, depth_vent, ref_station, const=None):
    """
    This function calculates the estimates the heatflux of a hydrothermal vent based on the rise height of the plume. The calculation is based on McDougall, T. J. (1990) and implemented as shown in the supplemental material of Wegener, G.(2024).

    Parameters
    ----------
    zmax : int or float
        Maximum rise height of the plume, as zmax = depth_vent - minimum_plume_depth.
    depth_vent : int or float
        Depth of the vent site.
    ref_station : pandas dataframe
        Dataframe containing the potential density anomaly and depth of a CTD station. This density profile is used to calculate the buoyancy frequency N^2.
    const : dictionary, optional
        Dictionary of physical constants needed for the calculation. Should be in the form of "const = {'rho0' : 1000, 'cp' :  3890, 'g' : 9.8, 'alpha' : 1.3 * 10**(-4)}". If none is provided, the default values are used, which are taken from supplemental material of Wegener, G.(2024). The default is None.

    Returns
    -------
    heat_flux : float
        The estimated heat flux of the vent site.
    Nsquared : float
        The buoyancy frequency of the stratification based on the reference station.
    """
    from hyvent.processing import sep_casts
    import numpy as np

    # default constants
    const = {'rho0' : 1000, 'cp' :  3890, 'g' : 9.8, 'alpha' : 1.3 * 10**(-4)} #, 'theta_bg' : 0}
    # paramers from Wegener, 2024 can be used for testing
    # N_gw = np.sqrt(0.25 * 10**(-7))
    # zmax_gw = 800

    # assign constants
    rho0 = const['rho0']
    cp = const['cp']
    g = const['g']
    alpha = const['alpha']
    #theta_bg = const['theta_bg']  #needed for calculation of Qv

    #get the first downcast of the reference station for calculating N
    N_cast = sep_casts(ref_station)[0]

    #get rho at the vent depth and rho at zmax depth
    rho_vent = N_cast['Sigma3'].iloc[(N_cast['DEPTH']-depth_vent).abs().argsort()[:1]]
    rho_vent = np.float64(rho_vent.iloc[0])
    rho_zmax = N_cast['Sigma3'].iloc[(N_cast['DEPTH']-(depth_vent-zmax)).abs().argsort()[:1]]
    rho_zmax = np.float64(rho_zmax.iloc[0])

    #calculate buoyancy frequency and heat flux
    Nsquared = (- g/rho0 * (rho_vent-rho_zmax)/-zmax).astype(float)
    N = np.sqrt(Nsquared)
    heat_flux = rho0 * cp * np.pi * N**3 *(zmax/5)**4 / g / alpha

    #calculate volume flux [not working]
    #Qv = heat_flux/(cp * rho0 * (theta_vent - theta_bg))

    return heat_flux, Nsquared

def get_fit_cast(data_fit, dens_var, min_dep, window_size, control_plot=False):
    """
    This function returns a list of all up- and down cast

    Parameters
    ----------
    data_fit : TYPE
        DESCRIPTION.
    dens_var : TYPE
        DESCRIPTION.
    min_dep : TYPE
        DESCRIPTION.
    control_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    from hyvent.processing import sep_casts
    from hyvent.misc import get_var
    import matplotlib.pyplot as plt

    cast_list = sep_casts(data_fit, window_size=window_size)

    for i, cast in enumerate(cast_list):
        cast_list[i] = cast[cast['DEPTH']>min_dep]

    #remove stationary part in station 028_01
    if data_fit['Station'].iloc[0] == '028_01':
        cast_list = cast_list[:8] + cast_list[17:]




    #plot all up and down casts, to choose one for the fit
    if control_plot == True:
        plt.figure()
        #cycle through colors:
        cm = plt.get_cmap('tab20')
        num_col = len(cast_list)
        plt.gca().set_prop_cycle('color', [cm(1.*i/num_col) for i in range(num_col)])

        for i, cast in enumerate(cast_list):
            plt.plot(cast[dens_var],cast['potemperature'],label=str(i),marker='.')

            plt.gca().invert_yaxis()
            plt.ylabel(get_var(dens_var)[0])
            plt.xlabel(get_var('potemperature')[0])
            plt.legend()
            plt.title(data_fit['Station'].iloc[0]+'; '+data_fit['SN'].iloc[0])
        plt.tight_layout()
        plt.show()


    # #take input from console for which cast should be used for fitting
    # plt.pause(0.1)
    # print('Type the number of the cast that should be used for fitting. Figure will be closed then:')
    # cast_no = int(input())
    # plt.close()

    # fit_cast = cast_list[cast_no]

    return(cast_list)

def calc_delta_densfit(data, dens_var, dens_cut, min_dep, fit_cast, fit_order=3, control_plot=False):

    from hyvent.misc import get_var
    from numpy.polynomial import polynomial as poly
    import matplotlib.pyplot as plt
    import pandas as pd

    #interpolate NaNs dan drop NaN rows
    if fit_cast[dens_var].isna().sum()>0:
        print('Interpolating '+str(fit_cast[dens_var].isna().sum())+' NaN values...')
        fit_cast[dens_var] = fit_cast[dens_var].interpolate()
        fit_cast = fit_cast.dropna(subset=[dens_var])

    #for station 028_01 all cast have a anomaly, therfore cut this part and only fit the rest
    fit_cast_part = fit_cast[(fit_cast[dens_var] < dens_cut)]
    coef = poly.polyfit(fit_cast_part[dens_var],fit_cast_part['potemperature'],fit_order)
    fit_cast['Fit'] = poly.polyval(fit_cast[dens_var],coef)

    # else:
    #     #do poly fit
    #     coef = poly.polyfit(fit_cast[dens_var],fit_cast['potemperature'],fit_order)
    #     fit_cast['Fit'] = poly.polyval(fit_cast[dens_var],coef)

    if control_plot == True:
        #plot data and fit
        plt.figure()
        plt.scatter(fit_cast['potemperature'],fit_cast[dens_var],label='Data',s=0.5,color=get_var('potemperature')[1])
        plt.plot(fit_cast['Fit'], fit_cast[dens_var],label='Fit',color='black')
        plt.gca().invert_yaxis()
        plt.ylabel(get_var(dens_var)[0])
        plt.xlabel(get_var('potemperature')[0])
        plt.legend()
        plt.title(fit_cast['Station'].iloc[0]+'; '+fit_cast['SN'].iloc[0])
        plt.tight_layout()
        plt.show()

    #merge with dataset
    fit_cast = fit_cast.sort_values(by='DEPTH')
    data['DEPTH'] = data['DEPTH'].interpolate()
    data = data.sort_values(by='DEPTH')
    data = pd.merge_asof(data, fit_cast[['DEPTH','Fit']], on='DEPTH', direction='nearest', tolerance=1)
    data = data.sort_index()

    if control_plot == True:
        #plot all data and fit result vs depth
        plt.figure()
        plt.tight_layout()
        plt.scatter(data['potemperature'][data['DEPTH']>min_dep],data['DEPTH'][data['DEPTH']>min_dep],label='Data',s=0.5,color = get_var('potemperature')[1])
        plt.plot(data['Fit'][data['DEPTH']>min_dep],data['DEPTH'][data['DEPTH']>min_dep],label='Fit evaluation',color='black')
        plt.gca().invert_yaxis()
        plt.ylabel(get_var('DEPTH')[0])
        plt.xlabel(get_var('potemperature')[0])
        plt.legend()
        plt.title(data['Station'].iloc[0]+'; '+data['SN'].iloc[0])
        plt.tight_layout()
        plt.show()

    # plt.pause(1)
    # print('Enter "y" to use this fit and continue:')

    # plt.close('all')

    data['Delta_potemperature'] = data['potemperature'] - data['Fit']
    del data['Fit']

    return data

