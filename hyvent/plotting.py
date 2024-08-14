#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:25:58 2024

@author: jonathan
"""

def plot_hist_he(data,lower_depth,upper_depth,bins,path_save='None'):
    import pandas as pd
    import matplotlib.pyplot as plt

    label_sel = 'Samples between '+str(upper_depth)+' - '+str(lower_depth)+'m'

    he_data = pd.concat(data.values(),ignore_index=True)
    plt.hist(he_data['delta3He'],bins=bins,label='All samples')
    plt.hist(he_data['delta3He'][(he_data['DepSM_mean']>upper_depth) & (he_data['DepSM_mean']<lower_depth)],bins=bins,label=label_sel,alpha=0.7)

    plt.xlabel('Delta 3He in %')
    plt.ylabel('Number of samples')
    plt.legend()

    if path_save != 'None':
        plt.savefig(path_save, dpi=300)

    plt.show()

def plot_hist_dorp(data,lower_depth,upper_depth,bins,ranges,path_save='None'):
    import pandas as pd
    import matplotlib.pyplot as plt

    label_sel = 'Samples between '+str(upper_depth)+' - '+str(lower_depth)+'m'

    dorp_data = pd.concat(data.values(),ignore_index=True)
    plt.hist(dorp_data['dORP'],bins=bins,range=ranges,label='All samples',log=True)
    plt.hist(dorp_data['dORP'][(dorp_data['DEPTH']>upper_depth) & (dorp_data['DEPTH']<lower_depth)],bins=bins,range=ranges,label=label_sel,log=True)

    plt.xlabel('dORP')
    plt.ylabel('Number of samples')
    plt.legend()


