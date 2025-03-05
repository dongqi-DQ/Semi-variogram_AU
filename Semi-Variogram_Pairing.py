#!/usr/bin/env python
# coding: utf-8

# In[1]:


import xarray as xr
import numpy as np
import os
from glob import glob
from tqdm import tqdm
import pandas as pd
from geopy.distance import geodesic 
import math
from sklearn.metrics.pairwise import haversine_distances
from math import radians
from multiprocess import Pool
import warnings


# In[2]:


from gstools import CovModel

class Stab(CovModel):
    def variogram(self, r):
        
        return self.nugget + self.sill * (1 - np.exp(-(3 * r) / self.len_scale))


# In[3]:


from scipy.optimize import curve_fit

# Define the piecewise function
def custom_curve(h, c, b, alpha):
    if np.isscalar(h):
        # Handle the scalar case
        if h == 0:
            return 0
        else:
            return c + b * (1 - np.exp(-3 * h / alpha))
    else:
        # Handle array inputs
        return np.where(h == 0, 0, c + b * (1 - np.exp(-3 * h / alpha)))

def fit_sci_curve(h_values, y_values,bins):

    # Use curve_fit to fit the custom function to the data
    # Initial guess for c, b, and alpha
    initial_guess = [0, 0.5, 5]  
    
    # Perform the curve fitting
    params, covariance = curve_fit(custom_curve, h_values, y_values, p0=initial_guess)
    
    # Extract the fitted parameters
    c_fitted, b_fitted, alpha_fitted = params
    
    y_fit = custom_curve(bins, c_fitted, b_fitted, alpha_fitted)
    return y_fit, alpha_fitted


# In[4]:


def calc_ratio_gamma(N11,N10):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        ratio = np.where(N10 + N11 != 0, 0.5 * (N10 / (N10 + N11)), np.nan)
        # if len(ratio.shape)==2:
        #     ratio[:,0] = 0
        # if len(ratio.shape)==1:
        #     ratio[0] = 0    
    return ratio


# In[5]:


# Define distance bins
bins = np.arange(5, 360, 10)

# Function to get pairwise distances
# including the station at the centre
def get_pairwise_distances(df):
    spec_stn = df["Spec_stn"].iloc[0]
    stn_list = df["Neighb_stn"].values
    stn_list = np.append(spec_stn, stn_list)
    stn_str = [str(ids) for ids in stn_list]
    distances =  df_dist.loc[stn_list, stn_str].values

    return distances

# Function to calculate N11 and N10 pairs
def calculate_pairs(df, bins):
    distances = get_pairwise_distances(df)
    n = len(df) + 1  # Include the station at the neighborhood center
    extreme = np.append(90, df['flag'].values)  # Append the center station's flag

    # Initialize arrays for N11 and N10
    N11 = np.zeros_like(bins, dtype=int)
    N10 = np.zeros_like(bins, dtype=int)

    # Precompute bin indices for all pairs
    bin_indices = np.digitize(distances, bins) - 1

    # Create masks for valid bin indices
    valid_mask = (bin_indices >= 0) & (bin_indices < len(bins))

    # Create masks for 1-1 and 1-0 pairs
    mask_11 = (extreme[:, None] > 0) & (extreme[None, :] > 0)  # 1-1 pairs
    mask_10 = ((extreme[:, None] > 0) & (extreme[None, :] == 0)) | ((extreme[:, None] == 0) & (extreme[None, :] > 0))  # 1-0 pairs

    # Iterate over all bins
    for bin_idx in range(len(bins)):
        # Find pairs that fall into this bin
        bin_mask = (bin_indices == bin_idx) & valid_mask

        # Count 1-1 and 1-0 pairs in this bin
        N11[bin_idx] = np.sum(mask_11 & bin_mask)
        N10[bin_idx] = np.sum(mask_10 & bin_mask)

    return N11/2, N10/2


# In[6]:


def calc_alpha(gamma_fit,bins_fit,bins):
    # gamma_fit = gamma[~np.isnan(gamma)]
    # bins_fit = bins[~np.isnan(bins)]
    y_fit = 0
    scale = 0
    try:
        _ = fit_model.fit_variogram(bins_fit, gamma_fit,  nugget=0, sill=0.5)#)#
        y_fit = fit_model.variogram(bins)
        scale = fit_model.len_scale
    except RuntimeError:
        try: 
            ## adjust the cape a little bit so the code can pass through
            _ = fit_model.fit_variogram(bins_fit, gamma_fit,  nugget=0, sill=0.51)#)#
            y_fit = fit_model.variogram(bins)
            scale = fit_model.len_scale
        except RuntimeError:
            y_fit,alpha = fit_sci_curve(bins_fit,gamma_fit,bins)
            # no_fit_days.append(days)
            scale = alpha
    return y_fit, scale


# In[7]:


def calc_pair_alpha(days):
    df = df_stations_p90[df_stations_p90["Day"]==days]
    if df["flag"].max() <=0:
        N11 = np.zeros_like(bins)
        N10 = np.zeros_like(bins)
        gamma = np.zeros_like(bins)
        y_fit = np.zeros_like(bins)
        scale = 0
    else:
        
        # Calculate N11 and N10 pairs
        N11, N10 = calculate_pairs(df, bins)
        gamma = calc_ratio_gamma(N11, N10)        
        gamma_fit = gamma.copy()
        bins_fit = bins.copy().astype(float)
        gamma_fit[N11+N10<=2] = np.nan
        bins_fit[N11+N10<=2] = np.nan
        ## make sure start with zeros
        # gamma_fit[0] = 0
        # bins_fit[0] =0
        ## remove NaNs 
        gamma_fit1 = gamma_fit[~np.isnan(gamma_fit)]
        bins_fit1 = bins_fit[~np.isnan(bins_fit)]
        ## in case all NaNs occured
        if len(gamma_fit1)<=2:
            y_fit = np.zeros_like(bins)
            scale = 0
        else:
            y_fit, scale = calc_alpha(gamma_fit1, bins_fit1, bins)
    local_bins_dict = {
        "Day": [days] * len(bins),
        "Date": [df_stations_p90["Date"][df_stations_p90["Day"]==days].values[0]] *len(bins),
        "Bins": bins.tolist(),
        "N11": N11.tolist(),
        "N10": N10.tolist(),
        "gamma": gamma.tolist(),
        "y_fit": y_fit.tolist()
    }
    return local_bins_dict, scale


# In[8]:


## my own directory
os.chdir("/g/data/k10/dl6968/Semi-variogram_AU/")


# In[10]:


fit_model = Stab(dim=2)#gs.Stable(dim=1)


# In[9]:


get_ipython().run_cell_magic('time', '', 'df_dist = pd.read_csv("./data/pairwise_distances.csv", index_col=0)\n')


# In[11]:


files = sorted(glob("./data/all_AU_p90/*_station_list_all_events.csv"))


# In[12]:


max_pool = 28
for file in tqdm(files, leave=True, position=0, total=len(files)):
    if not os.path.exists(file.replace("station","pair_scale")):
        df_stations_p90 = pd.read_csv(file)
        with Pool(max_pool) as p:
            pool_outputs = list(
                    p.imap(calc_pair_alpha,
                           np.unique(df_stations_p90["Day"])),
            )
        p.join()
        
        y_fit_arr = np.zeros((len(np.unique(df_stations_p90["Day"])), bins.shape[0]))
        N11_arr = np.zeros_like(y_fit_arr)
        N10_arr = np.zeros_like(y_fit_arr)
        gamma_arr = np.zeros_like(y_fit_arr)
        scale_arr = np.zeros(len(np.unique(df_stations_p90["Day"])))
        bins_dict = {"Day": [], "Date": [], "Bins": [], "N11": [], "N10": [], "gamma": [], "y_fit": []}
        for i in range(0, len(np.unique(df_stations_p90["Day"]))):
            scale_arr[i] = pool_outputs[i][1]
            for key in bins_dict:
                    bins_dict[key].extend(pool_outputs[i][0][key])
            
        scale_dict = {"extreme_dates": np.unique(df_stations_p90["Date"]), "scale": scale_arr}
        df_scale = pd.DataFrame.from_dict(scale_dict)
        df_scale.to_csv(file.replace("station","pair_scale"))
        
        df_bins =  pd.DataFrame.from_dict(bins_dict)
        df_bins.to_csv(file.replace("station","pair_bins"))
    


# In[ ]:




