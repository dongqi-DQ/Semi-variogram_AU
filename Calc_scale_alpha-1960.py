# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.2
#   kernelspec:
#     display_name: Python [conda env:myenv]
#     language: python
#     name: conda-env-myenv-py
# ---

# %%
import gstools as gs
from gstools import CovModel
import numpy as np
from IPython.display import clear_output
import pandas as pd

from tqdm import tqdm
import xarray as xr
from glob import glob
from multiprocess import Pool

from scipy.optimize import curve_fit
import os 
os.chdir("/g/data/k10/dl6968/Semi-variogram_AU/")


# %%
class Stab(CovModel):
    def variogram(self, r):
        
        return self.nugget + self.sill * (1 - np.exp(-(3 * r) / self.len_scale))


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
    initial_guess = [0.025, 0.5, 5]  
    
    # Perform the curve fitting
    params, covariance = curve_fit(custom_curve, h_values, y_values, p0=initial_guess)
    
    # Extract the fitted parameters
    c_fitted, b_fitted, alpha_fitted = params
    
    y_fit = custom_curve(bins, c_fitted, b_fitted, alpha_fitted)
    return y_fit, alpha_fitted


# %%
def calc_alpha(gamma_fit,bins_fit,bins):
    y_fit = np.zeros_like(bins)
    scale = 0
    try:
        fit_model = Stab(dim=2)  # Choose model type (Exponential, Spherical, Gaussian, etc.)
        _ = fit_model.fit_variogram(bins_fit, gamma_fit, nugget=0.025,sill=0.5-0.025)#)#
        y_fit = fit_model.variogram(bins)
        scale = fit_model.len_scale
        
    except RuntimeError:
        ## in case the variogram fit did not work
        ## this is usually a very small case
        ## mark it with -1 for later
        y_fit = np.zeros_like(bins)
        scale = -1
    return y_fit, scale


# %%
def gamma_to_alpha(args):
    N11, N10, gamma, bins = args
    gamma_fit = gamma.copy()
    bins_fit = bins.copy().astype(float)
    gamma_fit[N11+N10<=2] = np.nan
    bins_fit[N11+N10<=2] = np.nan
    
    ## remove NaNs 
    gamma_fit1 = gamma_fit[~np.isnan(gamma_fit)]
    bins_fit1 = bins_fit[~np.isnan(bins_fit)]
    
    ## in case all NaNs occured
    if len(gamma_fit1)<=2:
        y_fit = np.zeros_like(bins)
        scale = 0
    else:
        ## make sure this starts with zeros
        # gamma_fit1 = np.append(0, gamma_fit1)
        # bins_fit1 = np.append(0, bins_fit1)
        y_fit, scale = calc_alpha(gamma_fit1, bins_fit1, bins)
   
    return y_fit, scale


# %%
files = sorted(glob("./data/seasonal_p90_500km_1960/*pair_bins*.csv"))

# %%
for file in tqdm(files, position=0, leave=True):
    if not os.path.exists(file.replace("bins","scale_spherical")):
        df_bins = pd.read_csv(file, usecols=lambda col: not col.startswith("Unnamed"))
        # df_bins = df_bins[df_bins["Bins"]<=501].reset_index()
        args_list = [(df_bins[df_bins["Day"]==day]["N11"].values,df_bins[df_bins["Day"]==day]["N10"].values,
                  df_bins[df_bins["Day"]==day]["gamma"].values,df_bins[df_bins["Day"]==day]["Bins"].values,) for day in np.unique(df_bins["Day"])]
        ## multi-processing the extreme days
        max_pool =28
        with Pool(max_pool) as p:
            pool_outputs = list(
                    p.imap(gamma_to_alpha,
                           args_list),
            )
        p.join()
        ## save output
        y_fit_list = []
        scale_list = []
        for output in pool_outputs:
            y_fit_list.append(output[0])
            scale_list.append(output[1])
        ## save to csv
        df_bins["y_fit"] = np.concatenate(y_fit_list)
        df_bins.to_csv(file)
        scale_dict = {"extreme_dates": np.unique(df_bins["Date"]), "scale": scale_list}
        df_scale = pd.DataFrame.from_dict(scale_dict)
        df_scale.to_csv(file.replace("bins","scale_spherical"))
        
        clear_output(wait=True)

# %%

# %% [markdown]
# removed 49084, 78046, 78047
