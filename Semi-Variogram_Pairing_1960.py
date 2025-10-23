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
from IPython.display import clear_output
import warnings

# Suppress all warnings
warnings.filterwarnings('ignore', category=FutureWarning)


# %%
def calc_ratio_gamma(N11,N10):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        ratio = np.where(N10 + N11 != 0, 0.5 * (N10 / (N10 + N11)), np.nan)
    return ratio

# Function to get pairwise distances
# including the station at the centre
def get_pairwise_distances(df):
    spec_stn = df["Spec_stn"].iloc[0]
    stn_list = df["Neighb_stn"].values
    stn_list = np.append(spec_stn, stn_list)
    stn_str = [str(ids) for ids in stn_list]
    distances =  df_dist.loc[stn_list, stn_str].values

    return distances

def calculate_moving_pairs(args):
    df_nb, bins, days = args
    df = df_nb[df_nb['distance']<=500]
    if len(df[df["val"]>=0]) == 0:
        N11 = np.zeros_like(bins)
        N10 = np.zeros_like(bins)
        return N11, N10
    else:
        ## must have more than 10% stations had extremes
        if df["flag"].max() <=0 or len(df[df["flag"]>0])/len(df[df["val"]>=0])<0.1 or len(df)<20:
            N11 = np.zeros_like(bins)
            N10 = np.zeros_like(bins)
            return N11, N10
        else: 
            # Append the center station's flag
            in_out_arr = np.append(1, np.where(df_nb['distance'] > 500, 0, 1))
           
            ## inside is 1
            ## outside is 0
            pairwise_sums = in_out_arr[:, None] + in_out_arr[None, :]
            distances = get_pairwise_distances(df_nb)
            n = len(df_nb) + 1  # Include the station at the neighborhood center
            extreme = np.append(1, df_nb['flag'].values)  # Append the center station's flag
            extreme[extreme>0]=1
            pairwise_extreme = extreme[:, None] + extreme[None, :]
            
            # Initialize arrays for N11 and N10
            N11 = np.zeros_like(bins, dtype=int)
            N10 = np.zeros_like(bins, dtype=int)
            
            # Precompute bin indices for all pairs
            bin_indices = np.digitize(distances, bins) - 1
            
            # Create masks for valid bin indices
            valid_mask = (bin_indices >= 0) & (bin_indices < len(bins))
            
            # Create masks for 1-1 and 1-0 pairs
            mask_11 = (pairwise_extreme==2) & (pairwise_sums>0)  # 1-1 pairs for at least one station inside
            mask_10 = (pairwise_extreme==1) & (pairwise_sums>0)  # 1-0 pairs for at least one station inside
            pair_mask = pairwise_sums>0
            # Iterate over all bins
            for bin_idx in range(len(bins)):
                # Find pairs that fall into this bin
                bin_mask = (bin_indices == bin_idx) & valid_mask
            
                # Count 1-1 and 1-0 pairs in this bin
                N11[bin_idx] = np.sum(mask_11  & bin_mask) # & pair_mask
                N10[bin_idx] = np.sum(mask_10 & bin_mask)
            ## this is counted twice due to masking
            return N11/2, N10/2


# %%
## my own directory
os.chdir("/g/data/k10/dl6968/Semi-variogram_AU/")


# %%
# %%time
df_dist = pd.read_csv("./data/pairwise_distances.csv", index_col=0)


# %%
files = sorted(glob("./data/seasonal_p90_500km_1960/*_station_moving_list_all_events.csv"))

# %%
# Define distance bins
bins = np.arange(1,520,25)# np.arange(1, 1525, 25)


# %%
max_pool = 28
for file in files:
    if not os.path.exists(file.replace("station","pair_bins")):
        df_mv = pd.read_csv(file)
        # df_in = df_mv[df_mv["distance"]<=350]
        args_list = [ (df_mv[df_mv["Day"]==day], bins, day) for day in np.unique(df_mv["Day"]) ]
        with Pool(max_pool) as p:
            pool_outputs = list(tqdm(
                    p.imap(calculate_moving_pairs,
                           args_list),
                total=len(args_list),
                position=0, leave=True,  desc=file,
            )
            )
        p.join()
        
        bins_dict = {"Day": [], "Date": [], "Bins": [], "N11": [], "N10": [], "gamma": []}

        for i, output in enumerate(pool_outputs):
            bins_dict["Day"].append( [int(i)] * len(bins))
            bins_dict["Date"].append([df_mv[df_mv["Day"]==i]["Date"].iloc[0]] * len(bins))
            bins_dict["Bins"].append(bins)
            bins_dict["N11"].append(output[0])
            bins_dict["N10"].append(output[1])
            bins_dict["gamma"].append(calc_ratio_gamma(output[0],output[1]))

        out_dict = {}
        for keys in bins_dict.keys():
            out_dict[keys] = np.concatenate(bins_dict[keys])
        df_bins =  pd.DataFrame.from_dict(out_dict)
        df_bins.to_csv(file.replace("station","pair_bins"))
        
        clear_output(wait=True)

# %%
