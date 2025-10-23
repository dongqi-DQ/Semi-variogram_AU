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
from glob import glob
from multiprocess import Pool
from tqdm import tqdm
import os

# %%
inpath = "/g/data/w40/dl6968/BoM_daily_stations/percentiles/"
outpath = "/g/data/k10/dl6968/BoM_daily_station/prcp_pc_ts/"

files = sorted(glob(f"{inpath}*"))


# %%
def timestamp_daily(file):
    if not os.path.exists(file.replace(inpath,outpath)):
        with xr.open_dataset(file) as ds:
            ds_daily = ds.resample(time="D").mean()
            ds_daily.to_netcdf(file.replace(inpath,outpath))


# %%
max_pool = 28
with Pool(max_pool) as p:
    pool_outputs = list(
        tqdm(
            p.imap(timestamp_daily,
                   files[1500:]),
            total=len(files[1500:]),
            position=0, leave=True
        )
    )
    p.join()

# %%
