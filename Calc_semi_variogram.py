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
from multiprocess import Pool
from IPython.display import clear_output


# In[2]:


## my own directory
os.chdir("/g/data/k10/dl6968/Semi-variogram_AU/")


# In[25]:


### functions 
def find_neighbour_stations(df, df_lat, df_lon, df_id, center_lat, center_lon, search_radius=350):
    '''
    find stations within a given radius
    '''
    search_stations = []

    for i in range(0,len(df)):
        station  = (df_lat.iloc[i], df_lon.iloc[i])
        distance = geodesic((center_lat,center_lon), station).kilometers
    
        if distance <= search_radius:
                search_stations.append(df_id.iloc[i])

    local_df = df[df_id.isin(search_stations)]    
    
    return local_df

def haversine(lat1, lon1, lat2, lon2):
    '''
    calculate distance in km from lats and lons
    '''
    # Radius of the Earth in kilometers
    R = 6371.0
    
    # Convert latitude and longitude from degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)
    
    # Difference in latitudes and longitudes
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    # Haversine formula
    a = math.sin(dlat / 2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    
    # Distance in kilometers
    distance = R * c
    
    return distance
    
def distance_stn(station_id):
    '''
    calculate distance between two stations
    '''
    lat2 = df_search[df_search["ID"] == station_id]["Latitude"].values[0]
    lon2 = df_search[df_search["ID"] == station_id]["Longitude"].values[0]
                    
    local_distance = haversine(cent_lat, cent_lon, lat2, lon2)
    return local_distance

def bearing_stn(station_id):
    '''
    calculate bearing between a neighbour station and the center station
    '''
    lat2 = df_search[df_search["ID"] == station_id]["Latitude"].values[0]
    lon2 = df_search[df_search["ID"] == station_id]["Longitude"].values[0]
                    
    local_bearing = calculate_bearing(cent_lat, cent_lon, lat2, lon2)
    return local_bearing
    
def calculate_bearing(lat1, lon1, lat2, lon2):
    '''
    calculate bearing with given lats and lons
    '''
    # Convert latitude and longitude from degrees to radians
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)

    # Calculate the difference in longitude
    delta_lon = lon2 - lon1

    # Calculate the bearing using the formula
    x = math.sin(delta_lon) * math.cos(lat2)
    y = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(delta_lon)
    
    initial_bearing = math.atan2(x, y)

    # Convert from radians to degrees
    initial_bearing = math.degrees(initial_bearing)

    # Normalize the bearing to 0-360 degrees
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing

def calc_p90(ds,var_name=None):
    precip_da = ds

    if var_name!=None:
        precip_da = ds[var_name].load()
    filtered_da = precip_da.where(precip_da > 1, drop=True)
    # Compute the 99th percentile along the time dimension
    p90 = filtered_da.quantile(percentile, dim='time')
    
    # Find the datetimes where the rain rate is higher than the 90th percentile
    rain_rate_above_p90 = precip_da >= p90
    
    # Get the datetimes where rain rate is higher than the 90th percentile
    datetimes_above_p90 = ds['time'][rain_rate_above_p90]

    return p90, datetimes_above_p90

def find_val_extreme(station_id):
    '''
    identify whether an extreme day for the station at the semi-variogram center
    is also an extreme day for a neighbour station
    '''
    search_id = str(station_id).zfill(6)
    
    ds_search = xr.open_dataset(f"/g/data/k10/dl6968/BoM_daily_station/prcp_pc_ts/{search_id}.nc")
    precip_da = ds_search["prcp"]
    stn_p90 = df_pc[df_pc["ID"]==station_id][pc_str].values[0]
    if stn_p90>0:
    ## I've already filled in all the no data days with -1
        extreme_val = precip_da.sel(time=extreme_dates)
        extreme_pc = ds_search["percentile"].sel(time=extreme_dates)
        ## flag NaN days with -1, below p90 days 0, and take the values for above p90 days
        extreme_flag =  extreme_val.copy()
        extreme_flag = extreme_flag.where(extreme_val>=stn_p90, 0)
        extreme_flag = extreme_flag.where(extreme_val!=-1, -1)
        return extreme_val.values.flatten(), extreme_flag.values.flatten(), extreme_pc.values.flatten()
    else:
        extreme_val = np.zeros((len(extreme_dates), 1))  # Local array for this station
        extreme_flag = np.zeros((len(extreme_dates), 1))  # This is a flag with NaN as -1, <p90 as 0, and >=p90 as value
        # ds_pc = xr.open_dataset(f"/g/data/w40/dl6968/BoM_daily_stations/percentiles/{search_id}.nc")
        extreme_pc = np.zeros((len(extreme_dates), 1)) 
        extreme_pc[:] = -1
        extreme_val[:] = -1
        extreme_flag[:] = -1
        return extreme_val.flatten(), extreme_flag.flatten(), extreme_pc.flatten()
        


def process_avail_distance(dates):
    '''
    get distance between center and neighbor stations for plotting later
    also set filter for valid stations
    '''
    val_arr = df_flag[dates].values
    local_distance = df_search["Distance"].values
    
    ## need to have more than 20 stations in the neighbourhood
    count = len(np.argwhere(val_arr>=0))
    if count >= 20:
                
        return dates, local_distance
    else:
        return dates, None 



    
def process_date_stations(m):
    dates = extreme_dates[m]
    
    # Preallocate and avoid repetitive dictionary creation
    local_station_dict = {
        "Day": [], "Date": [], "Spec_stn": [], "cent_lat": [], "cent_lon": [],
        "Neighb_stn": [], "lat": [], "lon": [], "distance": [], "angle": [],
        "val": [], "flag": [], "pc": []
    }


    # Filter valid stations first
    valid_mask = df_val[dates] > -1
    valid_stations = df_val[valid_mask]

    # Retrieve data for all valid stations
    station_ids = valid_stations["ID"].values
    lats = df_search.set_index("ID").loc[station_ids, "Latitude"].values
    lons = df_search.set_index("ID").loc[station_ids, "Longitude"].values
    distances = distance_arr[m, valid_mask]
    local_bearings = df_search["Bearing"].values[valid_mask]
    values = valid_stations[dates].values
    flags = df_flag[dates][valid_mask].values
    pcs = df_pc_val[dates][valid_mask].values

    # Append data to the station dictionary
    local_station_dict["Day"] = [m] * len(valid_stations)
    local_station_dict["Date"] = [dates] * len(valid_stations)
    local_station_dict["Spec_stn"] = [spec_id] * len(valid_stations)
    local_station_dict["cent_lat"] = [cent_lat] * len(valid_stations)
    local_station_dict["cent_lon"] = [cent_lon] * len(valid_stations)
    local_station_dict["Neighb_stn"] = station_ids.tolist()
    local_station_dict["lat"] = lats.tolist()
    local_station_dict["lon"] = lons.tolist()
    local_station_dict["distance"] = distances.tolist()
    local_station_dict["angle"] = local_bearings.tolist()
    local_station_dict["val"] = values.tolist()
    local_station_dict["flag"] = flags.tolist()
    local_station_dict["pc"] = pcs.tolist()

    return local_station_dict


# In[ ]:





# In[13]:


##### main script starts here
## find available stations 

df = pd.read_csv("./data/BoM_daily_stations.csv")

## remove stations that do not have data 
exclude_stn = []
for stn_id in df["ID"]:
    bom_id = str(stn_id).zfill(6)
    # if not os.path.exists(f'/g/data/w40/dl6968/BoM_daily_stations/netcdf/{bom_id}.nc'):
    #     exclude_stn.append(stn_id)
    if not os.path.exists(f'/g/data/k10/dl6968/BoM_daily_station/prcp_pc_ts/{bom_id}.nc'):
        if stn_id not in exclude_stn:
            exclude_stn.append(stn_id)

## mannually remove some faulty stations
df = df[~df["ID"].isin(exclude_stn)& (df["End_Year"]>=1960) & (df['ID'] != 40592) & \
     (df['ID'] != 40593) & (df['ID'] != 58090)& (df['ID'] != 68002) & (df['ID'] != 64003) & (df['ID'] != 29051)\
   & (df['ID'] != 34050) & (df['ID']!=40646)& (df['ID']!=95009) & (df['ID']!=70041) & (df['ID']!=88089) \
  & (df['ID']!=68046) & (df['ID']!=40509) & (df['ID']!=68057)& (df['ID']!=63088) & (df['ID']!=68066) \
& (df['ID']!=43088)& (df['ID']!=86175)] 

# 
daily_lat = []
daily_lon = []
for i in range(0, len(df)):
    daily_lat.append(df["Latitude"].iloc[i])
    daily_lon.append(df["Longitude"].iloc[i])


# In[18]:


## for stations to be a neighbour station
df_neighbour = df[df["Years"] >=20]
## for stations to be a station to do semi-variogram
df_spec = df[df["Years"] >=20]
id_list = list(df_spec["ID"])
print(f"{len(id_list)} stations in the list")


# In[6]:


percentile = 0.90 ## this is for extreme days
pc_str = "P90" ## this is for the neighbours in the dataframe
df_pc = pd.read_csv("./data/BoM_stn_p90.csv")
max_pool = 28 ## number of CPUs for processing
max_radius = 350


# In[ ]:


for spec_id in id_list:
    csv_file = f"./data/all_AU_p90/{spec_id}_station_list_all_events.csv"
    ## do semi-variogram if file not exists
    if not os.path.exists(csv_file):
        print(f"Processing station {spec_id}")
        cent_lat, cent_lon = df[df["ID"]==spec_id]["Latitude"].values[0],  df[df["ID"]==spec_id]["Longitude"].values[0]
        df_search1 = find_neighbour_stations(df_neighbour, df_neighbour["Latitude"], df_neighbour["Longitude"], df_neighbour['ID'], cent_lat, cent_lon, max_radius)
        df_search = df_search1[(df_search1["ID"] != spec_id) & (df_search1["End_Year"]>=1960)].copy()
        distance = []
        bearings = []
        for stn_id in df_search["ID"]:
            distance.append(distance_stn(stn_id))
            bearings.append(bearing_stn(stn_id))
        df_search["Distance"] = distance
        df_search["Bearing"] = bearings
        
        ## if not enough neighbours
        ## skip the station
        if len(df_search)<20:
            continue
        else:
            ds_spec = xr.open_dataset("/g/data/k10/dl6968/BoM_daily_station/prcp_pc_ts/"+str(spec_id).zfill(6)+".nc")
            var = "prcp"
            ds_sel = ds_spec.sel(time=slice("1940-03-02", "2024-06-30"))
            prcp_p99, dt_p99 = calc_p90(ds_sel, var)
            ds_sel.close()
            ds_spec.close()
            extreme_dates = []
            for dates in dt_p99:
                yymmdd = str(dates.values)[:10]
                if yymmdd not in extreme_dates:
                    extreme_dates.append(yymmdd)
                    
            print("Look through neighbour stations")
            ## find values and flag for that extreme day
            ## search all stations within the 350 km (or any given radius)
            with Pool(max_pool) as p:
                pool_outputs = list(
                    tqdm(
                        p.imap(find_val_extreme,
                               df_search["ID"].values),
                        total=len(df_search["ID"]),
                        position=0, leave=True
                    )
                )
            p.join()
            
            bom_arr = np.zeros((len(extreme_dates), len(df_search) ))
            bom_flag = np.zeros((len(extreme_dates), len(df_search) ))
            bom_pc = np.zeros((len(extreme_dates), len(df_search) ))
            
            for istn in range(0, len(df_search)):
                bom_arr[:,istn] = pool_outputs[istn][0]
                bom_flag[:,istn] = pool_outputs[istn][1]
                bom_pc[:,istn] = pool_outputs[istn][2]
        
            bom_tmp = bom_arr.copy()
            bom_arr[np.isnan(bom_tmp)] = -1
            
            ## Assign the daily values for each extreme day and each neighbour station
            ## No values should be -1
            df_val = pd.DataFrame(bom_arr.T, index=df_search["ID"], columns=extreme_dates)
            
            # Reset the index and add the 'ID' column
            df_val.reset_index(inplace=True)
            df_val.rename(columns={'index': 'ID'}, inplace=True)
            del bom_tmp
            
            ## Assign the daily flags for each extreme day and each neighbour station
            ## No values should be -1
            ## Not an extreme is 0
            ## Keep the precipitation value if is an extreme
            bom_tmp = bom_flag.copy()
            bom_flag[np.isnan(bom_tmp)] = -1
            df_flag = pd.DataFrame(bom_flag.T, index=df_search["ID"], columns=extreme_dates)
                        
            # Reset the index and add the 'ID' column
            df_flag.reset_index(inplace=True)
            df_flag.rename(columns={'index': 'ID'}, inplace=True)
            
            del bom_tmp
            
            ## Assign the daily precentiles for each extreme day and each neighbour station
            ## only show the P90 value for the station if a station is available on the day
            bom_pc[np.isnan(bom_pc)] = -1
            
            df_pc_val = pd.DataFrame(bom_pc.T, index=df_search["ID"], columns=extreme_dates)
            
            # Reset the index and add the 'ID' column
            df_pc_val.reset_index(inplace=True)
            df_pc_val.rename(columns={'index': 'ID'}, inplace=True)
        ### station distance for plotting later
        
        distance_arr = np.zeros((len(extreme_dates), len(df_search)))
        
        print("Calc station distance")
        
        with Pool(max_pool) as p:
            pool_outputs = list(
                tqdm(
                    p.imap(process_avail_distance,
                           extreme_dates),
                    total=len(extreme_dates),
                    position=0, leave=True
                )
            )
        p.join()
        

        for i in range(0, len(extreme_dates)):
            date= pool_outputs[i][0]
            ## make sure the dates align
            if date == extreme_dates[i]:
                if pool_outputs[i][1] is None:
                    ## if results are None the data will be set to NaNs
                    distance_arr[i,:] = np.nan
                else:
                   
                    distance_arr[i,:] = pool_outputs[i][1]
        
        
        ### save to CSV
        # Main parallel loop
        station_dict = {"Day": [], "Date": [], "Spec_stn":[], "cent_lat": [], "cent_lon":[],"Neighb_stn": [],
                        "lat": [], "lon": [], "distance": [], "angle": [], "val": [], "flag": [],"pc":[]}

        
        print("Save to CSV")
        
        with Pool(max_pool) as p:
            pool_outputs = list(
                tqdm(
                    p.imap(process_date_stations,
                           range(0,len(extreme_dates))),
                    total=len(extreme_dates),
                    position=0, leave=True
                )
            )
        p.join()
        
        for results in pool_outputs:
            
            for key in station_dict:
                        station_dict[key].extend(results[key])
        
        df_stations =  pd.DataFrame.from_dict(station_dict)
        
        df_stations.to_csv(csv_file, index=False)
        
        clear_output(wait=True)


# In[ ]:




