# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 09:43:51 2016

@author: Gerrit
"""


import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import interpolate





directory = os.getcwd()
powercurve_dir = os.path.join(directory,'Databasepowercurves(February2016) - Processed')
windspeeddata_dir = os.path.join(directory,'windspeed_data')
plant_file = os.path.join(directory,'plant_list.csv')

plant_list = pd.read_csv(plant_file)

turbine_names = np.unique(plant_list['turbine_name'].values)


def wind_calc(wind_speed, turbine_name):    
    power_curve = curve_dict2[turbine_name]
    wind_power = np.array([])
    last_step = None
    for step in wind_speed:
        if step>len(power_curve)-1: # if wind speed is larger than the max value in power curve (30), output is zero
            wind_power = np.hstack((wind_power,0))
            last_step = step
            continue
        if last_step>10 and wind_power[-1]==0:
            wind_power = np.hstack((wind_power,0)) # if wind power last step was zero, and wind speed still > 10, leave it at zero (i.e. wait for wind to die down)
            last_step = step
            continue
        flo, cei = np.floor(step), np.ceil(step) # for linear interpolation of power curve
        last_step = (power_curve[cei]-power_curve[flo])*(step-flo)+power_curve[flo] # linear interpolation 
        wind_power = np.hstack((wind_power,last_step))
    return wind_power



# ------------- Fit a spline interpolation to each power curve -------------- #
curve_dict = {}
curve_dict2 = {}

for turbine_name in turbine_names:
    power_curve = pd.read_csv(os.path.join(powercurve_dir,turbine_name+'.csv'))
    x = power_curve['wind_speed'].values # wind speed (m/s)
    y = power_curve['power'].values # power (MW)
    tck = interpolate.splrep(x,y,k=1,s=0)
    curve_dict[turbine_name] = tck
    curve_dict2[turbine_name] = y

    power = interpolate.splev(x,tck,der=0)
    plt.plot(x,y)
    plt.plot(x,power)

#for i,row in plant_list.iterrows():
#    plant = row['plant_name']
#    turbine_name = row['turbine_name']
#    tck = curve_dict[turbine_name]
#
#    wind_data_df = pd.read_csv(os.path.join(windspeeddata_dir,plant+'.csv'))
#    wind_speed = wind_data_df['wind_speed'].values #TODO: index by position rather than name?
#    # TODO: need to look up hub height and adjust wind speed for hub height!
#    # TODO: also adjust for elevation of location --> higher elevation means lower output
#
#    wind_data_df['power'] = interpolate.splev(wind_speed,tck,der=0) #TODO: check naming is consistent with previous data
#
#    # adjustment: output is zero after cutout is reached until wind speed is back below cutout windspeed - 5m/s (20 m/s)
#    cutout_idxs = wind_data_df[wind_data_df['wind_speed']>25].index.values
#    for cutout_idx in cutout_idxs:
#        i = cutout_idx + 1
#        while wind_data_df.loc[i,'wind_speed'] >= 20 and i < len(wind_data_df):            
#            wind_data_df.loc[i,'power'] = 0
#            i = i + 1

# Ryan's script is WAYYY slower than my spline interpolation script (because we can apply that to the range at once, vectorization)

for i,row in plant_list.iterrows():
    plant = row['plant_name']
    turbine_name = row['turbine_name']

    wind_data_df = pd.read_csv(os.path.join(windspeeddata_dir,plant+'.csv'))
    wind_speed = wind_data_df['wind_speed'].values #TODO: index by position rather than name?
    # TODO: need to look up hub height and adjust wind speed for hub height!
    # TODO: also adjust for elevation of location --> higher elevation means lower output

    wind_data_df['power'] = wind_calc(wind_speed,turbine_name)


    
#    print plant
#    print turbine_name





# -------------- Code from wind simulator.py ------------------------------- #
# (code is in S:\Ryans_Documents\python\wind_sim)

#measure_height = 10
#
##Hub height as a function of DC nameplate (y = a*x^b)
#hh_a = 13.084
#hh_b = 0.2356
#
#stability_coefficient = 0.143
#
#
#summary_data = []
#
#wind_dir = os.listdir(r'S:\Ryans_Documents\python\wind_sim\Databasepowercurves\for_sim')
#
#lat_lon_elevtion = np.genfromtxt(r"S:\Ryans_Documents\python\solar_sim\lat_lon_elevation.csv", delimiter = ',')[:,:2]
#
#turbine_dic = {} 
#
#for turbine in wind_dir:
#    data = []
#    with open(r'S:/Ryans_Documents/python/wind_sim/Databasepowercurves/for_sim/'+turbine, 'rb') as f:
#        reader = csv.reader(f, delimiter = ',')
#        for row in reader:
#            data+=row
#    numbers = []     
#    for inst in data:
#        try:
#            numbers.append(float(inst))
#        except:
#            continue
#    turbine_dic[float(turbine[:-4])] = np.array(numbers[4:])/float(turbine[:-4])
#
#def wind_calc(dc_nameplate, wind_speed):
#    hub_height = hh_a*dc_nameplate**hh_b
#    u = wind_speed*(hub_height/measure_height)**stability_coefficient
#    keys = np.array(turbine_dic.keys())
#    
#    power_curve = turbine_dic[keys[np.argmin(np.abs(keys-dc_nameplate))]]
#    wind_power = np.array([])
#    last_step = None
#    for step in u:
#        if step>len(power_curve)-1: # if wind speed is larger than the max value in power curve (30), output is zero
#            wind_power = np.hstack((wind_power,0))
#            last_step = step
#            continue
#        if last_step>10 and wind_power[-1]==0:
#            wind_power = np.hstack((wind_power,0)) # if wind power last step was zero, and wind speed still > 10, leave it at zero (i.e. wait for wind to die down)
#            last_step = step
#            continue
#        flo, cei = np.floor(step), np.ceil(step) # for linear interpolation of power curve
#        last_step = (power_curve[cei]-power_curve[flo])*(step-flo)+power_curve[flo] # linear interpolation 
#        wind_power = np.hstack((wind_power,last_step))
#    return wind_power


