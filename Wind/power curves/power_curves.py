# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 09:43:51 2016

Script to test the wind calculation function
Uses made up wind_speed data files and plant_list files

@author: Gerrit
"""


import pandas as pd
import numpy as np
import os
import scipy.interpolate


# ------------- Defined functions ------------------------------------------- #

def wind_calc(wind_speed,density,measure_height,hub_height,tck,cutout, n_turbines=1, stability_coefficient = 0.143, hysteris=True):
    '''
    Calculates power output for array of wind_speed values, inputs:
        wind_speed: numpy array with wind speeds (m/s)
        density: numpy array with densitiies (kg/m3)
        measure_height: int, height at which wind speed and density is measured (m)
        hub_height: int, height of wind turbine hub (m)
        tck: spline interpolation of power curve
        cutout: int, cutout windspeed (m/s)
        stability_coefficient: wind shear factor for adjusting wind speed by hub height, default = 0.143
        hysteresis: boolean, whether to apply hysteresis after cutout, default True
    Returns numpy arrays adjusted wind speed and power output
    '''
    
    # Adjust wind speed for hub height, Reference Manual for SAM Wind Power Performance Model, p26-27
    adjusted_wind_speed = wind_speed * (float(hub_height)/measure_height)**stability_coefficient 

    # Adjust wind speed for wake losses, WIND Toolkit - Applied energy - Final, p5/p358
    adjusted_wind_speed = adjusted_wind_speed * (1-1/20.*(n_turbines-1)/7)

    # Apply power curve
    power = scipy.interpolate.splev(adjusted_wind_speed,tck,der=0) * n_turbines
    
    # Calculate density adjustment and adjust for air density, Reference Manual for SAM Wind Power Performance Model, p26-27
    reference_density = 1.225 # air density (kg/m3) at sea level at 15C 
#    density_adjustment = density / reference_density # to adjust power output for air density AFTER applying power curve
#    power = np.multiply(power,density_adjustment)
    
    if hysteris:
        # Adjust for hysteris: output is zero after cutout is reached until wind speed is back below (cutout windspeed - 5 m/s)
        cutout_idxs = np.where(adjusted_wind_speed > 25)[0]
        for cutout_idx in cutout_idxs:
            i = cutout_idx + 1
            while adjusted_wind_speed[i] >= (cutout -5) and i < len(adjusted_wind_speed):
                power[i] = 0
                i = i+1
#        print 'finished hysteris adjustment for %s events' % str(len(cutout_idxs))
    
    return (adjusted_wind_speed, power)

def wind_calc_ryan(wind_speed, density,measure_height,hub_height,power_curve,cutout, n_turbines=1, stability_coefficient = 0.143):
    '''
    Ryan's wind calculation function which doesn't use spline interpolation
    power_curve is just an array that starts with wind power at 0 m/s and gives it for each additional m/s wind
    '''
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


    
# ------------------ INPUTS -------------------------------------------------#

measure_height = 100 # measure height of the wind speed data

directory = os.getcwd()
powercurve_dir = os.path.join(directory,'Databasepowercurves(February2016) - Processed')
windspeeddata_dir = os.path.join(directory,'windspeed_data')
plant_file = os.path.join(directory,'plant_list.csv')
turbine_metadata_file = 'metadata.csv'

plant_list = pd.read_csv(plant_file)
turbine_names = np.unique(plant_list['turbine_name'].values)


# ------------- Fit a spline interpolation to each power curve -------------- #

turbine_metadata = pd.read_csv(os.path.join(powercurve_dir,turbine_metadata_file),index_col='turbine_name')
curve_dict = {}
for turbine_name in turbine_names:
    power_curve = pd.read_csv(os.path.join(powercurve_dir,turbine_name+'.csv'))
    x = power_curve['wind_speed'].values # wind speed (m/s)
    y = power_curve['power'].values # power (MW)
    tck = scipy.interpolate.splrep(x,y,k=1,s=0)
    cut_out_wind_speed = turbine_metadata.loc[turbine_name,'cut_out_wind_speed']
    curve_dict[turbine_name] = (tck,cut_out_wind_speed)


# ------------ Go through plant list and calculate wind power output -------- #
for i,row in plant_list.iterrows():
    plant = row['plant_name']
    turbine_name = row['turbine_name']
    hub_height = row['hub_height']
    tck,cut_out_wind_speed = curve_dict[turbine_name]

    wind_data_df = pd.read_csv(os.path.join(windspeeddata_dir,plant+'.csv'))
    wind_speed = wind_data_df['wind_speed'].values 
    density = wind_data_df['density'].values

    adjusted_wind_speed, power = wind_calc(wind_speed,density,measure_height,hub_height,tck,cut_out_wind_speed)

