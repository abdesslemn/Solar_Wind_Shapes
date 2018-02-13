# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 17:47:46 2017

@author: Gerrit
"""

import numpy as np
import site #import additional module for SAM simulation
site.addsitedir('S:\Gerrits_Documents\Renewable Profile Creation\PV\sdk-release\languages\python') # Use site.addsitedir() to set the path to the SAM SDK API. Set path to the python directory.
import sscapi # can only import this after directory is added to site module on where to find it

class PVSystem:
    """
    Class that holds NREL Software Simulation Core objects and data containers
    Uses SAM's Software Development Kit (SDK)
    Allows you to run PVWatts in python
    """
    def __init__(self,
                 weather_dataframe, # weather data from NSRDB (pd.DataFrame)
                 weather_meta_data, # meta data for the weather data file above (pd.DataFrame)
                 system_capacity_dc = 1, # DC system capacity (kW or MW)
                 dc_ac_ratio = 1.1, # Set DC/AC ratio (or inverter loading ratio). See https://sam.nrel.gov/sites/default/files/content/virtual_conf_july_2013/07-sam-virtual-conference-2013-woodcock.pdf
                 tilt = 0, # tilt of system in degrees (0 = horizontal)
                 azimuth = 180, # azimuth angle (in degrees) from north (180 = south facing)
                 inv_eff = 96, # inverter efficiency in percent
                 losses = 14.0757, # system losses in percent (soiling, shading, wiring, etc.)
                 array_type = 0, # specify fixed tilt system (0=Fixed, 1=Fixed Roof, 2=1 Axis Tracker, 3=Backtracted, 4=2 Axis Tracker)
                 gcr = 0.4, # ground coverage ratio
                 adjust_constant = 0 # constant loss adjustment
                 ):
        
        # Set up Software Simulation Core (SSC) Object
        ssc = sscapi.PySSC()

        # Set Up Data Containers
        dat = ssc.data_create() # container for all input data
        wfd = ssc.data_create() # weather file data container
        
        # Fill wfd container with weather data
        ssc.data_set_number(wfd, 'lat', weather_meta_data['Latitude'])
        ssc.data_set_number(wfd, 'lon', weather_meta_data['Longitude'])
        ssc.data_set_number(wfd, 'tz', weather_meta_data['Local Time Zone'])
        ssc.data_set_number(wfd, 'elev', weather_meta_data['Elevation'])
        ssc.data_set_array(wfd, 'year', weather_dataframe.index.year)
        ssc.data_set_array(wfd, 'month', weather_dataframe.index.month)
        ssc.data_set_array(wfd, 'day', weather_dataframe.index.day)
        ssc.data_set_array(wfd, 'hour', weather_dataframe.index.hour)
        ssc.data_set_array(wfd, 'minute', weather_dataframe.index.minute)
        ssc.data_set_array(wfd, 'dn', weather_dataframe['DNI'])
        ssc.data_set_array(wfd, 'df', weather_dataframe['DHI'])
        ssc.data_set_array(wfd, 'wspd', weather_dataframe['Wind Speed'])
        ssc.data_set_array(wfd, 'tdry', weather_dataframe['Temperature'])

        # Add wfd container to dat container
        ssc.data_set_table(dat, 'solar_resource_data', wfd)
        ssc.data_free(wfd)

        # Specify the system Configuration
        ssc.data_set_number(dat, 'system_capacity', system_capacity_dc)
        ssc.data_set_number(dat, 'dc_ac_ratio', dc_ac_ratio)
        ssc.data_set_number(dat, 'tilt', tilt)
        ssc.data_set_number(dat, 'azimuth', azimuth)
        ssc.data_set_number(dat, 'inv_eff', inv_eff)
        ssc.data_set_number(dat, 'losses', losses)
        ssc.data_set_number(dat, 'array_type', array_type)
        ssc.data_set_number(dat, 'gcr', gcr)
        ssc.data_set_number(dat, 'adjust:constant', adjust_constant)
        
        # Add the software simulation core object and the input data object as attributes to the PVSystem object
        self.ssc = ssc
        self.dat = dat

    def simulate(self):
        """
        Simulate PV generation for defined system
        Returns:
            output (np.array): Array of hourly generation in MW (ac)
        """
        
        # Create PVWatts module, execute, and save results in dataframe
        mod = self.ssc.module_create('pvwattsv5') # create a pvwatts module (pvwattsv5 is the current version of the module)
        self.ssc.module_exec(mod, self.dat)
        
        return np.array(self.ssc.data_get_array(self.dat, 'gen'))