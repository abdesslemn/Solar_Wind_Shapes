# -*- coding: utf-8 -*-
"""
Created on Mon Oct 06 14:06:53 2014
@author: zachary.ming & gerrit.demoor
Edited by Huai
"""

import os
import calendar
import numpy as np
import pandas as pd

import pvsim_sam
import load_data

    #functionality to create profiles for multiple location_data lists at once
    #functionality to aggregate profiles in this script at once, save them while creating profiles?

    #TODO: add a wait to the NSRD API ping
    #TODO: fix functionality for 30 min intervals
    #TODO: have a log of the inputs that were used? (so you know what timezone was used)
    #TODO: Limit MW per nsrdb so you get more diverse profiles? hard to force  to another site though, not sure how far lat,lon you have to go 
    
    # grid cells are 0.038 by 0.038 so there is a third number precision? (seems like nsrdb only uses 2 numbers of precision though)
    # i.e. downloading for lat 12.1292 is same as downloading for lat 12.13
    

# --------------------------- INPUTS ------------------------------------- ##

# key inputs
start_year = 1998 # Int, must be >= 1998
end_year = 1999 # Int, must be <= 2014
interval = 60 # Int, time interval in minutes, i.e., '30' is half hour intervals. Valid intervals are 30 & 60.
timezone = -8 # Int, timezone relative to UTC to output profiles. Note that NSRDB data will be reported in local timezone
leap_year = 'true' # Str, true will return leap day data if present, false will not. CASE SENSTIIVE!!!
include_date = True # include column with data index (triples file size)

# intputs for aggregation only
output_project_shapes = True # output individual project shapes
aggregate_by_zone = True # output zonal aggregation
aggregate_by_state = False # output statewide aggregation
normalize_agg = True # normalize profiles to installed capacity

# directory settings
rad_data_dir = 'NSRDB data' # this is the folder where the weather NSRDB data is stored
location_data_dir = 'location_data' # inputs folder
profiles_dir = 'generation_profiles' # outputs folder

# other
output_header = 'output (mw)'
years = range(start_year,end_year+1)
n_years = (end_year - start_year + 1)

 # ------------------------- HELPER FUNCTIONS ----------------------------- ##

def correct_timezone(array,timezone_adjustment):
    """
    shift profile by nperiods = timezone, to account for timezone
    Args:
        timezone_adjustment (int): number of periods to shift
        e.g. to go from UTC to PST, timezone_adjustment = -8
             to go from EST to PST, timezone_adjustment = -3
    """
    if timezone < 0:
        return np.hstack((array[-timezone_adjustment:],array[:-timezone_adjustment]))
    elif timezone > 0:
        return np.hstack((array[-timezone_adjustment:],array[:-timezone_adjustment]))
    else:
        return array

# ------------------------ SIMULATE PROFILES ----------------------------- ##

directory = os.getcwd()
location_data_path = os.path.join(directory,location_data_dir)
cases = os.listdir(location_data_path)

for case in cases:
    print '\n------------------------------------------------------'
    print 'processing case list ... ' + case
    print '------------------------------------------------------'

    # Create the necessary folders
    write_dir = case.split('.')[0]
    if not os.path.exists(os.path.join(directory,profiles_dir,write_dir)):
        os.makedirs(os.path.join(directory, profiles_dir,write_dir))

    # read in plant profile list
    projectlist = pd.read_csv(os.path.join(location_data_path,case))
    zones = np.unique(projectlist['zone'].values)
    states = np.unique(projectlist['state'].values)

    # initialize lists that track outputs
    zone_profiles = {} # cumulative output by zone
    state_profiles = {} # cumulative output by state
    shape_capacities = []
    shape_energies = []
    target_cfs = []
    cfs = [] # cf by project

    # loop over every project

    for i,row in projectlist.iterrows():

        # Extract the plant parameters
        project = row['project name']
        zone = row['zone']
        state = row['state']
        latitude = row['latitude'] # NSRDB lat and lon have 4 numbers of precision
        longitude = row['longitude']
        acsize = row['ac capacity (mw)']
        energy = row['energy (mwh)']
        dc_ac_ratio = row['inverter loading ratio']
        tilt = row['tilt']
        array_type = row['tracking type'] # (0=Fixed, 1=Fixed Roof, 2=1 Axis Tracker, 3=Backtracted, 4=2 Axis Tracker)

        print '\nprocessing project ... ' + project

        profiles_by_year = []
        # simulate each year individually
        for year in years:

            # load weather data from nsrdb file (download if necessary) including calculated sunangles (using sunangle library)
            weather_dataframe, weather_meta_data = load_data.get_nsrdb_data(year,latitude,longitude,rad_data_dir,interval,leap_year)

            # create an instance of a pv system with certain specifications from pvsim_sam library
            pvsystem = pvsim_sam.PVSystem(weather_dataframe, weather_meta_data, 
                                          system_capacity_dc = acsize * dc_ac_ratio, dc_ac_ratio = dc_ac_ratio, tilt = tilt, array_type = array_type)

            # simulate PV output on this particular pvsystem and weather data
            print 'simulating solar generation for ... ' + project + ' - ' + str(year)
            output = pvsystem.simulate()

            # shift profile by nperiods to match required timezone
            # note that SAM requires NSRDB data to be in lcoal timezone, and output profiles will also be in local timezone
            profile_timezone = weather_meta_data['Local Time Zone'] # note: this will return a series! (index [0] if you want number)
            timezone_adjustment = int(timezone - profile_timezone) # e.g. if profile is EST (-5) and prefered timezone is PTC (-8), shift is -8 + 5 = -3 
            output = correct_timezone(output,timezone_adjustment)
            
            # add leap year outputs using previous day generation (SAM skips the leap year days!)
            if calendar.isleap(year):
                profile_start = output[:((31+28) * 24 * 60/interval)]
                profile_add_day = output[((31+27) * 24 * 60/interval):((31+28) * 24 * 60/interval)]
                profile_end = output[((31+28) * 24 * 60/interval):]
                output = np.concatenate([profile_start, profile_add_day, profile_end])

            profiles_by_year.append(output)

        # Concatenate individual year profiles into one long profile
        profile = pd.DataFrame(data=np.hstack(tuple(profiles_by_year)),columns=[output_header])
        
        # Add date index if toggle is on
        if include_date:
            rng = pd.date_range(start=pd.datetime(start_year,1,1,0),end=pd.datetime(end_year,12,31,23),freq='H') #TODO: freq should be dependent on interval
            if leap_year == 'false': 
                rng = rng[~((rng.month == 2) & (rng.day == 29))] # remove leap days
            profile.index = rng
            profile.index.name = 'date (tz = ' + str(timezone) + ')'

        # output individual full (all years together) shapes if toggle is on
        if output_project_shapes:
            if start_year != end_year:
                filename = '_'.join([project,str(start_year),str(end_year)]) + '.csv'
            else:
                filename = '_'.join([project,str(start_year)]) + '.csv'
            filepath = os.path.join(profiles_dir,write_dir,filename)
            profile.to_csv(filepath,index=include_date)

        # calculate metrics and store for output
        shape_capacities.append(acsize)
        shape_energies.append(profile.sum()[output_header] / float(n_years) )
        cfs.append(profile.mean()[output_header] / acsize)
        target_cfs.append(energy / (acsize * 8760.0) ) # TODO: adjust for interval

        # add the project output to the running total of zonal outputs
        if zone in zone_profiles:
            zone_profiles[zone] += profile[output_header]
        else:
            zone_profiles[zone] = profile[output_header].copy() # copy, otherwise get's changed in place (double counting!)
            
        # add the project output to the running total of state outputs    
        if state in state_profiles:
            state_profiles[state] += profile[output_header]
        else:
            state_profiles[state] = profile[output_header].copy()


    # Construct and save a summary of sites file (projectlist + more metrics)
    projectlist['shape capacity (MW)'] = shape_capacities # redundant for now, but could be useful when we have option to scale for energy
    projectlist['shape energy (MWh)'] = shape_energies
    projectlist['target CF'] = target_cfs
    projectlist['shape CF'] = cfs
    filename = '_'.join(['summary_of_sites', case.split('.')[0], str(start_year), str(end_year)]) +'.csv'
    projectlist.to_csv(os.path.join(profiles_dir,write_dir,filename),index=False)

    # aggregate by zone if toggle is on (and normalize if selected)
    if aggregate_by_zone:
        for zone,profile in zone_profiles.items():
            print '\nAggregating output for zone ... %s' % zone

            capacity = projectlist.groupby(['zone'])['ac capacity (mw)'].sum()[zone]
            if normalize_agg:
                profile = profile / capacity

            filename = '_'.join(['aggregated',zone,str(start_year),str(end_year)]) +'.csv'
            filepath = os.path.join(profiles_dir,write_dir,filename)
            profile.to_csv(filepath,index=include_date,header=True)

    # aggregate by state if toggle is on (and normalize if selected)
    if aggregate_by_state:
        for state,profile in state_profiles.items():
            print '\nAggregating output for state ... %s' % state

            capacity = projectlist.groupby(['state'])['ac capacity (mw)'].sum()[state]
            if normalize_agg:
                profile = profile / capacity
            filename = '_'.join(['aggregated',state,str(start_year),str(end_year)]) +'.csv'
            filepath = os.path.join(profiles_dir,write_dir,filename)
            profile.to_csv(filepath,index=include_date,header=True)

profile.head(24)

