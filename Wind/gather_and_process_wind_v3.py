# -*- coding: utf-8 -*-
"""
Created on Mon Jun 03 14:26:05 2013
Edited Feb 2016

Script that gathers wind profiles based on site information and creates the RT, DA and HA profiles
HA and DA are provided by NREL in the wind toolkit
However, HA is weird so we use the h-1 output instead. # NOTE_TO_SELF: Is this still the case? Why?

All profiles are scaled by ENERGY, not capacity.  # NOTE_TO_SELF: What does this mean?
The script adds toolkit ids until the target energy is met. # NOTE_TO_SELF: Why do we not just scale?
When looking for toolkit IDs, the capacity in that id grid cell is used to determine the amount of turbines
e.g. grid cell with 16 MW wind with turbines of 0.5 MW has 32 turbines
(in WTK it would be 8 turbines max because WTK assumes 2 MW turbines)
Energy for those turbines is usually lower than WTK energy because lower hub height and worse turbine technology


CURRENT OBSERVATIONS:
energy budget is not completely met for regional aggregations even though we scale for energy
---> is because of shortages, checked for CIPB #NOTE_TO_SELF: What does this mean?
the peak is now way higher because we scale for energy so we keep adding sites to meet energy budget.
Question is whether this peak is realistic?
---> peak is now way more in line with installed capacity, checked for CIPB
(peak is 1580, WTK peak is 1050, installed cap is 1600)


@author: ryan + edits Gerrit + Manu
"""

import csv
import os
import numpy as np
import pandas as pd
import scipy.interpolate
from collections import defaultdict
from math import radians, cos, sin, asin, sqrt
import load_data

# TODO: use sscapi rather than our own power curves
# TODO: clean up script further (e.g. summary of sites should match location information (format and location in folder)
#       remove uncommented code, make file structure same as PV, toggle for normalization by capacity etc.)
#       have option for date column. Add timezone designation in it for clarity
#       leap days
# TODO: add aggregation by state and by zone same way as PV + toggles to turn on/off
# TODO: include script for adjusting shape (Elaine's method)
# TODO: include other wrappers to allow you to match cap factors by picking sites / tech
# (discussion w Doug, perhaps need geo library)
#       e.g. "get me wind shape of NY state as a whole that matches a 40% cap factor"

################################################################################
###############                     User Inputs                #################
################################################################################

# WARNING: make sure timezone assumptions are consistent! (see line 55, 319)

# User inputs
start_year = 2007
end_year = 2007
resolution = 60     # resolution of the WTK data in min
timezone = 0        # timezone vs. UTC (e.g. -8 = PST, -5 = EST). ASSUMES THAT WTK DATA IS IN UTC.
# DOUBLE CHECK THAT IT IS! #NOTE_TO_SELF: Is this checked yet?
MAX_DISTANCE = 50   # max distance in miles for grabbing toolkit sites
MIN_CAP_FACTOR = .20

# toggles 
include_date = True # include column with data index (triples file size)
output_project_shapes = True # output individual project shapes (vs. only aggregate shapes)
adjust_for_energy_shortage = False # True: adjust final aggregated profile for energy shortage
# (added 31/05/16 - slightly buggy still)
scale_to_energy = False # if true, wind farm shape is scaled to meet energy target.
# If False, shape will meet capacity target
exclude_certain_ids = False # to exclude certain WTK ids
# (e.g. for future wind data you want to look at the remaining IDs after running CA/common case)

# for debugging power curves
force_NREL_powercurves = False # True: use same powercurve as WTK,
# False: use user specified power curves from location data table
use_NREL_power = False # True: use the power output straight from NREL, so avoids using the power_calc functions.
# Should give similar answer to force_NREL_power_curves
output_WTK_RT = False # True: write out csvs of the direct wind toolkit power output
# (to compare with our calculated output that is based on WTK wind)

# parameters
MILES_TO_DEGREE = 69.           # degree latitude

# directory structure inputs
directory = os.getcwd()
profile_dir = os.path.join(directory,'generation_profiles')
location_data_dir = os.path.join(directory,'location_data')
wtk_data_dir = os.path.join(directory,'WTK data')
powercurve_dir = os.path.join(directory,'power curves','Databasepowercurves(February2016) - Processed')
used_ids_dir = os.path.join(profile_dir, "CA_existing_wind") # only used when you want to exclude certain ids
cases = os.listdir(location_data_dir)   # listdir: Returns  of all entries in the directory (except starting with .)

site_meta = pd.read_csv(os.path.join(wtk_data_dir, 'wtk_site_metadata.csv'))  # Reads all the meta data
#  NOTE_TO_SELF: Cleaner way please!
turbine_metadata = pd.read_csv(os.path.join(powercurve_dir, 'metadata.csv'), \
                               index_col='turbine_name')  # Names of all turbines in the dbase

################################################################################
###############               Defined Functions                #################
################################################################################

# ---------------------- Distance functions ----------------------------------#


def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))

    # 3957.5 mi is the radius of the Earth
    mi = 3957.5 * c
    return mi


def get_distance(lat, lon, meta_lon_lat):
    """
    Gives the array of distances from current point to WTK database points
    :param lat:
    :param lon:
    :param meta_lon_lat:
    :return: array of distances in miles
    """
    return np.array([haversine(lat, lon, mlat, mlon) for mlon, mlat in meta_lon_lat])


# ------------------- Conversion/Check functions  ----------------------------#


def correct_wind_time(array, nperiods):
    """
    shift profile by nperiods to account for timezone difference
    Args:
        nperiods (int): number of periods to shift
        e.g. to go from UTC to PST, timezone_adjustment = -8
             to go from EST to PST, timezone_adjustment = -3
    """
    if nperiods < 0:
        return np.hstack((array[-nperiods:],array[:-nperiods]))
    elif nperiods > 0:
        return np.hstack((array[-nperiods:],array[:-nperiods]))
    else:
        return array

# ---------------------- Power Curve Functions  ------------------------------#


def wind_calc(wind_speed, density, measure_hub_height_wtk, hub_height, tck, cutout, n_turbines=1, \
              stability_coefficient=0.143, hysteresis=True, dens_adj = True):
    """
    Calculates power output for array of wind_speed values, inputs:
        wind_speed: numpy array with wind speeds (m/s)
        density: numpy array with densities (kg/m3)
        measure_hub_height_wtk: int, height at which wind speed and density is measured (m)
        hub_height: int, height of wind turbine hub (m)
        tck: spline interpolation of power curve
        cutout: int, cutout windspeed (m/s)
        stability_coefficient: wind shear factor for adjusting wind speed by hub height, default = 0.143
        hysteresis: boolean, whether to apply hysteresis after cutout, default True
        dens_adj: boolean, whether to apply air density adjustment, default True
    Returns numpy arrays adjusted wind speed and power output
    """

    # Adjust wind speed for hub height, see Reference Manual for SAM Wind Power Performance Model, p26-27
    adjusted_wind_speed = wind_speed * (float(hub_height)/measure_hub_height_wtk)**stability_coefficient

    # Adjust wind speed for wake losses, see WIND Toolkit - Applied energy - Final, p5/p358
    n_turbines_wake = min(n_turbines,8) # limit turbines to 8 for wind adjustment calculation
    # TODO: allow more wake losses than 5% for more turbines? (e.g. 100s of 225 kW in Altamont)
    adjusted_wind_speed = adjusted_wind_speed * (1-1/20.*(n_turbines_wake-1)/7)

    # Adjust wind speed for air density, see WindPro / Energy Power Curve Air Density Correction, p4
    if dens_adj:
        reference_density = 1.225 # air density (kg/m3) at sea level at 15C
        x = np.array([0, 7.5, 12.5, 25])  # wind speed data points for air density exponent adjustment
        y = np.array([1./3, 1./3, 2./3, 2./3])  # exponent adjustment data points
        tck_air_density = scipy.interpolate.splrep(x, y, k=1, s=0)
        # B-spline representation of interpolation of degree k and smoothing s (larger is more)
        exponents = scipy.interpolate.splev(adjusted_wind_speed,tck_air_density,der=0)  # tck: 3 tuple
        # B spline
        adjusted_wind_speed = adjusted_wind_speed * (density / reference_density) ** (exponents)

    # Apply power curve
    power = scipy.interpolate.splev(adjusted_wind_speed,tck,der=0) * n_turbines

    # Force to zero when above cutout (interpolation calculates non zero wind speed for wind speeds just above cut out)
    power[adjusted_wind_speed>cutout] = 0

    if hysteresis:
        # Adjust for hysteresis: output is zero after cutout is reached until wind speed
        # is back below (cutout windspeed - 5) m/s
        cutout_idxs = np.where(adjusted_wind_speed > cutout)[0]
        for cutout_idx in cutout_idxs:
            i = cutout_idx + 1
            if i == len(adjusted_wind_speed):
                continue # cutout_idx was last element of wind_speed array
            else:
                while adjusted_wind_speed[i] >= (cutout -5):
                    power[i] = 0
                    i = i+1
                    if i == len(adjusted_wind_speed):
                        break # reached end of array
        # print 'finished hysteresis adjustment for %s events' % str(len(cutout_idxs))
    return (adjusted_wind_speed, power)


################################################################################
###########                  Main part of Script                ################
################################################################################

# ------ Set up list of years and dates --------------------------------------#
if resolution == 60:
    dates = pd.date_range(start=pd.datetime(start_year,1,1,0),end=pd.datetime(end_year,12,31,23),freq='H')
elif resolution == 5:
    dates = pd.date_range(start=pd.datetime(start_year,1,1,0),end=pd.datetime(end_year,12,31,23,55),freq='5Min')
else:
    sys.exit("Please check resolution")     # Throw exception and exit if resolution not 5 or 60

# ---- Create Dictionary of turbines mapped to a tuple of turbine power curves, cut out and max power ---- #
# NOTE_TO_SELF: Tuples are ordered whereas sets are not

turbine_names = np.unique(turbine_metadata.index.values)    # unique: Returns nd-array
curve_dict = {}
for turbine_name in turbine_names:
    power_curve = pd.read_csv(os.path.join(powercurve_dir,turbine_name+'.csv'))
    x = power_curve['wind_speed'].values # wind speed (m/s)
    y = power_curve['power'].values # power (MW)
    tck = scipy.interpolate.splrep(x,y,k=1,s=0) # spline inerpolation
    cut_out_wind_speed = turbine_metadata.loc[turbine_name,'cut_out_wind_speed']
    max_power = turbine_metadata.loc[turbine_name,'max_power']
    curve_dict[turbine_name] = (tck, cut_out_wind_speed, max_power)

# ------ Loop through cases of plant lists to create profiles ---------------- #

for case in cases:

    print '\n------------------------------------------------------'
    print 'processing case list ... ' + case
    print '------------------------------------------------------'

    # Create the necessary folders
    case_name = case.split('.')[0]   # Name before .csv (Actually the first name before .)
    if not os.path.exists(os.path.join(profile_dir,case_name)):
        os.makedirs(os.path.join(profile_dir, case_name))   # Makes directory if not existing

    # Read in the location information table
    projectlist = pd.read_csv(os.path.join(location_data_dir,case))     # Returns data frame of data in csv
    zones = np.unique(projectlist['zone'].values)
    states = np.unique(projectlist['state'].values)

    # Initialize lists/dicts that will get filled while looping through wind sites from project list
    if exclude_certain_ids: # unless you want to exclude ids from existing wind, overwrite list to be zero
        used_ids = pd.read_csv(os.path.join(used_ids_dir, "used_ids.csv"), header=None)[0].values.tolist()
    else:
        used_ids = []    # to keep track of which toolkit ids are used, avoid using them more than once
    wind_by_zone = defaultdict(dict)        # Returns a dictionary-like object
    production_shortage = []    # when not enough sites are found for each project
    summary_of_sites = pd.DataFrame(columns=['project name',
                                             'zone',
                                             'state',
                                             'latitude',
                                             'longitude',
                                             'target ac capacity (mw)',
                                             'target energy (mwh)',
                                             'share',
                                             'turbine name',
                                             'hub height (m)',
                                             'offshore',
                                             'shape capacity (mw)',
                                             'shape energy (mwh)',
                                             'target cf',
                                             'shape cf',
                                             'wind toolkit cf'
                                             ])
    # This creates an empty data frame with relevant column names

    # Start looping through wind farm projects in the project list

    for i, row in projectlist.iterrows():

        # Extract the project parameters
        project_name = row['project name']
        zone = row['zone']
        state = row['state']
        lat = row['latitude']
        lon = row['longitude']
        target_capacity = row['ac capacity (mw)']
        target_energy = row['energy (mwh)']
        share = row['share']        # NOTE_TO_SELF: What does share do?
        turbine_name = row['turbine name']
        hub_height = row['hub height (m)']
        offshore = row['offshore']  # flag to make sure we don't mix on and off-shore

        print '\nprocessing project ... ' + project_name
        print project_name, zone, state, lat, lon, target_capacity, target_energy, share

        # If no turbine name is given, use NREL power output data rather than applying power curve
        # (temporarily, only for this site)
        if pd.isnull(projectlist).loc[i,'turbine name']:
            use_NREL_power_temp = True
            print 'Using generic turbine shape from NREL'
        else:
            use_NREL_power_temp = False
            print 'Turbine type provided: ' + turbine_name
            tck, cutout, max_power = curve_dict[turbine_name]  # grab the power curve interpolation and cutout speed

        if not wind_by_zone.has_key(zone):
            wind_by_zone[zone] = np.zeros(len(dates))    # initializes a list as a value
            # NOTE_TO_SELF: Is dict the best class to use for wind_by_zone?

        if not (share and target_energy):
            continue

        # initialize results
        final_outputs = {}
        wind_production = np.zeros(len(dates))
        wind_production_WTK = np.zeros(len(dates))

        # list of state variables that we update through the for loop until we break out of it when hitting the target
        id_list, energy, capacity, enough_sites = [], 0, 0, False

        # Loop through list of toolkit ids ranked by distance from wind site until you hit energy target
        lon_lat_array = np.array(site_meta[['longitude', 'latitude']])   # Of all the 126,000 sites
        distance = get_distance(lat, lon, lon_lat_array)    # Array of distances for all 126,000 sites
        for toolkit_id in np.argsort(distance):     # np.argsort returns indices in a sorted order (ascending)
            # NOTE: This is because row_id is same as toolkit_id

            # Read in relevant toolkit id meta data
            site_cap_factor = site_meta.loc[toolkit_id,'capacity_factor']
            site_capacity = site_meta.loc[toolkit_id,'capacity']
            site_lat = site_meta.loc[toolkit_id,'latitude']
            site_lon = site_meta.loc[toolkit_id,'longitude']

            # n_turbines = np.around(site_capacity/2.) # NREL assumes 2 MW turbines
            # above is old code when trying to limit n_turbines to 8. Problem when there are 100s of small ones
            # now we assume capacity per grid cell is independent of turbine technology, only n_turbines varies

            try:  # 1,2,3
                NREL_turbine_type = 'NREL_' + str(int(site_meta.loc[toolkit_id, 'power_curve']))
            except:  # offshore # NOTE_TO_SELF: If try block has an error, this is passed to except
                NREL_turbine_type = 'NREL_' + str(site_meta.loc[toolkit_id, 'power_curve'])

                # NOTE_TO_SELF: Why is this needed?

            if force_NREL_powercurves:
                tck,cutout,max_power = curve_dict[NREL_turbine_type]

            n_turbines = site_capacity / max_power # float, not necessarily whole number

            # Do some checks and get out of loop or go to next iteration if needed
            if distance[toolkit_id]>MAX_DISTANCE:
                break # get out of for loop of wind toolkit ids as soon as you go over max distance
                # (in this case enough_sites flag will be false and there will be an energy shortage recorded)
            if toolkit_id in used_ids or site_cap_factor < MIN_CAP_FACTOR:
                continue # if toolkit id already used, skip this one and go to next
            if (NREL_turbine_type == 'NREL_offshore' and not offshore) or (NREL_turbine_type != 'NREL_offshore' and offshore):
                continue # if offshore status of toolkit id does not match project's, skip
                # (don't want to grab onshore profile for offshore or vice versa)

            # Read in or download toolkit data for the site of interest (if available)
            data = []
            for year in np.arange(start_year, end_year + 1):
                temp_data, metadata = load_data.get_wtk_data(year, toolkit_id, site_lat, site_lon, wtk_data_dir, interval= resolution)
                # NOTE_TO_SELF: If not using metadata, remove it?
                # WARNING: defaults to UTC = false, so if you want data in UTC, make it true!
                data.append(temp_data)
            data = pd.concat(data)      # Concatenates all data

            wind_speed = data['wind speed at 100m (m/s)'].values
            density = data['density at hub height (kg/m^3)'].values
            measure_hub_height_wtk = 100  # measure height of the wind speed data in WTK

            '''
             Wind profile power law:

             u / u_0 = (z / z_0) ^ alpha 

             where 
             u = wind velocity
             z = height
             _0 = reference
             alpha = stability coefficient = 0.143 for onshore and 0.11 for offshore 

            '''
            if offshore:
                stability_coefficient = 0.11
            else:
                stability_coefficient = 0.134

            if use_NREL_power or use_NREL_power_temp:
                power = data['power (MW)']
                adjusted_wind_speed = wind_speed
            else:
                (adjusted_wind_speed, power) = wind_calc(wind_speed, density, measure_hub_height_wtk, hub_height, tck, cutout,
                                                         n_turbines, stability_coefficient = stability_coefficient,
                                                         hysteresis=True, dens_adj=False)
            id_energy = power.mean() * 8766
            # NOTE_TO_SELF: This is better than adding because it takes care of average year (leap and non-leap)

            if scale_to_energy:
                # Update state variables and break out of loop if needed
                if (energy+id_energy) >= target_energy and energy:
                    enough_sites = True
                    print 'Reached required energy of ' + str(target_energy) + ' MWh'
                    if abs(target_energy - (energy + id_energy)) < abs(target_energy - energy):
                        # check whether adding the last toolkit id brings us closer to target or not. If yes, add it too
                        # NOTE_TO_SELF: Why is this check even needed?
                        # NOTE_TO_SELF: Please make this a function instead of repeating 4 times
                        id_list.append(toolkit_id)
                        energy += id_energy
                        capacity += site_capacity
                        wind_production += power
                        wind_production_WTK += data['power (MW)']
                        print 'added toolkit id ...' + str(toolkit_id)
                    break  # WE ARE DONE, ENERGY TARGET IS MET, BREAK OUT OF FOR LOOP
                else:  # the standard case: add another toolkit id to the wind site and go to next loop
                    id_list.append(toolkit_id)
                    energy += id_energy
                    capacity += site_capacity
                    wind_production += power
                    wind_production_WTK += data.variables['power'][:]
                    print 'added toolkit id ...' + str(toolkit_id)
                    print 'Reached ' + str(energy) + ' MW out of target ' + str(target_energy) + ' MW (' \
                          + str(energy * 1. / target_energy * 100) + '%)'

            else:   # scale to capacity
                # Update state variables and break out of loop if needed
                if (capacity+site_capacity) >= target_capacity and capacity:
                    enough_sites = True
                    print 'Reached required capacity of ' + str(target_capacity) + ' MW'
                    if abs(target_capacity-(capacity+site_capacity)) < abs(target_capacity - capacity):
                        # NOTE_TO_SELF: Is this what we want? (If target = 20 and site = 16, this takes only 1 site)
                        # check whether adding the last toolkit id brings us closer to target or not. If yes, add it too
                        id_list.append(toolkit_id)
                        energy += id_energy
                        capacity += site_capacity
                        wind_production += power
                        wind_production_WTK += data['power (MW)']
                        print 'added toolkit id ...' + str(toolkit_id)
                    break # WE ARE DONE, CAPACITY TARGET IS MET, BREAK OUT OF FOR LOOP
                else: # the standard case: add another toolkit id to the wind site and go to next loop
                    id_list.append(toolkit_id)
                    energy += id_energy
                    capacity += site_capacity
                    wind_production += power
                    wind_production_WTK += data['power (MW)']
                    print 'added toolkit id ...' + str(toolkit_id)
                    print 'Reached ' + str(capacity) + ' MW out of target ' + str(target_capacity) + ' MW (' \
                        + str(capacity * 1. / target_capacity * 100) + '%)'

        # --- end of loop through toolkit IDs (arrive here after break statement) ----- #

        used_ids += id_list  # update used _ids list (done after each project)

        # Scale up profile to target capacity (default) or energy
        if enough_sites:
            if scale_to_energy:
                embedded_scale_factor = (target_energy*1./energy)
            else:
                embedded_scale_factor = (target_capacity*1./capacity)
        else:
            embedded_scale_factor = 1.
            production_shortage.append([project_name, zone, lat, lon]+[(target_energy-energy)*share])

        # Correct output profiles for timezone (assumes WTK is in UTC; not necessarily true!)
        wind_production = correct_wind_time(wind_production, nperiods=timezone)*(share*embedded_scale_factor)
        wind_production_WTK = correct_wind_time(wind_production_WTK, nperiods=timezone)*(share*embedded_scale_factor)

        # Calculate and print capacity factor
        WTK_CF = np.mean(wind_production_WTK) / (capacity*embedded_scale_factor) # will be nan if no toolkit sites found
        calc_CF = np.mean(wind_production) / (capacity*embedded_scale_factor)
        print 'WTK cap factor for ' + project_name + ' is ... ' + str(WTK_CF)   # NOTE_TO_SELF: Why print string?
        print 'E3 calculated cap factor for ' + project_name + ' is ... ' + str(calc_CF)

        # Put data into data frames and put those dfs in final_outputs dict:
        final_outputs['E3'] = pd.DataFrame({'output (MW)': wind_production})
        final_outputs['WTK'] = pd.DataFrame({'output (MW)': wind_production_WTK})

        # Add generation profile to running total by zone
        wind_by_zone[zone] += wind_production

        target_CF = target_energy / (target_capacity * 8766.0)

        # add entry into summary of sites data frame
        summary_of_sites.loc[i] = [project_name, zone, state, lat, lon,
                                   target_capacity, target_energy, share, turbine_name, hub_height, offshore,
                                   capacity*embedded_scale_factor, energy*embedded_scale_factor,target_CF,calc_CF,WTK_CF]
        # NOTE_TO_SELF: Understand what loc does

        # write out outputs for each project in the projectlist (if toggle is on)
        if output_project_shapes:
            for _tag in ['E3', 'WTK']:
                df = final_outputs[_tag]
                if include_date:
                    df.index = dates
                    df.index.name = 'date (tz = ' + str(timezone) + ')'     # This is the name of index column

                tag = '_WTK' if _tag == 'WTK' else ''

                if start_year != end_year:
                    filename = '_'.join([project_name+tag,str(start_year),str(end_year)]) + '.csv'
                    # NOTE_TO_SELF: This is a nifty way to name files
                else:
                    filename = '_'.join([project_name+tag,str(start_year)]) + '.csv'
                filepath = os.path.join(profile_dir,case_name,filename)

                if _tag == 'WTK' and not output_WTK_RT:
                    continue # don't output WTK if toggle is off
                else:
                    df.to_csv(filepath, index=include_date)   # NOTE_TO_SELF: index = true writes row names


    # ------ adjust zonal output for energy shortages ---------------- #

    # Adjust zonal outputs for energy shortages (if any)
    if adjust_for_energy_shortage:
        if len(production_shortage) >= 1:
            print "adjusting outputs by zone for production shortage"
            for shortage in production_shortage:
                energy_shortage = shortage[4]
                zone = shortage[1]
                average_annual_energy_by_zone = np.mean(wind_by_zone[zone])*8766
                wind_by_zone[zone] = wind_by_zone[zone] * (1 + energy_shortage / average_annual_energy_by_zone) # scale up the whole zonal profile to make up for adjustment

    # ------------------- Write out outputs -------------------------- #

    # Write out aggregated outputs: aggregated by zone

    for zone in wind_by_zone.keys():
        df = pd.DataFrame({'output (MW)':wind_by_zone[zone]})
        if include_date:
            df.index = dates
            df.index.name = 'date (tz = ' + str(timezone) + ')'
        if start_year != end_year:
            filename = '_'.join(['aggregation',zone,str(start_year),str(end_year)]) + '.csv'
        else:
            filename = '_'.join(['aggregation',zone,str(start_year)]) + '.csv'

        filepath = os.path.join(profile_dir,case_name,filename)
        df.to_csv(filepath,index=include_date)

    # Write out aggregated outputs: full aggregation
    if start_year != end_year:
        filename = '_'.join(['aggregation_all_zones',str(start_year),str(end_year)]) + '.csv'
    else:
        filename = '_'.join(['aggregation_all_zones',str(start_year)]) + '.csv'

    sum_case_wind = np.sum([wind_by_zone[zone] for zone in wind_by_zone.keys()], axis = 0)
    df_all_hourly = pd.DataFrame({'output (MW)': sum_case_wind})

    if include_date:
        df_all_hourly.index = dates
        df_all_hourly.index.name = 'date (tz = ' + str(timezone) + ')'
    df_all_hourly.to_csv(os.path.join(profile_dir,case_name,filename), index=include_date)

    # NOTE_TO_SELF: What if !include_date??
            
    # Write out the summary of sites to for validation afterwards
    # (this is a file very similar to the location information input file in terms of layout)
    summary_of_sites.to_csv(os.path.join(profile_dir, case_name, 'summary_of_sites_' + case_name + '.csv'), index=False)    
    
    # Write out energy shortage outputs
    with open(os.path.join(profile_dir, case_name, 'wind_toolkit_target_energy_shortage.csv'), 'wb') as outfile:
        csv_writer = csv.writer(outfile, delimiter = ',')
        csv_writer.writerow(['Locations with energy shortage (MWh)', 'Zone', 'Latitude', 'Longitude', 'Energy shortage (MWh)'])
        for row in production_shortage:
            csv_writer.writerow(row)
    
    # Write out used ids        
    with open(os.path.join(profile_dir, case_name, 'used_ids.csv'), 'wb') as outfile:
        csv_writer = csv.writer(outfile, delimiter = ',')
        for row in used_ids:
            csv_writer.writerow([row])
            