# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 13:23:02 2016

@author: Gerrit
"""

import pandas as pd
import os

# for more info: https://nsrdb.nrel.gov/api-instructions

# FROM PVSIM: NECESSARY WEATHER DATA NEEDED
#        To do a PV simulation, the following parameters are needed:
#            dni (ndarray): direct normal irradiance (W/m^2)
#            ghi (ndarray): global horizontal irradiance (W/m^2)
#            dhi (ndarray): diffuse horizontal irradiance (W/m^2)
#            sun_zenith (ndarray): height of the sun from the horizon (degrees)
#            sun_azimuth (ndarray): angle clockwise from north (degrees)
#            temperature (ndarray): ambient temperature (C)
#            wind_speed (ndarray): wind speed at panel array (m/s)

# Declare all variables as strings. Spaces must be replaced with '+', i.e., change 'John Smith' to 'John+Smith'.

def get_nsrdb_data(year,
                   lat,
                   lon,
                   rad_data_dir,
                   interval=60,
                   leap_year = 'true',
                   utc = 'false',
                   attributes = 'ghi,dhi,dni,wind_speed_10m_nwp,surface_air_temperature_nwp,solar_zenith_angle'):
    '''
    Returns NSRDB data file as a Pandas DF for the year,lat,lon of interest, as well as the meta data file
    First checks whether there already exists a downloaded file. If not, queries the API and downloads file    
        year: string, year of interest
        lat: float, latiude of solar site
        lon: float, longitude of solar site
        rad_data_dir: string, directory where the radiation data will be saved
        interval: int, Set time interval in minutes, i.e., '30' is half hour intervals. Valid intervals are 30 & 60.
        leap_year: str, Set leap year to true or false. True will return leap day data if present, false will not.
        utc: str, Specify Coordinated Universal Time (UTC), 'true' will use UTC, 'false' will use the local time zone of the data.
            NOTE: In order to use the NSRDB data in SAM, you must specify UTC as 'false'. SAM requires the data to be in the
            local time zone.
        attributes: str, Set the attributes to extract (e.g., dhi, ghi, etc.), separated by commas.
    '''

    ## ------------------ LOGIN DATA + API KEY ---------------------- ##

    # You must request an NSRDB api key from the link above
    api_key = 'zUpLu4x7i0kTbKKN2UE1SnWOdeLLAjwCmUeSIVKb'
    # Your full name, use '+' instead of spaces.
    your_name = 'gerrit+de+moor'
    # Your reason for using the NSRDB.
    reason_for_use = 'elcc+analysis'
    # Your affiliation
    your_affiliation = 'e3'
    # Your email address
    your_email = 'gerrit@ethree.com'
    # Please join our mailing list so we can keep you up-to-date on new developments.
    mailing_list = 'false'
    
    # ----------------- ADJUST INTPUTS ------------------------------ ##
    
    lat = round(lat,2) # NSRDB data does not go beyodn 2 significant figures. i.e. data file for lat = 12.1292 is same as lat = 12.13
    lon = round(lon,2) # see above
    
    # Declare url string
    url = 'http://developer.nrel.gov/api/solar/nsrdb_0512_download.csv?wkt=POINT({lon}%20{lat})&names={year}&leap_day={leap}&interval={interval}&utc={utc}&full_name={name}&email={email}&affiliation={affiliation}&mailing_list={mailing_list}&reason={reason}&api_key={api}&attributes={attr}'.format(
        year=str(year), 
        lat=lat, # nsrdb precision is only up to 1/100 degrees (pixels are 0.038 x 0.038 degrees). NSRDB grabs closest pixel, so can be a few 1/100s off. 
        lon=lon, 
        leap=leap_year, 
        interval=str(interval), 
        utc=utc, 
        name=your_name, 
        email=your_email, 
        mailing_list=mailing_list, 
        affiliation=your_affiliation, 
        reason=reason_for_use, 
        api=api_key, 
        attr=attributes)
    
    # ----------------- DOWNLOAD DATA ------------------------------ ##    
    
    # Download metadata if not existing
    meta_filepath = os.path.join(os.getcwd(),rad_data_dir,'_'.join([str(lat),str(lon),'meta.csv']))
    if not os.path.isfile(meta_filepath):
         # if no nsrdb file is present, query nsrdb for the year, lat, lon of interest and save it to rad_data_dir folder
        info = pd.read_csv(url, nrows=1) # Return just the first 2 lines to get metadata:   
        info.to_csv(meta_filepath)      
    else:
        info = pd.read_csv(meta_filepath)    
    
    # Download actual radiation data if not existing
    nsrdb_filepath = os.path.join(os.getcwd(),rad_data_dir,'_'.join([str(lat),str(lon),str(year),'rad.csv']))
    if not os.path.isfile(nsrdb_filepath):  
        print 'no nsrdb data found for ... (' + str(lat) + ',' + str(lon) + ') year ' + str(year) +  ' ... downloading data from website'            
        # if no nsrdb file is present, query nsrdb for the year, lat, lon of interest and save it to rad_data_dir folder
        df = pd.read_csv(url, skiprows=2) # Return all but first 2 lines of csv to get data:
        df.to_csv(nsrdb_filepath)  
    else: 
        print 'found local nsrdb data for ... (' + str(lat) + ',' + str(lon) + ') year ' + str(year) +  '!'
        df = pd.read_csv(nsrdb_filepath, index_col=0)

    start_date = str(year) + "/01/01 00:00:00"
    end_date = str(year) + "/12/31 23:00:00"
    df = df.set_index(pd.date_range(start_date,end_date, freq=str(interval)+'Min')) # Set the time index in the pandas dataframe

    return df,info

