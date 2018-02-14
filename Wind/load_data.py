# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 13:23:02 2016

@author: Gerrit
Edits: Manu
"""

import pandas as pd
import os

# full overview of available wind toolkit data: https://www.nrel.gov/grid/wind-toolkit.html
# note that this focuses on Techno-Economic WIND Toolkit (2 TB 5 min data at 126,000 sites),
# not the Gridded Atmospheric WIND Toolkit (50 TB, Hourly data for all 2x2 km2 sites in US, Baja, Atlantic and Pacific

# WTK API guide: https://developer.nrel.gov/docs/wind/wind-toolkit/wind-toolkit-extract/
# TODO: add script to find closest wind toolkit location for specified lat-lon?

# Gridded Atmospheric WIND Toolkit API guide: https://github.com/NREL/hsds-examples (Not used here)


def get_wtk_data(year,
                 toolkit_id,
                 lat,
                 lon,
                 wtk_data_dir,
                 interval=60,   # TODO: use this to calculate hourly if needed
                 leap_year='true',
                 utc='false',
                 attributes='wind_speed,wind_direction,power,pressure,temperature,density'):
    '''
    Returns WTK data file as a Pandas DF for the year,lat,lon of interest, as well as the meta data file
    First checks whether there already exists a downloaded file. If not, queries the API and downloads file    
        year: string, year of interest
        toolkit_id: #NOTE_TO_SELF: What is this?
        lat: float, latiude of solar site
        lon: float, longitude of solar site
        wtk_data_dir: string, directory where the wind data will be saved
        interval: int, Set time interval in minutes, i.e., '30' is half hour intervals. Valid intervals are 5 & 60.
        leap_year: str, Set leap year to true or false. True will return leap day data if present, false will not.
        utc: str, Specify Coordinated Universal Time (UTC), 'true' will use UTC, 'false' will use the local time zone of the data.
            NOTE: In order to use the WTK data in SAM, you must specify UTC as 'false'. SAM requires the data to be in the
            local time zone.
        attributes: str, Set the attributes to extract (e.g., wind_speed), separated by commas.
    '''

    ## ------------------ LOGIN DATA + API KEY ---------------------- ##

    # You must request an WTK api key from the link above
    api_key = 'zUpLu4x7i0kTbKKN2UE1SnWOdeLLAjwCmUeSIVKb'
    # Your full name, use '+' instead of spaces.
    your_name = 'gerrit+de+moor'
    # Your reason for using the WTK.
    reason_for_use = 'elcc+analysis'
    # Your affiliation
    your_affiliation = 'e3'
    # Your email address
    your_email = 'gerrit@ethree.com'
    # Please join our mailing list so we can keep you up-to-date on new developments.
    mailing_list = 'false'
    
    # ----------------- ADJUST INPUTS ------------------------------ ##
        
    # Declare url string
    url = 'http://developer.nrel.gov/api/wind-toolkit/wind/wtk_download.csv?wkt=POINT({lon}%20{lat})&names={year}&' \
          'leap_day={leap}&utc={utc}&full_name={name}&email={email}&affiliation={affiliation}&mailing_' \
          'list={mailing_list}&reason={reason}&api_key={api}&attributes={attr}'\
        .format(year=str(year),
                lat=lat,
                lon=lon,
                leap=leap_year,
                utc=utc,
                name=your_name,
                email=your_email,
                mailing_list=mailing_list,
                affiliation=your_affiliation,
                reason=reason_for_use,
                api=api_key,
                attr=attributes
                )
    # NOTE_TO_SELF: why redefine year, lat, etc.?
#    print url # TODO: remove
    
    # ----------------- DOWNLOAD DATA ------------------------------ ##    
    
    # Download metadata if not existing
    meta_filepath = os.path.join(os.getcwd(),wtk_data_dir,'_'.join([str(toolkit_id),'meta.csv']))
    if not os.path.isfile(meta_filepath):
        # if no wtk file is present, query wtk for the year, lat, lon of interest and save it to wtk_data_dir folder
        info = pd.read_csv(url, nrows=2) # Return just the first 3 lines to get metadata:   
        info.to_csv(meta_filepath)      
    else:
        info = pd.read_csv(meta_filepath)    
    
    # Download actual radiation data if not existing
    start_date = str(year) + "/01/01 00:00:00"  # NOTE_TO_SELF: replace with pd.Timestamp()
    # start_date = pd.Timestamp(year)

    end_date = str(year) + "/12/31 23:55:00"

    time_int = str(interval) + 'min'

    min_filepath = os.path.join(os.getcwd(), wtk_data_dir, '5min')
    hour_filepath = os.path.join(os.getcwd(), wtk_data_dir, 'Hourly')
    if not os.path.exists(min_filepath):
        os.makedirs(min_filepath)   # Makes '5min' directory
    if not os.path.exists(hour_filepath):
        os.makedirs(hour_filepath)  # Makes 'Hourly' directory

    wtk_5min_filepath = os.path.join(os.getcwd(), wtk_data_dir, '5min',
                                     '_'.join([str(toolkit_id), str(year), '5min', 'wtk.csv']))
    wtk_hourly_filepath = os.path.join(os.getcwd(), wtk_data_dir, 'Hourly',
                                       '_'.join([str(toolkit_id), str(year), 'Hourly', 'wtk.csv']))

    if not os.path.isfile(wtk_5min_filepath):
        print 'no wtk data found locally for toolkit_id... ' + str(toolkit_id) + ', year ' + str(year) + ' ... downloading data from website'
        # if no wtk file is present, query wtk for the year, lat, lon of interest and save it to wtk_data_dir folder
        df = pd.read_csv(url, skiprows=3) # Return all but first 3 lines of csv to get data:
        df = df.set_index(pd.date_range(start_date,end_date, freq='5Min')) # Set the time index in the pandas dataframe
        df.to_csv(wtk_5min_filepath)
        df_hourly = df.resample('60min').mean()
        # https://stackoverflow.com/questions/35339139/where-is-the-documentation-on-pandas-freq-tags
        # if interval == 60:
        #    df.drop(['Minute'],axis=1,inplace=True)
        df_hourly.to_csv(wtk_hourly_filepath)
        df = df.resample(time_int).mean()
    else: 
        print 'found local wtk data for toolkit_id... ' + str(toolkit_id) + ', year ' + str(year) +  '!'
        df = pd.read_csv(wtk_5min_filepath, index_col=0, parse_dates=True)
#        df = df.set_index(pd.date_range(start_date,end_date, freq='5Min')) # Set the time index in the pandas dataframe
        df = df.resample(time_int).mean()

    return df, info
    
# background on requests module: http://docs.python-requests.org/en/master/user/quickstart/    
# import requests
# result = requests.request("GET", url)
# result.text (very slow because huge text string!)
    
    
#df, info = get_wtk_data(year=2009,toolkit_id=01,lat=39.90973623453719,lon=-104.23828125,wtk_data_dir="WTK data")

