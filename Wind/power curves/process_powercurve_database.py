# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:51:52 2016

This script goes through the power curves database from www.wind-power-program.com
and converts those .pow files in csv files:
    wind_speed: column with wind speed in m/s
    power: column with power output in MW
The script also creates a csv file "metadata.csv" with the metadata:
    turbine_name: turbine name, same as the .pow file name (index)
    cut_in_wind_speed: wind speed in m/s at which turbine will start producing power
    cut_out_wind_speed: wind speed in m/s at which the turbine will cut out its power production
    rotor_diameter: rotor diameter in m
    max_power: maximum power output in MW

@author: Gerrit
"""

import glob
import os
import pandas as pd

def list_files(path):
    """ 
    returns a list of names (with extension, without full path) of all files 
    in folder path (NOT including subdirectories)
    """    
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(name)
    return files
    # NOT USED    
    
def get_filepaths(directory):
    """
    This function will generate the file names in a directory AND its subdirectories
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


directory = os.getcwd()
out_folder = "Databasepowercurves(February2016) - Processed"
all_files = get_filepaths(directory)
all_pow_files = [pow_file for pow_file in all_files if pow_file[-3:] == 'pow']

metadata_dict = {}
for pow_file in all_pow_files:
    turbine_name = pow_file.split('\\')[-1][:-4]
    power_data = []
    with open(pow_file, 'rb') as infile:    
        for i,line in enumerate(infile):
            if i == 3:
                cut_out = float(line.strip().replace("\"",""))
            elif i == 1: 
                rotor_diameter = int(line.strip().replace("\"",""))
            elif i == 4:
                cut_in = float(line.strip().replace("\"",""))
            elif i > 4 and i < 35:
                power_data.append(float(line.strip().replace("\"",""))/1000) # convert to MW by dividing by 1000
    
    max_power = max(power_data)
    metadata_dict[turbine_name] = cut_in, cut_out, rotor_diameter, max_power
    data_dict = {'power': power_data, 'wind_speed': range(1,len(power_data)+1)}    
    pd.DataFrame(data=data_dict).to_csv(os.path.join(directory,out_folder,turbine_name+'.csv'),index=False)

metadata = pd.DataFrame.from_dict(metadata_dict,orient='index').sort_index()
metadata.to_csv(os.path.join(directory,out_folder,'metadata.csv'),index_label ='turbine_name', header=['cut_in_wind_speed','cut_out_wind_speed','rotor_diameter','max_power'])        
    
            