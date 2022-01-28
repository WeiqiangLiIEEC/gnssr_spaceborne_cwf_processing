from __future__ import print_function
import numpy as np
import pandas as pd
#import altair as alt
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import os
import requests   
from LIBRARY_Spaceborne_CW import (
    Integrate,
    Read_CW_File, 
    Create_MetaData_InterpolationFunctions,
    IntegrateCoherently,                      
    IntegrateIncoherently,
    Create_Interp1dFunction, 
    ApplyDelay, 
    ApplyCounterRotation,
    Find_RawIF_Track,
    Download_cWF_File,
    Download_cWF_File_List,
    generate_kml_scatter) 
import sys
if sys.version_info.major == 3:
     from urllib.parse import urlparse
else:
     from urlparse import urlparse    

username = 'username'     # Please replace with your username for the GOLD-RTR MINING Server
password = 'password'     # Please replace with password for the GOLD-RTR MINING Server
cWF_data_dir = 'cWF_Data/'

'''
Download raw IF track list files
'''
Download_cWF_File_List(username, password, cWF_data_dir)

'''
Searching for raw IF tracks crossing a target region
'''
list_tracks_target = 'test_list_tracks.txt'   # The text file including the raw IF tracks crossing the target region
Target_Region = [-82.7, 28.95, -82.6, 29.05]  # A rectangle target region defined by [Lon_min, Lat_min, Lon_max, Lat_max]
Start_time = '2018-09-19T00:00:00'     # Start of the time period for raw IF track searching
Stop_time  = '2018-09-20T00:00:00'     # End of the time period for raw IF track searching
Find_RawIF_Track(Target_Region, Start_time, Stop_time, cWF_data_dir + list_tracks_target, cWF_data_dir)

'''
Download complex waveform files in the list "list_tracks_target"
'''
fid_raw_list = open(cWF_data_dir + list_tracks_target, 'r')
lines1 = fid_raw_list.readlines()

for line in lines1:
    track_info = line
    m1 = track_info.split(",")
    url = m1[2][0:-1]
    cWF_file_name = os.path.basename(urlparse(url).path)
    if os.path.exists(cWF_data_dir + cWF_file_name):
        continue
    Download_cWF_File(url, username, password, cWF_data_dir)

'''
Processing the cWF file in the list"
'''
file_Delay    = 'None'
file_Rotation = 'None'
CohIntSamples = 10
IncohIntSamples = int(1000/CohIntSamples)

for line in lines1:
    track_info = line
    m1 = track_info.split(",")
    url = m1[2][0:-1]
    cWF_file_name = os.path.basename(urlparse(url).path)
    file_CW = cWF_data_dir + cWF_file_name
    print ('Processing complex waveform file: ' + cWF_file_name + '...')
    ds_MetaData, ds_cWF, UTC_start, SoW_start = Read_CW_File(file_CW)
    f_Lat_SP, f_Lon_SP, f_Alt_SP, f_incidence,  f_LanOrOcean = Create_MetaData_InterpolationFunctions(ds_MetaData)
    print ('Coherent integration time: %d ms, number of incoherent average: %d.'%(CohIntSamples, IncohIntSamples))
    WF_DW_CoherentInt, DW_Coherence, DW_PeakPhasors, WF_DW_IncoherentInt, Time_tag_CoherentInt, Time_tag_IncoherentInt = Integrate(ds_cWF,UTC_start,'DW',CohIntSamples, IncohIntSamples)
    WF_UP_CoherentInt, UP_Coherence, UP_PeakPhasors, WF_UP_IncoherentInt, Time_tag_CoherentInt, Time_tag_IncoherentInt = Integrate(ds_cWF,UTC_start,'UP',CohIntSamples, IncohIntSamples)
    Lat_coh = f_Lat_SP(Time_tag_CoherentInt)
    Lon_coh = f_Lon_SP(Time_tag_CoherentInt)
    Alt_coh = f_Alt_SP(Time_tag_CoherentInt)
    Incidence_coh = f_incidence(Time_tag_CoherentInt)
    Lat_inc = f_Lat_SP(Time_tag_IncoherentInt)
    Lon_inc = f_Lon_SP(Time_tag_IncoherentInt)
    Alt_inc = f_Alt_SP(Time_tag_IncoherentInt)
    Incidence_inc = f_incidence(Time_tag_IncoherentInt)
    
    kml = generate_kml_scatter('Coherence', 'Coherence Function', Lon_coh, Lat_coh, DW_Coherence, 0, 1.0)
    kml.save(cWF_data_dir + cWF_file_name[0:-2] + 'kml')
