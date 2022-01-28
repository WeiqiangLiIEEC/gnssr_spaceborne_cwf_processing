#!/usr/bin/env python
# coding: utf-8


#The use of the overlays in this library follows the reference:
# https://ocefpaf.github.io/python4oceanographers/blog/2014/03/10/gearth/
# More details could be found there


#=============================================
from __future__ import print_function
import numpy as np
import pandas as pd
import xarray as xr
import os
import matplotlib.pyplot as plt
import matplotlib.animation as  animation
import datetime

from scipy.interpolate import interp1d
from scipy.ndimage.interpolation import shift
import simplekml
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import scipy.io as sio

import matplotlib as mpl
from matplotlib import cm
import simplekml
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import scipy.io as sio
from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)
import requests
import sys
if sys.version_info.major == 3:
    from urllib.parse import urlparse
else:
    from urlparse import urlparse
cmap = cm.rainbow(range(256))
newcmp = ListedColormap(cmap)

#==========================================================
def Download_cWF_File_List(username, password, cWF_data_dir_local):
    """
    Download_cWF_File_List
        Download the lists of the complex waveform file from the GOLD-RTR server
    Arguments:
        username: username for the GOLD-RTR website
        password: password for the GOLD-RTR website
        cWF_data_dir_local:  local directory for the complex waveform files
    Return:
        1: the file is downloaded sucessfully
        0: unable to download the file
    """
    ulr_cygnss = 'https://www.ice.csic.es/research/gold_rtr_mining/data/Spaceborne_GNSSR_RawIFData_OnGroundProcessing/List_CYGNSS_RawIF_Track.txt'
    filename = cWF_data_dir_local + os.path.basename(urlparse(ulr_cygnss).path)
    try:
        print ("Downloading file: " + '"' + os.path.basename(urlparse(ulr_cygnss).path) + '"' + '...')
        r = requests.get(ulr_cygnss, auth=(username,password))
    except:
        print ("Cannot download CYGNSS complex waveform list file")
        return False

    if r.status_code == 200:
        open(filename, 'wb').write(r.content)
        print("File " +  '"' + os.path.basename(urlparse(ulr_cygnss).path) + '"' + ' downloaded!')
    else:
        print("Cannot download CYGNSS complex waveform list file")
        return False

    ulr_tds1 = 'https://www.ice.csic.es/research/gold_rtr_mining/data/Spaceborne_GNSSR_RawIFData_OnGroundProcessing/List_TDS1_RawIF_Track.txt'
    filename = cWF_data_dir_local + os.path.basename(urlparse(ulr_tds1).path)
    try:
        print ("Downloading file: " + '"' + os.path.basename(urlparse(ulr_tds1).path) + '"' + '...')
        r = requests.get(ulr_tds1, auth=(username,password))
    except:
        print ("Cannot download TDS-1 complex waveform list file")
        return False

    if r.status_code == 200:
        open(filename, 'wb').write(r.content)
        print("File " +  '"' + os.path.basename(urlparse(ulr_tds1).path) + '"' + ' downloaded!')
        return True
    else:
        print("Cannot download TDS-1 complex waveform list file")
        return False 
        
        
def isInsideRect(x1, y1, x2, y2, px, py):
    """
    isInsideRect
        To check the number of points in the rectangle
    Arguments: 
        Bottom-left corner of the rectangle (x1, y1)
        Top-right corner of the rectangle (x2, y2)
        Coordinates of the points (px, py)
    Return:
        The number of points in the rectangle
    """
    inRect = (px >= x1)*(px <= x2)*(py >= y1)*(py <= y2)
    return np.sum(inRect)

 
def Find_RawIF_Track(Region, Start_time, Stop_time, TrackListFile, cWF_data_dir_local):    
    """
    Find_RawIF_Track
        Find TDS-1 and CyGNSS Raw IF data tracks within the time period of [start_time, stop_time] and crossing the rectangle region [lon_min, lat_min, lon_max, lat_max]
    Arguments: 
        Region: A rectangle region defined by [lon_min, lat_min, lon_max, lat_max]
        Start_time: YYYY-MM-DDThh:mm:ss in UTC Time
        Stop_time:  YYYY-MM-DDThh:mm:ss in UTC Time
    Return:
        A list of raw IF tracks and their location in the GOLD-RTR server
    """
    f_TrackListFile = open(TrackListFile, 'w')
    dir_cygnss_cwf = 'https://www.ice.csic.es/research/gold_rtr_mining/data/Spaceborne_GNSSR_RawIFData_OnGroundProcessing/CYGNSS_CWF_Products/CYGNSS_Complex_Waveform_Release/'
    dir_tds1_cwf   = 'https://www.ice.csic.es/research/gold_rtr_mining/data/Spaceborne_GNSSR_RawIFData_OnGroundProcessing/TDS1_CWF_Products/TDS1_Complex_Waveform_Release/'
    start_time_obj = datetime.datetime.strptime(Start_time, '%Y-%m-%dT%H:%M:%S')
    stop_time_obj  = datetime.datetime.strptime(Stop_time,  '%Y-%m-%dT%H:%M:%S')
    # Find CYGNSS tracks
    n_track_CYGNSS = 0
    fid_raw_list_tds1 = open(cWF_data_dir_local+ 'List_CYGNSS_RawIF_Track.txt', 'r')
    lines1 = fid_raw_list_tds1.readlines()
    for line in lines1[1::]:
        m1 = line.split(",")
        track_time_str = m1[2]
        track_time_obj  = datetime.datetime.strptime(track_time_str,  '%Y-%m-%dT%H:%M:%S')
        if track_time_obj > stop_time_obj or track_time_obj < start_time_obj:
            continue
            
        Lat_0 = float(m1[5])
        Lon_0 = float(m1[6])
        Lat_1 = float(m1[7])
        Lon_1 = float(m1[8])            
        lat_itp = np.linspace(Lat_0, Lat_1, 1000)
        lon_itp = np.linspace(Lon_0, Lon_1, 1000)
        cross_area = isInsideRect(Region[0], Region[1], Region[2], Region[3], lon_itp, lat_itp)
        if cross_area < 10:
            continue
        n_track_CYGNSS = n_track_CYGNSS + 1
        print ('CYGNSS Raw IF collection: ' + m1[0], "GNSS Satellite: " + m1[1])
        f_TrackListFile.write(m1[0] + ',' + m1[1] + ',' + dir_cygnss_cwf + m1[10])
        
    # Find TDS-1 tracks
    n_track_TDS1 = 0
    fid_raw_list_tds1 = open(cWF_data_dir_local + 'List_TDS1_RawIF_Track.txt', 'r')
    lines1 = fid_raw_list_tds1.readlines()
    for line in lines1[1::]:
        m1 = line.split(",")
        track_time_str = m1[2][0:19]
        track_time_obj  = datetime.datetime.strptime(track_time_str,  '%Y-%m-%dT%H:%M:%S')
        if track_time_obj > stop_time_obj or track_time_obj < start_time_obj:
            continue
        Lat_0 = float(m1[5])
        Lon_0 = float(m1[6])
        Lat_1 = float(m1[7])
        Lon_1 = float(m1[8])            
        lat_itp = np.linspace(Lat_0, Lat_1, 1000)
        lon_itp = np.linspace(Lon_0, Lon_1, 1000)
        cross_area = isInsideRect(Region[0], Region[1], Region[2], Region[3], lon_itp, lat_itp)
        if cross_area < 10:
            continue
        n_track_TDS1 = n_track_TDS1 + 1      
        print ('TDS-1 Raw IF collection: ' + m1[0], "GNSS Satellite: " + m1[1])
        f_TrackListFile.write(m1[0] + ',' + m1[1] + ',' + dir_tds1_cwf + m1[10])
        
    if (n_track_CYGNSS + n_track_TDS1) > 0:
        print ('***********************************************************************************')
        print ("%d CYGNSS tracks and %d TDS-1 tracks found crossing the target region."%(n_track_CYGNSS, n_track_TDS1))
        print ("The list of the tracks has been saved to " + TrackListFile)
        print ('***********************************************************************************')
    else:
        print ('***********************************************************************************')
        print ("No track  found crossing the target region.")
        print ('***********************************************************************************')
 
#==========================================================
def Download_cWF_File(url, username, password, cWF_data_dir_local):
    """
    Download_cWF_File
        Download a complex waveform file from the GOLD-RTR server
    Arguments: 
        url: url of the complex waveform file
        username: username for the GOLD-RTR website
        password: password for the GOLD-RTR website
        cWF_data_dir_local:  local directory for the complex waveform files
    Return:
        1: the file is downloaded sucessfully
        0: unable to download the file
    """
    filename = cWF_data_dir_local + os.path.basename(urlparse(url).path)
    try:
        print ("Downloading file: " + '"' + os.path.basename(urlparse(url).path) + '"' + '...')
        r = requests.get(url, auth=(username,password))
    except:
        print ("Cannot download the complex waveform file")
        return False

    if r.status_code == 200:
        open(filename, 'wb').write(r.content)
        print("File " +  '"' + os.path.basename(urlparse(url).path) + '"' + ' downloaded!')
        return True
    else:
        print("Cannot download the complex waveform file")
        return False  
#==========================================================

def generate_kml_scatter(Name, asp, lon, lat, value, min_value, max_value):
    """
    generate_kml_scatter
    Arguments: 
    Name, label 
    asp, label
    lon,longitude array
    lat,latitude array
    value, value array
    min_value, minimum value to plot
    max_value, maximum value to plot
    Returns 
    kml, kml object
    """
    df=pd.DataFrame({
        'lon':lon,
        'lat':lat,
        'value':value
    })
    df=df.loc[(df['value'] >= min_value) &  (df['value'] <= max_value)]
    if len(df)==0:
        print('ERROR: no points to plot')
        kml = -1
        return kml
    lon = df['lon'].values
    lat = df['lat'].values
    value = df['value'].values
    if len(df)==0:
        print('ERROR: no points to plot')
        kml = -1
        return kml
    
    kml=simplekml.Kml()
    scale = (np.clip(value,min_value,max_value)-min_value)/(max_value-min_value)
    c = newcmp(scale)*255
    for i in range(len(lon)):
        if i ==0:
            pnt = kml.newpoint(name=Name, description = asp + '=' + '%5.2f'%value[i], coords=[(lon[i],lat[i])])
        else:
            pnt = kml.newpoint(name='', description = asp + '=' + '%5.2f'%value[i], coords=[(lon[i],lat[i])])

        pnt.style.labelstyle.color = simplekml.Color.black
        pnt.style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/pal2/icon15.png'
        pnt.style.iconstyle.scale = 0.5  # Icon thrice as big
        pnt.style.iconstyle.color = simplekml.Color.rgb(int(c[i,0]), int(c[i,1]), int(c[i,2]),  200)
      
    return kml
#==========================================================
def generate_kml_line(Name, asp, lon, lat, value, min_value, max_value):
    """
    generate_kml_scatter
    Arguments: 
    Name, label 
    asp, label
    lon,longitude array
    lat,latitude array
    value, value array
    min_value, minimum value to plot
    max_value, maximum value to plot
    Returns 
    kml, kml object
    """
    df=pd.DataFrame({
        'lon':lon,
        'lat':lat,
        'value':value
    })
    
    df=df.loc[(df['value'] >= min_value) &  (df['value'] <= max_value)]
    
    if len(df)==0:
        print('ERROR: no points to plot')
        kml = -1
        return kml
    lon = df['lon'].values
    lat = df['lat'].values
    value = df['value'].values
    scale_height  = 100000/(np.max(value)-np.min(value))
    height_value = (value-np.min(value))*scale_height
    scale = (np.clip(value,min_value,max_value)-min_value)/(max_value-min_value)
    c = newcmp(scale)*255
    kml=simplekml.Kml()
    for i in range(len(lon)-1):
        linestring = kml.newlinestring(name="A Line"+str(i))
        linestring.coords = [(lon[i],lat[i],value[i]), (lon[i+1],lat[i+1],height_value[i+1])]
        linestring.altitudemode = simplekml.AltitudeMode.relativetoground
        linestring.style.linestyle.color = simplekml.Color.rgb(int(c[i,0]), int(c[i,1]), int(c[i,2]),  200)
        
        linestring.extrude = 1
    

    return kml


#==========================================================
def CreateLegend(label):
    """
    CreateLegend
    Arguments: 
    label, label 
   
    Returns 
    saves a png file with a colorbar
    """
    fig = plt.figure(figsize=(1.0, 4.0), facecolor=None, frameon=False)
    ax0 = fig.add_axes([0.0, 0.05, 0.2, 0.9])
    fraction = 1  # .05
    norm = mpl.colors.Normalize(vmin=-3, vmax=99)
    cbar = ax0.figure.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap='rainbow'),
        ax=ax0, pad=.05, 
        label=label,
        fraction=fraction)
    cbar.set_label(label, color='r', labelpad=20)
    ax0.axis('off')
    if os.path.exist('legend.png'):
        os.remove('legend.png')
    plt.savefig('legend.png', 
        facecolor='red',
        transparent=False, 
        format='png')  # Change transparent to True if your colorbar is not on space :)
        
#==========================================================
def overlay_kml(kml,
                LonLatLimits,
                figs, colorbar=None, **kw):
    """
    overlay_kml
    Arguments: 
    kml a kml structure
    LonLatLimits Lon, Lat Limits
    figs,  a list of png figures to overlay
    colorbar, a png figure
   
    Returns 
    Creates a kml file
    """
    llcrnrlon  = LonLatLimits[0]
    llcrnrlat  = LonLatLimits[1]
    urcrnrlon  = LonLatLimits[2]
    urcrnrlat  = LonLatLimits[3]
    """TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw..."""

   
    altitude = kw.pop('altitude', 6e5)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        ground.name = kw.pop('name', 'overlay')
        ground.color = kw.pop('color', '9effffff')
        ground.atomauthor = kw.pop('author', 'ocefpaf')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'Matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.east = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.west = urcrnrlon

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1
    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    if os.path.exists(kmzfile):
        os.remove(kmzfile)
    print('kmzfile',kmzfile)
    kml.savekmz(kmzfile)

#==========================================================
def Integrate(ds_cWF,UTC_start, UD,CohIntSamples, IncohIntSamples):
    if UD == 'DW':
        WF = 1.*ds_cWF['wf_dw_i'].values + 1j*ds_cWF['wf_dw_q'].values
    elif UD == 'UP':
        WF = 1.*ds_cWF['wf_up_i'].values + 1j*ds_cWF['wf_up_q'].values
    else:
        print('ERROR: UD parameter not implemented')
        return -1
    WF_CoherentInt, Coherence, PeakPhasors, Time_tag_CoherentInt  =  IntegrateCoherently(WF, ds_cWF['Start_time'].values, CohIntSamples)
    WF_IncoherentInt, Time_tag_IncoherentInt = IntegrateIncoherently(WF_CoherentInt, Time_tag_CoherentInt, IncohIntSamples )
    return WF_CoherentInt, Coherence, PeakPhasors, WF_IncoherentInt, Time_tag_CoherentInt, Time_tag_IncoherentInt
#==========================================================
def ShiftCW(CW,Shift):
    """ Shift a Complex Waveform 
    Arguments:
    Complex Waveform CW
    Shift float real in samples
    Returns:
    Shifted CW
    """
    R = np.real(CW)
    I = np.imag(CW)
    R =  shift(R, Shift,order=1,mode='nearest')
    I =  shift(I, Shift,order=1,mode='nearest')
    return R+1j*I 
#==========================================================
def CounterRotateCW(CW,Shift):
    """ Counter rotate  a Complex Waveform 
    Arguments:
    Complex Waveform CW
    Shift float real in cycles
    Returns:
    Shifted wv 
    """
    CW = CW*np.exp(-1j*Shift)
    return CW
#==========================================================
def Read_CW_File(file_CW):
    """ Read file containing Complex Waveforms
    Arguments:
    Input file file_CW
    Returns
    MetaData xarray data set  ds_MetaData
    Comple Waveforms xarray data set  ds_cWF 
    Start UTC string UTC_start 
    Delay Resolution in meters DelayResolution 
    wavelength in meters wavelength
    """
    if os.path.isfile(file_CW):
        with xr.open_dataset(file_CW) as ds_all:
            UTC_start = ds_all.variables['UTC'].values
            SoW_start = ds_all.variables['SoW'].values
            print('Reading variables')
        with xr.open_dataset(file_CW, group='MetaData') as ds_MetaData:
            print('Reading MetaData dataset')
        with xr.open_dataset(file_CW, group='cWF') as ds_cWF:
            print('Reading cWF dataset')
        # ------------- Reserved for the second frequency -----------#
        #with xr.open_dataset(file_CW, group='cWF_2nd_Freq') as ds_cWF_2nd_Freq:
        #    print('Reading cWF dataset')
        #return ds_MetaData, ds_cWF, ds_cWF_2nd_Freq, UTC_start
        #print(ds_MetaData, ds_cWF, UTC_start, SoW_start
        return ds_MetaData, ds_cWF, UTC_start, SoW_start
    else:
        print('ERROR: file_CW ', file_CW, ' is not a file') 
        return -1
   
#==========================================================
def Create_MetaData_InterpolationFunctions(ds_MetaData):
    """ Create_MetaData_InterpolationFunctions"
    Arguments:
    Input Metadata xarray data set
    Return
    The following functions f_XX are created, 
    f_Lat_SP f_Lon_SP, f_Alt_SP, f_incidence,  f_LanOrOcean
    """
    MetaTime        =  ds_MetaData['MetaTime'].values
    Lat_SP          =  ds_MetaData['Lat_SP'].values
    Lon_SP          =  ds_MetaData['Lon_SP'].values
    Alt_SP          =  ds_MetaData['Alt_SP'].values
    incidence       =  ds_MetaData['incidence'].values
    LandOrOcean     =  ds_MetaData['LandOrOcean'].values
    f_Lat_SP = interp1d(MetaTime, Lat_SP, kind='cubic',fill_value='extrapolate')
    f_Lon_SP = interp1d(MetaTime, Lon_SP, kind='cubic',fill_value='extrapolate')
    f_Alt_SP = interp1d(MetaTime, Alt_SP, kind='cubic',fill_value='extrapolate')
    f_incidence = interp1d(MetaTime, incidence, kind='cubic')
    f_LanOrOcean = interp1d(MetaTime, LandOrOcean, kind='nearest')
    return f_Lat_SP, f_Lon_SP, f_Alt_SP, f_incidence,  f_LanOrOcean

#=============================================
def IntegrateCoherently(WF, time_tag, CohIntSamples):
    """ Compute the Coherence Function of  the peak of the complex waveforms
    Arguments:
    Complex waveforms WF
    Coherent Integration time, measured in msec CohIntSamples
    Returns
    A np array with the Coherently integrated waveforms, and
    A np array with the computed Coherent Function
    """
    #Integrate cWF to obtain WF integrated   
    N_waveforms_input = np.shape(WF)[0]
    N_waveforms_to_use = int(N_waveforms_input/CohIntSamples)*CohIntSamples
    N_waveforms_coherent = int(N_waveforms_input/CohIntSamples)
    N_lags             = np.shape(WF)[1]
    WF = np.copy(WF[0:N_waveforms_to_use,:])
    WF_time_tag = np.copy(time_tag[0:N_waveforms_to_use])
    
    #Integrate coherently
    WF = WF.reshape(N_waveforms_coherent,CohIntSamples,N_lags)
    WF_CoherentInt = WF.mean(axis=(1))
    WF_time_tag = WF_time_tag.reshape(N_waveforms_coherent,CohIntSamples)
    Time_tag_CoherentInt = WF_time_tag.mean(axis=(1))
    #Compute Coherence
    #pos_max_power_WF_CoherentInt = np.empty(N_waveforms_coherent, dtype=int)
    #for i in range(N_waveforms_coherent):
    #    pos_max_power_WF_CoherentInt[i] = np.argmax(np.abs(WF_CoherentInt[i,:]))
    #lag with the max power
    pwf = np.mean(np.abs(WF_CoherentInt)**2.0,axis=0)
    pos_max_power_WF = np.argmax(pwf)
    
    phasors = np.empty((N_waveforms_coherent,CohIntSamples), dtype= complex)
    PeakPhasors = np.empty((N_waveforms_coherent), dtype= complex)
    for i in range(N_waveforms_coherent):
        phasors[i,:] = WF[i,:,pos_max_power_WF]/np.abs(WF[i,:,pos_max_power_WF])
        PeakPhasors[i] = np.mean(WF[i,:,pos_max_power_WF])
    Coherence = np.abs(phasors.mean(axis=(1)))
    
    phasors = np.empty((N_waveforms_coherent, CohIntSamples), dtype= complex)
    #Compute phasors peak
    
    return WF_CoherentInt, Coherence, PeakPhasors, Time_tag_CoherentInt
#==========================================================
def IntegrateIncoherently(WF_CoherentInt, Time_tag_CoherentInt, IncohIntSamples ):

    """ Compute the Coherence Function of  the peak of the complex waveforms
    Arguments:
    Complex waveforms WF_CoherentInt
    Incoherent Integration time, measured in duration of coherent waveforms, IncoIntSamples
    Returns
    A np array with the integrated Waveforms, normalized to the floor 
    """
    #Integrate cWF to obtain WF integrated
    N_waveforms_coherent   = np.shape(WF_CoherentInt)[0]
    N_lags                 = np.shape(WF_CoherentInt)[1]
    N_waveforms_uncoherent = int(N_waveforms_coherent/IncohIntSamples)
    N_waveforms_to_use     = N_waveforms_uncoherent*IncohIntSamples
    WF_CoherentInt =  np.copy(WF_CoherentInt[0:N_waveforms_to_use,:])
    Time_tag_CoherentInt = np.copy(Time_tag_CoherentInt[0:N_waveforms_to_use])
    
    y = np.abs(WF_CoherentInt)**2
    y = y.reshape(N_waveforms_uncoherent,IncohIntSamples,N_lags)
    Time_tag_CoherentInt = Time_tag_CoherentInt.reshape(N_waveforms_uncoherent,IncohIntSamples)
    WF_IncoherentInt = y.mean(axis=(1))
    Time_tag_IncoherentInt = Time_tag_CoherentInt.mean(axis=(1))
    
    Floor =  np.median(WF_IncoherentInt[:,0:20])
    WF_IncoherentInt = WF_IncoherentInt/Floor
    
    return WF_IncoherentInt, Time_tag_IncoherentInt
#==========================================================
def Create_Interp1dFunction(File):
    """Create an interpolator function
    Arguments: 
    file name containg two columns, x and y; x should be sorted
    Returns:
    if succesful returns interpolator y=f(x)
    otherwise, returns -1
    """
    if os.path.isfile(File):
        my_csv = np.genfromtxt(File)
        f = interp1d(my_csv[:,0], my_csv[:,1])
        return f
    else:
        print('ERROR: file ', File, ' is not a file') 
        return -1
#==========================================================   
def ApplyDelay(wf, time_tag, f_Delay,DelayResolution):
    """Shift the set of CW
    Arguments:
    Complex Waveforms wf 
    Time tags associated to wf time_tag
    Function to associate to each WF a delay in meters f_Delay
    Delay resolution in meters DelayResolution
    Returns:
    Shifted wv
    """
    ShiftArray = f_Delay(time_tag)/DelayResolution
    N_waveforms = np.shape(wf)[0]
    for i in range(N_waveforms):
        wf[i,:] = ShiftCW(wf[i,:],ShiftArray[i])
    return wf
#==========================================================
def ApplyCounterRotation(wf, time_tag, f_Rotation, wavelength):
    """ Counter a set of CW
    Arguments:
    Complex Waveforms wf 
    Time tags associated to wf time_tag
    Function to associate to each WF a rotation  in meters
    wavelength   in meters wavelength
    Returns:
    counter rotated wv
    """
    RotationArray = f_Rotation(time_tag)/wavelength
   
    N_waveforms = np.shape(wf)[0]
    for i in range(N_waveforms):
        wf[i,:] =  CounterRotateCW(wf[i,:],RotationArray[i])
    return wf
#==========================================================
