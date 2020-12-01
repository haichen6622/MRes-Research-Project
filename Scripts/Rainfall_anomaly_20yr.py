#!/usr/bin/env python2
# -*- coding: utf-8 -*-

### Import all required packages ###
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt
import glob, os


### Please set working directory here ### -------------------
filepath_BSISO = "/Users/helenleft/Desktop/Data_processing/Data/BSISO/"
filepath_rain = "/Users/helenleft/Desktop/Data_processing/Data/TRMM_3B42_V7_Daily/"
filepath_plot = "/Users/helenleft/Desktop/Data_processing/Plots/Rainfall_anomaly/"

#############################################
### Read in the text file of BSISO phases ###
#############################################

#Change working directory
os.chdir(filepath_BSISO)
# Skip one line of header
BSISO_txt = np.genfromtxt("BSISO_25-90bpfil_pc_cdr_20yrs.txt", skip_header=1)
# Get phases and amplitude of BSISO
BSISO_phase = BSISO_txt[:,5]
BSISO_amp = BSISO_txt[:,6]

# Create a datetime array that is identical to BSISO data 
# Start with the first datetime
base = dt.date(1999, 1, 1)
# Iterate through the days, adding 1 day each time
BSISO_dates = np.array([base + dt.timedelta(days=x) for x in range(0, np.size(BSISO_phase))])

### Locate positions and extract time, phase and amp of BSISO during JJAS
# Extract month of JJAS
BSISO_months = [BSISO_dates[x].month for x in range(0, np.size(BSISO_dates))]
BSISO_months = np.array(BSISO_months)
JJAS_days = np.where(((BSISO_months==6) | (BSISO_months==7) | (BSISO_months==8)| (BSISO_months==9)))    

# Extract date, phase and amplitude data for JJAS days during 1999-2018
BSISO_dates_JJAS = BSISO_dates[JJAS_days[0]]
BSISO_amp_JJAS = BSISO_amp[JJAS_days[0]]
BSISO_phase_JJAS = BSISO_phase[JJAS_days[0]]

#################################
### Read in the rainfall data ###
#################################
#Change working directory
os.chdir(filepath_rain)

# Get the rainfall data file names,and sort by name
filenames = sorted(glob.glob('*.nc4'))

# Open first file to have a look what is in file
nc_fid = Dataset(filenames[0],'r') 
print(nc_fid.variables) 
# print(nc_fid.variables.keys()) # Have a look what is in file - variable names only

# Extract lat and lon from 1st file
lats = nc_fid.variables['lat'][:]
lons = nc_fid.variables['lon'][:]

# Set up an empty array to fill with all precipitation data
ppt = np.empty((np.size(filenames),np.size(lons),np.size(lats)))

# Set up an empty array to fill with times
rain_dates = [] 

### Loop through all the rainfall files and extract the rainfall data ###
count = 0 # Use to count the place of rainfall data in ppt array
for f in filenames:
    # Open file
    nc_fid = Dataset(f,'r') 
    
    # Read data and add to ppt array
    ppt[count,:,:] = nc_fid.variables['precipitation'][:] 
   
    # Retrieve date from file name and add to array rain_dates
    temp_date = dt.datetime.strptime(f[11:18+1], "%Y%m%d")
    rain_dates.append(temp_date)
    
    # Increment the count
    count = count+1 

# Check missing value in ppt and change to nan
ppt[np.where(ppt == -9999)] = np.nan
 
# Create an empty array to store precipitation data during JJAS
ppt_JJAS = np.empty((len(JJAS_days[0]),np.size(lons),np.size(lats)))

### Put precipitation during JJAS 1999-2018 into ppt_JJAS ###
count_JJAS = 0 # Use to count the place of rainfall data in ppt_JJAS array
for n in JJAS_days[0]:
    #print(n)
    ppt_JJAS[count_JJAS,:,:] = ppt[n,:,:]
    count_JJAS = count_JJAS + 1 # Increment the count_JJAS


##########################################################
### Calculate mean rainfall during active BSISO phases ###
##########################################################

# Create empty array to store mean rainfall during each BSISO phase
ppt_by_phase = np.empty((8,np.size(lons),np.size(lats)))

### Loop through each BSISO phase and extract rainfall data for each phase ###
for pha in np.arange(1,8+1):
    
    # Work out how many days fit criteria of active BSISO 
    num_days = 0
    for dd in np.arange(0,len(BSISO_dates_JJAS)):
        if BSISO_phase_JJAS[dd]==pha and BSISO_amp_JJAS[dd]>=1:
            num_days = num_days+1

    # Create temporary array to store rainfall data
    ppt_temp = np.empty((num_days,np.size(lons),np.size(lats))) 

    # Locate days during JJAS when BSISO_phase=mp and BSISO_amp>=1
    count = 0
    for dd in np.arange(0,len(BSISO_dates_JJAS)):
        
        # Add rainfall data to temporary array 
        if BSISO_phase_JJAS[dd]==pha and BSISO_amp_JJAS[dd]>=1: 
            ppt_temp[count,:,:] = ppt_JJAS[dd,:,:] 
            count = count + 1

    # Average rainfall data in temporary array and add to ppt_by_phase array
    ppt_by_phase[pha-1,:,:] = np.nanmean(ppt_temp, axis=0) 
        

#################################################################
### Create plots of rainfall anomaly during each BSISO phases ###
#################################################################    
# Compute mean rainfall for computation of rainfall anomaly by BSISO phase in plotting loop below
ppt_JJAS_mean = np.mean(ppt_JJAS, axis=0)

# Create plots of rainfall anomaly, 1 for each phase
for pha in np.arange(1,8+1):
    
    # Initiate new Figure
    fig = plt.figure(figsize=(6,8))

    # Set up axes and select map projection
    ax1 = plt.subplot(1,1,1,projection=ccrs.PlateCarree())
    
    # Add features
    ax1.add_feature(cfeature.COASTLINE)
    ax1.add_feature(cfeature.BORDERS)

    # Find anomally from long-term mean
    anoa =  ppt_by_phase[pha-1,:,:] - ppt_JJAS_mean
    
    # Permute the dimensions of anoa array
    data2plot = np.transpose(anoa)

    # Change values outside contour limits to min and max contours limits
    data2plot[data2plot > 5] = 5
    data2plot[data2plot < -5] = -5
    
    # Plot rainfall anomaly 
    im1 = ax1.contourf(lons, lats, data2plot, np.arange(-5,5+1,1), transform=ccrs.PlateCarree(), cmap=plt.cm.bwr)

    # Add a title
    ax1.set_title('Rainfall anomaly during June-September in 1999-2018:\n BSISO phase:' + str(pha), fontsize=14)
    
    # Add latitude and longitude values
    ax1.set_xticks(np.arange(65,95+1,5), crs=ccrs.PlateCarree())
    ax1.set_yticks(np.arange(5,40+1,5), crs=ccrs.PlateCarree())
    
    # # Label the axes
    ax1.set_xlabel("Longitude", fontsize=12)
    ax1.set_ylabel("Latitude", fontsize=12)
    
    # Add colorbar
    cbar = plt.colorbar(im1,shrink = 0.717,orientation='vertical') 
    # Add colorbar label
    cbar.ax.set_ylabel('mm/day')

    ### Change working directory and save plots
    os.chdir(filepath_plot)
    
    # Save plot
    filename_start = "BSISO_rainfall_anomaly_20yrs_p"
    filename_end = ".png"
    plt.savefig(filename_start + str(pha) + filename_end, dpi=600, bbox_inches='tight', pad_inches=0.5, format='png')
 





 
