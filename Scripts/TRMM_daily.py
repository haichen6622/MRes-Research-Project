#!/usr/bin/env python2
# -*- coding: utf-8 -*-

### Import all required packages ###
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob, os

### Please set working directory here ### ----------------
filepath_BSISO = "/Users/helenleft/Desktop/Data_processing/Data/BSISO/"
filepath_rain = "/Users/helenleft/Desktop/Data_processing/Data/TRMM_3B42_V7_2016_2018_JJAS_3HR/"
filepath_plot = "/Users/helenleft/Desktop/Data_processing/Plots/TRMM_Daily_pha/"
filepath_process = "/Users/helenleft/Desktop/Data_processing/Process_data/"

###################
### Choose site ###
###################
# Type site name Dharwad or Berambadi-------------------
site = "Dharwad"

# Define the latitude and longitude of the sites
if site == "Dharwad":
    lat_site = 15.50 
    lon_site = 74.99
    
elif site == "Berambadi":
    lat_site = 11.76
    lon_site = 76.59
            
else:
    print("No this site")

##############################
### Read in rainfall data  ###
##############################
# Change working directory
os.chdir(filepath_rain)

# Get the rainfall data file names, sorted by name
filenames = sorted(glob.glob('*.nc4'))

# Open first file
nc_fid = Dataset(filenames[0],'r') 
# Have a look what is in file - check the units
print(nc_fid.variables) 

# Get latitude and longitude of the region
lats = nc_fid['nlat'][:]
lons = nc_fid['nlon'][:]
   
# Locate the position of grid point that is nearest to the specified site
# Latitude
a = np.abs(lat_site - lats)
lat_ind = np.where(a == np.min(a))
print(lat_ind)  # print latitude value to check

# Longitude
b = np.abs(lon_site - lons)
lon_ind = np.where(b == np.min(b))
print(lon_ind) # print longitude value to check

# Create empty array to save the rainfall data of a site
pcp_all = []

# Get the time series of rainfall for a specified site
for i in range(0,len(filenames)): 
    # Open file
    data_TRMM = Dataset(filenames[i],'r') 
    
    # Get rainfall
    pcp_TRMM = data_TRMM['precipitation'][:]
    
    # Read in all the times but only the selected lat and lon
    pcp_temp = pcp_TRMM[lon_ind[0],lat_ind[0]][-1]
    
    # Combine rainfall of each year during JJAS to one time series  
    pcp_all.append(pcp_temp)

# Check missing value and replace them by nan 
mis_ind = np.where(pcp_all == -9999.9)[-1]
if not mis_ind:
    print("List is empty")
else:
    pcp_all[mis_ind] = np.nan

# Reshape the rainfall data
pcp_re = np.reshape(pcp_all, (np.size(pcp_all)/8,8)) 

### Save data ### 
# Change working directory
os.chdir(filepath_process)
# Save processed data for diurnal cycle analysis
np.save(filepath_process + "TRMM_3HR_" + site + ".npy", pcp_re)

###########################
### Read in BSISO data  ###
###########################
# Change directory 
os.chdir(filepath_BSISO)

# Read in the text file of BSISO phases, and skip one line of header
BSISO_data = np.genfromtxt("BSISO_25-90bpfil_pc_cdr_JJAS_1618.txt", skip_header=1)
  
# Get the phase and amplitude of BSISO 
BSISO_phase = BSISO_data[:,5]
BSISO_amp = BSISO_data[:,6]  

#######################################################################
### Calculate average daily rainfall during each active BSISO phase ###
#######################################################################

# Convert 3-hour rainfall intensity to daily rainfall
pcp_daily = np.nanmean(pcp_re,axis = 1)*24

# Calculate the average daily rainfall for all the days during JJAS 2016-2018
pcp_ave = np.nanmean(pcp_daily)

# Create empty array to save average daily rainfall of different BSISO phases
pcp_ave_phas = []

# Calculate average daily rainfall of each BSISO phases
for pha in np.arange(1,8+1,1):
    # Locate the date with strong BSISO activity
    BSISO_day = np.where((BSISO_phase == pha)&(BSISO_amp >= 1))
    # Extract daily rainfall during active BSISO days
    pcp_pha = pcp_daily[BSISO_day]
    # Calculate the average daily rainfall during active BSISO days
    pcp_ave_daily = np.nanmean(pcp_pha)
    # Put data into pcp_ave_phas
    pcp_ave_phas.append(pcp_ave_daily)

###########################################################
### Plot average daily rainfall during each BSISO phase ###
###########################################################
   
# Initiate new Figure
fig = plt.figure(figsize=(8,6))

# Set up axes
ax1 = fig.add_subplot(1,1,1)

# Plot data
ax1.bar(np.arange(1,8+1,1),pcp_ave_phas,)
ax1.bar(9,pcp_ave)

# Label the axes
ax1.set_xlabel("Phase", fontsize=12)
ax1.set_ylabel("Rainfall intensity (mm/day)", fontsize=12)

# Set the y-axis view limits
ax1.set_ylim(0,30)

# Set xticks
ax1.set_xticks(np.arange(1,9+1,1))
ax1.set_xticklabels(["1","2","3","4","5","6","7","8","Mean"])

# Add a title
ax1.set_title("Average daily rainfall at " + site + "\n during JJAS 2016-2018 with active BSISO activity",fontsize=14)

# Add legend
ax1.legend()

##################
### Save plots ###
##################

# Change working directory
os.chdir(filepath_plot)

# Save plot
filename_start = "Daily_rainfall_TRMM_8phas_2016_2018_"
filename_end = ".png"
plt.savefig(filename_start + site + filename_end, dpi=600, format='png')



