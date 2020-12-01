#!/usr/bin/env python2
# -*- coding: utf-8 -*-

### Import all required packages ###
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import os
from scipy.optimize import curve_fit


### Choose the site: Dharwad or Berambadi ### ------------
site  = "Berambadi"

### Please set working directory here ### ---------------
filepath_FluxTower = "/Users/helenleft/Desktop/"+ site + "/"
filepath_rain = "/Users/helenleft/Desktop/Data_processing/Data/IMERG/"
filepath_plot = "/Users/helenleft/Desktop/Data_processing/Plots/Drydowns/"
    
##########################################################
### Read in IMERG rainfall estimates at specified site ###
##########################################################
# Change working directory to load in rainfall data and soil moisture
os.chdir(filepath_rain)

# Read data
data_IMERG = np.load('GPM_IMERG_WesternGhats_2016.npz')
    
# Get rainfall
pcp_IMERG = data_IMERG['precip']
    
# Get latitude, longitude and time
lats = data_IMERG['lat']
lons = data_IMERG['lon']
time_dt_fmt = data_IMERG['time_dt_fmt']

# Define the location of the sites
if site == "Dharwad":
    lat_site = 15.50 
    lon_site = 74.99
    
elif site == "Berambadi":
    lat_site = 11.76
    lon_site = 76.59
            
else:
    print("No this site")
       
# Locate the grid point that is nearest to specified site
# Latitude
a = np.abs(lat_site - lats)
lat_ind = np.where(a == np.min(a))
print(lat_ind)  # print latitude value to check

# Longitude
b = np.abs(lon_site - lons)
lon_ind = np.where(b == np.min(b))
print(lon_ind) # print longitude value to check
    
# Read in all the times but only the selected lats and lons
pcp_30 = pcp_IMERG[:,lat_ind[0],lon_ind[0]]

# Converty 30-min rainfall to daily rainfall
pcp_30_re = np.reshape(pcp_30,(122,48))
pcp_daily = np.sum(pcp_30_re,axis = 1)/2

#######################################
### Read in soil moisture ###
#######################################
# Change working directory to get FluxTower data
os.chdir(filepath_FluxTower)

# Read in flux tower data
flt = np.genfromtxt(site + "_GapFilled.txt", skip_header=2)
len(flt)

# Extract soil moisture content
vwc_1 = flt[:,37]
vwc_2 = flt[:,39]

# Check missing value, replace -9999. with nan
miss_v = np.where(vwc_1==-9999)[-1]
vwc_1[miss_v] = np.nan

miss_v = np.where(vwc_2==-9999)[-1]
vwc_2[miss_v] = np.nan

# Create time series 
base_time = dt.datetime(2016,2,12,0,30)
timestamp = np.array([base_time + dt.timedelta(minutes = x) for x in np.arange(0, len(flt)*30,30)])
#print(timestamp[0]) # print to check

# Convert India local time to UTC 00  to match the TRMM data and BSISO index
timestamp_utc = timestamp -  dt.timedelta(minutes = 330)

# Extract data for 2016 JJAS 
# Find the index of start time and end time
start = np.where(timestamp_utc == dt.datetime(2016,6,1,0,00))[-1]
end = np.where(timestamp_utc == dt.datetime(2016,9,30,23,30))[-1]

# Extract time
t_16_JJAS = timestamp_utc[start[0]:end[0]+1]

# Extract soil moisture content 
vwc_1_JJAS = vwc_1[start[0]:end[0]+1]
vwc_2_JJAS = vwc_2[start[0]:end[0]+1]


# Convert 30-min soil moisture to daily average
vwc_1_JJAS_re = np.reshape(vwc_1_JJAS,(122,48))
vwc_1_JJAS_d = np.mean(vwc_1_JJAS_re,axis = 1)

vwc_2_JJAS_re = np.reshape(vwc_2_JJAS,(122,48))
vwc_2_JJAS_d = np.mean(vwc_2_JJAS_re,axis = 1)


##################################
### Filter soil drydown events ###
##################################

# Find out the dates of daily rainfall less than 5 mm
pcp_dry_ind = np.where(pcp_daily[0:122] < 5)[-1]
print(pcp_dry_ind)

# Create an empty list to save the index of drydown events 
dry_ind = list()

# Loop to screen out dry-down events
j=1
for i in pcp_dry_ind:

    if (i+1 in pcp_dry_ind)  == True:
        j = j+1
    
    else:
        if j>= 3:
            print(np.arange(i-j,i+1,1)) # print(np.arange(i-j+1,i+1,1))
            ind_temp = np.arange(i-j,i+1,1)
            
            dry_ind.append(ind_temp)
            
            j = 1
        else:
            j = 1


### Manually check and adjust the periods of drydown events ### 
# Eliminate the period during which soil moisture increase (e.g.due to irrigation) 
dry_ind_copy = dry_ind

if site == "Dharwad":
# For Dharwad soil moisture 1 -------------
    del dry_ind[6]
    dry_ind[5] = np.array([108,109,110,111])
    dry_ind[3] = np.array([84, 85, 86, 87])
    del dry_ind[2]

elif site == "Berambadi":
    # For Berambadi soil moisture 1 -------------
    dry_ind[1] = np.array([8,9,10,11,12,13])
    dry_ind[5] = np.array([57,58,59,60,61,62,63,64,65,66])
    dry_ind[6] = np.array([72,73,74,75])
    del dry_ind[7::]


#############################################
### Define the function for model fitting ###
#############################################
def func(x,a,b,c):
    return a*np.exp((-1/b) * x) + c


####################################
### Exponential fitting and plot ###
####################################

# Create number of days for x-axes label
x_date = np.arange(0,122,1)

# Create an array to save parameters from exponential fitting
para_ary = np.zeros((len(dry_ind),4))


# Initiate new figure
fig = plt.figure(figsize=(7,4))

# Set up axes
ax1 = fig.add_subplot(1,1,1)

# Carry out exponential fitting for each drydown event and plot the result 
for ii in range(0,len(dry_ind)):

    ### Check and delete the missing data ###
    y_sm = vwc_1_JJAS_d[dry_ind[ii]]
    ind_nan = np.where(np.isnan(y_sm))[-1]
    y_new = np.delete(y_sm,ind_nan) 

    x_d = np.arange(1,len(y_sm)+1,1)
    x_new = np.delete(x_d,ind_nan)
    
    # Exponential fitting 
    popt,pcov = curve_fit(func,x_new,y_new,method = 'trf',bounds = ((0,0,np.nanmin(vwc_1)),(np.inf,np.inf,np.min(y_new))))
    
    ### Calculate R2 ###
    # Residual sum of squares
    ss_res = np.sum((y_new - func(x_new, *popt)) ** 2)
    
    # Total sum of squares
    ss_tot = np.sum((y_new - np.mean(y_new)) ** 2)

    # r-squared
    r2 = 1 - (ss_res / ss_tot)
        
    # Save parameters
    para_ary[ii,0:3] = popt
    para_ary[ii,3] = r2
    
    # Get the date of x-axes for plotting
    x_axes = x_date[dry_ind[ii]]
    x_axes = np.delete(x_axes,ind_nan)
    
    # Plot the fitting results 
    ax1.plot(x_axes, func(x_new, *popt), ls = '--', color = 'dimgray', label= r'$fit: A=%5.2f, \tau=%5.2f, \theta_f=%5.2f, r^2=%5.2f $' % tuple(para_ary[ii,:])) # tuple(popt)

# Create dates in x-axes for in-situ soil moisture plot during dry-down events
x_date_sm = []
for k in range(0,len(dry_ind)):
    x_date_sm = np.concatenate((x_date_sm, dry_ind[k]), axis=None)  

# Plot in situ soil moisture
ax1.scatter(x_date_sm, vwc_1_JJAS_d[x_date_sm.astype(np.int64)], marker='x', color = 'orangered', label = 'In situ soil moisture')

# Add second y-axes for rainfall
ax2 = ax1.twinx()
# Plot daily rainfall
ax2.bar(x_date,pcp_daily[0:122],color = 'lightseagreen',label = "Rainfall")

# Add legend
ax1.legend(loc="best", bbox_to_anchor=(0.68, -0.15))
ax2.legend(loc="best", bbox_to_anchor=(1, -0.15))

# Set the y-axis view limits
#ax1.set_ylim(0,50)

# Set ticks and labels of x-axes
ax1.set_xticks([0,14,29,43,60,74,91,105,121])
ax1.set_xticklabels(["Jun-1","Jun-15","Jul-1","Jul-15","Aug-1","Aug-15","Sep-1","Sep-15","Sep-30"])

# Set label of x-axes and y-axes
ax1.set_xlabel("Date")
ax1.set_ylabel("Soil Moisture Content (%)")
ax2.set_ylabel("Rainfall (mm/day)")

# Add a title
ax1.set_title("Exponential fitting of drydowns at " + site + " during JJAS 2016",fontsize=14)

#plt.show()

#Change working directory
os.chdir(filepath_plot)

# Save plot
filename_start = "Drydown_JJAS_2016_"
filename_end = ".png"
plt.savefig(filename_start + site + filename_end, dpi=600, bbox_inches='tight', pad_inches=0.5, format='png')