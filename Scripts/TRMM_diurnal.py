#!/usr/bin/env python2
# -*- coding: utf-8 -*-

### Import all required packages ###
import numpy as np
import matplotlib.pyplot as plt
import os


### Please set working directory here ### -------------------
filepath_BSISO = "/Users/helenleft/Desktop/Data_processing/Data/BSISO/"
filepath_rain = "/Users/helenleft/Desktop/Data_processing/Process_data/"
filepath_plot = "/Users/helenleft/Desktop/Data_processing/Plots/Diurnal_cycle/"

###################
### Choose site ###
###################
# Type site name Dharwad or Berambadi ---------------------
site = "Berambadi"

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

# Load the rainfall data
pcp_re = np.load("TRMM_3HR_" + site + ".npy") 


################################################################
#### Carry out diurnal cycle analysis during JJAS 2016-2018 ####
################################################################

# Initiate new figure
fig = plt.figure(figsize=(8,6))

# Set up axes
ax1 = fig.add_subplot(1,1,1)


### Calculate and plot diurnal cycle for all the days during JJAS 2016-2018 ###
pcp_diur = np.nanmean(pcp_re,axis = 0)
# Change the result to local time
pcp_diur_mean = np.append(pcp_diur[7],pcp_diur[0:7]) 
# Plot diurnal cycle of rainfall for all days during JJAS 2016-2018
ax1.plot(pcp_diur_mean,'k--',label = "Mean" )


### Calculate and plot the diurnal cycle of rainfall during each phase with strong BSISO activity ###
# Change directory to get BSISO index
os.chdir(filepath_BSISO)

# Read in the text file of BSISO, and skip one line of header
BSISO_data = np.genfromtxt("BSISO_25-90bpfil_pc_cdr_JJAS_1618.txt", skip_header=1)
  
# Get the phase and amplitude of BSISO
BSISO_phase = BSISO_data[:,5]
BSISO_amp = BSISO_data[:,6]  

# Diurnal cycle analysis for each BSISO phase 
for pha in range(1,8+1):
    # Locate the data of strong BSISO activity
    BSISO_day = np.where((BSISO_phase == pha)&(BSISO_amp >= 1))
    # Calculate the total number of days with strong BSISO activity
    BSISO_day_sum = np.size(BSISO_day) 
    #print(BSISO_day)

    # Create an array to save the rainfall data during the BSISO phase
    pcp_BSISO = np.zeros((np.size(BSISO_day),8))

    # Loop and put rainfall data into pcp_BSISO array
    n=0
    for n in range(0,np.size(BSISO_day)):
        #print(n)
        rr = BSISO_day[0][n]
        #print("rr = ",rr)
        pcp_BSISO[n,:] = pcp_re[rr,:]

    # Calculate mean diurnal cycle 
    pcp_BSISO_diur = np.nanmean(pcp_BSISO,axis = 0)
    # Change the result to local time
    pcp_BSISO_diur_pha = np.append(pcp_BSISO_diur[7],pcp_BSISO_diur[0:7])

    # Plot the mean diurnal cycle of rainfall for each BSISO phase
    ax1.plot(pcp_BSISO_diur_pha,label = ("Phase " + str(pha) + " (" + str(BSISO_day_sum) + " days)") )

### Carry out diurnal cycle analysis for all days with weak BSISO activity during JJAS 2016-2018 ###
# Locate the data of weak BSISO activity
nBSISO_day = np.where(BSISO_amp < 1)
nBSISO_day_sum = np.size(nBSISO_day)

# Create an array to save the rainfall data of non-active BSISO days
pcp_nBSISO = np.zeros((nBSISO_day_sum,8))

# Put rainfall during non-active BSISO days into pcp_nBSISO
n=0
for n in range(0,np.size(BSISO_day)):
#print(n)
    rr = nBSISO_day[0][n]
    #print("rr = ",rr)
    pcp_nBSISO[n,:] = pcp_re[rr,:]

# Calculate diurnal cycle
pcp_nBSISO_diur = np.nanmean(pcp_nBSISO,axis = 0)
pcp_nBSISO_diur_new = np.append(pcp_nBSISO_diur[7],pcp_nBSISO_diur[0:7])

# Plot diurnal cycle during non-active BSISO days into pcp_nBSISO
ax1.plot(pcp_nBSISO_diur_new,linestyle = '-.',label = ('Weak ' + "(" + str(nBSISO_day_sum) + " days)") )

#######################################################################################
# Label the axes
ax1.set_xlabel("Time", fontsize=12)
ax1.set_ylabel("Rainfall intensity (mm/h)", fontsize=12)

# Set xticks
ax1.set_xticks(np.arange(8))
ax1.set_xticklabels(["2:30","5:30","8:30","11:30","14:30","17:30","20:30","23:30"])

# Set the y-axis view limits
ax1.set_ylim(0,2.5)

# Add a title
ax1.set_title("Diurnal cycle of rainfall at " + site + " during JJAS 2016-2018",fontsize=14)

# Add legend
ax1.legend()

##################
### Save plots ###
##################

#Change working directory
os.chdir(filepath_plot)

# Save plot
filename_start = "Diurnal_rainfall_TRMM_8phas_2016_2018_"
filename_end = ".png"
plt.savefig(filename_start + site + filename_end, dpi=600, format='png')





