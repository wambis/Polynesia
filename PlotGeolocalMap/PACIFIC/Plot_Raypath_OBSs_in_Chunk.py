#!/usr/bin/python
# -*- coding: utf-8 -*-
import obspy
import numpy as np
import os, sys
sys.path.append('/home/mw1685/libs/ucbpy3')
from glob import glob
from ndk_rec_file import NDKFile
from FigTools import plot_plates, plot_hotspots, plot_gc_clean, get_tectonics_path
import json, yaml
from yaml.loader import SafeLoader
import math
import statistics
from numpy import savetxt
from scipy import stats
from scipy.stats import norm
from check_2Dpolygons import check_2Dpolygons
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from obspy.imaging.mopad_wrapper import beach
from operator import truediv
import seaborn as sns
from math import log10, floor, ceil
from matplotlib.ticker import FormatStrFormatter
from matplotlib.lines import Line2D


#compute the percentiles for 2.5% and 97.5%

fname_evcatalogue = "/data/GLOBAL-EVENTS/ndk76-2020.ndk"
def Catalog(event_name):
    evcat       = evcatalogue.findname(event_name)
    #Set Parameter
    ref_time    = evcat.origin + evcat.ctime
    event_lat   = evcat.cposition[0]
    event_lon   = evcat.cposition[1]
    event_depth = evcat.cdepth
    event_CMT   = evcat.focal_mechanism
    M0          = evcat.M0
    return(ref_time, event_lat, event_lon, event_depth, event_CMT)







extension = "png"
#extension = "pdf"
dpi       = 500

#Chunk for station
#chunk       = np.loadtxt('ChunkSEM.90.90.centered.at.20S.110E')
chunk       = np.loadtxt('ChunkSEM.105.105.centered.at.20S.120E')
#got the Polygone from the Chunk points
poly_chunk  = list(zip(chunk[:, 0],chunk[:, 1]))
######################################################
#read the events catalog
evcatalogue = NDKFile(fname_evcatalogue)
print('reading event catalogue...')
print('found {:d} events in the catalogue'.format(evcatalogue.nevents))
#######################################################
#######################################################
#Open the config file
#Get the OBSs file
#grab all the data File
FILE_RAYPATH = "FILE_RAYPATHS_in_Chunk.dat" 
#Grab the regional station file and events file
#Grab the events file name in order to save the events
f_obs = open("STATIONS_OBSs_in_CHUNK.dat", 'r')
#Open a figure for a particular file
#fig, ax = plt.subplots(figsize=(14, 12), edgecolor='w')
#Define size list
figname = "Fig_Rapypath_OBS_in_Chunk_dot.png"
#Define a box
Fevent  = "EVENTS_GLADM25_in_CHUNK.dat"

#get the points of the chunk on the map
fig, ax = plt.subplots(figsize=(14, 16), edgecolor='w')
m       = Basemap(projection='nsper',lat_0 = -20.0, lon_0= -120., satellite_height= 2e7, resolution='c', ax = ax)
#m.plot(xk,yk,c='k',lw=5.0,linestyle='-') # NB Geodetic est bien pour placer les points sur la carte
#get the points of the chunk on the map
xk, yk  = m(chunk[:, 0], chunk[:, 1])
#Grab all the event-station pairs with the corresponding coordinates as described below
################ (network   ,  station)   :  (station_lat, station_lon ,station_elevation depth
INFOS_DICT = {('%s.%s'%(l.split()[0], l.split()[1])) : (l.split()[2], l.split()[3], l.split()[4], l.split()[5]) for l in f_obs}

###############

CHECK = set()
##################################
#get network and station information
#Assign to the variable the initial value
Fp            = open(FILE_RAYPATH, 'r') 
#readline 
Fp.readline() 
for l in Fp:
    event_id  = l.split()[0]
    event_lat = float(l.split()[1])
    event_lon = float(l.split()[2])
    stn_id    = l.split()[3]
    slat      = float(l.split()[5])
    slon      = float(l.split()[6])
############## Check ##############################
    if (float(slon) < 0.0)  :        slon  = slon % 360.
    if (float(event_lon) < 0.0)  :   event_lon  = event_lon % 360.
    #Check if the event-station contain a bad value
    #if check_2Dpolygons(float(event_lon), float(event_lat), poly_chunk) == 'IN':
    if(stn_id in INFOS_DICT.keys()):
        it_chk    = (event_id,  stn_id)
        #######
        if(it_chk not in CHECK):
            CHECK.add(it_chk) 
            #Get the station coordinate on the map
            xstn_lon, ystn_lat     = m(slon, slat)
            xevent_lon, yevent_lat = m(event_lon, event_lat)
            #plot the station on the map
            #m.plot(xstn_lon, ystn_lat, '^', markersize=30, c='g', markeredgecolor='k', alpha=0.95, axes=ax)
            m.plot(xstn_lon, ystn_lat, '^', markersize=30, c='g', markeredgecolor='k', alpha=1.0)
            ###########
            # Draw the great circle route
            m.drawgreatcircle(event_lon, event_lat, slon, slat, linewidth=0.4, color='k', alpha = 0.3)
            #plot the beachball
            ref_time, event_lat_1, event_lon_1, event_depth, event_CMT = Catalog(event_id)
            m.plot(xevent_lon, yevent_lat, 'o', markersize=15, c='k', markeredgecolor='r', alpha=0.7)
#            b = beach(event_CMT, xy=(xevent_lon, yevent_lat), width=50, size=150, linewidth=1.5, facecolor='r', alpha=0.95, axes=ax)
#            b.set_zorder(100)
#            ax.add_collection(b)
############################

print("NPATHS %d " %( len(CHECK) ))
print("******" * 50)
print("NSTATIONS %d " %( len(INFOS_DICT.keys()) ))
#########################
#get the limit of the axis
xlim, ylim = ax.get_xlim()
#Scaling factor
#yscl_factor  = ylim/120.0
yscl_factor  = ylim/100.0
#Set the axis-limit
ax.set_xlim(xlim - yscl_factor, ylim + yscl_factor)
ax.set_ylim(xlim - yscl_factor, ylim + yscl_factor)


#Add the legend to the plot
#ax.legend(handles=legend_elements, loc='upper right')

#ax.legend(handles=legend_elements, fontsize=18, loc='lower right')

#plt.tight_layout()
plot_plates(m, linewidth=5.0, linestyle='-', color='firebrick', alpha=0.5)
#Plot the check
plt.plot(xk,yk,c='k',lw=15.0,linestyle='--') # NB Geodetic est bien pour placer les points sur la carte
################################################################################
#Plot the  topography
#m.etopo(ax=ax)
#fig.savefig(figname)
plt.savefig(figname, bbox_inches="tight")
plt.close(fig)
