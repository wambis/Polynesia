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







#extension = "png"
extension = "pdf"
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
with open('configEvent.yaml') as Fym:
    #This allow to avoid the _io.TextIOWrapper error 
    Fig_params = yaml.load(Fym, Loader=SafeLoader) 
    #Fig_params = yaml.dump(Fig_params)
#Load the config file
#Grab the path of the data
path      = Fig_params['path']
#get the size
msize     = Fig_params['msize']
##############################
F_info    = Fig_params['Finfo']
#open the global station file
f_i = open(F_info, 'r')
f_i.readline()
#Get the OBSs file
#Grab the regional station file and events file
#Grab the events file name in order to save the events
#Open the file to save the events
#f_event = open(F_names, 'w') 
#open the global station file
f_event = open("EVENTS_GLADM25_in_CHUNK.dat", 'w')
#Open a figure for a particular file
#fig, ax = plt.subplots(figsize=(14, 12), edgecolor='w')
fig, ax = plt.subplots(figsize=(14, 16), edgecolor='w')
m       = Basemap(projection='nsper',lat_0 = -20.0, lon_0= -120., satellite_height= 2e7, resolution='c', ax = ax)
#Define size list
figname = "FigEvents_in_Chunk.png"
#Define a box
#get the points of the chunk on the map
xk, yk  = m(chunk[:, 0], chunk[:, 1])

#get the points of the chunk on the map
xk, yk  = m(chunk[:, 0], chunk[:, 1])
#plot the chunk in the Background of the map, if you plot at the end of the script it will be at the front of the map
#m.plot(xk,yk,c='k',lw=5.0,linestyle='-') # NB Geodetic est bien pour placer les points sur la carte
#Grab all the event-station pairs with the corresponding coordinates as described below
################ event_id       station      event_lat     event_lon      slat          slon       EPIC_DIST
#INFOS_DICT = {(l.split()[0], l.split()[3]): (l.split()[1], l.split()[2], l.split()[4], l.split()[5], l.split()[6]) for l in f_i}
################ event_id       station      event_lat #########################################
INFOS_DICT = {l.split()[0]: (l.split()[1], l.split()[2]) for l in f_i}
##################################
for event_name  in INFOS_DICT:
        #get event information
        #Event information 
        ref_time, event_lat_1, event_lon_1, event_depth, event_CMT = Catalog(event_name)
        #Assign to the variable the initial value
        event_lat = event_lat_1 
        event_lon = event_lon_1
        ############## Check ##############################
        if (float(event_lon) < 0.0)  :   event_lon  = event_lon % 360.
        #Check if the event-station contain a bad value
        if check_2Dpolygons(float(event_lon), float(event_lat), poly_chunk) == 'IN':
                #plot events on the map
                xevent_lon, yevent_lat = m(event_lon, event_lat)
                #plot the event on the map
                #m.plot(x_ev, y_ev, 'o', markersize=msize, c='lime', markeredgecolor='k')
                #Plots All stations on map
                #if(event_depth < 100.0):
                if(event_depth < 300.0):
                    b = beach(event_CMT, xy=(xevent_lon, yevent_lat), width=50, size=150, linewidth=1.5, facecolor='r', alpha=0.95, axes=ax)
                else:
                    b = beach(event_CMT, xy=(xevent_lon, yevent_lat), width=50, size=150, linewidth=1.5, facecolor='y', alpha=0.95, axes=ax)
                b.set_zorder(100)
                ax.add_collection(b)
                #m.drawgreatcircle(event_lon, event_lat, station_lon, station_lat, linewidth=20.0, color='k', alpha=0.6)
                ll="{:16s}  {:8.3f}  {:8.3f}  {:8.3f}".format(event_name, event_lat_1, event_lon_1, event_depth)
                f_event.write("%s\n"%(ll))
############################
#plt.tight_layout()
plot_plates(m, linewidth=5.0, linestyle='-', color='firebrick', alpha=0.5)
#Plot the check
#plt.plot(xk,yk,c='lime',lw=15.0,linestyle='--') # NB Geodetic est bien pour placer les points sur la carte
plt.plot(xk,yk,c='k',lw=15.0,linestyle='--') # NB Geodetic est bien pour placer les points sur la carte
################################################################################
#Plot the  topography
#m.etopo(ax=ax, alpha=0.85)
m.etopo(ax=ax)
#########################
#get the limit of the axis
xlim, ylim = ax.get_xlim()
#Scaling factor
#yscl_factor  = ylim/120.0
yscl_factor  = ylim/100.0
#Set the axis-limit
ax.set_xlim(xlim - yscl_factor, ylim + yscl_factor)
ax.set_ylim(xlim - yscl_factor, ylim + yscl_factor)
#fig.savefig(figname)
#plt.savefig('test_station.pdf')
plt.savefig(figname, bbox_inches="tight")
plt.close(fig)
