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
#evcatalogue = NDKFile(fname_evcatalogue)
#print('reading event catalogue...')
#print('found {:d} events in the catalogue'.format(evcatalogue.nevents))
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
#F_info    = Fig_params['Finfo']
F_stn     = Fig_params['FstnName']
#open the global station file
f_i = open(F_stn, 'r')
#Get the OBSs file
#Grab the regional station file and events file
#Grab the events file name in order to save the events
#open the global station file
f_station = open("STATIONS_GLADM25_in_CHUNK.dat", 'w')
#Open the OBSs stations
fobs = open("STATIONS_OBSs_in_CHUNK.dat", 'w')
#Open a figure for a particular file
#fig, ax = plt.subplots(figsize=(14, 12), edgecolor='w')
fig, ax = plt.subplots(figsize=(14, 16), edgecolor='w')
m       = Basemap(projection='nsper',lat_0 = -20.0, lon_0= -120., satellite_height= 2e7, resolution='c', ax = ax)
#Define size list
figname = "FigStations_in_Chunk.png"
#Define a box

#get the points of the chunk on the map
xk, yk  = m(chunk[:, 0], chunk[:, 1])
#plot the chunk in the Background of the map, if you plot at the end of the script it will be at the front of the map
#m.plot(xk,yk,c='k',lw=15.0,linestyle='--') # NB Geodetic est bien pour placer les points sur la carte

#m.plot(xk,yk,c='k',lw=5.0,linestyle='-') # NB Geodetic est bien pour placer les points sur la carte

#Grab all the event-station pairs with the corresponding coordinates as described below
################ (network   ,  station)   :  (station_lat, station_lon ,station_elevation depth
INFOS_DICT = {(l.split()[0], l.split()[1]) : (l.split()[2], l.split()[3], l.split()[4], l.split()[5]) for l in f_i}
##################################
for  netwk_station, i  in zip(INFOS_DICT, range(len(INFOS_DICT))):
        #get network and station information
        netwk = netwk_station[0]
        stn   = netwk_station[1]
        #Assign to the variable the initial value
        slat          = float(INFOS_DICT[netwk_station][0])
        slon          = float(INFOS_DICT[netwk_station][1])
        elevation     = float(INFOS_DICT[netwk_station][2])
        depth         = float(INFOS_DICT[netwk_station][3])
        slat_1        = slat 
        slon_1        = slon
        ############## Check ##############################
        if (float(slon) < 0.0)  :   slon  = slon % 360.
        #Check if the event-station contain a bad value
        if check_2Dpolygons(float(slon), float(slat), poly_chunk) == 'IN':
                #######
                ll="{:6s} {:6s} {:10.3f} {:10.3f} {:10.3f} {:10.3f}".format(netwk,stn, slat_1, slon_1 , elevation, depth)
                #Get the station coordinate on the map
                xstn_lon, ystn_lat = m(slon, slat)
                #plot the station on the map
                if((elevation + depth) < 0):
                    #m.plot(xstn_lon, ystn_lat, '^', markersize=30, c='g', markeredgecolor='k', alpha=0.95, axes=ax)
                    m.plot(xstn_lon, ystn_lat, '^', markersize=30, c='g', markeredgecolor='k', alpha=1.0)
                    fobs.write("%s\n"%(ll))
                else:
                    #m.plot(xstn_lon, ystn_lat, 'v', markersize=30, c='y', markeredgecolor='k',alpha=0.95, axes=ax)
                    m.plot(xstn_lon, ystn_lat, 'v', markersize=30, c='y', markeredgecolor='k',alpha=1.0)
                f_station.write("%s\n"%(ll))
############################
#Create a legend
legend_elements = [
    Line2D([0], [0], marker='v', color='w', label='Land-based',
           markerfacecolor='y', markersize=22, lw=12),
    Line2D([0], [0], marker='^', color='w', label='OBS',
           markerfacecolor='g', markersize=22, lw=12)
]

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
m.etopo(ax=ax)
#fig.savefig(figname)
plt.savefig(figname, bbox_inches="tight")
plt.close(fig)
