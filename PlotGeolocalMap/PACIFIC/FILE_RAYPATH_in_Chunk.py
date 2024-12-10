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

##Chunk for station
##chunk       = np.loadtxt('ChunkSEM.90.90.centered.at.20S.110E')
#chunk       = np.loadtxt('ChunkSEM.105.105.centered.at.20S.120E')
##got the Polygone from the Chunk points
#poly_chunk  = list(zip(chunk[:, 0],chunk[:, 1]))
#######################################################
##read the events catalog
#evcatalogue = NDKFile(fname_evcatalogue)
#print('reading event catalogue...')
#print('found {:d} events in the catalogue'.format(evcatalogue.nevents))
#######################################################
#######################################################
#Open the Raypths File
Fr    = open("FILE_RAYPATHS.dat", "w")
#Get the OBSs file
#f_obs = open("STATIONS_OBSs_in_CHUNK.dat", 'r')
f_stn = open("STATIONS_GLADM25_in_CHUNK.dat", 'r')
#grab all the data File
FILES = glob('%s/*.dat' %('ALLDATA')) 
#Grab the regional station file and events file
#Grab the events file name in order to save the events
#f_obs = open("STATIONS_OBSs_in_CHUNK.dat", 'r')
#Open a figure for a particular file
#fig, ax = plt.subplots(figsize=(14, 12), edgecolor='w')
#Define size list
#Define a box
Fevent  = "EVENTS_GLADM25_in_CHUNK.dat"
#get the points of the chunk on the map
#m.plot(xk,yk,c='k',lw=5.0,linestyle='-') # NB Geodetic est bien pour placer les points sur la carte

#Grab all the event-station pairs with the corresponding coordinates as described below
################ (network   ,  station)   :  (station_lat, station_lon ,station_elevation depth
INFOS_DICT = {('%s.%s'%(l.split()[0], l.split()[1])) : (l.split()[2], l.split()[3], l.split()[4], l.split()[5]) for l in f_stn}

###############
Fe              = open(Fevent,'r')
EVENTS_in_Chunk = {l.split()[0] for l in Fe}

####################
CHECK = set()
##################################
for  Fname, i  in zip(FILES, np.arange(len(FILES))):
        #get network and station information
        #Assign to the variable the initial value
        Fp            = open(Fname, 'r') 
        #readline 
        l_1 = Fp.readline() 
        if(i==0):
            Fr.write("%s\n" %(l_1.strip()))
        for l in Fp:
            event_id  = l.split()[0]
            stn_id    = l.split()[3]
        ############## Check ##############################
            if(event_id in EVENTS_in_Chunk):
                if(stn_id in INFOS_DICT):
                    it_chk = (event_id, stn_id)
                    if(it_chk not in CHECK):
                        Fr.write("%s\n" %(l.strip()))
                        CHECK.add(it_chk)
                #Get the station coordinate on the map
############################

