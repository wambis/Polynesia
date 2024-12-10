#!/usr/bin/python
# -*- coding: utf-8 -*-
import obspy
import numpy as np
import os, sys
sys.path.append('/home/wamba/ucbpy3')
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
#from check_2Dpolygons import check_2Dpolygons
from mpl_toolkits.basemap import Basemap
import matplotlib
#matplotlib.use('pdf')
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from obspy.imaging.mopad_wrapper import beach
from operator import truediv
#import seaborn as sns
from math import log10, floor, ceil
from matplotlib.ticker import FormatStrFormatter


#compute the percentiles for 2.5% and 97.5%
#Open a figure for a particular file
#fig, ax = plt.subplots(figsize=(14, 12), edgecolor='w')
fig, ax = plt.subplots(figsize=(14, 16), edgecolor='w')
m       = Basemap(projection='nsper',lat_0 = 9.5, lon_0= 4.2, satellite_height= 2e7, resolution='c', ax = ax)
#region = [8.5, 10.0, 3.5, 5.0]

## Define the Basemap with the Cameroon region
#mc = Basemap(
#    projection='cyl',  # Use 'cyl' projection for simplicity
#    llcrnrlat=2, urcrnrlat=13,  # Latitude bounds for Cameroon
#    llcrnrlon=8, urcrnrlon=16,  # Longitude bounds for Cameroon
#    resolution='i',  # Intermediate resolution
#    ax=ax
#)

m.drawcountries(color='black', linewidth=1.5)  # Country borders
m.drawcoastlines(color='gray', linewidth=1.0)  # Coastlines
####################################################
#mc.drawparallels(range(2, 14, 2), labels=[1, 0, 0, 0], fontsize=10)  # Parallels
#mc.drawmeridians(range(8, 17, 2), labels=[0, 0, 0, 1], fontsize=10)  # Meridians
#Define size list
figname = "Fig_Glob-Cameroon.png"
#get the points of the chunk on the map
#Plot the  topography
#m.etopo(ax=ax, alpha=0.85)
#m.etopo(ax=ax,alpha=0.85)
# Plot topography
m.etopo(alpha=1.0)
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
# Draw features on the map
# Save the figure
