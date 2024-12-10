
#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os, sys, math
import obspy
from obspy import read, read_inventory
import matplotlib.pyplot as plt 
from glob import glob
from scipy import stats
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from math import log10, floor, ceil
import numpy as np
import pandas as pd
from tabulate import tabulate
#for interpolation
from scipy.interpolate import InterpolatedUnivariateSpline
#######################################
from matplotlib.ticker import ScalarFormatter
from matplotlib.transforms import Affine2D
#######################
import xlrd
#################################################
from matplotlib.ticker import FormatStrFormatter
##############################################
import matplotlib.dates as mdates
#######################################
from matplotlib.colors import Normalize
import matplotlib.cm as cm
############################
import datetime as dtt 
import xarray as xr
import pandas as pd
from mhkit import dolfyn
from mhkit import dolfyn as dlfn
#################
from matplotlib.dates import DateFormatter
from matplotlib.dates import HourLocator
##############
import matplotlib.dates as mdates
from mhkit.dolfyn.adp import api
from pandas.core.common import flatten
##################################
from obspy.core import UTCDateTime

AphaDict = {0:'(a)', 1: '(b)', 2: '(c)', 3: '(d)', 4: '(e)',
            5: '(f)', 6: '(g)', 7:'(h)', 8:'(i)', 9: '(j)',10:'(k)',
            11:'(l)', 12: '(m)',13: '(n)',14:'(o)', 15:'(p)', 16: '(q)'}


import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
from cartopy import crs as ccrs
from cartopy import feature as cfeature

#File path to the GeoTIFF file
#tif_file = "output_GEBCOSubIceTopo.tif"
tif_file = "output_SRTM15Plus.tif"

#Open the GeoTIFF file
with rasterio.open(tif_file) as src:
    # Read the data
    topo_data = src.read(1)  # Reading the first band
    # Get geographic extent of the data
    extent = src.bounds

#Plot the data using Rasterio's built-in plot
plt.figure(figsize=(12, 8))
show(topo_data, cmap="terrain", title="Topography (GEBCO SubIceTopo)")

#Optionally plot with Cartopy
plt.figure(figsize=(12, 8))
ax = plt.axes(projection=ccrs.PlateCarree())

#Set map extent using the bounds of the raster
ax.set_extent([extent.left, extent.right, extent.bottom, extent.top], crs=ccrs.PlateCarree())

#Plot the topography data
plt.imshow(
    topo_data,
    extent=[extent.left, extent.right, extent.bottom, extent.top],
    transform=ccrs.PlateCarree(),
    cmap="terrain",
    origin="upper"
)

# Add map features
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=":")
ax.add_feature(cfeature.LAND, edgecolor="black", alpha=0.3)

# Add labels and title
plt.title("Topography with Cartopy (GEBCO SubIceTopo)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar(label="Elevation (m)")

# Show the plot
# Show the map
figname   = "Map_MtCameroon.png"
#Set the font size of yticks
plt.yticks(fontsize=13)
# Set the font size of xticks
plt.xticks(fontsize=14)
##Space between the subplots
#Aligned all the ylabels of the figure
#fig.align_ylabels()
#Save the figure
plt.savefig(figname, bbox_inches = 'tight', dpi = 300)
#fig.savefig(figname, bbox_inches = 'tight', dpi = 510)
#plt.savefig(figname)










##Grab the space between the figure
#fig_space = Fig_params['fig_space']
##Grab the figure size
#fig_size  = (Fig_params['fig_width'], Fig_params['fig_height'])
#
#
#if(pm3D_adcp_up == False or pm3D_adcp_down == False):
#    #Create  a normal 2D figure
#    fig, axs  = plt.subplots(nfigs, 1, sharex = False, figsize = fig_size)
#
#
#if(velocity_profile_up):
##Create the figure
# #Grab the figure size
# fig, axs  = plt.subplots(1, nfigs, sharex = False, figsize = fig_size)


#figname   = "Map_MtCameroon.png"
##Set the font size of yticks
#plt.yticks(fontsize=13)
## Set the font size of xticks
#plt.xticks(fontsize=14)
###Space between the subplots
##Aligned all the ylabels of the figure
#fig.align_ylabels()
##Save the figure
#fig.savefig(figname, bbox_inches = 'tight', dpi = 300)
##fig.savefig(figname, bbox_inches = 'tight', dpi = 510)
#fig.savefig(figname)
