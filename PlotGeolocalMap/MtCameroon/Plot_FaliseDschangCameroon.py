
#Plot the mseed seismic data files
import json, yaml
from yaml.loader import SafeLoader
import os, sys, math
import matplotlib.pyplot as plt 
from glob import glob
from scipy import stats
from scipy.stats import norm
import matplotlib.ticker as mticker
import numpy as np
import pygmt
import tempfile 
################################


# Set the GMT session to modern mode 
os.environ["GMT_SESSION_NAME"] = "modern"

# Define the region around Falaise de Dschang (approximate coordinates)
#region = [9.85, 10.05, 5.25, 5.45]  # [min_lon, max_lon, min_lat, max_lat]
region = [9.5, 11.0, 4.8, 6.3]  # [min_lon, max_lon, min_lat, max_lat]
# Create a PyGMT figure
fig = pygmt.Figure()

# Add the map, using a topographic relief dataset
fig.basemap(region=region, 
            projection="M6i", 
            frame=["a0.5", "WSne"])
######################################
#fig.grdimage("@earth_relief_01m", 
fig.grdimage("@earth_relief_15s",  #good for regional
             region=region, 
             shading=True, 
             #cmap="viridis")
             cmap="globe")

# Add a title

fig.text(x=10.2, y=5.5, text="Falaise of Dschang", font="18p,Helvetica-Bold,black")
#Add coastlines, rivers, and borders
fig.coast(region=region, resolution="10m", borders=[1, 2], rivers="1/0.25p,blue", shorelines="1/0.5p,black")

##Generate random coordinates for 50 seismic stations 
#np.random.seed(42)             #for reproducibility 
#lons = np.random.uniform(8.5, 10.0, 50) 
#lats = np.random.uniform(3.5, 5.0, 50)
#
## Plot the seismic stations 
#fig.plot(x=lons, y=lats, style="i0.5c", fill="red", pen="black", label="Seismic Stations")

# Add a legend 
fig.legend(position="JTR+jTR+o0.3c", box=True)

# Add a colorbar 
#fig.colorbar(frame=["a2000", "+lElevation (m)"])
fig.colorbar(position="JMR+o0.5c/0c+w12c/0.5c", frame=["a2000", "+lElevation (m)"])
#Show the plot
figname   = "Map_FaliseDscahngCameroon.png"
#Set the font size of yticks
#plt.yticks(fontsize=13)
## Set the font size of xticks
#Save the figure
#fig.savefig(figname, bbox_inches = 'tight', dpi = 300)
#fig.savefig(figname, bbox_inches = 'tight', dpi = 510)
fig.savefig(figname)
