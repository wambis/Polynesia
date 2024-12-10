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

#Define the region around Mount Cameroon (longitude and latitude bounds)
region = [8.5, 10.0, 3.5, 5.0]

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
fig.text(x=9.2, y=4.5, text="Mount Cameroon", font="18p,Helvetica-Bold,black")

#Add coastlines, rivers, and borders
fig.coast(region=region, resolution="10m", borders=[1, 2], rivers="1/0.25p,blue", shorelines="1/0.5p,black")

#Generate random coordinates for 50 seismic stations 
np.random.seed(42)             #for reproducibility 
lons = np.random.uniform(8.5, 10.0, 50) 
lats = np.random.uniform(3.5, 5.0, 50)

# Plot the seismic stations 
fig.plot(x=lons, y=lats, style="i0.5c", fill="red", pen="black", label="Seismic Stations")

# Add a legend 
fig.legend(position="JTR+jTR+o0.3c", box=True)

#fig.colorbar(frame=["a250", "x+lElevation", "y+lm"])
# Add a colorbar 
#fig.colorbar(frame=["a2000", "+lElevation (m)"])
fig.colorbar(position="JMR+o0.5c/0c+w12c/0.5c", frame=["a2000", "+lElevation (m)"])
#Show the plot




# Add grid lines
#fig.grid(region=region, spacing="0.5", frame=True)

# Show the plot


# Show the plot
# Show the map
figname   = "Map_MtCameroon.png"
#Set the font size of yticks
#plt.yticks(fontsize=13)
## Set the font size of xticks
#plt.xticks(fontsize=14)
##Space between the subplots
#Aligned all the ylabels of the figure
#fig.align_ylabels()
#Save the figure
#fig.savefig(figname, bbox_inches = 'tight', dpi = 300)
#fig.savefig(figname, bbox_inches = 'tight', dpi = 510)
fig.savefig(figname)
