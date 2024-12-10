import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Create a map projection
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

# Add background features (e.g., coastlines, country borders)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, edgecolor='black')

# Plot data on the map
# Replace this with your own data or use sample data for demonstration

# Create a legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Legend Label 1',
           markerfacecolor='blue', markersize=10),
    Line2D([0], [0], marker='s', color='w', label='Legend Label 2',
           markerfacecolor='red', markersize=10)
]

# Add the legend to the plot
ax.legend(handles=legend_elements, loc='lower right')

# Set a title for the map
ax.set_title('Map with Legend')

# Show the plot
#plt.show()
#Plot save
plt.savefig("TEST_fig.png", bbox_inches="tight")
