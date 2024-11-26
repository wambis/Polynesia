#!/usr/bin/python
# -*- coding: utf-8 -*-
import obspy
import numpy as np
import os, sys
sys.path.append('/home/mw1685/libs/ucbpy3')
from glob import glob
from ndk_rec_file import NDKFile
from FigTools import plot_plates, plot_hotspots, plot_gc_clean, get_tectonics_path
from glob import glob
from check_2Dpolygons import check_2Dpolygons
from mpl_toolkits.basemap import Basemap
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.mplot3d import Axes3D
import json, yaml
from yaml.loader import SafeLoader
import math
import statistics
from numpy import savetxt
from scipy import stats
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from operator import truediv
import seaborn as sns
from math import log10, floor, ceil
from matplotlib.ticker import FormatStrFormatter
from matplotlib.font_manager import FontProperties
from matplotlib import rc
#Enable LaTeX rendering in matplotlib 
#rc('text', usetex=True)
#from termcolor import colored


def EXTENSION(string):
    #split by the point
    Z_, ext_file = string.split(".") 
    return(ext_file)

def CheckCategory(path_file):
    Dirbasename          = os.path.basename(path_file) 
    event_name,ext_file  = Dirbasename.split(".")
    #Open a dictionary to keep the measurements 
    if(ext_file == "17_40"):
        tmp        = "17_40"
    elif(ext_file=="40_100#body_wave"):
        tmp        = "40_100#body_wave"
    elif(ext_file == "40_100#surface_wave"): 
        tmp        = "40_100#surface_wave"
    elif(ext_file == "90_250"): 
        tmp        = "90_250"
    return (event_name, tmp)



def BP_Type_func(ctype):
    #Open a dictionary to keep the measurements
    if ctype   == "17_40":
        bp         = "17-40 s"
        tmp        = "body_wave"
    elif ctype == "40_100#body_wave" :
        bp         = "40-100 s"
        tmp        = "body_wave"
    elif ctype == "40_100#surface_wave":
        bp         = "40-100 s"
        tmp        = "surface_wave"
    elif ctype == "90_250" :
        bp         = "90-250 s"
        tmp        = "surface_wave"
    return (tmp, bp)


def readwindow(json_file):
    if os.path.isfile(json_file):
        try:
            fp = open(json_file)
            windows    = json.load(fp)
        except ValueError:
            print("Check this window file: %s "%(json_file))
            sys.exit(1)
    return windows


def check_region_event(event_name):
    if event_name in REG_EVENTS:
        return True
    else:
        return False

###################################################
#compute the percentiles for 2.5% and 97.5%
#percs      = np.asarray([2.5, 97.5])
def Percentiles(data, percentiles):
    #print(data)
    #sys.exit()
    #calculate the values of percentiles
    P_L, P_H   =  np.percentile(data, percentiles, axis=0)
    #Trim the data according to the percentiles
    data_trim  = np.asarray([it for it in data if (it > P_L and it < P_H)]) 
    #Compute the max and the min of the percentage
    dmax       = max(data_trim)
    dmin       = min(data_trim)
    #return(P_L, P_H, dmin, dmax, selx)
    return(P_L, P_H, data_trim)
####################################################################
def AxisWrite1(ax, x, y, data, perc, fsize):
    #fsize is fontesize
    #perc, is the list of the percentiles
    ax.text(x, y/1.15, '%s-%s'%(perc[0], perc[1])+'%', horizontalalignment='center',verticalalignment='center',alpha=1.0, fontsize=fsize)
    ax.text(x, y/1.35, 'N= %d'%(data.size), horizontalalignment='center',verticalalignment='center',alpha=1.0, fontsize=fsize)
    ax.text(x, y/1.55, r'$\mu_{G} =$'+'%.1f'%(np.mean(data)), horizontalalignment='center',verticalalignment='center',alpha=1.0, fontsize=fsize)
    ax.text(x, y/1.85, r'$\sigma_{G} =$'+'%.1f'%(np.std(data)), horizontalalignment='center',verticalalignment='center',alpha=1.0, fontsize=fsize)
    ax.text(x, y/2.3,  'max = %.1f'%(max(data)), horizontalalignment='center',verticalalignment='center',alpha=1.0, fontsize=fsize)
    ax.text(x, y/2.9,  'min = %.1f'%(min(data)), horizontalalignment='center',verticalalignment='center',alpha=1.0, fontsize=fsize)
####################################################################
def AxisWrite2(ax, x, y, data, perc, fsize):
    ax.text(x, y/2.5,   '%s-%s'%(perc[0], perc[1])+'%', horizontalalignment='center',verticalalignment='center',alpha=1.0,fontsize=fsize)
    ax.text(x, y/1.85,  'N= %d' % (data.size), horizontalalignment='center',verticalalignment='center',alpha=1.0,  fontsize=fsize)
    ax.text(x, y/1.57,  r'$\mu_{R}=$'+'%.1f'%(np.mean(data)), horizontalalignment='center',verticalalignment='center',alpha=1.0,fontsize=fsize)
    ax.text(x, y/1.35,  r'$\sigma_{R}=$'+'%.1f'%(np.std(data)), horizontalalignment='center',verticalalignment='center',alpha=1.0,fontsize=fsize)
    ax.text(x, y/1.21,  'max = %.1f'%(max(data)), horizontalalignment='center',verticalalignment='center',alpha=1.0, fontsize=fsize)
    ax.text(x, y/1.1,   'min = %.1f'%(min(data)), horizontalalignment='center',verticalalignment='center',alpha=1.0, fontsize=fsize)
####################################################################
def parameters(data):  
    #compute the mean and the std
    mu          = np.mean(data)
    sigma       = np.std(data)
    ##################################
    #Compute the median
    med         = statistics.median(data)
    #Compute the counts and the bin for global and regional data
    counts, bins  = np.histogram(data)
    #compute the median of the Data
    #Sturge's rule for bins avaluations
    nbins   = int(1 + 3.322 * np.log10(data.size))
    ## OR #############################################################
    #nbins  = int(1+ np.log2(data_H.size)) 
    ###########################################################
    heights, edges = np.histogram(data,  bins=int(nbins))
    ################################################
    #Compute the bandwith #Set a width for each bar
    bwidth  = (max(data) - min(data))/nbins
    return(mu,sigma, med, nbins, bwidth)
####################################################################
def count_bin(data, nbins):
    fig_ , ax        = plt.subplots(1,1) 
    #N is the count in each bin, bins is the lower-limit of the bin
    N, bins, patches = ax.hist(data, bins = nbins)
    N_norm_max       = max(N/sum(N)) 
    #N_Limit          = (2.0  *  N_norm_max) + N_norm_max * (8.0/100.0) #good for binwidth = 0.5
    #N_Limit          =  N_norm_max + N_norm_max * (54.0/100.0)
    N_Limit          = 6.0 * N_norm_max + N_norm_max * (70.0/100.0) 
    plt.close(fig_)
    return(N, N_Limit, bins, patches)
####################################################################
def SYMETRIC_BIN(data):
    #compute the median of the Data
    #Sturge's rule for bins avaluations
    nbins   = int(1 + 3.322 * np.log10(data.size))
    ## OR #############################################################
    #nbins  = int(1+ np.log2(data.size))
    #Make the edge of the bin centered at -1/+1
    #nbins_symetric = np.arange(-nbins, nbins +1, 1)
    #Make the bins symetric At zero, so the central bin is centered at zero
    #nbins_symetric = np.arange(-nbins + 0.5, nbins + 0.5, 1)
    nbins_symetric = np.arange(-nbins + 0.5, nbins + 0.5, 0.185)
    return(nbins_symetric)
####################################################################
#def Plot_hist(ax, data, X, Pdf, Ymax, bins_sym, color,  label):
def Plot_hist(ax, data, X, Pdf, bins_sym, color,  label):
    #Args, ax-axis, data, X-coordinate of PDF, PDF, bins_sym (symetric), color, label
    #ax.hist(data, density=True, color = color, ec='k',lw=2.0, alpha=1.0, rwidth = 0.75, 
    N, bins, patches = ax.hist(data, density=True, color = 'grey', ec='k',lw=2.0, alpha=1.0, rwidth = 1.0,
                                    bins=bins_sym, histtype='bar')
    #Plot the histogram n, bins, patches = ax.hist(data, bins=30, edgecolor='black')
    for patch, left, right in zip(patches, bins[:-1], bins[1:]):
        #if left >= 0.4:
        if left >= 0.5:
            patch.set_facecolor('blue')
        #elif(right <= -0.4):
        elif(right <= -0.5):
            patch.set_facecolor('red')
    #Plot the Probability density function
    ax.plot(X, Pdf, color = color,linestyle='-', lw=5.0, alpha=1.0)
    ax.plot(X, Pdf, color = 'k',linestyle='-', lw=2.0, alpha=.9)
    #Scatter  
    #ax.scatter(0.0, Ymax * 0.97, marker ="o", edgecolor ="green",linewidths = 2, s = 50)
    #ax.legend([label])

def Scatter_G(ax, y, data,color) :
    #yp = y * 0.8
    #yp = y * 0.9
    yp  = 3.15
    mean = np.mean(data) 
    std  = np.std(data) 
    X1 = mean + std
    X2 = mean - std
    X_mid = (X1 + X2)/2.0
#    X = np.linspace(VL, VR,10) 
#    YY = np.asarray([yp,yp,yp,yp,yp,yp,yp,yp,yp,yp])
#    #plotting line within the given range
#    ax.scatter(float('%.1f'%(mean)), yp, marker ="o", color=color, edgecolor =color,linewidths = 3, s = 70)
#    ax.scatter(float('%.1f'%(mean)), yp, marker ="o", color=color, edgecolor = 'k',linewidths = 3, s = 72, alpha=.6)
    #############################################
    ax.scatter(X_mid, yp, marker ="o", color=color, edgecolor =color,linewidths = 3, s = 70)
    ax.scatter(X_mid, yp, marker ="o", color=color, edgecolor = 'k',linewidths = 3, s = 72, alpha=.6)
    #ax.plot(X, YY, color =color, lw=1.5)
    ax.plot([X1, X2], [yp, yp], color =color, lw=1.5)

def Scatter_R(ax, y, data, color) :
    #y= y * 0.8
    y= y * 0.82
    #y= y * 0.9
    #y= 3.15
    mean = np.mean(data) 
    std  = np.std(data) 
    X1 = mean + std
    X2 = mean - std
    X_mid = (X1 + X2)/2.0
#    X  = np.linspace(VL, VR,10) 
#    YY = np.asarray([y,y,y,y,y,y,y,y,y,y])
#    # plotting line within the given range
#    ax.scatter(float('%.1f'%(mean)), y, marker ="o", edgecolor =color, color=color, linewidths = 3, s = 70)
#    ax.scatter(float('%.1f'%(mean)), y, marker ="o", edgecolor ='k', color=color, linewidths = 3, s = 72,alpha=.6)
    #ax.scatter(0.0, y, marker ="o", edgecolor =color, color=color, linewidths = 3, s = 70)
    #ax.scatter(0.0, y, marker ="o", edgecolor ='k', color=color, linewidths = 3, s = 72,alpha=.6)
    #ax.scatter(X, YY, marker ="o", edgecolor ="green",linewidths =0.005)
    ############################################################
    ax.scatter(X_mid, y, marker ="o", color=color, edgecolor =color,linewidths = 3, s = 70)
    ax.scatter(X_mid, y, marker ="o", color=color, edgecolor = 'k',linewidths = 3, s = 72, alpha=.6)
    #ax.plot(X, YY , color =color, lw=1.5)
    ax.plot([X1, X2], [y, y], color =color, lw=1.5)

#Set the density to True if you want to normalize the Y-axis between [0,1]
density = True
################################
#Open the config file
#Fym = open("configG.yaml","r") 
with open('configs.yaml') as Fym:
    #This allow to avoid the _io.TextIOWrapper error 
    Fig_params = yaml.load(Fym, Loader=SafeLoader) 
    #Fig_params = yaml.dump(Fig_params)
#Load the config file
#Grab the path of the data
path      = Fig_params['path']
#Grab the regional station file and events file
fname_stn = Fig_params['FstnName']
fname_evt = Fig_params['FevtName']
#t1, t2 for text writen
t1         = Fig_params['t1']
t2         = Fig_params['t2']
############# Get the X-axis limit fixed by the user ##############
Xmin_fix       = Fig_params['Xmin']
Xmax_fix       = Fig_params['Xmax']
#get the fonte for text
fontsize_Title = Fig_params['fontsize_Title']
fontsize_Label = Fig_params['fontsize_Label']
fontsize_Txt   = Fig_params['fontsize_Txt']
fsize_abc      = Fig_params['fsize_abc']
#Grab the original percentiles
percs_O = np.asarray(Fig_params['percs_O'])
#Grab trim percentiles
percs   = np.asarray(Fig_params['percs'])
#Grab all the modes used for FWI
#get the type of the path to plot
path_average = Fig_params['path_average']
path_absolue = Fig_params['path_absolue']
#modes= ['17_40','40_100#body_wave','40_100#surface_wave','90_250']
modes = Fig_params['modes']
#############
YY_max = Fig_params['Ymax']
####
WaveTpes = Fig_params['WaveTpes']
#Create a figure label according to the different modes, and givent this possible list of labels
FigLab = {'17_40':['a','b','c'], '40_100_body_wave':['d','e','f'], '40_100_surface_wave':['g','h','i'], '90_250':['j','k','l']}
#X-axis visbility dependent on the mode
AxiVis = {'17_40':{'Ax': False, 'rc':'aquamarine', 'gc': 'deeppink'},   '40_100_body_wave':{'Ax': True,  'rc':  'aquamarine', 'gc':'deeppink'},
        '40_100_surface_wave':{'Ax':False, 'rc':'palegreen', 'gc': 'r'}, '90_250':{'Ax':True,      'rc': 'palegreen',   'gc': 'r'}}

#Set the locations of X-axis
locs = np.arange(Xmin_fix, Xmax_fix +1, 1)
#print(FigLab) 
#sys.exit()

#for m in modes:

#Regional station file
fs = open(fname_stn)
#REG_STATIONS = set(['%s.%s'%(l.split()[0], l.split()[1]) for l in fs])
################################### station name          slat_lat       stat_lon
REG_STATIONS = {'%s.%s'%(l.split()[0], l.split()[1]) : (l.split()[2], l.split()[3])  for l in fs}
##File for regional events
f    = open(fname_evt)
##grab regional events: event-name   event lat      event lon
REG_EVENTS = {l.split()[0].strip() : (l.split()[1], l.split()[2])  for l in f}
#Grab all the event used in FWI

#Chunk for station
msize = 15
chunk       = np.loadtxt('REGIONAL_EVENTS_and_STATIONS/ChunkSEM.105.105.centered.at.20S.120E')
#got the Polygone from the Chunk points
poly_chunk  = list(zip(chunk[:, 0],chunk[:, 1]))
#for mode in modes:
for File in WaveTpes:
    #Grab all the events for this mode (i.e: 40_100#body_wave) used in FWI
    mode   = File.split("_in_")[0]

    #events  = sorted([it for it in glob("%s/*"%(path)) if (os.path.isdir(it) and EXTENSION(it) == mode)])
    fp      = open(File, "r")
    fp.readline()
    #read all the lines
    lines   = fp.readlines()
    ########## Set Dictionary for Global Events ###################################
    #DicGlob     = {mode: {'Z': [], 'R': [], 'T': []}}
    # event_name     event_lat  event_lon  station     station_lat  station_lon  DT/T
    DicReg={(l.split()[0], l.split()[3]) : {'event_coord':(float(l.split()[1]), float(l.split()[2])), 'stn_coord':(float(l.split()[4]), float(l.split()[5]))} for l in lines}
    #print(DicReg)
    ##Event-station pair ## event_name      station
    DicRegPath   = {(l.split()[0], l.split()[3]): []  for l in lines}
    #category 
    #####################################
    for ll in lines:
        #Get the Name and the window type of the event
        event_name, event_lat, event_lon, stn_name,  slat,  slon,  DT = ll.split()
        #################################################
        path_i   = (event_name, stn_name)
        #Loop over the stations
        DicRegPath[path_i].append(float(DT))
    #####Plot the map and the histogram of each frequency band
    #Create a figure
    #fig, (ax1, ax2) = plt.subplots(2, 1,  figsize=(14, 14), edgecolor='w')
    fig, (ax1, ax2) = plt.subplots(2, 1,  figsize=(14, 16), edgecolor='w')
    #Set the basemap
    m       = Basemap(projection='nsper',lat_0 = -20.0, lon_0= -120., satellite_height= 2e7, resolution='c', ax = ax1)
    ##get the points of the chunk on the map
    xk, yk  = m(chunk[:, 0], chunk[:, 1])
    #set the figure name
    #figname = "Maps_Histograms_%s.png"%(File.replace(".dat", ""))
    figname = "Maps_Histograms_%s.pdf"%(File.replace(".dat", ""))
    #figname = "Maps_Histograms_%s_sorted.png"%(File.replace(".dat", ""))
    #open a list to append the data
    data_RF   = []
    paths_red = 0; paths_blue = 0; paths_gray = 0
    #Loop over the Dictionary
    for item in DicRegPath:
        #Open a figure for a particular file
        #to place the legend beneath the figure neatly
        #ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        #get the station, event-name and the DT value
        #get the station coordinates and the event coordinates
        station_lat, station_lon =  DicReg[item]['stn_coord']
        event_lat,   event_lon   =  DicReg[item]['event_coord']
        #print(event_lat, event_lon, station_lat, station_lon)
        station_lat = float(station_lat); station_lon = float(station_lon)
        event_lat   = float(event_lat) ;  event_lon   = float(event_lon)
        ############## Check ##############################
        DT_L        = DicRegPath[item]
        #Check if the event-station contain a bad value
        if(path_average== True):
            if(len(DT_L) > 1):
                DT  = np.mean(np.asarray(DT_L)) 
                data_RF.append(DT)
            else:
                DT = float(DT_L[0])
                data_RF.append(DT)
        elif(path_absolue==True):
            if(len(DT_L) > 1):
                DT = max(abs(x) for x in  np.asarray(DT_L))
                data_RF.append(DT)
            else:
                DT = float(DT_L[0])
                data_RF.append(DT)
        #plot events on the map
        x_ev, y_ev  = m(event_lon, event_lat)
        #plot the event on the map
        #m.plot(x_ev, y_ev, 'o', markersize=msize, c='lime', markeredgecolor='k', label='Events')
        #m.plot(x_ev, y_ev, 'o', markersize=msize, c='r', markeredgecolor='k')
        m.plot(x_ev, y_ev, 'o', markersize=msize, c='k')
        #Plots stations on the map
        x, y = m(station_lon, station_lat)
        #m.plot(x, y, '^', markersize = msize,     c='k', markeredgecolor='r')
        #m.plot(x, y, '^', markersize = msize,     c='k')
        m.plot(x, y, '^', markersize = 12,     c='k')
        #Draw the ray path
        #if(abs(DT) > 0.5):
        if(DT > -0.5 and   DT < 0.5):
            m.drawgreatcircle(event_lon, event_lat, station_lon, station_lat, linewidth=0.5, color='grey',alpha=0.2, zorder=1)
            paths_gray +=1
        elif(DT < -0.5):
            m.drawgreatcircle(event_lon, event_lat, station_lon, station_lat, linewidth=0.9, color='r',alpha=1.0, zorder =5)
            paths_red +=1
        elif(DT > 0.5):
            m.drawgreatcircle(event_lon, event_lat, station_lon, station_lat, linewidth=1.0, color='b',alpha=1.0, zorder=5)
            paths_blue +=1
    #Set the Tile of the map
    #get the DT valueto form the database
    #data_R = np.asarray([v for it in DicRegPath for v in DicRegPath[it] ])
    data_R = np.asarray(data_RF)
    #data_R = np.asarray([it for it in DicRegPath[p] for p in DicRegPath.keys()])
    #########Compute the Percentiles of the data and trim the latter ######################
    Q_LR, Q_HR, data_trim_R = Percentiles(data_R,percs)
    ################### Computer the parameters for the trim data ############################
    ##Get the std, mean, median, number of bin and the binwidth
    mu_R,sigma_R, med_R,  nb_,  bwidth_R  = parameters(data_trim_R)
    ##linspace between minimum and the maximum on the x-axis
    XX_R   = np.linspace(Q_LR, Q_HR, data_trim_R.size)
    ########################################################
    ####Compute the Probability Density Function (PDF) ##############################
    P_R    = norm.pdf(XX_R, mu_R, sigma_R)
    ###############################################
    nbins_R_sym =  SYMETRIC_BIN(data_trim_R)
    ####################################################
    N_R, N_Limit_R, bins_R, patches_R =  count_bin(data_R, nbins_R_sym)
    ##########################################
    YY_max_R      =float('%.2f'%(N_Limit_R))
    #################################################
    #################################################
    ax2.grid(visible=True, axis='both',  alpha = 0.9)
    #Args, axis, data, X-poistion of pdf,  PDF, Ymax, bins_sym (symetric), color, label
    #Scatter the heighest point of the histogram
    ################ PLOT THE REGIONAL DATA ###################################################
    Plot_hist(ax2, data_R, XX_R, P_R, nbins_R_sym, AxiVis[mode]['rc'], 'REG-M25')
    #Scatter the heighest point of the histogram
    Scatter_R(ax2, YY_max_R, data_trim_R,AxiVis[mode]['rc'])
    #Get y-limit of the ax1-axis
    ymin, ymax = ax2.get_ylim()
    ##########################################
    #Font of the y-axis
    ax2.yaxis.set_tick_params(labelsize=fontsize_Label)
    #Font x-axis
    ax2.xaxis.set_tick_params(labelsize=fontsize_Label)
    #set limit of the axis
    ax2.set_ylim(0.0, YY_max)
    ax2.set_xlim(Fig_params['Xmin'], Fig_params['Xmax'], emit=True)
    ### Create additional axis########################################
    ax22 =  ax2.twinx()
    #############################
    ax22.set_ylabel(WaveTpes[File], color = 'k', rotation='vertical',fontweight='bold', fontsize=fsize_abc, labelpad=5.0)
    ax22.set_yticks([])
    ######################

    ax1.set_title(WaveTpes[File], color = 'k', fontweight='bold', fontsize=fsize_abc, pad = 8.0)
    #args, axis, x-coordinate, y-coordinate, data (trim data), percentiles list, fontsize
    AxisWrite2(ax2,  t1, YY_max, data_trim_R, percs, fontsize_Txt)
    #args, axis, x-coordinate, y-coordinate, data (none-trim data), percentiles list, fontsize
    AxisWrite2(ax2, -t1, YY_max, data_R, percs_O, fontsize_Txt)
    #######################################################################
    ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
    ax2.xaxis.set_major_locator(mticker.FixedLocator(locs))
    ax2.set_xlabel(r'$\Delta T/T$'+'(%)', fontsize=fontsize_Label, color='k', labelpad=10, fontweight= 'bold')
    ############################ plt.tight_layout() ########################################
    #plot_plates(m, linewidth=5.0, linestyle='-', color='firebrick', zorder=1000, alpha=0.5)
    #plot_plates(m, linewidth=5.0, linestyle='-', color='firebrick', zorder=1000)
    plot_plates(m, linewidth=3.0, linestyle='-', color='k', zorder=1000, alpha=0.6)
    
    plot_hotspots(m, marker='^', ls='', markerfacecolor=(0,1,0.1),
                 markeredgecolor='k', zorder=1000, ms=12)
    #Plot the check
    ax1.plot(xk,yk,c='k',lw=5.0,linestyle='--') # NB Geodetic est bien pour placer les points sur la carte
    #plt.title(r'$\Delta T/T$'+' < 0.5 %',fontsize=24, loc='center',  color="k", pad =25)
    #ax1.title(r'$\Delta T/T$'+' < 0.5 % (Red)   and  ' + r'$\Delta T/T$'+' > 0.5 %  (Blue)' ,fontsize=24, loc='center',  color="k", pad =25)
    plt.title(r'$\Delta T/T$'+' < -0.5 % (red)   and  ' + r'$\Delta T/T$'+' > 0.5 %  (blue)' ,fontsize=24, loc='center',  color="k", pad =25)
    #plt.title(r'$\Delta T/T$'+' < -0.5 % ('+ r"$\textcolor{red}{red}$"+')'  +' and  ' + r'$\Delta T/T$'+' > 0.5 % ('+ r"$\textcolor{blue}{blue}$"+')' ,fontsize=24, loc='center',  color="k", pad =25)
    #r"\textcolor{red}{This is red text in LaTeX}"
    #plt.title(r'$\Delta T/T$'+' < -0.5 % '+ '('+colored("red", "red")+')' +' and  ' + r'$\Delta T/T$'+' > 0.5 % '+ '('+colored("blue", "blue")+')'  ,fontsize=24, loc='center',  color="k", pad =25)
    ################################################################################
    #fig.savefig(figname)
    plt.savefig(figname, bbox_inches="tight")
    plt.close(fig)
    #Write number of paths into a file
    fsname = "Path_%s.dat"%(File.replace(".dat", ""))
    fs     = open(fsname, "w") 
    fs.write("# number of path for each color\n") 
    fs.write("Number of grey paths: %d \n"%(paths_gray))
    fs.write("Number of red paths:  %d \n"%(paths_red))
    fs.write("Number of blue paths: %d \n"%(paths_red))
