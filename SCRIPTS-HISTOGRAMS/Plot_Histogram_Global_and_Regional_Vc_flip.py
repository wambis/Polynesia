#!/usr/bin/python
# -*- coding: utf-8 -*-
import obspy
import numpy as np
import os, sys
from glob import glob
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
    N_Limit          =  N_norm_max + N_norm_max * (54.0/100.0)
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
    nbins_symetric = np.arange(-nbins + 0.5, nbins + 0.5, 1)
    return(nbins_symetric)
####################################################################
#def Plot_hist(ax, data, X, Pdf, Ymax, bins_sym, color,  label):
def Plot_hist(ax, data, X, Pdf, bins_sym, color,  label):
    #Args, ax-axis, data, X-coordinate of PDF, PDF, bins_sym (symetric), color, label
    #ax.hist(data, density=True, color = color, ec='k',lw=2.0, alpha=1.0, rwidth = 0.75, 
    ax.hist(data, density=True, color = color, ec='k',lw=2.0, alpha=1.0, rwidth = 1.0, 
            bins=bins_sym, histtype='bar')
            #bins=bins_sym, histtype='bar', align='mid')
    #Plot the Probability density function
    ax.plot(X, Pdf, color = color,linestyle='-', lw=5.0, alpha=1.0)
    ax.plot(X, Pdf, color = 'k',linestyle='-', lw=2.0, alpha=.9)
    #Scatter  
    #ax.scatter(0.0, Ymax * 0.97, marker ="o", edgecolor ="green",linewidths = 2, s = 50)
    #ax.legend([label])

def Scatter_G(ax, y, data,color) :
    #yp = y * 0.985
    yp = y * 0.8
    mean = np.mean(data) 
    std  = np.std(data) 
    VR = mean + std
    VL = mean - std
    X = np.linspace(VL, VR,10) 
    YY = np.asarray([yp,yp,yp,yp,yp,yp,yp,yp,yp,yp])
    #plotting line within the given range
    ax.scatter(float('%.1f'%(mean)), yp, marker ="o", color=color, edgecolor =color,linewidths = 3, s = 70)
    ax.scatter(float('%.1f'%(mean)), yp, marker ="o", color=color, edgecolor = 'k',linewidths = 3, s = 72, alpha=.6)
    #############################################
    #ax.scatter(0.0, yp, marker ="o", color=color, edgecolor =color,linewidths = 3, s = 70)
    #ax.scatter(0.0, yp, marker ="o", color=color, edgecolor = 'k',linewidths = 3, s = 72, alpha=.6)
    ax.plot(X, YY, color =color, lw=1.5)

def Scatter_R(ax, y, data, color) :
    #y= y * 0.97
    y= y * 0.8
    mean = np.mean(data) 
    std  = np.std(data) 
    VR = mean + std
    VL = mean - std
    X  = np.linspace(VL, VR,10) 
    YY = np.asarray([y,y,y,y,y,y,y,y,y,y])
    # plotting line within the given range
    ax.scatter(float('%.1f'%(mean)), y, marker ="o", edgecolor =color, color=color, linewidths = 3, s = 70)
    ax.scatter(float('%.1f'%(mean)), y, marker ="o", edgecolor ='k', color=color, linewidths = 3, s = 72,alpha=.6)

    #ax.scatter(0.0, y, marker ="o", edgecolor =color, color=color, linewidths = 3, s = 70)
    #ax.scatter(0.0, y, marker ="o", edgecolor ='k', color=color, linewidths = 3, s = 72,alpha=.6)
    #ax.scatter(X, YY, marker ="o", edgecolor ="green",linewidths =0.005)
    ax.plot(X, YY , color =color, lw=1.5)

#Set the density to True if you want to normalize the Y-axis between [0,1]
density = True
################################
#Open the config file
#Fym = open("configG.yaml","r") 
with open('configZRT_Vc.yaml') as Fym:
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
#modes= ['17_40','40_100#body_wave','40_100#surface_wave','90_250']
modes = Fig_params['modes']
####
#Create a figure label according to the different modes, and givent this possible list of labels
FigLab = {'17_40':['a','b','c'], '40_100#body_wave':['d','e','f'], '40_100#surface_wave':['g','h','i'], '90_250':['j','k','l']}
#X-axis visbility dependent on the mode
AxiVis = {'17_40':{'Ax': False, 'rc':'aquamarine', 'gc': 'deeppink'},   '40_100#body_wave':{'Ax': True,  'rc':  'aquamarine', 'gc':'deeppink'},
        '40_100#surface_wave':{'Ax':False, 'rc':'palegreen', 'gc': 'r'}, '90_250':{'Ax':True,      'rc': 'palegreen',   'gc': 'r'}}

#print(FigLab) 
#sys.exit()

#for m in modes:

#Regional station file
fs = open(fname_stn)
REG_STATIONS = set(['%s.%s'%(l.split()[0], l.split()[1]) for l in fs])
##File for regional events
f    = open(fname_evt)
##grab regional events
REG_EVENTS = set([l.split()[0].strip() for l in f])
#Grab all the event used in FWI
#events= sorted([it for it in glob("%s/*"%(path)) if os.path.isdir(it)])
# Loop over the entire modes ##################################################
for mode in modes:
    #Grab all the events for this mode (i.e: 40_100#body_wave) used in FWI
    events  = sorted([it for it in glob("%s/*"%(path)) if (os.path.isdir(it) and EXTENSION(it) == mode)])
    ########## Set Dictionary for Global Events ###################################
    DicGlob = {mode: {'Z': [], 'R': [], 'T': []}}
    DicReg  = {mode: {'Z': [], 'R': [], 'T': []}}
    #category 
    wave_type, bandpass = BP_Type_func(mode)
    #grab the category(example 17_40s)
    category = bandpass.replace(' ','')
    #Check and reset the wave type
    if("body_wave" in mode or "17_40" in mode):
        wave_type = "Body wave"
    elif("surface_wave" in mode or "90_250" in mode):
        wave_type = "Surf wave"
    #Loop over the path contains in the variables events
    for pat in events:
        #Get the Name and the window type of the event
        C          = CheckCategory(pat) 
        event_name = C[0] 
        win_type   = C[1]; 
        #Grab the corresponding window file
        fwin       = "%s/windows.filter.json"%(pat)
        if os.path.isfile(fwin):
            windows= readwindow(fwin)
        #check if the it's body wave category in 17-40
        #################################################
        stations = windows.keys()
        #Loop over the stations
        for station in stations:
            #Loop over the component on the stations
            for component in windows[station].keys():
                #Number of windows on the component
                n           = len(windows[station][component])
                #Loop over the windows that belong to the considered component
                for cwind in windows[station][component]:
                    channel = cwind['channel_id']
                    comp    = channel.strip()[-1]
                    #Grab the weight of the window
                    weight  = float(cwind["window_weight"])
                    #Get the weighted cross-correlation time measurement from the window
                    #dlnA    = float(cwind["dlnA"])
                    center_index    = int(cwind["center_index"])
                    time_inv        = float(cwind["dt"]) 
                    #get the central time of the window
                    relative_starttime = float(cwind["relative_starttime"])
                    relative_endtime   = float(cwind["relative_endtime"])
                    cc_max             = float(cwind["cc_shift_in_seconds"])
                    #Compute the phase-velocity anomaly (which travel-time anomaly divided by the central time of the window)
                    value              = (100 * weight * cc_max)/(center_index * time_inv) 
                    #Save the data into the Global dictionary 
                    DicGlob[win_type][comp].append(value)
                    #Check if the event and the station are in the chunk (i.e. the considered region)
                    if(event_name in REG_EVENTS and station in REG_STATIONS):
                        DicReg[win_type][comp].append(value)


    ############################
    #set the components
    comps = ['Z', 'R', 'T']
    ############## #***************************************#########################
    #create a figures of 2 rows x 3 columns
    fig, axs = plt.subplots(2, 3, sharex=True, figsize=(26, 14)) 
    #Create the figure name
    #figname  = "ZRT_Histogram_%s_Vc.pdf"%(mode.replace("#","_"))
    figname  = "ZRT_Histogram_%s_Vc.png"%(mode.replace("#","_"))
    #Grid the figure
    plt.grid(visible=True, axis='both')
    ####Loop over the components
    for comp in comps:
        #Grad the Global and the Regional data
        data_G       = np.asarray(sorted(DicGlob[mode][comp]))
        data_R       = np.asarray(sorted(DicReg[mode][comp])) 
        ########Compute the Percentiles of the data and trim the latter ######################
        Q_LR, Q_HR, data_trim_R = Percentiles(data_R,percs)
        Q_LG, Q_HG, data_trim_G = Percentiles(data_G,percs)
        ################## Computer the parameters for the trim data ############################
        #Get the std, mean, median, number of bin and the binwidth
        mu_R,sigma_R, med_R,  nb_,  bwidth_R  = parameters(data_trim_R) 
        mu_G,sigma_G, med_G,  nb__, bwidth_G  = parameters(data_trim_G) 
        
        #linspace between minimum and the maximum on the x-axis
        XX_G   = np.linspace(Q_LG, Q_HG, data_trim_G.size)
        XX_R   = np.linspace(Q_LR, Q_HR, data_trim_R.size)
        #######################################################
        ###Compute the Probability Density Function (PDF) ##############################
        P_G    = norm.pdf(XX_G, mu_G, sigma_G)
        P_R    = norm.pdf(XX_R, mu_R, sigma_R)
        ###########################################
        #plot the gaussian if the 
        if density:
            #plot the histogram
            #To center the bins we needd to Compute nbins as an array in which will allow the center the bins around zero
            ###############################################################
            nbins_G_sym =  SYMETRIC_BIN(data_trim_G)
            nbins_R_sym =  SYMETRIC_BIN(data_trim_R)
            ####################################################
            N_G, N_Limit_G, bins_G, patches_G =  count_bin(data_G, nbins_G_sym)
            N_R, N_Limit_R, bins_R, patches_R =  count_bin(data_R, nbins_R_sym)
            ##########################################
            #sys.exit()
            YY_max_G      =float('%.2f'%(N_Limit_G))
            YY_max_R      =float('%.2f'%(N_Limit_R))
            YY_max        = max(YY_max_G, YY_max_R)
            #####################
            if(comp=="Z"):
                    ####################################
                    ax1 = axs[0,0] ; ax2 = axs[1,0]
                    ax1.grid(visible=True, axis='both',  alpha = 0.9)
                    ax2.grid(visible=True, axis='both',  alpha = 0.9)
                    #Args, axis, data, X-poistion of pdf,  PDF, Ymax, bins_sym (symetric), color, label
                    Plot_hist(ax1, data_G, XX_G, P_G, nbins_G_sym, AxiVis[mode]['gc'], 'M25')
                    #Scatter the heighest point of the histogram
                    Scatter_G(ax1, YY_max_G, data_trim_G,AxiVis[mode]['gc']) 
                    ################ PLOT THE REGIONAL DATA ###################################################
                    Plot_hist(ax2, data_R, XX_R, P_R, nbins_R_sym, AxiVis[mode]['rc'], 'REG-M25')
                    #Plot_hist(ax2, data_R, XX_R, P_R, nbins_G_sym, 'palegreen', 'REG-M25')
                    #Scatter the heighest point of the histogram
                    Scatter_R(ax2, YY_max_R, data_trim_R,AxiVis[mode]['rc'])
                    #Get y-limit of the ax1-axis
                    ymin, ymax = ax1.get_ylim()
                    #set limit of the axis
                    ax1.set_ylim(0.0, YY_max)
                    ax2.set_ylim(0.0, YY_max)
                    #### inverted the y-axis now ##########
                    ax2.invert_yaxis()
                    #Reduce the space between the two ax-axis##################
                    #plt.subplots_adjust(hspace=0.05)
                    plt.subplots_adjust(hspace=0.0)
                    ############## Write text on the ax-axis ######################################################
                    #args, axis, x-coordinate, y-coordinate, data (trim data), percentiles list, fontsize
                    AxisWrite1(ax1,  t1, ymax, data_trim_G, percs, fontsize_Txt)
                    #args, axis, x-coordinate, y-coordinate, data (none-trim data), percentiles list, fontsize
                    AxisWrite1(ax1, -t1, ymax, data_G, percs_O, fontsize_Txt)
                    ############## Write text on the ax2-axis ######################################################
                    #args, axis, x-coordinate, y-coordinate, data (trim data), percentiles list, fontsize
                    AxisWrite2(ax2,  t1, ymax, data_trim_R, percs, fontsize_Txt)
                    #args, axis, x-coordinate, y-coordinate, data (none-trim data), percentiles list, fontsize
                    AxisWrite2(ax2, -t1, ymax, data_R, percs_O, fontsize_Txt)
                    #######################################################################
                    #Write a lettter in Figure
                    #ax2.set_xlabel(r'$\Delta T/T$'+'(%)', fontsize=fontsize_Label, color='k', labelpad=10, fontweight= 'bold')
                    ax1.set_title(label= "%s,   Component:%s,    BP=%s "%(wave_type,comp,bandpass),alpha=1.0,fontweight='bold',fontsize=fontsize_Title, loc='center')
                    ax1.text(t2,ymax/8,'(%s)'%(FigLab[mode][0]),horizontalalignment='center',verticalalignment='center',alpha=1.0,fontweight='bold',fontsize=fsize_abc)
                    #Font of the y-axis
                    ax1.yaxis.set_tick_params(labelsize=fontsize_Label)
                    ax2.yaxis.set_tick_params(labelsize=fontsize_Label)
                    #Font x-axis
                    ax2.xaxis.set_tick_params(labelsize=fontsize_Label)
                    #Set the X-limit 
                    #x_ticks_labels = np.arange(Fig_params['Xmin'],Fig_params['Xmax'],5)
                    #ax1.set_xticklabels(x_ticks_labels)
                    #ax2.set_xticklabels(x_ticks_labels)
                    ax1.set_xlim(Fig_params['Xmin'], Fig_params['Xmax'], emit=True)
                    ax2.set_xlim(Fig_params['Xmin'], Fig_params['Xmax'], emit=True)
                    #################### Format Axis ##########################
                    ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
                    ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
                    ax1.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
                    ax2.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
                    #Axis visibility
                    #ax2.get_xaxis().set_visible(AxiVis[mode]['Ax'])    
                    if(AxiVis[mode]['Ax']) == False:
                        #Remove tick-lable on y-axis
                        ax2.set_xticklabels([])
                        #ax2.set_yticklabels([])
                    else:
                    #Write a lettter in Figure
                        ax2.set_xlabel(r'$\Delta T/T$'+'(%)', fontsize=fontsize_Label, color='k', labelpad=10, fontweight= 'bold')

                    #ax2.get_yaxis().set_visible(True)    
                    #Grid the figure
##################################################################################################
            if(comp=="R"):
                    ax1 = axs[0,1]; ax2 = axs[1,1]
                    #Grid the figure
                    ax1.grid(visible=True, axis='both',  alpha = 0.9)
                    ax2.grid(visible=True, axis='both',  alpha = 0.9)
                    #Args, axis, data, X-poistion of pdf,  PDF, Ymax, bins_sym (symetric), color, label
                    Plot_hist(ax1, data_G, XX_G, P_G, nbins_G_sym, AxiVis[mode]['gc'], 'M25')
                    #Scatter the heighest point of the histogram
                    #Scatter_G(ax1, YY_max_G, data_trim_G,'r') 
                    Scatter_G(ax1, YY_max, data_trim_G,AxiVis[mode]['gc']) 
                    ################ PLOT THE REGIONAL DATA ###################################################
                    Plot_hist(ax2, data_R, XX_R, P_R, nbins_R_sym, AxiVis[mode]['rc'], 'REG-M25')
                    #Plot_hist(ax2, data_R, XX_R, P_R, nbins_G_sym, 'palegreen', 'REG-M25')
                    #Scatter the heighest point of the histogram
                    #Scatter_R(ax2, YY_max_R, data_trim_R, 'palegreen')
                    Scatter_R(ax2, YY_max, data_trim_R, AxiVis[mode]['rc'])
                    #Get y-limit of the ax1-axis
                    ymin, ymax = ax1.get_ylim()
                    #set limit of the axis
                    ax1.set_ylim(0.0, YY_max)
                    ax2.set_ylim(0.0, YY_max)
                    #### inverted the y-axis now ##########
                    ax2.invert_yaxis()
                    #Reduce the space between the two ax-axis##################
                    #plt.subplots_adjust(hspace=0.05)
                    plt.subplots_adjust(hspace=0.0)
                    ############## Write text on the ax-axis ######################################################
                    #args, axis, x-coordinate, y-coordinate, data (trim data), percentiles list, fontsize
                    AxisWrite1(ax1,  t1, ymax, data_trim_G, percs, fontsize_Txt)
                    #args, axis, x-coordinate, y-coordinate, data (none-trim data), percentiles list, fontsize
                    AxisWrite1(ax1, -t1, ymax, data_G, percs_O, fontsize_Txt)
                    ############## Write text on the ax2-axis ######################################################
                    #args, axis, x-coordinate, y-coordinate, data (trim data), percentiles list, fontsize
                    AxisWrite2(ax2,  t1, ymax, data_trim_R, percs, fontsize_Txt)
                    #args, axis, x-coordinate, y-coordinate, data (none-trim data), percentiles list, fontsize
                    AxisWrite2(ax2, -t1, ymax, data_R, percs_O, fontsize_Txt)
                    #######################################################################
                    #Write a lettter in Figure
                    #ax2.set_xlabel(r'$\Delta T/T$'+'(%)', fontsize=fontsize_Label, color='k', labelpad=10, fontweight= 'bold')
                    ax1.set_title(label= "%s,   Component:%s,    BP=%s "%(wave_type,comp,bandpass),alpha=1.0,fontweight='bold',fontsize=fontsize_Title, loc='center')
                    #Write the figure label
                    ax1.text(t2,ymax/8,'(%s)'%(FigLab[mode][1]),horizontalalignment='center',verticalalignment='center',alpha=1.0,fontweight='bold',fontsize=fsize_abc)
                    #Axis visibility
                    #ax1.get_yaxis().set_visible(False)    
                    #ax2.get_yaxis().set_visible(False)    
                    #Font of the y-axis
                    ax1.yaxis.set_tick_params(labelsize=fontsize_Label)
                    #Font x-axis
                    ax2.xaxis.set_tick_params(labelsize=fontsize_Label)
                    #Axis visibility
                    #ax2.get_xaxis().set_visible(AxiVis[mode]['Ax'])    
                    #Set the X-limit 
                    #x_ticks_labels = np.arange(Fig_params['Xmin'],Fig_params['Xmax'],5)
                    #ax1.set_xticklabels(x_ticks_labels)
                    #ax2.set_xticklabels(x_ticks_labels)
                    ax1.set_xlim(Fig_params['Xmin'], Fig_params['Xmax'], emit=True)
                    ax2.set_xlim(Fig_params['Xmin'], Fig_params['Xmax'], emit=True)
                    #################### Format Axis ##########################
                    ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
                    ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
                    ax1.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
                    ax2.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
                    ############################################################
                    #Remove y-label ########################################
                    ax1.set_yticklabels([])
                    ax2.set_yticklabels([])
                    #################################################
                    if(AxiVis[mode]['Ax']) == False:
                        #Remove tick-lable on y-axis
                        ax2.set_xticklabels([])
                    else:
                        #Write a lettter in Figure
                        ax2.set_xlabel(r'$\Delta T/T$'+'(%)', fontsize=fontsize_Label, color='k', labelpad=10, fontweight= 'bold')
            if(comp=="T"):
                    ax1 = axs[0,2] ; ax2 = axs[1,2]
                    #Grid the figure
                    ax1.grid(visible=True, axis='both',  alpha = 0.9)
                    ax2.grid(visible=True, axis='both',  alpha = 0.9)
                    #Args, axis, data, X-poistion of pdf,  PDF, Ymax, bins_sym (symetric), color, label
                    Plot_hist(ax1, data_G, XX_G, P_G, nbins_G_sym, AxiVis[mode]['gc'], 'M25')
                    #Scatter the heighest point of the histogram
                    #Scatter_G(ax1, YY_max_G, data_trim_G,'r') 
                    Scatter_G(ax1, YY_max, data_trim_G,AxiVis[mode]['gc']) 
                    ################ PLOT THE REGIONAL DATA ###################################################
                    Plot_hist(ax2, data_R, XX_R, P_R, nbins_R_sym, AxiVis[mode]['rc'], 'REG-M25')
                    #Plot_hist(ax2, data_R, XX_R, P_R, nbins_G_sym, 'palegreen', 'REG-M25')
                    #Scatter the heighest point of the histogram
                    #Scatter_R(ax2, YY_max_R, data_trim_R,'palegreen')
                    Scatter_R(ax2, YY_max, data_trim_R,AxiVis[mode]['rc'])
                    #Get y-limit of the ax1-axis
                    ymin, ymax = ax1.get_ylim()
                    #set limit of the axis
                    ax1.set_ylim(0.0, YY_max)
                    ax2.set_ylim(0.0, YY_max)
                    #### inverted the y-axis now ##########
                    ax2.invert_yaxis()
                    #Reduce the space between the two ax-axis##################
                    #plt.subplots_adjust(hspace=0.05)
                    plt.subplots_adjust(hspace=0.0)
                    ############## Write text on the ax-axis ######################################################
                    #args, axis, x-coordinate, y-coordinate, data (trim data), percentiles list, fontsize
                    AxisWrite1(ax1,  t1, ymax, data_trim_G, percs, fontsize_Txt)
                    #args, axis, x-coordinate, y-coordinate, data (none-trim data), percentiles list, fontsize
                    AxisWrite1(ax1, -t1, ymax, data_G, percs_O, fontsize_Txt)
                    ############## Write text on the ax2-axis ######################################################
                    #args, axis, x-coordinate, y-coordinate, data (trim data), percentiles list, fontsize
                    AxisWrite2(ax2,  t1, ymax, data_trim_R, percs, fontsize_Txt)
                    #args, axis, x-coordinate, y-coordinate, data (none-trim data), percentiles list, fontsize
                    AxisWrite2(ax2, -t1, ymax, data_R, percs_O, fontsize_Txt)
                    #######################################################################
                    #Write a lettter in Figure
                    #ax2.set_xlabel(r'$\Delta T/T$'+'(%)', fontsize=fontsize_Label, color='k', labelpad=10, fontweight= 'bold')
                    ax1.set_title(label= "%s,   Component:%s,    BP=%s "%(wave_type,comp,bandpass),alpha=1.0,fontweight='bold',fontsize=fontsize_Title,loc='center')
                    #write the label on the figure
                    ax1.text(t2, ymax/8,'(%s)'%(FigLab[mode][2]),horizontalalignment='center',verticalalignment='center',alpha=1.0,fontweight='bold',fontsize=fsize_abc)
                    #Font of the y-axis
                    ax1.yaxis.set_tick_params(labelsize=fontsize_Label)
                    #Font x-axis
                    ax2.xaxis.set_tick_params(labelsize=fontsize_Label)
                    #Set the X-limit 
                    #x_ticks_labels = np.arange(Fig_params['Xmin'],Fig_params['Xmax'],5)
                    #ax1.set_xticklabels(x_ticks_labels)
                    #ax2.set_xticklabels(x_ticks_labels)
                    ax1.set_xlim(Fig_params['Xmin'], Fig_params['Xmax'], emit=True)
                    ax2.set_xlim(Fig_params['Xmin'], Fig_params['Xmax'], emit=True)
                    #Axis visibility
                    #ax2.get_xaxis().set_visible(AxiVis[mode]['Ax'])    
                    #################################
                    ax11 =  ax1.twinx()
                    ax22 =  ax2.twinx()
                    #########################
                    ax11.set_ylabel('GLOBAL', color = 'k', rotation='vertical', fontweight='bold', fontsize=fsize_abc, labelpad=2.0)
                    ax11.set_yticks([])
                    #############################
                    ax22.set_ylabel('REGIONAL', color = 'k', rotation='vertical',fontweight='bold', fontsize=fsize_abc, labelpad=2.0)
                    ax22.set_yticks([])
                    #################### Format Axis ##########################
                    ax1.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
                    ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
                    ax1.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
                    ax2.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
                    #############################################################
                    #Remove y-labels ########################################
                    ax1.set_yticklabels([])
                    ax2.set_yticklabels([])
                    #################################################
                    if(AxiVis[mode]['Ax']) == False:
                        #Remove tick-lable on x-axis
                        ax2.set_xticklabels([])
                    else:
                        #Write a lettter in Figure
                        ax2.set_xlabel(r'$\Delta T/T$'+'(%)', fontsize=fontsize_Label, color='k', labelpad=10, fontweight= 'bold')
    ################################################################################
    fig.savefig(figname)
    plt.close(fig)
