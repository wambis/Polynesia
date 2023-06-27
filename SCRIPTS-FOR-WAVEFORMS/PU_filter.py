#!/usr/bin/env python
import h5py
import re,yaml
import sys,os
import re, json
sys.path.append('/home/mw1685/libs/ucbpy3')
from ndk_rec_file import NDKFile
from UCBFilter import Filter
import argparse
from obspy.signal.invsim import cosine_sac_taper
from obspy.signal.util import _npts2nfft
import numpy as np
#from .rotate import rotate_stream
from pypaw import ProcASDF
#########################################
from pyasdf import ASDFDataSet
from pytomo3d import signal, adjoint
from glob import glob
#sys.path.append("/Volumes/wamba/Plot-Models-2021/ucbpy")
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.interpolate import InterpolatedUnivariateSpline
DTYPE = np.float32
np.set_printoptions(threshold=np.inf)
import obspy
#from scipy import signal1
import scipy 
from yaml.loader import SafeLoader
import obspy.signal.filter as fltobspy
from obspy.geodetics.base import gps2dist_azimuth
from obspy import read, read_inventory
from pytomo3d.window.io import get_json_content, WindowEncoder




def process_function(st, inv):
    #Attached response
    st.attach_response(inv)
    #processing the stream with the pytomo3d packages
    st = signal.process_stream(st,  inventory=inv, remove_response_flag = params["remove_response_flag"], water_level=int(params["water_level"]),
         pre_filt      = params["pre_filt"], filter_flag        = params["filter_flag"], starttime=starttime, endtime =endtime,
         resample_flag = params["resample_flag"], sampling_rate = float(params["sampling_rate"]),
         taper_type    = params["taper_type"], taper_percentage = float(params["taper_percentage"]), rotate_flag=params["rotate_flag"],
         sanity_check  = params["sanity_check"], event_latitude = event_latitude, event_longitude = event_longitude) 



def process_tream_func(stream, param_filt):
    st = signal.process_stream(stream, inventory=None, remove_response_flag=False, water_level=60,
    pre_filt      = param_filt["pre_filt"],   filter_flag      = param_filt["filter_flag"],  starttime=None, endtime=None,
    resample_flag = False, sampling_rate=1.0, taper_type       = param_filt["taper_type"],   taper_percentage = float(param_filt["taper_percentage"]),
    rotate_flag   = False, sanity_check = None, event_latitude = None, event_longitude = None)
    return st

##############################################
def LagCross(d, s):
     #Compute the cross correlation
    #corr = signal.correlate(obs_select - np.mean(obs_select), syn_select - np.mean(syn_select), mode="full")
    #corr = signal.correlate(d, s, mode="full")
    corr = scipy.signal.correlate(d, s, mode="full")
    #get the corresponding tau (i.e. lags)
    #lags = signal.correlation_lags(len(d), len(s), mode="full")
    lags = scipy.signal.correlation_lags(len(d), len(s), mode="full")
    corr /= np.max(corr)
    #grab sifted positions, get the corresponding tau (i.e. lags)
    correlation_max_index = np.argmax(corr)
    #get the lag corresponding to the maximum value of the cross-correlation
    tau  = lags[correlation_max_index]
    return tau

#Define the correlation funuction
def correlation(d, s):
       proc1   = (d -np.mean(d)) * (s -np.mean(s))
       proc2   = (d -np.mean(d)) * (d -np.mean(d))
       proc3   = (s -np.mean(s)) * (s -np.mean(s))
       c       = proc2.sum() * proc3.sum()
       corr    = proc1.sum()/np.sqrt(c)
       #corr    = "%.2f"%(corr)
       return corr


def Resample(data_in):
    #Resample to 1pt/s, set time_inv =1.0 then
    t_inv = 1.0
    #t0    = 500
    t0    = 0
    x_in  = data_in[:,0]
    #Reset the time
    x     = x_in - t0
   #######################
    y     = data_in[:,1]
    spl   = InterpolatedUnivariateSpline(x,y)
    x_new = np.arange(min(x), max(x), t_inv)
    y_new = spl(x_new)
    #define the dimension of the array
    dim   = (x_new.shape[0],2)
    #Define an empty 2-D array
    data_out = np.zeros(dim)
    #Fill the array
    data_out[:,0] = x_new 
    data_out[:,1] = y_new 
    return data_out


          

#def make_stream(data_in):
def make_stream(data):
    #tr = Trace(data=np.asarray(ds_UCB[:,1], np.float64))
    from obspy import Stream, Trace
    #Resample to 1pt/s, set time_inv =1.0 then
    t_inv = 1.0
    trace  = Trace()
    #########################
    stream = Stream()
    data_in = Resample(data)
    ########################
    t_     = data_in[:,0]
    ######set value to trace #################
    trace.times = t_
    trace.data  = data_in[:,1]
    trace.stats.sampling_rate = t_inv
    #Make SEMUCB Trace
    st = stream.append(trace)
    return st


###############################################
with open('config-mw.yaml') as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    params = yaml.load(Fym, Loader=SafeLoader)
############
#station test
net_stn_test  = params['net_stn']
#Grab the event id
event_id  = params['event_id']
#Grab windowns files
#get the json file
#Grab the frequency band
periodband= params['periodband']
#get the mode 
mode      = params['mode']
#get the color list in order to plot windows with different color
########################
par_file_name= params['par_file']
#Window file name
mode               = mode.split("#")[0]
#grab the parameter file
par_file  = open(par_file_name, "r")
param_filt= yaml.safe_load(par_file)
#make a list
DICT_STN_NETWK_UCB = [net_stn_test]
######################
fucb   = "UCB.G.KIP.Z.dat"
##Open the SEMUCB synthetic
ds_UCB = np.loadtxt(fucb) 
d = ds_UCB[:,1]
t = ds_UCB[:,0]
#Make a stream from data
syn_ucb_st = make_stream(ds_UCB)
####################################
ds_GLAD    = np.loadtxt("G.KIP.dat") 
syn_glad_st= make_stream(ds_GLAD)
##################################
st_ucb  = process_tream_func(syn_ucb_st, param_filt)
st_glad = process_tream_func(syn_glad_st, param_filt)
##############Grab the sampling rate##############################
smplintv_ucb   = st_ucb[0].stats.sampling_rate
smplintv_glad  = st_glad[0].stats.sampling_rate 
################get the time and the data from trace####################################
t_ucb       = st_ucb[0].times 
#shifted the data from UCB
t_ucb_shift = t_ucb -500
#t_ucb_shift = t -500
t_glad      = st_glad[0].times 
#grab filtered data 
data_ucb    = st_ucb[0].data
data_glad   = st_glad[0].data
##################################
##################Create a figure #####################
fig, axs = plt.subplots(1,figsize=(12, 4))
#figure name 
fig_name = "%s.png"%("FIG-PU-FILT")
################################## Plot normalized traces ################################
axs.plot(t_ucb_shift,      data_ucb/max(data_ucb), 'r', linewidth=1., alpha = 1.0,label='SEMUCB')
axs.plot(t_glad,           data_glad/max(data_glad),    'b', linewidth=1., alpha = 1.0, label='GLADM25')
#put the legend
axs.set_xlim(0, max(t_ucb_shift))
##################
axs.set_title('PU-Filter      %s  %s'%(event_id, net_stn_test) + '  '+ 'BP='+periodband.replace('_','-')+' s', loc = 'center',fontsize=14)
axs.set_ylabel('Displacement', color = 'k', fontsize=14, labelpad=10)
#Set axe for component
axs_c    = axs.twinx()
#########################
axs.legend(prop={'size': 9})
#Set the labels
axs.set_xlabel(r'Time (s) ', labelpad=0.5, fontsize=14)
axs_c.set_ylabel('Z', color = 'k', rotation='horizontal', fontsize=14, labelpad=10)
axs_c.set_yticks([])
#Save the figure
plt.savefig(fig_name,dpi=410,bbox_inches = 'tight', pad_inches = 0)
plt.close(fig)
