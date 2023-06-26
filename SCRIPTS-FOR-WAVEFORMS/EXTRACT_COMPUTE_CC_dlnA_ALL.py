import obspy
import numpy as np
import sys,os,math,json
import re,yaml
sys.path.append('/home/mw1685/libs/ucbpy3')
from ndk_rec_file import NDKFile
from glob import glob
import statistics
from scipy import stats
from scipy.stats import norm
import matplotlib.pyplot as plt
from operator import truediv
from math import log10, floor, ceil
from matplotlib.ticker import FormatStrFormatter
from scipy import signal
import obspy
import obspy.signal.filter as fltobspy
from obspy.geodetics.base import gps2dist_azimuth
from yaml.loader import SafeLoader
from pyasdf import ASDFDataSet
from pytomo3d.window.io import get_json_content, WindowEncoder
import time

#get the starting time
start_time = time.time()


def EXTENSION(string):
    #split by the point
    Z_, ext_file = string.split(".")
    return(ext_file)

def CheckCategory(path_file):
    Dirbasename          = os.path.basename(path_file) 
    event_name,ext_file  = Dirbasename.split(".")
    #Open a dictionary to keep the measurements 
    if "17_40" in ext_file:
        tmp        = "17_40"
    elif "#body_wave" in ext_file and "40_100" in  ext_file:
        tmp        = "40_100#body_wave"
    elif "#surface_wave" in ext_file and "40_100" in ext_file: 
        tmp        = "40_100#surface_wave"
    elif "90_250" in ext_file: 
        tmp        = "90_250"
    return (event_name, tmp)



def BP_Type_func(ctype):
    #Open a dictionary to keep the measurements
    if "17_40" in  ctype:
        bp         = "17-40 s"
        tmp        = "body_wave"
    elif "#body_wave" in ctype and "40_100" in  ctype:
        bp         = "40-100 s"
        tmp        = "body_wave"
    elif "#surface_wave" in ctype and "40_100" in ctype:
        bp         = "40-100 s"
        tmp        = "surface_wave"
    elif "90_250" in ctype:
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


#Function to open a ASDF file
def openASDF(ASDF_file):
    try:
        ds = ASDFDataSet(ASDF_file, mode="r")
    except:
        try:
            ds = ASDFDataSet(ASDF_file, mode="a")
        except:
            print("Check the input file: %s"%(ASDF_file))
            sys.exit(1)
    return ds

#Function to extract a trace
def ExtractTrace(st, comp):
    #tr      = st.select(component=comp)[0]
    tr      = st.select(id=comp)[0]
    #tr       = list(tr_)[0]
    #Convert data to micro-meter
    #tr.data = tr.data * 1000000
    Amax    = max(tr.data)
    #get the 10% of the Amax
    Amax_per_cent= (Amax * 10.0)/100.0
    #get the final Amax
    Amax_F  = math.ceil(Amax + Amax_per_cent)
    return (tr, Amax_F)

##Define the correlation funuction
def correlation(d, s):
       proc1   = (d -np.mean(d)) * (s -np.mean(s))
       proc2   = (d -np.mean(d)) * (d -np.mean(d))
       proc3   = (s -np.mean(s)) * (s -np.mean(s))
       c       = proc2.sum() * proc3.sum()
       corr    = proc1.sum()/np.sqrt(c)
       #corr    = "%.2f"%(corr)
       return corr

##Define the correlation funuction
#def correlation(d, s):
#       dxs   = (d) * (s)
#       corr    = dxs.sum()/np.sqrt((d ** 2).sum() * (s ** 2).sum())
#       #corr    = "%.2f"%(corr)
#       return corr

##def _xcorr_win(self, d, s):
#def Xcorr_win(d, s):
#        cc = np.correlate(d, s, mode="full")
#        time_shift = cc.argmax() - len(d) + 1
#        # Normalized cross correlation.
#        max_cc_value = cc.max() / np.sqrt((s ** 2).sum() * (d ** 2).sum())
#        return max_cc_value, time_shift
##########################################################
#get the rms
def VR(data, syn):
    diff2  = (data - syn)**2
    syn_   = syn    - np.mean(syn)
    data_  = data   - np.mean(data)
    diff2  = (data_ - syn_)**2
    data2  = (data_)**2
    V      = diff2.sum()/data2.sum()
    rms    = np.sqrt(V)
    return rms

def weigth_func(right_index, left_index, dt, period, CC):
    w  = (right_index - left_index) * dt / (period * CC) 
    return w

#Extract the window from a raw data
def get_window_from_trace(trace, window):
    starttime   = window['absolute_starttime']
    start_win   = window['relative_starttime']
    ended_win   = window['relative_endtime']
    #grab the index
    left_indx   = window['left_index']
    right_indx  = window['right_index']
    center_indx = window['center_index']
    #Grab the observed and the synthetics in corresponding window#########
    select_data = trace.data[left_indx : right_indx]
    select_time = trace.times()[left_indx    : right_indx]
    return(start_win, ended_win, select_data,select_time)

def LagCross(d, s):
    #Compute the cross correlation
    corr = signal.correlate(d, s, mode="full")
    #get the corresponding tau (i.e. lags)
    lags = signal.correlation_lags(len(d), len(s), mode="full")
    corr /= np.max(corr)
    #grab sifted positions, get the corresponding tau (i.e. lags)
    correlation_max_index = np.argmax(corr)
    #get the lag corresponding to the maximum value of the cross-correlation
    tau  = lags[correlation_max_index]
    return tau
##################################
def Delta_Amp(d, s, t1, t2): 
    #formula Asyn = \sqrt{[1/(t2 -t1) ] * \int_{t1}^{t2} (u_syn)^2 \,dt
    #formula Aobs = \sqrt{[1/(t2 -t1) ] * \int_{t1}^{t2} (u_obs)^2 \,dt
    d2   = d**2
    s2   = s**2
    Asyn = np.sqrt( (1.0/(t2 -t1)) * s2.sum() ) 
    Aobs = np.sqrt( (1.0/(t2 -t1)) * d2.sum() ) 
    dlnA = (Aobs - Asyn)/Asyn
    return dlnA

####################
def dlnA_win(d, s):
    return 0.5 * np.log(np.sum(d ** 2) / np.sum(s ** 2))


#####################################################################

with open('configExtract.yaml') as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)

#Grab the path of the data
path     = Fig_params['path_wid']
#Grab the regional station file and events file
path_obs = Fig_params['path_obs']
path_syn = Fig_params['path_syn']
#win_file = Fig_params['win_file']
#Grab the modes
modes    = Fig_params['modes']
#Grab the station file
FG       = Fig_params['FstnG']

F_info   = Fig_params['Finfo']
#Grab the global catalog
fname_evcatalogue = Fig_params['fname_evcatalogue']
###########################
path ="/mnt/tigress/TROMP/data/DATA_M25/window"
#Global stations:

#open the global station file 
f_i = open(F_info, 'r') 
#Grab all the event-station pairs with the corresponding coordinates as described below
################  event_id       station      event_lat     event_lon     slat            slon       EPIC_DIST
INFOS_DICT = {(l.split()[0], l.split()[3]): (l.split()[1], l.split()[2], l.split()[4], l.split()[5], l.split()[6]) for l in f_i}
#Read events catalog ###################################
#evcatalogue = NDKFile(fname_evcatalogue)
#print('reading event catalogue...')
#print('found {:d} events in the catalogue'.format(evcatalogue.nevents))
#ALL mode used in the full waveform inversion
#ll_   = "#event_id event_lat event_lon station netwk.stn.channel.comp slat slon EPIC_DIST relative_starttime relative_endtime shift(in s) CC_value(t) CC_value(0) wind_weight dlnA(t) dlnA(0)"
#

ll_="#event_id event_lat event_lon station channel_id        slat  slon    EPIC_DIST rtive_starttime rtive_endtime tau(in s) tau_new(in s) CC(t) CC_new(t) CC(0)  weight(t)  weight_new(t) weight(0) dlnA(t) dlnA_new(t) dlnA(0) rms(t) rms(0)"
####################################################
####################################################

#sys.exit()
for mode in modes:
    #Grab all the events for this mode (i.e: 40_100#body_wave) used in FWI
    f_name  = mode.replace("#","_") 
    fp      = open("%s.dat"%(f_name),"w")
    fp.write("%s\n"%(ll_))
    #Grab all the events 
    events  = [it for it in glob("%s/*"%(path)) if (os.path.isdir(it) and EXTENSION(it) == mode)] 
    #print("--- %s seconds ---" % (time.time() - start_time))
    #sys.exit()
    #for pat in [events[0]]:
    for pat in events:
        #Get the Name and the window type of the event
        C      = CheckCategory(pat)
        event_name = C[0]
        win_type   = C[1];
        fwin       = "%s/windows.filter.json"%(pat)
        if os.path.isfile(fwin):
            windows= readwindow(fwin)
        #get the category
        ctype              = win_type.split("#")[0]
        #Observed data file
        fobs_name          = "%s/%s.proc_obsd_%s.h5"%(path_obs,event_name,ctype)
        #Synthetic data file
        fsyn_name          = "%s/%s.proc_synt_%s.h5"%(path_syn,event_name,ctype)
        #print(fobs_name, fsyn_name) 
        #Open the observed and the synthtetic files
        obs_ds = openASDF("%s"%(fobs_name))
        syn_ds = openASDF("%s"%(fsyn_name))
        #####Grab the stations if the correponding JSON FILE#
        stations = sorted(windows.keys())
        #Get the event coordinate
        #########################################
        #Loop over the stations
        for station in stations:
            #try to grab informations for the corresponding station:
            try:
                tag_obs_name = obs_ds.waveforms[station].get_waveform_tags()[0]
                tag_syn_name = syn_ds.waveforms[station].get_waveform_tags()[0]
                #get the Stream of traces
                obs_st       = obs_ds.waveforms[station][tag_obs_name]
                syn_st       = syn_ds.waveforms[station][tag_syn_name]
                #get time serie
                obs_time     = obs_st[0].times()
                syn_time     = syn_st[0].times()
                #get the starttime
                obs_startime = str(obs_st[0].stats.starttime).split('.')[0]
                syn_startime = str(syn_st[0].stats.starttime).split('.')[0]
            except:
                continue
            item      = (event_name, station)
            #get the coordinates of the events and the stations as well as the epicentral distance
            event_lat = INFOS_DICT[item][0]
            event_lon = INFOS_DICT[item][1]
            slat      = INFOS_DICT[item][2]
            slon      = INFOS_DICT[item][3]
            EPIDIST   = INFOS_DICT[item][4]
            #make the id for the observed and the synthetic in order to easily grab the trace for these components
            #trace_obs_syn_id = {tr_obs.id : tr_syn.id for tr_obs, tr_syn in zip(obs_st, syn_st)}
            #Loop over the component on the stations
            for component in windows[station].keys():
                #set the synthetic component
                component_syn   = "%s.S3.MX%s"%(station, component.strip()[-1])
                #Extract the observed traces
                obs_tr, Amax_b  = ExtractTrace(obs_st, component)
                #Extract the synthtetic traces
                syn_tr, Amax_s  = ExtractTrace(syn_st, component_syn)
                #Loop over the windows that belong to the considered component
                for cwind in windows[station][component]:
                    channel_id  = cwind['channel_id']
                    comp    = channel_id.strip()[-1]
                    #Grab the weight of the window
                    weight  = float(cwind["window_weight"]) 
                    #Get the weighted cross-correlation time measurement from the window 
                    dlnA    = float(cwind["dlnA"]) 
                    left_index         = cwind['left_index']
                    right_index        = cwind['right_index']
                    time_inv           = cwind['dt']
                    #get the central time of the window
                    rtive_starttime    = float(cwind["relative_starttime"])
                    rtive_endtime      = float(cwind["relative_endtime"])
                    shift              = float(cwind["cc_shift_in_seconds"])
                    max_cc             = float(cwind["max_cc_value"])
                    min_period         = float(cwind["min_period"])
                    #get the windowing waveform of the observed and synthetic 
                    start_win_d, ended_win_d, obs_select, tobs_select = get_window_from_trace(obs_tr, cwind)
                    start_win_s, ended_win_s, syn_select, tsyn_select = get_window_from_trace(syn_tr, cwind)
                    #Compute the Amplitude pertubation
                    dlnA_0     = dlnA_win(obs_select, syn_select)
                    #compute the corralation between the data and the synthetics in the current model
                    max_cc_0   = correlation(obs_select, syn_select)
                    #Compute the rms for observed and synthetic
                    rms_0      = VR(obs_select, syn_select)
                    #Compute a weight at zero
                    weight_0  = weigth_func(right_index, left_index, time_inv, min_period, max_cc_0) 
                    ############################################
                    tau       = LagCross(obs_select, syn_select)
                    #Pick the commont interval between the observed and the synthetic
                    if(tau <= 0):
                        w_syn = syn_select[1-tau : ]
                        w_obs = obs_select[1 : len(obs_select) + tau]
                    if(tau >= 0):
                        w_syn = syn_select[1 : len(syn_select) - tau]
                        w_obs = obs_select[1+tau : ]
                    #Compute the  new rms and the correlation coefficient after shifting, this should be the given measurement
                    rms_tau   = VR(w_obs, w_syn)
                    #compute the correlation after shifting
                    max_cc_tau= correlation(w_obs, w_syn)
                    #Compute the Amplitude pertubation after shifting
                    dlnA_tau  = dlnA_win(w_obs, w_syn)
                    #compute the weight after shift
                    weight_tau  = weigth_func(right_index, left_index, time_inv, min_period, max_cc_tau) 
                    #Compute the new shift
                    shift_new = tau * time_inv 
                    #####################################
                    #print(max_cc, c_coef_0,tau*0.2, c_coef_0_, c_coef_tau, c_coef_tau_)
                    #sys.exit()
                    ##########################################
                    #ll="{:14s} {:5.2f}  {:5.2f}   {:12s}  {:20s}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}".format( 
                    #        event_name, float(event_lat), float(event_lon), station, channel_id, float(slat), float(slon), float(EPIDIST), relative_starttime, 
                    #        #relative_endtime, value, max_cc, max_cc_0, weight, dlnA, dlnA_0)
                    #        relative_endtime, value, max_cc, max_cc_0, weight, dlnA, dlnA_tau)

                    ll="{:14s} {:5.2f}  {:5.2f}  {:10s}  {:14s} {:8.2f} {:8.2f}  {:8.2f}  {:8.2f} {:8.2f} {:8.2f} {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f} {:8.2f} {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}  {:8.2f}".format(event_name, float(event_lat), float(event_lon), station, channel_id, float(slat), float(slon), float(EPIDIST), rtive_starttime, 
                            rtive_endtime, shift, shift_new, max_cc, max_cc_tau, max_cc_0, weight, weight_tau, weight_0, dlnA, dlnA_tau, dlnA_0, rms_tau, rms_0)
                    ###write into the file #############
                    fp.write('%s\n'%(ll))
