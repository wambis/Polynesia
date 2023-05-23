#!/usr/bin/env python
import argparse
import h5py
import re,yaml
import sys,os
#sys.path.append("/Volumes/wamba/Plot-Models-2021/ucbpy")
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
DTYPE = np.float32
np.set_printoptions(threshold=np.inf)
import obspy
from yaml.loader import SafeLoader
import obspy.signal.filter as fltobspy
from obspy.geodetics.base import gps2dist_azimuth
from obspy import read, read_inventory
######################################
from obspy.taup import TauPyModel
model   = TauPyModel(model="iasp91")
from array import array
from scipy import signal
#import pyflex
import glob
import os,math,json
from pyasdf import ASDFDataSet
from pytomo3d.window.io import get_json_content, WindowEncoder
import pylab as pl
import math
#
# See the website

def get_nbre_comps(windows, station):
    if not  isinstance(windows, dict):
        sys.exit("the argument windows should be dictinary")
    components =  win_param[station].keys()
    ncomp      =  len(components)
    components =  list(components)
    return (ncomp, components)

def get_nbre_window(windows, station, component):
    if not  isinstance(windows, dict):
        sys.exit("the argument windows should be dictinary")
    windows_out = win_param[station][component]
    if isinstance(windows_out, dict) :
        windows_out   = [windows_out]
    nwindow = len(windows_out)
    return (nwindow, windows_out)

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

def get_window_from_trace(trace, window):
    starttime = window['absolute_starttime']
    start_win = window['relative_starttime']
    ended_win = window['relative_endtime']
    #grab the index
    left_indx = window['left_index']
    right_indx= window['right_index']
    center_indx= window['center_index']
    #Grab the observed and the synthetics in corresponding window#########
    data      = trace.data[left_indx : right_indx]
    time      = trace.times()[left_indx    : right_indx]
    #tc        = start_win + (ended_win - start_win)/2.0
    tc        = start_win + (ended_win - start_win) *(2.0/3.0)
    return(start_win,tc, ended_win, data,time)

def ExtractTrace(st, id_channel):
    #print(st, id_channel) 
    tr      = st.select(id=id_channel)[0]
    #Convert data to micro-meter
    tr.data = tr.data * 1000000
    Amax    = max(tr.data)
    #get the 10% of the Amax
    Amax_per_cent= (Amax * 10.0)/100.0
    #get the final Amax
    Amax_F  = math.ceil(Amax + Amax_per_cent)
    return (tr, Amax_F)

def LagCross(d, s):
     #Compute the cross correlation
    #corr = signal.correlate(obs_select - np.mean(obs_select), syn_select - np.mean(syn_select), mode="full")
    corr = signal.correlate(d, s, mode="full")
    #get the corresponding tau (i.e. lags)
    lags = signal.correlation_lags(len(d), len(s), mode="full")
    corr /= np.max(corr)
    #grab sifted positions, get the corresponding tau (i.e. lags)
    correlation_max_index = np.argmax(corr)
    #get the lag corresponding to the maximum value of the cross-correlation
    tau  = lags[correlation_max_index]
    return tau


def ExtractTrace_NEW(st, comp):
    try:
        tr      = st.select(component=comp)[0]
    except:
        return None
    #Convert data to micro-meter
    #tr.data = tr.data * 1000000
    tr.data = tr.data 
    Amax    = max(tr.data)
    #get the 10% of the Amax
    Amax_per_cent= (Amax * 10.0)/100.0
    #get the final Amax
    Amax_F  = math.ceil(Amax + Amax_per_cent)
    return (tr, Amax_F)

def Average(lst):
    s=0.0
    for it in lst:
        s = s+float(it) 
    M = s/len(lst)

    return M

#Define the correlation funuction
def correlation(d, s):
       proc1   = (d -np.mean(d)) * (s -np.mean(s))
       proc2   = (d -np.mean(d)) * (d -np.mean(d))
       proc3   = (s -np.mean(s)) * (s -np.mean(s))
       c       = proc2.sum() * proc3.sum()
       corr    = proc1.sum()/np.sqrt(c)
       #corr    = "%.2f"%(corr)
       return corr

#def correlation(d, s):
#       dxs   = (d) * (s)
#       corr    = dxs.sum()/np.sqrt((d ** 2).sum() * (s ** 2).sum())
#       return corr
#def Xcorr_win(d, s):
#        cc = np.correlate(d, s, mode="full")
#        time_shift = cc.argmax() - len(d) + 1
#        # Normalized cross correlation.
#        max_cc_value = cc.max() / np.sqrt((s ** 2).sum() * (d ** 2).sum())
#        return max_cc_value, time_shift

def Xcorr_win(d, s):
        cc = np.correlate(d, s, mode="full")
        time_shift = cc.argmax() - len(d) + 1
        # Normalized cross correlation.
        max_cc_value = cc.max() / np.sqrt((s ** 2).sum() * (d ** 2).sum())
        return max_cc_value


def dlnA_win(d, s):
    return 0.5 * np.log(np.sum(d ** 2) / np.sum(s ** 2))


def VR(data, syn):
    diff2  = (data - syn)**2
    syn_   = syn    - np.mean(syn)
    data_  = data   - np.mean(data)
    diff2  = (data_ - syn_)**2
    data2  = (data_)**2
    V      = diff2.sum()/data2.sum()
    rms    = np.sqrt(V)
    return rms




#Set the parameter:
#set to true if you need the table output
TABLE      = True
#set to true if you need to plot non-picked seismogram
NON_PICKED = False
#NON_PICKED = True





###############################################
with open('configPlot.yaml') as Fym:
    #This allow to avoid the _io.TextIOWrapper error
    Fig_params = yaml.load(Fym, Loader=SafeLoader)

#Grab the path of the data
path_win  = Fig_params['path_wid']
#Grab the regional station file and events file
path_obs  = Fig_params['path_obs']
path_syn  = Fig_params['path_syn']
#Grab the event id
event_id  = Fig_params['event_id']
#Grab windowns files
#win_file = Fig_params['win_file']
#Grab the net.stn
net_stns  = Fig_params['net_stn']
#Grab the frequency band
net_stns  = Fig_params['net_stn']
#Grab the frequency
periodband= Fig_params['periodband']
#Axis-limit
x_lm_left = Fig_params['x_lm_left']
x_lm_right= Fig_params['x_lm_right']
#get the color list in order to plot windows with different color
fcolors   = Fig_params['fcolors']



#set the standard components with index
#CompDict      = {'Z':0,'R':1,'T':2}
CompDict      = {0:'Z', 1:'R', 2:'T'}
#Selec_Evens   = [l.split()[0] for l in Fs]  
Selec_Evens   = [event_id]  
#OBSs file name
fobs      = open("REGIONAL-EVENTS-LAND-STATIONS-SELECTED-in-Chunk.REG", "r")
LINES_OBS = fobs.readlines()

for l in LINES_OBS:
    category           = l.split()[0] 
    event_name         = l.split()[1] 
    station            = l.split()[2] 
    channel            = l.split()[3] 
    left_index         = int(l.split()[4]) 
    right_index        = int(l.split()[5]) 
    cc_shift_in_seconds= float(l.split()[-1]) 
    #Window file name
    fwin_name          = "%s/%s.%s/%s"%(path_win, event_name, category,"windows.filter.json")
    ctype              = category.split("#")[0]
    #Observed data file
    fobs_name          = "%s/%s.proc_obsd_%s.h5"%(path_obs,event_name,ctype)
    #Synthetic data file
    fsyn_name          = "%s/%s.proc_synt_%s.h5"%(path_syn,event_name,ctype)
    #fig_name           ="%s.%s.%s"%(event_name,station,channel)
    fig_name           ="%s.%s.%s.%s"%(event_name,station,channel[-1], ctype.replace('_','-'))
    BP                 = ctype
    if(os.path.isfile(fwin_name) and os.path.isfile(fobs_name) and os.path.isfile(fsyn_name)):
        #if(cc_shift_in_seconds < -7.0 or cc_shift_in_seconds > 7.0):
        #if(abs(cc_shift_in_seconds) > 5.0 and event_name in Selec_Evens):
        #if(cc_shift_in_seconds == 0.00 and event_name in Selec_Evens):
        if(event_name in Selec_Evens and station in net_stns and BP==periodband):
                #load corresponding parameter file
                #print("window selection config parameters:")
                f             = open(fwin_name)
                #grab the windows
                win_param     = json.load(f)
                f.close()
                #print(win_param) 
                #Grab the dictionary corresponding and channel
                window_dict        = win_param[station][channel][0]
                #grab all the stations
                #stations      = win_param.keys()
                relative_starttime = window_dict['relative_starttime']
                relative_endtime   = window_dict['relative_endtime']
                left_index         = window_dict['left_index']
                right_index        = window_dict['right_index']
                time_inv           = window_dict['dt']
                cc_shift_in_seconds= float(window_dict['cc_shift_in_seconds'])
                max_cc_value       = window_dict['max_cc_value']
                channel_id         = window_dict['channel_id']
                channel_id_2       = window_dict['channel_id_2']
                #stations           = [station]
                #try to grab the data and the synthetics
                obs_ds = openASDF("%s"%(fobs_name)) 
                syn_ds = openASDF("%s"%(fsyn_name))
                #Get event and event latitude and longitude
                event         = obs_ds.events[0]
                origin        = event.preferred_origin() or event.origins[0]
                event_lat     = origin.latitude
                event_lon     = origin.longitude
                resouce_id    = origin.resource_id
                event_id      = str(resouce_id).split('/')[-2]
                Mw            = '%.1f'%(event["magnitudes"][0]["mag"]) 
                #################################
                depth         = '%.1f'%(origin.depth/1000.0)
                depth_km      = '%s km'%(str(depth)) 
                #Open a table for all the stations 
                TAB_STN       = { }
                #Loop over the stations 
                for station in [station]:
                    basename   = station
                    stn        = station.split('.')[-1]
                    #try to grab informations:
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
                    #channel     = obs_st[0].stats.channel
                    network     = obs_st[0].stats.network
                    station_    = obs_st[0].stats.station
                    #Grab Statation latitude, longitude, depth and elevation
                    try:
                        slat    = obs_ds.waveforms[station].coordinates["latitude"]
                        slon    = obs_ds.waveforms[station].coordinates["longitude"]
                        Elv     = obs_ds.waveforms[station].coordinates["elevation_in_m"]
                        #depth  = obs_ds.waveforms[station].coordinates["local_depth_in_m"]
                    except:
                        continue
                    #compute de azimuth and the epicenter distance
                    D, az, baz  = gps2dist_azimuth(event_lat, event_lon, slat,slon)
                    D_km        = D/1000.
                    D_deg       ='%.1f'%(D_km/111.)
                    #Compute travel-tme
                    arrivals    = model.get_travel_times(source_depth_in_km =float(depth), distance_in_degree = float(D_deg))
                    #Dictionary for each phase and their travel time
                    Phases      = {p.name: p.time for p in arrivals}
                    ####################################################
                    #grab component from windows files
                    ncom,components= get_nbre_comps(win_param, station)
                    #Get the desired component from the json file 
                    #Create a figure
                    fig, axs = plt.subplots(2,figsize=(8, 8))
                    #Grid the figure
                    #plt.grid(visible=True, axis='both')
                    #space between the figure
                    fig.tight_layout()
                    #print(len(axs), list(axs)) 
                    #sys.exit()
                    #set the figure name
                    #figname  = station
                    figname  = '%s_%s'%(event_id, station)
                    sinfo    = '%s.%s'%(network,stn)
                    #Loop over the axis 
                    #for  icm in CompDict:
                    for  icm in range(len(axs)):
                            print(icm) 
                            comp        = channel_id[-1] 
                            #Get the synthetic component
                            #try to extract the componenent from the observed and synthetics traces
                            station_cmp = '%s.%s'%(station,comp)
                            #signed the station into the dictionary
                            TAB_STN[station_cmp]= {}
                            #Get all the windows from the components
                            nwin, windows   = get_nbre_window(win_param, station, channel_id)
                            obs_tr, Amax_b  = ExtractTrace(obs_st, channel_id) 
                            syn_tr, Amax_s  = ExtractTrace(syn_st, channel_id_2) 
                            Amax_F          = max(Amax_b, Amax_s)
                            Amax            = max(max(obs_tr.data), max(syn_tr.data))
                            Amax            = Amax - Amax * 9.0/1000.
                            #Amax           = Amax_F - Amax_F * 15.0/1000.
                            #Select the window desired
                            window_select   =[win for win in windows if(win['right_index']==right_index and win['left_index']== left_index and win['channel_id']==channel_id)]
                            window_select   = window_select[0]
                            #print(window_select) 
                            axs[icm].set_ylim(-abs(Amax_F), abs(Amax_F))
                            #sys.exit()
                            if(icm==0):
                                #axs[icm].set_title(sinfo+ '  '+ r'$\Delta$='+str(D_deg)+r'$^{\circ}$'+ ' Mw=%s  depth=%s km    %s'%(str(Mw), 
                                #        str(depth),event_name) +'  '+ 'BP='+BP.replace('_','-')+'s', loc = 'center',fontsize=10)

                                axs[icm].set_title(sinfo+ '  '+ r'$\Delta$='+str(D_deg)+r'$^{\circ}$'+ ' Mw=%s  depth=%s km'%(str(Mw),
                                        str(depth)) +'  '+ 'BP='+BP.replace('_','-')+' s', loc = 'center',fontsize=12, pad=0.5)
                                #Set phases on the plot
                                if('P' in Phases and 'S' in Phases):
                                        axs[icm].axvline(x=float(Phases['P']), color='k',linewidth=0.9, linestyle='dashdot', alpha=0.7)
                                        axs[icm].axvline(x=float(Phases['S']), color='k',linewidth=0.9, linestyle='dashdot', alpha=0.7)
                                if('PP' in Phases and 'SS' in Phases and float(D_deg) > 100.0):
                                        axs[icm].axvline(x=float(Phases['PP']), color='k',linewidth=0.9, linestyle='dashdot', alpha=0.7)
                                        axs[icm].axvline(x=float(Phases['SS']), color='k',linewidth=0.9, linestyle='dashdot', alpha=0.7)
                                #set axis limitation#########
                                ####axs[icm].set_ylim(-abs(Amax_F), abs(Amax_F))
                                #axs[icm].set_xlim(min(abs(obs_tr.times())), max(abs(obs_tr.times())))
                                axs[icm].set_xlim(x_lm_left, x_lm_right)
                                #get the minimum and the maximum
                                ymin, ymax= axs[icm].get_ylim() 
                                #Open an empty dictionary
                                comp_dict = { }
                                #Loop over the windows on a particular component
                                for iw, window in zip(range(nwin), windows):
                                        #print(icm)
                                        #New variable for correlation and variance
                                        corr_ = 'corr_%d'%(iw)
                                        var_  = 'var_%d'%(iw)
                                        #get the data window
                                        start_win_d,tc_d, ended_win_d, data_d, t_d = get_window_from_trace(obs_tr, window)
                                        #get the synthetic window
                                        start_win_s,tc_s, ended_win_s, data_s, t_s = get_window_from_trace(syn_tr, window)
                                        ##get the minimum and the maximum
                                        #ymin, ymax= axs[icm].get_ylim() 
                                        #compute the corralation between the data and the synthetics
                                        corr_0      = correlation(data_d, data_s)
                                        VAR_0       = VR(data_d, data_s)
                                        ##Put parameter into the dictionary
                                        comp_dict[corr_] = corr_0
                                        comp_dict[var_]  = VAR_0
                                        #plot the extracted window
                                        axs[icm].axvline(x=start_win_d, color='k',linewidth=1.5, linestyle='--', alpha=0.3)
                                        axs[icm].axvline(x=ended_win_d, color='k',linewidth=1.5, linestyle='--', alpha=0.3)
                                        #axs[icm].axvline(x=start_win_d, color=fcolors[iw],linewidth=1.5, linestyle='--', alpha=0.9)
                                        #axs[icm].axvline(x=ended_win_d, color=fcolors[iw],linewidth=1.5, linestyle='--', alpha=0.9)
                                        #print(iw)
                                        axs[icm].axvspan(start_win_d, ended_win_d, facecolor=fcolors[iw], alpha = 0.2)
                                        axs[icm].plot(t_d, data_d, 'r',linewidth=1.4)
                                        axs[icm].plot(t_s, data_s, 'b',linewidth=1.4)
                                        #Write correlation and the variance 
                                        axs[icm].text(tc_d,  ymax/2.0, '%.2f'%(float(corr_0)), horizontalalignment='center',verticalalignment='center',alpha=0.9,fontsize=4)
                                        axs[icm].text(tc_d, ymin/2.0,  '%.2f'%(float(VAR_0)), horizontalalignment='center',verticalalignment='center',alpha=0.9, fontsize=4)
                                #Plot the entire seismogram
                                axs_c    = axs[icm].twinx()
                                axs[icm].plot(obs_tr.times(), obs_tr.data, 'r', linewidth=1.2, alpha = 0.6,label='Observed')
                                axs[icm].plot(syn_tr.times(), syn_tr.data, 'b', linewidth=1.2, alpha = 0.6, label='GLADM25')
                                #axs[icm].legend(loc='upper right',prop={'size':5 })
                                axs[icm].legend(loc='upper right',prop={'size':10 })
                                #axs[icm].set_ylabel(r'displacement (${\mu}$m)')
                                #axs[icm].set_ylabel(r'Displacement (${\mu}$m)', labelpad=1.0)
                                axs[icm].set_ylabel(r'Displacement (${\mu}$m)', labelpad= -1.5, fontsize=8)
                                axs[icm].yaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
                                #axs[icm].set_xlim(x_lm_left, x_lm_right)
                                #print('point %.2f   %.2f'%(x_lm_left, x_lm_right))
                                if('P' in Phases and 'S' in Phases):
                                    axs[icm].text(float(Phases['P']) -24 , ymin/1.2,  'P', horizontalalignment='center',verticalalignment='center',alpha=0.9, fontsize=8)
                                    axs[icm].text(float(Phases['S']) -24, ymin/1.2,  'S', horizontalalignment='center',verticalalignment='center',alpha=0.9, fontsize=8)
                                if('PP' in Phases and 'SS' in Phases and float(D_deg) > 100.0):
                                    axs[icm].text(float(Phases['PP']) -24, ymin/1.2,  'PP', horizontalalignment='center',verticalalignment='center',alpha=0.9, fontsize=8)
                                    axs[icm].text(float(Phases['SS']) -24, ymin/1.2,  'SS', horizontalalignment='center',verticalalignment='center',alpha=0.9, fontsize=8)
                                #set the component on the axis
                                axs_c.set_ylabel(comp, color = 'k', rotation='vertical')
                                axs_c.set_yticks([])
                
                            if(icm==1):
                                #fig.tight_layout()
                                #get the data from seleceted window
                                start_win_d,tc_d, ended_win_d, data_select, td_select = get_window_from_trace(obs_tr, window_select) 
                                #get the synthetic from selected window
                                start_win_s,tc_s, ended_win_s, syn_select, ts_select = get_window_from_trace(syn_tr, window_select)
                                #Compute the cross correlation
                                tau   = LagCross(data_select, syn_select)
                                #get the starttime and the endtime of the selected window
                                relative_starttime = window_select['relative_starttime']
                                relative_endtime   = window_select['relative_endtime']
                                #Compute a new correlation
                                max_cc_0           = correlation(data_select, syn_select)
                                rms_0              = VR(data_select,syn_select)
                                #compute the amplitude pertubation
                                dlnA_0             = dlnA_win(data_select, syn_select)
                                ##############################
                                #Pick the commont interval between the observed and the synthetic 
                                if(tau <= 0):
                                    w_syn = syn_select[1-tau : ]
                                    w_obs = data_select[1 : len(data_select) + tau]
                                    t_obs = ts_select[1-tau  : ]
                                if(tau >= 0):
                                    w_syn = syn_select[1 : len(syn_select) - tau]
                                    w_obs = data_select[1+tau : ]
                                #Compute the  new rms and the correlation coefficient 
                                rms_tau          = VR(w_obs, w_syn)
                                max_cc_tau       = correlation(w_obs, w_syn)
                                #compute the amplitude pertubation tau
                                dlnA_tau         = dlnA_win(w_obs, w_syn)
                                #Shifted time of the synthetics
                                #compute the multiplicator factor after shifting
                                M_coef           = np.exp(dlnA_tau) 
                                #######################################################
                                Amax             = max(max(data_select), max(syn_select))
                                Amax_F           = Amax - Amax * 9.0/1000.
                                #print(Amax_F)
                                #Plot the window (i.e. a piece of seismogram) from the entire seismogram
                                #axs[icm].plot(ts_select, 200*data_select, 'r',linewidth=2,label='observed')
                                axs[icm].plot(ts_select, data_select, 'r',linewidth=2.3,label='Observed')
                                axs[icm].plot(ts_select, syn_select, 'b',linewidth=2.3, label='M25 synthetic')
                                ##Grid the plot
                                #axs[icm].grid()
                                #Now mutilplier the synthetic with the multiplicator factor after shifting 
                                syn_select        = M_coef * syn_select
                                if(tau < 0.0):
                                    axs[icm].plot(ts_select + tau*time_inv, syn_select, 'k--',linewidth=2., label = 'M25 advanced by %.1f s\nand scaled by %.1f'%(abs(tau*time_inv), M_coef))
                                elif(tau > 0.0):
                                    axs[icm].plot(ts_select + tau*time_inv, syn_select, 'k--',linewidth=2., label = 'M25 delayed by %.1f s\nand scaled by %.1f '%(abs(tau*time_inv), M_coef))
                                elif(tau == 0.0):
                                    axs[icm].plot(ts_select + tau*time_inv, syn_select, 'k--',linewidth=2., label = 'M25 nonshifted  %.1f s\nand scaled by %.1f'%(abs(tau*time_inv),M_coef))
                                #write the title on the figure
                                axs[icm].set_title('event-id: %s'%(event_id), loc = 'right',fontsize=9, pad=.05)
                                #Text to write on the figure
                                textstr   = "rms($\u03C4$)=%.2f\n\n rms(0)=%.2f"%(float(rms_tau),float(rms_0))
                                #textstr   = "rms($\delta t$)=%.2f\n\n rms(0)=%.2f"%(float(rms_tau),float(rms_0))
                                #axs[icm].set_title(sinfo+ '  '+ r'$\Delta$='+str(D_deg)+r'$^{\circ}$'+ ' Mw=%s  depth=%s km    %s'%(str(Mw), 
                                ################
                                #print(max_cc_tau, max_cc_0, tau*time_inv, fig_name)
                                #textstr_2 = "C($\delta t$)=%.2f\n\n C(0)=%.2f"%(float(max_cc_tau), float(max_cc_0))
                                textstr_2 = "C($\u03C4$)=%.2f\n\n C(0)=%.2f"%(float(max_cc_tau), float(max_cc_0))
                                #Get text for the amplitude
                                textstr_dlnA = "dlnA($\u03C4$)=%.2f\n\n dlnA(0)=%.2f"%(float(dlnA_tau), float(dlnA_0))
                                #textstr_dlnA = "dlnA($\delta t$)=%.2f\n\n dlnA(0)=%.2f"%(float(dlnA_tau), float(dlnA_0))
                                #put the labels on the axis
                                axs[icm].set_xlabel(r'Time (s) since '+obs_startime, labelpad=2.0, fontsize=9)
                                #write text on the figure
                                #axs[icm].text(min(ts_select) , abs(Amax_F)/2.0, textstr, fontsize=8)
                                axs[icm].text(min(ts_select) , abs(Amax_F), textstr, fontsize=12)
                                #write for the Cross-correlation
                                axs[icm].text(min(ts_select) , -abs(Amax_F)*1.5, textstr_2, fontsize=12)
                                #write for the amplitude
                                axs[icm].text(max(ts_select) * 0.95 , -abs(Amax_F)*1.5, textstr_dlnA, fontsize=12)
                                #axs[icm].text(max(ts_select) * 0.9 , -abs(Amax_F)*1.5, textstr_dlnA, fontsize=10)
                                #set axis limit
                                #axs[icm].set_ylim(-abs(Amax_F), abs(Amax_F))
                                axs[icm].set_ylim(-(abs(Amax_F) +abs(Amax_F)/1.5), (abs(Amax_F) +abs(Amax_F)/1.5))
                                #ax.text(0.5 , 0.5, textstr)
                                #set the component on the axis
                                axs[icm].set_yticks([])
                                #axs[icm].set_ylabel(r'displacement (${\mu}$m)', labelpad=2.0 )
                                #axs[icm].set_ylabel(r'displacement (${\mu}$m)')
                                #axs[icm].grid()
                                axs[icm].legend(prop={'size': 8}) 
                
                #set the figure name
                fname = fig_name+'.png'
                plt.savefig(fname,dpi=800)
                plt.close(fig)
