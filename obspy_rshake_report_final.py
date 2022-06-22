# -*- coding: utf-8 -*-
"""
Created on Sun May 22 10:53:50 2022

@author: Christopher Jeff Sanchez, Jeremy Alan Hilado
"""

# This python script creates a .csv report on the peak acceleration, velocity, 
# and displacement of an earthquake event, as well as its intensity. Currently 
# it can only process .mseed files from RSHAKE instruments. The overall 
# process is heavily based on the PRISM script used by USGS.

# To run, just place this .py file into the folder containing the .mseed files.
# Manually specify the name of the station and the network. You can edit them
# by changing the corresponding variables (station and network) at the bottom
# of this script.
# You can use any python IDE or run through a python terminal to 
# run the script.

# The names of the files in the parent folder must contain any channel
# from RSHAKE (EHZ, ENE, ENN, ENZ), otherwise it won't run.



from obspy import read_inventory
from obspy.core import read
from obspy.clients.fdsn import Client as RS_Client
from obspy.core import UTCDateTime, Stream
from obspy.core.inventory import Inventory, Network
#from obspy.signal.invsim import estimate_magnitude
from obspy.signal.trigger import recursive_sta_lta, trigger_onset, plot_trigger
from obspy.signal.util import _npts2nfft
from scipy import integrate
from pathlib import Path
import matplotlib.pyplot as plt
import os
import numpy as np
import warnings
import re
import pandas as pd
import eqAlarm
from eqAlarm import intensityConvert
import threading


#STA/LTA constants
INV_DIR = "inventories"
sta_window = 2
lta_window = 10
thresh_on = 1.5
thresh_off = 1


def convert_acc_to_vel_trace(tr): ### Converts accelaration waveform to velocity
    if tr.stats.units == "ACC":
        tr = improved_integration(tr)
        tr.stats.units = "VEL"
    elif tr.stats.units == "VEL":
        tr.stats.units = "VEL"
    else:
        print("Can't convert", tr.stats.units, "to VEL.")

    return tr

def convert_vel_to_dis_trace(tr): ### Converts velocity waveform to displacement
    if tr.stats.units == "VEL":
        tr = improved_integration(tr)
        tr.stats.units = "DIS"
    elif tr.stats.units == "DIS":
        tr.stats.units = "DIS"
    else:
        print("Can't convert", tr.stats.units, "to DIS.")

    return tr

def convert_counts_to_metric_trace(tr, metric_units, limit=True): ### Converts counts to metric
    if tr.stats.units == "COUNTS":
        tr = tr.copy()
        freq = tr.stats.sampling_rate
        #print(estimate_magnitude(inv[0][0][43].response ,max(abs(tr.data)), 0.1, 20))
        # Note: The following are applied before deconvolution
        # a bandpass filter
        # a quarter cosine taper, tapering 0.05 length of data
        lowcut = get_low_corner_freq(tr, plot=False)
        if limit:
            lowcut = lowcut if lowcut <= 0.5 else 0.5
        #print("used lowcut Hz:", lowcut)
        tr.remove_response(pre_filt=[lowcut, lowcut*1.1, 0.5*freq*0.9, 0.5*freq],
        #tr.remove_response(pre_filt=[0.1, 0.5, 49, 50],
                           output=metric_units, water_level=4.5, taper=True, taper_fraction=0.1)
        tr.stats.units = metric_units
    return tr, lowcut

def get_weighted_mean(x, weights):
    sum_weight = integrate.trapz(weights)
    weighted_x = x*weights

    return sum(weighted_x)/sum_weight

def get_weighted_variance(x, weights):
    sum_weight = integrate.trapz(weights)
    weighted_mean = sum(x*weights)/sum_weight
    var = sum(weights*(x - weighted_mean)*(x - weighted_mean))/sum_weight

    return var

def get_low_corner_freq(tr, event_onset = None, low_power_thresh = 0.015, use_end=True, plot=False, verbose=False):
    ###Gets low corner frequency based from baseline waveform prior or after event
    #use_end -> get baseline from after event if True, prior to event if False
    
    global event_type
    tr = tr.copy()
    if event_onset == None:
        # get event onset
        freq = tr.stats.sampling_rate
        sta_lta = recursive_sta_lta(tr, nsta=int(sta_window*freq), nlta=int(lta_window*freq))
        pick_pairs_ind = trigger_onset(sta_lta, thresh_on, thresh_off)
        event_onset = pick_pairs_ind[0][0]
        #print("event onset", str(event_onset))

    if use_end:
        noise = tr.data[-1*int(event_onset-1*tr.stats.sampling_rate):] # remove one second worth of samples
        #noise = tr.data[event_onset:event_onset+10000] # use signal
    else: # get noise from beginning
        noise = tr.data[:int(event_onset-1*tr.stats.sampling_rate)] # remove one second worth of samples
    noise = noise - np.mean(noise) #remove DC
    noise_resp = abs(np.fft.rfft(noise, n=_npts2nfft(len(noise)))) # get distance as magnitude
    noise_resp_x = np.fft.rfftfreq(_npts2nfft(len(noise)), 1/tr.stats.sampling_rate)
    nyquist_ind = len(noise_resp)//2
    #weighted_mean = get_weighted_mean(noise_resp_x[:nyquist_ind], noise_resp[:nyquist_ind])
    #weighted_var = get_weighted_variance(noise_resp_x[:nyquist_ind], noise_resp[:nyquist_ind])
    #print("Central Frequency:", weighted_mean, "Hz")
    #print("Average distance from center:", np.sqrt(weighted_var) , "Hz")
    noise_integral = integrate.cumtrapz(noise_resp[:nyquist_ind], noise_resp_x[:nyquist_ind], initial=0) # only up to nyquist freq
    noise_integral = noise_integral/noise_integral[-1] # normalize w max
    lowcut_ind = int((noise_integral > low_power_thresh).nonzero()[0][0])
    lowcut = noise_resp_x[lowcut_ind]
    if verbose:
        print("power% centered at lowcut:", noise_integral[lowcut_ind-1:lowcut_ind+2])
        print("lowcut_ind:",lowcut_ind)
        print("lowcut Hz:",lowcut)

    # get first occurence of greater than 1/2 threshold

    signal = tr.detrend("linear").data
    signal_resp = abs(np.fft.rfft(signal, n=_npts2nfft(len(signal)))) # get distance as magnitude
    signal_resp_x = np.fft.rfftfreq(_npts2nfft(len(signal)), 1/tr.stats.sampling_rate)
    final_signal = tr.copy()
    final_signal.remove_response(pre_filt=[lowcut, lowcut*1.1, 49, 50],
                       output="ACC", water_level=4.5, taper=True, taper_fraction=0.1)

    if plot:
        fig, axs = plt.subplots(nrows=3, ncols=2, constrained_layout=True)
        fig.set_size_inches(w=12, h=6)

        fig.suptitle("Event "+str(event_type)+" Noise Removal")

        # NOISE PLOT COLUMN
        axs[0][0].set_ylabel("Noise\nCounts")
        axs[0][0].plot(noise)
        axs[1][0].set_ylabel("Noise Response\nCounts")
        axs[1][0].plot(noise_resp_x,noise_resp)
        axs[2][0].set_ylabel("Noise Response\nNormalized Integral")
        axs[2][0].plot(noise_resp_x[:nyquist_ind],noise_integral)

        # SIGNAL PLOT COLUMN
        axs[0][1].set_ylabel("Signal\nCounts")
        axs[0][1].plot(tr.data)
        axs[1][1].set_ylabel("Signal Response\nCounts")
        axs[1][1].plot(signal_resp_x,signal_resp)
        axs[2][1].set_ylabel("Final Signal \nCounts")
        axs[2][1].plot(final_signal)

        for row in axs:
            for ax in row:
                ax.grid()
        plt.show()

    return lowcut

def improved_integration(tr):
    tr = tr.copy()
    tr.detrend("demean") # make sure no constant that will become linear function
    tr.integrate(method="cumtrapz")
    tr.detrend("linear") # (mx+b, ie the leakage due to cumtrapz)

    return tr

def get_inventory(inv_path, network, station, client_name="RASPISHAKE"):
    ### Gets inventory response file for RSHAKE instrument
    #print(inv_path)
    # create copy of latest inv, remove date so it can be used in remove_response
    if os.path.exists(inv_path):
        #print("reading inv")
        inv = read_inventory(inv_path)
    else:
        print("downloading inv")
        rs_client = RS_Client(client_name)
        inv = rs_client.get_stations(network="AM", station=station, level="RESP")
        latest_station_response = (inv[-1][-1]).copy()
        latest_station_response.start_date = None
        latest_station_response.end_date = None
        print(latest_station_response)
        for cha in latest_station_response:
            cha.start_date=None
            cha.end_date=None
        inv = Inventory(networks=[ \
                Network(code=network, stations=[latest_station_response])])
        Path(os.path.dirname(inv_path)).mkdir(parents=True, exist_ok=True)
        inv.write(inv_path, format="STATIONXML")
    return inv

def getIntensity(g): ### Gets intensity from PGA
    intensity_scale = {
        (0,0.00170): 'I',
        (0.00170,0.01400): 'II-III',
        (0.01400,0.03900): 'IV',
        (0.03900,0.09200): 'V',
        (0.09200,0.18000): 'VI',
        (0.18000,0.34000): 'VII',
        (0.34000,0.65000): 'VIII',
        (0.65000,1.24000): 'IX',
        (1.24000,5): 'X+'
    }
    for i in intensity_scale:
        if i[0] < g < i[1]:
            intensity = intensity_scale[i]
    return intensity

def obspyProcess(mseed_file, network, station, channel, plot_results = False, invOnline = True):

    # Obspy processing
    st = read(mseed_file, format="MSEED")
    metric_units = "ACC"

    inv_path = os.path.join(os.getcwd(), INV_DIR, network+"_"+station+"_"+channel+".xml")
    inv = get_inventory(inv_path, network, station)

    counts = st.select(channel=channel)[0]
    #freq = counts.stats.sampling_rate
    #plot_trigger(counts, recursive_sta_lta(counts, nsta=int(sta_window*freq), nlta=int(lta_window*freq)), 1.5, 1)
    counts.attach_response(inv)
    counts.stats.units = "COUNTS"
    counts_mean = counts.data.mean()
    counts_peak = max(abs(counts.data))
    acc, lowcut = convert_counts_to_metric_trace(counts, metric_units)
    #return
    acc_mean = acc.data.mean()
    acc_peak = max(abs(acc.data))
    vel = convert_acc_to_vel_trace(acc)
    vel_mean = vel.data.mean()
    vel_peak = max(abs(vel.data))
    dis = convert_vel_to_dis_trace(vel)
    dis_mean = dis.data.mean()
    dis_peak = max(abs(dis.data))
    
    # Plot results
    if plot_results:
        fig, axs = plt.subplots(nrows=4, ncols=3, constrained_layout=True)
        fig.set_size_inches(w=12, h=6)

        fig.suptitle(channel+" Channel in Metric Units")

        # OBSPY PLOT COLUMN
        axs[0][0].set_ylabel("OBSPY\nCounts")
        axs[0][0].plot(counts.data)
        x_text0 = len(counts.data)*0.65
        y_text0 = max(counts.data)*0.7
        axs[0][0].text(x_text0,y_text0,
                    "Mean = {:+.6f}counts\nPeak = {:.6f}counts".format(counts_mean,counts_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        _, _, counts_ymin, counts_ymax = axs[0][0].axis()

        axs[1][0].set_ylabel("OBSPY\nAcceleration")
        axs[1][0].plot(acc.data)
        x_text0 = len(acc.data)*0.65
        y_text0 = max(acc.data)*0.7
        axs[1][0].text(x_text0,y_text0,
                    "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(acc_mean,acc_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        _, _, acc_ymin, acc_ymax = axs[1][0].axis()

        axs[2][0].set_ylabel("OBSPY\nVelocity")
        axs[2][0].plot(vel.data)
        x_text1 = len(vel.data)*0.65
        y_text1 = max(vel.data)*0.7
        axs[2][0].text(x_text1,y_text1,
                    "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(vel_mean,vel_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        _, _, vel_ymin, vel_ymax = axs[2][0].axis()

        axs[3][0].set_ylabel("OBSPY\nDisplacement")
        axs[3][0].plot(dis.data)
        x_text2 = len(dis.data)*0.65
        y_text2 = max(dis.data)*0.7
        axs[3][0].text(x_text2,y_text2,
                    "Mean = {:+.6f}m\nPeak = {:.6f}m".format(dis_mean,dis_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        _, _, dis_ymin, dis_ymax = axs[3][0].axis()
        
        for row in axs:
            for ax in row:
                ax.grid()
        plt.show()
        
    intensity = (getIntensity(acc_peak/9.81))
    processing_info = [lowcut]
    instrument_info = [network, station, channel]
    obspy_peaks = [acc_peak, vel_peak, dis_peak]
    
    return (processing_info, instrument_info, obspy_peaks, intensity)

def threadMain(intensity,displacement):
    eqAlarm.alarm(intensity,displacement)
        
    
    
if __name__ == "__main__":
    
    # Set path
    parent_dir = Path(__file__).parent.resolve()
    #parent_dir = r"C:\Users\benit\Documents\Python Scripts\cosmos_downloads_v2"

    columns=["Station", "OBS-PeakAcc (m/s^2)", "OBS-PeakVel (m/s)", 
             "OBS-PeakDis (m)","Intensity"]
    rows = []
        
    doDownloadMseed =   False   #Set True if you still need to dload msd files
    doAlert         =   True    #Set True if you want to activate alarm
    
    #Set values
    
    network = "AM" ### <<< replace this with your network name ###
    station = 'R2A83' ### <<< replace this with your station name ###
    channel = ["ENE", "ENN", "ENZ"] ### contains all channels for RSHAKE 4D
    intThresh = 2 ### set intensity threshold for alarm
    
    # set data start/end times for dloading mseed (if mseed is not yet dloaded)
    start = UTCDateTime(2022, 5, 21, 21, 51, 0) # (YYYY, m, d, H, M, S)
    end = UTCDateTime(2022, 5, 21, 21, 55, 0) # (YYYY, m, d, H, M, S)

    
    if doDownloadMseed:
        rs = RS_Client('RASPISHAKE')
        # get waveforms and put them all into one Stream
        stream = Stream()
        for ch in channel:
            trace = rs.get_waveforms('AM', station, '00', ch, start, end)
            stream += trace
            
        #write mseed files for each channel
        for tr in stream: 
            tr.write(tr.id + ".mseed", format="MSEED") 
    
    
    mseedFiles = [file for file in os.listdir(parent_dir) 
                  if file.endswith(".mseed") and station in file]
    
    
    #Obspy processing
    for cha in channel:
        try:
            mseed_file = [file for file in mseedFiles if cha in file][0]
            processing_info, instrument_info, obspy_peaks, intensity = obspyProcess(os.path.join(parent_dir,mseed_file), network, station, cha)
            
            
            # collect output
            lowcut = processing_info[0]
            (acc, vel, dis) = obspy_peaks
            rows.append([mseed_file, acc, vel, dis, intensity])
        except:
            pass
    
    rows.append([""]*5)
        
    #saving report as csv
    df = pd.DataFrame(columns=columns,data=rows)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df)
    df.to_csv(os.path.join(parent_dir, "obspy_result.csv"), index=False)
    
    if doAlert:
        intensityListString = [i[4] for i in rows]
        intensityList = [x for x in [intensityConvert(i[4]) for i in rows] if x]
        displacementList = map(float,[x for x in [i[3] for i in rows] if x])
        if any(i>=intThresh for i in intensityList):
            intensity = intensityListString[np.argmax(intensityList)]
            displacement = max(displacementList)*100 #convert to centimeters
            # eqAlarm.alarm(intensity,displacement)
            
            thread = threading.Thread(target=threadMain(intensity,displacement))
            thread.start()
