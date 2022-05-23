from obspy import read_inventory
from obspy.core import read
from obspy.clients.fdsn import Client as RS_Client
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

INV_DIR = "inventories"
sta_window = 2
lta_window = 10
thresh_on = 1.5
thresh_off = 1

def get_PRISM_data(filename):
    with open(filename) as f:
        record = False
        start_count = False
        cnt = 0
        num_to_skip_at_start = 0 # for v1 files
        if ".V2c" in filename:
            num_to_skip_at_start = 5
            num_to_skip_at_start = 0
            if "BVDA" in filename:
                num_to_skip_at_start = 8
            num_to_skip_at_start = 0
        data = []
        for line in f:
            if start_count:
                cnt += 1
            if "End-of-data" in line:
                record = False
            if record and (cnt > num_to_skip_at_start):
                val = float(line.strip())
                data.append(val)
                start_count = False # no need to continuously increment
            #if "<SCNL>" in line:
            if "pts" in line:
                record = True
                start_count = True

        if len(data) == 0:
            warnings.warn("No data parsed")

        return np.array(data)

def get_event_info(prism_filename):
    with open(prism_filename) as f:
        target_line_number = 3
        line_number = 0
        target_line = ""

        # get line containing desired info
        for line in f:
            line_number += 1
            if line_number == target_line_number:
                target_line = line
                break

        # extract floats from line
        floats = (float(s) for s in re.findall(r"[-+]?(?:\d*\.\d+|\d+)",line))
        lat, lon, depth, mag = floats

    return lat, lon, depth, mag

def convert_acc_to_vel_trace(tr):
    if tr.stats.units == "ACC":
        tr = improved_integration(tr)
        tr.stats.units = "VEL"
    elif tr.stats.units == "VEL":
        tr.stats.units = "VEL"
    else:
        print("Can't convert", tr.stats.units, "to VEL.")

    return tr

def convert_vel_to_dis_trace(tr):
    if tr.stats.units == "VEL":
        tr = improved_integration(tr)
        tr.stats.units = "DIS"
    elif tr.stats.units == "DIS":
        tr.stats.units = "DIS"
    else:
        print("Can't convert", tr.stats.units, "to DIS.")

    return tr

def get_weighted_mean(x, weights):
    sum_weight = integrate.trapz(weights)
    weighted_x = x*weights

    return sum(weighted_x)/sum_weight

def get_weighted_variance(x, weights):
    sum_weight = integrate.trapz(weights)
    weighted_mean = sum(x*weights)/sum_weight
    var = sum(weights*(x - weighted_mean)*(x - weighted_mean))/sum_weight

    return var

def improved_integration(tr):
    tr = tr.copy()
    tr.detrend("demean") # make sure no constant that will become linear function
    tr.integrate(method="cumtrapz")
    tr.detrend("linear") # (mx+b, ie the leakage due to cumtrapz)

    return tr

def get_low_corner_freq(tr, event_onset = None, low_power_thresh = 0.015, noise_type="use_end", plot=False, save_plot=False, plot_info=None, verbose=False):
    global event_type
    tr = tr.copy()
    if event_onset == None:
        # get event onset
        freq = tr.stats.sampling_rate
        sta_lta = recursive_sta_lta(tr, nsta=int(sta_window*freq), nlta=int(lta_window*freq))
        pick_pairs_ind = trigger_onset(sta_lta, thresh_on, thresh_off)
        event_onset = pick_pairs_ind[0][0]
        #print("event onset", str(event_onset))

    if noise_type == "use_end":
        noise = tr.data[-1*int(event_onset-1*tr.stats.sampling_rate):] # remove one second worth of samples
        #noise = tr.data[event_onset:event_onset+10000] # use signal
        #noise = tr.data[-1*int(event_onset):int(-2*tr.stats.sampling_rate)] # remove two second worth of samples, at the right
    elif noise_type == "lowest_ave":
        tr = tr.copy()
        tr.detrend("linear") # remove dc
        window_size = event_onset # note that event onset is an index
        windowed_ave_arr = np.convolve(np.absolute(tr.data[window_size:]), np.ones(window_size), 'valid')/window_size
                            # convolve takes inner product of the two arrays and sum them,
                            # take abs to avoid cancelling of positive and negative amps
                            # skip first window size to avoid it being the minimum
                            # valid is to only perform when overlap is complete
        last_occur_min_index = np.argmin(windowed_ave_arr[::-1])
        last_occur_min_index = len(windowed_ave_arr) - last_occur_min_index - 1
                            # get first occur of min in reversed arr to get last occur
                            # then flip to normal indices
        #min_window_end_index = 2*window_size + np.argmin(windowed_ave_arr)
        noise_end_ind = 2*window_size + last_occur_min_index
        noise_start_ind = noise_end_ind - window_size
        noise = tr.data[noise_start_ind : noise_end_ind + 1]

        #print("---------")
        #print("w =", window_size)
        #print("len tr.data =", len(tr.data))
        #print("len window_ave_arr =", len(windowed_ave_arr))
        #print("last_occur_min_ind =", last_occur_min_index)
        #print("window_ave arr =", windowed_ave_arr)
        #print("argmin", np.argmin(windowed_ave_arr))
        #print("windowed_ave_arr[argmin+1]", windowed_ave_arr[np.argmin(windowed_ave_arr)+1])
        #print("windowed_ave_arr[argmin]", windowed_ave_arr[np.argmin(windowed_ave_arr)])
        #print("windowed_ave_arr[argmin-1]", windowed_ave_arr[np.argmin(windowed_ave_arr)-1])



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

        fig.suptitle("Noise Removal")

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
        if noise_type == "lowest_ave":
            axs[0][1].plot(np.array(range(2*window_size,2*window_size+len(windowed_ave_arr))), windowed_ave_arr) # also plot windowed ave
            #axs[0][1].plot(windowed_ave_arr) # also plot windowed ave
            y_min, y_max = axs[0][1].get_ylim()
            axs[0][1].vlines(noise_start_ind, y_min, y_max, color="r", lw=1)
            axs[0][1].vlines(noise_end_ind, y_min, y_max, color="b", lw=1)

        axs[1][1].set_ylabel("Signal Response\nCounts")
        axs[1][1].plot(signal_resp_x,signal_resp)
        axs[2][1].set_ylabel("Final Signal \nCounts")
        axs[2][1].plot(final_signal)

        for row in axs:
            for ax in row:
                ax.grid()
        #plt.show()
        station, mag =  plot_info
        fig_path = os.path.join(os.getcwd(), "results", station+"_results", "M"+str(mag)+"_noise.png")
        Path(os.path.dirname(fig_path)).mkdir(parents=True, exist_ok=True)
        plt.savefig(fig_path)


    return lowcut

def convert_counts_to_metric_trace(tr, metric_units, limit=True, plot=False, plot_info=None):
    if tr.stats.units == "COUNTS":
        tr = tr.copy()
        freq = tr.stats.sampling_rate
        #print(estimate_magnitude(inv[0][0][43].response ,max(abs(tr.data)), 0.1, 20))
        # Note: The following are applied before deconvolution
        # a bandpass filter
        # a quarter cosine taper, tapering 0.05 length of data
        lowcut = get_low_corner_freq(tr, noise_type="lowest_ave", plot=plot, plot_info=plot_info)
        if limit:
            lowcut = lowcut if lowcut <= 0.5 else 0.5
        #print("used lowcut Hz:", lowcut)
        tr.remove_response(pre_filt=[lowcut, lowcut*1.1, 0.5*freq*0.9, 0.5*freq],
        #tr.remove_response(pre_filt=[0.1, 0.5, 49, 50],
                           output=metric_units, water_level=4.5, taper=True, taper_fraction=0.1)
        tr.stats.units = metric_units
    return tr, lowcut

def get_inventory(inv_path, network, station, client_name="IRIS"):
    # create copy of latest inv, remove date so it can be used in remove_response
    if os.path.exists(inv_path):
        #print("reading inv")
        inv = read_inventory(inv_path)
    else:
        #print("downloading inv")
        rs_client = RS_Client(client_name)
        inv = rs_client.get_stations(network=network, station=station, level="RESP")
        latest_station_response = (inv[-1][-1]).copy()
        latest_station_response.start_date = None
        latest_station_response.end_date = None
        for cha in latest_station_response:
            cha.start_date=None
            cha.end_date=None
        inv = Inventory(networks=[ \
                Network(code=network, stations=[latest_station_response])])
        Path(os.path.dirname(inv_path)).mkdir(parents=True, exist_ok=True)
        inv.write(inv_path, format="STATIONXML")
    return inv

def compare(prism_files, mseed_file, network, station, channel, plot= True):
    # PRISM-processed-data import
    prism_v1_file, prism_acc_file, prism_vel_file, prism_dis_file = prism_files

    _, _, depth, mag = get_event_info(prism_v1_file)
    prism_v1 = get_PRISM_data(prism_v1_file)/100
    prism_v1_mean = prism_v1.mean()
    prism_v1_peak = max(abs(prism_v1))

    prism_acc = get_PRISM_data(prism_acc_file)/100 # original data is in cm
    prism_acc_mean = prism_acc.mean()
    prism_acc_peak = max(abs(prism_acc))

    prism_vel = get_PRISM_data(prism_vel_file)/100
    prism_vel_mean = prism_vel.mean()
    prism_vel_peak = max(abs(prism_vel))

    prism_dis = get_PRISM_data(prism_dis_file)/100
    prism_dis_mean = prism_dis.mean()
    prism_dis_peak = max(abs(prism_dis))


    # Obspy processing
    st = read(mseed_file, format="MSEED")
    metric_units = "ACC"

    inv_path = os.path.join(os.getcwd(), INV_DIR, network+"_"+station+".xml")
    inv = get_inventory(inv_path, network, station)

    counts = st.select(channel=channel)[0]
    #counts.data = counts.data - counts.data.mean()
    #freq = counts.stats.sampling_rate
    #plot_trigger(counts, recursive_sta_lta(counts, nsta=int(sta_window*freq), nlta=int(lta_window*freq)), 1.5, 1)
    counts.attach_response(inv)
    counts.stats.units = "COUNTS"
    counts_mean = counts.data.mean()
    counts_peak = max(abs(counts.data))
    acc, lowcut = convert_counts_to_metric_trace(counts, metric_units, limit=True, \
                                                plot=True, plot_info=(station, mag))
    acc_mean = acc.data.mean()
    acc_peak = max(abs(acc.data))
    vel = convert_acc_to_vel_trace(acc)
    vel_mean = vel.data.mean()
    vel_peak = max(abs(vel.data))
    dis = convert_vel_to_dis_trace(vel)
    dis_mean = dis.data.mean()
    dis_peak = max(abs(dis.data))


    # Absolute errors between prism and obspy
    accv1_ae = abs(acc.data - prism_v1)
    accv1_ae_mean = accv1_ae.mean()
    accv1_ae_peak = max(abs(accv1_ae))

    acc_ae = abs(acc.data - prism_acc)
    acc_ae_mean = acc_ae.mean()
    acc_ae_peak = max(abs(acc_ae))

    vel_ae = abs(vel.data - prism_vel)
    vel_ae_mean = vel_ae.mean()
    vel_ae_peak = max(abs(vel_ae))

    dis_ae = abs(dis.data - prism_dis)
    dis_ae_mean = dis_ae.mean()
    dis_ae_peak = max(abs(dis_ae))


    # Plot results
    if plot:
        fig, axs = plt.subplots(nrows=4, ncols=3, constrained_layout=True)
        fig.set_size_inches(w=17, h=9)

        fig.suptitle(station+" in Metric Units: M"+str(mag))

        # OBSPY PLOT COLUMN
        axs[0][0].set_ylabel("OBSPY\nCounts")
        axs[0][0].plot(counts.data)
        x_text0 = len(counts.data)*0.65
        y_text0 = max(counts.data)*0.7
        axs[0][0].text(x_text0,y_text0,
                    "Mean = {:+.6f}counts\nPeak = {:.6f}counts".format(counts_mean,counts_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[1][0].set_ylabel("OBSPY\nAcceleration")
        axs[1][0].plot(acc.data)
        x_text0 = len(acc.data)*0.65
        y_text0 = max(acc.data)*0.7
        axs[1][0].text(x_text0,y_text0,
                    "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(acc_mean,acc_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[2][0].set_ylabel("OBSPY\nVelocity")
        axs[2][0].plot(vel.data)
        x_text1 = len(vel.data)*0.65
        y_text1 = max(vel.data)*0.7
        axs[2][0].text(x_text1,y_text1,
                    "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(vel_mean,vel_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[3][0].set_ylabel("OBSPY\nDisplacement")
        axs[3][0].plot(dis.data)
        x_text2 = len(dis.data)*0.65
        y_text2 = max(dis.data)*0.7
        axs[3][0].text(x_text2,y_text2,
                    "Mean = {:+.6f}m\nPeak = {:.6f}m".format(dis_mean,dis_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        # PRISM PLOT COLUMN
        axs[0][1].set_ylabel("PRISM V1\n(scaled counts)")
        axs[0][1].plot(prism_v1)
        x_text0 = len(prism_v1)*0.65
        y_text0 = max(prism_v1)*0.7
        axs[0][1].text(x_text0,y_text0,
                    "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_v1_mean,prism_v1_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[1][1].set_ylabel("PRISM V2\nAcceleration")
        axs[1][1].plot(prism_acc)
        x_text0 = len(prism_acc)*0.65
        y_text0 = max(prism_acc)*0.7
        axs[1][1].text(x_text0,y_text0,
                    "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_acc_mean,prism_acc_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[2][1].set_ylabel("PRISM V2\nVelocity")
        axs[2][1].plot(prism_vel)
        x_text0 = len(prism_vel)*0.65
        y_text0 = max(prism_vel)*0.7
        axs[2][1].text(x_text0,y_text0,
                    "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(prism_vel_mean,prism_vel_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[3][1].set_ylabel("PRISM V2\nDisplacement")
        axs[3][1].plot(prism_dis)
        x_text0 = len(prism_dis)*0.65
        y_text0 = max(prism_dis)*0.7
        axs[3][1].text(x_text0,y_text0,
                    "Mean = {:+.6f}m\nPeak = {:.6f}m".format(prism_dis_mean,prism_dis_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        # ABSOLUTE ERRORS COLUMN
        axs[0][2].set_ylabel("OBSPY_MSDvsPRISM_V1\nAccel AbsErr")
        axs[0][2].plot(abs(accv1_ae))
        x_text2 = len(accv1_ae)*0.65
        y_text2 = max(accv1_ae)*0.7
        axs[0][2].text(x_text2,y_text2,
                    "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(accv1_ae_mean,accv1_ae_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[1][2].set_ylabel("OBSPYvsPRISM_V2\nAccel AbsErr")
        axs[1][2].plot(abs(acc_ae))
        x_text2 = len(acc_ae)*0.65
        y_text2 = max(acc_ae)*0.7
        axs[1][2].text(x_text2,y_text2,
                    "Mean = {:+.6f}m/2\nPeak = {:.6f}m/s2".format(acc_ae_mean,acc_ae_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[2][2].set_ylabel("OBSPYvsPRISM_V2\nVelocity AbsErr")
        axs[2][2].plot(abs(vel_ae))
        x_text2 = len(vel_ae)*0.65
        y_text2 = max(vel_ae)*0.7
        axs[2][2].text(x_text2,y_text2,
                    "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(vel_ae_mean,vel_ae_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        axs[3][2].set_ylabel("OBSPYvsPRISM_V2\nDisp AbsErr")
        axs[3][2].plot(abs(dis_ae))
        x_text2 = len(dis_ae)*0.65
        y_text2 = max(dis_ae)*0.7
        axs[3][2].text(x_text2,y_text2,
                    "Mean = {:+.6f}m\nPeak = {:.6f}m".format(dis_ae_mean,dis_ae_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))

        # get max vertical axis size
        cts_ymin, cts_ymax = 0, 0
        acc_ymin, acc_ymax = 0, 0
        vel_ymin, vel_ymax = 0, 0
        dis_ymin, dis_ymax = 0, 0
        for i in range(3):
            #_, _, cts_ymin_tmp, cts_ymax_tmp = axs[0][i].axis()
            #if cts_ymin_tmp < cts_ymin:
            #    cts_ymin = cts_ymin_tmp
            #if cts_ymax_tmp > cts_ymax:
            #    cts_ymax= cts_ymax_tmp

            _, _, acc_ymin_tmp, acc_ymax_tmp = axs[1][i].axis()
            if acc_ymin_tmp < acc_ymin:
                acc_ymin = acc_ymin_tmp
            if acc_ymax_tmp > acc_ymax:
                acc_ymax= acc_ymax_tmp

            _, _, vel_ymin_tmp, vel_ymax_tmp = axs[2][i].axis()
            if vel_ymin_tmp < vel_ymin:
                vel_ymin = vel_ymin_tmp
            if vel_ymax_tmp > vel_ymax:
                vel_ymax= vel_ymax_tmp

            _, _, dis_ymin_tmp, dis_ymax_tmp = axs[3][i].axis()
            if dis_ymin_tmp < dis_ymin:
                dis_ymin = dis_ymin_tmp
            if dis_ymax_tmp > dis_ymax:
                dis_ymax= dis_ymax_tmp


        # set vertical axis sizes
        for i in range(3):
            _, _, cts_ymin, cts_ymax= axs[0][i].axis() # don't set for counts
            axs[0][i].set_ylim([cts_ymin, cts_ymax])
            axs[1][i].set_ylim([acc_ymin, acc_ymax])
            axs[2][i].set_ylim([vel_ymin, vel_ymax])
            axs[3][i].set_ylim([dis_ymin, dis_ymax])

        # set row grids
        for row in axs:
            for ax in row:
                ax.grid()

        #plt.show()
        fig_path = os.path.join(os.getcwd(), "results", station+"_results", "M"+str(mag)+"_conv.png")
        Path(os.path.dirname(fig_path)).mkdir(parents=True, exist_ok=True)
        plt.savefig(fig_path)


    event_info = [depth, mag]
    processing_info = [lowcut]
    instrument_info = [network, station, channel]
    prism_peaks = [prism_acc_peak, prism_vel_peak, prism_dis_peak]
    obspy_peaks = [acc_peak, vel_peak, dis_peak]
    mean_abs_errors = [acc_ae_mean, vel_ae_mean, dis_ae_mean]

    return (event_info, processing_info, instrument_info, prism_peaks, obspy_peaks, mean_abs_errors)

if __name__ == "__main__":
    #print(compare(prism_files, mseed_file, network, station, channel, True ))
    # load prism and mseed
    network = "SB"
    channel = "HNZ"
    stations = ["WLA", "GVDA", "BVDA"]
    events = ["Event_1","Event_2","Event_3"]
    parent_dir = "/home/jeffsanchez/Downloads/cosmos_downloads_v2"

    columns=["Station", "Event #", "Mag", "--", \
             "OBS-PeakAcc", "PSM-PeakAcc", "PeakAccDiff", "%AccDiff","--", \
             "OBS-PeakVel", "PSM-PeakVel", "PeakVelDiff", "%VelDiff","--", \
             "OBS-PeakDis", "PSM-PeakDis", "PeakDisDiff", "%DisDiff", "--", \
             "lowcut"]
    rows = []

    for station in stations:
        station_dir = parent_dir + "/" + station
        for event in  events:
            station_event_dir = station_dir + "/" + event

            # get event tag and evid directory
            for dir_name in os.listdir(station_event_dir):
                if "evid" in dir_name:
                    event_tag = dir_name.split("_")[1]
                    evid_dir = dir_name
                    break

            mseed_file     = os.path.join(station_event_dir, evid_dir, "_".join((network, station, channel, "00", event_tag)) + ".msd")
            prism_v1_file = os.path.join(station_event_dir, "prismPROCESSED", "CI"+"."+event_tag, network+"."+station,"V1","_".join((network, station, channel, "00", event_tag)) + ".V1c")
            prism_acc_file = os.path.join(station_event_dir, "prismPROCESSED", "CI"+"."+event_tag, network+"."+station,"V2","_".join((network, station, channel, "00", event_tag)) + ".acc.V2c")
            prism_vel_file = os.path.join(station_event_dir, "prismPROCESSED", "CI"+"."+event_tag, network+"."+station,"V2","_".join((network, station, channel, "00", event_tag)) + ".vel.V2c")
            prism_dis_file = os.path.join(station_event_dir, "prismPROCESSED", "CI"+"."+event_tag, network+"."+station,"V2","_".join((network, station, channel, "00", event_tag)) + ".dis.V2c")
            prism_files = (prism_v1_file, prism_acc_file, prism_vel_file, prism_dis_file)

            #print(station)
            #print(event)
            (event_info, processing_info, instrument_info, prism_peaks, obspy_peaks, mean_abs_errors) = \
                compare(prism_files, mseed_file, network, station, channel, False )

            # collect output
            lowcut = processing_info[0]
            (acc, vel, dis) = obspy_peaks
            (p_acc, p_vel, p_dis) = prism_peaks
            (d_acc, d_vel, d_dis) = (abs(a-b) for a,b in zip(obspy_peaks, prism_peaks))
            (prcnt_acc, prcnt_vel, prcnt_dis) = (abs(a-b)/b for a,b in zip(obspy_peaks, prism_peaks))
            rows.append([station, event.split("_")[1], event_info[1], "", \
                         acc, p_acc, d_acc, prcnt_acc, "",
                         vel, p_vel, d_vel, prcnt_vel, "",
                         dis, p_dis, d_dis, prcnt_dis, "",
                         lowcut
            ])
        rows.append([""]*20)

    df = pd.DataFrame(columns=columns,data=rows)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(df)
    df.to_csv(parent_dir + "/comparison_result.csv", index=False)




