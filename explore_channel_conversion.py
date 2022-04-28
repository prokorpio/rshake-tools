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
import pandas as pd


INV_DIR = "inventories"
event_type = 1
if event_type == 1:
    mseed_file     = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_1/evid_15200401_13519/SB_WLA_HNZ_00_15200401.msd"
    prism_v1_file  = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_1/prismPROCESSED/CI.15200401/SB.WLA/V1/SB_WLA_HNZ_00_15200401.V1c"
    prism_acc_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_1/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.acc.V2c"
    prism_vel_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_1/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.vel.V2c"
    prism_dis_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_1/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.dis.V2c"
elif event_type == 2:
    mseed_file     = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_2/evid_14733516_4612/SB_WLA_HNZ_00_14733516.msd"
    prism_v1_file  = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_2/prismPROCESSED/CI.14733516/SB.WLA/V1/SB_WLA_HNZ_00_14733516.V1c"
    prism_acc_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_2/prismPROCESSED/CI.14733516/SB.WLA/V2/SB_WLA_HNZ_00_14733516.acc.V2c"
    prism_vel_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_2/prismPROCESSED/CI.14733516/SB.WLA/V2/SB_WLA_HNZ_00_14733516.vel.V2c"
    prism_dis_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_2/prismPROCESSED/CI.14733516/SB.WLA/V2/SB_WLA_HNZ_00_14733516.dis.V2c"
elif event_type == 3:
    mseed_file     = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/evid_10665149_4386/SB_WLA_HNZ_00_10665149.msd"
    prism_v1_file  = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V1/SB_WLA_HNZ_00_10665149.V1c"
    prism_acc_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.acc.V2c"
    prism_vel_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.vel.V2c"
    prism_dis_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.dis.V2c"
prism_files = (prism_v1_file, prism_acc_file, prism_vel_file, prism_dis_file)
network = "SB"
station = "WLA"
channel = "HNZ"

sta_window = 2
lta_window = 10
thresh_on = 1.5
thresh_off = 1

def select_event(event_type):
    if event_type == 1:
        mseed_file     = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_1/evid_15200401_13519/SB_WLA_HNZ_00_15200401.msd"
        prism_v1_file  = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_1/prismPROCESSED/CI.15200401/SB.WLA/V1/SB_WLA_HNZ_00_15200401.V1c"
        prism_acc_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_1/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.acc.V2c"
        prism_vel_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_1/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.vel.V2c"
        prism_dis_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_1/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.dis.V2c"
    elif event_type == 2:
        mseed_file     = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_2/evid_14733516_4612/SB_WLA_HNZ_00_14733516.msd"
        prism_v1_file  = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_2/prismPROCESSED/CI.14733516/SB.WLA/V1/SB_WLA_HNZ_00_14733516.V1c"
        prism_acc_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_2/prismPROCESSED/CI.14733516/SB.WLA/V2/SB_WLA_HNZ_00_14733516.acc.V2c"
        prism_vel_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_2/prismPROCESSED/CI.14733516/SB.WLA/V2/SB_WLA_HNZ_00_14733516.vel.V2c"
        prism_dis_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_2/prismPROCESSED/CI.14733516/SB.WLA/V2/SB_WLA_HNZ_00_14733516.dis.V2c"
    elif event_type == 3:
        mseed_file     = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_3/evid_10665149_4386/SB_WLA_HNZ_00_10665149.msd"
        prism_v1_file  = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V1/SB_WLA_HNZ_00_10665149.V1c"
        prism_acc_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.acc.V2c"
        prism_vel_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.vel.V2c"
        prism_dis_file = "C:/Users/benit/Documents/Python Scripts/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.dis.V2c"
        
        
    return (prism_v1_file, prism_acc_file, prism_vel_file, prism_dis_file), mseed_file

def get_PRISM_data(filename):
    with open(filename) as f:
        record = False
        start_count = False
        cnt = 0
        num_to_skip_at_start = 1 # for v1 files
        if ".V2c" in filename:
            num_to_skip_at_start = 5
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
            if "<SCNL>" in line:
                record = True
                start_count = True

        if len(data) == 0:
            warnings.warn("No data parsed")

        return np.array(data)

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

def convert_counts_to_metric_trace(tr, metric_units):
    if tr.stats.units == "COUNTS":
        tr = tr.copy()
        freq = tr.stats.sampling_rate
        #print(estimate_magnitude(inv[0][0][43].response ,max(abs(tr.data)), 0.1, 20))
        # Note: The following are applied before deconvolution
        # a bandpass filter
        # a quarter cosine taper, tapering 0.05 length of data
        lowcut = get_low_corner_freq(tr, plot=True)
        lowcut = lowcut if lowcut <= 0.5 else 0.5
        print("used lowcut Hz:", lowcut)
        tr.remove_response(pre_filt=[lowcut, lowcut*1.1, 0.5*freq*0.9, 0.5*freq],
        #tr.remove_response(pre_filt=[0.1, 0.5, 49, 50],
                           output=metric_units, water_level=4.5, taper=True, taper_fraction=0.1)
        tr.stats.units = metric_units
    return tr

def get_low_corner_freq(tr, event_onset = None, low_power_thresh = 0.015, use_end=False, plot=False):
    global event_type
    tr = tr.copy()
    if event_onset == None:
        # get event onset
        freq = tr.stats.sampling_rate
        sta_lta = recursive_sta_lta(tr, nsta=int(sta_window*freq), nlta=int(lta_window*freq))
        pick_pairs_ind = trigger_onset(sta_lta, thresh_on, thresh_off)
        event_onset = pick_pairs_ind[0][0]

    if use_end:
        noise = tr.data[-1*int(event_onset-1*tr.stats.sampling_rate):] # remove one second worth of samples
    else: # get noise from beginning
        noise = tr.data[:int(event_onset-1*tr.stats.sampling_rate)] # remove one second worth of samples
    noise = noise - np.mean(noise) #remove DC
    noise_resp = abs(np.fft.rfft(noise, n=_npts2nfft(len(noise)))) # get distance as magnitude
    noise_resp_x = np.fft.rfftfreq(_npts2nfft(len(noise)), 1/tr.stats.sampling_rate)
    nyquist_ind = len(noise_resp)//2
    noise_integral = integrate.cumtrapz(noise_resp[:nyquist_ind], noise_resp_x[:nyquist_ind], initial=0) # only up to nyquist freq
    noise_integral = noise_integral/noise_integral[-1] # normalize w max
    lowcut_ind = int((noise_integral > low_power_thresh).nonzero()[0][0])
    lowcut = noise_resp_x[lowcut_ind]
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

def get_inventory(inv_path, network, station, client_name="IRIS"):
    # create copy of latest inv, remove date so it can be used in remove_response
    if os.path.exists(inv_path):
        print("reading inv")
        inv = read_inventory(inv_path)
    else:
        print("downloading inv")
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

def compare(prism_files, mseed_file, network, station, channel, plot_results = True):
    # PRISM-processed-data import
    prism_v1_file, prism_acc_file, prism_vel_file, prism_dis_file = prism_files

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
    #freq = counts.stats.sampling_rate
    #plot_trigger(counts, recursive_sta_lta(counts, nsta=int(sta_window*freq), nlta=int(lta_window*freq)), 1.5, 1)
    counts.attach_response(inv)
    counts.stats.units = "COUNTS"
    counts_mean = counts.data.mean()
    counts_peak = max(abs(counts.data))
    acc = convert_counts_to_metric_trace(counts, metric_units)
#    return
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
#    print(dis_ae_mean)
    dis_ae_peak = max(abs(dis_ae))


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

        # PRISM PLOT COLUMN
        axs[0][1].set_ylabel("PRISM V1\n(scaled counts)")
        axs[0][1].plot(prism_v1)
        x_text0 = len(prism_v1)*0.65
        y_text0 = max(prism_v1)*0.7
        axs[0][1].text(x_text0,y_text0,
                    "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_v1_mean,prism_v1_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        axs[0][1].set_ylim([acc_ymin, acc_ymax])

        axs[1][1].set_ylabel("PRISM V2\nAcceleration")
        axs[1][1].plot(prism_acc)
        x_text0 = len(prism_acc)*0.65
        y_text0 = max(prism_acc)*0.7
        axs[1][1].text(x_text0,y_text0,
                    "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_acc_mean,prism_acc_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        axs[1][1].set_ylim([acc_ymin, acc_ymax])

        axs[2][1].set_ylabel("PRISM V2\nVelocity")
        axs[2][1].plot(prism_vel)
        x_text0 = len(prism_vel)*0.65
        y_text0 = max(prism_vel)*0.7
        axs[2][1].text(x_text0,y_text0,
                    "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(prism_vel_mean,prism_vel_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        axs[2][1].set_ylim([vel_ymin, vel_ymax])

        axs[3][1].set_ylabel("PRISM V2\nDisplacement")
        axs[3][1].plot(prism_dis)
        x_text0 = len(prism_dis)*0.65
        y_text0 = max(prism_dis)*0.7
        axs[3][1].text(x_text0,y_text0,
                    "Mean = {:+.6f}m\nPeak = {:.6f}m".format(prism_dis_mean,prism_dis_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        axs[3][1].set_ylim([dis_ymin, dis_ymax])

        # ABSOLUTE ERRORS COLUMN
        axs[0][2].set_ylabel("OBSPY_MSDvsPRISM_V1\nAccel AbsErr")
        axs[0][2].plot(abs(accv1_ae))
        x_text2 = len(accv1_ae)*0.65
        y_text2 = acc_ymax*0.7
        axs[0][2].text(x_text2,y_text2,
                    "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(accv1_ae_mean,accv1_ae_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        axs[0][2].set_ylim([acc_ymin, acc_ymax])

        axs[1][2].set_ylabel("OBSPYvsPRISM_V2\nAccel AbsErr")
        axs[1][2].plot(abs(acc_ae))
        x_text2 = len(acc_ae)*0.65
        y_text2 = acc_ymax*0.7
        axs[1][2].text(x_text2,y_text2,
                    "Mean = {:+.6f}m/2\nPeak = {:.6f}m/s2".format(acc_ae_mean,acc_ae_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        axs[1][2].set_ylim([acc_ymin, acc_ymax])

        axs[2][2].set_ylabel("OBSPYvsPRISM_V2\nVelocity AbsErr")
        axs[2][2].plot(abs(vel_ae))
        x_text2 = len(vel_ae)*0.65
        y_text2 = vel_ymax*0.7
        axs[2][2].text(x_text2,y_text2,
                    "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(vel_ae_mean,vel_ae_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        axs[2][2].set_ylim([vel_ymin, vel_ymax])

        axs[3][2].set_ylabel("OBSPYvsPRISM_V2\nDisp AbsErr")
        axs[3][2].plot(abs(dis_ae))
        x_text2 = len(dis_ae)*0.65
        y_text2 = dis_ymax*0.7
        axs[3][2].text(x_text2,y_text2,
                    "Mean = {:+.6f}m\nPeak = {:.6f}m".format(dis_ae_mean,dis_ae_peak),
                    ha="left", va="center",
                    bbox=dict(boxstyle="round", fc="white"))
        axs[3][2].set_ylim([dis_ymin, dis_ymax])

        for row in axs:
            for ax in row:
                ax.grid()
        plt.show()

    event_info = () # not yet implemented, parse from cosmos files
    instrument_info = [network, station, channel]
    prism_peaks = [prism_acc_peak, prism_vel_peak, prism_dis_peak]
    obspy_peaks = [acc_peak, vel_peak, dis_peak]
    mean_abs_errors = [acc_ae_mean, vel_ae_mean, dis_ae_mean]

    return (event_info, instrument_info, prism_peaks, obspy_peaks, mean_abs_errors)

if __name__ == "__main__":
#    print(compare(prism_files, mseed_file, network, station, channel, True ))

    numEvents = 3
    sampleDict = {}
    for event in range(1,numEvents+1):
        prism_files, mseed_file = select_event(event)
        allData = compare(prism_files, mseed_file, network, station, channel, False )
        #event_info, instrument_info, prism_peaks, obspy_peaks, mean_abs_errors = allData
        sampleDict[event] = allData
    
    count = 0
    prism_peaksAll = []
    obspy_peaksAll = []
    mean_abs_errorsAll = []
    peak_diffAll = []
    peak_diff_gAll = []
    for key in sampleDict:
        prism_peaks = np.array(sampleDict[key][2])
        obspy_peaks = np.array(sampleDict[key][3])
        mean_abs_errors = np.array( sampleDict[key][4])
        peak_diff = np.abs(np.subtract(prism_peaks,obspy_peaks))
        peak_diff_g = np.abs(np.subtract(prism_peaks,obspy_peaks)/9.81)

        prism_peaksAll.append(prism_peaks)
        obspy_peaksAll.append(obspy_peaks)
        mean_abs_errorsAll.append(mean_abs_errors)
        peak_diffAll.append(peak_diff)
        peak_diff_gAll.append(peak_diff_g)

    prism_peaksAll = np.array(prism_peaksAll)
    obspy_peaksAll = np.array(obspy_peaksAll)
    mean_abs_errorsAll = np.array(mean_abs_errorsAll)
    peak_diffAll = np.array(peak_diffAll)
    peak_diff_gAll = np.array(peak_diff_gAll)


    for i in range(0,numEvents):
        d = {
        'Event' : [1,2,3],
        'Prism PGA (m/s2)' : prism_peaksAll[:,i],
        'Obspy PGA (m/s2)' : obspy_peaksAll[:,i],
        'Peak Diff (m/s2)' : peak_diffAll[:,i],
        'Peak Difference (g)' : peak_diff_gAll[:,i],
        'Mean Abs Errors' : mean_abs_errorsAll[:,i]
        }

        df = pd.DataFrame(data=d)
        print(df)
        df.to_csv(r'C:\Users\benit\Documents\Python Scripts\rshake-tools-dev\Prism_Obspy_Table_' + str(i+1) + '.csv',index=False)
