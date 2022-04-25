from obspy import read_inventory
from obspy.core import read
from obspy.clients.fdsn import Client as RS_Client
from obspy.core.inventory import Inventory, Network
from pathlib import Path
import matplotlib.pyplot as plt
import os
import numpy as np
import warnings

INV_DIR = "inventories"
event_type = 3
prism_v1_file  = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V1/SB_WLA_HNZ_00_10665149.V1c"
prism_acc_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.acc.V2c"
prism_vel_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.vel.V2c"
prism_dis_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/prismPROCESSED/CI.10665149/SB.WLA/V2/SB_WLA_HNZ_00_10665149.dis.V2c"
prism_files = (prism_v1_file, prism_acc_file, prism_vel_file, prism_dis_file)
mseed_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Event_3/evid_10665149_4386/SB_WLA_HNZ_00_10665149.msd"
network = "SB"
station = "WLA"
channel = "HNZ"

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
        tr.remove_response(pre_filt=[0.1, 0.5, 0.95*freq, freq],
                           output=metric_units, water_level=4.5, taper=False)
        tr = correct_acceleration(tr)
        tr.stats.units = metric_units
    return tr

def correct_acceleration(tr):
    tr = tr.copy()

    # remove slew, assumption is stationary
    tr.detrend("linear")

    # integrate to velocity and detrend acc via velocity's diff'd trend
    # might not be necessary since linear detrend seems to be enough
    #vel = improved_integration(tr)
    #one_fit_coeffs, one_score, _,_,_ = np.polyfit(np.arange(tr.stats.npts),vel.data,1,full=True)
    #two_fit_coeffs, two_score, _,_,_ = np.polyfit(np.arange(tr.stats.npts),vel.data,2,full=True)
    #best_fit = np.poly1d(one_fit_coeffs if one_score > two_score else two_fit_coeffs) # make func
    #diffed_best_fit = np.gradient(best_fit(np.arange(tr.stats.npts)))
    #tr.data = tr.data - diffed_best_fit

    # Length is based on the the expectation the the data has 15sec pre-event buffer
    # p (cosine percentage), is set to be similar to prism implementation (they have no p)
    tr.taper(max_percentage=0.25, max_length=5, side="both", type="cosine", p=1)
    # no zero padding bc it's only necessary to accomodate cyclic conv property of
    # inverse FT, when integrating in f-domain

    if event_type == 3:
        freqmin= 0.5
        freqmax= 25
    elif event_type == 2:
        freqmin= 0.3
        freqmax= 35
    else:
        freqmin= 0.1
        freqmax= 40
    tr.filter("bandpass",freqmin=freqmin,freqmax=freqmax) # based on magnitude
    return tr

def improved_integration(tr):
    tr = tr.copy()
    tr.detrend("demean")
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
    counts.attach_response(inv)
    counts.stats.units = "COUNTS"
    counts_mean = counts.data.mean()
    counts_peak = max(abs(counts.data))
    acc = convert_counts_to_metric_trace(counts, metric_units)
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
    instrument_info = (network, station, channel)
    prism_peaks = (prism_acc_peak, prism_vel_peak, prism_dis_peak)
    obspy_peaks = (acc_peak, vel_peak, dis_peak)
    mean_abs_errors = (acc_ae_mean, vel_ae_mean, dis_ae_mean)

    return (event_info, instrument_info, prism_peaks, obspy_peaks, mean_abs_errors)

if __name__ == "__main__":
    print(compare(prism_files, mseed_file, network, station, channel, True ))



