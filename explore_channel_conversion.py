from obspy import read_inventory
from obspy.core import read
from obspy.clients.fdsn import Client as RS_Client
from obspy.core.inventory import Inventory, Network
from pathlib import Path
import matplotlib.pyplot as plt
import os
import numpy as np

#
def get_PRISM_data(filename):
    with open(filename) as f:
        record = False
        data = []
        for line in f:
            if "End-of-data" in line:
                record = False
            if record:
                val = float(line.strip())
                data.append(val)
            if "   50400 " in line:
                record = True

        return np.array(data)

# PRISM-processed-data import
prism_v1_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/prismPROCESSED/CI.15200401/SB.WLA/V1/SB_WLA_HNZ_00_15200401.V1c"
prism_acc_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.acc.V2c"
prism_vel_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.vel.V2c"
prism_dis_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.dis.V2c"

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
# create copy of latest inv, remove date so it can be used in remove_response
get_inventory = True
if get_inventory:
    network = "SB"
    station = "WLA"
    inv_dir = "inventories"
    inv_path = os.path.join(os.getcwd(), inv_dir, network+"_"+station+".xml")
    if os.path.exists(inv_path):
        print("reading inv")
        inv = read_inventory(inv_path)
    else:
        print("downloading inv")
        rs_client = RS_Client("IRIS")
        #rs_client = RS_Client("RASPISHAKE")
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

st = read("/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/evid_15200401_13519/SB_WLA_HNZ_00_15200401.mseed")
target_channel="HNZ"
counts = st.select(channel=target_channel)[0]
counts.stats.units = "COUNTS"
metric_units = "ACC"

def convert_acc_to_vel_trace(tr):
    tr = tr.copy() # make a deepcopy to avoid altering original
    if tr.stats.units == "ACC":
        tr.integrate(method="cumtrapz")
        tr.detrend("demean") # results to velocity relative to mean
        tr.stats.units = "VEL"
    elif tr.stats.units == "VEL":
        tr.stats.units = "VEL"
    else:
        print("Can't convert", tr.stats.units, "to VEL.")

    return tr

def convert_vel_to_dis_trace(tr):
    tr = tr.copy() # make a deepcopy to avoid altering original
    if tr.stats.units == "VEL":
        tr.integrate(method="cumtrapz")
        tr.detrend("demean") # results to dis relative to mean
        tr.stats.units = "DIS"
    elif tr.stats.units == "DIS":
        tr.stats.units = "DIS"
    else:
        print("Can't convert", tr.stats.units, "to VEL.")

    return tr

def convert_counts_to_metric_trace(tr, metric_units):
    if tr.stats.units == "COUNTS":
        tr = tr.copy()
        freq = tr.stats.sampling_rate
        tr.remove_response(inventory=inv, pre_filt=[0.1, 0.5, 0.95*freq, freq],
                           output=metric_units, water_level=4.5, taper=False)
        tr.stats.units = metric_units
    return tr

counts = counts # see above selection from channel
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
fig, axs = plt.subplots(nrows=4, ncols=3, constrained_layout=True)

fig.suptitle(target_channel+" Channel in Metric Units")

# OBSPY PLOT COLUMN
axs[0][0].set_ylabel("OBSPY\nCounts")
axs[0][0].plot(counts.data)
x_text0 = len(counts.data)*0.67
y_text0 = max(counts.data)*0.7
axs[0][0].text(x_text0,y_text0,
            "Mean = {:+.6f}counts\nPeak = {:.6f}counts".format(counts_mean,counts_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[1][0].set_ylabel("OBSPY\nAcceleration")
axs[1][0].plot(acc.data)
x_text0 = len(acc.data)*0.67
y_text0 = max(acc.data)*0.7
axs[1][0].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(acc_mean,acc_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[2][0].set_ylabel("OBSPY\nVelocity")
axs[2][0].plot(vel.data)
x_text1 = len(vel.data)*0.7
y_text1 = max(vel.data)*0.7
axs[2][0].text(x_text1,y_text1,
            "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(vel_mean,vel_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[3][0].set_ylabel("OBSPY\nDisplacement")
axs[3][0].plot(dis.data)
x_text2 = len(dis.data)*0.7
y_text2 = max(dis.data)*0.7
axs[3][0].text(x_text2,y_text2,
            "Mean = {:+.6f}m\nPeak = {:.6f}m".format(dis_mean,dis_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

# PRISM PLOT COLUMN
axs[0][1].set_ylabel("PRISM V1\n(scaled counts)")
axs[0][1].plot(prism_v1.data)
x_text0 = len(prism_v1.data)*0.67
y_text0 = max(prism_v1.data)*0.7
axs[0][1].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_v1_mean,prism_v1_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[1][1].set_ylabel("PRISM V2\nAcceleration")
axs[1][1].plot(prism_acc)
x_text0 = len(prism_acc)*0.67
y_text0 = max(prism_acc)*0.7
axs[1][1].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_acc_mean,prism_acc_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[2][1].set_ylabel("PRISM V2\nVelocity")
axs[2][1].plot(prism_vel)
x_text0 = len(prism_vel)*0.67
y_text0 = max(prism_vel)*0.7
axs[2][1].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_vel_mean,prism_vel_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[3][1].set_ylabel("PRISM V2\nDisplacement")
axs[3][1].plot(prism_dis)
x_text0 = len(prism_dis)*0.67
y_text0 = max(prism_dis)*0.7
axs[3][1].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_dis_mean,prism_dis_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

# ABSOLUTE ERRORS COLUMN
axs[0][2].set_ylabel("OBSPYvsPRISM_V1\nAccel AbsErr")
axs[0][2].plot(abs(accv1_ae))
x_text2 = len(accv1_ae)*0.7
y_text2 = max(accv1_ae)*0.7
axs[0][2].text(x_text2,y_text2,
            "Mean = {:+.6f}m\nPeak = {:.6f}m".format(accv1_ae_mean,accv1_ae_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[1][2].set_ylabel("OBSPYvsPRISM_V2\nAccel AbsErr")
axs[1][2].plot(abs(acc_ae))
x_text2 = len(acc_ae)*0.7
y_text2 = max(acc_ae)*0.7
axs[1][2].text(x_text2,y_text2,
            "Mean = {:+.6f}m\nPeak = {:.6f}m".format(acc_ae_mean,acc_ae_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[2][2].set_ylabel("OBSPYvsPRISM_V2\nVelocity AbsErr")
axs[2][2].plot(abs(vel_ae))
x_text2 = len(vel_ae)*0.7
y_text2 = max(vel_ae)*0.7
axs[2][2].text(x_text2,y_text2,
            "Mean = {:+.6f}m\nPeak = {:.6f}m".format(vel_ae_mean,vel_ae_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[3][2].set_ylabel("OBSPYvsPRISM_V2\nVelocity AbsErr")
axs[3][2].plot(abs(dis_ae))
x_text2 = len(dis_ae)*0.7
y_text2 = max(dis_ae)*0.7
axs[3][2].text(x_text2,y_text2,
            "Mean = {:+.6f}m\nPeak = {:.6f}m".format(dis_ae_mean,dis_ae_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))


for row in axs:
    for ax in row:
        ax.grid()
plt.show()



