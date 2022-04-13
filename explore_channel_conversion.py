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
a = get_PRISM_data("/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.acc.V2c")
print(a)

# PRISM-processed-data import
prism_acc_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.acc.V2c"
prism_vel_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.vel.V2c"
prism_dis_file = "/home/jeffsanchez/Downloads/cosmos_downloads/Mag8.2/prismPROCESSED/CI.15200401/SB.WLA/V2/SB_WLA_HNZ_00_15200401.dis.V2c"
prism_acc = get_PRISM_data(prism_acc_file)/100 # original data in cm
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

acc = convert_counts_to_metric_trace(counts, metric_units)
acc_mean = acc.data.mean()
acc_peak = max(abs(acc.data))
vel = convert_acc_to_vel_trace(acc)
vel_mean = vel.data.mean()
vel_peak = max(abs(vel.data))
dis = convert_vel_to_dis_trace(vel)
dis_mean = dis.data.mean()
dis_peak = max(abs(dis.data))


# Plot results
fig, axs = plt.subplots(nrows=3, ncols=2, constrained_layout=True)

fig.suptitle(target_channel+" Channel in Metric Units")

axs[0][0].set_ylabel("OBSPY Acceleration")
axs[0][0].plot(acc.data)
x_text0 = len(acc.data)*0.67
y_text0 = max(acc.data)*0.7
axs[0][0].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(acc_mean,acc_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[0][1].set_ylabel("PRISM Acceleration")
axs[0][1].plot(prism_acc)
x_text0 = len(prism_acc)*0.67
y_text0 = max(prism_acc)*0.7
axs[0][1].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_acc_mean,prism_acc_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[1][0].set_ylabel("Velocity")
axs[1][0].plot(vel.data)
x_text1 = len(vel.data)*0.7
y_text1 = max(vel.data)*0.7
axs[1][0].text(x_text1,y_text1,
            "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(vel_mean,vel_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[1][1].set_ylabel("PRISM Velocity")
axs[1][1].plot(prism_vel)
x_text0 = len(prism_vel)*0.67
y_text0 = max(prism_vel)*0.7
axs[1][1].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_vel_mean,prism_vel_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[2][0].set_ylabel("Displacement")
axs[2][0].plot(dis.data)
x_text2 = len(dis.data)*0.7
y_text2 = max(dis.data)*0.7
axs[2][0].text(x_text2,y_text2,
            "Mean = {:+.6f}m\nPeak = {:.6f}m".format(dis_mean,dis_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[2][1].set_ylabel("PRISM Displacement")
axs[2][1].plot(prism_dis)
x_text0 = len(prism_dis)*0.67
y_text0 = max(prism_dis)*0.7
axs[2][1].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(prism_dis_mean,prism_dis_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

for row in axs:
    for ax in row:
        ax.grid()
plt.show()



