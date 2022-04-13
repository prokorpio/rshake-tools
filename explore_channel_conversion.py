from obspy.core import read
import matplotlib.pyplot as plt

st = read("captures/AM_RE722_22-04-13T04:44:57/metric.mseed")
target_channel="ENE"
acc = st.select(channel=target_channel)[0]
acc.stats.units = "ACC"

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

def convert_vel_to_disp_trace(tr):
    tr = tr.copy() # make a deepcopy to avoid altering original
    if tr.stats.units == "VEL":
        tr.integrate(method="cumtrapz")
        tr.detrend("demean") # results to disp relative to mean
        tr.stats.units = "DISP"
    elif tr.stats.units == "DISP":
        tr.stats.units = "VEL"
    else:
        print("Can't convert", tr.stats.units, "to VEL.")

    return tr

acc_mean = acc.data.mean()
acc_peak = max(abs(acc.data))
vel = convert_acc_to_vel_trace(acc)
vel_mean = vel.data.mean()
vel_peak = max(abs(vel.data))
disp = convert_vel_to_disp_trace(vel)
disp_mean = disp.data.mean()
disp_peak = max(abs(disp.data))

fig, axs = plt.subplots(nrows=3, ncols=1, constrained_layout=True)

fig.suptitle(target_channel+" Channel in Metric Units")

axs[0].set_ylabel("Acceleration")
axs[0].plot(acc.data)
x_text0 = len(acc.data)*0.67
y_text0 = max(acc.data)*0.7
axs[0].text(x_text0,y_text0,
            "Mean = {:+.6f}m/s2\nPeak = {:.6f}m/s2".format(acc_mean,acc_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[1].set_ylabel("Velocity")
axs[1].plot(vel.data)
x_text1 = len(vel.data)*0.7
y_text1 = max(vel.data)*0.7
axs[1].text(x_text1,y_text1,
            "Mean = {:+.6f}m/s\nPeak = {:.6f}m/s".format(vel_mean,vel_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

axs[2].set_ylabel("Displacement")
axs[2].plot(disp.data)
x_text2 = len(disp.data)*0.7
y_text2 = max(disp.data)*0.7
axs[2].text(x_text2,y_text2,
            "Mean = {:+.6f}m\nPeak = {:.6f}m".format(disp_mean,disp_peak),
            ha="left", va="center",
            bbox=dict(boxstyle="round", fc="white"))

for ax in axs:
    ax.grid()
plt.show()



