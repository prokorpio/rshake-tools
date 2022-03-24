from obspy.clients.seedlink.easyseedlink import create_client
from obspy.realtime import RtTrace
from obspy.signal.trigger import trigger_onset
from collections import deque
import numpy as np
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2)

# derived from SC3's scautopick default config
sta_window = 2 # seconds
lta_window = 10 #80
ratio_window = 120 # must be > lta_window, and >= max_evt_len + 2*buffer
thres_on = 3
thres_off = 1.5
sensor_freq=100 # RS4D has 100Hz sensor
max_evt_len = 80*sensor_freq # Check longest possible event occurence (in sec)

rt_trace = RtTrace(max_length=ratio_window)

ratio = np.zeros(ratio_window*sensor_freq)
pick_pairs_ind = np.empty((0,2), int) # len is determined by diff of first & last
                              # elem should be within ratio_window*sensor_freq

def append_trace_to_realtime(tr):
    global ratio
    global pick_pairs_ind

    #print("Appending the following trace:")
    #print(tr)
    rt_trace.append(tr)

    print("RtTrace:")
    print(rt_trace)

    print("Calculating STA/LTA")
    sta_lta = rt_trace.copy().trigger("recstalta", sta=sta_window, lta=lta_window)
    len_new_samples = len(tr.data)
    print("Len new samples:", len_new_samples)
    # roll ratio
    ratio = np.roll(ratio,-len_new_samples)
    ratio[-len_new_samples:] = sta_lta[-len_new_samples:]
    print("Len ratio:", len(ratio))
    # roll pic_pairs
    if len(pick_pairs_ind) > 0:
        pick_pairs_ind = pick_pairs_ind - len_new_samples # shift indices
        pick_pairs_ind = pick_pairs_ind[(pick_pairs_ind>0).all(axis=1)] # remove negative indices

    print("Getting trigger times")
    # What happens if an event window is split between two batches?
        # tOff is taken to be at the last sample of the array.

    # get part of sta/lta data that may contain tOn,tOff pair
    # edge case is when a tOff is detected at first sample of new_samples
    new_sta_lta_left_buffered = ratio[-(max_evt_len + len_new_samples):]
    pick_pairs_ind_tmp = trigger_onset(new_sta_lta_left_buffered,
                        thres_on, thres_off, max_len=max_evt_len)

    if len(pick_pairs_ind_tmp) >0:
        print("temp picks")
        print(np.diff(pick_pairs_ind_tmp, axis = 1))
    else:
        print("temp picks")
        print("[]")

    # handle re-detected pairs
    last_sample_ind = len(new_sta_lta_left_buffered)-1
    for pair_ind in pick_pairs_ind_tmp:
        # skip if tOff was only considered bc its the end of batch
            # and not because it's actually <thres_off
        if (pair_ind[1] != last_sample_ind) \
            and new_sta_lta_left_buffered[pair_ind[1]]+1 > thres_off:

            # convert ind_tmp to ind on ratio
            pair_ind = pair_ind + len(ratio)- len(new_sta_lta_left_buffered)

            # check if pair is not yet recorded
            # Having at least one of the pair mems match any of the values
            # in the pair-list is considered complete match to avoid
            # overlapping triggers (happens on ends of new_sta_lta window)
            if not (pick_pairs_ind == pair_ind).any():
                pick_pairs_ind = np.vstack((pick_pairs_ind, pair_ind))

    print("Pick pairs len:")
    print(len(pick_pairs_ind))

    print()

    # plot data
    axs[0].plot(rt_trace.data)
    axs[1].plot(ratio)

    # plot triggers on ratio
    y_min, y_max = axs[1].get_ylim()
    axs[1].vlines(pick_pairs_ind[:,0], y_min, y_max, color='r', lw=2)
    axs[1].vlines(pick_pairs_ind[:,1], y_min, y_max, color='b', lw=2)

    fig.canvas.draw_idle()
    plt.pause(0.0001)
    axs[0].cla()
    axs[1].cla()

client = create_client("10.196.16.147", on_data=append_trace_to_realtime)
client.select_stream("AM", "RE722", "EHZ")
client.run()
