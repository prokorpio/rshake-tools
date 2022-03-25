#--------------------------------------------------------------------
# Filename: STALTA-trigger.py
#  Purpose: Realtime data capture from Seedlink via recursive STA/LTA
#   Author: Christopher Jeff Sanchez
#
# Copyright (C) 2022 C.J. Sanchez
#---------------------------------------------------------------------

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import UTCDateTime
from obspy.realtime import RtTrace
from obspy.clients.seedlink.easyseedlink import create_client
from obspy.signal.trigger import trigger_onset
from collections import deque
from pathlib import Path

# parse cmd line args
# create sub functions
# create automatic folder creation
# fix print statements
# capture all channels

fig, axs = plt.subplots(2)
directory = os.path.join(os.getcwd(), "captures")
target_mseed_dir = "mseed_files"
target_figures_dir = "png_files"

sta_window = 2 # seconds
lta_window = 15 # should've been 80
capture_buffer = 15 # seconds
stalta_window = 120 # must be >= settle_buff + max_evt_len + 2*capture_buffer
active_window = lta_window + stalta_window
thresh_on = 3
thresh_off = 1.5
sensor_freq=100 # RS4D has 100Hz sensor
max_evt_len = 30 # Check longest possible event occurence (in sec)

# lta_window is the required buffer to the left to compute stalta
rt_trace = RtTrace(max_length=active_window)

pick_pairs_ind = np.empty((0,2), int)
pick_pairs_to_save = np.empty((0,2), int)


def append_trace_to_realtime(tr):
    global ratio
    global pick_pairs_ind
    global pick_pairs_to_save

    freq = tr.stats.sampling_rate
    #print("Appending the following trace:")
    #print(tr)
    prev_trace_len = len(rt_trace.data)
    rt_trace.append(tr)

    print("RtTrace:")
    print(rt_trace)

    print("Calculating STA/LTA")
    sta_lta = rt_trace.copy().trigger("recstalta", \
                                        sta=sta_window, lta=lta_window)
    len_new_samples = len(tr.data)
    print("Len new samples:", len_new_samples)

    # Avoid spike at beginning, zero out misleading data
    nlta_window = int(lta_window*freq)
    if len(sta_lta) < nlta_window:
        sta_lta.data[:nlta_window] = 0
        # after the condition isn't met anymore, this is automatically done
        # see recstalta.c by Moritz Beyreuther

    # roll pic_pairs once sta_lta values starts rolling
    shift_amount = (prev_trace_len + len_new_samples) - active_window*freq
    if (len(pick_pairs_ind) > 0) and (shift_amount > 0):
        pick_pairs_ind = pick_pairs_ind - shift_amount # shift indices
        pick_pairs_ind = pick_pairs_ind[(pick_pairs_ind>=nlta_window).all(axis=1)] # remove negative indices

    print("Checking for data capture..")
    ncapture_buffer = capture_buffer*freq
    print("Pick pairs to save len:")
    print(len(pick_pairs_to_save))
    if (len(pick_pairs_to_save) > 0):
        if shift_amount > 0:
            pick_pairs_to_save = pick_pairs_to_save - shift_amount # shift indices

        for pair_ind in pick_pairs_to_save:
            print("testing:", pair_ind)
            print("with right dist of:", len(sta_lta) - pair_ind[1])
            if (len(sta_lta) - pair_ind[1] >= ncapture_buffer):
                pair_ind_buffered = pair_ind + [-ncapture_buffer, ncapture_buffer]
                pair_utc_buffered = rt_trace.stats.starttime.timestamp + pair_ind_buffered/freq
                [start, end] = [UTCDateTime(t) for t in pair_utc_buffered]

                print("Saving to:")
                print(start, end)
                mseed_to_save = rt_trace.slice(start,end)
                event_name = "_".join([rt_trace.stats.network, \
                                      rt_trace.stats.station,\
                            start.strftime("%y-%m-%dT%H:%M:%S")])
                target_dir = os.path.join(directory, event_name)
                Path(target_dir).mkdir(parents=True, exist_ok=True)
                mseed_path = os.path.join(target_dir, event_name +".mseed")
                plot_path = os.path.join(target_dir, event_name +".png")
                mseed_to_save.write(mseed_path, format="MSEED", reclen=512)
                mseed_to_save.plot(outfile=plot_path)
                pick_pairs_to_save = pick_pairs_to_save[1:]

    print("Getting trigger times")
    # What happens if an event window is split between two batches?
        # tOff is taken to be at the last sample of the array.

    # get part of sta/lta data that may contain tOn,tOff pair
    # edge case is when a tOff is detected at first sample of new_samples
    necessary_data_len = min((len(sta_lta),(max_evt_len*freq)+len_new_samples))
    new_sta_lta_left_buffered = sta_lta.data[-necessary_data_len:]
    pick_pairs_ind_tmp = trigger_onset(new_sta_lta_left_buffered,
                        thresh_on, thresh_off, max_len=max_evt_len*freq)

    if len(pick_pairs_ind_tmp) >0:
        print("temp picks")
        #print(np.diff(pick_pairs_ind_tmp, axis = 1))
        print(pick_pairs_ind_tmp)
    else:
        print("temp picks")
        print("[]")

    # handle re-detected pairs
    last_sample_ind = len(new_sta_lta_left_buffered)-1
    for pair_ind in pick_pairs_ind_tmp:
        # skip if tOff was only considered bc its the end of batch
            # and not because it's actually <thresh_off
        if (pair_ind[1] != last_sample_ind) \
            and new_sta_lta_left_buffered[pair_ind[1]]+1 > thresh_off:

            # convert ind_tmp to ind on sta_lta
            pair_ind = pair_ind + (len(sta_lta) - necessary_data_len)
            print("temp picks adjusted:")
            print(pair_ind)
            # check if pair is not yet recorded
            # Having at least one of the pair mems match any of the values
            # in the pair-list is considered complete match to avoid
            # overlapping triggers (happens on ends of new_sta_lta window)
            # It is important that necessary_data_len doesn't intersect
            # with left_edge_ind so that shorter versions of outgoing
            # pick-pairs aren't re-added to pic_pairs_ind
            if not (pick_pairs_ind == pair_ind).any():
                pick_pairs_ind = np.vstack((pick_pairs_ind, pair_ind))
                pick_pairs_to_save = np.vstack((pick_pairs_to_save, pair_ind))

    print("Pick pairs len:")
    print(len(pick_pairs_ind))


    print()

    # plot data
    axs[0].plot(rt_trace.data)
    axs[1].plot(sta_lta)

    # plot triggers on ratio
    y_min, y_max = axs[1].get_ylim()
    axs[1].vlines(pick_pairs_ind[:,0], y_min, y_max, color='r', lw=2)
    axs[1].vlines(pick_pairs_ind[:,1], y_min, y_max, color='b', lw=2)

    fig.canvas.draw_idle()
    plt.pause(0.0001)
    axs[0].cla()
    axs[1].cla()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Capture data from SeedLink server using STA/LTA")
    parser.add_argument("--sta_window", type=int, default=2,
        help="Length in seconds used for the short-term-average (STA)")
    parser.add_argument("--lta_window", type=int, default=15,
        help="Length in seconds used for the long-term-average (LTA)")
    parser.add_argument("--thresh_on", type=float, default=3.0,
        help="Upper threshold to create On-pick and activate data capture.")
    parser.add_argument("--thresh_off", type=float, default=1.5,
        help="Lower threshold to create Off-pick and finish data capture.")
    parser.add_argument("--stalta_window", type=int, default=120,
        help="Length in seconds where STA/LTA computation should have nonzero values")
    parser.add_argument("--capture_buffer", type=int, default=10,
        help="Buffer in seconds applied to left and right of data capture.")
    parser.add_argument("--max_evt_len", type=int, default=30,
        help="Length in seconds of maximum data capture before adding buffer.")
    parser.add_argument("--directory", type=str, default="captures",
        help="Target directory for outputs")

    args = parser.parse_args()
    print(args)
    print('done')
    #client = create_client("10.196.16.147", on_data=append_trace_to_realtime)
    #client.select_stream("AM", "RE722", "EHZ")
    #client.run()










