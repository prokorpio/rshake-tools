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

# parse cmd line args [DONE]
# create sub functions
# create automatic folder creation [DONE]
# fix print statements
# capture all channels
# handle dropped packets

#fig, axs = plt.subplots(2)
directory = os.path.join(os.getcwd(), "captures")

class PickWindow:
    """A window of rolling data where picks are computed as trigger
       for further data processing.
    """

    def __init__(self, sta_window, lta_window, stalta_window,
                 thresh_on, thresh_off, capture_buffer, max_evt_len):

        self.sta_window = sta_window
        self.lta_window = lta_window
        self.thresh_on = thresh_on
        self.thresh_off = thresh_off
        self.stalta_window = stalta_window # Must be >= settle_buff + max_evt_len + 2*capture_buffer
        self.capture_buffer = capture_buffer
        self.max_evt_len = max_evt_len

        self.active_window = lta_window + stalta_window

        self.rt_trace = RtTrace(max_length=self.active_window)
        self.pick_pairs_ind = np.empty((0,2), int)
        self.pick_pairs_is_processed = np.empty((0), bool)

        #fig, axs = plt.subplots(2)
        #self.fig = fig
        #self.axs = axs

    def roll_data(self, tr):
        print("Rolling data")
        freq = tr.stats.sampling_rate
        len_new_samples = len(tr.data)
        prev_trace_len = len(self.rt_trace.data)
        self.rt_trace.append(tr) # rolls automatically due to max_length

        shift_amount = (prev_trace_len + len_new_samples) - self.active_window*freq
        if (shift_amount > 0) and (len(self.pick_pairs_ind) > 0):
            self.pick_pairs_ind = self.pick_pairs_ind - shift_amount # shift indices via broadcast
            self.pick_pairs_ind = self.pick_pairs_ind[\
                    (self.pick_pairs_ind >= int(self.lta_window*freq)).all(axis=1)] # remove indices past stalta_window

    def calculate_new_picks(self, len_new_samples):
        print("Calculating new picks")
        print("    Calculating STA/LTA")
        freq = self.rt_trace.stats.sampling_rate
        sta_lta = self.rt_trace.copy().trigger("recstalta", \
                                            sta=self.sta_window, lta=self.lta_window)

        # Avoid spike at beginning, zero out misleading data
        if len(sta_lta) < int(self.lta_window*freq):
            sta_lta.data[:int(self.lta_window*freq)] = 0
            # After sta_lta becomes long than lta_window, this is automatically done
            # see recstalta.c by Moritz Beyreuther

        print("    Getting Pick Pairs")
        # Get part of sta/lta data that may contain tOn-tOff pick pair
        # Edge case is when a tOff is detected at first sample of new_samples
        data_len_to_pick = int(min((len(sta_lta),(self.max_evt_len*freq)+len_new_samples)))
        data_to_pick = sta_lta.data[-data_len_to_pick:]
        pick_pairs_ind_tmp = trigger_onset(data_to_pick,
                            self.thresh_on, self.thresh_off, max_len=self.max_evt_len*freq)

        # Handle re-detected pairs
        last_sample_ind = data_len_to_pick-1
        for pair_ind in pick_pairs_ind_tmp:

            # What happens if an event window is split between two batches?
            # tOff is taken to be at the last sample of the array despite
            # thresh_off condition not satisfied..
            # Skip such pairs.
            if (pair_ind[1] != last_sample_ind) \
                and data_to_pick[pair_ind[1]+1] < self.thresh_off:
                #print("temp tOff:", data_to_pick[pair_ind[1]+1])

                # Convert ind_tmp to ind on sta_lta
                pair_ind = pair_ind + (len(sta_lta) - data_len_to_pick)

                # Check if pair is not yet recorded
                # Having at least one of the pair mems match any of the values
                # in the pair-list is considered complete match to avoid
                # overlapping triggers (happens on ends of new_sta_lta window)
                # It is important that necessary_data_len doesn't intersect
                # with left_edge_ind so that shorter versions of outgoing
                # pick-pairs aren't re-added to pic_pairs_ind

                if not (self.pick_pairs_ind == pair_ind).any():
                    self.pick_pairs_ind = np.vstack((self.pick_pairs_ind, pair_ind))
                    self.pick_pairs_is_processed = np.append( \
                            self.pick_pairs_is_processed, False)

        print("Pick pairs len:")
        print(len(self.pick_pairs_ind))

        return sta_lta

    def get_processable_picks(self):
        print("Checking available picks")
        pick_pairs_utc = []
        if (len(self.pick_pairs_ind) > 0) and \
            not self.pick_pairs_is_processed.all():

            for i in np.where(~(self.pick_pairs_is_processed)):
                freq = self.rt_trace.stats.sampling_rate
                pair_ind = np.squeeze(self.pick_pairs_ind[i])
                ncapture_buffer = int(self.capture_buffer * freq)


                # condition for data capture
                if ((len(self.rt_trace.data) - 1) - pair_ind[1]) >= ncapture_buffer:
                    starttime = self.rt_trace.stats.starttime.timestamp
                    pair_ind_buffered = pair_ind + [-ncapture_buffer, ncapture_buffer]
                    pair_utc_buffered = starttime + pair_ind_buffered/freq
                    [start, end] = [UTCDateTime(t) for t in pair_utc_buffered]
                    pick_pairs_utc.append([start, end])
                    self.pick_pairs_is_processed[i] = True

        return pick_pairs_utc

def process_active_window(tr):
    #global ratio
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
        # After the condition isn't met anymore, this is automatically done
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

            # Check if pair is not yet recorded
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

    IP = "10.196.16.147"
    args = parser.parse_args()
    pickWindow = PickWindow(args.sta_window, args.lta_window, args.stalta_window,
                 args.thresh_on, args.thresh_off, args.capture_buffer, args.max_evt_len)
    fig, axs = plt.subplots(2)

    def process_data(tr):
        if tr.stats.channel != "EHZ":
            print(tr.stats.channel) # TODO: Buffer the other channels
            return

        pickWindow.roll_data(tr)
        sta_lta = pickWindow.calculate_new_picks(len(tr.data))
        pick_pairs_to_process = pickWindow.get_processable_picks()

        for pick_pair_utc in pick_pairs_to_process:
            print(pick_pair_utc)

        # plot data
        axs[0].plot(pickWindow.rt_trace.data)
        axs[1].plot(sta_lta)

        # plot triggers on ratio
        y_min, y_max = axs[1].get_ylim()
        axs[1].vlines(pickWindow.pick_pairs_ind[:,0], y_min, y_max, color='r', lw=2)
        axs[1].vlines(pickWindow.pick_pairs_ind[:,1], y_min, y_max, color='b', lw=2)

        fig.canvas.draw_idle()
        plt.pause(0.0001)
        axs[0].cla()
        axs[1].cla()

    client = create_client(IP, on_data=process_data)
    client.select_stream("AM", "RE722", "E??")
    client.run()
    print('done')










