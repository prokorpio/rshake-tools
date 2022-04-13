#--------------------------------------------------------------------
# Filename: STALTA-trigger.py
#  Purpose: Realtime data capture from Seedlink via recursive STA/LTA
#   Author: Christopher Jeff Sanchez
#
# Copyright (C) 2022 C.J. Sanchez
#---------------------------------------------------------------------

import os
import argparse
from collections import deque
from threading import Thread
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from obspy import read_inventory
from obspy.core import UTCDateTime
from obspy.realtime import RtTrace
from obspy.clients.seedlink.easyseedlink import create_client
from obspy.clients.seedlink.basic_client import Client as BasicSLClient
from obspy.clients.fdsn import Client as RS_Client
from obspy.signal.trigger import trigger_onset
from obspy.core.inventory import Inventory, Network

# TODO:
# parse cmd line args [DONE]
# create sub functions [DONE]
# create automatic folder creation [DONE]
# capture all channels [DONE]
# fix print statements
# handle dropped packets

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

def save(st, event_name, title, save_img=True, save_str=True):
    directory = os.path.join(os.getcwd(), "captures")
    target_dir = os.path.join(directory, event_name)
    Path(target_dir).mkdir(parents=True, exist_ok=True)

    if save_str:
        mseed_path = os.path.join(target_dir, title +".mseed")
        st.write(mseed_path, format="MSEED", reclen=512)

    if save_img:
        plot_path = os.path.join(target_dir, title +".png")
        st.plot(outfile=plot_path)

def convert_counts_to_metric(st, inv):
    stream = st.copy() # make a deepcopy to avoid altering original
    vel_channels = ["EHE", "EHN", "EHZ", "SHZ"]
    #vel_channels = ["EHE", "EHN", "EHZ", "SHZ", "HHZ"]
    acc_channels = ["ENE", "ENN", "ENZ"]

    for tr in stream:
        units = "COUNTS"
        if tr.stats.channel in vel_channels:
            units = "VEL"
        elif tr.stats.channel in acc_channels:
            units = "ACC"

        if units != "COUNTS":
            freq = tr.stats.sampling_rate
            tr.remove_response(inventory=inv, pre_filt=[0.1, 0.5, 0.95*freq, freq],
                               output=units, water_level=4.5, taper=False)

        tr.stats.units = units

    return stream

def convert_metric_to_disp(st):
    stream = st.copy() # make a deepcopy to avoid altering original
    for tr in stream:
        if tr.stats.units == "VEL":
            tr.integrate(method="cumtrapz")
            tr.stats.units = "DISP"
        elif tr.stats.units == "ACC":
            tr.integrate(method="cumtrapz")
            tr.integrate(method="cumtrapz")
            tr.stats.units = "DISP"
        elif tr.stats.units == "DISP":
            continue
        else:
            print("Can't convert", tr.stats.units, "to DISP.")

    return stream


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Capture data from SeedLink server using STA/LTA")
    parser.add_argument("--IP", type=str, default="rs.local", required=True,
        help="Host address of the raspberryshake unit")
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

    # get station information
    basic_sl_client = BasicSLClient(args.IP)
    channels = basic_sl_client.get_info(level="channel")
    for i, channel in enumerate(channels):
        if ("HZ" in channel[3]):
            network = channel[0]
            station = channel[1]
            basis_channel = channel[3]
            break
        elif i == len(channels)-1:
            raise RuntimeError("*HZ channel not available in rshake@"+args.IP)
    #network = "GE"
    #station = "WLF"
    #basis_channel = "HHZ"

    # create copy of latest inv, remove date so it can be used in remove_response
    inv_dir = "inventories"
    inv_path = os.path.join(os.getcwd(), inv_dir, network+"_"+station+".xml")
    if os.path.exists(inv_path):
        print("reading inv")
        inv = read_inventory(inv_path)
    else:
        print("downloading inv")
        #rs_client = RS_Client("IRIS")
        rs_client = RS_Client("RASPISHAKE")
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

    # select realtime stream
    rt_sl_client = create_client(args.IP) # client for realtime sl streaming
    rt_sl_client.select_stream(network, station, basis_channel)

    pickWindow = PickWindow(args.sta_window, args.lta_window, args.stalta_window,
                 args.thresh_on, args.thresh_off, args.capture_buffer, args.max_evt_len)

    fig, axs = plt.subplots(2) # for visual checking

    def download_convert_save(start, end):
        global basic_sl_client
        global network
        global station
        global inv

        event_name = "_".join([network, station,
            start.strftime("%y-%m-%dT%H:%M:%S")])
        print("THREAD:", event_name)
        #st = basic_sl_client.get_waveforms(network, station, "*", "HHZ", start, end)
        st = basic_sl_client.get_waveforms(network, station, "*", "*", start, end)
        print("THREAD:","Downloaded ST")
        save(st, event_name, "counts")
        print("THREAD:","Saved ST")
        st_metric  = convert_counts_to_metric(st, inv)
        print("THREAD:","Converted ST to metric")
        save(st_metric, event_name, "metric")
        print("THREAD:","Saved ST in metric")
        st_disp = convert_metric_to_disp(st_metric)
        print("THREAD:","Converted ST to DISP")
        save(st_disp, event_name, "DISP")
        print("THREAD:","Saved ST in DISP")

    counter_for_testing = 0
    processing_fns = [download_convert_save]
    def on_data_callback(tr):
        global counter_for_testing

        pickWindow.roll_data(tr)
        sta_lta = pickWindow.calculate_new_picks(len(tr.data))
        pick_pairs_to_process = pickWindow.get_processable_picks()

        #counter_for_testing += 1
        #print("COUNTER:", counter_for_testing)
        #if counter_for_testing == 10:
        #    #counter_for_testing = 0
        #    print("Processing stuff...")
        #    for fn in processing_fns:
        #        start = UTCDateTime()-5-1
        #        end = UTCDateTime()-1
        #        Thread(target=fn, args=(start, end)).start()

        for (start, end) in pick_pairs_to_process:
            print("Processing stuff...")
            for fn in processing_fns:
               Thread(target=fn, args=(start, end)).start()

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

    rt_sl_client.on_data = on_data_callback
    rt_sl_client.run()










