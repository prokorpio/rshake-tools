from RollingStream import RollingStream
from obspy.core import Stream
from helpers import (
    get_inventory,
    g_to_intensity,
    intensity_to_int,
    save_mseed,
)
from trace_conversion import (
    convert_counts_to_metric_trace,
    convert_acc_to_vel_trace,
    convert_vel_to_acc_trace,
    convert_vel_to_dis_trace,
)
import threading
import os


class Conversion(threading.Thread):
    """ Class to convert stream to disp on receipt of PICKS """

    def __init__(self, message_board, \
                 network, station, channels, sps, \
                 basis_channel, left_padding_duration, intensity_threshold, \
                 duration, \
                 inv_dir, captures_dir):

        # prepare thread info
        self.thread_name = self.__class__.__name__ + ":" + station
        threading.Thread.__init__(self, name=self.thread_name)

        # prepare station info
        self.network = network
        self.station = station
        self.channels = channels
        self.sps = sps
        self.inv = get_inventory(inv_dir, network, station)
        self.captures_dir = captures_dir

        # prepare source of realtime rolling data
        self.duration = duration
        self.channel_length = duration*sps
        self.rt_stream = RollingStream(self.channels, self.channel_length)

        # prepare conversion attribs
        self.basis_channel = basis_channel
        self.left_padding_duration = left_padding_duration
        self.intensity_threshold = intensity_threshold

        # prepare pub/sub comms
        self.src_traces_topic = station.upper() + "_" + "TRACES"
        self.src_picks_topic = "PICKS"
        self.dst_event_summary_topic = "EVENT_SUMMARY"
        self.dst_event_data_topic = "EVENT_DATA"
        self.message_board = message_board

    def convert_counts_to_metric(self, traces):
        """ Instrument response correction & bandpass filtering """

        vel_channels = ["EHE", "EHN", "EHZ", "SHZ"]
        acc_channels = ["ENE", "ENN", "ENZ"]

        metric_traces = Stream()
        for tr in traces:
            # set units to counts first
            tr.stats.units = "COUNTS"

            # get correct metric units
            if tr.stats.channel in vel_channels:
                units = "VEL"
            elif tr.stats.channel in acc_channels:
                units = "ACC"
            #else: raise Error/Warning

            # convert
            tr.attach_response(self.inv) # response file used in conversion
            metric_tr, _, _ = convert_counts_to_metric_trace(tr, units) # correction & filter
            metric_tr.stats.peak = max(abs(metric_tr.data))
            metric_traces.append(metric_tr)

        return metric_traces

    def convert_metric_to_acc(self, metric_traces):
        """ Convert traces with METRIC units into ACC units"""
        acc_traces = Stream()

        for tr in metric_traces:
            if tr.stats.units == "ACC":
                acc_traces.append(tr)
            elif tr.stats.units == "VEL":
                tr = convert_vel_to_acc_trace(tr)
                acc_traces.append(tr)

        return acc_traces

    def convert_metric_to_vel(self, metric_traces):
        """ Convert traces with METRIC units into VEL units"""
        vel_traces = Stream()

        for tr in metric_traces:
            if tr.stats.units == "ACC":
                tr = convert_acc_to_vel_trace(tr)
                vel_traces.append(tr)
            elif tr.stats.units == "VEL":
                vel_traces.append(tr)

        return vel_traces

    def convert_vel_to_dis(self, vel_traces):
        """ Convert traces with VEL units into DISP units"""
        dis_traces = Stream()

        for tr in vel_traces:
            if tr.stats.units == "VEL":
                tr = convert_vel_to_dis_trace(tr)
                dis_traces.append(tr)
            elif tr.stats.units == "DISP":
                dis_traces.append(tr)

        return dis_traces

    def set_stats_peak(self, traces):
        st = Stream()
        for tr in traces:
            tr.stats.peak = max(abs(tr.data))
            st.append(tr)
        return st

    def run(self):
        traces_queue = self.message_board.subscribe(self.src_traces_topic)
        picks_queue = self.message_board.subscribe(self.src_picks_topic)
        pick_msgs = []

        while True:
            # update rt_stream, assumes this wouldn't hog the loop
            for trace_msg in traces_queue.listen():
                print(self.thread_name, " RECEIVED: ", trace_msg['data'], " qsize: ", traces_queue.qsize(),sep='')

                additional_tr = trace_msg['data']
                _ = self.rt_stream.update_trace(additional_tr)

                if traces_queue.qsize() == 0:
                    break

            # process one pick, ensures rt_stream doesn't fill up
            tmp = list(picks_queue.listen(block=False))
            if len(tmp) > 0:
                pick_msgs.append(*tmp)
            if len(pick_msgs) > 0:
                message = pick_msgs.pop(0)
                print(self.thread_name, " RECEIVED: ", message['data'], " qsize: ", picks_queue.qsize(),sep='')

                pick = message['data']
                # NOTE: Only trigger base on picks from basis_channel
                #       to avoid picks from same event triggering multiple
                #       conversions. However, max pga is taken across all
                #       channels.
                if pick["channel"] == self.basis_channel:
                    # get acc traces
                    traces = self.rt_stream.get_traces_from(pick["on_time"]-self.left_padding_duration) # NOTE: this assumes the on_time is already in rt_stream
                    metric_st = self.convert_counts_to_metric(traces) # returns Stream object
                    metric_st = self.set_stats_peak(metric_st)
                    acc_st = self.convert_metric_to_acc(metric_st)
                    acc_st = self.set_stats_peak(acc_st)

                    # get max PGA, max intensity, and which channel
                    max_pga = 0
                    max_intensity_int = 0
                    max_intensity_str = ""
                    max_channel = "" # channel with max PGA
                    for acc_tr in acc_st:
                        if acc_tr.stats.peak > max_pga:
                            max_pga = acc_tr.stats.peak
                            max_intensity_str = g_to_intensity(max_pga/9.81)
                            max_intensity_int = intensity_to_int(max_intensity_str)
                            max_channel = acc_tr.stats.channel

                    # check against thresh, and get PGD if thresh is reached
                    if max_intensity_int >= self.intensity_threshold:
                        if max_channel == "EHZ":
                            vel_tr = metric_st.select(channel=max_channel)[0]
                        else:
                            acc_tr = acc_st.select(channel=max_channel)[0]
                            vel_tr = convert_acc_to_vel_trace(acc_tr)
                        dis_tr = convert_vel_to_dis_trace(vel_tr)
                        peak_dis = max(abs(dis_tr.data))*100 #convert to centimeters

                        # publish event summary
                        evt_summary = {
                                "PGA": max_pga,
                                "PGD": peak_dis, # in centimeters
                                "intensity": max_intensity_str,
                                "channel": max_channel,
                        }
                        self.message_board.publish(self.dst_event_summary_topic, evt_summary)

                        # convert the rest of the channels to disp
                        vel_st = self.convert_metric_to_vel(metric_st)
                        vel_st = self.set_stats_peak(vel_st)
                        dis_st = self.convert_vel_to_dis(vel_st)
                        dis_st = self.set_stats_peak(dis_st)

                        # save data and publish event data
                        event_name = "_".join([self.network, self.station,
                            pick["on_time"].strftime("%y-%m-%dT%H:%M:%S")])
                        dir_path = save_mseed(Stream(traces), event_name, "counts", self.captures_dir)
                        dir_path = save_mseed(metric_st, event_name, "metric", self.captures_dir)
                        dir_path = save_mseed(acc_st, event_name, "acc", self.captures_dir)
                        dir_path = save_mseed(vel_st, event_name, "vel", self.captures_dir)
                        dir_path = save_mseed(dis_st, event_name, "dis", self.captures_dir)

                        evt_data = {"path": dir_path}
                        self.message_board.publish(self.dst_event_data_topic, evt_data)
