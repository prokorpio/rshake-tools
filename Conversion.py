from RollingStream import RollingStream
from obspy.core import Stream
from helpers import (
    get_inventory,
    g_to_intensity,
    intensity_to_int,
    save_mseed,
    save_json,
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
                print(self.thread_name, " RECEIVED: PICK", message['data'], " qsize: ", picks_queue.qsize(),sep='')

                pick = message['data']
                # NOTE: Only trigger processing base on picks from basis_channel
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
                    max_pga_channel = "" # channel with max PGA
                    for acc_tr in acc_st:
                        if acc_tr.stats.peak > max_pga:
                            max_pga = acc_tr.stats.peak
                            max_intensity_str = g_to_intensity(max_pga/9.81)
                            max_intensity_int = intensity_to_int(max_intensity_str)
                            max_pga_channel = acc_tr.stats.channel

                    ## check against thresh
                    if max_intensity_int >= self.intensity_threshold:
                        # since threshold is reached, convert all to dis to get PGD
                        vel_st = self.convert_metric_to_vel(metric_st)
                        vel_st = self.set_stats_peak(vel_st)
                        dis_st = self.convert_vel_to_dis(vel_st)
                        dis_st = self.set_stats_peak(dis_st)

                        # get PGD
                        max_pgd = 0
                        for dis_tr in acc_st:
                            if dis_tr.stats.peak > max_pgd:
                                max_pgd = dis_tr.stats.peak
                                max_pgd_channel = dis_tr.stats.channel

                        # publish event summary
                        evt_summary = {
                            "PGA": max_pga,
                            "PGA_channel": max_pga_channel,
                            "intensity": max_intensity_str,
                            "PGD": max_pgd,
                            "PGD_channel": max_pgd_channel,
                        }
                        self.message_board.publish(self.dst_event_summary_topic, evt_summary)

                        ## create all event data
                        event_name = "_".join([self.network, self.station,
                            pick["on_time"].strftime("%y-%m-%dT%H:%M:%S")])
                        parent_dir = os.path.join(os.getcwd(), self.captures_dir)
                        dir_path = os.path.join(parent_dir, event_name)

                        # create summary report
                        report = {
                            "Event name": event_name,
                            "Detection time": str(pick["on_time"]),
                            "Intensity": max_intensity_str,
                            "PGA": max_pga,
                            "PGA_channel": max_pga_channel,
                            "PGD": max_pgd,
                            "PGD_channel": max_pgd_channel,
                            "Channel Peak Accel's (m/s*s)": \
                                {a.stats.channel:a.stats.peak for a in acc_st},
                            "Channel Peak Velo's (m/s)": \
                                {v.stats.channel:v.stats.peak for v in vel_st},
                            "Channel Peak Disp's (m)": \
                                {d.stats.channel:d.stats.peak for d in dis_st},
                        }

                        # write all to files
                        _ = save_json(report, "event_summary", dir_path)
                        _ = save_mseed(Stream(traces), "counts", dir_path)
                        _ = save_mseed(metric_st, "metric", dir_path)
                        _ = save_mseed(acc_st, "acc", dir_path)
                        _ = save_mseed(vel_st, "vel", dir_path)
                        _ = save_mseed(dis_st, "dis", dir_path)

                        # publish location of all data
                        evt_data = {"path": dir_path}
                        self.message_board.publish(self.dst_event_data_topic, evt_data)
