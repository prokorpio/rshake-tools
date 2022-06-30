import threading
from obspy.signal.trigger import trigger_onset
from obspy import UTCDateTime

from RollingStream import RollingStream

class DetectSpike(threading.Thread):
    """ Class to detect spikes from STATION_TRACES and publish to PICKS """

    def __init__(self, message_board, \
                 station, channels, sps, \
                 sta_sec, lta_sec, on_thresh, off_thresh, max_evt_dur,\
                 duration):
        """
            Parameters:
            - message_board: system PubSub board to receive mseed packets from
            - station: name of station whose traces are to be processed
            - channels: list of channel names to be processed
            - sps: station sampling rate / samples per second
            - sta_sec: length of data in seconds to get short-term-ave
            - lta_sec: length of data in seconds to get long-term-ave
            - on_thresh: sta/lta threshold to trigger an on-pick
            - off_thresh: sta/lta threshold to trigger an off-pick
            - max_evt_dur: assumed max length in seconds of on-to-off pick-pair
            - duration: max_length in seconds of traces to be processed
        """

        # prepare thread info
        self.thread_name = self.__class__.__name__ + ":" + station
        threading.Thread.__init__(self, name=self.thread_name)

        # prepare station info
        self.station = station
        self.channels = channels
        self.sps = sps

        # prepare source of realtime rolling data
        self.duration = duration
        self.channel_length = duration*sps
        self.rt_stream = RollingStream(self.channels, self.channel_length)

        # prepare STA/LTA params
        self.sta_sec, self.sta_len = sta_sec, sta_sec*sps
        self.lta_sec, self.lta_len = lta_sec, lta_sec*sps
        self.max_evt_dur, self.max_evt_len = max_evt_dur, max_evt_dur*sps
        self.on_thresh = on_thresh
        self.off_thresh = off_thresh

        # prepare pub/sub comms
        self.src_topic = station.upper() + "_" + "TRACES"
        self.dst_topic = "PICKS"
        self.message_board = message_board

    def _get_STALTA_pick_pairs(self, tr, prev_pick):
        """ Process a trace and return any new pick-pairs

            Parameters:
            - tr: trace to get picks from
            - prev_pick: channel and UTCDateTime of previously computed pick
        """

        # len(tr) must allow for calculation of pick-pair with dur=max_evt_dur
        if len(tr) < self.lta_len + self.max_evt_dur:
            return None

        sta_lta = tr.trigger("recstalta", sta=self.sta_sec, lta=self.lta_sec)

        # NOTE: Optimize by acting immediately on on-pick?
        # Only consider ratio of len=max_evt_len since any evt longer
        # than this won't be detected anyway
        data_len_to_pick = self.max_evt_len + 5 # 5 is just safety padding
        data_to_pick = sta_lta[-data_len_to_pick:]
        pick_pair_ind_tmp = trigger_onset(data_to_pick,
                            self.on_thresh, self.off_thresh, self.max_evt_len)

        # What happens if an event window is split between two batches?
        # tOff is taken to be at the last sample of the array despite
        # thresh_off condition not satisfied..
        # Skip such pairs.
        if len(pick_pair_ind_tmp) > 0:
            latest_pair = pick_pair_ind_tmp[-1] # only process latest
            last_sample_ind = data_len_to_pick - 1
            if (latest_pair[1] < last_sample_ind) \
                and (data_to_pick[latest_pair[1]+1] < self.off_thresh):

                # convert ind_tmp to ind on sta_lta
                latest_pair = latest_pair + (len(sta_lta) - data_len_to_pick)

                # convert indices to UTCDateTime, precision is for checking if new
                start_t = tr.stats.starttime
                on_pick_t = UTCDateTime(start_t + latest_pair[0]/self.sps, precision=0)
                off_pick_t = UTCDateTime(start_t + latest_pair[1]/self.sps, precision=0)

                # check if new
                prev_channel, prev_on_pick, prev_off_pick = prev_pick
                curr_channel = tr.stats.channel # ensure all channel can be picked
                if (curr_channel != prev_channel) \
                   and (on_pick_t != prev_on_pick) \
                   and (off_pick_t != prev_off_pick):

                    return (curr_channel, on_pick_t, off_pick_t)
                else:
                    return None

            else:
                return None


    def run(self):
        traces_queue = self.message_board.subscribe(self.src_topic) # subscribe to source of data
        prev_channel = ""
        prev_on_pick = UTCDateTime(0,precision=0)
        prev_off_pick = UTCDateTime(0,precision=0)

        # start
        while True:
            updated_tr = None

            # read all available data,
            # NOTE: this assumes it doesn't become longer than self.duration...
            for message in traces_queue.listen():
                #print(self.thread_name, " RECEIVED: ", message['data'], " qsize: ", traces_queue.qsize(),sep='')
                additional_tr = message['data']
                updated_tr = self.rt_stream.update_trace(additional_tr)

                if traces_queue.qsize() == 0:
                    break

            # get picks from newly updated trace
            if updated_tr is not None:
                pick = self._get_STALTA_pick_pairs(updated_tr, \
                            (prev_channel, prev_on_pick, prev_off_pick))
                if pick is not None:
                    #publish
                    pick_object = {
                        "channel": pick[0],
                        "on_time": pick[1], \
                        "off_time": pick[2],
                    }
                    self.message_board.publish(self.dst_topic, pick_object)

                    # record latest pick
                    prev_channel = pick_object['channel']
                    prev_on_pick = pick_object['on_time']
                    prev_off_pick = pick_object['off_time']
