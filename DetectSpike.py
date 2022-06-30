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

        # prepare station info and thread attribs
        self.station = station
        self.channels = channels
        self.sps = sps
        self.thread_name = self.__class__.__name__ + ":" + self.station
        threading.Thread.__init__(self, name=self.thread_name)

        # prepare source of realtime rolling data
        self.duration = duration
        self.channel_length = duration*sps
        self.rt_stream = RollingStream(self.channels, self.channel_length)

        # prepare STA/LTA params
        self.sta_sec, self.sta_len = sta_sec, sta_sec*self.sps
        self.lta_sec, self.lta_len = lta_sec, lta_sec*self.sps
        self.max_evt_dur, self.max_evt_len = max_evt_dur, max_evt_dur*self.sps
        self.on_thresh = on_thresh
        self.off_thresh = off_thresh

        # prepare pub/sub comms
        self.src_topic = self.station.upper() + "_" + "TRACES"
        self.dst_topic = "PICKS"
        self.message_board = message_board

    def _get_STALTA_pick_pairs(self, tr, prev_pick):
        """ Process a trace and return any new pick-pairs

            Parameters:
            - tr: trace to get picks from
            - prev_pick: UTCDateTime of previously computed pick-pair
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
            if (latest_pair[1] != last_sample_ind) \
                and (data_to_pick[latest_pair[1]+1] < self.off_thresh):

                # convert ind_tmp to ind on sta_lta
                latest_pair = latest_pair + (len(sta_lta) - data_len_to_pick)

                # convert indices to UTCDateTime, precision is for checking if new
                start_t = tr.stats.starttime
                on_pick_t = UTCDateTime(start_t + latest_pair[0]/self.sps, precision=0)
                off_pick_t = UTCDateTime(start_t + latest_pair[1]/self.sps, precision=0)

                # check if new
                prev_on_pick, prev_off_pick = prev_pick
                if (on_pick_t != prev_on_pick) and (off_pick_t != prev_off_pick):
                    return (on_pick_t, off_pick_t)
                else:
                    return None

            else:
                return None


    def run(self):
        # init vars
        queue = self.message_board.subscribe(self.src_topic) # subscribe to source of data
        prev_on_pick = UTCDateTime(0,precision=0)
        prev_off_pick = UTCDateTime(0,precision=0)

        # start
        while True:
            updated_tr = None

            # read all available data,
            # assumes it doesn't become longer than self.duration...
            for message in queue.listen():
                print(self.thread_name, " RECEIVED: ", message['data'], " qsize: ", queue.qsize(),sep='')
                additional_tr = message['data']
                updated_tr = self.rt_stream.update_trace(additional_tr)

                if queue.qsize() == 0:
                    break

            # get picks from newly updated trace
            if updated_tr is not None:
                pick_pair = self._get_STALTA_pick_pairs(updated_tr, (prev_on_pick, prev_off_pick))
                if pick_pair is not None:
                    #publish
                    pick_object = {
                        "channel": updated_tr.stats.channel,\
                        "on_time": pick_pair[0], \
                        "off_time": pick_pair[1],
                    }
                    self.message_board.publish(self.dst_topic, pick_object)

                    # record latest pick
                    prev_on_pick = pick_pair[0]
                    prev_off_pick = pick_pair[1]
