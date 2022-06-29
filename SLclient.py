import threading
from obspy.clients.seedlink import basic_client, easyseedlink

class SLclient(threading.Thread):
    """ Class to process and forward mseed from given station """

    def __init__(self, message_board, station_IP):
        """
            Parameters:
            - message_board: system PubSub board to publish mseed packets to

        """

        self.station_IP = station_IP
        self.network, self.station = self._get_info(station_IP)
        self.thread_name = self.__class__.__name__ + ":" + self.station
        threading.Thread.__init__(self, name=self.thread_name)

        self.topic = self.station.upper() + " " + "TRACES"
        self.message_board = message_board

        # setup realtime stream
        self.rt_sl_client = easyseedlink.create_client(station_IP)
        self.rt_sl_client.select_stream(self.network, self.station, "E??") # All channels
        self.rt_sl_client.on_data = self._data_pckt_rcvd_callback


    def _get_info(self,station_IP):
        client = basic_client.Client(station_IP)
        network, station = client.get_info(level="station")[0]
        return network, station

    def _data_pckt_rcvd_callback(self, tr):
        """ Add data to STATION_TRACES topic

            Parameters:
            - tr: trace Obspy object

            Note:
            - This fn is being called ~ once every 0.5 sec.
              Make sure each subscriber reads data from their
              queue at a faster rate to avoid queue overflow.

        """
        self.message_board.publish(self.topic, tr) # publish rate: ~ once every 0.5sec
                                                   # make sure each subscriber reads data
                                                   # from their queue at a faster rate

    def run(self):
        self.rt_sl_client.run()
