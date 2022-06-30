from pubsub.pubsub import PubSub
from SLclient import SLclient
from RollingStream import RollingStream
from DetectSpike import DetectSpike
from Conversion import Conversion
from obspy.core import read


def main():
    # TODO: get configs from configs file
    IP = "192.168.1.21" # rshake IP
    DUR = 30 # duration in sec (for plotting and other processes)
    STA_SEC = 2 # length of data in seconds to get short-term-ave
    LTA_SEC = 15 # length of data in seconds to get long-term-ave
    ON_THRESH = 3.0 # sta/lta threshold to trigger an on-pick
    OFF_THRESH = 1.5 # sta/lta threshold to trigger an off-pick / reset
    MAX_EVT_DUR = 20 # max duration of on-thresh to off-thresh in sta/lta
    BASIS_CHANNEL = "EHZ" # picks on this channel trigger conversion for all traces
    PADDING_DUR = 3 # duration of padding on the left of event onset
    INTENSITY_THRESHOLD = 2 # intensity greater than this triggers an alarm
    INV_DIR = "inventories" # folder to store station inventory files
    CAP_DIR = "captures" # folder to captured event mseed and png data

    # create communication link between threads
    message_board = PubSub(max_queue_in_a_channel=100)

    # create link to rshake
    SL_client = SLclient(message_board, IP)

    # create spike detector
    Picker = DetectSpike(
        message_board,
        SL_client.station, SL_client.channels, SL_client.sps,
        STA_SEC, LTA_SEC, ON_THRESH, OFF_THRESH, MAX_EVT_DUR,
        DUR
    )

    # create traces converter
    Converter = Conversion(
        message_board,
        SL_client.network, SL_client.station, SL_client.channels, SL_client.sps,
        BASIS_CHANNEL, PADDING_DUR, INTENSITY_THRESHOLD,
        DUR,
        INV_DIR, CAP_DIR
    )

    # subscribe to spike detector pick messages
    event_summary_queue = message_board.subscribe(Converter.dst_event_summary_topic)
    event_data_queue = message_board.subscribe(Converter.dst_event_data_topic)

    # start
    SL_client.start()
    Picker.start()
    Converter.start()
    while True:
        for summary in event_summary_queue.listen():
            print("MAIN: ", "RECEIVED: id=", summary['id'], " data=", summary['data'], " qsize=", event_summary_queue.qsize(), sep='')
            if event_summary_queue.qsize() == 0:
                break
        for dirpath in event_data_queue.listen():
            print("MAIN: ", "RECEIVED: id=", dirpath['id'], " data=", dirpath['data'], " qsize=", event_summary_queue.qsize(), sep='')
            directory = dirpath['data']['path']
            for unit in ["acc", "dis", "counts", "metric", "vel"]:
                st = read(directory + "/" + unit + ".mseed")
                st.plot(outfile = directory + "/" + unit + ".png")

            if event_data_queue.qsize() == 0:
                break

    # wait for all threads to exit gracefully
    #SL_client.join()
    #Picker.join()


if __name__ == "__main__":
    main()




