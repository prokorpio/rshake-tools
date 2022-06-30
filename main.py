from pubsub.pubsub import PubSub
from SLclient import SLclient
from RollingStream import RollingStream
from DetectSpike import DetectSpike


def main():
    # TODO: get configs
    IP = "192.168.1.21" # rshake IP
    DUR = 30 # duration in sec (for plotting and other processes)
    STA_SEC = 2
    LTA_SEC = 15
    ON_THRESH = 3.0
    OFF_THRESH = 1.5
    MAX_EVT_DUR = 20 # max duration of on-thresh to off-thresh in sta/lta

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

    # subscribe to spike detector pick messages
    message_queue = message_board.subscribe(Picker.dst_topic)

    # start
    SL_client.start()
    Picker.start()
    for message in message_queue.listen():
        print("RECEIVED: id=", message['id'], " data=", message['data'], " qsize=",message_queue.qsize(), sep='')

    # wait for all threads to exit gracefully
    SL_client.join()
    Picker.join()


if __name__ == "__main__":
    main()
