from pubsub.pubsub import PubSub
from SLclient import SLclient
from RollingStream import RollingStream

import threading

#TODO: get configs
# rshake IP
# buffer length : 30sec

def main():
    # get configs
    IP = "192.168.1.21"

    # create communication link between threads
    message_board = PubSub(max_queue_in_a_channel=100)

    # create link to rshake
    SL_client = SLclient(message_board, IP)

    # create receiving stream
    rt_stream = RollingStream(["EHZ", "ENN", "ENE", "ENZ"], channel_length=40*100)

    # subscribe
    message_queue = message_board.subscribe(SL_client.topic)

    # start
    SL_client.start()
    for message in message_queue.listen():
        tr = message['data']
        print("RECEIVED: id=", message['id'], " data=", message['data'], " qsize=",message_queue.qsize(), sep='')
        rt_stream.update_trace(tr)
        if len(rt_stream.traces[0]) > 3000:
            latest_traces = rt_stream.latest_traces(duration=30)
            for tr in latest_traces:
                print(tr)
            message_board.unsubscribe("RE722_TRACES", message_queue)
            break

    # wait for all threads to exit gracefully
    SL_client.join()


if __name__ == "__main__":
    main()
