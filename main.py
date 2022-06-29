from pubsub.pubsub import PubSub
from SLclient import SLclient

#TODO: get configs

def main():
    # get configs
    IP = "192.168.1.21"

    # create communication link between threads
    message_board = PubSub(max_queue_in_a_channel=100)

    # create link to rshake
    SL_client = SLclient(message_board, IP)
    SL_client.start()

    # create subscribers
    messages = message_board.subscribe(SL_client.topic)
    for message in messages.listen():
        print("RECEIVED: id=", message['id'], " data=", message['data'], " qsize=",messages.qsize(), sep='')

    # wait for all threads to exit gracefully
    SL_client.join()


if __name__ == "__main__":
    main()
