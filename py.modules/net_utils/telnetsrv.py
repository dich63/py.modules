#!/usr/bin/env python3

import socket

# Connect to the server with `telnet $HOSTNAME 5000`.

server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.setblocking(False)
server.bind((socket.gethostname(), 5000))
server.listen(5)

connections = []

while True:
    try:
        connection, address = server.accept()
        connection.setblocking(False)
        connections.append(connection)
    except BlockingIOError:
        pass

    for connection in connections:
        try:
            message = connection.recv(4096)
        except BlockingIOError:
            continue

        for connection in connections:
            connection.send(message)