# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 00:50:05 2021

@author: wwww
"""

import pika
import numpy as np
connection = pika.BlockingConnection(pika.ConnectionParameters('localhost'))
channel = connection.channel()
#channel.queue_declare(queue='hello')
channel.basic_publish(exchange='',
                      routing_key='hello',
                      body='np.array([1,2,3])')
print(" [x] Sent 'Hello World!'")
connection.close()



'''
channel = connection.channel()
method_frame, header_frame, body = channel.basic_get('test')
if method_frame:
    print(method_frame, header_frame, body)
    channel.basic_ack(method_frame.delivery_tag)
else:
    print('No message returned')
    
'''