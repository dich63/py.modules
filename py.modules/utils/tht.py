# -*- coding: utf-8 -*-
import code

import asyncio
from threading import Timer

async def hello_world():
    code.interact(local=globals());
    

def wait():
    code.interact(local=globals());

timer = Timer(1, wait)
timer.start()