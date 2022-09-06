# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 11:27:47 2021

@author: wwww
"""

from threading import Thread

class asyn_interact(Thread):
    def __init__(this,**kwd):
        Thread.__init__(this);
        this.kwd=kwd;
        this.g=globals()
        this.g=kwd
    def run(this):
        #this.g.update(this.kwd)
        import code;
        code.interact(local=this.kwd)
        #print(this.kwd);