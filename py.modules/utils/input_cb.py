# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 18:46:32 2022

@author: DICH
"""
import msvcrt

def input_cb(cb=lambda : None ):
    s='';
    c=0
    while True:
        while not msvcrt.kbhit():
            cb();
        c=msvcrt.getwch()
        msvcrt.putwch(c);
        if c=='\r':            
            break;
        s+=c;
        if c=="\x1a":            
            break;
        
    
    msvcrt.putwch('\n');    
    return s;
    