# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 14:26:55 2022

@author: wwww
"""


def check_url(u):
    u=u.lower();
    if u.find('http')==0:
        return ( u.find('https://')==0) or ( u.find('http://')==0)
    if u.find('ftp')==0:
        return ( u.find('ftps://')==0) or ( u.find('ftp://')==0)
    return False;


def uopen(u,m='rb'):
    
    if check_url(u):
        from urllib.request import urlopen
        return urlopen(u);
    else:
        return open(u,m)    
    
       
def readtxt(u,code='utf8'):
    return uopen(u,'rb').read().decode(code);