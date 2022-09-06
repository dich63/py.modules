# -*- coding: utf-8 -*-

from .uimport import httpimport
from .uimport.httpimport import add_remote_repo, remove_remote_repo

'''
from onfly.httpimport2 import add_remote_repo, remove_remote_repo;
import onfly.httpimport2 as httpimport;
'''
#import httpimport;
import logging;

httpimport.INSECURE=1
httpimport.RELOAD=1
url='http://127.0.0.1:8000/';

httpimport.logger.setLevel(logging.DEBUG)

def urlset(u):
    global url
    url=u;
    
def pkg(*lp):
    global url
    l=httpimport.logger.level;
    httpimport.logger.setLevel(logging.ERROR);    
    add_remote_repo(lp, url);    
    httpimport.logger.setLevel(l);


    
    

