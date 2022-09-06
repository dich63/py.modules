# -*- coding: utf-8 -*-

from httpimport import add_remote_repo, remove_remote_repo;
import httpimport;
import logging;

httpimport.INSECURE=1
httpimport.RELOAD=1

url='http://vhost37383.cpsite.ru/plugins/';
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


    
    

