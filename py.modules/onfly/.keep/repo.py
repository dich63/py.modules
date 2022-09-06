# -*- coding: utf-8 -*-

from httpimport import add_remote_repo, remove_remote_repo;
import httpimport;
httpimport.INSECURE=1

def pkg(*lp):
    url='http://127.0.0.1:8000/'
    add_remote_repo(lp, url);    
    

    
    

