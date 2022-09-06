#
import sys, importlib,requests


def load_module(name,code,g=None):
    spec = importlib.util.spec_from_loader(name, loader=None)
    module = importlib.util.module_from_spec(spec)
    exec(code, module.__dict__)
    sys.modules[name] = module        
    if  type(g) is dict:
        g[name] = module
    return module

def url2module(name,url,g=None):
    r =requests.get(url);
    if r.status_code==200:
        m=load_module(name,r.text,g);
        #m.__file__=url;
        m.__file__='Хер вам нету!!!';
        return m
    else:
        return None;

