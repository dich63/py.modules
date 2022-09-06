
from . import jsonclass,ndarray_marshal,error_marshal,uopen

import jsonrpc.marshals

jsonclass.compress_mode(1);

encode=jsonclass.encode_to_file
#decode=jsonclass.decode_from_file

def decode(fn,jslike=False):
    return jsonclass.decode(uopen.readtxt(fn),jslike);