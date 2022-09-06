# Enfold Enterprise Server
# Copyright(C), 2004-5, Enfold Systems, LLC - ALL RIGHTS RESERVED

# Enfold Systems, LLC
# 4617 Montrose Blvd., Suite C215
# Houston, Texas 77006 USA
# p. +1 713.942.2377 | f. +1 832.201.8856
# www.enfoldsystems.com
# info@enfoldsystems.com

import sys

if sys.platform in ('win32',):
    # Use win32process from pywin32
    import win32api
    import win32con
    import win32process
    import pywintypes
    import win32job

    def create_job(name=''):
        w=win32job;
        hjob=w.CreateJobObject(None,'')
        extended_info = w.QueryInformationJobObject(hjob, w.JobObjectExtendedLimitInformation)
        f=w.JOB_OBJECT_LIMIT_DIE_ON_UNHANDLED_EXCEPTION | w.JOB_OBJECT_LIMIT_BREAKAWAY_OK | w.JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        extended_info['BasicLimitInformation']['LimitFlags']=f;
        w.SetInformationJobObject(hjob, w.JobObjectExtendedLimitInformation, extended_info)
        return hjob;


    def assign_process_to_job(hjob,pid=0):
        hp=_get_handle_for_pid(pid)
        win32job.AssignProcessToJobObject(hjob,hp)


    def _get_handle_for_pid(pid, ro=True):
        if pid == 0:
            pHandle = win32process.GetCurrentProcess()
        else:
            flags = win32con.PROCESS_QUERY_INFORMATION
            if not ro:
                flags |= win32con.PROCESS_SET_INFORMATION
            try:
                pHandle = win32api.OpenProcess(flags, 0, pid)
            except pywintypes.error as e:
                raise ValueError(e)
        return pHandle

    def set_process_affinity_mask(pid, value):
        pHandle = _get_handle_for_pid(pid, False)
        current = win32process.GetProcessAffinityMask(pHandle)[0]
        try:
            win32process.SetProcessAffinityMask(pHandle, value)
        except win32process.error as e:
            raise ValueError(e)
        return current

    def get_process_affinity_mask(pid):
        pHandle = _get_handle_for_pid(pid)
        try:
            return win32process.GetProcessAffinityMask(pHandle)[0]
        except win32process.error as e:
            raise ValueError(e)

elif sys.platform in ('linux2'):
    from _affinity import set_process_affinity_mask, get_process_affinity_mask

    def create_job(name=''):
        return None
    
    def assign_process_to_job(hjob,pid=0):
        pass

else:
    def set_process_affinity_mask(pid, value):
        raise NotImplementedError

    def get_process_affinity_mask(pid):
        raise NotImplementedError
