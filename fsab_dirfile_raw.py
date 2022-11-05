import os
import numpy as np


class fsab_dirfile():
    """
    For analysing FSAB data if libgetdata is not available.
    
    """
    def __init__(self, path, reference='I0000'):
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        self.path = path
        
        #parse sweep data:
        self.sweep={}
        with open(os.path.join(path,'sweep'),'r') as file:
            lines=file.readlines()
        #ignore comments
        lines=[i for i in lines if not i.startswith('#')]
        #ignore metadata
        lines=[i for i in lines if not i.startswith('/')]
        #ignore empty
        lines=[i for i in lines if len(i)>1]
        for line in lines:
            line = line.strip().split(' ')
            fieldname,fieldtype,datatype=line[:3]
            data = np.array(line[3:],dtype=np.float)
            _,fiq,num = fieldname.split('_')
            kidnum=int(num)
            if not kidnum in self.sweep.keys():
                self.sweep[kidnum]={}
            if fiq == 'f':
                self.sweep[kidnum]['f'] = data
            if fiq == 'i':
                if not 'z' in self.sweep[kidnum].keys():
                    self.sweep[kidnum]['z'] = data.astype(np.cdouble)
                else:
                    self.sweep[kidnum]['z'] += data.astype(np.cdouble)
            if fiq == 'q':
                if not 'z' in self.sweep[kidnum].keys():
                    self.sweep[kidnum]['z'] = 1j*data.astype(np.cdouble)
                else:
                    self.sweep[kidnum]['z'] += 1j*data.astype(np.cdouble)
                
        self.numkids = len(self.sweep)
        print(self.numkids)
        
        #load tone freqs
        with open(os.path.join(path,'calibration'),'r') as file:
            lines=file.readlines()
            lines=[i for i in lines if i.startswith('_cal_tone_freq')]
            for line in lines:
                line = line.strip().split(' ')
                fieldname,fieldtype,datatype,data=line
                
                kidnum = int(fieldname[-4:])
                tonefreq = float(data)
                self.sweep[kidnum]['tone_freq'] = tonefreq
        
        self.start_time = np.loadtxt(os.path.join(self.path,'time_start.txt')).item()
        self.stop_time = np.loadtxt(os.path.join(self.path,'time_stop.txt')).item()
        
        print('Ready %s'%(self.path))
        
        
        
    def get_iq_data(self,kidnum):
        assert kidnum < self.numkids
        filename_i = os.path.join(self.path,'I%04d'%kidnum)
        filename_q = os.path.join(self.path,'Q%04d'%kidnum)
        
        i = np.fromfile(filename_i,dtype=np.float32)
        q = np.fromfile(filename_q,dtype=np.float32)
        
        return i + 1j*q
    
        
    def get_sync_data(self):
        filename = os.path.join(self.path,'Q1023')
        return np.fromfile(filename,dtype=np.float32)


        
        #with open(os.path.join(path,'format'),'r') as file:
            #lines=file.readlines()
            ##ignore comments
            #lines=[i for i in lines if not lines.startswith('#')]
            ##ignore metadata
            #lines=[i for i in lines if not lines.startswith('/')]
            ##ignore empty
            #lines=[i for i in lines if not i]
        
            
