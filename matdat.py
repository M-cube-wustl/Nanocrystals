import numpy as np
import matplotlib.pyplot as plt


def main():
    plt.close('all')
    # path to Results4Dong
    #path = r'C:\Users\jcavi\Box Sync\Research\Rohan\Projects\Nanocrystals\MATDAT\Results4Dong'
    path = #'your\path\here'
    
    datadir = path+'\\Exchange_Density\\individual_data'
    datapath = datadir+'\\mc_03-23-2020Parabola-k=0.03.matdat'
    # loading data
    simdata = NC(datapath)
    # plotting the 10 example individual trajectories
    plt.figure()
    plt.plot(simdata.time,simdata.indv.T)
    # plotting dG vs i
    plt.figure()
    plt.plot(simdata.swaps,simdata.dG)
    plt.show()
    print simdata.WTs.shape
    
    pass

    
    

    
class MATDAT:
    def __init__(self, fp):
        '''initializer for MATDAT class -- enter filepath to .matdat file'''
        f = open(fp,'r')
        s = f.read()
        s = s.replace('\t\t','').replace('\t\n','\n')
        ss = s.split('\n*\n')
        self.label = ss[0]
        self.data = {}
        for l in ss[1:-1]:
            (header, data) = self.split_data_header(l)
            if self.isNum(data):
                dataarr = self.str2arr(data)
            else:
                dataarr = data.split()
            self.data[header] = dataarr
  
    
    def rn(self,s):
        '''remove linebreaks'''
        return s.replace('\n','') 
        
    def str2arr(self,s):
        '''converts a tab and linebreak delimitted string into a numpy array'''
        ss = s.split('\n')[0:-1]
        out = []
        for i,line in enumerate(ss):
            sline = self.rn(line).split('\t')[0:]
            outi = []
            for j,d in enumerate(sline):
                outi.append(float(d))
            out.append(outi)
        return np.squeeze(np.asarray(out))
        
    def split_data_header(self,s):
        '''splits an entry in the MATDAT into the header and data'''
        ss = s.split('\n')
        header = ss[0]
        out = ''
        for line in ss[1:]:
            out += line+'\n'
        return (header,out)
        
        
    def isNum(self,s):
        '''check if the data is a number'''
        test = s.split('\n')[0].split('\t')[0]
        try:
            float(test)
            return True
        except ValueError:
            return False
            
class NC:
    def __init__(self,fp):
        '''wrapper of MATDAT specifically for these simulation datasets'''
        self.MD = MATDAT(fp)
        md = self.MD
        self.time  = md.data['time'] # simulation time (N_steps x 1 array)
        self.swaps = md.data['swap'] # swap number (20 x 1 array)
        self.dG    = md.data['dG'] # free energy change from ion swapping (20 x 1 array)
        self.p     = md.data['p'] # probability of ion swapping (20 x 1 array)
        self.ens   = md.data['ions_ens'] # ensemble ions swapped (N_steps x 1 array)
        self.indv  = md.data['ions_indv'] # individual ions swapped (10 x N_steps array)
        self.WTs   = md.data['waittimes'] # wait times of each NC (N_NC x 1 array)
        self.STs   = md.data['switchtimes'] # switching times of each NC (N_NC x 1 array)
        

    
if __name__ == '__main__':
    main()    
    
    
    
    
    
    
    
    
1    