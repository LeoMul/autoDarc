import numpy as np 

def readConfigFile(path_to_config_dot_dat):
    '''
    This routine will read the CONFIG.DAT file,
    which is perhaps a more reliable way of obtaining the data than in das,
    which may or may not be formatted in a particular way.
    
    The code is (probably) similar to that of autoDarc 
    https://github.com/LeoMul/autoDarc/blob/main/write_dstg2inp.py
    
    '''
    
    f = open(path_to_config_dot_dat,'r')

    firstline = f.readline().split()

    if len(firstline) == 2:
        korb = int(firstline[0])

    second_line = f.readline().split()
    
    num_orbs = int(len(second_line) / 2 )

    third_line = f.readline().split()
    
    num_csfs = int(third_line[0])
    
    princN = []
    orbital_L = [] #np.zeros(num_orbs,dtype=int)
    orbital_strings = [] 
    max_occupations = []
    for ii in range(0,len(second_line),2):
        princN.append(int(second_line[ii]))
        ll = int(second_line[ii+1])
        orbital_L.append(ll)
        max_occupations.append(2*(2*ll+1))


    f.close()


    #easiest way for me to do this:
    das_file_numpy = np.loadtxt(path_to_config_dot_dat,skiprows=5,dtype=int)[:,0:-1]
    #print(das_file_numpy)
    #temperory
    
    #print(orbital_strings)
    
    return das_file_numpy