#Writes a DSTG2.INP from the CONFIG.DAT 
#User input variables are nast and num pws.s
#This code is hilariously unpretty 
#It is meerely a convenience for me to get calculations going.


import numpy as np 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-l','--num_levels',  help='Number of close-coupling levels',type= int)
parser.add_argument('-p', '--num_partial_waves',  help='Number of N+1 partial waves',type=int)
args = parser.parse_args()

if not args.num_partial_waves:
    print("no partial waves requested")
    parser.print_help()
    exit()

if not args.num_levels:
    print("no cc levels requested")
    parser.print_help()
    exit()

angular_symbols = np.array(['s','p','d','f','g','h'])
angular_kappa   = np.array([ -1, -2, -3, -4, -5, -6])

f = open('CONFIG.DAT','r')

firstline = f.readline().split()

num_orbs = int(firstline[1])

second_line = f.readline().split()

princN = np.zeros(num_orbs,dtype=int)
orbL = np.zeros(num_orbs,dtype=int)
orbs = []

for ii in range(0,num_orbs,1):
    princN[ii] = int(second_line[2*ii])
    orbL[ii] =   int(second_line[2*ii+1])
    #print(orbL[ii])
    string = str(princN[ii])+angular_symbols[orbL[ii]]
    orbs.append(string)
    
num_csfs = int(f.readline().split()[0])

f.close()

#easiest way for me to do this:
csfs = np.loadtxt('CONFIG.DAT',skiprows=5,dtype=int)

def write_many_orbitals(core_orbitals,peel_orbitals,occupations,f):
    #first do core
    #f = open('DSTG2.INP','w')
    for ii in range(0,len(core_orbitals)):
        princ_n = int(core_orbitals[ii][0:len(core_orbitals[ii])-1])
        orbital_symbol = core_orbitals[ii][-1]
        index_of_angular_symbol = np.where(angular_symbols == orbital_symbol)[0]
        kappa = angular_kappa[index_of_angular_symbol][0]
        orbital_string = make_orbital_string(princ_n,kappa,[])
        f.write(orbital_string)
    for ii in range(0,len(peel_orbitals)):
        princ_n = int(peel_orbitals[ii][0:len(peel_orbitals[ii])-1])
        orbital_symbol = peel_orbitals[ii][-1]
        index_of_angular_symbol = np.where(angular_symbols == orbital_symbol)[0]
        kappa = angular_kappa[index_of_angular_symbol][0]
        orbital_string = make_orbital_string(princ_n,kappa,occupations[:,ii])
        f.write(orbital_string)

    return 0

def make_orbital_string(princ_n,kappa,occupations):

    orbital = '&ORB PRINC={} KAPPA={} '

    string = orbital.format(princ_n,kappa)

    if len(occupations) > 1:
        string += 'CSF= '
        counter = 0
        for ii in range(0,len(occupations)):
            counter += 1
            string += str(int(occupations[ii])) +' '
            if (counter > 5) and (ii < len(occupations)-1):
                 string +=', \n CSF={}* '.format(ii+1)
                 counter = 0

    string += '/ \n'
   # print(string)

    return string



def write_partial_wave():
    return 0

def write_dstg2(num_levels,num_partial_waves,core_orbitals,peel_orbitals,occupations):
    total_orbtials = len(core_orbitals) + len(peel_orbitals)

    num_csfs = np.shape(occupations)[0]

    f = open('DSTG2.INP','w')

    f.write('DSTG2: DARC\n')
    f.write('&PREINP NPW={} NP_PER_SUBWORLD=1 IDIMCHECK=0 IANGULAR=-1 /\n'.format(num_partial_waves))
    f.write('&DSTG2 NWM={} NMAN={} INAST={} NAST={}/\n'.format(total_orbtials,num_csfs,num_partial_waves,num_levels))

    write_many_orbitals(core_orbitals,peel_orbitals,occupations,f)

    f.write(' &ANGOPT /\n')
    f.write(' &JVALUE /\n')

    num_electrons = int(np.sum(occupations[0,:]))
    #print(num_electrons)
    
    if (num_electrons%2) == 0:
        jstart = 0.5
    else:
        jstart = 0.0
    j = jstart
    print(num_partial_waves)
    for ii in range(0,int(num_partial_waves/2)):
        f.write( '&SYM JTOT={}  NPTY= 1 /\n'.format(j))
        f.write( '&SYM JTOT={}  NPTY=-1 /\n'.format(j))
        j += 1
        
        
write_dstg2(args.num_levels,args.num_partial_waves,[],orbs,csfs)
