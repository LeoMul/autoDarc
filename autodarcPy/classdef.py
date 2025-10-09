import numpy as np 
import scipy as sci
from atomic_calc import overlap
from readConfig import readConfigFile

NUM_DARC_POINTS = 250 
FINESTRUC = 0.0072973525643

class autodarcCalc:

    def getGenOcc(self):
        
        nr_occs = readConfigFile('CONFIG.DAT')
        
        nr_av_occs = np.mean(nr_occs,axis=0)
        nvalence = len(nr_av_occs)
        
        
        #set to max occ at first 
        for ii in range(0,self.num_rel_orbs):
            self.genocc[ii] = 2.0 * np.abs(self.kappa[ii])
        
        counter = self.num_rel_orbs-1
        
        for jj in range(0,nvalence):
            #this loop technically goes backwards..
            nr_index = -1 * (jj+1)
            
            if self.orbL[nr_index] == 0:
                self.genocc[counter] = nr_av_occs[nr_index]
                counter -= 1 
            else:
                k1 = abs( self.kappa[counter]  )
                k2 = abs( self.kappa[counter-1])
                ktot = k1+k2
                self.genocc[counter  ] = k1 * nr_av_occs[nr_index]/ktot 
                self.genocc[counter-1] = k2 * nr_av_occs[nr_index]/ktot 

                
                
                counter -= 2 
            
            
        
        return 0 
    
    def setDarcGrid(self):
        
        self.numdarcpoints = NUM_DARC_POINTS 
        self.rnt = self.as_radial[1] 
        maxrho = np.log( self.as_radial[-1] /self.rnt ) 
            
        self.rhogrid = np.linspace(0,maxrho,NUM_DARC_POINTS)
        
        self.darcgrid = self.rnt * np.exp(self.rhogrid)
        

        return 0 
    
    def getRelativisticOrbitals(self):
        
        for ii in range(0,self.num_rel_orbs):
            self.calculateRelativisticOrbital(ii)
            
    
    def calculateRelativisticOrbital(self,ii):
        #for ii in range(0,self.num_rel_orbs):
        
        #I access the same array a lot of times here, does this matter?
        #performance wise?
        
        mapper = self.orbMap[ii]
        interp_scipy = sci.interpolate.CubicSpline(self.as_radial,self.as_orbitals[mapper,:])
        self.large_component[ii,:] = interp_scipy(self.darcgrid)
        derivative = interp_scipy.derivative()(self.darcgrid)
        self.small_component[ii,:] = FINESTRUC * ( derivative + self.kappa[ii] * self.large_component[ii,:] / self.darcgrid )
        
        correction_factor = 1.0 + 0.25 * FINESTRUC * FINESTRUC * ( self.nelec / self.darcgrid)
        self.small_component[ii,:]/= correction_factor
        
        ##for jj in range(0,self.numdarcpoints):
        ##    print(self.darcgrid[jj],self.large_component[ii,jj],self.small_component[ii,jj])
        ##
        ##    #print(self.princN_rel[ii],self.orbL_rel[ii],self.kappa[ii])
        
        norm_factor = overlap(
            self.large_component[ii,:],
            self.large_component[ii,:],
            self.small_component[ii,:],
            self.small_component[ii,:],
            self.darcgrid 
        )
        norm_factor = np.sqrt(norm_factor)
        
        self.large_component[ii,:] /= norm_factor 
        self.small_component[ii,:] /= norm_factor
        
        
        
        return 0 
    
    def writeTargetInp(self):
        f = open('TARGET.INP','w')
                
        #need to include zero here. 
        
        f.write('{:6} {:6}\n'.format(self.num_rel_orbs,self.numdarcpoints+1))
        
        for ii in range(0,self.num_rel_orbs):
            f.write(' {:11.5e}\n'.format(self.genocc[ii]))
        
        f.write(' {:11.5e}\n'.format(0.0))     
        for ii in range(0,self.numdarcpoints):
            f.write(' {:11.5e}\n'.format(self.darcgrid[ii]))     
        
        for ii in range(0,self.num_rel_orbs):
            f.write(' {:6} {:6}\n'.format(self.princN_rel[ii],self.kappa[ii]))
            f.write(' {:11.5e}\n'.format(0.0)) 
            for jj in range(0,self.numdarcpoints):
                f.write(' {:11.5e}\n'.format(self.large_component[ii,jj]))
            
            f.write(' {:6} {:6}\n'.format(self.princN_rel[ii],self.kappa[ii]))
            f.write(' {:11.5e}\n'.format(0.0)) 
            for jj in range(0,self.numdarcpoints):
                f.write(' {:11.5e}\n'.format(self.small_component[ii,jj]))
            
        
        return 0 

    
    def __init__(self,
                 nzed,
                 nelec,
                 chargearray,
                 as_radial,
                 potential,
                 as_orbitals,
                 orbitalParams,
                 princN,
                 orbL,
                 orbEnergyRyd,
                 orbitalFloat
                 ):
        self.nzed = nzed 
        self.nelec = nelec 
        self.chargearray = chargearray
        self.as_radial = as_radial
        self.potential = potential
        self.as_orbitals= as_orbitals
        self.orbitalParams= orbitalParams
        self.princN= princN
        self.orbL= orbL
        self.orbEnergyRyd=orbEnergyRyd
        self.orbitalFloat= orbitalFloat
        
        self.num_nr_orbs = len(orbL)
        num_rel_orbs = 0
        for ii in range(0,self.num_nr_orbs):
            num_rel_orbs += 1
            if orbL[ii] != 0:
                num_rel_orbs += 1
        
        
        
       # print(num_rel_orbs)
        self.num_rel_orbs = num_rel_orbs
        self.setDarcGrid()
        
        self.princN_rel = np.zeros(num_rel_orbs,dtype=int)
        self.kappa      = np.zeros(num_rel_orbs,dtype=int)
        self.orbL_rel   = np.zeros(num_rel_orbs,dtype=int)
        self.orbMap     = np.zeros(num_rel_orbs,dtype=int)
        self.genocc     = np.zeros(num_rel_orbs,dtype=float)
        
        counter = 0
        
        #ugly code
        for ii in range(0,self.num_nr_orbs):
            self.princN_rel[counter] = self.princN[ii]
            self.orbL_rel[counter]   = self.orbL[ii]
            self.kappa[counter]      =  -1 * self.orbL[ii] - 1
            self.orbMap[counter] = ii 
            counter += 1
            
            if orbL[ii] != 0:
                self.princN_rel[counter] = self.princN[ii] 
                self.orbL_rel[counter]   = self.orbL[ii]
                self.kappa[counter]      = self.orbL[ii]
                self.orbMap[counter] = ii 
                counter+=1 
            
        self.large_component = np.zeros([num_rel_orbs,self.numdarcpoints])
        self.small_component = np.zeros([num_rel_orbs,self.numdarcpoints])
        self.getGenOcc()
        for ii in range(0,self.num_rel_orbs):
            print('{:3} {:3} {:10.7f}'.format(
                self.princN_rel[ii],self.kappa[ii],self.genocc[ii]))