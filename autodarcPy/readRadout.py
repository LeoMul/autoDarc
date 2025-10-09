import fortranformat as ff 
import numpy as np 
from classdef import autodarcCalc 

def readRadout(path_to_radout):
    
    fileHeaderFormat = ff.FortranRecordReader('(3I5,4X,I9,I4,F5.1,5F7.3,A4)')
    orbiHeaderFormat = ff.FortranRecordReader('(3I5,I3,F5.1,2X,F12.6,I6,29X,A4)')
    dataRecordFormat = ff.FortranRecordReader('(I5,2(I4,2(1PE14.7)))')
    
    radout = open(path_to_radout,'r')
    
    header  = radout.readline()
    
    da = [0,0,0,0,0]
    
    key,num_orbs,mb,num_points,nelec,nzed,da[0],da[1],da[2],da[3],da[4],tlbl = fileHeaderFormat.read(header) 
    
    charge = np.zeros(num_points,dtype=float)
    radial = np.zeros(num_points,dtype=float)
    
    orbitals  = np.zeros([num_orbs,num_points],dtype=float)
    potential = np.zeros([num_orbs,num_points],dtype=float)
    orbitalParams  = np.zeros(num_orbs,dtype=float)
    princN = np.zeros(num_orbs,dtype=int)
    orbL = np.zeros(num_orbs,dtype=int)
    orbEnergyRyd = np.zeros(num_orbs,dtype=float)
    orbitalFloat = np.zeros(num_orbs,dtype=float)

    for i in range(0,num_points,2):
        line = radout.readline()
        key,dummy,radial[i],charge[i],dummy,radial[i+1],charge[i+1] = dataRecordFormat.read(line)
    #read(1,2) key,dummy,radial(i),charge(i),dummy,radial(i+1),charge(i+1)
    
    #print(radial[-1])
    
    for i in range(0,num_orbs):
        line = radout.readline()
        key,dummy,princN[i],orbL[i],orbitalFloat[i],orbEnergyRyd[i],dummy,tlbl = orbiHeaderFormat.read(line)
        for j in range(0,num_points,2):
            line = radout.readline()
            key,dummy,orbitals[i,j],potential[i,j],dummy,orbitals[i,j+1],potential[i,j+1] = dataRecordFormat.read(line)
        orbitalParams[i] = orbitals[i,0]
        orbitals[i,0]  = 0.0e0
    
    
    
    return autodarcCalc (nzed,nelec,charge,radial,potential,orbitals,orbitalParams,princN,orbL,orbEnergyRyd,orbitalFloat)