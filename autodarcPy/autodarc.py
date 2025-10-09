from readRadout import * 
adc = readRadout('radout')
adc.getRelativisticOrbitals() 
adc.writeTargetInp()