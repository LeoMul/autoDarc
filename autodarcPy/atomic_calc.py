import numpy as np  

#todo - put this in a class 
#      put orbitals into a class?

def overlap(large1,large2,small1,small2,radialgrid):
    
    #integrates an integral by trapezium rule - could be better
    
    array = large1 * large2 + small1 * small2 
    integral = np.trapz(array,radialgrid)

    #in this implemetation I don't have zero as a point, account for this with:
    integral += array[0] * 0.5 * radialgrid[0]
    #i.e adding the contribution from zero at the orgin.
    
    return integral 

