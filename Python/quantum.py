# -*- coding: utf-8 -*-
"""
A very rough class to calculate some parameters that we often use. Wanted to do it myself, 
but struggled with taking the square root of an imaginary matrix. Thus, the qutip package... :(
"""

import numpy
import scipy
import qutip


class State: # defines a state, that moslty contains the density matrix based on some phase value
    
    
    rho =[]
    
    
    
    
    def __init__(self, rho = None, dPhi = None):
        
        if rho is None and dPhi is None: # default state is negative correlated bell state
            
            
            self.rho = [[0j,0j,0j,0j] for x in range(0,4)]           
            self.rho = numpy.array(self.rho)
            self.rho[0,0] = 1+0j
            self.rho[0,3] = -1+0j
            self.rho[3,0] = -1+0j
            self.rho[3,3] = 1+0j
            
        elif dPhi is not None:    
            self.rho = [[0j,0j,0j,0j] for x in range(0,4)]
            self.rho = numpy.array(self.rho)
            self.rho[0,0] = 1+0j#numpy.cos(dPhi/4)**2+0j
            self.rho[0,3] =  1*numpy.exp(-1j*dPhi)+0j#numpy.cos(dPhi/4)*numpy.sin(dPhi/4)+0j
            self.rho[3,0] = 1*numpy.exp(-1j*dPhi)+0j#numpy.cos(dPhi/4)*numpy.sin(dPhi/4)+0j
            self.rho[3,3] =1+0j# numpy.sin(dPhi/4)**2+0j
             
             
        elif rho is not None:
            self.rho = rho 
            
        self.normalizeDensityMatrix()
            
            
    def printDensityMatrix(self): 
        print (self.rho)
        
    def normalizeDensityMatrix(self):
        self.rho = self.rho/sum([self.rho[0][0],self.rho[1][1],self.rho[2][2],self.rho[3][3]])
            
    
    def calcFidelity(self, state):
        """
        Got to check again whether this is the fidelity or the square root fidelity... 
        """
        
        obj1 = qutip.Qobj(self.rho)
        obj2 = qutip.Qobj(state.rho)
        
        F  =qutip.metrics.fidelity(obj1, obj2)
        return (F)
    
if __name__ == "__main__":
    """ 
    Example: define a state that is rotate by pi/4 and calculate the fidelity towards the negative bell state
    """
    p = State(dPhi =numpy.pi/4)
    BS = State()
    
    print (p.calcFidelity(BS))
            
            
            
            
            
            
            