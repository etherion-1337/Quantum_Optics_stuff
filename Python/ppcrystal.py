"""
Calculates the source parameters based on some inputs. Based on: 
    
https://pdfs.semanticscholar.org/e945/06ee76720c6e2a9cbd9658f35aecf6939a97.pdf

Everything here is collinear, type-0. I have tried to implement some emission angle dependency in a few different ways, but so far it has not worked.

The phasematching seems to work quite well, I think I am 8-9 nm off in the poling period. Should be 3.425 um, is 3.434 um for our case. This is
good enough.

"""

import numpy
import matplotlib.pyplot as plt

import quantum
import nonlinearCrystal
import helper
import time


class ppcrystal:

    
    
    def __init__(self, l = 10e3, material = "KTP3", poling = 3.4335, T = 40, alpha = 8.7e-6,beta = 11e-9, wlp = 405e-3, thetaS = 0, waistP =50):
        
        self.PPKTP = nonlinearCrystal.NonlinearCrystal(l = l, material = material, poling = poling) # oject from the nonlinear crystal class
        self.T = T # temperature of the crystal
        self.alpha = alpha 
        self.beta = beta
        self.wlp = wlp
        self.thetaS = thetaS
        self.waistP = waistP
        
        self.loadParams()
        self.calcNEff()
        self.calcResults()
        self.calcFWHM()
        
    def loadParams(self, wlsMin = 750):

        
        self.c = 3e8
        
        self.wls = numpy.linspace(0.75e3,2e3*self.wlp,101)*1e-3
        self.wli = 1/(1/self.wlp-1/self.wls)

    def calcNEff(self):         #https://pdfs.semanticscholar.org/e945/06ee76720c6e2a9cbd9658f35aecf6939a97.pdf

        self.nes = self.PPKTP.calcN(self.wls, theta=90)
        self.nos = self.PPKTP.calcN(self.wls, theta= 0)

        self.nei = self.PPKTP.calcN(self.wli, theta=90)
        self.noi = self.PPKTP.calcN(self.wli, theta= 0)

        self.nep = self.PPKTP.calcN(self.wlp, theta= 90)
        self.nop = self.PPKTP.calcN(self.wlp, theta= 0)


        if self.PPKTP.material == "KTP3":
            a1e= [9.9587,9.9228,-8.9603, 4.1010]
            a2e = [-1.1882e-2, 10.459e-2, -9.8136e-2, 3.1481e-2]
            
            a1o=  [6.2897, 6.3061, -6.0629, 2.6486]
            a2o = [-0.14445e-2, 2.2244e-2, -3.5770e-2, 1.3470e-2]
            
        
        self.n1ep = a1e[0]+ a1e[1]/self.wlp**1 + a1e[2]/self.wlp**2 + a1e[3]/self.wlp**3
        self.n2ep = a2e[0]+ a2e[1]/self.wlp**1 + a2e[2]/self.wlp**2 + a2e[3]/self.wlp**3
        
        self.n1op = a1o[0]+ a1o[1]/self.wlp**1 + a1o[2]/self.wlp**2 + a1o[3]/self.wlp**3
        self.n2op = a2o[0]+ a2o[1]/self.wlp**1 + a2o[2]/self.wlp**2 + a2o[3]/self.wlp**3
        
        self.n1es = a1e[0]+ a1e[1]/self.wls**1 + a1e[2]/self.wls**2 + a1e[3]/self.wls**3
        self.n2es = a2e[0]+ a2e[1]/self.wls**1 + a2e[2]/self.wls**2 + a2e[3]/self.wls**3
        
        self.n1os = a1o[0]+ a1o[1]/self.wls**1 + a1o[2]/self.wls**2 + a1o[3]/self.wls**3
        self.n2os = a2o[0]+ a2o[1]/self.wls**1 + a2o[2]/self.wls**2 + a2o[3]/self.wls**3
        
        self.n1ei = a1e[0]+ a1e[1]/self.wli**1 + a1e[2]/self.wli**2 + a1e[3]/self.wli**3
        self.n2ei = a2e[0]+ a2e[1]/self.wli**1 + a2e[2]/self.wls**2 + a2e[3]/self.wli**3
        
        self.n1oi = a1o[0]+ a1o[1]/self.wli**1 + a1o[2]/self.wli**2 + a1o[3]/self.wli**3
        self.n2oi = a2o[0]+ a2o[1]/self.wli**1 + a2o[2]/self.wli**2 + a2o[3]/self.wli**3

        self.dnes=(1e-6*self.n1es*(self.T-25)+1e-6*self.n2es*(self.T-25)**2)
        self.dnos=(1e-6*self.n1os*(self.T-25)+1e-6*self.n2os*(self.T-25)**2)
    
        self.dnei=(1e-6*self.n1ei*(self.T-25)+1e-6*self.n2ei*(self.T-25)**2)
        self.dnoi=( 1e-6*self.n1oi*(self.T-25)+1e-6*self.n2oi*(self.T-25)**2)
    
        self.dnop=(1e-6*self.n1op*(self.T-25)+1e-6*self.n2op*(self.T-25)**2)
        self.dnep=(1e-6*self.n1ep*(self.T-25)+1e-6*self.n2ep*(self.T-25)**2)   
    
        self.nes_T=(self.nes+self.dnes)
        self.nos_T=(self.nos+self.dnos)
        self.nei_T=(self.nei+self.dnei)
        self.noi_T=(self.noi+self.dnoi)
        self.nep_T=(self.nep+self.dnep)
        self.nop_T=(self.nop+self.dnop)

    def calcResults(self):

        self.L = self.PPKTP.l *(1+self.alpha*(self.T-25)+self.beta*(self.T-25))
    
        self.kp = 2*numpy.pi*self.nop_T/self.wlp
        self.ks =  2*numpy.pi*self.nos_T/self.wls
        self.ki =  2*numpy.pi*self.noi_T/self.wli
        
        self.dk = self.kp-self.ki-self.ks-2*numpy.pi/(self.PPKTP.poling*self.L/self.PPKTP.l)
    
        self.minDK = min(abs(self.dk))
        self.phasematchedSignalWL=self.wls[numpy.argmin(abs(self.dk))]
        self.phasematchedIdlerWL= self.wli[numpy.argmin(abs(self.dk))]    
        
        self.spectrumSignal=(numpy.sin(self.dk*self.L/2)/(self.dk*self.L/2))**2
        self.spectrumIdler=(numpy.sin(self.dk*self.L/2)/(self.dk*self.L/2))**2

        self.brightness=numpy.sum(self.spectrumSignal)
    
        self.spectrum = numpy.concatenate((self.spectrumSignal,self.spectrumIdler[::-1]))
        self.wl= numpy.concatenate ((self.wls, self.wli[::-1]))
   

        self.phi = numpy.zeros([ len(self.wls),len(self.thetaS)])
#        
#        # from https://arxiv.org/pdf/1404.1192.pdf
#
#        for i in range(len(self.thetaS)):
#
#            self.thetaI=numpy.arcsin(self.nos_T*self.ks*numpy.sin(self.thetaS[i])/(self.noi_T*self.ki))
#            self.qi = self.ki*numpy.sin(self.thetaI[i])
#            self.qs = self.ks*numpy.sin(self.thetaS[i])
#                
#            self.dkz = numpy.sqrt((self.ki*self.noi_T)**2 - self.qi**2) + numpy.sqrt((self.ks*self.nos_T)**2 - self.qs**2) - numpy.sqrt((self.kp*self.nop_T)**2 - (self.qs+self.qi)**2)
#            Eplus = numpy.exp(-self.waistP**2*(self.qs +self.qi)**2/4)
#                
#            LAMBDA = self.L*numpy.sqrt(self.kp-self.ks)*numpy.sqrt(self.ks)/(self.noi_T)/self.nos_T*Eplus*numpy.sinc((((self.dkz+2*numpy.pi/self.PPKTP.poling)*self.L))/2/numpy.pi)
#            self.phi[:,i] =  LAMBDA**2
            #self.thetaI=numpy.arcsin(self.nos_T*self.ks*numpy.sin(self.thetaS[i])/(self.noi_T*self.ki))
            #self.dkz=self.nop_T*self.kp-self.nos_T*self.ks*numpy.cos(self.thetaS[i])-self.noi_T*self.ki*numpy.cos(self.thetaI)
            #self.dky=-self.nos_T*self.ks*numpy.sin(self.thetaS[i])+self.noi_T*self.ki*numpy.sin(self.thetaI)
            #self.phi[:,i] =numpy.exp(-((0.1)**2)*(self.dky**2)/2)*(numpy.sin(0.5*self.dkz*self.L)/(0.5*self.dkz*self.L))**2
#    
    def calcFWHM(self):
        
        tmpSignal = self.wls[numpy.where(self.spectrumSignal>0.5*max(self.spectrumSignal))]
        tmpIdler = self.wli[numpy.where(self.spectrumIdler>0.5*max(self.spectrumIdler))]
    
        self.FWHMSignal = max(tmpSignal) - min(tmpSignal) # um
        self.FWHMIdler = max(tmpIdler) - min(tmpIdler) #um
        
        
    def plotSpectrum(self, color = 'k'):
        
        #plt.plot(self.wls, self.spectrumSignal, color)
        #plt.plot(self.wli, self.spectrumIdler, color)
        plt.plot(self.wl,self.spectrum)
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Intensity (arb. units)")
        plt.show()
        
    def plotAngle(self, color ='k'):
        #plt.plot(self.wls,1e3*self.phi, color)
        plt.contourf(((numpy.transpose(self.phi))), extent = ( min(self.wls), max(self.wls),min(self.thetaS), max(self.thetaS)))
        plt.colorbar()
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("Emission angle (mrad)")
        plt.show()
        
if __name__ =="__main__":
    
    ppktp = ppcrystal(T =50, thetaS = numpy.linspace(0, 0.01, 11))
    ppktp.plotSpectrum()
    ppktp.plotAngle()
    
    
    


