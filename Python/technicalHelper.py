# -*- coding: utf-8 -*-
import os
import helper
import numpy
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def thermistorTemperature(dacValue = 200, R0 = 10e3, R = []):

    T0 = 25+273.15;
    B = 3988; 
    
    if R == []:
        T=1/(1/T0+(1/B)*numpy.log(0.27*(3300./dacValue-1))) - 273.15;
    else: 
        T=1/(1/T0+(1/B)*numpy.log(R/R0)) - 273.15;

    return T
 
def dacToLaserCurrent(dacValue):
    return dacValue/4096.0*500

def laserCurrentToDac(current):
    return current/500.0*4096



def quantumEfficiency(wl):
    
    wlData = 1e-3*numpy.array([300,400,450,500,550,600,700,800,900,1000])
    etaData= [0.22, 0.53,0.73,0.83,0.87, 0.86, 0.74, 0.53, 0.32,0.1]
    
    pf = numpy.polyfit(wlData,etaData, 5)
    
    if wl < 0.300 or wl > 1.0000:
        return 0
    else:
        return numpy.polyval(pf,wl)
    
def dichroicMirrorFilterFunction(wl):

    eta = 1 - 0.5*numpy.exp(-(wl-0.810)**2/(0.01)**2)
    return eta


def pressureVSAltitude(h):
    
    hData = numpy.array([0,1, 2,3, 4, 5, 6, 7, 8,9,10, 15, 20, 25, 30,40,50,60,70,80])
    pressureData = numpy.array( [10.13,	8.988,	7.95	,7.012	,6.166	,5.405	,4.722	,4.111,	3.565,	3.08,	2.65,	1.211,0.5529,	0.2549	,0.1197,0.0287,0.007978,	0.002196,	0.00052,	0.00011]  )  
  
    def func(x, a,b):   
            return a * numpy.exp(-b*x);
    p0 = [10, 0.1]
    popt,pcov = curve_fit(func, hData,pressureData, p0 = p0)

    p = func(h, popt[0], popt[1])
    
    return p
    
    
def atmosphericTransmission(wl): #http://instrumentation.tamu.edu/aTmCam.html
    
    data = numpy.loadtxt("D:\\Dropbox\\python\\libraries\\data\\airTransmissionsSpectrum.csv", delimiter = ',')    
    res = numpy.interp(1e9*numpy.array(wl),data[:,0], data[:,1] )
    return res
        

def dichroicMirrorPhaseDifference(wl,centralWL = 810, name = "FF801"):
    
    if name == "FF801":
        filename = os.getcwd()+"\\data\\FF801Di02_phase.csv"
        
    data = numpy.loadtxt(filename, delimiter = ',')
    wlData = data[:,0]
    freq = 2*numpy.pi*3e8/wlData*1e9
    freq0 = 2*numpy.pi*3e8/wlData[helper.findNearestValue(wlData,810)[0]]*1e9

    plt.plot(freq-freq0)
    plt.show()
    
    plt.plot(1e-15*data[2000:,1]/2*(freq[2000:]-freq0)**2)
    plt.show()

    phaseData = 1e-15*data[:,1]/2*(freq-freq0)**2#https://www.rp-photonics.com/group_delay_dispersion.html 
    
    phi = []
    for i in range(len(wl)):
        phi.append(phaseData[helper.findNearestValue(wlData, wl[i])[0]])

    return phi


a = dichroicMirrorPhaseDifference(numpy.linspace(750, 810,101))
plt.plot( a[1000:2000])
    
    
    