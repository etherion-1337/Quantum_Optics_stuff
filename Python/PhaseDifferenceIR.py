import sys
import pylab
import numpy
import matplotlib.pyplot as plt
import matplotlib
import quantum
import phasematching
import nonlinearCrystal
import helper
import time



# prepare constants, initial guesses, wavelengths and angle ranges!
LYVO_pre = 0#0.79e3
LBBO1 = 9.0e3#19e3##14.66e3#8.5e32*9e3
LBBO2 =  9.1e3#19.15e3#15.49e3#9e32*9.1e3
LYVO =2*1e3;
LQ= 0e3#10e3 #@1.15e3#5e3;



c = 3e8;
nu_p = c/658e-9;
nu_i = nu_p/2-numpy.linspace(0.1,35,51)*1e12
nu_s = nu_p/2+numpy.linspace(0.1,35,51)*1e12


opticalAxis =45
pump = 0.658;
s = c/nu_s*1e6; #convert to microns!
i = c/nu_i*1e6;
wlArray = numpy.concatenate((numpy.array(list(reversed(s))),i))

res = 21
maxAngle = 1
alpha = numpy.linspace(-maxAngle,maxAngle,res)
gamma  = numpy.linspace(-maxAngle,maxAngle,res)
#halfAngle =alpha[len(alpha)/2:-1]

#################################################### PREPARE CRYSTALS #########
YVO_pre = nonlinearCrystal.NonlinearCrystal(l = LYVO_pre, material = "Silica", theta =90)

HWP_p = nonlinearCrystal.NonlinearCrystal(l = 1.08e3, material = "Q", theta = 90)
BBOI = nonlinearCrystal.NonlinearCrystal(l = LBBO1, material = "BBO", theta =opticalAxis)
BBOII = nonlinearCrystal.NonlinearCrystal(l = LBBO2, material = "BBO", theta =opticalAxis)

HWP1 =  nonlinearCrystal.NonlinearCrystal(l = 1.125e3, material = "MgF", theta =90)
HWP2 =  nonlinearCrystal.NonlinearCrystal(l = 0.892e3, material = "Q", theta = 90)

YVO = nonlinearCrystal.NonlinearCrystal(l = LYVO, material = "YVO", theta =90, oaOrientation="right")


################################################################################


"""
calculate the different effects of effective length and refractive indices
"""


#for opticalAxis in [28.77, 28.78, 28.79, 28.80, 28.81]:

 #   for LBBO in [2.0e3,3e3,4e3,5e3,6e3,8e3,10e3,12e3]:
######## PHASE FUNCTIONS
if 1:
    if 1:
        
            def calcSystemPhaseDifference(comp = 1, Lcomp = LYVO, alpha = 0, gamma = 0,yvo4AngleFactor=0, quartzAngleFactor = 0, BBOCompAngle = 4.7, pumpWL = 0.405):

                """ assumption: everything is just ordinarily polarized in the HWP. This adds a slight error (< 1 percent), but saves me the headache"""



                dPhiBBO1_1 = BBOI.calcPhase(pumpWL,pol = "e") #left
                dPhiBBO1_2 = BBOI.calcPhase(pumpWL,pol = "o") # right (the one that we align on)

                dPhiHWPpump_1 = HWP_p.calcPhase(pumpWL, pol = "e")/2+HWP_p.calcPhase(pumpWL, pol = "O")/2

                dPhiYVO_pre_2 = YVO_pre.calcPhase(pumpWL, pol = "o")

                dPhiHWP1_1 = (HWP1.calcPhase(s,pol = "o")+HWP1.calcPhase(i,pol = "o")+HWP1.calcPhase(s,pol = "e")+HWP1.calcPhase(i,pol = "e"))/2.
                dPhiHWP2_1 = (HWP2.calcPhase(s,pol = "o")+HWP2.calcPhase(i,pol = "o")+HWP2.calcPhase(s,pol = "e")+HWP2.calcPhase(i,pol = "e"))/2.

                dPhiHWP1_2 = (HWP1.calcPhase(s,pol = "o")+HWP1.calcPhase(i,pol = "o")+HWP1.calcPhase(s,pol = "o")+HWP1.calcPhase(i,pol = "o"))/2.
                dPhiHWP2_2 = (HWP2.calcPhase(s,pol = "o")+HWP2.calcPhase(i,pol = "o")+HWP2.calcPhase(s,pol = "o")+HWP2.calcPhase(i,pol = "o"))/2.


                dPhiBBOII_1 = BBOII.calcPhase(s, pol = "e") + BBOII.calcPhase(i, pol = "e")
                dPhiBBOII_2 = BBOII.calcPhase(s, pol = "o") + BBOII.calcPhase(i, pol = "o")

                dPhiYVO_1 = YVO.calcPhase(s,l = Lcomp, pol = "e") + YVO.calcPhase(i, l = Lcomp,pol = "e")
                dPhiYVO_2 = YVO.calcPhase(s, l = Lcomp,pol = "o") + YVO.calcPhase(i, l = Lcomp,pol = "o")

                
                                
                """ sum of all phases """

                dPhi = dPhiBBO1_1 - dPhiBBO1_2 - dPhiYVO_pre_2+dPhiHWPpump_1  + dPhiBBOII_1- dPhiBBOII_2 + dPhiYVO_1 - dPhiYVO_2 
            

                return dPhi
                
            def gradientDescent(Lstart = LYVO, alpha = 0):
                
                loopCounter = 0
                precision = 1e-7
                L =1.0*Lstart;
                gamma = 1
                grad = 0
                tiny = 1e-5*L
                grad = 0
                gradOld = 1
                
                
                
                while precision < abs(grad-gradOld):
                    loopCounter+=1
                    
                    gradOld = grad
                    dPhi1 = abs(calcSystemPhaseDifference(Lcomp=L, alpha = alpha))
                    dPhi2 = abs(calcSystemPhaseDifference(Lcomp=L+tiny, alpha = alpha))

                    
                    
                    S2 = sum((dPhi2-min(dPhi2))**2)
                    S1 = sum((dPhi1-min(dPhi1))**2)
                    grad = (S2-S1)/(tiny)
                    
                    
                    
                    
                    L -= grad*gamma
                    gamma = (1000*numpy.sqrt(abs(grad)))        
                    
                    
                   #print(L)
                    if loopCounter >200:
                        break
                return L
            
#############################
                
            print ("###########################################")
            d1 = (BBOI.calcWalkoff(0.658, opticalAxis)*BBOI.l)
            d2 = (BBOII.calcWalkoff(1.316, opticalAxis)*BBOII.l)
            
            print ("Walkoff splitter:   " + str(d1))
            print ("Walkoff combiner:   " + str(d2))
            
            print ("Walkoff difference: " + str(d1-d2))
            print ("Relative walkoff difference:  " + str((d1-d2)/d1) )
            print ("###########################################")
#############################            
            
            L = gradientDescent()
            
            print (L)
            
            pump = numpy.linspace(0.658,0.659,11)
            dPhi = []
            for j in range(len(pump)):
                nu_p = c/(pump[j]*1e-6);
                nu_i = nu_p/2-numpy.linspace(0.1,35,51)*1e12
                nu_s = nu_p/2+numpy.linspace(0.1,35,51)*1e12
                
                
                opticalAxis =28.81

                s = c/nu_s*1e6; #convert to microns!
                i = c/nu_i*1e6;
                
                dPhi.append(calcSystemPhaseDifference(Lcomp = L, alpha = 0, pumpWL = pump[j]));
                helper.setFonts(16)
                
            plt.contourf(s,pump, dPhi- min(dPhi[j]))    
            plt.colorbar()
            #plt.clabel("Phase difference (rad)")
            plt.ylabel("Pump wavelength (nm)")
            plt.xlabel("Signal wavelength (nm)")

                
            plt.figure()    
            plt.plot(s,dPhi[0]-(dPhi[0][0]), 'k')
            plt.plot(i,dPhi[0]-(dPhi[0][0]), 'k')
            plt.xlabel("Wavelength (nm)")
            plt.ylabel("Phase difference (rad)")


            visibility = []
            
            
            #numpy.savetxt("phase_pomo.txt", numpy.c_[s, dPhi])

            
