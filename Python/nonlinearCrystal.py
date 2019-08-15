"""
I have been using this for many many different configurations and types, thus all the unrelated stuff. In our case 
just used to calculate the Sellmeijer equations, the walkoff (maybe) and the phase difference.


"""

import numpy




def partFrac(a,b,wl):
        return a*wl**2 /(wl**2 - b**2);


def calcN(a,wl):

    
        if (not isinstance(wl, list)):
            wl = numpy.array([wl])

        nsq = numpy.ones(len(wl));
        N = (len(a)-1)/2
    
        for i in range(0,int(N)):
            nsq = nsq + partFrac(a[2*i+1], a[2*i+2], wl)    
        
        return numpy.sqrt(nsq)
     

class NonlinearCrystal:
    
    
    
    def __init__(self, l = 6e3, material = "BBO", theta = 28.8, oaOrientation = "up", poling = numpy.inf):
        self.l = l
        self.material = material
        self.theta = theta
        self.oAorientation = oaOrientation
        self.poling = poling
        self.calcNCoefficients();
        
        
    def calcNCoefficients(self):
        

        
        if self.material == "BBO":
            self.ao = [2.7359,0.01878,0.01822,0.01354];
            self.ae = [2.3753,0.01224,0.01667,0.01516];

        elif self.material == "YVO4" or self.material == "YVO":
            self.ae = [4.59905,0.110534,0.04813,0.0122676];
            self.ao  =[3.77879,0.069736,0.04724,0.0108133];
        elif self.material == "FusedSilica" or self.material == "Silica":
            self.ae = [1,0.6961663,  0.0684043,0.4079426 , 0.1162414, 0.89747949,9.896161];
            self.ao = [1,0.6961663,  0.0684043,0.4079426 , 0.1162414, 0.89747949,9.896161];

        elif self.material == "Quartz" or self.material == "Q":
            self.ae = [1,0.665721,  0.0600,0.503511 , 0.1060, 0.21479,0.1190 ,0.539173 ,8.792 ,1.8076613 ,19.70];
            self.ao = [1,0.663044,0.0600,0.517852 , 0.1060, 0.175912,0.1190,0.565380,8.844,1.675299,20.742];
        elif self.material == "MgF":
            self.ae  =[1,0.41344023,0.03684262,0.50497499 , 0.09076162,2.4904862,23.771995];
            self.ao  =[1,0.48755108,0.04338408,0.39875031 , 0.09461442,2.3120353,23.793604];
        elif self.material == "KTP":
            #self.ae = [3.29100,0.04140,0.03978,9.35522,31.45571]
            self.ae = [3.45018,0.04341,0.04597,16.98825,39.43799];
            self.ao = [4.59423,0.06206,0.04763,110.80672,86.12171];
        elif self.material == "KTP2":
            #self.ae = [3.29100,0.04140,0.03978,9.35522,31.45571]
            self.ae = [3.0333 ,0.04154,0.04547,0.01408];
            self.ao = [3.3134,0.05694,0.05658,0.1682];
            
        elif self.material == "KTP3":
            #self.ae = [3.29100,0.04140,0.03978,9.35522,31.45571]
            self.ae = [3.0333 ,0.04154,0.04547,0.01408];
            self.ao = [2.12725,1.184315, 5.14852e-2, 0.6603,100.00507,9.68956e-3]
        
        elif self.material == "LCR":
            self.ae = [1.6795,0.0048,0.0027] 
            self.ao = [1.5187, 0.0016, 0.0011]


            
            
        elif self.material == "CaCO3" or self.material == "Calcite":
            self.ae  =[1.35859695, 0.82427830,1.06689543e-2,0.14429128, 120];
            self.ao  =[1.73358749,0.96464345,1.94325203e-2, 1.82831454, 120];
        
        elif self.material == "BK7":
            self.ae = [1,1.03961212,numpy.sqrt(0.00600069867),0.231792344,numpy.sqrt(0.0200179144),1.01046945,numpy.sqrt(103.560653)]
            self.ao = [1,1.03961212,numpy.sqrt(0.00600069867),0.231792344,numpy.sqrt(0.0200179144),1.01046945,numpy.sqrt(103.560653)]
        else:
            print ("#######################################################################################")
            raise ValueError("NONLINEAR MATERIAL NOT IN DATABASE: CANNOT COMPUTE REFRACTIVE INDEX COEFFICIENTS")
            
            
    def calcN(self, wl, theta = 90):
        """
        theta is the angle with the optical axis!!!
        """
        
        
        #theta = theta+90
        
        
        if self.material == "Q" or self.material =="Quartz" or self.material =="MgF" or self.material =="BK7" or self.material == "FusedSilica" or self.material == "Silica":
            ne = calcN(self.ae, wl)
            no = calcN(self.ao, wl)
        elif self.material == "BBO" or self.material =="YVO4" or self.material == "YVO":
            ne = (self.ae[0] + self.ae[1] / (wl**2-self.ae[2]) - self.ae[3] * wl**2)**0.5;
            no = (self.ao[0] + self.ao[1] / (wl**2-self.ao[2]) - self.ao[3] * wl**2)**0.5;
        elif self.material =="KTP": ## kato
            ne = (self.ae[0] + self.ae[1] / (wl**2-self.ae[2]) + self.ae[3] / (wl**2-self.ae[4]))**0.5            
            no = (self.ao[0] + self.ao[1] / (wl**2-self.ao[2]) + self.ao[3] / (wl**2-self.ao[4]))**0.5  
        elif self.material =="KTP2": ##CHAUDHARY1
            ne = (self.ae[0] +  self.ae[1]/ (wl**2 - self.ae[2]) -  self.ae[3]*wl**2)**0.5
            no =  (self.ao[0] +  self.ao[1]/ (wl**2 - self.ao[2]) -  self.ao[3]*wl**2)**0.5

        elif self.material =="KTP3": ##Fradkin https://pdfs.semanticscholar.org/b6e1/bb6430a557912f0704e84f19416cc1f4653f.pdf
            ne = (self.ae[0] +  self.ae[1]/ (wl**2 - self.ae[2]) -  self.ae[3]*wl**2)**0.5
            no =  (self.ao[0] +  self.ao[1]/ (1- self.ao[2]/wl**2) + self.ao[3]/ (1- self.ao[4]/wl**2) -  self.ao[5]*wl**2)**0.5

        elif self.material == "LCR":
            ne = (self.ae[0]+self.ae[1]/wl**2+self.ae[2]/wl**4)
            no = (self.ao[0]+self.ao[1]/wl**2+self.ao[2]/wl**4)

        elif self.material == "CaCO3" or self.material == "Calcite":
            ne = (self.ae[0] + self.ae[1]*wl**2 / (wl**2-self.ae[2]) - self.ae[3] * wl**2/ (wl**2-self.ae[4]))**0.5;
            no = (self.ao[0] + self.ao[1]*wl**2 / (wl**2-self.ao[2]) - self.ao[3] * wl**2/ (wl**2-self.ao[4]))**0.5;
        else:
            print ("#######################################################################################")
            raise ValueError("NONLINEAR MATERIAL NOT IN DATABASE: CANNOT COMPUTE REFRACTIVE INDEX DATA")


        return (1/((numpy.cos((theta)/180.0*numpy.pi)/no)**2+((numpy.sin((theta)/180.0*numpy.pi)/ne)**2)))**0.5


    def calcPhase(self, wl, l =None, pol = "o",alpha= 0, gamma = 0):
        """
        Phase of a photon of wavelength wl in a crystal of length l (if not provided self.l) with the optical axis orientation of theta
        """
        if l == None:
            l = self.l;
            
        if pol == "e":
            if self.oAorientation == "up":
                n = self.calcN(wl, theta = self.theta-alpha)
            if self.oAorientation == "right":
                n = self.calcN(wl, theta = self.theta-gamma)
            
        else:
            n = self.calcN(wl, theta = 0)

#print rho*180/numpy.pi
        
       
        if (pol == "e"):
            thetaEff = numpy.pi/180.*(self.theta-alpha);
            rho = thetaEff-numpy.arctan(self.calcN(numpy.mean(wl))**2/self.calcN(numpy.mean(wl),thetaEff)**2*numpy.tan(thetaEff));
            kS = numpy.cos(rho) # correction factor, see altepter 2005
        else:
            kS = 1;
            rho = alpha/180.*numpy.pi;
            
         
            
#        print rho*180/numpy.pi
#        print kS
#        print 1/numpy.cos(rho)
#        
#        
#        print kS/numpy.cos(rho)
#        print 1/numpy.cos(gamma/180.*numpy.pi)
        #print kS/numpy.cos(rho)
        
        
        return 2*numpy.pi*l*kS/numpy.cos(alpha/180.*numpy.pi-rho)/numpy.cos(gamma/180.*numpy.pi)*(n/wl)    # correction factor, see altepter 2005
        #return 2*numpy.pi*L
    
    
    def calcPropagationAngle(self, wl, alpha = 0, pol = "o"):
        
            if pol == "e":
                
                alpha  = alpha/self.calcN(numpy.mean(wl), self.theta)
                

                


            else:
                alpha  = alpha/self.calcN(wl, 90)

            
            return alpha
        
        
    
    def calcPropagationAngleOutside(self, wl, alpha = 0, pol = "o"):
        
        if pol == "e":
            alpha  = alpha*self.calcN(wl, self.theta)
            
        else:
            alpha  = alpha*self.calcN(wl, 90)
            
        
        return alpha
    

    def calcIdlerEmissionAngle(self, wls, wli = None,  alpha = 0, theta = 28.80):
        
        if wli == None:
            
            wli = 1/(1/0.405- 1/wls)

        
        nos =  self.calcN(wls, 90);
        noi =  self.calcN(wli, 90);

        alpha = (alpha) / 180 * numpy.pi

        x = alpha**2/numpy.arcsin(nos*wls*numpy.sin(alpha)/(noi*wli))

        return (x)/numpy.pi*180 

    def calcWalkoff(self, wl = 0.405, theta = 45):
        
        neff= self.calcN(wl, theta)
        no = self.calcN(wl, 90)
        ne = self.calcN(wl, 0)
        rho = 1/2 * neff**2 * (1/no**2 - 1./ne**2) * numpy.sin(2*theta/180.0*numpy.pi);
        return rho
        
    
if __name__ == "__main__":
    
    cryst = NonlinearCrystal(material = "BBO")
    print(cryst.calcWalkoff(0.405, 28.8)*cryst.l)
    
    
    
    