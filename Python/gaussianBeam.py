import numpy
import matplotlib.pyplot as plt



class Lens:
    
    z = []
    f  =[]
    
    
    def __init__(self, z = 0, f = 10):
        self.f = f;
        self.z = z;

class GaussianBeam:
    
    elementList = [] # the list of lenses and other segments that are impacting the beam path
    segmentList = [] # different Gaussian beam segments


    def __init__(self, segmentList = [], elementList = []):
        
        self.elementList = elementList
        self.segmentList = segmentList
        
        
    def lensTransferFunction(self, segment, f):
        
        newSegment =GaussianBeamSegment(wl = segment.wl, z = segment.z[-1] )
        Q = segment.q[-1]
        newSegment.q = numpy.array([Q/(-Q/f+1)]);


        newSegment.w = numpy.sqrt(newSegment.wl/numpy.pi/-(1/newSegment.q).imag)  
        newSegment.R = 1/(1/newSegment.q).real
        return newSegment
    
    def transferBeam(self, res):
       
        for I in self.elementList:
            
            seg = self.segmentList[-1]
            
            seg.propagate(res,I.z-seg.z[0]);

            self.segmentList.append(self.lensTransferFunction(seg, I.f))
            
            
            
    def plotBeam(self):
        for I in self.segmentList:
            plt.plot(I.z, I.w, 'k')
            plt.plot(I.z, -I.w, 'k')
            #plt.show()
            #plt.hold(True)
        plt.show() 
    
    
            

class GaussianBeamSegment:
    """ this class contains a single segment of a gaussian beam. segments can only propagate in vacuum and
    are stitched together to define a GaussianBeam """
    

    wl = numpy.array([]) # wavelength
    q = numpy.array([]) # Gaussian beam parameter
    w0 = [] # waist size
    z0 = [] # origin of beam 
    
    z= numpy.array([]) # position array
    w = numpy.array([]) # beam waist size array
    R = numpy.array([]) # radius of curvature array
    
    
    
    
    def __init__(self, wl = 810e-9, z = 0, w = 50e-6, R = 1e10):
        """ Initialize a beam based one a position, a beam size at that position and the ROC.
        this is done to easily calculate the lens transfer"""
        
        self.wl = wl
        self.z = numpy.append(self.z,z)
        self.w = numpy.append(self.w,w)
        self.R= numpy.append(self.R,R)
        
        
    
    def calcInitialBeamParameters(self):
        self.w0 = self.w[0]/(1+(numpy.pi*self.w[0]**2/(self.wl*self.R[0]))**2)**0.5
        self.z0 = self.z[0] - self.R[0]/(1+(self.wl*self.R[0]/(numpy.pi*self.w[0]**2))**2)
        self.q = 1/(1/self.R - 1j*self.wl/numpy.pi/self.w**2)
        

        
        
        
    def propagate(self,dz, d):
        
        n = int(float(d)/dz)
        if n == 0:
            n =1
            
        for i in range(n):
            
            Q = self.q[i]
            self.q = numpy.append(self.q,(Q+dz))
            self.z = numpy.append(self.z,self.z[i]+dz)

        self.w = numpy.sqrt(self.wl/numpy.pi/-(1/self.q).imag)  
        self.R = 1/(1/self.q).real
            
            
        
    
    
    
if __name__ == "__main__":
    
    
    angle = 0.1;
    w = 660e-9/angle/numpy.pi;
    #w = 2.5e-6
    GS = GaussianBeamSegment(w =w, wl = 660e-9);
    GS.calcInitialBeamParameters()
    
    L1 =  Lens(z = 0.0075, f = 0.007 )
    L2 = Lens(z = 0.02, f = 0.2 )
    #L3 = Lens(z = 0.61, f = 1)


    GB= GaussianBeam()
    GB.elementList = [L1,L2]
    GB.segmentList.append(GS)
    GB.transferBeam(0.0001)
    GB.plotBeam()
    
    print("Waist size in crystal:")
    
    print(min(GB.segmentList[-2].w))


    
    