import matplotlib.pyplot as plt
 
from IPython.display import set_matplotlib_formats
import matplotlib
import numpy
import os



def saveFig(fig,filename= os.getcwd() + "\\out", mode = ".pdf"):
    
    fig.savefig(filename+mode, bbox_inches='tight')




def cm2inch(value):
    return value/2.54

def prepare(mode = "pdf"):
    set_matplotlib_formats('png', mode) # change to vector graphic pdf output
    matplotlib.rcParams.update({'errorbar.capsize': 3}) # errorbar cap sizes
    numpy.set_printoptions(threshold=numpy.nan) # not sure anymore, but important to avoid some errors ?!
    
def setFonts(size):  
    matplotlib.rc('font', family='sans-serif') 
    matplotlib.rc('font', serif='Helvetica Neue') 
    matplotlib.rc('text', usetex='false') 
    matplotlib.rcParams.update({'font.size': size})
    
def color(index):

    
    i = index%7;
    
    col =  [\
            (0.0, 0.0, 0.0), \
            (0.8, 0.3, 0.3), \
            (0.3, 0.7, 0.3), \
            (0.4, 0.4, 1.0), \
            (1.0, 1.0, 0.0), \
            (0.0, 1.0, 1.0), \
            (1.0, 0.0, 1.0), \
            ]
    return col[i]


def findNearestValue(arr, val):
    
    dist = numpy.inf
    for i in range(len(arr)):
        distNew = abs(arr[i]-val)
        if distNew < dist:
            dist = distNew
            index = i
            value = arr[i]
            
    return [int(index), value]
        
    
    
    
    