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
import math

#aspherics choice for pump collimation for a Ondax TO-660-35
def colli_focal(beam_size, beam_divergence): #returned value same unit as beam_size, beam_divergence is full divergence angle in degree.
    focal_length = (beam_size/2)/(math.tan(math.radians(beam_divergence/2)));
    return focal_length;

def NA_LD (LD_divergence): #in degree, full angle divergence
    LD_NA = math.sin(math.radians(LD_divergence/2));
    return LD_NA;

def f_hex (focal_length, lens_clear_aperture): #all units, including return, in mm
    lens_f_hex = focal_length/lens_clear_aperture;
    return lens_f_hex;

def NA_lens (focal_length, lens_clear_aperture): #all units, including return, in mm
   lens_NA = 1/(2*(f_hex(focal_length, lens_clear_aperture)));
   return lens_NA;

pump_full_divergence = 20; #in degree
target_colli_beam_size = 0.6; #in mm

pump_NA = NA_LD(pump_full_divergence);
pump_collimation_NA = 2*pump_NA; #NA of collimation lens taken to be twice of the LD's NA to collect everything from the LD.
print ("Pump LD's NA is %8.8f" %(pump_NA));
print ("Pump collimation lens's NA should be at least %5.5f" %(pump_collimation_NA));
print("###########################################")

pump_colli_focal_length = colli_focal(target_colli_beam_size, pump_full_divergence);
pump_colli_h_hex = f_hex(target_colli_beam_size, pump_full_divergence);
print ("Pump collimation lens's focal length is %5.5f" %(pump_colli_focal_length)); #larger focal length -> wider beam
print ("Pump collimation lens's f/# is %5.5f" %(pump_colli_h_hex));

lens_NA_pump = NA_lens(1.7, 3.7);
if lens_NA_pump <= pump_collimation_NA:
    print ("This pump collimation lens ONLY has NA of %5.5f, choose a lens with LARGER NA" %(lens_NA_pump));
else :
    print ("This pump collimation lens has NA of %5.5f, good enough" %(lens_NA_pump));

print ("###########################################")
#Coupling freespace Ondax TO-660-35 into SMF
def colli_couple_fiber_focal_length (wavelength, beam_waist, MFD):
    couple_focal_length = (1.7*beam_waist)*(math.pi*MFD)/(4*wavelength);
    return couple_focal_length;

#For fiber SM600 at Ondax TO-660
pump_wavelength = 0.00066 # in mm
pump_beam_waist = 0.6 # in mm
MFD = 0.0045 # in mm
pump_couple_fiber_lens_focal_length = colli_couple_fiber_focal_length(pump_wavelength,pump_beam_waist, MFD);
print("Pump collection lens from collimated beam should have focal length of %5.5f" %(pump_couple_fiber_lens_focal_length)); #in mm

#Focus pump from SMF into PPKTP (target waist = 150 um)

