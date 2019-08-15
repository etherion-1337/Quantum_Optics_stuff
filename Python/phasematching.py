



#################################################################
##############	Supporting physics module	#####################
#################################################################

#################################################################
#### DIMENSIONS USED:
####	-Wavelength in [nm], because lambda is protected lam is used
####	-Dimensions in [mm]
####	-Time 		in [s] or [fs]
####	-Angles 	in	[] (rad)
##################################################################

import math
from scipy import interpolate
import numpy as np



def Sellmeier(coeff=[0,0,0,0],lam = 785):
	
	return (coeff[0]+(coeff[1])/((lam*1e-3)**2-coeff[2])-coeff[3]*(lam*1e-3)**2)**0.5


def phasefunction(W=0.1,lpump=405,l_range=785,theta_range=0,coeff=[2.7359, 0.01878, 0.01822, 0.01354, 2.3753, 0.01224, 0.01667, 0.01516],cutangle=28.7991*np.pi/180,crystal_length=6,
	show_plot=False,return_grid=False):
	#if given a list of [lambdamin,lambdamax] make list accordingly
	if isinstance(l_range,list):
		lsignal=l_range[2]*np.array(list(range(round((1/l_range[2])*l_range[0]),round((1/l_range[2])*l_range[1]),1)))
		lamsize=len(lsignal)

	else:
		lsignal=np.array([l_range])
		lamsize=1
		
	#if given a list of [theta_min,theta_max] make list accordingly
	if isinstance(theta_range,list):
		theta_o=(np.pi/180)*(theta_range[2])*np.array(list(range(round(1/(theta_range[2])*theta_range[0]),round(1/(theta_range[2])*theta_range[1]),1)))
		thetasize=len(theta_o)
		
	else:
		theta_o=np.array([theta_range])
		thetasize=1
		

	c=coeff

	lp=lpump
	ls=lsignal
	li=ls*lp/(ls-lp);

	wp,ws,wi = (2*np.pi)/lp,(2*np.pi)/ls,(2*np.pi)/li

	nop = Sellmeier(coeff=c[0:4],lam=lp)
	nep = Sellmeier(coeff=c[4:8],lam=lp)
	nos = Sellmeier(coeff=c[0:4],lam=ls)
	nes = Sellmeier(coeff=c[4:8],lam=ls)
	noi = Sellmeier(coeff=c[0:4],lam=li)
	nei = Sellmeier(coeff=c[4:8],lam=li)

	nairs = 1+(0.05792105)/(238-(ls*1e-3)**-2)+(0.00167917)/(57.362-(ls*1e-3)**-2)
	thet_cut=cutangle
	npeff=np.sqrt(1/((np.cos(thet_cut)/nop)**2+(np.sin(thet_cut)/nep)**2))
	L=crystal_length



	wsmesh,theta_omesh=np.meshgrid(ws,theta_o)
	wimesh,theta_omesh=np.meshgrid(wi,theta_o)

	thet_s=theta_omesh
	thet_i=np.arcsin(nos*wsmesh*np.sin(thet_s)/(noi*wimesh))
	dkz=npeff*wp-nos*wsmesh*np.cos(thet_s)-noi*wimesh*np.cos(thet_i)
	dky=-nos*wsmesh*np.sin(thet_s)+noi*wimesh*np.sin(thet_i); print (dky, dkz)
	phi=np.exp(-((W*1e6)**2)*(dky**2)/2)*(np.sin(0.5*dkz*L*1e6)/(0.5*dkz*L*1e6))**2;# print(thet_s, thet_i)

	result = phi

	if show_plot:
		
		if lamsize == 1:
			plt.figure()
			plt.plot(theta_o*(180/np.pi),result)
			plt.xlabel('opening angle [*]')
			plt.show()

		elif thetasize==1:

			plt.figure()
			plt.plot(lsignal,result[0])
			plt.xlabel("wavelength [nm]")
			plt.show()

		else:
			plt.figure()
			CS=plt.contourf(lsignal,(180/np.pi)*theta_o,result,levels=np.linspace(0,1,20))
			plt.xlabel('wavelength [nm]')
			plt.ylabel('opening angle [*]')
			plt.show()
	
	if not show_plot:
		if not return_grid:
			if lamsize == 1 and thetasize == 1:
				return (result)
			else:
				return result
		else:
			return (lsignal,theta_o,result)


def phasefunction_vec(lpump=405,l_list=[785],theta_list=[0],coeff=[2.7359, 0.01878, 0.01822, 0.01354, 2.3753, 0.01224, 0.01667, 0.01516],cutangle=28.7991*np.pi/180,crystal_length=6):

	c=coeff

	lp=lpump
	ls=l_list
	li=ls*lp/(ls-lp);

	theta_o=theta_list

	wp,ws,wi = (2*np.pi)/lp,(2*np.pi)/ls,(2*np.pi)/li

	nop = Sellmeier(coeff=c[0:4],lam=lp)
	nep = Sellmeier(coeff=c[4:8],lam=lp)
	nos = Sellmeier(coeff=c[0:4],lam=ls)
	nes = Sellmeier(coeff=c[4:8],lam=ls)
	noi = Sellmeier(coeff=c[0:4],lam=li)
	nei = Sellmeier(coeff=c[4:8],lam=li)

	npeff=np.sqrt(1/((np.cos(cutangle)/nop)**2+(np.sin(cutangle)/nep)**2))
	L=crystal_length
	W=0.1

	thet_s=theta_o
	thet_i=np.arcsin(nos*ws*np.sin(thet_s)/(noi*wi))
	dkz=npeff*wp-nos*ws*np.cos(thet_s)-noi*wi*np.cos(thet_i)
	dky=-nos*ws*np.sin(thet_s)+noi*wi*np.sin(thet_i)
	phi=np.exp(-((W*1e6)**2)*(dky**2)/2)*(np.sin(0.5*dkz*L*1e6)/(0.5*dkz*L*1e6))**2

	return phi

	


