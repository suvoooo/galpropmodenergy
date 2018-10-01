#!/usr/bin/python 

# program that generates the input spectrum from PYTHIA for GALPROP 
# the code is similar to the dNdEinputforgal.py	
# but it's for a decaying DM of type ll\nu _R .. so the lepton spectrum from PYTHIA is multiplied with two
# right handed neutrino doesn't contribute to the electron spectrum ... 
# bins used for PYTHIA as E[i]=10**(i*0.001) ...so to get i from here we use, i = (10**3)*(log10(E)) 

import math 
import numpy 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+
#+                Definition of NFW profile 
#+
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# R_vir = virial radius
# M_vir = virial mass
# C_v = ratio of virial radius to scale radius 
# rho = density of the halo 
# BF = boost factor
Rsun = 8.33#kpc
SigmaV = 3.2e-26 #(cm^3/Sec)
Taurate = 3.0e-26
M_DM = 4000000 #1.0e3*(1.1**63) #1.0e3*(1.2**33)#MeV
fM_DM = M_DM/2000.
R_vir = 200 # Kpc
M_vir = 1e12 
pi = 3.14
factor = 3.81718e-8
rhomax = 1e18
alpha = 0.17
#BF = 1 #boost factor
#if BF == 1 :
#	print "Boost factor is 1"
print "Mass of the DM is %dGeV" %(M_DM/1000.)
print "! beware of type of DM !"
def rhosq(x,y,z) :
	R_sph = math.sqrt(x**2 + z**2 + y**2)
	if R_sph <= 0.1 :
		R_sph = 0.1
# the effect of cut-off radius has to be checked 
	else :
		R_sph = R_sph 
	C_v = 10.0 
	logt = math.log(1 + C_v) - (C_v/(1 + C_v))
	R_s = R_vir/C_v
	x = R_sph/R_s
	rho_s = M_vir/((4 * (pi) * (R_s)** 3)* logt) 
	rho = min(rhomax, rho_s / (x*((1 + x)**2))) *factor
	# corrected formula >>>> not x*(1 + x**2)
	return rho 

def flatrho(x,y,z):
	rho = 0.3 
	return rho

def einastorhosq(r,z) : 
	return math.exp(-2 * ( math.pow ((r/Rsun), alpha) -1 ) / alpha ) * factor 

#*****************************************************************
#*
#*                 read PYTHIA dN/dE *dE file
#*
#****************************************************************** 
PYpath = '/home/suvo/Downloads/decayDM/'
def readpyvals():
	pyfile =open(PYpath+"4point-tau-%ddN.dat"%(fM_DM),'r')
	pyvals = {}
	while True :
		stringline = pyfile.readline()
		if stringline =='':
			break 
		stringlist = stringline.split()
		E = float(stringlist.pop(0))*(10**3) 
		Flux = float(stringlist.pop(0))
		pyvals[E] = [Flux] 
	return pyvals 
pyvalsdict = readpyvals()
Epy = sorted(pyvalsdict.keys())
print pyvalsdict[Epy[3]][0]
print Epy[-16]
start_energy = input("write the start energy value of GALPROP :")
energy_factor = input("write the energy factor as in Galdef file :")
bin_number = input("write the number of energy bin used in GALPROP :")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +    Select Bin Boundaries
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Egalp = []
for i in range(bin_number):
	E = start_energy*(energy_factor**i)
	Egalp.append(E)
	if E > Epy[-1]: # this makes sure that galprop bins before and after the Mass of dark matter(=Emic[-1]) is selected, for example 800 GeV dark #matter, both the galprop bins before and after the dark matter mass will be included in the list Egalp. Usage of this is demonstrated later .  
		break 
print Egalp 
interpE = Egalp[-2]
#for i in range(len(Epy)):
#	if abs(Epy[-1] - Egalp[-2]) >=  abs(Egalp[-1] - Epy[-1]) :
#		interpE = Egalp[-1]
#	else :
		#interpE = Egalp[-2]	

print interpE
for i in range(len(Egalp)):
	if Egalp[i]== interpE:
#		print i  
		gali =i
Epy_last = start_energy*(energy_factor** (gali-0.5) ) 
print Epy_last

Epylast = []
for i in range(len(Epy)):
	if Epy[i] >= Epy_last :
		Epylast.append(Epy[i])  # this helps to select the last remaining energy bins
print len(Epylast) 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                Linear interpolation between PYTHIA bins and GALPROP bins
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def interpol(Egal,Epy_last):
	finalsum = 0. 
	intflux = 0.
	if Egal> Epy_last :
		 for lo in range(len(Epylast)):
			lo = lo +1 
			finalsum += pyvalsdict[Epy[-lo]][0]
			#print finalsum
		 return finalsum , True 
 
	for i in range(len(Epy)):
		E2 = Epy[i]
		if E2 > Egal :
			E1 = Epy[i-1]
			var = E2 /E1
			#print E2, E1
			index = ( math.log(pyvalsdict[E2][0]/(pyvalsdict[E1][0])) ) /(math.log(var))
			factor = (E2 -E1)
			#print index
			#if  index <= -40.: 		
			#	intflux = 0.
			#else :	
			intflux = pyvalsdict[E1][0] + ( ( pyvalsdict[E2][0] - pyvalsdict[E1][0] )*(Egal -E1) /(factor) ) 
			return intflux , False
        print "wrong" 

spclbinsize = (M_DM/2) - Epy_last
print spclbinsize
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +           Creating source file for GALPROP 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++




Ega={}
d_file = open(PYpath+'sourceDMtau.dat','w')
for i in range(-20,21,1):
	x = i/1
	for j in range(-20,21,1):
		y = j/1
		for k in range(-7,8,1):
			z = k/1 
			ltdm = False
			remain_bins = 500
			for l in range(bin_number):
				Ega[l] = start_energy*(energy_factor**l)
				#Ega[l-0.5] = start_energy*(energy_factor**(l-0.5))
				#Ega[l+0.5] = start_energy*(energy_factor**(l+0.5))
				#upbin = Ega[l-0.5]
				#lowbin = Ega[l+0.5]
				#binsize = lowbin - upbin
				ipy = (10**3)*math.log10(Ega[l])
				newbinsize = 10**( ( ipy +0.5 )*0.001 ) - 10**(( ipy -0.5 )*0.001 ) 
				
				if ltdm and remain_bins == 0 :
					flux=0.0
					d_file.write('{0:2}' .format(flux) + "\n")
					#data_file.write('{0:2f} {1:3} ' .format(Ega[l],source) + "\n")
#					print ltdm
					continue
				elif  ltdm : 
					pass 
				else:
#					print ltdm
	
					intpflx,ltdm=interpol(Ega[l], Epy_last)
					remain_bins= 1
				if ltdm and remain_bins > 0 :
						
					#print intpflx,ltdm
					#intpflux = interpol(Ega[l])
					sumflx = intpflx
					flux = 1.*(abs ( 1.0e3 * (Taurate) * sumflx * (math.pow((rhosq(x,y,z)/M_DM),1)/(spclbinsize)) ) )
					d_file.write('{0:2}' .format(flux) + "\n")
					remain_bins = remain_bins -1 
					
					#print Ega[l] , flux

				else :
					lastflx = intpflx
					flux = 1.*(abs ( 1.0e3 * (Taurate) * intpflx * (math.pow((rhosq(x,y,z)/M_DM),1)/(newbinsize)) ) ) 
					print newbinsize
					d_file.write('{0:2}' .format(flux) + "\n")
					#print Ega[l] , flux







d_file.close()		
#++++++++++++++++++++++++++++++++++++++++++
# + check interpolation
#++++++++++++++++++++++++++++++++++++++++++

#for i in range(bin_number):
#	Ega[i]=start_energy*(energy_factor**i)
#	interpol(Ega[i])
#	#print Ega[i], interpol(Ega[i])
#'''
# 
