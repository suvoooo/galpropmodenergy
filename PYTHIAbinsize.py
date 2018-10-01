#!/usr/bin/python 
#the program is written to multiply the PYTHIA binsize with the dN/dE obtained from pythia
# this reduces dN/dE only to dN and then the new dN is used as calculation for dN/dE for GALPROP 
# binsize 



import math 
PYpath = '/home/suvo/Downloads/decayDM/'
Mcdm = input("Mass of DM:")
fM_DM = Mcdm/2.
def readPYTHIAfile():
	pyfile = open(PYpath+"4point-el-%d.dat"%(fM_DM),'r')
	pyvals = {}
	while True :
		stringline = pyfile.readline()
		if stringline =='' :
			break 
		stringlist = stringline.split()
		E = float(stringlist.pop(0))
		Flux = float(stringlist.pop(0))
		pyvals[E]=[Flux]
	return pyvals 

pyvals = readPYTHIAfile()
Epy = sorted(pyvals.keys())
len(Epy)
z = {}
bind =[]
for i in range(1,4001):
	z[i] = 10**(i*0.001)
	z[i-0.5] = 10**((i-0.5)*0.001)
	z[i+0.5] = 10**((i+0.5)*0.001)
	dEn = z[i+0.5] - z[i-0.5]
	En = z[i]
	print dEn
	bind.append(dEn)
dN = []
for i in range(len(Epy)):
	dNdE = pyvals[Epy[i]][0] * bind[i]
	dN.append(dNdE)
d_file = open(PYpath+"4point-el-%ddN.dat"%(fM_DM),'w')
for i in range(len(Epy)):
	if dN[i] != 0. :	
		d_file.write('{0:2} {1:3}' .format(Epy[i] , dN[i])+ "\n")
	else :
		pass
d_file.close()
