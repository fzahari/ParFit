#!/usr/bin/env python

f=open("gplot0",'r')
lines0=f.readlines()
f.close()

f=open("gplot1",'r')
lines=f.readlines()
f.close()

import numpy as np
import matplotlib.pyplot as plt

X=[]
Y=[]
Z=[]
Z0=[]

for line in lines0:
	x,y,z=line.split()
	x,y,z=float(x),float(y),float(z)
	Z0.append(z)

for line in lines:
	x,y,z=line.split()
	x,y,z=float(x),float(y),float(z)
	X.append(x)
	Y.append(y)
	Z.append(z)

#plt.plot(X,Y,"o",X,Z)
plt.plot(X,Y,"o",label="MP2 energy")
plt.plot(X,Z0,label="MM3 energy (before optimization)")
plt.plot(X,Z,label="MM3 energy (after optimization)")
plt.xticks(range(0,370,30))
plt.title('MM3 fit to MP2/cc-pVTZ (Phenyl(dimethyl)phosphine Oxide)',fontsize=18,fontname="Times Roman",color="brown")
plt.xlabel("C=C-P-O dihedral angle (deg)",fontsize=14,fontname="Times Roman",color="black")
plt.ylabel("Relative energy (kcal/mol)",fontsize=14,fontname="Times Roman",color="black")
plt.legend(loc="lower right")
plt.xlim([0.,360.])
plt.ylim([-2.0,4.0])
plt.show()
