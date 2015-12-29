#!/usr/bin/env python
#
# This program plots the energy profiles of MM derived 
# energies after optimizing the parameters 
# using ParFit. 
# The program must be run from the directory containing
# the csv data file to plot. Change the filename and
# any variables, such as labels, and tiles.

# *****************************************************
# *  User: change file name to the csv file to plot.  *
# *  Note: file must be in current directory.         *
# *****************************************************
f=open("opt0_31.csv",'r')

lines=f.readlines()
f.close()

import numpy as np
import matplotlib.pyplot as plt

X=[]
Y=[]
Z=[]

for line in lines:
	x,y,z=line.split(',')
	x,y,z=float(x),float(y),float(z)
	X.append(x)
	Y.append(y)
	Z.append(z)

# *****************************************************
# *  Plot legend, range, and tick mark spacing.       *
# *****************************************************
plt.plot(X,Y,"o",label="MP2 energy")
plt.plot(X,Z,label="MM3 energy")
plt.xticks(range(0,360,30))
plt.legend(loc="lower right")

# *****************************************************
# *  Plot title and attributes.                       *
# *****************************************************
plt.title('Title: MM3 fit to MP2',fontsize=18,fontname="Times Roman",color="brown")

# *****************************************************
# *  Plot axis labels, attributes, and axis limits.   *
# *****************************************************
plt.xlabel("X-axis label: Dihedral Angle (deg)",fontsize=14,fontname="Times Roman",color="black")
plt.ylabel("Y-axis label: Relative energy (kcal/mol)",fontsize=14,fontname="Times Roman",color="black")
plt.xlim([0.,360.])
plt.ylim([-2.5,1.0])

plt.show()
