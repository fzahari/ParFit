#!/usr/bin/env python

import sys,os

dirs=sys.argv[1:]

if len(sys.argv)==1:
   print
   print "Help:"
   print
   print 'A) "./clean.py all" cleans up the entire content of the "../Data/ParFit" directory'
   print
   print 'B) "./clean.py dir1 dir2" ... " dirN" cleans up only the subdirectories "dir1", "dir2", ..., "dirN" of the ../Data/ParFit" directory'
   print
elif sys.argv[1]=="all":
   os.system("rm -rf ../Data/ParFit/*")
else:
   for d in dirs:
      os.system("rm -rf ../Data/ParFit/"+str(d)) 
