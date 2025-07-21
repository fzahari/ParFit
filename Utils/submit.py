#!/usr/bin/env python3

from os import system

comm1="~/tools/xgms2 -q ge -l acet_snap"
comm2=".log acet_snap"

for i in range(50):
	n=500+(i+1)*10
	system(comm1+str(n)+comm2+str(n))
