#!/usr/bin/env python3

from os import system

fnameb="opmmph-mp2-popt-dd-"

for n in range(0,185,5):
    sn0=str(n)
    if (n<10):
        sn="00"+sn0
    elif (n<100):
        sn="0"+sn0
    else:
        sn=sn0
    fname1=fnameb+sn+"-dlc.log"
    fname2=fnameb+sn+".log"
    system("mv "+fname1+" "+fname2)
