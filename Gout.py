#!/usr/bin/env python3

import sys
from Scan import DihAScan,BondScan,AnglScan
from _IO import ginp_inp

def pf_run(input_fname):
  
   scan_type,gopt_s_fnameb,tup,bes=ginp_inp(input_fname)

   if scan_type=="diha":
      bes=list(map(int,bes))
      ds=DihAScan(None,gopt_s_fnameb,None,"mm3",None,None,None,bes,tup)
      ds.read_gamess_outputs()
      ds.write_gouts_data()
   elif scan_type=="bond":
      bes=list(map(lambda x:10.*x,bes))
      bes=list(map(int,bes))
      ds=BondScan(None,gopt_s_fnameb,None,"mm3",None,None,None,bes,tup)
      ds.read_gamess_outputs()
      ds.write_gouts_data()
   elif scan_type=="angl":
      bes=list(map(lambda x:10.*x,bes))
      bes=list(map(int,bes))
      ds=AnglScan(None,gopt_s_fnameb,None,"mm3",None,None,None,bes,tup)
      ds.read_gamess_outputs()
      ds.write_gouts_data()

GO_input_fname="gout_inp"

lsa=len(sys.argv)
if lsa>2:
   print()
   print('Use: "./Gout.py name_of_Gout_input_file"')
   print()
elif lsa==2:
   GO_input_fname=h=sys.argv[1]

pf_run(GO_input_fname)
