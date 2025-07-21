#!/usr/bin/env python3

import sys
from Scan import DihAScan,BondScan,AnglScan
from _IO import ginp_inp

def pf_inp_run(input_fname):
  
   prog_type,scan_type,gopt_s_fnameb,tup,bes=ginp_inp(input_fname)

   if scan_type=="diha":
      bes=list(map(int,bes))
      ds=DihAScan(None,gopt_s_fnameb,None,"mm3",None,None,None,bes,tup)
      if prog_type=="nwchem":
         ds.write_nwchem_inputs()
      else:
         ds.write_gamess_inputs()
   elif scan_type=="bond":
      bes=list(map(lambda x:10.*x,bes))
      bes=list(map(int,bes))
      ds=BondScan(None,gopt_s_fnameb,None,"mm3",None,None,None,bes,tup)
      if prog_type=="nwchem":
         ds.write_nwchem_inputs()
      else:
         ds.write_gamess_inputs()
   elif scan_type=="angl":
      bes=list(map(lambda x:10.*x,bes))
      bes=list(map(int,bes))
      ds=AnglScan(None,gopt_s_fnameb,None,"mm3",None,None,None,bes,tup)
      if prog_type=="nwchem":
         ds.write_nwchem_inputs()
      else:
         ds.write_gamess_inputs()

GI_input_fname="ginp_inp"

lsa=len(sys.argv)
if lsa>2:
   print()
   print('Use: "./Ginp.py name_of_Ginp_input_file"')
   print()
elif lsa==2:
   GI_input_fname=h=sys.argv[1]

pf_inp_run(GI_input_fname)
