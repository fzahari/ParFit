#!/usr/bin/env python

from DihScan import DihScan
from IO import par_fit_inp

def pf_run(input_fname):
  
   gopt_s_fnameb,t1234,bes=ginp_inp(input_fname)

   ds=DihScan(gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,bes,t1234)
   ds.write_gamess_inputs()

default_input_fname="dih_scan_inp"

pf_run(default_input_fname)
