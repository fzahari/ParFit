#!/usr/bin/env python

import sys,os
import numpy

from scipy.optimize import minimize,basinhopping,anneal,fmin,fmin_powell,fmin_cg,fmin_tnc
from Scan import BondScan,AnglScan,DihAScan
from IO import par_fit_inp,read_add,write_add

from Ga import run_ga

def pf_run(PF_if):

   scan_type,gopt_type,gopt_s_fnameb,tup,bes,engine_path,mm,mode,alg,opt_lin,np,nc,step_int,csv=par_fit_inp(PF_if)
   #scan_type=scan_type.strip()   

   if engine_path=="":
      engine_path="../Engine"

   sdir=[]
   pref="../Data/ParFit/"+PF_if
   sdir.append(pref)
   if os.path.exists(pref):
      print
      print 'Warning: The directory',pref,'exists!'
      print
      sys.exit()
   os.mkdir(pref)   
   for gsf in gopt_s_fnameb:
      sd_gsf=pref+"/"+gsf
      sdir.append(sd_gsf)
      os.mkdir(sd_gsf)

   n=len(gopt_type)
   ds=[]
   if scan_type=="diha":
      for i in range(n):
         if not bes[i]==():
            bes[i]=map(int,bes[i])
         sds=DihAScan(sdir,gopt_s_fnameb[i],engine_path,mm,opt_lin,np,nc,bes[i],tup[i])
         ds.append(sds)
   elif scan_type=="bond":
      for i in range(n):
         if not bes[i]==():
            bes[i]=map(lambda x:10.*x,bes[i])
            bes[i]=map(int,bes[i])
         sds=BondScan(sdir,gopt_s_fnameb[i],engine_path,mm,opt_lin,np,nc,bes[i],tup[i])
         ds.append(sds)
   elif scan_type=="angl":
      for i in range(n):
         if not bes[i]==():
            bes[i]=map(lambda x:10.*x,bes[i])
            bes[i]=map(int,bes[i])
         sds=AnglScan(sdir,gopt_s_fnameb[i],engine_path,mm,opt_lin,np,nc,bes[i],tup[i])
         ds.append(sds)

   if not gopt_type[0]=="ginp":
      os.environ["ENGINE_DIR"]=engine_path+"/engine_dir"
      p,c,ol_templ,lines=read_add(mm,opt_lin,np,nc,1,scan_type)

   for i in  range(n):
      if gopt_type[i]=="full":
         ds[i].read_gamess_outputs()
      elif gopt_type[i]=="comp":
         ds[i].read_gouts_data()
      elif gopt_type[i]=="ginp":
         ds[i].write_gamess_inputs()
         quit()
      else:
         "Par_Fit: Wrong gopt_type!"

      ds[i].write_engine_inputs()

   def engine_rmse(p):

      n=len(ds)
      rmse=0.
      for i in range(n):
         write_add(sdir,p,c,mm,ol_templ,lines,1,engine_rmse.step,step_int)
         ds[i].run_scan(p,c,mm,ol_templ)
         rmse+=ds[i].calc_rmse(csv,i,engine_rmse.step,step_int)
      print engine_rmse.step,round(rmse/n,4),p

      engine_rmse.step+=1

      return round(rmse/n,4)
  
   engine_rmse.step=1

   def engine_rmse2(p):

      n=len(ds)
      rmse=0.
      for i in range(n):
         write_add(sdir,p,c,mm,ol_templ,lines,1,engine_rmse2.step,step_int)
         ds[i].run_scan(p,c,mm,ol_templ)
         rmse+=ds[i].calc_rmse(csv,i,engine_rmse2.step,step_int)

      engine_rmse2.step+=1

      return (round(rmse/n,4),)

   engine_rmse2.step=1

   if mode=="sense":
      eps=0.01
      np=len(p)
      ppM=numpy.zeros(np)
      ppP=numpy.zeros(np)
      p0=engine_rmse(p)
      for i in range(np): 
         for j in range(np):
            ppM[j]=p[j]
            ppP[j]=p[j]
         ppM[i]=ppM[i]-eps
         ppP[i]=ppP[i]+eps
         pM=engine_rmse(ppM)
         pP=engine_rmse(ppP)
         print i+1,(pM-p0)/eps,(pP-p0)/eps
   else:
      np=len(p)
      if alg=="ga":
         print "Warning: The genetic algorithm printout will not start immediately!"
         hof=run_ga(engine_rmse2,np,40)
         hof0=numpy.array(hof[0])
         write_add(sdir,hof0,c,mm,ol_templ,lines,1,"ga",None)
      elif alg=="fmin":
         #print fmin_powell(engine_rmse,p)
         print fmin(engine_rmse,p,ftol=0.2)
      elif alg=="hybr":
         print "Warning: The genetic algorithm printout will not start immediately!"
         hof=run_ga(engine_rmse2,np,40)
         hof0=numpy.array(hof[0])
         write_add(sdir,hof0,c,mm,ol_templ,lines,1,"ga",None)
         fmin(engine_rmse,hof0,ftol=0.2)
      else: 
         "'alg' is not a known algorithm!"

PF_input_fname="scan_inp"

lsa=len(sys.argv)
if lsa>2:
   print 
   print 'Use: "./ParFit.py name_of_ParFit_input_file"'
   print 
elif lsa==2:
   PF_input_fname=h=sys.argv[1]

pf_run(PF_input_fname)
