#!/usr/bin/env python

from os import environ
from numpy import zeros
from scipy.optimize import minimize,basinhopping,anneal,fmin,fmin_powell,fmin_cg,fmin_tnc
from DihScan import DihScan
from IO import par_fit_inp,read_add,write_add

from pyevolve import GSimpleGA
from pyevolve import G1DList
from pyevolve import Selectors
from pyevolve import Initializators, Mutators
from pyevolve import Consts
from pyevolve import DBAdapters

def run_ga(eval_func,n):
   
   genome = G1DList.G1DList(n)
   genome.setParams(rangemin=-3.0, rangemax=3.0)

   genome.initializator.set(Initializators.G1DListInitializatorReal)

   genome.mutator.set(Mutators.G1DListMutatorRealGaussian)

   genome.evaluator.set(eval_func)

   ga = GSimpleGA.GSimpleGA(genome)
   ga.setMinimax(Consts.minimaxType["minimize"])
   ga.setPopulationSize(20)
   ga.selector.set(Selectors.GRouletteWheel)
   ga.setGenerations(40)
   
   sqlite_adapter = DBAdapters.DBSQLite(identify="fit3b", resetDB=True)
   ga.setDBAdapter(sqlite_adapter)

   ga.evolve(freq_stats=1)

   bi=ga.bestIndividual()
   print bi
   return bi

default_input_fname="dih_scan_inp"

gopt_type,gopt_s_fnameb,t1234,bes,engine_path,mm,mode,alg,opt_lin,np,nc,csv=par_fit_inp(default_input_fname)
ds=DihScan(gopt_s_fnameb,engine_path,mm,bes,t1234)
if not gopt_type=="ginp":
   environ["ENGINE_DIR"]=engine_path+"engine_dir"
   p,c,ol_templ,lines=read_add(mm,opt_lin,np,nc,1)
#p[0]=p[2]=0.0
#write_add(p,c,mm,ol_templ,lines)

def engine_rmse(p):
    #print p
    write_add(p,c,mm,ol_templ,lines,0)
    ds.run_dih_scan(p,c,mm,ol_templ,lines)
    rmse=ds.calc_rmse(csv)
    print rmse
    return rmse

if gopt_type=="full":
    ds.read_gamess_outputs()
elif gopt_type=="comp":
    ds.read_gouts_data()
elif gopt_type=="ginp":
    ds.write_gamess_inputs()
    quit()
else:
    print "\nPar_Fit: Wrong run type. \nAccepted values: full, comp, and ginp. \nPlease check input file.\n\n"
ds.write_engine_inputs()

if mode=="sense":
   eps=0.01
   np=len(p)
   ppM=zeros(np)
   ppP=zeros(np)
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
      run_ga(engine_rmse,np)
   elif alg=="fmin":
      print fmin(engine_rmse,p)
   else: 
      "'alg' is not a known algorithm!"
