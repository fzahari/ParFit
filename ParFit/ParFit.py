#!/usr/bin/env python

from os import environ
from numpy import zeros,ndarray,random
from scipy.optimize import minimize,basinhopping,anneal,fmin,fmin_powell,fmin_cg,fmin_tnc
from DihScan import DihScan
from IO import par_fit_inp,read_add,write_add

from pyevolve import GSimpleGA
from pyevolve import G1DList
from pyevolve import Selectors
from pyevolve import Initializators, Mutators
from pyevolve import Consts
from pyevolve import DBAdapters

from deap import algorithms, base, creator, tools

def run_ga2(engine_rmse,np):

   creator.create("FitnessMax",base.Fitness,weights=(-1.0,))
   creator.create("Individual",ndarray,fitness=creator.FitnessMax)

   toolbox=base.Toolbox()
   toolbox.register("attr_float",random.random)
   toolbox.register("individual",tools.initRepeat,creator.Individual,toolbox.attr_float,n=np)
   toolbox.register("population",tools.initRepeat,list,toolbox.individual)
   toolbox.register("evaluate",engine_rmse)

   def cxComb(ind1,ind2):
      if random.randint(2)==0:
         tools.cxOnePoint(ind1,ind2)
      else:
         tools.cxBlend(ind1,ind2,0.0)
      return ind1,ind2

   toolbox.register("mate",cxComb)
   toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=3, indpb=1.)

   def checkBounds(min, max):
    def decorator(func):
        def wrapper(*args, **kargs):
            offspring = func(*args, **kargs)
            for child in offspring:
                for i in xrange(len(child)):
                    if child[i] > max:
                        child[i] = max
                    elif child[i] < min:
                        child[i] = min
            return offspring
        return wrapper
    return decorator

   MIN=-3.0
   MAX=3.0

   toolbox.decorate("mate", checkBounds(MIN, MAX))
   toolbox.decorate("mutate", checkBounds(MIN, MAX))

   toolbox.register("select",tools.selTournament,tournsize=3)

   pop=toolbox.population(n=50)
   algorithms.eaSimple(pop,toolbox,cxpb=0.3,mutpb=0.05,ngen=10)
   print(tools.selBest(pop,k=1)[0])

def run_ga(eval_func,n):
   
   genome = G1DList.G1DList(n)
   #genome.setParams(rangemin=-3.0, rangemax=3.0)
   genome.setParams(rangemin=-3.0, rangemax=3.0, bestrawscore=1.0, roundecimal=0)

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

def pf_run(input_fname):
  
   def engine_rmse(p):

      print p

      f=open("../Data/ParFit/step",'r')
      ls=f.readlines()
      f.close()
      step=int(ls[0])

      n=len(ds)
      rmse=0.
      for i in range(n):
         write_add(p,c,mm,ol_templ,lines,1,step,step_int)
         ds[i].run_dih_scan(p,c,mm,ol_templ)
         rmse+=ds[i].calc_rmse(csv,i,step,step_int)
      print step,rmse/n

      step+=1
      f=open("../Data/ParFit/step",'w')
      print >>f,step
      f.close()

      return rmse

   def engine_rmse2(p):

      print p

      f=open("../Data/ParFit/step",'r')
      ls=f.readlines()
      f.close()
      step=int(ls[0])

      n=len(ds)
      rmse=0.
      for i in range(n):
         write_add(p,c,mm,ol_templ,lines,1,step,step_int)
         ds[i].run_dih_scan(p,c,mm,ol_templ)
         rmse+=ds[i].calc_rmse(csv,i,step,step_int)
      print step,rmse/n

      step+=1
      f=open("../Data/ParFit/step",'w')
      print >>f,step
      f.close()

      return (rmse,)

   gopt_type,gopt_s_fnameb,t1234,bes,engine_path,mm,mode,alg,opt_lin,np,nc,step_int,csv=par_fit_inp(input_fname)

   n=len(gopt_type)
   ds=[]
   for i in range(n):
      sds=DihScan(gopt_s_fnameb[i],engine_path,mm,opt_lin,np,nc,bes[i],t1234[i])
      ds.append(sds)

   if not gopt_type=="ginp":
      environ["ENGINE_DIR"]=engine_path+"engine_dir"
      p,c,ol_templ,lines=read_add(mm,opt_lin,np,nc,1)

   f=open("../Data/ParFit/step",'w')
   print >>f,1
   f.close()
 
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
      elif alg=="ga2":
         run_ga2(engine_rmse2,np)
      elif alg=="fmin":
         #print fmin_powell(engine_rmse,p)
         print fmin(engine_rmse,p)
      else: 
         "'alg' is not a known algorithm!"
   print "inside"

default_input_fname="dih_scan_inp"

pf_run(default_input_fname)
