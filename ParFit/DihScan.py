#!/usr/bin/env python

from os import system,environ
from numpy import rad2deg,array,sqrt,pi,around
from IO import DihGOpt_Molecule,par_fit_inp,read_add,write_add
from _Engine import run_engine_timeout, pert_add_param
from GeomStr import default_mm3_type,default_mmff94_type

import copy_reg
import types

def _reduce_method(meth):
    return (getattr,(meth.__self__,meth.__func__.__name__))
copy_reg.pickle(types.MethodType,_reduce_method)

class DihScan(object):
    hartree2kcal_mol=627.509469
    def __init__(self,gopt_s_fnameb,engine_path,mm,da_range_tuple,dih_tuple=(),gopt_b_fnameb="mp2_base"):
        self.gopt_base_fnameb=gopt_b_fnameb
        self.engine_path=engine_path
        self.gopt_scan_fnameb=gopt_s_fnameb
        self._mm=mm
        assert len(da_range_tuple)==3
        self._da_rt=da_range_tuple
        if not dih_tuple==():
            assert len(dih_tuple)==4
        self._dt=tuple(map(lambda x:x-1,dih_tuple))
        self._ml=[]
        self._dad={}
        self._ged={}
        self._eed={}

    @property
    def mm(self):
        return self._mm

    @property
    def da_rt(self):
        return self._da_rt

    @property
    def dt(self):
        return self._dt

    @property
    def ml(self):
        return self._ml

    @property
    def dad(self):
        return self._dad

    @property
    def ged(self):
        return self._ged

    @property
    def eed(self):
        return self._eed

    def _gen_geometries(self):
        #fnameb=self.gopt_scan_fnameb
        #da_rt=self.da_rt
        #dgomol=DihGOpt_Molecule(dih_tuple=self.dt,mm=self.mm,pairs_tuple=self.pairs_tuple,fifties_tuple=self.fifties_tuple,name=fnameb)
        #print "inside _gen_geometries: ",fnameb,da_rt
        pass

    def write_gamess_inputs(self):
        #self._gen_geometries()
        fnameb=self.gopt_scan_fnameb
        da_rt=self.da_rt
        dgomol=DihGOpt_Molecule(dih_tuple=self.dt,mm=self.mm,name=fnameb)
        dgomol.read_ginp(fnameb)
        t1,t2,t3,t4=self.dt
        dar=dgomol.calc_dihedral(t1,t2,t3,t4)
        #dag=(dar/pi)*180.0
        #print dar, dag, int(around(dag,0))
        dgomol.dih_rot(-dar)
        #dar=dgomol.calc_dihedral(t1,t2,t3,t4)
        #dag=(dar/pi)*180.0
        #print dar, dag, int(around(dag,0))
        b,e,s=da_rt
        e+=1
        na=dgomol.na
        for n in range(b,e,s):
           sn0=str(n)
           if (n<10):
              sn="00"+sn0
           elif (n<100):
              sn="0"+sn0
           else:
              sn=sn0
           f=open("../Data/Gamess/"+fnameb+"-"+sn+".inp",'w')
           print >>f," $ZMAT DLC=.T. AUTO=.T. $END"
           print >>f," $ZMAT IFZMAT(1)=3,",t1+1,",",t2+1,",",t3+1,",",t4+1," FVALUE(1)=",float(n),"$END"
           print >>f," $CONTRL COORD=UNIQUE NZVAR=",3*na-6,"$END"
           for line in dgomol._ginp_templ:
              print >>f,line[:-1]
           dgomol.dih_rot(float(n))
           rl=dgomol.rl
           cl=dgomol.cl
           sl=dgomol.sl
           for i in range(na):
              x,y,z=rl[i]
              print >>f,"",sl[i],cl[i],x,y,z
           print >>f," $END"
           f.close()

    def read_gamess_outputs(self):
        self._ged={}
        b,e,s=self._da_rt
        e+=1
        for n in range(b,e,s):
            sn0=str(n)
            if (n<10):
                sn="00"+sn0
            elif (n<100):
                sn="0"+sn0
            else:
                sn=sn0
            fnameb=self.gopt_scan_fnameb+sn
            dgomol=DihGOpt_Molecule(dih_tuple=self._dt,mm=self._mm,name=fnameb)
            dgomol.read_gopt_log(fnameb)
            self._ml.append(dgomol)
            print "dogmol.da:",dgomol.da
            self._dad.update({dgomol.name:rad2deg(dgomol.da)})
            self._ged.update({dgomol.name:dgomol.e})

    def write_gouts_data(self):
        fnameb=self.gopt_scan_fnameb
        f=open("../Data/Gamess/"+fnameb+"scan",'w')
        t1,t2,t3,t4=self._dt
        print >>f,t1,t2,t3,t4
        b,e,s=self._da_rt
        print >>f,b,e,s
        l=len(self._ml)
        print >>f,l
        for i in range(l):
            m=self._ml[i]
            n=m._n
            na=m._na
            print >>f,n,na,self._dad[n],self._ged[n]
            for j in range(na):
                x,y,z=m.rl[j]
                print >>f,m.sl[j],x,y,z,m.cl[j]
        f.close()

    def read_gouts_data(self):
        fnameb=self.gopt_scan_fnameb
        f=open("../Data/Gamess/"+fnameb+"scan",'r')
        lines=f.readlines() 
        f.close()
        self._dt=map(int,lines[0][:-1].split())
        self._da_rt=map(int,lines[1][:-1].split())
        l=int(lines[2][:-1])
        lines=lines[3:]
        self._dad={}
        self._ged={}
        self._ml=[]
        for i in range(l):
            fnameb,na,d,e=lines[0][:-1].split()
            na=int(na)
            dgomol=DihGOpt_Molecule(dih_tuple=self._dt,mm=self._mm,name=fnameb)
            dgomol._sl=[]
            dgomol._rl=[]
            dgomol._cl=[]
            dgomol._tl=[]
            lines=lines[1:]
            for j in range(na):
                s,x,y,z,c=lines[j].split()
                dgomol._sl.append(s)
                dgomol._rl.append(array(map(float,[x,y,z]),'d'))
                dgomol._cl.append(float(c))
                if self._mm=="mm3":
                    dgomol._tl.append(default_mm3_type[s])
                elif self._mm=="mmff94":
                    dgomol._tl.append(default_mmff94_type[s])
                else:
                    print "DihGOpt_Molecule.read_gopt_log: Wrong MM-type!"
            dgomol._rl=array(dgomol._rl,'d')
            dgomol._na=len(dgomol._rl)
            dgomol.set_conn()
            self._ml.append(dgomol)
            d,e=map(float,[d,e])
            self._dad.update({fnameb:d})
            self._ged.update({fnameb:e})
            lines=lines[na:]

    def write_engine_inputs(self):
        ##dgomol=DihGOpt_Molecule(fifties_tuple=self.fifties_tuple,pairs_tuple=self.pairs_tuple,dih_tuple=self._dt)
        #dgomol=DihGOpt_Molecule()
        #dgomol.dt=(1,0,4,13)
        #fn_b="opmmm-mp2-opt"
        #dgomol.read_gopt_log(fn_b)
        #b,e,s=self._da_rt
        #e+=1
        #for n in range(b,e,s):
        #        sn0=str(n)
        #        if (n<10):
        #                sn="00"+sn0
        #        elif (n<100):
        #                sn="0"+sn0
        #        else:
        #                sn=sn0
	#	dgomol.write_inp_pcm(self.gopt_scan_fnameb+sn)
        for m in self._ml:
            fnameb=m.name
            m.write_inp_pcm(fnameb)

    def read_engine_outputs(self):
        self._eed={}
        for m in self._ml:
            fnameb=m.name
            m.read_out_pcm(fnameb)
            #self._dad.update({m.name:rad2deg(m.da)})
            self._eed.update({m.name:m.e})

    def run_dih_elem(self,m):
            fnameb=m.name
            coengine_name="coengine_"+fnameb
            f=open("../Data/Engine/"+coengine_name,'w')
            print >>f,"mode opt"
            print >>f,"infile  "+fnameb+"_inp.pcm"
            print >>f,"outfile "+fnameb+"_out.pcm"
            print >>f,"print 3"
            if self._mm=="mm3":
               print >>f,"forcefield mm3"
               print >>f,"addpar add_MM3.prm"
            elif self._mm=="mmff94":
               print >>f,"forcefield mmff94"
               print >>f,"addpar add_MMFF94.prm"
            else:
               print "DihScan.run_dih_scan: Wrong MM-name!"
            f.close()
            #
            timeout=1
            while True:
                status=run_engine_timeout(self.engine_path,coengine_name,timeout)
                if status:
                    break
                else:
                    p[0]+=0.001
                    #pert_add_param("add_MM3.prm")
                    write_add(p,c,mm,ol_templ,lines,1,None,None)
            #
            comm="rm ../Data/Engine/"+coengine_name
            system(comm)

    def run_dih_scan(self,p,c,mm,ol_templ,lines):
        from multiprocessing import Pool
        p=Pool()
        p.map(self.run_dih_elem,self._ml)

    def calc_rmse(self,csv,step,step_int):
        self.read_engine_outputs()
        m0=self._ml[0]
        n0=m0.name
        ge0,ee0=self._ged[n0],self._eed[n0]
        rmse=0.
        if csv=="csv_on" and (step-1)%step_int==0:
            f=open("../Data/ParFit/opt_"+str(step)+".csv",'w')
        for m in self._ml:
            name=m.name
            ged,eed=(self._ged[name]-ge0)*self.hartree2kcal_mol,self._eed[name]-ee0
            rmse+=(ged-eed)**2
            if csv=="csv_on" and (step-1)%step_int==0:
                print >>f,"%f,%f,%f"%(self._dad[name],ged,eed)
        if csv=="csv_on" and (step-1)%step_int==0:
            f.close()
        nm=len(self._ml)
        rmse=sqrt(rmse/nm)
        return rmse

if __name__=="__main__":
    ##ds.write_gamess_inputs()
    #ds.read_gamess_outputs()
    #ds.write_pcm()

    #print ds.gopt_base_fnameb
    #print ds.gopt_scan_fnameb
    #print ds.fifties_tuple
    #print ds.pairs_tuple
    #print ds.da_rt
    #print ds.dt

    #ds.read_gamess_outputs()
    #ds.write_gouts_data()
    #ds.read_gouts_data()
    #ds.write_engine_inputs()

    #engine_rmse(p)
    default_input_fname="dih_scan_inp"

    gopt_type,gopt_s_fnameb,t1234,bes,engine_path,mm,opt_lin,np,csv=par_fit_inp(default_input_fname)
    ds=DihScan(gopt_s_fnameb,engine_path,mm,bes,t1234)
    if not gopt_type=="ginp":
       environ["Engine_DIR"]=engine_path+"engine_dir"
       p,c,ol_templ,lines=read_add(mm,opt_lin,np,1)

    def engine_rmse(p):
       print p
       step=1
       step_int=10
       write_add(p,c,mm,ol_templ,lines,1,step,step_int)
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
       "Par_Fit: Wrong gopt_type!"
    ds.write_engine_inputs()
    engine_rmse(p)
