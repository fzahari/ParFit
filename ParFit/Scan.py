#!/usr/bin/env python

from os import system,environ
from numpy import rad2deg,array,sqrt,pi,around,deg2rad
from _IO import par_fit_inp,read_add,write_add
from _Engine import run_engine_timeout, pert_add_param
from GeomStr import Molecule,default_mm3_type,default_mmff94_type

import os

import copy_reg
import types

def _reduce_method(meth):
    return (getattr,(meth.__self__,meth.__func__.__name__))
copy_reg.pickle(types.MethodType,_reduce_method)

class ScanElem(Molecule):
    """\
    A class to handle Gamess and Engine input/output 
    on an individual element of a scan. 
    """
    input_pcm_head1="{PCM\nNA "
    input_pcm_head2="\nATOMTYPES "

    def __init__(self,a_list=(),tup=(),mm=None,t_list=[],name="Molecule"):
        super(ScanElem,self).__init__(a_list,mm,t_list,name)
        self._t=tuple(tup)
        self._v=None
        self._e=None
        self._ginp_templ1=None
        self._ginp_templ2=None

    @property
    def t(self):
        return self._t
    
    @t.setter
    def t(self,tup):
        self._t=tuple(tup)

    @property
    def v(self):
        if self._na>0:
            assert not self._t==()
            if len(self._t)==4:
               t1,t2,t3,t4=self._t
               self._v=self.calc_dihedral(t1,t2,t3,t4)
            elif len(self._t)==2:
               t1,t2=self._t
               self._v=self.calc_dist(t1,t2)
            elif len(self._t)==3:
               t1,t2,t3=self._t
               self._v=self.calc_angle(t1,t2,t3)
        return self._v

    @property
    def e(self):
        return self._e
 
    def read_gopt_log(self,fname_base):	
        gkey="EQUILIBRIUM GEOMETRY LOCATED"
        ekey="TOTAL ENERGY      ="
        #
        #fname="../Data/Gamess/"+fname_base+"-dlc.log"
        fname="../Data/Gamess/"+fname_base+".log"
        f=open(fname,'r')
        lines=f.readlines()
        f.close()
        #
        ln=len(lines)
        for i in range(ln):
            if gkey in lines[i]: break
        #
        #self._rl=[]
        #atom_list=[]
        for line in lines[i+4:]:
            if line[:-1]=="": break
            s,c,x,y,z=line[:-1].split()
            s=s.upper()
            self._sl.append(s)
            self._rl.append(array(map(float,[x,y,z]),'d'))
            self._cl.append(float(c))
            if self._mm=="mm3":
                self._tl.append(default_mm3_type[s])
            elif self._mm=="mmff94":
                self._tl.append(default_mmff94_type[s])
            else:
                print "ScanElem.read_gopt_log: Wrong MM-type!"
            #s,c,x,y,z=line[:-1].split()
            #symbl=s.upper()
            #coord=map(float,[x,y,z])
            #charg=float(c)
            #atom_list.append(Atom(symbl,coord,charg,self._mm))
        #self.a_list_update(atom_list)
        self._rl=array(self._rl,'d')
        self._na=len(self._rl)
        #
        for line in lines[i+5+self._na:]:
            if ekey in line: self._e=float(line.split()[3])
        #
        self.set_conn()
        self.bens_ring()
        self.bond_ord()

    def read_ginp(self,fname_base):	
        dkey="$DATA"
        ekey="$END"
        #
        fname="../Data/Gamess/"+fname_base+".inp"
        f=open(fname,'r')
        lines=f.readlines()
        f.close()
        #
        ln=len(lines)
        for i in range(ln):
            if dkey in lines[i].upper(): break
        #
        self._ginp_templ1=lines[:i+3]
        #
        #atom_list=[]
        for j in range(i+3,ln):
            if ekey in lines[j].upper(): break
            s,c,x,y,z=lines[j][:-1].split()
            s=s.upper()
            self._sl.append(s)
            self._rl.append(array(map(float,[x,y,z]),'d'))
            self._cl.append(float(c))
            if self._mm=="mm3":
                self._tl.append(default_mm3_type[s])
            elif self._mm=="mmff94":
                self._tl.append(default_mmff94_type[s])
            else:
                print "ScanElem.read_gopt_log: Wrong MM-type!"
            #s,c,x,y,z=lines[j][:-1].split()
            #symbl=s.upper()
            #coord=map(float,[x,y,z])
            #charg=float(c)
            #atom_list.append(Atom(symbl,coord,charg,self._mm))
        #self.a_list_update(atom_list)
        self._rl=array(self._rl,'d')
        self._na=len(self._rl)
        self._ginp_templ2=lines[j+1:]
        #
        #
        self.set_conn()
        self.bens_ring()
        self.bond_ord()

    def write_gopt_inp(self,fname_base):
        pass

    def write_inp_pcm(self,fname_base,styp):
        f=open("../Data/Engine/"+fname_base+"_inp.pcm",'w')
        if self._mm=="mm3":
            self.input_pcm_head2 += str(3)
        elif self._mm=="mmff94":
            self.input_pcm_head2 += str(7)
        else:
            print "ScanElem.write_inp_pcm: Wrong MM-type!"
        print >>f,self.input_pcm_head1+str(self._na)+self.input_pcm_head2
        for i in range(self._na):
            iconn=self._conn[i]
            if i in self._br:
                if self._mm=="mm3":
                    t=50
                elif self._mm=="mmff94":
                    t=37
                else:
                    print "ScanElem.write_inp_pcm: Wrong MM-type!"
            else:
                t=self.tl[i]
                if self._mm=="mm3":
                    if t==1 and len(iconn)==3: t=2
                    if t==7 and len(iconn)==2: t=6
                elif self._mm=="mmff94":
                    if t==1 and len(iconn)==3: t=2
                    if t==32 and len(iconn)==2: t=6
                else:
                    print "ScanElem.write_inp_pcm: Wrong MM-type!"
            x,y,z=self.rl[i]
            print >>f,"AT %d %d %10.4f %10.4f %10.4f B"%((i+1),t,x,y,z),
            for j in iconn:
                if (i,j) in self._db or (j,i) in self._db:
                    print >>f,j+1,"2",
                else:
                    print >>f,j+1,"1",
            print >>f
        if styp=="bond":
           print >>f,"}"
        elif styp=="angl":
           print >>f,"}"
        elif styp=="diha":
           t1,t2,t3,t4=self._t
           self._v=self.calc_dihedral(t1,t2,t3,t4)
           da=rad2deg(self._v)
           st1,st2,st3,st4=map(str,map(lambda x:x+1,self._t))
           tail="DD "+st1+" "+st2+" "+st3+" "+st4+" FROM %.1f TO %.1f BY 0.0\n}"%(da,da)
           print >>f,tail
        f.close() 

    def read_out_pcm(self,fname_base,styp):
        f=open("../Data/Engine/"+fname_base+"_out.pcm",'r')
        lines=f.readlines()
        f.close()
        self._na=int(lines[1].split()[1])
        self._rl=[]
        for line in lines[5:5+self._na]:
            self._rl.append(array(map(float,line.split()[3:6])))
        if styp=="bond":
           self._e=float(lines[2].split()[6])
        elif styp=="angl":
           self._e=float(lines[2].split()[6])
        elif styp=="diha":
           self._e=float(lines[2].split()[8])

    def diha_rot(self,diha_rot_val):
        if not diha_rot_val is None:
            super(ScanElem,self).diha_rot(self._t,diha_rot_val)

    def bond_tra(self,bond_tra_val):
        if not bond_tra_val is None:
            super(ScanElem,self).bond_tra(self._t,bond_tra_val)

    def angl_rot(self,angl_rot_val):
        if not angl_rot_val is None:
            super(ScanElem,self).angl_rot(self._t,angl_rot_val)

class Scan(object):
    """\
    A base class for doing a scan. 
    """
    hartree2kcal_mol=627.509469
    def __init__(self,sdir,gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,ran_tup,tup=(),gopt_b_fnameb="mp2_base",styp=None):
        self.sdir=sdir
        self.gopt_base_fnameb=gopt_b_fnameb
        self.engine_path=engine_path
        self.gopt_scan_fnameb=gopt_s_fnameb
        self._mm=mm
        self._opt_lin=opt_lin
        self._np=np
        self._nc=nc
        if not ran_tup==():
           assert len(ran_tup)==3
        self._rt=ran_tup
        #if not tup==():
        #   assert len(tup)==4
        if not tup==():
           self._t=tuple(map(lambda x:x-1,tup))
        else: 
           self._t=()
        self._ml=[]
        self._v={}
        self._ge={}
        self._ee={}
        self._styp=styp

    @property
    def mm(self):
        return self._mm

    @property
    def rt(self):
        return self._rt

    @property
    def t(self):
        return self._t

    @property
    def ml(self):
        return self._ml

    @property
    def v(self):
        return self._v

    @property
    def ge(self):
        return self._ge

    @property
    def ee(self):
        return self._ee

    def _gen_geometries(self):
        pass

    def read_gamess_outputs(self):
        self._ge={}
        b,e,s=self._rt
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
            se=ScanElem(tup=self._t,mm=self._mm,name=fnameb)
            se.read_gopt_log(fnameb)
            self._ml.append(se)
            self._v.update({se.name:rad2deg(se.v)})
            self._ge.update({se.name:se.e})

    def read_gouts_data(self):
        fnameb=self.gopt_scan_fnameb
        f=open("../Data/Gamess/"+fnameb+"scan",'r')
        lines=f.readlines() 
        f.close()
        self._t=map(int,lines[0][:-1].split())
        #self._t=map(lambda x:x-1,self._t)
        self._rt=map(int,map(float,lines[1][:-1].split()))
        l=int(lines[2][:-1])
        lines=lines[3:]
        self._v={}
        self._ged={}
        self._ml=[]
        for i in range(l):
            fnameb,na,v,e=lines[0][:-1].split()
            na=int(na)
            se=ScanElem(tup=self._t,mm=self._mm,name=fnameb)
            se._sl=[]
            se._rl=[]
            se._cl=[]
            se._tl=[]
            lines=lines[1:]
            #atom_list=[]
            for j in range(na):
                s,x,y,z,c=lines[j].split()
                s=s.upper()
                se._sl.append(s)
                se._rl.append(array(map(float,[x,y,z]),'d'))
                se._cl.append(float(c))
                if self._mm=="mm3":
                    se._tl.append(default_mm3_type[s])
                elif self._mm=="mmff94":
                    se._tl.append(default_mmff94_type[s])
                else:
                    print "ScanElem.read_gopt_log: Wrong MM-type!"
                #s,x,y,z,c=lines[j].split()
                #symbl=s.upper()
                #coord=map(float,[x,y,z])
                #charg=float(c)
                #atom_list.append(Atom(symbl,coord,charg,self._mm))
            #se.a_list_update(atom_list)
            se._rl=array(se._rl,'d')
            se._na=len(se._rl)
            se.set_conn()
            se.bens_ring()
            se.bond_ord()
            self._ml.append(se)
            v,e=map(float,[v,e])
            self._v.update({fnameb:v})
            self._ge.update({fnameb:e})
            lines=lines[na:]

    def write_engine_inputs(self):
        for m in self._ml:
            fnameb=m.name
            m.write_inp_pcm(fnameb,self._styp)

    def read_engine_outputs(self):
        self._ee={}
        for m in self._ml:
            fnameb=m.name
            m.read_out_pcm(fnameb,self._styp)
            #self._v.update({m.name:rad2deg(m.v)})
            self._ee.update({m.name:m.e})

    def run_dih_elem(self,m):
            fnameb=m.name
            coengine_name="coengine_"+fnameb
            #f=open("../Data/Engine/"+coengine_name,'w')
            f=open(coengine_name,'w')
            if self._styp=="diha":
               print >>f,"mode opt"
            else:
               print >>f,"mode single"
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
               print "Scan.run_scan: Wrong MM-name!"
            f.close()
            #
            timeout=1
            while True:
                status=run_engine_timeout(self.engine_path,coengine_name,timeout)
                if status:
                    break
                else:
                    p,c,ol_templ,lines=read_add(self._mm,self._opt_lin,self._np,self._nc,1,self._styp)
                    p[0]+=0.001
                    #pert_add_param("add_MM3.prm")
                    write_add(self.sdir,p,c,self._mm,ol_templ,lines,1,None,None)
            #
            #comm="rm ../Data/Engine/"+coengine_name
            comm="rm "+coengine_name
            system(comm)

    def run_scan(self,p,c,mm,ol_templ):
        os.chdir("../Data/Engine")
        from multiprocessing import Pool
        p=Pool()
        p.map(self.run_dih_elem,self._ml)
        p.terminate()
        #p.close()
        #p.join()
        #map(self.run_dih_elem,self._ml)
        os.chdir("../../ParFit")

    def calc_rmse(self,csv,mi,step,step_int):
        self.read_engine_outputs()
        m0=self._ml[0]
        n0=m0.name
        ge0,ee0=self._ge[n0],self._ee[n0]
        rmse=0.
        if csv=="csv_on" and (step-1)%step_int==0:
            f=open(self.sdir[mi+1]+"/opt"+"_"+str(step)+".csv",'w')
        for m in self._ml:
            name=m.name
            ge,ee=(self._ge[name]-ge0)*self.hartree2kcal_mol,self._ee[name]-ee0
            rmse+=(ge-ee)**2
            if csv=="csv_on" and (step-1)%step_int==0:
                print >>f,"%f,%f,%f"%(self._v[name],ge,ee)
        if csv=="csv_on" and (step-1)%step_int==0:
            f.close()
        nm=len(self._ml)
        rmse=sqrt(rmse/nm)
        return rmse

class BondScan(Scan):
    def __init__(self,sdir,gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,ran_tup,tup=(),gopt_b_fnameb="mp2_base"):
       if not tup==():
          assert len(tup)==2
       super(BondScan,self).__init__(sdir,gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,ran_tup,tup,gopt_b_fnameb,"bond")

    def write_gamess_inputs(self):
        #self._gen_geometries()
        fnameb=self.gopt_scan_fnameb
        rt=self.rt
        se=ScanElem(tup=self.t,mm=self.mm,name=fnameb)
        se.read_ginp(fnameb)
        t1,t2=self.t
        dar=se.calc_dist(t1,t2)
        b,e,s=rt
        se.bond_tra(0.1*b-dar)
        e+=1
        na=se.na
        for n in range(b,e,s):
           sn0=str(n)
           if (n<10):
              sn="00"+sn0
           elif (n<100):
              sn="0"+sn0
           else:
              sn=sn0
           f=open("../Data/Gamess/"+fnameb+"-"+sn+".inp",'w')
           #f=open("../Data/Gamess/"+fnameb+sn+".inp",'w')
           print >>f," $ZMAT DLC=.T. AUTO=.T. $END"
           print >>f," $ZMAT IFZMAT(1)=1,",t1+1,",",t2+1," FVALUE(1)=",0.1*float(n),"$END"
           print >>f," $CONTRL COORD=UNIQUE NZVAR=",3*na-6,"$END"
           for line in se._ginp_templ1:
              print >>f,line[:-1]
           if not n==b:
              se.bond_tra(0.1*s)
           rl=se.rl
           cl=se.cl
           sl=se.sl
           for i in range(na):
              x,y,z=rl[i]
              print >>f,"",sl[i],cl[i],x,y,z
           print >>f," $END"
           for line in se._ginp_templ2:
              print >>f,line[:-1]
           f.close()

    def write_gouts_data(self):
        fnameb=self.gopt_scan_fnameb
        f=open("../Data/Gamess/"+fnameb+"scan",'w')
        t1,t2=map(str,map(lambda x:x+1,self._t))
        print >>f,t1,t2
        b,e,s=map(lambda x:0.1*x,map(float,self._rt))
        print >>f,b,e,s
        l=len(self._ml)
        print >>f,l
        for i in range(l):
            m=self._ml[i]
            n=m._n
            na=m._na
            print >>f,n,na,self._v[n],self._ge[n]
            for j in range(na):
                x,y,z=m.rl[j]
                print >>f,m.sl[j],x,y,z,m.cl[j]
        f.close()

class AnglScan(Scan):
    def __init__(self,sdir,gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,ran_tup,tup=(),gopt_b_fnameb="mp2_base"):
       if not tup==():
          assert len(tup)==3
       super(AnglScan,self).__init__(sdir,gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,ran_tup,tup,gopt_b_fnameb,"angl")

    def write_gamess_inputs(self):
        #self._gen_geometries()
        fnameb=self.gopt_scan_fnameb
        rt=self.rt
        se=ScanElem(tup=self.t,mm=self.mm,name=fnameb)
        se.read_ginp(fnameb)
        t1,t2,t3=self.t
        dar=se.calc_angle(t1,t2,t3)
        b,e,s=rt
        se.angl_rot(deg2rad(0.1*b)-dar)
        e+=1
        na=se.na
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
           print >>f," $ZMAT IFZMAT(1)=2,",t1+1,",",t2+1,",",t3+1," FVALUE(1)=",0.1*float(n),"$END"
           print >>f," $CONTRL COORD=UNIQUE NZVAR=",3*na-6,"$END"
           for line in se._ginp_templ1:
              print >>f,line[:-1]
           if not n==b:
              se.angl_rot(deg2rad(0.1*s))
           rl=se.rl
           cl=se.cl
           sl=se.sl
           for i in range(na):
              x,y,z=rl[i]
              print >>f,"",sl[i],cl[i],x,y,z
           print >>f," $END"
           for line in se._ginp_templ2:
              print >>f,line[:-1]
           f.close()

    def write_gouts_data(self):
        fnameb=self.gopt_scan_fnameb
        f=open("../Data/Gamess/"+fnameb+"scan",'w')
        t1,t2,t3=map(str,map(lambda x:x+1,self._t))
        print >>f,t1,t2,t3
        b,e,s=map(lambda x:0.1*x,map(float,self._rt))
        print >>f,b,e,s
        l=len(self._ml)
        print >>f,l
        for i in range(l):
            m=self._ml[i]
            n=m._n
            na=m._na
            print >>f,n,na,self._v[n],self._ge[n]
            for j in range(na):
                x,y,z=m.rl[j]
                print >>f,m.sl[j],x,y,z,m.cl[j]
        f.close()

class DihAScan(Scan):
    def __init__(self,sdir,gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,ran_tup,tup=(),gopt_b_fnameb="mp2_base"):
       if not tup==():
          assert len(tup)==4
       super(DihAScan,self).__init__(sdir,gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,ran_tup,tup,gopt_b_fnameb,"diha")

    def write_gamess_inputs(self):
        #self._gen_geometries()
        fnameb=self.gopt_scan_fnameb
        rt=self.rt
        se=ScanElem(tup=self.t,mm=self.mm,name=fnameb)
        se.read_ginp(fnameb)
        t1,t2,t3,t4=self.t
        dar=se.calc_dihedral(t1,t2,t3,t4)
        b,e,s=rt
        se.diha_rot(deg2rad(b)-dar)
        e+=1
        na=se.na
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
           for line in se._ginp_templ1:
              print >>f,line[:-1]
           if not n==b:
              se.diha_rot(deg2rad(s))
           rl=se.rl
           cl=se.cl
           sl=se.sl
           for i in range(na):
              x,y,z=rl[i]
              print >>f,"",sl[i],cl[i],x,y,z
           print >>f," $END"
           for line in se._ginp_templ2:
              print >>f,line[:-1]
           f.close()

    def write_gouts_data(self):
        fnameb=self.gopt_scan_fnameb
        f=open("../Data/Gamess/"+fnameb+"scan",'w')
        t1,t2,t3,t4=map(str,map(lambda x:x+1,self._t))
        print >>f,t1,t2,t3,t4
        b,e,s=map(float,self._rt)
        print >>f,b,e,s
        l=len(self._ml)
        print >>f,l
        for i in range(l):
            m=self._ml[i]
            n=m._n
            na=m._na
            print >>f,n,na,self._v[n],self._ge[n]
            for j in range(na):
                x,y,z=m.rl[j]
                print >>f,m.sl[j],x,y,z,m.cl[j]
        f.close()

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

    gopt_type,gopt_s_fnameb,t1234,bes,engine_path,mm,mode,alg,opt_lin,np,nc,step_int,csv=par_fit_inp(default_input_fname)
    ds=DihScan(gopt_s_fnameb,engine_path,mm,opt_lin,np,nc,bes,t1234)
    if not gopt_type=="ginp":
       environ["ENGINE_DIR"]=engine_path+"engine_dir"
       p,c,ol_templ,lines=read_add(mm,opt_lin,np,nc,1)
  
    def engine_rmse(p):
       print p
       step=1
       step_int=10
       write_add(p,c,mm,ol_templ,lines,1,step,step_int)
       ds.run_scan(p,c,mm,ol_templ,lines)
       rmse=ds.calc_rmse(csv,step,step_int)
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

    gmol=ScanElem(mm="mmff94",name="opmmm")    
    print gmol.name
    #gmol.dt=(13,4,0,1)
    gmol.dt=(1,0,4,13)
    print gmol.dt
    fn_b="opmmm-mp2-opt"
    gmol.read_gopt_log(fn_b)
    print "tl:",gmol.tl
    print gmol.e
    print gmol.sl
    print gmol.rl
    print gmol.cl
    print gmol.da
    gmol.dih_rot(pi/2.0)
    print gmol.rl
    print gmol.da
    print gmol.write_inp_pcm("test")
    gmol=ScanElem(mm="mmff94")
    gmol.read_ginp("tempC2")    
    gmol.dt=(13,4,0,1)
    print gmol._ginp_templ
    print gmol._rl
    gmol.dih_rot(pi/2.0)
    print gmol._rl
