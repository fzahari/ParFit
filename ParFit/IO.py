#!/usr/bin/env python

from numpy import zeros,array,rad2deg,pi
from GeomStr import Molecule,default_mm3_type,default_mmff94_type

def ginp_inp(input_fname):
    f=open(input_fname,'r')
    lines=f.readlines()
    f.close()
    line1=lines[0][:-1].split(',')
    gopt_scan_fnameb,t1234,bes=line1
    gopt_scan_fnameb=gopt_scan_fnameb.strip()
    t1,t2,t3,t4=map(int,t1234.split())
    t1234=(t1,t2,t3,t4)
    b,e,s=map(int,bes.split())
    bes=(b,e,s)
    return gopt_scan_fnameb,t1234,bes

def par_fit_inp(input_fname):
    f=open(input_fname,'r')
    lines=f.readlines()
    f.close()
    line1=lines[0][:-1].split(',')
    if line1[0]=="mult": 
       n=int(line1[1])
       gopt_type=[]
       gopt_scan_fnameb=[]
       t1234=[]
       bes=[]
       for i in range(n):
          sgopt_type,sgopt_scan_fnameb,st1234,sbes=lines[i+1].split(',')
          gopt_type.append(sgopt_type.strip())
          gopt_scan_fnameb.append(sgopt_scan_fnameb.strip())
          t1,t2,t3,t4=map(int,st1234.split())    
          t1234.append((t1,t2,t3,t4))
          b,e,s=map(int,sbes.split())    
          bes.append((b,e,s))
    else:
       n=0
       gopt_type,gopt_scan_fnameb,t1234,bes=line1
       gopt_type=[gopt_type.strip(),]
       gopt_scan_fnameb=[gopt_scan_fnameb.strip(),]
       t1,t2,t3,t4=map(int,t1234.split())
       t1234=[(t1,t2,t3,t4),]
       b,e,s=map(int,bes.split())
       bes=[(b,e,s),]
    if gopt_type[0]=="ginp":
       return gopt_type,gopt_scan_fnameb,t1234,bes,None,None,None,None,None,None,None,None,None
    engine_path=lines[n+1][:-1]
    l2s=lines[n+2].split()
    mm=l2s[0].lower()
    if len(l2s)==2:
       mode=l2s[1].lower()
    else:
       mode="opt"
    np=0
    #plist=[]
    nc=0
    opt_lin={}
    alg=lines[n+3].strip()
    for line in lines[n+4:-1]:
	t=line.split()
        for i in range(3):
            t[i+1]=t[i+1].lower()
            if t[i+1][0]=='p': 
               #plab=t[i+1][1:]
               #if len(plab)>0:               
               #   if plab in plist:
               #      continue
               #   else:
               #      plist.append(plab)
               np+=1
            if t[i+1]=='c': 
               nc+=1
        opt_lin.update({int(t[0])-1:(t[1],t[2],t[3])})
    last_line=lines[-1].split()
    if len(last_line)==1:
       step_int=10
       csv=lines[-1].strip()
    elif len(last_line)==2:
       step_int,csv=lines[-1].split()
    else:
       print "Wrong csv line in the input file!"

    return gopt_type,gopt_scan_fnameb,t1234,bes,engine_path,mm,mode,alg,opt_lin,np,nc,step_int,csv

def read_add(mm,opt_lin,np,nc,fl):

    p=zeros(np)
    c=zeros(nc)
     
    if mm=="mm3":
        if fl==1:
           add_name="../Data/Engine/add_MM3.prm"
        else:
           add_name="../Data/Engine/add_MM3_PF.prm"
    elif mm=="mmff94":
        if fl==1:
           add_name="../Data/Engine/add_MMFF94.prm"
        else:
           add_name="../Data/Engine/add_MMFF94_PF.prm"
    else:
        print "read_add: Wrong MM-type!"

    f=open(add_name,'r')
    lines=f.readlines()
    f.close()
        
    i=j=0
    ol_keys=opt_lin.keys()
    ol_keys.sort()
    pdict={}
    ol_templ={}
    for ol_sk in ol_keys:
        s=""
        t=opt_lin[ol_sk]
        for k in range(3):
            if t[k][0]=='p':
                plab=t[k][1:]
                if len(plab)>0:
                   if pdict.has_key(plab):
                      ii=pdict[plab]
                   else:
                      pdict.update({plab:i})
                      ii=i
                      i+=1
                else:
                   ii=i
                   i+=1
                s+="p["+str(ii)+"],"
                #s+="p["+str(i)+"],"
                #i+=1
            elif t[k]=='c':
                s+="c["+str(j)+"],"
                j+=1
            else:
                print "Woops!" 
        v=lines[ol_sk].split()
        ol_vars=s[:-1]
        exec ol_vars+"="+v[5]+","+v[7]+","+v[9]
        v[1:5]=map(int,v[1:5])
        s1='"%s        %3i  %3i  %3i  %3i'%(v[0],v[1],v[2],v[3],v[4])
        s2='      %6.3f +1   %6.3f -2 %6.3f +3      "%('+ol_vars+')'
        ol_templ.update({ol_sk:s1+s2})

    return p,c,ol_templ,lines

def write_add(sdir,p,c,mm,ol_templ,lines,fl,step,step_int):
    if mm=="mm3":
       if fl==1:
          add_name="../Data/Engine/add_MM3.prm"
       else:
          add_name="../Data/Engine/add_MM3_PF.prm"
    elif mm=="mmff94":
       if fl==1:
          add_name="../Data/Engine/add_MMFF94.prm"
       else:
          add_name="../Data/Engine/add_MMFF94_PF.prm"
    else:
       print "write_add: Wrong MM-type!"

    f=open(add_name,'w')

    ol_keys=ol_templ.keys()
    ol_keys.sort()
    for ol_sk in ol_keys:
        exec 'lines['+str(ol_sk)+']='+ol_templ[ol_sk]

    for line in lines:
       print >>f,line[:-1]

    f.close()
   
    if step==None: return
    
    if step=="ga" or (step-1)%step_int==0:
       if mm=="mm3":
          add_name_arch=sdir[0]+"/add_MM3_"+str(step)+".prm"
       elif mm=="mmff94":
          add_name_arch=sdir[0]+"/add_MMFF94_"+str(step)+".prm"

       f=open(add_name_arch,'w')
       for line in lines:
          print >>f,line[:-1]
       f.close()

    return
    
class DihGOpt_Molecule(Molecule):
    """\
    A class to handle Gamess and Engine input/output. 
    For now, it works on (constraint dihedral angle)
    geometry optimization cases only.
    """
    input_pcm_head1="{PCM\nNA "
    input_pcm_head2="\nATOMTYPES "

    def __init__(self,a_list=(),dih_tuple=(),mm=None,t_list=[],name="Molecule"):
        super(DihGOpt_Molecule,self).__init__(a_list,mm,t_list,name)
        if not tuple(dih_tuple)==():
            assert len(dih_tuple)==4
        self._dt=tuple(dih_tuple)
        self._da=None
        self._e=None
        #self._ginp_templ=None

    @property
    def dt(self):
        return self._dt
    
    @dt.setter
    def dt(self,dih_tuple):
        assert len(dih_tuple)==4
        self._dt=tuple(dih_tuple)

    @property
    def da(self):
        if self._na>0:
            assert not self._dt==()
            t1,t2,t3,t4=self._dt
            print "da property1:",t1,t2,t3,t4
            self._da=self.calc_dihedral(t1,t2,t3,t4)
            print "da property2:",self._da
        return self._da

    @property
    def e(self):
        return self._e
 
    def read_gopt_log(self,fname_base):	
        gkey="EQUILIBRIUM GEOMETRY LOCATED"
        ekey="TOTAL ENERGY      ="
        #
        fname="../Data/Gamess/"+fname_base+"-dlc.log"
        f=open(fname,'r')
        lines=f.readlines()
        f.close()
        #
        ln=len(lines)
        for i in range(ln):
            if gkey in lines[i]: break
        #
        self._rl=[]
        for line in lines[i+4:]:
            if line[:-1]=="": break
            s,c,x,y,z=line[:-1].split()
            self._sl.append(s)
            self._rl.append(array(map(float,[x,y,z]),'d'))
            self._cl.append(float(c))
            if self._mm=="mm3":
                self._tl.append(default_mm3_type[s])
            elif self._mm=="mmff94":
                self._tl.append(default_mmff94_type[s])
            #else:
            #    print "DihGOpt_Molecule.read_gopt_log: Wrong MM-type!"
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
        self._ginp_templ=lines[:i+3]
        #
        for line in lines[i+3:]:
            if ekey in line.upper(): break
            s,c,x,y,z=line[:-1].split()
            self._sl.append(s)
            self._rl.append(array(map(float,[x,y,z]),'d'))
            self._cl.append(float(c))
            if self._mm=="mm3":
                self._tl.append(default_mm3_type[s])
            elif self._mm=="mmff94":
                self._tl.append(default_mmff94_type[s])
            #else:
            #    print "DihGOpt_Molecule.read_gopt_log: Wrong MM-type!"
        self._rl=array(self._rl,'d')
        self._na=len(self._rl)
        #
        #
        self.set_conn()
        self.bens_ring()
        self.bond_ord()

    def write_gopt_inp(self,fname_base):
        pass

    def write_inp_pcm(self,fname_base):
        f=open("../Data/Engine/"+fname_base+"_inp.pcm",'w')
        if self._mm=="mm3":
            self.input_pcm_head2=self.input_pcm_head2+str(3)
        elif self._mm=="mmff94":
            self.input_pcm_head2=self.input_pcm_head2+str(7)
        else:
            print "DihGOpt_Molecule.write_inp_pcm: Wrong MM-type!"
        print >>f,self.input_pcm_head1+str(self._na)+self.input_pcm_head2
        for i in range(self._na):
            iconn=self._conn[i]
            if i in self._br:
                if self._mm=="mm3":
                    t=50
                elif self._mm=="mmff94":
                    t=37
                else:
                    print "DihGOpt_Molecule.write_inp_pcm: Wrong MM-type!"
            else:
                t=self.tl[i]
                if self._mm=="mm3":
                    if t==1 and len(iconn)==3: t=2
                    if t==7 and len(iconn)==2: t=6
                elif self._mm=="mmff94":
                    if t==1 and len(iconn)==3: t=2
                    if t==32 and len(iconn)==2: t=6
                else:
                    print "DihGOpt_Molecule.write_inp_pcm: Wrong MM-type!"
            x,y,z=self.rl[i]
            print >>f,"AT %d %d %10.4f %10.4f %10.4f B"%((i+1),t,x,y,z),
            for j in iconn:
                if set((i,j)) in self._db:
                    print >>f,j+1,"2",
                else:
                    print >>f,j+1,"1",
            print >>f
        t1,t2,t3,t4=self._dt
        self._da=self.calc_dihedral(t1,t2,t3,t4)
        da=rad2deg(self._da)
        st1,st2,st3,st4=map(str,map(lambda x:x+1,self._dt))
        tail="DD "+st1+" "+st2+" "+st3+" "+st4+" FROM %.1f TO %.1f BY 0.0\n}"%(da,da)
        print >>f,tail
        f.close() 

    def read_out_pcm(self,fname_base):
        f=open("../Data/Engine/"+fname_base+"_out.pcm",'r')
        lines=f.readlines()
        f.close()
        self._na=int(lines[1].split()[1])
        self._rl=[]
        for line in lines[5:5+self._na]:
            self._rl.append(array(map(float,line.split()[3:6])))
        self._e=float(lines[2].split()[8])
 
    def dih_rot(self,dih_angle_rot):
        if not dih_angle_rot==None:
            super(DihGOpt_Molecule,self).dih_rot(self._dt,dih_angle_rot)

if __name__=="__main__":
    gmol=DihGOpt_Molecule(mm="mmff94",name="opmmm")    
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
    gmol=DihGOpt_Molecule(mm="mmff94")
    gmol.read_ginp("tempC2")    
    gmol.dt=(13,4,0,1)
    print gmol._ginp_templ
    print gmol._rl
    gmol.dih_rot(pi/2.0)
    print gmol._rl
