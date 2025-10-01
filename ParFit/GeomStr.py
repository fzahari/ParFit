#!/usr/bin/env python3

from numpy import array,pi
from ._GeomCalc import dist,angle,dihedral,norm_vec,uangl,rotu,tra,angle

keywords=["default_charge", "default_mm3_type", "default_mmff94_type", "cov_radii", "bond_ords"]

default_charge={}
default_mm3_type={}
default_mmff94_type={}
cov_radii={}
bond_ords={}

import os
f=open(os.path.join(os.path.dirname(__file__), "atomic.db"),'r')
lines=f.readlines()
f.close()

for line in lines:
   l=line[:-1]
   if l in keywords:
      l0=l
      exec(l0+"={}")
   ls=l.split()
   if len(ls)==2:
      a,b=ls
      if l0=="bond_ords":
         b=list(map(float,b.split(',')))
      else:
         b=float(b)
      exec(l0+"[a]=b")

class Atom(object):
    def __init__(self,symbol,r,charge=None,mm=None,mm_type=None):
        self._s=symbol
        assert len(r)==3
        self._r=array(r,'d')
        if charge is None:
            self._c=default_charge[self._s]
        else:
            self._c=charge
        self._mm=mm
        if self._mm=="mm3":
            if mm_type is None:
                self._t=default_mm3_type[self._s]
            else:
                self._t=mm_type
        elif self._mm=="mmff94":
            if mm_type is None:
                self._t=default_mmff94_type[self._s]
            else:
                self._t=mm_type
        else:
            print("Atom.__init__: Wrong MM-name!")

    @property
    def s(self):
        return self._s

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self,r):
        assert len(r)==3
        self._r=array(r,'d')

    @property
    def c(self):
        return self._c

    @c.setter
    def c(self,charge):
        self._c=charge

    @property
    def mm(self):
        return self._mm

    @mm.setter
    def mm(self,mm):
        self._mm=mm
        if self._mm=="mm3":
            self._t=default_mm3_type[self._s]
        elif self._mm=="mmff94":
            self._t=default_mmff94_type[self._s]
        else:
            print("Atom.mm: Wrong MM-name!")

    @property
    def t(self):
        return self._t
     
    @t.setter
    def t(self,mm_type):
        self._t=mm_type

class Molecule(object):
    def __init__(self, a_tuple=(), mm=None, t_list=[], name="Molecule"):
        if t_list is None:
            t_list = []
        self._sl=[]
        self._rl=[]
        self._cl=[]
        self._tl=t_list
        self._mm=mm
        self._br=[]
        self._db=[]
        self._tb=[]
        for a in a_tuple:
            #if not self._mm==None: a.mm=self._mm
            a.mm=self._mm
            self._sl.append(a.s)        
            self._rl.append(a.r)        
            #a.c
            self._cl.append(a.c)        
            #a.t
            self._tl.append(a.t)        
        if a_tuple:
            self._na=len(a_tuple)
            self._rl=array(self._rl,'d')
        else:
            self._na=None
        self._n=name
        self._conn=[]
   
    @property 
    def sl(self):
        return self._sl

    @property 
    def mm(self):
        return self._mm

    @mm.setter 
    def mm(self,mm):
        self._mm=mm

    @property 
    def rl(self):
        return self._rl

    @rl.setter
    def rl(self,rl_list):
        assert len(rl_list)==self._na
        self._rl=array(rl_list,'d')        
   
    @property 
    def cl(self):
        return self._cl

    @property 
    def tl(self):
        return self._tl

    @tl.setter 
    def tl(self,t_list):
        self._tl=t_list

    @property 
    def na(self):
        return self._na

    @property 
    def name(self):
        return self._n

    @property 
    def conn(self):
        return self._conn

    @property
    def br(self):
        return self._br

    @property
    def db(self):
        return self._db

    @property
    def tb(self):
        return self._tb

    def set_conn(self):
        wdv_dist_tol=0.45
        for i in range(self._na):
            self._conn.append([])
            for j in range(self._na):
                if (j!=i):
                    l1=self._sl[i]
                    l2=self._sl[j]
                    if self.calc_dist(i,j)<float(cov_radii[l1])+float(cov_radii[l2])+wdv_dist_tol: 
                        self._conn[i].append(j)
        return self._conn

    def calc_dist(self,i1,i2):
        assert i1<=self._na
        assert i2<=self._na
        return dist(self._rl[i1],self._rl[i2])

    def calc_angle(self,i1,i2,i3):
        assert i1<=self._na
        assert i2<=self._na
        assert i3<=self._na
        return angle(self._rl[i1],self._rl[i2],self._rl[i3])

    def calc_dihedral(self,i1,i2,i3,i4):
        assert i1<=self._na
        assert i2<=self._na
        assert i3<=self._na
        assert i4<=self._na
        return dihedral(self._rl[i1],self._rl[i2],self._rl[i3],self._rl[i4])

    def calc_bsplit(self,i1,i2):
        assert i1<=self._na
        assert i2<=self._na
        assert not self._conn==[]
        orig_nset=set([i1])
        orig_pset=set([i2])
        def bnd_spl_rec(nset,pset,conn):
            #print nset,pset
            if pset==set():
                return nset,pset
            else:	
                new_pset=set()
            for a1 in pset:
                rset=set(conn[a1])-nset
                #new_pset=new_pset|set(rset)
                new_pset=new_pset|rset
                #return bnd_spl_rec(nset|pset|new_pset,new_pset,conn)
            return bnd_spl_rec(nset|pset,new_pset,conn)
        nset,pset=bnd_spl_rec(orig_nset,orig_pset,self._conn)
        return list(nset-orig_nset-orig_pset)

    def bond_tra(self,bond_tup,bond_tra_val):
        assert len(bond_tup)==2
        t1,t2=bond_tup        
        u=norm_vec(self._rl[t1],self._rl[t2])
        bsl=self.calc_bsplit(t1,t2)
        self._rl[t2]=tra(self._rl[t2],u,bond_tra_val)
        for ai in bsl:
            self._rl[ai]=tra(self._rl[ai],u,bond_tra_val)

    def angl_rot(self,angl_tup,angl_rot_val):
        assert len(angl_tup)==3
        t1,t2,t3=angl_tup        
        u=uangl(self._rl[t1],self._rl[t2],self._rl[t3])
        bsl=self.calc_bsplit(t1,t2)
        for ai in bsl:
            delr=rotu(self._rl[ai]-self._rl[t2],u,angl_rot_val)
            self._rl[ai]=self._rl[t2]+delr

    def diha_rot(self,diha_tup,diha_rot_val):
        assert len(diha_tup)==4
        #self.set_conn()
        t1,t2,t3,t4=diha_tup
        u=norm_vec(self._rl[t2],self._rl[t3])
        bsl=self.calc_bsplit(t2,t3)
        for ai in bsl:
            delr=rotu(self._rl[ai]-self._rl[t2],u,diha_rot_val)
            self._rl[ai]=self._rl[t2]+delr
 
    def benz_ring(self):
        #mol=["C","C","C","H","H","H","H","H","H","C","C","C","C","C","C","C","C","C","C","C"]
        #conn=[[1,9],[0,2],[1,11],[0],[1],[2],[9],[10],[11],[0,10],[9,11],[10,2],[13,17],[12,14],[13,15],[14,16],[15,17],[16,12],[17,19],[18,0]]
        mol=self._sl
        #self.set_conn()
        conn=self._conn

        carbs=[]
        for i in range(len(mol)):
           if mol[i]=='C': carbs.append(i)

        brings=[]
        tcarbs1=carbs[:]
        while not tcarbs1==[]:
           c1=tcarbs1.pop(0)
           tcarbs2=tcarbs1[:]
           fbreak=False
           for a2 in conn[c1]:
              i2=0
              for tc2 in tcarbs2:
                 if a2==tc2:
                    c2=a2
                    tcarbs3=tcarbs2[:i2]+tcarbs2[i2+1:]
                    for a3 in conn[c2]:
                       i3=0
                       for tc3 in tcarbs3:
                          if a3==tc3:
                             c3=a3
                             tcarbs4=tcarbs3[:i3]+tcarbs3[i3+1:]
                             for a4 in conn[c3]:
                                i4=0
                                for tc4 in tcarbs4:
                                   if a4==tc4:
                                      c4=a4
                                      tcarbs5=tcarbs4[:i4]+tcarbs4[i4+1:]
                                      for a5 in conn[c4]:
                                         i5=0
                                         for tc5 in tcarbs5:
                                            if a5==tc5:
                                               c5=a5
                                               tcarbs6=tcarbs5[:i5]+tcarbs5[i5+1:]
                                               for a6 in conn[c5]:
                                                  i6=0
                                                  for tc6 in tcarbs6:
                                                     if a6==tc6:
                                                        c6=a6
                                                        tcarbs7=tcarbs6[:i6]+tcarbs6[i6+1:]
                                                        for a7 in conn[c6]:
                                                           if a7==c1:
                                                              brings+=[c1,c2,c3,c4,c5,c6]
                                                              tcarbs1=tcarbs7[:]
                                                              fbreak=True
                                                              break
                                                     if fbreak: break
                                                     i6+=1
                                                  if fbreak: break
                                            if fbreak: break
                                            i5+=1
                                         if fbreak: break
                                   if fbreak: break
                                   i4+=1
                                if fbreak: break
                          if fbreak: break
                          i3+=1      
                       if fbreak: break
                 if fbreak: break
                 i2+=1
              if fbreak: break
        self._br=brings
        return brings
 
    def bond_ord(self):
        nbrings=[]
        for i in range(self._na):
           inb=False
           for j in self._br:
              if i==j: inb=True; break
           if not inb: nbrings.append(i)

        #for i in range(self._na):
        #print nbrings,self._conn
        for i in nbrings:
           #print i, self._conn[i]
           for j in self._conn[i]:
              if j>i: 
                 if self._sl[i]=='P' and self._sl[j]=='O' or self._sl[i]=='O' and self._sl[j]=='P':
                    s,d,t=bond_ords['PO']
                    #if d>self.calc_dist(i,j)>t: print i,j,"double PO"; self._db.append((i,j))
                    #if t>self.calc_dist(i,j): print i,j,"triple PO"; self._tb.append((i,j))
                    if d>self.calc_dist(i,j)>t: self._db.append((i,j))
                    if t>self.calc_dist(i,j): self._tb.append((i,j))
                 elif self._sl[i]=='P' and self._sl[j]=='C' or self._sl[i]=='C' and self._sl[j]=='P':
                    s,d,t=bond_ords['PC']
                    #if d>self.calc_dist(i,j)>t: print i,j,"double PC"; self._db.append((i,j))
                    #if t>self.calc_dist(i,j): print i,j,"triple PC"; self._tb.append((i,j))
                    if d>self.calc_dist(i,j)>t: self._db.append((i,j))
                    if t>self.calc_dist(i,j): self._tb.append((i,j))
                 elif self._sl[i]=='O' and self._sl[j]=='C' or self._sl[i]=='C' and self._sl[j]=='O':
                    s,d,t=bond_ords['OC']
                    #if d>self.calc_dist(i,j)>t: print "double CO"; self._db.append((i,j))
                    #if t>self.calc_dist(i,j): print "triple CO"; self._tb.append((i,j))
                    if d>self.calc_dist(i,j)>t: self._db.append((i,j))
                    if t>self.calc_dist(i,j): self._tb.append((i,j))
                 elif self._sl[i]=='C' and self._sl[j]=='C':
                    s,d,t=bond_ords['CC']
                    #if d>self.calc_dist(i,j)>t: print "double OO"; self._db.append((i,j))
                    #if t>self.calc_dist(i,j): print "triple OO"; self._tb.append((i,j))
                    if d>self.calc_dist(i,j)>t: self._db.append((i,j))
                    if t>self.calc_dist(i,j): self._tb.append((i,j))
                 #else: print "Uknown bond type: ",self._sl[i],"-",self._sl[j]
        #print self._db
        return self._db
 
if __name__=="__main__":
    a=Atom('C',(1.,2.,3.),6.0,"mm3")
    print(a.s)
    print(a.r)
    a.r=[0.,0.,0.]
    print(a.r)
    print(a.c)
    a.c=5.0
    print(a.c)
    a.t=50
    print(a.t)
    a.t
    print(a.t)

    b=Atom('P',(1.,0.,0.),mm="mm3")
    c=Atom('H',(1.,1.,0.),mm="mm3")
    d=Atom('C',(1.,1.,-1.),mm="mm3")
    m=Molecule((a,b,c,d),"mmff94",name="Mol2")
    print(m.mm, m.tl)
    print(m.sl)
    print(m.rl)
    print(m.cl)
    print("tl: ", m.tl)
    print(m.na)
    print(m.name)
    print(m.conn)
    m.set_conn()
    print(m.conn)
    print(m.calc_dist(0,1))
    print(m.calc_angle(1,0,2))
    print(m.calc_dihedral(0,1,2,3))
    m.rl=((1.,1.,1.),(0.,0.,0.),(2.,2.,2.),(3.,3.,3.))
    #m.rl=((1.,1.,1.),(0.,0.,0.),(2.,2.,2.))
    print(m.rl)

    m=Molecule()
    print(m.sl)
    print(m.rl)
    print(m.cl)
    print(m.tl)
    print(m.na)
    print(m.name)
    print(m.conn)

    a=Atom('C',(0.,-1.,0.),mm="mmff94")
    #b=Atom('P',(1.,0.,0.),mm="mmff94")
    b=Atom('P',(0.,0.,0.),mm="mmff94")
    c=Atom('C',(2.,0.,0.),mm="mmff94")
    d=Atom('H',(3.,1.,0.),mm="mmff94")
    m=Molecule((a,b,c,d),"mm3")
    print("tl: ", m.tl)
    m.set_conn()
    print(m.conn)
    print(m.calc_bsplit(0,1))
    print("1", m.rl)
    #m.diha_rot((0,1,2,3),pi/2.)
    m.angl_rot((0,1,2),pi/2.)
    print("2", m.rl)
    m.bond_ord()
    print(m.db)
    print(m.benz_ring())
