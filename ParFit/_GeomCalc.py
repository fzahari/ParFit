#!/usr/bin/env python

#from numpy import rad2deg,arccos,cross,dot,allclose,cos,sin,matrix,array,pi
from numpy import arccos,cross,dot,allclose,cos,sin,matrix,array,pi
from numpy.linalg import norm

def norm_vec(r1,r2):
    #return r2-r1
    return (r2-r1)/norm(r2-r1)

def vangle(v1,v2):
    #return rad2deg(arccos(dot(v1/norm(v1),v2/norm(v2))))
    vdot=dot(v1/norm(v1),v2/norm(v2))
    if allclose(vdot,1.0): return 0.0 
    if allclose(vdot,-1.0): return pi/2. 
    return arccos(vdot)

def dist(r1,r2):
        return norm(r1-r2)

def angle(r1,r2,r3):
        return vangle(r1-r2,r3-r2)

def dihedral(r1,r2,r3,r4):
    v1 = norm_vec(r1, r2)
    v2 = norm_vec(r2, r3)
    v3 = norm_vec(r3, r4)
    v1xv2 = cross(v1,v2)
    v2xv3 = cross(v2,v3)
    sign=dot(v1xv2,v3)
    dih=vangle(v1xv2,v2xv3)
    #if sign<0.: dih=360.0-dih
    if sign<0.: dih=2.0*pi-dih
    if allclose(dih,2.0*pi): dih=0.0
    return dih

def rotu(r,u,th):

    try:
        if not allclose(norm(u),1.):
            raise Exception("u is not a unit vector!")
        
        dru=dot(r,u)
        nr=r-dru*u
        nnr=norm(nr)
        unr=nr/nnr
        v=cross(u,unr)

    except Exception, msg:
        print msg

    return cos(th)*nr+sin(th)*nnr*v+dru*u

def rotu2(r,u,th):

    try:
        if not allclose(norm(u),1.):
            raise Exception("u is not a unit vector!")

        c=cos(th)
        o=1.-c
        s=sin(th)

        rot=[[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
        rot=matrix(rot)

        rot[0,0]=c+u[0]**2*o;        rot[0,1]=u[0]*u[1]*o-u[2]*s; rot[0,2]=u[0]*u[2]*o+u[1]*s
        rot[1,0]=u[1]*u[0]*o+u[2]*s; rot[1,1]=c+u[1]**2*o;        rot[1,2]=u[1]*u[2]*o-u[0]*s
        rot[2,0]=u[2]*u[0]*o-u[1]*s; rot[2,1]=u[2]*u[1]*o+u[0]*s; rot[2,2]=c+u[2]**2*o

        r=matrix(r)
        r=r.T

        rot_r=rot*r

    except Exception, msg:
        print msg

    return array(rot_r.T)[0]

def tra(r,u,d):

   try:
      if not allclose(norm(u),1.):
         raise Exception("u is not a unit vector!")

      rt=r+d*u

   except Exception, msg:
      print msg

   return rt

def uangl(r1,r2,r3):
   v=cross(r1-r2,r3-r2) 
   return v/norm(v)

if __name__=="__main__":
    a=array([0.,0.,0.])
    b=array([1.,0.,0.])
    c=array([1.,1.,0.])
    d=array([1.,1.,1.])
    e=array([0.,1.,0.])
    #
    print dist(a,b)
    print angle(b,a,c)
    print dihedral(a,b,c,d)
    d=array([1.,1.,-1.])
    print dihedral(a,b,c,d)
    #
    #u=array([0.,0.,1.01])
    u=array([0.,0.,1.])
    r=array([2.,0.,0.])
    print rotu(r,u,pi/2.)
    # 
    print rotu2(r,u,pi/2)
    # 
    d /= norm(d)
    print rotu(r,d,pi)
    # 
    print rotu2(r,d,pi)
    print tra(r,u,3.)
    #
    print uangl(b,a,e)
