// Copyright (C) 2004-2016 Kevin E. Gilbert, dba Serena Software
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Kevin E. Gilbert
// Serena Software
// Box 3076
// Bloomington, IN 47402-3076
//
// gilbert@serenasoft.com

#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "angles.h"
#include "derivs.h"
#include "hess.h"
#include "utility.h"
#include "units.h"
#include "minim_values.h"
#include "opbend.h"
        
void eopbend()
{
      int i,ia,ib,ic,id;
      double e,angle,dot,cosine;
      double dt,dt2,dt3,dt4;
      double xia,yia,zia;
      double xib,yib,zib;
      double xic,yic,zic;
      double xid,yid,zid;
      double xad,yad,zad;
      double xbd,ybd,zbd,rbd2;
      double xcd,ycd,zcd;
      double xpd,ypd,zpd,rpd2;
      double xt,yt,zt,rt2,delta;

      energies.eopb = 0.0;

      if (minim_values.iprint)
      {
          fprintf(pcmoutfile,"\nOut of Plane Bending Terms\n");
          fprintf(pcmoutfile,"    At1   At2   At3   At4      Angle     Oopconst     Eoop\n");
      }

      for (i=0; i < angles.nang; i++)
      {
          if (angles.angin[i] == TRUE)
          {
              ia = angles.i13[i][0];
              ib = angles.i13[i][1];
              ic = angles.i13[i][2];
              id = angles.i13[i][3];
              if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use)
              {
                 xia = atom[ia].x;
                 yia = atom[ia].y;
                 zia = atom[ia].z;
                 xib = atom[ib].x;
                 yib = atom[ib].y;
                 zib = atom[ib].z;
                 xic = atom[ic].x;
                 yic = atom[ic].y;
                 zic = atom[ic].z;
                 xid = atom[id].x;
                 yid = atom[id].y;
                 zid = atom[id].z;
                 xad = xia - xid;
                 yad = yia - yid;
                 zad = zia - zid;
                 xbd = xib - xid;
                 ybd = yib - yid;
                 zbd = zib - zid;
                 xcd = xic - xid;
                 ycd = yic - yid;
                 zcd = zic - zid;
                 xt = yad*zcd - zad*ycd;
                 yt = zad*xcd - xad*zcd;
                 zt = xad*ycd - yad*xcd;
                 rt2 = xt*xt + yt*yt + zt*zt;
                 delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2;
                 xpd = xbd + xt*delta;
                 ypd = ybd + yt*delta;
                 zpd = zbd + zt*delta;
                 rbd2 = xbd*xbd + ybd*ybd + zbd*zbd;
                 rpd2 = xpd*xpd + ypd*ypd + zpd*zpd;
                 if (rbd2 != 0.0 && rpd2 != 0.0)
                 {
                    dot = xbd*xpd + ybd*ypd + zbd*zpd;
                    cosine = dot / sqrt(rbd2*rpd2);
                    if (cosine < -1.0)
                       cosine = -1.0;
                    if (cosine > 1.0)
                       cosine = 1.0;
                    angle = radian * acos(cosine);
                    dt = angle;
                    dt2 = dt * dt;
                    dt3 = dt2 * dt;
                    dt4 = dt2 * dt2;
                    e = units.angunit * angles.copb[i] * dt2
                        * (1.0+units.cang*dt+units.qang*dt2+units.pang*dt3+units.sang*dt4);
                    energies.eopb += e;
                    atom[ia].energy += e;
                    atom[ib].energy += e;
                    atom[ic].energy += e;
                    atom[id].energy += e;
                    if (minim_values.iprint)
                      fprintf(pcmoutfile,"Oop: %-4d- %-4d- %-4d- %-4d    %-8.3f    %-8.3f = %-8.4f\n",ia,ib,ic,id,
                        angle, angles.copb[i],e);
                 }
                 
              }
          }
      }
}

void eopbend1()
{
      int i,ia,ib,ic,id;
      double e,angle,dot,cosine;
      double dt,dt2,dt3,dt4;
      double deddt,termb,termp,term;
      double xia,yia,zia;
      double xib,yib,zib;
      double xic,yic,zic;
      double xid,yid,zid;
      double xad,yad,zad;
      double xbd,ybd,zbd,rbd2;
      double xcd,ycd,zcd;
      double xpd,ypd,zpd,rpd2;
      double xt,yt,zt,rt2,ptrt2;
      double xm,ym,zm,rm,delta,delta2;
      double dedxia,dedyia,dedzia;
      double dedxib,dedyib,dedzib;
      double dedxic,dedyic,dedzic;
      double dedxid,dedyid,dedzid;
      double dedxip,dedyip,dedzip;
      double dpdxia,dpdyia,dpdzia;
      double dpdxic,dpdyic,dpdzic;

      energies.eopb = 0.0;
      for (i=0; i <= natom; i++)
      {
          deriv.deopb[i][0] = 0.0;
          deriv.deopb[i][1] = 0.0;
          deriv.deopb[i][2] = 0.0;
      }

      for (i=0; i < angles.nang; i++)
      {
          if (angles.angin[i] == TRUE)
          {
              ia = angles.i13[i][0];
              ib = angles.i13[i][1];
              ic = angles.i13[i][2];
              id = angles.i13[i][3];
              if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use)
              {
                 xia = atom[ia].x;
                 yia = atom[ia].y;
                 zia = atom[ia].z;
                 xib = atom[ib].x;
                 yib = atom[ib].y;
                 zib = atom[ib].z;
                 xic = atom[ic].x;
                 yic = atom[ic].y;
                 zic = atom[ic].z;
                 xid = atom[id].x;
                 yid = atom[id].y;
                 zid = atom[id].z;
                 xad = xia - xid;
                 yad = yia - yid;
                 zad = zia - zid;
                 xbd = xib - xid;
                 ybd = yib - yid;
                 zbd = zib - zid;
                 xcd = xic - xid;
                 ycd = yic - yid;
                 zcd = zic - zid;
                 xt = yad*zcd - zad*ycd;
                 yt = zad*xcd - xad*zcd;
                 zt = xad*ycd - yad*xcd;
                 rt2 = xt*xt + yt*yt + zt*zt;
                 delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2;
                 xpd = xbd + xt*delta;
                 ypd = ybd + yt*delta;
                 zpd = zbd + zt*delta;
                 rbd2 = xbd*xbd + ybd*ybd + zbd*zbd;
                 rpd2 = xpd*xpd + ypd*ypd + zpd*zpd;
                 xm = ypd*zbd - zpd*ybd;
                 ym = zpd*xbd - xpd*zbd;
                 zm = xpd*ybd - ypd*xbd;
                 rm = sqrt(xm*xm + ym*ym + zm*zm);
                 if (rm != 0.0)
                 {
                    dot = xbd*xpd + ybd*ypd + zbd*zpd;
                    cosine = dot / sqrt(rbd2*rpd2);
                    if (cosine < -1.0)
                       cosine = -1.0;
                    if (cosine > 1.0)
                       cosine = 1.0;
                    angle = radian * acos(cosine);
                    dt = angle;
                    dt2 = dt * dt;
                    dt3 = dt2 * dt;
                    dt4 = dt2 * dt2;
                    e = units.angunit * angles.copb[i] * dt2
                        * (1.0+units.cang*dt+units.qang*dt2+units.pang*dt3+units.sang*dt4);

                    deddt = units.angunit * angles.copb[i] * dt
                         * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                             + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);
                    deddt = deddt * radian;
                    termb = -deddt / (rbd2*rm);
                    termp = deddt / (rpd2*rm);
                    dedxib = termb * (ybd*zm-zbd*ym);
                    dedyib = termb * (zbd*xm-xbd*zm);
                    dedzib = termb * (xbd*ym-ybd*xm);
                    dedxip = termp * (ypd*zm-zpd*ym);
                    dedyip = termp * (zpd*xm-xpd*zm);
                    dedzip = termp * (xpd*ym-ypd*xm);

                    delta2 = 2.0 * delta;
                    ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2;
                    term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd);
                    dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2;
                    term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd);
                    dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2;
                    term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd);
                    dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2;
                    term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad);
                    dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2;
                    term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad);
                    dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2;
                    term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad);
                    dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2;

                    dedxia = dpdxia;
                    dedyia = dpdyia;
                    dedzia = dpdzia;
                    dedxic = dpdxic;
                    dedyic = dpdyic;
                    dedzic = dpdzic;
                    dedxid = -dedxia - dedxib - dedxic;
                    dedyid = -dedyia - dedyib - dedyic;
                    dedzid = -dedzia - dedzib - dedzic;

                    energies.eopb += e;
                    deriv.deopb[ia][0] += dedxia;
                    deriv.deopb[ia][1] += dedyia;
                    deriv.deopb[ia][2] += dedzia;
                   
                    deriv.deopb[ib][0] += dedxib;
                    deriv.deopb[ib][1] += dedyib;
                    deriv.deopb[ib][2] += dedzib;
                    
                    deriv.deopb[ic][0] += dedxic;
                    deriv.deopb[ic][1] += dedyic;
                    deriv.deopb[ic][2] += dedzic;
                    
                    deriv.deopb[id][0] += dedxid;
                    deriv.deopb[id][1] += dedyid;
                    deriv.deopb[id][2] += dedzid;
                    virial.virx += xad*dedxia + xbd*dedxib + xcd*dedxic;
                    virial.viry += yad*dedyia + ybd*dedyib + ycd*dedyic;
                    virial.virz += zad*dedzia + zbd*dedzib + zcd*dedzic;
                 }
                 
              }
          }
      }
}

void eopbend2(int i)
{
    int j,k;
    int ia,ib,ic,id;
    double old,eps;
    double **d0; // d0[MAXATOM][3];

    d0 = dmatrix(0,natom,0,3);
    eps = 1.0e-7;

    for (k=0; k < angles.nang; k++)
    {
        if (angles.angin[k] == TRUE)
        {
            ia = angles.i13[k][0];
            ib = angles.i13[k][1];
            ic = angles.i13[k][2];
            id = angles.i13[k][3];
            if (ia == i || ib == i || ic == i || id == i)
            {
                eopbend2b(k);
                for (j = 0; j < 3; j++)
                {
                    d0[ia][j] = deriv.deopb[ia][j];
                    d0[ib][j] = deriv.deopb[ib][j];
                    d0[ic][j] = deriv.deopb[ic][j];
                    d0[id][j] = deriv.deopb[id][j];
                }
  // numerical x
                old = atom[i].x;
                atom[i].x += eps;
                eopbend2b(k);
                atom[i].x = old;
                for (j=0; j < 3; j++)
                {
                    hess.hessx[ia][j] += (deriv.deopb[ia][j]-d0[ia][j])/eps;
                    hess.hessx[ib][j] += (deriv.deopb[ib][j]-d0[ib][j])/eps;
                    hess.hessx[ic][j] += (deriv.deopb[ic][j]-d0[ic][j])/eps;
                    hess.hessx[id][j] += (deriv.deopb[id][j]-d0[id][j])/eps;
                }
  // numerical y
                old = atom[i].y;
                atom[i].y += eps;
                eopbend2b(k);
                atom[i].y = old;
                for (j=0; j < 3; j++)
                {
                    hess.hessy[ia][j] += (deriv.deopb[ia][j]-d0[ia][j])/eps;
                    hess.hessy[ib][j] += (deriv.deopb[ib][j]-d0[ib][j])/eps;
                    hess.hessy[ic][j] += (deriv.deopb[ic][j]-d0[ic][j])/eps;
                    hess.hessy[id][j] += (deriv.deopb[id][j]-d0[id][j])/eps;
                }
  // numerical z
                old = atom[i].z;
                atom[i].z += eps;
                eopbend2b(k);
                atom[i].z = old;
                for (j=0; j < 3; j++)
                {
                    hess.hessz[ia][j] += (deriv.deopb[ia][j]-d0[ia][j])/eps;
                    hess.hessz[ib][j] += (deriv.deopb[ib][j]-d0[ib][j])/eps;
                    hess.hessz[ic][j] += (deriv.deopb[ic][j]-d0[ic][j])/eps;
                    hess.hessz[id][j] += (deriv.deopb[id][j]-d0[id][j])/eps;
                }
            }
        }
    }
    free_dmatrix( d0 ,0,natom,0,3);
}

void eopbend2b(int k)
{
    int j,ia,ib,ic,id;
    double angle,dot,cosine;
    double dt,dt2,dt3,dt4;
    double deddt,termb,termp,term;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xid,yid,zid;
    double xad,yad,zad;
    double xbd,ybd,zbd,rbd2;
    double xcd,ycd,zcd;
    double xpd,ypd,zpd,rpd2;
    double xt,yt,zt,rt2,ptrt2;
    double xm,ym,zm,rm,delta,delta2;
    double dedxia,dedyia,dedzia;
    double dedxib,dedyib,dedzib;
    double dedxic,dedyic,dedzic;
    double dedxid,dedyid,dedzid;
    double dedxip,dedyip,dedzip;
    double dpdxia,dpdyia,dpdzia;
    double dpdxic,dpdyic,dpdzic;
    
    ia = angles.i13[k][0];
    ib = angles.i13[k][1];
    ic = angles.i13[k][2];
    id = angles.i13[k][3];
                    
      xia = atom[ia].x;
      yia = atom[ia].y;
      zia = atom[ia].z;
      xib = atom[ib].x;
      yib = atom[ib].y;
      zib = atom[ib].z;
      xic = atom[ic].x;
      yic = atom[ic].y;
      zic = atom[ic].z;
      xid = atom[id].x;
      yid = atom[id].y;
      zid = atom[id].z;

      for (j=0; j < 3; j++)
      {
          deriv.deopb[ia][j] = 0.0;
          deriv.deopb[ib][j] = 0.0;
          deriv.deopb[ic][j] = 0.0;
          deriv.deopb[id][j] = 0.0;
      }
          
      xad = xia - xid;
      yad = yia - yid;
      zad = zia - zid;
      xbd = xib - xid;
      ybd = yib - yid;
      zbd = zib - zid;
      xcd = xic - xid;
      ycd = yic - yid;
      zcd = zic - zid;
      xt = yad*zcd - zad*ycd;
      yt = zad*xcd - xad*zcd;
      zt = xad*ycd - yad*xcd;
      rt2 = xt*xt + yt*yt + zt*zt;
      delta = -(xt*xbd + yt*ybd + zt*zbd) / rt2;
      xpd = xbd + xt*delta;
      ypd = ybd + yt*delta;
      zpd = zbd + zt*delta;
      rbd2 = xbd*xbd + ybd*ybd + zbd*zbd;
      rpd2 = xpd*xpd + ypd*ypd + zpd*zpd;
      xm = ypd*zbd - zpd*ybd;
      ym = zpd*xbd - xpd*zbd;
      zm = xpd*ybd - ypd*xbd;
      rm = sqrt(xm*xm + ym*ym + zm*zm);
      if (rm != 0.0)
      {
         dot = xbd*xpd + ybd*ypd + zbd*zpd;
         cosine = dot / sqrt(rbd2*rpd2);
         if (cosine < -1.0)
            cosine = -1.0;
         if (cosine > 1.0)
            cosine = 1.0;
         angle = radian * acos(cosine);
         dt = angle;
         dt2 = dt * dt;
         dt3 = dt2 * dt;
         dt4 = dt2 * dt2;
         deddt = units.angunit * angles.copb[k] * dt
                   * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                   + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);
            
         deddt = deddt * radian;
         termb = -deddt / (rbd2*rm);
         termp = deddt / (rpd2*rm);
         dedxib = termb * (ybd*zm-zbd*ym);
         dedyib = termb * (zbd*xm-xbd*zm);
         dedzib = termb * (xbd*ym-ybd*xm);
         dedxip = termp * (ypd*zm-zpd*ym);
         dedyip = termp * (zpd*xm-xpd*zm);
         dedzip = termp * (xpd*ym-ypd*xm);

         delta2 = 2.0 * delta;
         ptrt2 = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2;
         term = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd);
         dpdxia = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2;
         term = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd);
         dpdyia = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2;
         term = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd);
         dpdzia = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2;
         term = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad);
         dpdxic = delta*(zad*dedyip-yad*dedzip) + term*ptrt2;
         term = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad);
         dpdyic = delta*(xad*dedzip-zad*dedxip) + term*ptrt2;
         term = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad);
         dpdzic = delta*(yad*dedxip-xad*dedyip) + term*ptrt2;

         dedxia = dpdxia;
         dedyia = dpdyia;
         dedzia = dpdzia;
         dedxic = dpdxic;
         dedyic = dpdyic;
         dedzic = dpdzic;
         dedxid = -dedxia - dedxib - dedxic;
         dedyid = -dedyia - dedyib - dedyic;
         dedzid = -dedzia - dedzib - dedzic;

         deriv.deopb[ia][0] = dedxia;
         deriv.deopb[ia][1] = dedyia;
         deriv.deopb[ia][2] = dedzia;
         deriv.deopb[ib][0] = dedxib;
         deriv.deopb[ib][1] = dedyib;
         deriv.deopb[ib][2] = dedzib;
         deriv.deopb[ic][0] = dedxic;
         deriv.deopb[ic][1] = dedyic;
         deriv.deopb[ic][2] = dedzic;
         deriv.deopb[id][0] = dedxid;
         deriv.deopb[id][1] = dedyid;
         deriv.deopb[id][2] = dedzid;
      }
}


