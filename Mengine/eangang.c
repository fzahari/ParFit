# Copyright (C) 2004-2016 Kevin E. Gilbert, dba Serena Software
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Kevin E. Gilbert
# Serena Software
# Box 3076
# Bloomington, IN 47402-3076
#
# gilbert@serenasoft.com

#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "angles.h"
#include "derivs.h"
#include "hess.h"
#include "utility.h"
#include "units.h"
#include "minim_values.h"
#include "angang.h"
#include "eangang.h"

void eangang()
{
      int i,k,l,ia,ib,ic,id,ie;
      double e,dt1,dt2,angle,dot,cosine, angle1;
      double xia,yia,zia;
      double xib,yib,zib;
      double xic,yic,zic;
      double xid,yid,zid;
      double xie,yie,zie;
      double xab,yab,zab,rab2;
      double xcb,ycb,zcb,rcb2;
      double xdb,ydb,zdb,rdb2;
      double xeb,yeb,zeb,reb2;

      energies.eangang = 0.0;
      if (minim_values.iprint)
      {
        fprintf(pcmoutfile,"\nAngle-Angle Interactions\n");
        fprintf(pcmoutfile,"             Angle 1                     Angle2          Thet1     Thet0      Thet2    Thet0   Const         Eaa\n");
      }

      for (i=0; i < angang.nangang; i++)
      {
          k = angang.iaa[i][0];
          l = angang.iaa[i][1];
          ia = angles.i13[k][0];
          ib = angles.i13[k][1];
          ic = angles.i13[k][2];
          id = angles.i13[l][0];
          ie = angles.i13[l][2];
          if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use || atom[ie].use)
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
            xie = atom[ie].x;
            yie = atom[ie].y;
            zie = atom[ie].z;
            xab = xia - xib;
            yab = yia - yib;
            zab = zia - zib;
            rab2 = xab*xab + yab*yab + zab*zab;
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
            xdb = xid - xib;
            ydb = yid - yib;
            zdb = zid - zib;
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb;
            xeb = xie - xib;
            yeb = yie - yib;
            zeb = zie - zib;
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb;
            if (rab2*rcb2*rdb2*reb2 != 0.0)
            {
               dot = xab*xcb + yab*ycb + zab*zcb;
               cosine = dot / sqrt(rab2*rcb2);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);
               angle1 = angle;
               dt1 = angle - angles.anat[k];
               dot = xdb*xeb + ydb*yeb + zdb*zeb;
               cosine = dot / sqrt(rdb2*reb2);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);
               dt2 = angle - angles.anat[l];
               e = units.aaunit * angang.fa[i] * dt1 * dt2;
               energies.eangang += e;
               atom[ia].energy += e;
               atom[ib].energy += e;
               atom[ic].energy += e;
               atom[id].energy += e;
               atom[ie].energy += e;
               if (minim_values.iprint)
                 fprintf(pcmoutfile,"AA: %2s(%-3d)- %2s(%-3d)- %2s(%-3d) %2s(%-3d)- %2s(%-3d)- %2s(%-3d): %-8.3f %-8.3f - %-8.3f %-8.3f %-8.3f = %8.4f\n",
                     atom[ia].name, ia, atom[ib].name, ib, atom[ic].name, ic, atom[id].name, id, atom[ib].name, ib, atom[ie].name, ie,angle1,
                     angles.anat[k], angle, angles.anat[l], angang.fa[i],e);
            }
          }        
      }

}

void eangang1()
{
      int i,k,l,ia,ib,ic,id,ie;
      double e,dt1,dt2,deddt1,deddt2,angle,dot,cosine;
      double xia,yia,zia;
      double xib,yib,zib;
      double xic,yic,zic;
      double xid,yid,zid;
      double xie,yie,zie;
      double xab,yab,zab,rab2;
      double xcb,ycb,zcb,rcb2;
      double xdb,ydb,zdb,rdb2;
      double xeb,yeb,zeb,reb2;
      double xp,yp,zp,rp;
      double xq,yq,zq,rq;
      double terma,termc,termd,terme;
      double dedxia,dedyia,dedzia;
      double dedxib,dedyib,dedzib;
      double dedxic,dedyic,dedzic;
      double dedxid,dedyid,dedzid;
      double dedxie,dedyie,dedzie;

      energies.eangang = 0.0;
      for (i=0; i <= natom; i++)
      {
          deriv.deaa[i][0] = 0.0;
          deriv.deaa[i][1] = 0.0;
          deriv.deaa[i][2] = 0.0;
      }

      for (i=0; i < angang.nangang; i++)
      {
          k = angang.iaa[i][0];
          l = angang.iaa[i][1];
          ia = angles.i13[k][0];
          ib = angles.i13[k][1];
          ic = angles.i13[k][2];
          id = angles.i13[l][0];
          ie = angles.i13[l][2];
          if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use || atom[ie].use)
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
            xie = atom[ie].x;
            yie = atom[ie].y;
            zie = atom[ie].z;
            xab = xia - xib;
            yab = yia - yib;
            zab = zia - zib;
            rab2 = xab*xab + yab*yab + zab*zab;
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
            xdb = xid - xib;
            ydb = yid - yib;
            zdb = zid - zib;
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb;
            xeb = xie - xib;
            yeb = yie - yib;
            zeb = zie - zib;
            reb2 = xeb*xeb + yeb*yeb + zeb*zeb;
            xp = ycb*zab - zcb*yab;
            yp = zcb*xab - xcb*zab;
            zp = xcb*yab - ycb*xab;
            rp = sqrt(xp*xp + yp*yp + zp*zp);
            xq = yeb*zdb - zeb*ydb;
            yq = zeb*xdb - xeb*zdb;
            zq = xeb*ydb - yeb*xdb;
            rq = sqrt(xq*xq + yq*yq + zq*zq);
            if (rp*rq != 0.0)
            {
               dot = xab*xcb + yab*ycb + zab*zcb;
               cosine = dot / sqrt(rab2*rcb2);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);
               dt1 = angle - angles.anat[k];
               dot = xdb*xeb + ydb*yeb + zdb*zeb;
               cosine = dot / sqrt(rdb2*reb2);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);
               dt2 = angle - angles.anat[l];
               e = units.aaunit * angang.fa[i] * dt1 * dt2;

               deddt1 = radian * e / dt1;
               deddt2 = radian * e / dt2;
               terma = -deddt1 / (rab2*rp);
               termc = deddt1 / (rcb2*rp);
               dedxia = terma * (yab*zp-zab*yp);
               dedyia = terma * (zab*xp-xab*zp);
               dedzia = terma * (xab*yp-yab*xp);
               dedxic = termc * (ycb*zp-zcb*yp);
               dedyic = termc * (zcb*xp-xcb*zp);
               dedzic = termc * (xcb*yp-ycb*xp);

               termd = -deddt2 / (rdb2*rq);
               terme = deddt2 / (reb2*rq);
               dedxid = termd * (ydb*zq-zdb*yq);
               dedyid = termd * (zdb*xq-xdb*zq);
               dedzid = termd * (xdb*yq-ydb*xq);
               dedxie = terme * (yeb*zq-zeb*yq);
               dedyie = terme * (zeb*xq-xeb*zq);
               dedzie = terme * (xeb*yq-yeb*xq);

               dedxib = -dedxia - dedxic - dedxid - dedxie;
               dedyib = -dedyia - dedyic - dedyid - dedyie;
               dedzib = -dedzia - dedzic - dedzid - dedzie;

               deriv.deaa[ia][0] += dedxia;
               deriv.deaa[ia][1] += dedyia;
               deriv.deaa[ia][2] += dedzia;
               
               deriv.deaa[ib][0] += dedxib;
               deriv.deaa[ib][1] += dedyib;
               deriv.deaa[ib][2] += dedzib;

               deriv.deaa[ic][0] += dedxic;
               deriv.deaa[ic][1] += dedyic;
               deriv.deaa[ic][2] += dedzic;

               deriv.deaa[id][0] += dedxid;
               deriv.deaa[id][1] += dedyid;
               deriv.deaa[id][2] += dedzid;
               
               deriv.deaa[ie][0] += dedxie;
               deriv.deaa[ie][1] += dedyie;
               deriv.deaa[ie][2] += dedzie;

               virial.virx += xab*dedxia + xcb*dedxic + xdb*dedxid + xeb*dedxie;
               virial.viry += yab*dedyia + ycb*dedyic + ydb*dedyid + yeb*dedyie;
               virial.virz += zab*dedzia + zcb*dedzic + zdb*dedzid + zeb*dedzie;

               energies.eangang += e;
            }
          }    
      }

}

void eangang2(int i)
{
    int j,k,iang,ij;
    int ia,ib,ic,id,ie;
    double old,eps;
    double **d0;     //  double d0[MAXATOM][3];

    d0 = dmatrix(0,natom, 0,3);
    eps = 1.0e-7;
    for (iang = 0; iang < angang.nangang; iang++)
    {
        j = angang.iaa[iang][0];
        k = angang.iaa[iang][1];
        ia = angles.i13[j][0];
        ib = angles.i13[j][1];
        ic = angles.i13[j][2];
        id = angles.i13[k][0];
        ie = angles.i13[k][2];

        if (ia == i || ib == i || ic == i || id == i || ie == i)
        {
            if (id == ia || id == ic)
            {
                id = ie;
                ie = 0;
            } else if (ie == ia || ie == ic)
                ie = 0;

            eangang2b(iang);
            for (ij=0; ij < 3; ij++)
            {
                d0[ia][ij] = deriv.deaa[ia][ij];
                d0[ib][ij] = deriv.deaa[ib][ij];
                d0[ic][ij] = deriv.deaa[ic][ij];
                d0[id][ij] = deriv.deaa[id][ij];
                if (ie != 0)
                   d0[ie][ij] = deriv.deaa[ie][ij];
            }
  //  x  components
            old = atom[i].x;
            atom[i].x += eps;
            eangang2b(iang);
            atom[i].x = old;
            for (ij=0; ij < 3; ij++)
            {
                hess.hessx[ia][ij] += (deriv.deaa[ia][ij]-d0[ia][ij])/eps;
                hess.hessx[ib][ij] += (deriv.deaa[ib][ij]-d0[ib][ij])/eps;
                hess.hessx[ic][ij] += (deriv.deaa[ic][ij]-d0[ic][ij])/eps;
                hess.hessx[id][ij] += (deriv.deaa[id][ij]-d0[id][ij])/eps;
                if (ie != 0)
                  hess.hessx[ie][ij] += (deriv.deaa[ie][ij]-d0[ie][ij])/eps;
            }
  //  y  components
            old = atom[i].y;
            atom[i].y += eps;
            eangang2b(iang);
            atom[i].y = old;
            for (ij=0; ij < 3; ij++)
            {
                hess.hessy[ia][ij] += (deriv.deaa[ia][ij]-d0[ia][ij])/eps;
                hess.hessy[ib][ij] += (deriv.deaa[ib][ij]-d0[ib][ij])/eps;
                hess.hessy[ic][ij] += (deriv.deaa[ic][ij]-d0[ic][ij])/eps;
                hess.hessy[id][ij] += (deriv.deaa[id][ij]-d0[id][ij])/eps;
                if (ie != 0)
                  hess.hessy[ie][ij] += (deriv.deaa[ie][ij]-d0[ie][ij])/eps;
            }
  //  z  components
            old = atom[i].z;
            atom[i].z += eps;
            eangang2b(iang);
            atom[i].z = old;
            for (ij=0; ij < 3; ij++)
            {
                hess.hessz[ia][ij] += (deriv.deaa[ia][ij]-d0[ia][ij])/eps;
                hess.hessz[ib][ij] += (deriv.deaa[ib][ij]-d0[ib][ij])/eps;
                hess.hessz[ic][ij] += (deriv.deaa[ic][ij]-d0[ic][ij])/eps;
                hess.hessz[id][ij] += (deriv.deaa[id][ij]-d0[id][ij])/eps;
                if (ie != 0)
                  hess.hessz[ie][ij] += (deriv.deaa[ie][ij]-d0[ie][ij])/eps;
            }
        }
    }
    free_dmatrix(d0, 0,natom, 0,3);
}

void eangang2b(int i)
{
    int j,k,ia,ib,ic,id,ie, ijj;
    double angle,dot,cosine,term;
    double dt1,dt2,deddt1,deddt2;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xid,yid,zid;
    double xie,yie,zie;
    double xab,yab,zab,rab2;
    double xcb,ycb,zcb,rcb2;
    double xdb,ydb,zdb,rdb2;
    double xeb,yeb,zeb,reb2;
    double xp,yp,zp,rp;
    double xq,yq,zq,rq;
    double terma,termc,termd,terme;
    double dedxia,dedyia,dedzia;
    double dedxib,dedyib,dedzib;
    double dedxic,dedyic,dedzic;
    double dedxid,dedyid,dedzid;
    double dedxie,dedyie,dedzie;

      j = angang.iaa[i][0];
      k = angang.iaa[i][1];
      ia = angles.i13[j][0];
      ib = angles.i13[j][1];
      ic = angles.i13[j][2];
      id = angles.i13[k][0];
      ie = angles.i13[k][2];
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
      xie = atom[ie].x;
      yie = atom[ie].y;
      zie = atom[ie].z;

      for (ijj=0; ijj < 3; ijj++)
      {
          deriv.deaa[ia][ijj] = 0.0;
          deriv.deaa[ib][ijj] = 0.0;
          deriv.deaa[ic][ijj] = 0.0;
          deriv.deaa[id][ijj] = 0.0;
          deriv.deaa[ie][ijj] = 0.0;
      }

      xab = xia - xib;
      yab = yia - yib;
      zab = zia - zib;
      rab2 = xab*xab + yab*yab + zab*zab;
      xcb = xic - xib;
      ycb = yic - yib;
      zcb = zic - zib;
      rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
      xdb = xid - xib;
      ydb = yid - yib;
      zdb = zid - zib;
      rdb2 = xdb*xdb + ydb*ydb + zdb*zdb;
      xeb = xie - xib;
      yeb = yie - yib;
      zeb = zie - zib;
      reb2 = xeb*xeb + yeb*yeb + zeb*zeb;
      xp = ycb*zab - zcb*yab;
      yp = zcb*xab - xcb*zab;
      zp = xcb*yab - ycb*xab;
      rp = sqrt(xp*xp + yp*yp + zp*zp);
      xq = yeb*zdb - zeb*ydb;
      yq = zeb*xdb - xeb*zdb;
      zq = xeb*ydb - yeb*xdb;
      rq = sqrt(xq*xq + yq*yq + zq*zq);
      if (rp*rq != 0.0)
      {
         dot = xab*xcb + yab*ycb + zab*zcb;
         cosine = dot / sqrt(rab2*rcb2);
         if (cosine < -1.0)
           cosine = -1.0;
         if (cosine > 1.0)
           cosine = 1.0;
         angle = radian * acos(cosine);
         dt1 = angle - angles.anat[j];
         dot = xdb*xeb + ydb*yeb + zdb*zeb;
         cosine = dot / sqrt(rdb2*reb2);
         if (cosine < -1.0)
           cosine = -1.0;
         if (cosine > 1.0)
           cosine = 1.0;
         angle = radian * acos(cosine);
         dt2 = angle - angles.anat[k];

         term = radian * units.aaunit * angang.fa[i];
         deddt1 = term * dt2;
         deddt2 = term * dt1;

         terma = -deddt1 / (rab2*rp);
         termc = deddt1 / (rcb2*rp);
         dedxia = terma * (yab*zp-zab*yp);
         dedyia = terma * (zab*xp-xab*zp);
         dedzia = terma * (xab*yp-yab*xp);
         dedxic = termc * (ycb*zp-zcb*yp);
         dedyic = termc * (zcb*xp-xcb*zp);
         dedzic = termc * (xcb*yp-ycb*xp);
         termd = -deddt2 / (rdb2*rq);
         terme = deddt2 / (reb2*rq);
         dedxid = termd * (ydb*zq-zdb*yq);
         dedyid = termd * (zdb*xq-xdb*zq);
         dedzid = termd * (xdb*yq-ydb*xq);
         dedxie = terme * (yeb*zq-zeb*yq);
         dedyie = terme * (zeb*xq-xeb*zq);
         dedzie = terme * (xeb*yq-yeb*xq);
         dedxib = -dedxia - dedxic - dedxid - dedxie;
         dedyib = -dedyia - dedyic - dedyid - dedyie;
         dedzib = -dedzia - dedzic - dedzid - dedzie;

         deriv.deaa[ia][0] += dedxia;
         deriv.deaa[ia][1] += dedyia;
         deriv.deaa[ia][2] += dedzia;
         deriv.deaa[ib][0] += dedxib;
         deriv.deaa[ib][1] += dedyib;
         deriv.deaa[ib][2] += dedzib;
         deriv.deaa[ic][0] += dedxic;
         deriv.deaa[ic][1] += dedyic;
         deriv.deaa[ic][2] += dedzic;
         deriv.deaa[id][0] += dedxid;
         deriv.deaa[id][1] += dedyid;
         deriv.deaa[id][2] += dedzid;
         deriv.deaa[ie][0] += dedxie;
         deriv.deaa[ie][1] += dedyie;
         deriv.deaa[ie][2] += dedzie;
      }
}
    
