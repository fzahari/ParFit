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
#include "torsions.h"
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "units.h"
#include "minim_values.h"
#include "strtor.h"

void estrtor()
{
      int i,ii,k,ia,ib,ic,id;
      double e,rcb,dr,rt2,ru2,rtru;
      double xt,yt,zt,xu,yu,zu;
      double cosine,v3,s3,phi3;
      double xia,yia,zia,xib,yib,zib;
      double xic,yic,zic,xid,yid,zid;
      double xba,yba,zba,xcb,ycb,zcb;
      double xdc,ydc,zdc;

      energies.estrtor = 0.0;
      if (minim_values.iprint && strtor.nstrtor > 0)
      {
          fprintf(pcmoutfile,"\nStretch-Torsion Terms\n");
          fprintf(pcmoutfile,"          At1     At2    At3     At4     Angle     Dr   StrTorConst    Phi3   Estb\n");
      }

      for (i = 0; i < strtor.nstrtor; i++)
      {
          ii = strtor.ist[i][0];
          ia = torsions.i14[ii][0];
          ib = torsions.i14[ii][1];
          ic = torsions.i14[ii][2];
          id = torsions.i14[ii][3];
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
            xba = xib - xia;
            yba = yib - yia;
            zba = zib - zia;
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            xdc = xid - xic;
            ydc = yid - yic;
            zdc = zid - zic;
            xt = yba*zcb - ycb*zba;
            yt = zba*xcb - zcb*xba;
            zt = xba*ycb - xcb*yba;
            xu = ycb*zdc - ydc*zcb;
            yu = zcb*xdc - zdc*xcb;
            zu = xcb*ydc - xdc*ycb;
            rt2 = xt*xt + yt*yt + zt*zt;
            ru2 = xu*xu + yu*yu + zu*zu;
            rtru = sqrt(rt2 * ru2);
            if (rtru != 0.0)
            {
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;
               v3 = torsions.v3[ii];
               s3 = torsions.ph3[ii];
               phi3 = 1.0 + s3*(4.0*pow(cosine,3.0) - 3.0*cosine);
               k = strtor.ist[i][1];
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               dr = rcb - bonds_ff.bl[k];
               e = units.storunit * strtor.stcon[i] * dr * phi3;
               energies.estrtor += e;
               atom[ia].energy += e;
               atom[ib].energy += e;
               atom[ic].energy += e;
               atom[id].energy += e;
               if (minim_values.iprint)
                 fprintf(pcmoutfile,"StrTor: %2s(%-3d)-%2s(%-3d)-%2s(%-3d)-%2s(%-3d) : %-8.3f %-8.3f %-8.3f %-8.3f = %-8.4f\n",
                    atom[ia].name, ia,atom[ib].name, ib,atom[ic].name, ic,atom[id].name, id,cosine,dr,
                    strtor.stcon[i],phi3, e);
            }
          }
      }
}

void estrtor1()
{
      int i,ii,k,ia,ib,ic,id;
      double e,sine,dedphi,dr,ddr;
      double ddrdx,ddrdy,ddrdz,term;
      double v3,s3,phi3,dphi3;
      double cosine,cosine2,cosine3;
      double rcb,rt2,ru2,rtru;
      double xt,yt,zt,xu,yu,zu;
      double xia,yia,zia,xib,yib,zib;
      double xic,yic,zic,xid,yid,zid;
      double xba,yba,zba,xcb,ycb,zcb;
      double xdc,ydc,zdc;
      double xca,yca,zca,xdb,ydb,zdb;
      double xtu,ytu,ztu;
      double dphidxt,dphidyt,dphidzt;
      double dphidxu,dphidyu,dphidzu;
      double dphidxia,dphidyia,dphidzia;
      double dphidxib,dphidyib,dphidzib;
      double dphidxic,dphidyic,dphidzic;
      double dphidxid,dphidyid,dphidzid;
      double dedxia,dedyia,dedzia;
      double dedxib,dedyib,dedzib;
      double dedxic,dedyic,dedzic;
      double dedxid,dedyid,dedzid;

      energies.estrtor = 0.0;
      for (i=0; i <= natom; i++)
      {
          deriv.destor[i][0] = 0.0;
          deriv.destor[i][1] = 0.0;
          deriv.destor[i][2] = 0.0;
      }

      for (i = 0; i < strtor.nstrtor; i++)
      {
          ii = strtor.ist[i][0];
          ia = torsions.i14[ii][0];
          ib = torsions.i14[ii][1];
          ic = torsions.i14[ii][2];
          id = torsions.i14[ii][3];
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
            xba = xib - xia;
            yba = yib - yia;
            zba = zib - zia;
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            xdc = xid - xic;
            ydc = yid - yic;
            zdc = zid - zic;
            xt = yba*zcb - ycb*zba;
            yt = zba*xcb - zcb*xba;
            zt = xba*ycb - xcb*yba;
            xu = ycb*zdc - ydc*zcb;
            yu = zcb*xdc - zdc*xcb;
            zu = xcb*ydc - xdc*ycb;
            xtu = yt*zu - yu*zt;
            ytu = zt*xu - zu*xt;
            ztu = xt*yu - xu*yt;
            rt2 = xt*xt + yt*yt + zt*zt;
            ru2 = xu*xu + yu*yu + zu*zu;
            rtru = sqrt(rt2 * ru2);
            if (rtru != 0.0)
            {
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);
               v3 = torsions.v3[ii];
               s3 = torsions.ph3[ii];
               cosine2 = cosine * cosine;
               cosine3 = cosine2 * cosine;
               phi3 = 1.0 + s3*(4.0*cosine3 - 3.0*cosine);
               dphi3 = s3 * (12.0*cosine2 - 3.0);

               k = strtor.ist[i][1];
               dr = rcb - bonds_ff.bl[k];
               term = units.storunit * strtor.stcon[i];
               e = term * dr * phi3;
               dedphi = -term * dr * sine * dphi3;
               ddr = term * phi3 / rcb;
               ddrdx = xcb * ddr;
               ddrdy = ycb * ddr;
               ddrdz = zcb * ddr;
               xca = xic - xia;
               yca = yic - yia;
               zca = zic - zia;
               xdb = xid - xib;
               ydb = yid - yib;
               zdb = zid - zib;
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb);
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb);
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb);
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb);
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb);
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb);

               dphidxia = zcb*dphidyt - ycb*dphidzt;
               dphidyia = xcb*dphidzt - zcb*dphidxt;
               dphidzia = ycb*dphidxt - xcb*dphidyt;
               dphidxib = yca*dphidzt - zca*dphidyt + zdc*dphidyu - ydc*dphidzu;
               dphidyib = zca*dphidxt - xca*dphidzt + xdc*dphidzu - zdc*dphidxu;
               dphidzib = xca*dphidyt - yca*dphidxt + ydc*dphidxu - xdc*dphidyu;
               dphidxic = zba*dphidyt - yba*dphidzt + ydb*dphidzu - zdb*dphidyu;
               dphidyic = xba*dphidzt - zba*dphidxt + zdb*dphidxu - xdb*dphidzu;
               dphidzic = yba*dphidxt - xba*dphidyt + xdb*dphidyu - ydb*dphidxu;
               dphidxid = zcb*dphidyu - ycb*dphidzu;
               dphidyid = xcb*dphidzu - zcb*dphidxu;
               dphidzid = ycb*dphidxu - xcb*dphidyu;
               dedxia = dedphi*dphidxia;
               dedyia = dedphi*dphidyia;
               dedzia = dedphi*dphidzia;
               dedxib = dedphi*dphidxib - ddrdx;
               dedyib = dedphi*dphidyib - ddrdy;
               dedzib = dedphi*dphidzib - ddrdz;
               dedxic = dedphi*dphidxic + ddrdx;
               dedyic = dedphi*dphidyic + ddrdy;
               dedzic = dedphi*dphidzic + ddrdz;
               dedxid = dedphi*dphidxid;
               dedyid = dedphi*dphidyid;
               dedzid = dedphi*dphidzid;
               
               energies.estrtor += e;
               deriv.destor[ia][0] += dedxia;
               deriv.destor[ia][1] += dedyia;
               deriv.destor[ia][2] += dedzia;

               deriv.destor[ib][0] += dedxib;
               deriv.destor[ib][1] += dedyib;
               deriv.destor[ib][2] += dedzib;

               deriv.destor[ic][0] += dedxic;
               deriv.destor[ic][1] += dedyic;
               deriv.destor[ic][2] += dedzic;

               deriv.destor[id][0] += dedxid;
               deriv.destor[id][1] += dedyid;
               deriv.destor[id][2] += dedzid;
               
               virial.virx = virial.virx - xba*dedxia + xdc*dedxid
                          + xcb*(dedxic+dedxid);
               virial.viry = virial.viry - yba*dedyia + ydc*dedyid
                          + ycb*(dedyic+dedyid);
               virial.virz = virial.virz - zba*dedzia + zdc*dedzid
                          + zcb*(dedzic+dedzid);
            }
          }
      }
}

void estrtor2(int iatom)
{
    int i,k,ia,ib,ic,id, ii;
    double dedphi,d2edphi2,sine;
    double v3,s3,phi3,dphi3,d2phi3;
    double cosine,cosine2,cosine3;
    double xia,yia,zia,xib,yib,zib;
    double xic,yic,zic,xid,yid,zid;
    double xba,yba,zba,xcb,ycb,zcb;
    double xdc,ydc,zdc;
    double xca,yca,zca,xdb,ydb,zdb;
    double xt,yt,zt,xu,yu,zu,xtu,ytu,ztu;
    double rt2,ru2,rtru,rcb;
    double dr,term,ddr,d2dr;
    double ddrdx,ddrdy,ddrdz;
    double d2drdxx,d2drdyy,d2drdzz;
    double d2drdxy,d2drdxz,d2drdyz;
    double dphidxt,dphidyt,dphidzt;
    double dphidxu,dphidyu,dphidzu;
    double dphidxia,dphidyia,dphidzia;
    double dphidxib,dphidyib,dphidzib;
    double dphidxic,dphidyic,dphidzic;
    double dphidxid,dphidyid,dphidzid;
    double xycb2,xzcb2,yzcb2;
    double rcbxt,rcbyt,rcbzt,rcbt2;
    double rcbxu,rcbyu,rcbzu,rcbu2;
    double dphidxibt,dphidyibt,dphidzibt;
    double dphidxibu,dphidyibu,dphidzibu;
    double dphidxict,dphidyict,dphidzict;
    double dphidxicu,dphidyicu,dphidzicu;
    double dxia,dyia,dzia,dxib,dyib,dzib;
    double dxic,dyic,dzic,dxid,dyid,dzid;
    double dxiaxia,dyiayia,dziazia,dxibxib,dyibyib,dzibzib;
    double dxicxic,dyicyic,dziczic,dxidxid,dyidyid,dzidzid;
    double dxiayia,dxiazia,dyiazia,dxibyib,dxibzib,dyibzib;
    double dxicyic,dxiczic,dyiczic,dxidyid,dxidzid,dyidzid;
    double dxiaxib,dxiayib,dxiazib,dyiaxib,dyiayib,dyiazib;
    double dziaxib,dziayib,dziazib,dxiaxic,dxiayic,dxiazic;
    double dyiaxic,dyiayic,dyiazic,dziaxic,dziayic,dziazic;
    double dxiaxid,dxiayid,dxiazid,dyiaxid,dyiayid,dyiazid;
    double dziaxid,dziayid,dziazid,dxibxic,dxibyic,dxibzic;
    double dyibxic,dyibyic,dyibzic,dzibxic,dzibyic,dzibzic;
    double dxibxid,dxibyid,dxibzid,dyibxid,dyibyid,dyibzid;
    double dzibxid,dzibyid,dzibzid,dxicxid,dxicyid,dxiczid;
    double dyicxid,dyicyid,dyiczid,dzicxid,dzicyid,dziczid;
    
      for (i = 0; i < strtor.nstrtor; i++)
      {
          ii = strtor.ist[i][0];
          ia = torsions.i14[ii][0];
          ib = torsions.i14[ii][1];
          ic = torsions.i14[ii][2];
          id = torsions.i14[ii][3];
          if (iatom == ia || iatom == ib || iatom == ic || iatom == id)
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
            xba = xib - xia;
            yba = yib - yia;
            zba = zib - zia;
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            xdc = xid - xic;
            ydc = yid - yic;
            zdc = zid - zic;
            xt = yba*zcb - ycb*zba;
            yt = zba*xcb - zcb*xba;
            zt = xba*ycb - xcb*yba;
            xu = ycb*zdc - ydc*zcb;
            yu = zcb*xdc - zdc*xcb;
            zu = xcb*ydc - xdc*ycb;
            xtu = yt*zu - yu*zt;
            ytu = zt*xu - zu*xt;
            ztu = xt*yu - xu*yt;
            rt2 = xt*xt + yt*yt + zt*zt;
            ru2 = xu*xu + yu*yu + zu*zu;
            rtru = sqrt(rt2 * ru2);
            if (rtru != 0.0)
            {
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;
               sine = (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru);
               v3 = torsions.v3[ii];
               s3 = torsions.ph3[ii];
               cosine2 = cosine * cosine;
               cosine3 = cosine2 * cosine;
               phi3 = 1.0 + s3*(4.0*cosine3 - 3.0*cosine);
               dphi3 = s3 * (12.0*cosine2 - 3.0);
               d2phi3 = -s3 * (36.0*cosine3 - 27.0*cosine);

               k = strtor.ist[i][1];
               dr = rcb - bonds_ff.bl[k];
               term = units.storunit * strtor.stcon[i];
               dedphi = -term * sine * dphi3;
               d2edphi2 = term * dr * d2phi3;
               ddr = 1.0 / rcb;
               ddrdx = xcb * ddr;
               ddrdy = ycb * ddr;
               ddrdz = zcb * ddr;
               d2dr = -term * phi3 / (rcb*rcb*rcb);
               d2drdxx = (xcb*xcb-rcb*rcb) * d2dr;
               d2drdyy = (ycb*ycb-rcb*rcb) * d2dr;
               d2drdzz = (zcb*zcb-rcb*rcb) * d2dr;
               d2drdxy = xcb * ycb * d2dr;
               d2drdxz = xcb * zcb * d2dr;
               d2drdyz = ycb * zcb * d2dr;

               xca = xic - xia;
               yca = yic - yia;
               zca = zic - zia;
               xdb = xid - xib;
               ydb = yid - yib;
               zdb = zid - zib;
               dphidxt = (yt*zcb - ycb*zt) / (rt2*rcb);
               dphidyt = (zt*xcb - zcb*xt) / (rt2*rcb);
               dphidzt = (xt*ycb - xcb*yt) / (rt2*rcb);
               dphidxu = -(yu*zcb - ycb*zu) / (ru2*rcb);
               dphidyu = -(zu*xcb - zcb*xu) / (ru2*rcb);
               dphidzu = -(xu*ycb - xcb*yu) / (ru2*rcb);

               xycb2 = xcb*xcb + ycb*ycb;
               xzcb2 = xcb*xcb + zcb*zcb;
               yzcb2 = ycb*ycb + zcb*zcb;
               rcbxt = -2.0 * rcb * dphidxt;
               rcbyt = -2.0 * rcb * dphidyt;
               rcbzt = -2.0 * rcb * dphidzt;
               rcbt2 = rcb * rt2;
               rcbxu = 2.0 * rcb * dphidxu;
               rcbyu = 2.0 * rcb * dphidyu;
               rcbzu = 2.0 * rcb * dphidzu;
               rcbu2 = rcb * ru2;
               dphidxibt = yca*dphidzt - zca*dphidyt;
               dphidxibu = zdc*dphidyu - ydc*dphidzu;
               dphidyibt = zca*dphidxt - xca*dphidzt;
               dphidyibu = xdc*dphidzu - zdc*dphidxu;
               dphidzibt = xca*dphidyt - yca*dphidxt;
               dphidzibu = ydc*dphidxu - xdc*dphidyu;
               dphidxict = zba*dphidyt - yba*dphidzt;
               dphidxicu = ydb*dphidzu - zdb*dphidyu;
               dphidyict = xba*dphidzt - zba*dphidxt;
               dphidyicu = zdb*dphidxu - xdb*dphidzu;
               dphidzict = yba*dphidxt - xba*dphidyt;
               dphidzicu = xdb*dphidyu - ydb*dphidxu;
               dphidxia = zcb*dphidyt - ycb*dphidzt;
               dphidyia = xcb*dphidzt - zcb*dphidxt;
               dphidzia = ycb*dphidxt - xcb*dphidyt;
               dphidxib = dphidxibt + dphidxibu;
               dphidyib = dphidyibt + dphidyibu;
               dphidzib = dphidzibt + dphidzibu;
               dphidxic = dphidxict + dphidxicu;
               dphidyic = dphidyict + dphidyicu;
               dphidzic = dphidzict + dphidzicu;
               dphidxid = zcb*dphidyu - ycb*dphidzu;
               dphidyid = xcb*dphidzu - zcb*dphidxu;
               dphidzid = ycb*dphidxu - xcb*dphidyu;
               dxia = dedphi * dphidxia;
               dyia = dedphi * dphidyia;
               dzia = dedphi * dphidzia;
               dxib = dedphi * dphidxib;
               dyib = dedphi * dphidyib;
               dzib = dedphi * dphidzib;
               dxic = dedphi * dphidxic;
               dyic = dedphi * dphidyic;
               dzic = dedphi * dphidzic;
               dxid = dedphi * dphidxid;
               dyid = dedphi * dphidyid;
               dzid = dedphi * dphidzid;
               dedphi = dedphi * dr;

         dxiaxia = rcbxt*dphidxia;
         dxiayia = rcbxt*dphidyia - zcb*rcb/rt2;
         dxiazia = rcbxt*dphidzia + ycb*rcb/rt2;
         dxiaxib = rcbxt*dphidxibt + xcb*(zca*ycb-yca*zcb)/rcbt2;
         dxiayib = rcbxt*dphidyibt + dphidzt + (xca*zcb*xcb+zca*yzcb2)/rcbt2;
         dxiazib = rcbxt*dphidzibt - dphidyt - (xca*ycb*xcb+yca*yzcb2)/rcbt2;
         dxiaxic = rcbxt*dphidxict + xcb*xt/rcbt2;
         dxiayic = rcbxt*dphidyict - dphidzt - (xba*zcb*xcb+zba*yzcb2)/rcbt2;
         dxiazic = rcbxt*dphidzict + dphidyt + (xba*ycb*xcb+yba*yzcb2)/rcbt2;
         dxiaxid = 0.0;
         dxiayid = 0.0;
         dxiazid = 0.0;
         dyiayia = rcbyt*dphidyia;
         dyiazia = rcbyt*dphidzia - xcb*rcb/rt2;
         dyiaxib = rcbyt*dphidxibt - dphidzt - (yca*zcb*ycb+zca*xzcb2)/rcbt2;
         dyiayib = rcbyt*dphidyibt + ycb*(xca*zcb-zca*xcb)/rcbt2;
         dyiazib = rcbyt*dphidzibt + dphidxt + (yca*xcb*ycb+xca*xzcb2)/rcbt2;
         dyiaxic = rcbyt*dphidxict + dphidzt + (yba*zcb*ycb+zba*xzcb2)/rcbt2;
         dyiayic = rcbyt*dphidyict + ycb*yt/rcbt2;
         dyiazic = rcbyt*dphidzict - dphidxt - (yba*xcb*ycb+xba*xzcb2)/rcbt2;
         dyiaxid = 0.0;
         dyiayid = 0.0;
         dyiazid = 0.0;
         dziazia = rcbzt*dphidzia;
         dziaxib = rcbzt*dphidxibt + dphidyt + (zca*ycb*zcb+yca*xycb2)/rcbt2;
         dziayib = rcbzt*dphidyibt - dphidxt - (zca*xcb*zcb+xca*xycb2)/rcbt2;
         dziazib = rcbzt*dphidzibt + zcb*(yca*xcb-xca*ycb)/rcbt2;
         dziaxic = rcbzt*dphidxict - dphidyt - (zba*ycb*zcb+yba*xycb2)/rcbt2;
         dziayic = rcbzt*dphidyict + dphidxt + (zba*xcb*zcb+xba*xycb2)/rcbt2;
         dziazic = rcbzt*dphidzict + zcb*zt/rcbt2;
         dziaxid = 0.0;
         dziayid = 0.0;
         dziazid = 0.0;
         dxibxic = -xcb*dphidxib/(rcb*rcb) - (yca*(zba*xcb+yt)-zca*(yba*xcb-zt))/rcbt2
            - 2.0*(yt*zba-yba*zt)*dphidxibt/rt2 - (zdc*(ydb*xcb+zu)-ydc*(zdb*xcb-yu))/rcbu2
            + 2.0*(yu*zdb-ydb*zu)*dphidxibu/ru2;
         dxibyic = -ycb*dphidxib/(rcb*rcb) + dphidzt + dphidzu - (yca*(zba*ycb-xt)+zca*(xba*xcb+zcb*zba))/rcbt2
            - 2.0*(zt*xba-zba*xt)*dphidxibt/rt2 + (zdc*(xdb*xcb+zcb*zdb)+ydc*(zdb*ycb+xu))/rcbu2
            + 2.0*(zu*xdb-zdb*xu)*dphidxibu/ru2;
         dxibxid = rcbxu*dphidxibu + xcb*xu/rcbu2;
         dxibyid = rcbyu*dphidxibu - dphidzu - (ydc*zcb*ycb+zdc*xzcb2)/rcbu2;
         dxibzid = rcbzu*dphidxibu + dphidyu + (zdc*ycb*zcb+ydc*xycb2)/rcbu2;
         dyibzib = ycb*dphidzib/(rcb*rcb) - (xca*(xca*xcb+zcb*zca)+yca*(ycb*xca+zt))/rcbt2
            - 2.0*(xt*zca-xca*zt)*dphidzibt/rt2 + (ydc*(xdc*ycb-zu)+xdc*(xdc*xcb+zcb*zdc))/rcbu2
            + 2.0*(xu*zdc-xdc*zu)*dphidzibu/ru2;
         dyibxic = -xcb*dphidyib/(rcb*rcb) - dphidzt - dphidzu + (xca*(zba*xcb+yt)+zca*(zba*zcb+ycb*yba))/rcbt2
            - 2.0*(yt*zba-yba*zt)*dphidyibt/rt2 - (zdc*(zdb*zcb+ycb*ydb)+xdc*(zdb*xcb-yu))/rcbu2
            + 2.0*(yu*zdb-ydb*zu)*dphidyibu/ru2;
         dyibyic = -ycb*dphidyib/(rcb*rcb) - (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
            - 2.0*(zt*xba-zba*xt)*dphidyibt/rt2 - (xdc*(zdb*ycb+xu)-zdc*(xdb*ycb-zu))/rcbu2
            + 2.0*(zu*xdb-zdb*xu)*dphidyibu/ru2;
         dyibxid = rcbxu*dphidyibu + dphidzu + (xdc*zcb*xcb+zdc*yzcb2)/rcbu2;
         dyibyid = rcbyu*dphidyibu + ycb*yu/rcbu2;
         dyibzid = rcbzu*dphidyibu - dphidxu - (zdc*xcb*zcb+xdc*xycb2)/rcbu2;
         dzibxic = -xcb*dphidzib/(rcb*rcb) + dphidyt + dphidyu - (xca*(yba*xcb-zt)+yca*(zba*zcb+ycb*yba))/rcbt2
            - 2.0*(yt*zba-yba*zt)*dphidzibt/rt2 + (ydc*(zdb*zcb+ycb*ydb)+xdc*(ydb*xcb+zu))/rcbu2
            + 2.0*(yu*zdb-ydb*zu)*dphidzibu/ru2;
         dzibzic = -zcb*dphidzib/(rcb*rcb) - (xca*(yba*zcb+xt)-yca*(xba*zcb-yt))/rcbt2
            - 2.0*(xt*yba-xba*yt)*dphidzibt/rt2 - (ydc*(xdb*zcb+yu)-xdc*(ydb*zcb-xu))/rcbu2
            + 2.0*(xu*ydb-xdb*yu)*dphidzibu/ru2;
         dzibxid = rcbxu*dphidzibu - dphidyu - (xdc*ycb*xcb+ydc*yzcb2)/rcbu2;
         dzibyid = rcbyu*dphidzibu + dphidxu + (ydc*xcb*ycb+xdc*xzcb2)/rcbu2;
         dzibzid = rcbzu*dphidzibu + zcb*zu/rcbu2;
         dxicxid = rcbxu*dphidxicu - xcb*(zdb*ycb-ydb*zcb)/rcbu2;
         dxicyid = rcbyu*dphidxicu + dphidzu + (ydb*zcb*ycb+zdb*xzcb2)/rcbu2;
         dxiczid = rcbzu*dphidxicu - dphidyu - (zdb*ycb*zcb+ydb*xycb2)/rcbu2;
         dyicxid = rcbxu*dphidyicu - dphidzu - (xdb*zcb*xcb+zdb*yzcb2)/rcbu2;
         dyicyid = rcbyu*dphidyicu - ycb*(xdb*zcb-zdb*xcb)/rcbu2;
         dyiczid = rcbzu*dphidyicu + dphidxu + (zdb*xcb*zcb+xdb*xycb2)/rcbu2;
         dzicxid = rcbxu*dphidzicu + dphidyu + (xdb*ycb*xcb+ydb*yzcb2)/rcbu2;
         dzicyid = rcbyu*dphidzicu - dphidxu - (ydb*xcb*ycb+xdb*xzcb2)/rcbu2;
         dziczid = rcbzu*dphidzicu - zcb*(ydb*xcb-xdb*ycb)/rcbu2;
         dxidxid = rcbxu*dphidxid;
         dxidyid = rcbxu*dphidyid + zcb*rcb/ru2;
         dxidzid = rcbxu*dphidzid - ycb*rcb/ru2;
         dyidyid = rcbyu*dphidyid;
         dyidzid = rcbyu*dphidzid + xcb*rcb/ru2;
         dzidzid = rcbzu*dphidzid;

         dxibxib = -dxiaxib - dxibxic - dxibxid;
         dxibyib = -dyiaxib - dxibyic - dxibyid;
         dxibzib = -dxiazib - dzibxic - dzibxid;
         dxibzic = -dziaxib - dxibzib - dxibzid;
         dyibyib = -dyiayib - dyibyic - dyibyid;
         dyibzic = -dziayib - dyibzib - dyibzid;
         dzibzib = -dziazib - dzibzic - dzibzid;
         dzibyic = -dyiazib - dyibzib - dzibyid;
         dxicxic = -dxiaxic - dxibxic - dxicxid;
         dxicyic = -dyiaxic - dyibxic - dxicyid;
         dxiczic = -dziaxic - dzibxic - dxiczid;
         dyicyic = -dyiayic - dyibyic - dyicyid;
         dyiczic = -dziayic - dzibyic - dyiczid;
         dziczic = -dziazic - dzibzic - dziczid;

         if (iatom == ia)
         {
            hess.hessx[ia][0] += dedphi*dxiaxia + d2edphi2*dphidxia*dphidxia;
            hess.hessy[ia][0] += dedphi*dxiayia + d2edphi2*dphidxia*dphidyia;
            hess.hessz[ia][0] += dedphi*dxiazia + d2edphi2*dphidxia*dphidzia;
            hess.hessx[ia][1] += dedphi*dxiayia + d2edphi2*dphidxia*dphidyia;
            hess.hessy[ia][1] += dedphi*dyiayia + d2edphi2*dphidyia*dphidyia;
            hess.hessz[ia][1] += dedphi*dyiazia + d2edphi2*dphidyia*dphidzia;
            hess.hessx[ia][2] += dedphi*dxiazia + d2edphi2*dphidxia*dphidzia;
            hess.hessy[ia][2] += dedphi*dyiazia + d2edphi2*dphidyia*dphidzia;
            hess.hessz[ia][2] += dedphi*dziazia + d2edphi2*dphidzia*dphidzia;
            hess.hessx[ib][0] += dedphi*dxiaxib - dxia*ddrdx + d2edphi2*dphidxia*dphidxib;
            hess.hessy[ib][0] += dedphi*dyiaxib - dyia*ddrdx + d2edphi2*dphidyia*dphidxib;
            hess.hessz[ib][0] += dedphi*dziaxib - dzia*ddrdx + d2edphi2*dphidzia*dphidxib;
            hess.hessx[ib][1] += dedphi*dxiayib - dxia*ddrdy + d2edphi2*dphidxia*dphidyib;
            hess.hessy[ib][1] += dedphi*dyiayib - dyia*ddrdy + d2edphi2*dphidyia*dphidyib;
            hess.hessz[ib][1] += dedphi*dziayib - dzia*ddrdy + d2edphi2*dphidzia*dphidyib;
            hess.hessx[ib][2] += dedphi*dxiazib - dxia*ddrdz + d2edphi2*dphidxia*dphidzib;
            hess.hessy[ib][2] += dedphi*dyiazib - dyia*ddrdz + d2edphi2*dphidyia*dphidzib;
            hess.hessz[ib][2] += dedphi*dziazib - dzia*ddrdz + d2edphi2*dphidzia*dphidzib;
            hess.hessx[ic][0] += dedphi*dxiaxic + dxia*ddrdx + d2edphi2*dphidxia*dphidxic;
            hess.hessy[ic][0] += dedphi*dyiaxic + dyia*ddrdx + d2edphi2*dphidyia*dphidxic;
            hess.hessz[ic][0] += dedphi*dziaxic + dzia*ddrdx + d2edphi2*dphidzia*dphidxic;
            hess.hessx[ic][1] += dedphi*dxiayic + dxia*ddrdy + d2edphi2*dphidxia*dphidyic;
            hess.hessy[ic][1] += dedphi*dyiayic + dyia*ddrdy + d2edphi2*dphidyia*dphidyic;
            hess.hessz[ic][1] += dedphi*dziayic + dzia*ddrdy + d2edphi2*dphidzia*dphidyic;
            hess.hessx[ic][2] += dedphi*dxiazic + dxia*ddrdz + d2edphi2*dphidxia*dphidzic;
            hess.hessy[ic][2] += dedphi*dyiazic + dyia*ddrdz + d2edphi2*dphidyia*dphidzic;
            hess.hessz[ic][2] += dedphi*dziazic + dzia*ddrdz + d2edphi2*dphidzia*dphidzic;
            hess.hessx[id][0] += dedphi*dxiaxid + d2edphi2*dphidxia*dphidxid;
            hess.hessy[id][0] += dedphi*dyiaxid + d2edphi2*dphidyia*dphidxid;
            hess.hessz[id][0] += dedphi*dziaxid + d2edphi2*dphidzia*dphidxid;
            hess.hessx[id][1] += dedphi*dxiayid + d2edphi2*dphidxia*dphidyid;
            hess.hessy[id][1] += dedphi*dyiayid + d2edphi2*dphidyia*dphidyid;
            hess.hessz[id][1] += dedphi*dziayid + d2edphi2*dphidzia*dphidyid;
            hess.hessx[id][2] += dedphi*dxiazid + d2edphi2*dphidxia*dphidzid;
            hess.hessy[id][2] += dedphi*dyiazid + d2edphi2*dphidyia*dphidzid;
            hess.hessz[id][2] += dedphi*dziazid + d2edphi2*dphidzia*dphidzid;
         } else if (iatom == ib)
         {
            hess.hessx[ib][0] += dedphi*dxibxib + d2edphi2*dphidxib*dphidxib - 2.0*dxib*ddrdx + d2drdxx;
            hess.hessy[ib][0] += dedphi*dxibyib + d2edphi2*dphidxib*dphidyib - dyib*ddrdx - dxib*ddrdy + d2drdxy;
            hess.hessz[ib][0] += dedphi*dxibzib + d2edphi2*dphidxib*dphidzib - dzib*ddrdx - dxib*ddrdz + d2drdxz;
            hess.hessx[ib][1] += dedphi*dxibyib + d2edphi2*dphidxib*dphidyib - dxib*ddrdy - dyib*ddrdx + d2drdxy;
            hess.hessy[ib][1] += dedphi*dyibyib + d2edphi2*dphidyib*dphidyib - 2.0*dyib*ddrdy + d2drdyy;
            hess.hessz[ib][1] += dedphi*dyibzib + d2edphi2*dphidyib*dphidzib - dzib*ddrdy - dyib*ddrdz + d2drdyz;
            hess.hessx[ib][2] += dedphi*dxibzib + d2edphi2*dphidxib*dphidzib - dxib*ddrdz - dzib*ddrdx + d2drdxz;
            hess.hessy[ib][2] += dedphi*dyibzib + d2edphi2*dphidyib*dphidzib - dyib*ddrdz - dzib*ddrdy + d2drdyz;
            hess.hessz[ib][2] += dedphi*dzibzib + d2edphi2*dphidzib*dphidzib - 2.0*dzib*ddrdz + d2drdzz;
            hess.hessx[ia][0] += dedphi*dxiaxib - dxia*ddrdx + d2edphi2*dphidxib*dphidxia;
            hess.hessy[ia][0] += dedphi*dxiayib - dxia*ddrdy + d2edphi2*dphidyib*dphidxia;
            hess.hessz[ia][0] += dedphi*dxiazib - dxia*ddrdz + d2edphi2*dphidzib*dphidxia;
            hess.hessx[ia][1] += dedphi*dyiaxib - dyia*ddrdx + d2edphi2*dphidxib*dphidyia;
            hess.hessy[ia][1] += dedphi*dyiayib - dyia*ddrdy + d2edphi2*dphidyib*dphidyia;
            hess.hessz[ia][1] += dedphi*dyiazib - dyia*ddrdz + d2edphi2*dphidzib*dphidyia;
            hess.hessx[ia][2] += dedphi*dziaxib - dzia*ddrdx + d2edphi2*dphidxib*dphidzia;
            hess.hessy[ia][2] += dedphi*dziayib - dzia*ddrdy + d2edphi2*dphidyib*dphidzia;
            hess.hessz[ia][2] += dedphi*dziazib - dzia*ddrdz + d2edphi2*dphidzib*dphidzia;
            hess.hessx[ic][0] += dedphi*dxibxic + d2edphi2*dphidxib*dphidxic + (dxib-dxic)*ddrdx - d2drdxx;
            hess.hessy[ic][0] += dedphi*dyibxic + d2edphi2*dphidyib*dphidxic + dyib*ddrdx - dxic*ddrdy - d2drdxy;
            hess.hessz[ic][0] += dedphi*dzibxic + d2edphi2*dphidzib*dphidxic + dzib*ddrdx - dxic*ddrdz - d2drdxz;
            hess.hessx[ic][1] += dedphi*dxibyic + d2edphi2*dphidxib*dphidyic + dxib*ddrdy - dyic*ddrdx - d2drdxy;
            hess.hessy[ic][1] += dedphi*dyibyic + d2edphi2*dphidyib*dphidyic + (dyib-dyic)*ddrdy - d2drdyy;
            hess.hessz[ic][1] += dedphi*dzibyic + d2edphi2*dphidzib*dphidyic + dzib*ddrdy - dyic*ddrdz - d2drdyz;
            hess.hessx[ic][2] += dedphi*dxibzic + d2edphi2*dphidxib*dphidzic + dxib*ddrdz - dzic*ddrdx - d2drdxz;
            hess.hessy[ic][2] += dedphi*dyibzic + d2edphi2*dphidyib*dphidzic + dyib*ddrdz - dzic*ddrdy - d2drdyz;
            hess.hessz[ic][2] += dedphi*dzibzic + d2edphi2*dphidzib*dphidzic + (dzib-dzic)*ddrdz - d2drdzz;
            hess.hessx[id][0] += dedphi*dxibxid - dxid*ddrdx + d2edphi2*dphidxib*dphidxid;
            hess.hessy[id][0] += dedphi*dyibxid - dxid*ddrdy + d2edphi2*dphidyib*dphidxid;
            hess.hessz[id][0] += dedphi*dzibxid - dxid*ddrdz + d2edphi2*dphidzib*dphidxid;
            hess.hessx[id][1] += dedphi*dxibyid - dyid*ddrdx + d2edphi2*dphidxib*dphidyid;
            hess.hessy[id][1] += dedphi*dyibyid - dyid*ddrdy + d2edphi2*dphidyib*dphidyid;
            hess.hessz[id][1] += dedphi*dzibyid - dyid*ddrdz + d2edphi2*dphidzib*dphidyid;
            hess.hessx[id][2] += dedphi*dxibzid - dzid*ddrdx + d2edphi2*dphidxib*dphidzid;
            hess.hessy[id][2] += dedphi*dyibzid - dzid*ddrdy + d2edphi2*dphidyib*dphidzid;
            hess.hessz[id][2] += dedphi*dzibzid - dzid*ddrdz + d2edphi2*dphidzib*dphidzid;
         } else if (iatom == ic)
         {
            hess.hessx[ic][0] += dedphi*dxicxic + d2edphi2*dphidxic*dphidxic + 2.0*dxic*ddrdx + d2drdxx;
            hess.hessy[ic][0] += dedphi*dxicyic + d2edphi2*dphidxic*dphidyic + dyic*ddrdx + dxic*ddrdy + d2drdxy;
            hess.hessz[ic][0] += dedphi*dxiczic + d2edphi2*dphidxic*dphidzic + dzic*ddrdx + dxic*ddrdz + d2drdxz;
            hess.hessx[ic][1] += dedphi*dxicyic + d2edphi2*dphidxic*dphidyic + dxic*ddrdy + dyic*ddrdx + d2drdxy;
            hess.hessy[ic][1] += dedphi*dyicyic + d2edphi2*dphidyic*dphidyic + 2.0*dyic*ddrdy + d2drdyy;
            hess.hessz[ic][1] += dedphi*dyiczic + d2edphi2*dphidyic*dphidzic + dzic*ddrdy + dyic*ddrdz + d2drdyz;
            hess.hessx[ic][2] += dedphi*dxiczic + d2edphi2*dphidxic*dphidzic + dxic*ddrdz + dzic*ddrdx + d2drdxz;
            hess.hessy[ic][2] += dedphi*dyiczic + d2edphi2*dphidyic*dphidzic + dyic*ddrdz + dzic*ddrdy + d2drdyz;
            hess.hessz[ic][2] += dedphi*dziczic + d2edphi2*dphidzic*dphidzic + 2.0*dzic*ddrdz + d2drdzz;
            hess.hessx[ia][0] += dedphi*dxiaxic + dxia*ddrdx + d2edphi2*dphidxic*dphidxia;
            hess.hessy[ia][0] += dedphi*dxiayic + dxia*ddrdy + d2edphi2*dphidyic*dphidxia;
            hess.hessz[ia][0] += dedphi*dxiazic + dxia*ddrdz + d2edphi2*dphidzic*dphidxia;
            hess.hessx[ia][1] += dedphi*dyiaxic + dyia*ddrdx + d2edphi2*dphidxic*dphidyia;
            hess.hessy[ia][1] += dedphi*dyiayic + dyia*ddrdy + d2edphi2*dphidyic*dphidyia;
            hess.hessz[ia][1] += dedphi*dyiazic + dyia*ddrdz + d2edphi2*dphidzic*dphidyia;
            hess.hessx[ia][2] += dedphi*dziaxic + dzia*ddrdx + d2edphi2*dphidxic*dphidzia;
            hess.hessy[ia][2] += dedphi*dziayic + dzia*ddrdy + d2edphi2*dphidyic*dphidzia;
            hess.hessz[ia][2] += dedphi*dziazic + dzia*ddrdz + d2edphi2*dphidzic*dphidzia;
            hess.hessx[ib][0] += dedphi*dxibxic + d2edphi2*dphidxic*dphidxib - (dxic-dxib)*ddrdx - d2drdxx;
            hess.hessy[ib][0] += dedphi*dxibyic + d2edphi2*dphidyic*dphidxib - dyic*ddrdx + dxib*ddrdy - d2drdxy;
            hess.hessz[ib][0] += dedphi*dxibzic + d2edphi2*dphidzic*dphidxib - dzic*ddrdx + dxib*ddrdz - d2drdxz;
            hess.hessx[ib][1] += dedphi*dyibxic + d2edphi2*dphidxic*dphidyib - dxic*ddrdy + dyib*ddrdx - d2drdxy;
            hess.hessy[ib][1] += dedphi*dyibyic + d2edphi2*dphidyic*dphidyib - (dyic-dyib)*ddrdy - d2drdyy;
            hess.hessz[ib][1] += dedphi*dyibzic + d2edphi2*dphidzic*dphidyib - dzic*ddrdy + dyib*ddrdz - d2drdyz;
            hess.hessx[ib][2] += dedphi*dzibxic + d2edphi2*dphidxic*dphidzib - dxic*ddrdz + dzib*ddrdx - d2drdxz;
            hess.hessy[ib][2] += dedphi*dzibyic + d2edphi2*dphidyic*dphidzib - dyic*ddrdz + dzib*ddrdy - d2drdyz;
            hess.hessz[ib][2] += dedphi*dzibzic + d2edphi2*dphidzic*dphidzib - (dzic-dzib)*ddrdz - d2drdzz;
            hess.hessx[id][0] += dedphi*dxicxid + dxid*ddrdx + d2edphi2*dphidxic*dphidxid;
            hess.hessy[id][0] += dedphi*dyicxid + dxid*ddrdy + d2edphi2*dphidyic*dphidxid;
            hess.hessz[id][0] += dedphi*dzicxid + dxid*ddrdz + d2edphi2*dphidzic*dphidxid;
            hess.hessx[id][1] += dedphi*dxicyid + dyid*ddrdx + d2edphi2*dphidxic*dphidyid;
            hess.hessy[id][1] += dedphi*dyicyid + dyid*ddrdy + d2edphi2*dphidyic*dphidyid;
            hess.hessz[id][1] += dedphi*dzicyid + dyid*ddrdz + d2edphi2*dphidzic*dphidyid;
            hess.hessx[id][2] += dedphi*dxiczid + dzid*ddrdx + d2edphi2*dphidxic*dphidzid;
            hess.hessy[id][2] += dedphi*dyiczid + dzid*ddrdy + d2edphi2*dphidyic*dphidzid;
            hess.hessz[id][2] += dedphi*dziczid + dzid*ddrdz + d2edphi2*dphidzic*dphidzid;
         }else if (iatom == id)
         {
            hess.hessx[id][0] += dedphi*dxidxid + d2edphi2*dphidxid*dphidxid;
            hess.hessy[id][0] += dedphi*dxidyid + d2edphi2*dphidxid*dphidyid;
            hess.hessz[id][0] += dedphi*dxidzid + d2edphi2*dphidxid*dphidzid;
            hess.hessx[id][1] += dedphi*dxidyid + d2edphi2*dphidxid*dphidyid;
            hess.hessy[id][1] += dedphi*dyidyid + d2edphi2*dphidyid*dphidyid;
            hess.hessz[id][1] += dedphi*dyidzid + d2edphi2*dphidyid*dphidzid;
            hess.hessx[id][2] += dedphi*dxidzid + d2edphi2*dphidxid*dphidzid;
            hess.hessy[id][2] += dedphi*dyidzid + d2edphi2*dphidyid*dphidzid;
            hess.hessz[id][2] += dedphi*dzidzid + d2edphi2*dphidzid*dphidzid;
            hess.hessx[ia][0] += dedphi*dxiaxid + d2edphi2*dphidxid*dphidxia;
            hess.hessy[ia][0] += dedphi*dxiayid + d2edphi2*dphidyid*dphidxia;
            hess.hessz[ia][0] += dedphi*dxiazid + d2edphi2*dphidzid*dphidxia;
            hess.hessx[ia][1] += dedphi*dyiaxid + d2edphi2*dphidxid*dphidyia;
            hess.hessy[ia][1] += dedphi*dyiayid + d2edphi2*dphidyid*dphidyia;
            hess.hessz[ia][1] += dedphi*dyiazid + d2edphi2*dphidzid*dphidyia;
            hess.hessx[ia][2] += dedphi*dziaxid + d2edphi2*dphidxid*dphidzia;
            hess.hessy[ia][2] += dedphi*dziayid + d2edphi2*dphidyid*dphidzia;
            hess.hessz[ia][2] += dedphi*dziazid + d2edphi2*dphidzid*dphidzia;
            hess.hessx[ib][0] += dedphi*dxibxid - dxid*ddrdx + d2edphi2*dphidxid*dphidxib;
            hess.hessy[ib][0] += dedphi*dxibyid - dyid*ddrdx + d2edphi2*dphidyid*dphidxib;
            hess.hessz[ib][0] += dedphi*dxibzid - dzid*ddrdx + d2edphi2*dphidzid*dphidxib;
            hess.hessx[ib][1] += dedphi*dyibxid - dxid*ddrdy + d2edphi2*dphidxid*dphidyib;
            hess.hessy[ib][1] += dedphi*dyibyid - dyid*ddrdy + d2edphi2*dphidyid*dphidyib;
            hess.hessz[ib][1] += dedphi*dyibzid - dzid*ddrdy + d2edphi2*dphidzid*dphidyib;
            hess.hessx[ib][2] += dedphi*dzibxid - dxid*ddrdz + d2edphi2*dphidxid*dphidzib;
            hess.hessy[ib][2] += dedphi*dzibyid - dyid*ddrdz + d2edphi2*dphidyid*dphidzib;
            hess.hessz[ib][2] += dedphi*dzibzid - dzid*ddrdz + d2edphi2*dphidzid*dphidzib;
            hess.hessx[ic][0] += dedphi*dxicxid + dxid*ddrdx + d2edphi2*dphidxid*dphidxic;
            hess.hessy[ic][0] += dedphi*dxicyid + dyid*ddrdx + d2edphi2*dphidyid*dphidxic;
            hess.hessz[ic][0] += dedphi*dxiczid + dzid*ddrdx + d2edphi2*dphidzid*dphidxic;
            hess.hessx[ic][1] += dedphi*dyicxid + dxid*ddrdy + d2edphi2*dphidxid*dphidyic;
            hess.hessy[ic][1] += dedphi*dyicyid + dyid*ddrdy + d2edphi2*dphidyid*dphidyic;
            hess.hessz[ic][1] += dedphi*dyiczid + dzid*ddrdy + d2edphi2*dphidzid*dphidyic;
            hess.hessx[ic][2] += dedphi*dzicxid + dxid*ddrdz + d2edphi2*dphidxid*dphidzic;
            hess.hessy[ic][2] += dedphi*dzicyid + dyid*ddrdz + d2edphi2*dphidyid*dphidzic;
            hess.hessz[ic][2] += dedphi*dziczid + dzid*ddrdz + d2edphi2*dphidzid*dphidzic;
         }
            }
          }
      }
}
