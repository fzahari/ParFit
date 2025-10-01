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
#include "derivs.h"
#include "hess.h"
#include "fix.h"
#include "units.h"
#include "minim_values.h"
#include "minim_control.h"

void etorsion()
{
    int i,ia, ib, ic, id;
    double e, rt2, ru2, rtru;
    double xt,yt,zt,xu,yu,zu;
    double v1,v2,v3,v4,v5,v6;
    double s1,s2,s3,s4,s5,s6;
    double cosine,cosine2,cosine3;
    double cosine4,cosine5,cosine6;
    double phi1,phi2,phi3;
    double phi4,phi5,phi6;
    double xia,yia,zia,xib,yib,zib;
    double xic,yic,zic,xid,yid,zid;
    double xba,yba,zba,xcb,ycb,zcb;
    double xdc,ydc,zdc;
    double xtu,ytu,ztu,rcb, angle,rdihed;
    double width,sine,dt2;
    int curang;
    
    energies.etor = 0.0;
    if (minim_values.iprint)
    {
      fprintf(pcmoutfile,"\nTorsion Terms\n");
      fprintf(pcmoutfile,"        At1      At2      At3      At4       Types          Angle     V1      V2        V3         Etor\n");
    }

    width = 0.0;
    
    for (i=0; i < torsions.ntor; i++)
    {
        ia = torsions.i14[i][0];
        ib = torsions.i14[i][1];
        ic = torsions.i14[i][2];
        id = torsions.i14[i][3];
        if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use )
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

               v1 = torsions.v1[i];
               s1 = torsions.ph1[i];
               v2 = torsions.v2[i];
               s2 = torsions.ph2[i];
               v3 = torsions.v3[i];
               s3 = torsions.ph3[i];
               v4 = torsions.v4[i];
               s4 = torsions.ph4[i];
               v5 = torsions.v5[i];
               s5 = torsions.ph5[i];
               v6 = torsions.v6[i];
               s6 = torsions.ph6[i];

//     compute the powers of the cosine of this angle

               cosine2 = cosine * cosine;
               cosine3 = cosine2 * cosine;
               cosine4 = cosine2 * cosine2;
               cosine5 = cosine3 * cosine2;
               cosine6 = cosine3 * cosine3;

               phi1 = 1.0 + s1*cosine;
               phi2 = 1.0 + s2*(2.0*cosine2 - 1.0);
               phi3 = 1.0 + s3*(4.0*cosine3 - 3.0*cosine);
               phi4 = 1.0 + s4*(8.0*(cosine4-cosine2) + 1.0);
               phi5 = 1.0 + s5*(16.0*cosine5 - 20.0*cosine3 + 5.0*cosine);
               phi6 = 1.0 + s6*(32.0*cosine6 - 48.0*cosine4 + 18.0*cosine2 - 1.0);

               e = units.torsunit * (v1*phi1 + v2*phi2 + v3*phi3
                                   + v4*phi4 + v5*phi5 + v6*phi6);
               energies.etor += e;
               atom[ia].energy += e;
               atom[ib].energy += e;
               atom[ic].energy += e;
               atom[id].energy += e;
               if(minim_values.iprint)
               {
                   rdihed = dihdrl(ia,ib,ic,id);
                   fprintf(pcmoutfile,"Tor:  %2s(%-3d)- %2s(%-3d)- %2s(%-3d)- %2s(%-3d)    %d %d %d %d     %8.2f    %-8.3f %-8.3f %-8.3f = %-8.4f\n",
                   atom[ia].name,ia,atom[ib].name,ib,atom[ic].name,ic,atom[id].name,id, atom[ia].type, atom[ib].type,
                    atom[ic].type, atom[id].type, rdihed,v1,v2,v3,e);
               }
            }
        }
     }
// fixed torsions
    for (i=0; i < fxtor.nfxtor; i++)
    {
        ia = fxtor.iatom[i][0];
        ib = fxtor.iatom[i][1];
        ic = fxtor.iatom[i][2];
        id = fxtor.iatom[i][3];
        if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use )
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
               cosine = MinFun(1.0, MaxFun(-1.0,cosine));
               angle = radian*acos(cosine);
               if (sine < 0.0) angle = -angle;
               
               curang = fxtor.curr_ang[i];
               if (curang > 360)
                  curang %= 360;
                              
               v1 = 2.0;
               s1 = 1.0;
               phi1 = s1*(curang - angle);
               if (phi1 > 180.0)
                  phi1 -= 360.0;
               else if (phi1 < -180.0)
                  phi1 += 360.0;
               dt2 = phi1*phi1;   
               e = units.torsunit * v1*dt2;
               energies.etor += e;
               atom[ia].energy += e;
               atom[ib].energy += e;
               atom[ic].energy += e;
               atom[id].energy += e;
               if(minim_values.iprint)
               {
                   rdihed = dihdrl(ia,ib,ic,id);
                   fprintf(pcmoutfile,"FxTor:  %2s(%-3d)- %2s(%-3d)- %2s(%-3d)- %2s(%-3d)  %8.2f    %-8.3f                   = %-8.4f\n",
                     atom[ia].name,ia,atom[ib].name,ib,atom[ic].name,ic,atom[id].name,id, rdihed,v1,e);
               }
            }
        }
     }


}

void etorsion1()
{
    int i,ia, ib, ic, id;
    int curang;
    double angle;
    double e, rt2, ru2, rtru,dedphi,sine,dt2;
    double xt,yt,zt,xu,yu,zu;
    double v1,v2,v3,v4,v5,v6;
    double s1,s2,s3,s4,s5,s6;
    double cosine,cosine2,cosine3;
    double cosine4,cosine5,cosine6;
    double phi1,phi2,phi3;
    double phi4,phi5,phi6;
    double xia,yia,zia,xib,yib,zib;
    double xic,yic,zic,xid,yid,zid;
    double xba,yba,zba,xcb,ycb,zcb;
    double xdc,ydc,zdc;
    double xca,yca,zca,xdb,ydb,zdb;
    double xtu,ytu,ztu,rcb;
    double dphi1,dphi2,dphi3;
    double dphi4,dphi5,dphi6;
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
    double dedxt,dedyt,dedzt, dedxu,dedyu,dedzu;
    double width;
    
    energies.etor = 0.0;
      for (i=0; i <= natom; i++)
      {
          deriv.detor[i][0] = 0.0;
          deriv.detor[i][1] = 0.0;
          deriv.detor[i][2] = 0.0;
      }

    width = 0.0;
    
    for (i=0; i < torsions.ntor; i++)
    {
        ia = torsions.i14[i][0];
        ib = torsions.i14[i][1];
        ic = torsions.i14[i][2];
        id = torsions.i14[i][3];
        if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use )
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

               v1 = torsions.v1[i];
               s1 = torsions.ph1[i];
               v2 = torsions.v2[i];
               s2 = torsions.ph2[i];
               v3 = torsions.v3[i];
               s3 = torsions.ph3[i];
               v4 = torsions.v4[i];
               s4 = torsions.ph4[i];
               v5 = torsions.v5[i];
               s5 = torsions.ph5[i];
               v6 = torsions.v6[i];
               s6 = torsions.ph6[i];
//     compute the powers of the cosine of this angle

               cosine2 = cosine * cosine;
               cosine3 = cosine2 * cosine;
               cosine4 = cosine2 * cosine2;
               cosine5 = cosine3 * cosine2;
               cosine6 = cosine3 * cosine3;

               phi1 = 1.0 + s1*cosine;
               phi2 = 1.0 + s2*(2.0*cosine2 - 1.0);
               phi3 = 1.0 + s3*(4.0*cosine3 - 3.0*cosine);
               phi4 = 1.0 + s4*(8.0*(cosine4-cosine2) + 1.0);
               phi5 = 1.0 + s5*(16.0*cosine5 - 20.0*cosine3 + 5.0*cosine);
               phi6 = 1.0 + s6*(32.0*cosine6 - 48.0*cosine4 + 18.0*cosine2 - 1.0);

               dphi1 = s1;
               dphi2 = s2 * (4.0*cosine);
               dphi3 = s3 * (12.0*cosine2 - 3.0);
               dphi4 = s4 * (32.0*cosine3 - 16.0*cosine);
               dphi5 = s5 * (80.0*cosine4 - 60.0*cosine2 + 5.0);
               dphi6 = s6 * (192.0*(cosine5-cosine3) + 36.0*cosine);

               e = units.torsunit * (v1*phi1 + v2*phi2 + v3*phi3
                                   + v4*phi4 + v5*phi5 + v6*phi6);
               dedphi = -sine * units.torsunit
                       * (v1*dphi1+v2*dphi2+v3*dphi3+v4*dphi4+v5*dphi5+v6*dphi6);

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
               dedxia = dedphi * dphidxia;
               dedyia = dedphi * dphidyia;
               dedzia = dedphi * dphidzia;
               dedxib = dedphi * dphidxib;
               dedyib = dedphi * dphidyib;
               dedzib = dedphi * dphidzib;
               dedxic = dedphi * dphidxic;
               dedyic = dedphi * dphidyic;
               dedzic = dedphi * dphidzic;
               dedxid = dedphi * dphidxid;
               dedyid = dedphi * dphidyid;
               dedzid = dedphi * dphidzid;

               energies.etor += e;
               deriv.detor[ia][0] += dedxia;
               deriv.detor[ia][1] += dedyia;
               deriv.detor[ia][2] += dedzia;
               
               deriv.detor[ib][0] += dedxib;
               deriv.detor[ib][1] += dedyib;
               deriv.detor[ib][2] += dedzib;

               deriv.detor[ic][0] += dedxic;
               deriv.detor[ic][1] += dedyic;
               deriv.detor[ic][2] += dedzic;

               deriv.detor[id][0] += dedxid;
               deriv.detor[id][1] += dedyid;
               deriv.detor[id][2] += dedzid;
               
               virial.virx = virial.virx - xba*dedxia + xdc*dedxid
                          + xcb*(dedxic+dedxid);
               virial.viry = virial.viry - yba*dedyia + ydc*dedyid
                          + ycb*(dedyic+dedyid);
               virial.virz = virial.virz - zba*dedzia + zdc*dedzid
                          + zcb*(dedzic+dedzid);
            }
        }
     }
// fixed torsions
    for (i=0; i < fxtor.nfxtor; i++)
    {
        ia = fxtor.iatom[i][0];
        ib = fxtor.iatom[i][1];
        ic = fxtor.iatom[i][2];
        id = fxtor.iatom[i][3];
        if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use )
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
               cosine = MinFun(1.0, MaxFun(-1.0,cosine));
               angle = radian*acos(cosine);
               if (sine < 0.0) angle = -angle;
               curang = fxtor.curr_ang[i];
               if (curang > 360)
                  curang %= 360;
                  
               v1 = 2.0;
               s1 = 1.0;
               phi1 = s1*(angle - curang);
               if (phi1 > 180.0)
                  phi1 -= 360.0;
               else if (phi1 < -180.0)
                  phi1 += 360.0;
               dt2 = phi1*phi1;   
               e = units.torsunit * v1*dt2;
               dedphi = 2.0*radian*units.torsunit*v1*phi1;

               xca = xic - xia;
               yca = yic - yia;
               zca = zic - zia;
               xdb = xid - xib;
               ydb = yid - yib;
               zdb = zid - zib;
               dedxt = dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
               dedyt = dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
               dedzt = dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
               dedxu = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
               dedyu = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
               dedzu = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);
               
               dedxia = zcb*dedyt - ycb*dedzt;
               dedyia = xcb*dedzt - zcb*dedxt;
               dedzia = ycb*dedxt - xcb*dedyt;
               dedxib = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
               dedyib = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
               dedzib = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
               dedxic = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
               dedyic = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
               dedzic = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
               dedxid = zcb*dedyu - ycb*dedzu;
               dedyid = xcb*dedzu - zcb*dedxu;
               dedzid = ycb*dedxu - xcb*dedyu;
               
               energies.etor += e;
               deriv.detor[ia][0] += dedxia;
               deriv.detor[ia][1] += dedyia;
               deriv.detor[ia][2] += dedzia;
               
               deriv.detor[ib][0] += dedxib;
               deriv.detor[ib][1] += dedyib;
               deriv.detor[ib][2] += dedzib;

               deriv.detor[ic][0] += dedxic;
               deriv.detor[ic][1] += dedyic;
               deriv.detor[ic][2] += dedzic;

               deriv.detor[id][0] += dedxid;
               deriv.detor[id][1] += dedyid;
               deriv.detor[id][2] += dedzid;
            }
        }
     }


}

void etorsion2(iatom)
{
    int i,ia,ib,ic,id;
    int curang;
    double agl, cosm, sinm;
    double dot, denom, denom2, denom3;
    double dedphi,d2edphi2,sine;
    double term1, term2, term3;
    double v1,v2,v3,v4,v5,v6;
    double s1,s2,s3,s4,s5,s6;
    double cosine,cosine2,cosine3;
    double cosine4,cosine5,cosine6;
    double xia,yia,zia,xib,yib,zib;
    double xic,yic,zic,xid,yid,zid;
    double xba,yba,zba,xcb,ycb,zcb;
    double xdc,ydc,zdc;
    double xca,yca,zca,xdb,ydb,zdb;
    double xt,yt,zt,xu,yu,zu,xtu,ytu,ztu;
    double rt2,ru2,rtru,rcb, rt2ru2;
    double phi1,phi2,phi3;
    double phi4,phi5,phi6;
    double dphi1,dphi2,dphi3;
    double dphi4,dphi5,dphi6;
    double d2phi1,d2phi2,d2phi3;
    double d2phi4,d2phi5,d2phi6;
    double dphidxt,dphidyt,dphidzt;
    double dphidxu,dphidyu,dphidzu;
    double dphidxt1,dphidyt1,dphidzt1;
    double dphidxu1,dphidyu1,dphidzu1;
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
    double yuzab, zuydb, ytzba, rtru32, rtru322;
    double zdcyu, xdczu, ydcxu, ycazt, zcaxt, xcayt;
    double width;

    width = 0.0;

    for (i=0; i < torsions.ntor; i++)
    {
        ia = torsions.i14[i][0];
        ib = torsions.i14[i][1];
        ic = torsions.i14[i][2];
        id = torsions.i14[i][3];
        if (ia == iatom || ib == iatom || ic == iatom || id == iatom)
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

               v1 = torsions.v1[i];
               s1 = torsions.ph1[i];
               v2 = torsions.v2[i];
               s2 = torsions.ph2[i];
               v3 = torsions.v3[i];
               s3 = torsions.ph3[i];
               v4 = torsions.v4[i];
               s4 = torsions.ph4[i];
               v5 = torsions.v5[i];
               s5 = torsions.ph5[i];
               v6 = torsions.v6[i];
               s6 = torsions.ph6[i];

//     compute the powers of the cosine of this angle

               cosine2 = cosine * cosine;
               cosine3 = cosine2 * cosine;
               cosine4 = cosine2 * cosine2;
               cosine5 = cosine3 * cosine2;
               cosine6 = cosine3 * cosine3;

               phi1 = 1.0 + s1*cosine;
               phi2 = 1.0 + s2*(2.0*cosine2 - 1.0);
               phi3 = 1.0 + s3*(4.0*cosine3 - 3.0*cosine);
               phi4 = 1.0 + s4*(8.0*(cosine4-cosine2) + 1.0);
               phi5 = 1.0 + s5*(16.0*cosine5 - 20.0*cosine3 + 5.0*cosine);
               phi6 = 1.0 + s6*(32.0*cosine6 - 48.0*cosine4 + 18.0*cosine2 - 1.0);

               dphi1 = s1;
               dphi2 = s2 * (4.0*cosine);
               dphi3 = s3 * (12.0*cosine2 - 3.0);
               dphi4 = s4 * (32.0*cosine3 - 16.0*cosine);
               dphi5 = s5 * (80.0*cosine4 - 60.0*cosine2 + 5.0);
               dphi6 = s6 * (192.0*(cosine5-cosine3) + 36.0*cosine);
            
               d2phi1 = -s1 * cosine;
               d2phi2 = -s2 * (8.0*cosine2 - 4.0);
               d2phi3 = -s3 * (36.0*cosine3 - 27.0*cosine);
               d2phi4 = -s4 * (128.0*(cosine4-cosine2) + 16.0);
               d2phi5 = -s5 * (400.0*cosine5 - 500.0*cosine + 125.0*cosine);
               d2phi6 = -s6 * (1152.0*cosine6 - 1728.0*cosine2 + 648.0*cosine2);

               dedphi = -sine * units.torsunit * (v1*dphi1 + v2*dphi2 + v3*dphi3
                             + v4*dphi4 + v5*dphi5 + v6*dphi6);
               d2edphi2 = units.torsunit * (v1*d2phi1 + v2*d2phi2 + v3*d2phi3
                             + v4*d2phi4 + v5*d2phi5 + v6*d2phi6);

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
               dyibyic = -ycb*dphidyib/(rcb*rcb)- (zca*(xba*ycb+zt)-xca*(zba*ycb-xt))/rcbt2
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
            hess.hessx[ib][0] += dedphi*dxiaxib + d2edphi2*dphidxia*dphidxib;
            hess.hessy[ib][0] += dedphi*dyiaxib + d2edphi2*dphidyia*dphidxib;
            hess.hessz[ib][0] += dedphi*dziaxib + d2edphi2*dphidzia*dphidxib;
            hess.hessx[ib][1] += dedphi*dxiayib + d2edphi2*dphidxia*dphidyib;
            hess.hessy[ib][1] += dedphi*dyiayib + d2edphi2*dphidyia*dphidyib;
            hess.hessz[ib][1] += dedphi*dziayib + d2edphi2*dphidzia*dphidyib;
            hess.hessx[ib][2] += dedphi*dxiazib + d2edphi2*dphidxia*dphidzib;
            hess.hessy[ib][2] += dedphi*dyiazib + d2edphi2*dphidyia*dphidzib;
            hess.hessz[ib][2] += dedphi*dziazib + d2edphi2*dphidzia*dphidzib;
            hess.hessx[ic][0] += dedphi*dxiaxic + d2edphi2*dphidxia*dphidxic;
            hess.hessy[ic][0] += dedphi*dyiaxic + d2edphi2*dphidyia*dphidxic;
            hess.hessz[ic][0] += dedphi*dziaxic + d2edphi2*dphidzia*dphidxic;
            hess.hessx[ic][1] += dedphi*dxiayic + d2edphi2*dphidxia*dphidyic;
            hess.hessy[ic][1] += dedphi*dyiayic + d2edphi2*dphidyia*dphidyic;
            hess.hessz[ic][1] += dedphi*dziayic + d2edphi2*dphidzia*dphidyic;
            hess.hessx[ic][2] += dedphi*dxiazic + d2edphi2*dphidxia*dphidzic;
            hess.hessy[ic][2] += dedphi*dyiazic + d2edphi2*dphidyia*dphidzic;
            hess.hessz[ic][2] += dedphi*dziazic + d2edphi2*dphidzia*dphidzic;
            hess.hessx[id][0] += dedphi*dxiaxid + d2edphi2*dphidxia*dphidxid;
            hess.hessy[id][0] += dedphi*dyiaxid + d2edphi2*dphidyia*dphidxid;
            hess.hessz[id][0] += dedphi*dziaxid + d2edphi2*dphidzia*dphidxid;
            hess.hessx[id][1] += dedphi*dxiayid + d2edphi2*dphidxia*dphidyid;
            hess.hessy[id][1] += dedphi*dyiayid + d2edphi2*dphidyia*dphidyid;
            hess.hessz[id][1] += dedphi*dziayid + d2edphi2*dphidzia*dphidyid;
            hess.hessx[id][2] += dedphi*dxiazid + d2edphi2*dphidxia*dphidzid;
            hess.hessy[id][2] += dedphi*dyiazid + d2edphi2*dphidyia*dphidzid;
            hess.hessz[id][2] += dedphi*dziazid + d2edphi2*dphidzia*dphidzid;
         }else if (iatom == ib)
         {
            hess.hessx[ib][0] += dedphi*dxibxib + d2edphi2*dphidxib*dphidxib;
            hess.hessy[ib][0] += dedphi*dxibyib + d2edphi2*dphidxib*dphidyib;
            hess.hessz[ib][0] += dedphi*dxibzib + d2edphi2*dphidxib*dphidzib;
            hess.hessx[ib][1] += dedphi*dxibyib + d2edphi2*dphidxib*dphidyib;
            hess.hessy[ib][1] += dedphi*dyibyib + d2edphi2*dphidyib*dphidyib;
            hess.hessz[ib][1] += dedphi*dyibzib + d2edphi2*dphidyib*dphidzib;
            hess.hessx[ib][2] += dedphi*dxibzib + d2edphi2*dphidxib*dphidzib;
            hess.hessy[ib][2] += dedphi*dyibzib + d2edphi2*dphidyib*dphidzib;
            hess.hessz[ib][2] += dedphi*dzibzib + d2edphi2*dphidzib*dphidzib;
            hess.hessx[ia][0] += dedphi*dxiaxib + d2edphi2*dphidxib*dphidxia;
            hess.hessy[ia][0] += dedphi*dxiayib + d2edphi2*dphidyib*dphidxia;
            hess.hessz[ia][0] += dedphi*dxiazib + d2edphi2*dphidzib*dphidxia;
            hess.hessx[ia][1] += dedphi*dyiaxib + d2edphi2*dphidxib*dphidyia;
            hess.hessy[ia][1] += dedphi*dyiayib + d2edphi2*dphidyib*dphidyia;
            hess.hessz[ia][1] += dedphi*dyiazib + d2edphi2*dphidzib*dphidyia;
            hess.hessx[ia][2] += dedphi*dziaxib + d2edphi2*dphidxib*dphidzia;
            hess.hessy[ia][2] += dedphi*dziayib + d2edphi2*dphidyib*dphidzia;
            hess.hessz[ia][2] += dedphi*dziazib + d2edphi2*dphidzib*dphidzia;
            hess.hessx[ic][0] += dedphi*dxibxic + d2edphi2*dphidxib*dphidxic;
            hess.hessy[ic][0] += dedphi*dyibxic + d2edphi2*dphidyib*dphidxic;
            hess.hessz[ic][0] += dedphi*dzibxic + d2edphi2*dphidzib*dphidxic;
            hess.hessx[ic][1] += dedphi*dxibyic + d2edphi2*dphidxib*dphidyic;
            hess.hessy[ic][1] += dedphi*dyibyic + d2edphi2*dphidyib*dphidyic;
            hess.hessz[ic][1] += dedphi*dzibyic + d2edphi2*dphidzib*dphidyic;
            hess.hessx[ic][2] += dedphi*dxibzic + d2edphi2*dphidxib*dphidzic;
            hess.hessy[ic][2] += dedphi*dyibzic + d2edphi2*dphidyib*dphidzic;
            hess.hessz[ic][2] += dedphi*dzibzic + d2edphi2*dphidzib*dphidzic;
            hess.hessx[id][0] += dedphi*dxibxid + d2edphi2*dphidxib*dphidxid;
            hess.hessy[id][0] += dedphi*dyibxid + d2edphi2*dphidyib*dphidxid;
            hess.hessz[id][0] += dedphi*dzibxid + d2edphi2*dphidzib*dphidxid;
            hess.hessx[id][1] += dedphi*dxibyid + d2edphi2*dphidxib*dphidyid;
            hess.hessy[id][1] += dedphi*dyibyid + d2edphi2*dphidyib*dphidyid;
            hess.hessz[id][1] += dedphi*dzibyid + d2edphi2*dphidzib*dphidyid;
            hess.hessx[id][2] += dedphi*dxibzid + d2edphi2*dphidxib*dphidzid;
            hess.hessy[id][2] += dedphi*dyibzid + d2edphi2*dphidyib*dphidzid;
            hess.hessz[id][2] += dedphi*dzibzid + d2edphi2*dphidzib*dphidzid;
         }else if (iatom == ic)
         {
            hess.hessx[ic][0] += dedphi*dxicxic + d2edphi2*dphidxic*dphidxic;
            hess.hessy[ic][0] += dedphi*dxicyic + d2edphi2*dphidxic*dphidyic;
            hess.hessz[ic][0] += dedphi*dxiczic + d2edphi2*dphidxic*dphidzic;
            hess.hessx[ic][1] += dedphi*dxicyic + d2edphi2*dphidxic*dphidyic;
            hess.hessy[ic][1] += dedphi*dyicyic + d2edphi2*dphidyic*dphidyic;
            hess.hessz[ic][1] += dedphi*dyiczic + d2edphi2*dphidyic*dphidzic;
            hess.hessx[ic][2] += dedphi*dxiczic + d2edphi2*dphidxic*dphidzic;
            hess.hessy[ic][2] += dedphi*dyiczic + d2edphi2*dphidyic*dphidzic;
            hess.hessz[ic][2] += dedphi*dziczic + d2edphi2*dphidzic*dphidzic;
            hess.hessx[ia][0] += dedphi*dxiaxic + d2edphi2*dphidxic*dphidxia;
            hess.hessy[ia][0] += dedphi*dxiayic + d2edphi2*dphidyic*dphidxia;
            hess.hessz[ia][0] += dedphi*dxiazic + d2edphi2*dphidzic*dphidxia;
            hess.hessx[ia][1] += dedphi*dyiaxic + d2edphi2*dphidxic*dphidyia;
            hess.hessy[ia][1] += dedphi*dyiayic + d2edphi2*dphidyic*dphidyia;
            hess.hessz[ia][1] += dedphi*dyiazic + d2edphi2*dphidzic*dphidyia;
            hess.hessx[ia][2] += dedphi*dziaxic + d2edphi2*dphidxic*dphidzia;
            hess.hessy[ia][2] += dedphi*dziayic + d2edphi2*dphidyic*dphidzia;
            hess.hessz[ia][2] += dedphi*dziazic + d2edphi2*dphidzic*dphidzia;
            hess.hessx[ib][0] += dedphi*dxibxic + d2edphi2*dphidxic*dphidxib;
            hess.hessy[ib][0] += dedphi*dxibyic + d2edphi2*dphidyic*dphidxib;
            hess.hessz[ib][0] += dedphi*dxibzic + d2edphi2*dphidzic*dphidxib;
            hess.hessx[ib][1] += dedphi*dyibxic + d2edphi2*dphidxic*dphidyib;
            hess.hessy[ib][1] += dedphi*dyibyic + d2edphi2*dphidyic*dphidyib;
            hess.hessz[ib][1] += dedphi*dyibzic + d2edphi2*dphidzic*dphidyib;
            hess.hessx[ib][2] += dedphi*dzibxic + d2edphi2*dphidxic*dphidzib;
            hess.hessy[ib][2] += dedphi*dzibyic + d2edphi2*dphidyic*dphidzib;
            hess.hessz[ib][2] += dedphi*dzibzic + d2edphi2*dphidzic*dphidzib;
            hess.hessx[id][0] += dedphi*dxicxid + d2edphi2*dphidxic*dphidxid;
            hess.hessy[id][0] += dedphi*dyicxid + d2edphi2*dphidyic*dphidxid;
            hess.hessz[id][0] += dedphi*dzicxid + d2edphi2*dphidzic*dphidxid;
            hess.hessx[id][1] += dedphi*dxicyid + d2edphi2*dphidxic*dphidyid;
            hess.hessy[id][1] += dedphi*dyicyid + d2edphi2*dphidyic*dphidyid;
            hess.hessz[id][1] += dedphi*dzicyid + d2edphi2*dphidzic*dphidyid;
            hess.hessx[id][2] += dedphi*dxiczid + d2edphi2*dphidxic*dphidzid;
            hess.hessy[id][2] += dedphi*dyiczid + d2edphi2*dphidyic*dphidzid;
            hess.hessz[id][2] += dedphi*dziczid + d2edphi2*dphidzic*dphidzid;
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
            hess.hessx[ib][0] += dedphi*dxibxid + d2edphi2*dphidxid*dphidxib;
            hess.hessy[ib][0] += dedphi*dxibyid + d2edphi2*dphidyid*dphidxib;
            hess.hessz[ib][0] += dedphi*dxibzid + d2edphi2*dphidzid*dphidxib;
            hess.hessx[ib][1] += dedphi*dyibxid + d2edphi2*dphidxid*dphidyib;
            hess.hessy[ib][1] += dedphi*dyibyid + d2edphi2*dphidyid*dphidyib;
            hess.hessz[ib][1] += dedphi*dyibzid + d2edphi2*dphidzid*dphidyib;
            hess.hessx[ib][2] += dedphi*dzibxid + d2edphi2*dphidxid*dphidzib;
            hess.hessy[ib][2] += dedphi*dzibyid + d2edphi2*dphidyid*dphidzib;
            hess.hessz[ib][2] += dedphi*dzibzid + d2edphi2*dphidzid*dphidzib;
            hess.hessx[ic][0] += dedphi*dxicxid + d2edphi2*dphidxid*dphidxic;
            hess.hessy[ic][0] += dedphi*dxicyid + d2edphi2*dphidyid*dphidxic;
            hess.hessz[ic][0] += dedphi*dxiczid + d2edphi2*dphidzid*dphidxic;
            hess.hessx[ic][1] += dedphi*dyicxid + d2edphi2*dphidxid*dphidyic;
            hess.hessy[ic][1] += dedphi*dyicyid + d2edphi2*dphidyid*dphidyic;
            hess.hessz[ic][1] += dedphi*dyiczid + d2edphi2*dphidzid*dphidyic;
            hess.hessx[ic][2] += dedphi*dzicxid + d2edphi2*dphidxid*dphidzic;
            hess.hessy[ic][2] += dedphi*dzicyid + d2edphi2*dphidyid*dphidzic;
            hess.hessz[ic][2] += dedphi*dziczid + d2edphi2*dphidzid*dphidzic;
         }
            }
        }
    }

// fixed torsion
    for (i=0; i < fxtor.nfxtor; i++)
    {
        ia = fxtor.iatom[i][0];
        ib = fxtor.iatom[i][1];
        ic = fxtor.iatom[i][2];
        id = fxtor.iatom[i][3];
        if (ia == iatom || ib == iatom || ic == iatom || id == iatom)
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
            rt2ru2 = rt2*ru2;
            rtru = sqrt(rt2 * ru2);
            if (rtru != 0.0)
            {
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               cosine = (xt*xu + yt*yu + zt*zu) / rtru;

               curang = fxtor.curr_ang[i];
               if (curang > 360)
                  curang %= 360;
                  
               if (curang > 180)
                  agl = -(360.0 - (double)curang -radian*acos(cosine));
               else if (fxtor.curr_ang[i] < 0)
                  agl = -(double)curang -radian*acos(cosine); 
               else
                  agl = (double)curang-radian*acos(cosine); 

               cosm = cos(agl/radian);
               sinm = sin(agl/radian);
               
               v1 = 2600.0;
               s1 = 1.0;
               phi1 = s1*(1.0-cosm);
 
               dedphi =   units.torsunit*v1*sinm;
               d2edphi2 = units.torsunit*v1*cosm;
               
               dot = xt*xu + yt*yu + zt*zu;
               denom2 = 1.0-((dot*dot)/(rt2*ru2));
               denom = sqrt(denom2);
//               if (denom < 0.00001)
//                  denom = 0.00001;
               denom3 = denom2*denom;
               
               xca = xic - xia;
               yca = yic - yia;
               zca = zic - zia;
               xdb = xid - xib;
               ydb = yid - yib;
               zdb = zid - zib;

               dphidxt = (yt*zcb - ycb*zt);
               dphidyt = (zt*xcb - zcb*xt);
               dphidzt = (xt*ycb - xcb*yt);
               dphidxu = (yu*zcb - ycb*zu);
               dphidyu = (zu*xcb - zcb*xu);
               dphidzu = (xu*ycb - xcb*yu);
               
               dphidxt1 = (yt*zdc - ydc*zt);
               dphidyt1 = (zt*xdc - zdc*xt);
               dphidzt1 = (xt*ydc - xdc*yt);
               dphidxu1 = (yu*zdc - ydc*zu);
               dphidyu1 = (zu*xdc - zdc*xu);
               dphidzu1 = (xu*ydc - xdc*yu);

               zdcyu = zdc*yu-ydc*zu;
               xdczu = xdc*zu-zdc*xu;
               ydcxu = ydc*xu-xdc*yu;
               ycazt = zt*yca-yt*zca;
               zcaxt = xt*zca-zt*xca;
               xcayt = yt*xca-xt*yca;
               rtru32 = rtru*rt2ru2;
               rtru322 =  2.0*rtru*rt2ru2;
               
// diagonal terms
               term1 =   dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2);
               term2 =  (3.0*4.0*dphidxt*dphidxt*dot)/(4.0*rt2*rt2*rtru) - (2.0*dphidxt*dphidxu)/(rt2*rtru)
                        -( (zcb*zcb+ycb*ycb)*dot)/(rt2*rtru);
               term3 = ( (2.0*dphidxt*dot*dot)/(rt2*rt2*ru2) - (2.0*dphidxu*dot)/rt2ru2 )*(term1);
               dxiaxia = (d2edphi2*term1*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//                                  
               term1 =   dphidyu/rtru - (dphidyt*dot*ru2)/(rtru*rt2ru2);
               
               term2 =   (3.0*4.0*dphidyt*dphidyt*dot*ru2*ru2)/(4.0*rt2ru2*rt2ru2*rtru) - (2.0*dphidyt*dphidyu*ru2)/(rt2ru2*rtru)
                        -( (zcb*zcb+xcb*xcb)*dot*ru2)/(rt2ru2*rtru);
               term3 = ( (2.0*dphidyt*dot*dot)/(rt2*rt2*ru2) - (2.0*dphidyu*dot)/rt2ru2 )*(term1);
               dyiayia = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =   dphidzu/rtru - (dphidzt*dot*ru2)/(rtru*rt2ru2);
               term2 =   (3.0*4.0*dphidzt*dphidzt*dot*ru2*ru2)/(4.0*rt2ru2*rt2ru2*rtru) - (2.0*dphidzt*dphidzu*ru2)/(rt2ru2*rtru)
                        -( (ycb*ycb+xcb*xcb)*dot*ru2)/(rt2ru2*rtru);
               term3 = ( (2.0*dphidzt*dot*dot)/(rt2*rt2*ru2) - (2.0*dphidzu*dot)/rt2ru2 )*(term1);
               dziazia = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru -(dot*(rt2*2.0*(zdcyu)+ru2*2.0*(ycazt)))/(rtru322);
               term2 =   -(2.0*(ydc*yca+zdc*zca))/rtru
                         -(dot*( 2.0*rt2*(ydc*ydc+zdc*zdc) + (8.0*(ycazt)*(zdcyu))
                         +(yca*yca+zca*zca)*ru2*2.0) )/(2.0*rt2ru2*rtru)
                         -( (yt*zdc-zt*ydc+zu*yca-yu*zca)*( rt2*2.0*(zdcyu)+ru2*2.0*(ycazt) ))/(rt2ru2*rtru)
                         +(3.0*dot*(rt2*2.0*(zdcyu)+ru2*2.0*(ycazt))
                         *( rt2*2.0*(zdcyu)+ru2*2.0*(ycazt)) )/(4.0*rt2ru2*rt2ru2*rtru);
               term3 = ( (2.0*(zdcyu)*dot*dot)/(rt2*ru2*ru2)
                        -(2.0*(yt*zdc-zt*ydc+zu*yca-yu*zca)*dot)/rt2ru2
                        +(2.0*(ycazt)*dot*dot)/(rt2*rt2*ru2) )*(term1);
               dxibxib = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = (zt*xdc-xt*zdc+xu*zca-zu*xca)/rtru - (dot*(rt2*2.0*(xdczu)+ru2*2.0*(zcaxt)))/(rtru322);
               term2 =   -(2.0*(zdc*zca+xdc*xca))/rtru
                         -( dot*(2.0*rt2*(xdc*xdc+zdc*zdc)+(8.0*(zcaxt)*(xdczu))+2.0*ru2*(xca*xca+zca*zca) ) )/(2.0*rt2ru2*rtru)
                         -( (zt*xdc-xt*zdc+xu*zca-zu*xca)*(rt2*2.0*(xdczu)+ru2*2.0*(zcaxt)) )/(rt2ru2*rtru)
                         +3.0*dot*(rt2*2.0*(xdczu)+ru2*2.0*(zcaxt))
                         *(rt2*2.0*(xdczu)+ru2*2.0*(zcaxt))/(4.0*rt2ru2*rt2ru2*rtru);
               term3 = ( (2.0*(xdczu)*dot*dot)/(rt2*ru2*ru2)
                        -(2.0*(zt*xdc-xt*zdc+xu*zca-zu*xca)*dot)/rt2ru2
                        +(2.0*(zcaxt)*dot*dot)/(rt2*rt2*ru2) )*(term1);
               dyibyib = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//               
               term1 = (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru -(dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt)))/(rtru322);
               term2 =   -(2.0*(ydc*yca+xdc*xca))/rtru
                         -(dot*( 2.0*rt2*(ydc*ydc+xdc*xdc) + (2.0*(2.0*yt*xca-2.0*xt*yca)*(2.0*ydc*xu-2.0*xdc*yu))
                         +(2.0*yca*yca+2.0*xca*xca)*ru2) )/(2.0*rt2ru2*rtru)
                         -( (xt*ydc - yt*xdc + yu*xca - xu*yca)*( (rt2*(2.0*ydc*xu-2.0*xdc*yu))+(ru2*(2.0*yt*xca-2.0*xt*yca)) ))/
                           (rt2ru2*rtru)
                         +(3.0*dot*(rt2*(2.0*ydc*xu-2.0*xdc*yu)+ru2*(2.0*yt*xca-2.0*xt*yca))
                         *( rt2*(2.0*ydc*xu-2.0*xdc*yu)+ru2*(2.0*yt*xca-2.0*xt*yca)) )/(4.0*rt2ru2*rt2ru2*rtru);
               term3 = ( (2.0*(ydcxu)*dot*dot)/(rt2*ru2*ru2)
                        -(2.0*(xt*ydc-yt*xdc+yu*xca-xu*yca)*dot)/rt2ru2
                        +(2.0*(xcayt)*dot*dot)/(rt2*rt2*ru2) )*(term1);
               dzibzib = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               yuzab = yu*zba - zu*yba;
               zuydb = zu*ydb - yu*zdb;
               ytzba = yt*zba - yba*zt;

               term1 = (zt*ydb-yt*zdb+zba*yu-zu*yba)/rtru - (dot* (rt2*(2.0*zuydb)+ru2*(2.0*ytzba)))/(2.0*rt2ru2*rtru);
               term2 = -(2.0*(yba*ydb+zba*zdb))/rtru
                       -(dot*(rt2*2.0*(ydb*ydb+zdb*zdb) + 2.0*(2.0*ytzba)*(2.0*zuydb) + ru2*2.0*(yba*yba+zba*zba)))/(2.0*rt2ru2*rtru)
                       -( (zt*ydb-yba*zu-yt*zdb+zba*yu)*(rt2*(2.0*zuydb) + ru2*(2.0*ytzba)))/(rt2ru2*rtru)
                       +( (3.0*dot*(rt2*(2.0*zuydb) + ru2*(2.0*ytzba)))*
                          ( rt2*(2.0*zuydb) + ru2*(2.0*ytzba)))/(4.0*rt2ru2*rt2ru2*rtru);

               term3 = ( (dot*dot*(2.0*zuydb))/(rt2*ru2*ru2)
                        -(2.0*dot*(zt*ydb-yba*zu-yt*zdb+zba*yu))/rt2ru2 + ((2.0*ytzba)*dot*dot)/(rt2*rt2*ru2))*( term1);  

               dxicxic = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = (xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*(zdb*xu-zu*xdb)+ru2*(zt*xba-xt*zba)))/(rtru*rt2ru2);
               term3 = ( dot*dot*2.0*(zdb*xu-zu*xdb)/(rt2*ru2*ru2) - 2.0*dot*(xt*zdb-zt*xdb+zu*xba-xu*zba)/rt2ru2 +
                         dot*dot*2.0*(xba*zt-xt*zba)/(rt2*rt2*ru2))*(term1);

               term2 = 3.0*dot*(rt2*2.0*(xu*zdb-zu*xdb)+ru2*2.0*(zt*xba-xt*zba))*(rt2*2.0*(xu*zdb-zu*xdb)+ru2*2.0*(zt*xba-xt*zba))/(4.0*rt2ru2*rt2ru2*rtru)
                      - ( dot*(rt2*2.0*(zdb*zdb+xdb*xdb)+8.0*(xba*zt-xt*zba)*(zdb*xu-zu*xdb)+ru2*2.0*(xba*xba+zba*zba) ))/(2.0*rt2ru2*rtru)
                      - (xt*zdb-zt*xdb+zu*xba-xu*zba)*(rt2*2.0*(xu*zdb-zu*xdb)+ru2*2.0*(zt*xba-xt*zba))/(rtru*rt2ru2)
                      - 2.0*(xba*xdb+zba*zdb)/rtru ;          
               dyicyic = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = (yt*xdb-xt*ydb+xu*yba-yu*xba)/rtru -(dot*(rt2*2.0*(xdb*yu-ydb*xu)+ru2*2.0*(xt*yba-yt*xba)))/(rtru322);
               term2 =   -(2.0*(ydb*yba+xdb*xba))/rtru
                         -(dot*( 2.0*rt2*(ydb*ydb+xdb*xdb) + (2.0*(2.0*xt*yba-2.0*yt*xba)*(2.0*xdb*yu-2.0*ydb*xu))
                         +(2.0*yba*yba+2.0*xba*xba)*ru2) )/(2.0*rt2ru2*rtru)
                         -( (yt*xdb - xt*ydb + xu*yba - yu*xba)*( (rt2*(2.0*xdb*yu-2.0*ydb*xu))+(ru2*(2.0*xt*yba-2.0*yt*xba)) ))/
                           (rt2ru2*rtru)
                         +(3.0*dot*(rt2*2.0*(xdb*yu-ydb*xu)+ru2*2.0*(xt*yba-yt*xba))
                         *( rt2*2.0*(xdb*yu-ydb*xu)+ru2*2.0*(xt*yba-yt*xba)) )/(4.0*rt2ru2*rt2ru2*rtru);
               term3 = ( (2.0*(xdb*yu-ydb*xu)*dot*dot)/(rt2*ru2*ru2) - (2.0*(yt*xdb - xt*ydb + xu*yba - yu*xba)*dot)/rt2ru2
                        +(2.0*(xt*yba-yt*xba)*dot*dot)/(rt2*rt2*ru2) )*(term1);
               dziczic = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);

//
               term1 = ( dphidxt/rtru - (dphidxu*dot*rt2)/(rtru*rt2ru2));
               term2 = ( (3.0*4.0*dphidxu*dphidxu*dot*rt2*rt2)/(4.0*rt2ru2*rt2ru2*rtru) - (2.0*dphidxt*dphidxu*rt2)/(rt2ru2*rtru)
                        -( (zcb*zcb+ycb*ycb)*dot*rt2)/(rt2ru2*rtru));
               term3 = ( (2.0*dphidxu*dot*dot)/(rt2*ru2*ru2) - (2.0*dphidxt*dot)/rt2ru2)*(term1);
               dxidxid = (d2edphi2*term1*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidyt/rtru - (dphidyu*dot*rt2)/(rtru*rt2ru2));
               term2 = ( (3.0*4.0*dphidyu*dphidyu*dot*rt2*rt2)/(4.0*rt2ru2*rt2ru2*rtru) - (2.0*dphidyt*dphidyu*rt2)/(rt2ru2*rtru)
                        -( (zcb*zcb+xcb*xcb)*dot*rt2)/(rt2ru2*rtru));
               term3 = ( (2.0*dphidyu*dot*dot)/(rt2*ru2*ru2) - (2.0*dphidyt*dot)/rt2ru2)*(term1);
               dyidyid = (d2edphi2*term1*term1)/denom2 + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidzt/rtru - (dphidzu*dot*rt2)/(rtru*rt2ru2));
               term2 = ( (3.0*4.0*dphidzu*dphidzu*dot*rt2*rt2)/(4.0*rt2ru2*rt2ru2*rtru) - (2.0*dphidzt*dphidzu*rt2)/(rt2ru2*rtru)
                        -( (ycb*ycb+xcb*xcb)*dot*rt2)/(rt2ru2*rtru));
               term3 = ( (2.0*dphidzu*dot*dot)/(rt2*ru2*ru2) - (2.0*dphidzt*dot)/rt2ru2)*(term1);
               dzidzid = (d2edphi2*term1*term1)/denom2 + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
// ========================================================
//  off diagonal terms
// ========================================================                                                            
               term1 = ( dphidxu/rtru - (dphidxt*dot*ru2)/(rtru32) )*
                       ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - (dot*(rt2*2.0*(zdcyu)+ru2*2.0*(ycazt)))/(rtru322) );
               term2 = (3.0*2.0*(zcb*yt-zt*ycb)*dot*ru2*(rt2*2.0*(zdcyu)+ru2*2.0*(ycazt)))/(4.0*rt2ru2*rt2ru2*rtru) 
                      - (dot*4.0*(zcb*yt-zt*ycb)*(zdcyu))/(2.0*rt2ru2*rtru)
                      - (2.0*(zcb*yt-zt*ycb)*ru2*(yt*zdc-zt*ydc+zu*yca-yu*zca))/(2.0*rt2ru2*rtru)
                      - (dot*ru2*2.0*(-yca*ycb+zca*(-zcb)))/(rtru322)
                      - (zcb*yu-zu*ycb)*(rt2*2.0*(zdcyu)+ru2*2.0*(ycazt))/(2.0*rt2ru2*rtru)
                      + (ycb*ydc + zcb*zdc)/rtru;
               term3 = ( 2.0*(zdcyu)*dot*dot/(rt2*ru2*ru2)-2.0*dot*(yt*zdc-zt*ydc+zu*yca-yu*zca)/rt2ru2
                        + 2.0*dot*dot*(ycazt)/(rt2*rt2*ru2))*
                       (dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));
               dxiaxib = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//                                                    
               term1 = ( dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2) )*
                       ( (zt*xdc-xt*zdc+xu*zca-zu*xca)/rtru - (dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zcaxt)))/(rtru322) );
               term3 = ( 2.0*dot*dot*(xdc*zu-xu*zdc)/(rt2*ru2*ru2)-2.0*dot*(zt*xdc-xt*zdc+xu*zca-zu*xca)/rt2ru2
                       + 2.0*dot*dot*(zcaxt)/(rt2*rt2*ru2))*
                       ( dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));
               term2 = (3.0*2.0*dphidxt*dot*ru2*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zcaxt)))/(4.0*rt2ru2*rt2ru2*rtru) 
                      - ( dot*4.0*dphidxt*(xdc*zu-xu*zdc))/rtru322
                      - ( dphidxt*2.0*ru2*(zt*xdc-xt*zdc+xu*zca-zu*xca))/rtru322
                      - ( dot*ru2*2.0*(xca*ycb+zt))/rtru322
                      - ( dphidxu*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zcaxt)))/rtru322
                      + ( xcb*ydc - 2.0*ycb*xdc)/rtru;                      
               dxiayib = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//
               term1 = ( dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2) )*
                       ( (xt*ydc-yt*xdc+xca*yu-xu*yca)/rtru - (dot*(rt2*2.0*(xu*ydc-yu*xdc)+ru2*2.0*(xcayt)))/(rtru322) );
               term3 = ( 2.0*dot*dot*(xu*ydc-yu*xdc)/(rt2*ru2*ru2)-2.0*dot*(xt*ydc-yt*xdc+xca*yu-xu*yca)/rt2ru2
                       + 2.0*dot*dot*(xcayt)/(rt2*rt2*ru2))*
                       (dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));                            
               term2 = (3.0*2.0*dphidxt*dot*ru2*(rt2*2.0*(xu*ydc-yu*xdc)+ru2*2.0*(xcayt)))/(4.0*rt2ru2*rt2ru2*rtru) 
                      - ( dot*4.0*dphidxt*(xu*ydc-yu*xdc))/rtru322
                      - ( dphidxt*2.0*ru2*(xt*ydc-yt*xdc+xca*yu-xu*yca))/rtru322
                      - ( dot*ru2*2.0*(xca*zcb-yt))/rtru322
                      - ( dphidxu*(rt2*2.0*(xu*ydc-yu*xdc)+ru2*2.0*(xcayt)))/rtru322
                      + ( xcb*zdc - 2.0*zcb*xdc)/rtru;                          
               dxiazib = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//                              
               term1 = ( dphidxu/rtru - (ru2*dot*dphidxt)/(rtru*rt2ru2))*
                       ( (zt*ydb-yt*zdb+zba*yu-zu*yba)/rtru - (dot*(rt2*2.0*(zu*ydb-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/(rtru322));
               term3 = ( 2.0*dot*dot*(zu*ydb-yu*zdb)/(rt2*ru2*ru2)-2.0*dot*(zt*ydb-yt*zdb+zba*yu-zu*yba)/rt2ru2
                       + 2.0*dot*dot*(zba*yt-zt*yba)/(rt2*rt2*ru2))*
                       (dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));                       
               term2 = (3.0*2.0*dphidxt*dot*ru2*(rt2*2.0*(zu*ydb-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/(4.0*rt2ru2*rt2ru2*rtru) 
                      - ( 4.0*dphidxt*dot*(zu*ydb-yu*zdb))/(2.0*rt2ru2*rtru)
                      - ( 2.0*dphidxt*ru2*(zt*ydb-yt*zdb+zba*yu-zu*yba))/(2.0*rt2ru2*rtru)
                      - ( dot*ru2*2.0*(zba*zcb+yba*ycb))/(2.0*rt2ru2*rtru)
                      - ( dphidxu*(rt2*2.0*(zu*ydb-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/(rtru322)
                      + (-ycb*ydb + zcb*(-zdb))/rtru;                                    
               dxiaxic = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//               
               term1 = ( dphidxu/rtru - (ru2*dot*dphidxt)/(rtru*rt2ru2))*
                       ( (xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(xu*zdb-zu*xdb)+ru2*2.0*(xba*zt-xt*zba)))/(rtru322));                        
               term3 = ( 2.0*dot*dot*(xu*zdb-zu*xdb)/(rt2*ru2*ru2)-2.0*dot*(xt*zdb-zt*xdb+zu*xba-xu*zba)/rt2ru2
                        + 2.0*dot*dot*(xba*zt-xt*zba)/(rt2*rt2*ru2))*
                       (dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));                       
               term2 = (3.0*2.0*dphidxt*dot*ru2*(rt2*2.0*(xu*zdb-zu*xdb)+ru2*2.0*(xba*zt-xt*zba)))/(4.0*rt2ru2*rt2ru2*rtru) 
                      - ( 4.0*dphidxt*dot*(xu*zdb-zu*xdb))/(2.0*rt2ru2*rtru)
                      - ( 2.0*dphidxt*ru2*(xt*zdb-zt*xdb+zu*xba-xu*zba))/(2.0*rt2ru2*rtru)
                      - ( dot*ru2*2.0*(-xba*ycb-zt))/(2.0*rt2ru2*rtru)
                      - ( dphidxu*(rt2*2.0*(xu*zdb-zu*xdb)+ru2*2.0*(xba*zt-xt*zba)) )/(rtru322)
                      + (xdc*ycb-xcb*ydc+xdb*ycb)/rtru;                                    
               dxiayic = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//                            
               term1 = ( dphidxu/rtru - (ru2*dot*dphidxt)/(rtru*rt2ru2))*
                       ( (yt*xdb-xt*ydb+xu*yba-yu*xba)/rtru - (dot*(rt2*2.0*(yu*xdb-xu*ydb)+ru2*2.0*(xt*yba-yt*xba)))/(rtru322));                        
               term3 = ( 2.0*dot*dot*(yu*xdb-xu*ydb)/(rt2*ru2*ru2)-2.0*dot*(yt*xdb-xt*ydb+xu*yba-yu*xba)/rt2ru2
                        + 2.0*dot*dot*(xt*yba-yt*xba)/(rt2*rt2*ru2))*
                       (dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));                      

               term2 = (3.0*2.0*dphidxt*dot*ru2*(rt2*2.0*(yu*xdb-xu*ydb)+ru2*2.0*(xt*yba-yt*xba)))/(4.0*rt2ru2*rt2ru2*rtru) 
                      - ( 4.0*dphidxt*dot*(yu*xdb-xu*ydb))/(2.0*rt2ru2*rtru)
                      - ( 2.0*dphidxt*ru2*(yt*xdb-xt*ydb+xu*yba-yu*xba))/(2.0*rt2ru2*rtru)
                      - ( dot*ru2*2.0*(yt-xba*zcb))/(2.0*rt2ru2*rtru)
                      - ( dphidxu*(rt2*2.0*(yu*xdb-xu*ydb)+ru2*2.0*(xt*yba-yt*xba)) )/(rtru322)
                      + (xdb*zcb+xdc*zcb-xcb*zdc)/rtru;                                    
               dxiazic = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//
               term1 = ( dphidxt/rtru - (rt2*dot*dphidxu)/(rtru*rt2ru2))*( dphidxu/rtru - (ru2*dot*dphidxt)/(rtru*rt2ru2));                        
               term3 = ( 2.0*dot*dot*dphidxu/(rt2*ru2*ru2) - 2.0*dot*dphidxt/rt2ru2)*(dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));
               term2 = (3.0*4.0*rt2*ru2*dot*dphidxt*dphidxu)/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidxu*dphidxu)/(rtru322)
                       - ( dot*4.0*dphidxt*dphidxu)/(rtru322)
                       - ( ru2*2.0*dphidxt*dphidxt)/(rtru322)
                       + (zcb*zcb+ycb*ycb)/rtru;                      
               dxiaxid = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//                                      
               term1 = ( dphidyt/rtru - (rt2*dot*dphidyu)/(rtru*rt2ru2))*( dphidxu/rtru - (ru2*dot*dphidxt)/(rtru*rt2ru2));                        
               term3 = ( 2.0*dot*dot*dphidyu/(rt2*ru2*ru2) - 2.0*dot*dphidyt/rt2ru2)*(dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));
               term2 = (3.0*4.0*rt2*ru2*dot*dphidxt*dphidyu)/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidxu*dphidyu)/(rtru322)
                       - ( dot*4.0*dphidxt*dphidyu)/(rtru322)
                       - ( ru2*2.0*dphidxt*dphidyt)/(rtru322)
                       - (xcb*ycb)/rtru;
               dxiayid = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//               
               term1 = ( dphidzt/rtru - (rt2*dot*dphidzu)/(rtru*rt2ru2))*( dphidxu/rtru - (ru2*dot*dphidxt)/(rtru*rt2ru2));                        
               term3 = ( 2.0*dot*dot*dphidzu/(rt2*ru2*ru2) - 2.0*dot*dphidzt/rt2ru2)*(dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));
               term2 = (3.0*4.0*rt2*ru2*dot*dphidxt*dphidzu)/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidxu*dphidzu)/(rtru322)
                       - ( dot*4.0*dphidxt*dphidzu)/(rtru322)
                       - ( ru2*2.0*dphidxt*dphidzt)/(rtru322)
                       - (xcb*zcb)/rtru;
               dxiazid = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//                                                  
               term1 = ( dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2) )*( dphidyu/rtru-(dot*ru2*dphidyt)/(rtru*rt2ru2) );                       
               term3 = ( 2.0*dot*dot*dphidyt/(rt2*rt2*ru2)-2.0*dot*dphidyu/rt2ru2)*(dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));                       
               term2 = (3.0*4.0*dphidxt*dot*ru2*ru2*dphidyt)/(4.0*rt2ru2*rt2ru2*rtru) 
                      - ( ru2*2.0*dphidyt*dphidxu)/rtru322
                      - ( ru2*2.0*dphidyu*dphidxt)/rtru322
                      + ( xcb*ycb*dot*ru2)/rtru32;                                  
               dxiayia = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//                                      
               term1 = ( dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2) )*( dphidzu/rtru-(dot*ru2*dphidzt)/(rtru*rt2ru2) );                       
               term3 = ( 2.0*dot*dot*dphidzt/(rt2*rt2*ru2)-2.0*dot*dphidzu/rt2ru2)*(dphidxu/rtru - (dphidxt*dot*ru2)/(rtru*rt2ru2));                       
               term2 = (3.0*4.0*dphidxt*dot*ru2*ru2*dphidzt)/(4.0*rt2ru2*rt2ru2*rtru) 
                      - ( ru2*2.0*dphidzt*dphidxu)/rtru322
                      - ( ru2*2.0*dphidzu*dphidxt)/rtru322
                      + ( xcb*zcb*dot*ru2)/rtru32;                                  
               dxiazia = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//
               term1 = ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322)*
                       ( (zt*xdc-xt*zdc+xu*zca-zu*xca)/rtru - dot*(rt2*2.0*(xdczu)+ru2*2.0*(zcaxt))/rtru322);
               term3 = ( dot*dot*2.0*(xdczu)/(rt2*ru2*ru2) - dot*2.0*(zt*xdc-xt*zdc+xu*zca-zu*xca)/rt2ru2
                         + dot*dot*2.0*(zcaxt)/(rt2*rt2*ru2))*
                       ((yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))*(rt2*2.0*(xdczu)+ru2*2.0*(zcaxt))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(-rt2*2.0*xdc*ydc + 4.0*(zcaxt)*(zdcyu) + 4.0*(yca*zt-yt*zca)*(xdczu) - 2.0*ru2*xca*yca) )/rtru322
                        - ( (zt*xdc-xt*zdc+xu*zca-zu*xca)*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)) )/rtru322
                        - ( (yt*zdc-zt*ydc+zu*yca-yu*zca)*(rt2*2.0*(xdczu)+ru2*2.0*(zcaxt)) )/rtru322
                        + (xdc*yca+xca*ydc)/rtru;
               dxibyib = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//
               term1 = ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322)*
                       ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))/rtru322);
               term3 = ( dot*dot*2.0*(ydcxu)/(rt2*ru2*ru2) - dot*2.0*(xt*ydc-yt*xdc+yu*xca-xu*yca)/rt2ru2
                         + dot*dot*2.0*(xcayt)/(rt2*rt2*ru2))*
                       ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);

               term2 =  3.0*dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(-rt2*2.0*xdc*zdc + 4.0*(yca*zt-yt*zca)*(ydcxu) + 4.0*(xcayt)*(zdcyu) - 2.0*ru2*xca*zca))/rtru322
                        - ( (xt*ydc-yt*xdc+yu*xca-xu*yca)*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)))/rtru322
                        - ( (yt*zdc-zt*ydc+zu*yca-yu*zca)*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt)))/rtru322
                        + (xca*zdc+xdc*zca)/rtru;
               dxibzib = (d2edphi2*term1)/(denom2) + (dedphi*term2)/denom - (dedphi*term3)/(2.0*denom3);
//
               term1 = ( (zt*ydb-yt*zdb+yu*zba-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/rtru322)*
                       ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term3 = ( dot*dot*2.0*(ydb*zu-yu*zdb)/(rt2*ru2*ru2) - dot*2.0*(zt*ydb-yt*zdb+yu*zba-zu*yba)/rt2ru2
                        + dot*dot*2.0*(zba*yt-zt*yba)/(rt2*rt2*ru2))*
                       ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322 );
               term2 =  3.0*dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(rt2*2.0*(-ydc*ydb-zdb*zdc)+4.0*(yca*zt-yt*zca)*(ydb*zu-yu*zdb)
                            +4.0*(zba*yt-zt*yba)*(zdcyu)+ru2*2.0*(-yba*yca-zca*zba)))/rtru322
                        - ( (yt*zdc-zt*ydc+zu*yca-yu*zca)*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/rtru322
                        - ( (zt*ydb-yt*zdb+yu*zba-zu*yba)*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)))/rtru322
                        + (yca*ydb+yba*ydc+zca*zdb+zba*zdc)/rtru;
               dxibxic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( (xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/rtru322)*
                       ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term3 = ( dot*dot*2.0*(zdb*xu-zu*xdb)/(rt2*ru2*ru2) - dot*2.0*(xt*zdb-zt*xdb+zu*xba-xu*zba)/rt2ru2
                         + dot*dot*2.0*(xba*zt-xt*zba)/(rt2*rt2*ru2))*
                       ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(rt2*2.0*(xdb*ydc+zu)+4.0*(xba*zt-xt*zba)*(zdcyu)+4.0*(yca*zt-yt*zca)*(zdb*xu-zu*xdb)+ru2*2.0*(xba*yca+zt)))/rtru322
                        - ( (xt*zdb-zt*xdb+zu*xba-xu*zba)*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)))/rtru322
                        - ( (yt*zdc-zt*ydc+zu*yca-yu*zca)*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba)))/rtru322
                        + (xba*ycb-xcb*yba-xdb*yca-xdc*ycb-xba*ydc+xcb*ydc)/rtru;
               dxibyic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 = ( (yt*xdb-xt*ydb+xu*yba-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322)*
                       ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322 );
               term3 = ( dot*dot*2.0*(xdb*yu-xu*ydb)/(rt2*ru2*ru2) - dot*2.0*(yt*xdb-xt*ydb+xu*yba-yu*xba)/rt2ru2
                       + dot*dot*2.0*(yba*xt-yt*xba)/(rt2*rt2*ru2))*
                       ( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(rt2*2.0*(xdb*zdc-yu)+4.0*(yba*xt-yt*xba)*(zdcyu)+
                           4.0*(yca*zt-yt*zca)*(xdb*yu-xu*ydb)+ru2*2.0*(xba*zca-yt)))/rtru322
                        - ( (yt*xdb-xt*ydb+xu*yba-yu*xba)*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)))/rtru322
                        - ( (yt*zdc-zt*ydc+zu*yca-yu*zca)*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/rtru322
                        + (xba*zcb-xcb*zba-xdb*zca-xdc*zcb-xba*zdc+xcb*zdc)/rtru;
               dxibzic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 = ( dphidxt/rtru - (dphidxu*dot*rt2)/(rtru*rt2ru2))*
                       ( (yt*zdc-zt*ydc+yca*zu-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term3 = ( dot*dot*2.0*dphidxu/(rt2*ru2*ru2) - 2.0*dot*dphidxt/rt2ru2)*
                       ( (yt*zdc-zt*ydc+yca*zu-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);

               term2 = 3.0*2.0*dphidxu*rt2*dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 2.0*dphidxu*rt2*(yt*zdc-zt*ydc+yca*zu-yu*zca))/rtru322
                       - ( dot*(rt2*2.0*(ycb*ydc+zcb*zdc)+4.0*dphidxu*(yca*zt-yt*zca)))/rtru322
                       - ( dphidxt*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)))/rtru322
                       + (-ycb*yca-zca*zcb)/rtru;
               dxibxid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidyt/rtru - (dphidyu*dot*rt2)/(rtru*rt2ru2))*
                       ( (yt*zdc-zt*ydc+yca*zu-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term3 = ( dot*dot*2.0*dphidyu/(rt2*ru2*ru2) - 2.0*dot*dphidyt/rt2ru2)*
                       ( (yt*zdc-zt*ydc+yca*zu-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term2 = 3.0*2.0*dphidyu*rt2*dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 2.0*dphidyu*rt2*(yt*zdc-zt*ydc+yca*zu-yu*zca))/rtru322
                       - ( dot*(rt2*2.0*(-xcb*ydc-zu)+4.0*dphidyu*(yca*zt-yt*zca)))/rtru322
                       - ( dphidyt*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)))/rtru322
                       + (xcb*yba+xcb*yca-xba*ycb)/rtru;
               dxibyid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 = ( dphidzt/rtru - (dphidzu*dot*rt2)/(rtru*rt2ru2))*
                       ( (yt*zdc-zt*ydc+yca*zu-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term3 = ( dot*dot*2.0*dphidzu/(rt2*ru2*ru2) - 2.0*dot*dphidzt/rt2ru2)*
                       ( (yt*zdc-zt*ydc+yca*zu-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term2 = 3.0*2.0*dphidzu*rt2*dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 2.0*dphidzu*rt2*(yt*zdc-zt*ydc+yca*zu-yu*zca))/rtru322
                       - ( dot*(rt2*2.0*(yu-xcb*zdc)+4.0*dphidzu*(yca*zt-yt*zca)))/rtru322
                       - ( dphidzt*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)))/rtru322
                       + (xcb*zba+xcb*zca-xba*zcb)/rtru;
               dxibzid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
//
               term1 = ( dphidzu/rtru - (dphidzt*dot*ru2)/(rtru*rt2ru2))*( dphidyu/rtru - (dphidyt*dot*ru2)/(rtru*rt2ru2));
               term3 = ( (2.0*dphidzt*dot*dot)/(rt2*rt2*ru2) - (2.0*dphidzu*dot)/rt2ru2)*
                         (dphidyu/rtru - (dphidyt*dot*ru2)/(rt2ru2*rtru) );
               term2 = ( (3.0*4.0*dphidzt*dphidyt*dot*ru2*ru2)/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 2.0*dphidyt*dphidzu*ru2)/rtru322
                       - ( 2.0*dphidzt*dphidyu*ru2)/rtru322
                       + ( (ycb*zcb)*dot*ru2)/(rt2ru2*rtru));
               dyiazia = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2))*( (yt*zdc-zt*ydc+zu*yca-yu*zca)/rtru -
                         dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term3 = ( dot*dot*2.0*(zdcyu)/(rt2*ru2*ru2) - dot*2.0*(yt*zdc-zt*ydc+zu*yca-yu*zca)/rt2ru2 +
                         dot*dot*2.0*(yca*zt-yt*zca)/(rt2*rt2*ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 =   3.0*2.0*dphidyt*dot*ru2*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( 4.0*dot*dphidyt*(zdcyu))/rtru322
                        - ( 2.0*dphidyt*ru2*(yt*zdc-zt*ydc+zu*yca-yu*zca))/rtru322
                        - ( dot*ru2*2.0*(xcb*yca-zt))/rtru322
                        - ( dphidyu*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)) )/rtru322
                        + (ycb*xdc-xcb*ydc-xcb*ydc)/rtru;
               dyiaxib = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                                     
//
               term1 = ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2))*( (zt*xdc-xt*zdc+xu*zca-zu*xca)/rtru -
                         dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zca*xt-zt*xca))/rtru322);
               term3 = ( dot*dot*2.0*(xdc*zu-xu*zdc)/(rt2*ru2*ru2) - dot*2.0*(zt*xdc-xt*zdc+xu*zca-zu*xca)/rt2ru2 +
                         dot*dot*2.0*(zca*xt-zt*xca)/(rt2*rt2*ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 =   3.0*2.0*dphidyt*dot*ru2*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zca*xt-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( 4.0*dot*dphidyt*(xdc*zu-xu*zdc))/rtru322
                        - ( 2.0*dphidyt*ru2*(zt*xdc-xt*zdc+xu*zca-zu*xca))/rtru322
                        + ( dot*ru2*2.0*(xca*xcb+zcb*zca))/rtru322
                        - ( dphidyu*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zca*xt-zt*xca)) )/rtru322
                        + (xcb*xdc+zcb*zdc)/rtru;
               dyiayib = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 = ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2))*( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru -
                         dot*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca))/rtru322);
               term3 = ( dot*dot*2.0*(ydc*xu-yu*xdc)/(rt2*ru2*ru2) - dot*2.0*(xt*ydc-yt*xdc+yu*xca-xu*yca)/rt2ru2 +
                         dot*dot*2.0*(xca*yt-xt*yca)/(rt2*rt2*ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 =   3.0*2.0*dphidyt*dot*ru2*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( 4.0*dot*dphidyt*(ydc*xu-yu*xdc))/rtru322
                        - ( 2.0*dphidyt*ru2*(xt*ydc-yt*xdc+yu*xca-xu*yca))/rtru322
                        - ( dot*ru2*2.0*(yca*zcb+xt))/rtru322
                        - ( dphidyu*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca)) )/rtru322
                        + (ycb*zdc-2.0*ydc*zcb)/rtru;
               dyiazib = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2))*
                       ((zt*ydb-yt*zdb+yu*zba-zu*yba)/rtru-dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/rtru322);
               term3 = ( dot*dot*2.0*(ydb*zu-yu*zdb)/(rt2*ru2*ru2) - dot*2.0*(zt*ydb-yt*zdb+yu*zba-zu*yba)/rt2ru2
                        + dot*dot*2.0*(zba*yt-zt*yba)/(rt2*rt2*ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 =  3.0*2.0*dphidyt*dot*ru2*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidyt*dot*(ydb*zu-yu*zdb))/rtru322
                       - ( 2.0*dphidyt*ru2*(zt*ydb-yt*zdb+yu*zba-zu*yba))/rtru322
                       - ( dot*ru2*2.0*(zt-xcb*yba))/rtru322
                       - ( dphidyu*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/rtru322
                       + (xcb*ydc+xcb*ydb-xdc*ycb)/rtru;
               dyiaxic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                                       
//
               term1 = ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2))*
                       ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru-dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/rtru322);                            
               term3 = ( dot*dot*2.0*(zdb*xu-zu*xdb)/(rt2*ru2*ru2) - dot*2.0*(xt*zdb-zt*xdb+zu*xba-xu*zba)/rt2ru2
                        + dot*dot*2.0*(xba*zt-xt*zba)/(rt2*rt2*ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 =  3.0*2.0*dphidyt*dot*ru2*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidyt*dot*(zdb*xu-zu*xdb))/rtru322
                       - ( 2.0*dphidyt*ru2*(xt*zdb-zt*xdb+zu*xba-xu*zba))/rtru322
                       - ( dot*ru2*2.0*(xba*xcb+zba*zcb))/rtru322
                       - ( dphidyu*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba)))/rtru322
                       - (xcb*xdb+zcb*zdb)/rtru;
               dyiayic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 = ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2))*
                       ((yt*xdb-xt*ydb+xu*yba-yu*xba)/rtru-dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);                            
               term3 = ( dot*dot*2.0*(xdb*yu-xu*ydb)/(rt2*ru2*ru2) - dot*2.0*(yt*xdb-xt*ydb+xu*yba-yu*xba)/rt2ru2
                        + dot*dot*2.0*(yba*xt-yt*xba)/(rt2*rt2*ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 =  3.0*2.0*dphidyt*dot*ru2*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidyt*dot*(xdb*yu-xu*ydb))/rtru322
                       - ( 2.0*dphidyt*ru2*(yt*xdb-xt*ydb+xu*yba-yu*xba))/rtru322
                       + ( dot*ru2*2.0*(yba*zcb+xt))/rtru322
                       - ( dphidyu*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/rtru322
                       + (ydc*zcb+ydb*zcb-ycb*zdc)/rtru;
               dyiazic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidxt/rtru - dot*dphidxu*rt2/(rtru*rt2ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term3 = ( 2.0*dot*dot*dphidxu/(rt2*ru2*ru2) - dot*2.0*dphidxt/rt2ru2)*
                       ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 = 3.0*4.0*dphidyt*rt2*dphidxu*dot*ru2/(4.0*rt2ru2*rt2ru2*rtru)
                       - (rt2*2.0*dphidxu*dphidyu)/(rtru322) 
                       - (4.0*dphidyt*dphidxu*dot)/(rtru322)  - (2.0*dphidxt*dphidyt*ru2)/(rtru322)
                       - xcb*ycb/rtru;
               dyiaxid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                          
//
               term1 = ( dphidyt/rtru - dot*dphidyu*rt2/(rtru*rt2ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term3 = ( 2.0*dot*dot*dphidyu/(rt2*ru2*ru2) - dot*2.0*dphidyt/rt2ru2)*
                       ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 = 3.0*4.0*dphidyt*rt2*dphidyu*dot*ru2/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidyu*dphidyu)/(rtru322)
                       - (4.0*dphidyt*dphidyu*dot)/(rtru322) - (2.0*dphidyt*dphidyt*ru2)/(rtru322)
                       + (xcb*xcb+zcb*zcb)/rtru;
               dyiayid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3); 
//
               term1 = ( dphidzt/rtru - dot*dphidzu*rt2/(rtru*rt2ru2))*(dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));                           
               term3 = ( 2.0*dot*dot*dphidzu/(rt2*ru2*ru2) - dot*2.0*dphidzt/rt2ru2)*
                       ( dphidyu/rtru - dot*ru2*dphidyt/(rtru*rt2ru2));
               term2 = 3.0*4.0*dphidyt*rt2*dphidzu*dot*ru2/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidzu*dphidyu)/(rtru322)
                       - (4.0*dphidyt*dphidzu*dot)/(rtru322) - (2.0*dphidzt*dphidyt*ru2)/(rtru322)
                       - zcb*ycb/rtru;
               dyiazid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3); 
//
               term1 = ( dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2))*
                       ( (yt*zdc-zt*ydc+yca*zu-yu*zca)/rtru - dot*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/rtru322);
               term3 = ( dot*dot*2.0*(zdcyu)/(rt2*ru2*ru2) - dot*2.0*(yt*zdc-zt*ydc+yca*zu-yu*zca)/rt2ru2 +
                         dot*dot*2.0*(yca*zt-yt*zca)/(rt2*rt2*ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term2 =  3.0*2.0*dphidzt*dot*ru2*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidzt*dot*(zdcyu))/rtru322
                       - ( 2.0*dphidzt*ru2*(yt*zdc-zt*ydc+yca*zu-yu*zca))/rtru322
                       - ( 2.0*dot*ru2*(xcb*zca+yt))/rtru322
                       - ( dphidzu*(rt2*2.0*(zdcyu)+ru2*2.0*(yca*zt-yt*zca)))/rtru322
                       + (xdc*zcb-xcb*zdc*2.0)/rtru;
               dziaxib = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 = ( dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2))*
                       ( (zt*xdc-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zca*xt-zt*xca))/rtru322);
               term3 = ( dot*dot*2.0*(xdc*zu-xu*zdc)/(rt2*ru2*ru2) - dot*2.0*(zt*xdc-xt*zdc+zca*xu-zu*xca)/rt2ru2 +
                         dot*dot*2.0*(zca*xt-zt*xca)/(rt2*rt2*ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));

               term2 =  3.0*2.0*dphidzt*dot*ru2*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zca*xt-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidzt*dot*(xdc*zu-xu*zdc))/rtru322
                       - ( 2.0*dphidzt*ru2*(zt*xdc-xt*zdc+zca*xu-zu*xca))/rtru322
                       - ( 2.0*dot*ru2*(ycb*zca-xt))/rtru322
                       - ( dphidzu*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(zca*xt-zt*xca)))/rtru322
                       + (ydc*zcb-ycb*zdc*2.0)/rtru;
               dziayib = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 = ( dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2))*
                       ( (xt*ydc-yt*xdc+xca*yu-xu*yca)/rtru - dot*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca))/rtru322);
               term3 = ( dot*dot*2.0*(ydc*xu-yu*xdc)/(rt2*ru2*ru2) - dot*2.0*(xt*ydc-yt*xdc+xca*yu-xu*yca)/rt2ru2 +
                         dot*dot*2.0*(xca*yt-xt*yca)/(rt2*rt2*ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term2 =  3.0*2.0*dphidzt*dot*ru2*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidzt*dot*(ydc*xu-yu*xdc))/rtru322
                       - ( 2.0*dphidzt*ru2*(xt*ydc-yt*xdc+xca*yu-xu*yca))/rtru322
                       + ( dot*ru2*2.0*(xcb*xca+yca*ycb))/rtru322
                       - (dphidzu*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca)))/rtru322
                       + (ycb*ydc+xcb*xdc)/rtru;
               dziazib = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2))*
                       ( (zt*ydb-yt*zdb+zba*yu-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/rtru322);
               term3 = ( dot*dot*2.0*(ydb*zu-yu*zdb)/(rt2*ru2*ru2) - dot*2.0*(zt*ydb-yt*zdb+zba*yu-zu*yba)/rt2ru2
                         + dot*dot*2.0*(zba*yt-zt*yba)/(rt2*rt2*ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term2 =  3.0*2.0*dphidzt*dot*ru2*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidzt*dot*(ydb*zu-yu*zdb))/rtru322
                       - ( 2.0*dphidzt*ru2*(zt*ydb-yt*zdb+zba*yu-zu*yba))/rtru322
                       - ( 2.0*dot*ru2*(-xcb*zba-yt))/rtru322
                       - ( dphidzu*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/rtru322
                       + (xcb*zdc+xcb*zdb-xdc*zcb)/rtru;
               dziaxic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//             
               term1 = ( dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2))*
                       ( (xt*zdb-zt*xdb+xba*zu-xu*zba)/rtru - dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/rtru322);
               term3 = ( dot*dot*2.0*(zdb*xu-zu*xdb)/(rt2*ru2*ru2) - dot*2.0*(xt*zdb-zt*xdb+xba*zu-xu*zba)/rt2ru2
                         + dot*dot*2.0*(xba*zt-xt*zba)/(rt2*rt2*ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));

               term2 =  3.0*2.0*dphidzt*dot*ru2*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidzt*dot*(zdb*xu-zu*xdb))/rtru322
                       - ( 2.0*dphidzt*ru2*(xt*zdb-zt*xdb+xba*zu-xu*zba))/rtru322
                       - ( 2.0*dot*ru2*(xt-ycb*zba))/rtru322
                       - ( dphidzu*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba)))/rtru322
                       + (ycb*zdc+ycb*zdb-ydc*zcb)/rtru;
               dziayic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                      
               term1 = ( dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2))*
                       ( (yt*xdb-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);
               term3 = ( dot*dot*2.0*(xdb*yu-xu*ydb)/(rt2*ru2*ru2) - dot*2.0*(yt*xdb-xt*ydb+yba*xu-yu*xba)/rt2ru2
                         + dot*dot*2.0*(yba*xt-yt*xba)/(rt2*rt2*ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term2 =  3.0*2.0*dphidzt*dot*ru2*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 4.0*dphidzt*dot*(xdb*yu-xu*ydb))/rtru322
                       - ( 2.0*dphidzt*ru2*(yt*xdb-xt*ydb+yba*xu-yu*xba))/rtru322
                       - ( 2.0*dot*ru2*(yba*ycb+xba*xcb))/rtru322
                       - ( dphidzu*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/rtru322
                       - (ycb*ydb+xcb*xdb)/rtru;
               dziazic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1  =  (dphidxt/rtru - dot*rt2*dphidxu/(rtru*rt2ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term3  =  ( dot*dot*2.0*dphidxu/(rt2*ru2*ru2) - 2.0*dot*dphidxt/rt2ru2)*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));

               term2 =  3.0*4.0*dphidzt*rt2*dphidxu*dot*ru2/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( rt2*2.0*dphidxu*dphidzu)/rtru322
                        - ( 4.0*dphidzt*dphidxu*dot)/rtru322
                        - ( 2.0*dphidxt*dphidzt*ru2)/rtru322
                        - xcb*zcb/rtru;
               dziaxid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1  =  (dphidyt/rtru - dot*rt2*dphidyu/(rtru*rt2ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term3  =  ( dot*dot*2.0*dphidyu/(rt2*ru2*ru2) - 2.0*dot*dphidyt/rt2ru2)*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term2 =  3.0*4.0*dphidzt*rt2*dphidyu*dot*ru2/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( rt2*2.0*dphidyu*dphidzu)/rtru322
                        - ( 4.0*dphidzt*dphidyu*dot)/rtru322
                        - ( 2.0*dphidyt*dphidzt*ru2)/rtru322
                        - ycb*zcb/rtru;
               dziayid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1  =  (dphidzt/rtru - dot*rt2*dphidzu/(rtru*rt2ru2))*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term3  =  ( dot*dot*2.0*dphidzu/(rt2*ru2*ru2) - 2.0*dot*dphidzt/rt2ru2)*(dphidzu/rtru - dot*ru2*dphidzt/(rtru*rt2ru2));
               term2 =  3.0*4.0*dphidzt*rt2*dphidzu*dot*ru2/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( rt2*2.0*dphidzu*dphidzu)/rtru322
                        - ( 4.0*dphidzt*dphidzu*dot)/rtru322
                        - ( 2.0*dphidzt*dphidzt*ru2)/rtru322
                        + ( xcb*xcb+ycb*ycb)/rtru;
               dziazid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                                        
               term1 = ( (xt*ydc-yt*xdc+xca*yu-xu*yca)/rtru - dot*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca))/rtru322)*
                       ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term3 = ( dot*dot*2.0*(ydc*xu-yu*xdc)/(rt2*ru2*ru2) - dot*2.0*(xt*ydc-yt*xdc+xca*yu-xu*yca)/rt2ru2 +
                         dot*dot*2.0*(xca*yt-xt*yca)/(rt2*rt2*ru2))*
                       ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term2 = 3.0*dot*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca))*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( dot*(2.0*rt2*ydc*(-zdc)+4.0*(xt*zca-zt*xca)*(ydc*xu-yu*xdc)+4.0*(xca*yt-xt*yca)*(xdc*zu-xu*zdc)-2.0*yca*zca*ru2))/rtru322
                       - ( (xdc*zt-xt*zdc+zca*xu-zu*xca)*(rt2*2.0*(ydc*xu-yu*xdc)+ru2*2.0*(xca*yt-xt*yca)))/rtru322
                       - ( (xt*ydc-yt*xdc+xca*yu-xu*yca)*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca)))/rtru322
                       + (ydc*zca+yca*zdc)/rtru;
               dyibzib = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( (ydb*zt-yt*zdb+zba*yu-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-zdb*yu)+ru2*2.0*(yt*zba-zt*yba))/rtru322)*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term3 =  ( dot*dot*2.0*(ydb*zu-zdb*yu)/(rt2*ru2*ru2) - dot*2.0*(ydb*zt-yt*zdb+zba*yu-zu*yba)/rt2ru2
                        + dot*dot*2.0*(yt*zba-zt*yba)/(rt2*rt2*ru2))*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);                         
               term2 =  3.0*dot*(rt2*2.0*(ydb*zu-zdb*yu)+ru2*2.0*(yt*zba-zt*yba))*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(rt2*2.0*(xdc*ydb-zu)+4.0*(zca*xt-zt*xca)*(ydb*zu-yu*zdb)
                           +4.0*(zba*yt-zt*yba)*(xdc*zu-xu*zdc)+ru2*2.0*(xca*yba-zt)))/rtru322
                        - ( (xdc*zt-xt*zdc+zca*xu-zu*xca)*(rt2*2.0*(ydb*zu-zdb*yu)+ru2*2.0*(yt*zba-zt*yba)))/rtru322
                        - ( (ydb*zt-yt*zdb+zba*yu-zu*yba)*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca)) )/rtru322
                        + (xcb*yba-xdc*yba-xba*ycb+xdc*ycb-xca*ydb-xcb*ydc)/rtru;
               dyibxic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( (zdb*xt-zt*xdb+xba*zu-xu*zba)/rtru - dot*(rt2*2.0*(zdb*xu-xdb*zu)+ru2*2.0*(zt*xba-xt*zba))/rtru322)*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term3 =  ( dot*dot*2.0*(zdb*xu-xdb*zu)/(rt2*ru2*ru2) - dot*2.0*(zdb*xt-zt*xdb+xba*zu-xu*zba)/rt2ru2
                         + dot*dot*2.0*(zt*xba-xt*zba)/(rt2*rt2*ru2))*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(zdb*xu-xdb*zu)+ru2*2.0*(zt*xba-xt*zba))*
                                (rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(rt2*2.0*(-xdb*xdc-zdc*zdb)+4.0*(zt*xba-xt*zba)*(xdc*zu-xu*zdc)
                            +4.0*(xt*zca-zt*xca)*(zdb*xu-xdb*zu)+ ru2*2.0*(-xba*xca-zba*zca)))/rtru322
                        - ( (xdc*zt-xt*zdc+zca*xu-zu*xca)*(rt2*2.0*(zdb*xu-xdb*zu)+ru2*2.0*(zt*xba-xt*zba)))/rtru322
                        - ( (zdb*xt-zt*xdb+xba*zu-xu*zba)*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca)))/rtru322
                        + (xba*xdc+xca*xdb+zba*zdc+zca*zdb)/rtru;
               dyibyic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( (xdb*yt-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-ydb*xu)+ru2*2.0*(xt*yba-yt*xba))/rtru322)*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term3 =  ( dot*dot*2.0*(xdb*yu-ydb*xu)/(rt2*ru2*ru2) - dot*2.0*(xdb*yt-xt*ydb+yba*xu-yu*xba)/rt2ru2
                          + dot*dot*2.0*(xt*yba-yt*xba)/(rt2*rt2*ru2))*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(xdb*yu-ydb*xu)+ru2*2.0*(xt*yba-yt*xba))*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(rt2*2.0*(ydb*zdc+xu) + 4.0*(xt*zca-zt*xca)*(xdb*yu-ydb*xu)
                          + 4.0*(xt*yba-yt*xba)*(xdc*zu-xu*zdc) + ru2*2.0*(yba*zca+xt)))/rtru322
                        - ( (xdc*zt-xt*zdc+zca*xu-zu*xca)*(rt2*2.0*(xdb*yu-ydb*xu)+ru2*2.0*(xt*yba-yt*xba)))/rtru322
                        - ( (xdb*yt-xt*ydb+yba*xu-yu*xba)*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca)))/rtru322
                        + (yba*zcb-ycb*zba-ydb*zca-ydc*zcb-yba*zdc+ycb*zdc)/rtru;
               dyibzic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( dphidxt/rtru - dot*rt2*dphidxu/(rtru*rt2ru2))*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term3 =  ( dot*dot*2.0*(zcb*yu-zu*ycb)/(rt2*ru2*ru2) - 2.0*dphidxt*dot/rt2ru2)*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);

               term2 =  3.0*2.0*rt2*(zcb*yu-zu*ycb)*dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( 2.0*rt2*(zcb*yu-zu*ycb)*(xdc*zt-xt*zdc+zca*xu-zu*xca))/rtru322
                        - ( dot*(rt2*2.0*(zu-xdc*ycb)+4.0*(zcb*yu-zu*ycb)*(zca*xt-zt*xca)))/rtru322
                        - ( dphidxt*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca)))/rtru322
                        + (xba*ycb+xca*ycb-xcb*yba)/rtru;
               dyibxid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( dphidyt/rtru - dot*rt2*dphidyu/(rtru*rt2ru2))*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term3 =  ( dot*dot*2.0*dphidyu/(rt2*ru2*ru2) - 2.0*dphidyt*dot/rt2ru2)*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term2 =  3.0*2.0*rt2*dphidyu*dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( 2.0*rt2*dphidyu*(xdc*zt-xt*zdc+zca*xu-zu*xca))/rtru322
                        - ( dot*(rt2*2.0*(xcb*xdc+zcb*zdc)+4.0*dphidyu*(zca*xt-zt*xca)))/rtru322
                        - ( dphidyt*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca)))/rtru322
                        + (-xca*xcb-zcb*zca)/rtru;
               dyibyid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                         
               term1 =  ( dphidzt/rtru - dot*rt2*dphidzu/(rtru*rt2ru2))*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term3 =  ( dot*dot*2.0*dphidzu/(rt2*ru2*ru2) - 2.0*dphidzt*dot/rt2ru2)*
                        ( (xdc*zt-xt*zdc+zca*xu-zu*xca)/rtru - dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/rtru322);
               term2 =  3.0*2.0*rt2*dphidzu*dot*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( 2.0*rt2*dphidzu*(xdc*zt-xt*zdc+zca*xu-zu*xca))/rtru322
                        - ( dot*(rt2*2.0*(-ycb*zdc-xu)+4.0*dphidzu*(zca*xt-zt*xca)))/rtru322
                        - ( dphidzt*(rt2*2.0*(xdc*zu-xu*zdc)+ru2*2.0*(xt*zca-zt*xca)))/rtru322
                        + (ycb*zba+ycb*zca-yba*zcb)/rtru;
               dyibzid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( (zt*ydb-yt*zdb+zba*yu-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/rtru322)*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term3 =  ( dot*dot*2.0*(ydb*zu-yu*zdb)/(rt2*ru2*ru2) - dot*2.0*(zt*ydb-yt*zdb+zba*yu-zu*yba)/rt2ru2
                        + dot*dot*2.0*(yt*zba-zt*yba)/(rt2*rt2*ru2))*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( dot*(rt2*2.0*(xdc*zdb+yu) + 4.0*(xca*yt-xt*yca)*(ydb*zu-yu*zdb)
                                                     + 4.0*(zba*yt-zt*yba)*(ydc*xu-yu*xdc) + ru2*2.0*(xca*zba+yt)))/rtru322
                       - ( (xt*ydc-yt*xdc+yu*xca-xu*yca)*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/rtru322
                       - ( (zt*ydb-yt*zdb+zba*yu-zu*yba)*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))) /rtru322
                       + (xcb*zba-xdc*zba-xba*zcb+xdc*zcb-xca*zdb-xcb*zdc)/rtru;
               dzibxic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( (xt*zdb-zt*xdb+xba*zu-xu*zba)/rtru - dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/rtru322)*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term3 =  ( dot*dot*2.0*(zdb*xu-zu*xdb)/(rt2*ru2*ru2) - dot*2.0*(xt*zdb-zt*xdb+xba*zu-xu*zba)/rt2ru2
                        + dot*dot*2.0*(zt*xba-xt*zba)/(rt2*rt2*ru2))*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( dot*(rt2*2.0*(ydc*zdb-xu) + 4.0*(xba*zt-xt*zba)*(ydcxu)
                                                     + 4.0*(xcayt)*(zdb*xu-zu*xdb) + ru2*2.0*(yca*zba-xt)))/rtru322
                       - ( (xt*zdb-zt*xdb+xba*zu-xu*zba)*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt)))/rtru322
                       - ( (xt*ydc-yt*xdc+yu*xca-xu*yca)*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))) /rtru322
                       + (ycb*zba-ydc*zba-yba*zcb+ydc*zcb-yca*zdb-ycb*zdc)/rtru;
               dzibyic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( (yt*xdb-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322)*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term3 =  ( dot*dot*2.0*(xdb*yu-xu*ydb)/(rt2*ru2*ru2) - dot*2.0*(yt*xdb-xt*ydb+yba*xu-yu*xba)/rt2ru2
                        + dot*dot*2.0*(xt*yba-yt*xba)/(rt2*rt2*ru2))*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( dot*(rt2*2.0*(-xdc*xdb-ydb*ydc) + 4.0*(xcayt)*(xdb*yu-xu*ydb)
                                                          + 4.0*(yba*xt-yt*xba)*(ydcxu) - ru2*2.0*(xca*xba+yba*yca)))/rtru322
                       - ( (yt*xdb-xt*ydb+yba*xu-yu*xba)*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt)))/rtru322
                       - ( (xt*ydc-yt*xdc+yu*xca-xu*yca)*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/rtru322
                       + (xca*xdb+xba*xdc+yca*ydb+yba*ydc)/rtru;
               dzibzic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//               
               term1 = ( dphidxt/rtru - dot*rt2*dphidxu/(rtru*rt2ru2))*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term3 = (dphidxu*2.0*dot*dot/(rt2*ru2*ru2) - 2.0*dphidxt*dot/rt2ru2)*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term2 =  3.0*2.0*rt2*dphidxu*dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidxu*(xt*ydc-yt*xdc+yu*xca-xu*yca))/rtru322
                       - ( dot*(rt2*(2.0*(-xdc*zcb)-2.0*yu)+4.0*dphidxu*(xca*yt-xt*yca)))/rtru322
                       - ( dphidxt*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt)))/rtru322
                       + (xca*zcb + xba*zcb-xcb*zba)/rtru;
               dzibxid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidyt/rtru - dot*rt2*dphidyu/(rtru*rt2ru2))*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term3 = (dphidyu*2.0*dot*dot/(rt2*ru2*ru2) - 2.0*dphidyt*dot/rt2ru2)*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);

               term2 =  3.0*2.0*rt2*dphidyu*dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidyu*(xt*ydc-yt*xdc+yu*xca-xu*yca))/rtru322
                       - ( dot*(rt2*2.0*(-ydc*zcb+xu)+4.0*dphidyu*(xca*yt-xt*yca)))/rtru322
                       - ( dphidyt*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt)))/rtru322
                       + (yba*zcb + yca*zcb-ycb*zba)/rtru;
               dzibyid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 = ( dphidzt/rtru - dot*rt2*dphidzu/(rtru*rt2ru2))*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term3 = (dphidzu*2.0*dot*dot/(rt2*ru2*ru2) - 2.0*dphidzt*dot/rt2ru2)*
                        ( (xt*ydc-yt*xdc+yu*xca-xu*yca)/rtru - dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt) )/rtru322);
               term2 =  3.0*2.0*rt2*dphidzu*dot*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidzu*(xt*ydc-yt*xdc+yu*xca-xu*yca))/rtru322
                       - ( dot*(rt2*2.0*(xdc*xcb+ycb*ydc)+4.0*dphidzu*(xca*yt-xt*yca)))/rtru322
                       - ( dphidzt*(rt2*2.0*(ydcxu)+ru2*2.0*(xcayt)))/rtru322
                       - (xca*xcb+yca*ycb)/rtru;
               dzibzid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  ( (zt*ydb-yt*zdb+zba*yu-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/rtru322)*
                        ( (xt*zdb-zt*xdb+xba*zu-xu*zba)/rtru - dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/rtru322);
               term3 =  ( dot*dot*2.0*(zdb*xu-zu*xdb)/(rt2*ru2*ru2) - dot*2.0*(xt*zdb-zt*xdb+xba*zu-xu*zba)/rt2ru2
                        + dot*dot*2.0*(xba*zt-xt*zba)/(rt2*rt2*ru2))*
                        ((zt*ydb-yt*zdb+zba*yu-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(4.0*(xba*zt-xt*zba)*(ydb*zu-yu*zdb)-rt2*2.0*ydb*xdb+4.0*(zba*yt-zt*yba)*(zdb*xu-zu*xdb)-ru2*2.0*xba*yba))/rtru322
                        - ( (xt*zdb-zt*xdb+xba*zu-xu*zba)*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/rtru322
                        - ( (zt*ydb-yt*zdb+zba*yu-zu*yba)*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(xba*zt-xt*zba)))/rtru322
                        + (xba*ydb+yba*xdb)/rtru;
               dxicyic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                         
               term1 =  ( (zt*ydb-yt*zdb+zba*yu-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/rtru322)*
                        ( (yt*xdb-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);
               term3 =  ( dot*dot*2.0*(xdb*yu-xu*ydb)/(rt2*ru2*ru2) - dot*2.0*(yt*xdb-xt*ydb+yba*xu-yu*xba)/rt2ru2
                        + dot*dot*2.0*(yba*xt-yt*xba)/(rt2*rt2*ru2))*
                        ((zt*ydb-yt*zdb+zba*yu-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))/rtru322);
               term2 =  3.0*dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba))*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( dot*(4.0*(yba*xt-yt*xba)*(ydb*zu-yu*zdb)-rt2*2.0*zdb*xdb+4.0*(zba*yt-zt*yba)*(xdb*yu-xu*ydb)-ru2*2.0*xba*zba))/rtru322
                        - ( (yt*xdb-xt*ydb+yba*xu-yu*xba)*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(zba*yt-zt*yba)))/rtru322
                        - ( (zt*ydb-yt*zdb+zba*yu-zu*yba)*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/rtru322
                        + (xba*zdb+zba*xdb)/rtru;
               dxiczic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  (dphidxt/rtru - rt2*dot*dphidxu/(rtru*rt2ru2))*
                        ( (zt*ydb-yt*zdb+yu*zba-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/rtru322);
               term3 =  ( 2.0*dphidxu*dot*dot/(rt2*ru2*ru2) - 2.0*dphidxt*dot/rt2ru2)*
                        ( (zt*ydb-yt*zdb+yu*zba-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/rtru322);
               term2 =  3.0*2.0*rt2*dphidxu*dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidxu*(zt*ydb-yt*zdb+yu*zba-zu*yba))/rtru322
                       - ( dot*(rt2*2.0*(-ydb*ycb-zdb*zcb)+4.0*dphidxu*(zba*yt-zt*yba)))/rtru322
                       - ( dphidxt*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba)))/rtru322
                       + (zba*zcb+yba*ycb)/rtru;
               dxicxid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 =  (dphidyt/rtru - rt2*dot*dphidyu/(rtru*rt2ru2))*
                        ( (zt*ydb-yt*zdb+yu*zba-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/rtru322);
               term3 =  ( 2.0*dphidyu*dot*dot/(rt2*ru2*ru2) - 2.0*dphidyt*dot/rt2ru2)*
                        ( (zt*ydb-yt*zdb+yu*zba-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/rtru322);
               term2 =  3.0*2.0*rt2*dphidyu*dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidyu*(zt*ydb-yt*zdb+yu*zba-zu*yba))/rtru322
                       - ( dot*(rt2*2.0*(xcb*ydb+zu)+4.0*dphidyu*(zba*yt-zt*yba)))/rtru322
                       - ( dphidyt*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba)))/rtru322
                       + (xba*ycb-2.0*yba*xcb)/rtru;
               dxicyid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                      
               term1 =  (dphidzt/rtru - rt2*dot*dphidzu/(rtru*rt2ru2))*
                        ( (zt*ydb-yt*zdb+yu*zba-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/rtru322);
               term3 =  ( 2.0*dphidzu*dot*dot/(rt2*ru2*ru2) - 2.0*dphidzt*dot/rt2ru2)*
                        ( (zt*ydb-yt*zdb+yu*zba-zu*yba)/rtru - dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/rtru322);
               term2 =  3.0*2.0*rt2*dphidzu*dot*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidzu*(zt*ydb-yt*zdb+yu*zba-zu*yba))/rtru322
                       - ( dot*(rt2*2.0*(zdb*xcb-yu)+4.0*dphidzu*(zba*yt-zt*yba)))/rtru322
                       - ( dphidzt*(rt2*2.0*(ydb*zu-yu*zdb)+ru2*2.0*(yt*zba-zt*yba)))/rtru322
                       + (xba*zcb-2.0*zba*xcb)/rtru;
               dxiczid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                         
               term1 = ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322)*
                       ((yt*xdb-xt*ydb+xu*yba-yu*xba)/rtru - (dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(xt*yba-yt*xba)))/rtru322 );            
               term3 = ( dot*dot*2.0*(xdb*yu-xu*ydb)/(rt2*ru2*ru2) - dot*2.0*(yt*xdb-xt*ydb+xu*yba-yu*xba)/rt2ru2
                       + dot*dot*2.0*(yba*xt-yt*xba)/(rt2*rt2*ru2))*
                       ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322);
               term2 = 3.0*dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(xt*yba-yt*xba))*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba))/(4.0*rt2ru2*rt2ru2*rtru)
                      - ( dot*(-rt2*2.0*zdb*ydb+4.0*(xba*zt-xt*zba)*(xdb*yu-xu*ydb)-ru2*2.0*yba*zba
                                              + 4.0*(yba*xt-yt*xba)*(zdb*xu-zu*xdb)))/rtru322
                      - ( (xt*zdb-zt*xdb+zu*xba-xu*zba)*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(xt*yba-yt*xba)) )/rtru322
                      - ( (xdb*yt-xt*ydb+yba*xu-yu*xba)*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)) )/rtru322
                      + (yba*zdb+zba*ydb)/rtru;
               dyiczic = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//
               term1 = ( dphidxt/rtru - rt2*dot*dphidxu/(rtru*rt2ru2))*
                       ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322);
               term3 = ( 2.0*dot*dot*dphidxu/(rt2*ru2*ru2) - 2.0*dphidxt*dot/rt2ru2)*
                       ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322);
               term2 = ( 3.0*2.0*rt2*dphidxu*dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)) )/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*(zcb*yu-zu*ycb)*(xt*zdb-zt*xdb+xba*zu-zba*xu))/rtru322
                       - ( dot*(rt2*2.0*(xdb*ycb-zu)+4.0*(xba*zt-xt*zba)*(zcb*yu-zu*ycb)))/rtru322
                       - ( dphidxt*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322
                       + (xcb*yba-2.0*xba*ycb)/rtru ;
               dyicxid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3); 
//                                       
               term1 = ( dphidyt/rtru - rt2*dot*dphidyu/(rtru*rt2ru2))*
                       ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322);
               term3 = ( 2.0*dot*dot*dphidyu/(rt2*ru2*ru2) - 2.0*dphidyt*dot/rt2ru2)*
                       ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322);

               term2 = ( 3.0*2.0*rt2*dphidyu*dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)) )/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidyu*(xt*zdb-zt*xdb+xba*zu-zba*xu))/rtru322
                       - ( dot*(rt2*2.0*(-xdb*xcb-zdb*zcb)+4.0*(xba*zt-xt*zba)*(xcb*zu-xu*zcb)))/rtru322
                       - ( dphidyt*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322
                       + (xcb*xba+zba*zcb)/rtru ;
               dyicyid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3); 
//               
               term1 = ( dphidzt/rtru - rt2*dot*dphidzu/(rtru*rt2ru2))*
                       ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322);
               term3 = ( 2.0*dot*dot*dphidzu/(rt2*ru2*ru2) - 2.0*dphidzt*dot/rt2ru2)*
                       ((xt*zdb-zt*xdb+zu*xba-xu*zba)/rtru - (dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322);
               term2 = ( 3.0*2.0*rt2*dphidzu*dot*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)) )/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*(ycb*xu-yu*xcb)*(xt*zdb-zt*xdb+xba*zu-zba*xu))/rtru322
                       - ( dot*(rt2*2.0*(zdb*ycb+xu)+4.0*(xba*zt-xt*zba)*(ycb*xu-yu*xcb)))/rtru322
                       - ( dphidzt*(rt2*2.0*(zdb*xu-zu*xdb)+ru2*2.0*(zt*xba-xt*zba)))/rtru322
                       + (zcb*yba-2.0*zba*ycb)/rtru ;
               dyiczid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3); 
//
               term1 = ( dphidxt/rtru - rt2*dphidxu*dot/(rtru*rt2ru2))*(dphidyt/rtru - rt2*dphidyu*dot/(rtru*rt2ru2));
               term3 = ( 2.0*dphidyu*dot*dot/(rt2*ru2*ru2) - 2.0*dphidyt*dot/rt2ru2)*
                       ( dphidxt/rtru - rt2*dphidxu*dot/(rtru*rt2ru2));
               term2 =  3.0*4.0*rt2*rt2*dphidxu*dphidyu*dot/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( 2.0*dphidyt*rt2*dphidxu)/rtru322
                        - ( 2.0*dphidxt*rt2*dphidyu)/rtru322
                        - ( -xcb*ycb*rt2*dot)/rtru32;
               dxidyid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                    
               term1 = ( dphidxt/rtru - rt2*dphidxu*dot/(rtru*rt2ru2))*(dphidzt/rtru - rt2*dphidzu*dot/(rtru*rt2ru2));
               term3 = ( 2.0*dphidzu*dot*dot/(rt2*ru2*ru2) - 2.0*dphidzt*dot/rt2ru2)*(dphidxt/rtru - rt2*dphidxu*dot/(rtru*rt2ru2));
               term2 =  3.0*4.0*rt2*rt2*dphidxu*dphidzu*dot/(4.0*rt2ru2*rt2ru2*rtru)
                        - ( 2.0*dphidzt*rt2*dphidxu)/rtru322
                        - ( 2.0*dphidxt*rt2*dphidzu)/rtru322
                        - ( -xcb*zcb*rt2*dot)/(rtru*rt2ru2);
               dxidzid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);               
//
               term1 = ( dphidzt/rtru - rt2*dphidzu*dot/(rtru*rt2ru2))*(dphidyt/rtru - rt2*dphidyu*dot/(rtru*rt2ru2));
               term3 = ( 2.0*dphidzu*dot*dot/(rt2*ru2*ru2) - 2.0*dphidzt*dot/rt2ru2)*(dphidyt/rtru - rt2*dphidyu*dot/(rtru*rt2ru2));

               term2 =  3.0*4.0*rt2*rt2*dphidzu*dphidyu*dot/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( 2.0*dphidyt*rt2*dphidzu)/rtru322
                       - ( 2.0*dphidzt*rt2*dphidyu)/rtru322
                       - ( -dot*rt2*ycb*zcb)/(rt2ru2*rtru);
               dyidzid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);                         
//
               term1 =  ( dphidxt/rtru - rt2*dot*dphidxu/(rtru*rt2ru2))*
                        ( (xdb*yt-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);

               term3 = ( (2.0*dphidxu*dot*dot)/(rt2*ru2*ru2) - (2.0*dphidxt*dot)/rt2ru2)*
                        ( (xdb*yt-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);

               term2 = ( (3.0*2.0*dphidxu*dot*rt2*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*dphidxu*(xdb*yt-xt*ydb+yba*xu-yu*xba))/rtru322
                       - ( dot*(rt2*2.0*(xdb*zcb+yu)+4.0*dphidxu*(yba*xt-yt*xba)))/rtru322
                       - ( dphidxt*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/rtru322
                       + (xcb*zba-2.0*xba*zcb)/(rtru));
               dzicxid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                                      
               term1 =  ( dphidyt/rtru - rt2*dot*dphidyu/(rtru*rt2ru2))*
                        ( (xdb*yt-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);

               term3 = ( (2.0*dphidyu*dot*dot)/(rt2*ru2*ru2) - (2.0*dphidyt*dot)/rt2ru2)*
                        ( (xdb*yt-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);

               term2 = ( (3.0*2.0*dphidyu*dot*rt2*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*(xdb*yt-xt*ydb+yba*xu-yu*xba)*dphidyu)/rtru322
                       - ( dot*(rt2*2.0*(ydb*zcb-xu)+4.0*(yba*xt-yt*xba)*dphidyu))/rtru322
                       - ( (xcb*zt-xt*zcb)*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/rtru322
                       + (ycb*zba-2.0*yba*zcb)/(rtru));
               dzicyid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);
//                                     
               term1 =  ( dphidzt/rtru - rt2*dot*dphidzu/(rtru*rt2ru2))*
                        ( (xdb*yt-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);

               term3 = ( (2.0*dphidzu*dot*dot)/(rt2*ru2*ru2) - (2.0*dphidzt*dot)/rt2ru2)*
                        ( (xdb*yt-xt*ydb+yba*xu-yu*xba)/rtru - dot*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba))/rtru322);

               term2 = ( (3.0*2.0*dphidzu*dot*rt2*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/(4.0*rt2ru2*rt2ru2*rtru)
                       - ( rt2*2.0*(xdb*yt-xt*ydb+yba*xu-yu*xba)*dphidzu)/rtru322
                       - ( dot*(rt2*2.0*(-xcb*xdb-ycb*ydb)+4.0*(yba*xt-yt*xba)*dphidzu))/rtru322
                       - ( dphidzt*(rt2*2.0*(xdb*yu-xu*ydb)+ru2*2.0*(yba*xt-yt*xba)))/rtru322
                       + (xba*xcb+yba*ycb)/(rtru));
               dziczid = (d2edphi2*term1)/(denom2) + dedphi*term2/denom - dedphi*term3/(2.0*denom3);              



         if (iatom == ia)
         {
            hess.hessx[ia][0] += dxiaxia ;
            hess.hessy[ia][0] += dxiayia ;
            hess.hessz[ia][0] += dxiazia ;
            hess.hessx[ia][1] += dxiayia ;
            hess.hessy[ia][1] += dyiayia ;
            hess.hessz[ia][1] += dyiazia ;
            hess.hessx[ia][2] += dxiazia ;
            hess.hessy[ia][2] += dyiazia ;
            hess.hessz[ia][2] += dziazia ;
            hess.hessx[ib][0] += dxiaxib ;
            hess.hessy[ib][0] += dyiaxib ;
            hess.hessz[ib][0] += dziaxib ;
            hess.hessx[ib][1] += dxiayib ;
            hess.hessy[ib][1] += dyiayib ;
            hess.hessz[ib][1] += dziayib ;
            hess.hessx[ib][2] += dxiazib ;
            hess.hessy[ib][2] += dyiazib ;
            hess.hessz[ib][2] += dziazib ;
            hess.hessx[ic][0] += dxiaxic ;
            hess.hessy[ic][0] += dyiaxic ;
            hess.hessz[ic][0] += dziaxic ;
            hess.hessx[ic][1] += dxiayic ;
            hess.hessy[ic][1] += dyiayic ;
            hess.hessz[ic][1] += dziayic ;
            hess.hessx[ic][2] += dxiazic ;
            hess.hessy[ic][2] += dyiazic ;
            hess.hessz[ic][2] += dziazic ;
            hess.hessx[id][0] += dxiaxid ;
            hess.hessy[id][0] += dyiaxid ;
            hess.hessz[id][0] += dziaxid ;
            hess.hessx[id][1] += dxiayid ;
            hess.hessy[id][1] += dyiayid ;
            hess.hessz[id][1] += dziayid ;
            hess.hessx[id][2] += dxiazid ;
            hess.hessy[id][2] += dyiazid ;
            hess.hessz[id][2] += dziazid ;
         }else if (iatom == ib)
         {
            hess.hessx[ib][0] += dxibxib ;
            hess.hessy[ib][0] += dxibyib ;
            hess.hessz[ib][0] += dxibzib ;
            hess.hessx[ib][1] += dxibyib ;
            hess.hessy[ib][1] += dyibyib ;
            hess.hessz[ib][1] += dyibzib ;
            hess.hessx[ib][2] += dxibzib ;
            hess.hessy[ib][2] += dyibzib ;
            hess.hessz[ib][2] += dzibzib ;
            hess.hessx[ia][0] += dxiaxib ;
            hess.hessy[ia][0] += dxiayib ;
            hess.hessz[ia][0] += dxiazib ;
            hess.hessx[ia][1] += dyiaxib ;
            hess.hessy[ia][1] += dyiayib ;
            hess.hessz[ia][1] += dyiazib ;
            hess.hessx[ia][2] += dziaxib ;
            hess.hessy[ia][2] += dziayib ;
            hess.hessz[ia][2] += dziazib ;
            hess.hessx[ic][0] += dxibxic ;
            hess.hessy[ic][0] += dyibxic ;
            hess.hessz[ic][0] += dzibxic ;
            hess.hessx[ic][1] += dxibyic ;
            hess.hessy[ic][1] += dyibyic ;
            hess.hessz[ic][1] += dzibyic ;
            hess.hessx[ic][2] += dxibzic ;
            hess.hessy[ic][2] += dyibzic ;
            hess.hessz[ic][2] += dzibzic ;
            hess.hessx[id][0] += dxibxid ;
            hess.hessy[id][0] += dyibxid ;
            hess.hessz[id][0] += dzibxid ;
            hess.hessx[id][1] += dxibyid ;
            hess.hessy[id][1] += dyibyid ;
            hess.hessz[id][1] += dzibyid ;
            hess.hessx[id][2] += dxibzid ;
            hess.hessy[id][2] += dyibzid ;
            hess.hessz[id][2] += dzibzid ;
         }else if (iatom == ic)
         {
            hess.hessx[ic][0] += dxicxic ;
            hess.hessy[ic][0] += dxicyic ;
            hess.hessz[ic][0] += dxiczic ;
            hess.hessx[ic][1] += dxicyic ;
            hess.hessy[ic][1] += dyicyic ;
            hess.hessz[ic][1] += dyiczic ;
            hess.hessx[ic][2] += dxiczic ;
            hess.hessy[ic][2] += dyiczic ;
            hess.hessz[ic][2] += dziczic ;
            hess.hessx[ia][0] += dxiaxic ;
            hess.hessy[ia][0] += dxiayic ;
            hess.hessz[ia][0] += dxiazic ;
            hess.hessx[ia][1] += dyiaxic ;
            hess.hessy[ia][1] += dyiayic ;
            hess.hessz[ia][1] += dyiazic ;
            hess.hessx[ia][2] += dziaxic ;
            hess.hessy[ia][2] += dziayic ;
            hess.hessz[ia][2] += dziazic ;
            hess.hessx[ib][0] += dxibxic ;
            hess.hessy[ib][0] += dxibyic ;
            hess.hessz[ib][0] += dxibzic ;
            hess.hessx[ib][1] += dyibxic ;
            hess.hessy[ib][1] += dyibyic ;
            hess.hessz[ib][1] += dyibzic ;
            hess.hessx[ib][2] += dzibxic ;
            hess.hessy[ib][2] += dzibyic ;
            hess.hessz[ib][2] += dzibzic ;
            hess.hessx[id][0] += dxicxid ;
            hess.hessy[id][0] += dyicxid ;
            hess.hessz[id][0] += dzicxid ;
            hess.hessx[id][1] += dxicyid ;
            hess.hessy[id][1] += dyicyid ;
            hess.hessz[id][1] += dzicyid ;
            hess.hessx[id][2] += dxiczid ;
            hess.hessy[id][2] += dyiczid ;
            hess.hessz[id][2] += dziczid ;
         }else if (iatom == id)
         {
            hess.hessx[id][0] += dxidxid ;
            hess.hessy[id][0] += dxidyid ;
            hess.hessz[id][0] += dxidzid ;
            hess.hessx[id][1] += dxidyid ;
            hess.hessy[id][1] += dyidyid ;
            hess.hessz[id][1] += dyidzid ;
            hess.hessx[id][2] += dxidzid ;
            hess.hessy[id][2] += dyidzid ;
            hess.hessz[id][2] += dzidzid ;
            hess.hessx[ia][0] += dxiaxid ;
            hess.hessy[ia][0] += dxiayid ;
            hess.hessz[ia][0] += dxiazid ;
            hess.hessx[ia][1] += dyiaxid ;
            hess.hessy[ia][1] += dyiayid ;
            hess.hessz[ia][1] += dyiazid ;
            hess.hessx[ia][2] += dziaxid ;
            hess.hessy[ia][2] += dziayid ;
            hess.hessz[ia][2] += dziazid ;
            hess.hessx[ib][0] += dxibxid ;
            hess.hessy[ib][0] += dxibyid ;
            hess.hessz[ib][0] += dxibzid ;
            hess.hessx[ib][1] += dyibxid ;
            hess.hessy[ib][1] += dyibyid ;
            hess.hessz[ib][1] += dyibzid ;
            hess.hessx[ib][2] += dzibxid ;
            hess.hessy[ib][2] += dzibyid ;
            hess.hessz[ib][2] += dzibzid ;
            hess.hessx[ic][0] += dxicxid ;
            hess.hessy[ic][0] += dxicyid ;
            hess.hessz[ic][0] += dxiczid ;
            hess.hessx[ic][1] += dyicxid ;
            hess.hessy[ic][1] += dyicyid ;
            hess.hessz[ic][1] += dyiczid ;
         }
         }
        }
    }
}
