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
#include "gmmx.h"
#include "field.h"
#include "fix.h"
#include "units.h"
#include "minim_values.h"
#include "eangle.h"
        
void eangle()
{
    int i, ia, ib, ic, id;
    double e, angle, dot, cosine, factor;
    double dt, dt2, dt3, dt4;
    double xia, yia, zia;
    double xib, yib, zib;
    double xic, yic, zic;
    double xid, yid, zid;
    double xab, yab, zab, rab2;
    double xcb, ycb, zcb, rcb2;
    double xad, yad, zad;
    double xbd, ybd, zbd;
    double xcd,ycd,zcd;
    double xip,yip,zip;
    double xap,yap,zap,rap2;
    double xcp,ycp,zcp,rcp2;
    double xt,yt,zt,rt2,delta;



    if (minim_values.iprint)
    {
        fprintf(pcmoutfile,"\nAngle Terms \n");
        fprintf(pcmoutfile,"             At1      At2      At3     Angle     Thet0   Tconst      Ebend\n");
    }
    gmmx.hybrida = FALSE;

    for (i=0; i < angles.nang; i++)
    {
        ia = angles.i13[i][0];
        ib = angles.i13[i][1];
        ic = angles.i13[i][2];
        id = angles.i13[i][3];
        if (atom[ia].use || atom[ib].use || atom[ic].use)
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
            if ( (angles.angin[i] == FALSE) || field.type == MMFF94 )
            {
                xab = xia - xib;
                yab = yia - yib;
                zab = zia - zib;
                rab2 = xab*xab + yab*yab + zab*zab;
                xcb = xic - xib;
                ycb = yic - yib;
                zcb = zic - zib;
                rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
                if (rab2 != 0.0 && rcb2 != 0.00)
                {
                    dot = xab*xcb + yab*ycb + zab*zcb;
                    cosine = dot / sqrt(rab2*rcb2);
                    if (cosine > 1.0)
                       cosine = 1.0;
                    if (cosine < -1.0)
                       cosine = -1.0;
                    angle = radian*acos(cosine);

                    if (angle < 0.75*angles.anat[i])
                      if ( atom[ib].mmx_type != 80) 
                         gmmx.hybrida = TRUE;

                    if (angles.angtype[i] == HARMONIC)
                    {
                        dt = angle - angles.anat[i];
                        dt2 = dt*dt;
                        dt3 = dt2*dt;
                        dt4 = dt2*dt2;
                        e = units.angunit*angles.acon[i]*dt2
                          *(1.0 + units.cang*dt + units.qang*dt2 + units.pang*dt3 + units.sang*dt4);
                    } else
                    {
                        factor = angle/radian;
                        e = units.angunit*angles.fcon[i]*
                          (angles.c0[i] + angles.c1[i]*cosine + angles.c2[i]*cos(2*factor)
                          + angles.c3[i]*cos(3*factor) + angles.c4[i]*cos(4*factor)
                          + angles.c5[i]*cos(5*factor) + angles.c6[i]*cos(6*factor));
                    }
                    energies.ebend += e;
                    atom[ia].energy += e;
                    atom[ib].energy += e;
                    atom[ic].energy += e;
                    if (minim_values.iprint)
                    {
                       if (angles.angtype[i] == HARMONIC)
                          fprintf(pcmoutfile,"Angle 1-3: %2s(%-3d)- %2s(%-3d)- %2s(%-3d)  %-8.3f  %-8.3f %-8.4f = %-8.4f\n",
                             atom[ia].name,ia,atom[ib].name, ib,atom[ic].name, ic, angle, angles.anat[i],angles.acon[i], e);
                       else
                          fprintf(pcmoutfile,"Angle 1-3: %2s(%-3d)- %2s(%-3d)- %2s(%-3d)  %-8.3f         %-8.4f = %-8.4f\n",
                             atom[ia].name,ia,atom[ib].name, ib,atom[ic].name, ic, angle, angles.fcon[i], e);
                    }
                }  
            } else
            {
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
                xip = xib + xt*delta;
                yip = yib + yt*delta;
                zip = zib + zt*delta;
                xap = xia - xip;
                yap = yia - yip;
                zap = zia - zip;
                rap2 = xap*xap + yap*yap + zap*zap;
                xcp = xic - xip;
                ycp = yic - yip;
                zcp = zic - zip;
                rcp2 = xcp*xcp + ycp*ycp + zcp*zcp;
                if (rap2 != 0.0 && rcp2 != 0.0)
                {
                    dot = xap*xcp + yap*ycp + zap*zcp;
                    cosine = dot/ sqrt(rap2*rcp2);
                    if (cosine < -1.0)
                      cosine = -1.0;
                    if (cosine > 1.0)
                      cosine = 1.0;
                    angle = radian*acos(cosine);
                    dt = angle - angles.anat[i];
                    dt2 = dt*dt;
                    dt3 = dt2*dt;
                    dt4 = dt2*dt2;
                    e = units.angunit*angles.acon[i]*dt2
                      *(1.0 + units.cang*dt + units.qang*dt2 + units.pang*dt3 + units.sang*dt4);
                    energies.ebend += e;
                   atom[ia].energy += e;
                   atom[ib].energy += e;
                   atom[ic].energy += e;
                    if (minim_values.iprint)
                     fprintf(pcmoutfile,"Angle InP: %2s(%-3d)- %2s(%-3d)- %2s(%-3d)  %-8.3f  %-8.3f %-8.4f = %-8.4f\n",
                       atom[ia].name, ia, atom[ib].name, ib, atom[ic].name, ic, angle, angles.anat[i], angles.acon[i], e);

                }
            }
        }         
    }
     
}

void eangle1()
{
    int i, ia, ib, ic, id;
    double temp1, temp2;
    double e, angle, dot, cosine, factor, sine;
    double dt, dt2, dt3, dt4;
    double deddt,terma,termc,term;
    double xia, yia, zia;
    double xib, yib, zib;
    double xic, yic, zic;
    double xid, yid, zid;
    double xab, yab, zab, rab2;
    double xcb, ycb, zcb, rcb2;
    double xp,yp,zp,rp;
    double xad, yad, zad;
    double xbd, ybd, zbd;
    double xcd,ycd,zcd;
    double xip,yip,zip;
    double xap,yap,zap,rap2;
    double xcp,ycp,zcp,rcp2;
    double xt,yt,zt,rt2,ptrt2;
    double xm,ym,zm,rm,delta,delta2;
    double dedxia,dedyia,dedzia;
    double dedxib,dedyib,dedzib;
    double dedxic,dedyic,dedzic;
    double dedxid,dedyid,dedzid;
    double dedxip,dedyip,dedzip;
    double dpdxia,dpdyia,dpdzia;
    double dpdxic,dpdyic,dpdzic;
    double dtemp;

    energies.ebend = 0.0F;
      for (i=0; i <= natom; i++)
      {
          deriv.dea[i][0] = 0.0;
          deriv.dea[i][1] = 0.0;
          deriv.dea[i][2] = 0.0;
      }
    gmmx.hybrida = FALSE;
    for (i=0; i < angles.nang; i++)
    {
        ia = angles.i13[i][0];
        ib = angles.i13[i][1];
        ic = angles.i13[i][2];
        id = angles.i13[i][3];
        if (atom[ia].use || atom[ib].use || atom[ic].use)
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
            if ( (angles.angin[i] == FALSE) || field.type == MMFF94 )
            {
                xab = xia - xib;
                yab = yia - yib;
                zab = zia - zib;
                rab2 = xab*xab + yab*yab + zab*zab;
                xcb = xic - xib;
                ycb = yic - yib;
                zcb = zic - zib;
                rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
                xp = ycb*zab - zcb*yab;
                yp = zcb*xab - xcb*zab;
                zp = xcb*yab - ycb*xab;
                rp = sqrt(xp*xp + yp*yp + zp*zp);
                if (rp != 0.0 )
                {
                    dot = xab*xcb + yab*ycb + zab*zcb;
                    cosine = dot / sqrt(rab2*rcb2);
                    if (cosine > 1.0)
                       cosine = 1.0;
                    if (cosine < -1.0)
                       cosine = -1.0;
                    angle = radian*acos(cosine);
                    if (angle < 0.75*angles.anat[i])
                      if ( atom[ib].mmx_type != 80) 
                         gmmx.hybrida = TRUE;
                    if (angles.angtype[i] == HARMONIC)
                    {
                        dt = angle - angles.anat[i];
                        dt2 = dt*dt;
                        dt3 = dt2*dt;
                        dt4 = dt2*dt2;

                        e = units.angunit*angles.acon[i]*dt2
                          *(1.0 + units.cang*dt + units.qang*dt2 + units.pang*dt3 + units.sang*dt4);
                        deddt = units.angunit*angles.acon[i]* dt * radian
                            * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                                 + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);
                    } else
                    {
                        factor = angle/radian;
                        sine = sin(factor);
                        e = units.angunit*angles.fcon[i]*
                          (angles.c0[i] + angles.c1[i]*cosine + angles.c2[i]*cos(2*factor)
                          + angles.c3[i]*cos(3*factor) + angles.c4[i]*cos(4*factor)
                          + angles.c5[i]*cos(5*factor) + angles.c6[i]*cos(6*factor));
                        deddt = -units.angunit*angles.fcon[i]*
                          (angles.c1[i]*sine + 2*angles.c2[i]*sin(2*factor)
                          + 3*angles.c3[i]*sin(3*factor) + 4*angles.c4[i]*sin(4*factor)
                          + 5*angles.c5[i]*sin(5*factor) + 6*angles.c6[i]*sin(6*factor));
                    }

/*                   terma = -deddt / (rab2*rp);
                   termc = deddt / (rcb2*rp);
                   dedxia = terma * (yab*zp-zab*yp);
                   dedyia = terma * (zab*xp-xab*zp);
                   dedzia = terma * (xab*yp-yab*xp);
                   dedxic = termc * (ycb*zp-zcb*yp);
                   dedyic = termc * (zcb*xp-xcb*zp);
                   dedzic = termc * (xcb*yp-ycb*xp);
                   dedxib = -dedxia - dedxic;
                   dedyib = -dedyia - dedyic;
                   dedzib = -dedzia - dedzic;     */
         
                   dtemp = dot*dot/(rab2*rcb2);
                   if (fabs(dtemp) >= 1.0)
                     dtemp = 0.9999*dtemp/(fabs(dtemp));
                   terma = sqrt(1.0-dtemp);
                   termc = sqrt(rab2*rcb2);
                   dedxia = -deddt*( xcb/termc - (xab*cosine)/rab2)/terma;
                   dedyia = -deddt*( ycb/termc - (yab*cosine)/rab2)/terma;
                   dedzia = -deddt*( zcb/termc - (zab*cosine)/rab2)/terma;

                   dedxic = -deddt*( xab/termc - (xcb*cosine)/rcb2)/terma;
                   dedyic = -deddt*( yab/termc - (ycb*cosine)/rcb2)/terma;
                   dedzic = -deddt*( zab/termc - (zcb*cosine)/rcb2)/terma;

                   dedxib = -deddt*( (-xab-xcb)/termc + (dot*(xcb*rab2 + xab*rcb2))/(rab2*rcb2*termc))/terma;
                   dedyib = -deddt*( (-yab-ycb)/termc + (dot*(ycb*rab2 + yab*rcb2))/(rab2*rcb2*termc))/terma;
                   dedzib = -deddt*( (-zab-zcb)/termc + (dot*(zcb*rab2 + zab*rcb2))/(rab2*rcb2*termc))/terma;  
                  
                              
                   deriv.dea[ia][0] += dedxia;
                   deriv.dea[ia][1] += dedyia;
                   deriv.dea[ia][2] += dedzia;

                   deriv.dea[ib][0] += dedxib;
                   deriv.dea[ib][1] += dedyib;
                   deriv.dea[ib][2] += dedzib;

                   deriv.dea[ic][0] += dedxic;
                   deriv.dea[ic][1] += dedyic;
                   deriv.dea[ic][2] += dedzic;
                   
                   virial.virx += xab*dedxia + xcb*dedxic;
                   virial.viry += yab*dedyia + ycb*dedyic;
                   virial.virz += zab*dedzia + zcb*dedzic;

                   energies.ebend += e;
                }  
            } else
            {
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
                xip = xib + xt*delta;
                yip = yib + yt*delta;
                zip = zib + zt*delta;
                xap = xia - xip;
                yap = yia - yip;
                zap = zia - zip;
                rap2 = xap*xap + yap*yap + zap*zap;
                xcp = xic - xip;
                ycp = yic - yip;
                zcp = zic - zip;
                rcp2 = xcp*xcp + ycp*ycp + zcp*zcp;
                xm = ycp*zap - zcp*yap;
                ym = zcp*xap - xcp*zap;
                zm = xcp*yap - ycp*xap;
                rm = sqrt(xm*xm + ym*ym + zm*zm);
                if (rm != 0.0)
                {
                    dot = xap*xcp + yap*ycp + zap*zcp;
                    cosine = dot/ sqrt(rap2*rcp2);
                    if (cosine < -1.0)
                      cosine = 1.0;
                    if (cosine > 1.0)
                      cosine = 1.0;
                    temp1 = acos(cosine);
                    temp2 = radian;
                    angle = temp1*temp2;
                    dt = angle - angles.anat[i];
                    dt2 = dt*dt;
                    dt3 = dt2*dt;
                    dt4 = dt2*dt2;
                    e = units.angunit*angles.acon[i]*dt2
                      *(1.0 + units.cang*dt + units.qang*dt2 + units.pang*dt3 + units.sang*dt4);

                    deddt = units.angunit *angles.acon[i]* dt
                          * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                           + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);

                  deddt = deddt * radian;
                  terma = -deddt / (rap2*rm);
                  termc = deddt / (rcp2*rm);
                  dedxia = terma * (yap*zm-zap*ym);
                  dedyia = terma * (zap*xm-xap*zm);
                  dedzia = terma * (xap*ym-yap*xm);
                  dedxic = termc * (ycp*zm-zcp*ym);
                  dedyic = termc * (zcp*xm-xcp*zm);
                  dedzic = termc * (xcp*ym-ycp*xm);
                  dedxip = -dedxia - dedxic;
                  dedyip = -dedyia - dedyic;
                  dedzip = -dedzia - dedzic;

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

                  dedxia = dedxia + dpdxia;
                  dedyia = dedyia + dpdyia;
                  dedzia = dedzia + dpdzia;
                  dedxib = dedxip;
                  dedyib = dedyip;
                  dedzib = dedzip;
                  dedxic = dedxic + dpdxic;
                  dedyic = dedyic + dpdyic;
                  dedzic = dedzic + dpdzic;
                  dedxid = -dedxia - dedxib - dedxic;
                  dedyid = -dedyia - dedyib - dedyic;
                  dedzid = -dedzia - dedzib - dedzic;
                  
                   energies.ebend += e;
                   deriv.dea[ia][0] += dedxia;
                   deriv.dea[ia][1] += dedyia;
                   deriv.dea[ia][2] += dedzia;

                   deriv.dea[ib][0] += dedxib;
                   deriv.dea[ib][1] += dedyib;
                   deriv.dea[ib][2] += dedzib;

                   deriv.dea[ic][0] += dedxic;
                   deriv.dea[ic][1] += dedyic;
                   deriv.dea[ic][2] += dedzic;

                   deriv.dea[id][0] += dedxid;
                   deriv.dea[id][1] += dedyid;
                   deriv.dea[id][2] += dedzid;

                  virial.virx += xad*dedxia + xbd*dedxib + xcd*dedxic;
                  virial.viry += yad*dedyia + ybd*dedyib + ycd*dedyic;
                  virial.virz += zad*dedzia + zbd*dedzib + zcd*dedzic;
                }
            }
        }         
    }
 // fixed angle code 
    for (i=0; i < fixangle.nfxangle; i++)
    {
        ia = fixangle.ifxangle[i][0];
        ib = fixangle.ifxangle[i][1];
        ic = fixangle.ifxangle[i][2];
        id = 0;
        if (atom[ia].use || atom[ib].use || atom[ic].use)
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
                xab = xia - xib;
                yab = yia - yib;
                zab = zia - zib;
                rab2 = xab*xab + yab*yab + zab*zab;
                xcb = xic - xib;
                ycb = yic - yib;
                zcb = zic - zib;
                rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
                xp = ycb*zab - zcb*yab;
                yp = zcb*xab - xcb*zab;
                zp = xcb*yab - ycb*xab;
                rp = sqrt(xp*xp + yp*yp + zp*zp);
                if (rp != 0.0 )
                {
                    dot = xab*xcb + yab*ycb + zab*zcb;
                    cosine = dot / sqrt(rab2*rcb2);
                    if (cosine > 1.0)
                       cosine = 1.0;
                    if (cosine < -1.0)
                       cosine = -1.0;
                    angle = radian*acos(cosine);
                    if (angle < 0.75*angles.anat[i])
                      if ( atom[ib].mmx_type != 80) 
                         gmmx.hybrida = TRUE;

                    dt = angle - fixangle.anat[i];
                    dt2 = dt*dt;
                    dt3 = dt2*dt;
                    dt4 = dt2*dt2;

                    e = units.angunit*fixangle.acon[i]*dt2;
                    deddt = units.angunit*fixangle.acon[i]* dt*2.0;
                    deddt = deddt * radian;

                   terma = sqrt(1.00-( (dot*dot)/(rab2*rcb2)));
                   termc = sqrt(rab2*rcb2);
                   dedxia = -deddt*( xcb/termc - (xab*cosine)/rab2)/terma;
                   dedyia = -deddt*( ycb/termc - (yab*cosine)/rab2)/terma;
                   dedzia = -deddt*( zcb/termc - (zab*cosine)/rab2)/terma;

                   dedxic = -deddt*( xab/termc - (xcb*cosine)/rcb2)/terma;
                   dedyic = -deddt*( yab/termc - (ycb*cosine)/rcb2)/terma;
                   dedzic = -deddt*( zab/termc - (zcb*cosine)/rcb2)/terma;

                   dedxib = -deddt*( (-xab-xcb)/termc + (dot*(xcb*rab2 + xab*rcb2))/(rab2*rcb2*termc))/terma;
                   dedyib = -deddt*( (-yab-ycb)/termc + (dot*(ycb*rab2 + yab*rcb2))/(rab2*rcb2*termc))/terma;
                   dedzib = -deddt*( (-zab-zcb)/termc + (dot*(zcb*rab2 + zab*rcb2))/(rab2*rcb2*termc))/terma;
                  
                              
                   deriv.dea[ia][0] += dedxia;
                   deriv.dea[ia][1] += dedyia;
                   deriv.dea[ia][2] += dedzia;

                   deriv.dea[ib][0] += dedxib;
                   deriv.dea[ib][1] += dedyib;
                   deriv.dea[ib][2] += dedzib;

                   deriv.dea[ic][0] += dedxic;
                   deriv.dea[ic][1] += dedyic;
                   deriv.dea[ic][2] += dedzic;

                   energies.ebend += e;
                }  
        }         
    }
    //  end of fixed angles
}


void eangle2(int i)
{
    int j,k,ia,ib,ic,id;
    double old,eps;
    double **d0;   //  d0[MAXATOM][3];

    d0 = dmatrix(0,natom,0,3);

    eangle2a(i);

    eps = 1.0e-7;
    for (k=0; k < angles.nang; k++)
    {
        if ( (angles.angin[k]  == TRUE) && (field.type != MMFF94) )
        {
             ia = angles.i13[k][0];
             ib = angles.i13[k][1];
             ic = angles.i13[k][2];
             id = angles.i13[k][3];
             if (i == ia || i == ib || i == ic || i == id)
             {
                 eangle2b(k);
                 for (j=0; j < 3; j++)
                 {
                     d0[ia][j] = deriv.dea[ia][j];
                     d0[ib][j] = deriv.dea[ib][j];
                     d0[ic][j] = deriv.dea[ic][j];
                     d0[id][j] = deriv.dea[id][j];
                 }
          // numerical x derivative
                  old = atom[i].x;
                  atom[i].x += eps;
                  eangle2b(k);
                  atom[i].x = old;
                  for (j=0; j < 3; j++)
                  {
                      hess.hessx[ia][j] += (deriv.dea[ia][j]-d0[ia][j])/eps;
                      hess.hessx[ib][j] += (deriv.dea[ib][j]-d0[ib][j])/eps;
                      hess.hessx[ic][j] += (deriv.dea[ic][j]-d0[ic][j])/eps;
                      hess.hessx[id][j] += (deriv.dea[id][j]-d0[id][j])/eps;
                  }
          // numerical y derivative
                  old = atom[i].y;
                  atom[i].y += eps;
                  eangle2b(k);
                  atom[i].y = old;
                  for (j=0; j < 3; j++)
                  {
                      hess.hessy[ia][j] += (deriv.dea[ia][j]-d0[ia][j])/eps;
                      hess.hessy[ib][j] += (deriv.dea[ib][j]-d0[ib][j])/eps;
                      hess.hessy[ic][j] += (deriv.dea[ic][j]-d0[ic][j])/eps;
                      hess.hessy[id][j] += (deriv.dea[id][j]-d0[id][j])/eps;
                  }
          // numerical z derivative
                  old = atom[i].z;
                  atom[i].z += eps;
                  eangle2b(k);
                  atom[i].z = old;
                  for (j=0; j < 3; j++)
                  {
                      hess.hessz[ia][j] += (deriv.dea[ia][j]-d0[ia][j])/eps;
                      hess.hessz[ib][j] += (deriv.dea[ib][j]-d0[ib][j])/eps;
                      hess.hessz[ic][j] += (deriv.dea[ic][j]-d0[ic][j])/eps;
                      hess.hessz[id][j] += (deriv.dea[id][j]-d0[id][j])/eps;
                  }
             }
        }
    }
    // fixed angle code 
    if (fixangle.nfxangle > 0)
       eangle2a_fixed(i);
       
    free_dmatrix(d0 ,0,natom,0,3);
}
                     
void eangle2a(int iatom)
{
    int i,ia,ib,ic;
    double angle,dot,cosine,factor,sine;
    double dt,dt2,dt3,dt4;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xab,yab,zab,rab;
    double xcb,ycb,zcb,rcb;
    double xp,yp,zp,rp,rp2;
    double xrab,yrab,zrab,rab2;
    double xrcb,yrcb,zrcb,rcb2;
    double xabp,yabp,zabp;
    double xcbp,ycbp,zcbp;
    double deddt,d2eddt2,terma,termc;
    double ddtdxia,ddtdyia,ddtdzia;
    double ddtdxib,ddtdyib,ddtdzib;
    double ddtdxic,ddtdyic,ddtdzic;
    double dxiaxia,dxiayia,dxiazia;
    double dxibxib,dxibyib,dxibzib;
    double dxicxic,dxicyic,dxiczic;
    double dyiayia,dyiazia,dziazia;
    double dyibyib,dyibzib,dzibzib;
    double dyicyic,dyiczic,dziczic;
    double dxibxia,dxibyia,dxibzia;
    double dyibxia,dyibyia,dyibzia;
    double dzibxia,dzibyia,dzibzia;
    double dxibxic,dxibyic,dxibzic;
    double dyibxic,dyibyic,dyibzic;
    double dzibxic,dzibyic,dzibzic;
    double dxiaxic,dxiayic,dxiazic;
    double dyiaxic,dyiayic,dyiazic;
    double dziaxic,dziayic,dziazic;

    for (i=0; i < angles.nang; i++)
    {
        ia = angles.i13[i][0];
        ib = angles.i13[i][1];
        ic = angles.i13[i][2];
        if (iatom == ia || iatom == ib || iatom == ic)
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

            if ( (angles.angin[i] == FALSE) || field.type == MMFF94 )
            {
               xab = xia - xib;
               yab = yia - yib;
               zab = zia - zib;
               rab = sqrt(xab*xab + yab*yab + zab*zab);
               xcb = xic - xib;
               ycb = yic - yib;
               zcb = zic - zib;
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               xp = ycb*zab - zcb*yab;
               yp = zcb*xab - xcb*zab;
               zp = xcb*yab - ycb*xab;
               rp = sqrt(xp*xp + yp*yp + zp*zp);
               if (rp != 0.0)
               {
                  dot = xab*xcb + yab*ycb + zab*zcb;
                  cosine = dot / (rab*rcb);
                  if (cosine < -1.0)
                     cosine = -1.0;
                  if (cosine > 1.0)
                     cosine = 1.0;
                  angle = radian * acos(cosine);

                  if (angles.angtype[i] == HARMONIC)
                  {
                      dt = angle - angles.anat[i];
                      dt2 = dt * dt;
                      dt3 = dt2 * dt;
                      dt4 = dt3 * dt;
                      deddt = units.angunit * angles.acon[i] * dt 
                           * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                               + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);
                      d2eddt2 = units.angunit * angles.acon[i]
                            * (2.0 + 6.0*units.cang*dt + 12.0*units.qang*dt2
                                + 20.0*units.pang*dt3 + 30.0*units.sang*dt4);
                        terma = -radian/ (rab*rab*rp);
                        termc = radian / (rcb*rcb*rp);
                  } else
                  {
                        factor = angle/radian;
                        sine = sin(factor);
                        deddt = -units.angunit*angles.fcon[i]*
                          (angles.c1[i]*sine + 2*angles.c2[i]*sin(2*factor)
                          + 3*angles.c3[i]*sin(3*factor) + 4*angles.c4[i]*sin(4*factor)
                          + 5*angles.c5[i]*sin(5*factor) + 6*angles.c6[i]*sin(6*factor));
                        d2eddt2 = -units.angunit*angles.fcon[i]*
                          (angles.c1[i]*cosine + 4*angles.c2[i]*cos(2*factor)
                          + 9*angles.c3[i]*cos(3*factor) + 16*angles.c4[i]*cos(4*factor)
                          + 25*angles.c5[i]*cos(5*factor) + 36*angles.c6[i]*cos(6*factor));
                        terma = -1.0 / (rab*rab*rp);
                        termc = 1.0 / (rcb*rcb*rp);
                  }
                  ddtdxia = terma * (yab*zp-zab*yp);
                  ddtdyia = terma * (zab*xp-xab*zp);
                  ddtdzia = terma * (xab*yp-yab*xp);
                  ddtdxic = termc * (ycb*zp-zcb*yp);
                  ddtdyic = termc * (zcb*xp-xcb*zp);
                  ddtdzic = termc * (xcb*yp-ycb*xp);
                  ddtdxib = -ddtdxia - ddtdxic;
                  ddtdyib = -ddtdyia - ddtdyic;
                  ddtdzib = -ddtdzia - ddtdzic;
                  rab2 = 2.0 / (rab*rab);
                  xrab = xab * rab2;
                  yrab = yab * rab2;
                  zrab = zab * rab2;
                  rcb2 = 2.0 / (rcb*rcb);
                  xrcb = xcb * rcb2;
                  yrcb = ycb * rcb2;
                  zrcb = zcb * rcb2;
                  rp2 = 1.0 / (rp*rp);
                  xabp = (yab*zp-zab*yp) * rp2;
                  yabp = (zab*xp-xab*zp) * rp2;
                  zabp = (xab*yp-yab*xp) * rp2;
                  xcbp = (ycb*zp-zcb*yp) * rp2;
                  ycbp = (zcb*xp-xcb*zp) * rp2;
                  zcbp = (xcb*yp-ycb*xp) * rp2;

                  dxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab);
                  dxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab);
                  dxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab);
                  dyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab);
                  dyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab);
                  dziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab);
                  dxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb);
                  dxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb);
                  dxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb);
                  dyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb);
                  dyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb);
                  dziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb);
                  dxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp;
                  dxiayic = -terma*xab*yab - ddtdxia*yabp;
                  dxiazic = -terma*xab*zab - ddtdxia*zabp;
                  dyiaxic = -terma*xab*yab - ddtdyia*xabp;
                  dyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp;
                  dyiazic = -terma*yab*zab - ddtdyia*zabp;
                  dziaxic = -terma*xab*zab - ddtdzia*xabp;
                  dziayic = -terma*yab*zab - ddtdzia*yabp;
                  dziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp;

                  dxibxia = -dxiaxia - dxiaxic;
                  dxibyia = -dxiayia - dyiaxic;
                  dxibzia = -dxiazia - dziaxic;
                  dyibxia = -dxiayia - dxiayic;
                  dyibyia = -dyiayia - dyiayic;
                  dyibzia = -dyiazia - dziayic;
                  dzibxia = -dxiazia - dxiazic;
                  dzibyia = -dyiazia - dyiazic;
                  dzibzia = -dziazia - dziazic;
                  dxibxic = -dxicxic - dxiaxic;
                  dxibyic = -dxicyic - dxiayic;
                  dxibzic = -dxiczic - dxiazic;
                  dyibxic = -dxicyic - dyiaxic;
                  dyibyic = -dyicyic - dyiayic;
                  dyibzic = -dyiczic - dyiazic;
                  dzibxic = -dxiczic - dziaxic;
                  dzibyic = -dyiczic - dziayic;
                  dzibzic = -dziczic - dziazic;
                  dxibxib = -dxibxia - dxibxic;
                  dxibyib = -dxibyia - dxibyic;
                  dxibzib = -dxibzia - dxibzic;
                  dyibyib = -dyibyia - dyibyic;
                  dyibzib = -dyibzia - dyibzic;
                  dzibzib = -dzibzia - dzibzic;

                  if (ia == iatom)
                  {
                     hess.hessx[ia][0] += deddt*dxiaxia + d2eddt2*ddtdxia*ddtdxia;
                     hess.hessx[ia][1] += deddt*dxiayia + d2eddt2*ddtdxia*ddtdyia;
                     hess.hessx[ia][2] += deddt*dxiazia + d2eddt2*ddtdxia*ddtdzia;
                     hess.hessy[ia][0] += deddt*dxiayia + d2eddt2*ddtdyia*ddtdxia;
                     hess.hessy[ia][1] += deddt*dyiayia + d2eddt2*ddtdyia*ddtdyia;
                     hess.hessy[ia][2] += deddt*dyiazia + d2eddt2*ddtdyia*ddtdzia;
                     hess.hessz[ia][0] += deddt*dxiazia + d2eddt2*ddtdzia*ddtdxia;
                     hess.hessz[ia][1] += deddt*dyiazia + d2eddt2*ddtdzia*ddtdyia;
                     hess.hessz[ia][2] += deddt*dziazia + d2eddt2*ddtdzia*ddtdzia;
                     hess.hessx[ib][0] += deddt*dxibxia + d2eddt2*ddtdxia*ddtdxib;
                     hess.hessx[ib][1] += deddt*dyibxia + d2eddt2*ddtdxia*ddtdyib;
                     hess.hessx[ib][2] += deddt*dzibxia + d2eddt2*ddtdxia*ddtdzib;
                     hess.hessy[ib][0] += deddt*dxibyia + d2eddt2*ddtdyia*ddtdxib;
                     hess.hessy[ib][1] += deddt*dyibyia + d2eddt2*ddtdyia*ddtdyib;
                     hess.hessy[ib][2] += deddt*dzibyia + d2eddt2*ddtdyia*ddtdzib;
                     hess.hessz[ib][0] += deddt*dxibzia + d2eddt2*ddtdzia*ddtdxib;
                     hess.hessz[ib][1] += deddt*dyibzia + d2eddt2*ddtdzia*ddtdyib;
                     hess.hessz[ib][2] += deddt*dzibzia + d2eddt2*ddtdzia*ddtdzib;
                     hess.hessx[ic][0] += deddt*dxiaxic + d2eddt2*ddtdxia*ddtdxic;
                     hess.hessx[ic][1] += deddt*dxiayic + d2eddt2*ddtdxia*ddtdyic;
                     hess.hessx[ic][2] += deddt*dxiazic + d2eddt2*ddtdxia*ddtdzic;
                     hess.hessy[ic][0] += deddt*dyiaxic + d2eddt2*ddtdyia*ddtdxic;
                     hess.hessy[ic][1] += deddt*dyiayic + d2eddt2*ddtdyia*ddtdyic;
                     hess.hessy[ic][2] += deddt*dyiazic + d2eddt2*ddtdyia*ddtdzic;
                     hess.hessz[ic][0] += deddt*dziaxic + d2eddt2*ddtdzia*ddtdxic;
                     hess.hessz[ic][1] += deddt*dziayic + d2eddt2*ddtdzia*ddtdyic;
                     hess.hessz[ic][2] += deddt*dziazic + d2eddt2*ddtdzia*ddtdzic;
                  } else if (ib == iatom)
                  {
                     hess.hessx[ib][0] += deddt*dxibxib + d2eddt2*ddtdxib*ddtdxib;
                     hess.hessx[ib][1] += deddt*dxibyib + d2eddt2*ddtdxib*ddtdyib;
                     hess.hessx[ib][2] += deddt*dxibzib + d2eddt2*ddtdxib*ddtdzib;
                     hess.hessy[ib][0] += deddt*dxibyib + d2eddt2*ddtdyib*ddtdxib;
                     hess.hessy[ib][1] += deddt*dyibyib + d2eddt2*ddtdyib*ddtdyib;
                     hess.hessy[ib][2] += deddt*dyibzib + d2eddt2*ddtdyib*ddtdzib;
                     hess.hessz[ib][0] += deddt*dxibzib + d2eddt2*ddtdzib*ddtdxib;
                     hess.hessz[ib][1] += deddt*dyibzib + d2eddt2*ddtdzib*ddtdyib;
                     hess.hessz[ib][2] += deddt*dzibzib + d2eddt2*ddtdzib*ddtdzib;
                     hess.hessx[ia][0] += deddt*dxibxia + d2eddt2*ddtdxib*ddtdxia;
                     hess.hessx[ia][1] += deddt*dxibyia + d2eddt2*ddtdxib*ddtdyia;
                     hess.hessx[ia][2] += deddt*dxibzia + d2eddt2*ddtdxib*ddtdzia;
                     hess.hessy[ia][0] += deddt*dyibxia + d2eddt2*ddtdyib*ddtdxia;
                     hess.hessy[ia][1] += deddt*dyibyia + d2eddt2*ddtdyib*ddtdyia;
                     hess.hessy[ia][2] += deddt*dyibzia + d2eddt2*ddtdyib*ddtdzia;
                     hess.hessz[ia][0] += deddt*dzibxia + d2eddt2*ddtdzib*ddtdxia;
                     hess.hessz[ia][1] += deddt*dzibyia + d2eddt2*ddtdzib*ddtdyia;
                     hess.hessz[ia][2] += deddt*dzibzia + d2eddt2*ddtdzib*ddtdzia;
                     hess.hessx[ic][0] += deddt*dxibxic + d2eddt2*ddtdxib*ddtdxic;
                     hess.hessx[ic][1] += deddt*dxibyic + d2eddt2*ddtdxib*ddtdyic;
                     hess.hessx[ic][2] += deddt*dxibzic + d2eddt2*ddtdxib*ddtdzic;
                     hess.hessy[ic][0] += deddt*dyibxic + d2eddt2*ddtdyib*ddtdxic;
                     hess.hessy[ic][1] += deddt*dyibyic + d2eddt2*ddtdyib*ddtdyic;
                     hess.hessy[ic][2] += deddt*dyibzic + d2eddt2*ddtdyib*ddtdzic;
                     hess.hessz[ic][0] += deddt*dzibxic + d2eddt2*ddtdzib*ddtdxic;
                     hess.hessz[ic][1] += deddt*dzibyic + d2eddt2*ddtdzib*ddtdyic;
                     hess.hessz[ic][2] += deddt*dzibzic + d2eddt2*ddtdzib*ddtdzic;
                  }else if (ic == iatom)
                  {
                     hess.hessx[ic][0] += deddt*dxicxic + d2eddt2*ddtdxic*ddtdxic;
                     hess.hessx[ic][1] += deddt*dxicyic + d2eddt2*ddtdxic*ddtdyic;
                     hess.hessx[ic][2] += deddt*dxiczic + d2eddt2*ddtdxic*ddtdzic;
                     hess.hessy[ic][0] += deddt*dxicyic + d2eddt2*ddtdyic*ddtdxic;
                     hess.hessy[ic][1] += deddt*dyicyic + d2eddt2*ddtdyic*ddtdyic;
                     hess.hessy[ic][2] += deddt*dyiczic + d2eddt2*ddtdyic*ddtdzic;
                     hess.hessz[ic][0] += deddt*dxiczic + d2eddt2*ddtdzic*ddtdxic;
                     hess.hessz[ic][1] += deddt*dyiczic + d2eddt2*ddtdzic*ddtdyic;
                     hess.hessz[ic][2] += deddt*dziczic + d2eddt2*ddtdzic*ddtdzic;
                     hess.hessx[ib][0] += deddt*dxibxic + d2eddt2*ddtdxic*ddtdxib;
                     hess.hessx[ib][1] += deddt*dyibxic + d2eddt2*ddtdxic*ddtdyib;
                     hess.hessx[ib][2] += deddt*dzibxic + d2eddt2*ddtdxic*ddtdzib;
                     hess.hessy[ib][0] += deddt*dxibyic + d2eddt2*ddtdyic*ddtdxib;
                     hess.hessy[ib][1] += deddt*dyibyic + d2eddt2*ddtdyic*ddtdyib;
                     hess.hessy[ib][2] += deddt*dzibyic + d2eddt2*ddtdyic*ddtdzib;
                     hess.hessz[ib][0] += deddt*dxibzic + d2eddt2*ddtdzic*ddtdxib;
                     hess.hessz[ib][1] += deddt*dyibzic + d2eddt2*ddtdzic*ddtdyib;
                     hess.hessz[ib][2] += deddt*dzibzic + d2eddt2*ddtdzic*ddtdzib;
                     hess.hessx[ia][0] += deddt*dxiaxic + d2eddt2*ddtdxic*ddtdxia;
                     hess.hessx[ia][1] += deddt*dyiaxic + d2eddt2*ddtdxic*ddtdyia;
                     hess.hessx[ia][2] += deddt*dziaxic + d2eddt2*ddtdxic*ddtdzia;
                     hess.hessy[ia][0] += deddt*dxiayic + d2eddt2*ddtdyic*ddtdxia;
                     hess.hessy[ia][1] += deddt*dyiayic + d2eddt2*ddtdyic*ddtdyia;
                     hess.hessy[ia][2] += deddt*dziayic + d2eddt2*ddtdyic*ddtdzia;
                     hess.hessz[ia][0] += deddt*dxiazic + d2eddt2*ddtdzic*ddtdxia;
                     hess.hessz[ia][1] += deddt*dyiazic + d2eddt2*ddtdzic*ddtdyia;
                     hess.hessz[ia][2] += deddt*dziazic + d2eddt2*ddtdzic*ddtdzia;
                  }
               }
            }
        }
    }
}

void eangle2b(int i)
{
    int ia,ib,ic,id;
    double angle,dot,cosine;
    double dt,dt2,dt3,dt4;
    double deddt,terma,termc,term;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xid,yid,zid;
    double xad,yad,zad;
    double xbd,ybd,zbd;
    double xcd,ycd,zcd;
    double xip,yip,zip;
    double xap,yap,zap,rap2;
    double xcp,ycp,zcp,rcp2;
    double xt,yt,zt,rt2,ptrt2;
    double xm,ym,zm,rm,delta,delta2;
    double dedxia,dedyia,dedzia;
    double dedxib,dedyib,dedzib;
    double dedxic,dedyic,dedzic;
    double dedxid,dedyid,dedzid;
    double dedxip,dedyip,dedzip;
    double dpdxia,dpdyia,dpdzia;
    double dpdxic,dpdyic,dpdzic;

    ia = angles.i13[i][0];
    ib = angles.i13[i][1];
    ic = angles.i13[i][2];
    id = angles.i13[i][3];
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

    deriv.dea[ia][0] = 0.0;
    deriv.dea[ia][1] = 0.0;
    deriv.dea[ia][2] = 0.0;
    deriv.dea[ib][0] = 0.0;
    deriv.dea[ib][1] = 0.0;
    deriv.dea[ib][2] = 0.0;
    deriv.dea[ic][0] = 0.0;
    deriv.dea[ic][1] = 0.0;
    deriv.dea[ic][2] = 0.0;
    deriv.dea[id][0] = 0.0;
    deriv.dea[id][1] = 0.0;
    deriv.dea[id][2] = 0.0;

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
      xip = xib + xt*delta;
      yip = yib + yt*delta;
      zip = zib + zt*delta;
      xap = xia - xip;
      yap = yia - yip;
      zap = zia - zip;
      rap2 = xap*xap + yap*yap + zap*zap;
      xcp = xic - xip;
      ycp = yic - yip;
      zcp = zic - zip;
      rcp2 = xcp*xcp + ycp*ycp + zcp*zcp;
      xm = ycp*zap - zcp*yap;
      ym = zcp*xap - xcp*zap;
      zm = xcp*yap - ycp*xap;
      rm = sqrt(xm*xm + ym*ym + zm*zm);
      if (rm != 0.0)
      {
         dot = xap*xcp + yap*ycp + zap*zcp;
         cosine = dot / sqrt(rap2*rcp2);
         if (cosine < -1.0)
           cosine = -1.0;
         if (cosine > 1.0)
           cosine = 1.0;
         angle = radian * acos(cosine);

         dt = angle - angles.anat[i];
         dt2 = dt * dt;
         dt3 = dt2 * dt;
         dt4 = dt2 * dt2;
         deddt = units.angunit * angles.acon[i] * dt
                 * (2.0 + 3.0*units.cang*dt + 4.0*units.qang*dt2
                    + 5.0*units.pang*dt3 + 6.0*units.sang*dt4);

         deddt = deddt * radian;
         terma = -deddt / (rap2*rm);
         termc = deddt / (rcp2*rm);
         dedxia = terma * (yap*zm-zap*ym);
         dedyia = terma * (zap*xm-xap*zm);
         dedzia = terma * (xap*ym-yap*xm);
         dedxic = termc * (ycp*zm-zcp*ym);
         dedyic = termc * (zcp*xm-xcp*zm);
         dedzic = termc * (xcp*ym-ycp*xm);
         dedxip = -dedxia - dedxic;
         dedyip = -dedyia - dedyic;
         dedzip = -dedzia - dedzic;

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
         dedxia = dedxia + dpdxia;
         dedyia = dedyia + dpdyia;
         dedzia = dedzia + dpdzia;
         dedxib = dedxip;
         dedyib = dedyip;
         dedzib = dedzip;
         dedxic = dedxic + dpdxic;
         dedyic = dedyic + dpdyic;
         dedzic = dedzic + dpdzic;
         dedxid = -dedxia - dedxib - dedxic;
         dedyid = -dedyia - dedyib - dedyic;
         dedzid = -dedzia - dedzib - dedzic;

         deriv.dea[ia][0] = dedxia;
         deriv.dea[ia][1] = dedyia;
         deriv.dea[ia][2] = dedzia;
         deriv.dea[ib][0] = dedxib;
         deriv.dea[ib][1] = dedyib;
         deriv.dea[ib][2] = dedzib;
         deriv.dea[ic][0] = dedxic;
         deriv.dea[ic][1] = dedyic;
         deriv.dea[ic][2] = dedzic;
         deriv.dea[id][0] = dedxid;
         deriv.dea[id][1] = dedyid;
         deriv.dea[id][2] = dedzid;
      }
}

void eangle2a_fixed(int iatom)
{
    int i,ia,ib,ic;
    double angle,dot,cosine;
    double dt,dt2,dt3,dt4;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xab,yab,zab,rab;
    double xcb,ycb,zcb,rcb;
    double xp,yp,zp,rp,rp2;
    double xrab,yrab,zrab,rab2;
    double xrcb,yrcb,zrcb,rcb2;
    double xabp,yabp,zabp;
    double xcbp,ycbp,zcbp;
    double deddt,d2eddt2,terma,termc;
    double ddtdxia,ddtdyia,ddtdzia;
    double ddtdxib,ddtdyib,ddtdzib;
    double ddtdxic,ddtdyic,ddtdzic;
    double dxiaxia,dxiayia,dxiazia;
    double dxibxib,dxibyib,dxibzib;
    double dxicxic,dxicyic,dxiczic;
    double dyiayia,dyiazia,dziazia;
    double dyibyib,dyibzib,dzibzib;
    double dyicyic,dyiczic,dziczic;
    double dxibxia,dxibyia,dxibzia;
    double dyibxia,dyibyia,dyibzia;
    double dzibxia,dzibyia,dzibzia;
    double dxibxic,dxibyic,dxibzic;
    double dyibxic,dyibyic,dyibzic;
    double dzibxic,dzibyic,dzibzic;
    double dxiaxic,dxiayic,dxiazic;
    double dyiaxic,dyiayic,dyiazic;
    double dziaxic,dziayic,dziazic;

    for (i=0; i < fixangle.nfxangle; i++)
    {
        ia = fixangle.ifxangle[i][0];
        ib = fixangle.ifxangle[i][1];
        ic = fixangle.ifxangle[i][2];
        if (iatom == ia || iatom == ib || iatom == ic)
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

               xab = xia - xib;
               yab = yia - yib;
               zab = zia - zib;
               rab = sqrt(xab*xab + yab*yab + zab*zab);
               xcb = xic - xib;
               ycb = yic - yib;
               zcb = zic - zib;
               rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
               xp = ycb*zab - zcb*yab;
               yp = zcb*xab - xcb*zab;
               zp = xcb*yab - ycb*xab;
               rp = sqrt(xp*xp + yp*yp + zp*zp);
               if (rp != 0.0)
               {
                  dot = xab*xcb + yab*ycb + zab*zcb;
                  cosine = dot / (rab*rcb);
                  if (cosine < -1.0)
                     cosine = -1.0;
                  if (cosine > 1.0)
                     cosine = 1.0;
                  angle = radian * acos(cosine);
                
                  dt = angle - fixangle.anat[i];
                  dt2 = dt * dt;
                  dt3 = dt2 * dt;
                  dt4 = dt3 * dt;
                  deddt = units.angunit * fixangle.acon[i] * dt*2.0;
                  d2eddt2 = units.angunit * fixangle.acon[i]*2.0;
        
                  terma = -radian / (rab*rab*rp);
                  termc = radian / (rcb*rcb*rp);
                  ddtdxia = terma * (yab*zp-zab*yp);
                  ddtdyia = terma * (zab*xp-xab*zp);
                  ddtdzia = terma * (xab*yp-yab*xp);
                  ddtdxic = termc * (ycb*zp-zcb*yp);
                  ddtdyic = termc * (zcb*xp-xcb*zp);
                  ddtdzic = termc * (xcb*yp-ycb*xp);
                  ddtdxib = -ddtdxia - ddtdxic;
                  ddtdyib = -ddtdyia - ddtdyic;
                  ddtdzib = -ddtdzia - ddtdzic;
                  rab2 = 2.0 / (rab*rab);
                  xrab = xab * rab2;
                  yrab = yab * rab2;
                  zrab = zab * rab2;
                  rcb2 = 2.0 / (rcb*rcb);
                  xrcb = xcb * rcb2;
                  yrcb = ycb * rcb2;
                  zrcb = zcb * rcb2;
                  rp2 = 1.0 / (rp*rp);
                  xabp = (yab*zp-zab*yp) * rp2;
                  yabp = (zab*xp-xab*zp) * rp2;
                  zabp = (xab*yp-yab*xp) * rp2;
                  xcbp = (ycb*zp-zcb*yp) * rp2;
                  ycbp = (zcb*xp-xcb*zp) * rp2;
                  zcbp = (xcb*yp-ycb*xp) * rp2;

                  dxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab);
                  dxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab);
                  dxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab);
                  dyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab);
                  dyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab);
                  dziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab);
                  dxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb);
                  dxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb);
                  dxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb);
                  dyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb);
                  dyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb);
                  dziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb);
                  dxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp;
                  dxiayic = -terma*xab*yab - ddtdxia*yabp;
                  dxiazic = -terma*xab*zab - ddtdxia*zabp;
                  dyiaxic = -terma*xab*yab - ddtdyia*xabp;
                  dyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp;
                  dyiazic = -terma*yab*zab - ddtdyia*zabp;
                  dziaxic = -terma*xab*zab - ddtdzia*xabp;
                  dziayic = -terma*yab*zab - ddtdzia*yabp;
                  dziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp;

                  dxibxia = -dxiaxia - dxiaxic;
                  dxibyia = -dxiayia - dyiaxic;
                  dxibzia = -dxiazia - dziaxic;
                  dyibxia = -dxiayia - dxiayic;
                  dyibyia = -dyiayia - dyiayic;
                  dyibzia = -dyiazia - dziayic;
                  dzibxia = -dxiazia - dxiazic;
                  dzibyia = -dyiazia - dyiazic;
                  dzibzia = -dziazia - dziazic;
                  dxibxic = -dxicxic - dxiaxic;
                  dxibyic = -dxicyic - dxiayic;
                  dxibzic = -dxiczic - dxiazic;
                  dyibxic = -dxicyic - dyiaxic;
                  dyibyic = -dyicyic - dyiayic;
                  dyibzic = -dyiczic - dyiazic;
                  dzibxic = -dxiczic - dziaxic;
                  dzibyic = -dyiczic - dziayic;
                  dzibzic = -dziczic - dziazic;
                  dxibxib = -dxibxia - dxibxic;
                  dxibyib = -dxibyia - dxibyic;
                  dxibzib = -dxibzia - dxibzic;
                  dyibyib = -dyibyia - dyibyic;
                  dyibzib = -dyibzia - dyibzic;
                  dzibzib = -dzibzia - dzibzic;

                  if (ia == iatom)
                  {
                     hess.hessx[ia][0] += deddt*dxiaxia + d2eddt2*ddtdxia*ddtdxia;
                     hess.hessx[ia][1] += deddt*dxiayia + d2eddt2*ddtdxia*ddtdyia;
                     hess.hessx[ia][2] += deddt*dxiazia + d2eddt2*ddtdxia*ddtdzia;
                     hess.hessy[ia][0] += deddt*dxiayia + d2eddt2*ddtdyia*ddtdxia;
                     hess.hessy[ia][1] += deddt*dyiayia + d2eddt2*ddtdyia*ddtdyia;
                     hess.hessy[ia][2] += deddt*dyiazia + d2eddt2*ddtdyia*ddtdzia;
                     hess.hessz[ia][0] += deddt*dxiazia + d2eddt2*ddtdzia*ddtdxia;
                     hess.hessz[ia][1] += deddt*dyiazia + d2eddt2*ddtdzia*ddtdyia;
                     hess.hessz[ia][2] += deddt*dziazia + d2eddt2*ddtdzia*ddtdzia;
                     hess.hessx[ib][0] += deddt*dxibxia + d2eddt2*ddtdxia*ddtdxib;
                     hess.hessx[ib][1] += deddt*dyibxia + d2eddt2*ddtdxia*ddtdyib;
                     hess.hessx[ib][2] += deddt*dzibxia + d2eddt2*ddtdxia*ddtdzib;
                     hess.hessy[ib][0] += deddt*dxibyia + d2eddt2*ddtdyia*ddtdxib;
                     hess.hessy[ib][1] += deddt*dyibyia + d2eddt2*ddtdyia*ddtdyib;
                     hess.hessy[ib][2] += deddt*dzibyia + d2eddt2*ddtdyia*ddtdzib;
                     hess.hessz[ib][0] += deddt*dxibzia + d2eddt2*ddtdzia*ddtdxib;
                     hess.hessz[ib][1] += deddt*dyibzia + d2eddt2*ddtdzia*ddtdyib;
                     hess.hessz[ib][2] += deddt*dzibzia + d2eddt2*ddtdzia*ddtdzib;
                     hess.hessx[ic][0] += deddt*dxiaxic + d2eddt2*ddtdxia*ddtdxic;
                     hess.hessx[ic][1] += deddt*dxiayic + d2eddt2*ddtdxia*ddtdyic;
                     hess.hessx[ic][2] += deddt*dxiazic + d2eddt2*ddtdxia*ddtdzic;
                     hess.hessy[ic][0] += deddt*dyiaxic + d2eddt2*ddtdyia*ddtdxic;
                     hess.hessy[ic][1] += deddt*dyiayic + d2eddt2*ddtdyia*ddtdyic;
                     hess.hessy[ic][2] += deddt*dyiazic + d2eddt2*ddtdyia*ddtdzic;
                     hess.hessz[ic][0] += deddt*dziaxic + d2eddt2*ddtdzia*ddtdxic;
                     hess.hessz[ic][1] += deddt*dziayic + d2eddt2*ddtdzia*ddtdyic;
                     hess.hessz[ic][2] += deddt*dziazic + d2eddt2*ddtdzia*ddtdzic;
                  } else if (ib == iatom)
                  {
                     hess.hessx[ib][0] += deddt*dxibxib + d2eddt2*ddtdxib*ddtdxib;
                     hess.hessx[ib][1] += deddt*dxibyib + d2eddt2*ddtdxib*ddtdyib;
                     hess.hessx[ib][2] += deddt*dxibzib + d2eddt2*ddtdxib*ddtdzib;
                     hess.hessy[ib][0] += deddt*dxibyib + d2eddt2*ddtdyib*ddtdxib;
                     hess.hessy[ib][1] += deddt*dyibyib + d2eddt2*ddtdyib*ddtdyib;
                     hess.hessy[ib][2] += deddt*dyibzib + d2eddt2*ddtdyib*ddtdzib;
                     hess.hessz[ib][0] += deddt*dxibzib + d2eddt2*ddtdzib*ddtdxib;
                     hess.hessz[ib][1] += deddt*dyibzib + d2eddt2*ddtdzib*ddtdyib;
                     hess.hessz[ib][2] += deddt*dzibzib + d2eddt2*ddtdzib*ddtdzib;
                     hess.hessx[ia][0] += deddt*dxibxia + d2eddt2*ddtdxib*ddtdxia;
                     hess.hessx[ia][1] += deddt*dxibyia + d2eddt2*ddtdxib*ddtdyia;
                     hess.hessx[ia][2] += deddt*dxibzia + d2eddt2*ddtdxib*ddtdzia;
                     hess.hessy[ia][0] += deddt*dyibxia + d2eddt2*ddtdyib*ddtdxia;
                     hess.hessy[ia][1] += deddt*dyibyia + d2eddt2*ddtdyib*ddtdyia;
                     hess.hessy[ia][2] += deddt*dyibzia + d2eddt2*ddtdyib*ddtdzia;
                     hess.hessz[ia][0] += deddt*dzibxia + d2eddt2*ddtdzib*ddtdxia;
                     hess.hessz[ia][1] += deddt*dzibyia + d2eddt2*ddtdzib*ddtdyia;
                     hess.hessz[ia][2] += deddt*dzibzia + d2eddt2*ddtdzib*ddtdzia;
                     hess.hessx[ic][0] += deddt*dxibxic + d2eddt2*ddtdxib*ddtdxic;
                     hess.hessx[ic][1] += deddt*dxibyic + d2eddt2*ddtdxib*ddtdyic;
                     hess.hessx[ic][2] += deddt*dxibzic + d2eddt2*ddtdxib*ddtdzic;
                     hess.hessy[ic][0] += deddt*dyibxic + d2eddt2*ddtdyib*ddtdxic;
                     hess.hessy[ic][1] += deddt*dyibyic + d2eddt2*ddtdyib*ddtdyic;
                     hess.hessy[ic][2] += deddt*dyibzic + d2eddt2*ddtdyib*ddtdzic;
                     hess.hessz[ic][0] += deddt*dzibxic + d2eddt2*ddtdzib*ddtdxic;
                     hess.hessz[ic][1] += deddt*dzibyic + d2eddt2*ddtdzib*ddtdyic;
                     hess.hessz[ic][2] += deddt*dzibzic + d2eddt2*ddtdzib*ddtdzic;
                  }else if (ic == iatom)
                  {
                     hess.hessx[ic][0] += deddt*dxicxic + d2eddt2*ddtdxic*ddtdxic;
                     hess.hessx[ic][1] += deddt*dxicyic + d2eddt2*ddtdxic*ddtdyic;
                     hess.hessx[ic][2] += deddt*dxiczic + d2eddt2*ddtdxic*ddtdzic;
                     hess.hessy[ic][0] += deddt*dxicyic + d2eddt2*ddtdyic*ddtdxic;
                     hess.hessy[ic][1] += deddt*dyicyic + d2eddt2*ddtdyic*ddtdyic;
                     hess.hessy[ic][2] += deddt*dyiczic + d2eddt2*ddtdyic*ddtdzic;
                     hess.hessz[ic][0] += deddt*dxiczic + d2eddt2*ddtdzic*ddtdxic;
                     hess.hessz[ic][1] += deddt*dyiczic + d2eddt2*ddtdzic*ddtdyic;
                     hess.hessz[ic][2] += deddt*dziczic + d2eddt2*ddtdzic*ddtdzic;
                     hess.hessx[ib][0] += deddt*dxibxic + d2eddt2*ddtdxic*ddtdxib;
                     hess.hessx[ib][1] += deddt*dyibxic + d2eddt2*ddtdxic*ddtdyib;
                     hess.hessx[ib][2] += deddt*dzibxic + d2eddt2*ddtdxic*ddtdzib;
                     hess.hessy[ib][0] += deddt*dxibyic + d2eddt2*ddtdyic*ddtdxib;
                     hess.hessy[ib][1] += deddt*dyibyic + d2eddt2*ddtdyic*ddtdyib;
                     hess.hessy[ib][2] += deddt*dzibyic + d2eddt2*ddtdyic*ddtdzib;
                     hess.hessz[ib][0] += deddt*dxibzic + d2eddt2*ddtdzic*ddtdxib;
                     hess.hessz[ib][1] += deddt*dyibzic + d2eddt2*ddtdzic*ddtdyib;
                     hess.hessz[ib][2] += deddt*dzibzic + d2eddt2*ddtdzic*ddtdzib;
                     hess.hessx[ia][0] += deddt*dxiaxic + d2eddt2*ddtdxic*ddtdxia;
                     hess.hessx[ia][1] += deddt*dyiaxic + d2eddt2*ddtdxic*ddtdyia;
                     hess.hessx[ia][2] += deddt*dziaxic + d2eddt2*ddtdxic*ddtdzia;
                     hess.hessy[ia][0] += deddt*dxiayic + d2eddt2*ddtdyic*ddtdxia;
                     hess.hessy[ia][1] += deddt*dyiayic + d2eddt2*ddtdyic*ddtdyia;
                     hess.hessy[ia][2] += deddt*dziayic + d2eddt2*ddtdyic*ddtdzia;
                     hess.hessz[ia][0] += deddt*dxiazic + d2eddt2*ddtdzic*ddtdxia;
                     hess.hessz[ia][1] += deddt*dyiazic + d2eddt2*ddtdzic*ddtdyia;
                     hess.hessz[ia][2] += deddt*dziazic + d2eddt2*ddtdzic*ddtdzia;
                  }
               }
        }
    }
}

