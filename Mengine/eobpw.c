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
#include "units.h"
#include "minim_values.h"
#include "opbend.h"
        
void eopbend_wilson()
{
      int i,ia,ib,ic,id;
      double e,angle,dot,cosine;
      double dt,dt2;
      double xia,yia,zia;
      double xib,yib,zib;
      double xic,yic,zic;
      double xid,yid,zid;
      double xab,yab,zab,rab2,rab;
      double xbc,ybc,zbc,rbc2,rbc;
      double xbd,ybd,zbd,rbd2,rbd;
      double sine;

      energies.eopb = 0.0;

      if (minim_values.iprint)
      {
          fprintf(pcmoutfile,"\nWilson Out of Plane Bending Terms\n");
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
                 
                 xab = xia - xib;
                 yab = yia - yib;
                 zab = zia - zib;
                 rab2 = (xab*xab+yab*yab+zab*zab);
                 rab = sqrt(rab2);
                 
                 xbc = xic - xib;
                 ybc = yic - yib;
                 zbc = zic - zib;
                 rbc2 = xbc*xbc+ybc*ybc+zbc*zbc;
                 rbc = sqrt(rbc2);
                 
                 if (rab2 != 0.0 && rbc2 != 0.0)
                 {
                    dot = xab*xbc + yab*ybc + zab*zbc;
                    cosine = dot /(rab*rbc);
                    if (cosine < -1.0)
                       cosine = -1.0;
                    if (cosine > 1.0)
                       cosine = 1.0;
                    sine = sqrt(1.00-cosine*cosine);

                    xbd = atom[id].x - atom[ib].x;
                    ybd = atom[id].y - atom[ib].y;
                    zbd = atom[id].z - atom[ib].z;
                    rbd2 = xbd*xbd+ybd*ybd+zbd*zbd;
                    rbd = sqrt(rbd2);
                    
                    cosine = ( xbd*(yab*zbc-zab*ybc) + ybd*(zab*xbc-xab*zbc) + zbd*(xab*ybc-yab*xbc))/((rab*rbc*rbd)*sine);

                    angle = radian * asin(cosine);
                    dt = angle;
                    dt2 = dt * dt;
                    e = units.angunit * angles.copb[i] * dt2;
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

void eopbend_wilson1()
{
      int i,ia,ib,ic,id;
      double e,angle,dot,cosine;
      double dt,dt2;
      double xia,yia,zia;
      double xib,yib,zib;
      double xic,yic,zic;
      double xid,yid,zid;
      double xab,yab,zab,rab2,rab;
      double xbc,ybc,zbc,rbc2,rbc;
      double xbd,ybd,zbd,rbd2,rbd;
      double sine;
      double deddt, terma, denom,termc;
      double dedxia, dedxib, dedxic, dedxid;
      double dedyia, dedyib, dedyic, dedyid;
      double dedzia, dedzib, dedzic, dedzid;

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
                 
                 xab = xia - xib;
                 yab = yia - yib;
                 zab = zia - zib;
                 rab2 = (xab*xab+yab*yab+zab*zab);
                 rab = sqrt(rab2);
                 
                 xbc = xic - xib;
                 ybc = yic - yib;
                 zbc = zic - zib;
                 rbc2 = xbc*xbc+ybc*ybc+zbc*zbc;
                 rbc = sqrt(rbc2);
                 
                 if (rab2 != 0.0 && rbc2 != 0.0)
                 {
                    dot = xab*xbc + yab*ybc + zab*zbc;
                    cosine = dot /(rab*rbc);
                    if (cosine < -1.0)
                       cosine = -1.0;
                    if (cosine > 1.0)
                       cosine = 1.0;
                    sine = sqrt(1.00-cosine*cosine);

                    xbd = atom[id].x - atom[ib].x;
                    ybd = atom[id].y - atom[ib].y;
                    zbd = atom[id].z - atom[ib].z;
                    rbd2 = xbd*xbd+ybd*ybd+zbd*zbd;
                    rbd = sqrt(rbd2);
                    
                    cosine = ( xbd*(yab*zbc-zab*ybc) + ybd*(zab*xbc-xab*zbc) + zbd*(xab*ybc-yab*xbc))/((rab*rbc*rbd)*sine);

                    angle = radian * asin(cosine);
                    dt = angle;
                    dt2 = dt * dt;
                    e = units.angunit * angles.copb[i] * dt2;
                    energies.eopb += e;

                    deddt = 2.0*units.angunit*angles.copb[i]*dt;
                    deddt *= radian;
                    
                    terma = xbd*(yab*zbc-zab*ybc) + ybd*(zab*xbc-xab*zbc) + zbd*(xab*ybc - yab*xbc);
                    termc = rab*rbc*rbd*sine;
                    denom = sqrt(1.00 - (terma*terma)/(rab2*rbd2*rbc2*sine*sine));
                    
                    dedxia = (-(terma*( (2.0*xab*dot*dot)/(rab2*rab2*rbc2) - (2.0*xbc*dot)/(rab2*rbc2))/(2.0*termc*sine*sine) )
                             - (xab*terma)/(rab2*termc) + (ybc*zbd-ybd*zbc)/termc ) / denom;
                    dedyia = (-(terma*( (2.0*yab*dot*dot)/(rab2*rab2*rbc2) - (2.0*ybc*dot)/(rab2*rbc2))/(2.0*termc*sine*sine) )
                             - (yab*terma)/(rab2*termc) + (zbc*xbd-zbd*xbc)/termc ) / denom;
                    dedzia = (-(terma*( (2.0*zab*dot*dot)/(rab2*rab2*rbc2) - (2.0*zbc*dot)/(rab2*rbc2))/(2.0*termc*sine*sine) )
                             - (zab*terma)/(rab2*termc) + (xbc*ybd-xbd*ybc)/termc ) / denom;
                             
                    dedxib = ( (xbd*terma)/(termc*rbd2) + (ybc*zab-yab*zbc+ybd*(zic-zia)+zbd*(yia-yic))/termc
                              +(xbc*terma)/(termc*rbc2) + (xab*terma)/(termc*rab2)
                              -( (-(2.0*xbc*dot*dot)/(rab2*rbc2*rbc2) - (2.0*(2.0*xib-xia-xic)*dot)/(rab2*rbc2)
                                 -(2.0*xab*dot*dot)/(rab2*rab2*rbc2))*terma)/(2.0*termc*sine*sine)      )/ denom;
                                 
                    dedyib = ( (ybd*terma)/(termc*rbd2) + (xab*zbc-zab*xbc+xbd*(zab-zbc)+zbd*(xic-xia))/termc
                              +(ybc*terma)/(termc*rbc2) + (yab*terma)/(termc*rab2)
                              -( (-(2.0*ybc*dot*dot)/(rab2*rbc2*rbc2) - (2.0*(2.0*yib-yia-yic)*dot)/(rab2*rbc2)
                                 -(2.0*yab*dot*dot)/(rab2*rab2*rbc2))*terma)/(2.0*termc*sine*sine)      )/ denom;
 
                    dedzib = ( (zbd*terma)/(termc*rbd2) + (xbc*yab-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic))/termc
                              +(zbc*terma)/(termc*rbc2) + (zab*terma)/(termc*rab2)
                              -( (-(2.0*zbc*dot*dot)/(rab2*rbc2*rbc2) - (2.0*(2.0*zib-zia-zic)*dot)/(rab2*rbc2)
                                 -(2.0*zab*dot*dot)/(rab2*rab2*rbc2))*terma)/(2.0*termc*sine*sine)      )/ denom;
                    
                    dedxic = (-(terma*( (2.0*xbc*dot*dot)/(rab2*rbc2*rbc2) - (2.0*xab*dot)/(rab2*rbc2) )/(2.0*termc*sine*sine) )
                             - (xbc*terma)/(rbc2*termc) + (zab*ybd-yab*zbd)/termc ) / denom;
                    dedyic = (-(terma*( (2.0*ybc*dot*dot)/(rab2*rbc2*rbc2) - (2.0*yab*dot)/(rab2*rbc2) )/(2.0*termc*sine*sine) )
                             - (ybc*terma)/(rbc2*termc) + (xab*zbd-zab*xbd)/termc ) / denom;
                    dedzic = (-(terma*( (2.0*zbc*dot*dot)/(rab2*rbc2*rbc2) - (2.0*zab*dot)/(rab2*rbc2) )/(2.0*termc*sine*sine) )
                             - (zbc*terma)/(rbc2*termc) + (yab*xbd-xab*ybd)/termc ) / denom;

                    dedxid = ( (yab*zbc-zab*ybc)/termc - (xbd*terma)/(rbd2*termc)  )/ denom;
                    dedyid = ( (zab*xbc-xab*zbc)/termc - (ybd*terma)/(rbd2*termc)  )/ denom;
                    dedzid = ( (xab*ybc-yab*xbc)/termc - (zbd*terma)/(rbd2*termc)  )/ denom;

                    deriv.deopb[ia][0] += deddt*dedxia;
                    deriv.deopb[ia][1] += deddt*dedyia;
                    deriv.deopb[ia][2] += deddt*dedzia;

                    deriv.deopb[ib][0] += deddt*dedxib;
                    deriv.deopb[ib][1] += deddt*dedyib;
                    deriv.deopb[ib][2] += deddt*dedzib;
                    
                    deriv.deopb[ic][0] += deddt*dedxic;
                    deriv.deopb[ic][1] += deddt*dedyic;
                    deriv.deopb[ic][2] += deddt*dedzic;
                    
                    deriv.deopb[id][0] += deddt*dedxid;
                    deriv.deopb[id][1] += deddt*dedyid;
                    deriv.deopb[id][2] += deddt*dedzid;

                 }
                 
              }
          }
      }
}


double term1aa(double rab,double rcross,double terma,double term1xa,double termc,
              double termss2, double termcrab2)
{
     double term1;
     
     term1 = rcross/termc - (term1xa*terma)/termss2 - (rab*terma)/(termcrab2);
     return(term1);
}
        
void eopbend_wilson2(int iatom)
{
      int i,ia,ib,ic,id;
      double e,angle,dot,cosine;
      double dt,dt2;
      double xia,yia,zia;
      double xib,yib,zib;
      double xic,yic,zic; 
      double xid,yid,zid;
      double xab,yab,zab,rab2,rab;
      double xbc,ybc,zbc,rbc2,rbc;
      double xbd,ybd,zbd,rbd2,rbd;
      double sine, rab2rbc2, term1, term2, term3, term4, term5, term6;
      double term1xa, term1ya, term1za;
      double term1xc, term1yc, term1zc;
      double deddt, terma, denom,termc;
      double dedxiaxia, dedxiaxib, dedxiaxic, dedxiaxid, dedxiayia, dedxiayib, dedxiayic, dedxiayid, dedxiazia, dedxiazib, dedxiazic,dedxiazid;
      double dedxibxib, dedxibxic, dedxibxid, dedxibyia, dedxibyib, dedxibyic, dedxibyid, dedxibzia, dedxibzib, dedxibzic, dedxibzid;
      double dedxicxic, dedxicxid, dedxicyia, dedxicyib, dedxicyic, dedxicyid, dedxiczia, dedxiczib, dedxiczic, dedxiczid;
      double dedxidxid, dedxidyia, dedxidyib, dedxidyic, dedxidyid, dedxidzia, dedxidzib, dedxidzic, dedxidzid;
      double dedyiayia, dedyiayib, dedyiayic, dedyiayid, dedyiazia, dedyiazib, dedyiazic, dedyiazid;
      double dedyibyib, dedyibyic, dedyibyid, dedyibzia, dedyibzib, dedyibzic, dedyibzid;
      double dedyicyic, dedyicyid, dedyiczia, dedyiczib, dedyiczic, dedyiczid;
      double dedyidyid, dedyidzia, dedyidzib, dedyidzic, dedyidzid;
      double dedziazia, dedziazib, dedziazic, dedziazid;
      double dedzibzib, dedzibzic, dedzibzid;
      double dedziczic, dedziczid;
      double dedzidzid;
      double termss, denom2,denom3,dot2,terma2,termss2;
      double termc2, terms4, ddtdenom, ddt1denom, ddt2denom;
      double xab2,yab2,zab2,xbc2,ybc2,zbc2,xbd2,ybd2,zbd2;
      double yabzbc,yabxbd, ybczbd,ybdzab, xbdzbc, xabzbd, xabybc, xbcybd, xbczab;
      double term1xaa, term1zca, termcrab2, termcrbc2, xiabc, yiabc, ziabc;
      double rrab22rbc2, rrab2rbc22;
      double zabybc, xabzbc, yabxbc;
      double termcrbd22;

      for (i=0; i < angles.nang; i++)
      {
          if (angles.angin[i] == TRUE)
          {
              ia = angles.i13[i][0];
              ib = angles.i13[i][1];
              ic = angles.i13[i][2];
              id = angles.i13[i][3];
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
                 
                 xab = xia - xib;
                 yab = yia - yib;
                 zab = zia - zib;
                 rab2 = (xab*xab+yab*yab+zab*zab);
                 rab = sqrt(rab2);
                 
                 xbc = xic - xib;
                 ybc = yic - yib;
                 zbc = zic - zib;
                 rbc2 = xbc*xbc+ybc*ybc+zbc*zbc;
                 rbc = sqrt(rbc2);
                 
                 if (rab2 != 0.0 && rbc2 != 0.0)
                 {
                    dot = xab*xbc + yab*ybc + zab*zbc;
                    dot2 = dot*dot;
                    cosine = dot /(rab*rbc);
                    if (cosine < -1.0)
                       cosine = -1.0;
                    if (cosine > 1.0)
                       cosine = 1.0;
                    sine = sqrt(1.00-cosine*cosine);

                    xbd = atom[id].x - atom[ib].x;
                    ybd = atom[id].y - atom[ib].y;
                    zbd = atom[id].z - atom[ib].z;
                    rbd2 = xbd*xbd+ybd*ybd+zbd*zbd;
                    rbd = sqrt(rbd2);
                    
                    xab2 = xab*xab;
                    yab2 = yab*yab;
                    zab2 = zab*zab;
                    xbc2 = xbc*xbc;
                    ybc2 = ybc*ybc;
                    zbc2 = zbc*zbc;
                    xbd2 = xbd*xbd;
                    ybd2 = ybd*ybd;
                    zbd2 = zbd*zbd;
                    
                    cosine = ( xbd*(yab*zbc-zab*ybc) + ybd*(zab*xbc-xab*zbc) + zbd*(xab*ybc-yab*xbc))/((rab*rbc*rbd)*sine);

                    angle = asin(cosine);
                    dt = angle*radian;
                    dt2 = dt * dt;
                    e = units.angunit * angles.copb[i] * dt2;

                    deddt = 2.0*units.angunit*angles.copb[i]*radian*radian;
                    
                    terma = xbd*(yab*zbc-zab*ybc) + ybd*(zab*xbc-xab*zbc) + zbd*(xab*ybc - yab*xbc);
                    terma2 = terma*terma;
                    termc = rab*rbc*rbd*sine;
                    termc2 = termc*termc;
                    denom = sqrt(1.00 - (terma2)/termc2);
                    denom2 = (1.00 - (terma2)/termc2);
                    denom3 = denom2*denom;
                    ddtdenom = deddt/denom2;
                    ddt1denom = (deddt*angle)/denom;
                    ddt2denom = (0.5*deddt*angle)/denom3;
                    rab2rbc2 = rab2*rbc2;
                    rrab22rbc2 = rab2*rab2rbc2;
                    rrab2rbc22 = rbc2*rab2rbc2;
                    termcrbd22 = termc*rbd2*rbd2;
                    
                    termss = termc*sine*sine;
                    terms4 = 4.0*termss*sine*sine;
                    termss2 = (2.0*termss);
                    term1xa = (2.0*xab*dot2)/(rrab22rbc2) - (2.0*xbc*dot)/rab2rbc2;
                    term1ya = (2.0*yab*dot2)/(rrab22rbc2) - (2.0*ybc*dot)/rab2rbc2;
                    term1za = (2.0*zab*dot2)/(rrab22rbc2) - (2.0*zbc*dot)/rab2rbc2;
                    term1xc = (2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xab*dot)/rab2rbc2;
                    term1yc = (2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yab*dot)/rab2rbc2;
                    term1zc = (2.0*zbc*dot2)/(rrab2rbc22) - (2.0*zab*dot)/rab2rbc2;
                    term1xaa = term1xa*terma;
                    term1zca = term1zc*terma;
                    termcrab2 = termc*rab2;
                    termcrbc2 = termc*rbc2;
                    
                    ybczbd = ybc*zbd-ybd*zbc;
                    ybdzab = ybd*zab-yab*zbd;
                    yabzbc = yab*zbc-ybc*zab;
                    yabxbd = yab*xbd-xab*ybd;
                    xbdzbc = xbd*zbc-xbc*zbd;
                    xabzbd = xab*zbd-zab*xbd;
                    xabybc = xab*ybc-yab*xbc;
                    xbcybd = xbc*ybd-ybc*xbd;
                    xbczab = xbc*zab-xab*zbc;
                    xabzbc = (xab*zbc-zab*xbc+xbd*(zia-zic)+zbd*(xic-xia));
                    yabxbc = (xbc*yab-ybc*xab+xbd*(ybc-yab)+ybd*(xia-xic));
                    zabybc = (ybc*zab-yab*zbc+ybd*(zic-zia)+zbd*(yia-yic));
                    xiabc = (2.0*xib-xia-xic);
                    yiabc = (2.0*yib-yia-yic);
                    ziabc = (2.0*zib-zia-zic);
                             
                    term1 = term1aa(xab, ybczbd, terma, term1xa, termc, termss2, termcrab2);
                    term3 = (8.0*xab*xbc*dot)/(rrab22rbc2)- (2.0*xbc2)/rab2rbc2
                           -(8.0*xab2*dot2)/(rab2*rrab22rbc2) + (2.0*dot2)/(rrab22rbc2);
                    term2 =  - (term1xa*(ybczbd))/(termss)
                             - (2.0*xab*(ybczbd))/(termcrab2)
                             + (3.0*terma*term1xa*term1xa)/(terms4)
                             - (terma*term3)/termss2
                             + (xab*term1xaa)/(termss*rab2)
                             + (3.0*xab2*terma)/(rab2*termcrab2) - terma/(termcrab2);
                    term4 = (- (2.0*terma*(ybczbd))/(termc2)
                             + (term1xaa*terma)/(termc*termss)
                             + (2.0*xab*terma2)/(termc2*rab2) )
                           *(  (ybczbd)/termc - (term1xaa)/termss2
                             - (xab*terma)/(termcrab2));
                    dedxiaxia = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);
                    term3 = (4.0*xab*xbc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*xbc2*dot)/(rrab2rbc22) - (2.0*xiabc*xbc)/rab2rbc2
                             +(4.0*xab*xiabc*dot)/(rrab22rbc2) - (4.0*xab*xbc*dot)/(rrab22rbc2)  + (2.0*dot)/rab2rbc2
                             +(8.0*xab2*dot2)/(rrab22rbc2*rab2) - (2.0*dot2)/(rrab22rbc2);
                    term1 = ( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    term2 =  (xbd*(ybczbd))/(termc*rbd2)
                            -(xbd*term1xaa)/(2.0*termss*rbd2)
                            -(xab*xbd*terma)/(rab2*rbd2*termc)
                            -(term1xa*zabybc)/termss2
                            -(xab*zabybc)/(termcrab2)
                            -(term5*(ybczbd))/termss2
                            +(xbc*(ybczbd))/(termcrbc2)
                            +(xab*(ybczbd))/(termcrab2)
                            +(3.0*term5*term1xaa)/(terms4)
                            -(term3*terma)/termss2
                            +(xab*term5*terma)/(2.0*termss*rab2)
                            -(xbc*term1xaa)/(2.0*termss*rbc2)
                            -(xab*term1xaa)/(2.0*termss*rab2)
                            -(xab*xbc*terma)/(termcrab2*rbc2)
                            -(3.0*xab2*terma)/(termcrab2*rab2)
                            +terma/(termcrab2);
                    term4 = (-(2.0*xbd*terma2)/(termc2*rbd2) - (2.0*zabybc*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*xbc*terma2)/(termc2*rbc2) - (2.0*xab*terma2)/(termc2*rab2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiaxib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xbc2*dot)/(rrab2rbc22) - (4.0*xab*xbc*dot2)/(rab2rbc2*rab2rbc2) - (2.0*xab*xbc)/rab2rbc2
                            +(4.0*xab2*dot)/(rrab22rbc2) - 2.0*dot/rab2rbc2;
                    term1 = ( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    term2 = -(term1xa*(ybdzab))/termss2
                            -(xab*(ybdzab))/(termcrab2)
                            -(term1xc*(ybczbd))/termss2
                            -(xbc*(ybczbd))/(rbc2*termc)
                            +(3.0*term1xa*term1xc*terma)/(terms4)
                            -(term5*terma)/termss2
                            +(xab*term1xc*terma)/(2.0*rab2*termss)
                            +(xbc*term1xaa)/(2.0*rbc2*termss)
                            +(xab*xbc*terma)/(rab2*rbc2*termc);
                    term4 = (-(2.0*(ybdzab)*terma)/(termc2) + (term1xc*terma2)/(termc*termss)
                             +(2.0*xbc*terma2)/(termc2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));        
                    dedxiaxic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc)*((ybczbd)/termc
                             -(term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    term2 =  -(xbd*(ybczbd))/(termc*rbd2)
                             +(xbd*term1xaa)/(2.0*termss*rbd2)
                             +(xab*xbd*terma)/(rab2*rbd2*termc)
                             -((yabzbc)*term1xa)/termss2
                             -(xab*(yabzbc))/(termcrab2);
                    term4 = ( (2.0*xbd*terma2)/(termc2*rbd2) - (2.0*(yabzbc)*terma)/(termc2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2)); 
                    dedxiaxid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xbc*yab*dot)/(rrab22rbc2) + (4.0*xab*ybc*dot)/(rrab22rbc2) - (8.0*xab*yab*dot2)/(rrab22rbc2)
                                     -(2.0*xbc*ybc)/rab2rbc2;
                    term1 = ( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    term2 =  -(term1xa*(xbdzbc))/termss2
                             -(xab*(xbdzbc))/(termcrab2)
                             -(term1ya*(ybczbd))/termss2
                             -(yab*(ybczbd))/(termcrab2)
                             +(3.0*term1xa*term1ya*terma)/(terms4)
                             -(term5*terma)/termss2
                             +(yab*term1xaa)/(2.0*rab2*termss)
                             +(xab*term1ya*terma)/(2.0*rab2*termss)
                             +(3.0*xab*yab*terma)/(rab2*termcrab2);
                    term4 = (-(2.0*(xbdzbc)*terma)/(termc2) + (term1ya*terma2)/(termc*termss)
                             +(2.0*yab*terma2)/(termc2*rab2) )
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiayia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term3 = (4.0*xab*ybc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*xbc*ybc*dot)/(rrab2rbc22) - (2.0*yiabc*xbc)/rab2rbc2
                             +(4.0*xab*yiabc*dot)/(rrab22rbc2) - (4.0*yab*xbc*dot)/(rrab22rbc2)
                             +(8.0*xab*yab*dot2)/(rrab22rbc2*rab2);
                    term1 = ( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2))
                           *( (ybd*terma)/(termc*rbd2) + (xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    term2 =  (ybd*(ybczbd))/(termc*rbd2)
                            -(ybd*term1xaa)/(2.0*termss*rbd2)
                            -(xab*ybd*terma)/(rab2*rbd2*termc)
                            +(zbc-zbd)/termc
                            -(term1xa*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/termss2
                            -(xab*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/(termcrab2)
                            -(term5*(ybczbd))/termss2
                            +(ybc*(ybczbd))/(termcrbc2)
                            +(yab*(ybczbd))/(termcrab2)
                            +(3.0*term5*term1xaa)/(terms4)
                            -(term3*terma)/termss2
                            -(ybc*term1xaa)/(2.0*termss*rbc2)
                            -(yab*term1xaa)/(2.0*termss*rab2)
                            +(xab*term5*terma)/(2.0*termss*rab2)
                            -(xab*ybc*terma)/(termcrab2*rbc2)
                            -(3.0*xab*yab*terma)/(termcrab2*rab2);
                    term4 = (-(2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*ybc*terma2)/(termc2*rbc2) - (2.0*yab*terma2)/(termc2*rab2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiayib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xbc*ybc*dot)/(rrab2rbc22) - (4.0*xab*ybc*dot2)/(rab2rbc2*rab2rbc2) - (2.0*xbc*yab)/rab2rbc2
                            +(4.0*xab*yab*dot)/(rrab22rbc2);
                    term1 = ( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    term2 =  (zbd)/termc
                            -(term1xa*(xabzbd))/termss2
                            -(xab*(xabzbd))/(termcrab2)
                            -(term1yc*(ybczbd))/termss2
                            -(ybc*(ybczbd))/(rbc2*termc)
                            +(3.0*term1yc*term1xaa)/(terms4)
                            -(term5*terma)/termss2
                            +(xab*term1yc*terma)/(2.0*rab2*termss)
                            +(ybc*term1xaa)/(2.0*termss*rbc2)
                            +(xab*ybc*terma)/(termcrab2*rbc2);
                    term4 = (-(2.0*(xabzbd)*terma)/(termc2) + (term1yc*terma2)/(termc*termss)
                             +(2.0*ybc*terma2)/(termc2*rbc2) )
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiayic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*( (ybczbd)/termc
                             -(term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    term2 =  -(ybd*(ybczbd))/(termc*rbd2)
                             +(ybd*term1xaa)/(2.0*termss*rbd2)
                             +(xab*ybd*terma)/(rab2*rbd2*termc)
                             -((xbczab)*term1xa)/termss2
                             -(zbc)/termc
                             -(xab*(xbczab))/(termcrab2);
                    term4 = ( (2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xbczab)*terma)/(termc2) )
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiayid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xab*zbc*dot)/(rrab22rbc2) + (4.0*xbc*zab*dot)/(rrab22rbc2) - (2.0*xbc*zbc)/rab2rbc2
                             -(8.0*xab*zab*dot2)/(rab2*rrab22rbc2);
                    term1 = ( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2))
                           *( (xbcybd)/termc - (term1za*terma)/termss2 - (zab*terma)/(termcrab2));
                    term2 =  -(term1xa*(xbcybd))/termss2
                             -(xab*(xbcybd))/(termcrab2)
                             -(zab*(ybczbd))/(termcrab2)
                             -(term1za*(ybczbd))/termss2
                             +(zab*term1xaa)/(2.0*rab2*termss)
                             +(3.0*xab*zab*terma)/(rab2*termcrab2)
                             -(term5*terma)/termss2
                             +(3.0*term1xa*term1za*terma)/(terms4)
                             +(xab*term1za*terma)/(2.0*rab2*termss);
                    term4 = (-(2.0*(xbcybd)*terma)/(termc2) + (2.0*zab*terma2)/(termc2*rab2)
                             +(term1za*terma2)/(termc*termss))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiazia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                                
                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);                    
                    term3 =  (4.0*xab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*xbc*zbc*dot)/(rrab2rbc22)
                            -(2.0*xbc*ziabc)/rab2rbc2 - (4.0*xbc*zab*dot)/(rrab22rbc2)
                            +(4.0*xab*ziabc*dot)/(rrab22rbc2) + (8.0*xab*zab*dot2)/(rrab22rbc2*rab2);
                    term1 = ( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2))
                           *( (zbd*terma)/(termc*rbd2) + (xbc*yab-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic))/termc
                             +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term5*terma)/termss2);
                    term2 =  (zbd*(ybczbd))/(termc*rbd2)
                            -(zbd*term1xaa)/(2.0*termss*rbd2)
                            -(xab*zbd*terma)/(rab2*rbd2*termc)
                            -(term1xa*(xbc*yab-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic)))/termss2                                    
                            -(xab*(xbc*yab-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic)))/(termcrab2)
                            +(ybd-ybc)/termc
                            +(zbc*(ybczbd))/(termcrbc2)
                            +(zab*(ybczbd))/(termcrab2)
                            -(term5*(ybczbd))/termss2
                            -(zbc*term1xaa)/(2.0*termss*rbc2)
                            -(zab*term1xaa)/(2.0*termss*rab2)
                            -(xab*zbc*terma)/(termcrab2*rbc2)
                            -(3.0*xab*zab*terma)/(termcrab2*rab2) 
                            -(term3*terma)/termss2
                            +(3.0*term5*term1xaa)/(terms4)
                            +(xab*term5*terma)/(2.0*termss*rab2);
                    term4 = (-(2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xbc*yab-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic))*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiazib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xbc*zbc*dot)/(rrab2rbc22) - (4.0*xab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (2.0*xbc*zab)/rab2rbc2
                             + (4.0*xab*zab*dot)/(rrab22rbc2);
                    term1 = ( (yabxbd)/termc - (term1zca)/termss2 - (zbc*terma)/(termcrbc2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    term2 =  -(term1xa*(yabxbd))/termss2
                             -(xab*(yabxbd))/(termcrab2)
                             -ybd/termc
                             -(term1zc*(ybczbd))/termss2
                             -(zbc*(ybczbd))/(termcrbc2)
                             +(3.0*term1zc*term1xaa)/(terms4)
                             -(term5*terma)/termss2
                             +(xab*term1zca)/(2.0*rab2*termss)
                             +(zbc*term1xaa)/(2.0*rbc2*termss)
                             +(xab*zbc*terma)/(rab2*rbc2*termc);
                    term4 = (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                             +(2.0*zbc*terma2)/(rbc2*termc2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiazic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*( (ybczbd)/termc
                             -(term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    term2 =  -(zbd*(ybczbd))/(rbd2*termc)
                             +(term1xa*zbd*terma)/(2.0*termss*rbd2)
                             +(xab*zbd*terma)/(rab2*rbd2*termc)
                             -(term1xa*(xabybc))/termss2
                             -(xab*(xabybc))/(termcrab2)
                             + ybc/termc;
                    term4 = ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                           *( (ybczbd)/termc - (term1xaa)/termss2 - (xab*terma)/(termcrab2));
                    dedxiazid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term3 = -(8.0*xbc2*dot2)/(rrab2rbc22*rbc2) - (8.0*xbc*dot*xiabc)/(rrab2rbc22)                         
                            -(8.0*xab*xbc*dot2)/(rab2rbc2*rab2rbc2) + (2.0*dot2)/(rrab2rbc22)
                            -(2.0*xiabc*xiabc)/rab2rbc2
                            -(8.0*xab*xiabc*dot)/(rrab22rbc2) - (8.0*xab2*dot2)/(rrab22rbc2*rab2)
                            -(4.0*dot)/rab2rbc2 + (2.0*dot2)/(rrab22rbc2);
                    term1 =  (xbd*terma)/(termc*rbd2) + zabybc/termc
                            -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2);
                    term2 =  (3.0*xbd2*terma)/(termcrbd22)
                            +(2.0*xbd*zabybc)/(termc*rbd2)
                            -(xbd*term5*terma)/(termss*rbd2)
                            +(2.0*xbc*xbd*terma)/(termc*rbd2*rbc2)
                            +(2.0*xab*xbd*terma)/(termc*rbd2*rab2)
                            -terma/(termc*rbd2)
                            -(term5*zabybc)/(termc*sine)
                            +(2.0*xbc*zabybc)/(termcrbc2)
                            +(2.0*xab*zabybc)/(termcrab2)                                    
                            +(3.0*term5*term5*terma)/(terms4)
                            -(term3*terma)/termss2
                            -(xbc*term5*terma)/(termss*rbc2)
                            -(xab*term5*terma)/(termss*rab2)
                            +(3.0*xbc2*terma)/(termcrbc2*rbc2)
                            +(2.0*xab*xbc*terma)/(termcrab2*rbc2)
                            -terma/(termcrbc2)
                            +(3.0*xab2*terma)/(termcrab2*rab2)
                            -terma/(termcrab2);
                    term4 =  (-(2.0*xbd*terma2)/(termc2*rbd2) - (2.0*zabybc*terma)/(termc2)
                              +(term5*terma2)/(termc*termss) - (2.0*xbc*terma2)/(termc2*rbc2) - (2.0*xab*terma2)/(termc2*rab2))
                            *( (xbd*terma)/(termc*rbd2) + zabybc/termc - (term5*terma)/(termc*2.0*sine*sine)
                              +(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2)) ;                                                                                       
                    dedxibxib = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4; 

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term3 = (8.0*xbc2*dot2)/(rrab2rbc22*rbc2) - (4.0*xab*xbc*dot)/(rrab2rbc22)
                           +(4.0*xbc*dot*xiabc)/(rrab2rbc22) + (4.0*xab*xbc*dot2)/(rab2rbc2*rab2rbc2)
                           -(2.0*dot2)/(rrab2rbc22) - (2.0*xab*xiabc)/rab2rbc2 - (4.0*xab2*dot)/(rrab22rbc2)
                           +(2.0*dot)/rab2rbc2;
                    term1 = ( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    term2 =  (xbd*(ybdzab))/(termc*rbd2)
                            -(xbd*term1xc*terma)/(2.0*termss*rbd2)
                            -(xbc*xbd*terma)/(termc*rbd2*rbc2)
                            -(term5*(ybdzab))/(termc*2.0*sine*sine)
                            +(xbc*(ybdzab))/(termcrbc2)
                            +(xab*(ybdzab))/(termcrab2)
                            -(term1xc*zabybc)/(termc*2.0*sine*sine)
                            -(xbc*zabybc)/(termcrbc2)
                            +(3.0*term1xc*term5*terma)/(terms4)
                            -(term3*terma)/termss2
                            -(xbc*term1xc*terma)/(2.0*termss*rbc2)
                            -(xab*term1xc*terma)/(2.0*termss*rab2)
                            +(xbc*term5*terma)/(2.0*termss*rbc2)
                            -(3.0*xbc2*terma)/(termcrbc2*rbc2)
                            -(xab*xbc*terma)/(termcrab2*rbc2)
                            +terma/(termcrbc2);
                    term4 = (-(2.0*(ybdzab)*terma)/(termc2) + (term1xc*terma2)/(termc*termss)
                             +(2.0*xbc*terma2)/(termc2*rbc2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibxic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term1 = (-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc)*( (xbd*terma)/(termc*rbd2)
                             +zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    term2 = -(3.0*xbd2*terma)/(termcrbd22)
                            +(xbd*(yabzbc))/(termc*rbd2)
                            -(xbd*zabybc)/(termc*rbd2)
                            +(xbd*term5*terma)/(2.0*termss*rbd2)
                            -(xbd*xbc*terma)/(termcrbc2*rbd2)
                            -(xab*xbd*terma)/(termcrab2*rbd2)
                            +terma/(termc*rbd2)
                            -(term5*(yabzbc))/termss2
                            +(xbc*(yabzbc))/(termcrbc2)
                            +(xab*(yabzbc))/(termcrab2);
                    term4 = ( (2.0*xbd*terma2)/(termc2*rbd2) - (2.0*(yabzbc)*terma)/(termc2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibxid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term3 = (8.0*xab*yab*dot2)/(rrab22rbc2*rab2) - (4.0*ybc*xab*dot)/(rrab22rbc2)
                            +(4.0*yab*dot*xiabc)/(rrab22rbc2) + (4.0*yab*xbc*dot2)/(rab2rbc2*rab2rbc2)
                            -(2.0*ybc*xiabc)/(rab2rbc2) - (4.0*xbc*ybc*dot)/(rrab2rbc22);
                    term1 = ( (xbd*terma)/(termc*rbd2) + (ybc*zab-yab*zbc+ybd*(zic-zia)+zbd*(yia-yic) )/termc
                             -(term5*terma)/termss2+(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    term2 =  (xbd*(xbdzbc))/(termc*rbd2)
                            -(xbd*term1ya*terma)/(2.0*termss*rbd2)
                            -(xbd*yab*terma)/(termcrab2*rbd2)
                            +(zbd-zbc)/termc
                            -(term5*(xbdzbc))/termss2
                            +(xbc*(xbdzbc))/(termcrbc2)
                            +(xab*(xbdzbc))/(termcrab2)
                            -(term1ya*zabybc)/termss2
                            -(yab*zabybc)/(termcrab2)
                            +(3.0*term1ya*term5*terma)/(terms4)
                            -(term3*terma)/termss2
                            +(yab*term5*terma)/(2.0*termss*rab2)
                            -(xbc*term1ya*terma)/(2.0*termss*rbc2)
                            -(xab*term1ya*terma)/(2.0*termss*rab2)
                            -(3.0*xab*yab*terma)/(termcrab2*rab2)
                            -(yab*xbc*terma)/(termcrbc2*rab2);
                    term4 = (-(2.0*(xbdzbc)*terma)/(termc2) + (term1ya*terma2)/(termc*termss)
                             +(2.0*yab*terma2)/(termc2*rab2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibyia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                                     
                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term6 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                                      
                    term3 = -(8.0*xbd*ybc*dot2)/(rrab2rbc22*rbc2) - (4.0*xbd*dot*(yiabc))/(rrab2rbc22)
                            -(4.0*ybc*dot*xiabc)/(rrab2rbc22) - (4.0*xbd*yab*dot2)/(rab2rbc2*rab2rbc2)
                            -(4.0*xab*ybc*dot2)/(rab2rbc2*rab2rbc2) - (2.0*xiabc*yiabc)/rab2rbc2
                            -(4.0*yab*dot*xiabc)/(rrab22rbc2) - (4.0*xab*yiabc*dot)/(rrab22rbc2)
                            -(8.0*xab*yab*dot2)/(rrab22rbc2*rab2);
                    term1 = ( (xbd*terma)/(termc*rbd2) + (ybc*zab-yab*zbc+ybd*(zic-zia)+zbd*(yia-yic) )/termc
                             -(term5*terma)/termss2+(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2))
                           *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             -(term6*terma)/termss2+(ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    term2 =  (3.0*xbd*ybd*terma)/(termcrbd22)
                            +(xbd*xabzbc)/(termc*rbd2)
                            +(ybd*zabybc)/(termc*rbd2)
                            -(ybd*term5*terma)/(2.0*termss*rbd2)
                            -(xbd*term6*terma)/(2.0*termss*rbd2)
                            +(xbd*ybc*terma)/(termcrbc2*rbd2)
                            +(xbd*ybd*terma)/(termcrbc2*rbd2)
                            +(xbd*yab*terma)/(termcrab2*rbd2)
                            +(xab*ybd*terma)/(termcrab2*rbd2)
                            -(term5*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/termss2
                            +(xbc*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/(termcrbc2)
                            +(xab*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/(termcrab2)
                            -(term6*zabybc)/termss2
                            +(ybc*zabybc)/(termcrbc2)
                            +(yab*zabybc)/(termcrab2)
                            +(3.0*term5*term6*terma)/(terms4)
                            -(term3*terma)/termss2
                            -(ybc*term5*terma)/(2.0*termss*rbc2)
                            -(yab*term5*terma)/(2.0*termss*rab2)
                            -(xbc*term6*terma)/(2.0*termss*rbc2)
                            -(xab*term6*terma)/(2.0*termss*rab2)
                            +(3.0*xbc*ybc*terma)/(termcrbc2*rbc2)
                            +(xbc*yab*terma)/(termcrab2*rbc2)
                            +(xab*ybc*terma)/(termcrab2*rbc2)
                            +(3.0*xab*yab*terma)/(termcrab2*rab2);
                    term4 = (-(2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))*terma)/(termc2)
                             +(term6*terma2)/(termc*termss) - (2.0*ybc*terma2)/(termc2*rbc2) - (2.0*yab*terma2)/(termc2*rab2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbd*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibyib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term3 = (8.0*xbc*ybc*dot2)/(rrab2rbc22*rbc2) - (4.0*xbc*yab*dot)/(rrab2rbc22)
                            +(4.0*ybc*dot*xiabc)/(rrab2rbc22) + (4.0*xab*ybc*dot2)/(rab2rbc2*rab2rbc2)
                            -(2.0*yab*xiabc)/(rab2rbc2) - (4.0*xab*yab*dot)/(rrab22rbc2);
                    term1 = ( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2))
                           *( (xbd*terma)/(termc*rbd2) + (ybc*zab-yab*zbc+ybd*(zic-zia)+zbd*(yia-yic) )/termc
                             -(term5*terma)/termss2+(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2) );                         
                    term2 =  (xbd*(xabzbd))/(termc*rbd2)
                            -(xbd*term1yc*terma)/(2.0*termss*rbd2)
                            -(xbd*ybc*terma)/(termcrbc2*rbd2)
                            +(zab-zbd)/termc
                            -(term5*(xabzbd))/termss2
                            +(xbc*(xabzbd))/(termcrbc2)
                            +(xab*(xabzbd))/(termcrab2)
                            -(term1yc*zabybc)/termss2
                            -(ybc*zabybc)/(termcrbc2)
                            +(3.0*term1yc*term5*terma)/(terms4)
                            -(term3*terma)/termss2
                            -(xbc*term1yc*terma)/(2.0*termss*rbc2)
                            -(xab*term1yc*terma)/(2.0*termss*rab2)
                            +(ybc*term5*terma)/(2.0*termss*rbc2)
                            -(3.0*xbc*ybc*terma)/(termcrbc2*rbc2)
                            -(xab*ybc*terma)/(termcrbc2*rab2);
                    term4 = (-(2.0*(xabzbd)*terma)/(termc2) + (term1yc*terma2)/(termc*termss)
                             +(2.0*ybc*terma2)/(termc2*rbc2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibyic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term1 = (-(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*( (xbd*terma)/(termc*rbd2)
                             +zabybc/termc - (term5*terma)/termss2
                             +(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    term2 =  -(3.0*xbd*ybd*terma)/(termcrbd22)
                             +(xbd*(xbczab))/(termc*rbd2)
                             -(ybd*zabybc)/(termc*rbd2)
                             +(ybd*term5*terma)/(2.0*termss*rbd2)
                             -(xbc*ybd*terma)/(termcrbc2*rbd2)
                             -(xab*ybd*terma)/(termcrab2*rbd2)
                             -((xbczab)*term5)/termss2
                             +(xbc*(xbczab))/(termcrbc2)
                             +(zbc-zab)/termc
                             +(xab*(xbczab))/(termcrab2);
                    term4 = ( (2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xbczab)*terma)/(termc2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibyid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term3 = (8.0*xab*zab*dot2)/(rrab22rbc2*rab2) - (4.0*zbc*xab*dot)/(rrab22rbc2)
                            +(4.0*zab*dot*xiabc)/(rrab22rbc2) + (4.0*zab*xbc*dot2)/(rab2rbc2*rab2rbc2)
                            -(2.0*zbc*xiabc)/(rab2rbc2) - (4.0*xbc*zbc*dot)/(rrab2rbc22);
                    term1 = ( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2+(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2))
                           *( (xbcybd)/termc - (term1za*terma)/termss2 - (zab*terma)/(termcrab2));
                    term2 =  (xbd*(xbcybd))/(termc*rbd2)
                            -(xbd*zab*terma)/(termcrab2*rbd2)
                            -(xbd*term1za*terma)/(2.0*termss*rbd2)
                            -(term5*(xbcybd))/termss2
                            +(xbc*(xbcybd))/(termcrbc2)
                            +(xab*(xbcybd))/(termcrab2)
                            +(ybc-ybd)/termc
                            -(zab*zabybc)/(termcrab2)
                            -(term1za*zabybc)/termss2
                            +(zab*term5*terma)/(2.0*termss*rab2)
                            -(xbc*zab*terma)/(termcrbc2*rab2)
                            -(3.0*xab*zab*terma)/(termcrab2*rab2)
                            -(term3*terma)/termss2                                    
                            +(3.0*term1za*term5*terma)/(terms4)                                    
                            -(xbc*term1za*terma)/(2.0*termss*rbc2)
                            -(xab*term1za*terma)/(2.0*termss*rab2);
                    term4 = (-(2.0*(xbcybd)*terma)/(termc2) + (term1za*terma2)/(termc*termss)
                             +(2.0*zab*terma2)/(termc2*rab2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibzia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term6 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);                    
                    term3 = -(8.0*xbc*zbc*dot2)/(rrab2rbc22*rbc2) - (4.0*xbc*dot*(ziabc))/(rrab2rbc22)
                            -(4.0*zbc*dot*xiabc)/(rrab2rbc22) - (4.0*xbc*zab*dot2)/(rab2rbc2*rab2rbc2)
                            -(4.0*xab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (2.0*xiabc*ziabc)/rab2rbc2
                            -(4.0*zab*dot*xiabc)/(rrab22rbc2) - (4.0*xab*ziabc*dot)/(rrab22rbc2)
                            -(8.0*xab*zab*dot2)/(rrab22rbc2*rab2);
                    term1 = ( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2+(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2))
                           *( (zbd*terma)/(termc*rbd2) + (xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                             +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2)-(term6*terma)/termss2);
                    term2 =   (3.0*xbd*zbd*terma)/(termcrbd22)
                             +(xbd*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termc*rbd2)
                             +(zbd*zabybc)/(termc*rbd2)
                             +(xbd*zbc*terma)/(termcrbc2*rbd2)
                             +(xbd*zab*terma)/(termcrab2*rbd2)
                             -(xbd*term6*terma)/(2.0*termss*rbd2)                                    
                             -(zbd*term5*terma)/(2.0*termss*rbd2)                                    
                             +(xbc*zbd*terma)/(termcrbc2*rbd2)
                             +(xab*zbd*terma)/(termcrab2*rbd2)                                    
                             -(term5*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/termss2
                             +(xbc*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrbc2)
                             +(xab*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrab2)
                             +(zbc*zabybc)/(termcrbc2)
                             +(zab*zabybc)/(termcrab2)                                    
                             -(term6*zabybc)/termss2
                             -(zbc*term5*terma)/(2.0*termss*rbc2)
                             -(zab*term5*terma)/(2.0*termss*rab2)
                             +(3.0*xbc*zbc*terma)/(termcrbc2*rbc2)
                             +(xbc*zab*terma)/(termcrab2*rbc2)
                             +(xab*zbc*terma)/(termcrab2*rbc2)
                             +(3.0*xab*zab*terma)/(termcrab2*rab2)
                             -(term3*terma)/termss2                                    
                             +(3.0*term5*term6*terma)/(terms4)                                    
                             -(xbc*term6*terma)/(2.0*termss*rbc2)
                             -(xab*term6*terma)/(2.0*termss*rab2);
                    term4 = (-(2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))*terma)/(termc2)
                             +(term6*terma2)/(termc*termss) - (2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibzib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term3 =  (8.0*xbc*zbc*dot2)/(rrab2rbc22*rbc2) - (4.0*xbc*zab*dot)/(rrab2rbc22)
                            +(4.0*zbc*dot*xiabc)/(rrab2rbc22) + (4.0*xab*zbc*dot2)/(rab2rbc2*rab2rbc2)
                            -(2.0*zab*xiabc)/(rab2rbc2) - (4.0*xab*zab*dot)/(rrab22rbc2);
                    term1 = ( (yabxbd)/termc - (term1zca)/termss2 - (zbc*terma)/(termcrbc2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2+(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    term2 =   (xbd*(yabxbd))/(termc*rbd2)
                             -(xbd*term1zca)/(2.0*termss*rbd2)
                             -(xbd*zbc*terma)/(termcrbc2*rbd2)
                             -(term5*(yabxbd))/termss2
                             +(xbc*(yabxbd))/(termcrbc2)
                             +(xab*(yabxbd))/(termcrab2)
                             +(ybd-yab)/termc
                             -(term1zc*zabybc)/termss2
                             -(zbc*zabybc)/(termcrbc2)
                             +(3.0*term1zc*term5*terma)/(terms4)
                             -(term3*terma)/termss2
                             -(xbc*term1zca)/(2.0*termss*rbc2)
                             -(xab*term1zca)/(2.0*termss*rab2)
                             +(zbc*term5*terma)/(2.0*termss*rbc2)
                             -(3.0*xbc*zbc*terma)/(termcrbc2*rbc2)
                             -(xab*zbc*terma)/(termcrbc2*rab2);
                    term4 =  (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                              +(2.0*zbc*terma2)/(termc2*rbc2))
                            *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                              -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibzic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*xbc*dot2)/(rrab2rbc22) - (2.0*xiabc*dot)/rab2rbc2 - (2.0*xab*dot2)/(rrab22rbc2);                    
                    term1 = (-(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*( (xbd*terma)/(termc*rbd2)
                             +zabybc/termc - (term5*terma)/termss2
                             +(xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    term2 =  -(3.0*xbd*zbd*terma)/(termcrbd22)
                             +(xbd*(xabybc))/(termc*rbd2)
                             -(zbd*zabybc)/(termc*rbd2)
                             +(zbd*term5*terma)/(2.0*termss*rbd2)
                             -(xbc*zbd*terma)/(termcrbc2*rbd2)
                             -(xab*zbd*terma)/(termcrab2*rbd2)
                             -((xabybc)*term5)/termss2
                             +(xbc*(xabybc))/(termcrbc2)
                             +(yab-ybc)/termc
                             +(xab*(xabybc))/(termcrab2);
                    term4 = ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                           *( (xbd*terma)/(termc*rbd2) + zabybc/termc
                             -(term5*terma)/termss2 + (xbc*terma)/(termcrbc2) + (xab*terma)/(termcrab2));
                    dedxibzid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                                   
                    term3 = (8.0*xab*xbc*dot)/(rrab2rbc22)- (2.0*xab2)/rab2rbc2 - (8.0*xbc2*dot2)/(rbc2*rrab2rbc22)
                             + (2.0*dot2)/(rrab2rbc22);
                    term1 = term1aa(xbc, ybdzab, terma, term1xc, termc, termss2, termcrbc2);
                    term2 = - (term1xc*(ybdzab))/(termss)
                            - (2.0*xbc*(ybdzab))/(termcrbc2)
                            + (3.0*terma*term1xc*term1xc)/(terms4)
                            - (terma*term3)/termss2
                            + (xbc*term1xc*terma)/(termss*rbc2)
                            + (3.0*xbc2*terma)/(rbc2*rbc2*termc)
                            - terma/(rbc2*termc);
                    term4 = (-(2.0*terma*(ybdzab))/(termc2)
                             +(term1xc*terma2)/(termc*termss)
                             +(2.0*xbc*terma2)/(rbc2*termc2) )
                           *((ybdzab)/termc - (term1xc*terma)/termss2
                            -(xbc*terma)/(rbc2*termc));                             
                    dedxicxic = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4; 


                    term1 = (-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc)*( (ybdzab)/termc
                             -(term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2));
                    term2 = -(xbd*(ybdzab))/(termc*rbd2)
                            +(xbd*term1xc*terma)/(2.0*termss*rbd2)
                            +(xbc*xbd*terma)/(rbc2*rbd2*termc)
                            -((yabzbc)*term1xc)/termss2
                            -(xbc*(yabzbc))/(rbc2*termc);
                    term4 = ( (2.0*xbd*terma2)/(termc2*rbd2) - (2.0*(yabzbc)*terma)/(termc2) )
                           *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(rbc2*termc));
                    dedxicxid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xbc*ybc*dot)/(rrab2rbc22) - (4.0*xbc*yab*dot2)/(rab2rbc2*rab2rbc2) - (2.0*xab*ybc)/(rab2rbc2)
                              + (4.0*xab*yab*dot)/(rrab22rbc2);
                    term1 = ( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));        
                    term2 = -(zbd)/termc
                            -(term1xc*(xbdzbc))/termss2
                            -(xbc*(xbdzbc))/(rbc2*termc)
                            -(term1ya*(ybdzab))/termss2
                            -(yab*(ybdzab))/(termcrab2)
                            +(3.0*term1xc*term1ya*terma)/(terms4)
                            -(term5*terma)/termss2
                            +(yab*term1xc*terma)/(2.0*rab2*termss)
                            +(xbc*term1ya*terma)/(2.0*rbc2*termss)
                            +(xbc*yab*terma)/(rab2*rbc2*termc);
                    term4 = (-(2.0*(xbdzbc)*terma)/(termc2) + (term1ya*terma2)/(termc*termss)
                             +(2.0*yab*terma2)/(termc2*rab2))
                           *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(rbc2*termc));
                    dedxicyia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term3 = (4.0*yab*xbc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*xab*ybc*dot)/(rrab2rbc22) - (2.0*yiabc*xab)/rab2rbc2
                             +(4.0*xbc*yiabc*dot)/(rrab22rbc2) - (4.0*xab*yab*dot)/(rrab22rbc2)
                             +(8.0*xbc*ybc*dot2)/(rrab2rbc22*rbc2);
                    term1 = ( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2))
                           *( (ybd*terma)/(termc*rbd2) + (xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    term2 =  (ybd*(ybdzab))/(termc*rbd2)
                            -(ybd*term1xc*terma)/(2.0*termss*rbd2)
                            -(xbc*ybd*terma)/(rbc2*rbd2*termc)
                            +(zbd-zab)/termc
                            -(term1xc*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/termss2
                            -(xbc*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/(termcrbc2)
                            -(term5*(ybdzab))/termss2
                            +(ybc*(ybdzab))/(termcrbc2)
                            +(yab*(ybdzab))/(termcrab2)
                            +(3.0*term5*term1xc*terma)/(terms4)
                            -(term3*terma)/termss2
                            -(ybc*term1xc*terma)/(2.0*termss*rbc2)
                            -(yab*term1xc*terma)/(2.0*termss*rab2)
                            +(xbc*term5*terma)/(termss*rbc2)
                            -(yab*xbc*terma)/(termcrab2*rbc2)
                            -(3.0*xbc*ybc*terma)/(termcrbc2*rbc2);
                    term4 = (-(2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*ybc*terma2)/(termc2*rbc2) - (2.0*yab*terma2)/(termc2*rab2))
                           *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(rbc2*termc));
                    dedxicyib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xbc*yab*dot)/(rrab2rbc22) - (8.0*xbc*ybc*dot2)/(rrab2rbc22*rbc2) + (4.0*xab*ybc*dot)/(rrab2rbc22)
                              - (2.0*xab*yab)/(rab2rbc2);
                    term1 = ( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2))
                           *( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2));
                    term2 =   -(term1xc*(xabzbd))/termss2
                              -(xbc*(xabzbd))/(termcrbc2)
                              -(term1yc*(ybdzab))/termss2
                              -(ybc*(ybdzab))/(termcrbc2)
                              +(3.0*term1xc*term1yc*terma)/(terms4)
                              -(term5*terma)/termss2
                              +(ybc*term1xc*terma)/(2.0*rbc2*termss)
                              +(xbc*term1yc*terma)/(2.0*termss*rbc2)
                              +(3.0*xbc*ybc*terma)/(termcrbc2*rbc2);
                    term4 = (-(2.0*(xabzbd)*terma)/(termc2) + (term1yc*terma2)/(termc*termss)
                             +(2.0*ybc*terma2)/(termc2*rbc2))
                           *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2));
                    dedxicyic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*( (ybdzab)/termc
                             -(term1yc*terma)/termss2 - (xbc*terma)/(termcrbc2));
                    term2 =  -(ybd*(ybdzab))/(termc*rbd2)
                             +(ybd*term1xc*terma)/(2.0*termss*rbd2)
                             +(xbc*ybd*terma)/(termcrbc2*rbd2)
                             -((xbczab)*term1xc)/termss2 
                             -(xbc*(xbczab))/(termcrbc2)
                             +zab/termc;
                    term4 = ( (2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xbczab)*terma)/(termc2))
                           *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(rbc2*termc));
                    dedxicyid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xbc*zbc*dot)/(rrab2rbc22) + (4.0*xab*zab*dot)/(rrab22rbc2) - (2.0*xab*zbc)/rab2rbc2
                             -(4.0*xbc*zab*dot2)/(rab2rbc2*rab2rbc2);
                    term1 = ( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2))
                           *( (xbcybd)/termc - (term1za*terma)/termss2 - (zab*terma)/(termcrab2));
                    term2 =  -(term1xc*(xbcybd))/termss2
                             -(xbc*(xbcybd))/(termcrbc2)
                             +ybd/termc
                             -(zab*(ybdzab))/(termcrab2)
                             -(term1za*(ybdzab))/termss2
                             +(zab*term1xc*terma)/(2.0*rab2*termss)
                             -(term5*terma)/termss2
                             +(xbc*zab*terma)/(rab2*rbc2*termc)
                             +(3.0*term1xc*term1za*terma)/(terms4)
                             +(xbc*term1za*terma)/(2.0*rbc2*termss);
                    term4 = (-(2.0*(xbcybd)*terma)/(termc2) + (2.0*zab*terma2)/(termc2*rab2)
                             +(term1za*terma2)/(termc*termss))
                           *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(rbc2*termc));
                    dedxiczia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);                    
                    term3 = (4.0*zab*xbc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*xab*zbc*dot)/(rrab2rbc22) - (2.0*ziabc*xab)/rab2rbc2
                             +(4.0*xbc*ziabc*dot)/(rrab22rbc2) - (4.0*zab*xab*dot)/(rrab22rbc2)
                             +(8.0*xbc*zbc*dot2)/(rrab22rbc2*rab2);
                    term1 =  ( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2))
                            *( (zbd*terma)/(termc*rbd2) +  (xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                              +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term5*terma)/termss2);
                    term2 =  (zbd*(ybdzab))/(termc*rbd2)
                            -(zbd*term1xc*terma)/(2.0*termss*rbd2)
                            -(xbc*zbd*terma)/(rbc2*rbd2*termc)
                            -(term1xc*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/termss2
                            -(xbc*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrbc2)
                            +(yab-ybd)/termc
                            +(zbc*(ybdzab))/(termcrbc2)
                            +(zab*(ybdzab))/(termcrab2)
                            -(term5*(ybdzab))/termss2
                            -(zbc*term1xc*terma)/(2.0*termss*rbc2)
                            -(zab*term1xc*terma)/(2.0*termss*rab2)
                            -(term3*terma)/termss2
                            -(3.0*xbc*zbc*terma)/(termcrbc2*rbc2)
                            -(zab*xbc*terma)/(termcrab2*rbc2)
                            +(3.0*term5*term1xc*terma)/(terms4)
                            +(xbc*term5*terma)/(termss*rbc2);
                    term4 = (-(2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2))
                           *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(rbc2*termc));
                    dedxiczib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*xbc*zab*dot)/(rrab2rbc22) + (4.0*xab*zbc*dot)/(rrab2rbc22) - (2.0*xab*zab)/rab2rbc2
                             - (8.0*xbc*zbc*dot2)/(rbc2*rrab2rbc22);
                    term1 = ( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2))
                           *( (yabxbd)/termc - (term1zca)/termss2 - (zbc*terma)/(termcrbc2));
                    term2 =  -(term1xc*(yabxbd))/termss2
                             -(xbc*(yabxbd))/(rbc2*termc)
                             -(term1zc*(ybdzab))/termss2
                             -(zbc*(ybdzab))/(termcrbc2)
                             +(3.0*term1zc*term1xc*terma)/(terms4)
                             -(term5*terma)/(2.0*termss*rbc2)
                             +(xbc*term1zca)/(2.0*rbc2*termss)
                             +(zbc*term1xc*terma)/(2.0*rbc2*termss)
                             +(3.0*xbc*zbc*terma)/(rbc2*rbc2*termc);
                    term4 = (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                             +(2.0*zbc*terma2)/(rbc2*termc2))
                           *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(rbc2*termc));
                    dedxiczic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*( (ybdzab)/termc
                             -(term1xc*terma)/termss2 - (xbc*terma)/(termcrbc2));
                    term2 =   -(zbd*(ybdzab))/(rbd2*termc)
                              +(term1xc*zbd*terma)/(2.0*termss*rbd2)
                              +(xbc*zbd*terma)/(rbc2*rbd2*termc)
                              -(term1xc*(xabybc))/termss2
                              -(xbc*(xabybc))/(rbc2*termc)
                              -yab/termc;
                    term4 =  ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                            *( (ybdzab)/termc - (term1xc*terma)/termss2 - (xbc*terma)/(rbc2*termc));
                    dedxiczid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = -(xbd*terma)/(termc*rbd2) + (yabzbc)/termc;
                    term2 =  (3.0*xbd2*terma)/(termcrbd22) - (2.0*xbd*(yabzbc))/(termc*rbd2)
                            -(terma)/(termc*rbd2);
                    term4 = ((2.0*xbd*terma2)/(termc2*rbd2) - (2.0*(yabzbc)*terma)/(termc2))
                            *(-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc );                                 
                    dedxidxid = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;


                    term1 = (-(xbd*terma)/(termc*rbd2)+(yabzbc)/termc)*( (xbdzbc)/termc
                             -(term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    term2 = -(xbd*(xbdzbc))/(termc*rbd2)
                            +(xbd*term1ya*terma)/(2.0*termss*rbd2)
                            +(xbd*yab*terma)/(termcrab2*rbd2)
                            -(term1ya*(yabzbc))/termss2
                            +(zbc)/termc
                            -(yab*(yabzbc))/(termcrab2);
                    term4 = (-(2.0*(xbdzbc)*terma)/(termc2) + (term1ya*terma2)/(termc*termss)
                             +(2.0*yab*terma2)/(termc2*rab2))
                           *(  -(xbd*terma)/(termc*rbd2) + (yabzbc)/termc );
                    dedxidyia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term1 = (-(xbd*terma)/(termc*rbd2)+(yabzbc)/termc)*( (ybd*terma)/(termc*rbd2)
                             +(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))/termc - (term5*terma)/termss2
                             +(ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    term2 = -(3.0*xbd*ybd*terma)/(termcrbd22)
                            +(ybd*(yabzbc))/(termc*rbd2)
                            -(xbd*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/(termc*rbd2)
                            +(xbd*term5*terma)/(2.0*termss*rbd2)
                            -(xbd*ybc*terma)/(termcrbc2*rbd2)
                            -(xbd*yab*terma)/(termcrab2*rbd2)
                            -(term5*(yabzbc))/termss2
                            +(ybc*(yabzbc))/(termcrbc2)
                            +(zab-zbc)/termc
                            +(yab*(yabzbc))/(termcrab2);
                    term4 = (-(2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*ybc*terma2)/(termc2*rbc2) - (2.0*yab*terma2)/(termc2*rab2))
                           *(-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc );
                    dedxidyib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(xbd*terma)/(termc*rbd2)+(yabzbc)/termc)*( (xabzbd)/termc
                             -(term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2));
                    term2 = -(xbd*(xabzbd))/(termc*rbd2)
                            +(xbd*term1yc*terma)/(2.0*termss*rbd2)
                            +(xbd*ybc*terma)/(termcrbc2*rbd2)
                            -(term1yc*(yabzbc))/termss2
                            -(ybc*(yabzbc))/(rbc2*termc)
                            -zab/termc;
                    term4 = (-(2.0*(xabzbd)*terma)/(termc2) + (term1yc*terma2)/(termc*termss)
                             +(2.0*ybc*terma2)/(rbc2*termc2))
                           *(-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc );       
                    dedxidyic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*( -(xbd*terma)/(termc*rbd2)
                             +(yabzbc)/termc);
                    term2 =  (3.0*xbd*ybd*terma)/(termcrbd22)
                            -(xbd*(xbczab))/(termc*rbd2)
                            -(ybd*(yabzbc))/(termc*rbd2);
                    term4 = ( (2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xbczab)*terma)/(termc2))
                           *(-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc);
                    dedxidyid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(xbd*terma)/(termc*rbd2)+(yabzbc)/termc)*( (xbcybd)/termc
                             -(zab*terma)/(termcrab2) - (term1za*terma)/termss2);
                    term2 =  -(xbd*(xbcybd))/(termc*rbd2)
                             +(xbd*zab*terma)/(termcrab2*rbd2)
                             +(xbd*term1za*terma)/(2.0*termss*rbd2)
                             -(ybc)/termc
                             -(zab*(yabzbc))/(termcrab2)
                             -(term1za*(yabzbc))/termss2;
                    term4 =  (-(2.0*(xbd*ybc-ybd*xbc)*terma)/(termc2) + (term1za*terma2)/(termc*termss)
                              +(2.0*zab*terma2)/(termc2*rab2))
                            *(-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc );
                    dedxidzia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);                    
                    term1 = (-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc)*( (zbd*terma)/(termc*rbd2)
                             +(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc + (zbc*terma)/(termcrbc2)
                             +(zab*terma)/(termcrab2) - (term5*terma)/termss2);
                    term2 =  -(3.0*xbd*zbd*terma)/(termcrbd22)
                             -(xbd*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termc*rbd2)
                             +(zbd*(yabzbc))/(termc*rbd2)
                             -(xbd*zbc*terma)/(termcrbc2*rbd2)
                             -(xbd*zab*terma)/(termcrab2*rbd2)
                             +(xbd*term5*terma)/(2.0*termss*rbd2)
                             +(zbc*(yabzbc))/(termcrbc2)
                             +(ybc-yab)/termc
                             +(zab*(yabzbc))/(termcrab2)
                             -(term5*(yabzbc))/termss2;
                    term4 = (-(2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2))
                           *(-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc );
                    dedxidzib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(xbd*terma)/(termc*rbd2)+(yabzbc)/termc)*( (yabxbd)/termc
                             -(term1zca)/termss2 - (zbc*terma)/(termcrbc2));
                    term2 =  -(xbd*(yabxbd))/(termc*rbd2)
                             +(xbd*term1zca)/(2.0*termss*rbd2)
                             +(xbd*zbc*terma)/(termcrbc2*rbd2)
                             -(term1zc*(yabzbc))/termss2
                             -(zbc*(yabzbc))/(rbc2*termc)
                             +yab/termc;
                    term4 = (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                             +(2.0*zbc*terma2)/(rbc2*termc2))
                           *(-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc );
                    dedxidzic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*(-(xbd*terma)/(termc*rbd2)
                             +(yabzbc)/termc);
                    term2 =   (3.0*xbd*zbd*terma)/(termcrbd22)
                             -(xbd*(ybc*xab-yab*xbc))/(termc*rbd2)
                             -(zbd*(yabzbc))/(termc*rbd2);
                    term4 =  ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(ybc*xab-yab*xbc)*terma)/(termc2))
                            *(-(xbd*terma)/(termc*rbd2) + (yabzbc)/termc );
                    dedxidzid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                                     
                    term3 = (8.0*yab*ybc*dot)/(rrab22rbc2)- (2.0*ybc2)/rab2rbc2 - (8.0*yab2*dot2)/(rab2*rrab22rbc2)
                             + (2.0*dot2)/(rrab22rbc2);
                    term1 = term1aa(yab, xbdzbc, terma, term1ya, termc, termss2, termcrab2);
                    term2 = -(term1ya*(xbdzbc))/(termss)
                            -(2.0*yab*(xbdzbc))/(termcrab2)
                            +(3.0*terma*term1ya*term1ya)/(terms4)
                            -(terma*term3)/termss2
                            +(yab*term1ya*terma)/(termss*rab2)
                            +(3.0*yab2*terma)/(rab2*termcrab2)
                            -terma/(termcrab2);
                    term4 = (-(2.0*terma*(xbdzbc))/(termc2)
                             +(term1ya*terma2)/(termc*termss)
                             +(2.0*yab*terma2)/(termc2*rab2) )
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2
                             -(yab*terma)/(termcrab2));
                    dedyiayia = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term3 = (4.0*yab*ybc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*ybc2*dot2)/(rrab2rbc22) - (2.0*yiabc*ybc)/rab2rbc2
                             +(4.0*yab*yiabc*dot)/(rrab22rbc2) - (4.0*yab*ybc*dot)/(rrab22rbc2) +(2.0*dot)/rab2rbc2
                             +(8.0*yab2*dot2)/(rrab22rbc2*rab2) - (2.0*dot2)/(rrab2rbc22);             
                    term1 = ( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2))
                           *( (ybd*terma)/(termc*rbd2) + (xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    term2 =  (ybd*(xbdzbc))/(termc*rbd2)
                            -(ybd*term1ya*terma)/(2.0*termss*rbd2)
                            -(yab*ybd*terma)/(rab2*rbd2*termc) 
                            -(term5*(xbdzbc))/termss2
                            +(ybc*(xbdzbc))/(termcrbc2)
                            +(yab*(xbdzbc))/(termcrab2)
                            -(term1ya*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/termss2
                            -(yab*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia)))/(termcrab2)
                            +(3.0*term5*term1ya*terma)/(terms4)
                            -(term3*terma)/termss2
                            +(yab*term5*terma)/(2.0*termss*rab2)
                            -(ybc*term1ya*terma)/(2.0*termss*rbc2)
                            -(yab*term1ya*terma)/(2.0*termss*rab2)
                            -(yab*ybc*terma)/(termcrab2*rbc2)
                            -(3.0*yab2*terma)/(termcrab2*rab2)
                            +terma/(termcrab2);
                    term4 = (-(2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xab*zbc-xbc*zab+xbd*(zia-zic)+zbd*(xic-xia))*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*ybc*terma2)/(termc2*rbc2)
                             -(2.0*yab*terma2)/(termc2*rab2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    dedyiayib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*ybc2*dot)/(rrab2rbc22) + (4.0*yab2*dot)/(rrab22rbc2) - (4.0*yab*ybc*dot2)/(rab2rbc2*rab2rbc2)
                            - (2.0*yab*ybc)/rab2rbc2;
                    term1 = ( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    term2 =  -(term1ya*(xabzbd))/termss2
                             -(yab*(xabzbd))/(termcrab2)
                             -(term1yc*(xbdzbc))/termss2
                             -(ybc*(xbdzbc))/(rbc2*termc)
                             +(3.0*term1yc*term1ya*terma)/(terms4)
                             -(term5*terma)/termss2
                             +(yab*term1yc*terma)/(2.0*termss*rab2)
                             +(ybc*term1ya*terma)/(2.0*termss*rbc2)
                             +(yab*ybc*terma)/(rab2*rbc2*termc);
                    term4 = (-(2.0*(xabzbd)*terma)/(termc2) + (term1yc*terma2)/(termc*termss)
                             +(2.0*ybc*terma2)/(rbc2*termc2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    dedyiayic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                                     
                    term1 =  (-(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*( (xbdzbc)/termc
                              -(term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    term2 =  -(ybd*(xbdzbc))/(termc*rbd2)
                             +(ybd*term1ya*terma)/(2.0*termss*rbd2)
                             +(yab*ybd*terma)/(rab2*rbd2*termc)
                             -(term1ya*(xbczab))/termss2
                             -(yab*(xbczab))/(termcrab2);
                    term4 = ( (2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xbczab)*terma)/(termc2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    dedyiayid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*ybc*zab*dot)/(rrab22rbc2) + (4.0*yab*zbc*dot)/(rrab22rbc2) - (8.0*yab*zab*dot2)/(rrab22rbc2*rab2)
                            - (2.0*zbc*ybc)/rab2rbc2;
                    term1 = ( (xbdzbc)/termc - (term1ya*terma)/termss2 -(yab*terma)/(termcrab2))
                           *( (xbcybd)/termc - (zab*terma)/(termcrab2) - (term1za*terma)/termss2);
                    term2 =  -(term1ya*(ybd*xbc-xbd*ybc))/termss2
                             -(yab*(ybd*xbc-xbd*ybc))/(termcrab2)
                             -(zab*(xbdzbc))/(termcrab2)
                             -(term1za*(xbdzbc))/termss2
                             +(zab*term1ya*terma)/(2.0*termss*rab2)
                             +(3.0*yab*zab*terma)/(termcrab2*rab2)
                             -(term5*terma)/termss2
                             +(3.0*term1ya*term1za*terma)/(terms4)
                             +(yab*term1za*terma)/(2.0*termss*rab2);
                    term4 = (-(2.0*(ybd*xbc-xbd*ybc)*terma)/(termc2) + (2.0*zab*terma2)/(termc2*rab2)
                             +(term1za*terma2)/(termc*termss))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    dedyiazia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);                    
                    term3 = (4.0*yab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*ybc*zbc*dot)/(rrab2rbc22) - (2.0*ziabc*ybc)/rab2rbc2
                             +(4.0*yab*ziabc*dot)/(rrab22rbc2) - (4.0*zab*ybc*dot)/(rrab22rbc2)
                             +(8.0*yab*zab*dot2)/(rrab22rbc2*rab2);
                    term1 = ( (xbdzbc)/termc - (term1ya*terma)/termss2 -(yab*terma)/(termcrab2))
                           *( (zbd*terma)/(termc*rbd2) + (xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                             +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term5*terma)/termss2);
                    term2 =   (zbd*(xbdzbc))/(termc*rbd2)
                             -(term1ya*zbd*terma)/(2.0*termss*rbd2)
                             -(yab*zbd*terma)/(termcrab2*rbd2)
                             -(term1ya*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/termss2
                             -(yab*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrab2)
                             +(xbc-xbd)/termc
                             +(zbc*(xbdzbc))/(termcrbc2)
                             +(zab*(xbdzbc))/(termcrab2)
                             -(term5*(xbdzbc))/termss2
                             -(zbc*term1ya*terma)/(2.0*termss*rbc2)
                             -(zab*term1ya*terma)/(2.0*termss*rab2)
                             -(yab*zbc*terma)/(termcrab2*rbc2)
                             -(3.0*yab*zab*terma)/(termcrab2*rab2)
                             -(term3*terma)/termss2
                             +(3.0*term1ya*term5*terma)/(terms4)
                             +(yab*term5*terma)/(2.0*termss*rab2);
                   term4 = (-(2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))*terma)/(termc2)
                            -(2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2)
                            +(term5*terma2)/(termc*termss))
                          *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    dedyiazib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*ybc*zbc*dot)/(rrab2rbc22) + (4.0*yab*zab*dot)/(rrab22rbc2) - (4.0*yab*zbc*dot2)/(rab2rbc2*rab2rbc2)
                            - (2.0*zab*ybc)/rab2rbc2;
                    term1 = ( (yabxbd)/termc - (term1zca)/termss2 - (zbc*terma)/(termcrbc2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    term2 =  -(term1ya*(yabxbd))/termss2
                             -(yab*(yabxbd))/(termcrab2)
                             + xbd/termc
                             -(term1zc*(xbdzbc))/termss2
                             -(zbc*(xbdzbc))/(rbc2*termc)
                             +(3.0*term1zc*term1ya*terma)/(terms4)
                             -(term5*terma)/termss2
                             +(yab*term1zca)/(2.0*termss*rab2)
                             +(zbc*term1ya*terma)/(2.0*termss*rbc2)
                             +(yab*zbc*terma)/(rab2*rbc2*termc);
                    term4 = (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                             +(2.0*zbc*terma2)/(rbc2*termc2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    dedyiazic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(zbd*terma)/(termc*rbd2) + (xab*ybc-xbc*yab)/termc)*( (xbdzbc)/termc
                             -(term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    term2 =  -(zbd*(xbdzbc))/(termc*rbd2)
                             +(zbd*term1ya*terma)/(2.0*termss*rbd2)
                             +(yab*zbd*terma)/(rab2*rbd2*termc)
                             -(term1ya*(xabybc))/termss2
                             -(yab*(xabybc))/(termcrab2)
                             - xbc/termc;
                    term4 = ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                           *( (xbdzbc)/termc - (term1ya*terma)/termss2 - (yab*terma)/(termcrab2));
                    dedyiazid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                 
                    term3 = -(8.0*ybc2*dot2)/(rrab2rbc22*rbc2)
                            -(8.0*yiabc*ybc*dot)/(rrab2rbc22) + (2.0*dot2)/(rrab2rbc22)
                            -(8.0*yab*ybc*dot2)/(rab2rbc2*rab2rbc2) - (8.0*yab*yiabc)/(rab2rbc2)
                            -(4.0*dot)/rab2rbc2 - (8.0*yab2*dot2)/(rrab22rbc2) + (2.0*dot2)/(rrab22rbc2)
                            -(2.0*yiabc*yiabc)/rab2rbc2;
                    term1 = (ybd*terma)/(termc*rbd2) + xabzbc/termc
                           -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2);
                    term2 =  (3.0*ybd2*terma)/(termcrbd22)
                            +(2.0*ybd*xabzbc)/(termc*rbd2)
                            -(term5*ybd*terma)/(termss*rbd2)
                            +(2.0*ybc*ybd*terma)/(termc*rbd2*rbc2)
                            +(2.0*yab*ybd*terma)/(termc*rbd2*rab2)
                            -terma/(termc*rbd2)
                            -(term5*xabzbc)/(termss)
                            +(2.0*ybc*xabzbc)/(termcrbc2)
                            +(2.0*yab*xabzbc)/(termcrab2)
                            +(3.0*term5*term5*terma)/(terms4)
                            -(term3*terma)/termss2 
                            -(ybc*term5*terma)/(termss*rbc2)
                            -(yab*term5*terma)/(termss*rab2)
                            +(3.0*ybc2*terma)/(termcrbc2*rbc2)
                            +(2.0*yab*ybc*terma)/(termcrab2*rbc2)
                            +(3.0*yab2*terma)/(termcrab2*rab2)
                            -terma/(termcrab2) - terma/(termcrbc2);
                    term4 = ( -(2.0*ybd*terma2)/(termc2*rbd2)
                              -(2.0*xabzbc*terma)/(termc2)
                              +(term5*terma2)/(termc*termss) - (2.0*ybc*terma2)/(termc2*rbc2)
                              -(2.0*yab*terma2)/(termc2*rab2))
                           *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             +(ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2) - (term5*terma)/termss2);                                                      
                    dedyibyib = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;
                                    
                    term3 = (8.0*ybc2*dot2)/(rrab2rbc22*rbc2) - (4.0*yab*ybc*dot)/(rrab2rbc22)
                           +(4.0*ybc*dot*yiabc)/(rrab2rbc22) + (4.0*yab*ybc*dot2)/(rab2rbc2*rab2rbc2)
                           -(2.0*dot2)/(rrab2rbc22) - (2.0*yab*yiabc)/rab2rbc2 - (4.0*yab2*dot)/(rrab22rbc2)
                           +(2.0*dot)/rab2rbc2;
                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term1 = ( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2))
                           *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));                            
                    term2 =  (ybd*(xabzbd))/(termc*rbd2)
                            -(ybd*term1yc*terma)/(2.0*termss*rbd2)
                            -(ybc*ybd*terma)/(termcrbc2*rbd2)
                            -(term5*(xabzbd))/termss2
                            +(ybc*(xabzbd))/(termcrbc2)
                            +(yab*(xabzbd))/(termcrab2)
                            -(term1yc*xabzbc)/termss2
                            -(ybc*xabzbc)/(termcrbc2)
                            +(3.0*term1yc*term5*terma)/(terms4)
                            -(term3*terma)/termss2
                            -(ybc*term1yc*terma)/(2.0*termss*rbc2)
                            -(yab*term1yc*terma)/(2.0*termss*rab2)
                            +(ybc*term5*terma)/(2.0*termss*rbc2)
                            -(3.0*ybc2*terma)/(termcrbc2*rbc2)
                            -(yab*ybc*terma)/(termcrbc2*rab2)
                            +(terma)/(termcrbc2);
                    term4 = (-(2.0*(xabzbd)*terma)/(termc2) + (term1yc*terma2)/(termc*termss)
                             +(2.0*ybc*terma2)/(termc2*rbc2))
                           *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             +(ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2) - (term5*terma)/termss2);
                    dedyibyic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term1 = (-(ybd*terma)/(termc*rbd2) + (xbczab)/termc )*( (ybd*terma)/(termc*rbd2)
                             +xabzbc/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    term2 =  -(3.0*ybd2*terma)/(termcrbd22)
                             +(ybd*(xbczab))/(termc*rbd2)
                             -(ybd*xabzbc)/(termc*rbd2)
                             +(ybd*term5*terma)/(2.0*termss*rbd2)
                             -(ybc*ybd*terma)/(termcrbc2*rbd2)
                             -(yab*ybd*terma)/(termcrab2*rbd2)
                             +terma/(termc*rbd2)
                             -(term5*(xbczab))/termss2
                             +(ybc*(xbczab))/(termcrbc2)
                             +(yab*(xbczab))/(termcrab2);
                    term4 = ( (2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xbczab))/(termc2))
                           *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             +(ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2) - (term5*terma)/termss2);                            
                    dedyibyid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term3 = (4.0*zab*ybc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*zbc*ybc*dot)/(rrab2rbc22) - (2.0*yiabc*zbc)/rab2rbc2
                             +(4.0*zab*yiabc*dot)/(rrab22rbc2) - (4.0*yab*zbc*dot)/(rrab22rbc2)
                             +(8.0*zab*yab*dot2)/(rrab22rbc2*rab2);
                    term1 = ( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2))
                           *( (xbcybd)/termc - (zab*terma)/(termcrab2) - (term1za*terma)/termss2);
                    term2 =   (ybd*(xbcybd))/(termc*rbd2)
                             -(ybd*zab*terma)/(termcrab2*rbd2)
                             -(ybd*term1za*terma)/(2.0*termss*rbd2)
                             -(term5*(xbcybd))/termss2
                             +(ybc*(xbcybd))/(termcrbc2)
                             +(yab*(xbcybd))/(termcrab2)
                             +(xbd-xbc)/termc
                             -(zab*xabzbc)/(termcrab2)
                             -(term1za*xabzbc)/termss2
                             +(zab*term5*terma)/(2.0*termss*rab2)
                             -(ybc*zab*terma)/(termcrab2*rbc2)
                             -(3.0*yab*zab*terma)/(termcrab2*rab2)
                             -(term3*terma)/termss2
                             +(3.0*term5*term1za*terma)/(terms4)
                             -(ybc*term1za*terma)/(2.0*termss*rbc2)
                             -(yab*term1za*terma)/(2.0*termss*rab2);
                    term4 =  (-(2.0*(xbcybd)*terma)/(termc2) + (2.0*zab*terma2)/(termc2*rab2)
                              +(term1za*terma2)/(termc*termss))
                            *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                              -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    dedyibzia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term6 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);                    
                    term3 = -(8.0*ybc*zbc*dot2)/(rrab2rbc22*rbc2) - (4.0*ybc*ziabc*dot)/(rrab2rbc22)
                            -(4.0*zbc*dot*yiabc)/(rrab2rbc22) -(4.0*ybc*zab*dot2)/(rab2rbc2*rab2rbc2)
                            -(4.0*yab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (2.0*yiabc*ziabc)/rab2rbc2
                            -(4.0*zab*dot*yiabc)/(rrab22rbc2) - (4.0*yab*dot*ziabc)/(rrab22rbc2)
                            -(8.0*yab*zab*dot2)/(rrab22rbc2*rab2);
                    term1 = ( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2))
                           *( (zbd*terma)/(termc*rbd2) + (yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                             +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term6*terma)/termss2);
                    term2 =   (3.0*ybd*zbd*terma)/(termcrbd22)
                             +(ybd*(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termc*rbd2)
                             +(zbd*xabzbc)/(termc*rbd2)
                             +(ybd*zbc*terma)/(termc*rbd2*rbc2)
                             +(ybd*zab*terma)/(termc*rbd2*rab2)
                             -(ybd*term6*terma)/(2.0*termss*rbd2)
                             -(zbd*term5*terma)/(2.0*termss*rbd2)
                             +(ybc*zbd*terma)/(termcrbc2*rbd2)
                             +(yab*zbd*terma)/(termcrab2*rbd2)
                             -(term5*(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/termss2
                             +(ybc*(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrbc2)
                             +(yab*(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrab2)
                             +(zbc*xabzbc)/(termcrbc2)
                             +(zab*xabzbc)/(termcrab2)
                             -(term6*xabzbc)/termss2
                             -(zbc*term5*terma)/(2.0*termss*rbc2)
                             -(zab*term5*terma)/(2.0*termss*rab2)
                             +(3.0*ybc*zbc*terma)/(termcrbc2*rbc2)
                             +(ybc*zab*terma)/(termcrbc2*rab2)
                             +(yab*zbc*terma)/(termcrbc2*rab2)
                             +(3.0*yab*zab*terma)/(termcrab2*rab2)
                             -(term3*terma)/termss2
                             +(3.0*term5*term6*terma)/(terms4)
                             -(ybc*term6*terma)/(2.0*termss*rbc2)
                             -(yab*term6*terma)/(2.0*termss*rab2);
                    term4 =  (-(2.0*zbd*terma2)/(termc2) - (2.0*(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))*terma)/(termc2)
                              -(2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2)
                              +(term6*terma2)/(termc*termss))
                            *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                              -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    dedyibzib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                                    
                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term3 = (4.0*yab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*yab*zab*dot)/(rrab2rbc22) - (2.0*yiabc*zab)/rab2rbc2
                             +(4.0*zbc*yiabc*dot)/(rrab2rbc22) - (4.0*zab*ybc*dot)/(rrab2rbc22)
                             +(8.0*zbc*ybc*dot2)/(rrab2rbc22*rbc2);
                    term1 = ( (xbd*yab-ybd*xab)/termc - (term1zca)/termss2 - (zbc*terma)/(termcrbc2))
                           *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    term2 = (ybd*(yabxbd))/(termc*rbd2)
                           -(ybd*term1zca)/(2.0*termss*rbd2)
                           -(ybd*zbc*terma)/(termcrbc2*rbd2)
                           -(term5*(yabxbd))/termss2
                           +(ybc*(yabxbd))/(termcrbc2)
                           +(yab*(yabxbd))/(termcrab2)
                           +(xab-xbd)/termc
                           -(term1zc*xabzbc)/termss2
                           -(zbc*xabzbc)/(termcrbc2)
                           +(3.0*term5*term1zca)/(terms4)
                           -(term3*terma)/termss2
                           -(ybc*term1zca)/(2.0*termss*rbc2)
                           -(yab*term1zca)/(2.0*termss*rab2)                                   
                           +(zbc*term5*terma)/(2.0*termss*rbc2)
                           -(3.0*ybc*zbc*terma)/(termcrbc2*rbc2)                                    
                           -(zbc*yab*terma)/(termcrab2*rbc2);
                    term4 =  (-(2.0*(yabxbd)*terma)/(termc2) + (2.0*zbc*terma2)/(termc2*rbc2)
                              +(term1zca*terma)/(termc*termss))
                            *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                              -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    dedyibzic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*ybc*dot2)/(rrab2rbc22) - (2.0*yiabc*dot)/rab2rbc2 - (2.0*yab*dot2)/(rrab22rbc2);                    
                    term1 = (-(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*((ybd*terma)/(termc*rbd2)
                             +xabzbc/termc - (term5*terma)/termss2
                             +(ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    term2 =  -(3.0*ybd*zbd*terma)/(termcrbd22)
                             +(ybd*(xabybc))/(termc*rbd2)
                             -(zbd*xabzbc)/(termc*rbd2)
                             +(term5*zbd*terma)/(2.0*termss*rbd2)
                             -(ybc*zbd*terma)/(termcrbc2*rbd2)
                             -(yab*zbd*terma)/(termcrab2*rbd2)
                             -(term5*(xabybc))/termss2
                             +(ybc*(xabybc))/(termcrbc2)
                             +(yab*(xabybc))/(termcrab2)
                             +(xbc-xab)/termc;
                    term4 = ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                           *( (ybd*terma)/(termc*rbd2) + xabzbc/termc
                             -(term5*terma)/termss2 + (ybc*terma)/(termcrbc2) + (yab*terma)/(termcrab2));
                    dedyibzid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                                     
                    term3 = (8.0*yab*ybc*dot)/(rrab2rbc22)- (2.0*yab2)/rab2rbc2 - (8.0*ybc2*dot2)/(rbc2*rrab2rbc22)
                             + (2.0*dot2)/(rrab2rbc22);
                    term1 = term1aa(ybc, xabzbd, terma, term1yc, termc, termss2, termcrbc2);
                    term2 = -(term1yc*(xabzbd))/(termss)
                            -(2.0*ybc*(xabzbd))/(termcrbc2)
                            +(3.0*terma*term1yc*term1yc)/(terms4)
                            -(terma*term3)/termss2
                            +(ybc*term1yc*terma)/(termss*rbc2)
                            +(3.0*ybc2*terma)/(rbc2*rbc2*termc) - terma/(rbc2*termc);
                    term4 = (-(2.0*terma*(xabzbd))/(termc2)
                             +(term1yc*terma2)/(termc*termss)
                             +(2.0*ybc*terma2)/(rbc2*termc2) )
                           *( (xabzbd)/termc - (term1yc*terma)/termss2
                             -(ybc*terma)/(rbc2*termc));
                    dedyicyic = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*( (xabzbd)/termc
                             -(term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2));
                    term2 = -(ybd*(xabzbd))/(termc*rbd2)
                            +(ybd*term1yc*terma)/(2.0*termss*rbd2)
                            +(ybc*ybd*terma)/(rbc2*rbd2*termc)
                            -(term1yc*(xbczab))/termss2
                            -(ybc*(xbczab))/(rbc2*termc);
                    term4 = ( (2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xbczab)*terma)/(termc2))
                            *( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrab2));
                    dedyicyid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*ybc*zbc*dot)/(rrab2rbc22) + (4.0*yab*zab*dot)/(rrab22rbc2) - (4.0*ybc*zab*dot2)/(rab2rbc2*rab2rbc2)
                            - (2.0*zbc*yab)/rab2rbc2;
                    term1 = ( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2))
                           *( (xbcybd)/termc - (zab*terma)/(termcrbc2) - (term1za*terma)/termss2);
                    term2 =  -(term1yc*(ybd*xbc-xbd*ybc))/termss2
                             -(ybc*(ybd*xbc-xbd*ybc))/(rbc2*termc)
                             - xbd/(termc)
                             -(zab*(xabzbd))/(termcrab2)
                             -(term1za*(xabzbd))/termss2
                             +(zab*term1yc*terma)/(2.0*termss*rab2)
                             +(ybc*zab*terma)/(termcrab2*rab2)
                             -(term5*terma)/termss2
                             +(ybc*zab*terma)/(termcrab2*rbc2)
                             +(3.0*term1yc*term1za*terma)/(terms4)
                             +(ybc*term1za*terma)/(2.0*termss*rab2);
                    term4 =  (-(2.0*(ybd*xbc-xbd*ybc)*terma)/(termc2) + (2.0*zab*terma2)/(termc2*rab2)
                              +(term1za*terma2)/(termc*termss))
                            *( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(rbc2*termc));
                    dedyiczia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);
                    term3 =  (8.0*ybc*zbc*dot2)/(rrab2rbc22*rbc2) + (4.0*ybc*ziabc*dot)/(rrab2rbc22)
                             -(4.0*yab*zbc*dot)/(rrab2rbc22) + (4.0*zab*ybc*dot2)/(rab2rbc2*rab2rbc2)
                             -(4.0*yab*zab*dot)/(rrab2rbc22) - (2.0*ziabc*yab)/rab2rbc2;
                    term1 = ( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2))
                           *( (zbd*terma)/(termc*rbd2) + (xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                             +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term5*terma)/termss2);
                    term2 =   (zbd*(xabzbd))/(termc*rbd2)
                             -(term1yc*zbd*terma)/(2.0*termss*rbd2)
                             -(ybc*zbd*terma)/(termcrbc2*rbd2)
                             -(term1yc*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/termss2
                             -(ybc*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrab2)
                             +(xbd-xab)/termc
                             +(zbc*(xabzbd))/(termcrbc2)
                             +(zab*(xabzbd))/(termcrab2)
                             -(term5*(xabzbd))/termss2
                             -(zbc*term1yc*terma)/(2.0*termss*rbc2)
                             -(zab*term1yc*terma)/(2.0*termss*rab2)
                             -(term3*terma)/termss2
                             -(3.0*ybc*zbc*terma)/(termcrbc2*rbc2)
                             -(ybc*zab*terma)/(termcrab2*rbc2)                                    
                             +(3.0*term1yc*term5*terma)/(terms4)                                    
                             +(ybc*term5*terma)/(2.0*termss*rab2);
                    term4 = (-(2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))*terma)/(termc2)
                             -(2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2)
                             +(term5*terma2)/(termc*termss))
                          *(  (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2));          
                    dedyiczib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = (4.0*ybc*zab*dot)/(rrab2rbc22) + (4.0*yab*zbc*dot)/(rrab2rbc22) - (8.0*ybc*zbc*dot2)/(rrab2rbc22*rbc2)
                            - (2.0*zab*yab)/rab2rbc2;
                    term1 = ( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2))
                           *( (yabxbd)/termc - (term1zca)/termss2 - (zbc*terma)/(termcrbc2));
                    term2 =  -(term1yc*(yabxbd))/termss2
                             -(ybc*(yabxbd))/(rbc2*termc)
                             -(term1zc*(xabzbd))/termss2
                             -(zbc*(xabzbd))/(rbc2*termc)
                             +(3.0*term1zc*term1yc*terma)/(terms4)
                             -(term5*terma)/termss2
                             +(ybc*term1zca)/(2.0*termss*rbc2)
                             +(zbc*term1yc*terma)/(2.0*termss*rbc2)
                             +(3.0*ybc*zbc*terma)/(rbc2*rbc2*termc);
                    term4 = (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                             +(2.0*zbc*terma2)/(rbc2*termc2))
                           *( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(rbc2*termc));
                    dedyiczic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = ( -(zbd*terma)/(termc*rbd2) + (xab*ybc-xbc*yab)/termc)*( (xabzbd)/termc
                              -(term1yc*terma)/termss2 - (ybc*terma)/(termcrbc2));
                    term2 =   -(zbd*(xabzbd))/(termc*rbd2)
                              +(zbd*term1yc*terma)/(2.0*termss*rbd2)
                              +(ybc*zbd*terma)/(rbc2*rbd2*termc)
                              -(term1yc*(xabybc))/termss2
                              -(ybc*(xabybc))/(rbc2*termc) + xab/termc;
                    term4 =  ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                            *( (xabzbd)/termc - (term1yc*terma)/termss2 - (ybc*terma)/(rbc2*termc));
                    dedyiczid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = -(ybd*terma)/(termc*rbd2) + (xbczab)/termc;
                    term2 = (3.0*ybd2*terma)/(termcrbd22) - (2.0*ybd*(xbczab))/(termc*rbd2)
                             -(terma)/(termc*rbd2);
                    term4 = ( (2.0*ybd*terma2)/(termc2*rbd2) - (2.0*(xbczab)*terma)/(termc2))
                           *(-(ybd*terma)/(termc*rbd2) + (xbczab)/termc);
                    dedyidyid = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*((xbcybd)/termc
                             -(zab*terma)/(termcrab2) - (term1za*terma)/termss2);
                    term2 =  -(ybd*(xbcybd))/(termc*rbd2) + (ybd*zab*terma)/(rab2*rbd2*termc)
                             +(ybd*term1za*terma)/(2.0*termss*rbd2) + xbc/termc
                             -(zab*(xbczab))/(termcrab2)
                             -((xbczab)*term1za)/termss2;
                    term4 = (-(2.0*(xbcybd)*terma)/(termc2) + (2.0*zab*terma2)/(termc2*rab2)
                             +(term1za*terma2)/(termc*termss))
                           *(-(ybd*terma)/(termc*rbd2) + (xbczab)/termc ) ;      
                    dedyidzia = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                               
                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);                    
                    term1 = ( -(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*( (zbd*terma)/(termc*rbd2)
                             +(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc + (zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2)
                             -(term5*terma)/termss2 );
                    term2 =  -(3.0*ybd*zbd*terma)/(termcrbd22)
                             -(ybd*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termc*rbd2)
                             +(zbd*(xbczab))/(termc*rbd2)
                             -(ybd*zbc*terma)/(termcrbc2*rbd2)
                             -(ybd*zab*terma)/(termcrab2*rbd2)
                             +(ybd*term5*terma)/(2.0*termss*rbd2)
                             +(zbc*(xbczab))/(termcrbc2)
                             +(zab*(xbczab))/(termcrab2)
                             +(xab-xbc)/termc
                             -(term5*(xbczab))/termss2;
                    term4 =  (-(2.0*zbd*terma2)/(termc2*rbd2)
                              -(2.0*terma*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termc2)
                              -(2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2)
                              +(term5*terma2)/(termc*termss*sine))
                            *(-(ybd*terma)/(termc*rbd2) + (xbczab)/termc );
                    dedyidzib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = ( -(ybd*terma)/(termc*rbd2) + (xbczab)/termc)*((yabxbd)/termc
                              -(term1zca)/termss2 - (zbc*terma)/(termcrbc2) );
                    term2 =  -(ybd*(yabxbd))/(termc*rbd2)
                             +(ybd*term1zca)/(2.0*termss*rbd2)
                             +(ybd*zbc*terma)/(rbc2*rbd2*termc)
                             -((xbczab)*term1zc)/termss2
                             -(zbc*(xbczab))/(termcrbc2)
                             -xab/termc ;
                    term4 = (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                             +(2.0*zbc*terma2)/(termc2*rbc2))
                           *(-(ybd*terma)/(termc*rbd2) + (xbczab)/termc );
                    dedyidzic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = ( -(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*(-(ybd*terma)/(termc*rbd2)
                              +(xbczab)/termc);
                    term2 =  (3.0*ybd*zbd*terma)/(termcrbd22) - (ybd*(xabybc))/(termc*rbd2)
                            -(zbd*(xbczab))/(termc*rbd2);
                    term4 = ((2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                           *(-(ybd*terma)/(termc*rbd2) + (xbczab)/termc );
                    dedyidzid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
 
                    term3 = (8.0*zab*zbc*dot)/(rrab22rbc2)- (2.0*zbc2)/rab2rbc2 - (8.0*zab2*dot2)/(rab2*rrab22rbc2)
                             + (2.0*dot2)/(rrab22rbc2);
                    term1 = term1aa(zab, xbcybd, terma, term1za, termc, termss2, termcrab2);
                    term2 = -(2.0*zab*(xbcybd))/(termcrab2)
                            -(term1za*(xbcybd))/(termss)
                            -(terma)/(termcrab2)
                            +(3.0*zab2*terma)/(rab2*termcrab2)
                            +(zab*term1za*terma)/(termss*rab2)
                            +(3.0*terma*term1za*term1za)/(terms4)
                            -(term3*terma)/termss2;
                    term4 = ( -(2.0*terma*(xbcybd))/(termc2)
                              +(term1za*terma2)/(termc*termss)
                              +(2.0*zab*terma2)/(termc2*rab2) )
                           *(  (xbcybd)/termc - (term1za*terma)/termss2
                              -(zab*terma)/(termcrab2));
                    dedziazia = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;
                    
                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2);                    
                    term3 = (4.0*zab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (4.0*zbc2*dot)/(rrab2rbc22) - (2.0*ziabc*zbc)/rab2rbc2
                             +(4.0*zab*ziabc*dot)/(rrab22rbc2) - (4.0*zab*zbc*dot)/(rrab22rbc2)
                             +(8.0*zab2*dot2)/(rrab22rbc2*rab2) + (2.0*dot)/rab2rbc2 - (2.0*dot2)/(rrab22rbc2);
                    term1 = ( (zbd*terma)/(termc*rbd2) + (xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                             +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term5*terma)/termss2)
                           *( (xbcybd)/termc - (zab*terma)/(termcrab2) - (term1za*terma)/termss2);
                    term2 =  (zbd*(xbcybd))/(termc*rbd2)
                            -(zab*zbd*terma)/(rab2*rbd2*termc)
                            -(zbd*term1za*terma)/(2.0*termss*rbd2)
                            +(zbc*(xbcybd))/(termcrbc2)
                            -(zab*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrab2)
                            +(zab*(xbcybd))/(termcrab2)
                            -(term5*(xbcybd))/termss2
                            -(term1za*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/termss2
                            -(zab*zbc*terma)/(termcrab2*rbc2)
                            +(terma)/(termcrab2)
                            -(3.0*zab2*terma)/(termcrab2*rab2)
                            +(zab*term5*terma)/(2.0*termss*rab2)
                            -(zbc*term1za*terma)/(2.0*termss*rbc2)
                            -(zab*term1za*terma)/(2.0*termss*rab2)
                            +(3.0*term5*term1za*terma)/(terms4)
                            -(term3*terma)/termss2;
                    term4 = (-(2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xbc*yab-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))*terma)/(termc2)
                             +(term5*terma2)/(termc*termss) - (2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2))
                           *( (xbcybd)/termc - (term1za*terma)/termss2
                             -(zab*terma)/(termcrab2));
                    dedziazib = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term3 = (4.0*zbc2*dot)/(rab2rbc2*rbd2) - (4.0*zab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (2.0*zab*zbc)/rab2rbc2
                            - (2.0*dot)/rab2rbc2 + (4.0*zab2*dot)/(rrab22rbc2);
                    term1 = ( (yabxbd)/termc - (term1zca)/termss2 - (zbc*terma)/(termcrbc2))
                           *( (xbcybd)/termc - (zab*terma)/(termcrab2) - (term1za*terma)/termss2);
                    term2 =  -((xbcybd)*term1zc)/termss2
                             -(zbc*(xbcybd))/(termcrbc2)
                             -(zab*(yabxbd))/(termcrab2)
                             -(term1za*(yabxbd))/termss2
                             +(zab*term1zca)/(2.0*termss*rab2)
                             -(term3*terma)/termss2
                             +(zab*zbc*terma)/(termcrab2*rbc2)
                             +(3.0*term1zc*term1za*terma)/(terms4)
                             +(zbc*term1za*terma)/(2.0*termss*rbc2);
                    term4 = (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                             +(2.0*zbc*terma2)/(termc2*rbc2))
                           *( (xbcybd)/termc - (zab*terma)/(termcrab2)  - (term1za*terma)/termss2);
                    dedziazic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*( (xbcybd)/termc
                             -(zab*terma)/(termcrab2) - (term1za*terma)/termss2 );
                    term2 = -(zbd*(xbcybd))/(termc*rbd2)
                            +(zab*zbd*terma)/(termcrab2*rbd2)
                            +(zbd*term1za*terma)/(2.0*termss*rbd2)
                            -(zab*(xabybc))/(termcrab2)
                            -((xabybc)*term1za)/termss2;
                    term4 = ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                           *( (xbcybd)/termc - (zab*terma)/(termcrab2) - (term1za*terma)/termss2);
                    dedziazid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2); 
                    term3 = -(8.0*zbc2*dot2)/(rrab2rbc22*rbc2)
                            -(8.0*ziabc*zbc*dot)/(rrab2rbc22) + (2.0*dot2)/(rrab2rbc22)
                            -(8.0*zab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (8.0*zab*ziabc)/(rab2rbc2)
                            -(4.0*dot)/rab2rbc2 - (8.0*zab2*dot2)/(rrab22rbc2) + (2.0*dot2)/(rrab22rbc2)
                            -(2.0*ziabc*ziabc)/rab2rbc2;                            
                    term1 =  (zbd*terma)/(termc*rbd2) + yabxbc/termc
                            +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term5*terma)/termss2;                            
                    term2 =  (3.0*zbd2*terma)/(termcrbd22)
                            +(2.0*zbd*(yab*xbc-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic)))/(termc*rbd2)
                            -terma/(termc*rbd2)
                            +(2.0*zbc*zbd*terma)/(termc*rbd2*rbc2)
                            +(2.0*zab*zbd*terma)/(termc*rbd2*rab2)
                            -(term5*zbd*terma)/(termss*rbd2)
                            +(2.0*zbc*(yab*xbc-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic)))/(termcrbc2)
                            +(2.0*zab*(yab*xbc-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic)))/(termcrab2)
                            -(term5*(yab*xbc-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic)))/(termss)                                    
                            +(3.0*zbc2*terma)/(termcrbc2*rbc2)
                            -terma/(termcrbc2)
                            +(2.0*zab*zbc*terma)/(termcrab2*rbc2)
                            -terma/(termcrab2)
                            +(3.0*zab2*terma)/(termcrab2*rab2)                                                                      
                            -(term5*zbc*terma)/(termcrbc2*sine)
                            -(zab*term5*terma)/(termss*rab2)
                            +(3.0*term5*term5*terma)/(terms4)
                            -(term3*terma)/termss2;
                    term4 = (-(2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(yab*xbc-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic))*terma)/(termc2)
                             -(2.0*zbc*terma2)/(termc2*rbc2) - (2.0*zab*terma2)/(termc2*rab2)
                             +(term5*terma2)/(termc*termss))
                           *( (zbd*terma)/(termc*rbd2) + (yab*xbc-xab*ybc+xbd*(ybc-yab)+ybd*(xia-xic))/termc
                             +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term5*terma));
                    dedzibzib = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term3 = (8.0*zbc2*dot2)/(rrab2rbc22*rbc2) - (4.0*zab*zbc*dot)/(rrab2rbc22)
                            +(4.0*ziabc*zbc*dot)/(rrab2rbc22) - (2.0*dot2)/(rrab2rbc22)
                            +(4.0*zab*zbc*dot2)/(rab2rbc2*rab2rbc2) - (2.0*zab*ziabc)/(rab2rbc2)
                            +(2.0*dot)/rab2rbc2 - (4.0*zab2*dot)/(rrab22rbc2);
                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2); 
                    term1 = ( (yabxbd)/termc - (term1zca)/termss2 - (zbc*terma)/(termcrbc2))
                           *( (zbd*terma)/(termc*rbd2) + (yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                              + (zbc*terma)/(termcrbc2) +(zab*terma)/(termcrab2) - (term5*terma)/termss2 );
                    term2 =  (zbd*(yabxbd))/(termc*rbd2)
                            -(zbd*term1zca)/(2.0*termss*rbd2)
                            -(zbc*zbd*terma)/(termcrbc2*rbd2)
                            -(term1zc*(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/termss2
                            +(zbc*(yabxbd))/(termcrbc2)
                            -(zbc*(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termcrbc2)
                            +(zab*(yabxbd))/(termcrab2)
                            -(term5*(yabxbd))/termss2
                            -(zbc*term1zca)/(2.0*termss*rbc2)
                            -(zab*term1zca)/(2.0*termss*rab2)
                            -(term3*terma)/termss2
                            -(3.0*zbc2*terma)/(termcrbc2*rbc2)
                            +terma/(termcrbc2)
                            -(zab*zbc*terma)/(termcrab2*rbc2)
                            +(3.0*term1zc*term5*terma)/(terms4)
                            +(zbc*term5*terma)/(2.0*termss*rbc2);
                    term4 = (-(2.0*(yabxbd)*terma)/(termc2) + (term1zca*terma)/(termc*termss)
                             +(2.0*zbc*terma2)/(termc2*rbc2))
                           *( (zbd*terma)/(termc*rbd2) + (yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                             +(zbc*terma)/(termcrbc2) + (zab*terma)/(termcrab2) - (term5*terma)/termss2);
                    dedzibzic = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term5 = -(2.0*zbc*dot2)/(rrab2rbc22) - (2.0*ziabc*dot)/rab2rbc2 - (2.0*zab*dot2)/(rrab22rbc2); 
                    term1 = ( -(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*( (zbd*terma)/(termc*rbd2)
                              +(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc + (zbc*terma)/(termcrbc2)
                              +(zab*terma)/(termcrab2) - (term5*terma)/termss2 );
                    term2 = -(3.0*zbd2*terma)/(termcrbd22)
                            +(zbd*(xabybc))/(termc*rbd2)
                            -(zbd*(yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic)))/(termc*rbd2)
                            +terma/(termc*rbd2)
                            -(zbc*zbd*terma)/(termcrbc2*rbd2)
                            -(zab*zbd*terma)/(termcrab2*rbd2)
                            +(term5*zbd*terma)/(2.0*termss*rbd2)
                            +(zbc*(xabybc))/(termcrbc2)
                            +(zab*(xabybc))/(termcrab2)
                            -(term5*(xabybc))/termss2;
                    term4 = ((2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                           *((zbd*terma)/(termc*rbd2) + (yab*xbc-xab*ybc+xbd*(yic-yia)+ybd*(xia-xic))/termc
                            +(zbc*terma)/(termcrbc2) +(zab*terma)/(termcrab2) - (term5*terma)/termss2);
                    dedzibzid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = term1aa(zbc, yabxbd, terma, term1zc, termc, termss2, termcrbc2);
                    term3 = (8.0*zab*zbc*dot)/(rrab2rbc22)- (2.0*zab2)/rab2rbc2 - (8.0*zbc2*dot2)/(rbc2*rrab2rbc22)
                             + (2.0*dot2)/(rrab2rbc22);
                    term2 = - (term1zc*(yabxbd))/(termss)
                            - (2.0*zbc*(yabxbd))/(termcrbc2)
                            + (3.0*terma*term1zc*term1zc)/(terms4)
                            - (terma*term3)/termss2
                            + (zbc*term1zca)/(termss*rbc2)
                            + (3.0*zbc2*terma)/(rbc2*rbc2*termc)
                            - terma/(rbc2*termc);
                    term4 = ( - (2.0*terma*(yabxbd))/(termc2)
                              + (term1zca*terma)/(termc*termss)
                              + (2.0*zbc*terma2)/(rbc2*termc2) )
                           *(   (yabxbd)/termc - (term1zca)/termss2
                              - (zbc*terma)/(rbc2*termc));                                   
                    dedziczic = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;

                    term1 = (-(zbd*terma)/(termc*rbd2) + (xabybc)/termc)*( (yabxbd)/termc
                             -(term1zca)/termss2 - (zbc*terma)/(termcrbc2));
                    term2 =  -(zbd*(yabxbd))/(rbd2*termc) + (term1zc*zbd*terma)/(2.0*termss*rbd2)
                             +(zbc*zbd*terma)/(rbc2*rbd2*termc) - (term1zc*(xabybc))/termss2
                             -(zbc*(xabybc))/(termcrbc2);
                    term4 = ( (2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                           *( (yabxbd)/termc - (term1zca)/termss2 - (zbc*terma)/(rbc2*termc));
                    dedziczid = ddtdenom*term1 + ddt1denom*term2 - ddt2denom*term4;
                    
                    term1 =  -(zbd*terma)/(termc*rbd2) + (xabybc)/termc;
                    term2 =   (3.0*zbd2*terma)/(termcrbd22) - (2.0*(xabybc)*zbd)/(termc*rbd2)
                             -(terma)/(termc*rbd2);
                    term4 =   ((2.0*zbd*terma2)/(termc2*rbd2) - (2.0*(xabybc)*terma)/(termc2))
                             *( -(zbd*terma)/(termc*rbd2) + (xabybc)/termc );
                    dedzidzid = ddtdenom*term1*term1 + ddt1denom*term2 - ddt2denom*term4;
                                  
         if (iatom == ia)
         {
            hess.hessx[ia][0] += dedxiaxia ;
            hess.hessx[ia][1] += dedxiayia ;
            hess.hessx[ia][2] += dedxiazia ;
            
            hess.hessy[ia][0] += dedxiayia ;            
            hess.hessy[ia][1] += dedyiayia ;
            hess.hessy[ia][2] += dedyiazia ;
            
            hess.hessz[ia][0] += dedxiazia ;
            hess.hessz[ia][1] += dedyiazia ;            
            hess.hessz[ia][2] += dedziazia ;
                       
            hess.hessx[ib][0] += dedxiaxib ;
            hess.hessx[ib][1] += dedxiayib ;
            hess.hessx[ib][2] += dedxiazib ;
            
            hess.hessy[ib][0] += dedxibyia ;            
            hess.hessy[ib][1] += dedyiayib ;
            hess.hessy[ib][2] += dedyiazib ;
            
            hess.hessz[ib][0] += dedxibzia ;
            hess.hessz[ib][1] += dedyibzia ;            
            hess.hessz[ib][2] += dedziazib ;
                      
            hess.hessx[ic][0] += dedxiaxic ;
            hess.hessy[ic][0] += dedxicyia ;
            hess.hessz[ic][0] += dedxiczia ;
                        
            hess.hessx[ic][1] += dedxiayic ;
            hess.hessy[ic][1] += dedyiayic ;
            hess.hessz[ic][1] += dedyiczia ;
            
            hess.hessx[ic][2] += dedxiazic ;
            hess.hessy[ic][2] += dedyiazic ;
            hess.hessz[ic][2] += dedziazic ;
            
            hess.hessx[id][0] += dedxiaxid ;
            hess.hessy[id][0] += dedxidyia ;
            hess.hessz[id][0] += dedxidzia ;
            hess.hessx[id][1] += dedxiayid ;
            hess.hessy[id][1] += dedyiayid ;
            hess.hessz[id][1] += dedyidzia ;
            hess.hessx[id][2] += dedxiazid ;
            hess.hessy[id][2] += dedyiazid ;
            hess.hessz[id][2] += dedziazid ;
         }else if (iatom == ib)
         {
            hess.hessx[ib][0] += dedxibxib ;
            hess.hessy[ib][0] += dedxibyib ;
            hess.hessz[ib][0] += dedxibzib ;           
            hess.hessx[ib][1] += dedxibyib ;
            hess.hessy[ib][1] += dedyibyib ;
            hess.hessz[ib][1] += dedyibzib ;           
            hess.hessx[ib][2] += dedxibzib ;
            hess.hessy[ib][2] += dedyibzib ;
            hess.hessz[ib][2] += dedzibzib ;
                       
            hess.hessx[ia][0] += dedxiaxib ;
            hess.hessy[ia][0] += dedxiayib ;
            hess.hessz[ia][0] += dedxiazib ;            
            hess.hessx[ia][1] += dedxibyia ;
            hess.hessy[ia][1] += dedyiayib ;
            hess.hessz[ia][1] += dedyiazib ;            
            hess.hessx[ia][2] += dedxibzia ;
            hess.hessy[ia][2] += dedyibzia ;
            hess.hessz[ia][2] += dedziazib ;
                      
            hess.hessx[ic][0] += dedxibxic ;
            hess.hessy[ic][0] += dedxicyib ;
            hess.hessz[ic][0] += dedxiczib ;
                      
            hess.hessx[ic][1] += dedxibyic ;
            hess.hessy[ic][1] += dedyibyic ;
            hess.hessz[ic][1] += dedyiczib ;
                     
            hess.hessx[ic][2] += dedxibzic ;
            hess.hessy[ic][2] += dedyibzic ;
            hess.hessz[ic][2] += dedzibzic ;
                     
            hess.hessx[id][0] += dedxibxid ;
            hess.hessy[id][0] += dedxidyib ;
            hess.hessz[id][0] += dedxidzib ;
                       
            hess.hessx[id][1] += dedxibyid ;
            hess.hessy[id][1] += dedyibyid ;
            hess.hessz[id][1] += dedyidzib ;
                      
            hess.hessx[id][2] += dedxibzid ;
            hess.hessy[id][2] += dedyibzid ;
            hess.hessz[id][2] += dedzibzid ;
         }else if (iatom == ic)
         {
            hess.hessx[ic][0] += dedxicxic ;
            hess.hessy[ic][0] += dedxicyic ;
            hess.hessz[ic][0] += dedxiczic ;
            hess.hessx[ic][1] += dedxicyic ;
            hess.hessy[ic][1] += dedyicyic ;
            hess.hessz[ic][1] += dedyiczic ;
            hess.hessx[ic][2] += dedxiczic ;
            hess.hessy[ic][2] += dedyiczic ;
            hess.hessz[ic][2] += dedziczic ;
            hess.hessx[ia][0] += dedxiaxic ;
            hess.hessy[ia][0] += dedxiayic ;
            hess.hessz[ia][0] += dedxiazic ;
            hess.hessx[ia][1] += dedxicyia ;
            hess.hessy[ia][1] += dedyiayic ;
            hess.hessz[ia][1] += dedyiazic ;
            hess.hessx[ia][2] += dedxiczia ;
            hess.hessy[ia][2] += dedyiczia ;
            hess.hessz[ia][2] += dedziazic ;
            hess.hessx[ib][0] += dedxibxic ;
            hess.hessy[ib][0] += dedxibyic ;
            hess.hessz[ib][0] += dedxibzic ;
            hess.hessx[ib][1] += dedxicyib ;
            hess.hessy[ib][1] += dedyibyic ;
            hess.hessz[ib][1] += dedyibzic ;            
            hess.hessx[ib][2] += dedxiczib ;
            hess.hessy[ib][2] += dedyiczib ;
            hess.hessz[ib][2] += dedzibzic ;            
            hess.hessx[id][0] += dedxicxid ;
            hess.hessy[id][0] += dedxidyic ;
            hess.hessz[id][0] += dedxidzic ;            
            hess.hessx[id][1] += dedxicyid ;
            hess.hessy[id][1] += dedyicyid ;
            hess.hessz[id][1] += dedyidzic ;            
            hess.hessx[id][2] += dedxiczid ;
            hess.hessy[id][2] += dedyiczid ;
            hess.hessz[id][2] += dedziczid ;
         }else if (iatom == id)
         {
            hess.hessx[id][0] += dedxidxid ;
            hess.hessy[id][0] += dedxidyid ;
            hess.hessz[id][0] += dedxidzid ;
            
            hess.hessx[id][1] += dedxidyid ;
            hess.hessy[id][1] += dedyidyid ;
            hess.hessz[id][1] += dedyidzid ;
            
            hess.hessx[id][2] += dedxidzid ;
            hess.hessy[id][2] += dedyidzid ;
            hess.hessz[id][2] += dedzidzid ;
            
            hess.hessx[ia][0] += dedxiaxid ;
            hess.hessy[ia][0] += dedxiayid ;
            hess.hessz[ia][0] += dedxiazid ;
            hess.hessx[ia][1] += dedxidyia ;
            hess.hessy[ia][1] += dedyiayid ;
            hess.hessz[ia][1] += dedyiazid ;
            hess.hessx[ia][2] += dedxidzia ;
            hess.hessy[ia][2] += dedyidzia ;
            hess.hessz[ia][2] += dedziazid ;
            
            hess.hessx[ib][0] += dedxibxid ;
            hess.hessy[ib][0] += dedxibyid ;
            hess.hessz[ib][0] += dedxibzid ;
            hess.hessx[ib][1] += dedxidyib ;
            hess.hessy[ib][1] += dedyibyid ;
            hess.hessz[ib][1] += dedyibzid ;
            hess.hessx[ib][2] += dedxidzib ;
            hess.hessy[ib][2] += dedyidzib ;
            hess.hessz[ib][2] += dedzibzid ;
            
            hess.hessx[ic][0] += dedxicxid ;
            hess.hessy[ic][0] += dedxicyid ;
            hess.hessz[ic][0] += dedxiczid ;
            hess.hessx[ic][1] += dedxidyic ;
            hess.hessy[ic][1] += dedyicyid ;
            hess.hessz[ic][1] += dedyiczid ;
            hess.hessx[ic][2] += dedxidzic ;
            hess.hessy[ic][2] += dedyidzic ;
            hess.hessz[ic][2] += dedziczid ;
         }


                 }
                 
              }
          }
      }
}

