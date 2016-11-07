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
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "cutoffs.h"
#include "pot.h"
#include "units.h"
#include "minim_values.h"
#include "dipole.h"

void edipole()
{
    int i, j, ia, ib, ic, id;
    double xi,yi,zi,xk,yk,zk;
    double xq,yq,zq,xr,yr,zr;
    double f,fi,fik;
    double e,ri2,rk2,rirkr3,r2;
    double doti,dotk,dotp;



    energies.eu = 0.0;
    f = 14.388/units.dielec;
    if (minim_values.iprint)
    {
        fprintf(pcmoutfile,"\nBond Dipole Terms\n");
        fprintf(pcmoutfile,"     At1     At2          At3     At4      Bmom-i     Bmom-j     Rik           Edd\n");
    }

    for (i=0; i < bonds_ff.nbnd -1 ; i++)
    {
      if (fabs(bonds_ff.bmom[i]) > 0.01)
        {
          ia = bonds_ff.idpl[i][0];
          ib = bonds_ff.idpl[i][1];
            xi = atom[ib].x - atom[ia].x;
            yi = atom[ib].y - atom[ia].y;
            zi = atom[ib].z - atom[ia].z;
            ri2 = xi*xi + yi*yi + zi*zi;
            xq = atom[ia].x + 0.5*xi;
            yq = atom[ia].y + 0.5*yi;
            zq = atom[ia].z + 0.5*zi;
            fi = f*bonds_ff.bmom[i];
            for (j= i+1; j < bonds_ff.nbnd; j++)
            {
	      if(fabs(bonds_ff.bmom[j]) > 0.01)
                {
                  ic = bonds_ff.idpl[j][0];
                  id = bonds_ff.idpl[j][1];
                  if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use)
                  {
                    if (ic != ia && ic != ib && id != ia && id != ib)
                    {
                      xk = atom[id].x - atom[ic].x;
                      yk = atom[id].y - atom[ic].y;
                      zk = atom[id].z - atom[ic].z;
                      xr = xq - atom[ic].x - 0.5*xk;
                      yr = yq - atom[ic].y - 0.5*yk;
                      zr = zq - atom[ic].z - 0.5*zk;
                      r2 = xr*xr + yr*yr + zr*zr;
                      rk2 = xk*xk + yk*yk + zk*zk;
                      rirkr3 = sqrt(ri2*rk2*r2) * r2;
                      if (rirkr3 < 1.0) rirkr3 = 1.0;
                        dotp = xi*xk + yi*yk + zi*zk;
                        doti = xi*xr + yi*yr + zi*zr;
                        dotk = xk*xr + yk*yr + zk*zr;
                        fik = fi*bonds_ff.bmom[j];
                        e = fik*(dotp-3.00*doti*dotk/r2)/rirkr3;
                        energies.eu += e;
                        atom[ia].energy += e;
                        atom[ib].energy += e;
                        if (minim_values.iprint)
                            fprintf(pcmoutfile,"DD: %2s(%-3d)- %2s(%-3d) -- %2s(%-3d)- %2s(%-3d) :  %-8.3f  %-8.3f   %-8.3f  = %8.4f\n",
                            atom[ia].name, ia,atom[ib].name,ib,atom[ic].name,ic,atom[id].name,id,
                             bonds_ff.bmom[i], bonds_ff.bmom[j], rirkr3, e);
                    }
                  }
                }
            }
        }
    }
}
// =============================
void edipole1()
{
    int i, j, ia, ib, ic, id;
    double xi,yi,zi,xk,yk,zk;
    double xq,yq,zq,xr,yr,zr;
    double f,fi,fik;
    double e,ri2,rk2,rirkr3;
    double doti,dotk,dotp;
    double r2;
    double si1,si2,sk1,sk2;
    double term,dedr,dedrirk;
    double deddoti,deddotk,deddotp;
    double termx,termy,termz;
    double dedrirkri2,dedrirkrk2;
    double termxi,termyi,termzi;
    double termxk,termyk,termzk;
    double dedxi1,dedyi1,dedzi1;
    double dedxi2,dedyi2,dedzi2;
    double dedxk1,dedyk1,dedzk1;
    double dedxk2,dedyk2,dedzk2;



    energies.eu = 0.0;
    f = 14.388/units.dielec;
      for (i=0; i <= natom; i++)
      {
          deriv.deqq[i][0] = 0.0;
          deriv.deqq[i][1] = 0.0;
          deriv.deqq[i][2] = 0.0;
      }

    for (i=0; i < bonds_ff.nbnd -1 ; i++)
    {
      if (fabs(bonds_ff.bmom[i]) > 0.01)
        {
            ia = bonds_ff.idpl[i][0];
            ib = bonds_ff.idpl[i][1];
            si1 = 0.5;
            si2 = 0.5;
            xi = atom[ib].x - atom[ia].x;
            yi = atom[ib].y - atom[ia].y;
            zi = atom[ib].z - atom[ia].z;
            ri2 = xi*xi + yi*yi + zi*zi;
            xq = atom[ia].x + 0.5*xi;
            yq = atom[ia].y + 0.5*yi;
            zq = atom[ia].z + 0.5*zi;
            fi = f*bonds_ff.bmom[i];
            for (j= i+1; j < bonds_ff.nbnd; j++)
            {
	      if(fabs(bonds_ff.bmom[j]) > 0.01)
                {
                  ic = bonds_ff.idpl[j][0];
                  id = bonds_ff.idpl[j][1];
                  if (atom[ia].use || atom[ib].use || atom[ic].use || atom[id].use)
                  {
                    if (ic != ia && ic != ib && id != ia && id != ib)
                    {
                      sk1 = 0.5;
                      sk2 = 0.5;
                      xk = atom[id].x - atom[ic].x;
                      yk = atom[id].y - atom[ic].y;
                      zk = atom[id].z - atom[ic].z;
                      xr = xq - atom[ic].x - 0.5*xk;
                      yr = yq - atom[ic].y - 0.5*yk;
                      zr = zq - atom[ic].z - 0.5*zk;
                      r2 = xr*xr + yr*yr + zr*zr;
                      rk2 = xk*xk + yk*yk + zk*zk;
                      rirkr3 = sqrt(ri2*rk2*r2) * r2;
                      if (rirkr3 < 1.0) rirkr3 = 1.0;
                      dotp = xi*xk + yi*yk + zi*zk;
                      doti = xi*xr + yi*yr + zi*zr;
                      dotk = xk*xr + yk*yr + zk*zr;
                      fik = fi*bonds_ff.bmom[j];
                      e = fik*(dotp-3.00*doti*dotk/r2)/rirkr3;

                       term = -fik / (rirkr3*r2);
                       deddotp = -term * r2;
                       deddoti = term * 3.0*dotk;
                       deddotk = term * 3.0*doti;
                       dedr = term * (3.0*dotp-15.0*doti*dotk/r2);
                       dedrirk = -e;
                       dedrirkri2 = dedrirk / ri2;
                       dedrirkrk2 = dedrirk / rk2;

                       termx = dedr*xr + deddoti*xi + deddotk*xk;
                       termy = dedr*yr + deddoti*yi + deddotk*yk;
                       termz = dedr*zr + deddoti*zi + deddotk*zk;
                       termxi = dedrirkri2*xi + deddotp*xk + deddoti*xr;
                       termyi = dedrirkri2*yi + deddotp*yk + deddoti*yr;
                       termzi = dedrirkri2*zi + deddotp*zk + deddoti*zr;
                       termxk = dedrirkrk2*xk + deddotp*xi + deddotk*xr;
                       termyk = dedrirkrk2*yk + deddotp*yi + deddotk*yr;
                       termzk = dedrirkrk2*zk + deddotp*zi + deddotk*zr;

                       dedxi1 = si1*termx - termxi;
                       dedyi1 = si1*termy - termyi;
                       dedzi1 = si1*termz - termzi;
                       dedxi2 = si2*termx + termxi;
                       dedyi2 = si2*termy + termyi;
                       dedzi2 = si2*termz + termzi;
                       dedxk1 = -sk1*termx - termxk;
                       dedyk1 = -sk1*termy - termyk;
                       dedzk1 = -sk1*termz - termzk;
                       dedxk2 = -sk2*termx + termxk;
                       dedyk2 = -sk2*termy + termyk;
                       dedzk2 = -sk2*termz + termzk;

                      energies.eu += e;
                      deriv.deqq[ia][0] += dedxi1;
                      deriv.deqq[ia][1] += dedyi1;
                      deriv.deqq[ia][2] += dedzi1;

                      deriv.deqq[ib][0] += dedxi2;
                      deriv.deqq[ib][1] += dedyi2;
                      deriv.deqq[ib][2] += dedzi2;

                      deriv.deqq[ic][0] += dedxk1;
                      deriv.deqq[ic][1] += dedyk1;
                      deriv.deqq[ic][2] += dedzk1;

                      deriv.deqq[id][0] += dedxk2;
                      deriv.deqq[id][1] += dedyk2;
                      deriv.deqq[id][2] += dedzk2;

                      virial.virx += xr*(dedxi1+dedxi2);
                      virial.viry += yr*(dedyi1+dedyi2);
                      virial.virz += zr*(dedzi1+dedzi2);
                    }
                  }
                }
            }
        }
    }
}
// =====================
void edipole2(int i)
{
    int ia,ib, ic, id;
    int idipole,kdipole;
    double f,fi,fik;
    double xi,yi,zi,xk,yk,zk;
    double xq,yq,zq,xr,yr,zr;
    double r2,ri2,rk2,rirkr3;
    double doti,dotk,dotp;
    double si1,si2,sk1,sk2;
    double term,dedr,dedrirk;
    double deddoti,deddotk,deddotp;
    double termx,termy,termz;
    double termxi,termyi,termzi;
    double termxk,termyk,termzk;
    double anum,r2inv,ri2inv,dotik,xrr2,yrr2,zrr2;
    double xiri2,yiri2,ziri2;
    double xkrk2,ykrk2,zkrk2;
    double xixr,xiyr,xizr,yixr,yiyr,yizr,zixr,ziyr,zizr;
    double xkxr,xkyr,xkzr,ykxr,ykyr,ykzr,zkxr,zkyr,zkzr;
    double xixk,xiyk,xizk,yixk,yiyk,yizk,zixk,ziyk,zizk;
    double xrxr,xryr,xrzr,yryr,yrzr,zrzr;
    double xidotk,yidotk,zidotk,xkdoti,ykdoti,zkdoti;
    double factor,factori,factork,part,partik;
    double dtdxi1,dtdyi1,dtdzi1,dtdxi2,dtdyi2,dtdzi2;
    double dtxdxi1,dtxidxi1,dtxkdxi1,dtxdxi2,dtxidxi2,dtxkdxi2;
    double dtydxi1,dtyidxi1,dtykdxi1,dtydxi2,dtyidxi2,dtykdxi2;
    double dtzdxi1,dtzidxi1,dtzkdxi1,dtzdxi2,dtzidxi2,dtzkdxi2;
    double dtxdyi1,dtxidyi1,dtxkdyi1,dtxdyi2,dtxidyi2,dtxkdyi2;
    double dtydyi1,dtyidyi1,dtykdyi1,dtydyi2,dtyidyi2,dtykdyi2;
    double dtzdyi1,dtzidyi1,dtzkdyi1,dtzdyi2,dtzidyi2,dtzkdyi2;
    double dtxdzi1,dtxidzi1,dtxkdzi1,dtxdzi2,dtxidzi2,dtxkdzi2;
    double dtydzi1,dtyidzi1,dtykdzi1,dtydzi2,dtyidzi2,dtykdzi2;
    double dtzdzi1,dtzidzi1,dtzkdzi1,dtzdzi2,dtzidzi2,dtzkdzi2;
    
    f = 14.388/units.dielec;
    dtdxi1 = dtdyi1 = dtdzi1 = dtdxi2 = dtdyi2 = dtdzi2 = 0.0;
    dtxdxi1 = dtxidxi1 = dtxkdxi1 = dtxdxi2 = dtxidxi2 = dtxkdxi2 = 0.0;
    dtydxi1 = dtyidxi1 = dtykdxi1 = dtydxi2 = dtyidxi2 = dtykdxi2 = 0.0;
    dtzdxi1 = dtzidxi1 = dtzkdxi1 = dtzdxi2 = dtzidxi2 = dtzkdxi2 = 0.0;
    dtxdyi1 = dtxidyi1 = dtxkdyi1 = dtxdyi2 = dtxidyi2 = dtxkdyi2 = 0.0;
    dtydyi1 = dtyidyi1 = dtykdyi1 = dtydyi2 = dtyidyi2 = dtykdyi2 = 0.0;
    dtzdyi1 = dtzidyi1 = dtzkdyi1 = dtzdyi2 = dtzidyi2 = dtzkdyi2 = 0.0;
    dtxdzi1 = dtxidzi1 = dtxkdzi1 = dtxdzi2 = dtxidzi2 = dtxkdzi2 = 0.0;
    dtydzi1 = dtyidzi1 = dtykdzi1 = dtydzi2 = dtyidzi2 = dtykdzi2 = 0.0;
    dtzdzi1 = dtzidzi1 = dtzkdzi1 = dtzdzi2 = dtzidzi2 = dtzkdzi2 = 0.0;

    for (idipole = 0; idipole < bonds_ff.nbnd; idipole++)
    {
        ia = bonds_ff.idpl[idipole][0];
        ib = bonds_ff.idpl[idipole][1];
        if ( (ia == i || ib == i) && bonds_ff.bmom[idipole] != 0.0)
        {
            si1 = 0.5;
            si2 = 0.5;
            xi = atom[ib].x - atom[ia].x;
            yi = atom[ib].y - atom[ia].y;
            zi = atom[ib].z - atom[ia].z;
            ri2 = xi*xi + yi*yi + zi*zi;
            xq = atom[ia].x + xi*si2;
            yq = atom[ia].y + yi*si2;
            zq = atom[ia].z + zi*si2;
            fi = f*bonds_ff.bmom[idipole];
            for (kdipole = 0; kdipole < bonds_ff.nbnd; kdipole++)
            {
                ic = bonds_ff.idpl[kdipole][0];
                id = bonds_ff.idpl[kdipole][1];
                if (ic != ia && ic != ib && id != ia && id != ib && fabs(bonds_ff.bmom[kdipole]) >  0.01 )
                {
                   sk1 = 0.5;
                   sk2 = 0.5;
                   xk = atom[id].x - atom[ic].x;
                   yk = atom[id].y - atom[ic].y;
                   zk = atom[id].z - atom[ic].z;
                   xr = xq - atom[ic].x - xk*sk2; 
                   yr = yq - atom[ic].y - yk*sk2; 
                   zr = zq - atom[ic].z - zk*sk2; 
                   r2 = xr*xr + yr*yr + zr*zr;
                   rk2 = xk*xk + yk*yk + zk*zk;
                   rirkr3 = sqrt(ri2*rk2*r2) * r2;
                   if (rirkr3 < 1.0) rirkr3 = 1.0;
                   dotp = xi*xk + yi*yk + zi*zk;
                   doti = xi*xr + yi*yr + zi*zr;
                   dotk = xk*xr + yk*yr + zk*zr;
                   fik = fi * bonds_ff.bmom[kdipole];
                  dotik = doti * dotk;
                  anum = dotp*r2 - 3.0*dotik;
                  r2inv = 15.0 / r2;
                  ri2inv = 1.0 / ri2;
                  xrr2 = xr / r2;
                  yrr2 = yr / r2;
                  zrr2 = zr / r2;
                  xiri2 = xi / ri2;
                  yiri2 = yi / ri2;
                  ziri2 = zi / ri2;
                  xkrk2 = xk / rk2;
                  ykrk2 = yk / rk2;
                  zkrk2 = zk / rk2;
                  xixr = xi * xr;
                  xiyr = xi * yr;
                  xizr = xi * zr;
                  yixr = yi * xr;
                  yiyr = yi * yr;
                  yizr = yi * zr;
                  zixr = zi * xr;
                  ziyr = zi * yr;
                  zizr = zi * zr;
                  xkxr = xk * xr;
                  xkyr = xk * yr;
                  xkzr = xk * zr;
                  ykxr = yk * xr;
                  ykyr = yk * yr;
                  ykzr = yk * zr;
                  zkxr = zk * xr;
                  zkyr = zk * yr;
                  zkzr = zk * zr;
                  xixk = xi * xk;
                  xiyk = xi * yk;
                  xizk = xi * zk;
                  yixk = yi * xk;
                  yiyk = yi * yk;
                  yizk = yi * zk;
                  zixk = zi * xk;
                  ziyk = zi * yk;
                  zizk = zi * zk;
                  xrxr = 3.0 * xr * xr;
                  xryr = 3.0 * xr * yr;
                  xrzr = 3.0 * xr * zr;
                  yryr = 3.0 * yr * yr;
                  yrzr = 3.0 * yr * zr;
                  zrzr = 3.0 * zr * zr;
                  xidotk = xi * dotk;
                  yidotk = yi * dotk;
                  zidotk = zi * dotk;
                  xkdoti = xk * doti;
                  ykdoti = yk * doti;
                  zkdoti = zk * doti;
                  
                  term = -fik / (rirkr3*r2);
                  deddotp = -term * r2;
                  deddoti = term * 3.0*dotk;
                  deddotk = term * 3.0*doti;
                  dedr = term * (3.0*dotp-15.0*dotik/r2);
                  dedrirk = term * anum;

                  termx = dedr*xr + deddoti*xi + deddotk*xk;
                  termy = dedr*yr + deddoti*yi + deddotk*yk;
                  termz = dedr*zr + deddoti*zi + deddotk*zk;
                  termxi = dedrirk*xiri2 + deddotp*xk + deddoti*xr;
                  termyi = dedrirk*yiri2 + deddotp*yk + deddoti*yr;
                  termzi = dedrirk*ziri2 + deddotp*zk + deddoti*zr;
                  termxk = dedrirk*xkrk2 + deddotp*xi + deddotk*xr;
                  termyk = dedrirk*ykrk2 + deddotp*yi + deddotk*yr;
                  termzk = dedrirk*zkrk2 + deddotp*zi + deddotk*zr;

                  if (i == ia)
                  {
                     dtdxi1 = -5.0*si1*xrr2 + xiri2;
                     part = si1*xkdoti - dotk*xr + si1*xidotk - 2.0*si1*dotik*xrr2;
                     partik = -xk*r2 + 2.0*si1*dotp*xr - 3.0*si1*xkdoti + 3.0*xr*dotk- 3.0*si1*xidotk;
                     factor = 3.0*si1*dotp - 6.0*xkxr + 6.0*si1*xixk - 3.0*dotk - r2inv*(xr*part+si1*dotik);
                     factori = 3.0*si1*dotk + si1*xkxr + xiri2*partik - anum*(ri2inv-2.0*xiri2*xiri2);
                     factork = r2 + 3.0*si1*doti + si1*xixr - xrxr + xkrk2*partik;
                     dtxdxi1 = dtdxi1*termx + term*factor;
                     dtxidxi1 = dtdxi1*termxi + term*factori;
                     dtxkdxi1 = dtdxi1*termxk + term*factork;
                     factor = -3.0*xkyr - 3.0*ykxr + 3.0*si1*xiyk + 3.0*si1*yixk - r2inv*yr*part;
                     factori = -2.0*si1*ykxr + 3.0*si1*xkyr  + yiri2*partik + 2.0*anum*yiri2*xiri2;
                     factork = -2.0*si1*yixr - xryr + 3.0*si1*xiyr + ykrk2*partik;
                     dtydxi1 = dtdxi1*termy + term*factor;
                     dtyidxi1 = dtdxi1*termyi + term*factori;
                     dtykdxi1 = dtdxi1*termyk + term*factork;
                     factor = -3.0*xkzr - 3.0*zkxr + 3.0*si1*xizk + 3.0*si1*zixk - r2inv*zr*part;
                     factori = -2.0*si1*zkxr + 3.0*si1*xkzr + ziri2*partik + 2.0*anum*ziri2*xiri2;
                     factork = -2.0*si1*zixr - xrzr + 3.0*si1*xizr + zkrk2*partik;
                     dtzdxi1 = dtdxi1*termz + term*factor;
                     dtzidxi1 = dtdxi1*termzi + term*factori;
                     dtzkdxi1 = dtdxi1*termzk + term*factork;

                     dtdyi1 = -5.0*si1*yrr2 + yiri2;
                     part = si1*ykdoti - dotk*yr + si1*yidotk - 2.0*si1*dotik*yrr2;
                     partik = -yk*r2 + 2.0*si1*dotp*yr - 3.0*si1*ykdoti + 3.0*yr*dotk - 3.0*si1*yidotk;
                     factor = -3.0*ykxr - 3.0*xkyr + 3.0*si1*yixk + 3.0*si1*xiyk - r2inv*xr*part;
                     factori = -2.0*si1*xkyr + 3.0*si1*ykxr + xiri2*partik + 2.0*anum*xiri2*yiri2;
                     factork = -2.0*si1*xiyr - xryr + 3.0*si1*yixr + xkrk2*partik;
                     dtxdyi1 = dtdyi1*termx + term*factor;
                     dtxidyi1 = dtdyi1*termxi + term*factori;
                     dtxkdyi1 = dtdyi1*termxk + term*factork;
                     factor = 3.0*si1*dotp - 6.0*ykyr + 6.0*si1*yiyk - 3.0*dotk - r2inv*(yr*part+si1*dotik);
                     factori = 3.0*si1*dotk + si1*ykyr + yiri2*partik - anum*(ri2inv-2.0*yiri2*yiri2);
                     factork = r2 + 3.0*si1*doti + si1*yiyr - yryr + ykrk2*partik;
                     dtydyi1 = dtdyi1*termy + term*factor;
                     dtyidyi1 = dtdyi1*termyi + term*factori;
                     dtykdyi1 = dtdyi1*termyk + term*factork;
                     factor = -3.0*ykzr - 3.0*zkyr + 3.0*si1*yizk + 3.0*si1*ziyk - r2inv*zr*part;
                     factori = -2.0*si1*zkyr + 3.0*si1*ykzr + ziri2*partik + 2.0*anum*ziri2*yiri2;
                     factork = -2.0*si1*ziyr - yrzr + 3.0*si1*yizr + zkrk2*partik;
                     dtzdyi1 = dtdyi1*termz + term*factor;
                     dtzidyi1 = dtdyi1*termzi + term*factori;
                     dtzkdyi1 = dtdyi1*termzk + term*factork;

                     dtdzi1 = -5.0*si1*zrr2 + ziri2;
                     part = si1*zkdoti - dotk*zr + si1*zidotk - 2.0*si1*dotik*zrr2;
                     partik = -zk*r2 + 2.0*si1*dotp*zr - 3.0*si1*zkdoti + 3.0*zr*dotk  - 3.0*si1*zidotk;
                     factor = -3.0*zkxr - 3.0*xkzr + 3.0*si1*zixk  + 3.0*si1*xizk - r2inv*xr*part;
                     factori = -2.0*si1*xkzr + 3.0*si1*zkxr  + xiri2*partik + 2.0*anum*xiri2*ziri2;
                     factork = -2.0*si1*xizr - xrzr + 3.0*si1*zixr  + xkrk2*partik;
                     dtxdzi1 = dtdzi1*termx + term*factor;
                     dtxidzi1 = dtdzi1*termxi + term*factori;
                     dtxkdzi1 = dtdzi1*termxk + term*factork;
                     factor = -3.0*zkyr - 3.0*ykzr + 3.0*si1*ziyk+ 3.0*si1*yizk - r2inv*yr*part;
                     factori = -2.0*si1*ykzr + 3.0*si1*zkyr + yiri2*partik + 2.0*anum*yiri2*ziri2;
                     factork = -2.0*si1*yizr - yrzr + 3.0*si1*ziyr + ykrk2*partik;
                     dtydzi1 = dtdzi1*termy + term*factor;
                     dtyidzi1 = dtdzi1*termyi + term*factori;
                     dtykdzi1 = dtdzi1*termyk + term*factork;
                     factor = 3.0*si1*dotp - 6.0*zkzr + 6.0*si1*zizk - 3.0*dotk - r2inv*(zr*part+si1*dotik);
                     factori = 3.0*si1*dotk + si1*zkzr + ziri2*partik - anum*(ri2inv-2.0*ziri2*ziri2);
                     factork = r2 + 3.0*si1*doti + si1*zizr - zrzr + zkrk2*partik;
                     dtzdzi1 = dtdzi1*termz + term*factor;
                     dtzidzi1 = dtdzi1*termzi + term*factori;
                     dtzkdzi1 = dtdzi1*termzk + term*factork;

                  } else if (i == ib)
                  {
                     dtdxi2 = -5.0*si2*xrr2 - xiri2;
                     part = si2*xkdoti + dotk*xr + si2*xidotk - 2.0*si2*dotik*xrr2;
                     partik = xk*r2 + 2.0*si2*dotp*xr - 3.0*si2*xkdoti - 3.0*xr*dotk - 3.0*si2*xidotk;
                     factor = 3.0*si2*dotp + 6.0*xkxr+ 6.0*si2*xixk + 3.0*dotk - r2inv*(xr*part+si2*dotik);
                     factori = 3.0*si2*dotk + si2*xkxr + xiri2*partik + anum*(ri2inv-2.0*xiri2*xiri2);
                     factork = -r2 + 3.0*si2*doti + si2*xixr + xrxr + xkrk2*partik;
                     dtxdxi2 = dtdxi2*termx + term*factor;
                     dtxidxi2 = dtdxi2*termxi + term*factori;
                     dtxkdxi2 = dtdxi2*termxk + term*factork;
                     factor = 3.0*xkyr + 3.0*ykxr + 3.0*si2*xiyk + 3.0*si2*yixk - r2inv*yr*part;
                     factori = -2.0*si2*ykxr + 3.0*si2*xkyr + yiri2*partik - 2.0*anum*yiri2*xiri2;
                     factork = -2.0*si2*yixr + xryr + 3.0*si2*xiyr  + ykrk2*partik;
                     dtydxi2 = dtdxi2*termy + term*factor;
                     dtyidxi2 = dtdxi2*termyi + term*factori;
                     dtykdxi2 = dtdxi2*termyk + term*factork;
                     factor = 3.0*xkzr + 3.0*zkxr + 3.0*si2*xizk + 3.0*si2*zixk - r2inv*zr*part;
                     factori = -2.0*si2*zkxr + 3.0*si2*xkzr  + ziri2*partik - 2.0*anum*ziri2*xiri2;
                     factork = -2.0*si2*zixr + xrzr + 3.0*si2*xizr + zkrk2*partik;
                     dtzdxi2 = dtdxi2*termz + term*factor;
                     dtzidxi2 = dtdxi2*termzi + term*factori;
                     dtzkdxi2 = dtdxi2*termzk + term*factork;

                     dtdyi2 = -5.0*si2*yrr2 - yiri2;
                     part = si2*ykdoti + dotk*yr + si2*yidotk - 2.0*si2*dotik*yrr2;
                     partik = yk*r2 + 2.0*si2*dotp*yr - 3.0*si2*ykdoti - 3.0*yr*dotk - 3.0*si2*yidotk;
                     factor = 3.0*ykxr + 3.0*xkyr + 3.0*si2*yixk + 3.0*si2*xiyk - r2inv*xr*part;
                     factori = -2.0*si2*xkyr + 3.0*si2*ykxr + xiri2*partik - 2.0*anum*xiri2*yiri2;
                     factork = -2.0*si2*xiyr + xryr + 3.0*si2*yixr + xkrk2*partik;
                     dtxdyi2 = dtdyi2*termx + term*factor;
                     dtxidyi2 = dtdyi2*termxi + term*factori;
                     dtxkdyi2 = dtdyi2*termxk + term*factork;
                     factor = 3.0*si2*dotp + 6.0*ykyr + 6.0*si2*yiyk + 3.0*dotk - r2inv*(yr*part+si2*dotik);
                     factori = 3.0*si2*dotk + si2*ykyr + yiri2*partik + anum*(ri2inv-2.0*yiri2*yiri2);
                     factork = -r2 + 3.0*si2*doti + si2*yiyr + yryr + ykrk2*partik;
                     dtydyi2 = dtdyi2*termy + term*factor;
                     dtyidyi2 = dtdyi2*termyi + term*factori;
                     dtykdyi2 = dtdyi2*termyk + term*factork;
                     factor = 3.0*ykzr + 3.0*zkyr + 3.0*si2*yizk + 3.0*si2*ziyk - r2inv*zr*part;
                     factori = -2.0*si2*zkyr + 3.0*si2*ykzr + ziri2*partik - 2.0*anum*ziri2*yiri2;
                     factork = -2.0*si2*ziyr + yrzr + 3.0*si2*yizr  + zkrk2*partik;
                     dtzdyi2 = dtdyi2*termz + term*factor;
                     dtzidyi2 = dtdyi2*termzi + term*factori;
                     dtzkdyi2 = dtdyi2*termzk + term*factork;

                     dtdzi2 = -5.0*si2*zrr2 - ziri2;
                     part = si2*zkdoti + dotk*zr + si2*zidotk - 2.0*si2*dotik*zrr2;
                     partik = zk*r2 + 2.0*si2*dotp*zr - 3.0*si2*zkdoti - 3.0*zr*dotk - 3.0*si2*zidotk;
                     factor = 3.0*zkxr + 3.0*xkzr + 3.0*si2*zixk + 3.0*si2*xizk - r2inv*xr*part;
                     factori = -2.0*si2*xkzr + 3.0*si2*zkxr + xiri2*partik - 2.0*anum*xiri2*ziri2;
                     factork = -2.0*si2*xizr + xrzr + 3.0*si2*zixr + xkrk2*partik;
                     dtxdzi2 = dtdzi2*termx + term*factor;
                     dtxidzi2 = dtdzi2*termxi + term*factori;
                     dtxkdzi2 = dtdzi2*termxk + term*factork;
                     factor = 3.0*zkyr + 3.0*ykzr + 3.0*si2*ziyk + 3.0*si2*yizk - r2inv*yr*part;
                     factori = -2.0*si2*ykzr + 3.0*si2*zkyr + yiri2*partik - 2.0*anum*yiri2*ziri2;
                     factork = -2.0*si2*yizr + yrzr + 3.0*si2*ziyr  + ykrk2*partik;
                     dtydzi2 = dtdzi2*termy + term*factor;
                     dtyidzi2 = dtdzi2*termyi + term*factori;
                     dtykdzi2 = dtdzi2*termyk + term*factork;
                     factor = 3.0*si2*dotp + 6.0*zkzr + 6.0*si2*zizk + 3.0*dotk  - r2inv*(zr*part+si2*dotik);
                     factori = 3.0*si2*dotk + si2*zkzr + ziri2*partik + anum*(ri2inv-2.0*ziri2*ziri2);
                     factork = -r2 + 3.0*si2*doti + si2*zizr + zrzr + zkrk2*partik;
                     dtzdzi2 = dtdzi2*termz + term*factor;
                     dtzidzi2 = dtdzi2*termzi + term*factori;
                     dtzkdzi2 = dtdzi2*termzk + term*factork;
                  }

                  if (i == ia)
                  {
                     hess.hessx[ia][0] += si1*dtxdxi1 - dtxidxi1;
                     hess.hessx[ia][1] += si1*dtydxi1 - dtyidxi1;
                     hess.hessx[ia][2] += si1*dtzdxi1 - dtzidxi1;
                     hess.hessx[ib][0] += si2*dtxdxi1 + dtxidxi1;
                     hess.hessx[ib][1] += si2*dtydxi1 + dtyidxi1;
                     hess.hessx[ib][2] += si2*dtzdxi1 + dtzidxi1;
                     hess.hessx[ic][0] += -sk1*dtxdxi1 - dtxkdxi1;
                     hess.hessx[ic][1] += -sk1*dtydxi1 - dtykdxi1;
                     hess.hessx[ic][2] += -sk1*dtzdxi1 - dtzkdxi1;
                     hess.hessx[id][0] += -sk2*dtxdxi1 + dtxkdxi1;
                     hess.hessx[id][1] += -sk2*dtydxi1 + dtykdxi1;
                     hess.hessx[id][2] += -sk2*dtzdxi1 + dtzkdxi1;
                     hess.hessy[ia][0] += si1*dtxdyi1 - dtxidyi1;
                     hess.hessy[ia][1] += si1*dtydyi1 - dtyidyi1;
                     hess.hessy[ia][2] += si1*dtzdyi1 - dtzidyi1;
                     hess.hessy[ib][0] += si2*dtxdyi1 + dtxidyi1;
                     hess.hessy[ib][2] += si2*dtzdyi1 + dtzidyi1;
                     hess.hessy[ib][1] += si2*dtydyi1 + dtyidyi1;
                     hess.hessy[ic][0] += -sk1*dtxdyi1 - dtxkdyi1;
                     hess.hessy[ic][1] += -sk1*dtydyi1 - dtykdyi1;
                     hess.hessy[ic][2] += -sk1*dtzdyi1 - dtzkdyi1;
                     hess.hessy[id][0] += -sk2*dtxdyi1 + dtxkdyi1;
                     hess.hessy[id][1] += -sk2*dtydyi1 + dtykdyi1;
                     hess.hessy[id][2] += -sk2*dtzdyi1 + dtzkdyi1;
                     hess.hessz[ia][0] += si1*dtxdzi1 - dtxidzi1;
                     hess.hessz[ia][1] += si1*dtydzi1 - dtyidzi1;
                     hess.hessz[ia][2] += si1*dtzdzi1 - dtzidzi1;
                     hess.hessz[ib][0] += si2*dtxdzi1 + dtxidzi1;
                     hess.hessz[ib][1] += si2*dtydzi1 + dtyidzi1;
                     hess.hessz[ib][2] += si2*dtzdzi1 + dtzidzi1;
                     hess.hessz[ic][0] += -sk1*dtxdzi1 - dtxkdzi1;
                     hess.hessz[ic][1] += -sk1*dtydzi1 - dtykdzi1;
                     hess.hessz[ic][2] += -sk1*dtzdzi1 - dtzkdzi1;
                     hess.hessz[id][0] += -sk2*dtxdzi1 + dtxkdzi1;
                     hess.hessz[id][1] += -sk2*dtydzi1 + dtykdzi1;
                     hess.hessz[id][2] += -sk2*dtzdzi1 + dtzkdzi1;
                  } else if (i == ib)
                  {
                     hess.hessx[ia][0] += si1*dtxdxi2 - dtxidxi2;
                     hess.hessx[ia][1] += si1*dtydxi2 - dtyidxi2;
                     hess.hessx[ia][2] += si1*dtzdxi2 - dtzidxi2;
                     hess.hessx[ib][0] += si2*dtxdxi2 + dtxidxi2;
                     hess.hessx[ib][1] += si2*dtydxi2 + dtyidxi2;
                     hess.hessx[ib][2] += si2*dtzdxi2 + dtzidxi2;
                     hess.hessx[ic][0] += -sk1*dtxdxi2 - dtxkdxi2;
                     hess.hessx[ic][1] += -sk1*dtydxi2 - dtykdxi2;
                     hess.hessx[ic][2] += -sk1*dtzdxi2 - dtzkdxi2;
                     hess.hessx[id][0] += -sk2*dtxdxi2 + dtxkdxi2;
                     hess.hessx[id][1] += -sk2*dtydxi2 + dtykdxi2;
                     hess.hessx[id][2] += -sk2*dtzdxi2 + dtzkdxi2;
                     hess.hessy[ia][0] += si1*dtxdyi2 - dtxidyi2;
                     hess.hessy[ia][1] += si1*dtydyi2 - dtyidyi2;
                     hess.hessy[ia][2] += si1*dtzdyi2 - dtzidyi2;
                     hess.hessy[ib][0] += si2*dtxdyi2 + dtxidyi2;
                     hess.hessy[ib][1] += si2*dtydyi2 + dtyidyi2;
                     hess.hessy[ib][2] += si2*dtzdyi2 + dtzidyi2;
                     hess.hessy[ic][0] += -sk1*dtxdyi2 - dtxkdyi2;
                     hess.hessy[ic][1] += -sk1*dtydyi2 - dtykdyi2;
                     hess.hessy[ic][2] += -sk1*dtzdyi2 - dtzkdyi2;
                     hess.hessy[id][0] += -sk2*dtxdyi2 + dtxkdyi2;
                     hess.hessy[id][1] += -sk2*dtydyi2 + dtykdyi2;
                     hess.hessy[id][2] += -sk2*dtzdyi2 + dtzkdyi2;
                     hess.hessz[ia][0] += si1*dtxdzi2 - dtxidzi2;
                     hess.hessz[ia][1] += si1*dtydzi2 - dtyidzi2;
                     hess.hessz[ia][2] += si1*dtzdzi2 - dtzidzi2;
                     hess.hessz[ib][0] += si2*dtxdzi2 + dtxidzi2;
                     hess.hessz[ib][1] += si2*dtydzi2 + dtyidzi2;
                     hess.hessz[ib][2] += si2*dtzdzi2 + dtzidzi2;
                     hess.hessz[ic][0] += -sk1*dtxdzi2 - dtxkdzi2;
                     hess.hessz[ic][1] += -sk1*dtydzi2 - dtykdzi2;
                     hess.hessz[ic][2] += -sk1*dtzdzi2 - dtzkdzi2;
                     hess.hessz[id][0] += -sk2*dtxdzi2 + dtxkdzi2;
                     hess.hessz[id][1] += -sk2*dtydzi2 + dtykdzi2;
                     hess.hessz[id][2] += -sk2*dtzdzi2 + dtzkdzi2;
                  }
                }
            }
        }
    }

}
