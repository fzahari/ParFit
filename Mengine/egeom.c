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
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "field.h"
#include "fix.h"
#include "units.h"
#include "minim_values.h"

void egeom()
{
   int i, it,kt,itype,ktype;
   double xr, yr, zr, rik, rik2, bcorr;
   double dt, dt2, e;

   energies.estr = 0.0F;

   if (minim_values.iprint)
   {
       fprintf(pcmoutfile,"\nFixed Distance Terms \n");
       fprintf(pcmoutfile,"        At1       At2     R       BLen    Bconst      Eb\n");
   }

// fixed distances
   for (i=0; i < fixdis.nfxstr; i++)
   {
      it = fixdis.ifxstr[i][0];
      kt = fixdis.ifxstr[i][1];
      itype = atom[it].type;
      ktype = atom[kt].type;
      bcorr = 0.0;
      if ( atom[it].use || atom[kt].use )
      {
         xr = atom[it].x - atom[kt].x;
         yr = atom[it].y - atom[kt].y;
         zr = atom[it].z - atom[kt].z;
         rik2 = xr*xr + yr*yr + zr*zr;
         rik = sqrt(rik2);
         dt = rik - fixdis.fxstrd[i];
         
         dt2 = dt*dt;
         e = units.bndunit *fixdis.fxstrc[i]*dt2;

         energies.estr += e;
         atom[it].energy += e;
         atom[kt].energy += e;
         if (minim_values.iprint)
           fprintf(pcmoutfile,"Fixed Dist: %2s(%-3d) - %2s(%-3d) %-8.3f %-8.3f %-8.3f = %-8.4f\n",atom[it].name,it,atom[kt].name
              ,kt, rik, fixdis.fxstrd[i], fixdis.fxstrc[i],e);
      }
   }  
}
void egeom1()
{
/* compute stretching energy and first derivatives */
   int i, it, kt;
   double xr, yr, zr, rik, rik2;
   double dt, dt2, e, deddt;
   double de,dedx,dedy,dedz;
 
   energies.estr = 0.0F;
      for (i=0; i <= natom; i++)
      {
          deriv.deb[i][0] = 0.0;
          deriv.deb[i][1] = 0.0;
          deriv.deb[i][2] = 0.0;
      }
//  fixed distances
   for (i=0; i < fixdis.nfxstr; i++)
   {
      it = fixdis.ifxstr[i][0];
      kt = fixdis.ifxstr[i][1];
      if ( atom[it].use || atom[kt].use )
      {
         xr = atom[it].x - atom[kt].x;
         yr = atom[it].y - atom[kt].y;
         zr = atom[it].z - atom[kt].z;
         rik2 = xr*xr + yr*yr + zr*zr;
         rik = sqrt(rik2);
         dt = rik - fixdis.fxstrd[i];
         
         dt2 = dt*dt;
         e = units.bndunit *fixdis.fxstrc[i]*dt2;
         deddt = 2.0 * units.bndunit * fixdis.fxstrc[i] * dt;
         if (fabs(rik) < 0.1)
            de = 0.0;
         else
            de = deddt/rik;
         dedx = de*xr;
         dedy = de*yr;
         dedz = de*zr;
         energies.estr += e;
         deriv.deb[it][0] += dedx;
         deriv.deb[it][1] += dedy;
         deriv.deb[it][2] += dedz;
         deriv.deb[kt][0] -= dedx;
         deriv.deb[kt][1] -= dedy;
         deriv.deb[kt][2] -= dedz;         

         virial.virx += xr*dedx;
         virial.viry += yr*dedy;
         virial.virz += zr*dedz;
      }
   }
}

void egeom2(int ia)
{
    int j,k,m,ibond;
    double  dt,dt2;
    double xr,yr,zr,rik,rik2,deddt,d2eddt2;
    double de,term,termx,termy,termz,d2e[3][3];

// fixed distances
    for (m=0; m < MAXIAT; m++)
    {
        if (atom[ia].iat[m]!= 0 && atom[ia].bo[m] != 9)
        {
            ibond = find_fixed_bond(ia, atom[ia].iat[m]);
            if (ibond > -1)
            {
              if (fixdis.ifxstr[ibond][0] == ia)
               k = fixdis.ifxstr[ibond][1];
              else
               k = fixdis.ifxstr[ibond][0];

              xr = atom[ia].x - atom[k].x;
              yr = atom[ia].y - atom[k].y;
              zr = atom[ia].z - atom[k].z;
              rik2 = xr*xr + yr*yr + zr*zr;
              rik = sqrt(rik2);
              dt = rik - fixdis.fxstrd[ibond];
              dt2 = dt * dt;

              deddt = 2.0 * units.bndunit * fixdis.fxstrc[ibond] * dt;
              d2eddt2 = 2.0 * units.bndunit * fixdis.fxstrc[ibond];
               
	      if (fabs(rik2) < 0.1)
              {
                de = 0.0;
                term = 0.0;
              }else
              {
                de = deddt / rik;
                term = (d2eddt2-de) / rik2;
              }

              termx = term * xr;
              termy = term * yr;
              termz = term * zr;
              d2e[0][0] = termx*xr + de;
              d2e[1][0] = termx*yr;
              d2e[2][0] = termx*zr;
              d2e[0][1] = d2e[1][0];
              d2e[1][1] = termy*yr + de;
              d2e[2][1] = termy*zr;
              d2e[0][2] = d2e[2][0];
              d2e[1][2] = d2e[2][1];
              d2e[2][2] = termz*zr + de;
    
              for (j=0; j < 3; j++)
              {
                 hess.hessx[ia][j] += d2e[j][0];
                 hess.hessy[ia][j] += d2e[j][1];
                 hess.hessz[ia][j] += d2e[j][2];
                 hess.hessx[k][j] -= d2e[j][0];
                 hess.hessy[k][j] -= d2e[j][1];
                 hess.hessz[k][j] -= d2e[j][2];
              }
            }
        }
    }
}
