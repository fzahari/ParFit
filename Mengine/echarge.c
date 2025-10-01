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
#include "field.h"
#include "cutoffs.h"
#include "pot.h"
#include "attached.h"
#include "units.h"
#include "minim_values.h"
#include "skip.h"
#include "gastype.h"
#include "echarge.h"
       
void echarge()
{
   int i, j, k;
   double xi, yi, zi, xk,yk,zk;
   double xr, yr, zr;
   double e, rik, rik2,cc;
   double rdielc;
   double cutoff;

   cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
   energies.eu = 0.0;
   cc = 4.80298*4.80298*14.39418;
   rdielc = 1.0/units.dielec;

   if (minim_values.iprint)
   {
       fprintf(pcmoutfile,"\nCharge-Charge Interactions : Diele = %7.3f\n",units.dielec);
       fprintf(pcmoutfile,"      At1      At2      q1        q2         Rik          Eqq\n");
   }
   for (i=1; i <= natom; i++)
      skip[i] = 0;

   for (i=1; i < natom; i++)
   {
      if (atom[i].charge != 0.0 )
      {
         for (k=0; k < MAXIAT; k++)
         {
            if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
              skip[atom[i].iat[k]] = i;
         }
         for(k=0; k < attached.n13[i]; k++)
            skip[attached.i13[k][i]] = i;
         for(k=0; k < attached.n14[i]; k++)
            skip[attached.i14[k][i]] = -i;
         
         for(j=i+1; j <= natom; j++)
         {
           if (atom[i].use || atom[j].use)
           {
            if ( skip[j] != i)
            {
             if (atom[j].charge != 0.0)
             {
                xi = atom[i].x;
                yi = atom[i].y;
                zi = atom[i].z;
                  
                xk = atom[j].x;
                yk = atom[j].y;
                zk = atom[j].z;
                
                xr = xi - xk;
                yr = yi - yk;
                zr = zi - zk;
                rik2 =  xr*xr + yr*yr + zr*zr;
                if (rik2 < cutoff)
                {
                   rik = sqrt(rik2);
                   e = cc*rdielc*atom[i].charge*atom[j].charge/rik;
                   if (skip[j] == -i)
                     e /= units.chgscale;  
                   energies.eu += e;
                   atom[i].energy += e;
                   atom[j].energy += e;
                 
                  if (minim_values.iprint)
                     fprintf(pcmoutfile,"QQ: %2s(%-3d)- %2s(%-3d)   %-8.3f  %-8.3f   %-8.3f = %8.4f\n",
                      atom[i].name, i, atom[j].name, j,atom[i].charge, atom[j].charge,rik,e);
                }
             }
            }
           }
         }
      }
   }
}
//  ==================
void echarge1()
{
   int i,j, k;
   double xi, yi, zi, fik;
   double xr, yr, zr;
   double xk,yk,zk,rdn;
   double e, rik, cc, rik2;
   double de, dedx, dedy, dedz;
   double rdielc;
   double cutoff;

   cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
   cc = 4.80298*4.80298*14.39418;
   rdielc = 1.0/units.dielec;
   rdn = 1.0;
   if (field.type == MMX || field.type == MM2 )
     rdn = .915;
   else if (field.type == MM3)
     rdn = .923;

   for (i=1; i <= natom; i++)
      skip[i] = 0;

   for (i=1; i < natom; i++)
   {
       if (atom[i].charge != 0.0 )
       {
         for (k=0; k < MAXIAT; k++)
         {
            if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
              skip[atom[i].iat[k]] = i;
         }
         for(k=0; k < attached.n13[i]; k++)
            skip[attached.i13[k][i]] = i;
         for(k=0; k < attached.n14[i]; k++)
            skip[attached.i14[k][i]] = -i;

           for(j=i+1; j <= natom; j++)
           {
             if (atom[i].use || atom[j].use)
             {
               if (atom[j].charge != 0.0 &&  skip[j] != i  )
               {
                  xi = atom[i].x;
                  yi = atom[i].y;
                  zi = atom[i].z;

                  xk = atom[j].x;
                  yk = atom[j].y;
                  zk = atom[j].z;

                xr = xi - xk;
                yr = yi - yk;
                zr = zi - zk;
                rik2 = xr*xr +yr*yr + zr*zr;
                if (rik2 < cutoff)
                {
                   rik = sqrt( rik2);
                   fik = cc*rdielc*atom[i].charge*atom[j].charge;
                   if (skip[j] == -i)
                     fik /= units.chgscale;
                   
                   e = fik/rik;

                   de = -fik/rik2;

                   de = de/rik;
                   dedx = de*xr;
                   dedy = de*yr;
                   dedz = de*zr;
                  
                   deriv.deqq[i][0] += dedx;
                   deriv.deqq[i][1] += dedy;
                   deriv.deqq[i][2] += dedz;
         
                   deriv.deqq[j][0] -= dedx;
                   deriv.deqq[j][1] -= dedy;
                   deriv.deqq[j][2] -= dedz;
                  virial.virx += xr*dedx;
                  virial.viry += yr*dedy;
                  virial.virz += zr*dedz;
                   energies.eu += e;
                }
               }
             }
           }
       }
   }
}

void echarge2(int i)
{
    int j,k;
    double fik,de,d2e,d2edx,d2edy,d2edz;
    double xi,yi,zi,xr,yr,zr,term[3][3];
    double xk,yk,zk,rdn;
    double r,r2;
    double cc, rdielc;
    double cutoff;
    
    if (fabs(atom[i].charge) < 0.01)
       return;

    cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
       
    cc = 4.80298*4.80298*14.39418;
    rdielc = 1.0/units.dielec;
    rdn = 1.0;
    if (field.type == MMX || field.type == MM2 )
     rdn = .915;
    else if (field.type == MM3)
     rdn = .923;
    
    for (j=1; j <= natom; j++)
      skip[j] = 0;
      
    skip[i] = i;
    for (k=0; k < MAXIAT; k++)
    {
        if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
          skip[atom[i].iat[k]] = i;
    }
    for(k=0; k < attached.n13[i]; k++)
       skip[attached.i13[k][i]] = i;
    for(k=0; k < attached.n14[i]; k++)
       skip[attached.i14[k][i]] = -i;
    
    xi = atom[i].x;
    yi = atom[i].y;
    zi = atom[i].z;

    
    for (k=1; k <= natom; k++)
    {
        if (i != k && atom[k].charge != 0.0  &&  (skip[k] != i) )
        {
               xk = atom[k].x;
               yk = atom[k].y;
               zk = atom[k].z;

                xr = xi - xk;
                yr = yi - yk;
                zr = zi - zk;
            r2 = xr*xr + yr*yr + zr*zr;
            if (r2 < cutoff)
            {
              r = sqrt(r2);
              fik = cc*rdielc*atom[i].charge*atom[k].charge;
              if (skip[k] == -i)
                 fik /= units.chgscale;
              de = -fik/r2;
              d2e = -2.0*de/r;
              de = de / r;
              d2e = (d2e-de) / r2;
              d2edx = d2e * xr;
              d2edy = d2e * yr;
              d2edz = d2e * zr;
              term[0][0]= d2edx*xr + de;
              term[1][0] = d2edx*yr;
              term[2][0] = d2edx*zr;
              term[0][1] = term[1][0];
              term[1][1] = d2edy*yr + de;
              term[2][1] = d2edy*zr;
              term[0][2] = term[2][0];
              term[1][2] = term[2][1];
              term[2][2] = d2edz*zr + de;
              for (j=0; j < 3; j++)
              {
                  hess.hessx[i][j] += term[j][0];
                  hess.hessy[i][j] += term[j][1];
                  hess.hessz[i][j] += term[j][2];
                  hess.hessx[k][j] -= term[j][0];
                  hess.hessy[k][j] -= term[j][1];
                  hess.hessz[k][j] -= term[j][2];
              }
            }
        }
    }
}
