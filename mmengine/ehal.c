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

#define EXTERN extern

#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"

#include "energies.h"
#include "derivs.h"
#include "hess.h"
#include "torsions.h"
#include "nonbond.h"
#include "field.h"
#include "cutoffs.h"
#include "attached.h"
#include "pot.h"
#include "solv.h"
#include "units.h"
#include "minim_values.h"
#include "angles.h"
        
int ishbond(int, int);
int iscoord_bond(int, int);
        
EXTERN int *skip;

void image(double *, double *, double *, int);
void ehal(void);
void ehal1(void);
void ehal2(int);
      
void ehal()
{
      int i,ia, ib, ita, itb, j, k;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double e,rv,eps;
      double rik,rik2,rik7;
      double sigma, kappa, kappa7;
      double rv7,rho, tau;
      double cutoff;
      
      cutoff = cutoffs.vdwcut*cutoffs.vdwcut;
      energies.evdw = 0.0;
      energies.e14  = 0.0;
      sigma = 1.12;
      kappa = 1.07;
      kappa7 = pow(kappa,7.0);

      if (minim_values.iprint)
      {
          fprintf(pcmoutfile,"\nVDW Terms - MMFF Buffered 14-7 Potential\n");
          fprintf(pcmoutfile,"     At1   At2    Rik    Radius     Eps       Evdw\n");
      }
      for (i=1; i <= natom; i++)
         skip[i] = 0;
      for (i=1; i < natom; i++)
      {
         for (j=0; j < MAXIAT; j++)
             if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
               skip[atom[i].iat[j]] = i;
         for (j=0; j < attached.n13[i]; j++)
               skip[attached.i13[j][i]] = i;
         for(k=0; k < attached.n14[i]; k++)
            skip[attached.i14[k][i]] = -i;
         if (high_coord.ncoord > 0)
         {
            for (k=0; k < high_coord.ncoord; k++)
            {
                for (j=0; j < 3; j++)
                {
                   if (high_coord.i13[k][j] == i)
                   {
                      skip[high_coord.i13[k][0]] = i;
                      skip[high_coord.i13[k][1]] = i;
                      skip[high_coord.i13[k][2]] = i;
                      break;
                   }
                }
            }
         }

         for(j=i+1; j <= natom; j++)
         {
             if (skip[j] != i)
             {
               ia = i;
               ib = j;
               if (atom[ia].use || atom[ib].use)
               {
                 ita = atom[ia].type;
                 itb = atom[ib].type;
         
                 xi = atom[ia].x;
                 yi = atom[ia].y;
                 zi = atom[ia].z;
                 xk = atom[ib].x;
                 yk = atom[ib].y;
                 zk = atom[ib].z;

                 xr = xi - xk;
                 yr = yi - yk;
                 zr = zi - zk;
                 rik2 = xr*xr + yr*yr + zr*zr;
                 if (rik2 < cutoff)
                 {
                     rv = nonbond.vrad[ita][itb];
                     eps = nonbond.veps[ita][itb];
                     if (skip[j] == -i)
                       eps /= units.v14scale;
                    rik = sqrt(rik2);
                    rv7 = pow(rv,7.0);
                    rik7 = pow(rik,7.0);
                    rho = rik7 + (sigma-1.00)*rv7;
                    tau = rik + (kappa-1.00)*rv;

                    e = eps*kappa7*(rv7/pow(tau,7.0))*(sigma*rv7/rho-2.00);
                    if (skip[j] == -i)
                      energies.e14 += e;
                    else
                      energies.evdw += e;
                    atom[ia].energy += e;
                    atom[ib].energy += e;
                    if (minim_values.iprint)
                      fprintf(pcmoutfile,"VDW: %-4d - %-4d %-8.3f %-8.3f %-8.3f = %8.4f\n",ia, ib, sqrt(rik2),rv,eps,e);
                 }
               }
             }
         }
      }
}

void ehal1()
{
      int i,ia, ib, j, k, ita, itb;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double e,rv,eps;
      double rik,rik2;
      double rik6,rik7;
      double sigma, kappa, kappa7;
      double rv7,rho, tau, rv14;
      double dedx,dedy,dedz,de;
      double cutoff;
      
      cutoff = cutoffs.vdwcut*cutoffs.vdwcut;
      
      energies.evdw = 0.0;
      energies.e14  = 0.0;
      
      sigma = 1.12;
      kappa = 1.07;
      kappa7 = pow(kappa,7.0);
      
      for (i=1; i <= natom; i++)
      {
          deriv.devdw[i][0] = 0.0;
          deriv.devdw[i][1] = 0.0;
          deriv.devdw[i][2] = 0.0;
          deriv.de14[i][0] = 0.0;
          deriv.de14[i][1] = 0.0;
          deriv.de14[i][2] = 0.0;
      }

      for (i=1; i <= natom; i++)
         skip[i] = 0;

      for (i=1; i < natom; i++)
      {
         for (j=0; j < MAXIAT; j++)
             if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
               skip[atom[i].iat[j]] = i;
         for (j=0; j < attached.n13[i]; j++)
               skip[attached.i13[j][i]] = i;
         for(k=0; k < attached.n14[i]; k++)
            skip[attached.i14[k][i]] = -i;
         if (high_coord.ncoord > 0)
         {
            for (k=0; k < high_coord.ncoord; k++)
            {
                for (j=0; j < 3; j++)
                {
                   if (high_coord.i13[k][j] == i)
                   {
                      skip[high_coord.i13[k][0]] = i;
                      skip[high_coord.i13[k][1]] = i;
                      skip[high_coord.i13[k][2]] = i;
                      break;
                   }
                }
            }
         }
         ia = i;
         ita = atom[i].type;

         for(j=i+1; j <= natom; j++)
         {
            if (atom[i].use || atom[j].use)
            {
              if (skip[j] != i)
              {
                 ib = j;
                 itb = atom[j].type;
          
                 xi = atom[ia].x;
                 yi = atom[ia].y;
                 zi = atom[ia].z;
                 xk = atom[ib].x;
                 yk = atom[ib].y;
                 zk = atom[ib].z;

                 xr = xi - xk;
                 yr = yi - yk;
                 zr = zi - zk;
                 rik2 = xr*xr + yr*yr + zr*zr;
                 if (rik2 < cutoff)
                 {   
                     rv = nonbond.vrad[ita][itb];
                     eps = nonbond.veps[ita][itb];
                    if (skip[j] == -i)
                        eps /= units.v14scale;
                    rik = sqrt(rik2);
                    rv7 = pow(rv,7.0);
                    rv14 = rv7*rv7;
                    rik6 = rik2*rik2*rik2;
                    rik7 = rik6*rik;
                    rho = rik7 + (sigma-1.00)*rv7;
                    tau = rik + (kappa-1.00)*rv;

                    e = eps*kappa7*(rv7/pow(tau,7.0))*(sigma*rv7/rho-2.00);
         
                    de = -7.00*eps*kappa7*(rv7/pow(tau,8.00))*(sigma*rv7/rho-2.00)
                         -7.00*eps*kappa7*sigma*rv14*rik6/(rho*rho*pow(tau,7.00));

                    de /= rik;
                    dedx = de*xr;
                    dedy = de*yr;
                    dedz = de*zr;
         
                    if (skip[j] == -i)
                    {
                      energies.e14 += e;
                      deriv.de14[ia][0] += dedx;
                      deriv.de14[ia][1] += dedy;
                      deriv.de14[ia][2] += dedz;
                      deriv.de14[ib][0] -= dedx;
                      deriv.de14[ib][1] -= dedy;
                      deriv.de14[ib][2] -= dedz;
                    } else
                    {
                      energies.evdw += e;
                      deriv.devdw[ia][0] += dedx;
                      deriv.devdw[ia][1] += dedy;
                      deriv.devdw[ia][2] += dedz;
                      deriv.devdw[ib][0] -= dedx;
                      deriv.devdw[ib][1] -= dedy;
                      deriv.devdw[ib][2] -= dedz;
                    }
                    virial.virx += xr*dedx;
                    virial.viry += yr*dedy;
                    virial.virz += zr*dedz;
                 }
              }
            }
         }
      }
}

void ehal2(int iatom)
{
      int j,ia, ib, k, ita, itb;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double rv,eps;
      double rik,rik2, rik5;
      double rik6,rik7, rik12;
      double sigma, kappa, kappa7;
      double rv7,rho, tau, rv14;
      double de,d2e;
      double d2edx,d2edy,d2edz,term[3][3];
      double cutoff;
      
      cutoff = cutoffs.vdwcut*cutoffs.vdwcut;
     
      sigma = 1.12;
      kappa = 1.07;
      kappa7 = pow(kappa,7.0);

    for (j=1; j <= natom; j++)
      skip[j] = 0;

    skip[iatom] = iatom;
    for (k=0; k < MAXIAT; k++)
    {
        if (atom[iatom].iat[k] != 0 && atom[iatom].bo[k] != 9)
          skip[atom[iatom].iat[k]] = iatom;
    }
    for(k=0; k < attached.n13[iatom]; k++)
       skip[attached.i13[k][iatom]] = iatom;
    for(k=0; k < attached.n14[iatom]; k++)
       skip[attached.i14[k][iatom]] = -iatom;
          if (high_coord.ncoord > 0)
          {
             for (k=0; k < high_coord.ncoord; k++)
             {
                 for (j=0; j < 3; j++)
                 {
                     if (high_coord.i13[k][j] == iatom )
                     {
                      skip[high_coord.i13[k][0]] = iatom;
                      skip[high_coord.i13[k][1]] = iatom;
                      skip[high_coord.i13[k][2]] = iatom;
                      break;
                     }
                 }
             }
          }

    ia = iatom;
    ita = atom[ia].type;    
    xi = atom[ia].x;
    yi = atom[ia].y;
    zi = atom[ia].z;
    
   for (j=1; j <= natom; j++)
   {
       ib = j;
       if (skip[j] != iatom)
       {
          itb = atom[ib].type;

         xk = atom[ib].x;
         yk = atom[ib].y;
         zk = atom[ib].z;

         xr = xi - xk;
         yr = yi - yk;
         zr = zi - zk;
         rik2 = xr*xr + yr*yr + zr*zr;
         if (rik2 < cutoff)
         {
             rv = nonbond.vrad[ita][itb];
             eps = nonbond.veps[ita][itb];
             if (skip[j] == -iatom)
              eps /= units.v14scale;
           
            rik = sqrt(rik2);
            rv7 = pow(rv,7.0);
            rv14 = rv7*rv7;
            rik6 = rik2*rik2*rik2;
            rik5 = rik6/rik;
            rik7 = rik6*rik;
            rik12 = rik6*rik6;
            rho = rik7 + (sigma-1.00)*rv7;
            tau = rik + (kappa-1.00)*rv;

            de = -7.00*eps*kappa7*(rv7/pow(tau,8.00))*(sigma*rv7/rho-2.00)
                 -7.00*eps*kappa7*sigma*rv14*rik6/(rho*rho*pow(tau,7.00));

            d2e = 56.0*eps*kappa7*(rv7/pow(tau,9.00))*(sigma*rv7/rho-2.00)
                 + 98.0*eps*kappa7*sigma*rv14*rik6/(rho*rho*pow(tau,8.00))
                 + 98.0*eps*kappa7*sigma*rv14*rik12/(rho*rho*rho*pow(tau,7.00))
                 - 42.0*eps*kappa7*sigma*rv14*rik5/(rho*rho*pow(tau,7.00));

           de /= rik;
           d2e = (d2e-de)/rik2;
           d2edx = d2e*xr;
           d2edy = d2e*yr;
           d2edz = d2e*zr;

           term[0][0] = d2edx*xr + de;
           term[1][0] = d2edx*yr;
           term[2][0] = d2edx*zr;
           term[0][1] = term[1][0];
           term[1][1] = d2edy*yr + de;
           term[2][1] = d2edy*zr;
           term[0][2] = term[2][0];
           term[1][2] = term[2][1];
           term[2][2] = d2edz*zr + de;

           if (ia == iatom)
           {
                 for (k=0; k < 3; k++)
                 {
                     hess.hessx[ia][k] += term[k][0];
                     hess.hessy[ia][k] += term[k][1];
                     hess.hessz[ia][k] += term[k][2];
                     hess.hessx[ib][k] -= term[k][0];
                     hess.hessy[ib][k] -= term[k][1];
                     hess.hessz[ib][k] -= term[k][2];
                 }
           } else
           {
                 for (k=0; k < 3; k++)
                 {
                     hess.hessx[ia][k] -= term[k][0];
                     hess.hessy[ia][k] -= term[k][1];
                     hess.hessz[ia][k] -= term[k][2];
                     hess.hessx[ib][k] += term[k][0];
                     hess.hessy[ib][k] += term[k][1];
                     hess.hessz[ib][k] += term[k][2];
                 }
           }
         }      
      }
   }
}

