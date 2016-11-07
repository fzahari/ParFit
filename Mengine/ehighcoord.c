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

#include "energies.h"
#include "angles.h"
#include "derivs.h"
#include "hess.h"
#include "nonbond.h"
#include "utility.h"
#include "field.h"
#include "cutoffs.h"
#include "attached.h"
#include "pot.h"
#include "solv.h"
#include "units.h"
#include "minim_values.h"

void ehigh_coord(void);
void ehigh_coord1(void);
void ehigh_coord2(int);
void image(double *, double *, double *, int);
int isbond(int,int);
int isangle(int,int);

void ehigh_coord(void)
{
   int    i,ia, ib, ita, itb;
   double xi, yi, zi, xk,yk,zk;
   double xr, yr, zr;
   double e, rik2;
   double p,p2,p6,p12,rv,eps;
   double expcut,expterm,expmerge;
   double pif,posterm;

   

   expcut = 2.0;
   expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));
   if (minim_values.iprint)
   {
          fprintf(pcmoutfile,"\nPOS Metal Atom Angles (At1 - Metal - At2\n");
          fprintf(pcmoutfile,"                 At1        At2    Rik    Radius     Eps   Posterm       Evdw\n");
   }

   for (i=0; i < high_coord.ncoord; i++)
   {
       ia = high_coord.i13[i][0];
       ib = high_coord.i13[i][2];
      
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
        rv = nonbond.vrad[ita][itb];
        eps = nonbond.veps[ita][itb];
        pif = (double)nonbond.ipif[ita][itb];
	posterm = 1.0;
	if (isbond(ia,ib)) posterm = units.pos_12scale;
        if (isangle(ia,ib)) posterm = units.pos_13scale;
	   
        p2 = (rv*rv)/rik2;
        p6 = p2*p2*p2;
        if (p2 < 14.0)
        {
            p = sqrt(p2);
            expterm = units.aterm*exp(-units.bterm/p);
            e = posterm*eps *(expterm- pif*units.cterm*p6);
        } else
        {
            p12 = p6*p6;
            e = posterm*expmerge*eps*p12;
        }
        energies.evdw += e;
        atom[ia].energy += e;
        atom[ib].energy += e;
        if (minim_values.iprint)
        {
           if ( e > 0.1 || e < -0.1)            
               fprintf(pcmoutfile,"High Coord 1,3: %2s(%-3d) - %2s(%-3d) %-8.3f %-8.3f %-8.3f  %-4.2f = %8.4f\n",
                             atom[ia].name, ia, atom[ib].name, ib, sqrt(rik2),rv,eps,posterm,e);
        }
   }
}
/* ========================================================= */
void ehigh_coord1(void)
{
   int     i,ia, ib, ita, itb;
   double xi, yi, zi, xk,yk,zk;
   double xr, yr, zr;
   double e, rik2;
   double p,p2,p6,p12,rv,eps;
   double expcut,expterm,expmerge;
   double dedx,dedy,dedz,de;
   double rik,rvterm;
   double pif,posterm;

   expcut = 2.0;
   expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));

   for (i=0; i < high_coord.ncoord; i++)
   {
      ia = high_coord.i13[i][0];
      ib = high_coord.i13[i][2];
      
      ita = atom[ia].type;
      itb = atom[ib].type;
      if (atom[ia].use || atom[ib].use)
      {
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
        rv = nonbond.vrad[ita][itb];
        eps = nonbond.veps[ita][itb];
        pif = (double)nonbond.ipif[ita][itb];
		posterm = 1.0;
		if (isbond(ia,ib)) posterm = units.pos_12scale;
		if (isangle(ia,ib)) posterm = units.pos_13scale;
		p2 = (rv*rv)/rik2;
        p6 = p2*p2*p2;
        rik = sqrt(rik2);
        if (p2 < 14.0)
        {
            p = sqrt(p2);
            rvterm = -units.bterm / rv;
            expterm = units.aterm*exp(-units.bterm/p);
            e = posterm*eps *(expterm- pif*units.cterm*p6);
            de = posterm*eps * (rvterm*expterm+6.0*pif*units.cterm*p6/rik);
        } else
        {
            p12 = p6*p6;
            e = posterm*expmerge*eps*p12;
            de = posterm*(-12.0)* expmerge * eps * p12/rik;
        }
        de /= rik;
        dedx = de * xr;
        dedy = de * yr;
        dedz = de * zr;
        energies.evdw += e;
        deriv.devdw[ia][0] += dedx;
        deriv.devdw[ia][1] += dedy;
        deriv.devdw[ia][2] += dedz;
        deriv.devdw[ib][0] -= dedx;
        deriv.devdw[ib][1] -= dedy;
        deriv.devdw[ib][2] -= dedz;
        atom[ia].energy += e;
        atom[ib].energy += e;
      }
   }
}

/* ========================================================= */
void ehigh_coord2(int iatom)
{
   int     i,ia, ib, ita, itb, k;
   double xi, yi, zi, xk,yk,zk;
   double xr, yr, zr;
   double de,d2e,p,p2,p6,p12,eps,rv;
   double rik,rik2;
   double d2edx,d2edy,d2edz,term[3][3];
   double expcut,expterm,expmerge;
   double rvterm,rvterm2;
   double pif,posterm;

   
   expcut = 2.0;
   expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));

   for (i=0; i < high_coord.ncoord; i++)
   {
      ia = high_coord.i13[i][0];
      ib = high_coord.i13[i][2];
      
      ita = atom[ia].type;
      itb = atom[ib].type;
      if (ia == iatom || ib == iatom )
      {
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
        rv = nonbond.vrad[ita][itb];
        eps = nonbond.veps[ita][itb];
        pif = (double)nonbond.ipif[ita][itb];
		posterm = 1.0;
		if (isbond(ia,ib)) posterm = units.pos_12scale;
		if (isangle(ia,ib)) posterm = units.pos_13scale;
		p2 = (rv*rv)/rik2;
        p6 = p2*p2*p2;
        rik = sqrt(rik2);
        if (p2 < 14.0)
        {
            p = sqrt(p2);
            rvterm = -units.bterm / rv;
            rvterm2 = rvterm*rvterm;
            expterm = units.aterm*exp(-units.bterm/p);
            de = posterm*eps * (rvterm*expterm+6.0*pif*units.cterm*p6/rik);
            d2e = posterm*eps * (rvterm2*expterm-42.0*pif*units.cterm*p6/rik2);
        } else
        {
            p12 = p6*p6;
            de = posterm*(-12.0)* expmerge * eps * p12/rik;
            d2e = posterm*(156.0) *expmerge *eps*p12/rik2;
        }
        de = de / rik;
        d2e = (d2e-de) / rik2;
        d2edx = d2e * xr;
        d2edy = d2e * yr;
        d2edz = d2e * zr;
        term[0][0] = d2edx*xr + de;
        term[1][0] = d2edx*yr;
        term[2][0] = d2edx*zr;
        term[0][1] = term[1][0];
        term[1][1] = d2edy*yr + de;
        term[2][1] = d2edy*zr;
        term[0][2] = term[2][0];
        term[1][2] = term[2][1];
        term[2][2] = d2edz*zr + de;
        for (k=0; k < 3; k++)
        {
            hess.hessx[ia][k] += term[k][0];
            hess.hessy[ia][k] += term[k][1];
            hess.hessz[ia][k] += term[k][2];                     
            hess.hessx[ib][k] -= term[k][0];
            hess.hessy[ib][k] -= term[k][1];
            hess.hessz[ib][k] -= term[k][2];
        }                   
      }
   }
}










