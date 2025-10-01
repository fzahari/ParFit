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
#include "derivs.h"
#include "hess.h"
#include "field.h"
#include "cutoffs.h"
#include "units.h"
#include "minim_values.h"
#include "coordb.h"        
        
void ecoord_bond()
{
      int i,ia, ib, ic;
      int ia1,ia2;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk, xl,yl,zl;
      double x1,y1,z1,x2,y2,z2,x3,y3,z3,rpk,cospi;
      double cospk;
      double e,p,p2,p6,p12;
      double rik,rik2, embe;
      double expcut,expcut2,expterm,expmerge;
      double disik, disik2;
      double cutoff;

      cutoff = cutoffs.vdwcut*cutoffs.vdwcut;
      p6 = xl = yl = zl = 0.0;
      if (minim_values.iprint)
      {
          fprintf(pcmoutfile,"\nCoord Bonds - Buckingham Potential\n");
          fprintf(pcmoutfile,"                 At1        At2    Rik    Radius     Eps          Evdw\n");
      }
      expcut = 2.0;
      expcut2 = expcut * expcut;
      expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));

      for (i=0; i < cbonds.nbnd; i++)
      {
        ia = cbonds.icbonds[i][0];  // metal
        ib = cbonds.icbonds[i][1];  // lp or pi atom
        ic = cbonds.icbonds[i][2];

        if (atom[ia].use || atom[ib].use || atom[ic].use)
        {

          xi = atom[ia].x;
          yi = atom[ia].y;
          zi = atom[ia].z;

          xk = atom[ib].x;
          yk = atom[ib].y;
          zk = atom[ib].z;

          if (ic > 0)
          {
              xl = atom[ic].x;
              yl = atom[ic].y;
              zl = atom[ic].z;
          }
          xr = xi - xk;
          yr = yi - yk;
          zr = zi - zk;

          rik2 = (xr*xr + yr*yr + zr*zr);
          if (rik2 < cutoff)
          {
              if (ic == 0)  // metal to pi atom
              {
                  if (atom[ib].type == 4)
                    cospk = 1.0;
                  else
                  {
//                    coscalc(ib, cbonds.ia1[i], cbonds.ia2[i], ia, &cospk);
                        ia1 = cbonds.ia1[i];
                        ia2 = cbonds.ia2[i];
                        rik = sqrt(rik2);
                        x1 = atom[ia1].x - atom[ib].x;
                        y1 = atom[ia1].y - atom[ib].y;
                        z1 = atom[ia1].z - atom[ib].z;
                        x2 = atom[ia2].x - atom[ib].x;
                        y2 = atom[ia2].y - atom[ib].y;
                        z2 = atom[ia2].z - atom[ib].z;
                        x3 = y1*z2 - z1*y2;
                        y3 = z1*x2 - x1*z2;
                        z3 = x1*y2 - y1*x2;
                        rpk = sqrt( x3*x3 + y3*y3 + z3*z3 );
                        cospi = (x3*xr + y3*yr + z3*zr)/(rpk*rik);
                        cospk = (cospi*cospi);
                  }
                  embe = cospk*cbonds.lpd[i];
                  p2 = cbonds.vrad[i]*cbonds.vrad[i]/rik2;
                  p6 = p2*p2*p2;
                  if (p2 < 14.0)
                  {
                    p = sqrt(p2);
                    expterm = 15.0*units.aterm*exp(-units.bterm/p);
                    e = cbonds.veps[i] *(expterm-embe*p2);   // *units.cterm
                  } else
                  {
                    p12 = p6*p6;
                    e = expmerge*cbonds.veps[i]*p12;
                  }
              } else  // metal to lp-donor
              {
                  disik2 = (xi-xl)*(xi-xl) + (yi-yl)*(yi-yl) + (zi-zl)*(zi-zl);
                  disik = sqrt(disik2);
                  
                  embe = (2.0*sqrt(rik2)/disik)*cbonds.lpd[i];
                  p2 = cbonds.vrad[i]*cbonds.vrad[i]/rik2;
                  if (p2 < 14.0)
                  {
                    p = sqrt(p2);
                    expterm = 15.0*units.aterm*exp(-units.bterm/p);
                    e = cbonds.veps[i] *(expterm-embe*p2);   // *units.cterm
                  } else
                  {
                    p12 = p6*p6;
                    e = expmerge*cbonds.veps[i]*p12;
                  }
              }  
              energies.evdw += e;
              atom[ia].energy += e;
              atom[ib].energy += e;
              if (minim_values.iprint)
              {
                 if (cbonds.icbonds[i][2] == 0 )
                   fprintf(pcmoutfile,"CBND: %2s(%-3d) - %2s(%-3d) %-8.3f %-8.3f %-8.3f = %8.4f\n",
                      atom[ia].name, ia, atom[ib].name, ib, sqrt(rik2),cbonds.vrad[i],cbonds.veps[i],e);
                 else
                   fprintf(pcmoutfile,"CBND: %2s(%-3d) - %2s(%-3d) %-8.3f %-8.3f %-8.3f = %8.4f\n",
                      atom[ia].name, ia, atom[ic].name, ic, sqrt(rik2),cbonds.vrad[i],cbonds.veps[i],e);
              }
          }
        }
      }
}

void ecoord_bond1()
{
      int i,ia, ib, ic, ia1,ia2;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk, xl,yl,zl;
      double x1,y1,z1,x2,y2,z2,x3,y3,z3,rpk,rpk2;
      double xa1,ya1,za1,xa2,ya2,za2;
      double cospk,cospi,vrad,eps;
      double e,p,p2,p6,p12;
      double rik,rik2, embe, rvterm;
      double expcut,expcut2,expterm,expmerge;
      double dedx,dedy,dedz,de;
      double cutoff, disik, disik2;

      cutoff = cutoffs.vdwcut*cutoffs.vdwcut;

      expcut = 2.0;
      expcut2 = expcut * expcut;
      expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));
      expterm = rvterm = 0.0;
      xl = yl = zl = 0.0;
      disik = disik2 =  0.0;
      p6 = cospi = za2 = ya2 = xa2 = za1 = ya1 = xa1 = 0.0;
      rpk2 = rpk = x3 = y3 = z3 = x2 = y2 = z2 = x1 = y1 = z1 = 0.0;
      ia1 = ia2 = 0;

      for (i=0; i < cbonds.nbnd; i++)
      {
        ia = cbonds.icbonds[i][0];
        ib = cbonds.icbonds[i][1];
        ic = cbonds.icbonds[i][2];

        if (atom[ia].use || atom[ib].use || atom[ic].use)
        {
          xi = atom[ia].x;
          yi = atom[ia].y;
          zi = atom[ia].z;

          xk = atom[ib].x;
          yk = atom[ib].y;
          zk = atom[ib].z;
          if (ic > 0)
          {
              xl = atom[ic].x;
              yl = atom[ic].y;
              zl = atom[ic].z;
          }

          xr = xi - xk;
          yr = yi - yk;
          zr = zi - zk;

          rik2 = (xr*xr + yr*yr + zr*zr);
          if (rik2 < cutoff)
          {
              rik = sqrt(rik2);
              vrad = cbonds.vrad[i];
              eps = cbonds.veps[i];
              rvterm = -units.bterm/cbonds.vrad[i];
              if (ic == 0)  // metal to pi atom
              {
                  if (atom[ib].type == 4)
                    cospk = 1.0;
                  else
                  {
//                    coscalc(ib, cbonds.ia1[i], cbonds.ia2[i], ia, &cospk);
                        ia1 = cbonds.ia1[i];
                        ia2 = cbonds.ia2[i];
                        x1 = atom[ia1].x - atom[ib].x;
                        y1 = atom[ia1].y - atom[ib].y;
                        z1 = atom[ia1].z - atom[ib].z;
                        x2 = atom[ia2].x - atom[ib].x;
                        y2 = atom[ia2].y - atom[ib].y;
                        z2 = atom[ia2].z - atom[ib].z;
                        xa1 = atom[ia1].x;
                        ya1 = atom[ia1].y;
                        za1 = atom[ia1].z;
                        xa2 = atom[ia2].x;
                        ya2 = atom[ia2].y;
                        za2 = atom[ia2].z;
                        
                        x3 = y1*z2 - z1*y2;
                        y3 = z1*x2 - x1*z2;
                        z3 = x1*y2 - y1*x2;
                        rpk2 = x3*x3 + y3*y3 + z3*z3;
                        rpk = sqrt(rpk2);

                        cospi = (x3*xr + y3*yr + z3*zr)/(rpk*rik);
                        cospk = (cospi*cospi);
                  }
                embe = cospk*cbonds.lpd[i];
                p2 = vrad*vrad/rik2;
                p6 = p2*p2*p2;
                if (p2 < 14.0)
                {
                  p = sqrt(p2);
                  expterm = 15.0*units.aterm*exp(-units.bterm/p);
                  e = eps*(expterm-embe*p2);   // *units.cterm
                  de = eps*(rvterm*expterm+2.0*embe*p2/rik);  // *units.cterm
                } else
                {
                  p12 = p6*p6;
                  e = expmerge*eps*p12;
                  de = -12.0*expmerge*eps*p12/rik;
                }
              } else  // metal to lp-donor
              {
                  disik2 = (xi-xl)*(xi-xl) + (yi-yl)*(yi-yl) + (zi-zl)*(zi-zl);
                  disik = sqrt(disik2);
                  
                  embe = (2.0*sqrt(rik2)/disik)*cbonds.lpd[i];
                  p2 = cbonds.vrad[i]*cbonds.vrad[i]/rik2;
                  if (p2 < 14.0)
                  {
                    p = sqrt(p2);
                    expterm = 15.0*units.aterm*exp(-units.bterm/p);
                    e = cbonds.veps[i] *(expterm-embe*p2);   // *units.cterm
                    de = cbonds.veps[i] *(rvterm*expterm+2.0*embe*p2/rik);  // *units.cterm 
                  } else
                  {
                    p12 = p6*p6;
                    e = expmerge*cbonds.veps[i]*p12;
                    de = -12.0*expmerge*cbonds.veps[i]*p12/rik;
                  }
              }  
 
              energies.evdw += e;
              de /= rik;
              dedx = de * xr;
              dedy = de * yr;
              dedz = de * zr;
              if (ic == 0)
              {
                  if (atom[ib].type == 4)
                  {
                        deriv.devdw[ia][0] += dedx;
                        deriv.devdw[ia][1] += dedy;
                        deriv.devdw[ia][2] += dedz;
                        deriv.devdw[ib][0] -= dedx;
                        deriv.devdw[ib][1] -= dedy;
                        deriv.devdw[ib][2] -= dedz;
                  } else
                  {
                        deriv.devdw[ia][0] += eps*(rvterm*expterm*xr/rik + 4.0*embe*p2*xr/rik2 - 2.0*cbonds.lpd[i]*p2*x3*cospi/(rpk*rik));
                        deriv.devdw[ia][1] += eps*(rvterm*expterm*yr/rik + 4.0*embe*p2*yr/rik2 - 2.0*cbonds.lpd[i]*p2*y3*cospi/(rpk*rik));
                        deriv.devdw[ia][2] += eps*(rvterm*expterm*zr/rik + 4.0*embe*p2*zr/rik2 - 2.0*cbonds.lpd[i]*p2*z3*cospi/(rpk*rik));
                 
                        deriv.devdw[ib][0] += eps*(-rvterm*expterm*xr/rik - 4.0*embe*p2*xr/rik2
                         - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(yr*(za2-za1)+y2*z1-y1*z2+zr*(ya1-ya2))/(rpk2*rik2))
                         + (cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*2.0*(z3*(ya1-ya2)+y3*(za2-za1))/(rpk2*rpk2*rik2)) );
                        deriv.devdw[ib][1] += eps*(-rvterm*expterm*yr/rik - 4.0*embe*p2*yr/rik2
                          - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(zr*(xa2-xa1)-x2*z1+x1*z2+xr*(za1-za2))/(rpk2*rik2))
                          + (cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*2.0*(x3*(za1-za2)+z3*(xa2-xa1))/(rpk2*rpk2*rik2)) );
                        deriv.devdw[ib][2] += eps*(-rvterm*expterm*zr/rik - 4.0*embe*p2*zr/rik2
                         - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*(ya2-ya1)+x2*y1-x1*y2+yr*(xa1-xa2))/(rpk2*rik2))
                          + (cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*2.0*(y3*(xa1-xa2)+x3*(ya2-ya1))/(rpk2*rpk2*rik2)) );
                   
                        deriv.devdw[ia1][0] += eps*( cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*(2.0*(z3*y2-y3*z2))/(rpk2*rpk2*rik2)
                          - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(y2*zr-z2*yr)/(rpk2*rik2) ));
                        deriv.devdw[ia1][1] += eps*( cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*(2.0*(x3*z2-z3*x2))/(rpk2*rpk2*rik2)
                           - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(z2*xr-x2*zr)/(rpk2*rik2) ));
                        deriv.devdw[ia1][2] += eps*( cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*(2.0*(y3*x2-x3*y2))/(rpk2*rpk2*rik2)
                           - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(x2*yr-y2*xr)/(rpk2*rik2) ));
                    
                        deriv.devdw[ia2][0] += eps*( cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*(2.0*(y3*z1-z3*y1))/(rpk2*rpk2*rik2)
                          - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(z1*yr-zr*y1)/(rpk2*rik2) ));
                                      
                        deriv.devdw[ia2][1] += eps*( cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*(2.0*(z3*x1-x3*z1))/(rpk2*rpk2*rik2)
                          - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(x1*zr-xr*z1)/(rpk2*rik2) ));
                    
                        deriv.devdw[ia2][2] += eps*( cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(xr*x3+yr*y3+zr*z3)*(2.0*(x3*y1-y3*x1))/(rpk2*rpk2*rik2)
                          - (2.0*cbonds.lpd[i]*p2*(xr*x3+yr*y3+zr*z3)*(y1*xr-yr*x1)/(rpk2*rik2) ));
                  }
              }else
              {
                 deriv.devdw[ia][0] += cbonds.veps[i]*(rvterm*expterm*xr/rik + embe*p2*(xi-xl)/disik2 + embe*p2*xr/rik2);
                 deriv.devdw[ia][1] += cbonds.veps[i]*(rvterm*expterm*yr/rik + embe*p2*(yi-yl)/disik2 + embe*p2*yr/rik2);
                 deriv.devdw[ia][2] += cbonds.veps[i]*(rvterm*expterm*zr/rik + embe*p2*(zi-zl)/disik2 + embe*p2*zr/rik2);
                 
                 deriv.devdw[ib][0] -= cbonds.veps[i]*(rvterm*expterm*xr/rik + embe*p2*xr/rik2);
                 deriv.devdw[ib][1] -= cbonds.veps[i]*(rvterm*expterm*yr/rik + embe*p2*yr/rik2);
                 deriv.devdw[ib][2] -= cbonds.veps[i]*(rvterm*expterm*zr/rik + embe*p2*zr/rik2);
                 
                 deriv.devdw[ic][0] -= cbonds.veps[i]*embe*p2*(xi-xl)/disik2;
                 deriv.devdw[ic][1] -= cbonds.veps[i]*embe*p2*(yi-yl)/disik2;
                 deriv.devdw[ic][2] -= cbonds.veps[i]*embe*p2*(zi-zl)/disik2;
              }
          }
        }
      }
}

void ecoord_bond2(int iatom)
{
   int i,ia, ib, ic, j;
   double xi,yi,zi,xr,yr,zr;
   double xk,yk,zk,xl,yl,zl;
   double cospk;
   double e,p,p2,p6,p12;
   double rik,rik2, embe, rvterm, rvterm2;
   double expcut,expcut2,expterm,expmerge;
   double de, d2e;
   double d2edx,d2edy,d2edz,term[3][3];
   double disik, disik2;
   double cutoff, term1, term2, term3;
   double dedxiaxia,dedxiayia,dedxiazia,dedxiaxib,dedxiayib,dedxiazib,dedxiaxic,dedxiayic,dedxiazic;
   double dedxibyia,dedxibzia,dedxibxib,dedxibyib,dedxibzib,dedxibxic,dedxibyic,dedxibzic;
   double dedxicyia,dedxiczia,dedxicyib,dedxiczib,dedxicxic,dedxicyic,dedxiczic;
   double dedyiayia,dedyiazia,dedyiayib,dedyiazib,dedyiayic,dedyiazic;
   double dedyibzia,dedyibyib,dedyibzib,dedyibyic,dedyibzic;
   double dedyiczia,dedyiczib,dedyicyic,dedyiczic;
   double dedziazia,dedziazib,dedziazic;
   double dedzibzib,dedzibzic,dedziczic;
   
   cutoff = cutoffs.vdwcut*cutoffs.vdwcut;

   expcut = 2.0;
   expcut2 = expcut * expcut;
   expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));
   expterm = rvterm = 0.0;
   xl = yl = zl = 0.0;

   for (i=0; i < cbonds.nbnd; i++)
   {
       ia = cbonds.icbonds[i][0];
       ib = cbonds.icbonds[i][1];
       ic = cbonds.icbonds[i][2];
       if (ia == iatom || ib == iatom || ic == iatom)
       {
          
          xi = atom[ia].x;
          yi = atom[ia].y;
          zi = atom[ia].z;

          xk = atom[ib].x;
          yk = atom[ib].y;
          zk = atom[ib].z;

          if (ic > 0)
          {
             xl = atom[ic].x;
             yl = atom[ic].y;
             zl = atom[ic].z;
          }
          xr = xi - xk;
          yr = yi - yk;
          zr = zi - zk;

          rik2 = (xr*xr + yr*yr + zr*zr);
          if (rik2 < cutoff)
          {
               rik = sqrt(rik2);
               if (ic == 0)
               {
                  if (atom[ib].type == 4)
                    cospk = 1.0;
                  else
                    coscalc(ib, cbonds.ia1[i], cbonds.ia2[i], ia, &cospk);
                 embe = cospk*cbonds.lpd[i];
                 if (rik2 < 49.0)
                 {
                   p2 = cbonds.vrad[i]*cbonds.vrad[i]/rik2;
                   p6 = p2*p2*p2;
                   if (p2 < 14.0)
                   {
                    p = sqrt(p2);
                    rvterm = -units.bterm/cbonds.vrad[i];
                    rvterm2 = rvterm*rvterm;
                    expterm = 15.0*units.aterm*exp(-units.bterm/p);
                    e = cbonds.veps[i] *(expterm-embe*p6);                  // *units.cterm
                    de = cbonds.veps[i] *(rvterm*expterm+2.0*embe*p2/rik);  // *units.cterm 
                    d2e = cbonds.veps[i] * (rvterm2*expterm-6.0*embe*p2/rik2); // units.cterm
                   } else
                   {
                    p12 = p6*p6;
                    e = expmerge*cbonds.veps[i]*p12;
                    de = -12.0*expmerge*cbonds.veps[i]*p12/rik;
                    d2e = 156.0 *expmerge *cbonds.veps[i]*p12/rik;
                   }

                   de /= rik;
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
            
                   for(j=0; j < 3; j++)
                   {
                        hess.hessx[ia][j] += term[j][0];
                        hess.hessy[ia][j] += term[j][1];
                        hess.hessz[ia][j] += term[j][2];
                        hess.hessx[ib][j] -= term[j][0];
                        hess.hessy[ib][j] -= term[j][1];
                        hess.hessz[ib][j] -= term[j][2];
                   }
                 }
               } else
               {
                  disik2 = (xi-xl)*(xi-xl) + (yi-yl)*(yi-yl) + (zi-zl)*(zi-zl);
                  disik = sqrt(disik2);
                  
                  embe = (2.0*sqrt(rik2)/disik)*cbonds.lpd[i];
                  p2 = cbonds.vrad[i]*cbonds.vrad[i]/rik2;
                  
                   if (rik2 < 49.0)
                   {
                     p2 = cbonds.vrad[i]*cbonds.vrad[i]/rik2;
                     p6 = p2*p2*p2;
                     if (p2 < 14.0)
                     {
                       p = sqrt(p2);
                       rvterm = -units.bterm/cbonds.vrad[i];
                       rvterm2 = rvterm*rvterm;
                       expterm = 15.0*units.aterm*exp(-units.bterm/p);
                       e = cbonds.veps[i] *(expterm-embe*p2);   // *units.cterm
                       de = cbonds.veps[i] *(rvterm*expterm+2.0*embe*p2/rik);  // *units.cterm 
                       d2e = cbonds.veps[i] * (rvterm2*expterm-6.0*embe*p2/rik2); // units.cterm
                     } else
                     {
                       p12 = p6*p6;
                       e = expmerge*cbonds.veps[i]*p12;
                       de = -12.0*expmerge*cbonds.veps[i]*p12/rik;
                       d2e = 156.0 *expmerge *cbonds.veps[i]*p12/rik;
                     }
                     term1 = cbonds.veps[i] *rvterm*expterm/rik2;
                     term2 = cbonds.veps[i] *rvterm*rvterm*expterm/rik2;
                     term3 = cbonds.veps[i] *embe*p2;
                     
                     dedxiaxia = -term1*xr*xr/rik +term1*rik + term2*xr*xr -3.0*term3*(xi-xl)*(xi-xl)/(disik2*disik2)
                                 -2.0*term3*xr*(xi-xl)/(disik2*rik2) + term3/disik2 - 3.0*term3*xr*xr/(rik2*rik2) + term3/rik2;
                     dedxiayia = -term1*xr*yr/rik + term2*xr*yr -3.0*term3*(xi-xl)*(yi-yl)/(disik2*disik2) - term3*yr*(xi-xl)/rik2*disik2
                                 -term3*xr*(yi-yl)/(rik2*disik2) - 3.0*term3*xr*yr/(rik2*rik2);
                     dedxiazia = -term1*xr*zr/rik + term2*xr*zr -3.0*term3*(xi-xl)*(zi-zl)/(disik2*disik2) - term3*zr*(xi-xl)/rik2*disik2
                                 -term3*xr*(zi-zl)/(rik2*disik2) - 3.0*term3*xr*zr/(rik2*rik2);

                     dedxiaxib = term1*xr*xr/rik - term1*rik - term2*xr*xr +3.0*term3*xr*xr/(rik2*rik2) + term3*xr*(xi-xl)/(rik2*disik2) - term3/rik2;
                     dedxiayib = term1*xr*yr/rik - term2*xr*yr +3.0*term3*xr*yr/(rik2*rik2) + term3*yr*(xi-xl)/(rik2*disik2);
                     dedxiazib = term1*xr*zr/rik - term2*xr*zr +3.0*term3*xr*zr/(rik2*rik2) + term3*zr*(xi-xl)/(rik2*disik2);

                     dedxiaxic = 3.0*term3*(xi-xl)*(xi-xl)/(disik2*disik2) + term3*xr*(xi-xl)/(disik2*rik2) - term3/disik2;
                     dedxiayic = 3.0*term3*(xi-xl)*(yi-yl)/(disik2*disik2) + term3*xr*(yi-yl)/(disik2*rik2);                                 
                     dedxiazic = 3.0*term3*(xi-xl)*(zi-zl)/(disik2*disik2) + term3*xr*(zi-zl)/(disik2*rik2);                                 

                     dedxibyia = term1*xr*yr/rik - term2*xr*yr + term3*xr*(yi-yl)/(disik2*rik2) + 3.0*term3*xr*yr/(rik2*rik2);
                     dedxibzia = term1*xr*zr/rik - term2*xr*zr + term3*xr*(zi-zl)/(disik2*rik2) + 3.0*term3*xr*zr/(rik2*rik2);
                     dedxibxib = -term1*xr*xr/rik + term1*rik + term2*xr*xr - 3.0*term3*xr*xr/(rik2*rik2) + term3/rik2;
                     dedxibyib = -term1*xr*yr/rik + term2*xr*yr - 3.0*term3*xr*yr/(rik2*rik2);
                     dedxibzib = -term1*xr*zr/rik + term2*xr*zr - 3.0*term3*xr*zr/(rik2*rik2);
                     dedxibxic = -term3*xr*(xi-xl)/(rik2*disik2);
                     dedxibyic = -term3*xr*(yi-yl)/(rik2*disik2);
                     dedxibzic = -term3*xr*(zi-zl)/(rik2*disik2);
                     
                     dedxicyia = 3.0*term3*(xi-xl)*(yi-yl)/(disik2*disik2) + term3*(xi-xl)*yr/(rik2*disik2);
                     dedxiczia = 3.0*term3*(xi-xl)*(zi-zl)/(disik2*disik2) + term3*(xi-xl)*zr/(rik2*disik2);
                     dedxicyib = -term3*(xi-xl)*yr/(rik2*disik2);
                     dedxiczib = -term3*(xi-xl)*zr/(rik2*disik2);
                     dedxicxic = -3.0*term3*(xi-xl)*(xi-xl)/(disik2*disik2) + term3/(disik2);
                     dedxicyic = -3.0*term3*(xi-xl)*(yi-yl)/(disik2*disik2);
                     dedxiczic = -3.0*term3*(xi-xl)*(zi-zl)/(disik2*disik2);

                     dedyiayia = -term1*yr*yr/rik +term1*rik + term2*yr*yr -3.0*term3*(yi-yl)*(yi-yl)/(disik2*disik2)
                                 -2.0*term3*yr*(yi-yl)/(disik2*rik2) + term3/disik2 - 3.0*term3*yr*yr/(rik2*rik2) + term3/rik2;
                     dedyiazia = -term1*yr*zr/rik + term2*yr*zr -3.0*term3*(yi-yl)*(zi-zl)/(disik2*disik2) - term3*zr*(yi-yl)/rik2*disik2
                                 -term3*yr*(zi-zl)/(rik2*disik2) - 3.0*term3*yr*zr/(rik2*rik2);
                     dedyiayib = term1*yr*yr/rik - term1*rik - term2*yr*yr +3.0*term3*yr*yr/(rik2*rik2) + term3*yr*(yi-yl)/(rik2*disik2) - term3/rik2;
                     dedyiazib = term1*yr*zr/rik - term2*yr*zr +3.0*term3*yr*zr/(rik2*rik2) + term3*zr*(yi-yl)/(rik2*disik2);
                     dedyiayic = 3.0*term3*(yi-yl)*(yi-yl)/(disik2*disik2) + term3*yr*(yi-yl)/(disik2*rik2) - term3/disik2;
                     dedyiazic = 3.0*term3*(yi-yl)*(zi-zl)/(disik2*disik2) + term3*yr*(zi-zl)/(disik2*rik2);
                     
                     dedyibzia = term1*yr*zr/rik - term2*yr*zr + term3*yr*(zi-zl)/(disik2*rik2) + 3.0*term3*yr*zr/(rik2*rik2);
                     dedyibyib = -term1*yr*yr/rik + term1*rik + term2*yr*yr - 3.0*term3*yr*yr/(rik2*rik2) + term3/rik2;
                     dedyibzib = -term1*yr*zr/rik + term2*yr*zr - 3.0*term3*yr*zr/(rik2*rik2);
                     dedyibyic = -term3*yr*(yi-yl)/(rik2*disik2);
                     dedyibzic = -term3*yr*(zi-zl)/(rik2*disik2);

                     dedyiczia = 3.0*term3*(yi-yl)*(zi-zl)/(disik2*disik2) + term3*(yi-yl)*zr/(rik2*disik2);
                     dedyiczib = -term3*(yi-yl)*zr/(rik2*disik2);
                     dedyicyic = -3.0*term3*(yi-yl)*(yi-yl)/(disik2*disik2) + term3/disik2;
                     dedyiczic = -3.0*term3*(yi-yl)*(zi-zl)/(disik2*disik2);

                     dedziazia = -term1*zr*zr/rik +term1*rik + term2*zr*zr -3.0*term3*(zi-zl)*(zi-zl)/(disik2*disik2)
                                 -2.0*term3*zr*(zi-zl)/(disik2*rik2) + term3/disik2 - 3.0*term3*zr*zr/(rik2*rik2) + term3/rik2;               
                     dedziazib = term1*zr*zr/rik - term1*rik - term2*zr*zr +3.0*term3*zr*zr/(rik2*rik2) + term3*zr*(zi-zl)/(rik2*disik2) - term3/rik2;
                     dedziazic = 3.0*term3*(zi-zl)*(zi-zl)/(disik2*disik2) + term3*zr*(zi-zl)/(disik2*rik2) - term3/disik2;

                     dedzibzib = -term1*zr*zr/rik + term1*rik + term2*zr*zr - 3.0*term3*zr*zr/(rik2*rik2) + term3/rik2;
                     dedzibzic = -term3*zr*(zi-zl)/(rik2*disik2);
                     dedziczic = -3.0*term3*(zi-zl)*(zi-zl)/(disik2*disik2) + term3/(disik2);
                                  
                       if (iatom == ia)
                       {
                          hess.hessx[ia][0] += dedxiaxia ;
                          hess.hessy[ia][0] += dedxiayia ;
                          hess.hessz[ia][0] += dedxiazia ;
                          hess.hessx[ia][1] += dedxiayia ;
                          hess.hessy[ia][1] += dedyiayia ;
                          hess.hessz[ia][1] += dedyiazia ;
                          hess.hessx[ia][2] += dedxiazia ;
                          hess.hessy[ia][2] += dedyiazia ;
                          hess.hessz[ia][2] += dedziazia ;
                          hess.hessx[ib][0] += dedxiaxib ;
                          hess.hessy[ib][0] += dedxibyia ;
                          hess.hessz[ib][0] += dedxibzia ;
                          hess.hessx[ib][1] += dedxiayib ;
                          hess.hessy[ib][1] += dedyiayib ;
                          hess.hessz[ib][1] += dedyibzia ;
                          hess.hessx[ib][2] += dedxiazib ;
                          hess.hessy[ib][2] += dedyiazib ;
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
                       }
                   }
               }  
          }
       }
   }
}

void coscalc(int i,int ia1,int ia2,int k,double *cospi)
{
        double rik, rpk, x1, x2, x3, xn, y1, y2, y3, yn, z1, z2, z3, zn;

        /* calculate 1-(1-cosine)**2 between p ao on atom i and vector i,k */

        x1 = atom[ia1].x - atom[i].x;
        y1 = atom[ia1].y - atom[i].y;
        z1 = atom[ia1].z - atom[i].z;
        x2 = atom[ia2].x - atom[i].x;
        y2 = atom[ia2].y - atom[i].y;
        z2 = atom[ia2].z - atom[i].z;
        x3 = y1*z2 - z1*y2;
        y3 = z1*x2 - x1*z2;
        z3 = x1*y2 - y1*x2;
        rpk = sqrt( x3*x3 + y3*y3 + z3*z3 );
        xn = atom[k].x - atom[i].x;
        yn = atom[k].y - atom[i].y;
        zn = atom[k].z - atom[i].z;
        rik = sqrt( xn*xn + yn*yn + zn*zn );
        *cospi = fabs( (x3*xn + y3*yn + z3*zn)/(rpk*rik) );
        /*jjg  stmt below for cos **2 dependence
         *         cospi=cospi*cospi */
        *cospi = 1. - ((1. - *cospi)*(1. - *cospi));
        return;
} 

