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
#include "bonds_ff.h"
#include "angles.h"
#include "field.h"
#include "cutoffs.h"
#include "energies.h"
#include "attached.h"
#include "derivs.h"
#include "hess.h"
#include "minim_values.h"
#include "skip.h"
#include "solv.h"

void esolv()
{
    int i,j;
    double e,ai,ri,rb, elocal;
    double term,probe,term1;
    double dwater, cc,f,cutoff;
    double xi,yi,zi,xr,yr,zr,rik2;
    double rb2,fgb;

    probe = 1.40;
    dwater = 78.3;
    cc = 4.80298*4.80298*14.39418;
    cutoff = cutoffs.chrgcut*cutoffs.chrgcut;

   if (minim_values.iprint)
       fprintf(pcmoutfile,"\nSolvation Energy\n");

   born();
   if (minim_values.iprint)
   {
       fprintf(pcmoutfile,"\n Born Radii of Atoms\n");
       for (i=1; i <= natom; i+= 5)
         fprintf(pcmoutfile,"%-3d %-8.4f   %-3d %-8.4f   %-3d %-8.4f   %-3d %-8.4f %-3d %-8.4f\n",i,solvent.rborn[i],
          i+1,solvent.rborn[i+1],i+2,solvent.rborn[i+2],i+3,solvent.rborn[i+3],i+4,solvent.rborn[i+4]);
   }
//
    if (minim_values.iprint)
       fprintf(pcmoutfile,"\nSolvation Surface Term\n");
    term = 4.0*PI;
    for (i=1; i <= natom; i++)
    {
      if (atom[i].mmx_type != 20)
      {
        ai = solvent.asolv[i];
        ri = solvent.rsolv[i];
        rb = solvent.rborn[i];
        if (rb != 0.0)
        {
            term1 = (ri/rb)*(ri/rb);
            e = ai*term*(ri+probe)*(ri+probe)*term1*term1*term1;
            energies.esolv += e;
            if (minim_values.iprint)
              fprintf(pcmoutfile,"Atom %-3d  %-8.4f\n",i,e);
        }
      }
    }
// polarization term
    f = -cc*(1.0 - 1.0/dwater);
    if (minim_values.iprint)
      fprintf(pcmoutfile,"\nSolvation Polar Term  and Total Energy\n"); 

    for (i=1; i <= natom; i++)
    {
      elocal = 0.0;
      if (fabs(atom[i].charge) < 0.010 && atom[i].mmx_type != 20)
      {
         for(j=i; j <= natom; j++)
         {
           if (atom[i].use || atom[j].use)
           {
              if (fabs(atom[j].charge) < 0.010 && atom[j].mmx_type != 20)
              {
                 xi = atom[i].x;
                 yi = atom[i].y;
                 zi = atom[i].z;

                 xr = xi - atom[j].x;
                 yr = yi - atom[j].y;
                 zr = zi - atom[j].z;
                 rik2 =  xr*xr + yr*yr + zr*zr;
                 if (rik2 < cutoff)
                 {
                    rb2 = solvent.rborn[i]*solvent.rborn[j];
                    fgb = sqrt(rik2 + rb2*exp((-0.25*rik2)/rb2));
                    e = (f*atom[i].charge*atom[j].charge)/fgb;
                    if (i == j)
                    {
                        e *= 0.5;
                        elocal += e;
                        energies.esolv += e;
                    } else
                    {
                        elocal += e;
                        energies.esolv += e;
                    }
                }
              }
           }
         }
         if (minim_values.iprint)
            fprintf(pcmoutfile,"Atom  %-3d   %-8.4f %-8.4f\n",i,elocal,energies.esolv);
      }
    }    
}
// =========================================
void esolv1()
{
    int i,j;
    double e,de,derb,ai,ri,rb;
    double term,probe,term1;
    double dwater, cc,f,cutoff;
    double xi,yi,zi,xr,yr,zr,rik2,r;
    double rb2,fgb,fgb2,expterm;
    double rbi,rbk,dedx,dedy,dedz,drbi,drbk;
    double temp;

    probe = 1.40;
    dwater = 78.3;
    cc = 4.80298*4.80298*14.39418;
    cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
   for (i=1; i <= natom; i++)
   {
       deriv.drb[i] = 0.0;
       deriv.desolv[i][0] = 0.0;
       deriv.desolv[i][1] = 0.0;
       deriv.desolv[i][2] = 0.0;
   }

    born();
    term = 4.0*PI;
    for (i=1; i <= natom; i++)
    {
      if (atom[i].mmx_type != 20)
      {
        ai = solvent.asolv[i];
        ri = solvent.rsolv[i];
        rb = solvent.rborn[i];
        if (rb != 0.0)
        {
            term1 = (ri/rb)*(ri/rb);
            e = ai*term*(ri+probe)*(ri+probe)*term1*term1*term1;
            energies.esolv += e;
            deriv.drb[i] -= 6.0*e/rb;
            temp = deriv.drb[i];
        }
      }
    }
// polarization term
    f = -cc*(1.0 - 1.0/dwater);

    for (i=1; i <= natom; i++)
    {
      if (atom[i].charge != 0.0 && atom[i].mmx_type != 20)
      {         
         xi = atom[i].x;
         yi = atom[i].y;
         zi = atom[i].z;
         rbi = solvent.rborn[i];
         for(j=i; j <= natom; j++)
         {
           if (atom[i].use || atom[j].use)
           {
             if (atom[j].charge != 0.0 && atom[j].mmx_type != 20)
             {
                 
                 xr = xi - atom[j].x;
                 yr = yi - atom[j].y;
                 zr = zi - atom[j].z;
                 rik2 =  xr*xr + yr*yr + zr*zr;
                 if (rik2 < cutoff)
                 {
                    r = sqrt(rik2);
                    rbk = solvent.rborn[j];
                    rb2 = rbi*rbk;
                    expterm = exp(-0.25*rik2/rb2);
                    fgb2 = rik2 + rb2*expterm;
                    fgb = sqrt(fgb2);
                    e = f*atom[i].charge*atom[j].charge/fgb;
                    de = -e*(r-0.25*r*expterm)/fgb2;
                    temp = 0.5 + 0.125*rik2/rb2;
                    derb = -e*expterm*temp/fgb2;
                    if (i == j)
                    {
                        e *= 0.5;
                        energies.esolv += e;
                        drbi = derb*rbk;
                        deriv.drb[i] += drbi;
                        temp = deriv.drb[i];
                    } else
                    {
                        energies.esolv += e;
                        de /= r;
                        dedx = de*xr;
                        dedy = de*yr;
                        dedz = de*zr;
                        deriv.desolv[i][0] += dedx;
                        deriv.desolv[i][1] += dedy;
                        deriv.desolv[i][2] += dedz;
                        deriv.desolv[j][0] -= dedx;
                        deriv.desolv[j][1] -= dedy;
                        deriv.desolv[j][2] -= dedz;
                        
                        drbi = derb*rbk;
                        drbk = derb*rbi;
                        deriv.drb[i] += drbi;
                        deriv.drb[j] += drbk;
                    }
                }
             }
           }
         }
      }
    }
//
    born1();
}
/* ===================================================== */
void born1()
{
    int i,k;
    double xi,yi,zi,cc,rvdw;
    double xr,yr,zr;
    double de;
    double r,r2,r6;
    double p5inv,pip5;
    double gpi,vk,ratio;
    double ccf,cosq,dccf;
    double sinq,term,theta;
    double rb2,ri,rk,sk,sk2;
    double lik,lik2,lik3;
    double uik,uik2,uik3;
    double dlik,duik;
    double t1,t2,t3;
    double dedx,dedy,dedz;
    double temp;

    cc = 4.80298*4.80298*14.39418;
    if (solvent.type == STILL)
    {
        for (i=1; i <= natom; i++)
          skip[i] = 0;
        p5inv = 1.0/solvent.p5;
        pip5 = PI*solvent.p5;
        for (i=1; i <= natom; i++)
        {
          if (atom[i].mmx_type != 20)
          {
            xi = atom[i].x;
            yi = atom[i].y;
            zi = atom[i].z;
            gpi = 2.0*solvent.rborn[i]*solvent.rborn[i]/cc;
            skip[i] = i;
            for (k=0; k < MAXIAT; k++)
            {
               if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
                 skip[atom[i].iat[k]] = i;
            }
            for(k=0; k < attached.n13[i]; k++)
               skip[attached.i13[k][i]] = i;
            for (k = 1; k <= natom; k++)
            {
                if (skip[k] != i  && atom[k].mmx_type != 20)
                {
                    xr = atom[k].x - xi;
                    yr = atom[k].y - yi;
                    zr = atom[k].z - zi;
                    vk = solvent.vsolv[k];
                    r2 = xr*xr + yr*yr + zr*zr;
                    r = sqrt(r2);
                    r6 = r2*r2*r2;
                    rvdw = solvent.rsolv[i] + solvent.rsolv[k];
                    ratio = r2/(rvdw*rvdw);
                    if (ratio > p5inv)
                    {
                        ccf = 1.0;
                        dccf = 0.0;
                    } else
                    {
                        theta = ratio*pip5;
                        cosq = cos(theta);
                        term = 0.5*(1.0-cosq);
                        ccf = term*term;
                        sinq = sin(theta);
                        dccf = 2.0*term*sinq*pip5*ratio;
                    }
                    de = deriv.drb[i]*solvent.p4*gpi*vk*(4.0*ccf-dccf)/r6;
                    temp = deriv.drb[i];
                    dedx = de*xr;
                    dedy = de*yr;
                    dedz = de*zr;
                    deriv.desolv[i][0] += dedx;
                    deriv.desolv[i][1] += dedy;
                    deriv.desolv[i][2] += dedz;
                    deriv.desolv[k][0] -= dedx;
                    deriv.desolv[k][1] -= dedy;
                    deriv.desolv[k][2] -= dedz;
                    
                }
            }
          }
        }
    } else if (solvent.type == HCT)
    {
        for (i=1; i <= natom; i++)
        {
          if (atom[i].mmx_type != 20)
          {
            xi = atom[i].x;
            yi = atom[i].y;
            zi = atom[i].z;
            ri = solvent.rsolv[i] + solvent.doffset;
            rb2 = solvent.rborn[i]*solvent.rborn[i];
            for (k=1; k <= natom; k++)
            {
                if (i != k  && atom[k].mmx_type != 20)
                {
                    xr = atom[k].x - xi;
                    yr = atom[k].y - yi;
                    zr = atom[k].z - zi;
                    r2 = xr*xr + yr*yr + zr*zr;
                    r = sqrt(r2);
                    rk = solvent.rsolv[k]+solvent.doffset;
                    sk = rk * solvent.shct[k];
                    sk2 = sk*sk;
                    if (ri < (r+sk) )
                    {
                        lik = 1.0/ MaxFun(ri,r-sk);
                        uik = 1.0/(r+sk);
                        lik2 = lik*lik;
                        uik2 = uik*uik;
                        lik3 = lik * lik2;
                        uik3 = uik * uik2;
                        dlik = 1.0;
                        if (ri >= (r-sk)) dlik = 0.0;
                        duik = 1.0;
                        t1 = 0.5*lik2 + 0.25*sk2*lik3/r - 0.25*((lik/r)+(lik3*r));
                        t2 = -0.5*uik2 - 0.25*sk2*(uik3/r) + 0.25*((uik/r)+uik3*r);
                        t3 = 0.125*(1.0+sk2/r2)*(lik2-uik2) + 0.25*log(uik/lik)/r2;
                        de = deriv.drb[i] * rb2 * (dlik*t1+duik*t2+t3) / r;
                        dedx = de*xr;
                        dedy = de*yr;
                        dedz = de*zr;
                        deriv.desolv[i][0] += dedx;
                        deriv.desolv[i][1] += dedy;
                        deriv.desolv[i][2] += dedz;
                        deriv.desolv[k][0] -= dedx;
                        deriv.desolv[k][1] -= dedy;
                        deriv.desolv[k][2] -= dedz;
                    }
                }
            }
          }
        }
    }
}
/*    ============================================== */
void esolv2(int ia)
{
    int jj,k;
    double probe,cc,cutoff;
    double fi,fik,de,d2e;
    double d2edx,d2edy,d2edz;
    double xi,yi,zi,xr,yr,zr;
    double r,r2;
    double dwater,rb2;
    double expterm;
    double fgb,fgb2,dfgb,dfgb2,d2fgb;
    double term[3][3];

    probe = 1.40;
    dwater = 78.3;
    cc = 4.80298*4.80298*14.39418;

    if (fabs(atom[ia].charge) < 0.01 || atom[ia].mmx_type == 20)
       return;
        
    cutoff = cutoffs.chrgcut*cutoffs.chrgcut;
    fi = -cc*(1.0-1.0/dwater)*atom[ia].charge;
    xi = atom[ia].x;
    yi = atom[ia].y;
    zi = atom[ia].z;

    for (k = 1; k <= natom; k++)
    {
        if (ia != k && atom[k].charge != 0.0 && atom[k].mmx_type != 20)
        {
            xr = xi - atom[k].x;
            yr = yi - atom[k].y;
            zr = zi - atom[k].z;
            r2 = xr*xr + yr*yr + zr*zr;
            if (r2 < cutoff)
            {
                r = sqrt(r2);
                fik = fi*atom[k].charge;
                rb2 = solvent.rborn[ia]*solvent.rborn[k];
                expterm = exp(-0.25*r2/rb2);
                fgb2 = r2 + rb2*expterm;
                fgb = sqrt(fgb2);
                dfgb = (1.0-0.25*expterm)*r/fgb;
                dfgb2 = dfgb*dfgb;
                d2fgb = -dfgb2/fgb + dfgb/r + 0.125*(r2/rb2)*expterm/fgb;
                de = -fik*dfgb/fgb2;
                d2e = -fik*(d2fgb-2.0*dfgb2/fgb)/fgb2;
                de /= r;
                d2e = (d2e-de)/r;
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
                for (jj=0; jj < 3; jj++)
                {
                    hess.hessx[ia][jj] += term[jj][0];
                    hess.hessy[ia][jj] += term[jj][1];
                    hess.hessz[ia][jj] += term[jj][2];
                    hess.hessx[k][jj] -= term[jj][0];
                    hess.hessy[k][jj] -= term[jj][1];
                    hess.hessz[k][jj] -= term[jj][2];
                }
            }
        }
    }
}

                
                
                
            
    
