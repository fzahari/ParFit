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

// mode
#define NONE      0
#define NEWTON    1
#define TNCG      2
#define DTNCG     3
#define AUTO      4
// method
#define DIAG      5
#define SSOR      6
#define ICCG      7
#define BLOCK     8


#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "derivs.h"
#include "utility.h"
#include "minvar.h"
#include "solve.h"

int stable;
double *f, *f_diag;
int *c_index, *c_value;

void solve(int mode,int method,int negtest,int nvar,double *p,double *x,
           double *g, double *h, int *h_init, int *h_stop, int *h_index,
           double *h_diag,int *icycle,int *iter_cg,int *fg_call,double (*fgvalue) ())
{
    int i,j,k, iter,maxiter,maxhess;
    double *m,  *d, *si, *r,  *q;  // 
    double alpha,beta,delta,sigma;
    double f_sigma;
    double *x_sigma;
    double *g_sigma;
    double g_norm,g_rms,eps,converge;
    double hj,gg,dq,rr,dd,rs,rs_new,r_norm;
    int maxvar = 3*natom;
   if (natom < 300)
      maxhess = (3*natom*(3*natom-1))/2;
   else if (natom < 800)
      maxhess = (3*natom*(3*natom-1))/4;
   else
      maxhess = (3*natom*(3*natom-1))/20;

    m = dvector(0,maxvar);
    r = dvector(0,maxvar);
    si = dvector(0,maxvar);
    d = dvector(0,maxvar);
    q = dvector(0,maxvar);
    x_sigma = dvector(0,maxvar);
    g_sigma = dvector(0,maxvar);
    f_diag = dvector(0,maxvar);
    f = dvector(0,maxhess);
    c_index = ivector(0,maxhess);
    c_value = ivector(0,maxhess);

      if (mode != DTNCG &&  method != NONE)
      {
         for (i=0; i < nvar; i++)
            m[i] = 1.0 / sqrt(fabs(h_diag[i]));
            
         for (i=0; i < nvar; i++)
         {
            g[i] = g[i] * m[i];
            h_diag[i] *=  m[i] * m[i];
            for (j = h_init[i]; j < h_stop[i]; j++)
            {
               k = h_index[j];
               h[j] *= m[i] * m[k];
            }
         }
      }

    iter = 0;
    maxiter = 0;
    eps = 0.0;
    gg = 0.0;
    for (i=0; i < nvar; i++)
    {
        p[i] = 0.0;
        r[i] = -g[i];
        gg += g[i]*g[i];
    }
    g_norm = sqrt(gg);
    precond(method, iter,nvar, si, r, h, h_init, h_stop, h_index, h_diag);
    rs = 0.0;
    for (i=0; i < nvar; i++)
    {
        d[i] = si[i];
        rs += r[i]*si[i];
    }
    
    if (mode == NEWTON)
    {
        eps = 1.0e-10;
        maxiter = nvar;
    } else if (mode == TNCG || mode == DTNCG)
    {
        delta = 1.0;
        eps = delta/(double) *icycle;
        g_rms = g_norm/sqrt((double) nvar);
        eps = MinFun(eps,g_rms);
        converge = 1.0;
        maxiter = (int) (10.0*sqrt((double)nvar));
    }
    iter = 1;
    
// evaluate the matrix vector product
    while (TRUE)
    {
        if (mode == TNCG || mode == NEWTON)
        {
            for (i=0; i < nvar; i++)
               q[i] = 0.0;
            for (i=0; i < nvar; i++)
            {
               q[i] += h_diag[i]*d[i];
               for(j= h_init[i]; j < h_stop[i]; j++)
               {
                   k = h_index[j];
                   hj = h[j];
                   q[i] += hj*d[k];
                   q[k] += hj*d[i];
               }
            }
        }else if (mode == DTNCG)
        {
            dd = 0.0;
            for (i=0; i < nvar; i++)
               dd += d[i]*d[i];
               
            sigma = 1.0e-7 / sqrt(dd);
            
            for (i=0; i < nvar; i++)
               x_sigma[i] = x[i] + sigma*d[i];
               
            fg_call = fg_call + 1;
            f_sigma = fgvalue (x_sigma,g_sigma);
            
            for (i=0; i < nvar; i++)
               q[i] = (g_sigma[i]-g[i]) / sigma;
        }

// test for direction of negative curvature
         dq = 0.0;
         for (i=0; i < nvar; i++)
            dq += d[i]*q[i];
            
         if (negtest)
         {
            if (dq <= 0.0)
            {
               if (iter == 1)
               {
                  for (i=0; i < nvar; i++)
                     p[i] = d[i];
               }
               goto L_10;
            }
         }

// test truncated Newton termination criterion
        alpha = rs / dq;
         rr = 0.0;
         for (i=0; i < nvar; i++)
         {
            p[i] += alpha*d[i];
            r[i] -= alpha*q[i];
            rr += r[i]*r[i];
         }
         r_norm = sqrt(rr);
         if (r_norm/g_norm <= eps)
            goto L_10;

// precondition
          precond(method,iter,nvar, si, r,h,h_init,h_stop,
                  h_index,h_diag);

// update direction
          rs_new = 0.0;
          for (i=0; i < nvar; i++)
            rs_new += r[i]*si[i];
            
          beta = rs_new/rs;
          rs = rs_new;
          
          for (i=0; i < nvar; i++)
            d[i] = si[i] + beta*d[i];
            
// check for overlimit next iteration
          if (iter >= maxiter)
             goto L_10;

          iter++;
    }

// termination

L_10:
      if (mode != DTNCG && method != NONE)
      {
         for (i=0; i < nvar; i++)
         {
            p[i] *= m[i];
            g[i] /= m[i];
         }
      }
      iter_cg += iter;
    free_dvector( m ,0,maxvar);
    free_dvector( r ,0,maxvar);
    free_dvector( si ,0,maxvar);
    free_dvector( d ,0,maxvar);
    free_dvector( q ,0,maxvar);
    free_dvector( x_sigma ,0,maxvar);
    free_dvector( g_sigma ,0,maxvar);
    free_dvector( f_diag ,0,maxvar);
    free_dvector( f ,0,maxhess);
    free_ivector( c_index ,0,maxhess);
    free_ivector( c_value ,0,maxhess);
}

void precond(int method, int iter, int nvar, double *s, double *r,
              double *h, int *h_init, int *h_stop,int *h_index,
              double *h_diag)
{
    int i,j,k,ii,kk,iii,kkk,nblock,ix,iy,iz,icount;
    int *c_init, *c_stop;                // c_init[maxvar],c_stop[maxvar];
    double a[6],b[3];
    double omega,factor,*diag;            // diag[maxvar];
    double maxalpha, alpha, f_i, f_k;
    int  maxvar;

   maxvar = 3*natom;
   c_init = ivector(0,maxvar);
   c_stop = ivector(0,maxvar);
   diag = dvector(0,maxvar);

   
// no preconditioning
    if (method == NONE)
    {
        for (i=0; i < nvar; i++)
           s[i] = r[i];
    }
// diagonal
    if (method == DIAG)
    {
        for (i=0; i < nvar; i++)
           s[i] = r[i]/ fabs(h_diag[i]);
    }
// block diagonal
    if (method == BLOCK)
    {
        nblock = 3;
        for (i=1; i <= natom; i++)
        {
            iz = (i-1)*3+2;
            iy = (i-1)*3+1;
            ix = (i-1)*3;
            a[0] = h_diag[ix];
            if (h_index[h_init[ix]] == iy) 
               a[1] = h[h_init[ix]];
            else
               a[1] = 0.0;

            if (h_index[h_init[ix]+1] == iz)
               a[2] = h[h_init[ix]+1];
            else
               a[2] = 0.0;

            a[3] = h_diag[iy];
            if (h_index[h_init[iy]] == iz)
               a[4] = h[h_init[iy]];
            else
               a[4] = 0.0;

            a[5] = h_diag[iz];
            b[0] = r[ix];
            b[1] = r[iy];
            b[2] = r[iz];
            cholesky (nblock,a,b);
            s[ix] = b[0];
            s[iy] = b[1];
            s[iz] = b[2];
        }
    }

// SSOR
      if (method == SSOR)
      {
         omega = 1.0;
         factor = 2.0 - omega;
         for (i=0; i < nvar; i++)
         {
            s[i] = r[i] * factor;
            diag[i] = h_diag[i] / omega;
         }
         for (i=0; i < nvar; i++)
         {
            s[i] = s[i] / diag[i];
            for (j=h_init[i]; j< h_stop[i]; j++)
            {
               k = h_index[j];
               s[k] = s[k] - h[j]*s[i];
            }
         }
         for (i=nvar-1; i >= 0; i--)
         {
            s[i] = s[i] * diag[i];
            for (j=h_init[i]; j< h_stop[i]; j++)
            {
               k = h_index[j];
               s[i] = s[i] - h[j]*s[k];
            }
            s[i] = s[i] / diag[i];
         }
      }
// incomplete cholesky preconditioning
      if (method == ICCG && iter == 0)
      {
         column (nvar,h_init,h_stop,h_index,
                     c_init,c_stop,c_index,c_value);
         stable = TRUE;
         icount = 0;
         maxalpha = 2.1;
         alpha = -0.001;
L_10:
         if (alpha <= 0.0)
            alpha = alpha + 0.001;
         else
            alpha = 2.0 * alpha;

         if (alpha > maxalpha)
         {
            stable = FALSE;
            fprintf(pcmoutfile,"Precond - Incomplete Cholesky is Unstable\n");
         } else
         {
            factor = 1.0 + alpha;
            for (i=0; i < nvar; i++)
            {
               f_diag[i] = factor * h_diag[i];
               for(j = c_init[i]; j<= c_stop[i]; j++)
               {
                  k = c_index[j];
                  f_i = f[c_value[j]];
                  f_diag[i] = f_diag[i] - (f_i*f_i*f_diag[k]);                  
                  icount++;
               }
               if (f_diag[i] <= 0.0)  goto L_10;
               if (f_diag[i] < 1.0e-7)  f_diag[i] = 1.0e-7;
               f_diag[i] = 1.0 / f_diag[i];
               for( j = h_init[i]; j < h_stop[i]; j++)
               {
                  k = h_index[j];
                  f[j] = h[j];
                  ii = c_init[i];
                  kk = c_init[k];
                 while (ii <= c_stop[i] && kk <= c_stop[k])
                  {
                     iii = c_index[ii];
                     kkk = c_index[kk];
                     if (iii < kkk)
                        ii = ii + 1;
                     else if (kkk < iii) 
                        kk = kk + 1;
                     else
                     {
                        f_i = f[c_value[ii]];
                        f_k = f[c_value[kk]];
                        f[j] = f[j] - (f_i*f_k*f_diag[iii]);
                        ii = ii + 1;
                        kk = kk + 1;
                        icount = icount + 1;
                     }
                  }
               }
            }
         }
      }
//     incomplete cholesky preconditioning  (solution phase)

      if (method == ICCG)
      {
         if (stable)
         {
            for (i=0; i < nvar; i++)
               s[i] = r[i];
  
            for (i=0; i < nvar; i++)
            {
               s[i] = s[i] * f_diag[i];
              for( j = h_init[i]; j < h_stop[i]; j++)
               {
                  k = h_index[j];
                  s[k] = s[k] - f[j]*s[i];
               }
            }
            for (i = nvar-1; i >= 0; i--)
            {
               s[i] = s[i] / f_diag[i];
               for( j = h_init[i]; j < h_stop[i]; j++)
               {
                  k = h_index[j];
                  s[i] = s[i] - f[j]*s[k];
               }
               s[i] = s[i] * f_diag[i];
            }
         } else
         {
            for (i=0; i < nvar; i++)
               s[i] = r[i] / fabs(h_diag[i]);
         }
      }
   free_ivector(c_init ,0,maxvar);
   free_ivector(c_stop ,0,maxvar);
   free_dvector(diag ,0,maxvar);
}

void cholesky(int nvar, double *a, double *b)
{
    int i,j,k,ii,ij,ik,im,jk,jm,ki,kk;
    double r,sa,t;

      ii = 0;
      for (i=0; i < nvar; i++)
      {
         im = i - 1;
         if (i != 0)
         {
            ij = i;
            for( j = 0; j <  im; j++)
            {
               r = a[ij];
               if (j != 0)
               {
                  ik = i;
                  jk = j;
                  jm = j - 1;
                  for ( k = 0; k < jm; k++)
                  {
                     r = r - a[ik]*a[jk];
                     ik = nvar-1 - k + ik;
                     jk = nvar-1 - k + jk;
                  }
               }
               a[ij] = r;
               ij = nvar-1 - j + ij;
            }
         }
         r = a[ii];
         if (i != 0)
         {
            kk = 0;
            ik = i;
            for (k = 0; k < im; k++)
            {
               sa = a[ik];
               t = sa * a[kk];
               a[ik] = t;
               r = r - sa*t;
               ik = nvar-1 - k + ik;
               kk = nvar-1 - k +1  + kk;
            }
         }
         a[ii] = 1.0/ r;
         ii = nvar-1 - i + 1 + ii;
      }

//     solve linear equations; first solve Ly = b for y

      for (i=0; i < nvar; i++)
      {
         if (i != 0)
         {
            ik = i;
            im = i - 1;
            r = b[i];
            for (k = 0; k < im; k++)
            {
               r = r - b[k]*a[ik];
               ik = nvar-1 - k + ik;
            }
            b[i] = r;
         }
      }

//     now, solve (D)(L transpose)(x) = y for x

      ii = (nvar*(nvar+1)/2)-1;
      for ( j = 0; j < nvar; j++)
      {
         i = nvar-1  - j;
         r = b[i] * a[ii];
         if (j != 0)
         {
            im = i + 1;
            ki = ii + 1;
            for( k = im; k < nvar; k++)
            {
               r = r - a[ki]*b[k];
               ki++;
            }
         }
         b[i] = r;
         ii = ii - j - 2;
      }
}    
