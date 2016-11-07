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
#include "derivs.h"
#include "utility.h"
#include "minvar.h"
#include "search.h"

#define  Success     0
#define  ReSearch    1
#define  WideAngle   2
#define  BadIntpln   3
#define  IntplnErr   4
#define  blank       5
#define  Failure     6

// #define min(a,b)  (((a) < (b) ? (a) : (b))

void search(int nvar, double *f, double *g, double *x, double *p,
double f_move,double *angle,int *ncalls,double (*fgvalue)(),int *status)
{
  int i, intpln;
  int maxvar;
  double *s; //   s[maxvar];
  double s_norm,g_norm,cosang;
  double step,parab,cube,cubstp,sss,ttt;
  double f_0,f_1,f_a,f_b,f_c,sg_0,sg_1,sg_a,sg_b,sg_c;
  
  maxvar = 3*natom;
  s = dvector(0,maxvar);
  
  if (fabs(minvar.cappa) < 0.0001)   minvar.cappa = 0.1;
  if (minvar.stpmin == 0.0)          minvar.stpmin = 1.0e-20;
  if (fabs(minvar.stpmax) < 0.0001)  minvar.stpmax = 2.0;
  if (fabs(minvar.angmax) < 0.1)     minvar.angmax = 180.0;
  if (minvar.intmax == 0)            minvar.intmax = 5;
       
    for (i=0; i < nvar; i++)
      s[i] = p[i];

     g_norm = 0.0;
     s_norm = 0.0;
     for (i=0; i < nvar; i++)
     {
         g_norm += g[i]*g[i];
         s_norm += s[i]*s[i];
     }
     g_norm = sqrt(g_norm);
     s_norm = sqrt(s_norm);
     if (isnan(g_norm) != 0)
     { 
// Hay commented out       printf("Found nan in search\n");
       *status = Failure;
       free_dvector(s,0,maxvar);
       return;
     }

     f_0 = *f;
     sg_0 = 0.0;
     for (i=0; i < nvar; i++)
     {
         s[i] = s[i] / s_norm;
         sg_0 += s[i]*g[i];
     }
      cosang = -sg_0 / g_norm;
      if (cosang < -1.00)
         cosang = -1.0;
      if (cosang > 1.00)
         cosang = 1.0;
      *angle = radian * acos(cosang);
      if (*angle > minvar.angmax)
      {
         *status = WideAngle;
         free_dvector(s,0,maxvar);
         return;
      }
      step = 2.0 * fabs(f_move/sg_0);
      step = MinFun(step,s_norm);
      if (step > minvar.stpmax)  step = minvar.stpmax;
      if (step < minvar.stpmin)  step = minvar.stpmin;

L_10:
      intpln = 0;
      f_b = f_0;
      sg_b = sg_0;

L_20:
      f_a = f_b;
      sg_a = sg_b;
      for (i=0; i < nvar; i++)
         x[i] = x[i] + step*s[i];
         
//  get new function and gradient

      ncalls++;
      f_b = fgvalue(x,g);
      sg_b = 0.0;

      for (i=0; i < nvar; i++)
         sg_b += s[i]*g[i];

      if (fabs(sg_b/sg_0) <= minvar.cappa && f_b < f_a)
      {
         *f = f_b;
         *status =  Success;
         free_dvector(s,0,maxvar);
         return;
      }
      
      if (sg_a*sg_b < 0.0 || f_a < f_b) goto L_30;
      if ( fabs(sg_a-sg_b) < 0.0000001)
      {
         *f = f_b;
         *status =  Success;
         free_dvector(s,0,maxvar);
         return;
      }
         step = 2.0 * step;
         if (sg_b > sg_a)
         {
            parab = (f_a-f_b) / (sg_b-sg_a);
            if (parab > 2.0*step)  parab = 2.0 * step;
            if (parab < 0.5*step)  parab = 0.5 * step;
            step = parab;
         }
         if (step > minvar.stpmax)  step = minvar.stpmax;
         goto L_20;
 
L_30:
      intpln++;
      sss = 3.0*(f_b-f_a)/step - sg_a - sg_b;
      ttt = sss*sss - sg_a*sg_b;
      if (ttt < 0.0)
      {
         *f = f_b;
         *status = IntplnErr;
         free_dvector(s,0,maxvar);
         return;
      }
      ttt = sqrt(ttt);
      cube = step * (sg_b+ttt+sss)/(sg_b-sg_a+2.0*ttt);
      if (cube < 0.0 || cube > step)
      {
         *f = f_b;
         *status = IntplnErr;
          free_dvector(s,0,maxvar);
        return;
      }
      for (i=0; i < nvar; i++)
         x[i] -= cube*s[i];
// get new function and gradient and test
      ncalls++;
      f_c = fgvalue(x,g);
      sg_c = 0.0;
      for (i=0; i < nvar; i++)
          sg_c += s[i]*g[i];

      if (fabs(sg_c/sg_0) <= minvar.cappa)
      {
         *f = f_c;
         *status = Success;
         free_dvector(s,0,maxvar);
         return;
      }
       if (f_c <= f_a || f_c <= f_b)
       {
         cubstp = MinFun(fabs(cube),fabs(step-cube));
         if (cubstp >= minvar.stpmin && intpln < minvar.intmax)
         {
            if (sg_a*sg_b < 0.0)
            {
               if (sg_a*sg_c < 0.0)
               {
                  f_b = f_c;
                  sg_b = sg_c;
                  step = step - cube;
               }else
               {
                  f_a = f_c;
                  sg_a = sg_c;
                  step = cube;
                  for(i=0; i < nvar; i++)
                     x[i] = x[i] + cube*s[i];
                  
               }
            } else
            {
               if (sg_a*sg_c < 0.0 || f_a <= f_c)
               {
                  f_b = f_c;
                  sg_b = sg_c;
                  step = step - cube;
               } else
               {
                  f_a = f_c;
                  sg_a = sg_c;
                  step = cube;
                  for(i=0; i < nvar; i++)
                     x[i] = x[i] + cube*s[i];
                 
               }
            }
            goto L_30;
         }
       }

//  failed interpolation
       if (f_a < f_b && f_a < f_c)
       {
          f_1 = f_a;
          sg_1 = sg_a;
          for (i=0; i < nvar; i++)
             x[i] += (cube-step)*s[i];
       } else if (f_b < f_a && f_b < f_c)
       {
          f_1 = f_b;
          sg_1 = sg_b;
          for (i=0; i < nvar; i++)
             x[i] += cube*s[i];
       } else
       {
          f_1 = f_c;
          sg_1 = sg_c;
       }
// try to restart
       if (f_1 > f_0)
       {
           ncalls++;
           *f = fgvalue(x,g);
           *status = IntplnErr;
           free_dvector(s,0,maxvar);
           return;
       }
       f_0 = f_1;
       sg_0 = sg_1;
       if (sg_1 > 0.0)
       {
           for (i=0; i < nvar; i++)
              s[i] = -s[i];
            sg_0 = -sg_1;
       }
       if (cube > (step-cube) )
          step = cube/10.0;
       else
          step = (step-cube)/10.0;

       if (step < minvar.stpmin) step = minvar.stpmin;

       if (*status == ReSearch)
       {
           ncalls++;
           *f = fgvalue(x,g);
           *status = BadIntpln;
            free_dvector(s,0,maxvar);
           return;
       } else
       {
           *status = ReSearch;
           goto L_10;
       }     
}

