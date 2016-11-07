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
#include "diagc.h"

void tred2(double **a,int n, double *d, double *e)
{

    int i,j,k,l, nn;
    double scale, hh, h, g, f;

    nn = n-1;
    for (i= nn; i >= 1; i--)
    {
        l = i - 1;
        h = scale = 0.0;
        if (l > 0)
        {
            for (k=0; k <= l; k++)
                scale += fabs(a[i][k]);
            if (fabs(scale) < 0.001)
                e[i] = a[i][l];
            else
            {
                for (k=0; k <= l; k++)
                {
                    a[i][k] /= scale;
                    h += a[i][k]*a[i][k];
                }
                f = a[i][l];
                g = f>0 ? -sqrt(h) : sqrt(h);
                e[i] = scale*g;
                h -= f*g;
                a[i][l] = f-g;
                f = 0.0;
                for (j=0; j <= l; j++)
                {
                    a[j][i] = a[i][j]/h;
                    g = 0.0;
                    for (k=0; k <= j; k++)
                        g += a[j][k]*a[i][k];
                    for (k = j+1; k <= l; k++)
                        g += a[k][j]*a[i][k];
                    e[j] = g/h;
                    f += e[j]*a[i][j];
                }
                hh = f/(h+h);
                for (j=0; j <= l; j++)
                {
                    f = a[i][j];
                    e[j] = g = e[j] - hh*f;
                    for (k=0; k <= j; k++)
                      a[j][k] -= (f*e[k]+g*a[i][k]);
                }
            }
        }else
           e[i] = a[i][l];
        d[i] = h;
    }
    d[0] = 0.0;
    e[0] = 0.0;
    for (i=0; i < n; i++)
    {
        l = i - 1;
        if (i > 0)
        {
            for (j=0; j <= l; j++)
            {
                g = 0.0;
                for (k=0; k <= l; k++)
                    g += a[i][k]*a[k][j];
                for (k=0; k <= l; k++)
                    a[k][j] -= g*a[k][i];
            }
        }
        d[i] = a[i][i];
        a[i][i] = 1.0;
        for (j=0; j <= l; j++)
        {
            a[j][i] = a[i][j] = 0.0;
        }
    }
}
/* =========================================== */
void tqli(double *d, double *e, int n, double **z)
{
    int m,l,iter,i,k,j;
    double s,r,p,g,f,dd,c,b;

    for (i=1; i < n; i++)
       e[i-1] = e[i];
    e[n-1] = 0.0;

    for (l=0; l < n; l++)
    {
        iter = 0;
        do
        {
            for (m=l; m < n-1; m++)
            {
                dd = fabs(d[m]) + fabs(d[m+1]);
                if (fabs(e[m])+dd == dd) break;
            }
            if (m!= l)
            {
                if (iter++ == 30)
                {
                    fprintf(pcmoutfile,"Error in tqli\n");
                    message_alert("Error in Tqli","Error");
                    return;
                }
                g = (d[l+1]-d[l])/(2.0*e[l]);
                r = sqrt((g*g)+1.0);
                g = d[m]-d[l]+e[l]/(g+sign(r,g));
                s = c = 1.0;
                p = 0.0;
                for (i=m-1; i >= l; i--)
                {
                    f = s*e[i];
                    b = c*e[i];
                    if (fabs(f) >= fabs(g))
                    {
                        c = g/f;
                        r = sqrt((c*c)+1.0);
                        e[i+1] = f*r;
                        c *= (s=1.0/r);
                    } else
                    {
                        s = f/g;
                        r = sqrt((s*s)+1.0);
                        e[i+1] = g*r;
                        s *= (c=1.0/r);
                    }
                    g = d[i+1]-p;
                    r = (d[i]-g)*s+2.0*c*b;
                    p = s*r;
                    d[i+1] = g+p;
                    g = c*r-b;
                    for (k=0; k < n; k++)
                    {
                        f = z[k][i+1];
                        z[k][i+1] = s*z[k][i]+c*f;
                        z[k][i] = c*z[k][i]-s*f;
                    }
                }
                d[l] = d[l]-p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m!= l);
    }
// sort d and z
    for (i=0; i < n; i++)
    {
        p  = d[k=i];
        for (j=i+1; j < n; j++)
            if (d[j] >= p) p = d[k=j];
        if (k != i)
        {
            d[k] = d[i];
            d[i] = p;
            for (j=0; j < n; j++)
            {
                p = z[j][i];
                z[j][i] = z[j][k];
                z[j][k] = p;
            }
        }
    }  
}
