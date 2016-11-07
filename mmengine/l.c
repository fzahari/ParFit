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
#include "utility.h"
#include "bonds_ff.h"
#include "torsions.h"
#include "pimatx.h"
#include "pistuf.h"
#include "pidata.h"
#include "pibonds.h"
#include "pisave.h"
#include "l.h"
#include "minim_values.h"
     
void  otable()
{
        int i, ik, k, lt, ng, no, nread;
        float a1, a2, b1, b2, c1, c2, rik, v2, v3, v4;
        double a, a3, a4, a5, a6, a7, aa, b, b3, b4, b5, b6, b7, bb, ddldl, 
         dpipi, dsgsb, dsgsg, exa, exb, onepi, onesg, p, p10, p11, p2, p3, 
         p4, p5, p6, p7, p8, p9, s0, s1, s1a, s1b, s1c, s2, s2a, s2b, 
         s2c, s5, sadsg, sasb, spi, ssg, sx, sy, t, t1, t2, twopi, twosg, 
         v1, x1, x2, xs, xx, xy, xy2, xz2, ys, yx, yy, yz, yz2, zeta;

        no = pistuf.norb - 1;
        nread = 0;
	sadsg = 0.0;
        for( i = 0; i < no; i++ )
        {
                penint.pen[i][i] = 0;
                ik = i + 1;
                for( k = ik; k < pistuf.norb; k++ )
                {
                        lt = 1;
                        zeta = (pimatx.zt[i] + pimatx.zt[k])/2;
                        rik = pibonds.riks[nread];
                        a1 = pibonds.a1s[nread];
                        b1 = pibonds.b1s[nread];
                        c1 = pibonds.c1s[nread];
                        a2 = pibonds.a2s[nread];
                        b2 = pibonds.b2s[nread];
                        c2 = pibonds.c2s[nread];
                        nread += 1;
                        p = zeta*rik;
                        t = (pimatx.zt[i] - pimatx.zt[k])/(pimatx.zt[i] + pimatx.zt[k]);
                        p2 = p*p;
                        p3 = p2*p;
                        p4 = p2*p2;
                        p5 = p2*p3;
                        p6 = p3*p3;
                        p7 = p4*p3;
                        p8 = p4*p4;
                        p9 = p5*p4;
                        p10 = p5*p5;
                        p11 = p5*p6;
                        if( t )
                        {
                                if( fabs( t ) < .04 )
                                {
                                        lt = 2;
                                }
                                else
                                {
                                        v1 = (t + 1./t)/2;
                                        v2 = v1*v1;
                                        v3 = v1*v2;
                                        v4 = v2*v2;
                                        a = (1. + t)*p;
                                        aa = a*a;
                                        a3 = a*aa;
                                        a4 = aa*aa;
                                        a5 = aa*a3;
                                        a6 = a3*a3;
                                        a7 = a3*a4;
                                        b = (1. - t)*p;
                                        bb = b*b;
                                        b3 = b*bb;
                                        b4 = bb*bb;
                                        b5 = bb*b3;
                                        b6 = b3*b3;
                                        b7 = b3*b4;
                                        exa = exp( -a );
                                        exb = exp( -b );
                                        t1 = sqrt((1. - t*t));
                                        t1 = 1./(t1*t*(p*p*p));
                                        x1 = (1. - v1)*(1. - v1);
                                        x2 = (1. + v1)*(1. + v1);
                                        s1 = (24.*x2*(1. + a) + 12.*aa*(1. + v1) + 2.*
                                         a3)*exa;
                                        s2 = (24.*x1*(1. + b) + 12.*bb*(1. - v1) + 2.*
                                         b3)*exb;
                                        spi = fabs( t1*(x2*s2 - x1*s1) );
                                        s1 = (48.*x2*(1. + a + .5*aa) + (10. + 12.*v1)*
                                         a3 + 2.*a4)*exa;
                                        s2 = (48.*x1*(1. + b + .5*bb) + (10. - 12.*v1)*
                                         b3 + 2.*b4)*exb;
                                        ssg = fabs( t1*(x2*s2 - x1*s1) );
                                        s1 = (.0625*(8. - v1 - 27.*v2 - 30.*v3 - 10.*v4) + 
                                         (a/32.)*(11. - 19.*v1 - 44.*v2 - 20.*v3) + (aa/
                                         16.)*(1. - 5.*v1 - 4.*v2) - (v1*a3/24.))*exp( -2.*
                                         a );
                                        s2 = (.0625*(8. + v1 - 27.*v2 + 30.*v3 - 10.*v4) + 
                                         (b/32.)*(11. + 19.*v1 - 44.*v2 + 20.*v3) + (bb/
                                         16.)*(1. + 5.*v1 - 4.*v2) + (v1*b3/24.))*exp( -2.*
                                         b );
                                        sasb = (zeta/p)*(1. - (1. - v1)*x1*s1 - (1. + 
                                         v1)*x2*s2);
                                        ng = 0;
                                        while( True )
                                        {
                                                s0 = zeta/((1. - t)*(1. - t)*p3);
                                                t1 = (1. - v1)*(1. - v1)*(1. - v1)*(1. - v1);
                                                t2 = (1. + v1)*(1. + v1)*(1. + v1);
                                                s1 = (.03125*(21. + 44.*v1 + 35.*v2 + 10.*
                                                 v3)*(1. + 2.*a) + (a*a/24.)*(29. + 56.*v1 + 
                                                 40.*v2 + 10.*v3) + ((a*a*a)/12.)*(8. + 12.*
                                                 v1 + 5.*v2) + ((a*a*a*a)/24.)*(5. + 4.*v1) + 
                                                 ((a*a*a*a*a)/36.))*exp( -2.*a );
                                                s2a = 0.03125*(11. + 7.*v1 - 39.*v2 + 35.*
                                                 v3 - 10.*v4)*(1. + 2.*b) + (b*b/24.)*(14. + 
                                                 13.*v1 - 51.*v2 + 40.*v3 - 10.*v4);
                                                s2b = ((b*b*b)/12.)*(3. + 6.*v1 - 12.*v2 + 
                                                 5.*v3) + ((b*b*b*b)/24.)*(1. + 5.*v1 - 4.*
                                                 v2) + (v1*(b*b*b*b*b)/36.);
                                                s2 = (s2a + s2b)*exp( -2.*b );
                                                s5 = s0*(1. - t1*s1 - t2*s2);
                                                ng += 1;
                                                v1 = -v1;
                                                v3 = -v3;
                                                t = -t;
                                                p10 = a;
                                                a = b;
                                                b = p10;
                                                if (ng == 1)
                                                 sadsg = s5;
                                                else if( ng == 2 )
                                                {
                                                   dsgsb = s5;
                                                   break;
                                                }
                                        }
                                        s0 = 6.*zeta/((1. + t)*(1. + t)*(1. - t)*(1. - t)*p5);
                                        t1 = (1. - v1)*(1. - v1)*(1. - v1)*(1. - v1);
                                        t2 = (1. + v1)*(1. + v1)*(1. + v1)*(1. + v1);
                                        exa = exp( -2.*a );
                                        exb = exp( -2.*b );
                                        s1a = .03125*(16. + 29.*v1 + 20.*v2 + 5.*v3)*(1. + 
                                         2.*a) + (aa/144.)*(139. + 246.*v1 + 165.*v2 + 
                                         40.*v3);
                                        s1c = (a3/72.)*(43. + 72.*v1 + 45.*v2 + 10.*v3) + 
                                         (a4/216.)*(57. + 88.*v1 + 50.*v2 + 10.*v3);
                                        s1b = (a5/432.)*(39. + 52.*v1 + 20.*v2) + (a6/
                                         216.)*(5. + 4.*v1) + (a7/324.);
                                        s1 = (s1a + s1b + s1c)*exa;
                                        s2a = .03125*(16. - 29.*v1 + 20.*v2 - 5.*v3)*(1. + 
                                         2.*b) + (bb/144.)*(139. - 246.*v1 + 165.*v2 - 
                                         40.*v3);
                                        s2c = (b3/72.)*(43. - 72.*v1 + 45.*v2 - 10.*v3) + 
                                         (b4/216.)*(57. - 88.*v1 + 50.*v2 - 10.*v3);
                                        s2b = (b5/432.)*(39. - 52.*v1 + 20.*v2) + (b6/
                                         216.)*(5. - 4.*v1) + (b7/324.);
                                        s2 = (s2a + s2b + s2c)*exb;
                                        dsgsg = s0*(1. - t1*s1 - t2*s2);
                                        s0 = 2.*s0/3.;
                                        s1a = .03125*(16. + 29.*v1 + 20.*v2 + 5.*v3)*(1. + 
                                         a + a) + (aa/96.)*(91. + 159.*v1 + 105.*v2 + 
                                         25.*v3) + (a3/48.)*(27. + 43.*v1 + 25.*v2 + 5.*
                                         v3);
                                        s1b = (a4/48.)*(11. + 14.*v1 + 5.*v2) + (a5/288.)*
                                         (17. + 12.*v1) + a6/144.;
                                        s2a = .03125*(16. - 29.*v1 + 20.*v2 - 5.*v3)*(1. + 
                                         b + b) + (bb/96.)*(91. - 159.*v1 + 105.*v2 - 
                                         25.*v3) + (b3/48.)*(27. - 43.*v1 + 25.*v2 - 5.*
                                         v3);
                                        s2b = (b4/48.)*(11. - 14.*v1 + 5.*v2) + (b5/288.)*
                                         (17. - 12.*v1) + b6/144.;
                                        s1 = (s1a + s1b)*exa;
                                        s2 = (s2a + s2b)*exb;
                                        dpipi = s0*(1. - t1*s1 - t2*s2);
                                        s0 = zeta/((1. + t)*(1. + t)*(1. - t)*(1. - t)* p5);
                                        s1 = (.03125*(16. + 29.*v1 + 20.*v2 + 5.*v3)*(1. + 
                                         2.*a) + (aa/48.)*(43. + 72.*v1 + 45.*v2 + 10.*
                                         v3) + (a3/24.)*(11. + 14.*v1 + 5.*v2) + (a4/24.)*
                                         (3. + 2.*v1) + (a5/72.))*exa;
                                        s2 = (.03125*(16. - 29.*v1 + 20.*v2 - 5.*v3)*(1. + 
                                         2.*b) + (bb/48.)*(43. - 72.*v1 + 45.*v2 - 10.*
                                         v3) + (b3/24.)*(11. - 14.*v1 + 5.*v2) + (b4/24.)*
                                         (3. - 2.*v1) + (b5/72.))*exb;
                                        ddldl = s0*(1. - t1*s1 - t2*s2);
                                        goto L_90;
                                }
                        }
                        s1 = exp( -p );
                        spi = fabs( (1. + p + .4*p2 + p3/15.)*s1 );
                        ssg = fabs( (-1. - p - .2*p2 + 2.*p3/15. + p4/15.)*s1 );
                        s1 = 1. + 419.*p/256. + 163.*p2/128. + 119.*p3/192. + 
                         5.*p4/24. + 0.05*p5 + p6/120. + p7/1260.;
                        s2 = exp( -2.*p );
                        sasb = (zeta/p)*(1. - s1*s2);
                        s0 = 1. + 2.*p + 2.*p2 + 4.*p3/3. + 2.*p4/3.;
                        s1 = (s0 + 191.*p5/720. + 31.*p6/360. + 19.*p7/840. + 
                         17.*p8/3780. + p9/1890.)*s2;
                        sadsg = (zeta/p3)*(1. - s1);
                        dsgsb = sadsg;
                        s1 = (s0 + 511.*p5/1920. + 253.*p6/2880. + 2237.*p7/90720. + 
                         71.*p8/11340. + p9/630. + 13.*p10/34020. + p11/17010.)*s2;
                        dsgsg = (6.*zeta/p5)*(1. - s1);
                        s1 = (s0 + 1027.*p5/3840. + 521.*p6/5760. + 67.*p7/2520. + 
                         67.*p8/10080. + 19.*p9/15120. + p10/7560.)*s2;
                        dpipi = (4.*zeta/p5)*(1. - s1);
                        s1 = (s0 + 253.*p5/960. + 119.*p6/1440. + 11.*p7/560. + 
                         p8/315. + p9/3780.)*s2;
                        ddldl = (zeta/p5)*(1. - s1);
L_90:
                        s0 = fabs( spi );
                        spi *= b1*b2 + c1*c2;
                        ssg *= a1*a2;
                        pimatx.hc[i][k] = spi + 0.773*ssg;
                        pimatx.f[i][k] = s0;
                        yz2 = 6.75*ddldl;
                        xy2 = 6.75*dpipi;
                        xz2 = xy2;
                        xx = sasb + 3.*sadsg + 3.*dsgsb + 9.*dsgsg;
                        xy = sasb - 1.5*sadsg + 3.*dsgsb - 4.5*dsgsg;
                        yx = xy + 4.5*sadsg - 4.5*dsgsb;
                        yy = sasb - 1.5*sadsg - 1.5*dsgsb + 2.25*dsgsg + yz2;
                        yz = yy - 2.*yz2;
                        s1 = a1*a1*a2*a2*xx + (b2*b2 + c2*c2)*a1*a1*xy + (b1*b1 + 
                         c1*c1)*a2*a2*yx + (b1*b1*b2*b2 + c1*c1*c2*c2)*yy + (c1*
                         c1*b2*b2 + b1*b1*c2*c2)*yz;
                        pimatx.g[i][k] = 27.204*(s1 + 4.*((b1*b2 + c1*c2)*a1*
                         a2*xy2 + b1*c1*b2*c2*yz2));
                        pimatx.g[k][i] = pimatx.g[i][k];
                        ys = sasb - 1.5*dsgsb;
                        sy = sasb - 1.5*sadsg;
                        xs = sasb + 3.0*dsgsb;
                        sx = sasb + 3.0*sadsg;
                        if( lt == 2 )
                        {
                                p = .5*rik*(pimatx.zt[i] + pimatx.zt[k]);
                                t = p;
                        }
                        else
                        {
                                p = rik*pimatx.zt[i];
                                t = rik*pimatx.zt[k];
                        }
                        s1 = 1./rik*(1. - (1. + 1.5*p + p*p + (p*p*p)/3.)*exp( -2.* p ));
                        s2 = 1./(rik*p*p)*(1. - (1. + 2.*p + 2.*p*p + 4.*(p*p*p)/
                         3. + 2.*(p*p*p*p)/3. + 2.*(p*p*p*p*p)/9.)*exp( -2.*p ));
                        s1a = pimatx.qp[k]*(s1 - 1.5*s2);
                        s1b = pimatx.qp[k]*(s1 + 3.*s2);
                        t1 = 1./rik*(1. - (1. + 1.5*t + t*t + (t*t*t)/3.)*exp( -2.*
                         t ));
                        t2 = 1./(rik*t*t)*(1. - (1. + 2.*t + 2.*t*t + 4.*(t*t*t)/
                         3. + 2.*(t*t*t*t)/3. + 2.*(t*t*t*t*t)/9.)*exp( -2.*t ));
                        s2a = pimatx.qp[i]*(t1 - 1.5*t2);
                        s2b = pimatx.qp[i]*(t1 + 3.0*t2);
                        s1 = 1;
                        s2 = 1;
                        if( pimatx.qp[k] >= 5 )
                                s1 = 2;
                        if( pimatx.qp[k] > 5 )
                                s2 = 2;
                        t1 = 1;
                        t2 = 1;
                        if( pimatx.qp[i] >= 5 )
                                t1 = 2;
                        if( pimatx.qp[i] > 5 )
                                t2 = 2;
                        onepi = s1a - (s1*ys + s2*yz + yy + yx);
                        onesg = s1b - (s1*xs + s2*xy + xy + xx);
                        twopi = s2a - (t1*sy + t2*yz + yy + xy);
                        twosg = s2b - (t1*sx + t2*yx + yx + xx);
                        penint.pen[i][k] = 27.204*(onesg*a1*a1 + onepi*(1. - a1*a1));
                        penint.pen[k][i] = 27.204*(twosg*a2*a2 + twopi*(1. - a2*a2));
                }
        }
        for( i = 0; i < pistuf.norb; i++ )
        {
                pimatx.g[i][i] = 27.204*pimatx.zt[i]*(93./256. + 81./2880.);
        }
        return;
}
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
void diag(int norb,float  **aa)
{
        int i, j, k, l, m;
        float b, sigm;
        double **a, **s;
        double  alm, com1, com2, com3, com4, cst, omga, snt, um, v, vsubf, vsubi;


        a = dmatrix(0, norb, 0, norb);
        s = dmatrix(0, norb, 0, norb);
        
        for( i = 0; i < norb; i++ )
        {
           for( j = 0; j < norb; j++ )
           {
               a[i][j] = aa[i][j];
           }
        }        
        vsubi = 0.;
        b = 0;
        for( i = 0; i < norb; i++ )
        {
            for( k = 0; k < norb; k++ )
            {
                 if( i == k )
                 {
                        s[i][k] = 1.0;
                 }
                 else
                 {
                        s[i][k] = 0.0;
                        vsubi += a[i][k]*a[i][k];
                 }
            }
        }
        vsubi = sqrt( vsubi );
        vsubf = vsubi/1.0e8;
        sigm = pistuf.norb;
        v = vsubi;
        while( True )
        {
           v /= sigm;
           while( True )
           {
              m = 1;
              while( True )
              {
                 l = 0;
                 while( True )
                 {
                    if( fabs( a[l][m] ) >= v )
                    {
                        b = 1.0;
                        alm = -a[l][m];
                        um = .5*(a[l][l] - a[m][m]);
                        omga = alm/(sqrt( alm*alm + um*um ));
                        if( um < 0.0 )
                            omga = -omga;
                        snt = omga/(sqrt( 2.*(1. + (sqrt( 1. - omga*omga ))) ));
                        cst = sqrt( 1. - snt*snt );
                        i = 0;
                        while( True )
                        {
                            com1 = a[i][l];
                            com2 = a[i][m];
                            a[i][l] = com1*cst - com2*snt;
                            a[i][m] = com1*snt + com2*cst;
                            com1 = s[i][l];
                            com2 = s[i][m];
                            s[i][l] = com1*cst - com2*snt;
                            s[i][m] = com1*snt + com2*cst;
                            if( i == norb-1 )
                                 break;
                            i += 1;
                         }
                         com1 = a[l][l];
                         com2 = a[m][m];
                         com3 = a[l][m];
                         com4 = a[m][l];
                         a[l][l] = com1*cst - com4*snt;
                         a[m][m] = com2*cst + com3*snt;
                         a[l][m] = com3*cst - com2*snt;
                         a[m][l] = a[l][m];
                         i = 0;
                         while( True )
                         {
                             a[l][i] = a[i][l];
                             a[m][i] = a[i][m];
                             if( i == norb-1 )
                                break;
                             i += 1;
                          }
                    }
                    if( !(l - m + 1) )
                         break;
                    l += 1;
                 }
                 if( m == norb-1 )
                     break;
                 m += 1;
              }
              if( b != 1.0 )
                  break;
              b = 0;
           }
           if( v <= vsubf )
              break;
        }
        for( i = 0; i < norb; i++ )
        {
            for( j = 0; j < norb; j++ )
            {
                aa[i][j] = a[i][j];
                pimatx.hc[i][j] = s[i][j];
            }
        }

        free_dmatrix(a, 0, norb, 0, norb);
        free_dmatrix(s, 0, norb, 0, norb);        
        return;
}
// =============================================

void border(int icall)
{
   int i,k,ii, nread,nn;
   float  rik, a1, b1, c1, a2, b2, c2;


        /*  CALCULATE P & P*B, STORE IN G-MATRIX (ZERO OUT DIAGONAL). 
         *  NOTE:   THE FACTOR -2.6915 CONVERTS ETHYLENE BETA TO UNITY. */
         
        nn = pistuf.norb - 1;
        for( i = 0; i < pistuf.norb; i++ )
        {
                for( k = i; k < pistuf.norb; k++ )
                {
                        if( i == k )
                        {
                                pimatx.g[i][i] = densit.den[i][i];
                        }
                        else
                        {
                                pimatx.g[i][k] = densit.den[i][k];
                                pimatx.g[k][i] = densit.den[i][k]*pimatx.hc[i][k]/(-2.6915);
                        }
                }
        }
        /*  DENSITY MATRIX IN upper HALF MATRIX FORM IS TRANSTERED TO R(I,K) TO PRINT 
         *  WHOLE MATRIX  */

        for( i = 0; i < pistuf.norb; i++ )
        {
                for( k = i; k < pistuf.norb; k++ )
                {
                        mtape.pensave[i][k] = pimatx.g[i][k];
                        mtape.pensave[k][i] = pimatx.g[i][k];
                }
        }
        if (pistuf.jprint) printmatrix(mtape.pensave,"Density Matrix");
        for( i = 0; i < pistuf.norb; i++ )
        {
                pimatx.eedd[i] = pimatx.g[i][i];
                pimatx.g[i][i] = 0.;
        }
        /* save density matrix in pensave */
        if( pistuf.nplane == 1 )
        {
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( k = 0; k < pistuf.norb; k++ )
                        {
                                mtape.pensave[i][k] = pimatx.g[i][k];
                        }
                }
        }
        /*  IF THIS IS FIRST CALL TO BORDER SAVE P & P*B AND F-MATRIX ON JSAVE.  */
        if( icall == 1 )
        {
                /*  IF NPLANE=1 SAVE DIR.COS. FOR PLANAR VESCF CALC'N. 
                 *  ALSO SAVE 'PLANAR' F-MATRIX IN GAMMA (RETRIEVE IN STMT 40)  */
                if( pistuf.nplane == 1 )
                {
                        for( i = 0; i < pistuf.norb; i++ )
                        {
                                for( k = 0; k < pistuf.norb; k++ )
                                {
                                        pimatx.gamma[i][k] = pisave.fksave[i][k];
                                }
                        }
                        nread = 0;
                        for( i = 0; i < nn; i++ )
                        {
                                ii = i + 1;
                                for( k = ii; k < pistuf.norb; k++ )
                                {
                                        rik = pibonds.riks[nread];
                                        a1 = pibonds.a1s[nread];
                                        b1 = pibonds.b1s[nread];
                                        c1 = pibonds.c1s[nread];
                                        a2 = pibonds.a2s[nread];
                                        b2 = pibonds.b2s[nread];
                                        c2 = pibonds.c2s[nread];
                                        if( rik <= 3.78 )
                                        {
                                                /*  MODIFY DIR. COS. TO MAKE BONDED PAIRS OF PI ORBITALS PARALLEL.  */
                                                a1 = 0.;
                                                b1 = 1.;
                                                c1 = 0.;
                                                a2 = 0.;
                                                b2 = 1.;
                                                c2 = 0.;
                                                pibonds.a1s[nread] = a1;
                                                pibonds.b1s[nread] = b1;
                                                pibonds.c1s[nread] = c1;
                                                pibonds.a2s[nread] = a2;
                                                pibonds.b2s[nread] = b2;
                                                pibonds.c2s[nread] = c2;
                                        }
                                        nread += 1;
                                }
                        }
                }

                /*  NOTE:  G CONTAINS P/P*B MATRIX, F CONTAINS F-MATRIX FOR NEXT VESCF SE */
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( k = 0; k < pistuf.norb; k++ )
                        {
                                pisave.fjsave[i][k] = pimatx.f[i][k];
                                pisave.gjsave[i][k] = pimatx.g[i][k];
                        }
                }
                /*  IF NPLANE=2, NOTHING IS DONE WITH KSAVE SINCE ITS P*B IS NEVER ALTER-
                 *  ED.  IF NPLANE=0, A NEW P*B IS WRITTEN ON KSAVE.  IF NPLANE=1, THE 
                 *  'PLANAR' F-MATRIX IS READ OFF KSAVE AND THE 'PLANAR' DIR. COS. ARE 
                 *  WRITTEN ON.  */
                if( pistuf.nplane >= 1 )
                {
                        if( pistuf.nplane <= 1 )
                        {
                                /*  TRANSFER 'PLANAR' F-MATRIX FROM GAMMA TO F  */
                                for( i = 0; i < pistuf.norb; i++ )
                                {
                                        for( k = 0; k < pistuf.norb; k++ )
                                        {
                                                pimatx.f[i][k] = pimatx.gamma[i][k];
                                        }
                                }
                        }
                        goto L_80;
                }
        }

        /*  FOR NPLANE=0 & 1(WITH ICALL=2)  WRITE P*B ON KSAVE.  */
        for( i = 0; i < pistuf.norb; i++ )
        {
                for( k = 0; k < pistuf.norb; k++ )
                {
                        pisave.gksave[i][k] = pimatx.g[i][k];
                }
        }
        /*  FOR NPLANE=1 (REACHED ONLY IF ICALL=2) ADD THE F-MATRIX.  */
        if( pistuf.nplane == 1 )
        {
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( k = 0; k < pistuf.norb; k++ )
                        {
                                pisave.fksave[i][k] = pimatx.f[i][k];
                        }
                }
        }
L_80:
        return;
} 
/* --------------------------------------- */
void  piden()
{
     int i, il;
     
        /*  if ndc=4 adds pi charge density to each pi atom i */
        for( i = 0; i < pistuf.norb; i++ )
        {
                il = pistuf.listpi[i];
                atom[il].charge = atom[il].sigma_charge + pimatx.q[i] - pimatx.eedd[i];
        }
        return;
}

void piseq(int nt,int iset)
{
   int i, k;


        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         *
         *  THIS SUBROUTINE CALLS THE VESCF SUBROUTINES AND FROM THE RESULTING
         *  BOND-ORDER (P) AND BOND-ORDER*BETA (P*B) MATRIX, CALCULATES TORSIONAL
         *  CONSTANTS AND STRETCHING CONSTANTS IN PICON.  THE READING,
         *  WRITING AND PUNCHING OF MATRICES IS ALSO DONE THRU THIS SUBROUTINE.
         *  IF NT=-1 IT READS MATRICES FROM CARDS, WRITES THEM ON SCRATCH-TAPES
         *  AND PRINTS THEM.  (IF ANY REQUIRED MATRICES ARE MISSING IT GENERATES
         *  A HUCKEL MATRIX IN ITS PLACE).  IF NT=-2 IT ONLY WRITES MATRICES
         *  ON PAPER AND PUNCHES THEM.
         *
         * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         * */

        /*  AN ENTRY POINT FOR READING MATRICES FROM CARDS AND WRITING THEM
         *  ON JSAVE & KSAVE. */
        if( nt == -1 )
        {
                piread( iset );

        }
        else if( !nt )
        {
                for( i = 0; i < pistuf.norb; i++ )
                {
                        for( k = 0; k < pistuf.norb; k++ )
                        {
                                pimatx.f[i][k] = pisave.fjsave[i][k];
                        }
                }
                /*       CALCULATE DIRECTION COSINES & WRITE ON JSAVE. */
                dircos();
                /*       CALCULATE PI-ELECTRON DISTRIBUTION. */
                pimatx.npbmx = pistuf.nplane;
                ed123.nscft++;
                vescf();
                /*       CALCULATE BOND ORDERS AND WRITE THEM ON JSAVE (ALONG WITH THE LAST
                 *       F-MATRIX GENERATED BY VESCF).  IF NPLANE=0 WRITE MOST RECENT P*B
                 *       MATRIX ON KSAVE. */
                border(1);
                /*       SELECT METHOD OF CALCULATING TORSIONAL CONSTANTS.  IF NPLANE=0,
                 *       NOTHING IS DONE SINCE P*B IS ALREADY ON KSAVE.  IF NPLANE=1 A 'PLANAR
                 *       VESCF MUST BE DONE TO GENERATE P*B MATRIX. */
                if( pistuf.nplane == 1 )
                {
                        /*         CALCULATE 'PLANAR' PI-ELECTRON DISTRIBUTION.  NOTE:  THE 'PLANAR'
                         *         F-MATRIX WAS ALREADY READ OFF KSAVE IN BORDER AND 'PLANAR' DIRECTION
                         *         COSINES WRITTEN ON. */
                        pimatx.npbmx = 0;
                        vescf();
                        /*         CALCULATE 'PLANAR' P*B MATRIX, WRITE IT ON KSAVE ALONG WITH THE
                         *         'PLANAR' F-MATRIX JUST GENERATED. */
                        border(2);
                }
                /*       MODIFY TORSIONAL AND STRETCHING CONSTANTS. */
                picon();

        }
        else if( nt == -2 )
        {
                /*       FOR WRITING MATRICES ON PAPER AND (IF NT=-2) PUNCHING THEM ON CARDS. */
                pirite();
        }
        else
        {
                goto L_1;
        }
        return;
L_1:
        ;
}

