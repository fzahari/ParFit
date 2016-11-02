#include "pcwin.h"
#include "pcmod.h"
#include "bonds_ff.h"
#include "angles.h"
#include "field.h"
#include "attached.h"
#include "vdwk.h"
#include "solv.h"
#include "skip.h"

#define STILL 1
#define HCT   2

// choices are gbsa models: Analytical Still and HCT
// allocate memory in get_mem

void ksolv()
{
    int i,j,k, jji,it,kt;
    int ia,ib,ic,id;
    double ri,ri2,rk,cc;
    double r,r2,r4,rab,rbc;
    double cosine,factor;
    double h,ratio,term;
        char hatext[60];

    solvent.doffset = -0.09;
    solvent.p1 = 0.073;
    solvent.p2 = 0.921;
    solvent.p3 = 6.211;
    solvent.p4 = 15.236;
    solvent.p5 = 1.254;

    if (solvent.type == HCT)
    {
        for (i=1; i <= natom; i++)
        {
            solvent.asolv[i] = 0.0054;
            solvent.shct[i] = 0.0;
            if (atom[i].atomnum == 1) solvent.shct[i] = 0.85;
            if (atom[i].mmx_type == 20) solvent.shct[i] = 0.85;
            if (atom[i].atomnum == 6) solvent.shct[i] = 0.72;
            if (atom[i].atomnum == 7) solvent.shct[i] = 0.79;
            if (atom[i].atomnum == 8) solvent.shct[i] = 0.85;
            if (atom[i].atomnum == 9) solvent.shct[i] = 0.88;
            if (atom[i].atomnum == 15) solvent.shct[i] = 0.86;
            if (atom[i].atomnum == 16) solvent.shct[i] = 0.96;
            if (atom[i].atomnum == 26) solvent.shct[i] = 0.88;
        }
    } else if (solvent.type == STILL)
    {
        for (i=1; i <= natom; i++)
        {
            solvent.asolv[i] = 0.0049;
        }
    }
//  assign standard radii
    for (i=1; i <= natom; i++)
    {
        solvent.rsolv[i] = 0.0;
        if (atom[i].atomnum == 1)
        {
            solvent.rsolv[i] = 1.25;
            if (atom[atom[i].iat[0]].atomnum == 7) solvent.rsolv[i] = 1.15;
            if (atom[atom[i].iat[0]].atomnum == 8) solvent.rsolv[i] = 1.05;
        } else if (atom[i].mmx_type == 20)
        {
            solvent.rsolv[i] = 0.8;
        } else if (atom[i].atomnum == 3)
        {
            solvent.rsolv[i] = 1.432;
        } else if (atom[i].atomnum == 6)
        {
            solvent.rsolv[i] = 1.90;
            jji = 0;
            for (j=0; j < 4; j++)
            {
                if (atom[i].iat[j] != 0)
                   jji++;
            }
            if (jji == 3) solvent.rsolv[i] = 1.875;
            if (jji == 2) solvent.rsolv[i] = 1.825;
        } else if (atom[i].atomnum == 7)
        {
            solvent.rsolv[i] = 1.7063;
            jji = 0;
            for (j=0; j < 4; j++)
            {
                if (atom[i].iat[j] != 0)
                   jji++;
            }
            if (jji == 4) solvent.rsolv[i] = 1.625;
            if (jji == 1) solvent.rsolv[i] = 1.60;
        } else if (atom[i].atomnum == 8)
        {
            solvent.rsolv[i] = 1.535;
            jji = 0;
            for (j=0; j < 4; j++)
            {
                if (atom[i].iat[j] != 0)
                   jji++;
            }
            if (jji == 1) solvent.rsolv[i] = 1.48;
        } else if (atom[i].atomnum == 9)
        {
            solvent.rsolv[i] = 1.47;
        } else if (atom[i].atomnum == 10)
        {
            solvent.rsolv[i] = 1.39;
        } else if (atom[i].atomnum == 11)
        {
            solvent.rsolv[i] = 1.992;
        } else if (atom[i].atomnum == 12)
        {
            solvent.rsolv[i] = 1.70;
        } else if (atom[i].atomnum == 14)
        {
            solvent.rsolv[i] = 1.80;
        } else if (atom[i].atomnum == 15)
        {
            solvent.rsolv[i] = 1.87;
        } else if (atom[i].atomnum == 16)
        {
            solvent.rsolv[i] = 1.775;
        } else if (atom[i].atomnum == 17)
        {
            solvent.rsolv[i] = 1.735;
        } else if (atom[i].atomnum == 18)
        {
            solvent.rsolv[i] = 1.70;
        } else if (atom[i].atomnum == 19)
        {
            solvent.rsolv[i] = 2.123;
        } else if (atom[i].atomnum == 20)
        {
            solvent.rsolv[i] = 1.817;
        } else if (atom[i].atomnum == 35)
        {
            solvent.rsolv[i] = 1.90;
        } else if (atom[i].atomnum == 36)
        {
            solvent.rsolv[i] = 1.812;
        } else if (atom[i].atomnum == 37)
        {
            solvent.rsolv[i] = 2.26;
        } else if (atom[i].atomnum == 53)
        {
            solvent.rsolv[i] = 2.10;
        } else if (atom[i].atomnum == 54)
        {
            solvent.rsolv[i] = 1.967;
        } else if (atom[i].atomnum == 55)
        {
            solvent.rsolv[i] = 2.507;
        } else if (atom[i].atomnum == 56)
        {
            solvent.rsolv[i] = 2.188;
        }else
        {
                    solvent.rsolv[i] = vdw1.rad[atom[i].type]*1.00; 
            sprintf(hatext,"Atom %d  radius = %f",i, solvent.rsolv[i]);
            message_alert(hatext,"Error");
            break;
        }
    }
//  atomic volumes for Still
    if (solvent.type == STILL)
    {
        for (i=1; i <= natom; i++)
        {
            if (atom[i].mmx_type != 20)
            {
               solvent.vsolv[i] = (4.0*PI/3.0)*solvent.rsolv[i]*solvent.rsolv[i]*solvent.rsolv[i];
               ri = solvent.rsolv[i];
               ri2 = ri*ri;
               for (j=0; j < MAXIAT; j++)
               {
                  if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].mmx_type != 20)
                  {
                    k = find_bond(i,atom[i].iat[j]);
                    rk = solvent.rsolv[atom[i].iat[j]];
                    r = 1.01 * bonds_ff.bl[k];
                    ratio = (rk*rk - ri2 - r*r)/(2.0*ri*r);
                    h = ri*(1.0+ratio);
                    term = (PI/3.0)*h*h*(3.0*ri-h);
                    solvent.vsolv[i] -= term;
                  }
               }
            }
        }

// 1,2 and 1,3 polarization
        cc = 4.80298*4.80298*14.39418;

        for (i=1; i <= natom; i++)
        {
            if (atom[i].mmx_type != 20)
                 solvent.gpol[i] = -0.5*cc/(solvent.rsolv[i]+solvent.doffset+solvent.p1);
        }
        
       for (i=0; i < bonds_ff.nbnd; i++)
       {
          it = bonds_ff.i12[i][0];
          kt = bonds_ff.i12[i][1];
          if (atom[it].mmx_type != 20 && atom[kt].mmx_type != 20)
          {
             r = bonds_ff.bl[i];
             r4 = r*r*r*r;
             solvent.gpol[it] += solvent.p2*solvent.vsolv[kt]/r4;
             solvent.gpol[kt] += solvent.p2*solvent.vsolv[it]/r4;
          }
       }

       for (i=0; i < angles.nang; i++)
       {
           ia = angles.i13[i][0];
           ib = angles.i13[i][1];
           ic = angles.i13[i][2];
           if (atom[ia].mmx_type != 20 && atom[ib].mmx_type != 20 && atom[ic].mmx_type != 20)
           {
             factor = 1.0;
             for (j=0; j < MAXIAT; j++)
             {
               if (atom[ia].iat[j] != 0)
               {
                   id = atom[ia].iat[j];
                   if (id == ic)
                      factor = 0.0;
                   else if (id != ib)
                   {
                       for (k=0; k < MAXIAT; k++)
                       {
                           if (atom[ic].iat[k] != 0)
                           {
                               if (atom[ic].iat[k] == id)
                                   factor = 0.5;
                           }
                       }
                   }
               }
             }
     
             id = find_bond(ia,ib);
             rab = bonds_ff.bl[id];
             id = find_bond(ib,ic);
             rbc = bonds_ff.bl[id];
             cosine = cos(angles.anat[i]/radian);
             r2 = rab*rab + rbc*rbc - 2.0*rab*rbc*cosine;
             r4 = r2*r2;
             solvent.gpol[ia] += factor*solvent.p3*solvent.vsolv[ic]/r4;
             solvent.gpol[ic] += factor*solvent.p3*solvent.vsolv[ia]/r4;
           }
       }
    }
}
/* =================================================== */
void born()
{
    int i,k;
    double cc;
    double ratio;
    double xi,yi,zi,ri;
    double rk,sk,sk2,sum;
    double lik,lik2,uik,uik2;
    double xr,yr,zr,rvdw;
    double r,r2,r4;
    double gpi,pip5,p5inv;
    double theta,term,ccf;

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
               gpi = solvent.gpol[i];
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
                  if (skip[k] != i && atom[k].mmx_type != 20)
                  {
                    xr = atom[k].x - xi;
                    yr = atom[k].y - yi;
                    zr = atom[k].z - zi;
                    r2 = xr*xr + yr*yr + zr*zr;
                    r4 = r2*r2;
                    rvdw = solvent.rsolv[i] + solvent.rsolv[k];
                    ratio = r2/(rvdw*rvdw);
                    if (ratio > p5inv)
                    {
                        ccf = 1.0;
                    } else
                    {
                        theta = ratio*pip5;
                        term = 0.5*(1.0-cos(theta));
                        ccf = term*term;
                    }
                    gpi += solvent.p4*ccf*solvent.vsolv[k]/r4;
                  }
               }
               solvent.rborn[i] = -0.5*cc/gpi;
            }
        }
    }else if (solvent.type == HCT)
    {
        for (i=1; i <= natom; i++)
        {
            if (atom[i].mmx_type != 20)
            {
               xi = atom[i].x;
               yi = atom[i].y;
               zi = atom[i].z;
               ri = solvent.rsolv[i] + solvent.doffset;
               sum = 1.0/ri;
               for (k=1; k <= natom; k++)
               {
                   if (i != k && atom[k].mmx_type != 20)
                   {
                      xr = atom[k].x - xi;
                      yr = atom[k].y - yi;
                      zr = atom[k].z - zi;
                      r2 = xr*xr + yr*yr + zr*zr;
                      r = sqrt(r2);
                      rk = solvent.rsolv[k]+solvent.doffset;
                      sk = rk * solvent.shct[k];
                      sk2 = sk*sk;
                      if (ri < (r+sk))
                      {
                        lik = 1.0/ MaxFun(ri,r-sk);
                        uik = 1.0/(r+sk);
                        lik2 = lik*lik;
                        uik2 = uik*uik;
                        term = lik - uik + 0.25*r*(uik2-lik2) + (0.5/r)*log(uik/lik)
                              + (0.25*sk2/r)*(lik2-uik2);
                        sum -= 0.5*term;
                      }
                   }
               }
               theta = 1.0/sum;
               solvent.rborn[i] = MaxFun(ri,theta);
            }
        }
    }
}
                   
            

              
      
