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
#include "units.h"
#include "minim_values.h"
#include "mm3hbond.h"
#include "angles.h"
        
int isbond(int,int);
int isangle(int,int);
int is_hbondmm3(int,int);
float get_bond(int,int);
        
EXTERN int *skip;

void image(double *, double *, double *, int);
// =========================================
int is_hbondmm3(int ia,int ib)
{
    int i, ita,itb;
    
    if (atom[ia].atomnum == 1)
    {
       ita = atom[ia].type;
       itb = atom[ib].type;
    } else if (atom[ib].atomnum == 1)
    {
       ita = atom[ib].type;
       itb = atom[ia].type;
   } else
       return FALSE;

    for (i=0; i < mm3hbond.npair; i++)
    {
        if (ita == mm3hbond.htype[i] && itb == mm3hbond.otype[i])
           return TRUE;
    } 
    return FALSE;
}
// =====================================================================
void ebuckmm3()
{
      int i,ia, ib, ita, itb, j,k;
      int i1,i2,i3;
      float ideal;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double xab,yab,zab,rab2,rab;
      double xcb,ycb,zcb,rcb2;
      double dot,cosine,fterm;
      double e,p,p2,p6,p12,rv,eps;
      double rik2;
      double expcut,expterm,expmerge;
      double rdn;
      double cutoff, pif;

      cutoff = cutoffs.vdwcut*cutoffs.vdwcut;

      if (minim_values.iprint)
      {
          fprintf(pcmoutfile,"\nVDW Terms - MM3 Buckingham Potential\n");
          fprintf(pcmoutfile,"          At1        At2    Rik    Radius     Eps          Evdw\n");
      }

      rdn = 1.0;
      if (field.type == MMX || field.type == MM2 )
        rdn = .915;
      else if (field.type == MM3)
        rdn = .923;
        
      expcut = 2.0;
      expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));

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
         ita = atom[ia].type;

         for(j=i+1; j <= natom; j++)
         {
           if (atom[i].use || atom[j].use)
           {
              if (skip[j] != i)
              {
                 ib = j;
                 itb = atom[ib].type;
                         
                 if (atom[ia].atomnum == 1)
                 {
                     // compute reduced distance
                     xi = rdn*(atom[ia].x-atom[atom[ia].iat[0]].x) + atom[atom[ia].iat[0]].x;
                     yi = rdn*(atom[ia].y-atom[atom[ia].iat[0]].y) + atom[atom[ia].iat[0]].y;
                     zi = rdn*(atom[ia].z-atom[atom[ia].iat[0]].z) + atom[atom[ia].iat[0]].z;
                 } else
                 {
                     xi = atom[ia].x;
                     yi = atom[ia].y;
                     zi = atom[ia].z;
                 }
                 if (atom[ib].atomnum == 1)
                 {
                     // compute reduced distance
                     xk = rdn*(atom[ib].x-atom[atom[ib].iat[0]].x) + atom[atom[ib].iat[0]].x;
                     yk = rdn*(atom[ib].y-atom[atom[ib].iat[0]].y) + atom[atom[ib].iat[0]].y;
                     zk = rdn*(atom[ib].z-atom[atom[ib].iat[0]].z) + atom[atom[ib].iat[0]].z;
                 } else
                 {
                     xk = atom[ib].x;
                     yk = atom[ib].y;
                     zk = atom[ib].z;
                 }
                 xr = xi - xk;
                 yr = yi - yk;
                 zr = zi - zk;
                 rik2 = xr*xr + yr*yr + zr*zr;
                 if (rik2 < cutoff)
                 {
                    fterm = 1.0;
                    if (is_hbondmm3(ia,ib) && skip[j] != -i)
                    {
                        rv = nonbond.vrad[ita][itb];
                        eps = nonbond.veps[ita][itb]/units.dielec;
                        pif = (double)nonbond.ipif[ita][itb];
                        if (atom[ia].atomnum == 1)
                        {
                            i1 = ia;
                            i2 = atom[ia].iat[0];
                            i3 = ib;
                        } else
                        {
                            i1 = ib;
                            i2 = atom[ib].iat[0];
                            i3 = ia;
                        }
                        xab = atom[i1].x - atom[i2].x;
                        yab = atom[i1].y - atom[i2].y;
                        zab = atom[i1].z - atom[i2].z;
                        rab2 = xab*xab + yab*yab + zab*zab;
                        rab = sqrt(rab2);
                        xcb = atom[i3].x - atom[i2].x;
                        ycb = atom[i3].y - atom[i2].y;
                        zcb = atom[i3].z - atom[i2].z;
                        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
                        if (rcb2 < 0.01) rcb2 = 0.01;
                        dot = xab*xcb + yab*ycb + zab*zcb;
                        cosine = dot/sqrt(rab2*rcb2);
                        ideal = get_bond(i1,i2);
                        fterm = cosine*(rab/ideal);
                         
                    } else
                    {
                        rv = nonbond.vrad[ita][itb];
                        eps = nonbond.veps[ita][itb];
                        pif = (double)nonbond.ipif[ita][itb];
                    }
                    if (skip[j] == -i)
                    {
                        rv = nonbond.vrad14[ita][itb];
                        eps = nonbond.veps14[ita][itb];
                        eps =eps/units.v14scale;
                        pif = 1.0;
                    }
                     p2 = (rv*rv)/rik2;
                     p6 = p2*p2*p2;
//                     if (p2 < 14.0 )
					 if (p2 < 9.5 )  // changed 6/20/2011
                     {
                        p = sqrt(p2);
                        expterm = units.aterm*exp(-units.bterm/p);
                        e = eps *(expterm- pif*units.cterm*p6*fterm);
                     } else
                     {
                        p12 = p6*p6;
                        e = expmerge*eps*p12;
                     }
                     if (e > 10000.)
                         fprintf(pcmoutfile,"Error in evdw %d %d\n",ia,ib);
                         
                     if (skip[j] == -i)
                       energies.e14 += e;
                     else
                       energies.evdw += e;
                     atom[ia].energy += e;
                     atom[ib].energy += e;
                     if (minim_values.iprint)
                    {
                           if ( is_hbondmm3(ia,ib))
                             fprintf(pcmoutfile,"VDW_MM3: %2s(%-3d) - %2s(%-3d) %-8.3f %-8.3f %-8.3f = %8.4f   HB\n",
                                 atom[ia].name, ia, atom[ib].name, ib, sqrt(rik2),rv,eps,e);
                           else
                             fprintf(pcmoutfile,"VDW_MM3: %2s(%-3d) - %2s(%-3d) %-8.3f %-8.3f %-8.3f = %8.4f\n",
                                 atom[ia].name, ia, atom[ib].name, ib, sqrt(rik2),rv,eps,e);
                     }
                 }
              }
           }
         }
      }
}
// ============================================================
void ebuckmm31()
{
      int i,ia, ib, ita, itb, j,k, iv, kv;
      int i1,i2,i3, use_hb;
      float ideal;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double xab,yab,zab,rab2,rab;
      double xcb,ycb,zcb,rcb2;
      double sine,ratio,deddr,deddt;
      double xp,yp,zp,rp,dot,cosine,fterm;
      double e,p,p2,p6,p12,rv,eps;
      double rik,rik2;
      double expcut,expterm,expmerge;
      double rdn,rvterm;
      double redi, redk;
      double dedx,dedy,dedz,de;
      double cutoff,pif;
      double term,terma,termc,dedxia,dedyia,dedzia;
      double dedxib,dedyib,dedzib, dedxic,dedyic,dedzic;

// Hay - initialize variables to remove compiler warning
      i1 = i2 = i3 = 0;
      rp = zp = yp = xp = 0.0;
      deddr = deddt = 0.0;
      rcb2 = zcb = ycb = xcb = 0.0;
      rab2 = zab = yab = xab = 0.0;

      cutoff = cutoffs.vdwcut*cutoffs.vdwcut;
      rdn = 1.0;
      if (field.type == MMX || field.type == MM2 )
        rdn = .915;
      else if (field.type == MM3)
        rdn = .923;
        
      expcut = 2.0;
      expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));
// non-bond term
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
         redi = 1.0;
         ita = atom[i].type;
         if (atom[ia].atomnum == 1)
         {
            // compute reduced distance
             iv = atom[ia].iat[0];
             redi = 1.0-rdn;
             xi = rdn*(atom[ia].x-atom[iv].x) + atom[iv].x;
             yi = rdn*(atom[ia].y-atom[iv].y) + atom[iv].y;
             zi = rdn*(atom[ia].z-atom[iv].z) + atom[iv].z;
         } else
         {
             iv = ia;
             xi = atom[ia].x;
             yi = atom[ia].y;
             zi = atom[ia].z;
         }

         for(j=i+1; j <= natom; j++)
         {
           if (atom[i].use || atom[j].use)
           {
              if (skip[j] != i)
              {
                 ib = j;
                 itb = atom[j].type;
                   
                 if (atom[ib].atomnum == 1)
                 {
                     // compute reduced distance
                     kv = atom[ib].iat[0];
                     redk = 1.0-rdn;
                     xk = rdn*(atom[ib].x-atom[kv].x) + atom[kv].x;
                     yk = rdn*(atom[ib].y-atom[kv].y) + atom[kv].y;
                     zk = rdn*(atom[ib].z-atom[kv].z) + atom[kv].z;
                 } else
                 {
                     kv = ib;
                     redk = 1.0;
                     xk = atom[ib].x;
                     yk = atom[ib].y;
                     zk = atom[ib].z;
                 }
                 xr = xi - xk;
                 yr = yi - yk;
                 zr = zi - zk;
                 rik2 = xr*xr + yr*yr + zr*zr;

                 if (rik2 < cutoff)
                 {
                    fterm = 1.0;
                    use_hb = FALSE;
                    if (is_hbondmm3(ia,ib) && skip[j] != -i)
                    {
                        use_hb = TRUE;
                        rv = nonbond.vrad[ita][itb];
                        eps = nonbond.veps[ita][itb]/units.dielec;
                        pif = (double)nonbond.ipif[ita][itb];
                        if (atom[ia].atomnum == 1)
                        {
                            i1 = ia;
                            i2 = atom[ia].iat[0];
                            i3 = ib;
                        } else
                        {
                            i1 = ib;
                            i2 = atom[ib].iat[0];
                            i3 = ia;
                        }
                        xab = atom[i1].x - atom[i2].x;
                        yab = atom[i1].y - atom[i2].y;
                        zab = atom[i1].z - atom[i2].z;
                        rab2 = xab*xab + yab*yab + zab*zab;
                        rab = sqrt(rab2);
                        xcb = atom[i3].x - atom[i2].x;
                        ycb = atom[i3].y - atom[i2].y;
                        zcb = atom[i3].z - atom[i2].z;
                        rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
                        if (rcb2 < 0.01) rcb2 = 0.01;
                        xp = ycb*zab - zcb*yab;
                        yp = zcb*xab - xcb*zab;
                        zp = xcb*yab - ycb*xab;
                        rp = sqrt(xp*xp + yp*yp + zp*zp);
                        dot = xab*xcb + yab*ycb + zab*zcb;
                        cosine = dot / sqrt(rab2*rcb2);
                        sine = sqrt(fabs(1.0-cosine*cosine));
                        ideal = get_bond(i1,i2);
                        ratio = rab / ideal;
                        fterm = cosine * ratio;
                        deddt = -sine * ratio;
                        deddr = cosine / (rab*ideal);
                    } else
                    {
                        rv = nonbond.vrad[ita][itb];
                        eps = nonbond.veps[ita][itb];
                        pif = (double)nonbond.ipif[ita][itb];
                    }
                     if (skip[j] == -i)
                     {
                        rv = nonbond.vrad14[ita][itb];
                        eps = nonbond.veps14[ita][itb];
                        eps /= units.v14scale;
                        pif = 1.0;
                     }
                     p2 = (rv*rv)/rik2;
                     p6 = p2*p2*p2;
                     rik = sqrt(rik2);
					 //                     if (p2 < 14.0 )
					 if (p2 < 9.5 )  // changed 6/20/2011
                     {
                       p = sqrt(p2);
                       rvterm = -units.bterm / rv;
                       expterm = units.aterm*exp(-units.bterm/p);
                       e = eps *(expterm-pif*units.cterm*p6*fterm);
                       de = eps * (rvterm*expterm+6.0*pif*units.cterm*p6*fterm/rik);
                     } else
                     {
                        p12 = p6*p6;
                        e = expmerge*eps*p12;
                        de = -12.0* expmerge * eps * p12/rik;
                     }

                     de /= rik;
                     dedx = de * xr;
                     dedy = de * yr;
                     dedz = de * zr;

                     if (skip[j] == -i)
                     {
                       energies.e14 += e;
                       if (iv == ia)
                       {
                         deriv.de14[ia][0] += dedx;
                         deriv.de14[ia][1] += dedy;
                         deriv.de14[ia][2] += dedz;
                       } else
                       {
                         deriv.de14[ia][0] += dedx*rdn;
                         deriv.de14[ia][1] += dedy*rdn;
                         deriv.de14[ia][2] += dedz*rdn;
                         deriv.de14[iv][0] += dedx*redi;
                         deriv.de14[iv][1] += dedy*redi;
                         deriv.de14[iv][2] += dedz*redi;
                       }
                       if (kv == ib)
                       {
                         deriv.de14[ib][0] -= dedx;
                         deriv.de14[ib][1] -= dedy;
                         deriv.de14[ib][2] -= dedz;
                       } else
                       {
                         deriv.de14[ib][0] -= dedx*rdn;
                         deriv.de14[ib][1] -= dedy*rdn;
                         deriv.de14[ib][2] -= dedz*rdn;
                         deriv.de14[kv][0] -= dedx*redk;
                         deriv.de14[kv][1] -= dedy*redk;
                         deriv.de14[kv][2] -= dedz*redk;
                       }
                     } else
                     {
                       energies.evdw += e;
                       if (iv == ia)
                       {
                         deriv.devdw[ia][0] += dedx;
                         deriv.devdw[ia][1] += dedy;
                         deriv.devdw[ia][2] += dedz;
                       } else
                       {
                         deriv.devdw[ia][0] += dedx*rdn;
                         deriv.devdw[ia][1] += dedy*rdn;
                         deriv.devdw[ia][2] += dedz*rdn;
                         deriv.devdw[iv][0] += dedx*redi;
                         deriv.devdw[iv][1] += dedy*redi;
                         deriv.devdw[iv][2] += dedz*redi;
                       }
                       if (kv == ib)
                       {
                         deriv.devdw[ib][0] -= dedx;
                         deriv.devdw[ib][1] -= dedy;
                         deriv.devdw[ib][2] -= dedz;
                       } else
                       {
                         deriv.devdw[ib][0] -= dedx*rdn;
                         deriv.devdw[ib][1] -= dedy*rdn;
                         deriv.devdw[ib][2] -= dedz*rdn;
                         deriv.devdw[kv][0] -= dedx*redk;
                         deriv.devdw[kv][1] -= dedy*redk;
                         deriv.devdw[kv][2] -= dedz*redk;
                       }
                     }
                     if (use_hb)
                     {
                         term = eps*units.cterm*p6;
                         deddt *= term;
                         deddr *= term;
                         if (rp >0.000001)
                            terma = deddt/(rab2*rp);
                         else
                            terma = deddt/(0.000001*rab2);
                         if (rp > 0.000001)
                             termc = -deddt/(rcb2*rp);
                         else
                             termc = -deddt/(rcb2*0.000001);
                        dedxia = terma * (yab*zp-zab*yp) - deddr*xab;
                        dedyia = terma * (zab*xp-xab*zp) - deddr*yab;
                        dedzia = terma * (xab*yp-yab*xp) - deddr*zab;
                        dedxic = termc * (ycb*zp-zcb*yp);
                        dedyic = termc * (zcb*xp-xcb*zp);
                        dedzic = termc * (xcb*yp-ycb*xp);
                        dedxib = -dedxia - dedxic;
                        dedyib = -dedyia - dedyic;
                        dedzib = -dedzia - dedzic;
                        
                        deriv.devdw[i1][0] += dedxia;
                        deriv.devdw[i1][1] += dedyia;
                        deriv.devdw[i1][2] += dedzia;
                        deriv.devdw[i2][0] += dedxib;
                        deriv.devdw[i2][1] += dedyib;
                        deriv.devdw[i2][2] += dedzib;
                        deriv.devdw[i3][0] += dedxic;
                        deriv.devdw[i3][1] += dedyic;
                        deriv.devdw[i3][2] += dedzic;
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
// =============================================
void ebuckmm32(int iatom)
{
      int i,ia, ib, iv, kv, j, k, ita, itb;
      int i1,i2,i3;
      float ideal;
      double xi,yi,zi,xr,yr,zr;
      double xk,yk,zk;
      double xab,yab,zab,rab2,rab;
      double xcb,ycb,zcb,rcb2;
      double ratio;
      double dot,cosine,fterm;
      double rdn;
      double de,d2e,p,p2,p6,p12,eps,rv;
      double redi,rediv,redk,redkv;
      double redi2,rediv2,rediiv;
      double redik,redivk,redikv;
      double rik,rik2;
      double d2edx,d2edy,d2edz,term[3][3];
      double expcut,expterm,expmerge;
      double rvterm,rvterm2;
      double cutoff,pif;
      
      cutoff = cutoffs.vdwcut*cutoffs.vdwcut;

      rdn = 1.0;
      if (field.type == MMX || field.type == MM2 )
        rdn = .915;
      else if (field.type == MM3)
        rdn = .923;
        
      expcut = 2.0;
      expmerge = (units.aterm*exp(-units.bterm/expcut) - units.cterm*(pow(expcut,6.0)))/ (pow(expcut,12.0));

    for (j=1; j <= natom; j++)
      skip[j] = 0;

    for (i=1; i <= natom; i++)
    {
        ia = i;
        iv = i;
        if (atom[ia].atomnum == 1)
            iv = atom[ia].iat[0];
        if (ia == iatom || iv == iatom)
        {
          skip[ia] = ia;
          ita = atom[ia].type;
          for (k=0; k < MAXIAT; k++)
          {
              if (atom[ia].iat[k] != 0 && atom[ia].bo[k] != 9)
                   skip[atom[ia].iat[k]] = ia;
          }
          for(k=0; k < attached.n13[ia]; k++)  // iatom
            skip[attached.i13[k][ia]] = ia;
          for(k=0; k < attached.n14[ia]; k++)  // iatom
             skip[attached.i14[k][ia]] = -ia;
          if (high_coord.ncoord > 0)
          {
             for (k=0; k < high_coord.ncoord; k++)
             {
                 for (j=0; j < 3; j++)
                 {
                     if (high_coord.i13[k][j] == ia )
                     {
                      skip[high_coord.i13[k][0]] = ia;
                      skip[high_coord.i13[k][1]] = ia;
                      skip[high_coord.i13[k][2]] = ia;
                      break;
                     }
                 }
             }
          }
             
     
         if (atom[ia].atomnum == 1)
         {
           redi = rdn;
           // compute reduced distance
           xi = rdn*(atom[ia].x-atom[iv].x) + atom[iv].x;
           yi = rdn*(atom[ia].y-atom[iv].y) + atom[iv].y;
           zi = rdn*(atom[ia].z-atom[iv].z) + atom[iv].z;
         } else
         {
            redi = 1.0;
            xi = atom[ia].x;
            yi = atom[ia].y;
            zi = atom[ia].z;
         }
        
         rediv = 1.0-redi;
         redi2 = redi*redi;
         rediv2 = rediv*rediv;
         rediiv = redi*rediv;
     

     for (j=1; j <= natom; j++)
     {
         ib = j;               
        if (skip[j] != ia )
        {
            itb = atom[ib].type;
            
         if (atom[ib].atomnum == 1)
         {
             // compute reduced distance
             redk = rdn;
             kv = atom[ib].iat[0];
             xk = rdn*(atom[ib].x-atom[kv].x) + atom[kv].x;
             yk = rdn*(atom[ib].y-atom[kv].y) + atom[kv].y;
             zk = rdn*(atom[ib].z-atom[kv].z) + atom[kv].z;
         } else
         {
             redk = 1.0;
             kv = ib;
             xk = atom[ib].x;
             yk = atom[ib].y;
             zk = atom[ib].z;
         }
            
         xr = xi - xk;
         yr = yi - yk;
         zr = zi - zk;
         rik2 = xr*xr + yr*yr + zr*zr;
         if (rik2 < cutoff)
         {
              fterm = 1.0;
              if (is_hbondmm3(ia,ib) && skip[j] != -i)
              {
                  rv = nonbond.vrad[ita][itb];
                  eps = nonbond.veps[ita][itb]/units.dielec;
                  pif = (double)nonbond.ipif[ita][itb];
                  if (atom[ia].atomnum == 1)
                  {
                        i1 = ia;
                        i2 = atom[ia].iat[0];
                        i3 = ib;
                   } else
                   {
                        i1 = ib;
                        i2 = atom[ib].iat[0];
                        i3 = ia;
                   }
                   xab = atom[i1].x - atom[i2].x;
                   yab = atom[i1].y - atom[i2].y;
                   zab = atom[i1].z - atom[i2].z;
                   rab2 = xab*xab + yab*yab + zab*zab;
                   rab = sqrt(rab2);
                   xcb = atom[i3].x - atom[i2].x;
                   ycb = atom[i3].y - atom[i2].y;
                   zcb = atom[i3].z - atom[i2].z;
                   rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
                   if (rcb2 < 0.01) rcb2 = 0.01;
                   dot = xab*xcb + yab*ycb + zab*zcb;
                   cosine = dot / sqrt(rab2*rcb2);
                   ideal = get_bond(i1,i2);
                   ratio = rab / ideal;
                   fterm = cosine * ratio;
              } else
              {
                    rv = nonbond.vrad[ita][itb];
                    eps = nonbond.veps[ita][itb];
                    pif = (double)nonbond.ipif[ita][itb];
              }
//
             if (skip[j] == -iatom)
             {
                rv = nonbond.vrad14[ita][itb];
                eps = nonbond.veps14[ita][itb];
                eps /= units.v14scale;
                pif = 1.0;
             }
             p2 = (rv*rv)/rik2;
             p6 = p2*p2*p2;
             rik = sqrt(rik2);
			 //                     if (p2 < 14.0 )
			 if (p2 < 9.5 )  // changed 6/20/2011
             {
               p = sqrt(p2);
               rvterm = -units.bterm / rv;
               rvterm2 = rvterm*rvterm;
               expterm = units.aterm*exp(-units.bterm/p);
               de = eps * (rvterm*expterm+6.0*pif*units.cterm*p6*fterm/rik);
               d2e = eps * (rvterm2*expterm-42.0*pif*units.cterm*p6*fterm/rik2);
             } else
             {
                p12 = p6*p6;
                de = -12.0* expmerge * eps * p12/rik;
                d2e = 156.0 *expmerge *eps*p12/rik2;
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
             if (ia == iatom )
             {
                if (ia == iv && ib == kv)
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
                } else if (ib == kv)
                {
                   for (k=0; k < 3; k++)
                   {
                     hess.hessx[ia][k] += term[k][0]*redi2;
                     hess.hessy[ia][k] += term[k][1]*redi2;
                     hess.hessz[ia][k] += term[k][2]*redi2;                     
                     hess.hessx[ib][k] -= term[k][0]*redi;
                     hess.hessy[ib][k] -= term[k][1]*redi;
                     hess.hessz[ib][k] -= term[k][2]*redi;                     
                     hess.hessx[iv][k] += term[k][0]*rediiv;
                     hess.hessy[iv][k] += term[k][1]*rediiv;
                     hess.hessz[iv][k] += term[k][2]*rediiv;
                 }
              } else if (ia == iv)
              {
                 redkv = 1.0 - redk;
                 for (k=0; k < 3; k++)
                 {
                     hess.hessx[ia][k] += term[k][0];
                     hess.hessy[ia][k] += term[k][1];
                     hess.hessz[ia][k] += term[k][2];                    
                     hess.hessx[ib][k] -= term[k][0]*redk;
                     hess.hessy[ib][k] -= term[k][1]*redk;
                     hess.hessz[ib][k] -= term[k][2]*redk;
                     hess.hessx[kv][k] -= term[k][0]*redkv;
                     hess.hessy[kv][k] -= term[k][1]*redkv;
                     hess.hessz[kv][k] -= term[k][2]*redkv;
                 }
               } else
               {
                  redkv = 1.0 -redk;
                  redik = redi*redk;
                  redikv = redi*redkv;
                  for (k=0; k < 3; k++)
                  {
                     hess.hessx[ia][k] += term[k][0]*redi2;
                     hess.hessy[ia][k] += term[k][1]*redi2;
                     hess.hessz[ia][k] += term[k][2]*redi2;
                     hess.hessx[ib][k] -= term[k][0]*redik;
                     hess.hessy[ib][k] -= term[k][1]*redik;
                     hess.hessz[ib][k] -= term[k][2]*redik;                                          
                     hess.hessx[iv][k] += term[k][0]*rediiv;
                     hess.hessy[iv][k] += term[k][1]*rediiv;
                     hess.hessz[iv][k] += term[k][2]*rediiv;                     
                     hess.hessx[kv][k] -= term[k][0]*redikv;
                     hess.hessy[kv][k] -= term[k][1]*redikv;
                     hess.hessz[kv][k] -= term[k][2]*redikv;
                 }
               }
             }else if (iv == iatom)
             {
                 if (ib == kv)
                 {
                   for (k=0; k < 3; k++)
                   {
                     hess.hessx[ia][k] += term[k][0]*rediiv;
                     hess.hessy[ia][k] += term[k][1]*rediiv;
                     hess.hessz[ia][k] += term[k][2]*rediiv;                    
                     hess.hessx[ib][k] -= term[k][0]*rediv;
                     hess.hessy[ib][k] -= term[k][1]*rediv;
                     hess.hessz[ib][k] -= term[k][2]*rediv;                     
                     hess.hessx[iv][k] += term[k][0]*rediv2;
                     hess.hessy[iv][k] += term[k][1]*rediv2;
                     hess.hessz[iv][k] += term[k][2]*rediv2;
                  }
                } else 
                {
                     redivk = rediv*redk;
                  for (k=0; k < 3; k++)
                  {
                     hess.hessx[ia][k] += term[k][0]*rediiv;
                     hess.hessy[ia][k] += term[k][1]*rediiv;
                     hess.hessz[ia][k] += term[k][2]*rediiv;                   
                     hess.hessx[ib][k] -= term[k][0]*redivk;
                     hess.hessy[ib][k] -= term[k][1]*redivk;
                     hess.hessz[ib][k] -= term[k][2]*redivk;
                     hess.hessx[iv][k] += term[k][0]*rediv2;
                     hess.hessy[iv][k] += term[k][1]*rediv2;
                     hess.hessz[iv][k] += term[k][2]*rediv2;
                     hess.hessx[kv][k] -= term[k][0]*rediv2;
                     hess.hessy[kv][k] -= term[k][1]*rediv2;
                     hess.hessz[kv][k] -= term[k][2]*rediv2;
                  }
               }
             }
           }
         }
       }
      }
    }
}
             

