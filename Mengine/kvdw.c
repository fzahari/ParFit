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
#include "nonbond.h"
#include "pot.h"
#include "field.h"
#include "atom_k.h"
#include "vdwk.h"
#include "units.h"
#include "mm3hbond.h"
#include "ebond.h"
               
void kvdw()
{
   int i, j, k;
   int it, ita, itb, found, ianum,ibnum;
   int dapair;
   char pa[4],pb[4],pt[7];
   long int hbmask;
   float rad_ia, rad_ib;
   float eps1, eps2, seps1, seps2;
   float rii,rjj,gammaij, rij, epsij;
   
// mm3 setup
   if (field.type == MM3)
   {
        mm3hbond.npair = 0;
        for (i=0; i < vdwpr_k.nvdwpr; i++)
        {
            if (atom_k.number[vdwpr_k.ia1[i]] == 1)  // looking for hydrogens
            {
                mm3hbond.htype[mm3hbond.npair] = vdwpr_k.ia1[i];
                mm3hbond.otype[mm3hbond.npair] = vdwpr_k.ia2[i];
                mm3hbond.npair++;
            } else if (atom_k.number[vdwpr_k.ia2[i]] == 1) 
            {
                if (!( (vdwpr_k.ia2[i] == 5 && vdwpr_k.ia1[i] == 1) || (vdwpr_k.ia2[i] == 36 && vdwpr_k.ia1[i] == 1) ))
                {
                        mm3hbond.htype[mm3hbond.npair] = vdwpr_k.ia2[i];
                        mm3hbond.otype[mm3hbond.npair] = vdwpr_k.ia1[i];
                        mm3hbond.npair++;
                }
            }
        }
   }
//                
   nonbond.npair = 0;
   hbmask = 1 << HBOND_MASK;
   for (i=0; i <= nonbond.maxnbtype; i++)
   {
       for (j=0; j <= nonbond.maxnbtype; j++)
           nonbond.iNBtype[j][i] = 0;
   }
   
   for (i=1; i < natom; i++)
   {
       ita = atom[i].type;
       ianum = atom[i].atomnum;
       for (j= i+1; j <= natom; j++)
       {
           itb = atom[j].type;
           ibnum = atom[j].atomnum;
           if (ita < itb)
              it = 1000*ita + itb;
           else
              it = 1000*itb + ita;
              
           found = FALSE;
           if (nonbond.iNBtype[ita][itb] != 0)
              found = TRUE;

           if (found == FALSE)
           {
               if (field.type == MMFF94)
               {
                        dapair = FALSE;
                        rii = vdw1.a[ita]*pow(vdw1.alpha[ita],0.25);
                        rjj = vdw1.a[itb]*pow(vdw1.alpha[itb],0.25);
                        if ((atom[i].atomnum == 1 && ita != 5) || (atom[j].atomnum ==1  && itb != 5) )  // polar hydrogens
                        {
                                rij = 0.5*(rii+rjj);
                                if ( (strcmp(vdw1.da[ita],"A") == 0) && (strcmp(vdw1.da[itb],"D") == 0) )
                                {
                                        dapair = TRUE;
                                }
                                if ( (strcmp(vdw1.da[itb],"A") == 0) && (strcmp(vdw1.da[ita],"D") == 0) )
                                {
                                        dapair = TRUE;
                                }
                        } else
                        {
                                gammaij = (rii-rjj)/(rii+rjj);
                                rij = 0.5*(rii+rjj)*(1.00+0.2*(1.00- exp(-12.0*gammaij*gammaij)));
                        }
                        epsij = ((181.16*vdw1.g[ita]*vdw1.g[itb]*vdw1.alpha[ita]*vdw1.alpha[itb])/
                           ( sqrt(vdw1.alpha[ita]/vdw1.n[ita]) + sqrt(vdw1.alpha[itb]/vdw1.n[itb])) ) / pow(rij,6.0);

                        if (dapair == TRUE)
                        {
                                rij *= 0.8;
                                epsij *= 0.5;
                        }

                   nonbond.vrad[ita][itb] = rij;
                   nonbond.veps[ita][itb] = epsij;
                   nonbond.ipif[ita][itb] = 1.0;
                   nonbond.vrad[itb][ita] = rij;
                   nonbond.veps[itb][ita] = epsij;
                   nonbond.ipif[itb][ita] = 1.0;

                   nonbond.vrad14[ita][itb] = nonbond.vrad[ita][itb];
                   nonbond.vrad14[itb][ita] = nonbond.vrad[itb][ita];
                   nonbond.veps14[ita][itb] = nonbond.veps[ita][itb];
                   nonbond.veps14[itb][ita] = nonbond.veps[itb][ita];

               } else
               {
                  rad_ia = vdw1.rad[ita];
                  rad_ib = vdw1.rad[itb];
                  if (strcmp(field.radiustype,"SIGMA") == 0)
                  {
                      rad_ia *= 1.122462048;
                      rad_ib *= 1.122462048;
                  }  
                  if (strcmp(field.radiussize,"DIAMETER") == 0)
                  {
                      rad_ia /= 2.0;
                      rad_ib /= 2.0;
                  }
                  eps1 = vdw1.eps[ita];
                  eps2 = vdw1.eps[itb];
       
                  if (strcmp(field.radiusrule,"ARITHMETIC") == 0)
                  {
                     nonbond.vrad[ita][itb] = rad_ia + rad_ib;
                     nonbond.vrad[itb][ita] = rad_ia + rad_ib;
                  } else if (strcmp(field.radiusrule,"GEOMETRIC") == 0)
                  {
                     nonbond.vrad[ita][itb] = 2.0*( sqrt(rad_ia)*sqrt(rad_ib));
                     nonbond.vrad[itb][ita] = 2.0*( sqrt(rad_ia)*sqrt(rad_ib));
                  } else if (strcmp(field.radiusrule,"CUBIC-MEAN") == 0)
                  {
                     if ((atom[i].atomnum == 1 && ita != 5) || (atom[j].atomnum ==1  && itb != 5) )  // polar hydrogens
                     {
                        nonbond.vrad[ita][itb] = rad_ia + rad_ib;
                        nonbond.vrad[itb][ita] = rad_ia + rad_ib;
                     }else
                     {
                        nonbond.vrad[ita][itb] = 2.0*( ((rad_ia*rad_ia*rad_ia) + (rad_ib*rad_ib*rad_ib))/(rad_ia*rad_ia + rad_ib*rad_ib) );
                        nonbond.vrad[itb][ita] = 2.0*( ((rad_ia*rad_ia*rad_ia) + (rad_ib*rad_ib*rad_ib))/(rad_ia*rad_ia + rad_ib*rad_ib) );
                     }
                  } else
                  {
                     nonbond.vrad[ita][itb] = rad_ia + rad_ib;
                     nonbond.vrad[itb][ita] = rad_ia + rad_ib;
                  }
      
                  if (strcmp(field.epsrule,"ARITHMETIC") == 0)
                  {
                      nonbond.veps[ita][itb] = 0.5*(eps1 + eps2);
                      nonbond.veps[itb][ita] = 0.5*(eps1 + eps2);
                  } else if (strcmp(field.epsrule,"GEOMETRIC") == 0)
                  {
                      nonbond.veps[ita][itb] = sqrt( eps1*eps2 );
                      nonbond.veps[itb][ita] = sqrt( eps1*eps2 );
                  } else if (strcmp(field.epsrule,"HARMONIC") == 0)
                 {
                      nonbond.veps[ita][itb] = 2.0* (eps1*eps2)/(eps1+eps2);
                      nonbond.veps[itb][ita] = 2.0* (eps1*eps2)/(eps1+eps2);
                  } else if (strcmp(field.epsrule,"HHG") == 0)
                  {
                     if ((atom[i].atomnum == 1 && ita != 5) || (atom[j].atomnum ==1  && itb != 5) )  // polar hydrogens
                     {
                       nonbond.veps[ita][itb] = 0.5*(eps1 + eps2);
                       nonbond.veps[itb][ita] = 0.5*(eps1 + eps2);
                     } else
                     {
                        seps1 = sqrt(eps1);
                        seps2 = sqrt(eps2);
                        nonbond.veps[ita][itb] = 4.0*(eps1*eps2)/( (seps1+seps2)*(seps1+seps2));
                        nonbond.veps[itb][ita] = 4.0*(eps1*eps2)/( (seps1+seps2)*(seps1+seps2));
                     }
                  }else
                  {
                      nonbond.veps[ita][itb] = sqrt( eps1*eps2 );
                      nonbond.veps[itb][ita] = sqrt( eps1*eps2 );
                  }

                  if (field.type == MM3)
                  {
                      nonbond.ipif[ita][itb] = 1;
                      nonbond.ipif[itb][ita] = 1;
                  } else if (field.type == MMX)
                  {
                      nonbond.ipif[ita][itb] = 1;
                      nonbond.ipif[itb][ita] = 1;
                      if( (ita == 2 || ita == 4 || ita == 40 || ita == 48) && 
                          (itb == 2 || itb == 4 || itb == 40 || itb == 48) )
                          {
                               nonbond.ipif[ita][itb] = 0;
                               nonbond.ipif[itb][ita] = 0;
                          }
                  }
                   nonbond.vrad14[ita][itb] = nonbond.vrad[ita][itb];
                   nonbond.vrad14[itb][ita] = nonbond.vrad[itb][ita];
                   nonbond.veps14[ita][itb] = nonbond.veps[ita][itb];
                   nonbond.veps14[itb][ita] = nonbond.veps[itb][ita];
               } 
  // look for any special vdw pairs
               numeral(atom[i].type,pa,3);
               numeral(atom[j].type,pb,3);
      
              if (atom[i].type < atom[j].type)
              {
                strcpy(pt,pa);
                strcat(pt,pb);
              } else
              {
                strcpy(pt,pb);
                strcat(pt,pa);
              }

// look for vdw pairs
                for (k= 0; k < vdwpr_k.nvdwpr ; k++)
                {
                   if (strcmp(vdwpr_k.kv[k],pt) == 0)
                   {
                      nonbond.vrad[ita][itb] = vdwpr_k.radius[k];
                      nonbond.vrad[itb][ita] = vdwpr_k.radius[k];
                      nonbond.veps[ita][itb] = vdwpr_k.eps[k];
                      nonbond.veps[itb][ita] = vdwpr_k.eps[k];
                      if (ianum != 1 && ibnum != 1 )
                      {
                         nonbond.veps[ita][itb] /= units.dielec;
                         nonbond.veps[itb][ita] /= units.dielec;
                      }
                      if (ita == 1 && itb == 5)
                      {
                            nonbond.vrad14[ita][itb] = nonbond.vrad[ita][itb];
                            nonbond.vrad14[itb][ita] = nonbond.vrad[itb][ita];
                            nonbond.veps14[itb][ita] = nonbond.veps[itb][ita];
                            nonbond.veps14[ita][itb] = nonbond.veps[ita][itb];
                      }
                      break;
                   }
                }

             nonbond.iNBtype[ita][itb] = it;
             nonbond.iNBtype[itb][ita] = it;
             
             if (pot.use_hbond)
             {
                 if ( atom[i].flags & hbmask || atom[i].mmx_type == 20)
                 {
                     for (k=0; k < MAXIAT; k++)
                     {
                         if (atom[j].iat[k] != 0 && atom[j].bo[k] != 9)
                         {
                             if (atom[atom[j].iat[k]].flags & hbmask || atom[atom[j].iat[k]].mmx_type == 20)
                             {
                                nonbond.vrad[ita][itb] -= 0.25;
                                nonbond.vrad[itb][ita] -= 0.25;
                                break;
                             }
                         }
                     }
                 } else if (atom[j].flags & hbmask || atom[j].mmx_type == 20)
                 {
                     for (k=0; k < MAXIAT; k++)
                     {
                         if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
                         {
                             if (atom[atom[i].iat[k]].flags & hbmask || atom[atom[i].iat[k]].mmx_type == 20)
                             {
                                nonbond.vrad[ita][itb] -= 0.25;
                                nonbond.vrad[itb][ita] -= 0.25;
                                break;
                             }
                         }
                     }
                 }
             }
             nonbond.npair++;
           }
       }
   }
   //
   for (i=1; i <= natom ; i++)
   {
     ita = atom[i].type;
     if (fabs(vdw1.rad[ita]) < 0.1)
       {
	 rii = vdw1.a[ita]*pow(vdw1.alpha[ita],0.25);
	 vdw1.rad[ita] = 0.5*rii;
       }
   }
}

