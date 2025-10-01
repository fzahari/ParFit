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
#include "torsions.h"
#include "rings.h"
#include "field.h"
#include "atom_k.h"
#include "fix.h"
#include "angles.h"
#include "pistuf.h"
#include "pisave.h"
#include "tork.h"
#include "tsbond.h"

        
void ktorsion()
{
    int i, j, ierr, itor, ii, ki, k;
    int ia, ib, ic, id;
    int ita, itb, itc, itd;
    int cl_a, cl_b, cl_c, cl_d;
    int use_ring4, use_ring5;
    long int mask;
    char izero[4];
    char pa[4],pb[4],pc[4],pd[4],pt[13], kv1[13];
    char k1[13], k2[13], k3[13], k4[13], k5[13], k6[13], k7[13], k8[13], k9[13], k10[13],
              k11[13], k12[13], k13[13], k14[13], k15[13];
    float v1, v2, v3, rangle;

    mask = (1L << 0);
    itor = 0;
    use_ring4 = FALSE;
    if (rings.nring4 > 0 && torkn1.ntor4 > 0)
       use_ring4 = TRUE;

    use_ring5 = FALSE;
    if (rings.nring5 > 0 && torkn1.ntor5 > 0)
       use_ring5 = TRUE;

    for (i=0; i < torsions.ntor; i++)
    {
        ia = torsions.i14[i][0];
        ib = torsions.i14[i][1];
        ic = torsions.i14[i][2];
        id = torsions.i14[i][3];

        ita = atom[ia].type;
        itb = atom[ib].type;
        itc = atom[ic].type;
        itd = atom[id].type;

        cl_a = atom[ia].tclass;
        cl_b = atom[ib].tclass;
        cl_c = atom[ic].tclass;
        cl_d = atom[id].tclass;

        if (field.type == MMX)
        {
            if (ita == 40)
               ita = 2;
            if (itb == 40)
               itb = 2;
            if (itc == 40)
               itc = 2;
            if (itd == 40)
               itd = 2;

            if (ita == 56 && itb == 56 && itc == 56 && itd == 56)
            {
                ita = 1; itb = 1; itc = 1; itd = 1;
            }
        }
        numeral(ita,pa,3);
        numeral(itb,pb,3);
        numeral(itc,pc,3);
        numeral(itd,pd,3);
        
        if (itb < itc )
           four(pa,pb, pc, pd, pt);
        else if (itb > itc)
           four(pd,pc, pb, pa, pt);
        else if (ita < itd)
           four(pa,pb, pc, pd, pt);
        else
           four(pd,pc, pb, pa, pt);
        strcpy(izero,"  0");

        if (field.type == MMFF94)
        {

           numeral(atom_k.tclass1[ita],pa,3);
           numeral(atom_k.tclass[itb],pb,3);
           numeral(atom_k.tclass[itc],pc,3);
           numeral(atom_k.tclass1[itd],pd,3);
           
        
           if (atom_k.tclass[itb] < atom_k.tclass[itc] )
           {
              four(pa, pb, pc, izero,k2);
              four(izero ,pb, pc, pd,k3);
           } else if (atom_k.tclass[itb] > atom_k.tclass[itc])
           {
              four(izero,pc, pb, pa, k2);
              four(pd,pc, pb, izero, k3);
           } else if (atom_k.tclass1[ita] < atom_k.tclass1[itd])
           {
              four(pa,pb, pc, izero, k2);
              four(izero,pb, pc, pd, k3);
           } else
           {
              four(izero,pc, pb, pa, k2);
              four(pd,pc, pb, izero, k3);
           }
        } else
        {
           four( izero, pc, pb, pa, k1 );
           four( pd, pc, pb, izero, k2 );
           four( izero, pb, pc, pd, k3 );
        }
        
        numeral(ita,pa,3);
        numeral(itb,pb,3);
        numeral(itc,pc,3);
        numeral(itd,pd,3);

        four( pa, pb, pc, izero, k4 );
        four( izero, pc, pb, izero, k5 );
        four( izero, pb, pc, izero, k6 );
        four( pd, pc, izero, izero, k7 );
        four( izero, izero, pc,pd, k8 );

        four( izero, izero, pb, pa, k9 );
        four( pa, pb, izero, izero, k10 );
        four( izero, izero, pc, izero, k11 );
        four( izero, pc, izero, izero, k12 );
        four( izero, izero, pb, izero, k13 );
        four( izero, pb, izero, izero, k14 );

        numeral(cl_a,pa,3);
        numeral(cl_b,pb,3);
        numeral(cl_c,pc,3);
        numeral(cl_d,pd,3);
        if (itb < itc )
           four(pa,pb, pc, pd, k15);
        else if (itb > itc)
           four(pd,pc, pb, pa, k15);
        else if (ita < itd)
           four(pa,pb, pc, pd, k15);
        else
           four(pd,pc, pb, pa, k15);

        ierr = FALSE;
//   check for pi bonds
        if ((atom[ib].flags & mask) && (atom[ic].flags & mask) )
        {
           if (field.type == MMX)
           {
              for (j=0; j < torknp.npitor; j++)
              {
                strcpy(kv1,torknp.kv[j]);
                if ( strcmp(kv1,pt) == 0 || strcmp(kv1,k1) == 0 || strcmp(kv1,k2) == 0|| strcmp(kv1,k3) == 0
                  || strcmp(kv1,k4) == 0 || strcmp(kv1,k5) == 0 || strcmp(kv1,k6) == 0 || strcmp(kv1,k7) == 0
                  || strcmp(kv1,k8) == 0 )
                {
                    torsions.v1[i] = torknp.tv1[j];
                    torsions.ph1[i] = torknp.ph1[j];
                    torsions.v2[i] = torknp.tv2[j];
                    torsions.ph2[i] = torknp.ph2[j];
                    torsions.v3[i] = torknp.tv3[j];
                    torsions.ph3[i] = torknp.ph3[j];
                    ierr = TRUE;
                    ii = -1;
                    ki = -1;
                    for (k=0; k < pistuf.norb; k++)
                    {
                        if (pistuf.listpi[k] == ib && ib < ic)
                          ii = k;
                        if (pistuf.listpi[k] == ib && ib > ic)
                          ki = k;
                        if (pistuf.listpi[k] == ic && ib < ic)
                          ki= k;
                        if (pistuf.listpi[k] == ic && ib > ic)
                          ii = k;
                    }     
                    if (itor > pistuf.ntpi)
                        message_alert("Error in pi torsion","Error");   
                    pitorsave.itsave[itor] = ii;
                    pitorsave.ktsave[itor] = ki;
                    pitorsave.ncsave[itor] = i;
                    pitorsave.tcsave[itor] = torsions.v2[i];
                    itor++;
                    goto L_10;
                }
              }
            } else if (field.type == MM3)
            {
                for(j=0; j < torkn1.ntor; j++)
                {
                    strcpy(kv1,torkn1.kv[j]);
                    if (strcmp(pt,kv1) == 0)
                    {
                       torsions.v1[i] = torkn1.tv1[j];
                       torsions.ph1[i] = torkn1.phase1[j];     
                       torsions.v2[i] = torkn1.tv2[j];
                       torsions.ph2[i] = torkn1.phase2[j];     
                       torsions.v3[i] = torkn1.tv3[j];
                       torsions.ph3[i] = torkn1.phase3[j];     
                       ierr = TRUE;
                       ii = -1;
                       ki = -1;
                       for (k=0; k < pistuf.norb; k++)
                       {
                         if (pistuf.listpi[k] == ib && ib < ic)
                          ii = k;
                         if (pistuf.listpi[k] == ib && ib > ic)
                          ki = k;
                         if (pistuf.listpi[k] == ic && ib < ic)
                          ki= k;
                         if (pistuf.listpi[k] == ic && ib > ic)
                          ii = k;
                       }     
                       if (itor > pistuf.ntpi)
                         message_alert("Error in pi torsion","Error");   
                       pitorsave.itsave[itor] = ii;
                       pitorsave.ktsave[itor] = ki;
                       pitorsave.ncsave[itor] = i;
                       pitorsave.tcsave[itor] = torsions.v2[i];
                       itor++;
                       goto L_10;
                    }
                 }
            }
            // get here with no pi parameters   
                    ii = -1;
                    ki = -1;
                    for (k=0; k < pistuf.norb; k++)
                    {
                        if (pistuf.listpi[k] == ib && ib < ic)
                          ii = k;
                        if (pistuf.listpi[k] == ib && ib > ic)
                          ki = k;
                        if (pistuf.listpi[k] == ic && ib < ic)
                          ki= k;
                        if (pistuf.listpi[k] == ic && ib > ic)
                          ii = k;
                    }
                    if (ii < 0 || ii > natom)
                      message_alert("Error assigning Pi Torsions","Error");
                    if (ki < 0 || ki > natom)
                      message_alert("Error assigning Pi Torsions","Error");
                         
                    pitorsave.itsave[itor] = ii;
                    pitorsave.ktsave[itor] = ki;
                    pitorsave.ncsave[itor] = i;
                    pitorsave.tcsave[itor] = 0.0;
                    itor++;            
        }
//   check for linear angles in high coordinate atoms
        if (high_coord.ncoord > 0)
        {
            for (k = 0; k < high_coord.ncoord; k++)
            {
                if (high_coord.i13[k][1] == ib) // high coordinate atom in center
                {
                    if (high_coord.i13[k][0] == ia && high_coord.i13[k][2] == ic)
                    {
                        angle(ia,ib,ic,&rangle);
                        if (rangle > 175 || rangle < -175)
                        {
                           torsions.v1[i] = 0.0;
                           torsions.v2[i] = 0.0;
                           torsions.v3[i] = 0.0;
                           ierr = TRUE;
                           goto L_10;
                        }
                    } else if (high_coord.i13[k][2] == ia && high_coord.i13[k][0] == ic)
                    {
                        angle(ia,ib,ic,&rangle);
                        if (rangle > 175 || rangle < -175)
                        {
                           torsions.v1[i] = 0.0;
                           torsions.v2[i] = 0.0;
                           torsions.v3[i] = 0.0;
                           ierr = TRUE;
                           goto L_10;
                        }
                    }
                } else if ( high_coord.i13[k][1] == ic)
                {
                    if (high_coord.i13[k][0] == ib && high_coord.i13[k][2] == id)
                    {
                        angle(ib,ic,id,&rangle);
                        if (rangle > 175 || rangle < -175)
                        {
                           torsions.v1[i] = 0.0;
                           torsions.v2[i] = 0.0;
                           torsions.v3[i] = 0.0;
                           ierr = TRUE;
                           goto L_10;
                        }
                    } else if (high_coord.i13[k][2] == ib&& high_coord.i13[k][0] == id)
                    {
                        angle(ib,ic,id,&rangle);
                        if (rangle > 175 || rangle < -175)
                        {
                           torsions.v1[i] = 0.0;
                           torsions.v2[i] = 0.0;
                           torsions.v3[i] = 0.0;
                           ierr = TRUE;
                           goto L_10;
                        }
                    }
                }
            }
        }
                    
//   check for transition state bonds
        if (field.type == MMX)
        {
            // check for TS atoms
            if( (itb >= 49 && itb <= 55) || itb == 43 || itb == 45 ||
                (itc >= 49 && itc <= 55) || itc == 43 || itc == 45 || 
                (itb == 29 && itc == 29) || (ita==45 && itb==29 && itc==29 && itd==50)||
                (itb==29 && itc == 49) || (itb==29 && itc ==50)||
                (itb==49 && itc == 29) || (itb==50 && itc == 29))
               {
                  transomeg( &v1, &v2, &v3, &ierr, &i, pt, k1,k2,k3,k4,k5,k6,k7,k8,k9,k10 );
                  if( ierr == TRUE )
                  {
                      torsions.v1[i] =  v1;
                      torsions.v2[i] =  v2;
                      torsions.v3[i] =  v3;
                      torsions.ph1[i] = 1;     
                      torsions.ph2[i] = -1;     
                      torsions.ph3[i] = 1;
                      goto L_10;
                  }
               }
        }
//  check for four membered rings
        if ( isbond(ia,id) )
        {
          for(j=0; j < torkn1.ntor4; j++)
          {
             strcpy(kv1,torkn1.kv4[j]);
             if (strcmp(pt,kv1) == 0 || strcmp(k1,kv1) == 0 || strcmp(k2,kv1) == 0 || strcmp(k3,kv1) == 0 ||      
                 strcmp(k4,kv1) == 0 || strcmp(k5,kv1) == 0 || strcmp(k6,kv1) == 0 || strcmp(k7,kv1) == 0 ||      
                 strcmp(k8,kv1) == 0 || strcmp(k9,kv1) == 0 || strcmp(k10,kv1) == 0 || strcmp(k11,kv1) == 0 ||      
                 strcmp(k12,kv1) == 0 || strcmp(k13,kv1) == 0 || strcmp(k14,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv41[j];
               torsions.ph1[i] = torkn1.phase41[j];     
               torsions.v2[i] =  torkn1.tv42[j];
               torsions.ph2[i] = torkn1.phase42[j];     
               torsions.v3[i] =  torkn1.tv43[j];
               torsions.ph3[i] = torkn1.phase43[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
        }
//   delocalized torsions 
        if ( is_delocalbond(ib,ic) )
        {
          for(j=0; j < torkn1.ntordel; j++)
          {
             strcpy(kv1,torkn1.kvdel[j]);
             if (strcmp(pt,kv1) == 0 || strcmp(k1,kv1) == 0 || strcmp(k2,kv1) == 0 || strcmp(k3,kv1) == 0 ||      
                 strcmp(k4,kv1) == 0 || strcmp(k5,kv1) == 0 || strcmp(k6,kv1) == 0 || strcmp(k7,kv1) == 0 ||      
                 strcmp(k8,kv1) == 0 || strcmp(k9,kv1) == 0 || strcmp(k10,kv1) == 0 || strcmp(k11,kv1) == 0 ||      
                 strcmp(k12,kv1) == 0 || strcmp(k13,kv1) == 0 || strcmp(k14,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tvdel1[j];
               torsions.ph1[i] = torkn1.phasedel1[j];     
               torsions.v2[i] =  torkn1.tvdel2[j];
               torsions.ph2[i] = torkn1.phasedel2[j];     
               torsions.v3[i] =  torkn1.tvdel3[j];
               torsions.ph3[i] = torkn1.phasedel3[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
        }
//   check for five membered rings
         if (is_ring54(ia,ib,ic,id) )
         {
          for(j=0; j < torkn1.ntor5; j++)
          {
             strcpy(kv1,torkn1.kv5[j]);
             if (strcmp(pt,kv1) == 0 || strcmp(k1,kv1) == 0 || strcmp(k2,kv1) == 0 || strcmp(k3,kv1) == 0 ||      
                 strcmp(k4,kv1) == 0 || strcmp(k5,kv1) == 0 || strcmp(k6,kv1) == 0 || strcmp(k7,kv1) == 0 ||      
                 strcmp(k8,kv1) == 0 || strcmp(k9,kv1) == 0 || strcmp(k10,kv1) == 0 || strcmp(k11,kv1) == 0 ||      
                 strcmp(k12,kv1) == 0 || strcmp(k13,kv1) == 0 || strcmp(k14,kv1) == 0 )
             {
               torsions.v1[i] =  torkn1.tv51[j];
               torsions.ph1[i] = torkn1.phase51[j];     
               torsions.v2[i] =  torkn1.tv52[j];
               torsions.ph2[i] = torkn1.phase52[j];     
               torsions.v3[i] =  torkn1.tv53[j];
               torsions.ph3[i] = torkn1.phase53[j];
               ierr = TRUE;
               goto L_10;
             }  
          }
         }
//   check regular parameters
         if (field.type == MMFF94)
         {
             for(j=0; j < torkn1.ntor; j++)   // specific constants
             {
                 strcpy(kv1,torkn1.kv[j]);
                 if (strcmp(pt,kv1) == 0)
                 {
                     torsions.v1[i] = torkn1.tv1[j];
                     torsions.ph1[i] = torkn1.phase1[j];     
                     torsions.v2[i] = torkn1.tv2[j];
                     torsions.ph2[i] = torkn1.phase2[j];     
                     torsions.v3[i] = torkn1.tv3[j];
                     torsions.ph3[i] = torkn1.phase3[j];     
                     torsions.v4[i] = torkn1.tv4[j];
                     torsions.ph4[i] = torkn1.phase4[j];     
                     torsions.v5[i] = torkn1.tv5[j];
                     torsions.ph5[i] = torkn1.phase5[j];     
                     torsions.v6[i] = torkn1.tv6[j];
                     torsions.ph6[i] = torkn1.phase6[j];
                     goto L_10;
                 }
             }
             for(j=0; j < torkn1.ntor; j++)  //  class 2
             {
                 strcpy(kv1,torkn1.kv[j]);
                 if (strcmp(k2,kv1) == 0 || strcmp(k3,kv1) == 0 )
                 {
                     torsions.v1[i] = torkn1.tv1[j];
                     torsions.ph1[i] = torkn1.phase1[j];     
                     torsions.v2[i] = torkn1.tv2[j];
                     torsions.ph2[i] = torkn1.phase2[j];     
                     torsions.v3[i] = torkn1.tv3[j];
                     torsions.ph3[i] = torkn1.phase3[j];     
                     torsions.v4[i] = torkn1.tv4[j];
                     torsions.ph4[i] = torkn1.phase4[j];     
                     torsions.v5[i] = torkn1.tv5[j];
                     torsions.ph5[i] = torkn1.phase5[j];     
                     torsions.v6[i] = torkn1.tv6[j];
                     torsions.ph6[i] = torkn1.phase6[j];
                     goto L_10;
                 }
             }
             for(j=0; j < torkn1.ntor; j++)  // wild cards
             {
                 strcpy(kv1,torkn1.kv[j]);
                 if (strcmp(k5,kv1) == 0 || strcmp(k6,kv1) == 0 || strcmp(k7,kv1) == 0 || strcmp(k8,kv1) == 0 )
                 {
                     torsions.v1[i] = torkn1.tv1[j];
                     torsions.ph1[i] = torkn1.phase1[j];     
                     torsions.v2[i] = torkn1.tv2[j];
                     torsions.ph2[i] = torkn1.phase2[j];     
                     torsions.v3[i] = torkn1.tv3[j];
                     torsions.ph3[i] = torkn1.phase3[j];     
                     torsions.v4[i] = torkn1.tv4[j];
                     torsions.ph4[i] = torkn1.phase4[j];     
                     torsions.v5[i] = torkn1.tv5[j];
                     torsions.ph5[i] = torkn1.phase5[j];     
                     torsions.v6[i] = torkn1.tv6[j];
                     torsions.ph6[i] = torkn1.phase6[j];
                     goto L_10;
                 }
             }
         }              
// check regular specific parameters
         for(j=0; j < torkn1.ntor; j++)
         {
             strcpy(kv1,torkn1.kv[j]);
             if (strcmp(pt,kv1) == 0)
             {
               torsions.v1[i] = torkn1.tv1[j];
               torsions.ph1[i] = torkn1.phase1[j];     
               torsions.v2[i] = torkn1.tv2[j];
               torsions.ph2[i] = torkn1.phase2[j];     
               torsions.v3[i] = torkn1.tv3[j];
               torsions.ph3[i] = torkn1.phase3[j];     
               torsions.v4[i] = torkn1.tv4[j];
               torsions.ph4[i] = torkn1.phase4[j];     
               torsions.v5[i] = torkn1.tv5[j];
               torsions.ph5[i] = torkn1.phase5[j];     
               torsions.v6[i] = torkn1.tv6[j];
               torsions.ph6[i] = torkn1.phase6[j];
// check for pi torsions = mm3 does not define pitorsion parameter so we need to pick it up here
               if (itb == 2 && itc == 2 && field.type == MMX)
               {
                   for (k=0; k < MAXIAT; k++)
                   {
                       if (atom[ib].iat[k] == ic)
                       {
                           if (atom[ib].bo[k] == 1 && is_delocalbond(ib,ic) )
                           {
                               torsions.v1[i] = 0.5;
                               torsions.v2[i] = 1.0;
                               torsions.v3[i] = 0.0;
                               ierr = TRUE;
                               goto L_10;
                               break;
                           }
                       }
                   }
               }
               goto L_10;
               break;
             }  
          }
// check regular generalized parameters
         for(j=0; j < torkn1.ntor; j++)
         {
             strcpy(kv1,torkn1.kv[j]);
             if (strcmp(pt,kv1) == 0 || strcmp(k1,kv1) == 0 || strcmp(k2,kv1) == 0 || strcmp(k3,kv1) == 0 ||      
                 strcmp(k4,kv1) == 0 || strcmp(k5,kv1) == 0 || strcmp(k6,kv1) == 0 || strcmp(k7,kv1) == 0 ||      
                 strcmp(k8,kv1) == 0 || strcmp(k9,kv1) == 0 || strcmp(k10,kv1) == 0 || strcmp(k11,kv1) == 0 ||      
                 strcmp(k12,kv1) == 0 || strcmp(k13,kv1) == 0 || strcmp(k14,kv1) == 0 || strcmp(k15,kv1) == 0)
             {
               torsions.v1[i] = torkn1.tv1[j];
               torsions.ph1[i] = torkn1.phase1[j];     
               torsions.v2[i] = torkn1.tv2[j];
               torsions.ph2[i] = torkn1.phase2[j];     
               torsions.v3[i] = torkn1.tv3[j];
               torsions.ph3[i] = torkn1.phase3[j];     
               torsions.v4[i] = torkn1.tv4[j];
               torsions.ph4[i] = torkn1.phase4[j];     
               torsions.v5[i] = torkn1.tv5[j];
               torsions.ph5[i] = torkn1.phase5[j];     
               torsions.v6[i] = torkn1.tv6[j];
               torsions.ph6[i] = torkn1.phase6[j];
// check for pi torsions = mm3 does not define pitorsion parameter so we need to pick it up here
               if (itb == 2 && itc == 2 && field.type == MMX)
               {
                   for (k=0; k < MAXIAT; k++)
                   {
                       if (atom[ib].iat[k] == ic)
                       {
                           if (atom[ib].bo[k] == 1 && is_delocalbond(ib,ic) )
                           {
                               torsions.v1[i] = 0.5;
                               torsions.v2[i] = 1.0;
                               torsions.v3[i] = 0.0;
                               ierr = TRUE;
                               goto L_10;
                               break;
                           }
                       }
                   }
               }
               goto L_10;
               break;
             }  
          }
//  missing parameters
          if (ierr == FALSE)
          {
	    //              Missing_constants = TRUE;
              fprintf(errfile,"Torsion constants missing: angle: %d %d %d %d  types: %d %d %d %d\n",ia,ib,ic,id,
                      ita, itb, itc, itd);
          }
L_10:
    continue;
    }

    if (itor < pistuf.ntpi)
    {
        message_alert("itor < pistuf.ntpi","ERROR");
    }
    if (allene.nallene > 0)
    {
        for (i = 0; i < allene.nallene; i++)
        {
            torsions.v1[allene.ntor[i]] = 0;
            torsions.v2[allene.ntor[i]] = -11.5;
            torsions.v3[allene.ntor[i]] = 0;
        }
    }
}
/*  --------------------------------------------   */
void four(char *pa, char *pb, char *pc, char *pd, char *pt)
{
        strcpy(pt,pa);
        strcat(pt,pb);
        strcat(pt,pc);
        strcat(pt,pd);
} 
void  transomeg(float *v1, float *v2, float *v3, int *ierr, int *ktcon_ent, char *kk, char *k1, char *k2, char *k3,
         char *k4, char *k5, char *k6, char *k7, char *k8, char *k9, char *k10)
{
         int iadj, iii, it, jjj,  kkk, nbond;

        /* do generalized first
         * will record if v1 v2 v3 are changed or if itrans flag set to 2
         * now do specific ones for overide of generalized */
        *v1 = 0.0F;
        *v2 = 0.0F;
        *v3 = 0.0F;
        nbond = 0L;
        iadj = 0L;
        for( iii = 0L; iii < 15; iii++ )
        {
                if (ts_bondorder.fbnd[iii] > 0.0F )
                        nbond += 1;
                if( ts_bondorder.fbnd[iii] > 0.0F )
                {
                        for( jjj = 1L; jjj <= natom; jjj++ )
                        {
                                if( atom[jjj].type == 49L )
                                {
                                        for( kkk = 0L; kkk < MAXIAT; kkk++ )
                                        {
                                                if( atom[jjj].iat[kkk] != 0L && atom[jjj].bo[kkk] != 9L )
                                                {
                                                        it = atom[atom[jjj].iat[kkk]].type;
                                                        if( it == 29L || it == 30L || it == 50L || 
                                                            it == 53L || it == 42L )
                                                                goto L_65434;
                                                }
                                        }
                                }
                        }
                goto L_65433;
L_65434:
                iadj = it;
                }
L_65433:
                ;
        }
        /* special constants for one bond additions to olefins */
        if( nbond < 2L )
        {
                if( ts_bondorder.fbnd[0L] > 0.0F )
                {
                        if( strcmp(k2," 49 49 50  0") == 0 || strcmp(k4," 49 49 50  0") == 0 )
                        {
                                if( iadj == 50L )
                                {
                                        *v2 = -10.;
                                        if( strcmp(k2," 49 49 50  0") == 0 )
                                        {
                                                //*ktcon_ent = k2;
                                                *ierr = 1;
                                                return;
                                        }
                                        if( strcmp(k4," 49 49 50  0") == 0  )
                                        {
                                                //*ktcon_ent = k4;
                                                *ierr = 1;
                                                return;
                                        }
                                }
                        }
                        else if( strcmp(k1,"  0 29 49 49") == 0 || strcmp(k3,"  0 29 49 49") == 0 ||  strcmp(k1,"  0 29 50 50") == 0 || strcmp(k3,"  0 29 50 50") == 0 )
                        {
                                if( iadj == 29L )
                                {
                                        *v2 = -10.;
                                        if( strcmp(k1,"  0 29 49 49") == 0 || strcmp(k1,"  0 29 50 50") == 0 )
                                        {
                                                //*ktcon_ent = k1;
                                                *ierr = 1;
                                                return;
                                        }
                                        if( strcmp(k3,"  0 29 49 49") == 0  || strcmp(k3,"  0 29 50 50") == 0 )
                                        {
                                                //*ktcon_ent = k3;
                                                *ierr = 1;
                                                return;
                                        }
                                }
                        }
                        else if( strcmp(k1,"  0 30 49 49") == 0 || strcmp(k3,"  0 30 49 49") == 0 || strcmp(k1,"  0 30 50 50") == 0 || strcmp(k3,"  0 30 50 50") == 0 )
                        {
                                if( iadj == 30L )
                                {
                                        *v2 = -10.;
                                        if( strcmp(k1,"  0 30 49 49") == 0 ||  strcmp(k1,"  0 30 50 50") == 0 )
                                        {
                                                //*ktcon_ent = k1;
                                                *ierr = 1;
                                                return;
                                        }
                                        if( strcmp(k3,"  0 30 49 49") == 0 || strcmp(k3,"  0 30 50 50") == 0 )
                                        {
                                                //*ktcon_ent = k3;
                                                *ierr = 1;
                                                return;
                                        }
                                }
                        }
                        else if( strcmp(k1,"  0 48 49 49") == 0 || strcmp(k3,"  0 48 49 49") == 0 || strcmp(k1,"  0 48 50 50") == 0 || strcmp(k3,"  0 48 50 50") == 0 )
                        {
                                if( iadj == 48L )
                                {
                                        *v2 = -10.;
                                        if( strcmp(k1,"  0 48 49 49") == 0 || strcmp(k1,"  0 48 50 50") == 0 )
                                        {
                                                //*ktcon_ent = k1;
                                                *ierr = 1;
                                                return;
                                        }
                                        if( strcmp(k3,"  0 48 49 49") == 0 || strcmp(k3,"  0 48 50 50") == 0 )
                                        {
                                                //*ktcon_ent = k3;
                                                *ierr = 1;
                                                return;
                                        }

                                }
                        }
                        else if( strcmp(k1,"  0 42 49 49") == 0 || strcmp(k3,"  0 42 49 49") == 0 ||  strcmp(k1,"  0 42 50 50") == 0 || strcmp(k3,"  0 42 50 50") == 0 )
                        {
                                if( iadj == 42L )
                                {
                                        *v2 = -10.;
                                        if( strcmp(k1,"  0 42 49 49") == 0 ||  strcmp(k1,"  0 42 50 50") == 0 )
                                        {
                                                //*ktcon_ent = k1;
                                                *ierr = 1;
                                                return;
                                        }
                                        if( strcmp(k3,"  0 42 49 49") == 0  || strcmp(k3,"  0 42 50 50") == 0 )
                                        {
                                                //*ktcon_ent = k3;
                                                *ierr = 1;
                                                return;
                                        }
                                }
                        }
                }
        }
        if( nbond < 2L && iadj == 29L )
        {
                /* 8/9/87 to match Houk's torsion */
                if( strcmp(kk,"  1 49 49 29") == 0 )
                {
                        //*ktcon_ent = kk;
                        *v1 = -.55;
                        *v2 = .55;
                        *v3 = .34;
                        *ierr = 1;
                        return;
                }
        }
// start of regular transition state constants
        if( strcmp(kk," 50 29 49 49") == 0 )
        {
/*                *v2 = -(1.0 - ts_bondorder.fbnd[0])*5.;
                *v3 = ts_bondorder.fbnd[0]*5.;*/
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 49 29 50 50") == 0 )
        {
/*                *v2 = -(1.0 - ts_bondorder.fbnd[1])*5.;
                *v3 = ts_bondorder.fbnd[1]*5.;*/
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 49 49 50 50") == 0 )
        {
                *v2 = 30.*(1. - ts_bondorder.fbnd[0]);
                if( ts_bondorder.fbnd[0L] > 0.0F && ts_bondorder.fbnd[1L] > 0.0F )
                {
                        if( ts_bondorder.fbnd[1] > ts_bondorder.fbnd[0] )
                                *v2 = 30.*(1. - ts_bondorder.fbnd[1]);
                }
//              //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 49 49 50 53") == 0 )
        {
                *v2 = 30.*(1. - ts_bondorder.fbnd[0]);
                if( ts_bondorder.fbnd[0L] > 0.0F && ts_bondorder.fbnd[6L] > 0.0F )
                {
                        if( ts_bondorder.fbnd[6] > ts_bondorder.fbnd[0] )
                                *v2 = 30.*(1. - ts_bondorder.fbnd[6]);
                }
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 45 49 50 53") == 0 )
        {
                if( ts_bondorder.fbnd[3L] > 0.0F && ts_bondorder.fbnd[6L] > 0.0F )
                {
                        *v2 = 30.*(1. - ts_bondorder.fbnd[3]);
                        if( ts_bondorder.fbnd[6] > ts_bondorder.fbnd[3] )
                                *v2 = 30.*(1. - ts_bondorder.fbnd[6]);
                }
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 45 49 50 43") == 0 )
        {
                if( ts_bondorder.fbnd[3L] > 0.0F && ts_bondorder.fbnd[4L] > 0.0F )
                {
                        *v2 = 30.*(1. - ts_bondorder.fbnd[3]);
                        if( ts_bondorder.fbnd[4] > ts_bondorder.fbnd[3] )
                                *v2 = 30.*(1. - ts_bondorder.fbnd[4]);
                }
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 51 49 50 53") == 0 )
        {
                if( ts_bondorder.fbnd[5L] > 0.0F && ts_bondorder.fbnd[6L] > 0.0F )
                {
                        /* add 3/18/88 */
                        *v2 = 30.*(1. - ts_bondorder.fbnd[5]);
                        if( ts_bondorder.fbnd[6] > ts_bondorder.fbnd[5] )
                                *v2 = 30.*(1. - ts_bondorder.fbnd[6]);
                }
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 49 49 50 55") == 0 )
        {
                if( ts_bondorder.fbnd[0L] > 0.0F && ts_bondorder.fbnd[10L] > 0.0F )
                {
                        *v2 = 30.*(1. - ts_bondorder.fbnd[0]);
                        if( ts_bondorder.fbnd[10] > ts_bondorder.fbnd[0] )
                                *v2 = 30.*(1. - ts_bondorder.fbnd[10]);
                }       
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 45 49 50 55") == 0 )
        {
                if( ts_bondorder.fbnd[3L] > 0.0F && ts_bondorder.fbnd[10L] > 0.0F )
                {
                        *v2 = 30.*(1. - ts_bondorder.fbnd[3]);
                        if( ts_bondorder.fbnd[10] > ts_bondorder.fbnd[3] )
                                *v2 = 30.*(1. - ts_bondorder.fbnd[10]);
                }
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 51 49 50 55") == 0 )
        {
                if( ts_bondorder.fbnd[5L] > 0.0F && ts_bondorder.fbnd[10L] > 0.0F )
                {
                        *v2 = 30.*(1. - ts_bondorder.fbnd[5]);
                        if( ts_bondorder.fbnd[10] > ts_bondorder.fbnd[5] )
                                *v2 = 30.*(1. - ts_bondorder.fbnd[10]);
                }
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 45  1 29 50") == 0 || strcmp(kk," 45  1 29 29") == 0 )
        {
                *v2 = -30.*(1. - ts_bondorder.fbnd[2]);
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 49  1 45 29") == 0 )
        {
                *v2 = 0.;
                *v3 = -.5;
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
        else if( strcmp(kk," 45 49 50 50") == 0 )
        {
                *v2 = 30.*(1. - ts_bondorder.fbnd[1]);
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }       
        else if( strcmp(kk," 29 29 49 45") == 0 )
        {
                *v2 = -30.*(1. - ts_bondorder.fbnd[3]);
                //*ktcon_ent = kk;
                *ierr = 1;
                return;
        }
// general 
        if( strcmp(k1,"  0  2 49 49") == 0 || strcmp(k1,"  0  2 50 50") == 0 )
        {
                *v1 = .2;
                *v2 = -3.5;
                //*ktcon_ent = k1;
                *ierr = 1;
                return;
        }
        else if( strcmp(k3,"  0  2 49 49") == 0 || strcmp(k3,"  0  2 50 50") == 0 )
        {
                *v1 = .2;
                *v2 = -3.5;
                //*ktcon_ent = k3;
                *ierr = 1;
                return;
        }
        else if( strcmp(k1,"  0  3 49 49") == 0 || strcmp(k1,"  0  3 50 50") == 0 )
        {
                *v1 = .5;
                *v2 = -5.;
                //*ktcon_ent = k1;
                *ierr = 1;
                return;
        }
        else if( strcmp(k3,"  0  3 49 49") == 0 || strcmp(k3,"  0  3 50 50") == 0 )
        {
                *v1 = .5;
                *v2 = -5.;
                //*ktcon_ent = k3;
                *ierr = 1;
                return;
        }
        else if( strcmp(k1,"  0  2 49 50") == 0 || strcmp(k1,"  0  2 50 49") == 0 )
        {
                *v2 = 2.;
                //*ktcon_ent = k1;
                *ierr = 1;
                return;
        }
        else if( strcmp(k3,"  0  2 49 50") == 0 || strcmp(k3,"  0  2 50 49") == 0 )
        {
                *v2 = 2.;
                //*ktcon_ent = k3;
                *ierr = 1;
                return;
        }
        else if( strcmp(k1,"  0  3 49 50") == 0 || strcmp(k1,"  0  3 50 49") == 0 )
        {
                *v2 = 2.;
                //*ktcon_ent = k1;
                *ierr = 1;
                return;
        }
        else if( strcmp(k3,"  0  3 49 50") == 0 || strcmp(k3,"  0  3 50 49") == 0 )
        {
                *v2 = 3.;
                //*ktcon_ent = k3;
                *ierr = 1;
                return;
        }
        else if( strcmp(k1,"  0  2 51 29") == 0 )
        {
                *v2 = -10. + ts_bondorder.fbnd[5]*30.;
                //*ktcon_ent = k1;
                *ierr = 1;
                return;
        }
        else if( strcmp(k3,"  0  2 51 29") == 0 )
        {
                *v2 = -10. + ts_bondorder.fbnd[5]*30.;
                //*ktcon_ent = k3;
                *ierr = 1;
                return;
        }
//jjg added next 6 stmts....
        else if( strcmp(k2, " 49 49 29  0") == 0 || strcmp(k2," 50 50 29  0")  == 0 )
        {
                *v2 = 0.;
         //       *ktcon_ent = k2;
                *ierr = 1;
                return;
        }
        else if( strcmp(k4, " 49 49 29  0") == 0 || strcmp(k4," 50 50 29  0")  == 0 )
        {
                *v2 = 0.;
         //       *ktcon_ent = k4;
                *ierr = 1;
                return;
        }
        else if( strcmp(k5, "  0 49 29  0")  == 0 )
        {
             if(ts_bondorder.fbnd[0]!=0.)
               { *v2 = 5.0*ts_bondorder.fbnd[0];
        //        *ktcon_ent = k5;
                *ierr = 1;
                }
                return;
        }
        else if( strcmp(k5, "  0 50 29  0") == 0 )
        {
             if(ts_bondorder.fbnd[1]!=0.)
               { *v2 = 5.0*ts_bondorder.fbnd[1];
      //          *ktcon_ent = k5;
                *ierr = 1;
                }
                return;
        }
        else if( strcmp(k6, "  0 49 29  0") == 0 )
        {
             if(ts_bondorder.fbnd[0]!=0.)
               { *v2 = 5.0*ts_bondorder.fbnd[0];
    //            *ktcon_ent = k6;
                *ierr = 1;
                }
                return;
        }
        else if( strcmp(k6, "  0 50 29  0" ) == 0 )
        {
             if(ts_bondorder.fbnd[1]!=0.)
               { *v2 = 5.0*ts_bondorder.fbnd[1];
  //              *ktcon_ent = k6;
                *ierr = 1;
                }
                return;
        }
        //end jjg
//
        if( strcmp(k6,"  0 49 50  0") == 0 )
        {
                *v2 = 0.;
                //*ktcon_ent = k6;
                *ierr = 1;
                return;
        }
        else if( strcmp(k6,"  0 50 49  0") == 0 )
        {
                *v2 = 0.;
                //*ktcon_ent = k6;
                *ierr = 1;
                return;
        }
        else if( strcmp(k5,"  0 49 49  0") == 0 || strcmp(k6,"  0 49 49  0") == 0 )
        {
 /*               *v1 = .01 + .2*ts_bondorder.fbnd[0];
                *v2 = .01 + .27*ts_bondorder.fbnd[0];*/
                *v3 = .01 + .5*ts_bondorder.fbnd[0];
                //*ktcon_ent = k5;
                *ierr = 1;
                return;
        }
        else if( strcmp(k5,"  0 50 50  0") == 0 || strcmp(k6,"  0 50 50  0") == 0 )
        {
/*                *v1 = .01 + .2*ts_bondorder.fbnd[1];
                *v2 = .01 + .27*ts_bondorder.fbnd[1]; */
                *v3 = .01 + .5*ts_bondorder.fbnd[1];
                //*ktcon_ent = k5;
                *ierr = 1;
                return;
        }
        else if( strcmp(k5,"  0 29 29  0") == 0 || strcmp(k6,"  0 29 29  0") == 0 )
        {
                *v2 = 8.0;
                //*ktcon_ent = k5;
                *ierr = 1;
                return;
        }
        return;
}

