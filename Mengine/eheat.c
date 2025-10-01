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
#include "pimatx.h"
#include "rings.h"
#include "utility.h"
#include "torsions.h"
#include "gmmx.h"
#include "pistuf.h"
#include "eheat.h"
#include "pisave.h"

static float hfo = 0.0;
static float si = 0.0;  
static int iunk = 0;  
static int ilw;

float get_heat_of_formation()
{
  return hfo;
}

void  eheat()
{
        int i, iatj, iatom, ico, ii, ia,ib,ic, nb,
         iit, ilogary, imet, isp, it, iti, itj, itk, itrans, 
         j,  k, kt, l, m, mask,  nbr, nc, ncl, nnplus,
         neob, nf, nh, nio, nj, nk, nni, nnl, noppi, nox, 
         nrad, ns, nsi, numoh;
        int ifound, center, icount;
        int jatm, katm;
        int  *nstruc, *nbond;
        float be, bent, ent,tent,ssls,stot, pfc, sbe, sicprn,  sstr, str, 
          strless;
        char pt[13],pa[4],pb[4],pc[4],pd[4];

        /*  calculates heats of formation and entropies
         *      hf = energy + sum bond enthalpy + sum structural enthalpy
         *             + partition function contribution */

        /*     nn & ee arrays refer to names and values of bond enthalpy increments.
         *     mm & ss arrays refer to structural names and values.
         *
         *jjg npara is # bond parameters; ilogary is flag that pi atoms are to be
         *    calculated by ppp;  */

        nstruc = ivector(0,ehpara.nevhf);
        nbond = ivector(0,ehpara.neheat);

        ilogary = 0;
        imet = 0;
        iunk = 0;
        be = 0.;
        bent = 0.;
        ilw = 0;
        sicprn = 0.;
/*  mask for pi atoms */
    mask = 1L << PI_MASK ;

        for( i = 0; i < ehpara.neheat; i++ )
                nbond[i] = 0;
        for( i = 0; i < ehpara.nevhf; i++ )
                nstruc[i] = 0;

       if ( !gmmx.run )    
        fprintf(pcmoutfile,"Heat of Formation, Strain Energies and Entropies at 300 k (units are kcal or eu.)\n");
        /*-----------------------------------------------------------------------
         *  if pi bonds involved for which lo whitehead paramaterization is inadequate:
         *  find bond types reset ilw flag to empirical scheme based soley on bos. */
        noppi = 0;
        if( pistuf.norb > 0 && ilw == 0 )
        {
                noppi = 1;
                for( i = 0; i < pistuf.norb; i++ )
                {
                        iatom = pistuf.listpi[i];
                        iti = atom[iatom].type;
                        for( j = 0; j < MAXIAT; j++ )
                        {
                                if( atom[iatom].iat[j] != 0 )
                                {
                                        iatj = atom[iatom].iat[j];
                                        if( ( atom[iatj].flags & mask ) )
                                        {
                                                itj = atom[iatj].type;
                                                if( iti <= itj )
                                                        itk = iti*100 + itj;
                                                if( iti > itj )
                                                        itk = itj*100 + iti;
                                                if( itk == 630 || itk == 930 || itk == 241 || itk == 741
                                                     || itk == 4142 || itk == 242 || itk == 342 || itk == 246
                                                     || itk ==  348 || itk == 366 )
                                                {
                                                        noppi = 0;
                                                        ilw = 1;
                                                }
                                        }
                                }
                        }
                }
        }
      if ( !gmmx.run )
      {    
        fprintf(pcmoutfile,"Bond Enthalpy (be) and Entropy:\n"); 
        fprintf(pcmoutfile,"#   Bond or Structure          Normal     Total      Strainless Tot    Tot S contrib.\n");
      }
/*  find the enthalpy increments for each bond in the molecule */

        for( i = 1; i <= natom; i++ )
        {
            it = atom[i].type;
            if( it != 20 )
            {
               for( ii = 0; ii < MAXIAT; ii++ )
               {
                   if( atom[i].iat[ii] != 0 && atom[i].bo[ii] != 9 )
                   {
                       k = atom[i].iat[ii];
                       kt = atom[k].type;
                       if( k <= i )
                       {
                           if( kt != 20 )
                           {
                                if( it <= kt )
                                {
                                    kt += it*100;
                                }else
                                {
                                    kt = kt*100 + it;
                                }
                                if( ( atom[i].flags & mask ) &&  ( atom[k].flags & mask ) )
                                    ilogary = 1;
                                else
                                {
                                   /*  i-k is a bond and kt is the packed bond type.  search the list
                                    *  of bond types (nn(l) as listed in data stmt) to find enthalpy incrm. */
                                   ifound = FALSE;
                                   for( l = 0; l < ehpara.neheat; l++ )
                                   {
                                     nnl = ehpara.nn[l];
                                     if( nnl != 0 && nnl == kt )
                                     {
                                        nbond[l] += 1;
                                        ifound = TRUE;
                                     }
                                   }
                                   if (ifound == FALSE)
                                   {
                                      if ( !gmmx.run )  
                                         fprintf(pcmoutfile,"Error - bond %d %d does not have programmed enthalpy increments\n",i,k);
                                      iunk = 1;
                                   }
                                }
                           }
                       }
                   }
               }
            }
        }
        /*  identify structural features.  these are coded by numbers-see data */
        for( i = 1; i <= natom; i++ )
        {
           nc = 0;
           nh = 0;
           ns = 0;
           nox = 0;
           nni = 0;
           nf = 0;
           ncl = 0;
           nbr = 0;
           nio = 0;
           nsi = 0;
           ico = 0;
           isp = 0;
           neob = 0;
           nrad = 0;
           numoh = 0;
           ifound = FALSE;
           it = atom[i].type;
//
           if (it == 1 || it == 2 || it == 4 || it == 6 || it == 8 || it == 22 ||
                    it == 30 || it == 56 || it == 57)
           {                            
                nnplus = FALSE;
                for( j = 0; j < MAXIAT; j++ )
                {
                     if( atom[i].iat[j] != 0 && atom[i].bo[j] != 9 )
                     {
                          ii = atom[i].iat[j];
                          if (atom[ii].type == 41)
                             nnplus = TRUE;
                          if( !(( atom[i].flags & mask ) && ( atom[ii].flags & mask )) )
                          {
                               iit = atom[ii].type;
                               if( iit == 56 || iit == 57)
                               {
                                  nc += 1;
                                  if (iit == 57) neob += 1;
                               }
                               if (iit == 22)
                                  nc += 1;
                               if( iit == 36 )
                                   iit = 5;
                               if( iit == 21 )
                                    numoh += 1;
                               if( iit == 29 )
                                    nrad += 1;
                               if( iit >= 1 && iit <= 4 )
                               {
                                   nc += 1;
                                   if( iit == 3 )
                                       ico = 1;
                                   if( iit == 4 )
                                       isp = 1;
                                   if( iit == 2)
                                       neob += 1;
                               }
                               if( iit == 5 || iit == 36)
                                     nh += 1;

                               if( iit == 6 )
                                     nox += 1;
                               if( iit == 8 )
                                     nni += 1;
                               if( iit == 11 )
                                     nf += 1;
                               if( iit == 12 )
                                     ncl += 1;
                               if( iit == 13 )
                                     nbr += 1;
                               if( iit == 14 )
                                     nio += 1;
                               if( iit == 15 )
                                     ns += 1;
                               if( iit == 19 )
                                     nsi += 1;
                          }
                     }
                }
                
                /*  the following tests are for all structural features based on csp3
                 *  environment.  np is structural prmtr index. */
              center = -1;
              if( it == 8 )  // sec and tert amines
              {
                 if (nc == 2)
                 {
                     center = 8;
                     numeral(0,pa,3);
                     numeral(1,pb,3);
                     numeral(1,pc,3);
                     numeral(23,pd,3);
                     strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 } else if (nc == 3)
                 {
                     center = 8;
                     numeral(0,pa,3);
                     numeral(1,pb,3);
                     numeral(1,pc,3);
                     numeral(1,pd,3);
                     strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 }
              }
              if( it == 30 ) // sec and tert cations
              {
                  if (nc == 2)
                  {
                     center = 30;
                     numeral(0,pa,3);
                     numeral(1,pb,3);
                     numeral(1,pc,3);
                     numeral(5,pd,3);
                     strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                  } else if (nc == 3)
                  {
                     center = 30;
                     numeral(0,pa,3);
                     numeral(1,pb,3);
                     numeral(1,pc,3);
                     numeral(1,pd,3);
                     strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                  }
              }
              if( it == 2 || it == 57)
              {
                   if (nc == 3)
                   {
                     center = 2;
                     numeral(0,pa,3);
                     numeral(1,pb,3);
                     numeral(1,pc,3);
                     numeral(2,pd,3);
                     strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                   }
                   if (nnplus == TRUE)
                   {
                     center = 2;
                     numeral(0,pa,3);
                     numeral(2,pb,3);
                     numeral(2,pc,3);
                     numeral(41,pd,3);
                     strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                   }
              }
              if (it == 1 || it == 56 ) // || it == 22)
              {
                    if (nox == 1 && nh == 3)  // o-methyl
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if ( nox == 1 && nc == 2)  // o-sec carbon ether & alcohol
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(5,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( (nox == 1) && (nc == 3) ) // o-tert
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nox == 2 && nh == 2 )   // geminal
                    { 
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nox == 2 && nc == 1 && nh == 1 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(5,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nox == 2 && nc == 2 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nox == 3 && nh == 1 )
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(6,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nox == 1 && ncl == 1 )  // general anomeric
                    {
                       center = 1;
                       numeral(0,pa,3);
                       numeral(0,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                     /*  gen. geminal for anomeric */
                    }else if( nox >= 1 && nni >= 1 ) // general anomeric
                    {
                       center = 1;
                       numeral(0,pa,3);
                       numeral(0,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nox >= 1 && nf >= 1 ) // general anomeric
                    {
                       center = 1;
                       numeral(0,pa,3);
                       numeral(0,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nox >= 1 && ns >= 1 ) // general anomeric
                    {
                       center = 1;
                       numeral(0,pa,3);
                       numeral(0,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                        /*  amines */
                    }else if( nni >= 1 && nf >= 1 )  // amines  general anomeric
                    {
                       center = 1;
                       numeral(0,pa,3);
                       numeral(0,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( (nni == 1) && (nh == 3) ) // general anomeric
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(8,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( (nni >= 1) && (nc == 2) )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(5,pc,3);
                       numeral(8,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( (nni == 1) && (nc == 3) )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(8,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nf >= 2 ) // florides
                    {
                       center = 1;
                       numeral(0,pa,3);
                       numeral(0,pb,3);
                       numeral(11,pc,3);
                       numeral(11,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nf == 1 && nh == 3 )
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(11,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nf == 1 && nc == 2 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(5,pc,3);
                       numeral(11,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nf == 1 && nc == 3 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(11,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( ncl ==1 && nh == 3 )  // chlorides
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(12,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( ncl == 1 && nc == 2 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(5,pc,3);
                       numeral(12,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( ncl == 1 && nc == 3 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(12,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nbr ==1 && nh == 3 )  // bromides
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(13,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nbr == 1 && nc == 2 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(5,pc,3);
                       numeral(13,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nbr == 1 && nc == 3 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(13,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nio ==1 && nh == 3 )  // iodides
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(14,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nio == 1 && nc == 2 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(5,pc,3);
                       numeral(14,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nio == 1 && nc == 3 )
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(14,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if ( ns >= 2)
                    {
                       center = 1;
                       numeral(0,pa,3);
                       numeral(0,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if (ns >=1 && nf >= 1)
                    {
                       center = 1;
                       numeral(0,pa,3);
                       numeral(0,pb,3);
                       numeral(6,pc,3);
                       numeral(6,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if (ns == 1 && nh == 3)
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(15,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if ( ns >=1 && nc == 2)
                    {
                        if (ns == 1)
                        {
                            center = 1;
                            numeral(1,pa,3);
                            numeral(1,pb,3);
                            numeral(5,pc,3);
                            numeral(15,pd,3);
                            strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                        } else
                        {
                            center = 1;
                            numeral(1,pa,3);
                            numeral(1,pb,3);
                            numeral(1,pc,3);
                            numeral(15,pd,3);
                            strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                        }
                    }else if (ns == 1 && nc == 3)
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(15,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if (nsi ==1 && nh == 3)
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(19,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if (nrad == 1 && nh == 3)
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(29,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if (nc == 3 && nh == 1)  // np 1
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(5,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if (nc == 4)       // np 2
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(1,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nc == 1 && nh == 3 && atom[atom[i].iat[0]].type == 22 ) /*  c(sp3)-methyl */
                    {
                       center = 1;
                       numeral(5,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(22,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nc == 1 && nh == 3 && isp == 0 ) /*  c(sp3)-methyl */
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(5,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( nc == 4 && neob == 1 )           /*  c(sp2)-tert-c(sp3) */
                    {
                       center = 1;
                       numeral(1,pa,3);
                       numeral(1,pb,3);
                       numeral(1,pc,3);
                       numeral(2,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if( (nc == 1 && nh == 3) && isp != 0 ) /*  c(sp)-methyl */
                    {
                       center = 1;
                       numeral(4,pa,3);
                       numeral(5,pb,3);
                       numeral(5,pc,3);
                       numeral(5,pd,3);
                       strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    }else if (ico != 0)
                    {
                        if (nc ==1 && nh == 3)  // me carbonyl
                        {
                           center = 301;
                           numeral(1,pa,3);
                           numeral(5,pb,3);
                           numeral(5,pc,3);
                           numeral(5,pd,3);
                           strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                        }else if (nc == 3 && nh == 1) // sec carbonyl
                        {
                           center = 301;
                           numeral(1,pa,3);
                           numeral(1,pb,3);
                           numeral(1,pc,3);
                           numeral(5,pd,3);
                           strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                        }else if (nc == 4)  // tert carbonyl
                        {
                           center = 301;
                           numeral(1,pa,3);
                           numeral(1,pb,3);
                           numeral(1,pc,3);
                           numeral(1,pd,3);
                           strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                        }
                    }
              }
              if (center != -1)
              {
                 for (l=0; l < ehpara.nevhf; l++)
                 {
                     if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                        nstruc[l] += 1;
                 }                      
              }
           }
        }  // end of i loop
   // finished atom based structural types
   // correct for small ringss
   // correct for sp2 in cyclopropanes
   // find specific functional groups
                /*jjg trigonal centers in cyclopropane-to correct for external angle strain-
                 *    but find cyclopropenes to ignore trig c corr. & sub 8.0 for c3enes.
                 *    for c=n in 3-ring sub 8.0/2. */

         if (rings.nring3 != 0)
         {
             center = 0;
             numeral(0,pa,3);
             numeral(22,pb,3);
             numeral(22,pc,3);
             numeral(22,pd,3);
             strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
             for (l=0; l < ehpara.nevhf; l++)
             {
                 if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                     nstruc[l] += rings.nring3;
             }                      
         }
         if (rings.nring4 != 0)
         {
             // read ring4 correction
             center = 0;
             numeral(56,pa,3);
             numeral(56,pb,3);
             numeral(56,pc,3);
             numeral(56,pd,3);
             strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
             for (l=0; l < ehpara.nevhf; l++)
             {
                 if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                     nstruc[l] += rings.nring4;
             }                      
         }
         if (rings.nring5 != 0)
         {
             // read ring5 correction
             center = 0;
             numeral(1,pa,3);
             numeral(1,pb,3);
             numeral(1,pc,3);
             numeral(1,pd,3);
             strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
             for (l=0; l < ehpara.nevhf; l++)
             {
                 if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                     nstruc[l] += rings.nring5;
             }                      
         }
   // correct for sp2 in c3 and cyclopropene
         if (rings.nring3 != 0)
         {
             for (i=0; i < rings.nring3; i++)
             {
                 nj = 0;
                 nk = 0;
                 for(j=0; j < 3; j++)
                 {
                     if (atom[rings.r13[i][j]].type == 2)
                       nj++;
                     if (atom[rings.r13[i][j]].type == 37)
                       nk++;
                 }
                 if (nk > 0)  // c=n in three membered ring
                 {
                    center = 0;
                    numeral(0,pa,3);
                    numeral(2,pb,3);
                    numeral(22,pc,3);
                    numeral(37,pd,3);
                    strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                    for (l=0; l < ehpara.nevhf; l++)
                    {
                         if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                              nstruc[l] += 1;
                    }                      
                 }
                 if (nj > 0 && nk == 0)
                 {
                     if (nj == 1)  // exocyclic
                     {
                        center = 0;
                        numeral(0,pa,3);
                        numeral(2,pb,3);
                        numeral(22,pc,3);
                        numeral(22,pd,3);
                        strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                        for (l=0; l < ehpara.nevhf; l++)
                        {
                            if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                              nstruc[l] += 1;
                        }                      
                     } else if (nj == 2 || nj == 3)
                     {
                           ia = rings.r13[i][0];
                           ib = rings.r13[i][1];
                           ic = rings.r13[i][2];
                           nk = 0;
                           for (j=0; j < MAXIAT; j++)
                           {
                               if (atom[ia].iat[j] != 0)
                               {
                                   if (atom[ia].bo[j] == 2)
                                   {
                                       if (atom[ia].iat[j] == ib || atom[ia].iat[j] == ic)
                                         nk++;
                                   }
                               }
                               if (atom[ib].iat[j] != 0)
                               {
                                   if (atom[ib].bo[j] == 2)
                                   {
                                       if (atom[ia].iat[j] == ic)
                                         nk++;
                                   }
                               }
                           }
                           if (nj == 2 && nk == 1) // in ring
                           {
                               center = 0;
                               numeral(0,pa,3);
                               numeral(2,pb,3);
                               numeral(2,pc,3);
                               numeral(22,pd,3);
                               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                               for (l=0; l < ehpara.nevhf; l++)
                               {
                                   if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                     nstruc[l] += 1;
                               }                      
                           } else if (nj == 2 && nk == 0)
                           {
                               center = 0;
                               numeral(0,pa,3);
                               numeral(2,pb,3);
                               numeral(22,pc,3);
                               numeral(22,pd,3);
                               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                               for (l=0; l < ehpara.nevhf; l++)
                               {
                                   if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                     nstruc[l] += 2;
                               }                      
                           } else if (nj == 3 && nk == 1)
                           {
                               center = 0;
                               numeral(0,pa,3);
                               numeral(2,pb,3);
                               numeral(2,pc,3);
                               numeral(22,pd,3);
                               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                               for (l=0; l < ehpara.nevhf; l++)
                               {
                                   if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                     nstruc[l] += 1;
                               }                      
                               center = 0;
                               numeral(0,pa,3);
                               numeral(2,pb,3);
                               numeral(22,pc,3);
                               numeral(22,pd,3);
                               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                               for (l=0; l < ehpara.nevhf; l++)
                               {
                                   if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                     nstruc[l] += 1;
                               }                      
                           } else if (nj == 3 && nk == 0)
                           {
                               center = 0;
                               numeral(0,pa,3);
                               numeral(2,pb,3);
                               numeral(22,pc,3);
                               numeral(22,pd,3);
                               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                               for (l=0; l < ehpara.nevhf; l++)
                               {
                                   if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                     nstruc[l] += 3;
                               }                      
                           }
                     }
                 }
             }
         }
//  look for c=N+
         
//  bridgehead carbons in bicyclobutane
         if (rings.nring3 > 1)
         {
             ifound = FALSE;
             icount = 0;
             for (i=0; i < rings.nring3-1; i++)
             {
                 ia = rings.r13[i][0];
                 ib = rings.r13[i][1];
                 ic = rings.r13[i][2];
                 for (j=i+1; j < rings.nring3; j++)
                 {
                     for (k=0; k < 3; k++)
                     {
                         if (ia == rings.r13[j][k])
                         {
                             icount++;
                             ifound = TRUE;
                         }
                         if (ib == rings.r13[j][k])
                         {
                             icount++;
                             ifound = TRUE;
                         }
                         if (ic == rings.r13[j][k])
                         {
                             icount++;
                             ifound = TRUE;
                         }
                     }
                 }
             }
             if (ifound == TRUE)
             {
                  center = 22;
                  numeral(22,pa,3);
                  numeral(22,pb,3);
                  numeral(22,pc,3);
                  numeral(22,pd,3);
                  strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                  for (l=0; l < ehpara.nevhf; l++)
                  {
                       if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                           nstruc[l] += icount;
                  }                     
             }
         }              
  /*  o-c-c-o, n-c-c-n , o-c-c-n and c-c-c=c */
        for (i=1; i <= natom; i++)
        {
            if (atom[i].type == 1 || atom[i].type == 56)  // 1-1-2-2
            {
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 &&
                       (atom[atom[i].iat[j]].type == 1 || atom[atom[i].iat[j]].type == 56))
                    {
                        jatm = atom[i].iat[j];
                        for (k=0; k < MAXIAT; k++)
                        {
                            if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9)
                            {
                                if (atom[jatm].iat[k] != i && atom[atom[jatm].iat[k]].type == 2)
                                {
                                    katm = atom[jatm].iat[k];
                                    for (m=0; m < MAXIAT; m++)
                                    {
                                        if (atom[katm].iat[m] != 0 && atom[katm].bo[m] != 9)
                                        {
                                            if (atom[atom[katm].iat[m]].type == 2) //  && i < atom[katm].iat[m])
                                            {
                                                center = 0;
                                                numeral(1,pa,3);
                                                numeral(1,pb,3);
                                                numeral(2,pc,3);
                                                numeral(2,pd,3);
                                                strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                                                for (l=0; l < ehpara.nevhf; l++)
                                                {
                                                     if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                                         nstruc[l] += 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }                                                       
            if (atom[i].type == 6)  // 6-1-1-6 & 6-1-1-8
            {
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].type == 1)
                    {
                        jatm = atom[i].iat[j];
                        for (k=0; k < MAXIAT; k++)
                        {
                            if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9)
                            {
                                if (atom[jatm].iat[k] != i && atom[atom[jatm].iat[k]].type == 1)
                                {
                                    katm = atom[jatm].iat[k];
                                    for (m=0; m < MAXIAT; m++)
                                    {
                                        if (atom[katm].iat[m] != 0 && atom[katm].bo[m] != 9)
                                        {
                                            if (atom[atom[katm].iat[m]].type == 6 && i < atom[katm].iat[m])
                                            {
                                                center = 0;
                                                numeral(6,pa,3);
                                                numeral(1,pb,3);
                                                numeral(1,pc,3);
                                                numeral(6,pd,3);
                                                strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                                                for (l=0; l < ehpara.nevhf; l++)
                                                {
                                                     if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                                         nstruc[l] += 1;
                                                }
                                            }                    
                                            if (atom[atom[katm].iat[m]].type == 8 && i < atom[katm].iat[m])
                                            {
                                                center = 0;
                                                numeral(6,pa,3);
                                                numeral(1,pb,3);
                                                numeral(1,pc,3);
                                                numeral(8,pd,3);
                                                strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                                                for (l=0; l < ehpara.nevhf; l++)
                                                {
                                                     if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                                         nstruc[l] += 1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } //
            if (atom[i].type == 8)  // 8-1-1-8
            {
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].type == 1)
                    {
                        jatm = atom[i].iat[j];
                        for (k=0; k < MAXIAT; k++)
                        {
                            if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9)
                            {
                                if (atom[jatm].iat[k] != i && atom[atom[jatm].iat[k]].type == 1)
                                {
                                    katm = atom[jatm].iat[k];
                                    for (m=0; m < MAXIAT; m++)
                                    {
                                        if (atom[katm].iat[m] != 0 && atom[katm].iat[m] != 9)
                                        {
                                            if (atom[atom[katm].iat[m]].type == 8 && i < atom[katm].iat[m])
                                            {
                                                center = 0;
                                                numeral(8,pa,3);
                                                numeral(1,pb,3);
                                                numeral(1,pc,3);
                                                numeral(8,pd,3);
                                                strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                                                for (l=0; l < ehpara.nevhf; l++)
                                                {
                                                     if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                                         nstruc[l] += 1;
                                                }
                                            }                    
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } //
            if (atom[i].type == 15)   // 15-1-1-15
            {
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].type == 1)
                    {
                        jatm = atom[i].iat[j];
                        for (k=0; k < MAXIAT; k++)
                        {
                            if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9)
                            {
                                if (atom[jatm].iat[k] != i && atom[atom[jatm].iat[k]].type == 1)
                                {
                                    katm = atom[jatm].iat[k];
                                    for (m=0; m < MAXIAT; m++)
                                    {
                                        if (atom[katm].iat[m] != 0 && atom[katm].iat[m] != 9)
                                        {
                                            if (atom[atom[katm].iat[m]].type == 15 && i < atom[katm].iat[m])
                                            {
                                                center = 0;
                                                numeral(15,pa,3);
                                                numeral(1,pb,3);
                                                numeral(1,pc,3);
                                                numeral(15,pd,3);
                                                strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                                                for (l=0; l < ehpara.nevhf; l++)
                                                {
                                                     if (center == ehpara.cent[l] && strcmp(pt,ehpara.kheat[l]) == 0)
                                                         nstruc[l] += 1;
                                                }
                                            }                    
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } //
        }
         be = 0.0;
         sbe = 0.0;
         /*  print increment sums for bonds */
         for( i = 0; i < ehpara.neheat; i++ )
         {
                nb = nbond[i];
                if( nb != 0 )
                {
                        if( nb > 1000 )
                        {
                                /*  unknown bonds appear here */
                                it = ehpara.nn[i]/100;
                                kt = ehpara.nn[i] - it*100;
                                nb -= 1000;
                               if ( !gmmx.run )
                                fprintf(pcmoutfile,"Bond %d Type %d - %d is Unknown\n",nb,it,kt);
                                iunk = 1;
                        }
                        else
                        {
                                str = ehpara.ee[i];
                                ent = ehpara.et[i];
                                strless = ehpara.estr[i];
                                sstr = nb*str;
                                tent = nb*ent;
                                ssls = nb*strless;
                                if( ent == 999. )
                                        tent = 0.;
                                bent += tent;
                                be += sstr;
                                sbe += ssls;
                               if ( !gmmx.run )   
                                fprintf(pcmoutfile,"%-4d %s %11.3f  %9.3f  %9.3f %9.3f   %9.3f\n",nb,ehpara.cc[i],str,sstr,strless,ssls,tent); 
                        }
                }
         }
//
        /*  entropy contributions from rings 3 through 8-no symmetry corrections */
        stot = bent + rings.nring3*31.5 + rings.nring4*29.0;
        /*  print increment sums for structural parameters */
        for( i = 0; i < ehpara.nevhf; i++ )
        {
                ns = nstruc[i];
                str = ehpara.ss[i];
                if( ns != 0 )
                {
                        sstr = ns*str;
                        be += sstr;
                       if ( !gmmx.run )   
                        fprintf(pcmoutfile,"%-4d %s %11.3f  %9.3f \n",ns,ehpara.mm[i],str,sstr);
                }
        }
       if ( !gmmx.run )
       {    
        fprintf(pcmoutfile,"\n                             BE=%9.3f  S=%9.3f\n\n",be,stot);
        fprintf(pcmoutfile," 3 & 4 Ring corrections to entropy are included w/o symmetry corrections.\n");
        fprintf(pcmoutfile,"   for each 5-ring add 26 eu.; for each 6 &7-ring add 16 eu.\n");
        fprintf(pcmoutfile,"   for each 8-ring add 14 eu.; for higher rings add 12 eu. each.\n");
        fprintf(pcmoutfile,"   there are no symmetry corrections to the entropy.\n\n");
       }
        pfc = 2.4;
        hfo = energies.total + be + pfc;
        si = energies.total + be - sbe;
        itrans = 0;
        if( ilogary != 1 )
        {
          if ( !gmmx.run )
          {    
                if( ilogary == 0 )
                {
                        fprintf(pcmoutfile," Heat of Formation Calculation\n");
                        fprintf(pcmoutfile," Partition Function Contribution (PFC):\n");
                        fprintf(pcmoutfile," Conformational Population Increment (POP): %7.3f\n",0.0);
                        fprintf(pcmoutfile," Torsional Contribution (TOR)             : %7.3f\n",0.0);
                        fprintf(pcmoutfile," Translation/Rotation Term (T/R)          : %7.3f\n",2.4);
                        fprintf(pcmoutfile,"                                PFC  =    : %7.3f\n",pfc);
                }
                fprintf(pcmoutfile,"\n Heat of Formation (hf0) = energy + be + pfc  : %16.3f\n",hfo);
           }
        }
       if ( !gmmx.run )   
        fprintf(pcmoutfile,"\n Strain Energy (energy+environment corrs.)=   : %16.3f\n",si);
        if( ilogary == 1 )
                eheatpi( &be, &sbe, &ssls, &sstr, &imet );
        if( iunk == 1 )
        {
               if ( !gmmx.run )    
                fprintf(pcmoutfile,"CAUTION - delta hf is not correct because of missing parameters\n");
        }
        free_ivector(nbond , 0,ehpara.neheat);
        free_ivector(nstruc, 0, ehpara.nevhf);
        return;
} 
/* ------------------------------ */
void eheatpi(float *be,float *sbe,float *ssls,float *sstr,int *imet)
{
        int i, iatm1, iatm2, ihydc, ii, iip, ikl, iik,
         ilp, imult, inum, ip, it, itiat, itin, itkn, k, kbo, kl,  
         kmatom, knum, kp, kt, l, nb, ncbnyl, nclpi, nimmo, 
         nnl, noazo, nocarb, nohyd, noimin, nonit, nooxy, nopi, nosp, 
         npt[60], ntyp2, ntyp23, ntyp28, mask, nodi, notri,
         kk, jj;
        int *nep, *npbond;
        int ia,ib,ic,id;
        double ang,dihed;
        float alim, blim, ccrho, del2, delsig, disp, ealt, ecarb, ecomp, 
         epitor, eppi, esig, pbe, pfc, pibe, rho, sls;
        float *epp;

        /* iuse=1 read pi parameters, no calc
         * iuse=2 find pi parameters & output final heat of formation
         * iuse=2 is called only after complete calc in eheat
         * ilw=1 use empirical scheme to calc.
         * ilw=0 use lo whitehead - like scheme (incomplete rel to ilw=1) */
        /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         * */
	nclpi = kp = ip = 0;
	eppi = 0.0;
        mask = 1L << 0;

        nep = ivector(0,200);
        npbond = ivector(0,200);
        epp = vector(0,200);
        
        alim = .4000;
        blim = .84;
        pbe = 0.;
        *sbe = 0.;
        pibe = 0.;
        epitor = 0.;
        ccrho = 0.;
        esig = 0.;
        ecarb = 0.;
        ecomp = 0.;
        ealt = 0.;
        imult = 0;

        for( i = 0; i < epiheat.npihf; i++ )
        {
                npbond[i] = 0;
                epp[i] = .0;
                nep[i] = 0;
                epp[i] = epiheat.eep[i];
        }
        for( i = 0; i < 60; i++ )
        {
                npt[i] = 0;
        }
        for( i = 0; i < pistuf.norb; i++ )
        {
                ilp = atom[pistuf.listpi[i]].type;
                npt[ilp] += 1;
        }
        kmatom = 0;
//
                nimmo = 0;
                noazo = 0;
                ihydc = 1;
                if( ilw == 0 )
                        goto L_99956;
                for( i = 0; i < pistuf.norb; i++ )
                {
                        inum = pistuf.listpi[i];
                        for( ii = 0; ii < MAXIAT; ii++ )
                        {
                                if( atom[inum].iat[ii] == 0 || atom[inum].bo[ii] == 9 )
                                        break;
                                knum = atom[inum].iat[ii];
                                if( !(!( atom[knum].flags & mask ) || knum > inum) )
                                {
                                        itin = atom[inum].type;
                                        itkn = atom[knum].type;
                                        kbo = atom[inum].bo[ii];
                                        if( kbo >= 2 )
                                        {
                                                imult += 1;
                                                if( itin == 2 && itkn == 2 )
                                                        nep[0] += 1;
                                                if( itin == 57 && itkn == 57 )
                                                        nep[0] += 1;
                                                if( itin == 4 && itkn == 4 )
                                                        nep[1] += 1;
                                                if( (itin == 2 && itkn == 4) || (itin == 4 && itkn == 2) )
                                                        nep[2] += 1;
                                                if( (itin == 2 && itkn == 3) || (itin == 3 && itkn == 2) )
                                                        nep[3] += 1;
                                                if( (itin == 3 && itkn == 7) || (itin == 7 && itkn == 3) )
                                                        nep[4] += 1;
                                                if( (itin == 2 && itkn == 37) || (itin == 37 && itkn == 2) )
                                                        nep[9] += 1;
                                                if( (itin == 41 && itkn == 2) || (itin == 2 && itkn == 41) )
                                                        nep[17] += 1;
                                                if( (itin == 4 && itkn == 10) || (itin == 10 && itkn == 4) )
                                                        nep[18] += 1;
                                                if( (itin == 7 && itkn == 41) || (itin == 41 && itkn == 7) ||
                                                 (itin == 41 && itkn == 66) || (itin == 66 && itkn == 41))
                                                        nep[19] += 1;
                                                if( itin == 37 && itkn == 37 )
                                                        nep[20] += 1;
                                                if( (itin == 2 && itkn == 46) || (itin == 46 && itkn == 2) )
                                                        nep[24] += 1;
                                        }
                                        else
                                        {
                                              if( itin == 2 && itkn == 2 )
                                                        nep[0] += 1;
                                                if( itin == 57 && itkn == 57 )
                                                        nep[0] += 1;
                                                if( itin == 4 && itkn == 4 )
                                                        nep[1] += 1;
                                                if( (itin == 2 && itkn == 4) || (itin == 4 && itkn == 2) )
                                                        nep[2] += 1;
                                                if( (itin == 2 && itkn == 3) || (itin == 3 && itkn == 2) )
                                                        nep[3] += 1;
                                                if( (itin == 3 && itkn == 7) || (itin == 7 && itkn == 3) )
                                                        nep[4] += 1;
                                                if( (itin == 2 && itkn == 37) || (itin == 37 && itkn == 2) )
                                                        nep[9] += 1;
                                                if( (itin == 41 && itkn == 2) || (itin == 2 && itkn == 41) )
                                                        nep[17] += 1;
                                                if( (itin == 4 && itkn == 10) || (itin == 10 && itkn == 4) )
                                                        nep[18] += 1;
                                                if( (itin == 7 && itkn == 41) || (itin == 41 && itkn == 7) ||
                                                 (itin == 41 && itkn == 66) || (itin == 66 && itkn == 41))
                                                        nep[19] += 1;
                                                if( itin == 37 && itkn == 37 )
                                                        nep[20] += 1;
                                                if( (itin == 2 && itkn == 46) || (itin == 46 && itkn == 2) )
                                                        nep[24] += 1;
//
                                                if( (itin == 2 && itkn == 9) || (itin == 9 && itkn == 2) )
                                                        nep[5] += 1;
                                                if( (itin == 2 && itkn == 29) || (itin == 29 && itkn == 2) )
                                                        nep[6] += 1;
                                                if( (itin == 2 && itkn == 30) || (itin == 30 && itkn == 2) )
                                                        nep[7] += 1;
                                                if( (itin == 2 && itkn == 48) || (itin == 48 && itkn == 2) )
                                                        nep[8] += 1;
                                                if( (itin == 2 && itkn == 6) || (itin == 6 && itkn == 2) )
                                                        nep[10] += 1;
                                                if( (itin == 6 && itkn == 30) || (itin == 30 && itkn == 6) )
                                                        nep[11] += 1;
                                                if( (itin == 9 && itkn == 30) || (itin == 30 && itkn == 9) )
                                                        nep[12] += 1;
                                                if( (itin == 3 && itkn == 9) || (itin == 9 && itkn == 3) )
                                                        nep[13] += 1;
                                                if( (itin == 3 && itkn == 6) || (itin == 6 && itkn == 3) )
                                                        nep[14] += 1;
                                                if( (itin == 3 && itkn == 48) || (itin == 48 && itkn == 3) )
                                                        nep[15] += 1;
                                                if( (itin == 42 && itkn == 2) || (itin == 2 && itkn == 42) )
                                                        nep[16] += 1;
                                                if( itin == 37 && itkn == 37 )
                                                {
                                                        nep[20] += 1;
                                                        noazo += 1;
                                                }
                                                if( (itin == 3 && itkn == 37) || (itin == 37 && itkn == 3) )
                                                        nep[21] += 1;
                                                if( (itin == 41 && itkn == 42) || (itin == 42 && itkn == 41) ||
                                                    (itin == 41 && itkn == 66) || (itin == 66 && itkn == 41) )
                                                        nep[22] += 1;
                                                if( (itin == 3 && itkn == 42) || (itin == 42 && itkn == 3) )
                                                        nep[23] += 1;
                                        }
                                }
                        }
                }
                if( imult > 0 )
                {
                      nep[0] /= 2;
                      nep[0] /= 2;
                      nep[1] /= 2;
                      nep[2] /= 2;
                      nep[3] /= 2;
                      nep[4] /= 2;
                      nep[9] /= 2;
                      nep[18] /= 2;
                      nep[19] /= 2;
                      nep[20] /= 2;
                      nep[24] /= 2;
                }
//       kmatom = 1;

        nclpi = nep[0];
        if( kmatom == 1 )
        {
                nep[1] = (npt[3] - npt[9])/2;
                if( nep[9] > npt[36] )
                        nep[9] /= 2;
                if( nep[17] > npt[40] )
                        nep[17] /= 2;
                nclpi = pimatx.npiel - nep[1] - nep[4] - nep[9] - nep[18]/2 - 
                 nep[19] - npt[5] - npt[8] - npt[28]/2 - npt[41] - npt[45] - 
                 npt[47];
                nep[0] = nclpi;
                if( (npt[36] >= 2 || npt[40] >= 2) || npt[36] + npt[40] >= 2 )
                {
                       if ( !gmmx.run )    
                        fprintf(pcmoutfile," Caution check for correct numbers of multiple bonds in lists below\n");
                }
        }
        for( i = 0; i < epiheat.npihf; i++ )
        {
                epp[i] *= nep[i];
        }
        goto L_140;
L_99956:
        ihydc = 0;

L_140:
        ;
        for( i = 1; i <= natom; i++ )
        {
                it = atom[i].type;
                if (it == 57)
                   it = 2;
                if( ( atom[i].flags & mask ) )
                {
                        if( ihydc == 0 )
                        {
                                nopi = 0;
                                nohyd = 0;
                                nocarb = 0;
                                nooxy = 0;
                                nonit = 0;
                                noimin = 0;
                                ncbnyl = 0;
                                nosp = 0;
                                ntyp2 = 0;
                                ntyp28 = 0;
                                ntyp23 = 0;
                                nodi = 0;
                                notri = 0;
                                for( kl = 0; kl < MAXIAT; kl++ )
                                {
                                        if( atom[i].iat[kl] != 0 && atom[i].bo[kl] != 9 )
                                        {
                                                ikl = atom[i].iat[kl];
                                                itiat = atom[ikl].type;
                                                if (itiat == 57)  // cyclobutene correction
                                                    itiat = 2;    
                                                if( ( atom[ikl].flags & mask ) )
                                                {
                                                        nopi += 1;
                                                        if( itiat == 2 )
                                                        {
                                                                ntyp2 += 1;
                                                                jj = 0;
                                                                for (kk = 0; kk < MAXIAT; kk++)
                                                                {
                                                                        if(atom[atom[ikl].iat[kk]].type == 2 || atom[atom[ikl].iat[kk]].type == 57 )
                                                                                jj += 1;
                                                                }
                                                                if (jj == 1)
                                                                        nodi += 1;
                                                                if (jj == 2)
                                                                        notri += 1;
                                                        }
                                                        if( itiat == 3 )
                                                                ncbnyl += 1;
                                                        if( itiat == 4 )
                                                                nosp += 1;
                                                        if( itiat == 6 )
                                                                nooxy += 1;
                                                        if( itiat == 9 )
                                                                nonit += 1;
                                                        if( itiat == 37 )
                                                                noimin += 1;
                                                }
                                                if( itiat == 1 )
                                                        nocarb += 1;
                                                if( itiat == 5 )
                                                        nohyd += 1;
                                                if( itiat == 28 )
                                                        ntyp28 += 1;
                                                if (itiat == 23)
                                                        ntyp23 += 1;
                                        }
                                }
                                if( it == 2 || it == 30 || it == 29 || it == 48 || it == 57 )
                                {
                                        if( nopi == 1 )
                                        {
                                                ecarb += 92.2;
                                                if( nocarb == 2 )
                                                {
                                                        ecarb -= 7.0;
                                                }
                                                else if( nohyd == 1 && nocarb == 1 )
                                                {
                                                        ecarb -= 3.5;
                                                }
                                        }
                                        else if( nopi == 2 )
                                        {
                                                ecarb += 141.95;
                                                if( nocarb == 1 )
                                                        ecarb -= 1.8;
                                                if( nooxy == 1 )
                                                        ecarb += 12.;
                                        }
                                        else
                                        {
                                                ecarb += 189.8;
                                                if( nonit == 1 )
                                                        ecarb += 21.;
                                                if( nooxy == 1 )
                                                        ecarb += 13.4;
                                        }
                                        if( it == 30 )
                                                ecarb += 222.;
                                        if( it == 48 )
                                                ecarb -= 230.;
                                }
                                if( it == 3 )
                                {
                                        if( nopi == 1 )
                                        {
                                                if( nohyd == 2 )
                                                        ecarb += 98.5;
                                                if( nohyd == 1 )
                                                        ecarb += 85.2 - 1.7;
                                                if( nohyd == 0 )
                                                        ecarb += 75.2;
                                        }
                                        else if( nopi == 2 )
                                        {
                                                ecarb += 145.2;
                                                if( nohyd == 0 )
                                                        ecarb -= 10.0;
                                                if( nohyd == 1 )
                                                        ecarb -= 1.7;
                                                if( nooxy == 1 )
                                                        ecarb += -23. - 1.5 + 33.8 + 17.;
                                                if( nonit == 1 )
                                                        ecarb += -3.9 + 28. + 17.9;
                                        }
                                        else
                                        {
                                                ecarb += 200.5;
                                        }
                                }
                                if( it == 4 )
                                {
                                        if( nopi == 1 )
                                        {
                                                ecarb += 119.;
                                                if( nohyd == 0 )
                                                        ecarb -= 3.4;
                                        }
                                        if( nopi == 2 )
                                        {
                                                ecarb += 179.6;
                                                if( nosp == 2 )
                                                        ecarb -= 3.85;
                                        }
                                }
                                if( it == 6 )
                                {
                                        ecarb -= 327.9;
                                        if( ntyp2 == 1 && ntyp28 == 1 )
                                                ecarb += 3.5;
                                        if( nopi == 2 )
                                                ecarb += 75.;
                                }
                                if( it == 7 )
                                        ecarb += 25.7;
                                if( it == 9 )
                                {
                                        ecarb -= 273.;
                                        if( nocarb == 2 )
                                                ecarb -= 1.3;
                                        /* for pyrrole type */
                                        if( nopi == 2 && ntyp23 == 1)
                                        {
                                                ecarb += 173.0F;                
                                                if ( nodi == 1 && notri == 1)
                                                        ecarb -= 31.0F;
                                                if ( nodi == 0 && notri == 2)
                                                        ecarb -= 58.0F;

                                        } else if (nopi == 2 && ntyp23 == 0)    // N - subst pyrrole
                                        {
                                                ecarb += 167.50F;
                                        }
                                        /* for aniline */
                                        if( nopi == 1 && ntyp2 == 1 )
                                                ecarb += 0.4;           // old value = +2.4
                                        /* for ene amides */
                                        if( ntyp2 == 1 && ncbnyl == 1 )
                                                ecarb -= 27.;
                                }
                                if( it == 10 )
                                        ecarb += 173.;
                                if( it == 37 )
                                {
                                        if( nopi == 2 )
                                        {
                                                ecarb += 128.6 + 19. + 5.7;
                                                ecarb += -8.5*noimin;
                                                ecarb += -36.*nonit;
                                        }
                                        else
                                        {
                                                ecarb += 71.8 + 16. + 5.;
                                                ecarb += -3.5*noimin;
                                                if( nohyd == 1 )
                                                        ecarb += 3.;
                                        }
                                }
                                if( it == 42 )
                                {
                                        ecarb -= 217.5;
                                        if( ncbnyl == 1 )
                                                ecarb -= 110.5;
                                }
                                if( it == 48 )
                                {
                                        ecarb += 12.*ncbnyl;
                                        ecarb += -45.*nosp;
                                }
                        }
                        for( ii = 0; ii < MAXIAT; ii++ )
                        {
                            if( atom[i].iat[ii] != 0 && atom[i].bo[ii] != 9 )
                            {
                               k = atom[i].iat[ii];
                               if( ( atom[k].flags & mask ) )
                               {
                                   if( k <= i )
                                   {
                                      kt = atom[k].type;
                                      if( ihydc == 0 )
                                      {
                                           if( kt == 29 || kt == 30 || kt == 48 || kt == 57 )
                                                  kt = 2;
                                           if( it == 29 || it == 30 || it == 48 || it == 57 )
                                                  it = 2;
                                      }
                                      if( it > kt )
                                      {
                                           kt = kt*100 + it;
                                      }else
                                      {
                                           kt += it*100;
                                      }
                                      for (iik=0; iik < pistuf.norb; iik++)
                                      {
                                          if (pistuf.listpi[iik] == i)
                                             ip = iik;
                                          if (pistuf.listpi[iik] == k)
                                             kp = iik;
                                      }
                                      rho = pimatx.g[kp][ip];
                                      /*  if hydrocarbon or hydrocarbon radical or cation calc pi compression only */
                                      if( ihydc == 0 )
                                      {
                                          disp = sqrt( (atom[i].x - atom[k].x)* (atom[i].x - atom[k].x) + 
                                                                             (atom[i].y - atom[k].y)* (atom[i].y - atom[k].y) + 
                                                                             (atom[i].z - atom[k].z)* (atom[i].z - atom[k].z) );
                                          esig += 90.2;
                                          if( kt == 202 || kt == 203 )
                                          {
                                               ecomp += 517.*(disp - 1.51)*(disp - 1.51);
                                          }else if( kt == 204 )
                                          {
                                               ecomp += 517.*(disp - 1.47)*(disp - 1.47);
                                          }else if( kt == 404 )
                                          {
                                               ecomp += 517.*(disp - 1.42)*(disp - 1.42);
                                          }else if( ((kt == 206 || kt == 306) || kt == 237) || kt == 337 )
                                          {
                                               ecomp += 517.*(disp - 1.39)*(disp - 1.39);
                                          }else if( kt == 307 )
                                          {
                                               ecomp += 517.*(disp - 1.31)*(disp - 1.31);
                                          }else if( kt == 209 || kt == 309 )
                                          {
                                               ecomp += 517.*(disp - 1.43)*(disp - 1.43);
                                          }else if( kt == 410 )
                                          {
                                               esig += 70.;
                                               ecomp += 1000.*(disp - 1.24)*(disp - 1.24);
                                          }else if( kt == 3737 || kt == 937 )
                                          {
                                               esig -= 25.;
                                          }else
                                          {
                                              if ( !gmmx.run )    
                                               fprintf(pcmoutfile," Error  pi-bond %d - %d does not have programmed enthalphy increments\n",i,k);
                                               iunk = 1;
                                          }
                                      }else
                                      {
                                        /*  if i-k is part of a delocalized system.  search the pi-bond list,
                                         *  nnp, to find enthalpy increment.
                                         *
                                         *     search thru pi-bond list */
                                           for( l = 0; l < epiheat.npihf; l++ )
                                           {
                                               nnl = epiheat.nnp[l];
                                               if( nnl == 0 )
                                                   goto L_184;
                                               if( nnl == kt )
                                                   goto L_99957;
                                           }
                                          if ( !gmmx.run )    
                                           fprintf(pcmoutfile," Error  pi-bond %d - %d does not have programmed enthalphy increments\n",i,k);
                                           iunk = 1;
                                           continue;
                                           /*    the bond type is located */
L_99957:
                                           npbond[l] += 1;

//Hay debug
                                           printf("l = %3d\n",l);
                                           printf("npbond[l] = %3d\n",npbond[l]);


                                           /*  use epp(l) to store pi hf contributions calc'd as below */
                                           if( kt == 202 || kt == 404 || kt == 204 ||  kt == 203 || kt == 307 ||
                                               kt == 229 || kt == 237 || kt == 246 || kt == 241 || kt == 410 ||
                                               kt == 741 || kt == 3737 || kt == 4166)
                                           {
                                               if( kt == 202 )
                                               {
                                                   delsig = sqrt( (atom[i].x - atom[k].x)*(atom[i].x - atom[k].x) +
                                                                  (atom[i].y - atom[k].y)*(atom[i].y - atom[k].y) + 
                                                                  (atom[i].z - atom[k].z)*(atom[i].z - atom[k].z) ) - 1.51;
                                                   del2 = delsig*delsig;
                                                   esig += 400.*del2;
                                                   ccrho += rho;
                                               } else
                                               {
                                                   if( rho < alim )
                                                       epp[l] += rho*rho*epiheat.aa[l];
                                                   if( rho >= alim && rho <= blim )
                                                       epp[l] += epiheat.bb[l];
                                                   if( rho > blim )
                                                       epp[l] += -epiheat.ccc[l]*(1 - rho)*(1-rho);
                                               }
                                               for( iip = 0; iip < pistuf.ntpi; iip++ )
                                               {
                                                    if( (torsions.i14[pitorsave.ncsave[iip]][1] == kp &&
                                                         torsions.i14[pitorsave.ncsave[iip]][2] == ip) ||
                                                        (torsions.i14[pitorsave.ncsave[iip]][1] == ip &&
                                                         torsions.i14[pitorsave.ncsave[iip]][2] == kp) )
                                                    {
                                                         ia = torsions.i14[pitorsave.ncsave[iip]][0];
                                                         ib = torsions.i14[pitorsave.ncsave[iip]][1];
                                                         ic = torsions.i14[pitorsave.ncsave[iip]][2];
                                                         id = torsions.i14[pitorsave.ncsave[iip]][3];
                                                         dihed = dihdrl(ia,ib,ic,id);
                                                         ang = .01745329F*fabs( dihed );                           
                                                         epitor += (1.0-cos(2.0*ang))*100.0*.02*rho;
                                                    }
                                               }
                                               continue;
                                           }else if( kt == 230 || kt == 248 || kt == 206 || kt == 630 || 
                                                      kt == 930 || kt == 309 || kt == 306 || kt == 348 || kt == 242 || 
                                                      kt == 337 || kt == 209 || kt == 4142 || kt == 342 || kt == 366 || kt == 4166 )
                                           {
                                                 epp[l] += rho*epiheat.bb[l];
                                                 for( iip = 0; iip < pistuf.ntpi; iip++ )
                                                 {
                                                      if( (torsions.i14[pitorsave.ncsave[iip]][1] == kp &&
                                                           torsions.i14[pitorsave.ncsave[iip]][2] == ip) ||
                                                           (torsions.i14[pitorsave.ncsave[iip]][1] == ip &&
                                                           torsions.i14[pitorsave.ncsave[iip]][2]) )
                                                       {
                                                          ia = torsions.i14[pitorsave.ncsave[iip]][0];
                                                          ib = torsions.i14[pitorsave.ncsave[iip]][1];
                                                          ic = torsions.i14[pitorsave.ncsave[iip]][2];
                                                          id = torsions.i14[pitorsave.ncsave[iip]][3];
                                                          dihed = dihdrl(ia,ib,ic,id);
                                                          ang = .01745329F*fabs( dihed );                          
                                                          epitor += (1.0-cos(2.0*ang))*100.0* .005*rho;
                                                       }
                                                 }
                                                 continue;
                                           }
                                           /*     bond is not in nnp list.  add it but start npbond counter at 1001. */
L_184:
                                           npbond[l] = 1001;
                                           epiheat.nnp[l] = kt;
                                      }
                                   }
                               }
                            }
                        }
                }
        }
        if( ihydc == 0 )
        {
                eppi = 23.061*elowht.epi;
                epp[0] = 0.;
        }
        else
        {
                epp[0] += (nclpi - ccrho)*35.;
               if ( !gmmx.run )
               {    
                fprintf(pcmoutfile," Parameteres for bonds in pi calculations\n");
                fprintf(pcmoutfile," type i - type j      description    contributions    aa    bb   cc\n");
               }
                for( i = 0; i < epiheat.npihf; i++ )
                {
                        if( npbond[i] >= 1 && npbond[i] <= 1000 )
                        {
                                iatm1 = epiheat.nnp[i]/100;
                                iatm2 = epiheat.nnp[i] - iatm1*100;
                               if ( !gmmx.run )    
                                fprintf(pcmoutfile," %d   %d  %s  %8.3f   %8.3f   %8.3f  %8.3f\n",iatm1,iatm2,
                                                epiheat.ccp[i],epiheat.eep[i],epiheat.aa[i],epiheat.bb[i],epiheat.ccc[i]);
                        }
                }
        }
        if ( !gmmx.run )
        {    
         fprintf(pcmoutfile,"The derived quantities above are for all sigma and nondelocalized\n");
         fprintf(pcmoutfile,"pi systems.  Note the following which includes delocalized pi systems.\n");
         }
        if( ihydc != 0 )
        {
                if ( !gmmx.run )
                {    
                 fprintf(pcmoutfile," PI bond enthalpy sums from bond orders\n"); 
                 fprintf(pcmoutfile," #   Bond or Structure         Sums from pi bos\n");
                } 
                for( i = 0; i < epiheat.npihf; i++ )
                {
                        nb = npbond[i];
                        if( nb != 0 )
                        {
                                if( nb > 1001 )
                                {
                                        it = epiheat.nnp[i]/100;
                                        kt = epiheat.nnp[i] - it*100;
                                        nb -= 100;
                                }
                                else
                                {
                                        sls = epp[i];
                                        if( i == 21 )
                                                sls += -noazo*10.;
                                        pibe += sls;
                                       if ( !gmmx.run )    
                                        fprintf(pcmoutfile,"  %d  %s  %11.3f\n",nb,epiheat.ccp[i],sls);
                                }
                        }
                }
               if ( !gmmx.run )   
                fprintf(pcmoutfile,"\n                   pibe =  %9.3f\n",pibe);
        }
        pfc = 2.4;
       if ( !gmmx.run )    
        fprintf(pcmoutfile,"\npartition function contribution (pfc)          = %9.3f\n",pfc); 
         
        if( ihydc == 0 )
        {
                hfo = energies.total + *be + pfc - esig - eppi + ecarb +  ecomp;
                if ( !gmmx.run )
                {    
                 fprintf(pcmoutfile,"strainless  sigma energy (esig)                = %9.3f\n",-esig);
                 fprintf(pcmoutfile,"sigma compression energy (ecomp)               = %9.3f\n",ecomp);
                 fprintf(pcmoutfile,"heavy atom contribution to delta hf (ecarb)    = %9.3f\n",ecarb);
                 fprintf(pcmoutfile,"pople pi energy (eppi)                         = %9.3f\n",-eppi);
                 fprintf(pcmoutfile,"delta hf (energy+be+pfc+esig+ecomp+ecarb+eppi) = %9.3f\n\n",hfo);
                }
        }
        else
        {
                hfo = energies.total + *be + pibe + pfc + esig;
                if ( !gmmx.run )
                {    
                 fprintf(pcmoutfile,"strainless  sigma energy (esig)                = %9.3f\n",esig);
                 fprintf(pcmoutfile,"csp2-csp2 sigma compression energy (ecomp)     = %9.3f\n",ecomp);
                 fprintf(pcmoutfile,"delta h f (energy+be+pibe+pfc+ecomp)           = %9.3f\n",hfo);
                }
        }
       if ( !gmmx.run ) 
       {  
         fprintf(pcmoutfile,"warning: entropy contribution from pi system is unknown\n");
        if( piheat.tepi > 0.002306 )
        {
           fprintf(pcmoutfile,"Caution, the energy difference in the last two cycles was %10.5f kcal\n",piheat.tepi);
        }
        else
        {
                fprintf(pcmoutfile,"The pi vescf calculation converged.  The energy difference between\n");
        fprintf(pcmoutfile,"the last two cycles was %10.5f kcal\n",piheat.tepi);
        }
       }

//Hay debug
        printf("nep[0] = %3d\n",nep[0]);    
        printf("nep[1] = %3d\n",nep[1]);    
        printf("nep[2] = %3d\n",nep[2]);    
        printf("nep[3] = %3d\n",nep[3]);    
        printf("nep[4] = %3d\n",nep[4]);    


        free_ivector( nep ,0,200);
        free_ivector( npbond, 0,200);
        free_vector( epp ,0,200);

        return;
}
