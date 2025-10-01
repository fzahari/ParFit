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
#include "bonds_ff.h"
#include "rings.h"
#include "field.h"
#include "fix.h"
#include "pistuf.h"
#include "pisave.h"
#include "tsbond.h"
#include "files.h"
#include "bondk.h"
#include "ebond.h"

int kbond()
{
    long int pi_mask;
    int i, j,k, iit, kit, it, ierr,nbpi, ii, kk, i1, i2;
    int igeni, igenk, ifound;
    int ia, ib, iend, kend, iy;
    int numi, numj, iatt_electro[6], jatt_electro[6];
    float jatt_const[6],iatt_const[6];
    float bconst, blength;
    float sslope, tslope, tslope2;
    double dbl,dks,bl0;
    char hatext[80];
    char pa[4],pb[4], pt[7], pt_gen[7];

    pi_mask = 1L << PI_MASK;   /*  this mask is for pi atoms - now uses bit 0 */
    nbpi = 0;
    if( natom != 0 )
    {

      for( i = 0; i < bonds_ff.nbnd; i++ )
      {
         ia = bonds_ff.i12[i][0];
         ib = bonds_ff.i12[i][1];
         iit = atom[bonds_ff.i12[i][0]].type;
         kit = atom[bonds_ff.i12[i][1]].type;

         igeni = atom[bonds_ff.i12[i][0]].tclass;
         igenk = atom[bonds_ff.i12[i][1]].tclass;

// metal general class
         if (iit > 300)
           igeni = 300;
         if (kit > 300)
           igenk = 300;

         if( field.type == MMX)
         {
             if (iit == 40 && ( kit != 2 && kit != 3 && kit != 4 && kit != 40) )
                iit = 2;
             if (kit == 40 && ( iit != 2 && iit != 3 && iit != 4 && iit != 40) )
                kit = 2;
         }
         numeral(iit, pa, 3);
         numeral(kit, pb, 3);
         strcpy(pt,pa);
         strcat(pt,pb);
         if( iit <= kit )
         {
            strcpy(pt,pa);
            strcat(pt,pb);
         }else
         {
            strcpy(pt,pb);
            strcat(pt,pa);
         }

         numeral(igeni, pa, 3);
         numeral(igenk, pb, 3);
         if( igeni <= igenk )
         {
            strcpy(pt_gen,pa);
            strcat(pt_gen,pb);
         }else
         {
            strcpy(pt_gen,pb);
            strcat(pt_gen,pa);
         }

       ierr = FALSE;
// first check for transition state constants
       if (field.type == MMX )
       {
           if ( (iit == 37 || iit == 43|| iit == 45 || iit == 49 || iit == 50 || iit == 51 || iit == 52 || iit == 53 ||
                 iit == 54 || iit == 55 || iit == 29) && ( kit == 43 || kit == 45 || kit == 49 || kit == 50
                 || kit == 51 || kit == 52 || kit == 53 || kit == 54 || kit == 55 || kit == 29) )
           {
               // check to see if bond orders are set
               ifound = FALSE;
               for (ii = 0; ii < 11; ii++)
               {
                   if (ts_bondorder.fbnd[ii] != 0.0)
                      ifound = TRUE;
               }
               if (ifound == FALSE)
              {
                  message_alert("TS bond orders are not set","ERROR");
                  return FALSE;
              }
                  
              transbond(&ierr, pt, &bconst, &blength);
              if (ierr == TRUE)
              {
                  bonds_ff.bk[i] = bconst;
                  bonds_ff.bl[i] = blength;
              }
           }
           if (ierr == TRUE)
              goto L_10;   
           if ( (iit == 1 && kit == 45 ) || (kit == 1 && iit == 45 ) )
           { 
              transbond(&ierr, pt, &bconst, &blength);
              if (ierr == TRUE)
              {
                  bonds_ff.bk[i] = bconst;
                  bonds_ff.bl[i] = blength;
              }
           }
       }
       if (ierr == TRUE)
         goto L_10;   
//  check for pi bonds
      if (bondpk.use_pibond == TRUE)
      {
          if ( (atom[ia].flags & pi_mask) && (atom[ib].flags & pi_mask) )
          {
            find_pibond(&ierr, pt, &bconst, &blength, &sslope, &tslope, &tslope2);
             if (ierr == TRUE)
             {
                 bonds_ff.bk[i] = bconst;
                 bonds_ff.bl[i] = blength;
                 ii = -1;
                 kk = -1;
                 ia = bonds_ff.i12[i][0];
                 ib = bonds_ff.i12[i][1];
                 for (j=0; j < pistuf.norb; j++)
                 {
                     if (pistuf.listpi[j] == ia && ia < ib)
                       ii = j;
                     if (pistuf.listpi[j] == ia && ia > ib)
                       kk = j;
                     if (pistuf.listpi[j] == ib && ia < ib)
                       kk = j;
                     if (pistuf.listpi[j] == ib && ia > ib)
                       ii = j;
                 }
                 pibondsave.ibsave[nbpi] = ii;
                 pibondsave.kbsave[nbpi] = kk;
                 pibondsave.nbsave[nbpi] = i;
                 pibondsave.bksave[nbpi] = bconst;
                 pibondsave.blsave[nbpi] = blength;
                 pibondsave.pksave[nbpi] = sslope;
                 pibondsave.plsave[nbpi] = tslope;
                 pibondsave.plsave2[nbpi] = tslope2;
                 nbpi++;
             } else
             {
	       set_missing_constants(TRUE);
                 sprintf(hatext,"Error pi bond %s not found. Please check parameter file",pt);
                 message_alert(hatext,"Error");
             }
          }
      }
      if (ierr == TRUE)
         goto L_10;
//  check 3 membered ring
      if (rings.nring3 > 0 && bondk1.nbnd3 > 0 && ( is_ring31(ia) && is_ring31(ib) ) )
      {
        for (j=0; j < bondk1.nbnd3; j++)
        {
          if ( (strcmp(bondk1.kb3[j], pt) == 0)  || (strcmp(bondk1.kb3[j],pt_gen) == 0) )
          {
              bonds_ff.bk[i] = bondk1.s3[j];
              bonds_ff.bl[i] = bondk1.t3[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              break;
          }
         }
         if (ierr != TRUE)
         {
             set_missing_constants(TRUE);
             sprintf(hatext,"Bond constants missing for cyclopropane bond %d - %d of type %d-%d\n",bonds_ff.i12[i][0],bonds_ff.i12[i][1],
             atom[bonds_ff.i12[i][0]].type, atom[bonds_ff.i12[i][1]].type);
             fprintf(pcmoutfile,"%s\n",hatext);
            fprintf(errfile,"%s\n",hatext);
         }
      }
      if (ierr == TRUE)
       goto L_10;
//  check 4 membered ring
      if (rings.nring4 > 0 && bondk1.nbnd4 > 0 && ( is_ring41(ia) && is_ring41(ib) ) )
      {
         for (j=0; j < bondk1.nbnd4; j++)
         {
          if ( (strcmp(bondk1.kb4[j],pt) == 0)  || (strcmp(bondk1.kb4[j],pt_gen) == 0) )
          {
              bonds_ff.bk[i] = bondk1.s4[j];
              bonds_ff.bl[i] = bondk1.t4[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              break;
          }
         }
         if (ierr != TRUE)
         {
            set_missing_constants(TRUE);
            sprintf(hatext,"Bond constants missing for cyclobutane bond %d - %d of type %d-%d\n",bonds_ff.i12[i][0],bonds_ff.i12[i][1],
            atom[bonds_ff.i12[i][0]].type, atom[bonds_ff.i12[i][1]].type);
            fprintf(pcmoutfile,"%s\n",hatext);
            fprintf(errfile,"%s\n",hatext);
         }
      }
      if (ierr == TRUE)
       goto L_10;
//  check 5 membered ring
      if (rings.nring5 > 0 && bondk1.nbnd5 > 0 && (is_ring51(ia) && is_ring51(ib)))
      {
         for (j=0; j < bondk1.nbnd5; j++)
         {
          if (strcmp(bondk1.kb5[j],pt) == 0 )
          {
              bonds_ff.bk[i] = bondk1.s5[j];
              bonds_ff.bl[i] = bondk1.t5[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              break;
          }
         }
         // don't fail on missing param - check for regular param first
     }
      if (ierr == TRUE)
       goto L_10;
//  delocalized bonds in MMFF94
      if (field.type == MMFF94 && bondk1.ndeloc > 0)
      {
          if ( is_delocalbond(ia, ib) )
          {
             for (k = 0; k < MAXIAT; k++)
             {
                 if (atom[ia].iat[k] == ib)
                 {
                    if (atom[ia].bo[k] == 1) // single bond
                    {
                      for (j=0; j < bondk1.ndeloc; j++)
                      {
                          if (strcmp(bondk1.kbdel[j],pt) == 0)
                          {
                              bonds_ff.bk[i] = bondk1.sdel[j];
                              bonds_ff.bl[i] = bondk1.tdel[j];
                              bonds_ff.index[i] = 1;
                              ierr = TRUE;
                              break;
                          }
                      }
                    }
                 }
              }
          }
          if (ierr == TRUE)
            goto L_10;
      }
//  regular parameters
      for (j=0; j < bondk1.nbnd; j++)
      {
          if ( (strcmp(bondk1.kb[j], pt) == 0) )
          {
              bonds_ff.bk[i] = bondk1.s[j];
              bonds_ff.bl[i] = bondk1.t[j];
              bonds_ff.index[i] = 0;
              ierr = TRUE;
              if ( ((strcmp(pt,"  2  2") == 0) || (strcmp(pt,"  2  3") == 0))&& field.type == MMX)
              {
                  for (k=0; k < MAXIAT; k++)
                  {
                      if (atom[ia].iat[k] == ib)
                      {
                          if (atom[ia].bo[k] == 1)  // butadiene
                          {
                              bonds_ff.bl[i] = 1.48;
                              break;
                          }
                      }
                  }
              }
              break;
          }
      }
      if (ierr == TRUE)
         goto L_10;      
//  check for metal bonds - have not found specific metal bonds - now looking for general mmx type bonds
      if ( iit >= 300 || kit >= 300  )
      {
          find_metalbond(&ierr,bonds_ff.i12[i][0],bonds_ff.i12[i][1], pt_gen, &bconst, &blength);
          if (ierr == TRUE)
          {
              bonds_ff.bk[i] = bconst;
              bonds_ff.bl[i] = blength;
          }
      }
      if (ierr == TRUE)
      {
          fprintf(pcmoutfile,"Using generalized metal bonds for MMX force field - %d %d of types %d %d\n",ia,ib,iit,kit);
          goto L_10;
      }
      if (ierr != TRUE)
      {
          set_missing_constants(TRUE);
          sprintf(hatext,"Bond constants missing for bond %d - %d of type %d - %d\n",bonds_ff.i12[i][0],bonds_ff.i12[i][1],
           atom[bonds_ff.i12[i][0]].type, atom[bonds_ff.i12[i][1]].type);
          fprintf(pcmoutfile,"%s\n",hatext);
          fprintf(errfile,"%s\n",hatext);
      }
L_10:
      continue;
          
    }
  }
  if (field.type == MM3) // need check for electronegativity correction
  {
      for( i = 0; i < bonds_ff.nbnd; i++ )
      {
         ia = bonds_ff.i12[i][0];
         ib = bonds_ff.i12[i][1];
         iit = atom[bonds_ff.i12[i][0]].tclass;
         kit = atom[bonds_ff.i12[i][1]].tclass;
         if (iit > kit)
         {
             iend = ib;
             kend = ia;
             it = 1000*kit+iit;
         }else
         {
             iend = ia;
             kend = ib;
             it = 1000*iit+kit;
         }
         numi = 0;
         for (j=0; j < 6; j++)
         {
           iatt_electro[j] = 0;
           iatt_const[j] =  0.0;
           jatt_electro[j] = 0;
           jatt_const[j] = 0.0;
         }
         for (j=0; j < electroneg.nelecti; j++)
         {
             if (it == electroneg.ibond[j])
             {
                 for (k=0; k < MAXIAT; k++)
                 {
                     if (atom[iend].iat[k] != kend && atom[atom[iend].iat[k]].tclass == electroneg.iattach[j])
                     {
                         iatt_electro[numi]++;
                         iatt_const[numi] = electroneg.icorr[j];
                         numi++;
                         if (numi > 5) message_alert("error in kbond numi > 5","Error");
                     }
                 }
                 if (iit == kit)
                 {
                     for (k=0; k < MAXIAT; k++)
                     {
                         if (atom[kend].iat[k] != iend && atom[atom[kend].iat[k]].tclass == electroneg.iattach[j])
                         {
                             iatt_electro[numi]++;
                             iatt_const[numi] = electroneg.icorr[j];
                             numi++;
                             if (numi > 5) message_alert("error in kbond numi > 5","Error");
                         }
                     }
                 }
             }
         }                 
         // loop over ib   "      "     "             "
         numj = 0;
         jatt_electro[0] = jatt_electro[1] = jatt_electro[2] = 0;
         jatt_const[0] = jatt_const[1] = jatt_const[2] = 0.0;
         for (j=0; j < electroneg.nelectj; j++)
         {
             if (it == electroneg.jbond[j])
             {
                 for (k=0; k < MAXIAT; k++)
                 {
                     if (atom[kend].iat[k] != iend && atom[atom[kend].iat[k]].tclass == electroneg.jattach[j])
                     {
                         jatt_electro[numj]++;
                         jatt_const[numj] = electroneg.jcorr[j];
                         numj++;
                         if (numj > 5) message_alert("error in kbond numj > 5","Error");
                     }
                 }
             }
         }
         // correct bond length
         dbl = 0.0;
         iy = -1;
         for (j=0; j < numi; j++)
         {
             iy++;
             dbl += iatt_const[j]*powi(0.62,iy);
         }
         for (j=0; j < numj; j++)
         {
             iy++;
             dbl += jatt_const[j]*powi(0.62,iy);
         }
         bonds_ff.bl[i] += dbl;
          //adjust force const
         if (dbl != 0.0)
         {
             if (it == 1005 || it == 5022 || it == 5056)
             {
                 bl0 = bonds_ff.bl[i] - dbl;
                 dks = deltaks(dbl,bl0);
                 bonds_ff.bk[i] += dks;
             }
         }
  
      }
  }
  if (get_missing_constants())
  {
      fprintf(pcmoutfile,"Error assigning constants in: %s\n",Struct_Title);
      fprintf(errfile,"%s\n",hatext);
      return FALSE;
  }

  if (fixdis.FixBond)  // fixed bonds replace bond constant with fixed constant
  {
      for (i=0; i < fixdis.nfxstr; i++)
      {
          ia = fixdis.ifxstr[i][0];
          ib = fixdis.ifxstr[i][1];
          if (is_bond(ia,ib))
          {
              for (j=0; j < bonds_ff.nbnd; j++)
              {
                  i1 = bonds_ff.i12[j][0];
                  i2 = bonds_ff.i12[j][1];
                  if ( ((ia == i1) && (ib == i2)) || ((ia == i2) && (ib == i1)) )
                  {
                      bonds_ff.bl[j] = fixdis.fxstrd[i];
                      bonds_ff.bk[j] = fixdis.fxstrc[i];
                  }
              }
          }
      }
  }
          
  return TRUE;     
}
// ==============================================
double deltaks(double dbl,double bl0)
{
    double a,b,pi,clite,amu,conv,delta;
    a = 1.3982;
    b = -0.0001023;
    pi = acos(-1.0);
    clite = 2.9979e10;
    amu = 1.660531e-24;
    conv = 4.0*pi*pi*clite*clite*1.008*amu*1.0e-5;
    delta = conv*dbl*2.0*(bl0-a)/(b*b);
    return(delta);
}
/* ------------------------------------------  */
void find_metalbond(int *ierr,int at1, int at2, char *pt, float *bconst, float *blength)
{
    int j, k, itm, ito, nusi, ntm, mi,imi;
    long int mask5, mask6, mask7, mask8;
    float a13;

    *blength = 0.0;

      for (j=0; j < bondk1.nbnd; j++)
      {
          if (strcmp(bondk1.kb[j], pt) == 0)
          {
              *bconst = bondk1.s[j];
              *blength = bondk1.t[j];
              *ierr = TRUE;
              break;
          }
      }
    if (strcmp(pt,"300300") == 0) 
    {
        if (atom[at1].type == atom[at2].type)
        {
            *blength += 2*atom[at1].radius - 2.52;
        } else
        {
            *blength += atom[at1].radius + atom[at2].radius -2.52;
        }
        for (k=0; k < MAXIAT; k++)
        {
            if (atom[at1].iat[k] == at2)
            {
                if (atom[at1].bo[k] > 1 && atom[at1].bo[k] != 9)
                   *blength += -atom[at1].bo[k]*0.19;
            }
        }
    }
    if ( (atom[at1].type >= 300 && atom[at2].type < 300 ) || (atom[at2].type >= 300 && atom[at1].type < 300)  )
    {
        if (atom[at1].type >= 300)
        {
            itm = at1;
            ito = at2;
        }else
        {
            itm = at2;
            ito = at1;
        }
        if (fabs(*blength) < 0.1)
           *blength += atom[itm].vdw_radius - 1.26;
        /*  flags 5 and 6 are used to mark electron count
         *    1,0 for saturated:  0,1 gt 18 electron and 1,1 for unsaturated
         *
         *  flags 7 and 8 are used for high spin, low spin and sq planar
         *    1,0 for low spin:  0,1 for sq planar : 1,1, for high spin
         *
         * find out if unsaturated and then the spin state
         *
         *  nusi - 0 : coord sat'd   (also high spin)
         *  nusi - 1 : >18 electrons
         *  nusi - 2 : low spin
         *  nusi - 3 : square planar */
          mask5 = 1L << SATMET_MASK;
          mask6 = 1L << GT18e_MASK;
          mask7 = 1L << LOWSPIN_MASK;
          mask8 = 1L << SQPLAN_MASK;
          nusi = 0;
          if( ( atom[itm].flags & mask6 ) &&  !( atom[itm].flags & mask5 ) )
          {
              nusi = 1;
          }else if( ( atom[itm].flags & mask5 ) && ( atom[itm].flags & mask6 ) )
          {
                if( ( atom[itm].flags & mask7 ) && !( atom[itm].flags & mask8 ) )
                {
                     nusi = 2;
                }else if(( atom[itm].flags & mask8 ) && !( atom[itm].flags & mask7 ) )
                {
                     nusi = 3;
                }
          }

          if( nusi == 2 || nusi == 3 )
          {
             /* unsaturated carbon */
             if( atom[ito].type == 2 )
             {
                *blength -= 0.15;
                            /* unsaturated oxygen or sulfur */
             }else if( atom[ito].type == 6 || atom[ito].type == 15 )
             {
               *blength -= 0.1;
                           /* unsaturated halogens */
             }else if( atom[ito].type >= 11 && atom[ito].type <= 15 )
             {
               *blength -= 0.2;
             }
             *bconst *= 1.5;
          }
          if( nusi == 1 )
              *blength += .15;

          /* sq planar complexes
           * found a metal-halogen search other bonds to metal for trans effects */
          if( nusi == 3 && (atom[ito].type== 12 || atom[ito].type == 13) )
          {
              for( mi = 0; mi < MAXIAT; mi++ )
              {
                  imi = atom[itm].iat[mi];
                  if( atom[itm].iat[mi] != 0 )
                  {
                      ntm = atom[imi].type;
                      if( ntm == 1 || ntm == 2 || ntm == 4 || ntm == 5 || ntm == 25 || ntm == 30 )
                      {
                          angle( ito, itm, imi,&a13 );
                          if( a13 > 150 || a13 <  -150. )
                          {
                              *blength += .20;
                              fprintf(pcmoutfile,"Trans influence (+0.20 a) on metal-cl or br bond for %d - %d\n",itm,ito); 
                          }
                      }
                  }
              }
          }

          /* more sq planar tests */
          if( nusi == 3 && (atom[ito].type == 1 || atom[ito].type == 2 || atom[ito].type ==  5 ||
                            atom[ito].type == 25 || atom[ito].type == 30) )
          {
              for( mi = 0; mi < MAXIAT; mi++ )
              {
                 imi = atom[itm].iat[mi];
                 if( imi != 0 )
                 {
                     ntm = atom[imi].type;
                     if( ntm == 12 || ntm == 13 )
                     {
                        angle( ito, itm, imi,&a13);
                        if( a13 > 150 || a13 < -150. )
                        {
                           *blength -= .1;
                           fprintf(pcmoutfile,"trans influence (-0.1 a) on metal-ligand bond: %d - %d\n",itm,ito);
                        }
                     }
                 }
              }
          }
          /* metal to nitrile */
          if( atom[ito].type == 10 )
          {
             for( mi = 0; mi < MAXIAT; mi++ )
             {
                 if( atom[ito].iat[mi] != 0 )
                 {
                    ntm = atom[atom[ito].iat[mi]].type;
                    if( ntm == 4 || ntm == 10 )
                    {
                        *blength = 1.91;
                        *bconst = 2.0;
                    }
                 }
             }
          }
          /* metal to amine */
          if( atom[ito].type == 7 )
          {
            for( mi = 0; mi < MAXIAT; mi++ )
            {
                if( atom[ito].iat[mi] != 0 )
                {
                   if( atom[atom[ito].iat[mi]].type == 3 )
                   {
                      *blength = 1.70;
                      *bconst = 2.0;
                   }
                }
            }
          }
         /* metal to ammonium */
         if( atom[ito].type == 41 )
         {
           for( mi = 0; mi < MAXIAT; mi++ )
           {
              if( atom[ito].iat[mi] != 0 )
              {
                 ntm = atom[atom[ito].iat[mi]].type;
                 if( ntm == 7 )
                     *blength = 1.67;
                 if( ntm == 37 || ntm == 41 )
                     *blength = 1.73;
              }
           }
         }
    }
        
}
/* ------------------------------------------  */
void find_pibond (int *ierr, char *pt, float *bconst, float *blength, float *sslope, float *tslope,float *tslope2)
{
    int m;
    *ierr = FALSE;
    for (m = 0; m < bondpk.npibond; m++)
    {
       if (strcmp(bondpk.kb[m],pt) == 0)
       {
          if (bondpk.bk[m] != 0.0 && bondpk.bl[m] != 0.0)
          {
             *bconst = bondpk.bk[m];
             *blength = bondpk.bl[m];
             *sslope = bondpk.sslop[m];
             *tslope = bondpk.tslop[m];
             *tslope2 = bondpk.tslop2[m];
             *ierr = TRUE;
             break;
           }
       }
    }
}
/* ------------------------------------------  */
int find_bond(int ia, int ib)
{
    int i;
    int iz;
    
    iz = -1;
    for (i=0; i < bonds_ff.nbnd; i++)
    {
        if (ia == bonds_ff.i12[i][0] && ib == bonds_ff.i12[i][1])
           return(i);
        if (ib == bonds_ff.i12[i][0] && ia == bonds_ff.i12[i][1])
           return(i);
    }
    return(iz);
}

/* ------------------------------------------  */
void  angle(int k,int i,int imi,float *a13)
{
        float disik, disim, dx1,dy1,dz1, dx2,dy2,dz2;

        /*  returns angle between k-i-imi */

     dx1 = atom[i].x - atom[k].x;
     dy1 = atom[i].y - atom[k].y;
     dz1 = atom[i].z - atom[k].z;
     dx2 = atom[i].x - atom[imi].x;
     dy2 = atom[i].y - atom[imi].y;
     dz2 = atom[i].z - atom[imi].z;

     disik = sqrt( dx1*dx1 + dy1*dy1 + dz1*dz1);
     disim = sqrt( dx2*dx2 + dy2*dy2 + dz2*dz2);

        *a13 = angs( dx1, dy1, dz1, dx2, dy2, dz2, disik, disim );

        *a13 *= 57.2958;

        return;
} 
/* ------------------------------ */
double angs(float x1,float y1,float z1,float x2,float y2,float z2,float d1,float d2)
{
        double angs_v, cosa;

        /*     *** calculates angles for subroutine ebend *** */
        cosa = (x1*x2 + y1*y2 + z1*z2)/(d1*d2);
        if( fabs( cosa ) > 1.0 )
                cosa /= fabs( cosa );
        angs_v = arcos( cosa );
        return( angs_v );
} 
/* ----------------------- */
double arcos(float x)
{
        double arcos_v, y;

        y = x;
        if( y > 1.0 )
                y = 1.0;
        if( y < -1.0 )
                y = -1.0;
        arcos_v = acos( y );
        return( arcos_v );
} 
void transbond(int *ierr, char *pt,float *bconst,float *blength)
{
        int iadj, iatjk, iii, it1st, it2nd, itjk, jjj, kkk, nbond;
        float b4953, bmomk, s4950, t4950;

	s4950 = t4950 = 0.0;
// C* C*
        if( strcmp(pt," 49 49") == 0)
        {
                *bconst = 4.4*ts_bondorder.fbnd[0];
                *blength = 1.53 - 0.60*log( ts_bondorder.fbnd[0] );
        }
// C# C#
        else if(  strcmp(pt," 50 50") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[1];
                *blength = 1.53 - 0.60*log( ts_bondorder.fbnd[1] );
        }
// C* C#
        else if(  strcmp(pt," 49 50") == 0 )
        {
                if( ts_bondorder.fbnd[0] > 0.0F && ts_bondorder.fbnd[1] > 0.0F )
                {
                        s4950 = 9.6 - ts_bondorder.fbnd[0]*ts_bondorder.fbnd[1]*5.2;
                        t4950 = 1.337 + (ts_bondorder.fbnd[0] + ts_bondorder.fbnd[1])*.1;
                }
                if( ts_bondorder.fbnd[0] > 0.0F && ts_bondorder.fbnd[4] > 0.0F )
                {
                        s4950 = 9.6 - ts_bondorder.fbnd[0]*ts_bondorder.fbnd[4]*5.2;
                        t4950 = 1.337 + (ts_bondorder.fbnd[0] + ts_bondorder.fbnd[4])*.1;
                }
                if( ts_bondorder.fbnd[0] > 0.0F && ts_bondorder.fbnd[6] > 0.0F )
                {
                        s4950 = 9.6 - ts_bondorder.fbnd[0]*ts_bondorder.fbnd[6]*5.2;
                        t4950 = 1.337 + (ts_bondorder.fbnd[0] + ts_bondorder.fbnd[6])*.1;
                }
                if( ts_bondorder.fbnd[3] > 0.0F && ts_bondorder.fbnd[4] > 0.0F )
                {
                        s4950 = 9.6 - ts_bondorder.fbnd[3]*ts_bondorder.fbnd[4]*5.2;
                        t4950 = 1.337 + (ts_bondorder.fbnd[3] + ts_bondorder.fbnd[4])*.1;
                }
                if( ts_bondorder.fbnd[1] > 0.0F && ts_bondorder.fbnd[5] > 0.0F )
                {
                        s4950 = 9.6 - ts_bondorder.fbnd[1]*ts_bondorder.fbnd[5]*5.2;
                        t4950 = 1.337 + (ts_bondorder.fbnd[1] + ts_bondorder.fbnd[5])*.1;
                }
                if( ts_bondorder.fbnd[5] > 0.0F && ts_bondorder.fbnd[6] > 0.0F )
                {
                        s4950 = 9.6 - ts_bondorder.fbnd[5]*ts_bondorder.fbnd[6]*5.2;
                        t4950 = 1.337 + (ts_bondorder.fbnd[5] + ts_bondorder.fbnd[6])*.1;
                }
                /* newc new on 8/9/87 */
                if( ts_bondorder.fbnd[3] > 0.0F && ts_bondorder.fbnd[6] > 0.0F )
                {
                        s4950 = 9.6 - ts_bondorder.fbnd[3]*ts_bondorder.fbnd[6]*5.2;
                        t4950 = 1.337 + (ts_bondorder.fbnd[3] + ts_bondorder.fbnd[6])*.1;
                }
                if( ts_bondorder.fbnd[1] > 0.0F && ts_bondorder.fbnd[3] > 0.0F )
                {
                        s4950 = 9.6 - ts_bondorder.fbnd[1]*ts_bondorder.fbnd[3]*5.2;
                        t4950 = 1.337 + (ts_bondorder.fbnd[1] + ts_bondorder.fbnd[3])*.1;
                }
                *bconst = s4950;
                *blength = t4950;
// C(sp3) H*
        } else if(  strcmp(pt,"  1 45") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[2];
                *blength = 1.10 - 0.60*log( ts_bondorder.fbnd[2] );
// B# H*
        }       else if(  strcmp(pt," 43 45") == 0 )
        {
                *bconst = 3.5*(1. - ts_bondorder.fbnd[3]);
                *blength = 1.05 - 0.60*log( 1. - ts_bondorder.fbnd[3] );
// H* C*
        }       else if(  strcmp(pt," 45 49") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[3];
                *blength = 1.10 - 0.60*log( ts_bondorder.fbnd[3] );
// B# C#
        }       else if(  strcmp(pt," 43 50") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[4];
                *blength = 1.56 - 0.6*log( ts_bondorder.fbnd[4] );
// C# C$ 
        } else if(  strcmp(pt," 49 51") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[5];
                *blength = 1.50 - 0.60*log( ts_bondorder.fbnd[5] );
// C# O#
        } else if(  strcmp(pt," 50 53") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[6];
                *blength = 1.40 - 0.60*log( ts_bondorder.fbnd[6] );
// C* C% (SN2 pentavalent C)
        } else if(  strcmp(pt," 49 52") == 0)
        {
                *bconst = 4.4*ts_bondorder.fbnd[7];
                *blength = 1.53 - 0.60*log( ts_bondorder.fbnd[7] );
// C* O#
        } else if(  strcmp(pt," 49 53") == 0 )
        {
                *bconst = 10.-(5.*ts_bondorder.fbnd[0]*ts_bondorder.fbnd[11]);
                *blength = 1.21 +0.13*(ts_bondorder.fbnd[0]+ts_bondorder.fbnd[11]);
                b4953 = 0.;
                bmomk = 0.;
// C# C%
        } else if(  strcmp(pt," 50 52") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[7];
                *blength = 1.53 - 0.60*log( ts_bondorder.fbnd[7] );
// C% I%
        } else if(  strcmp(pt," 52 54") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[8];
                *blength = 2.2 - 0.60*log( ts_bondorder.fbnd[8] );
// C% O#  
        } else if(  strcmp(pt," 52 53") == 0 )
        {
                *bconst = 5.0*ts_bondorder.fbnd[9];
                *blength = 1.407 - .6*log( ts_bondorder.fbnd[9] );
// C# N#
        } else if(  strcmp(pt," 50 55") == 0 )
        {
                bmomk = 0.;
                *bconst = 5.1*ts_bondorder.fbnd[10];
                *blength = 1.44 - 0.60*log( ts_bondorder.fbnd[10] );
// B# O#  new relative to pcm386
        } else if(  strcmp(pt," 43 53") == 0) // B# O#
        {
                bmomk = 0.;
                *bconst = 6.2*ts_bondorder.fbnd[11];
                *blength = 1.53 - 0.60*log( ts_bondorder.fbnd[11] );
        }  
        /*   Use Pauling's Rule for distance as a function bond order. BO must be
         *   greater than zero-Pauling constant is 2X that for BO 1 to 3. */
        if( ts_bondorder.fbnd[0] > 0.000 )
        {
                nbond = 0;
                iadj = 0;
                for( iii = 0; iii < 15; iii++ )
                {
                        if (ts_bondorder.fbnd[iii] > 0.0F)
                                nbond += 1;
                        if( ts_bondorder.fbnd[iii] > 0.0F )
                        {
                                for( jjj = 1; jjj <= natom; jjj++ )
                                {
                                        if( atom[jjj].type == 49 )
                                        {
                                                for( kkk = 0; kkk < MAXIAT; kkk++ )
                                                {
                                                        iatjk = atom[jjj].iat[kkk];
                                                        if( iatjk != 0)
                                                        {
                                                                itjk = atom[iatjk].type;
                                                                if( itjk == 29 || itjk == 30 || 
                                                                 itjk == 50 || itjk == 53 || itjk == 42)
                                                                        goto L_54323;
                                                        }
                                                }
                                        }
                                }
                        goto L_54322;
L_54323:
                        iadj = itjk;
                        }
L_54322:
                        ;
                }
                if( nbond < 2 )
                {
// C* C#              
                        if(  strcmp(pt," 49 50") == 0 && ts_bondorder.fbnd[0] > 0.0F )
                        {
                                if( iadj == 50 )
                                {
                                        *bconst = 9.6 - ts_bondorder.fbnd[0]*5.2;
                                        *blength = 1.337 + ts_bondorder.fbnd[0]*.19;
                                }
                        }
// C* O#              
                        if(  strcmp(pt," 49 53") == 0 && ts_bondorder.fbnd[0] > 0.0F )
                        {
                                if( iadj == 53 )
                                {
                                        *bconst = 10. - ts_bondorder.fbnd[0]*5.;
                                        *blength = 1.21 + ts_bondorder.fbnd[0]*.20;
                                }
                        }
                }
        }
        sscanf(pt,"%d %d",&it1st,&it2nd);
        if( it2nd == 49 )
        {
// C. C*
                if( it1st == 29 )
                {
                        *bconst = 9.6 - ts_bondorder.fbnd[0]*5.2;
                        *blength = 1.337 + ts_bondorder.fbnd[0]*.17;
// C+ C*
                } else if( it1st == 30 )
                {
                        *bconst = 9.6 - ts_bondorder.fbnd[0]*5.2;
                        *blength = 1.337 + ts_bondorder.fbnd[0]*.17;
// C- C*
                } else if( it1st == 48 )
                {
                        *bconst = 9.6 - ts_bondorder.fbnd[0]*5.2;
                        *blength = 1.337 + ts_bondorder.fbnd[0]*.17;
// C(sp2) C*
                } else if( it1st == 2 )
                {
                        *bconst = 5.2;
                        *blength = 1.46 + .08*ts_bondorder.fbnd[0];
// C(CO) C*
                } else if( it1st == 3 )
                {
                        *bconst = 5.2;
                        *blength = 1.46 + .07*ts_bondorder.fbnd[0];
// C(sp) C*
                } else if( it1st == 4 )
                {
                        *bconst = 5.2;
                        *blength = 1.40 + .06*ts_bondorder.fbnd[0];
// N C*
                } else if( it1st == 8 )
                {
                        *bconst = 9.0 - 4.*ts_bondorder.fbnd[0];
                        *blength = 1.36 + .09*ts_bondorder.fbnd[0];
                }
        } else if( it2nd == 50 )
        {
                if( it1st == 2 )
                {
                        *bconst = 5.2;
                        *blength = 1.46 + .08*ts_bondorder.fbnd[1];
                } else if( it1st == 3 )
                {
                        *bconst = 5.2;
                        *blength = 1.46 + .07*ts_bondorder.fbnd[1];
                } else if( it1st == 4 )
                {
                        *bconst = 5.2;
                        *blength = 1.40 + .06*ts_bondorder.fbnd[1];
                } else if( it1st == 8 )
                {
                        *bconst = 9.0 - 4.*ts_bondorder.fbnd[1];
                        *bconst = 1.36 + .09*ts_bondorder.fbnd[1];
                } else if( it1st == 29 )
                {
                        *bconst = 9.6 - ts_bondorder.fbnd[0]*5.2;
                        *blength = 1.337 + ts_bondorder.fbnd[0]*.17;
                } else if( it1st == 30 )
                {
                        *bconst = 9.6 - ts_bondorder.fbnd[0]*5.2;
                        *blength = 1.337 + ts_bondorder.fbnd[0]*.17;
                } else if( it1st == 48 )
                {
                        *bconst = 9.6 - ts_bondorder.fbnd[0]*5.2;
                        *blength = 1.337 + ts_bondorder.fbnd[0]*.17;
                }
        }
        /* assume all double bonds are delocalized with type 29 (radical) carbons */
// C. C.
        if(  strcmp(pt," 29 29") == 0 )
        {
                *bconst = 7.0;
                *blength = 1.40;
// N(sp2) C$
        } else if(  strcmp(pt," 37 51") == 0 )
        {
                *bconst = 18.5 - 8.2*ts_bondorder.fbnd[5];
                *blength = 1.15 + .15*ts_bondorder.fbnd[5];
// C. O#
        } else if(  strcmp(pt," 29 53") == 0 )
        {
                *bconst = 10.8 - 5.44*ts_bondorder.fbnd[6];
                *blength = 1.208 + .20*ts_bondorder.fbnd[6];
                bmomk = b4953;
// N(sp3) O#
        } else if(  strcmp(pt,"  8 53") == 0 )
        {
                *bconst = 7.5;
                *blength = 1.27 + .20*ts_bondorder.fbnd[6];
// N(sp2) O#
        } else if(  strcmp(pt," 37 53") == 0 )
        {
                *bconst = 7.8 - 2.8*ts_bondorder.fbnd[5];
                *blength = 1.20 + .24*ts_bondorder.fbnd[5];
// C. N#
        } else if(  strcmp(pt," 29 55") == 0 )
        {
                *bconst = 10.3 - 5.2*ts_bondorder.fbnd[10];
                *blength = 1.27 + .17*ts_bondorder.fbnd[10];
                bmomk = 0.;
// C* N#
        } else if(  strcmp(pt," 49 55") == 0 )
        {
                *bconst = 10. - 0.5*ts_bondorder.fbnd[0]*ts_bondorder.fbnd[10];
                *blength = 1.27 + .1*ts_bondorder.fbnd[0]*ts_bondorder.fbnd[10];
                bmomk = 0.;
// N(sp2) N#
        } else if(  strcmp(pt," 37 55") == 0 )
        {
                *bconst = 7.8;
                *blength = 1.40;
// C(cyclopropane) H*
        } else if(  strcmp(pt," 22 45") == 0 )
        {
                *bconst = 4.4*ts_bondorder.fbnd[2];
                *blength = 1.10 - 0.60*log( ts_bondorder.fbnd[2] );
        }
        if( fabs(*blength) < 0.1 )
        {
                *ierr = 0;
        } else
        {
                bmomk = 0.;
                *ierr = 1;
        }
        return;
} 
/* ------------------------------------------  */
int is_delocalbond(int ia, int ib)
{
    // test for delocalized bond
    // if ia & ib are aromatic and part of 6 membered ring
    // if not aromatic is bond order 1

    long int aromatic_mask;
    int i, j;
    int jdbl, kdbl;

    aromatic_mask = (1L << AROMATIC_MASK);

    if ( (atom[ia].flags & aromatic_mask) && (atom[ib].flags & aromatic_mask)
         && (is_ring62(ia,ib) != FALSE) )
         return(FALSE);

    if ( (atom[ia].flags & aromatic_mask) && (atom[ib].flags & aromatic_mask)
         && (is_ring52(ia,ib) != FALSE) )
         return(FALSE);

    if ( (atom[ia].flags & aromatic_mask) && (atom[ib].flags & aromatic_mask) )
         return(TRUE);
         
// use mmff atom types since only called for mmff
   if (atom[ia].type == 37 && atom[ib].type == 37 && (is_ring62(ia,ib) != FALSE))
      return FALSE;
   if (atom[ia].type == 37 && atom[ib].type == 37 && (is_ring52(ia,ib) != FALSE))
      return FALSE;
//   if (atom[ia].type == 2 && atom[ib].type == 2 && (is_ring62(ia,ib) != FALSE))
//      return FALSE;
//   if (atom[ia].type == 2 && atom[ib].type == 2 && (is_ring52(ia,ib) != FALSE))
//      return FALSE;

   for (i=0; i < MAXIAT; i++)
   {
       if (atom[ia].iat[i] == ib)
       {
           if (atom[ia].bo[i] == 1)
           {
             jdbl = FALSE;
             kdbl = FALSE;
             for (j=0; j < MAXIAT; j++)
             {
                 if (atom[ia].bo[j] >= 2 && atom[ia].bo[j] != 9)
                   jdbl = TRUE;
                 if (atom[ib].bo[j] >= 2 && atom[ib].bo[j] != 9)
                   kdbl = TRUE;
             }
             if (jdbl == TRUE && kdbl == TRUE)
                return (TRUE);
             else if (jdbl == TRUE && (atom[ib].flags & aromatic_mask) )
                return (TRUE);
             else if (kdbl == TRUE && (atom[ia].flags & aromatic_mask) )
                return (TRUE);            
             else
                return(FALSE);
           }
       }
   }
   return(FALSE);
}
#define IS_ODD(j)       ((j) & 1 )

double sign(double a,double b )         /* floating polong int sign transfer */
{
        return( b < 0.0 ? -fabs(a) : fabs(a) );
}

double powi(double x, int n )     /* returns:  x^n */
{
   double p;       /* holds partial product */

    if( x == 0 )
        return( 0. );

    if( n < 0 ){    /* test for negative exponent */
        n = -n;
        x = 1/x;
        }

    p = IS_ODD(n) ? x : 1;  /* test & set zero power */

    while( n >>= 1 ){       /* now do the other powers */
        x *= x;                 /* sq previous power of x */
        if( IS_ODD(n) ) /* if low order bit set */
                p *= x;         /*      then, multiply partial product by latest power of x */
        }

    return( p );
}
/* ================================ */
float get_bond(int ia, int ib)
{
    int i, ita,itb;
    char pa[4],pb[4], pt[7];

    ita = atom[ia].type;
    itb = atom[ib].type;
    numeral(ita,pa,3);
    numeral(itb,pb,3);
    if (ita <= itb)
    {
        strcpy(pt,pa);
        strcat(pt,pb);
    } else
    {
        strcpy(pt,pb);
        strcat(pt,pa);
    }
    for (i=0; i < bondk1.nbnd; i++)
    {
        if ( (strcmp(bondk1.kb[i],pt) == 0))
        {
            return (bondk1.t[i]);
        }
    }
    return 1.53;  // default to c-c single
}

