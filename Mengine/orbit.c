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
#include "pot.h"
#include "pimatx.h"
#include "field.h"
#include "pistuf.h"
#include "pisave.h"
#include "pidata.h"
#include "pibonds.h"
#include "piatomk.h"
#include "orbit.h"

/*   listpi  - list of pi atom numbers                maxatom
     pibond  - index into bonds list, at1 and at2     maxpibond x 3
     pitor   - index into tors list, and ib           maxpitor x 2
     piplane - three atoms making plane about pi atom maxpi x 3
     dynamically allocate memory here - deallocate before quitting  */

void free_orbit()
{
    float morb2;
    
     free_matrix(pimatx.hc ,0, pistuf.norb, 0, pistuf.norb);
     free_matrix(pimatx.v  ,0, pistuf.norb, 0, pistuf.norb);
     free_matrix(pimatx.f  ,0, pistuf.norb, 0, pistuf.norb);
     free_matrix(pimatx.g  ,0, pistuf.norb, 0, pistuf.norb);
     free_matrix(pimatx.gamma ,0, pistuf.norb, 0, pistuf.norb);
     free_vector(pimatx.en  ,0,pistuf.norb);
     free_vector(pimatx.zor ,0,pistuf.norb);
     free_vector(pimatx.zt  ,0,pistuf.norb);
     free_vector(pimatx.zz  ,0,pistuf.norb);
     free_vector(pimatx.emz ,0,pistuf.norb);
     free_vector(pimatx.w   ,0,pistuf.norb);
     free_vector(pimatx.w1  ,0,pistuf.norb);
     free_vector(pimatx.w2  ,0,pistuf.norb);
     free_vector(pimatx.w3  ,0,pistuf.norb);
     free_vector(pimatx.q   ,0,pistuf.norb);
     free_ivector(pimatx.qp ,0,pistuf.norb);
     free_vector(pimatx.eedd,0,pistuf.norb);
     free_vector(pimatx.edd ,0,pistuf.norb);

     free_vector(local.da ,0, pistuf.norb);
     free_vector(local.adnmo ,0, pistuf.norb);
     free_vector(local.adnmt ,0, pistuf.norb);
     free_vector(local.d12 ,0, pistuf.norb);
     free_vector(local.d13 ,0, pistuf.norb);
     free_vector(local.fn ,0, pistuf.norb);
     free_vector(local.f22 ,0, (pistuf.norb*(pistuf.norb+1)) );
     free_vector(local.f27 ,0, (pistuf.norb*(pistuf.norb+1)) );
     free_vector(local.f29 ,0, (pistuf.norb*(pistuf.norb+1)) );
     free_vector(local.f1 ,0, (pistuf.norb*(pistuf.norb+1)) );
     free_vector(local.fa1 ,0, (pistuf.norb*(pistuf.norb+1)) );

    free_vector(pitorsave.tcsave ,0,pistuf.ntpi);
    free_ivector(pistuf.listpi,  0, pistuf.norb);
    free_imatrix(pistuf.piplane, 0, pistuf.norb, 0, 3);
    free_imatrix(pistuf.pibond,  0, pistuf.nbpi, 0, 3);    
    //    free_imatrix(pistuf.pitor,   0, pistuf.ntpi, 0, 2);
    free_ivector(pibondsave.ibsave ,0,pistuf.nbpi);
    free_ivector(pibondsave.kbsave ,0,pistuf.nbpi);
    free_ivector(pibondsave.nbsave ,0,pistuf.nbpi);
    free_vector(pibondsave.bksave ,0,pistuf.nbpi);
    free_vector(pibondsave.blsave ,0,pistuf.nbpi);
    free_vector(pibondsave.pksave ,0,pistuf.nbpi);
    free_vector(pibondsave.plsave ,0,pistuf.nbpi);
    free_vector(pibondsave.plsave2 ,0,pistuf.nbpi);
    free_ivector(pitorsave.itsave ,0,pistuf.ntpi);
    free_ivector(pitorsave.ktsave ,0,pistuf.ntpi);
    free_ivector(pitorsave.ncsave ,0,pistuf.ntpi);
    morb2 = pistuf.norb*pistuf.norb;
    free_vector(pibonds.riks ,0, morb2);
    free_vector(pibonds.a1s ,0, morb2);
    free_vector(pibonds.b1s ,0, morb2);
    free_vector(pibonds.c1s ,0, morb2);
    free_vector(pibonds.a2s ,0, morb2);
    free_vector(pibonds.b2s ,0, morb2);
    free_vector(pibonds.c2s ,0, morb2);
    free_matrix(pisave.fjsave ,0, pistuf.norb, 0, pistuf.norb);
    free_matrix(pisave.fksave ,0, pistuf.norb, 0, pistuf.norb);
    free_matrix(pisave.gjsave ,0, pistuf.norb, 0, pistuf.norb);
    free_matrix(pisave.gksave ,0, pistuf.norb, 0, pistuf.norb);
    free_matrix(penint.pen ,0, pistuf.norb, 0, pistuf.norb);
    free_vector(penint.ed ,0, pistuf.norb);
    free_matrix(mtape.vmtape ,0, pistuf.norb, 0, pistuf.norb);
    free_matrix(mtape.vbtape ,0, pistuf.norb, 0, pistuf.norb);
    free_matrix(mtape.pensave ,0, pistuf.norb, 0, pistuf.norb);
    free_matrix(densit.den ,0, pistuf.norb+1, 0, pistuf.norb);
    free_matrix(densit.dena ,0, pistuf.norb+1, 0, pistuf.norb);
    free_matrix(densit.denb ,0, pistuf.norb+1, 0, pistuf.norb);

    free_vector(ed123.ed1, 0, pistuf.norb);
    free_vector(ed123.ed2, 0, pistuf.norb);
    free_vector(ed123.ed3, 0, pistuf.norb);
    free_vector(ed123.ed4, 0, pistuf.norb);
    free_vector(ed123.ed5, 0, pistuf.norb);
     
    pistuf.norb = 0;
    pistuf.nbpi = 0;
    pistuf.ntpi = 0;

}

void orbital()
{
    int i,ia,ib,j, found;
    int iorb;
    long int mask;
    float morb2;

    mask = (1L << 0);
    
     pot.use_picalc = FALSE;
     if (field.type == MMFF94 || field.type == AMBER)
       return;
     pistuf.norb = 0;
     pistuf.nbpi = 0;
     pistuf.ntpi = 0;
// count pi atoms
     for (i=1; i <= natom; i++)
     {
         if (atom[i].flags & mask)  // marked as pi atom check type is ok
         {
             found = FALSE;
             for (j=0; j < piatomk.npiatom; j++)
             {
                 if (atom[i].type == piatomk.kat[j])
                     found = TRUE;
             }
             if (found == TRUE)
                pistuf.norb++;
             else
                atom[i].flags &= ~mask;
         }
     }
     if (pistuf.norb == 0 )
       return;

     uvnpi.ifirst = TRUE;
     uvnpi.kfirst = TRUE;     
     
     if (pistuf.norb > 0)
       pot.use_picalc = TRUE;
       
     pistuf.listpi = ivector(0, pistuf.norb);
     pistuf.piplane = imatrix(0, pistuf.norb, 0, 3);
     for (i=0; i < pistuf.norb; i++)
        pistuf.listpi[i] = 0;

     pimatx.hc = matrix(0, pistuf.norb, 0, pistuf.norb);
     pimatx.v  = matrix(0, pistuf.norb, 0, pistuf.norb);
     pimatx.f  = matrix(0, pistuf.norb, 0, pistuf.norb);
     pimatx.g  = matrix(0, pistuf.norb, 0, pistuf.norb);
     pimatx.gamma = matrix(0, pistuf.norb, 0, pistuf.norb);
     pimatx.en = vector(0,pistuf.norb);
     pimatx.zor = vector(0,pistuf.norb);
     pimatx.zt = vector(0,pistuf.norb);
     pimatx.zz = vector(0,pistuf.norb);
     pimatx.emz = vector(0,pistuf.norb);
     pimatx.w  = vector(0,pistuf.norb);
     pimatx.w1 = vector(0,pistuf.norb);
     pimatx.w2 = vector(0,pistuf.norb);
     pimatx.w3 = vector(0,pistuf.norb);
     pimatx.q = vector(0,pistuf.norb);
     pimatx.qp = ivector(0,pistuf.norb);
     pimatx.eedd = vector(0,pistuf.norb);
     pimatx.edd = vector(0,pistuf.norb);
         local.da = vector(0, pistuf.norb);
         local.adnmo = vector(0, pistuf.norb);
         local.adnmt = vector(0, pistuf.norb);
         local.d12 = vector(0, pistuf.norb);
         local.d13 = vector(0, pistuf.norb);
         local.fn = vector(0, pistuf.norb);
         
         local.f22 = vector(0, (pistuf.norb*(pistuf.norb+1)) );
         local.f27 = vector(0, (pistuf.norb*(pistuf.norb+1)) );
         local.f29 = vector(0, (pistuf.norb*(pistuf.norb+1)) );
         local.f1 = vector(0, (pistuf.norb*(pistuf.norb+1)) );
         local.fa1 = vector(0, (pistuf.norb*(pistuf.norb+1)) );
     
     morb2 = pistuf.norb*pistuf.norb;
     pibonds.riks = vector(0, morb2);
     pibonds.a1s = vector(0, morb2);
     pibonds.b1s = vector(0, morb2);
     pibonds.c1s = vector(0, morb2);
     pibonds.a2s = vector(0, morb2);
     pibonds.b2s = vector(0, morb2);
     pibonds.c2s = vector(0, morb2);
     pisave.fjsave = matrix(0, pistuf.norb, 0, pistuf.norb);
     pisave.fksave = matrix(0, pistuf.norb, 0, pistuf.norb);
     pisave.gjsave = matrix(0, pistuf.norb, 0, pistuf.norb);
     pisave.gksave = matrix(0, pistuf.norb, 0, pistuf.norb);
     mtape.vmtape = matrix(0, pistuf.norb, 0, pistuf.norb);
     mtape.vbtape = matrix(0, pistuf.norb, 0, pistuf.norb);
     mtape.pensave = matrix(0, pistuf.norb, 0, pistuf.norb);
     penint.pen = matrix(0, pistuf.norb, 0, pistuf.norb);
     penint.ed = vector(0, pistuf.norb);

     densit.den = matrix(0, pistuf.norb+1, 0, pistuf.norb);
     densit.dena = matrix(0, pistuf.norb+1, 0, pistuf.norb);
     densit.denb = matrix(0, pistuf.norb+1, 0, pistuf.norb);

     ed123.nscft = 0;
     ed123.ed1 = vector(0, pistuf.norb);
     ed123.ed2 = vector(0, pistuf.norb);
     ed123.ed3 = vector(0, pistuf.norb);
     ed123.ed4 = vector(0, pistuf.norb);
     ed123.ed5 = vector(0, pistuf.norb);

     for (i=0; i < pistuf.norb; i++)
     {
         ed123.ed1[i] = 0.0;
         ed123.ed2[i] = 0.0;
         ed123.ed3[i] = 0.0;
         ed123.ed4[i] = 0.0;
         ed123.ed5[i] = 0.0;         
         for (j=0; j < pistuf.norb; j++)
            pimatx.f[i][j] = 0.0;
     }
     iorb = 0;
     for (i=1; i <= natom; i++ )
     {
         if (atom[i].flags & mask )
         {
             pistuf.listpi[iorb] = i;
             iorb++;
         }
     }
     piplane( pistuf.piplane);

// find pi bonds
     pistuf.nbpi = 0;
     for (i=0; i < bonds_ff.nbnd; i++)
     {
         ia = bonds_ff.i12[i][0];
         ib = bonds_ff.i12[i][1];
         if ( (atom[ia].flags & mask) && (atom[ib].flags & mask) )
           pistuf.nbpi++;
     }
     if (pistuf.nbpi > 0)
     {
         iorb = 0;
         pistuf.pibond = imatrix(0, pistuf.nbpi, 0, 3);
         pibondsave.ibsave = ivector(0,pistuf.nbpi);
         pibondsave.kbsave = ivector(0,pistuf.nbpi);
         pibondsave.nbsave = ivector(0,pistuf.nbpi);
         pibondsave.bksave = vector(0,pistuf.nbpi);
         pibondsave.blsave = vector(0,pistuf.nbpi);
         pibondsave.pksave = vector(0,pistuf.nbpi);
         pibondsave.plsave = vector(0,pistuf.nbpi);
         pibondsave.plsave2 = vector(0,pistuf.nbpi);
         
         for (i = 0; i < bonds_ff.nbnd; i++)
         {
             ia = bonds_ff.i12[i][0];
             ib = bonds_ff.i12[i][1];
             if ( (atom[ia].flags & mask) && (atom[ib].flags & mask) )
             {
                 pistuf.pibond[iorb][0] = i;
                 pistuf.pibond[iorb][1] = ia;
                 pistuf.pibond[iorb][2] = ib;
                 iorb++;
             }
         }
     }
// find pi torsions
     pistuf.ntpi = 0;
     for (i=0; i < torsions.ntor; i++)
     {
         ia = torsions.i14[i][1];
         ib = torsions.i14[i][2];
         if ( (atom[ia].flags & mask) && (atom[ib].flags & mask) )
           pistuf.ntpi++;
     }
     if (pistuf.ntpi > 0)
     {
         iorb = 0;
	 // pistuf.pitor = imatrix(0, pistuf.ntpi, 0, 2);
         pitorsave.itsave = ivector(0,pistuf.ntpi);
         pitorsave.ktsave = ivector(0,pistuf.ntpi);
         pitorsave.ncsave = ivector(0,pistuf.ntpi);
         pitorsave.tcsave = vector(0,pistuf.ntpi);

         for (i=0; i < pistuf.ntpi; i++)
           pitorsave.tcsave[i] = 0.0;
           
	 /* for (i=0; i < torsions.ntor; i++)
         {
             ia = torsions.i14[i][1];
             ib = torsions.i14[i][2];
             if ( (atom[ia].flags & mask) && (atom[ib].flags & mask) )
             {
                 pistuf.pitor[iorb][0] = i;
                 pistuf.pitor[iorb][1] = ia;
                 iorb++;
             }
	     }*/
     }
     // find groups of pi atoms
     find_pigroups();            
// resort listpi
    if (pigroup.nfrag > 1)
    {
       for (i=0; i < pistuf.norb; i++)
          pistuf.listpi[i] = pigroup.fraglist[i];
    }
}
// =============================================
//  assign piatom parameters 
void korbit()
{
    int i, nfill, ia, it, j, iclass;

    nfill = 0;
    for (i=0; i < pistuf.norb; i++)
    {
       ia = pistuf.listpi[i];
       it = atom[ia].type;
       iclass = 0;
       if (field.type == MM3)
          iclass = atom[ia].tclass;

       for (j=0; j < piatomk.npiatom; j++)
       {
           if (it == piatomk.kat[j] || iclass == piatomk.kat[j])
           {
               pimatx.zor[i] = piatomk.zor[j]/2.0;
               if (it == 2 || it == 3 || it == 4 || it == 29 || it == 30 || it == 48)
                   pimatx.zor[i] = zetac(ia,it)/2.0;
               pimatx.zz[i]  = piatomk.zz[j];
               pimatx.emz[i] = piatomk.emz[j];
               pimatx.w1[i]  = piatomk.w1[j];
               pimatx.w2[i]  = piatomk.w2[j];
               pimatx.w3[i]  = piatomk.w3[j];
               pimatx.q[i]   = piatomk.q[j];
               pimatx.qp[i]  = piatomk.qp[j];
               nfill += pimatx.q[i];
               if (it == 30 || it == 48 || it == 42)
                  pimatx.q[i] = 1.0;
               if (it == 66  )
                  pimatx.q[i] = 1.0 ;
               pigroup.nelecs[pigroup.grlist[i]-1] += pimatx.q[i];
               break;
           }
       }
    }
    pimatx.npiel = nfill/2;
}
// ====================================
float zetac(int i,int it)
{
    int ii, iit, l, nc, nh, nonpi, mask;
    float zetac_v;
        /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         *
         *  CALCULATES ORBITAL EXPONENTS (ZOR) FOR CARBON ATOMS WHICH ARE 
         *  ATTACHED TO H'S AND OTHER C'S ONLY. 
         *      SEE ALLINGER, JACS, 87, 2081 (1965). 
         *
         * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         * */


        zetac_v = 3.250;
        nonpi = 1;
        nc = 0;
        nh = 0;
        mask = 1L << 0;
        for( l = 0; l < MAXIAT; l++ )
        {
            ii = atom[i].iat[l];
            if(  ii && atom[i].bo[l] != 9 )
            {
               if( !( atom[ii].flags & mask) )
               {
                  nonpi += 1;
                  iit = atom[ii].type;
                  if( iit == 1 )
                  {
                    nc += 1;
                  }else if( !(iit >= 1 && iit <= 4) )
                  {
                     if( iit == 5 )
                        nh += 1;
                  }
                }
             }
        }
        if( it == 4 )
        {
           if( nh == 1 )
           {
       /*  ACETYLENE TYPE CSP  (1H)  */
              zetac_v = 3.056;
           }else if( nc == 1 )
           {
       /*  2-butyne type csp */
              zetac_v = 2.946;
            }
        }else if( nonpi == 2 )
        {
           if( nh == 1 )
           {
       /* benzene */
              zetac_v = 3.183;
            }else if( nc == 1 )
            {
       /* alkyl benzene */
               zetac_v = 3.073;
             }
         }else if( nonpi == 3 )
         {
             if( nh == 2 )
             {
        /* ethylene */
                zetac_v = 3.132;
             }else if( nc == 2 )
             {
        /* dialkyl sp2 carbon */
                zetac_v = 2.938;
              }else if( nc == 1 && nh == 1 )
              {
        /* mono h, mono alkyl sp2 carbon */
                 zetac_v = 3.022;
              }
        }
        return( zetac_v );
}
// ========================================
void piplane(int **piplane)
{
    int i,ia, ib, icnt, ig,il, it, j, loop, kt, k;

    for (i=0; i < pistuf.norb; i++)
    {
        il = pistuf.listpi[i];
        it = atom[il].type;
        if (atom[il].atomnum == 6)
          it = 2;
        else if (atom[il].atomnum == 7)
          it = 9;
        else if (atom[il].atomnum == 8)
          it = 7;
        icnt = 0;
        for (j=0; j < MAXIAT; j++)
        {
            if (atom[il].iat[j] != 0 && atom[il].bo[j] != 9)
            {
                piplane[i][icnt] = atom[il].iat[j];
                icnt++;
                if (icnt > 3)
                  break;
            }
        }
        if (icnt == 1)
        {
            piplane[i][1] = il;
            if (it == 2 || it == 3 || it == 6 || it == 8 || it == 9)
            {
                ia = piplane[i][0];
                if (atom[ia].iat[0] != il && atom[ia].iat[0] != 0)
                   piplane[i][2] = atom[ia].iat[0];
                else if (atom[ia].iat[1] != il && atom[ia].iat[1] != 0)
                   piplane[i][2] = atom[ia].iat[1];
                else
                {
                   message_alert("Error in planes 1","Error");
                   return;
                }
            }
            else if (it == 4)  // acetylene
            {
                ia = piplane[i][0];
                loop = 0;
L_10: 
                kt = atom[ia].type;
                if (it == 1 || it == 5 || it == 10)
                {
                    ia = piplane[i][1];
                    loop++;
                    if (loop >= 2)
                    {
                        message_alert("Error in planes 2","Error");
                        return;
                    }
                    goto L_10;
                }else if (kt == 2 || kt == 3 || kt == 6 || kt == 8 || kt == 9)
                {
                   piplane[i][2] = atom[ia].iat[0];
                   if (piplane[i][2] == il)
                     piplane[i][2] = atom[ia].iat[1];
                   break;
                } else if (kt == 4)
                {
                    ia = piplane[i][1];
                    loop++;
                     if (loop >= 2)
                    {
                        message_alert("Error in planes 3","Error");
                        return;
                    }
                    goto L_10;
                } else if (kt == 7)
                {                    
                    message_alert("Error in planes 4","Error");
                    return;
                }
            } else if (it == 7)  // carbonyl oxygen
            {
                ia = il;  // oxygen
                ib = atom[ia].iat[0];
                piplane[i][0] = il;
                piplane[i][1] = ib;
                for(k=0; k < MAXIAT; k++)
                {
                    if (atom[ib].iat[k] != 0 && atom[ib].iat[k] != il)
                    {
                        piplane[i][2] = atom[ib].iat[k];
                        break;
                    }
                }
            } else if (it == 10)
            {
                ia = piplane[i][1];
            } else
            {
                    message_alert("Error in planes 5","Error");
                    return;
            }
            ib = atom[ia].iat[0];
            if (ib == il)
               ib = atom[ia].iat[1];
            kt = atom[ib].type;
            if (it == 2 || it == 3 || it == 6 || it == 8 || it == 9)
            {
                ig = atom[ib].iat[0];
                if (ig == ia)
                    ig = atom[ib].iat[1];
                piplane[i][2] = ig;
            }
        } else if (icnt == 2)
           piplane[i][2] = il;
    }             
}
// ===========================================
void find_pigroups()
{
    int i,found, ia,ib, ncount,j,k, repeat,icount;
    int *used;

    pigroup.nfrag = 0;
    icount = 0;
    used = ivector(0,pistuf.norb);
//
    for (i=0; i < MAXORB/4; i++)
    {
        pigroup.nao[i] = 0;
        pigroup.nelecs[i] = 0;
    }
    for (i=0; i < pistuf.norb; i++)
    {
        used[i] = FALSE;
        pigroup.fraglist[i] = 0;
        pigroup.grlist[i] = 0;
    }

//  start
     ia = pistuf.listpi[0];
     pigroup.fraglist[icount] = ia;
     pigroup.grlist[icount] = 1;
     icount++;
     pigroup.nao[0] = 1;
     used[0] = TRUE;
     
// find all attached to ia
     for (j=1; j < pistuf.norb; j++)
     {
         ib = pistuf.listpi[j];
         if (!used[j])
         {
             if (is_bond(ia,ib))
             {
                 pigroup.fraglist[icount] = ib;
                 pigroup.grlist[icount] = pigroup.nfrag+1;
                 icount++;
                 used[j] = TRUE;
                 pigroup.nao[0]++;
             }
         }
     }
// loop through list until all in group are found
LOOP1:
     found = FALSE;
     ncount = pigroup.nao[pigroup.nfrag];
     for (i=0; i < ncount; i++)
     {
         for (j=0; j < pistuf.norb; j++)
         {
             ia = pigroup.fraglist[j];
             for (k=0; k < pistuf.norb; k++)
             {
                 ib = pistuf.listpi[k];
                 if (ia != ib && !used[k])
                 {
                     if (is_bond(ia,ib))
                     {
                         pigroup.fraglist[icount] = ib;
                         pigroup.grlist[icount] = pigroup.nfrag+1;
                         icount++;
                         used[k] = TRUE;
                         pigroup.nao[pigroup.nfrag]++;
                         found = TRUE;
                     }
                 }
             }
         }
     }
     if (found == TRUE)
        goto LOOP1;
// any left
     repeat = FALSE;
     for (i=0; i < pistuf.norb; i++)
     {
         if (!used[i])
         {
             pigroup.nfrag++;
             ia = pistuf.listpi[i];
             pigroup.fraglist[icount] = ia;
             pigroup.grlist[icount] = pigroup.nfrag+1;
             icount++;
             used[i] = TRUE;
             pigroup.nao[pigroup.nfrag] = 1;
             repeat = TRUE;
             break;
         }
     }
     if (repeat == TRUE)
        goto LOOP1;
//
     free_ivector(used, 0, pistuf.norb);
    pigroup.nfrag++;
}       
// ================================
void initial_picalc()
{
    int i,j;
    pistuf.nplane = 1;
    pistuf.iuv = 1;
    pistuf.ihuck = 0;
    uvnpi.ifirst = TRUE;
     uvnpi.kfirst = TRUE;     
    for (i=0; i < pistuf.norb; i++)
    {
         ed123.ed1[i] = 0.0;
         ed123.ed2[i] = 0.0;
         ed123.ed3[i] = 0.0;
         ed123.ed4[i] = 0.0;
         ed123.ed5[i] = 0.0;         
         for (j=0; j < pistuf.norb; j++)
            pimatx.f[i][j] = 0.0;
    }
    piseq(-1,0);
}        
                      
