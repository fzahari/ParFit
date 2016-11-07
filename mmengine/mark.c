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
#include "pot.h"
#include "field.h"
#include "piatomk.h"
#include "mark.h"

      
// ==================================================
int is_piatom(int ia)
{
    int i, it;

    if (field.type == MMX )
    {
        it = atom[ia].mmx_type;
        for (i=0; i < piatomk.npiatom; i++)
        {
           if (it == piatomk.kat[i])
              return TRUE;
        }
    } else if (field.type == MM3)
    {
        it = atom[ia].mm3_type;
        for (i=0; i < piatomk.npiatom; i++)
        {
           if (it == piatomk.kat[i])
              return TRUE;
        }
    } else
    {
        it = atom[ia].mmx_type;

        if( it ==  2 || it ==  3 || it ==  4 || it ==  6 || it == 7 || it ==  9 || it == 10 || it == 42 ||
           it == 37 || it == 41 || it == 46 || it == 47 || it == 48 || it == 29 || 
           it == 30 || it == 15 || it == 38 || it == 66 || it == 57 )
        return TRUE;
    }
    return FALSE;
}  
/* -------------------------------------- */
void  pimarkselection()
{
  return;
}
// ============================================
void bbchk()
{
        int back, i, i1, i2, i3, 
         iat2k, ibib, ibond1, ii, ijj, ip1, ip2, ip3, iatk, 
         ip4, ipf1, ipf1n, ipf4, ipf4n, it, j, jii, jjldx, k, k1, 
         k2, kdx, kk, kt, l, ldx, lev, maxlev, mnlev, 
         mnp1, mnp2, mxnode, nb, node, path[200][2];
        int found;
        long int mask;
        int *iara, **iarb;  


        /*     iara flag =1(aromatic atom) =0(non-aromatic atom)
         *     iarb flag =1(aromatic bond) =2(biphenyl-like bond)
         *     =3(butadiene-like bond) =0(other)
         *     call this subroutine from ce1(link) by making statement #316
         *     the call statement.
         *     use iarb flags in kbond and komega where appropriate
         *     one can potentially use this subr. to set logary flags
         *     for use in model
         * */

        iara = ivector(0,natom+1);
        iarb = imatrix(0,3,0,natom+1);
        
        mask = 1L << PI_MASK;
        for( i = 1; i <= natom; i++ )
        {
                atom[i].flags &= ~mask;
        }
        for( i = 0; i <= natom; i++ )
        {
                iara[i] = 0;
                for( j = 0; j < 3; j++ )
                {
                        iarb[j][i] = 0;
                }
        }

        nb = 0;
        for( i = 1; i <= natom; i++ )
        {
                it = atom[i].mmx_type;
                if (is_piatom(i) )
                {
                        for( j = 0; j < MAXIAT; j++ )
                        {
                                if( atom[i].iat[j] && i < atom[i].iat[j] )
                                {
                                        i1 = i;
                                        i2 = atom[i].iat[j];
                                        i3 = atom[i].bo[j];
                                        it = atom[i1].mmx_type;
                                        kt = atom[i2].mmx_type;
                                        /*  NOW CHECK FOR AROMATIC RINGS INVOLVING ATOMS I1 AND I2 */
                                        path[0][0] = i1;
                                        path[1][0] = i2;
                                        maxlev = 20;
                                        mnlev = 2;
                                        /* changed from maxiat to 3 */
                                        mxnode = 3;
                                        mnp1 = mxnode + 1;
                                        mnp2 = mxnode + 2;
                                        for( lev = 0; lev < 6; lev++ )
                                        {
                                                path[lev][1] = 0;
                                        }
                                        for( lev = 2; lev < 6; lev++ )
                                        {
                                                path[lev][0] = 0;
                                        }
                                        lev = 2;
                                        while( True )
                                        {
                                                if( lev <= maxlev )
                                                {
                                                        kdx = path[lev][0];
                                                        node = path[lev][1] + 1;
                                                        while( node <= mxnode )
                                                        {
                                                                ldx = atom[kdx].iat[node];
                                                                if( ldx > 0 )
                                                                {
                                                                        jjldx = 0;
                                                                        /* changed maxiat to 3 */
                                                                        for( ijj = 0; ijj < 3; ijj++ )
                                                                        {
                                                                                if( atom[ldx].iat[ijj] && atom[ldx].bo[ijj] )
                                                                                        jjldx += 1;
                                                                        }
                                                                        if( jjldx == 1 )
                                                                                ldx = 0;
                                                                }
                                                                it = 0;
                                                                if( ldx )
                                                                        it = atom[ldx].mmx_type;
                                                                if( it == 2 || it == 37 || it == 40 || it == 41 )
                                                                {
                                                                        ldx = ldx;
                                                                }
                                                                else
                                                                {
                                                                        ldx = 0;
                                                                }
                                                                if( !ldx )
                                                                {
                                                                        node += 1;
                                                                }
                                                                else
                                                                {
                                                                        back = lev - 1;
                                                                        while( back > 0 && ldx != path[back][0] )
                                                                        {
                                                                                back -= 1;
                                                                        }
                                                                        if( !back )
                                                                        {
                                                                                path[lev][1] = node;
                                                                                lev += 1;
                                                                                path[lev][0] = ldx;
                                                                                node = mnp2;
                                                                        }
                                                                        else if( back != 1 )
                                                                        {
                                                                                node += 1;
                                                                                /*            IF(LEV.GE.5) THEN */
                                                                        }
                                                                        else if( lev == 6 )
                                                                        {
                                                                                /*       AROMATIC RING FOUND */
                                                                                goto L_70;
                                                                        }
                                                                        else
                                                                        {
                                                                                node += 1;
                                                                        }
                                                                }
                                                        }
                                                        if( node == mnp1 )
                                                        {
                                                                /*       IF(LEV.EQ.MNLEV) GOTO 100  ! ACYCLIC OR NON-AROM. SYSTEM */
                                                                if( lev == mnlev )
                                                                        goto L_181;
                                                                path[lev][0] = 0;
                                                                path[lev][1] = 0;
                                                                lev -= 1;
                                                        }
                                                }
                                                else
                                                {
                                                        path[lev][0] = 0;
                                                        path[lev][1] = 0;
                                                        lev -= 1;
                                                }
                                        }

                                        /*c     AROMATIC RING FOUND */
L_70:
                                        ;
                                        for( k = 0; k <= lev; k++ )
                                        {
                                                k1 = path[k][0];
                                                k2 = path[k + 1][0];
                                                iara[k1] = 1;
                                                iara[k2] = 1;
                                                iarb[0][nb] = k1;
                                                iarb[1][nb] = k2;
                                                iarb[2][nb] = 1;
                                                nb += 1;
                                                if( atom[k1].mmx_type != 40 )
                                                        atom[k1].flags |= mask;
                                                if( atom[k2].mmx_type != 40 )
                                                        atom[k2].flags |= mask;
                                        }
                                        k1 = path[0][0];
                                        k2 = path[lev][0];
                                        iarb[0][nb] = k1;
                                        iarb[1][nb] = k2;
                                        iarb[2][nb] = 1;
                                        nb += 1;
                                }
L_181:
                                ;
                        }
                }
        }

        for( i = 1; i <= natom; i++ )
        {
                it = atom[i].mmx_type;
                if (is_piatom(i) )
                {
                        for( j = 0; j < MAXIAT; j++ )
                        {
                                if( (atom[i].iat[j] && atom[i].bo[j] != 9) && 
                                 i < atom[i].iat[j] )
                                {
                                        /*  search iarb array to see if this bond has been found before */
                                        for( k = 0; k < nb; k++ )
                                        {
                                                if( iarb[0][k] == i && iarb[1][k] == atom[i].iat[j] )
                                                        goto L_182;
                                        }
                                        i1 = i;
                                        i2 = atom[i].iat[j];
                                        i3 = atom[i].bo[j];
                                        it = atom[i1].mmx_type;
                                        kt = atom[i2].mmx_type;
                                        if( it == 2 && kt == 2 && i3 == 1 )
                                        {

                                                /*    HAVE NOW IDENTIFIED A TYPE 2-TYPE 2 SINGLE BOND NOT IN AN
                                                 *   AROMATIC RING */
                                                if( iara[i1] == 1 && iara[i2] == 1 )
                                                {
                                                        /*         BIPENYL-LIKE BOND (2-2 SINGLE BOND BETWEEN 2 AROM RINGS) */
                                                        iarb[2][k] = 2;
                                                }
                                                else
                                                {
                                                        iarb[2][k] = 3;
                                                }
                                        }
                                }
L_182:
                                ;
                        }
                }
        }
         /*         CHECK FOR BUTADIENE LIKE PI ATOM SEQUENCES   * */
        for( i = 1; i <= natom; i++ )
        {
                for( j = 0; j < MAXIAT; j++ )
                {
                        if( atom[i].iat[j] && atom[i].bo[j] != 9 )
                        {
                                i1 = i;
                                i2 = atom[i].iat[j];
                                i3 = atom[i].bo[j];
                                ip2 = atom[i1].mmx_type;
                                ip3 = atom[i2].mmx_type;
                                ipf1 = 0;
                                ipf4 = 0;
                                ipf1n = 0;
                                ipf4n = 0;
                                if( i3 == 1 )
                                {
                                    /* includes sulfur */
                                   if( is_piatom(i1) && is_piatom(i2) )
                                   {
                                      for( k = 0; k < MAXIAT; k++ )
                                      {
                                          if( atom[i1].iat[k] )
                                          {
                                             ip1 = atom[atom[i1].iat[k]].mmx_type;
                                             if( ip1 )
                                             {
                                                 if( atom[i1].iat[k] != i2 )
                                                 {
                                                    iatk = atom[i1].iat[k];
                                                    if( ibond(iatk,i1) )
                                                    {
                                                        /* sulfur */
                                                        if(is_piatom(iatk) && (ibond( iatk, i1 ) > 1) )
                                                           goto L_184;
                                                    }
                                                 }
                                             }
                                          }
                                      }
                                      goto L_183;
L_184:
                                      ipf1 = 1;
                                      ipf1n = atom[i1].iat[k];
L_183:
                                      for( kk = 0; kk < MAXIAT; kk++ )
                                      {
                                          iat2k = atom[i2].iat[kk];
                                          if( iat2k && atom[i2].bo[kk] != 9 )
                                          {
                                              ip4 = atom[iat2k].mmx_type;
                                              if( ip4 )
                                              {
                                                 if( iat2k != i1 )
                                                 {
                                                    ibond1 = ibond( iat2k, i2 );
                                                    ibib = ibond1;
                                                    if( is_piatom(iat2k) && (ibib > 1) )
                                                           goto L_186;
                                                 }
                                              }
                                          }
                                      }
                                      goto L_185;
L_186:
                                      ipf4 = 1;
                                      ipf4n = atom[i2].iat[kk];
                                   }
L_185:
                                   if( (ipf1 == 1 && ipf4 == 1) && i3 == 1 )
                                   {
                                                atom[ipf1n].flags |= mask;
                                                atom[i1].flags |= mask;
                                                atom[i2].flags |= mask;
                                                atom[ipf4n].flags |= mask;
                                        }
                                }
                        }
                }
        }

        /* Check for conjugated c+, c. c- */
        for( i = 1; i <= natom; i++ )
        {
             if( atom[i].mmx_type == 29 || atom[i].mmx_type == 30 || atom[i].mmx_type == 48 )
             {
                 for( k = 0; k < MAXIAT; k++ )
                 {
                     j = atom[i].iat[k];
                     if( j != 0 && atom[i].bo[k] != 9 )
                     {
                         ip4 = atom[j].mmx_type;
                         /* sulfur */
                         if( is_piatom(j) )
                         {
                            for( ii = 0; ii < MAXIAT; ii++ )
                            {
                                jii = atom[j].iat[ii];
                                if( jii && atom[j].bo[ii] != 9 )
                                {
                                    if( jii != i )
                                    {
                                        ip4 = atom[jii].mmx_type;
                                        if( is_piatom(jii) )
                                        {
                                            atom[i].flags |= mask;
                                            atom[j].flags |= mask;
                                        }
                                    }
                                }
                            }
                         }
                     }
                 }
             }
        }


        /*       CHECK FOR 1ST ATOM SUBSTITUENTS ON PI FLAGGED CONJUGATED ATOMS         * */
L_187:
        found = FALSE;        
        for( l = 1; l <= natom; l++ )
        {
            if( atom[l].flags & mask )
            {
                for( k = 0; k < MAXIAT; k++ )
                {
                    if( atom[l].iat[k] && atom[l].bo[k] <= 3 )
                    {
                         ip1 = atom[atom[l].iat[k]].mmx_type;
                         if( ip1 == 0 )
                            break;
                         /* sulfur */
                        if( is_piatom(atom[l].iat[k]) && ip1 != 40 && !(atom[atom[l].iat[k]].flags & mask) )
                        {
                             atom[atom[l].iat[k]].flags |= mask;
                             found = TRUE;
                        }
                    }
                }
            }

        }
        if (found == TRUE) goto L_187;
//
        // unmark allenes
        for (i=1; i <= natom; i++)
        {
            if (atom[i].mmx_type == 4)
            {
                if (atom[i].iat[0] != 0 && atom[i].iat[1] != 0)
                {
                    if (atom[i].bo[0] == 2 && atom[i].bo[1] == 2)
                    {
                        atom[i].flags &= ~mask;
                        atom[atom[i].iat[0]].flags &= ~mask;
                        atom[atom[i].iat[1]].flags &= ~mask;
                    }
                }
            }
        }

        free_imatrix(iarb, 0,3, 0, natom+1);
        free_ivector(iara, 0, natom+1);
        return;
}
/* =====================================  */
int ibond(int i1,int j1)
{
        int i, ib, ibond_v, ii, j;


        if( (!i1 || !j1) || i1 == j1 )
        {
                ibond_v = 0;
        }
        else
        {

                i = i1;
                j = j1;
                if( i > j )
                {
                        ii = i;
                        i = j;
                        j = ii;
                }
                for( ib = 0; ib < MAXIAT; ib++ )
                {
                        if( atom[i].iat[ib] == j )
                                goto L_1;
                }
                ibond_v = 0;
                goto L_2;
L_1:
                ibond_v = atom[i].bo[ib];
        }
L_2:
        return( ibond_v );
} 
/* =============================================  */
int checkbond(int i, int j)
{
        int loop;
        for (loop = 0; loop < MAXIAT; loop++)
        {
                if (atom[i].iat[loop] == j)
                        return loop;
        }
        return (-1);
}

