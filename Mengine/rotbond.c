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
#include "utility.h"
#include "rotbond.h"
        
void rotb(int ia1,int n1,int n2,int ia4, float degree)
{
        int issnum;
        long int mask;
        int i;
        float aa, cos1, cos2, denom, rot1[4][4], sin1, 
         sin2, trans1[4][4], x1, x2, xform1[4][4], xform2[4][4], y1, y2, 
         z1, z2, xtmp, ytmp, ztmp;

        /* subroutine to rotate about a bond. the atoms of the bond are n1 and n2,
         * the degree to rotate is idegree, and issnum is the substructure number
         * we get from chain to tell matxform which atoms to rotate
         * we assume chain has been called to set up the atoms to rotate
         * we use xaxis and xyplan to generate the rotation matrices to put
         * n1 at the origin and n2 on the x axis
         * */

        issnum = 1;
        mask = (1L << issnum);
        chain( n1, n2, issnum );

        x1 = atom[n1].x;
        y1 = atom[n1].y;
        z1 = atom[n1].z;
        x2 = atom[n2].x - x1;
        y2 = atom[n2].y - y1;
        z2 = atom[n2].z - z1;
        denom = sqrt( x2*x2 + y2*y2 );
        if( fabs( denom ) < .00001 )
        {
                cos1 = 1.0;
                sin1 = 0.0;
        }
        else
        {
                sin1 = -y2/denom;
                cos1 = x2/denom;
        }
        x2 = x2*cos1 - y2*sin1;
        denom = sqrt( x2*x2 + z2*z2 );
        if( fabs( denom ) < .00001 )
        {
                cos2 = 1.0;
                sin2 = 0.0;
        }
        else
        {
                sin2 = z2/denom;
                cos2 = x2/denom;
        }
        aa = degree - dihdrl( ia1, n1, n2, ia4 );

        i = 3;
        xtmp = -x1;
        ytmp = -y1;
        ztmp = -z1;
        mattrans( trans1, &xtmp, &ytmp, &ztmp );
        matrotsc( rot1, &i, &sin1, &cos1 );
        matriply( trans1, rot1, xform1 );

        i = 2;
        matrotsc( rot1, &i, &sin2, &cos2 );
        matriply( xform1, rot1, xform2 );

        i = 1;
        matrotat( rot1, &i, &aa );
        matriply( xform2, rot1, xform1 );

        i = 2;
        xtmp = -sin2;
        matrotsc( rot1, &i, &xtmp, &cos2 );
        matriply( xform1, rot1, xform2 );

        i = 3;
        xtmp = -sin1;
        matrotsc( rot1, &i, &xtmp, &cos1 );
        matriply( xform2, rot1, xform1 );

        mattrans( trans1, &x1, &y1, &z1 );
        matriply( xform1, trans1, xform2 );

        matxform( xform2, issnum );

        for( i = 1; i <= natom; i++ )
        {
                atom[i].substr[0] &= ~mask ;
        }
        return;
} 
/* ------------------------------------------ */
void chain(int ia,int ib,int issnum)
{
        int newmem;
        int i, ij, j;
        long int mask;

        /**** this routine starts at atom ia and fills
         *     iarray with the atoms attached directly
         *     or indirectly to ia but not attached either
         *     directly or indirectly to ib.  iarray is
         *     filled such that a non-zero value for the
         *     i-th member means that atom i is either
         *     directly or indirectly attached to ia. *** */
        /**** zero the output array *** */

        mask = 1L << issnum ;
        
        /**** find the atoms adjacent to ia other than atom ib *** */
        
        for( j = 0; j < MAXIAT; j++ )
        {
           if( atom[ia].iat[j])
           {
               if( atom[ia].iat[j] != ib )
               {
                  atom[atom[ia].iat[j]].substr[0] = 0 ;
                  atom[atom[ia].iat[j]].substr[0] |= mask ;
               }
           }
        }

        /**** repetitively loop through the output array
         *     looking for adjacent atoms until no new
         *     attachments are found *** */
        while( TRUE )
        {
           newmem = FALSE;
           for( i = 1; i <= natom; i++ )
           {
              if( atom[i].substr[0] & mask)   // atom not previously marked
              {
                 for( j = 0; j < MAXIAT; j++ )
                 {
                      if( atom[i].iat[j] != 0 &&  atom[i].iat[j] != ia && !(atom[atom[i].iat[j]].substr[0] & mask))
                      {
                         /**** a new atom has been found *** */
                           atom[atom[i].iat[j]].substr[0] = 0 ;
                           atom[atom[i].iat[j]].substr[0] |= mask ;
                           newmem = TRUE;
                      }
                 }
              }
           }
           /**** check for newly added atoms *** */
           if( !newmem )
              break;
        }
        if( ib )
        {
                if( atom[ib].substr[0] & mask)
                {
                        atom[ib].substr[0] &= ~mask;
                        for( ij = 0; ij < MAXIAT; ij++ )
                        {
                                if( atom[ib].iat[ij] )
                                {
//                                        if( atom[ib].bo[ij] != 9 )
                                        {
                                                if( atom[ib].iat[ij] != ia )
                                                   atom[atom[ib].iat[ij]].substr[0] &= ~mask;
                                        }
                                }
                        }
                }
        }
        return;
} 

/*     -------------------------------------- */
void matident(float matrix[][4])
{
        int i, j;

        for( i = 0; i < 4; i++ )
        {
            for( j = 0; j < 4; j++ )
            {
               if( i == j )
                  matrix[j][i] = 1.;
               else
                  matrix[j][i] = 0.;
            }
        }
        return;
} 
/*     -------------------------------------- */
void matrotsc(float matrix[][4],int *iaxis,float *sine,float *cose)
{

        matident( matrix );
        if( *iaxis == 1 ){
                matrix[1][1] = *cose;
                matrix[2][1] = *sine;
                matrix[1][2] = -*sine;
                matrix[2][2] = *cose;
                }
        else if( *iaxis == 2 ){
                matrix[0][0] = *cose;
                matrix[2][0] = -*sine;
                matrix[0][2] = *sine;
                matrix[2][2] = *cose;
                }
        else if( *iaxis == 3 ){
                matrix[0][0] = *cose;
                matrix[1][0] = *sine;
                matrix[0][1] = -*sine;
                matrix[1][1] = *cose;
                }
        return;
} 
/*     -------------------------------------- */
void matrotat(float matrix[][4],int *iaxis,float *angle)
{
        float cose, sine;

        matident( matrix );
        sine = sin( -*angle*0.017453 );
        cose = cos( -*angle*0.017453 );
        if( *iaxis == 1 ){
                matrix[1][1] = cose;
                matrix[2][1] = sine;
                matrix[1][2] = -sine;
                matrix[2][2] = cose;
                }
        else if( *iaxis == 2 ){
                matrix[0][0] = cose;
                matrix[2][0] = -sine;
                matrix[0][2] = sine;
                matrix[2][2] = cose;
                }
        else if( *iaxis == 3 ){
                matrix[0][0] = cose;
                matrix[1][0] = sine;
                matrix[0][1] = -sine;
                matrix[1][1] = cose;
                }
        return;
} 
/*     -------------------------------------- */
void mattrans(float matrix[][4],float *x,float *y,float *z)
{

        matident( matrix );
        matrix[0][3] = *x;
        matrix[1][3] = *y;
        matrix[2][3] = *z;
        return;
} 
/*     -------------------------------------- */
void matriply(float m1[][4],float m2[][4],float m3[][4])
{
        int i, j, k;

        for( i = 0; i < 4; i++ ){
                for( j = 0; j < 4; j++ ){
                        m3[j][i] = 0.;
                        for( k = 0; k < 4; k++ ){
                                m3[j][i] = m3[j][i] + m1[k][i]*m2[j][k];
                                }
                        }
                }
        return;
} 
/*     -------------------------------------- */
void matxform(float matrix[][4],int issnum)
{
    int i, k;
    long int mask;
    float x1[3];

    mask = 0;
    if (issnum != -1)
    mask = 1L << issnum;

    for( i = 1; i <= natom; i++ )
    {
       if ( (atom[i].substr[0] & mask) || issnum == -1)
       {
           for( k = 0; k < 3; k++ )
           {
                 x1[k] = atom[i].x*matrix[k][0] + atom[i].y*matrix[k][1] + atom[i].z*
                     matrix[k][2] + matrix[k][3];
            }
            atom[i].x = x1[0];
            atom[i].y = x1[1];
            atom[i].z = x1[2];
       }
    }
    return;
}

