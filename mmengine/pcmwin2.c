#include "pcwin.h"
#include "pcmod.h"
#include "field.h"
#include "atom_k.h"
#include "minim_control.h"
#include "pcmwin2.h"
                
void hadd(void)
{
    type();
    if (field.type == MM3)
    {
        nhadd();
        hcoord();
	//        hdel(1);
    } else if (field.type == MMFF94)
    {
        nhadd();
        hcoord();
	//    hdel(1);
    } else
    {
        set_atomtypes(MMX);
        nhadd();
        hcoord();
    }
    type();
    generate_bonds();
        
}

void nhadd()
{
        int i, i1, ik, im, it, itads, j, jk, nhadds, newatom;
        int jji, bo,ncount;
        static int iadd[MAXVDW]={2,2,2,2,0,3,1,6,5,1,
                                  0,0,0,0,3,6,3,2,2,0,
                                  0,2,0,0,6,5,2,0,5,5,
                                  0,0,0,3,0,0,6,1,1,2,
                                  2,4,2,5,0,6,0,6,2,2,
                                  3,3,3,0,6,2,2,2,0,0,
                                  0,0,0,0,0,1,6,0,0,0,
                                  0,0,0,0,0};

        /*     this routine calculates the #'s of h to be added
         *     to carbon and packs it symbolically into iat(i,j) table
         *     as well as into the connection table via the appropriate
         *     manipulation instructions.  hydrogen coordinates
         *     are not added by this routine. lp are also added to oxygen. */
        /*     0 to disregard : 1 and 4 for lp's only
         *     2 and 5 for h's only and 3 and 6 for h and lp
         *     add 3 for trigonal atoms */

   ncount = natom;
   for( i = 1 ; i <= ncount; i++ )
   {
        if( atom[i].mmx_type > 0 && (atom[i].mmx_type < 300) && atom[i].atomnum != 1 ) 
        {
           jji = 0;
           bo = 0;
           for (j=0; j < MAXIAT; j++)
           {
              if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
              {
                 jji ++;
                 bo += atom[i].bo[j];
              }
           }
           nhadds = atom_k.ligands[atom[i].type] - jji;

           if (atom[i].mmx_type >= 11 && atom[i].mmx_type <= 14)
               nhadds = 0;
           if (atom[i].mmx_type == 41 || atom[i].mmx_type == 46)
             nhadds = atom_k.ligands[atom[i].type] - bo;
           if (atom[i].mmx_type == 49 || atom[i].mmx_type == 50 || atom[i].mmx_type == 51)
             nhadds = atom_k.ligands[atom[i].type] - bo;
           
/*  add atom and bond it to atom i      */
           if( nhadds > 0 )  
           {
              for( im = 1; im <= nhadds; im++ )
              {
                  newatom = make_atom(5,0.F,0.F,0.F,"H");
                  set_atomdata(newatom,5,5,5,34,29);
                  atom[newatom].molecule = atom[i].molecule;
                  atom[newatom].residue = atom[i].residue;
                  atom[newatom].biotype = atom[i].biotype;
                  if (newatom == -1)
                  {
                      message_alert("Error in nhadd newatom = -1","Error");
                      hdel(0);
                      return;
                  }
                  make_bond(i,newatom,1);
                  for (j=0; j < MAXSS; j++)
                      atom[newatom].substr[j] = atom[i].substr[j];   
              }
           }
        }
   }

/*       ***  check for lp's ***        */
   if (field.type == MMX)
   {     
      ncount = natom;
      for( i = 1 ; i <= ncount; i++ )
      {
           it = iadd[atom[i].mmx_type-1];         
           if( it >= 1 )
           {
               itads = atom[i].mmx_type ;
               if (it == 1 || it == 3 || it == 4 || it == 6 ) 
               {
                   nhadds = 4;
                   for( i1 = 0; i1 < MAXIAT; i1++ )
                   {
                       if( atom[i].bo[i1] != 0 && atom[i].bo[i1] != 9 )
                           nhadds = nhadds - atom[i].bo[i1];
                   }
                   if( itads == 6 || itads == 15 || itads == 34 || itads == 42 || itads == 53 || itads == 48 )
                   {
                        for( ik = 0; ik < MAXIAT; ik++ )
                        {
                            if( atom[i].iat[ik] != 0 && atom[i].bo[ik] != 9 )
                            {
                                jk = atom[atom[i].iat[ik]].mmx_type;
                                if( jk == 2 || jk == 3 || jk == 4 || jk == 9
                                 || jk == 10 || jk == 30 || jk == 37 || jk == 40 )
                                {
                                    nhadds = nhadds - 1;
                                    break;
                                }
                            }
                        }
                   }
                                
/*         correction for sulfoxide */
                   if( itads == 17 )
                       nhadds = nhadds + 1;

                   if (itads == 66)
                       nhadds = 2;

                   if( nhadds > 0 ) 
                   {
                      for( i1 = 1; i1 <= nhadds; i1++ )
                      {
                         newatom = make_atom(20, 0.F,0.F,0.F,"");
                         set_atomdata(newatom,20,0,0,0,0);
                         atom[newatom].molecule = atom[i].molecule;
                         atom[newatom].residue = atom[i].residue;
                         atom[newatom].biotype = atom[i].biotype;
                         make_bond(i,newatom,1);
                         for (j=0; j < MAXSS; j++)
                             atom[newatom].substr[j] = atom[i].substr[j];    
                      }
                   }
               }
           }
        }
   }
   return;
} 

void  hcoord()
{
        int i, i1, i3, i4, i5, i6, i7, i8, 
         i9, iadj[10], ii, ii1, ii2, ii3, iii, 
         ij, ijk[3], it, j, j1, j2, j3, j4, j5, j8,  
         j9, k1, k2, l1, na, nhadds, ni;
        float adj, dis, xpnt, ypnt, zpnt;

        getvec();
        for( i = 1; i <= natom; i++ )
        {
           it = atom[i].mmx_type;
           if( it < 300)
           {
              for( j = 0; j < 10; j++ )
                   iadj[j] = 0;

              for( j = 0; j < 3; j++ )
                   ijk[j] = 0;

              nhadds = 0;
              l1 = 0;
              for( j = 0; j < MAXIAT; j++ )
              {
                  if( atom[i].iat[j] != 0 && atom[i].bo[j] !=  9 )
                  {
                      if( atom[atom[i].iat[j]].mmx_type == 5 || atom[atom[i].iat[j]].mmx_type == 20 )
                           nhadds++ ;
                      else
                      {
                           iadj[l1] = atom[i].iat[j];
                           l1++;
                      }
                  }
              }
              if( nhadds != 0 )
              {
/*       **** calculate the coordinates of the center we
 *       are going to move to the x-axis */
                 xpnt = atom[i].x;
                 ypnt = atom[i].y;
                 zpnt = atom[i].z;
                 for( ii = 1; ii <= (natom + 4); ii++ )
                 {
                       atom[ii].x = atom[ii].x - xpnt;
                       atom[ii].y = atom[ii].y - ypnt;
                       atom[ii].z = atom[ii].z - zpnt;
                 }
                 xpnt = 0.F;
                 ypnt = 0.F;
                 zpnt = 0.F;
                 adj = 0.F;
                 for( i1 = 0; i1 < 10; i1++ )
                 {
                     if( iadj[i1] != 0 )
                     {
                         na = iadj[i1];
                         adj += 1.0001F ;
                         xpnt += atom[na].x;
                         ypnt += atom[na].y;
                         zpnt += atom[na].z;
                     }
                 }
                 if( adj > 0.01F )
                 {
                      xpnt /= adj;
                      ypnt /= adj;
                      zpnt /= adj;
                      xaxis( &xpnt, &ypnt, &zpnt );
                      if( xpnt > 0.F )
                      {
                         for( i3 = 1; i3 <= (natom + 4); i3++ )
                                atom[i3].x = -atom[i3].x;
                      }
/*       ****  after moving, check the # of h added to the carbon i */
                      if( nhadds == 1 )
                      {
/*       **** there is only 1 h added to the carbon i */
                         for( i4 = 0; i4 < MAXIAT; i4++ )
                         {
                            if( atom[i].iat[i4] != 0 && atom[i].bo[i4] != 9 )
                            {
                                i5 = atom[i].iat[i4];
                                if( atom[i5].mmx_type == 5 || atom[i5].mmx_type == 20 )
                                {
/*           i5 is the atom to attach */
                                    atom[i5].x = 1.10F;
                                    if( atom[i5].mmx_type == 20 )
                                          atom[i5].x = 0.6F;
                                    atom[i5].y = 0.00F;
                                    atom[i5].z = 0.00F;
                                    if( atom[i].atomnum == 16 )
                                    {
                                        atom[i5].x = 0.F;
                                        atom[i5].y = 0.3796F;
                                        atom[i5].z = 1.043F;
                                    }
                                }
                            }
                         }
                      }
/*       **** there are 2 h added to the carbon i */
                      if( nhadds == 2 )
                      {
                          if( it == 2 || it == 3  || it == 29 || it == 30 || it == 9 || it == 41 || it == 46
                                      || it == 48 || it == 7 || it == 37 || it == 39
                                      || it == 42 || it == 38 || it == 40 || it == 44 || it == 66)
                          {
                               for( i8 = 0; i8 < 10; i8++ )
                               {
                                   if( iadj[i8] != 0 )
                                   {
                                       ii1 = iadj[i8];
                                       for( ii2 = 0; ii2 < MAXIAT; ii2++ )
                                       {
                                           if( atom[ii1].iat[ii2] != 0 && atom[ii1].iat[ii2] != i &&
                                            atom[ii1].bo[ii2] != 9 )
                                           {
                                                  ii3 = atom[ii1].iat[ii2];
                                                  xpnt = atom[ii3].x;
                                                  ypnt = atom[ii3].y;
                                                  zpnt = atom[ii3].z;
                                                  xyplan( &xpnt, &ypnt, &zpnt );
                                                  ij = 0;
                                                  for( i9 = 0; i9 < MAXIAT; i9++ )
                                                  {
                                                      if( atom[i].iat[i9] != 0 && atom[i].bo[i9] != 9 )
                                                      {
                                                           iii = atom[i].iat[i9];
                                                           if( atom[iii].mmx_type == 5 ||  atom[iii].mmx_type == 20 )
                                                           {
                                                                ijk[ij] = iii;
                                                                ij++;
                                                                if( ij == 2 )
                                                                {
              /*           *** =ch2 *** */
                                                                      atom[ijk[0]].x = 0.5582F;
                                                                      atom[ijk[1]].x = 0.5582F;
                                                                      atom[ijk[0]].y = 0.9478F;
                                                                      atom[ijk[1]].y = -0.9478F;
                                                                      atom[ijk[0]].z = 0.00F;
                                                                      atom[ijk[1]].z = 0.00F;
                                                                      if( it == 7 || it == 38 || it ==  42 || it == 39 || it == 66)
                                                                      {
                                                                          atom[ijk[0]].x = 0.30F;
                                                                          atom[ijk[1]].x = 0.30F;
                                                                          atom[ijk[0]].y = 0.52F;
                                                                          atom[ijk[1]].y = -0.52F;
                                                                      }
                                                                      if( atom[i].mmx_type == 37 )
                                                                      {
                                                                         atom[ijk[1]].x = 0.35F;
                                                                         atom[ijk[1]].y = -0.5F;
                                                                      }
                                                                }
                                                           }
                                                      }
                                                  }
                                                  break;
                                           }
                                       }
                                   }
                               }
                          } else
                          {
/*       ** there are 2 atoms attached to carbon i,
 *       one of which is "ni" and move these two atoms to xy-plane */
                              for( i6 = 0; i6 < 10; i6++ )
                              {
                                  if( iadj[i6] != 0 )
                                  {
                                      ni = iadj[i6];
                                      xpnt = atom[ni].x;
                                      ypnt = atom[ni].y;
                                      zpnt = atom[ni].z;
                                      xyplan( &xpnt, &ypnt, &zpnt );
                                  }
                              }
                              ij = 0;
                              for( i7 = 0; i7 < MAXIAT; i7++ )
                              {
                                  if( atom[i].iat[i7] != 0 && atom[i].bo[i7] != 9 )
                                  {
                                      iii = atom[i].iat[i7];
                                      if( atom[iii].mmx_type == 5 || atom[iii].mmx_type == 20 )
                                      {
                                          ijk[ij] = iii;
                                          ij++ ;
                                          if( ij == 2 )
                                          {
                  /*       *** -ch2- type 1 atom *** */
                                              atom[ijk[0]].x = 0.6758F;
                                              atom[ijk[1]].x = 0.6758F;
                                              atom[ijk[0]].y = 0.00F;
                                              atom[ijk[1]].y = 0.00F;
                                              atom[ijk[0]].z = 0.8807F;
                                              atom[ijk[1]].z = -0.8807F;
                                              if( atom[i].mmx_type == 25 ||  atom[i].mmx_type == 46 ||  atom[i].mmx_type == 67)
                                              {
                                                  atom[ijk[0]].x = 0.3796F;
                                                  atom[ijk[1]].x = 0.3796F;
                                                  atom[ijk[0]].y = 1.043F;
                                                  atom[ijk[1]].y = -0.5215F;
                                                  atom[ijk[0]].z = 0.0000F;
                                                  atom[ijk[1]].z = 0.9032F;
                                              }
                                              if( atom[i].mmx_type == 22 ||  atom[i].mmx_type == 30 )
                                              {
                 /*         *** -ch2- type 22 and 3o atoms *** */
                                                   atom[ijk[0]].x = 0.2F;
                                                   atom[ijk[1]].x = 0.2F;
                                                   atom[ijk[0]].y = 0.0F;
                                                   atom[ijk[1]].y = 0.0F;
                                                   atom[ijk[0]].z = 1.0F;
                                                   atom[ijk[1]].z = -1.0F;
                                              }
                 /*       *** -o(lp)2- *** */
                                              if (atom[i].mmx_type == 6 && (atom[iadj[0]].mmx_type == 3 ||
                                                  atom[iadj[0]].mmx_type == 2 || atom[iadj[0]].mmx_type == 30 ||
                                                  atom[iadj[0]].mmx_type == 37 || atom[iadj[0]].mmx_type == 40) )
                                              {
                                                   atom[ijk[0]].x = 0.450F;
                                                   atom[ijk[1]].x = 0.2488F;
                                                   atom[ijk[0]].z = 0.0F;
                                                   atom[ijk[1]].z = 0.0F;
                                                   atom[ijk[0]].y = 0.9500F;
                                                   atom[ijk[1]].y = -0.5460F;
                                              }
                                              else if( atom[i].mmx_type == 6 || atom[i].mmx_type == 15 || atom[i].mmx_type == 
                                                  42 || atom[i].mmx_type == 34 || atom[i].mmx_type == 35 || atom[i].mmx_type == 53 )
                                              {
                                                   atom[ijk[0]].x = 0.2488F;
                                                   atom[ijk[1]].x = 0.2488F;
                                                   atom[ijk[0]].y = 0.0F;
                                                   atom[ijk[1]].y = 0.0F;
                                                   atom[ijk[0]].z = 0.5460F;
                                                   atom[ijk[1]].z = -0.5460F;
                                              }                                                 
                                          }
                                      }
                                  }
                              }
                          }
                      }
           /*  end of nhadds = 2
                                         *         **** there are 3 h added to the carbon i *** */
                      if( nhadds == 3 )
                      {
                         for( j1 = 0; j1 < 10; j1++ )
                         {
                             if( iadj[j1] != 0 )
                             {
                                j2 = iadj[j1];
                                if( atom[j2].atomnum != 1 && atom[j2].mmx_type != 20 )
                                     na = j2;
                                if( atom[j2].mmx_type == 2 || atom[j2].mmx_type == 3 )
                                {
                                    for( j3 = 0; j3 < MAXIAT; j3++ )
                                    {
                                         if( (atom[j2].iat[j3] != 0 && atom[j2].iat[j3] !=  i) &&
                                         atom[j2].bo[j3] !=  9 )
                                         {
                                             j4 = atom[j2].iat[j3];
                                             if( atom[j4].mmx_type == 2 || (atom[j4].mmx_type == 7 && atom[j4].mmx_type ==  38) )
                                             {
                                                 xpnt = atom[j4].x;
                                                 ypnt = atom[j4].y;
                                                 zpnt = atom[j4].z;
                                                 xyplan( &xpnt, &ypnt, &zpnt );
                                                 if( atom[j4].y < 0.F )
                                                 {
                                                     for( j5 = 1; j5 <= (natom + 4); j5++ )
                                                         atom[j5].y = -atom[j5].y;
                                                 }
                                             }
                                         }
                                    }
                                }else
                                {
 /*  end of type 2 and type 3 adjacent* */
                                   for( j8 = 0; j8 < MAXIAT; j8++ )
                                   {
                                       if( (atom[j2].iat[j8] != 0 && atom[j2].iat[j8] != i)
                                        && atom[j2].bo[j8] != 9 )
                                        {
                                            j9 = atom[j2].iat[j8];
//                                            if( atom[j9].mmx_type != 5 &&  atom[j9].mmx_type != 20 )
                                            if(atom[j9].mmx_type != 20 )
                                            {
                                                xpnt = atom[j9].x;
                                                ypnt = atom[j9].y;
                                                zpnt = atom[j9].z;
                                                xyplan( &xpnt, &ypnt, &zpnt );
                                                if( atom[j9].y > 0.F )
                                                {
                                                    for( k1 = 1; k1 <= (natom + 4); k1++ )
                                                          atom[k1].y = -atom[k1].y;
                                                }
                                            }
                                        }
                                   }
                                }
/* do type 4 here */
                                ij = 0;
                                for( k2 = 0; k2 < MAXIAT; k2++ )
                                {
                                    if( !(atom[i].iat[k2] == 0 || atom[i].bo[k2] == 9) )
                                    {
                                        iii = atom[i].iat[k2];
                                        if( !(atom[iii].mmx_type != 5 &&  atom[iii].mmx_type != 20) )
                                        {
                                            ijk[ij] = iii;
                                            ij++;
                                            if( ij == 3 )
                                            {
                                                dis = 0.F;
                                              /*       *** -ch3 ***
                                               *       *** add the hydrogens displaced along x enough
                                               *       to give a normal length c-c bond *** */
                                                atom[ijk[0]].x = 0.3796F + dis;
                                                atom[ijk[1]].x = 0.3796F + dis;
                                                atom[ijk[2]].x = 0.3796F + dis;
                                                atom[ijk[0]].y = 1.043F;
                                                atom[ijk[1]].y = -0.5215F;
                                                atom[ijk[2]].y = -0.5215F;
                                                atom[ijk[0]].z = 0.00F;
                                                atom[ijk[1]].z = 0.9032F;
                                                atom[ijk[2]].z = -0.9032F;
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
        }
        /*     *** reorient the structure in starting orientation and
         *     turn off atom markers *** */
        revec();
        check_methane();
        return;
}


void revec()
{
        int i;
        float xpnt, ypnt, zpnt;
        /*     *** this routine used the unit vectors and origin saved
         *     by 'subroutine getvec' to return the molecule to
         *     its original orientation *** */

        if( natom <= MAXATOM - 4 )
        {
            for( i = 1; i <= (natom + 4); i++ )
            {
                 atom[i].x = atom[i].x - atom[natom + 4].x;
                 atom[i].y = atom[i].y - atom[natom + 4].y;
                 atom[i].z = atom[i].z - atom[natom + 4].z;
            }
            xpnt = atom[natom + 1].x;
            ypnt = atom[natom + 1].y;
            zpnt = atom[natom + 1].z;
            xaxis( &xpnt, &ypnt, &zpnt );
            xpnt = atom[natom + 2].x;
            ypnt = atom[natom + 2].y;
            zpnt = atom[natom + 2].z;
            xyplan( &xpnt, &ypnt, &zpnt );
            if( atom[natom + 1].x < 0. )
            {
               for( i = 1; i <= (natom + 4); i++ )
                     atom[i].x = -atom[i].x;
            }
            if( atom[natom + 2].y < 0. )
            {
               for( i = 1; i <= (natom + 4); i++ )
                     atom[i].y = -atom[i].y;
            }
            if( atom[natom + 3].z < 0. )
            {
               for( i = 1; i <= (natom + 4); i++ )
                     atom[i].z = -atom[i].z;
            }
        }
        return;
} 
/*     ------------------------ */
void  xyplan(float *xr, float *yr, float *zr)
{
        int i;
        float cos3, denom, savex, sin3;

        /*     this subroutine rotates the molecule to place atom (xr,yr,*zr)
         *     in the x,y-plane.  note0 y-coordinate will be +. */

        denom = sqrt( *yr**yr + *zr**zr );
        if( fabs( denom ) >= .00001F )
        {
           sin3 = *zr/denom;
           cos3 = *yr/denom;
           for( i = 1; i <= (natom + 4); i++ )
           {
               savex = atom[i].y;
               atom[i].y = atom[i].y*cos3 + atom[i].z*sin3;
               atom[i].z = atom[i].z*cos3 - savex*sin3;
           }
           savex = *yr;
           *yr = *yr*cos3 + *zr*sin3;
           *zr = *zr*cos3 - savex*sin3;
        }
        return;
}
/*     ------------------------ */
void xaxis(float *xpnt,float *ypnt,float *zpnt)
{
        int i;
        float cos1, cos2, denom, save, sin1, sin2;

        /*     this subroutine rotates the molecule about z-axis and then y-axis to
         *     place point with coordinates (i*xpnt,iypnt,izpnt) along x-axis. */

        denom = sqrt( *xpnt**xpnt + *ypnt**ypnt );
        if( fabs( denom ) < .00001F )
        {
                sin1 = 0.F;
                cos1 = 1.F;
        }else
        {
                sin1 = *ypnt/denom;
                cos1 = *xpnt/denom;
                for( i = 1; i <= (natom + 4); i++ )
                {
                        save = atom[i].x;
                        atom[i].x = atom[i].x*cos1 + atom[i].y*sin1;
                        atom[i].y = atom[i].y*cos1 - save*sin1;
                }
                save = *xpnt;
                *xpnt = *xpnt*cos1 + *ypnt*sin1;
                *ypnt = *ypnt*cos1 - save*sin1;
        }
        denom = sqrt( *xpnt**xpnt + *zpnt**zpnt );
        if( fabs( denom ) >= .00001F )
        {
                sin2 = *zpnt/denom;
                cos2 = *xpnt/denom;
                for( i = 1; i <= (natom + 4); i++ )
                {
                        save = atom[i].x;
                        atom[i].x = atom[i].x*cos2 + atom[i].z*sin2;
                        atom[i].z = atom[i].z*cos2 - save*sin2;
                }
                save = *xpnt;
                *xpnt = *xpnt*cos2 + *zpnt*sin2;
                *zpnt = *zpnt*cos2 - save*sin2;
        }
        return;
}
/* ==============================================  */
void getvec()
{

        /*     *** this routine saves unit vectors and origin in the
         *     n+1 - n+4 slots of x, y, and z *** */

        if(natom > MAXATOM - 4 )
                return;
        atom[natom + 1].x = 1.F;
        atom[natom + 1].y = 0.F;
        atom[natom + 1].z = 0.F;
        atom[natom + 2].x = 0.F;
        atom[natom + 2].y = 1.F;
        atom[natom + 2].z = 0.F;
        atom[natom + 3].x = 0.F;
        atom[natom + 3].y = 0.F;
        atom[natom + 3].z = 1.F;
        atom[natom + 4].x = 0.F;
        atom[natom + 4].y = 0.F;
        atom[natom + 4].z = 0.F;
        return;
} 
/* ==============================================  */
void hdel(int lptest)
{
/*  lptest = 0 to remove H and lone pair  */
/*  lptest = 1 to remove lone pairs only  */
/*  lptest = 2 to remove lp & reset various groups for amber  */
   int i, iatt, iatt1, iatyp, it, it1,ntemp;
        /*     *** delete hydrogens *** */
   ntemp = natom;
   for( i = 1; i <= ntemp; i++ )
   {
       if( atom[i].mmx_type != 0 )
       {
          iatyp = atom[i].mmx_type;
          if( (lptest == 0 && (iatyp == 5 || iatyp == 20 )) || (lptest >= 1 && iatyp == 20) ) 
          {
              iatt = atom[i].iat[0];
              iatt1 = atom[i].iat[1];
              if( iatt != 0 || iatt1 != 0)  
              {
                 it = atom[iatt].mmx_type;
                 it1 = atom[iatt1].mmx_type;
                 if (it >= 11 && it <= 14)   // halogens
                    goto L_10;
                 if( (it < 300) && (it1 < 300))   // attached metals
                    deleteatom(i);
               }
           }
         }
L_10:
         continue;
   }
     reseq();
     generate_bonds();
     return;
} 
/* --------------------------  */
void reseq()
{
    int i, ia, ib, iplus, ixx, j, jincr, k;

        /*     ***this subroutine will repack the atom
         *     ***connection table, making sure to fill in
         *     **all "holes" caused by deleting atoms or bonds.
         *     ***the necessary common blocks involve arrays
         *     of the connection table (contb) */

    for( i = 1; i < MAXATOM; i++ )
    {
        if( atom[i].type != 0 )
        {
           for( ia = 0; ia < (MAXIAT - 1); ia++ )
           {
               if( atom[i].iat[ia] == 0 )
               {
                                        /*            *** for every zero entry in the array iat
                                         *             for atom i, want to move down all remaining
                                         *            entries of iat by one slot; after doing so
                                         *             want to store that zero entry in the last
                                         *            slot--iat *** */
                   for( ib = ia + 1; ib < MAXIAT; ib++ )
                   {
                       atom[i].iat[ib - 1] = atom[i].iat[ib];
                       atom[i].bo[ib - 1] = atom[i].bo[ib];
                   }
                   atom[i].iat[MAXIAT - 1] = 0;
               }
           }
        }
    }

/*     *** now want to renumber atoms to leave
 *     no "holes" caused by the deleted atoms *** */
     jincr = 0;
/*  i points to first zero entry
 *  kincr points to next real entry */
     for( i = 1; i < (MAXATOM - 1); i++ )
     {
        iplus = 0;
        if( atom[i].type == 0 )
        {
            for( jincr = i; jincr < MAXATOM; jincr++ )
            {
                if( atom[jincr].type != 0 )
                {
                    iplus = jincr;
                    break;
                }
            }
            if (iplus == 0)
                 break;
            for( j = 0; j < MAXIAT; j++ )
            {
                atom[i].iat[j] = atom[iplus].iat[j];
                atom[i].bo[j] = atom[iplus].bo[j];
            }
            for ( j=0; j < MAXSSCLASS; j++)
                  atom[i].substr[j] = atom[iplus].substr[j];

            atom[i].x = atom[iplus].x;
            atom[i].y = atom[iplus].y;
            atom[i].z = atom[iplus].z;
            atom[i].type = atom[iplus].type;
            atom[i].tclass = atom[iplus].tclass;
            atom[i].type = atom[iplus].type;
            atom[i].mmx_type = atom[iplus].mmx_type;
            atom[i].mm3_type = atom[iplus].mm3_type;
            atom[i].amber_type = atom[iplus].amber_type;
            atom[i].mmff_type = atom[iplus].mmff_type;
            atom[i].atomnum = atom[iplus].atomnum;
            atom[i].serno = atom[iplus].serno;
            atom[i].molecule = atom[iplus].molecule;
            atom[i].residue = atom[iplus].residue;
            atom[i].biotype = atom[iplus].biotype;
            atom[i].formal_charge = atom[iplus].formal_charge;
            atom[i].atomwt = atom[iplus].atomwt;
            atom[i].use = atom[iplus].use;
            atom[i].color = atom[iplus].color;
            atom[i].chrg_color = atom[iplus].chrg_color;
            strcpy(atom[i].name,atom[iplus].name);
            atom[i].charge = atom[iplus].charge;
            atom[i].flags = atom[iplus].flags;
            atom[i].radius = atom[iplus].radius;
            atom[i].vdw_radius = atom[iplus].vdw_radius;
            atom[i].energy = atom[iplus].energy;
        /*           *** delete the atom whose statistics were just
         *           stored in the i-th atom's connection table *** */
            for( ixx = 0; ixx < MAXIAT; ixx++ )
            {
                atom[iplus].iat[ixx] = 0;
                atom[iplus].bo[ixx] = 0;
            }
            atom[iplus].x = 0.0F;
            atom[iplus].y = 0.0F;
            atom[iplus].z = 0.0F;
            atom[iplus].type = 0;
            atom[iplus].tclass = 0;
            atom[iplus].mmx_type = 0;
            atom[iplus].mm3_type = 0;
            atom[iplus].amber_type = 0;
            atom[iplus].mmff_type = 0;
            atom[iplus].atomnum = 0;
            atom[iplus].serno = 0;
            atom[iplus].molecule = 0;
            atom[iplus].residue = 0;
            atom[iplus].biotype = 0;
            atom[iplus].formal_charge = 0;
            atom[iplus].atomwt = 0.0F;
            atom[iplus].use = 0;             
            atom[iplus].color = 0;
            atom[iplus].chrg_color = 7;
            atom[iplus].charge = 0.0F;
            atom[iplus].flags = 0;
            atom[iplus].radius = 0.0F;
            atom[iplus].energy = 0.0F;
            strcpy(atom[iplus].name,"");
            for ( j=0; j < MAXSSCLASS; j++)
                atom[iplus].substr[j] = 0;

            for( j = 0; j < MAXIAT; j++ )
            {
                if( atom[i].iat[j] != 0 )
                {
                    for( k = 0; k < MAXIAT; k++ )
                    {
                        if( atom[atom[i].iat[j]].iat[k] == iplus )
                             atom[atom[i].iat[j]].iat[k] = i;
                    }
                }
            }
        }
     }
/*     ***now want to do the same resequencing for iat array */

/* recalc number of atoms */
     natom = 0;
     for (i = 1; i < MAXATOM; i++)
     {
        if (atom[i].type != 0)
           natom++;
     }
     return;
}

void check_methane()
{
        int i, j, it, nhadds, i5;
        float xpnt, ypnt, zpnt;

        for( i = 1; i <= natom; i++ )
        {
           it = atom[i].type;
           if( it < 300)
           {
              nhadds = 0;
              for( j = 0; j < MAXIAT; j++ )
              {
                 if( atom[i].iat[j] != 0 && atom[i].bo[j] !=  9 )
                 {
                     if( atom[atom[i].iat[j]].type == 5 || atom[atom[i].iat[j]].type == 20 )
                         nhadds++ ;
                 }
              }
              if (nhadds == 4)
              {
  // move i to origin
                 xpnt = atom[i].x;
                 ypnt = atom[i].y;
                 zpnt = atom[i].z;
                 for (j=1; j <= natom; j++)
                 {
                     atom[j].x -= xpnt;
                     atom[j].y -= ypnt;
                     atom[j].z -= zpnt;
                 }
 // add hydrogens
                 i5 = atom[i].iat[0];
                 atom[i5].x = 1.113F;
                 atom[i5].y = 0.0F;
                 atom[i5].z = 0.0F;

                 i5 = atom[i].iat[1];
                 atom[i5].x = -0.372F;
                 atom[i5].y = -0.524F;
                 atom[i5].z =  0.909F;
        
                 i5 = atom[i].iat[2];
                 atom[i5].x = -0.372F;
                 atom[i5].y = -0.524F;
                 atom[i5].z = -0.909F;

                 i5 = atom[i].iat[3];
                 atom[i5].x = -0.373F;
                 atom[i5].y =  1.049F;
                 atom[i5].z =  0.0F;
 // move back
                 for (j=1; j <= natom; j++)
                 {
                     atom[j].x += xpnt;
                     atom[j].y += ypnt;
                     atom[j].z += zpnt;
                 }
              }
              if (nhadds == 3 && (it == 30  || it == 29))
              {
  // move i to origin
                 xpnt = atom[i].x;
                 ypnt = atom[i].y;
                 zpnt = atom[i].z;
                 for (j=1; j <= natom; j++)
                 {
                     atom[j].x -= xpnt;
                     atom[j].y -= ypnt;
                     atom[j].z -= zpnt;
                 }
 // add hydrogens
                 i5 = atom[i].iat[0];
                 atom[i5].x = 1.113F;
                 atom[i5].y = 0.0F;
                 atom[i5].z = 0.0F;

                 i5 = atom[i].iat[1];
                 atom[i5].x = -0.5582F;
                 atom[i5].y = 0.9478F;
                 atom[i5].z =  0.0F;
        
                 i5 = atom[i].iat[2];
                 atom[i5].x = -0.5582F;
                 atom[i5].y = -0.9478F;
                 atom[i5].z = 0.0F;

 // move back
                 for (j=1; j <= natom; j++)
                 {
                     atom[j].x += xpnt;
                     atom[j].y += ypnt;
                     atom[j].z += zpnt;
                 }
              }
           }
        }
}
