#include "pcwin.h"
#include "pcmod.h"

#include "angles.h"
#include "rings.h"
#include "field.h"
#include "bonds_ff.h"
#include "tsbond.h"
#include "eangle.h"
#include "angk.h"

void kangle(void);
void set_missing_constants(int);        

void kangle()
{
    long int pi_mask, aromatic_mask;
    int i, j, k, izero, ierr, nh, nRc, nligand, jji;
    int ia, ib, ic, ie;
    int iaa, ibb, icc;
    int cl_a, cl_b, cl_c;
    int nbo1, nb1, nb2;
    char pa[4],pb[4],pc[4],pt[10], pt1[10], pt2[10],pt3[10], pt4[10], pt5[10], pt6[10];
    char zero[4];
    float bkan, at;

    pi_mask = (1L << PI_MASK);
    aromatic_mask = (1 << AROMATIC_MASK);

        if( angles.nang != 0 )
        {
           for( i = 0; i < angles.nang; i++ )
           {
               ia = angles.i13[i][0];
               ib = angles.i13[i][1];
               ic = angles.i13[i][2];
               iaa = atom[ia].type;
               ibb = atom[ib].type;
               icc = atom[ic].type;
               cl_a = atom[ia].tclass;
               cl_b = atom[ib].tclass;
               cl_c = atom[ic].tclass;
               
/* get bond index each bond for MMFF94  */
               nb1 = find_bond(ia,ib);
               nb2 = find_bond(ib,ic);
                              
/* special metal cases */
               if( iaa >= 300)
                    cl_a = 300;
                       
               if( ibb >= 300)
                    cl_b = 300;
                        
               if( icc >= 300)
                    cl_c = 300;

/* special MMX cases  */
               if ( field.type == MMX)
               {                                               
                    if( iaa == 40 )
                         iaa = 2;
                       
                    if( ibb == 40 )
                         ibb = 2;
                        
                    if( icc == 40 )
                         icc = 2;
                         
                }
                numeral(iaa,pa,3);
                numeral(ibb,pb,3);
                numeral(icc,pc,3);
                if( iaa <= icc )
                {
                   strcpy(pt,pa);
                   strcat(pt,pb);
                   strcat(pt,pc);
                }
                if( iaa > icc )
                {
                   strcpy(pt,pc);
                   strcat(pt,pb);
                   strcat(pt,pa);
                }
                izero = 0;
                strcpy(zero,"  0");
                strcpy(pt1,pa); strcat(pt1,pb); strcat(pt1,zero);
                strcpy(pt2,pc); strcat(pt2,pb); strcat(pt2,zero);
                strcpy(pt3,zero); strcat(pt3,pb); strcat(pt3,pc);
                strcpy(pt4,zero); strcat(pt4,pb); strcat(pt4,pa);
                strcpy(pt5,zero); strcat(pt5,pb); strcat(pt5,zero);

                numeral(cl_a,pa,3);
                numeral(cl_b,pb,3);
                numeral(cl_c,pc,3);
                if (cl_a < cl_c)
                {
                   strcpy(pt6,pa); strcat(pt6,pb); strcat(pt6,pc);
                }else
                {
                   strcpy(pt6,pc); strcat(pt6,pb); strcat(pt6,pa);
                }
                angles.acon[i] = 0.0;
                angles.anat[i] = 0.0;
                angles.angtype[i] = HARMONIC;

                ierr = FALSE;
/*     check TS constants  */
                if (field.type == MMX)
                {
                    if (is_TS(iaa,ibb,icc))
                    {
                        /* get TS parameters */
                        nRc = is_ring43(ia, ib, ic);
                        if (nRc == TRUE)
                        {
                            ts_four(&ierr, &bkan, &at, pt, pt1, pt3, pt5);
                            if (ierr == TRUE)
                            {
                                angles.anat[i] = at;
                                angles.acon[i] = bkan;
                                angles.index[i] = 0;
                            }
                        } else
                        {
                            transtheta( &ierr, &bkan, &at, pt, pt1, pt3, pt5);
                            if (ierr == TRUE)
                            {
                               angles.anat[i] = at;
                               angles.acon[i] = bkan;
                               angles.index[i] = 0;
                             }
                        }                          
                        if (ierr == TRUE)
                          goto L_10;
                        else
                        {
			  set_missing_constants(TRUE);
                           fprintf(errfile,"TS Angle parameters missing for angle %d-%d-%d, types: %d-%d-%d\n",
                             ia, ib, ic, iaa, ibb, icc);
                        }
                    }
                }
/*     check pi angles */
                if (angkp1.npiang > 0)
                {
                    if (atom[ib].flags & pi_mask)
                    {
                        for (j=0; j < angkp1.npiang; j++)
                        {
                            if ( (strcmp(angkp1.kpa[j],pt) == 0)  || (strcmp(angkp1.kpa[j],pt1) == 0) ||
                                (strcmp(angkp1.kpa[j],pt2) == 0) || (strcmp(angkp1.kpa[j],pt3) == 0) ||
                                (strcmp(angkp1.kpa[j],pt5) == 0) || (strcmp(angkp1.kpa[j],pt5) == 0) )
                            {
                                angles.acon[i] = angkp1.pacon[j];
                                angles.anat[i] = angkp1.panat[j][0];
                                angles.index[i] = 0;
                                ierr = TRUE;
                                break;
                            }
                        }
                        if (ierr == TRUE)
                           goto L_10;
                        else
                        {
                           set_missing_constants(TRUE);
                           fprintf(errfile,"PI Angle parameters missing for angle %d-%d-%d, types: %d-%d-%d\n",
                             ia, ib, ic, iaa, ibb, icc);
                        }
                    }
                }
/*  check delocalized angles in MMFF94 cyclopropanes */
                if (field.type == MMFF94 && angk1.ndel3 > 0)
                {
                    if (isbond(ia,ic) )
                    {
                        nbo1 = bonds_ff.index[nb1] + bonds_ff.index[nb2];
                        if (nbo1 == 1 || nbo1 == 2)
                        { 
                           for(j=0; j < angk1.ndel3; j++)
                           {
                               if (strcmp(angk1.kdel3[j],pt) == 0)
                               {
                                   angles.acon[i] = angk1.condel3[j];
                                   angles.anat[i] = angk1.angdel3[j][0];
                                   if (nbo1 == 1)
                                      angles.index[i] = 5;
                                   else if (nbo1 == 2)
                                      angles.index[i] = 6;
                                   ierr = TRUE;
                                   break;
                               }
                           }
                        }
                    
                        if (ierr == TRUE)
                          goto L_10;
                    }
                } 
/*     check three membered rings */
                if ( (rings.nring3 > 0) && (angk1.nang3 > 0) )
                {
                    if (isbond(ia,ic))
                    {
/*                         search three membered ring parameters */
                      for (j = 0; j < angk1.nang3; j++)
                      {
                        if ( (strcmp(angk1.ktype3[j],pt) == 0) || (strcmp(angk1.ktype3[j],pt1) == 0) ||
                            (strcmp(angk1.ktype3[j],pt2) == 0) || (strcmp(angk1.ktype3[j],pt3) == 0) ||           
                            (strcmp(angk1.ktype3[j],pt4) == 0) || (strcmp(angk1.ktype3[j],pt5) == 0))
                        {
			  if ( fabs(angk1.ang3[j][1]) < 0.01 && fabs(angk1.ang3[j][2]) < 0.01)
                           {
                             angles.acon[i] = angk1.con3[j];
                             angles.anat[i] = angk1.ang3[j][0];
                             angles.index[i] = 3;
                             ierr = TRUE;
                             break;
                           } else
                           {
                              nh = 0;
                              for (k=0; k < MAXIAT; k++)
                              {
                                ie = atom[ib].iat[k];
                                if ( (atom[ie].type == 5 || atom[ie].type == 36) && (ie != ia) && (ie != ic))
                                  nh++;
                              }
                              angles.acon[i] = angk1.con3[j];
                              angles.anat[i] = angk1.ang3[j][nh];
                              angles.index[i] = 3;
                              ierr = TRUE;
                              break;
                           }
                        }
                     }
                     if (ierr == TRUE)
                       goto L_10;
                     else
                     {
                            set_missing_constants(TRUE);
                           fprintf(errfile,"Angle3 parameters missing for angle %d-%d-%d, types: %d-%d-%d\n",
                             ia, ib, ic, iaa, ibb, icc);
                     }
                    }
                }
/*  check delocalized angles in MMFF94 cyclobutanes */
                if (field.type == MMFF94 && angk1.ndel4 > 0)
                {
                    nRc = is_ring43(ia, ib, ic);
                    if (nRc == TRUE)
                    {
                        nbo1 = bonds_ff.index[nb1] + bonds_ff.index[nb2];
                        if (nbo1 == 1 || nbo1 == 2)
                        {
                           for(j=0; j < angk1.ndel4; j++)
                           {
                               if (strcmp(angk1.kdel4[j],pt) == 0)
                               {
                                   angles.acon[i] = angk1.condel4[j];
                                   angles.anat[i] = angk1.angdel4[j][0];
                                   if (nbo1 == 1)
                                      angles.index[i] = 7;
                                   else if (nbo1 == 2)
                                      angles.index[i] = 8;
                                   ierr = TRUE;
                                   break;
                               }
                            }
                        }
                        if (ierr == TRUE)
                          goto L_10;
                    }
                } 
/*     check four memebered rings */
                if ( (rings.nring4 > 0) && (angk1.nang4 > 0) )
                {
                    nRc = is_ring43(ia, ib, ic);
                    if (nRc == TRUE)
                    {
 /*                        search four membered ring parameters */
                      for (j = 0; j < angk1.nang4; j++)
                      {
                        if ( (strcmp(angk1.ktype4[j],pt) == 0) || (strcmp(angk1.ktype4[j],pt1) == 0) ||
                            (strcmp(angk1.ktype4[j],pt2) == 0) || (strcmp(angk1.ktype4[j],pt3) == 0) ||           
                            (strcmp(angk1.ktype4[j],pt4) == 0) || (strcmp(angk1.ktype4[j],pt5) == 0))
                        {
			  if ( fabs(angk1.ang4[j][1]) < 0.01 && fabs(angk1.ang4[j][2]) < 0.01)
                           {
                             angles.acon[i] = angk1.con4[j];
                             angles.anat[i] = angk1.ang4[j][0];
                             angles.index[i] = 4;
                             ierr = TRUE;
                             break;
                           } else
                           {
                              nh = 0;
                              for (k=0; k < MAXIAT; k++)
                              {
                                ie = atom[ib].iat[k];
                                if ( (atom[ie].type == 5 || atom[ie].type == 36) && (ie != ia) && (ie != ic))
                                  nh++;
                              }
                              angles.acon[i] = angk1.con4[j];
                              angles.anat[i] = angk1.ang4[j][nh];
                              angles.index[i] = 4;
                              ierr = TRUE;
                              break;
                           }
                        }
                     }
                     if (ierr == TRUE)
                        goto L_10;
                     else
                     {
                           set_missing_constants(TRUE);                         
                           fprintf(errfile,"Angle4 parameters missing for angle %d-%d-%d, types: %d-%d-%d\n",
                             ia, ib, ic, iaa, ibb, icc);
                     }
                    }
                }                    
/*     check five membered rings */
                if ( (rings.nring5 > 0) && (angk1.nang5 > 0) )
                {
                    nRc = is_ring53(ia,ib, ic);
                    if (nRc == TRUE)
                    {
 /*                        search five memebered ring parameters */
                      for (j = 0; j < angk1.nang5; j++)
                      {
                        if ( (strcmp(angk1.ktype5[j],pt) == 0) || (strcmp(angk1.ktype5[j],pt1) == 0) ||
                            (strcmp(angk1.ktype5[j],pt2) == 0) || (strcmp(angk1.ktype5[j],pt3) == 0) ||           
                            (strcmp(angk1.ktype5[j],pt4) == 0) || (strcmp(angk1.ktype5[j],pt5) == 0) )
                        {
			  if ( fabs(angk1.ang5[j][1]) < 0.01 && fabs(angk1.ang5[j][2]) < 0.01)
                           {
                             angles.acon[i] = angk1.con5[j];
                             angles.anat[i] = angk1.ang5[j][0];
                             angles.index[i] = 0;
                             ierr = TRUE;
                             break;
                           } else
                           {
                              nh = 0;
                              for (k=0; k < MAXIAT; k++)
                              {
                                ie = atom[ib].iat[k];
                                if ( (atom[ie].type == 5 || atom[ie].type == 36) && (ie != ia) && (ie != ic))
                                  nh++;
                              }
                              angles.acon[i] = angk1.con5[j];
                              angles.anat[i] = angk1.ang5[j][nh];
                              angles.index[i] = 0;
                              ierr = TRUE;
                              break;
                           }
                        }
                     }
                     if (ierr == TRUE)
                        goto L_10;
                     // don't fail on missing 5 ring parameters
                   }
                }
/*  check delocalized angles in MMFF94  */
                if (field.type == MMFF94 && angk1.ndel > 0)
                {
                   if ( bonds_ff.index[nb1] == 1 || bonds_ff.index[nb2] == 1 )
                   {
                      nbo1 = bonds_ff.index[nb1] + bonds_ff.index[nb2];
                      for(j=0; j < angk1.ndel; j++)
                      {
                         if (strcmp(angk1.kdel[j],pt) == 0)
                         {
                             angles.acon[i] = angk1.condel[j];
                             angles.anat[i] = angk1.angdel[j][0];
                             if (nbo1 == 1)
                                angles.index[i] = 1;
                             else if (nbo1 == 2)
                                angles.index[i] = 2;
                             ierr = TRUE;
                             break;
                         }
                      }
                      if (ierr == TRUE)
                        goto L_10;
                   }
                } 
// look for fourier angles
                if (angf.nfang > 0)
                {
                    for (j=0; j < angf.nfang; j++)
                    {
                        if (strcmp(angf.kftype[j],pt) == 0)
                        {
                            angles.fcon[i] = angf.fcon[j];
                            angles.c0[i] = angf.fc0[j];
                            angles.c1[i] = angf.fc1[j];
                            angles.c2[i] = angf.fc2[j];
                            angles.c3[i] = angf.fc3[j];
                            angles.c4[i] = angf.fc4[j];
                            angles.c5[i] = angf.fc5[j];
                            angles.c6[i] = angf.fc6[j];
                            angles.angtype[i] = FOURIER;
                            ierr = TRUE;
                            break;
                        }
                    }
                    if (ierr == TRUE)
                       goto L_10;
                }
                        
/*     check regular parameters  for specific angles*/
                for (j = 0; j < angk1.nang; j++)
                {
                    if ( (strcmp(angk1.ktype[j],pt) == 0) )
                    {
		      if ( fabs(angk1.ang[j][1]) < 0.01 && fabs(angk1.ang[j][2]) < 0.01)
                        {
                          angles.acon[i] = angk1.con[j];
                          angles.anat[i] = angk1.ang[j][0];
                          angles.index[i] = 0;
                          ierr = TRUE;
                          break;
                        } else
                        {
                            nh = 0;
                            for (k=0; k < MAXIAT; k++)
                            {
                                ie = atom[ib].iat[k];
                                if ( (atom[ie].atomnum == 1) && (ie != ia) && (ie != ic))
                                  nh++;
                            }
                            angles.acon[i] = angk1.con[j];
                            angles.anat[i] = angk1.ang[j][nh];
                            angles.index[i] = 0;
                            ierr = TRUE;
                            break;
                        }
                    }
                }
                if (ierr == TRUE)
                   goto L_10;
// did not find any specific look for generalized                         
                for (j = 0; j < angk1.nang; j++)
                {
                    if ( (strcmp(angk1.ktype[j],pt) == 0) || (strcmp(angk1.ktype[j],pt1) == 0)  ||
                         (strcmp(angk1.ktype[j],pt2) == 0) || (strcmp(angk1.ktype[j],pt3) == 0) ||           
                         (strcmp(angk1.ktype[j],pt4) == 0) || (strcmp(angk1.ktype[j],pt5) == 0) ||
                         (strcmp(angk1.ktype[j],pt6) == 0))
                    {
		      if ( fabs(angk1.ang[j][1]) < 0.01 && fabs(angk1.ang[j][2]) < 0.01)
                        {
                          angles.acon[i] = angk1.con[j];
                          angles.anat[i] = angk1.ang[j][0];
                          angles.index[i] = 0;
                          ierr = TRUE;
                          break;
                        } else
                        {
                            nh = 0;
                            for (k=0; k < MAXIAT; k++)
                            {
                                ie = atom[ib].iat[k];
                                if ( (atom[ie].atomnum == 1) && (ie != ia) && (ie != ic))
                                  nh++;
                            }
                            angles.acon[i] = angk1.con[j];
                            angles.anat[i] = angk1.ang[j][nh];
                            angles.index[i] = 0;
                            ierr = TRUE;
                            break;
                        }
                    }
                }                         
                if (ierr != TRUE)
                {
                      set_missing_constants(TRUE);
                      fprintf(errfile,"Specific Angle parameters missing for angle %d-%d-%d, types: %d-%d-%d\n",
                         ia, ib, ic, iaa, ibb, icc);
                }
L_10:
                if (field.type == MMX && ierr != TRUE)
                {
                    if (ibb >= 300)
                    {
                        nligand = 0;
                        jji = 0;
                        for ( k=0; k < MAXIAT; k++)
                        {
                            if (atom[ib].iat[k] != 0 &&  atom[ib].bo[k] != 9)
                              jji++;
                            if (atom[ib].bo[k] == 9)
                              nligand++;
                        }
                        if (jji == 2 && nligand == 0)
                        {
                            angles.acon[i] = .1;
                            angles.anat[i] = 180.;
                        }
                        if ((jji == 3 && nligand == 0) || (jji == 2 && nligand == 1) )
                        {
                            angles.acon[i] = .1;
                            angles.anat[i] = 120.;
                        }    
                        if ((jji == 4 && nligand == 0) || (jji == 3 && nligand == 1) || (jji == 2 && nligand == 2))
                        {
                            if (jji == 4 && (atom[ib].type >= 300 && (atom[ib].flags & (1L << SATMET_MASK) && 
                                atom[ib].flags & (1L << GT18e_MASK) && atom[ib].flags & (1L <<  SQPLAN_MASK) &&
                                !(atom[ib].flags & (1L << LOWSPIN_MASK)))) )
                            { 
                              angles.acon[i] = .0001;
                              angles.anat[i] = 90.;
                            } else
                            {
                              angles.acon[i] = .1;
                              angles.anat[i] = 109.5;
                            }
                        }    
                            
                    }
                }
           }
        }
}

/* --------------------------------- */
long kang(int i,int j,int k)
{
        long int kang_v;
        static long i10000 = 10000;
        static long i100 = 100;

        kang_v = i10000*i + i100*j + k;
        return( kang_v );
} 

/* -------------------- */
long ipack(int i,int j,int k,int l)
{
        long int ipack_v;
        static long itbig = 200000000;
        static long itmid = 2000000;
        static long ithou = 1000;

        ipack_v = itbig*i + itmid*j + ithou*k + l;
        return( ipack_v );
}
/* ---------------------------------------------- */
int is_TS(int ia, int ib, int ic)
{

    if ( ia == 49 || ia == 50 || ia == 43 || ia == 45 || ia == 53 ||
         ia == 52 || ia == 51 || ia == 54 || ia == 55 )
         return TRUE;
         
    if ( ib == 49 || ib == 50 || ib == 43 || ib == 45 || ib == 53 ||
         ib == 52 || ib == 51 || ib == 54 || ib == 55 )
         return TRUE;

    if ( ic == 49 || ic == 50 || ic == 43 || ic == 45 || ic == 53 ||
         ic == 52 || ic == 51 || ic == 54 || ic == 55 )
         return TRUE;
         
    return FALSE;
}
/* ---------------------------------------------- */
void transtheta(int *ierr,float *bkan,float *at,char *pt,char *pt1,char *pt3,char *pt5)
{
        int  iadj, iii, itajk, jjj, kkk, nbond;

        /* look for transition state parameters first */
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
                                                if( atom[jjj].iat[kkk] != 0 &&  atom[jjj].bo[kkk] != 9 )
                                                {
                                                        itajk = atom[atom[jjj].iat[kkk]].type;
                                                        if( itajk == 29 || itajk == 30 || itajk == 50 || 
                                                            itajk == 53 || itajk == 42 )
                                                                goto L_291;
                                                }
                                        }
                                }
                        }
                        continue;
L_291:
                        iadj = itajk;
                }
        }
        if( (strcmp(pt," 45 49 50") == 0) && ts_bondorder.fbnd[4] <= 0.001 )
        {
                *bkan = .05 + .3*ts_bondorder.fbnd[3];
                *at = 90. + 19.5*ts_bondorder.fbnd[3];
                *ierr = TRUE;
        }
        else if( strcmp(pt,"  1 45 49" ) == 0)
        {
                *bkan = .35*ts_bondorder.fbnd[2];
                *at = 109.5;
                *ierr = TRUE;
        }
        else if( strcmp(pt,"  3 49 49") == 0 )
        {
                if( nbond < 2 && iadj == 30 )
                {
                        *bkan = .01 + .9*ts_bondorder.fbnd[0];
                        *at = 90. - ts_bondorder.fbnd[0]*30.;
                        *ierr = TRUE;
                }
        }
        else if( strcmp(pt1," 49 49  0") == 0 || strcmp(pt3,"  0 49 49") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[0];
                *at = 90. + ts_bondorder.fbnd[0]*19.5;
                *ierr = TRUE;
        }
        else if( strcmp(pt1," 50 50  0") == 0 || strcmp(pt3,"  0 50 50") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[1];
                *at = 90. + ts_bondorder.fbnd[1]*19.5;
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 49 45") == 0 || strcmp(pt1," 45 49  0") == 0 )
        {
                *bkan = .4;
                *at = 90. + ts_bondorder.fbnd[3]*19.5;
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 43 45") == 0 || strcmp(pt1," 45 43  0") == 0 )
        {
                *bkan = 0.001 + .35*(1. - ts_bondorder.fbnd[4]);
                *at = 120. - (ts_bondorder.fbnd[4]*30.);
                *ierr = TRUE;
        }
        else if( strcmp(pt," 49  8 50") == 0 )
        {
                *bkan = .4;
                *at = 120. - 10.5*ts_bondorder.fbnd[0];
                *ierr = TRUE;
        }
        else if( strcmp(pt," 49  8 53") == 0 )
        {
                *bkan = .4;
                *at = 120. - 10.5*ts_bondorder.fbnd[0];
                *ierr = TRUE;
                /* new from 7-1-86 */
        }
        else if( strcmp(pt," 20 53 52") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[9];
                *at = 90. + 19.5*ts_bondorder.fbnd[9];
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 53 52") == 0 || strcmp(pt1," 52 53  0") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[9];
                *at = 180. - 70.5*ts_bondorder.fbnd[9];
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 53 50") == 0 || strcmp(pt1," 50 53  0") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[6];
                *at = 90. + 19.5*ts_bondorder.fbnd[6];
                *ierr = TRUE;
        }
        else if( strcmp(pt5,"  0 53  0") == 0 )
        {
                *bkan = .4;
                *at = 120. - 10.5*ts_bondorder.fbnd[6];
                *ierr = TRUE;
        }
        else if( ts_bondorder.fbnd[1] > 0. && strcmp(pt," 20  8 49") == 0 )
        {
                *bkan = .01 + .4*(1 - ts_bondorder.fbnd[0] + ts_bondorder.fbnd[1]);
                *at = 90. + 19.5*(1 - ts_bondorder.fbnd[0] + ts_bondorder.fbnd[1]);
                *ierr = TRUE;
        }
        else if( ts_bondorder.fbnd[6] > 0. && strcmp(pt," 20  8 49") == 0 )
        {
                *bkan = .01 + .4*(1 - ts_bondorder.fbnd[0] + ts_bondorder.fbnd[6]);
                *at = 90. + 19.5*(1 - ts_bondorder.fbnd[0] + ts_bondorder.fbnd[6]);
                *ierr = TRUE;
        }
        else if( strcmp(pt," 20  8 50") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[1]*ts_bondorder.fbnd[0];
                *at = 90. + 19.5*ts_bondorder.fbnd[1]*ts_bondorder.fbnd[0];
                *ierr = TRUE;
        }
        else if( strcmp(pt," 20  8 53") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[6]*ts_bondorder.fbnd[0];
                *at = 90. + 19.5*ts_bondorder.fbnd[6]*ts_bondorder.fbnd[0];
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0  8 49") == 0 || strcmp(pt1," 49  8  0") == 0 )
        {
                *bkan = .4;
                *at = 120. - 10.5*ts_bondorder.fbnd[0];
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0  8 50") == 0 || strcmp(pt1," 50  8  0") == 0 )
        {
                *bkan = .4;
                *at = 120. - 10.5*ts_bondorder.fbnd[1];
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0  8 53") == 0 || strcmp(pt1," 53  8  0") == 0 )
        {
                *bkan = .4;
                *at = 120. - 10.5*ts_bondorder.fbnd[6];
                *ierr = TRUE;
        }
        else if( strcmp(pt," 51 37 53") == 0 )
        {
                *bkan = .1 + .3*ts_bondorder.fbnd[5];
                *at = 180. - 60.*ts_bondorder.fbnd[5];
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 51 49") == 0 || strcmp(pt1," 49 51  0") == 0 )
        {
                *bkan = .15 + .3*ts_bondorder.fbnd[5];
                *at = 90. + 30.*ts_bondorder.fbnd[5];
                *ierr = TRUE;
        }
        else if( strcmp(pt5,"  0 51  0") == 0 )
        {
                *bkan = .15 + .3*ts_bondorder.fbnd[5];
                *at = 180. - 60.*ts_bondorder.fbnd[5];
                *ierr = TRUE;
        }
        else if( strcmp(pt," 49 52 54") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[7]*ts_bondorder.fbnd[8];
                *at = 180.;
                *ierr = TRUE;
        }
          else if (strcmp(pt," 43 50 49") == 0)
          {
             *bkan=0.05+.35*ts_bondorder.fbnd[4];
             *at=60.+(ts_bondorder.fbnd[4]*ts_bondorder.fbnd[4])*25.0;
//             *at=60.+(ts_bondorder.fbnd[4]*ts_bondorder.fbnd[4])*49.5;
             *ierr = TRUE;
          }
        else if( strcmp(pt," 49 52 53") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[7]*ts_bondorder.fbnd[9];
                *at = 180.;
                *ierr = TRUE;
        }
        else if( strcmp(pt," 53 52 54") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[8]*ts_bondorder.fbnd[9];
                *at = 180.;
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 52 49") == 0 || strcmp(pt1," 49 52  0") == 0 )
        {
                if( ts_bondorder.fbnd[9] > 0.0001 )
                {
                        *bkan = .3*(ts_bondorder.fbnd[7] + ts_bondorder.fbnd[9]);
                        *at = 90. + 20.*(ts_bondorder.fbnd[7] - ts_bondorder.fbnd[9]);
                        *ierr = TRUE;
                }
                else
                {
                        *at = 90. + 20.*(ts_bondorder.fbnd[7] - ts_bondorder.fbnd[8]);
                        *bkan = .3*(ts_bondorder.fbnd[7] + ts_bondorder.fbnd[8]);
                        *ierr = TRUE;
                }
        }
        else if( strcmp(pt3,"  0 52 54") == 0 || strcmp(pt1," 54 52  0") == 0 )
        {
                if( ts_bondorder.fbnd[9] > 0.0001 )
                {
                        *bkan = .3*(ts_bondorder.fbnd[8] + ts_bondorder.fbnd[9]);
                        *at = 90. + 20.*(ts_bondorder.fbnd[8] - ts_bondorder.fbnd[9]);
                        *ierr = TRUE;
                }
                else
                {
                        *at = 90. + 20.*(ts_bondorder.fbnd[8] - ts_bondorder.fbnd[7]);
                        *bkan = .3*(ts_bondorder.fbnd[7] + ts_bondorder.fbnd[8]);
                        *ierr = TRUE;
                }
        }
        else if( strcmp(pt3,"  0 52 53") == 0 || strcmp(pt1," 53 52  0") == 0 )
        {
                if( ts_bondorder.fbnd[8] > 0.0001 )
                {
                        *at = 90. + 20.*(ts_bondorder.fbnd[9] - ts_bondorder.fbnd[8]);
                        *bkan = .3*(ts_bondorder.fbnd[8] + ts_bondorder.fbnd[9]);
                        *ierr = TRUE;
                }
                else
                {
                        *bkan = .3*(ts_bondorder.fbnd[7] + ts_bondorder.fbnd[9]);
                        *at = 90. + 20.*(ts_bondorder.fbnd[9] - ts_bondorder.fbnd[7]);
                        *ierr = TRUE;
                }
        }
        else if( strcmp(pt3,"  0 49 52") == 0 || strcmp(pt1," 52 49  0") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[7];
                *at = 90. + 19.5*ts_bondorder.fbnd[7];
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 53 52") == 0 || strcmp(pt1," 52 53  0") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[9];
                *at = 90. + 19.5*ts_bondorder.fbnd[9];
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 49 51") == 0 || strcmp(pt1," 51 49  0") == 0 )
        {
                *bkan = .3 + ts_bondorder.fbnd[5]*.15;
                *at = 90. + ts_bondorder.fbnd[5]*19.5;
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 50 53") == 0 || strcmp(pt1," 53 5  00") == 0 )
        {
                *bkan = .3 + ts_bondorder.fbnd[6]*.15;
                *at = 90. + ts_bondorder.fbnd[6]*19.5;
                *ierr = TRUE;
        }
        else if ( strcmp(pt, " 50 43 53") == 0) 
        {
                *bkan = .2;
                *at = 110. + 20.*fabs(ts_bondorder.fbnd[4]*ts_bondorder.fbnd[11]);
                *ierr = TRUE;
        }
        else if( (strcmp(pt3, " 0 43 53") == 0 || strcmp (pt1, " 50 43  0") == 0) && ts_bondorder.fbnd[4] !=0.)  
        {
                *bkan = .4 *ts_bondorder.fbnd[4];
                *at = 90. + 30.*ts_bondorder.fbnd[4];
                *ierr = TRUE;
        }
        else if( (strcmp(pt3, "  0 43 53") == 0  || strcmp(pt1, " 53 43  0") == 0) && ts_bondorder.fbnd[11] !=0.)
        {
                *bkan = .6 *ts_bondorder.fbnd[11];
                *at = 90. + 30.*ts_bondorder.fbnd[11];
                *ierr = TRUE;
        }
        else if ( strcmp (pt5, " 0 43  0") == 0 )
        {
                *bkan = .4;
                *at = 120.;
                *ierr = TRUE;
        }
        else if( strcmp(pt5,"  0 50  0") == 0 )
        {
                *bkan = .40;
                *at = 115.;
                if( ts_bondorder.fbnd[1] > 0 )
                        *at = 120. - ts_bondorder.fbnd[1]*10.5;
                if( ts_bondorder.fbnd[4] > 0 )
                        *at = 120. - ts_bondorder.fbnd[4]*10.5;
                if( ts_bondorder.fbnd[6] > 0 )
                        *at = 120. - ts_bondorder.fbnd[6]*10.5;
                *ierr = TRUE;
        }
        else if( strcmp(pt5,"  0 49  0") == 0 )
        {
                *bkan = .40;
                *at = 115.;
                if( ts_bondorder.fbnd[0] > 0. )
                        *at = 120. - ts_bondorder.fbnd[0]*10.5;
                *ierr = TRUE;
        }
        else if( strcmp(pt5,"  0 29  0") == 0 )
        {
                *bkan = .4;
                *at = 120.;
                *ierr = TRUE;
        }
        else if( strcmp(pt5,"  0 37  0") == 0 )
        {
                *bkan = .4;
                *at = 120.;
                *ierr = TRUE;
        }
        else if( strcmp(pt5,"  0 52  0") == 0 )
        {
                *bkan = .4;
                *at = 120.;
                if( ts_bondorder.fbnd[7] > 0. )
                        *at = 120. - 10.5*(fabs( ts_bondorder.fbnd[8] - ts_bondorder.fbnd[8]) );
                if( ts_bondorder.fbnd[9] > 0. )
                        *at = 120. - 10.5*(fabs( ts_bondorder.fbnd[9] - ts_bondorder.fbnd[9]) );
                /* new 3/28/88 */
                *ierr = TRUE;
        }
        else if( strcmp(pt3,"  0 55 50") == 0 || strcmp(pt1," 55 50  0") == 0 )
        {
                *bkan = .01 + .4*ts_bondorder.fbnd[10];
                *at = 90. + 19.5*ts_bondorder.fbnd[10];
                *ierr = TRUE;
        }
        else if( strcmp(pt5,"  0 55  0") == 0 )
        {
                *bkan = .4;
                *at = 120. - 10.5*ts_bondorder.fbnd[10];
                *ierr = TRUE;
        }
        return;
}
/*    ------------------------------------------  */
void ts_four(int *ierr, float *bkan, float *at, char *pt, char *pt1,char *pt3,char *pt5)
{
         if (strcmp(pt," 45 43 50") == 0)
         {
            *bkan = .3501-.35*ts_bondorder.fbnd[4];
            *at = 120.-ts_bondorder.fbnd[4]*19.5;
//            *at = 109.-ts_bondorder.fbnd[4]*19.5;
            *ierr = TRUE;
          }else if (strcmp(pt," 43 50 49") == 0)
          {
             *bkan=0.05+.35*ts_bondorder.fbnd[4];
             *at=60.+(ts_bondorder.fbnd[4]*ts_bondorder.fbnd[4])*49.5;
             *ierr = TRUE;
          } else if (strcmp(pt," 43 45 49") == 0)
          {
             *bkan=0.05;
             *at=77.;
//             *at=60.;
             *ierr = TRUE;
         }else if (strcmp(pt," 45 49 50") == 0)
          {
             if(ts_bondorder.fbnd[4] > 0.001) 
             {
                 *bkan=.05+.3*ts_bondorder.fbnd[4];
                 *at=104.0;
//                 *at=109.5;
                 *ierr = TRUE;
             }
          }else if (strcmp(pt5,"  0 49  0") == 0)
          {
             *bkan=.35;
             *at=120.-ts_bondorder.fbnd[0]*10.5;
             *ierr = TRUE;
          }else if (strcmp(pt5,"  0 50  0") == 0)
          {
             *bkan=.35;
             *at=120.-ts_bondorder.fbnd[1]*10.5;
             *ierr = TRUE;
          }else if ( strcmp(pt1," 45  1  0") == 0 || strcmp(pt3,"  0  1 45") == 0 )
          {
             *bkan=.4+3.*(1.-ts_bondorder.fbnd[2]);
             *at=90.+ts_bondorder.fbnd[2]*19.5;
             *ierr = TRUE;
          }
}
            
