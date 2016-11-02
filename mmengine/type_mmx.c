#include "pcwin.h"
#include "pcmod.h"
#include "field.h"
#include "utility.h"
#include "type_mmx.h"

void type_mmx()
{
   int i, j, ij, jji, jjk, jjbo, jj_bo, iatype,nh;
   int k, l, jatm, katm, latm, jbo,kbo,lbo, ismetal,matm;
   int nnit,nnh,nc,no,nox,nitatom;
   int ktype,jtype, jjik;
   int  adjn, nplus, adjn1, noxide, icurr;
   int mmxtype,mm3type,mmfftype,ambertype,oplstype, ia, ib, ia1, ib1;
   int ndouble, ntriple;
   int icycl3, icycl4, icycl5, icycl6;
   int full_ring, non_pi;
   int array[7];
   int nit[6];
   long int aromatic_mask, mask6, type_mask,pi_mask;

     
   aromatic_mask = (1 << AROMATIC_MASK);
   pi_mask = (1L << PI_MASK);
   mask6 = (1L << RING6);
   type_mask = (1L << NO_RETYPE);
   ambertype = 0;
   oplstype = 0;
   ia = ib = jatm = katm = latm = matm = 0;

   for (i=1; i <= natom; i++)
   {
      mmxtype = atom[i].mmx_type;
      mm3type = atom[i].mm3_type;
      mmfftype = atom[i].mmff_type;
      if (atom[i].atomnum != 0 && !(atom[i].flags & type_mask) )
      {
          if (mmxtype == 100 || mmxtype == 101 || mmxtype == 102 || mmxtype == 103 ) // user defined type 
          {
              goto L_20;
          }
          if (atom[i].atomnum == 1) // hydrogens
          {
              if (atom[i].mmx_type == 60)  // dummy atoms
              {
                  atom[i].flags |= (1L << DUMMY); // mark dummy atom
                  mmxtype = 60;
                  mm3type = 200;
                  mmfftype = 5;
                  ambertype = 33;
                  oplstype = 29;
                  goto L_10;
              }
              if (atom[i].mmx_type == 45)  // TS H 
              {
                  mmxtype = 45;
                  mm3type = 0;
                  mmfftype = 0;
                  ambertype = 33;
                 goto L_10;
              }
              if (atom[i].mmx_type == 36) // Deuterium
              {
                  mmxtype = 36;
                  mm3type = 36;
                  mmfftype = 5;
                  ambertype = 33;
                   goto L_10;
              }                
              jji = 0;
              for (j=0; j < MAXIAT; j++)
              {
                  if (atom[i].iat[j] != 0 && (atom[i].bo[j] == 1 || atom[i].bo[j] == 9))
                     jji++;
              }
              if (jji == 2 ) // hydrogen bonded to two atoms or bonded to one and coordinated to another
              {
                  mmxtype = 70;
                  goto L_10;
              }
              if (atom[atom[i].iat[0]].atomnum == 6)   // carbon
              {
                  mmxtype = 5;
                  mm3type = 5;
                  mmfftype = 5;
                  ambertype = 34;
                  oplstype = 29;
                  jatm = atom[i].iat[0];
                  for (j=0; j < 3; j++)   // Acetylene
                  {
                      if (atom[jatm].bo[j] == 3)
                      {
                        mm3type = 124;
                        oplstype = 33;
                        goto L_10;
                      }
                  }
                  for (j=0; j < MAXIAT; j++)  // ammonium
                  {
                      if (atom[jatm].iat[j] != 0 && atom[atom[jatm].iat[j]].atomnum == 7)
                      {
                          jjik = 0;
                          katm = atom[jatm].iat[j];
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[katm].iat[j] != 0)
                                jjik++;
                          }
                          if (jjik == 4) // H-C-N+
                          {
                              ambertype = 38;
                              goto L_10;
                          }
                      }
                  }
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[jatm].bo[j] == 2 && !(atom[jatm].flags && aromatic_mask) )
                      {
                          if (atom[atom[jatm].iat[j]].atomnum == 6)
                             oplstype = 30;
                          else if (atom[atom[jatm].iat[j]].atomnum == 7)
                             oplstype = 31;
                          goto L_10;
                      }
                      if (atom[jatm].bo[j] == 2 && atom[atom[jatm].iat[j]].atomnum == 8 )
                      {
                          oplstype = 32;
                          goto L_10;
                      }
                  }
                  jji = 0;
                  nnit = 0;
                  nitatom = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                     if (atom[jatm].iat[j] != 0 && atom[atom[jatm].iat[j]].bo[j] != 9 && atom[atom[jatm].iat[j]].mmx_type != 20)
                      jji++;
                     if (atom[atom[jatm].iat[j]].atomnum == 7 || atom[atom[jatm].iat[j]].atomnum == 8 ||
                        atom[atom[jatm].iat[j]].atomnum == 16)
                     {
                         nnit++;
                         nitatom = atom[jatm].iat[j];
                     }
                  }
                  if ( jji == 3 && nnit == 2)
                  {
                     ambertype = 40; // H5
                     oplstype = 43;
                  } else if ( jji == 4 && nnit == 2)
                  {
                     ambertype = 36; // H2
                     oplstype = 29;         
                  } else if ( jji == 3 && nnit == 1)
                  {
                     ambertype = 39; // H4
                     oplstype = 42;
                     jjk = 0;
                     for (j=0; j < MAXIAT; j++)
                     {
                         if (atom[nitatom].iat[j] != 0 && atom[atom[nitatom].iat[j]].mmx_type != 20)
                            jjk++;
                     }
                     if (jjk == 2)
                        oplstype = 34;
                  }else if ( jji == 4 && nnit == 1)
                  {
                     ambertype = 35; // H1
                     oplstype = 29;
                  }else if ( jji == 3 && nnit == 0)
                  {
                     ambertype = 33; // HA
                     oplstype = 34;
                  }else if ( jji == 4 && nnit == 0)
                  {
                     ambertype = 34; // HC
                     oplstype = 29;
                  }else if ( jji == 4 && nnit == 3)
                  {
                     ambertype = 37; // H3
                     oplstype = 29;
                  }
                  goto L_10;
              } else if (atom[atom[i].iat[0]].atomnum == 7)  // nitrogen
              {
                  jji = 0;
                  jatm = atom[i].iat[0];
                  iatype = atom[jatm].mmx_type;
                  noxide = FALSE;
                  ismetal = FALSE;
                  mmxtype = 23;
                  mm3type = 23;
                  mmfftype = 23;
                  ambertype = 29;
                  oplstype = 38;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[jatm].iat[j] != 0 && atom[atom[jatm].iat[j]].mmx_type != 20 )
                      {
                         jji += atom[jatm].bo[j];
                         if (atom[atom[jatm].iat[j]].mmx_type >= 300)
                           ismetal = TRUE;
                         if (atom[atom[jatm].iat[j]].atomnum == 8)
                         {
                             katm = atom[jatm].iat[j];
                             jjk = 0;
                             for (k=0; k < MAXIAT; k++)
                             {
                                 if (atom[katm].iat[k] != 0 && atom[katm].bo[k] != 9)
                                    jjk++;
                             }
                             if (jjk == 1)
                                noxide = TRUE;
                         }
                      }
                  }
                  if (jji == 4)   //  N+
                  {
                      mmxtype = 24;
                      mm3type = 48;
                      if (ismetal == TRUE)
                      {
                          mmxtype = 23;
                          mm3type = 23;
                      }
                      mmfftype = 36;
                      oplstype = 40;
                      if (noxide == TRUE)
                         mmfftype = 23;
                      goto L_10;
                  }
                  if (atom[jatm].flags & aromatic_mask)
                  {
                      mmxtype = 23;
                      mm3type = 23;
                      mmfftype = 23;
                      oplstype = 38;
                      if (is_cyclo5(jatm,array) )
                      {
                          jjk = 0;
                          for (j=0; j < 5; j++)
                          {
                              if (atom[array[j]].atomnum == 7)
                                 jjk++;
                          }
                          if (jjk == 2)
                          {
                              jj_bo = 0;
                              for (j=0; j < 5; j++)
                              {
                                  if (atom[array[j]].atomnum == 7 && array[j] != jatm)
                                  {
                                      katm = array[j];
                                      nh = 0;
                                      for(k=0; k < MAXIAT; k++)
                                      {
                                          if (atom[katm].iat[k] != 0 && atom[katm].bo[k] != 9 && atom[atom[katm].iat[k]].mmx_type != 20)
                                          {
                                             jj_bo += atom[katm].bo[k];
                                             if (atom[atom[katm].iat[k]].atomnum == 1)
                                               nh++;
                                          }
                                      }
                                      if (jj_bo == 4 && nh > 0)
                                      {
                                          mmfftype = 36;
                                          goto L_10;
                                      }
                                  }
                              }
                          }
                      }
 //                     goto L_10;
                  }
                  if (atom[atom[i].iat[0]].mmff_type == 56 || atom[atom[i].iat[0]].mmff_type == 55)
                  {
                       mmfftype = 36;
                       goto L_10;
                  }
                  if (atom[atom[i].iat[0]].mmff_type == 62)
                  {
                      mmfftype = 23;
                      goto L_10;
                  }
                  if (atom[atom[i].iat[0]].mmff_type == 81)
                  {
                      mmfftype = 36;
                      goto L_10;
                  }
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[jatm].iat[j] != 0 && atom[atom[jatm].iat[j]].mmx_type != 20)
                      {
                         if (atom[jatm].bo[j] == 2 && (atom[atom[jatm].iat[j]].atomnum == 6 || atom[atom[jatm].iat[j]].atomnum == 7))  // imine
                         {
                             mmxtype = 23;
                             mm3type = 23;
                             mmfftype = 27;
                             goto L_10;
                         }
                         if (atom[jatm].bo[j] == 2 && atom[atom[jatm].iat[j]].atomnum == 16 )  // h-n=s
                         {
                             mmxtype = 23;
                             mm3type = 23;
                             mmfftype = 28;
                             goto L_10;
                         }
                         if (atom[jatm].bo[j] == 2 && atom[atom[jatm].iat[j]].atomnum == 15 )  // h-n=p
                         {
                             mmxtype = 23;
                             mm3type = 23;
                             mmfftype = 28;
                             goto L_10;
                         }
                         if (atom[atom[jatm].iat[j]].atomnum == 16 && atom[atom[jatm].iat[j]].mmff_type == 18) // thioamide
                         {
                             mmxtype = 23;
                             mm3type = 23;
                             mmfftype = 28;
                             goto L_10;
                         }
                         if (atom[atom[jatm].iat[j]].atomnum == 15 )  // phosphonamide
                         {
                             katm = atom[jatm].iat[j];
                             for (k=0; k < MAXIAT; k++)
                             {
                                 if (atom[katm].iat[k] != 0 && atom[katm].bo[k] >= 2)
                                 {
                                     mmxtype = 23;
                                     mm3type = 28;
                                     mmfftype = 28;
                                     goto L_10;
                                  }
                             }                             
                         }
                         if (atom[atom[jatm].iat[j]].atomnum == 6 || atom[atom[jatm].iat[j]].atomnum == 7 || atom[atom[jatm].iat[j]].atomnum == 16)  // amide and enamine and sulfonamide
                         {
                             katm = atom[jatm].iat[j];
                             for (k=0; k < MAXIAT; k++)
                             {
                                 if (atom[katm].iat[k] != 0 && atom[katm].iat[k] != jatm)
                                 {
                                     if (atom[katm].bo[k] == 3 && atom[atom[katm].iat[k]].atomnum == 7)
                                     {
                                         mmxtype = 23;
                                         mm3type = 23;
                                         if (atom[jatm].mmff_type != 8)
                                            mmfftype = 28;
                                         else
                                            mmfftype = 23;
                                         goto L_10;
                                     } else if (atom[katm].bo[k] == 3)
                                     {
                                         mmfftype = 28;
                                         goto L_10;
                                     }
                                     if (atom[katm].bo[k] == 2)
                                     {
                                         mmfftype = 28;
                                         if (atom[atom[katm].iat[k]].atomnum == 8 || atom[atom[katm].iat[k]].atomnum == 16)  // amide
                                         {
                                             mmxtype = 23;
                                             mm3type = 28;
                                             mmfftype = 28;
                                             oplstype = 39;
                                             goto L_10;
                                         } else if (atom[atom[katm].iat[k]].atomnum == 6 || atom[atom[katm].iat[k]].atomnum == 7)
                                         {
                                             mmxtype = 23;
                                             mm3type = 28;
                                             if (atom[jatm].mmff_type != 8)
                                               mmfftype = 28;
                                             else
                                               mmfftype = 23;
                                             if (atom[atom[i].iat[0]].mmff_type == 56 || atom[atom[i].iat[0]].mmff_type == 55)
                                                mmfftype = 36;
                                             goto L_10;
                                         }
                                     }
                                 }
                             }
                         }
                      }
                  }
                  mmxtype = 23;  // amine 
                  mm3type = 23;
                  mmfftype = 23;
                  goto L_10;
              } else if (atom[atom[i].iat[0]].atomnum == 8) // oxygen
              {
                  mmxtype = 21;
                  mm3type = 21;
                  mmfftype = 21;
                  ambertype = 31;
                  oplstype = 36;
                  jatm = atom[i].iat[0];
                  if (atom[jatm].iat[0] != i)
                    katm = atom[jatm].iat[0];
                  else
                    katm = atom[jatm].iat[1];
                  jji = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[jatm].iat[j] != 0 && atom[atom[jatm].iat[j]].atomnum == 1)
                        jji++;
                  }
                  if (atom[katm].atomnum == 1 && jji == 2) // water
                  {
                      mmxtype = 21;
                      mm3type = 21;
                      mmfftype = 31;
                      ambertype = 30;
                      oplstype = 41;
                      goto L_10;
                  } else if (atom[katm].atomnum == 1 && jji == 3) // h3o+
                  {
                      mmxtype = 24;
                      mm3type = 21;
                      mmfftype = 50;
                      ambertype = 31;
                      oplstype = 36;
                      goto L_10;
                  }
                  if (atom[katm].atomnum == 15) // h-o-p
                  {
                      mmxtype = 24;
                      mm3type = 24;
                      mmfftype = 24;
                      ambertype = 31;
                      oplstype = 36;
                      goto L_10;
                  }
                  if (atom[katm].atomnum == 16) // h-o-s
                  {
                      mmxtype = 24;
                      mm3type = 24;
                      mmfftype = 33;
                      ambertype = 31;
                      oplstype = 36;
                      goto L_10;
                  }
                  jji = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] != 9 && atom[atom[jatm].iat[j]].mmx_type != 20)
                         jji++;
                  }
                  if (jji == 3 || jji == 4)  // O+
                  {
                      mmxtype = 24;
                      mm3type = 21;
                      mmfftype = 50;
                      ambertype = 31;
                      goto L_10;
                  }                      
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[atom[katm].iat[j]].atomnum == 8 && atom[katm].iat[j] != jatm)
                      {
                          if (atom[katm].bo[j] == 2)   // carboxyl
                          {
                              mmxtype = 24;
                              mm3type = 24;
                              mmfftype = 24;
                              ambertype = 31;
                              goto L_10;
                          }
                      }
                      if (atom[katm].atomnum == 6 && (atom[atom[katm].iat[j]].atomnum == 6 || atom[atom[katm].iat[j]].atomnum == 7)
                          && atom[katm].iat[j] != jatm)  // enol
                      {
                          if (atom[katm].bo[j] == 2)
                          {
                              mmxtype = 28;
                              mm3type = 73;
                              mmfftype = 29;
                              ambertype = 31;
                              goto L_10;
                          }
                      }
                  }
                  if (jji == 2)  // OH
                  {
                      mmxtype = 21;
                      mm3type = 21;
                      mmfftype = 21;
                      for (j=0; j < MAXIAT; j++)  // H-O=C
                      {
                          if (atom[jatm].iat[j] == katm && atom[jatm].bo[j] == 2)
                          {
                              mmxtype = 24;
                              mmfftype = 52;
                              goto L_10;
                          }
                      }
                      goto L_10;
                  }
                  goto L_10;
              } else if (atom[atom[i].iat[0]].atomnum == 5) // boron
              {
                  mmxtype = 23;
                  if (atom[atom[i].iat[0]].mmx_type == 43)   // ts hydroboration
                    mmxtype = 5;
                  mm3type = 5;
                  mmfftype = 71;
                  ambertype = 32;
                  oplstype = 37;
                  goto L_10;
              } else if (atom[atom[i].iat[0]].atomnum == 9 || atom[atom[i].iat[0]].atomnum == 17 ||
                          atom[atom[i].iat[0]].atomnum == 35 || atom[atom[i].iat[0]].atomnum == 53) // halogens
              {
                  mmxtype = 21;
                  mm3type = 5;
                  mmfftype = 71;
                  ambertype = 32;
                  oplstype = 37;
                  goto L_10;
              } else if (atom[atom[i].iat[0]].atomnum == 16) // sulfur
              {
                  mmxtype = 21;
                  mm3type = 5;
                  mmfftype = 71;
                  ambertype = 32;
                  oplstype = 37;
                  goto L_10;
              } else if (atom[atom[i].iat[0]].atomnum == 15) // phosphorous
              {
                  mmxtype = 23;
                  mm3type = 5;
                  mmfftype = 71;
                  ambertype = 0;
                  oplstype = 36;
                  goto L_10;
              } else  //  bridging hydrogens
              {
                  jji = 0;
                  mm3type = 5;
                  for (j=0; j < MAXIAT; j++)
                     if (atom[i].iat[j] != 0)
                        jji++;
                  if (jji == 2)
                    mmxtype = 70;
                  else
                  {
                    mmxtype = 5;
                    mm3type = 5;
                    mmfftype = 5;
                    ambertype = 31;
                  }
                  goto L_10;
              }
          } else if (atom[i].atomnum == 5) // boron
          {
              if (atom[i].mmx_type == 27)  // four coordinate boron 
              {
                 mmxtype = 27;
                 mm3type = 27;
                 mmfftype = 0;
                 goto L_10;
              }
              if (atom[i].mmx_type == 43)  // Transition state boron 
              {
                 mmxtype = 43;
                 mm3type = 0;
                 mmfftype = 0;
                 goto L_10;
              }
              jji = 0;
              for (j=0; j < MAXIAT; j++)
              {
                  if (atom[i].iat[j] != 0)
                     jji++;
              }
              if (jji == 4)
              {
                 mmxtype = 27;
                 mm3type = 27;
                 mmfftype = 0;
              } else
              {
                 mmxtype = 26;
                 mm3type = 26;
                 mmfftype = 0;
              }
              goto L_10;
// =========================== Carbon ===============================
          } else if (atom[i].atomnum == 6) // carbon
          {
              ndouble = 0;
              ntriple = 0;
              jji = 0;
              icycl3 = find_rsize(3,i);
              icycl4 = find_rsize(4,i);
              icycl5 = find_rsize(5,i);
              icycl6 = find_rsize(6,i);
              
              for (j=0; j < MAXIAT; j++)
              {
                  if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                  {
                      jji++;
                      if (atom[i].bo[j] == 3)
                         ntriple++;
                      if (atom[i].bo[j] == 2)
                         ndouble++;
                  }
              }
              //  check here for types that can not be done by rules but must be set by user
              //  and thus should not be changed
              if (mmxtype == 29)   // C radical
              {
                  mmxtype = 29;
                  mm3type = 29;
                  mmfftype = 0;
                  ambertype = 0;
                  oplstype = 0;
                  goto L_10;
              } else if (mmxtype == 30) // C cation
              {
                  mmxtype = 30;
                  mm3type = 30;
                  mmfftype = 0;
                  ambertype = 0;
                  oplstype = 0;
                  goto L_10;
              } else if (mmxtype == 48) // C anion
              {
                  mmxtype = 48;
                  mm3type = 0;
                  mmfftype = 0;
                  ambertype = 0;
                  oplstype = 0;
                  goto L_10;
              } else if (mmxtype == 40) // aromatic carbon
              {
                  mmxtype = 40;
                  mm3type = 50;
                  if (icycl5 >= 1)
                    mm3type = 113;
                  mmfftype = 37;
                  ambertype = 3;
                  oplstype = 7;
                  goto L_10;
              } else if (mmxtype == 49 || mmxtype == 50 || mmxtype == 51 || mmxtype == 52)  // TS atoms
              {
                  mm3type = 0;
                  mmfftype = 0;
                  ambertype = 0;
                  oplstype = 0;
                 goto L_10;
              }
//   check rings
              if (icycl5 >= 1 && icycl6 >= 1 && atom[i].flags & aromatic_mask)  // CB and CN
              {
                  get_rsize(i,5,0,array);
                  for (j=0; j < 5; j++)
                  {
                      if ( !(atom[array[j]].flags & aromatic_mask))
                         goto L_NOTPURINE;
                  }
                  nnh = 0;
                  nnit = 0;
                  for (j=0; j < 5; j++)
                  {
                      if (atom[array[j]].atomnum == 7)
                      {
                          nnit++;
                          for (k=0; k < 4; k++)
                          {
                              if (atom[atom[array[j]].iat[k]].atomnum == 1)
                                nnh++;
                          }
                      }
                  }
                  if (nnit == 2 && nnh == 2)
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 78;
                      ambertype  = 11;
                      oplstype = 11;
                      goto L_10;
                  }
//
                  nnh = 0;
                  nnit = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                   if (atom[atom[i].iat[j]].atomnum == 7)
                   {
                       nnit++;
                       for (k=0; k < MAXIAT; k++)
                       {
                           if (atom[atom[atom[i].iat[j]].iat[k]].atomnum == 1)
                           {
                               nnh++;
                           }
                       }
                   }
                  }
              }
L_NOTPURINE:
              if (icycl5 >= 1 && atom[i].flags & aromatic_mask)
              {
                  mmxtype = 2;
                  mm3type = 2;
                  for (ij=0; ij < icycl5; ij++)
                  {
                      // do ambertypes first
                      nh = 0;
                      nnit = 0;
                      nnh = 0;
                      for (j=0; j < MAXIAT; j++)
                      {
                          if (atom[atom[i].iat[j]].atomnum == 1)
                              nh++;
                          else if (atom[atom[i].iat[j]].atomnum == 7)
                          {
                              nnit++;
                              for(k=0; k < MAXIAT; k++)
                              {
                                 if (atom[atom[atom[i].iat[j]].iat[k]].atomnum == 1)
                                    nnh++;
                              }
                          }
                      }
                      if (nh ==1 && nnit == 2 && nnh == 1)
                      {
                         ambertype = 8;
                         oplstype = 15;
                      }else if (nh == 1 && nnh == 1)
                      {
                         ambertype = 7;
                         oplstype = 12;
                      }else if (nh == 1 && nnit == 1)
                      {
                         ambertype = 6;
                         oplstype = 13;
                      }else if (nh == 1 && nnit == 2 && nnh == 0)
                      {
                         ambertype = 12;
                         oplstype = 18;
                      } else if (nh == 0 && nnit == 1)
                      {
                         ambertype = 5;
                         oplstype = 12;
                      } else if (nh == 0 && nnit == 0)
                      {
                         ambertype = 10;
                         oplstype = 17;
                      }
                      // do mmff types
                      get_rsize(i,5,ij,array);
                      if (aromatic_5(ij,array))
                      {
                          icurr = -1;
                          nplus = FALSE;
                          for (j=0; j < 5; j++)
                          {
                              if (array[j] == i)
                                icurr = j;
                              if (atom[array[j]].atomnum == 7)
                              {
                                  jjk = 0;
                                  jj_bo = 0;
                                  for (k=0; k < MAXIAT; k++)
                                  {
                                      if (atom[array[j]].iat[k] != 0 && atom[array[j]].bo[k] != 9 && atom[atom[array[j]].iat[k]].mmx_type != 20)
                                      {
                                          jjk++;
                                          jj_bo += atom[array[j]].bo[k];
                                      }
                                  }
                                  if (jjk == 2 && jj_bo == 2) // divalent N anion
                                  {
                                      mmfftype = 78;
                                      goto L_10;
                                  }
                                  if (jj_bo == 4)
                                     nplus = TRUE;
                              }
                          }

                          // check alpha
                          ia = (icurr+4)%5;
                          ib = (icurr+6)%5;
                          if (atom[array[ia]].atomnum == 7 && atom[array[ib]].atomnum == 7)  // n=c-n
                          {
                              jatm = array[ia];
                              katm = array[ib];
                              jbo = 0;
                              kbo = 0;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] != 9 && atom[atom[jatm].iat[j]].mmx_type != 20)
                                     jbo += atom[jatm].bo[j];
                                  if (atom[katm].iat[j] != 0 && atom[katm].bo[j] != 9 && atom[atom[katm].iat[j]].mmx_type != 20)
                                     kbo += atom[katm].bo[j];
                              }
                              if ( jbo == 4 && kbo == 3)
                              {
                                  adjn = FALSE;
                                  for (j=0; j < MAXIAT; j++)
                                  {
                                      if (atom[katm].iat[j] != 0 && atom[katm].bo[j] == 2)
                                         adjn = TRUE;
                                  }
                                  if (adjn == FALSE)
                                  {
                                     mmfftype = 80;
                                     goto L_10;
                                  }
                              }else if (jbo == 3 && kbo == 4)
                              {
                                  adjn = FALSE;
                                  for (j=0; j < MAXIAT; j++)
                                  {
                                      if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] == 2)
                                         adjn = TRUE;
                                  }
                                  if (adjn == FALSE)
                                  {
                                     mmfftype = 80;
                                     goto L_10;
                                  }
                              }
                          }
                          if (atom[array[ia]].atomnum == 7 && atom[array[ib]].atomnum == 6)  // n=c-c
                          {
                              oplstype = 23;
//                              goto L_10;
                          } else if (atom[array[ia]].atomnum == 6 && atom[array[ib]].atomnum == 7)  // n=c-c
                          {
                              oplstype = 23;
 //                             goto L_10;
                          }
                          if (nplus == TRUE) // found N+ in ring
                          {
                              noxide = FALSE;
                              jjk = 0;
                              for (j=0; j < 5; j++)
                              {
                                  if (atom[array[j]].atomnum == 7)
                                  {
                                      nit[jjk] = j;
                                      jjk++;
                                  }
                              }
                              if (jjk >= 2)
                              {
                                  jatm = array[nit[0]];
                                  katm = array[nit[1]];
                                  jbo = 0;
                                  kbo = 0;
                                  adjn = FALSE;
                                  adjn1 = FALSE;
                                  for (k=0; k < MAXIAT; k++)
                                  {
                                       if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9 && atom[atom[jatm].iat[k]].mmx_type != 20)
                                              jbo += atom[jatm].bo[k];
                                       if (atom[katm].iat[k] != 0 && atom[katm].bo[k] != 9 && atom[atom[katm].iat[k]].mmx_type != 20)
                                              kbo += atom[katm].bo[k];
                                       if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] == 2)
                                              adjn = TRUE;
                                       if (atom[katm].iat[k] != 0 && atom[katm].bo[k] == 2)
                                              adjn1 = TRUE;
                                  }
                                  if ( (jbo == 4 && kbo == 3 && adjn1 == FALSE) ||
                                       (jbo == 3 && kbo == 4 && adjn == FALSE)  )
                                  {
                                      noxide = FALSE;
                                      if (jbo == 4)
                                        ia1 = jatm;
                                      else
                                        ia1 = katm;
                                      for (k=0; k < MAXIAT; k++)
                                      {
                                          if (atom[ia1].iat[k] != 0 && atom[atom[ia1].iat[k]].atomnum == 8)
                                          {
                                              jjk = 0;
                                              ib1 = atom[ia1].iat[k];
                                              for (l=0; l < MAXIAT; l++)
                                              {
                                                  if (atom[ib1].iat[l] != 0 && atom[ib1].bo[l] != 9)
                                                    jjk++;
                                              }
                                              if (jjk == 1)
                                                 noxide = TRUE;
                                          }
                                      }
                                      if (noxide == FALSE)
                                      {             
                                         mmfftype = 78;
                                         goto L_10;
                                      }
                                  }
                              }
                          }

                          if (atom[array[ia]].atomnum != 6 && atom[array[ib]].atomnum != 6)  // x=c-x
                          {
                              jatm = array[ia];
                              katm = array[ib];
                              jbo = 0;
                              kbo = 0;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] != 9 && atom[atom[jatm].iat[j]].mmx_type != 20)
                                     jbo += atom[jatm].bo[j];
                                  if (atom[katm].iat[j] != 0 && atom[katm].bo[j] != 9 && atom[atom[katm].iat[j]].mmx_type != 20)
                                     kbo += atom[katm].bo[j];
                              }
                              if ( (jbo == 4 && kbo == 2) || (jbo == 2 && kbo == 4))
                              {
                                  for (j=0; j < MAXIAT; j++)
                                  {
                                      if (atom[i].iat[j] != 0 && atom[i].iat[j] != jatm && atom[i].iat[j] != katm)
                                      {
                                          if (atom[atom[i].iat[j]].atomnum == 7)
                                          {
                                               mmfftype = 80;
                                               goto L_10;
                                          }
                                      }
                                  }
                              }
                          }
                          if (atom[array[ia]].atomnum == 7)  // alpha n
                          {
                              jatm = array[ia];
                              adjn = FALSE;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] ==2)
                                     adjn = TRUE;
                              }
                              if (adjn == FALSE)
                              {
                                  mmfftype = 63;
                                  goto L_10;
                              }
                          }
                          if (atom[array[ib]].atomnum == 7)  // alpha n
                          {
                              jatm = array[ib];
                              adjn = FALSE;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] == 2)
                                     adjn = TRUE;
                              }
                              if (adjn == FALSE)
                              {
                                  mmfftype = 63;
                                  goto L_10;
                              }
                          }
                          if (atom[array[ia]].atomnum == 8 || atom[array[ib]].atomnum == 8 ||  // alpha o or s
                              atom[array[ia]].atomnum == 16 || atom[array[ib]].atomnum == 16 )
                          {
                              mmfftype = 63;
                              oplstype = 14;
                              goto L_10;
                          }
                          // check beta
                          ia = (icurr+3)%5;
                          ib = (icurr+7)%5;
                          if (atom[array[ia]].atomnum == 7)  // c=x-n
                          {
                              jatm = array[ia];
                              adjn = FALSE;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] == 2)
                                     adjn = TRUE;
                              }
                              if (adjn == FALSE)
                              {
                                  mmfftype = 64;
                                  goto L_10;
                              }
                          }
                          if (atom[array[ib]].atomnum == 7)  // c=x-n
                          {
                              jatm = array[ib];
                              adjn = FALSE;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] == 2)
                                     adjn = TRUE;
                              }
                              if (adjn == FALSE)
                              {
                                  mmfftype = 64;
                                  goto L_10;
                              }
                          }
                          if (atom[array[ia]].atomnum == 8 || atom[array[ib]].atomnum == 8 ||  // beta o or s
                              atom[array[ia]].atomnum == 16 || atom[array[ib]].atomnum == 16 )
                          {
                              mmfftype = 64;
                              oplstype = 22;
                              nnit = 0;
                              for (j=0; j < 5; j++)
                              {
                                  if (atom[array[j]].atomnum != 6)
                                    nnit++;
                              }
                              if (nnit > 1)
                                 oplstype = 7;
                              goto L_10;
                          }
                      }
                  }
              }
              if (icycl6 >= 1 && atom[i].flags & aromatic_mask)
              {
                  mmxtype = 2;
                  mm3type = 2;
                  mmfftype = 37;
                  for (k = 0; k < icycl6; k++)
                  {
                      get_rsize(i,6,k,array);
                      {
                          full_ring = TRUE;
                          non_pi = FALSE;
                          for (j=0; j < 6; j++)  // check all ring atoms are pi
                          {
                              if ( !(atom[array[j]].flags & pi_mask) )
                                 non_pi = TRUE;
                          }
                          if (non_pi == TRUE)
                             full_ring = FALSE;
                      }
                      if (full_ring == TRUE)
                         break;
                  }
                  if (full_ring == FALSE) mmfftype = 2;
 //
                  ambertype = 3;
                  oplstype = 7;
                  nh = 0;
                  nnit = 0;
                  no = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                    if (atom[atom[i].iat[j]].atomnum == 1)
                      nh++;
                    else if (atom[atom[i].iat[j]].atomnum == 7)
                      nnit++;
                    else if (atom[atom[i].iat[j]].atomnum == 8)
                    {
                      no++;
                      if (atom[i].bo[j] >= 2) mmfftype = 3;
                    }
                  }
                  if (nh ==1 && nnit == 2)
                  {
                     ambertype = 13;
                     oplstype = 19;
                     goto L_10;
                  }
                  if (no == 1)
                  {
                    ambertype = 2; // tyrosine ring
                    oplstype =  7;
                    goto L_10;
                  }
                  if (jji == 3) // aromatic carbon from pdb file
                  {
                     ambertype = 3;
                     oplstype =  7;
                     goto L_10;
                  }
              }
              if (icycl3 >= 1 )
              {
                  mmxtype = 22;
                  mm3type = 22;  // type 38 for cyclopropene, type 67 for cyclopropanone
                  mmfftype = 22;
                  if (ndouble >= 1)
                  {
                      for (j=0; j < MAXIAT; j++)
                      {
                          if (atom[i].iat[j] != 0)
                          {
                              if (atom[i].bo[j] == 2 && atom[atom[i].iat[j]].atomnum == 8)
                              {
                                  mmxtype = 3;
                                  mm3type = 67;
                                  mmfftype = 3;
                                  goto L_10;
                              }
                          }
                      }
                     mm3type = 38;
                  }
                  mmfftype = 22;
                  goto L_10;
              }
//   start of check on bondorders
              if (ntriple == 1)
              {
                  mmxtype = 4;
                  mm3type = 4;
                  mmfftype = 4;
                  ambertype = 0;
                  oplstype =  24;
                  if (jji == 1 && atom[atom[i].iat[0]].atomnum == 7)
                  {
                     mmfftype = 60;  // isonitrile carbon
                     oplstype =  0;
                  }
                  if ( (atom[atom[i].iat[0]].mmx_type >= 300) &&
                       atom[atom[i].iat[1]].mmx_type == 46)
                  {
                         mmxtype = 63;
                         goto L_10;
                  }
                  if ( (atom[atom[i].iat[1]].mmx_type >= 300) &&
                       atom[atom[i].iat[0]].mmx_type == 46)
                  {
                         mmxtype = 63;
                         goto L_10;
                  }
                  for (j=0; j < MAXIAT; j++)
                  {
                      if ( (atom[atom[i].iat[j]].mmx_type >= 300) && atom[i].bo[j] == 3)
                      {
                         mmxtype = 62;
                         goto L_10;
                      }
                  }
                  if (atom[atom[i].iat[0]].atomnum == 1 || atom[atom[i].iat[1]].atomnum == 1)
                    oplstype = 24;
                  goto L_10;
              }
              if (ndouble == 2)   // allenes and ketenes
              {
                  mmxtype = 4;
                  mm3type = 68;
                  mmfftype = 4;
                  ambertype = 0;
                  for (j=0; j < MAXIAT; j++)  // metal carbene
                  {
                      if ( atom[atom[i].iat[j]].mmx_type >= 300)
                         mmxtype = 61;
                  }
                  goto L_10;
              }
              if (ndouble == 1)
              {
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0)
                      {
                          if (atom[i].bo[j] == 2)
                          {
                              jatm = atom[i].iat[j];
                              break;
                          }
                      }
                  }
                  if (atom[jatm].atomnum == 15) // c=p
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 3;
                      oplstype = 10;
                      goto L_10;
                  }
                  if (atom[jatm].atomnum == 16) // c=s
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 3;
                      ambertype = 3;
                      oplstype = 10;
                      for (k=0; k < MAXIAT; k++)
                      {
                          if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
                          {
                              if (atom[atom[i].iat[k]].atomnum == 16 && atom[i].iat[k] != jatm)
                              {
                                  katm = atom[i].iat[k];
                                  jjk = 0;
                                  for (l=0; l < MAXIAT; l++)
                                  {
                                      if (atom[katm].iat[l] != 0 && atom[katm].bo[l] != 9)
                                        jjk++;
                                  }
                                  if (jjk == 1)  // thiocarboxylate
                                  {
                                      mmfftype = 41;
                                      goto L_10;
                                  }
                              }
                          }
                      }    
                      goto L_10;
                  }
                  if (atom[jatm].atomnum == 8)  // C=O
                  {
                      mmxtype = 3;
                      mm3type = 3;
                      mmfftype = 3;
                      ambertype = 2;
                      oplstype = 5;
                      if (is_ring31(i))          // cyclopropanone
                          mm3type = 67;
                      else if (is_ring41(i))   // cyclobutanone
                          mm3type = 58;
                      else                     // carboxylate
                      {
                          for(k=0; k < MAXIAT; k++)
                          {
                              if (atom[i].iat[k] != 0)
                              {
                                  if (atom[atom[i].iat[k]].atomnum == 8 && atom[i].iat[k] != jatm)
                                  {
                                      katm = atom[i].iat[k]; 
                                      jjk = 0;
                                      for (l=0; l < MAXIAT; l++)
                                      {
                                          if (atom[katm].iat[l] != 0 && atom[atom[katm].iat[l]].mmx_type != 20)
                                            jjk++;
                                      }
                                      if (jjk == 1)
                                      {
                                           mmfftype = 41;
                                           oplstype = 5;
                                           break;
                                      }
                                  }
                              }
                          }
                      }
                      goto L_10;
                  } else if (atom[jatm].atomnum == 7) // C=N
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 3;
                      ambertype = 3;
                      oplstype = 10;
                      if (jji == 3)
                      {
                          if (atom[i].iat[0] == jatm)
                          {
                              katm = atom[i].iat[1];
                              latm = atom[i].iat[2];
                          } else if (atom[i].iat[1] == jatm)
                          {
                              katm = atom[i].iat[0];
                              latm = atom[i].iat[2];
                          } else if (atom[i].iat[2] == jatm)
                          {
                              katm = atom[i].iat[0];
                              latm = atom[i].iat[1];
                          }
                          if (atom[jatm].atomnum == 7 && atom[katm].atomnum == 7 && atom[latm].atomnum == 7)
                          {
                              jbo = 0;
                              kbo = 0;
                              lbo = 0;
                              adjn = FALSE;
                              oplstype = 8;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] != 9 && atom[atom[jatm].iat[j]].mmx_type != 20)
                                       jbo += atom[jatm].bo[j];
                                  if (atom[katm].iat[j] != 0 && atom[katm].bo[j] != 9 && atom[atom[katm].iat[j]].mmx_type != 20)
                                  {
                                       kbo += atom[katm].bo[j];
                                       if (atom[katm].bo[j] == 2)
                                           adjn = TRUE;
                                  }
                                  if (atom[latm].iat[j] != 0 && atom[latm].bo[j] != 9)
                                  {
                                       lbo += atom[latm].bo[j];
                                       if (atom[latm].bo[j] == 2)
                                           adjn = TRUE;
                                  }
                              }
                              if (jbo == 4 && kbo == 3 && lbo == 3 && adjn == FALSE)
                              {  
                                  mmfftype = 57;
                                  oplstype = 8;
                                  goto L_10;
                              }
                          }
                          if (atom[jatm].atomnum == 7 && atom[katm].atomnum == 7)
                          {
                              jbo = 0;
                              kbo = 0;
                              adjn = FALSE;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] != 9 && atom[atom[jatm].iat[j]].mmx_type != 20)
                                       jbo += atom[jatm].bo[j];
                                  if (atom[katm].iat[j] != 0 && atom[katm].bo[j] != 9 && atom[atom[katm].iat[j]].mmx_type != 20)
                                  {
                                       kbo += atom[katm].bo[j];
                                       if (atom[katm].bo[j] == 2)
                                           adjn = TRUE;
                                  }
                              }
                              if (jbo == 4 && kbo == 3 && adjn == FALSE)
                              {  
                                  mmfftype = 57;
                                  oplstype = 8;
                                  goto L_10;
                              }
                              if (jbo == 3 && kbo == 3 && adjn == FALSE) // cytosine c4
                              {
                                  if (atom[latm].atomnum != 7)
                                  {
                                     oplstype = 9;
                                     goto L_10;
                                  }
                              }
                          }
                          if (atom[jatm].atomnum == 7 && atom[latm].atomnum == 7)
                          {
                              jbo = 0;
                              lbo = 0;
                              adjn = FALSE;
                              for (j=0; j < MAXIAT; j++)
                              {
                                  if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] != 9 && atom[atom[jatm].iat[j]].mmx_type != 20)
                                       jbo += atom[jatm].bo[j];
                                  if (atom[latm].iat[j] != 0 && atom[latm].bo[j] != 9 && atom[atom[latm].iat[j]].mmx_type != 20)
                                  {
                                       lbo += atom[latm].bo[j];
                                       if (atom[latm].bo[j] == 2)
                                           adjn = TRUE;
                                  }
                              }
                              if (jbo == 4 && lbo == 3 && adjn == FALSE)
                              {  
                                  mmfftype = 57;
                                  oplstype = 8;
                                  goto L_10;
                              }
                              if (jbo == 3 && lbo == 3 && adjn == FALSE) // cytosine c4
                              {
                                  if (atom[katm].atomnum != 7)
                                  {
                                     oplstype = 9;
                                     goto L_10;
                                  }
                              }
                          }
                      }
                      goto L_10;
                  } else if (atom[jatm].atomnum == 6) // C=C
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 2;
                      ambertype = 4;
                      oplstype = 10;
                      if (is_ring31(i))    // cyclopropene
                      {
                          mm3type = 38;
                          goto L_10;
                      }
                      if (is_ring41(i))  // cyclobutene
                      {
                          mmxtype = 57;
                          mm3type = 57;
                          mmfftype = 30;
                          goto L_10;
                      }

                      goto L_10;
                  } else   // default for cases not dealt with yet
                  {
                      mmxtype = 2;
                      mm3type = 2;
                      mmfftype = 2;
                      ambertype = 4;
                      oplstype = 10;
                      goto L_10;
                  }
                  goto L_10;
              }
      // get here with only single bonds to carbon
              mmxtype = 1;
              mm3type = 1;
              mmfftype = 1;
              ambertype = 1;
              oplstype = 1;
              if (is_ring31(i))
              {
                  mmxtype = 22;
                  mm3type = 22;
                  mmfftype = 22;
                  oplstype = 28;
                  goto L_10;
              } else if (is_ring41(i))
              {
                  mmxtype = 56;
                  mm3type = 56;
                  mmfftype = 20;
                  goto L_10;
              } else  // look for adjacent triple bonds for opls
              {
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[atom[i].iat[j]].atomnum == 6)
                      {
                          jatm = atom[i].iat[j];
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] == 3)
                              {
                                  oplstype = 25;
                                  nh = 0;
                                  for (l = 0; l < MAXIAT; l++)
                                  {
                                      if (atom[i].iat[l] != 0 && atom[atom[i].iat[l]].atomnum == 1)
                                          nh++;
                                  }
                                  if (nh == 0) oplstype = 27;
                                  if (nh == 1) oplstype = 26;
                                  if (nh == 2) oplstype = 25;
                                  goto L_10;
                              }
                          }
                      }
                  }
              }
              // acetal type
              nh = 0;
              for (j=0; j < MAXIAT; j++)
              {
                  if (atom[i].iat[j] != 0 && atom[atom[i].iat[j]].atomnum == 8 && atom[i].bo[j] != 2)
                    nh++;
              }
              if (nh == 2 || nh == 3)
                  oplstype = 21;
              goto L_10;
          } else if (atom[i].atomnum == 7) // nitrogen
          {
              jji = 0;
              nc = 0;
              nh = 0;
              jjbo = 0;
              ndouble = 0;
              ntriple = 0;
              nplus = FALSE;
              noxide = FALSE;
              for (j=0; j < MAXIAT; j++)
              {
                if (atom[i].iat[j] != 0 && atom[atom[i].iat[j]].mmx_type != 20 && atom[atom[i].iat[j]].mmx_type < 300)
                  jji++;
                if (atom[atom[i].iat[j]].atomnum == 1)
                  nh++;
                if (atom[atom[i].iat[j]].atomnum == 6)
                  nc++;
              }
              icycl5 = find_rsize(5,i);
              icycl6 = find_rsize(6,i);
// amber types
              if (icycl5 >= 1 && nh == 0 && (nc == 2 || nc == 1))
              {
                  ambertype = 16; // NB
                  oplstype = 59;
                  goto L_1001;
              }else if (icycl5 >=1 && nh <= 1 && (nc == 2 || nc == 1))
              {
                // need to check for proline
                for (j=0; j < MAXIAT; j++)
                {
                    if (atom[atom[i].iat[j]].atomnum == 6)
                    {
                        jatm = atom[i].iat[j];
                        for (k=0; k < MAXIAT; k++)
                        {
                            if (atom[atom[jatm].iat[k]].atomnum == 6)
                            {
                                katm = atom[jatm].iat[k];
                                for (l=0; l < MAXIAT; l++)
                                {
                                    if (atom[katm].iat[l] != 0 && atom[katm].bo[l] == 2 && atom[atom[katm].iat[l]].atomnum == 8)
                                    {
                                        ambertype = 14;
                                        oplstype = 53;
                                        goto L_1001;
                                    }
                                }
                            }
                        }
                    }
                }       
                ambertype = 15; // NA
                oplstype = 58;
                goto L_1001;
              } else if (icycl5 >= 1 && nh == 0 && nc == 3)
              {
                  for (j=0; j < MAXIAT; j++)
                  {
                    if (atom[atom[i].iat[j]].atomnum == 6)
                    {
                        jatm = atom[i].iat[j];
                        for (k=0; k < MAXIAT; k++)
                        {
                            if (atom[atom[jatm].iat[k]].atomnum == 8 && atom[jatm].bo[k] == 2)
                            {
                                ambertype = 14;
                                oplstype = 55;
                                goto L_1001;
                            }
                        }
                    }
                  }
                  ambertype = 18; // N*
                  oplstype = 61;
                  goto L_1001;
              }
              if (icycl6 >= 1)
              {
                  if (nh == 1)
                  {
                    ambertype = 15; // NA
                    oplstype = 58;
                  } else if (nh == 0 && nc == 3)
                  {
                    ambertype = 18; // N*;
                    oplstype = 61;
                  }else if (nh == 0 && nc == 2)
                  {
                     ambertype = 17; // NC
                     oplstype = 60;
                  }
                  goto L_1001;
              }
              if (jji == 4)
              {
                 ambertype = 20; // N3
                 oplstype = 56;
                 goto L_1001;
              }
            // only amides and guanadinium ions left
              for (j=0; j< MAXIAT; j++)
              {
                 if (atom[atom[i].iat[j]].atomnum == 6)
                 {
                    jatm = atom[i].iat[j];
                    for (k=0; k < MAXIAT; k++)
                    {
                        if (atom[atom[jatm].iat[k]].atomnum == 8 && atom[jatm].bo[k] == 2)
                        {
                            ambertype = 14; // N amide
                            oplstype = 55;
                            goto L_1001;
                        }
                        if (atom[atom[jatm].iat[k]].atomnum == 7 && atom[jatm].bo[k] == 2)
                        {
                            ambertype = 19; // guanadinium
                            oplstype = 57;
                            goto L_1001;
                        }
                        if (atom[atom[jatm].iat[k]].atomnum == 6 && atom[jatm].bo[k] == 2)
                        {
                            ambertype = 19; // guanadinium
                            oplstype = 54;
                            goto L_1001;
                        }
                    }
                 }
              }
              ambertype = 20; // N2
              oplstype = 53;
L_1001:                 
//
              if (icycl5 >= 1 && atom[i].flags & aromatic_mask)
              {
                  jji = 0;
                  jjbo = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].mmx_type != 20 &&
                         atom[atom[i].iat[j]].mmx_type < 300 )
                      {
                          jji++;
                          jjbo += atom[i].bo[j];
                          if (atom[atom[i].iat[j]].atomnum == 8)
                          {
                              katm = atom[i].iat[j];
                              jjk = 0;
                              for (k=0; k < MAXIAT; k++)
                              {
                                  if (atom[katm].iat[k] != 0 && atom[katm].bo[k] != 9)
                                     jjk++;
                              }
                              if (jjk == 1)
                                 noxide = TRUE;
                          }
                      }
                  }
                  get_rsize(i,5,0,array);
                  for(j=0; j < 5; j++)
                  {
                      if (atom[array[j]].atomnum == 7)
                      {
                          jjk = 0;
                          jj_bo = 0;
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[array[j]].iat[k] != 0 && atom[array[j]].bo[k] != 9)
                              {
                                  jjk++;
                                  jj_bo += atom[array[j]].bo[k];
                              }
                          }
                          if (jji == 2 && jj_bo == 2)  // found divalent anion
                          {
                              mmfftype = 76;
                              goto L_10;
                          }
                      }
                  }
                  if (jji == 2 && jjbo == 2)  // divalent anion
                  {
 //                     mmxtype = 77;
                      mmfftype = 76;
                      goto L_10;
                  }
                  if (jjbo == 4 && noxide == TRUE)
                  {
                      mmxtype = 41;
                      mm3type = 39;
                      mmfftype = 82;
                      goto L_10;
                  }
                  if (jjbo == 4)
                  {
                      mmxtype = 41;
                      mm3type = 39;
                      mmfftype = 81;
                      goto L_10;
                  }
                  // assume only one ring need to fix this
                  get_rsize(i,5,0,array);
                  ndouble = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[i].bo[j] == 2)
                        ndouble++;
                  }
                  if (ndouble == 1)  // check for alpha or beta N,O,S
                  {
                      mmxtype = 37;
                      mm3type = 37;
                      icurr  = -1;
                      for (j=0; j < 5; j++)
                      {
                          if (array[j] == i)
                             icurr = j;
                      }
                      if (icurr == 0)
                      {
                          ia = array[1];
                          ib = array[4];
                      } else if (icurr == 4)
                      {
                          ia = array[0];
                          ib = array[3];
                      } else
                      {
                          ia = array[icurr-1];
                          ib = array[icurr+1];
                      }
                      
                      if (atom[ia].atomnum == 7)  // alpha
                      {
                         adjn = FALSE;
                         for (k=0; k < MAXIAT; k++)
                         {
                             if (atom[ia].iat[k] != 0 && atom[ia].bo[k] == 2)
                                adjn = TRUE;
                         }
                         if (adjn == FALSE)
                         {
                             mmfftype = 65;
                             goto L_10;
                         }
                      }
                      if (atom[ib].atomnum == 7)
                      {
                          adjn = FALSE;
                          for (k=0; k < MAXIAT; k++)
                          {
                             if (atom[ib].iat[k] != 0 && atom[ib].bo[k] == 2)
                                      adjn = TRUE;
                          }
                          if (adjn == FALSE)
                          {
                              mmfftype = 65;
                              goto L_10;
                          }
                      }
                      if (atom[ia].atomnum == 8 || atom[ia].atomnum == 16)
                      {
                              mmfftype = 65;
                              goto L_10;
                      }
                      if (atom[ib].atomnum == 8 || atom[ib].atomnum == 16)
                      {
                              mmfftype = 65;
                              goto L_10;
                      }
                      // now check beta
                      if (icurr == 0)
                      {
                          ia = array[2]; ib = array[3];
                      } else if (icurr == 1)
                      {
                          ia = array[3]; ib = array[4];
                      } else if (icurr == 2)
                      {
                          ia = array[4]; ib = array[0];
                      } else if (icurr == 3)
                      {
                          ia = array[0]; ib = array[1];
                      } else if (icurr == 4)
                      {
                          ia = array[1]; ib = array[2];
                      }
                      if (atom[ia].atomnum == 7)
                      {
                          jjk = 0;
                          jj_bo = 0;
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[ia].iat[k] != 0 && atom[ia].bo[k] != 9)
                              {
                                  jjk++;
                                  jj_bo += atom[ia].bo[k];
                              }
                          }
                          if (jjk == 3 && jj_bo == 3)
                          {
                              mmfftype = 66;
                              goto L_10;
                          }                             
                      }
                      if (atom[ib].atomnum == 7)
                      {
                          jjk = 0;
                          jj_bo = 0;
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[ib].iat[k] != 0 && atom[ib].bo[k] != 9)
                              {
                                  jjk++;
                                  jj_bo += atom[ib].bo[k];
                              }
                          }
                          if (jjk == 3 && jj_bo == 3)
                          {
                              mmfftype = 66;
                              goto L_10;
                          }                             
                      }
                      if (atom[ia].atomnum == 8 || atom[ia].atomnum == 16)
                      {
                          mmfftype = 66;
                          goto L_10;
                      }
                      if (atom[ib].atomnum == 8 || atom[ib].atomnum == 16)
                      {
                          mmfftype = 66;
                          goto L_10;
                      }
                  }
                  // single bonds only
                  icurr = -1;
                  nplus = FALSE;
                  for (j=0; j < 5; j++)
                  {
                      if (array[j] == i)
                         icurr = j;
                      if (atom[array[j]].atomnum == 7)
                         nplus = TRUE;
                  }
                  if (nplus == TRUE)
                  {
                      if (icurr == 0)
                      {
                          ia = array[2]; ib = array[3];
                      } else if (icurr == 1)
                      {
                          ia = array[3]; ib = array[4];
                      } else if (icurr == 2)
                      {
                          ia = array[4]; ib = array[0];
                      } else if (icurr == 3)
                      {
                          ia = array[0]; ib = array[1];
                      } else if (icurr == 4)
                      {
                          ia = array[1]; ib = array[2];
                      }
                      if (atom[ia].atomnum == 7)
                      {
                          jj_bo = 0;
                          noxide = FALSE;
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[ia].iat[k] != 0 && atom[ia].bo[k] != 9)
                                 jj_bo += atom[ia].bo[k];
                              if (atom[atom[ia].iat[k]].atomnum == 8)
                              {
                                  jjk = 0;
                                  katm = atom[ia].iat[k];
                                  for (l=0; l < MAXIAT; l++)
                                  {
                                      if (atom[katm].iat[l] != 0 && atom[katm].bo[l] != 9)
                                         jjk++;
                                  }
                                  if (jjk == 1)
                                     noxide = TRUE;
                              }
                          }
                          if (jj_bo == 4 && noxide == FALSE)
                          {
                              mmfftype = 81;  // histidine ??
                              goto L_10;
                          }
                      }
                      if (atom[ib].atomnum == 7)
                      {
                          jj_bo = 0;
                          noxide = FALSE;
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[ib].iat[k] != 0 && atom[ib].bo[k] != 9)
                                 jj_bo += atom[ib].bo[k];
                              if (atom[atom[ib].iat[k]].atomnum == 8)
                              {
                                  jjk = 0;
                                  katm = atom[ib].iat[k];
                                  for (l=0; l < MAXIAT; l++)
                                  {
                                      if (atom[katm].iat[l] != 0 && atom[katm].bo[l] != 9)
                                         jjk++;
                                  }
                                  if (jjk == 1)
                                     noxide = TRUE;
                              }
                          }
                          if (jj_bo == 4 && noxide == FALSE)
                          {
                              mmfftype = 81;  // histidine ??
                              goto L_10;
                          }
                      }
                  }
                  // need test for pyrol here
                  ndouble = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                      {
                         katm = atom[i].iat[j];
                         for (k=0; k < MAXIAT; k++)
                         {
                           if (atom[katm].iat[k] != 0 && atom[katm].bo[k] == 2)
                              ndouble++;
                         }
                     }
                  }
                  if (ndouble >= 1)
                  {
                     mmxtype = 9;
                     mm3type = 9;
                  }
                  // failed tests use general type
                  mmfftype = 39;
                  goto L_10;
              }
              if (icycl6 >= 1 && atom[i].flags & aromatic_mask)
              {
                  mmxtype = 37;
                  mm3type = 37;
                  mmfftype = 38;
                  nnit = 0;
                  get_rsize(i,6,0,array);
                  for (j=0; j < 6; j++)
                  {
                      if (atom[array[j]].atomnum == 7)
                        nnit++;
                  }
                  jji = 0;
                  jjbo = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].mmx_type != 20 &&
                      atom[atom[i].iat[j]].mmx_type < 300 )
                      {
                          jji++;
                          jjbo += atom[i].bo[j];
                          if (atom[atom[i].iat[j]].atomnum == 8)
                          {
                              katm = atom[i].iat[j];
                              jjk =0;
                              for (k=0; k < MAXIAT; k++)
                              {
                                  if (atom[katm].iat[k] != 0 && atom[katm].bo[k] != 9 && atom[atom[katm].iat[k]].mmx_type != 20)
                                     jjk++;
                              }
                              if (jjk == 1)
                                 noxide = TRUE;
                          }
                      }
                  }
                  if (jjbo == 4 && noxide == TRUE)
                  {
                      mmxtype = 41;
                      mm3type = 143;
                      mmfftype = 69;
                  } else if (jjbo == 4)
                  {
                      mmxtype = 41;
                      mm3type = 111;
                      mmfftype = 58;
                  }
                  goto L_10;
              }
// non cyclic systems
              jji = 0;
              jjbo = 0;
              ndouble = 0;
              ntriple = 0;           
              for (j=0; j < MAXIAT; j++)
              {
                  if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].mmx_type != 20 && atom[atom[i].iat[j]].mmx_type < 300)
                  {
                      jji++;
                      jjbo += atom[i].bo[j];
                      if (atom[i].bo[j] == 2)
                         ndouble++;
                      if (atom[i].bo[j] == 3)
                         ntriple++;
                  }
              }
              if (ntriple == 1)   // nitriles and isonitriles
              {
                  mmxtype = 10;
                  mm3type = 10;
                  if (jji == 1)
                     mmfftype = 42;
                  else if (jji == 2)
                  {
                     mmfftype = 61;
                     mmxtype = 68;
                     // nitrile oxides
                     if (atom[atom[i].iat[0]].atomnum == 8 || atom[atom[i].iat[1]].atomnum == 8)
                        mmxtype = 41;
                  }
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[atom[i].iat[j]].mmx_type >= 300)
                          mmxtype = 65;
                  }
                  goto L_10;
              }
        // two double bonds to N - azides
              if (ndouble == 2 && jji == 2)
              {
                  mmxtype = 37;
                  mm3type = 46;
                  mmfftype = 53;
                  goto L_10;
              } else if (ndouble == 2 && jji == 3) // nitrate
              {
                  mmxtype = 41;
                  mm3type = 45;
                  mmfftype = 45;
                  goto L_10;
              }
          // single double bond to N
              if (ndouble == 1)
              {
                  mmxtype = 37;
                  mm3type = 37;
                  mmfftype = 9;
                  if (jji == 3 && jjbo == 4) //=N+
                  {
                      mmxtype = 41;
                      mm3type = 39;
//                      goto L_10;
                  }
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[i].bo[j] == 2)
                      {
                          jatm = atom[i].iat[j];
                          break;
                      }
                  }
                  if (jji == 1 && atom[jatm].atomnum == 7) // N=N=X  azido and diazo
                  {
                      mmfftype = 47;
                      goto L_10;
                  }
                  if (atom[jatm].atomnum == 8)  // N=O
                  {
                      if (jjbo == 4)   // nitro
                      {
                          mmxtype = 41;
                          mm3type = 46;
                          mmfftype = 45;
                          goto L_10;
                      } else if (jjbo == 3)  // nitroso
                      {
                          mmxtype = 37;
                          mm3type = 0;
                          mmfftype = 46;
                          goto L_10;
                      }
                  }else if (atom[jatm].atomnum == 16) // N=S=O
                  {
                      for (k=0; k < MAXIAT; k++)
                      {
                          if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9 && atom[jatm].iat[k] != i)
                          {
                              if (atom[atom[jatm].iat[k]].atomnum == 8 && atom[jatm].bo[k] == 2)
                              {
                                    mmfftype = 48;
                                    goto L_10;
                              }
                          }
                      }
                  } else if (atom[jatm].atomnum == 7) // N=N
                  {
                      mmfftype = 9;
                      if (jjbo == 4)
                      {
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[i].iat[k] != 0 && atom[i].iat[k] != jatm)
                              {
                                  if (atom[atom[i].iat[k]].atomnum == 8) // n=n-o  azoxy
                                  {
                                      mmxtype = 37;
                                      mm3type = 37;
                                      mmfftype = 67;
                                      goto L_10;
                                  }
                              }
                          }
                      }
                      goto L_10;                     
                  } else if (atom[jatm].atomnum == 6) // N=C
                  {
                      if (jjbo == 4) // c=n+
                      {
                          mmxtype = 41;
                          for (j=0; j < MAXIAT; j++)
                          {
                              if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                              {
                                  if (atom[atom[i].iat[j]].atomnum == 8) // c=n+-o
                                  {
                                      katm = atom[i].iat[j];
                                      jjk = 0;
                                      for (k=0; k < MAXIAT; k++)
                                      {
                                          if (atom[katm].iat[k] != 0 && atom[katm].bo[k] != 9 && atom[atom[katm].iat[k]].mmx_type != 20)
                                             jjk++;
                                      }
                                      if (jjk == 1)
                                        mmfftype = 67;
                                      else if (jjk == 2)
                                        mmfftype = 54;
                                      goto L_10;
                                  }
                              }
                          }
                          mmfftype = 54;
                          goto L_10;
                      }
//                    n=c
                      jjk = 0;
                      for (k=0; k < MAXIAT; k++)
                      {
                          if (atom[jatm].iat[k] != 0)
                             jjk++;
                      }
                      if (jjk == 3)
                      {
                          if (atom[jatm].iat[0] == i)
                          {
                              katm = atom[jatm].iat[1];
                              latm = atom[jatm].iat[2];
                          } else if (atom[jatm].iat[1] == i)
                          {
                              katm = atom[jatm].iat[0];
                              latm = atom[jatm].iat[2];
                          } else if (atom[jatm].iat[2] == i)
                          {
                              katm = atom[jatm].iat[0];
                              latm = atom[jatm].iat[1];
                          }
                          jbo = 0;
                          kbo = 0;
                          lbo = 0;
                          adjn = FALSE;
                          if (atom[katm].atomnum == 7 && atom[latm].atomnum == 7)
                          {
                             for (k=0; k < MAXIAT; k++)
                             {
                                 if (atom[i].iat[k] != 0 && atom[i].bo[k] != 0 && atom[atom[i].iat[k]].mmx_type != 20)
                                    jbo += atom[i].bo[k];
                                 if (atom[katm].iat[k] != 0 && atom[katm].bo[k] != 0 && atom[atom[katm].iat[k]].mmx_type != 20)
                                 {
                                    kbo += atom[katm].bo[k];
                                    if (atom[katm].bo[k] == 2)
                                       adjn = TRUE;
                                 }
                                 if (atom[latm].iat[k] != 0 && atom[latm].bo[k] != 0 && atom[atom[latm].iat[k]].mmx_type != 20)
                                 {
                                    lbo += atom[latm].bo[k];
                                    if (atom[latm].bo[k] == 2)
                                       adjn = TRUE;
                                 }
                             }
                             if (jbo == 4 && kbo == 3 && lbo == 3 && adjn == FALSE)
                             {
                                mmfftype = 56;
                                goto L_10;
                             }
                          }
                          if (atom[katm].atomnum == 7)
                          {
                              jbo = 0;
                              kbo = 0;
                              adjn = FALSE;
                              for (k=0; k < MAXIAT; k++)
                              {
                                 if (atom[i].iat[k] != 0 && atom[i].bo[k] != 0 && atom[atom[i].iat[k]].mmx_type != 20)
                                    jbo += atom[i].bo[k];
                                 if (atom[katm].iat[k] != 0 && atom[katm].bo[k] != 0 && atom[atom[katm].iat[k]].mmx_type != 20)
                                 {
                                    kbo += atom[katm].bo[k];
                                    if (atom[katm].bo[k] == 2)
                                       adjn = TRUE;
                                 }
                              }
                              if (jbo == 4 && kbo == 3 && adjn == FALSE)
                              {
                                  mmfftype = 55;
                                  goto L_10;
                              }
                          }
                          if (atom[latm].atomnum == 7)
                          {
                              jbo = 0;
                              lbo = 0;
                              adjn = FALSE;
                              for (k=0; k < MAXIAT; k++)
                              {
                                 if (atom[i].iat[k] != 0 && atom[i].bo[k] != 0)
                                    jbo += atom[i].bo[k];
                                 if (atom[latm].iat[k] != 0 && atom[latm].bo[k] != 0)
                                 {
                                    lbo += atom[latm].bo[k];
                                    if (atom[latm].bo[k] == 2)
                                       adjn = TRUE;
                                 }
                              }
                              if (jbo == 4 && lbo == 3 && adjn == FALSE)
                              {
                                  mmfftype = 55;
                                  goto L_10;
                              }
                          }         
                      }
                      mmxtype = 37;
                      mm3type = 37;
                      mmfftype = 9;
                      goto L_10;                 
                  } // got here with not assignment check other side for mmff
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[i].iat[j] != jatm && atom[i].bo[j] != 9)
                      {
                          katm = atom[i].iat[j];
                          if (atom[katm].atomnum == 16)
                          {
                              for (k=0; k < MAXIAT; k++)
                              {
                                  if (atom[katm].iat[k] != 0 && atom[katm].iat[k] != i && atom[katm].bo[k] == 2)
                                  {
                                      if (atom[atom[katm].iat[k]].atomnum == 8)
                                      {
                                          mmfftype = 43;
                                          goto L_10;
                                      }
                                  }
                              }
                          }
                      }
                  }
                  goto L_10;
              }
        // goto here with only single bonds to nitrogen
               mmxtype = 8;
               mm3type = 8;
               mmfftype = 8;
               if (jji == 2)     // divalent anionic N
               {
                   for (j=0; j < MAXIAT; j++)
                   {
                       if (atom[atom[i].iat[j]].mmx_type > 300) // attached metal
                       {
                           for (k=0; k < MAXIAT; k++)
                           {
                               if (atom[i].iat[k] != 0 && atom[i].iat[k] != atom[i].iat[j])
                               {
                                   jatm = atom[i].iat[k];
                                   for (l=0; l < MAXIAT; l++)
                                   {
                                       if (atom[jatm].iat[l] != 0)
                                       {
                                           if (atom[jatm].bo[l] >= 2)
                                           {
                                               mmxtype = 9;
                                               goto L_10;
                                           }
                                       }
                                   }
                               }
                           }
                       }
                   }
//                   mmxtype = 77;
                   mmfftype = 62;
                   goto L_10;
               }
               if (jji == 4 && jjbo == 4)
               {
                   mmxtype = 41;
                   mm3type = 39;
                   mmfftype = 34;
                   oplstype = 56;
                   for (j=0; j < MAXIAT; j++)
                   {
                       if (atom[i].iat[j] != 0 && (atom[atom[i].iat[j]].mmx_type >= 300))
                       {
                           mmxtype = 8;
                           mm3type = 8;
                       }
                       if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                       {
                           if (atom[atom[i].iat[j]].atomnum == 8)
                              mmfftype = 68;
                       }
                   }
                   goto L_10;
               }
               ndouble = 0;
               latm = 0;
               for (j=0; j < MAXIAT; j++)
               {
                   if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                   {
                       katm = atom[i].iat[j];
                       for (k=0; k < MAXIAT; k++)
                       {
                           if (atom[katm].iat[k] != 0 && atom[katm].bo[k] == 2)
                           {
                               if (ndouble == 0)
                               {
                                  ndouble++;
                                  latm = katm;  // atom with double bond
                               } else
                               {
                                   // check atomic number and chose higher atomic number
                                   if (atom[latm].atomnum < atom[katm].atomnum)
                                      latm = katm;
                                   ndouble++;
                               }
                           }
                           if (atom[katm].iat[k] != 0 && atom[katm].bo[k] == 3)
                           {
                              ndouble++;
                              latm = katm;
                           }
                       }
                   }
               }
               if (ndouble >= 1)
               {
                   mmxtype = 9;
                   mm3type = 9;
                   mmfftype = 10;
                   if (atom[latm].atomnum == 6) // n-c=
                   {
                       jjk = 0;
                       for (k=0; k < MAXIAT; k++)
                       {
                           if (atom[latm].iat[k] != 0 && atom[latm].iat[k] != i && atom[latm].bo[k] == 2)
                           {
                               if (atom[atom[latm].iat[k]].atomnum == 8)  // amide
                               {
                                   mmxtype = 9;
                                   mm3type = 9;
                                   mmfftype = 10;
                                   goto L_10;
                               } else if (atom[atom[latm].iat[k]].atomnum == 6)  // eneamine
                               {
                                   mmxtype = 9;
                                   mmfftype = 40;
                                   goto L_10;
                               } else if (atom[atom[latm].iat[k]].atomnum == 15)  // n-p=
                               {
                                   mmxtype = 9;
                                   mm3type = 9;
                                   mmfftype = 10;
                                   goto L_10;
                               } else if (atom[atom[latm].iat[k]].atomnum == 16)  // thioamide
                               {
                                   mmxtype = 9;
                                   mm3type = 9;
                                   mmfftype = 10;
                                   goto L_10;
                               } else if (atom[atom[latm].iat[k]].atomnum == 7)  // eneamine, guanadine
                               {
                                   matm = atom[latm].iat[k];
                                   jjk++;
                               }
                           }
                       }
                       if (jjk == 1)  // n-c=n
                       {
                           jjbo = 0;
                           for (k=0; k < MAXIAT; k++)
                           {
                               if (atom[matm].iat[k] != 0 && atom[atom[matm].iat[k]].mmx_type != 20)
                                 jjbo += atom[matm].bo[k];
                           }
                           if (jjbo >= 4)
                              mmfftype = 55;
                           else
                              mmfftype = 10;
                           goto L_10;
                       } else if (jjk == 2) // n-c=n(n)
                       {
                           mmfftype = 56;
                           goto L_10;
                       }
                   } else if (atom[latm].atomnum == 7) // n-n=
                   {
                       mmfftype = 10;
                       goto L_10;
                   } else if (atom[latm].atomnum == 16) // n-s=
                   {
                      jjk = 0;
                      for (k=0; k < MAXIAT; k++)
                      {
                         if (atom[latm].iat[k] != 0 && atom[latm].bo[k] != 9)
                         {
                            if (atom[atom[latm].iat[k]].atomnum == 8 && atom[latm].bo[k] == 2)
                                     jjk++;
                               }
                           }
                           if (jjk >= 2)
                           {
                              mmxtype = 9;
                              mm3type = 155;
                              mmfftype = 43;
                              goto L_10;
                           }
                   } else if (atom[latm].atomnum == 5) // n-b=
                   {
                      mmxtype = 9;
                      mm3type = 9;
                      mmfftype = 9;
                      goto L_10;
                   } else if (atom[latm].atomnum == 15) // n-p=
                   {
                       jjk = 0;
                       for (l=0; l < MAXIAT; l++)
                       {
                          if (atom[latm].iat[l] != 0 && atom[latm].bo[l] != 9)
                             jjk++;
                       }
                /*       if (jjk > 3)
                       {
                          mmxtype = 8;
                          mm3type = 8;
                          mmfftype = 8;
                          goto L_10;
                       } */
                   }
               }
//
               for (j=0; j < MAXIAT; j++)
               {
                   if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                   {             
                       if (atom[atom[i].iat[j]].atomnum == 6) // carbon
                       {
                           for (k=0; k < MAXIAT; k++)
                           {
                               if (atom[atom[i].iat[j]].iat[k] != 0 && atom[atom[i].iat[j]].iat[k] != i)
                               {
                                   if (atom[atom[i].iat[j]].bo[k] == 3 && atom[atom[atom[i].iat[j]].iat[k]].atomnum == 7)
                                   {
                                       mmfftype = 43;
                                       goto L_10;
                                   }
                               }
                           }
                       } else if (atom[atom[i].iat[j]].atomnum == 7) // nitrogen - n-n bond  - n-n-c=o
                       {
                           jatm = atom[i].iat[j];
                           for (k=0; k < MAXIAT; k++)
                           {
                               if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i)
                               {
                                   katm = atom[jatm].iat[k];
                                   if (atom[katm].atomnum == 6)
                                   {
                                       for (l=0; l < MAXIAT; l++)
                                       {
                                           if (atom[katm].iat[l] != 0 && atom[katm].iat[l] != jatm)
                                           {
                                               if (atom[atom[katm].iat[l]].atomnum == 8 && atom[katm].bo[l] == 2)
                                               {
                                                   mmxtype = 9;
                                                   goto L_10;
                                               }
                                           }
                                       }
                                   }
                               }
                           }
                       }
                   }
               }
               if (is_cyclo5(i,array) && atom[i].flags & aromatic_mask)  // heterocycle
               {
                   for (j=1; j < 5; j++)
                   {
                       if (atom[array[j]].atomnum != 6)
                       {
                          mmxtype = 8;
                          mm3type = 8;
                          mmfftype = 79;
                          goto L_10;
                       }
                   }
               } 
               for (j=0; j < MAXIAT; j++)  // guanadinium
               {
                   if (atom[i].iat[j] != 0 && atom[atom[i].iat[j]].atomnum == 6)
                   {
                       jatm = atom[i].iat[j];
                       jjk = 0;
                       for (k=0; k < MAXIAT; k++)
                       {
                           if (atom[jatm].iat[k] != 0)
                             jjk++;
                       }
                       if (jjk == 3 && (atom[atom[jatm].iat[0]].atomnum == 7 &&
                             atom[atom[jatm].iat[1]].atomnum == 7 && atom[atom[jatm].iat[2]].atomnum == 7) )
                       {
                           mmfftype = 56;
                           oplstype = 57;
                           goto L_10;
                       }
                   }
               }
              //  mmxtypes that should not be changed - set by user
              if (atom[i].mmx_type == 41)  // N+
              {
                  mmxtype = 41;
                  mm3type = 39;
                  mmfftype = atom[i].mmff_type;
                  if (mmfftype == 0 &&  (mmfftype != 34 || mmfftype != 54 || mmfftype != 55 || mmfftype != 56
                       || mmfftype != 58 || mmfftype != 69 || mmfftype != 80 || mmfftype != 81) )
                     mmfftype = 34;
                  goto L_10;
              } else if (atom[i].mmx_type == 55)  // TS Nitrogen
              {
                  mmxtype = 55;
                  mm3type = 0;
                  mmfftype = 0;
                  goto L_10;
              }
               goto L_10;           
          } else if (atom[i].atomnum == 8) // oxygen
          {
/* ========================  Oxygen  ==================================== */
              if (atom[i].mmx_type == 53)  // TS Oxygen
              {
                  mmxtype = 53;
                  mm3type = 0;
                  mmfftype = 0;
                  ambertype = 0;
                  goto L_10;
              }
              jjbo = 0;
              for (j=0; j < MAXIAT; j++)
              {
                  if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].mmx_type != 20
                     && (atom[atom[i].iat[j]].mmx_type < 300))
                    jjbo += atom[i].bo[j];
              }
              if (jjbo == 3)  // o+
              {
                  mmxtype = 46;
                  mmfftype = 49;
                  ambertype = 0;
                  mm3type = 0;
                  oplstype = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[i].bo[j] == 2)
                      {
                          mmfftype = 51;
                          goto L_10;
                      }
                  }
                  goto L_10;
              }
              mmxtype = atom[i].mmx_type;
              if (mmxtype == 66)
              {
                mm3type = 47;
                mmfftype = 32;
                ambertype = 25;
                oplstype = 50;
              }
              if (mmxtype == 6 || mmxtype == 7 || mmxtype == 66)
              {
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].bo[j] == 2)  // x=0
                      {
                          mmxtype = 7;
                          mm3type = 7;
                          mmfftype = 7;
                          ambertype = 25;
                          oplstype = 49;
                          jatm = atom[i].iat[j];
                          if (atom[jatm].atomnum == 15)   // p=o
                          {
                              for (k=0; k < MAXIAT; k++)
                              {
                                  jjk = 0;
                                  if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i && atom[atom[jatm].iat[k]].atomnum == 8)
                                  {
                                      katm = atom[jatm].iat[k];
                                      for (l=0; l < MAXIAT; l++)
                                      {
                                          if (atom[katm].iat[l] != 0 && atom[katm].bo[l] != 9 && atom[atom[katm].iat[l]].mmx_type != 20)
                                              jjk++;
                                      }
                                      if (jjk == 1 && atom[katm].mmx_type == 42)
									  {
                                         mmxtype = 66;
										 mm3type = 47;
									  }
                                  }
                              }
                              mmfftype = 32;
                              oplstype = 50;
                              goto L_10;
                          }
                          if (atom[jatm].atomnum == 16)   // s=o
                          {
                              for (k=0; k < MAXIAT; k++)
                              {
                                  jjk = 0;
                                  jj_bo = 0;
                                  if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i && atom[atom[jatm].iat[k]].atomnum == 8)
                                  {
                                      katm = atom[jatm].iat[k];
                                      for (l=0; l < MAXIAT; l++)
                                      {
                                          if (atom[katm].iat[l] != 0 && atom[katm].bo[l] != 9 && atom[atom[katm].iat[l]].mmx_type != 20)
                                          {
                                              jjk++;
                                              jj_bo += atom[katm].bo[l];
                                          }
                                      }
                                      if (jjk == 1 && jj_bo == 1 && atom[katm].mmx_type == 42)
                                      {
                                         mmxtype = 66;
										 mm3type = 47;
                                         mmfftype = 32;
                                         oplstype = 50;
                                         goto L_10;
                                      }
                                  }
                              }
                          }
                          if (atom[jatm].atomnum == 6)    // c=o
                          {
                              mmfftype = 7;
                              for (k=0; k < MAXIAT; k++)
                              {
                                  if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i && atom[jatm].bo[k] != 9)
                                  {
                                      if (atom[atom[jatm].iat[k]].atomnum == 8)
                                      {
                                          katm = atom[jatm].iat[k];
                                          jjk = 0;
                                          for (l=0; l < MAXIAT; l++)
                                          {
                                              if (atom[katm].iat[l] != 0 && atom[katm].bo[l] != 9 && atom[atom[katm].iat[l]].mmx_type != 20)
                                                jjk++;
                                          }
                                          if (jjk == 1)   // carboxylate
                                          {
                                              mmxtype = 66;
											  mm3type = 47;
                                              mmfftype = 32;
                                              ambertype = 25;
                                              oplstype = 50;
                                              goto L_10;
                                          }
                                          if (jjk == 2)  // ester and acids
                                          {
                                              if ( atom[atom[katm].iat[0]].atomnum != 1 && atom[atom[katm].iat[1]].atomnum != 1) // ester
                                              {
                                                  mm3type = 78;
                                                  ambertype = 25;
                                                  oplstype = 49;
                                                  goto L_10;
                                              } else if ( atom[atom[katm].iat[0]].atomnum == 1 || atom[atom[katm].iat[1]].atomnum == 1) // acid
                                              {
                                                  mm3type = 77;
                                                  ambertype = 24;
                                                  oplstype = 49;
                                                  goto L_10;
                                              }
                                          }
                                      }
                                  }
                              }
                              for (k=0; k < MAXIAT; k++)
                              {
                                  if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i && atom[jatm].bo[k] != 9)
                                  {
                                      if (atom[atom[jatm].iat[k]].atomnum == 7)  // amide
                                      {
                                          ambertype = 24;
                                          oplstype = 49;
                                          mm3type = 7;  // should be 79 but database only has 58-79
                                          goto L_10;
                                      }
                                  }
                              }
                              for (k=0; k < MAXIAT; k++)
                              {
                                  if (atom[jatm].iat[k] != 0 && atom[jatm].iat[k] != i && atom[jatm].bo[k] != 9)
                                  {
                                      if (atom[atom[jatm].iat[k]].atomnum == 6) // vinyl ketone
                                      {
                                          katm = atom[jatm].iat[k];
                                          for (l=0; l < MAXIAT; l++)
                                          {
                                              if (atom[katm].iat[l] != 0 && atom[katm].bo[l] == 2)
                                              {
                                                  mm3type = 81;
                                                  goto L_10;
                                              }
                                          }
                                      }
                                  }
                              }
                          }           
                                          
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[atom[jatm].iat[k]].atomnum == 8 && atom[jatm].iat[k] != i && atom[jatm].bo[k] == 1)
                              {
                                  jji = 0;
                                  latm = atom[jatm].iat[k];
                                  for (l=0; l < MAXIAT; l++)
                                  {
                                      if (atom[latm].iat[l] != 0 && atom[latm].bo[l] != 9 && atom[atom[latm].iat[l]].mmx_type != 20)
                                        jji++;
                                  }
                                  if (jji == 1)
                                  {
                                      mmxtype = 66;
                                      mm3type =  47;
                                      mmfftype = 32;
                                      ambertype = 24;
                                      oplstype = 50;
                                      goto L_10;
                                  }
                              } else if (atom[atom[jatm].iat[k]].atomnum == 8 && atom[jatm].iat[k] != i && atom[jatm].bo[k] == 2)
                              {
                                  mmfftype = 32;
                                  oplstype = 50;
                                  goto L_10;
                              } else if (atom[atom[jatm].iat[k]].atomnum == 16 && atom[jatm].iat[k] != i && atom[jatm].bo[k] == 2)
                              {
                                  jjk = 0;
                                  latm = atom[jatm].iat[k];
                                  for (l=0; l < MAXIAT; l++)
                                  {
                                      if (atom[latm].iat[l] != 0 && atom[latm].bo[l] != 9 && atom[atom[latm].iat[l]].mmx_type != 20)
                                        jjk++;
                                  }
                                  if (jjk == 1)
                                  {
                                      mmxtype = 66;
                                      mm3type = 47;
                                      oplstype = 50;
                                  }
                                  mmfftype = 32;
                                  goto L_10;
                              } else if (atom[atom[jatm].iat[k]].atomnum == 7 && atom[jatm].iat[k] != i && atom[jatm].bo[k] == 2)
                              {
                                  mmfftype = 32;
                                  goto L_10;
                              }
                          }    
                          goto L_10;
                      } else if (atom[i].bo[j] == 3)
                      {
                          mmxtype = 46;
                          mmfftype = 53;
                          goto L_10;
                      }
                  } // got here with only single bonds - don't reset 66
                  jji = 0;
                  nh = 0;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].mmx_type != 20 &&
                          atom[atom[i].iat[j]].mmx_type < 300)
                          jji++;
                      if (atom[atom[i].iat[j]].atomnum == 1)
                          nh++;
                  }
                  if (jji == 1) // only one bond
                  {
                      oplstype = 52;
                      if (atom[atom[i].iat[0]].atomnum == 6)
                      {
                          jatm = atom[i].iat[0];
                          for (j=0; j < MAXIAT; j++)
                          {
                              if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] == 2)
                              {
                                  
                                  if (atom[atom[jatm].iat[j]].atomnum == 8) // c=c
                                  {
                                     oplstype = 50;
                                     goto L_10;
                                  }
                                  if (atom[atom[jatm].iat[j]].atomnum == 6) // c=c
                                  {
                                     mmfftype = 35;
                                     goto L_10;
                                  }
                              }
                          }
                      }
                      if (atom[atom[i].iat[0]].atomnum == 7)  //  n-o
                      {
                          mmxtype = 66;
                          mm3type = 69;
                          mmfftype = 32;
                          oplstype = 50;
                          goto L_10;
                      }
                      if (atom[atom[i].iat[0]].atomnum == 15 || atom[atom[i].iat[0]].atomnum == 16)  //  p-o
                      {
                          mmxtype = 66;
                          mm3type = 7;
                          mmfftype = 32;
                          oplstype = 50;
                          goto L_10;
                      }
                      
                  }
                  // test for epoxides
                  if (is_ring31(i))
                  {
                      mm3type = 49;
                      mmfftype = 6;
                      goto L_10;
                  } 
                  if (is_ring51(i))
                  {
                     if (atom[i].flags & aromatic_mask)
                     {
                         mmfftype = 59;
                         mm3type = 41;
                         ambertype = 23;
                         oplstype = 47;
                         goto L_10;
                     }
                  }
                  if (jji == 2)
                  {
                      if (nh == 2) // mmff water
                      {
                          mmfftype = 70;
                          ambertype = 21;
                          oplstype = 51;
                          goto L_10;
                      } else if (nh == 1) // ROH and RCOOH
                      {
                          mm3type = 6;
                          mmfftype = 6;
                          ambertype = 22;
                          oplstype = 44;
                          jatm = atom[i].iat[0];
                          for (j=0; j < MAXIAT; j++)
                          {
                              if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] == 2)
                              {
                                  if (atom[atom[jatm].iat[j]].atomnum == 8)
                                  {
                                      mm3type = 75;
                                      oplstype = 46;
                                      break;
                                  }
                              }
                          }
                          goto L_10;
                      } else if (atom[atom[i].iat[0]].atomnum == 15 || atom[atom[i].iat[1]].atomnum == 15)
                      {
                          mm3type = 159;
                          goto L_10;
                      } else
                      {
                          mmfftype = 6;
                          ambertype = 23;   // ether
                          oplstype = 47;
                          jatm = atom[i].iat[0];  // acetal
                          katm = atom[i].iat[1];
                          nox = 0;
                          for (j=0; j < MAXIAT; j++)
                          {
                              if (atom[jatm].iat[j] != 0 && atom[atom[jatm].iat[j]].atomnum == 8 && atom[jatm].bo[j] != 2)
                                 nox++;
                          }
                          if (nox == 2 || nox == 3)
                             oplstype = 48;
                          nox = 0;
                          for (j=0; j < MAXIAT; j++)
                          {
                              if (atom[katm].iat[j] != 0 && atom[atom[katm].iat[j]].atomnum == 8 && atom[katm].bo[j] != 2)
                                 nox++;
                          }
                          if (nox == 2 || nox == 3)
                             oplstype = 48;
                             
                      }
                      jatm = atom[i].iat[0];
                      jtype = 0;
                      katm = atom[i].iat[1];
                      ktype = 0;
                      for (j=0; j < MAXIAT; j++)
                      {
                          if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] == 2)
                          {
                              if (atom[atom[jatm].iat[j]].atomnum == 8)          // C=O
                              {
                                 jtype = 3;
                                 break;
                              } else if (atom[atom[jatm].iat[j]].atomnum == 6)  // C=C
                              {
                                  jtype = 2;
                                  break;
                              } else if (atom[atom[jatm].iat[j]].atomnum == 15) // P=O
                              {
                                  jtype = 4;
                                  break;
                              }
                          }
                      }
                      for (j=0; j < MAXIAT; j++)
                      {
                          if (atom[katm].iat[j] != 0 && atom[katm].bo[j] == 2)
                          {
                              if (atom[atom[katm].iat[j]].atomnum == 8)
                              {
                                 ktype = 3;
                                 break;
                              }else if (atom[atom[katm].iat[j]].atomnum == 6)
                              {
                                  ktype = 2;
                                  break;
                              }else if (atom[atom[katm].iat[j]].atomnum == 15)
                              {
                                  ktype = 4;
                                  break;
                              }
                          }
                      }
                      if (jtype == 3 && ktype == 3) // anhydrides
                         mm3type = 148;
                      else if (jtype == 3 || ktype == 3)  // carboxyl & ester
                         mm3type = 75;
                      else if (jtype == 4 || ktype == 4)  // phosphate ester
                         mm3type = 159;
                      else if (jtype == 2 && ktype == 2) // furan type ??
                         mm3type = 41;
                      else if (jtype == 2 || ktype == 2) // vinyl type ??
                         mm3type = 41;
                      else
                         mm3type = 6;
                      goto L_10;
                  }
                  if (mmxtype != 66)
                  {
                     mmxtype = 6;
                     mm3type = 6;
                     mmfftype = 6;
                     ambertype = 22;
                     oplstype = 47;
                  }
                  goto L_10;
              } else if (mmxtype == 42)
              {
                  mmfftype = 32;
                  ambertype = 25;
                  oplstype = 52;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[atom[i].iat[j]].atomnum == 16 || atom[atom[i].iat[j]].atomnum == 15)
                      {
                          mmxtype = 66;
                          mm3type = 47;
                          mmfftype = 32;
                          oplstype = 50;
                          goto L_10;
                      }
                      if (atom[atom[i].iat[j]].atomnum == 7)
                      {
                          mmxtype = 42;
                          mm3type = 69;
                          jj_bo = 0;
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[atom[i].iat[j]].iat[k] != 0 && atom[atom[i].iat[j]].bo[k] != 9)
                                  jj_bo += atom[atom[i].iat[j]].bo[k];
                          }
                          if (jj_bo >= 4)
                          {
                              mmxtype = 66;
                              mm3type = 69;
                              mmfftype = 32;
                              oplstype = 50;
                          } else
                            mmfftype = 35;
                         goto L_10; 
                      }
                      if (atom[atom[i].iat[j]].atomnum == 6 )
                      {
                          jatm = atom[i].iat[j];
                          mmfftype = 35;
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[jatm].iat[k] != i && atom[atom[jatm].iat[k]].atomnum == 8 && atom[jatm].bo[k] == 2)
                              {
                                  mmxtype = 66;
                                  mm3type = 47;
                                  mmfftype = 32;
                                  oplstype = 50;
                                  goto L_10;
                              }
                          }
                      }
                  }
                  goto L_10;
              } else
                 goto L_10;
          } else if (atom[i].atomnum == 13) // aluminum
          {
              mmxtype = 44;
              mm3type = 0;
              jji = 0;
              for (j=0; j < MAXIAT; j++)
                   if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                      jji++;
              if (jji == 4)
                 mmxtype = 58;
              goto L_10;
          } else if (atom[i].atomnum == 15) // phosphorus
          {
              mmxtype = 25;
              mm3type = 25;
              mmfftype = 25;
              ambertype = 28;
              oplstype = 66;
              jji = 0;
              for (j=0; j < MAXIAT; j++)
              {
                   if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                      jji++;
                   if (atom[atom[i].iat[j]].atomnum == 6 && atom[i].bo[j] == 2)
                   {
                      mmxtype = 67;
                      mmfftype = 75;
                      oplstype = 67;
                      goto L_10;
                   }
                   if (atom[atom[i].iat[j]].atomnum == 8 && atom[i].bo[j] == 2)
                      mm3type = 153;
              }
              if (jji >= 5)
              {
                 mmxtype = 47;
                 mm3type = 60;
                 oplstype = 67;
              }else if (jji == 3)
              {
                  mmxtype = 25;
                  mmfftype = 26;
              }
              goto L_10;
          } else if (atom[i].atomnum == 16) // sulfur
          {
              mmxtype = 15;
              mm3type = 15;
              mmfftype = 15;
              ambertype = 26;
              oplstype = 69;
              if (atom[atom[i].iat[0]].atomnum == 1 || atom[atom[i].iat[1]].atomnum == 1)
              {
                  ambertype = 27;
                  oplstype = 68;
              }
              jji = 0;
              ndouble = 0;
              for (j=0; j < MAXIAT; j++)
              {
                  if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9 && atom[atom[i].iat[j]].mmx_type != 20 &&
                  (atom[atom[i].iat[j]].mmx_type < 300))
                  {
                      jji++;
                      if (atom[i].bo[j] == 2)
                         ndouble++;
                  }
              }
              if (jji == 1 && ndouble == 1)
              {
                  if (atom[atom[i].iat[0]].atomnum == 16)  // s=s-
                  {
                      mmfftype = 72;
                      goto L_10;
                  }
                  if (atom[atom[i].iat[0]].atomnum == 6)  // s=c
                  {
                      mmxtype = 38;
                      mm3type = 74;
                      if (atom[i].bo[0] == 2)
                      {
                         mmfftype = 16;
                         jatm = atom[i].iat[0];
                          for (k=0; k < MAXIAT; k++)
                          {
                              if (atom[jatm].iat[k] != 0 && atom[jatm].bo[k] != 9)
                              {
                                  if (atom[atom[jatm].iat[k]].atomnum == 16 && atom[jatm].iat[k] != i)
                                  {
                                      katm = atom[jatm].iat[k];
                                      jj_bo = 0;
                                      for (l=0; l < MAXIAT; l++)
                                      {
                                          if (atom[katm].iat[l] != 0 && atom[katm].bo[l] != 9)
                                              jj_bo++;
                                      }
                                      if (jj_bo == 1)  //  s=c-s    i-jatm-katm
                                      {
                                          mmfftype = 72;
                                          goto L_10;
                                      }
                                  }
                              }
                          }
                          goto L_10;
                      }else if (atom[i].bo[0] == 1)
                      {
                          mmfftype = 72;
                          goto L_10;
                      }
                      goto L_10;
                  }
                  if (atom[atom[i].iat[0]].atomnum == 15 ) // s=p
                  {
                      mmxtype = 38;
                      mm3type = 74;
                      mmfftype = 72;
                      goto L_10;
                  }                  
              }
              if (jji == 2)
              {
                  mmxtype = 15;
                  mm3type = 15;
                  mmfftype = 15;
                  if (is_cyclo5(i,array))  // thiophene
                  {
                      if (atom[i].flags & aromatic_mask)
                      {
                          mmfftype = 44;
                          mm3type = 42;
                          goto L_10;
                      }
                  }
                  if (ndouble == 2)
                      mmfftype = 74;
                  goto L_10;
              }
              if (jji == 3)
              {
                  mmxtype = 17;
                  mm3type = 17;
                  mmfftype = 17;
                  if (ndouble == 2)
                    mmfftype = 73;
                  if (ndouble == 3)
                    mmfftype = 18;
                  goto L_10;
              }
              if (jji == 4)
              {
                  mmxtype = 18;
                  mm3type = 18;
                  for (j=0; j < MAXIAT; j++)
                  {
                      if (atom[atom[i].iat[j]].atomnum == 7 ) // sulfamide
                         mm3type = 154;
                  }
                  mmfftype = 18;
                  goto L_10;
              }
              goto L_10;
          } else if (atom[i].atomnum == 17) // chlorine
          {
              mm3type = 12;
              mmxtype = 12;
              mmfftype = 12;
              oplstype = 63;
              jji = 0;
              for (j=0; j < MAXIAT; j++)
                   if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                      jji++;
              if (jji == 2)
                 mmxtype = 74;  // bridging chlorine
              else if (jji == 4) // perchlorate
                 mmfftype = 77;
              else if (jji == 0)
                 mmfftype = 90;
              goto L_10;         
          } else if (atom[i].atomnum == 34) // selenium
          {
              mm3type = 34;
              mmxtype = 34;
              for (j=0; j < MAXIAT; j++)
              {
                  if (atom[i].bo[j] == 2)
                     mmxtype = 39;
              }
              goto L_10;
          } else if (atom[i].atomnum == 9)  // Florine
          {
              mmxtype = 11;
              mm3type = 11;
              mmfftype = 11;
              oplstype = 62;
              jji = 0;
              for (j=0; j < MAXIAT; j++)
                 if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                      jji++;
              if (jji == 0)
                mmfftype = 89;
              goto L_10;
          } else if (atom[i].atomnum == 35)   // Bromine
          {
              mmxtype = 13;
              mm3type = 13;
              mmfftype = 13;
              oplstype = 64;
              jji = 0;
              for (j=0; j < MAXIAT; j++)
                 if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                      jji++;
              if (jji == 0)
                mmfftype = 91;
              goto L_10;
          } else if (atom[i].atomnum == 53)    // Iodine
          {
              if (atom[i].mmx_type == 54)  // Sn2 I 
              {
                 mmxtype = 54;
                 mm3type = 0;
                 mmfftype = 0;
                 goto L_10;
              }
              mmxtype = 14;
              mm3type = 14;
              mmfftype = 14;
              oplstype = 65;
              goto L_10;
          } else if (atom[i].atomnum == 14 )  // Silicon
          {
              mmxtype = 19;
              mm3type = 19;
              mmfftype = 19;
              goto L_10;
          } else if (atom[i].atomnum == 50 ) // Tin
          {
              mmxtype = 32;
              mm3type = 32;
              mmfftype = 0;
              goto L_10;
          } else if (mmxtype >= 300) // metal atom - try to assign MMFF type
          {
              mm3type = mmxtype;
              ambertype = mmxtype;
              mmfftype = mmxtype;
              oplstype = mmxtype;
          }
L_10:
          set_atomtype(i,mmxtype,mm3type,mmfftype,ambertype,oplstype);
L_20:
          continue;  // do nothing
      }
   }
   // finished normal now type amber
   set_atomtypes(field.type);
}
/* --------------------------------------------------------- */
int is_cyclo5(int ia, int *array)
{
    int i,j,k,l;
    int jatm,katm,latm;

    for (i=0; i < 6; i++)
       array[i] = 0;
       
    for(i=0; i < MAXIAT; i++)
    {
        if (atom[ia].iat[i] != 0)
        {
            jatm = atom[ia].iat[i];
            for(j=0; j < MAXIAT; j++)
            {
                if (atom[jatm].iat[j] != 0 && atom[jatm].iat[j] != ia && atom[jatm].bo[j] != 9)
                {
                    katm = atom[jatm].iat[j];
                    if ( !isbond(ia,katm))  // three membered ring
                    {
                      for(k=0; k < MAXIAT; k++)
                      {
                        if (atom[katm].iat[k] != 0 && atom[katm].iat[k] != jatm && atom[katm].iat[k] != ia && atom[katm].bo[k] != 9)
                        {
                            latm = atom[katm].iat[k];
                            if ( !isbond(ia,latm))  // four membered ring
                            {
                              for (l=0; l < MAXIAT; l++)
                              {
                                if (atom[latm].iat[l] != 0 && atom[latm].iat[l] != katm && atom[latm].iat[l] != jatm &&
                                    atom[latm].iat[l] != ia && atom[latm].bo[l] != 9)
                                    {
                                        if (isbond(ia,atom[latm].iat[l]))
                                        {
                                            array[0] = ia;
                                            array[1] = jatm;
                                            array[2] = katm;
                                            array[3] = latm;
                                            array[4] = atom[latm].iat[l];
                                            return TRUE;
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
    return FALSE;
}
/* --------------------------------------------------------- */                                           
int is_cyclo6(int ia, int *array)
{
    int i,j,k,l,m;
    int jatm,katm,latm,matm,xatm;

    for (i=0; i < 6; i++)
       array[i] = 0;
       
    for(i=0; i < MAXIAT; i++)
    {
        if (atom[ia].iat[i] != 0)
        {
            jatm = atom[ia].iat[i];
            for(j=0; j < MAXIAT; j++)
            {
                if (atom[jatm].iat[j] != 0 && atom[jatm].iat[j] != ia && atom[jatm].bo[j] != 9)
                {
                    katm = atom[jatm].iat[j];
                    for(k=0; k < MAXIAT; k++)
                    {
                        if (atom[katm].iat[k] != 0 && atom[katm].iat[k] != jatm && atom[katm].iat[k] != ia && atom[katm].bo[k] != 9)
                        {
                            latm = atom[katm].iat[k];
                            if (!isbond(ia,latm))  // bail out for four memebered rings
                            {
                              for (l=0; l < MAXIAT; l++)
                              {
                                if (atom[latm].iat[l] != 0 && atom[latm].iat[l] != katm && atom[latm].iat[l] != jatm &&
                                    atom[latm].iat[l] != ia && atom[latm].bo[l] != 9)
                                {
                                    matm = atom[latm].iat[l];
                                    if (!(isbond(ia,matm)))  // bail out for five membered ring
                                    {
                                      for (m=0; m < MAXIAT; m++)
                                      {
                                        xatm = atom[matm].iat[m];
                                        if (xatm != 0 && xatm != matm && xatm != latm && xatm != katm && xatm != jatm && xatm != ia)
                                        { 
                                            if (isbond(ia,xatm))
                                            {
                                               array[0] = ia;
                                               array[1] = jatm;
                                               array[2] = katm;
                                               array[3] = latm;
                                               array[4] = matm;
                                               array[5] = xatm;
                                               return TRUE;
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
        }
    }
    return FALSE;
}
/* --------------------------------------------------------- */
/* ==================================================  */
// look for ionic types input with bonds and adjust to make ionic
//
void adjust_mmfftypes()
{
    int i,j, iatt;

    for (i=1; i <= natom; i++)
    {
        if (atom[i].mmff_type == 89 ) // F-
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                    deletebond(i,atom[i].iat[j]);
            }
        } else if (atom[i].mmff_type == 90 ) // CL-
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                    deletebond(i,atom[i].iat[j]);
            }
        } else if (atom[i].mmff_type == 91 ) // Br-
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                    deletebond(i,atom[i].iat[j]);
            }
        } else if (atom[i].mmff_type == 92 ) // Li+
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                {
                    iatt = atom[i].iat[j];
                    deletebond(i,iatt);
                    if (atom[iatt].atomnum == 8)
                       set_atomtype(iatt,42,0,92,0,0);
                }
            }
        } else if (atom[i].mmff_type == 93 ) // Na+
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                {
                    iatt = atom[i].iat[j];
                    deletebond(i,iatt);
                    if (atom[iatt].atomnum == 8)
                       set_atomtype(iatt,42,0,93,0,0);
                }
            }
        } else if (atom[i].mmff_type == 94 ) // K+
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                {
                    iatt = atom[i].iat[j];
                    deletebond(i,iatt);
                    if (atom[iatt].atomnum == 8)
                       set_atomtype(iatt,42,0,94,0,0);
                }
            }
        } else if (atom[i].mmff_type == 95 ) // Zn+2
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                {
                    iatt = atom[i].iat[j];
                    deletebond(i,iatt);
                    if (atom[iatt].atomnum == 8)
                       set_atomtype(iatt,42,0,95,0,0);
                }
            }
        } else if (atom[i].mmff_type == 96 ) // Ca+2
        {
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0)
                {
                    iatt = atom[i].iat[j];
                    deletebond(i,iatt);
                    if (atom[iatt].atomnum == 8)
                       set_atomtype(iatt,42,0,96,0,0);
                }
            }
        }
    }
}    
