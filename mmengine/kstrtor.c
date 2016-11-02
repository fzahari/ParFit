#include "pcwin.h"
#include "pcmod.h"
#include "torsions.h"
#include "pot.h"
#include "field.h"
#include "crossterm.h"
#include "ebond.h"
       
struct t_strtor {
    int nstrtor,ist[MAXTOR][2];
    float stcon[MAXTOR];
    }   strtor;

void kstrtor()
{
      int i,j,iz, ia,ita, USE;
      int ib,ic,itb,itc;
      char pa[4],pb[4],pt[7];

      strtor.nstrtor = 0;
      for (i = 0; i < torsions.ntor; i++)
      {
          USE = TRUE;
          ib = torsions.i14[i][1];
          ic = torsions.i14[i][2];
          itb = atom[ib].type;
          itc = atom[ic].type;
          if (itb == 2 || itb == 57)
          {
            ia = torsions.i14[i][0];
            ita = atom[ia].type;
            if (ita == 2 || ita == 57 || ita == 3 || ita == 9)
                 USE = FALSE;
          }
          if (itc == 2 || itc == 57)
          {
            ia = torsions.i14[i][3];
            ita = atom[ia].type;
            if (ita == 2 || ita == 57 || ita == 3 || ita == 9)
                 USE = FALSE;
          }
          if (itb == 3)
          {
             ia = torsions.i14[i][0];
             ita = atom[ia].type;
             if (ita == 7) USE = FALSE;
          }
          if (itc == 3)
          {
             ia = torsions.i14[i][3];
             ita = atom[ia].type;
             if (ita == 7) USE = FALSE;
          }
          
          numeral(itb,pa,3);
          numeral(itc,pb,3);
          if (itb < itc )
          {
            strcpy(pt,pa);
            strcat(pt,pb);
          } else
          {
            strcpy(pt,pb);
            strcat(pt,pa);
          }

          for (j=0; j < crossterm_k.nstrtor; j++)
          {
              if (strcmp(pt,crossterm_k.str_tor[j]) == 0 && USE == TRUE)
              {
                strtor.stcon[strtor.nstrtor] = crossterm_k.str_torcon[j];
                strtor.ist[strtor.nstrtor][0] = i;  // torsion number
                iz = find_bond(ib, ic);
                strtor.ist[strtor.nstrtor][1] = iz;
                strtor.nstrtor++;
                break;
              }
          }
      }
      if (strtor.nstrtor == 0)
         pot.use_strtor = FALSE;
}           

