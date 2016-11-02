#include "pcwin.h"
#include "pcmod.h"
#include "torsions.h"
#include "minim_control.h"
#include "tork.h"


int is_allene(int i, int j, int k, int l)
{
    if (atom[i].mmx_type == 2 && atom[j].mmx_type == 4 && atom[k].mmx_type == 2)
       return FALSE;
    if (atom[j].mmx_type == 2 && atom[k].mmx_type == 4 && atom[l].mmx_type == 2)
       return FALSE;
    return TRUE;
}
// ==================================================
int is_linear(int ia)
{
    if (atom[ia].mmx_type == 4 || atom[ia].mmx_type == 10)
       return TRUE;
    if (atom[ia].mmff_type == 53)
       return TRUE;
    return FALSE;
}       
//   ==================================      
int istorsion(int ia1, int ia2)
{
    int i;

    for (i=0; i < torsions.ntor; i++)
    {
        if (torsions.i14[i][0] == ia1 && torsions.i14[i][3] == ia2)
          return TRUE;
        if (torsions.i14[i][0] == ia2 && torsions.i14[i][3] == ia1)
          return TRUE;
    }
    return FALSE;
}
//   ==================================         
void get_torsions()
{
    int i, j, k, l,m;
    int ia1, ia2, nRc;
    int katm, latm;
    
    i = j = k = l = m = 0;
    torsions.ntor = 0;
    allene.nallene = 0;
    for (i=1; i <=natom; i++)
    {
        ia1 = i;
        for(j=0; j< MAXIAT; j++)
        {
            if (atom[i].iat[j] != 0 && isbond(i,atom[i].iat[j]) && atom[atom[i].iat[j]].type < 300 )
            {
                ia2 = atom[i].iat[j];
                if (is_linear(ia2) == FALSE)
                {
                  for(k=0; k<MAXIAT; k++)
                  {
                    if (atom[ia2].iat[k] != 0 && isbond(atom[i].iat[j],atom[atom[i].iat[j]].iat[k]) && atom[atom[i].iat[j]].type < 300 )
                    {
                        katm = atom[atom[i].iat[j]].iat[k];
                    if (isbond(ia2,katm) && atom[katm].type < 300 )
                    {
                        if (i != katm)
                        {
                            for(l=0; l<MAXIAT; l++)
                            {
                                if (isbond(katm,atom[katm].iat[l]) && !is_linear(katm) )
                                {
                                  latm = atom[katm].iat[l];
                                  nRc = is_allene(i,atom[i].iat[j],katm,latm);
                                  if ( (latm != atom[i].iat[j]) && (latm >i) && nRc)
                                  {
                                      torsions.i14[torsions.ntor][0] = i;
                                      torsions.i14[torsions.ntor][1] = atom[i].iat[j];
                                      torsions.i14[torsions.ntor][2] = katm;
                                      torsions.i14[torsions.ntor][3] = latm;
                                      torsions.ntor++;
                                      if (torsions.ntor > MAXTOR)
                                      {
                                          message_alert("Too many torsions. PCMODEL will now exit","Torsion Setup");
                                          fprintf(pcmoutfile,"Too many torsions in torsion setup.\n");
                                          exit(0);
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
    // look for allenes
    for (i=1; i <= natom; i++)
    {
        if (atom[i].mmx_type == 4)
        {
            if (atom[i].bo[0] == 2 && atom[i].bo[1] == 2) // got an cummulene
            {
                if (atom[atom[i].iat[0]].mmx_type == 2 && atom[atom[i].iat[1]].mmx_type == 2) // got an allene
                {
                    katm = atom[i].iat[0];
                    latm = atom[i].iat[1];
                    if (atom[katm].iat[0] == i)
                    {
                        j = atom[katm].iat[1];
                        k = atom[katm].iat[2];
                    }else if (atom[katm].iat[1] == i)
                    {
                        j = atom[katm].iat[0];
                        k = atom[katm].iat[2];
                    }else if (atom[katm].iat[2] == i)
                    {
                        j = atom[katm].iat[1];
                        k = atom[katm].iat[2];
                    }
                    
                    if (atom[latm].iat[0] == i)
                    {
                        l = atom[latm].iat[1];
                        m = atom[latm].iat[2];
                    }else if (atom[katm].iat[1] == i)
                    {
                        l = atom[latm].iat[0];
                        m = atom[latm].iat[2];
                    }else if (atom[katm].iat[2] == i)
                    {
                        l = atom[latm].iat[1];
                        m = atom[latm].iat[2];
                    }
                    torsions.i14[torsions.ntor][0] = j;
                    torsions.i14[torsions.ntor][1] = katm;
                    torsions.i14[torsions.ntor][2] = latm;
                    torsions.i14[torsions.ntor][3] = l;
                    allene.ntor[allene.nallene] = torsions.ntor;
                    allene.nallene++;
                    torsions.ntor++;
                    
                    torsions.i14[torsions.ntor][0] = j;
                    torsions.i14[torsions.ntor][1] = katm;
                    torsions.i14[torsions.ntor][2] = latm;
                    torsions.i14[torsions.ntor][3] = m;
                    allene.ntor[allene.nallene] = torsions.ntor;
                    allene.nallene++;
                    torsions.ntor++;
                    
                    torsions.i14[torsions.ntor][0] = k;
                    torsions.i14[torsions.ntor][1] = katm;
                    torsions.i14[torsions.ntor][2] = latm;
                    torsions.i14[torsions.ntor][3] = l;
                    allene.ntor[allene.nallene] = torsions.ntor;
                    allene.nallene++;
                    torsions.ntor++;
                        
                    torsions.i14[torsions.ntor][0] = k;
                    torsions.i14[torsions.ntor][1] = katm;
                    torsions.i14[torsions.ntor][2] = latm;
                    torsions.i14[torsions.ntor][3] = m;
                    allene.ntor[allene.nallene] = torsions.ntor;
                    allene.nallene++;
                    torsions.ntor++;
                    if (torsions.ntor > MAXTOR)
                    {
                        message_alert("Too many allene torsions. PCMODEL will now exit","Torsion Setup");
                        fprintf(pcmoutfile,"Too many allene torsions in torsion setup.\n");
                    }
                    if (allene.nallene > 10)
                    {
                        message_alert("Too many allene torsions. PCMODEL will now exit","Torsion Setup");
                        fprintf(pcmoutfile,"Too many allene torsions in torsion setup.\n");
                    }
                }
            }
        }
    }                  
}
// ======================================================================
void max_torsions()
{
    int i, j, k, l;
    int ia1, ia2;
    int katm, latm;
    
    torsions.ntor = 0;
    for (i=1; i <=natom; i++)
    {
        ia1 = i;
        for(j=0; j< MAXIAT; j++)
        {
            if (isbond(i,atom[i].iat[j]) )
            {
                ia2 = atom[i].iat[j];
                for(k=0; k<MAXIAT; k++)
                {
                    if (isbond(atom[i].iat[j],atom[atom[i].iat[j]].iat[k]) )
                    {
                        katm = atom[atom[i].iat[j]].iat[k];
                        if (i != katm)
                        {
                            for(l=0; l<MAXIAT; l++)
                            {
                                if (isbond(katm,atom[katm].iat[l]) )
                                {
                                  latm = atom[katm].iat[l];
                                  if ( (latm != atom[i].iat[j]) && (latm >i) )
                                  {
                                      torsions.ntor++;
                                  }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // find allenes
    for (i=0; i <= natom; i++)
    {
        if (atom[i].mmx_type == 4)
        {
            if (atom[i].bo[0] == 2 && atom[i].bo[1] == 2)
            {
                if (atom[atom[i].iat[0]].mmx_type == 2 && atom[atom[i].iat[1]].mmx_type == 2)
                {
                    torsions.ntor += 4;
                }
            }
        }
    }                   
}


