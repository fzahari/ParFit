#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"
#include "energies.h"
#include "torsions.h"
#include "gmmx.h"
#include "pot.h"
#include "field.h"
#include "nonbond.h"
#include "pistuf.h"
#include "minim_values.h"
#include "files.h"
#include "pcmfile.h"
#include "minim_control.h"
#include "dipmom.h"
#include "optimize.h"
#include "eheat.h"
#include "gmmx_run.h"

#include <time.h>
#include <errno.h>

#define NL gettoken(); if (pcmfile.head == 1000) goto L_30;
static char gmmx_comment[64];
static int *jji;
static int *iest, *irest;
static int npibond, nbond, **ibstereo;
int restart_job;
static double *xold,*yold,*zold;

static int n15, nbad15, **k15;
static int const_fail;

int input_type;

int nerr;


void read_gmmxinp()
{
   char filename[90], line[121];
   char dummy1[20],dummy2[20];
   static char szarWorkBuf[36];      
   int i,j, oldfield, sec_cycle,nheav, nhyd, nret;
   FILE *wfile, *infile;

// put up dialog box to get input name
         
       initialize();
       strcpy(filename,Openbox.fname);
     
       wfile = fopen_path(Openbox.path,Openbox.fname,"r");
       fscanf(wfile,"%s\n",gmmx_data.finame);
          strcpy(szarWorkBuf, gmmx_data.finame);

         for (i = 0; i < 36; i++)
         {
                gmmx_data.jobname[i] = szarWorkBuf[i];
                if (szarWorkBuf[i] == '.')
                {
                    gmmx_data.jobname[i] = '\0';
                        break;

                }
         }

       fscanf(wfile,"%s\n",gmmx_data.foname);
       
        sec_cycle = FALSE;
        pcmfile.nocaps = False;
L_30:
        FetchRecord(wfile,pcmfile.string);
        *(pcmfile.string + strlen(pcmfile.string) ) = '\0';
        if (feof(wfile))
                goto L_170;
        pcmfile.head = 1;
L_40:
        NL;
        if (strcmp(pcmfile.token, "METHOD") == 0) 
        {
                NL;
                gmmx_data.method = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "SEED") == 0) 
        {
                NL;
                gmmx_data.iseed = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "FF") == 0) 
        {
                NL;
                minim_control.field = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "LNANT") == 0) 
        {
                NL;
                  gmmx_data.lnant = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "ECUT") == 0) 
        {
                NL;
                  gmmx_data.ecut = atoi(pcmfile.token);
                goto L_40;
        } else if (strcmp(pcmfile.token,"HYBRID") == 0)
        {
                NL;
                gmmx_data.hybrid = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "HEAT") == 0) 
        {
                NL;
                  gmmx_data.heat = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "DOPI") == 0) 
        {
                NL;
                  gmmx_data.nopi = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "CHIG") == 0) 
        {
                NL;
                gmmx_data.chig = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "BAD") == 0) 
        {
                NL;
                NL;
                  gmmx_data.bad15 = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "COMP") == 0) 
        {
                NL;
                  gmmx_data.comp = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "HBOND") == 0) 
        {
                NL;
                  gmmx_data.hbond = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "COMPMET") == 0) 
        {
                NL;
                if (sec_cycle == FALSE)
                  gmmx_data.comp_method = atoi(pcmfile.token);
                else
                  gmmx_data.comp_method2 = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "STOP") == 0) 
        {
                NL;
                  gmmx_data.kstop = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "MIN") == 0) 
        {
                NL;
                  gmmx_data.kmin = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "DUP") == 0) 
        {
                NL;
                  gmmx_data.kdup = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "MAX") == 0) 
        {
                NL;
                  gmmx_data.max_search = atoi(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "EWINDOW") == 0) 
        {
                NL;
                if (sec_cycle == FALSE)
                  gmmx_data.ewindow = atof(pcmfile.token);
                else
                  gmmx_data.ewindow2 = atof(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "BOLTZTEMP") == 0) 
        {
                NL;
                  gmmx_data.boltz = atof(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "MIN_CONTACT") == 0) 
        {
                NL;
                  gmmx_data.min_contact = atof(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "ECUTOFF") == 0) 
        {
                NL;
//                if (sec_cycle == FALSE)
                  gmmx_data.ecutoff = atof(pcmfile.token);
//                else
//                  gmmx_data.ecutoff2 = atof(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "CRMS") == 0) 
        {
                NL;
//                if (sec_cycle == FALSE)
                  gmmx_data.crmsa = atof(pcmfile.token);
//                else
//                  gmmx_data.crmsa2 = atof(pcmfile.token);
                goto L_40;
        }  else if (strcmp(pcmfile.token, "NUMRINGS") == 0) 
        {
                NL;
                gmmx_data.nrings = atoi(pcmfile.token);
                for (i=0; i < gmmx_data.nrings; i++)
                {
                    FetchRecord(wfile,line);
                    sscanf(line,"%s %d\n",dummy1, &gmmx_data.rng_size[i]);
                    for(j=0; j < gmmx_data.rng_size[i]; j++)
                       fscanf(wfile,"%d ", &gmmx_data.ring_atoms[j][i]);
                    FetchRecord(wfile,line);
                    sscanf(line,"%s %d %s %d %d %d %d\n",dummy1, &gmmx_data.ring_resolution[i], dummy2, &gmmx_data.clo_bond[0][i],
                       &gmmx_data.clo_bond[1][i], &gmmx_data.clo_bond[2][i], &gmmx_data.clo_bond[3][i]);

                    FetchRecord(wfile,line);
                    sscanf(line,"%s %f %s %f %f\n",dummy1, &gmmx_data.clo_distance[i], dummy2, &gmmx_data.clo_ang[0][i],
                       &gmmx_data.clo_ang[1][i]);

                    FetchRecord(wfile,line);
                    sscanf(line,"%s %d\n",dummy1, &gmmx_data.nring_bonds[i]);
                    for (j=0; j < gmmx_data.nring_bonds[i]; j++)
                    {
                       FetchRecord(wfile,line);
                       sscanf(line,"%s %d %d %d\n", dummy1, &gmmx_data.ring_bond_data[i][j][0], &gmmx_data.ring_bond_data[i][j][1],
                          &gmmx_data.ring_bond_data[i][j][2]);
                    }
                }
                goto L_30;
        }  else if (strcmp(pcmfile.token, "NUMBONDS") == 0) 
        {
                NL;
                gmmx_data.nbonds = atoi(pcmfile.token);
                for (i=0; i < gmmx_data.nbonds; i++)
                {
                    FetchRecord(wfile,line);
                    sscanf(line,"%s %d %d %d %d %d\n", dummy1, &gmmx_data.bond_data[i][0], &gmmx_data.bond_data[i][1],
                      &gmmx_data.bond_data[i][2],&gmmx_data.bond_data[i][3],&gmmx_data.bond_data[i][4]);

/*Hay debug
                    printf("Number of bonds = %3d\n", gmmx_data.nbonds);
                    printf("ith jth and incr  %3d %3d %3d\n", gmmx_data.bond_data[1][0],
                    gmmx_data.bond_data[1][1], gmmx_data.bond_data[1][2]); */

                }
                goto L_30;
        } else if (strcmp(pcmfile.token,"COMPA") == 0)
        {
            NL;
            gmmx_data.ncompatoms = atoi(pcmfile.token);
            for (i=0; i < gmmx_data.ncompatoms; i += 20)
            {
                    FetchRecord(wfile,line);
                sscanf(line,"%d %d %d %d %d %d %d %d %d %d%d %d %d %d %d%d %d %d %d %d",
                  &gmmx_data.comp_atoms[i],&gmmx_data.comp_atoms[i+1],&gmmx_data.comp_atoms[i+2],&gmmx_data.comp_atoms[i+3],&gmmx_data.comp_atoms[i+4],
                  &gmmx_data.comp_atoms[i+5],&gmmx_data.comp_atoms[i+6],&gmmx_data.comp_atoms[i+7],&gmmx_data.comp_atoms[i+8],&gmmx_data.comp_atoms[i+9],
                  &gmmx_data.comp_atoms[i+10],&gmmx_data.comp_atoms[i+11],&gmmx_data.comp_atoms[i+12],&gmmx_data.comp_atoms[i+13],&gmmx_data.comp_atoms[i+14],
                  &gmmx_data.comp_atoms[i+15],&gmmx_data.comp_atoms[i+16],&gmmx_data.comp_atoms[i+17],&gmmx_data.comp_atoms[i+18],&gmmx_data.comp_atoms[i+19]);
            }
            goto L_30;
        }  else if (strcmp(pcmfile.token, "COMMENT") == 0)
        { 
                sscanf(pcmfile.string,"Comment %64c",gmmx_comment);
                for(i=0; i < 60; i++)
                {
                    if (gmmx_comment[i] == '\n')
                    {
                        gmmx_comment[i] = '\0';
                        break;
                    }
                }
                goto L_30;
        }  else if (strcmp(pcmfile.token, "INCLUDE") == 0)
        { 
                NL;
                gmmx_data.include_ftype = atoi(pcmfile.token);
                gmmx_data.include_file = TRUE;
                fscanf(wfile,"%s\n",gmmx_data.include_name);
                for(i=0; i < 36; i++)
                {
                    if (gmmx_data.include_name[i] == '\n')
                    {
                        gmmx_data.include_name[i] = '\0';
                        break;
                    }
                }
                goto L_30;
        }  else if (strcmp(pcmfile.token, "SECOND_CYCLE") == 0)
        { 
                sec_cycle = TRUE;
                fscanf(wfile,"%s\n",gmmx_data.foname);
                fscanf(wfile,"%s\n",gmmx_data.foname2);
                goto L_30;
        }  else if (strcmp(pcmfile.token, "Q_DIS") == 0) 
        {
                NL;
                  gmmx_data.ndist = atoi(pcmfile.token);
       if (gmmx_data.ndist  >> 0)
       {
          for(j=0; j < gmmx_data.ndist; j++)
          {
            fscanf(wfile,"%d %d ",&gmmx_data.dist_atoms[j][0], &gmmx_data.dist_atoms[j][1]);
          }
        }
                goto L_40;
        }  else if (strcmp(pcmfile.token, "Q_ANG") == 0) 
        {
                NL;
                  gmmx_data.nang = atoi(pcmfile.token);
       if (gmmx_data.nang  >> 0)
       {
          for(j=0; j < gmmx_data.nang; j++)
          {
            fscanf(wfile,"%d %d %d ",&gmmx_data.ang_atoms[j][0], &gmmx_data.ang_atoms[j][1], &gmmx_data.ang_atoms[j][2]);
          }
        }
                goto L_40;
        }  else if (strcmp(pcmfile.token, "Q_TOR") == 0) 
        {
                NL;
                  gmmx_data.ndihed = atoi(pcmfile.token);
       if (gmmx_data.ndihed  >> 0)
       {
          for(j=0; j < gmmx_data.ndihed; j++)
          {
            fscanf(wfile,"%d %d %d %d ",&gmmx_data.dihed_atoms[j][0], &gmmx_data.dihed_atoms[j][1], &gmmx_data.dihed_atoms[j][2], &gmmx_data.dihed_atoms[j][3]);
          }
        }
                goto L_40;
        }  else if (strcmp(pcmfile.token, "Q_PMR") == 0) 
        {
                NL;
                  gmmx_data.npmr = atoi(pcmfile.token);
       if (gmmx_data.npmr  >> 0)
       {
          for(j=0; j < gmmx_data.npmr; j++)
          {
            fscanf(wfile,"%d %d ",&gmmx_data.pmr_atoms[j][0], &gmmx_data.pmr_atoms[j][1]);
          }
        }
                goto L_40;
        }

L_170:
     fclose(wfile);
     
    strcpy(Openbox.fname,gmmx_data.finame);
    oldfield = field.type;
    if (field.type != MMX)
        field.type = MMX;
// check filetype
    infile = fopen(gmmx_data.finame,"r");
    if (infile == NULL)
    {
        printf("Unable to open - %s - Please check filename\n",gmmx_data.finame);
        exit(0);
    }
    input_type = 0;
    while ( FetchRecord(infile,line) )
    {
        if (strncasecmp(line,"{PCM",4) == 0)
        {
            input_type = FTYPE_PCM;
            break;
        } else
        {
            input_type = FTYPE_SDF;
            break;
        }
    }
    fclose(infile);
    if (input_type == FTYPE_PCM)
       pcmfin(1,0);
       
    find_rings();
    generate_bonds();
    nret = setup_calculation();
    if (nret == TRUE)
    {
        minimize();
    }
    end_calculation();
    
    nheav = 0;
    nhyd = 0;
    for (i=1; i <= natom; i++)
    {
        if (atom[i].atomnum != 1 && atom[i].atomnum != 2)
            nheav++;
        else
            nhyd++;
    }

    gmmx_data.kstop = 5;
    gmmx_data.kmin = 5*nheav;
    if (gmmx_data.kdup == 0)
    {
       gmmx_data.kdup = 5*nheav/2;
       if (gmmx_data.kdup > 50)
          gmmx_data.kdup = 50;
    }
}    
void run_gmmx()
{
    int nret;
    char astring[80];

//    FILE *rmsfile;
    // assume force field is set by pcmodel
    // generate file names
    // setup calculation
    // if include file read and minimize
    // call montc to do search
    // run 2nd cycle
    // process queries and output results
    
    if (restart_job == FALSE)
    {
       strcpy(gmmx_data.jobname,"scratch");
       strcpy(gmmx.tmpname, gmmx_data.jobname);
       strcat(gmmx.tmpname, ".tmp");
    
       strcpy(gmmx.rmsname, gmmx_data.jobname);
       strcat(gmmx.rmsname, ".rms");
    
       strcpy(gmmx.pkmname, gmmx_data.jobname);
       strcat(gmmx.pkmname, ".pkm");
    
       strcpy(gmmx.rstname, gmmx_data.jobname);
       strcat(gmmx.rstname, ".rst");
    }

     
    if(gmmx_data.nopi)
    {
//       pimarkselection();
    } else
    {
       pireset();    
    }
    
    gmmx.run = TRUE;
    gmmx_abort = FALSE;
    nret = setup_calculation();
    if (nret == FALSE)
    {
            printf("Error setting up calculation\n");
            return;
    }

    xold = dvector(0,natom+1);
    yold = dvector(0,natom+1);
    zold = dvector(0,natom+1);

    gmmx.job_done = FALSE;
    montc();
//    printf("done in montc\n");
    
    free_dvector(xold,0,natom+1);
    free_dvector(yold,0,natom+1);
    free_dvector(zold,0,natom+1);
    
    end_calculation();

    // put code here to run 2nd cycle
    if (gmmx.job_done == TRUE)
    {
        if (gmmx_abort == STOP)
        {
           sort_file(gmmx.tmpname,gmmx_data.foname);
           remove(gmmx.rmsname);
           remove(gmmx.tmpname);
           sprintf(astring,"Gmmx job %s is finished. The lowest energy structure will be displayed", gmmx_data.jobname);
           message_alert(astring,"GMMX Run");

           strcpy(Openbox.fname,gmmx_data.foname);
    
           initialize();
           pcmfin(1,0);
     
           if(gmmx_data.nopi)
             pimarkselection();  
           else
             pireset();
           mmxsub(0);
           gmmx.run = FALSE;
        }else if (gmmx_abort == PAUSE)
        {
            write_restart();
            sprintf(astring,"Restart file %s writen",gmmx.rstname);
            message_alert(astring,"GMMX Restart");
        } else
        {
            if (gmmx.nfound <= 1)
            {
                strcpy(Savebox.fname,gmmx_data.foname2);
                pcmfout(1);
                gmmx.run = FALSE;
                return;
            }
           sort_file(gmmx.tmpname,gmmx_data.foname);
           remove(gmmx.rmsname);
           remove(gmmx.tmpname);
           strcpy(Openbox.fname,gmmx_data.foname);
           gmmx2();
        }
    }
    gmmx.run = FALSE;
}
/* ====================================================== */
void montc()
{
    int i, j, jj, nrc, kntcount,k, nsf;
    int ia1, ia2, jj1, jj2, ia0, ia3;
    clock_t start_time;
    double etot;
    double **oldxyz;

    if (restart_job == FALSE)
    {
       gmmx.namvt = 0;
       gmmx.nrbond = 0;
       gmmx.neminim = 0;
       gmmx.nfound = 0;
       gmmx.nsaved = 0;
       gmmx.nsflow = 1;
       gmmx.ndupl = 0;
       gmmx.maxdupl = 0;
       gmmx.nlast = 0;
       gmmx.hybrida = FALSE;
       gmmx.ebendold = 1000.1;
       for (i=0; i < MAXFND; i++)
         gmmx.nsfound[i] = 0;
    }
    kntcount = 0;
    oldxyz = dmatrix(0,natom+1, 0,3);

//    start_time = clock();
     start_time = 12;
     srand(start_time);
//  count number of attachments to each atom
    jji = ivector(1,natom+1);
    for (i=1; i <= natom; i++)
    {
        jj = 0;
        for (j=0; j < MAXIAT; j++)
          if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
            jj++;
        jji[i] = jj;
    }

//  check chirality
    hster(0, jji);
//  check epimers
    if (gmmx_data.chig)
        nrc = lepimer(0);
//  check double bond stereochemistry
    nrc = bstereo(0);
//  check for bad15
    nrc = bad15(0);

    if (restart_job == FALSE)
    {
//  count number of movable atoms
       for (i=1; i <= natom; i++)
          if (jji[i] > 1)
            gmmx.namvt++;
         
//  count number of rotatable bonds and store in rot bond array
       for (i=0; i < gmmx_data.nrings; i++)
       {
          deletebond(gmmx_data.clo_bond[1][i], gmmx_data.clo_bond[2][i]);
          for (j=0; j < gmmx_data.rng_size[i]-1; j++)
          {
            if (gmmx_data.ring_atoms[j][i] != gmmx_data.clo_bond[1][i] &&
             gmmx_data.ring_atoms[j+1][i] != gmmx_data.clo_bond[2][i])
            {
                ia1 = gmmx_data.ring_atoms[j][i];
                ia2 = gmmx_data.ring_atoms[j+1][i];
                ia0 = 0;
                ia3 = 0;
                jj1 = 0;
                jj2 = 0;
                for (k=0; k < MAXIAT; k++)
                {
                    if (atom[ia1].iat[k] != 0 && atom[ia1].bo[k] != 9)
                    {
                         jj1++;
                         if (ia0 == 0 && atom[ia1].iat[k] != ia2)
                           ia0 = atom[ia1].iat[k]; 
                    }
                }
                for (k=0; k < MAXIAT; k++)
                {
                    if (atom[ia2].iat[k] != 0 && atom[ia2].bo[k] != 9)
                    {
                        jj2++;
                        if (ia3 == 0 && atom[ia2].iat[k] != ia1)
                           ia3 = atom[ia2].iat[k];
                    }
                }
                if (jj1 > 1 && jj2 > 1 )
                {  
                    gmmx_data.bond_data[gmmx_data.nbonds][0] = ia1;
                    gmmx_data.bond_data[gmmx_data.nbonds][1] = ia2;
                    gmmx_data.bond_data[gmmx_data.nbonds][3] = ia0;
                    gmmx_data.bond_data[gmmx_data.nbonds][4] = ia3;
                    if (gmmx_data.rng_size[i] <=7 )
                       gmmx_data.bond_data[gmmx_data.nbonds][2] = 30;
                    else
                      gmmx_data.bond_data[gmmx_data.nbonds][2] = 60;
                    gmmx_data.nbonds++;
                }
            }
          }
          if (gmmx_data.ring_atoms[0][i] != gmmx_data.clo_bond[1][i] &&
             gmmx_data.ring_atoms[gmmx_data.rng_size[i]-1][i] != gmmx_data.clo_bond[2][i])
          {
                ia1 = gmmx_data.ring_atoms[0][i];
                ia2 = gmmx_data.ring_atoms[gmmx_data.rng_size[i]-1][i];
                ia0 = 0;
                ia3 = 0;
                jj1 = 0;
                jj2 = 0;
                for (k=0; k < MAXIAT; k++)
                {
                    if (atom[ia1].iat[k] != 0 && atom[ia1].bo[k] != 9)
                    {
                        jj1++;
                        if (ia0 == 0 && atom[ia1].iat[k] != ia2)
                           ia0 = atom[ia1].iat[k];
                    }
                }
                for (k=0; k < MAXIAT; k++)
                {
                    if (atom[ia2].iat[k] != 0 && atom[ia2].bo[k] != 9)
                    {
                         jj2++;
                         if (ia3 == 0 && atom[ia2].iat[k] != ia1)
                            ia3 = atom[ia2].iat[k];
                    }
                }
                if (jj1 > 1 && jj2 > 1 )
                {  
                    gmmx_data.bond_data[gmmx_data.nbonds][0] = ia1;
                    gmmx_data.bond_data[gmmx_data.nbonds][1] = ia2;
                    gmmx_data.bond_data[gmmx_data.nbonds][3] = ia0;
                    gmmx_data.bond_data[gmmx_data.nbonds][4] = ia3;
                    if (gmmx_data.rng_size[i] <=7 )
                       gmmx_data.bond_data[gmmx_data.nbonds][2] = 30;
                    else
                      gmmx_data.bond_data[gmmx_data.nbonds][2] = 60;
                    gmmx_data.nbonds++;
                }
            }
            make_bond(gmmx_data.clo_bond[1][i], gmmx_data.clo_bond[2][i],1);
       }
       gmmx.nrbond = gmmx_data.nbonds;
       if (gmmx_data.include_file == TRUE)
       {
         include_min();
       } else
       {
        if (pistuf.norb > 1) initial_picalc();
         minimize();

// Hay debug
//         eheat();

         optimize_data.final_energy = energies.total;
         optimize_data.final_heat = get_heat_of_formation();
         
         gmmx.eminim = energies.total;
         gmmx.ecurrent = energies.total;
         gmmx.econf[gmmx.nfound] = gmmx.eminim; // first conformation found
         gmmx.hfconf[gmmx.nfound] = get_heat_of_formation();
         gmmx.dpmconf[gmmx.nfound] = dipolemom.total;
         sprintf(Struct_Title,"%i %9.3f %9.3f",gmmx.nfound + 1,gmmx.ecurrent,get_heat_of_formation());

         if (gmmx_data.include_file == FALSE)
         {
           files.append = FALSE;
           strcpy(Savebox.path,Openbox.path);
           strcpy(Savebox.fname,gmmx.tmpname);
           pcmfout(1);
         }
         gmmx_data.include_file = FALSE;
         gmmx.neminim++;
         gmmx.nfound++;
         gmmx.isitdupl = UNIQUE;
          nrc = rmsck(0,&nsf);
       }
    } else
    {
       strcpy(Savebox.path,Openbox.path);
       strcpy(Savebox.fname,gmmx.tmpname);
    }

    const_fail = 0;
    while (TRUE)
    {
        for(i=1; i <=natom; i++)
        {
            oldxyz[i][0] = atom[i].x;
            oldxyz[i][1] = atom[i].y;
            oldxyz[i][2] = atom[i].z;
        }           
        // generate new structure
        if (gmmx_data.method == 1)
        {
            gmmx.smethod = BONDS;
            mvbond();
        }else if ( gmmx_data.method == 2)
        {
            gmmx.smethod = CART;
            mvatom();
        }else
        {
            if (kntcount%2 == 0)
            {
              gmmx.smethod = BONDS;
              mvbond();
            }else
            {
              gmmx.smethod = CART;
              mvatom();
            }
            kntcount++; 
        }        
        // check on bad interactions
        nrc = check_structure();
        if (nrc == FALSE)
        { /* get old coordinates and start again */
            kntcount--;
            const_fail++;
            if (const_fail > 10000)
            {
               kntcount++;
//             printf("> 10000 constraint failures\n");
            }
            for(i=1; i <=natom; i++)
            {
                atom[i].x = oldxyz[i][0];
                atom[i].y = oldxyz[i][1];
                atom[i].z = oldxyz[i][2];
            }
            goto L_20;
        }         
        const_fail = 0;
        gmmx.hybrida = FALSE;
        pcmfout(2);
        initial_picalc();
        minimize();

// Hay debug
//        eheat();

        optimize_data.final_energy = energies.total;
        optimize_data.final_heat = get_heat_of_formation();
        gmmx.ecurrent = energies.total;
        gmmx.neminim++;

        if (gmmx.hybrida == TRUE && gmmx_data.hybrid == TRUE)
        {
            gmmx.isitdupl = BADHYBRID;
            for(i=1; i <=natom; i++)
            {
                atom[i].x = oldxyz[i][0];
                atom[i].y = oldxyz[i][1];
                atom[i].z = oldxyz[i][2];
            }
            goto L_101;
        }

        if (gmmx.ecurrent >= 1000000.0)
        {
            for(i=1; i <=natom; i++)
            {
                atom[i].x = oldxyz[i][0];
                atom[i].y = oldxyz[i][1];
                atom[i].z = oldxyz[i][2];
            }
            gmmx.isitdupl = BADENERGY;
            goto L_101;
        }         
        
        // check for bad energy or epimerization
        nrc = check_id();
        if (nrc == TRUE)
        {
            if (gmmx.ecurrent < (gmmx.eminim+gmmx_data.ewindow) )
            {
               sprintf(Struct_Title,"%i %9.3f %9.3f",gmmx.nfound + 1,gmmx.ecurrent,get_heat_of_formation());
               files.append = TRUE;
               optimize_data.final_energy = energies.total;
               optimize_data.final_heat = get_heat_of_formation();
               pcmfout(1);
               if (gmmx.ecurrent < -600.0)
               {
                   gradient();
                   printf("error - exit in gmmx_run energy < -600\n");
                   fprintf(pcmoutfile,"error - exit in gmmx_run energy < -600\n");
                   fclose(pcmoutfile);
                   exit(0);
               }
               gmmx.nlastk = gmmx.nfound;
               gmmx.nfound++;               
               gmmx.isitdupl = UNIQUE;
               if (gmmx.ecurrent < gmmx.eminim)
               {
                   gmmx.eminim = gmmx.ecurrent;
                   etot = energy();
                   gmmx.nsflow = gmmx.nsaved;
               }
            }else
            {
                gmmx.isitdupl = UNIQUEHIGH;
                   for(i=1; i <=natom; i++)
                   {
                      atom[i].x = oldxyz[i][0];
                      atom[i].y = oldxyz[i][1];
                      atom[i].z = oldxyz[i][2];
                   }         
            }
        } else
        {
                   for(i=1; i <=natom; i++)
                   {
                      atom[i].x = oldxyz[i][0];
                      atom[i].y = oldxyz[i][1];
                      atom[i].z = oldxyz[i][2];
                   }         
        }
L_101:
        write_gmmx();
L_20:
        nrc = check_stop();
        if (nrc == FALSE)
            break;
    }            
    // release all memory
    nrc = bad15(-1);
    if (gmmx_data.chig)
        nrc = lepimer(-1);
    nrc = bstereo(-1);
    hster(-1, jji);
    free_ivector(jji, 1, natom+1);
    free_dmatrix(oldxyz, 0,natom+1,0,3);   
}
/* ========================================  */
int bstereo(int mode)
{
    int i,j, jjl, jjk;
    int ia, ib,ic,id;
    double rdihed;
    
    if (mode == 0)
    {
        generate_bonds();
        npibond = 0;
        for (i=0; i < bonds.numbonds; i++)
        {
            if (bonds.bondorder[i] == 2)
            {
                jjl = 0;
                jjk = 0;
                if (atom[bonds.ia1[i]].atomnum == 6 && atom[bonds.ia2[i]].atomnum == 6)
                {
                    for (j=0; j < MAXIAT; j++)
                    {
                        if ( atom[bonds.ia1[i]].iat[j] != 0 && atom[atom[bonds.ia1[i]].iat[j]].atomnum != 1)
                           jjl++;
                        if (atom[bonds.ia2[i]].iat[j] != 0 && atom[atom[bonds.ia2[i]].iat[j]].atomnum != 1)
                           jjk++;
                    }
                    if (jjl >= 2 && jjk >= 2)
                         npibond++;
                }
            }
        }
        if (npibond > 0)
        {
           ibstereo = imatrix(0,npibond, 0,5);
           nbond = 0;
           for (i=0; i < bonds.numbonds; i++)
           {
               ia = 0;
               id = 0;
               if (bonds.bondorder[i] == 2)
               {
                   if (atom[bonds.ia1[i]].atomnum == 6 && atom[bonds.ia2[i]].atomnum == 6)
                   {
                       ib = bonds.ia1[i];
                       ic = bonds.ia2[i];
                       for (j=0; j < MAXIAT; j++)
                       {
                           if (atom[ib].iat[j] != 0 && atom[ib].bo[j] != 9)
                           {
                               if (atom[ib].iat[j] != ic && atom[atom[ib].iat[j]].atomnum != 1)
                               {
                                   ia = atom[ib].iat[j];
                                   break;
                               }
                           }
                       }
                       for (j=0; j < MAXIAT; j++)
                       {
                           if (atom[ic].iat[j] != 0 && atom[ic].bo[j] != 9)
                           {
                               if (atom[ic].iat[j] != ib && atom[atom[ic].iat[j]].atomnum != 1)
                               {
                                   id = atom[ic].iat[j];
                                   break;
                               }
                           }
                       }
                       if (ia != 0 && id != 0)
                       {
                           ibstereo[nbond][0] = ia;
                           ibstereo[nbond][1] = ib;
                           ibstereo[nbond][2] = ic;
                           ibstereo[nbond][3] = id;
                           rdihed = dihdrl(ia,ib,ic,id);
                           ibstereo[nbond][4] = 1;
                           if (rdihed > 90.0 || rdihed < -90.0)
                              ibstereo[nbond][4] = 2;
                           nbond++;
                       }
                   }
               }    
           }
        }      
    } else if (mode == 1)
    {
        // check bond stereochemistry
        for (i=0; i < nbond; i++)
        {
            ia = ibstereo[i][0];
            ib = ibstereo[i][1];
            ic = ibstereo[i][2];
            id = ibstereo[i][3];
            rdihed = dihdrl(ia,ib,ic,id);
            j = 1;
            if (rdihed > 90.0 || rdihed < -90.0)
               j = 2;
            if (j != ibstereo[i][4])
               return FALSE;
        }            
    } else if (mode == -1)
    {
        if (npibond > 0)
           free_imatrix(ibstereo ,0,npibond, 0,5);
    }
    
    return TRUE;
}
/* ========================================  */
void hster(int mode, int *jji)
{
    // mode = 0 allocate memory and setup irest array
    // mode = -1 free memory
    // mode = 1 check stereochemistry of hydrogens
    int i, irsa1, irsa2, irsa3, irsa4, irsd;
    double zsrd, xtmp, ytmp, ztmp;

    if (mode == -1)
    {
        free_ivector( irest, 1, natom+1);
        return;
    }
    if (mode == 0)
    {
        irest = ivector(1, natom+1);

        for (i=1; i <= natom; i++)
        {
            if ( atom[i].atomnum == 6 && jji[i] == 4)
            {
                irsa1 = atom[i].iat[0];
                irsa2 = atom[i].iat[1];
                irsa3 = atom[i].iat[2];
                irsa4 = atom[i].iat[3];
                if (irsa1 == 0 || irsa2 == 0 || irsa3 == 0 || irsa4 == 0 )
                  break;
                zsrd = dihdrl(irsa1, irsa2, irsa3, irsa4);
                irsd = (int) zsrd;
                irest[i] = irsd;
            } else
                 irest[i] = 0;
        }
    } else if (mode == 1)
    {
        for (i=1; i <= natom; i++)
        {
            if (irest[i] > 0 )
            {
                irsa1 = atom[i].iat[0];
                irsa2 = atom[i].iat[1];
                irsa3 = atom[i].iat[2];
                irsa4 = atom[i].iat[3];
                if (irsa1 == 0 || irsa2 == 0 || irsa3 == 0 || irsa4 == 0 )
                  break;
                zsrd = dihdrl(irsa1, irsa2, irsa3, irsa4);
                irsd = (int) zsrd;
                if (irsd != irest[i])
                {
                    // swap 3 and 4
                    xtmp = atom[irsa3].x;
                    ytmp = atom[irsa3].y;
                    ztmp = atom[irsa3].z;
                    atom[irsa3].x = atom[irsa4].x;
                    atom[irsa3].x = atom[irsa4].y;
                    atom[irsa3].x = atom[irsa4].z;
                    atom[irsa4].x = xtmp;
                    atom[irsa4].x = ytmp;
                    atom[irsa4].x = ztmp;
                }
            }
        }
    }
}  
/* ========================================  */
int lepimer(int mode)
{
    // mode = 0 allocate memory and check epimer
    // mode = -1 free memory
    // mode = 1 check stereochemistry of C, N+ and S+
    // ADDED P to this list on 7-14-2015

    int i, j, type, typemm3, irsa1, irsa2, irsa3, irsa4, isrd, temp, icount;
    double zsrd;

    if (mode == -1)
    {
        free_ivector(iest, 1, natom+1);
        return FALSE;
    } else if (mode == 0)
    {
        iest = ivector(1, natom+1);

        for (i=1; i <= natom; i++)
        {
            iest[i] = 0;
            type = atom[i].mmx_type;
            typemm3 = atom[i].mm3_type;
            if (type == 1 || type == 16 || type == 17 || type == 41 || type == 22 || typemm3 == 153)
            {
                icount = 0;
                for (j=0; j < 4; j++)  // one of first 3 is H or LP, must have 2 H's
                {
                    if (atom[atom[i].iat[j]].atomnum == 1 || atom[atom[i].iat[j]].mmx_type == 20 || atom[atom[i].iat[j]].mmx_type == 0)
                    {
                        icount++;
                        if (icount > 1)
                        {
                           iest[i] = 0;
                           goto L_10;
                        }
                    }
                }
                irsa1 = atom[i].iat[0];
                irsa2 = atom[i].iat[1];
                irsa3 = atom[i].iat[2];
                irsa4 = atom[i].iat[3];
                if (irsa1 == 0 || irsa2 == 0 || irsa3 == 0 || irsa4 == 0)
                {
                    iest[i] = 0;
                    goto L_10;
                }
                zsrd = dihdrl(irsa1, irsa2, irsa3, irsa4);
                isrd = (int) zsrd;
                if (isrd > 5)
                {
                   iest[i] = 1;
                }else if (isrd < -5)
                {
                   iest[i] = -1;
                }
            }
L_10:
            continue;
        }
        return FALSE;              
    } else if (mode == 1)
    {
        for (i=1; i <= natom; i++)
        {
            if (iest[i] != 0 )
            {
                irsa1 = atom[i].iat[0];
                irsa2 = atom[i].iat[1];
                irsa3 = atom[i].iat[2];
                irsa4 = atom[i].iat[3];
                zsrd = dihdrl(irsa1, irsa2, irsa3, irsa4);
                isrd = (int) zsrd;
                if (isrd > 0 )
                  temp = 1;
                else
                  temp = -1;
                if (iest[i] != temp)
                {
                   return FALSE;   
                }             
            }
        }
        return TRUE;
    }
    return FALSE;
}
/*  ============================================== */
void mvbond()
{
    int i, irand, ia1, ia2, ires, ia0, ia3, icount,j, irot;
    float rnb, rotdeg, randtmp,rate;
    double  drr, rdist;

    icount = 0;
 
    if (gmmx_data.nrings > 0)
    {
        for (i=0; i < gmmx_data.nrings; i++)
        {
            // delete closure bond
            deletebond(gmmx_data.clo_bond[1][i], gmmx_data.clo_bond[2][i]);
        }
    }
    icount = 0;
    // now rotate all bonds
    rnb = 0.80;
    if (gmmx_data.nbonds <= 4)
      rnb = 1.0;
    else if (gmmx_data.nbonds > 4 && gmmx_data.nbonds < 10)
      rnb = 0.60;

     if (gmmx.neminim > 0)
     {
      rate = (float)gmmx.nfound/(float)gmmx.neminim;
      if (rate < 0.05)
        rnb = .1;
      else if (rate < 0.07)
        rnb = .2;
      else if (rate < 0.10)
        rnb = .4;
     }
    if (const_fail > 25000)
      rnb = 0.2;

L_1:
// save coordinates
    for (i=1; i <= natom; i++)
    {
        xold[i] = atom[i].x;
        yold[i] = atom[i].y;
        zold[i] = atom[i].z;
    }

    for (i=0; i < gmmx_data.nbonds; i++)
    {
        irand = rand();
        randtmp = (float)irand/(float)RAND_MAX;
        if (randtmp < rnb)
        {
           ia1 = gmmx_data.bond_data[i][0];
           ia2 = gmmx_data.bond_data[i][1];
           ia0 = gmmx_data.bond_data[i][3];
           ia3 = gmmx_data.bond_data[i][4];
           ires = gmmx_data.bond_data[i][2];

//Hay debug
//           printf("rotate on bond %3d %3d\n",ia1,ia2);


           if (ires == 0)
             irot = 0;
           else
             irot = 360/ires;

              irand = rand();
              randtmp = (float)irand/(float)RAND_MAX;
           // irand = (int)(randtmp*(float)irot)+1.0;
           rotdeg = randtmp*ires*irot;
           if (rotdeg > 360)
              rotdeg -= 360;
           //
           irand = rand();
           randtmp = (float)irand/(float)RAND_MAX;
           if (randtmp > 0.5)
             rotdeg = -rotdeg;
           if (rotdeg != 0)
           {
              rotdeg += dihdrl(ia0,ia1,ia2,ia3);
              rotb(ia0, ia1, ia2, ia3, rotdeg);
           }
        }
    }
    if (gmmx_data.nrings > 0)
    {
        // check closure distance
        for (i=0; i < gmmx_data.nrings; i++)
        {
            ia1 = gmmx_data.clo_bond[1][i];
            ia2 = gmmx_data.clo_bond[2][i];
            rdist = distance(ia1, ia2);
            drr = fabs(rdist-gmmx_data.clo_distance[i]);
            if (drr > 0.90)
            {
                    for (j=1; j <= natom; j++)
                    {
                        atom[j].x = xold[j];
                        atom[j].y = yold[j];
                        atom[j].z = zold[j];
                    }
                    icount++;
                    goto L_1;
            }
        }
        // check closure angle
     }
    
    for(i=0; i < gmmx_data.nrings; i++)
      make_bond(gmmx_data.clo_bond[1][i],gmmx_data.clo_bond[2][i],1);
}
/*  ============================================== */
void mvatom()
{
    int i, j, nheav, itest, itest1, irand;
    int moved;
    double rnb, pn, amov[3], temp,rate;    

    moved = FALSE;
    nheav = 0;
    for (i= 1; i <= natom; i++)
    {
        if (atom[i].type != 0 && atom[i].atomnum > 1)
         nheav++;
    }

    rnb = 0.50;
    if (nheav <= 6)
      rnb = 1.0;
    if (const_fail > 1000)
      rnb = 0.90;
    rate = 0.9;
    pn = 0.7;
    if (gmmx.neminim > 0)
    {
      rate = (double)gmmx.nfound/(double)gmmx.neminim;
      if (rate < 0.05)
        pn = 3.0;
      else if (rate < 0.07)
        pn = 2.0;
      else if (rate < 0.10)
        pn = 1.5;
     }
 L_1:
    for (i=1; i <= natom; i++)
    {
        if (atom[i].use && jji[i] > 1)
        {
            itest = findcd(i); // check constrained distances
            if (itest)
            {
                itest1 = findtor(i); // check constrained torsions
                if (itest1)
                {
                    irand = rand();
                    temp = (double)irand/(double)RAND_MAX;
                    if (temp <= rnb)
                    {
                        for (j=0; j < 3; j++)
                        {
                            irand = rand();
                            temp = (double)irand/(double)RAND_MAX;
                            if (temp > 0.50)
                               pn = -pn;
                            irand = rand();
                            temp = (double)irand/(double)RAND_MAX;
                            if (temp < 0.5)
                               temp = 0.5;
                            amov[j] = pn*temp;
                        }
                        atom[i].x += amov[0];
                        atom[i].y += amov[1];
                        atom[i].z += amov[2];
                        moved = TRUE;
                        for (j=0; j < MAXIAT; j++)
                        {
                          if (atom[i].iat[j] != 0 && jji[atom[i].iat[j]] == 1)
                          {
                              atom[atom[i].iat[j]].x += amov[0];
                              atom[atom[i].iat[j]].y += amov[1];
                              atom[atom[i].iat[j]].z += amov[2];
                          }
                        }
                    }
                }
            }
        }
    }
    if (moved == FALSE) 
      goto L_1;
}
/* =========================================  */
int check_id()
{
    int nrc, nsf, ndupl_old;

    nsf = 0;
    ndupl_old = gmmx.ndupl;
    nrc = rmsck(0, &nsf);
    if (nrc == FALSE)
    {
        // structure rejected as identical
        gmmx.nsfound[nsf]++;
        if (nsf != gmmx.nlast)
        {
            if ( (gmmx.ecurrent < (gmmx.eminim+gmmx_data.ewindow)) &&
                 (gmmx.nlast != gmmx.nlastk) )
                      gmmx.ndupl++;
            if (gmmx.ndupl > gmmx.maxdupl) 
                      gmmx.maxdupl = gmmx.ndupl;
            gmmx.nlast = nsf;
        } 
        if (nsf == gmmx.nsflow && nsf != gmmx.nlast)
           gmmx.neminim++;       
        return FALSE;
    }
    nrc = bstereo(1);
    if (nrc == FALSE)
    {
        gmmx.nsflow = gmmx.nsfmin;
        gmmx.ndupl = ndupl_old;
        gmmx.isitdupl = BADDBOND;
        return FALSE;
    }
    if (gmmx_data.chig)  // check for epimerization
    {
       nrc = lepimer(1);
       if (nrc == FALSE)
       {
         gmmx.nsflow = gmmx.nsfmin;
         gmmx.ndupl = ndupl_old;
         gmmx.isitdupl = BADEPIMER;
         return FALSE;
       }
    }
    nrc = findcd(0);  // check for constrained distances
    if (nrc == FALSE)
    {
        gmmx.nsflow = gmmx.nsfmin;
        gmmx.ndupl = ndupl_old;        
        gmmx.isitdupl = BADCONDIST;        
        return FALSE;
    }
    
    nrc = findtor(0); // check for constrained torsions
    if (nrc == FALSE)
    {
        gmmx.nsflow = gmmx.nsfmin;
        gmmx.ndupl = ndupl_old;
        gmmx.isitdupl = BADCONTOR;
        return FALSE;
    }

    return TRUE;
}
   
/* =========================================  */
int check_structure()
{
    int nrc;

    if (gmmx.smethod == CART)
    {
        nrc = bad15(1); 
        if (nrc == FALSE)
        { 
//        printf("bad 15\n");
          return FALSE;
        }
        else
        {
           return TRUE;
        }
    }

    if (gmmx_data.chig)  // check for epimerization
    {
       nrc = lepimer(1);
       if (nrc == FALSE)
       {
//       printf("bad epimer\n");
         return FALSE;
       }
    }

    nrc = bstereo(1);  // check for double bond isomerization
    if (nrc == FALSE)
    {
//     printf("bad stereo\n");
       return FALSE;
    }

    nrc = findcd(0);  // check for constrained distances
    if (nrc == FALSE)
    {
//     printf("bad cd\n");
       return FALSE;
    }
       
    nrc = findtor(0); // check for constrained torsions
    if (nrc == FALSE)
    {
//     printf("bad ctor\n");
       return FALSE;
    }

     nrc = bad15(1);   // bad 1,5 CC interactions
     if (nrc == FALSE)
    {
//     printf("bad 15\n");
       return FALSE;
    }

    nrc = close_contact(1);  // bad vdw interactions
    if (nrc == FALSE)
    {
//     printf("bad 14\n");
       return FALSE;
    }
    
    return TRUE;  // passed all tests
}
/* =========================================  */

int findcd(int atom)
{
    // atom = 0 check all distances
    // atom > 0 look for specific distance containing this atom
    return TRUE;
}

int findtor(int atom)
{
    return TRUE;
}

//PCModel 9 version of function
/* =========================================  */
int close_contact(int mode)
{
  int i, j;
  double dx,dy,dz;

  for (i=1; i < natom; i++)
  {
    for (j=i+1; j <= natom; j++)
      {
        dx = fabs(atom[i].x - atom[j].x);
        dy = fabs(atom[i].y - atom[j].y);
        dz = fabs(atom[i].z - atom[j].z);
        if ( (dx + dy + dz) < 0.6)
        {
            return FALSE;
        }
      }
  }
    return TRUE;
}

/* =========================================  */
int bad15(int mode)
{
    // mode = 0 setup and allocate memory
    // mode = -1 release memory
    // mode = 1 check for bad 15

    int i, j, ia1, ia2, ia3, ia4, it1, it2, k, found, count;
    double rdist;

    if (mode == -1)
    {
        n15 = 0;
        free_imatrix(k15, 0,4*natom, 0,2);
        return TRUE;
    }else if (mode == 0)
    {
        n15 = 0;
        nbad15 = 0;
        k15 = imatrix(0, 4*natom, 0, 2);

        for (i = 0; i < torsions.ntor; i++)
        {
            ia1 = torsions.i14[i][0];
            ia2 = torsions.i14[i][1];
            ia3 = torsions.i14[i][2];
            ia4 = torsions.i14[i][3];
            
            if (atom[ia1].atomnum == 6)
            {
                for (j = 0; j < MAXIAT; j++)
                {
                    if (atom[atom[ia4].iat[j]].atomnum == 6 && atom[ia4].iat[j] != ia3)
                    {
                        if ( torsions.i14[i][0] < atom[ia4].iat[j])
                        {
                            it1 = torsions.i14[i][0];
                            it2 = atom[ia4].iat[j];
                        }else
                        {
                            it2 = torsions.i14[i][0];
                            it1 = atom[ia4].iat[j];
                        }
                        // check to see if already in list
                        found = FALSE;
                        for (k=0; k < n15; k++)
                        {
                            if (  it1 == k15[k][0] && it2 == k15[k][1])
                              found = TRUE;
                        }
                        if (found == FALSE)
                        {  
                           k15[n15][0] = it1;
                           k15[n15][1] = it2;
                           rdist = distance(it1, it2);
                           if (rdist < gmmx_data.bad15_cutoff)
                              nbad15++;
                           n15++;
                        }
                    }
                }
            }
            if (atom[ia4].atomnum == 6)
            {
                for (j = 0; j < MAXIAT; j++)
                {
                    if (atom[atom[ia1].iat[j]].atomnum == 6 && atom[ia1].iat[j] != ia2)
                    {
                        if ( torsions.i14[i][3] < atom[ia1].iat[j])
                        {
                            it1 = torsions.i14[i][3];
                            it2 = atom[ia1].iat[j];
                        }else
                        {
                            it2 = torsions.i14[i][3];
                            it1 = atom[ia1].iat[j];
                        }
                        // check to see if already in list
                        found = FALSE;
                        for (k=0; k < n15; k++)
                        {
                            if (  it1 == k15[k][0] && it2 == k15[k][1])
                              found = TRUE;
                        }
                        if (found == FALSE)
                        {  
                           k15[n15][0] = it1;
                           k15[n15][1] = it2;
                           rdist = distance(it1, it2);
                           if (rdist < gmmx_data.bad15_cutoff)
                              nbad15++;
                           n15++;
                        }
                    }
                }
            }
            if (nbad15 < 5) nbad15 = 5;
        }
    } else if (mode == 1)
    {
        count = 0;
        for (i=0; i < n15; i++)
        {
            ia1 = k15[i][0];
            ia2 = k15[i][1];
            rdist = distance(ia1, ia2);
            if (rdist < gmmx_data.bad15_cutoff)
            {
                count++;
                if (count > nbad15 + 1)
                {
                   return FALSE;
                }
            }
        }
    }
    return TRUE;
}
/* =========================================  */
int check_stop()
{
    if (gmmx_abort == TRUE)
       return FALSE;
    else if (gmmx.neminim >= gmmx_data.max_search)
    {
       gmmx.job_done = TRUE;
       return FALSE;
    }else if (gmmx.neminim > gmmx_data.kmin && gmmx.ndupl > gmmx_data.kdup)
    {
       gmmx.job_done = TRUE;
       return FALSE;
    }else if (gmmx.neminim > gmmx_data.kmin && gmmx.nsfound[gmmx.nsflow] > gmmx_data.kmin)
    {
       gmmx.job_done = TRUE;
       return FALSE;
    }else if (gmmx_abort == STOP)
    {
       gmmx.job_done = TRUE;
       return FALSE;
    }else if (gmmx_abort == PAUSE)
    {
       gmmx.job_done = TRUE;
       return FALSE;
    } else
       return TRUE;
}
/* =========================================  */
int rmsck(int mode, int *nsfound)
{
    // mode 0 compare as is - if new write to disk

    int i, k, nrc, current_pos, icount;
    long int nread;
    float denergy;
    FILE *rmsfile;
    double ecurrent, xtmp, ytmp, ztmp;
    double **xyz1, **xyz2;
//    double xyz1[150][3], xyz2[150][3];
    char errstring[50];
    
    if (gmmx.nsaved > 0)
    {
       *nsfound = -1;
       current_pos = 0;
       rmsfile = NULL;
       rmsfile = fopen_path(Openbox.path,gmmx.rmsname,"rb");
       rewind(rmsfile);
       icount= 0;
       if (rmsfile == NULL)
       {
           message_alert("Error opening rmsfile","GMMX Error");
           return FALSE;
       }

       xyz1 = dmatrix(1,natom+1,0,3);
       xyz2 = dmatrix(1,natom+1,0,3);

       // copy current coordinates into array1
       for (i=1; i <= natom; i++)
       {
           xyz1[i][0] = atom[i].x;
           xyz1[i][1] = atom[i].y;
           xyz1[i][2] = atom[i].z;
       }
       
       for (i=0; i < gmmx.nsaved; i++)
       {
           denergy = fabs(gmmx.ecurrent-gmmx.econf[i]);
           if( denergy < 0.1)  // ( denergy < gmmx_data.ecutoff)
           {
               icount = i;
               if (icount > 0)
               {
                  nread = ((3*natom + 1)*(icount-current_pos)*sizeof(xtmp));
                  fseek(rmsfile,nread,SEEK_CUR);
               } 
               
               nrc = fread( &ecurrent, sizeof(ecurrent),1,rmsfile);
               for(k=1; k <= natom; k++)
               {
                       nrc = fread( &xtmp, sizeof(xtmp), 1, rmsfile);
                       nrc = fread( &ytmp, sizeof(xtmp), 1, rmsfile);
                       nrc = fread( &ztmp, sizeof(xtmp), 1, rmsfile);
                       xyz2[k][0] = xtmp;
                       xyz2[k][1] = ytmp;
                       xyz2[k][2] = ztmp;
               }
               icount++;
               
               current_pos=icount;
               nrc = fast_compare(xyz1, xyz2);
               if (nrc == TRUE)
               {
                   *nsfound = i;
                   break;
               }
           }
       }

       free_dmatrix(xyz1, 1,natom+1,0,3);
       free_dmatrix(xyz2, 1,natom+1,0,3);
       
       fclose(rmsfile);
       if (*nsfound > -1)
          return FALSE;
    }
    // get here with new structure
    rmsfile = NULL;
    if (gmmx.nsaved == 0)
       rmsfile = fopen_path(Openbox.path,gmmx.rmsname,"wb");
    else
       rmsfile = fopen_path(Openbox.path,gmmx.rmsname,"ab+");

    if (rmsfile == NULL)
    {
        sprintf(errstring,"Error %s\n",strerror(errno) );
        return FALSE;
    }
    nrc = fwrite(&gmmx.ecurrent, sizeof(gmmx.ecurrent), 1, rmsfile);
    for (i=1; i <= natom; i++)
    {
        nrc = fwrite(&atom[i].x, sizeof(atom[i].x), 1, rmsfile);
        nrc = fwrite(&atom[i].y, sizeof(atom[i].y), 1, rmsfile);
        nrc = fwrite(&atom[i].z, sizeof(atom[i].z), 1, rmsfile);
    }
    fclose(rmsfile);
    gmmx.econf[gmmx.nsaved] = gmmx.ecurrent;
    gmmx.hfconf[gmmx.nsaved] = get_heat_of_formation();
    gmmx.dpmconf[gmmx.nsaved] = dipolemom.total;
    gmmx.nsfound[gmmx.nsaved] = 1;
    if (gmmx.ecurrent < gmmx.eminim+gmmx_data.ewindow)
       gmmx.ndupl = 0;
    gmmx.nsaved++;
    *nsfound = gmmx.nsaved;
    //    gmmx.nlast = gmmx.nsaved;  
    return TRUE;
}        
/* ====================================================== */
int fast_compare(double **xyz1, double **xyz2)
//int fast_compare(double xyz1[150][3], double xyz2[150][3])
{
        int i,  maxca, nheavy, j, iatom, k;
        double a, crms,  xav, yav, zav, drms;
        double **xb;
        int **katoms;

        /*     *** COMP COMMAND *** */
        /*     XB(250),YB(250),ZB(250) COMPARISON STRUCTURE
         *     X(250),Y(250),Z(250) ORIGINAL STRUCTURE
         *     N=#OF NON-HYDROGEN ATOMS ^GN IN RMSCK
         *     DRMS=RMS COMP OF STRUCTURES
         *        MSR = 0 NORMAL COMPARE
         *        MSR = 1 NUMERICAL ISOMERS COMPARE
         *        MSR = 2 REFLECTION ISOMER COMPARE
         *        MSR = 3 EQUIVALENT ISOMER COMPARE
         *        MSR = 4 SEQUENTIAL ISOMER COMPARE
         *        NSRMD = 5 COMPARE SPECIFIC ATOMS ONLY
         * */

        katoms = imatrix(1,natom+1,0,2);
        crms = gmmx_data.crmsa;
        for( i = 1; i <= natom; i++ )
        {
            katoms[i][0] = 0;
            katoms[i][1] = 0;
        }
        
        // atom 1
/* center  structures  */
        xav = 0.0;
        yav = 0.0;
        zav = 0.0;
        a = 0.0;
        for( i = 1; i <= natom; i++ )
        {
            if (atom[i].atomnum != 1)
            {
                xav += xyz2[i][0];
                yav += xyz2[i][1];
                zav += xyz2[i][2];
                a += 1.0;
            }
        }

        /*     *** NOW MOVE THE COMPARISON STRUCTURE CENTER OF MASS
         *           TO THE COORDINATE SYSTEM ORIGIN *** */
        xav /= a;
        yav /= a;
        zav /= a;
        for( i = 1; i <= natom; i++ )
        {
                xyz2[i][0] -= xav;
                xyz2[i][1] -= yav;
                xyz2[i][2] -= zav;
        }
/* ============================   */     
        /*     *** AVERAGE MAIN STRUCTURE COORDINATES OF COMMON ATOMS *** */
        xav = 0.0;
        yav = 0.0;
        zav = 0.0;
        a = 0.0;
        for(i = 1; i <= natom; i++ )
        {
            if (atom[i].atomnum != 1)
            {
                xav += xyz1[i][0];
                yav += xyz1[i][1];
                zav += xyz1[i][2];
                a += 1.0;
            }
        }
        /*     *** NOW MOVE THE MAIN STRUCTURE CENTER OF MASS
         *           TO THE COORDINATE SYSTEM ORIGIN *** */
        xav /= a;
        yav /= a;
        zav /= a;
        for(i = 1; i <= natom; i++ )
        {
                xyz1[i][0] -= xav;
                xyz1[i][1] -= yav;
                xyz1[i][2] -= zav;
        }
/* ================================================= */        

        if( gmmx_data.comp == 5 )   /* COMPARE SPECIFIC ATOMS ONLY */
        {
            // data for number of compare atoms needs to be done
            for( i = 1; i <= natom; i++ )
            {
                katoms[i][0] = 0;
                katoms[i][1] = 0;
            }
            maxca = 1;
            for(i=0; i < gmmx_data.ncompatoms; i++)
            {
                katoms[maxca][0] = gmmx_data.comp_atoms[i];
                katoms[maxca][1] = gmmx_data.comp_atoms[i];
                maxca++;
            }
            maxca--;
            if (gmmx_data.comp_method == 0)
               fndumt(maxca,katoms, xyz1, xyz2);
            else
               slow_compare(maxca,katoms,xyz1,xyz2);
            if (nerr == -1)
            {
                drms = 10000.0;
                free_imatrix(katoms, 1,natom+1,0,2);
               return FALSE;
            }
            drms = compute_rms(maxca, katoms, xyz1, xyz2);
            if (drms > gmmx_data.crmsa)
            {
                free_imatrix(katoms, 1,natom+1,0,2);
                return FALSE;
            }else
            {
                gmmx.isitdupl = IDENT;
                free_imatrix(katoms, 1,natom+1,0,2);
                return TRUE;
            }
            
        } else if( gmmx_data.comp == 3 ) // numerical isomers
        {
            nheavy = 0;
            xb = dmatrix(0,natom+1,0,3);
            for(i=1; i <= natom; i++)    // save xyz2
            {
                xb[i][0] = xyz2[i][0];
                xb[i][1] = xyz2[i][1];
                xb[i][2] = xyz2[i][2];
            }
            maxca = 1;
            for( i = 1; i <= natom; i++ )
            {
                if (atom[i].mmx_type != 5 && atom[i].mmx_type != 20)
                {
                
                    katoms[maxca][0] = i;
                    maxca++;
                    nheavy++;
                }
                
            }
            maxca -= 1;
            
            for (i=0; i < nheavy; i++)   // 1 to N
            {
                for (j=1; j <= nheavy; j++)
                {
                    iatom = i+j;
                    if (iatom > nheavy)
                       iatom -= nheavy;
                    katoms[j][1] = katoms[iatom][0];
                }
                for(k=1; k <= nheavy; k++)  // copy xb into xyz2
                {
                    xyz2[k][0] = xb[k][0];
                    xyz2[k][1] = xb[k][1];
                    xyz2[k][2] = xb[k][2];
                }
                if (gmmx_data.comp_method == 0)
                {
                   fndumt(maxca,katoms, xyz1, xyz2);
                } else
                {
                   slow_compare(maxca,katoms,xyz1,xyz2);
                }
                if (nerr == -1)
                {
                   printf("nerr == -1\n");
                   drms = 10000.0;
                   free_imatrix(katoms, 1,natom+1,0,2);
                   free_dmatrix(xb,0,natom+1,0,3);
                   return FALSE;
                }
                drms = compute_rms(maxca, katoms, xyz1, xyz2);
                if (drms < gmmx_data.crmsa)  // structure is identical
                {
                    if (i== 0)
                       gmmx.isitdupl = IDENT;
                    else
                       gmmx.isitdupl = NUMISOMER;
                    free_imatrix(katoms,1, natom+1,0,2);
                    free_dmatrix(xb,0,natom+1,0,3);
                    return TRUE;
                }
                // check enantiomer 
                   for(k=1; k <= nheavy; k++)  // copy xb into xyz2
                   {
                       xyz2[k][0] = xb[k][0];
                       xyz2[k][1] = xb[k][1];
                       xyz2[k][2] = -xb[k][2];
                   }
                   if (gmmx_data.comp_method == 0)
                   {
                       fndumt(maxca,katoms, xyz1, xyz2);
                   } else
                   {
                       slow_compare(maxca,katoms,xyz1,xyz2);
                   }
                   if (nerr == -1)
                   {
                      printf("nerr == -1\n");
                      drms = 10000.0;
                      free_imatrix(katoms, 1,natom+1,0,2);
                      free_dmatrix(xb,0,natom+1,0,3);
                      return FALSE;
                   }
                   drms = compute_rms(maxca, katoms, xyz1, xyz2);
                   if (drms < gmmx_data.crmsa)  // structure is identical
                   {
                       if (i== 0)
                         gmmx.isitdupl = IDENT;
                       else
                         gmmx.isitdupl = NUMISOMER;
                       free_imatrix(katoms,1, natom+1,0,2);
                       free_dmatrix(xb,0,natom+1,0,3);
                       return TRUE;
                   }
            }
            
            for (i = nheavy-1; i >= 0; i--)   // n to 1
            {
                for (j= 1; j <= nheavy; j++)
                {
                    iatom = nheavy+i-j;
                    if (iatom > nheavy)
                       iatom -= nheavy;
                    if (iatom < 1)
                        iatom += nheavy;
                    katoms[j][1] = katoms[iatom][0];
                }
                for(k=1; k <= nheavy; k++)  // copy xb into xyz2
                {
                    xyz2[k][0] = xb[k][0];
                    xyz2[k][1] = xb[k][1];
                    xyz2[k][2] = xb[k][2];
                }
                if (gmmx_data.comp_method == 0)
                {
                   fndumt(maxca,katoms, xyz1, xyz2);
                } else
                {
                   slow_compare(maxca,katoms,xyz1,xyz2);
                }
                if (nerr == -1)
                {
                   printf("nerr == -1\n");
                   drms = 10000.0;
                   free_imatrix(katoms, 1,natom+1,0,2);
                   free_dmatrix(xb,0,natom+1,0,3);
                   return FALSE;
                }
                drms = compute_rms(maxca, katoms, xyz1, xyz2);
                if (drms < gmmx_data.crmsa)  // structure is identical
                {
                    if (i== 0)
                       gmmx.isitdupl = IDENT;
                    else
                       gmmx.isitdupl = NUMISOMER;
                    free_imatrix(katoms,1, natom+1,0,2);
                    free_dmatrix(xb,0,natom+1,0,3);
                    return TRUE;
                }
                // check enantiomer                   
                for(k=1; k <= nheavy; k++)  // copy xb into xyz2
                {
                    xyz2[k][0] = xb[k][0];
                    xyz2[k][1] = xb[k][1];
                    xyz2[k][2] = -xb[k][2];
                }
                if (gmmx_data.comp_method == 0)
                {
                   fndumt(maxca,katoms, xyz1, xyz2);
                } else
                {
                   slow_compare(maxca,katoms,xyz1,xyz2);
                }
                if (nerr == -1)
                {
                   printf("nerr == -1\n");
                   drms = 10000.0;
                   free_imatrix(katoms, 1,natom+1,0,2);
                   free_dmatrix(xb,0,natom+1,0,3);
                   return FALSE;
                }
                drms = compute_rms(maxca, katoms, xyz1, xyz2);
                if (drms < gmmx_data.crmsa)  // structure is identical
                {
                    if (i== 0)
                       gmmx.isitdupl = IDENT;
                    else
                       gmmx.isitdupl = NUMISOMER;
                    free_imatrix(katoms,1, natom+1,0,2);
                    free_dmatrix(xb,0,natom+1,0,3);
                    return TRUE;
                }                   
            }
            free_imatrix(katoms,1, natom+1,0,2);
            free_dmatrix(xb,0,natom+1,0,3);
            return FALSE;          
        } else if( gmmx_data.comp == 4 ) // reflection isomers
        {
        } else // compare all heavy atoms
        {
            maxca = 1;
            for( i = 1; i <= natom; i++ )
            {
                if (atom[i].type != 5)
                {
                    katoms[maxca][0] = i;
                    katoms[maxca][1] = i;
                    maxca++;
                }
            }
            maxca = 0;
            for (i=1; i <= natom; i++)
            {
                if (katoms[i][0] != 0 && katoms[i][1] != 0)
                   maxca++;
            }
            
            if (gmmx_data.comp_method == 0)
            {
               fndumt(maxca,katoms, xyz1, xyz2);
               slow_compare(maxca,katoms,xyz1,xyz2);
            }else
            {
               slow_compare(maxca,katoms,xyz1,xyz2);
            }
            if (nerr == -1)
            {
                drms = 10000.0;
                free_imatrix(katoms, 1,natom+1,0,2);
                return FALSE;
            }
            drms = compute_rms(maxca, katoms, xyz1, xyz2);          
            if (drms > gmmx_data.crmsa)
            {
                gmmx.isitdupl = UNIQUE;
                free_imatrix(katoms, 1,natom+1,0,2);
                return FALSE;
            }else
            {
                gmmx.isitdupl = IDENT;
                free_imatrix(katoms, 1,natom+1,0,2);
                return TRUE;
            }
        }
        return FALSE;
}

//double compute_rms(int maxca, int **katoms, double xyz1[150][3], double xyz2[150][3])
double compute_rms(int maxca, int **katoms, double **xyz1, double **xyz2)
{
    float drms, dave,d, d1, d2, d3, ds;
    double xtmp,ytmp,ztmp;
    int i, ia, ib;
             
    drms = 0.0;
    dave = 0.0;
    d = 0.0;
    d1 = 0.0;
    d2 = 0.0;
    d3 = 0.0;
    for( i = 1; i <= maxca; i++ )
    {
        ia = katoms[i][0];
        ib = katoms[i][1];
        xtmp = xyz1[ia][0];
        ytmp = xyz1[ia][1];
        ztmp = xyz1[ia][2];
        xtmp = xyz2[ib][0];
        ytmp = xyz2[ib][1];
        ztmp = xyz2[ib][2];
        
        d1 = xyz1[katoms[i][0]][0] - xyz2[katoms[i][1]][0];
        d2 = xyz1[katoms[i][0]][1] - xyz2[katoms[i][1]][1];
        d3 = xyz1[katoms[i][0]][2] - xyz2[katoms[i][1]][2];
        d = d1*d1 + d2*d2 + d3*d3;
        drms += d;
        if( d < 0.0 )
            d = -d;
        if( d < 0.00001 || d > 10000.0 )
            ds = 0.0;
        else
        {
            ds = sqrt( d );
            if( fabs( ds ) > gmmx_data.crmsa )
            {
               drms = ds;
               return drms;
            }
        }
        dave += ds;
    }
        
     if( drms < 0.00001 || drms > 10000. )
         drms = 0.0;
     else
         drms = sqrt( drms/(float)( maxca ) );

     return drms;
}
/* -------------------------------------------- */
void fndumt(int npnt,int ** katoms,double **xyz1,double **xyz2)
{
        int  i, j,  k,  n,  norder;
        double b[3][3], d[3], dsq, fac, rmatrx[3][3], rms, rtrmat[3][3], 
         save, savex, savey, savez, umat[3][3], z[3][3];
         
        /*    *** TRANSLATE POINTS XYZ1 AND XYZ2 TO SUPERIMPOSE CENTROIDS
         *       *** ON THE ORIGIN FOR PROPER ROTATION */
        rms = 0.0;
        /*    *** DETERMINE THE R MATRIX (RMATRX) *** */

        for( j = 0; j < 3; j++ )
        {
            for( i = 0; i < 3; i++ )
            {
                        rmatrx[i][j] = 0.0;
            }
        }
        for( n = 1; n <= npnt;  n++ )
        {

                rmatrx[0][0] += 0.1*xyz1[katoms[n][0]][0]*xyz2[katoms[n][1]][0];
                rmatrx[1][0] += 0.1*xyz1[katoms[n][0]][1]*xyz2[katoms[n][1]][0];
                rmatrx[2][0] += 0.1*xyz1[katoms[n][0]][2]*xyz2[katoms[n][1]][0];
                rmatrx[0][1] += 0.1*xyz1[katoms[n][0]][0]*xyz2[katoms[n][1]][1];
                rmatrx[1][1] += 0.1*xyz1[katoms[n][0]][1]*xyz2[katoms[n][1]][1];
                rmatrx[2][1] += 0.1*xyz1[katoms[n][0]][2]*xyz2[katoms[n][1]][1];
                rmatrx[0][2] += 0.1*xyz1[katoms[n][0]][0]*xyz2[katoms[n][1]][2];
                rmatrx[1][2] += 0.1*xyz1[katoms[n][0]][1]*xyz2[katoms[n][1]][2];
                rmatrx[2][2] += 0.1*xyz1[katoms[n][0]][2]*xyz2[katoms[n][1]][2];
        }

        /*    *** DETERMINE THE PRODUCT OF THE TRANSPOSED R MATRIX
         *        WITH THE R MATRIX (RTRMAT) *** */

        for( j = 0; j < 3; j++ )
        {
            for( i = 0; i < 3; i++ )
            {
                rtrmat[j][i] = rmatrx[0][i]*rmatrx[0][j] + rmatrx[1][i]*
                         rmatrx[1][j] + rmatrx[2][i]*rmatrx[2][j];
            }
        }

        /*    *** DETERMINE THE THREE EIGENVALUES (D) AND MUTUALLY ORTHOGONAL
         *        EIGENVECTORS (Z) FROM RTRMAT *** */
        norder = 3;
        nerr = 0;
        eigen( norder, rtrmat, d, z );
        if( nerr == -1 )
                return;
        /*       *** ORDER THE EIGENVALUES IN DESCENDING ORDER. */

        for( i = 0; i < 2; i++ )
        {
            for( j = i + 1; j < 3; j++ )
            {
                if( d[i] < d[j] )
                {
                   save = d[i];
                   d[i] = d[j];
                   d[j] = save;
                   for( k = 0; k < 3; k++ )
                   {
                        save = z[k][i];
                        z[k][i] = z[k][j];
                        z[k][j] = save;
                   }
                }
            }
        }
        
        /*    *** SET Z(3)=Z(1)XZ(2) FOR RIGHT-HANDED SYSTEM *** */

        z[0][2] = z[1][0]*z[2][1] - z[2][0]*z[1][1];
        z[1][2] = z[2][0]*z[0][1] - z[0][0]*z[2][1];
        z[2][2] = z[0][0]*z[1][1] - z[1][0]*z[0][1];

        /*    *** DETERMINE THE NORMALIZED PRODUCTS OF THE R MATRIX WITH THE
         *        THREE EIGENVECTORS Z => VECTORS B *** */

        for( j = 0; j < 2; j++ )
        {
            for( i = 0; i < 3; i++ )
            {
                b[i][j] = rmatrx[i][0]*z[0][j] + rmatrx[i][1]*z[1][j] + 
                         rmatrx[i][2]*z[2][j];
            }
            fac = sqrt( b[0][j]*b[0][j] + b[1][j]*b[1][j] + b[2][j]* b[2][j] );
            if( fabs(fac) < 0.01 )
            {
               rms = -1.0;
               return;
            }
            for( i = 0; i < 3; i++ )
            {
                b[i][j] /= fac;
            }
        }

        /*    *** SET B(3)=B(1)XB(2) *** */

        b[0][2] = b[1][0]*b[2][1] - b[2][0]*b[1][1];
        b[1][2] = b[2][0]*b[0][1] - b[0][0]*b[2][1];
        b[2][2] = b[0][0]*b[1][1] - b[1][0]*b[0][1];

        /*    *** FORM THE OUTPUT MATRIX UMAT *** */

        for( j = 0; j < 3; j++ )
        {
           for( i = 0; i < 3; i++ )
           {
               umat[i][j] = b[i][0]*z[j][0] + b[i][1]*z[j][1] + b[i][2]*z[j][2];
           }
        }

        /*    *** ROTATE XYZ2 ACCORDINGLY AND CALCULATE RMS DEVIATION (RMS) *** */

        if( rms < 0.0 )
                return;
        dsq = 0.0;
        for( n = 1; n <= npnt; n++ )
        {
                savex = xyz2[katoms[n][1]][0];
                savey = xyz2[katoms[n][1]][1];
                savez = xyz2[katoms[n][1]][2];
                xyz2[katoms[n][1]][0] = umat[0][0]*savex + umat[0][1]*savey + umat[0][2]*savez;
                xyz2[katoms[n][1]][1] = umat[1][0]*savex + umat[1][1]*savey + umat[1][2]*savez;
                xyz2[katoms[n][1]][2] = umat[2][0]*savex + umat[2][1]*savey + umat[2][2]*savez;
        }
        return;
} 
/* ----------------------------------------- */
void eigen(int norder,double a[3][3],double *d, double v[3][3])
{
        int  i,  ip,  iq,  j,  nrot;
        double c, g, h, s, sm, t, tau, theta, tresh;
        double *b, *z;

        b = dvector(0,norder);
        z = dvector(0,norder);
        nerr = 0;
        for( ip = 0; ip < norder; ip++ )
        {
           for( iq = 0; iq < norder; iq++ )
           {
               v[ip][iq] = 0.;
           }
           v[ip][ip] = 1.;
        }
        for( ip = 0; ip < norder; ip++ )
        {
                b[ip] = d[ip] = a[ip][ip];
                z[ip] = 0.;
        }
        nrot = 0;
        for( i = 0; i < 50; i++ )
        {
            sm = 0.;
            for( ip = 0; ip < norder - 1; ip++ )
            {
                for( iq = ip + 1; iq < norder; iq++ )
                {
                     sm += fabs( a[ip][iq] );
                }
            }
            if( sm == 0. )
            {
                free_dvector(b,0,norder);
                free_dvector(z,0,norder);               
                return;
            }
            if( i < 3 )
                tresh = 0.2*sm/(norder*norder);
            else
                tresh = 0.;

            for( ip = 0; ip < norder - 1; ip++ )
            {
                for( iq = ip + 1; iq < norder; iq++ )
                {
                    g = 100.*fabs( a[ip][iq] );
                    if( ((i > 3) && (fabs( d[ip] ) + g == fabs( d[ip] ))) && 
                                 (fabs( d[iq] ) + g == fabs( d[iq] )) )
                    {
                        a[ip][iq] = 0.;
                    } else if( fabs( a[ip][iq] ) > tresh )
                    {
                         h = d[iq] - d[ip];
                         if( fabs( h ) + g == fabs( h ) )
                         {
                              t = a[ip][iq]/h;
                         }else
                         {
                              theta = 0.5*h/a[ip][iq];
                              t = 1./(fabs( theta ) + sqrt( 1. + theta*theta ));
                              if( theta < 0. )
                                  t = -t;
                         }
                         c = 1.0/sqrt( 1.0 + t*t );
                         s = t*c;
                         tau = s/(1. + c);
                         h = t*a[ip][iq];
                         z[ip] -= h;
                         z[iq] += h;
                         d[ip] -= h;
                         d[iq] += h;
                         a[ip][iq] = 0.;
                         for( j = 0; j < ip - 1;  j++ )
                         {
                             g = a[j][ip];
                             h = a[j][iq];
                             a[j][ip] = g-s*(h+g*tau);
                             a[j][iq] = h+s*(g-h*tau);
                         }
                         for( j = ip + 1; j < iq - 1; j++ )
                         {
                             g = a[ip][j];
                             h = a[j][iq];
                             a[ip][j] = g-s*(h+g*tau);
                             a[j][iq] = h+s*(g-h*tau);
                          }
                          for( j = iq + 1; j < norder; j++ )
                          {
                             g = a[ip][j];
                             h = a[iq][j];
                             a[ip][j] = g-s*(h+g*tau);
                             a[iq][j] = h+s*(g-h*tau);
                          }
                          for( j = 0; j < norder; j++ )
                          {
                             g = v[j][ip];
                             h = v[j][iq];
                             v[j][ip] = g-s*(h+g*tau);
                             v[j][iq] = h+s*(g-h*tau);
                          }
                          nrot++;
                    }
                }
            }
            for( ip = 0; ip < norder;  ip++ )
            {
                b[ip] += z[ip];
                d[ip] = b[ip];
                z[ip] = 0.;
            }
        }
        nerr = -1;
        free_dvector(b,0,norder);
        free_dvector(z,0,norder);               
        return;
}
/* =======================================================  */
/* ===================================================== */
void sort_file(char *infile, char *outfile)
{
    int i, icount, first, cur_struct,  nconf;
    char inputline[90],dummy1[20], dummy2[20];
    float energy, elowest, eminim;
    float heat;
    FILE *rfile;


    // read input file count structures and
    // read energies
    // write output file starting at lowest energy
    rfile = fopen_path(Openbox.path,infile,"r");
    if (rfile == NULL) 
    {
         return;
    }

    icount = 0;
    nconf = 0;

    while ( FetchRecord(rfile,inputline))
    {
        if (strncasecmp(inputline,"{PCM",4) == 0)
        {
                sscanf(inputline,"%s %s %f %f",dummy1, dummy2, &energy, &heat);
                gmmx.econf[icount] = energy;
                gmmx.hfconf[icount] = heat;
                icount++;
        }
    }
    fclose(rfile);

   strcpy(Openbox.fname,infile);
   strcpy(Savebox.fname,outfile);
   first = TRUE;

    eminim = 1000.0;
    cur_struct = 0;
    while (TRUE)
    {
        elowest = 1000.0;
        for (i=0; i < icount; i++)
        {
            if (gmmx.econf[i] < elowest)
            {
               cur_struct = i;
               elowest = gmmx.econf[i];
            }
        }
        if (elowest < eminim)
           eminim = elowest;
           
        if (elowest == 1000.0)
           break;
        if (elowest < eminim + gmmx_data.ewindow)
        {
           if (first == TRUE)
           {
              files.append = FALSE;
              first = FALSE;
           } else
              files.append = TRUE;
           initialize(); 
           pcmfin(cur_struct+1,0);
           nconf ++;
           sprintf(Struct_Title,"%i %9.3f %9.3f",nconf,gmmx.econf[cur_struct],gmmx.hfconf[cur_struct]);
           pcmfout(1);
           gmmx.econf[cur_struct] = 10000.0;
        } else
             gmmx.econf[cur_struct] = 10000.0;     
    } 
}
/* ================================================  */
void slow_compare(int npts, int **katoms, double **xyz1, double **xyz2)
//void slow_compare(int npts, int **katoms, double xyz1[150][3], double xyz2[150][3])
{
        int nobetter[3];
        int i,j, k, l;
        float a, c, co, cosi, d, dmin, save, si, sign, sine, x1, y1, z1;

	si = co = 0.0;
        nobetter[0] = FALSE;
        nobetter[1] = FALSE;
        nobetter[2] = FALSE;
        a = 1.74532925;
        
        while( TRUE )
        {
          a *= 0.1;
          dmin = 1000.0;
          if( a < 0.001 )
             break;
          for( k = 0; k < 5; k++ )
          {
              /*       *** Y ROTATION *** */
             nobetter[0] = TRUE;
             for( l = 0; l < 2; l++ )
             {
                sign = 1.0;
                if( l == 1 )
                   sign = -1.0;
                for( j = 0; j < 19; j++ )
                {
                   d = 0.0;
                   c = a*(float)(j)*sign;
                   sine = sin(c);
                   cosi = cos(c);
                   for( i = 1; i <= npts; i++ )
                   {
                      x1 = xyz2[katoms[i][1]][0]*cosi + xyz2[katoms[i][1]][2]* sine;
                      z1 = xyz2[katoms[i][1]][2]*cosi - xyz2[katoms[i][1]][0]* sine;

                      d +=  (xyz1[katoms[i][0]][0] - x1)*(xyz1[katoms[i][0]][0] - x1) 
                          + (xyz1[katoms[i][0]][2] - z1)*(xyz1[katoms[i][0]][2] - z1) 
                          + (xyz1[katoms[i][0]][1] - xyz2[katoms[i][1]][1])*
                            (xyz1[katoms[i][0]][1] - xyz2[katoms[i][1]][1]);
                    }
                    if( d < dmin )
                    {
                       nobetter[0] = FALSE;
                       dmin = d;
                       si = sine;
                       co = cosi;
                    }
               }
            }
            if( !nobetter[0] )
            {
               for( i = 1; i <= npts; i++ )
               {
                     save = xyz2[katoms[i][1]][0];
                     xyz2[katoms[i][1]][0] = save*co + xyz2[katoms[i][1]][2]*si;
                     xyz2[katoms[i][1]][2] = xyz2[katoms[i][1]][2]*co - save*si;
               }
            }
            /*       *** X ROTATION *** */
            nobetter[1] = TRUE;
            for( l = 0; l < 2; l++ )
            {
                sign = 1.0;
                if( l == 1 )
                   sign = -1.0;
                for( j = 0; j < 19; j++ )
                {
                   d = 0.0;
                   c = a*(float)(j)*sign;
                   sine = sin(c);
                   cosi = cos(c);
                   for( i = 1; i <= npts; i++ )
                   {
                      z1 = xyz2[katoms[i][1]][2]*cosi + xyz2[katoms[i][1]][1]*sine;
                      y1 = xyz2[katoms[i][1]][1]*cosi - xyz2[katoms[i][1]][2]*sine;

                      d += (xyz1[katoms[i][0]][1] - y1)*(xyz1[katoms[i][0]][1] - y1) 
                         + (xyz1[katoms[i][0]][2] - z1)*(xyz1[katoms[i][0]][2] - z1) 
                         + (xyz1[katoms[i][0]][0] - xyz2[katoms[i][1]][0])*
                           (xyz1[katoms[i][0]][0] - xyz2[katoms[i][1]][0]);
                   }
                   if( d < dmin )
                   {
                      nobetter[1] = FALSE;
                      dmin = d;
                      co = cosi;
                      si = sine;
                    }
                }
             }
             if( !nobetter[1] )
             {
                for( i = 1; i <= npts; i++ )
                {
                     save = xyz2[katoms[i][1]][2];
                     xyz2[katoms[i][1]][2] = save*co + xyz2[katoms[i][1]][1]*si;
                     xyz2[katoms[i][1]][1] = xyz2[katoms[i][1]][1]*co - save*si;
                }
             }
             /*       *** Z ROTATION *** */
             nobetter[2] = TRUE;
             for( l = 0; l < 2; l++ )
             {
                sign = 1.0;
                if( l == 1 )
                  sign = -1.0;
                for( j = 0; j < 19; j++ )
                {
                   d = 0.0;
                   c = a*(float)(j)*sign;
                   sine = sin( c );
                   cosi = cos( c );
                   for( i = 1; i  <= npts; i++ )
                   {
                      y1 = xyz2[katoms[i][1]][1]*cosi + xyz2[katoms[i][1]][0]* sine;
                      x1 = xyz2[katoms[i][1]][0]*cosi - xyz2[katoms[i][1]][1]* sine;

                      d += (xyz1[katoms[i][0]][0] - x1)*(xyz1[katoms[i][0]][0] - x1) 
                         + (xyz1[katoms[i][0]][1] - y1)*(xyz1[katoms[i][0]][1] - y1) 
                         + (xyz1[katoms[i][0]][2] - xyz2[katoms[i][1]][2])*
                           (xyz1[katoms[i][0]][2] - xyz2[katoms[i][1]][2]);
                   }
                   if( d < dmin )
                   {
                      nobetter[2] = FALSE;
                      dmin = d;
                      co = cosi;
                      si = sine;
                    }
                 }
            }
            if( !nobetter[2] )
            {
              for( i = 1; i <= npts; i++ )
              {
                    save = xyz2[katoms[i][1]][1];
                    xyz2[katoms[i][1]][1] = save*co + xyz2[katoms[i][1]][0]*si;
                    xyz2[katoms[i][1]][0] = xyz2[katoms[i][1]][0]*co - save*si;
              }
           }
           if( (nobetter[0] && nobetter[1]) && nobetter[2] )
               break;
        }
    }
        return;
} 
void include_min()
{
    int iftype,i, nrc;

     files.append = FALSE;
     strcpy(Savebox.path,Openbox.path);
     strcpy(Savebox.fname,gmmx.tmpname);
     strcpy(Openbox.fname,gmmx_data.include_name);
     Openbox.ftype = gmmx_data.include_ftype;
     iftype = gmmx_data.include_ftype;
     files.nfiles = 0;
     check_numfile( iftype);

     gmmx_abort = FALSE;

/*  do calculations */
     initialize();
     if (iftype == FTYPE_PCM)
         pcmfin(1,0);
                           
     
     if(gmmx_data.nopi)
       pimarkselection();  
     else
       pireset();    
                   
//  check chirality
     hster(0, jji);
//  check epimers
    if (gmmx_data.chig)  // check for epimerization
       nrc = lepimer(0);
     nrc = bstereo(0);
//  check double bond stereochemistry

     files.batch = True;
     gmmx.nfound = 0;
     gmmx.neminim = 0;
     gmmx.eminim = 1000;
     gmmx.nsaved = 0;
     gmmx.nsflow = 1;
     gmmx.maxdupl = 0;
     
     if (iftype == FTYPE_PCM)
     {
         rdfile(0);
         initialize();
     }
     for (i = 1; i <= files.nfiles; i++)
     {
         initialize();
         if (iftype == FTYPE_PCM)
             rdfile(i);

/* write structure number information  */
        files.icurrent = i;
//  do calculation
     
        if(gmmx_data.nopi)
          pimarkselection();  
        else
          pireset();
          
        set_active();
        generate_bonds();
         
        if (pistuf.norb > 1) initial_picalc();
        minimize();
        dipolemom.xdipole = 0.0;
        dipolemom.ydipole = 0.0;
        dipolemom.zdipole = 0.0;
        if (pot.use_charge) charge_dipole();
        if (pot.use_dipole) dipole_dipole();
        dipolemom.total = sqrt(dipolemom.xdipole*dipolemom.xdipole +
                             dipolemom.ydipole*dipolemom.ydipole +
                             dipolemom.zdipole*dipolemom.zdipole);
        
        gmmx.ecurrent = energies.total;
        gmmx.neminim++;
        nrc = check_id();
        if (nrc == TRUE)
        {
            if (gmmx.ecurrent < (gmmx.eminim+gmmx_data.ewindow2) )
            {
               sprintf(Struct_Title,"%i %9.3f %9.3f",gmmx.nfound + 1,gmmx.ecurrent,get_heat_of_formation());
               if (i == 1)
                  files.append = FALSE;
               else
                  files.append = TRUE;
               optimize_data.final_energy = energies.total;
               optimize_data.final_heat = get_heat_of_formation();
               pcmfout(1);
               gmmx.nfound++;
               gmmx.isitdupl = UNIQUE;
               if (gmmx.ecurrent < gmmx.eminim)
               {
                   gmmx.eminim = gmmx.ecurrent;
                   gmmx.nsflow = gmmx.nsaved;
               }
            }else
            {
                gmmx.isitdupl = UNIQUEHIGH;
            }
        }
        write_gmmx();
        nrc = check_stop();
        if (nrc == FALSE)
            break;
     }
     if (iftype == FTYPE_PCM)
         rdfile(-1);
}
/* ===================================================  */
void write_restart()
{
    char filename[256];
    FILE *wfile;
    int i,j;

    strcpy(filename,gmmx_data.jobname);
    strcat(filename,".rst");

    wfile = fopen_path(Savebox.path,filename,"w");
    // write gmmx_data block
    fprintf(wfile,"%s\n",gmmx_data.jobname);
    fprintf(wfile,"%s\n",gmmx_data.finame);
    fprintf(wfile,"%s\n",gmmx_data.foname);
    fprintf(wfile,"%s\n",gmmx_data.foname2);
    if (strlen(gmmx_data.include_name) == 0)
       fprintf(wfile,"none\n");
    else
      fprintf(wfile,"%s\n",gmmx_data.include_name);
    fprintf(wfile,"%d %d %d %d %d %d\n",gmmx_data.restart,gmmx_data.include_file,
      gmmx_data.include_ftype,gmmx_data.method, gmmx_data.iseed, gmmx_data.its);
    fprintf(wfile,"%d %d %d %d %d %d\n",gmmx_data.lnant,gmmx_data.heat,
      gmmx_data.nopi,gmmx_data.chig, gmmx_data.bad15, gmmx_data.hbond);
    fprintf(wfile,"%d %d %d %d %d\n",gmmx_data.ecut,gmmx_data.kstop,
      gmmx_data.kdup,gmmx_data.kmin, gmmx_data.max_search);
//
    fprintf(wfile,"%d %d\n",gmmx_data.comp, gmmx_data.ncompatoms);
    for (i=0; i < gmmx_data.ncompatoms; i += 20)
       fprintf(wfile,"%d %d %d %d %d %d %d %d %d %d%d %d %d %d %d%d %d %d %d %d\n",
       gmmx_data.comp_atoms[i],gmmx_data.comp_atoms[i+1],gmmx_data.comp_atoms[i+2],gmmx_data.comp_atoms[i+3],gmmx_data.comp_atoms[i+4],
       gmmx_data.comp_atoms[i+5],gmmx_data.comp_atoms[i+6],gmmx_data.comp_atoms[i+7],gmmx_data.comp_atoms[i+8],gmmx_data.comp_atoms[i+9],
       gmmx_data.comp_atoms[i+10],gmmx_data.comp_atoms[i+11],gmmx_data.comp_atoms[i+12],gmmx_data.comp_atoms[i+13],gmmx_data.comp_atoms[i+14],
       gmmx_data.comp_atoms[i+15],gmmx_data.comp_atoms[i+16],gmmx_data.comp_atoms[i+17],gmmx_data.comp_atoms[i+18],gmmx_data.comp_atoms[i+19]);
//
     fprintf(wfile,"%d %d %d %d\n",gmmx_data.nsrms,gmmx_data.nsrmsa,gmmx_data.comp_method,gmmx_data.comp_method2);
     fprintf(wfile,"%d %d %d %d\n",gmmx_data.qpmr,gmmx_data.qdist,gmmx_data.qang,gmmx_data.qdihed);
     fprintf(wfile,"%d %d %d %d\n",gmmx_data.npmr,gmmx_data.ndist,gmmx_data.nang,gmmx_data.ndihed);
//
     if (gmmx_data.npmr > 0)
     {
         for (i=0; i < gmmx_data.npmr; i++)
            fprintf(wfile,"%d %d\n",gmmx_data.pmr_atoms[i][0], gmmx_data.pmr_atoms[i][1]);
     }
//
     if (gmmx_data.ndist > 0)
     {
         for (i=0; i < gmmx_data.ndist; i++)
            fprintf(wfile,"%d %d\n",gmmx_data.dist_atoms[i][0], gmmx_data.dist_atoms[i][1]);
     }
//
     if (gmmx_data.nang > 0)
     {
         for (i=0; i < gmmx_data.nang; i++)
            fprintf(wfile,"%d %d %d\n",gmmx_data.ang_atoms[i][0], gmmx_data.ang_atoms[i][1]
               , gmmx_data.ang_atoms[i][2]);
     }
//
     if (gmmx_data.ndihed > 0)
     {
         for (i=0; i < gmmx_data.ndihed; i++)
            fprintf(wfile,"%d %d %d %d\n",gmmx_data.dihed_atoms[i][0], gmmx_data.dihed_atoms[i][1]
               , gmmx_data.dihed_atoms[i][2], gmmx_data.dihed_atoms[i][3]);
     }
//
     fprintf(wfile,"%d %d %d %d %d\n",gmmx_data.nrings, gmmx_data.rng_size[0],
        gmmx_data.rng_size[1],gmmx_data.rng_size[2],gmmx_data.rng_size[3]);
     if (gmmx_data.nrings > 0)
     {
         for (i=0; i < gmmx_data.nrings; i++)
         {
             for (j=0; j < gmmx_data.rng_size[i]; j++)
                fprintf(wfile,"%d ",gmmx_data.ring_atoms[j][i]);
             fprintf(wfile,"\n");
             fprintf(wfile,"%d %d %d %d %d\n",gmmx_data.ring_resolution[i], gmmx_data.clo_bond[0][i],
                gmmx_data.clo_bond[1][i],gmmx_data.clo_bond[2][i],gmmx_data.clo_bond[3][i]);

             fprintf(wfile,"%5.2f %6.2f %6.2f %d %d %d %d %d %d\n",gmmx_data.clo_distance[i], gmmx_data.clo_ang[0][i],
                gmmx_data.clo_ang[1][i], gmmx_data.clo_bond[4][i], gmmx_data.clo_bond[0][i],gmmx_data.clo_bond[1][i],gmmx_data.clo_bond[2][i],
                gmmx_data.clo_bond[3][i], gmmx_data.clo_bond[5][i]);

             fprintf(wfile,"%d\n",gmmx_data.nring_bonds[i]);
             for (j=0; j < gmmx_data.nring_bonds[i]; j++)
                fprintf(wfile,"%d %d %d\n", gmmx_data.ring_bond_data[i][j][0], gmmx_data.ring_bond_data[i][j][1],
                   gmmx_data.ring_bond_data[i][j][2]);
         }
     }
//
     fprintf(wfile,"%d\n",gmmx_data.nbonds);
     for (i=0; i < gmmx_data.nbonds; i++)
        fprintf(wfile,"%d %d %d %d %d\n", gmmx_data.bond_data[i][0], gmmx_data.bond_data[i][1],
           gmmx_data.bond_data[i][2],gmmx_data.bond_data[i][3],gmmx_data.bond_data[i][4]);
//
     fprintf(wfile,"%f %f %f %f\n",gmmx_data.ewindow,gmmx_data.boltz,gmmx_data.min_contact,
        gmmx_data.ecutoff);
     fprintf(wfile,"%f %f %f %f\n",gmmx_data.bad15_cutoff,gmmx_data.ewindow2,gmmx_data.ermsa,
        gmmx_data.crmsa);
//
//  end of gmmx_data
//
//  start of gmmx block
//
    fprintf(wfile,"%s\n",gmmx.tmpname);
    fprintf(wfile,"%s\n",gmmx.rmsname);
    fprintf(wfile,"%s\n",gmmx.pkmname);
    fprintf(wfile,"%s\n",gmmx.rstname);
    fprintf(wfile,"%d %d %d %d %d\n",gmmx.run,gmmx.namvt,
       gmmx.nrbond,gmmx.smethod,gmmx.hybrida);
    fprintf(wfile,"%d %d %d %d %d\n",gmmx.neminim,gmmx.nfound,
       gmmx.nsfmin,gmmx.nsflow,gmmx.ndupl);
    fprintf(wfile,"%d %d %d %d %d\n",gmmx.maxdupl,gmmx.nsaved,
       gmmx.nlast,gmmx.nlastk,gmmx.job_done);
    fprintf(wfile,"%f %g %g %g\n",gmmx.ebendold,gmmx.eminim,gmmx.ecurrent,
       gmmx.elast);
    for (i=0; i < gmmx.nfound; i++)
       fprintf(wfile,"%d %f %f %f\n",gmmx.nsfound[i],gmmx.econf[i],
         gmmx.hfconf[i],gmmx.dpmconf[i]);
//
     fclose(wfile);
     files.append = FALSE;
     strcpy(Savebox.fname,gmmx_data.finame);
     pcmfout(1);
}                
/* ===================================================  */
int read_restart()
{
    char line[101], dummy[21],dummy1[21],dummy2[21],dummy3[21];
    FILE *wfile;
    int i,j;

    wfile = fopen_path(Openbox.path,Openbox.fname,"r");
    if (wfile == NULL)
    {
        message_alert("Could not read restart file","GMMX Restart Error");
        return FALSE;
    }
    // write gmmx_data block
    fscanf(wfile,"%s",gmmx_data.jobname);
    fscanf(wfile,"%s",gmmx_data.finame);
    fscanf(wfile,"%s",gmmx_data.foname);
    fscanf(wfile,"%s",gmmx_data.foname2);
    fscanf(wfile,"%s",gmmx_data.include_name);
    fscanf(wfile,"%d %d %d %d %d %d",&gmmx_data.restart,&gmmx_data.include_file,
      &gmmx_data.include_ftype,&gmmx_data.method, &gmmx_data.iseed, &gmmx_data.its);
    fscanf(wfile,"%d %d %d %d %d %d",&gmmx_data.lnant,&gmmx_data.heat,
      &gmmx_data.nopi,&gmmx_data.chig, &gmmx_data.bad15, &gmmx_data.hbond);
    fscanf(wfile,"%d %d %d %d %d",&gmmx_data.ecut,&gmmx_data.kstop,
      &gmmx_data.kdup,&gmmx_data.kmin, &gmmx_data.max_search);
//
    fscanf(wfile,"%d %d",&gmmx_data.comp, &gmmx_data.ncompatoms);
    for (i=0; i < gmmx_data.ncompatoms; i += 20)
       fscanf(wfile,"%d %d %d %d %d %d %d %d %d %d%d %d %d %d %d%d %d %d %d %d",
       &gmmx_data.comp_atoms[i],&gmmx_data.comp_atoms[i+1],&gmmx_data.comp_atoms[i+2],&gmmx_data.comp_atoms[i+3],&gmmx_data.comp_atoms[i+4],
       &gmmx_data.comp_atoms[i+5],&gmmx_data.comp_atoms[i+6],&gmmx_data.comp_atoms[i+7],&gmmx_data.comp_atoms[i+8],&gmmx_data.comp_atoms[i+9],
       &gmmx_data.comp_atoms[i+10],&gmmx_data.comp_atoms[i+11],&gmmx_data.comp_atoms[i+12],&gmmx_data.comp_atoms[i+13],&gmmx_data.comp_atoms[i+14],
       &gmmx_data.comp_atoms[i+15],&gmmx_data.comp_atoms[i+16],&gmmx_data.comp_atoms[i+17],&gmmx_data.comp_atoms[i+18],&gmmx_data.comp_atoms[i+19]);
//
     fscanf(wfile,"%d %d %d %d",&gmmx_data.nsrms,&gmmx_data.nsrmsa,&gmmx_data.comp_method,&gmmx_data.comp_method2);
     fscanf(wfile,"%d %d %d %d",&gmmx_data.qpmr,&gmmx_data.qdist,&gmmx_data.qang,&gmmx_data.qdihed);
     fscanf(wfile,"%d %d %d %d",&gmmx_data.npmr,&gmmx_data.ndist,&gmmx_data.nang,&gmmx_data.ndihed);
//
     if (gmmx_data.npmr > 0)
     {
         for (i=0; i < gmmx_data.npmr; i++)
            fscanf(wfile,"%d %d",&gmmx_data.pmr_atoms[i][0], &gmmx_data.pmr_atoms[i][1]);
     }
//
     if (gmmx_data.ndist > 0)
     {
         for (i=0; i < gmmx_data.ndist; i++)
            fscanf(wfile,"%d %d",&gmmx_data.dist_atoms[i][0], &gmmx_data.dist_atoms[i][1]);
     }
//
     if (gmmx_data.nang > 0)
     {
         for (i=0; i < gmmx_data.nang; i++)
            fscanf(wfile,"%d %d %d",&gmmx_data.ang_atoms[i][0], &gmmx_data.ang_atoms[i][1]
               , &gmmx_data.ang_atoms[i][2]);
     }
//
     if (gmmx_data.ndihed > 0)
     {
         for (i=0; i < gmmx_data.ndihed; i++)
            fscanf(wfile,"%d %d %d %d",&gmmx_data.dihed_atoms[i][0], &gmmx_data.dihed_atoms[i][1]
               , &gmmx_data.dihed_atoms[i][2], &gmmx_data.dihed_atoms[i][3]);
     }
//
     fscanf(wfile,"%d %d %d %d %d",&gmmx_data.nrings, &gmmx_data.rng_size[0],
        &gmmx_data.rng_size[1],&gmmx_data.rng_size[2],&gmmx_data.rng_size[3]);
     if (gmmx_data.nrings > 0)
     {
         for (i=0; i < gmmx_data.nrings; i++)
         {
             for (j=0; j < gmmx_data.rng_size[i]; j++)
                fscanf(wfile,"%d ",&gmmx_data.ring_atoms[j][i]);
             FetchRecord(wfile,line);
             sscanf(line,"%d %d %d %d %d",&gmmx_data.ring_resolution[i], &gmmx_data.clo_bond[0][i],
                &gmmx_data.clo_bond[1][i],&gmmx_data.clo_bond[2][i],&gmmx_data.clo_bond[3][i]);
         
             FetchRecord(wfile,line);
             sscanf(line,"%f %f %f %d %d %d %d %d %d",&gmmx_data.clo_distance[i], &gmmx_data.clo_ang[0][i],
                &gmmx_data.clo_ang[1][i], &gmmx_data.clo_bond[4][i], &gmmx_data.clo_bond[0][i],&gmmx_data.clo_bond[1][i],&gmmx_data.clo_bond[2][i],
                &gmmx_data.clo_bond[3][i], &gmmx_data.clo_bond[5][i]);

             FetchRecord(wfile,line);
             sscanf(line,"%d",&gmmx_data.nring_bonds[i]);
             for (j=0; j < gmmx_data.nring_bonds[i]; j++)
                fscanf(wfile,"%d %d %d", &gmmx_data.ring_bond_data[i][j][0], &gmmx_data.ring_bond_data[i][j][1],
                   &gmmx_data.ring_bond_data[i][j][2]);
         }
     }
//
     fscanf(wfile,"%d",&gmmx_data.nbonds);
     for (i=0; i < gmmx_data.nbonds; i++)
        fscanf(wfile,"%d %d %d %d %d", &gmmx_data.bond_data[i][0], &gmmx_data.bond_data[i][1],
           &gmmx_data.bond_data[i][2],&gmmx_data.bond_data[i][3],&gmmx_data.bond_data[i][4]);
//
     fscanf(wfile,"%f %f %f %f",&gmmx_data.ewindow,&gmmx_data.boltz,&gmmx_data.min_contact,
        &gmmx_data.ecutoff);
     fscanf(wfile,"%f %f %f %f",&gmmx_data.bad15_cutoff,&gmmx_data.ewindow2,&gmmx_data.ermsa,
        &gmmx_data.crmsa);
//
//  end of gmmx_data
//
//  start of gmmx block
//
    FetchRecord(wfile,line);
    FetchRecord(wfile,line);
     sscanf(line,"%s",gmmx.tmpname);
    FetchRecord(wfile,line);
     sscanf(line,"%s",gmmx.rmsname);
    FetchRecord(wfile,line);
     sscanf(line,"%s",gmmx.pkmname);
    FetchRecord(wfile,line);
     sscanf(line,"%s",gmmx.rstname);

    FetchRecord(wfile,line);
    sscanf(line,"%d %d %d %d %d",&gmmx.run,&gmmx.namvt,
       &gmmx.nrbond,&gmmx.smethod,&gmmx.hybrida);

    FetchRecord(wfile,line);
    sscanf(line,"%d %d %d %d %d",&gmmx.neminim,&gmmx.nfound,
       &gmmx.nsfmin,&gmmx.nsflow,&gmmx.ndupl);

    FetchRecord(wfile,line);
    sscanf(line,"%d %d %d %d %d",&gmmx.maxdupl,&gmmx.nsaved,
       &gmmx.nlast,&gmmx.nlastk,&gmmx.job_done);
       
    FetchRecord(wfile,line);
    sscanf(line,"%s %s %s %s",dummy, dummy1, dummy2,dummy3);
    gmmx.ebendold = atof(dummy);
    gmmx.eminim = atof(dummy1);
    gmmx.ecurrent = atof(dummy2);
    gmmx.elast = atof(dummy3);
    
    for (i=0; i < gmmx.nfound; i++)
       fscanf(wfile,"%d %f %f %f",&gmmx.nsfound[i],&gmmx.econf[i],
         &gmmx.hfconf[i],&gmmx.dpmconf[i]);
//
     fclose(wfile);
     return TRUE;
}
/* ==================================================  */             
void quatfit(int n1,int **katoms, double **xyz1, double **xyz2)
{
        int i, j, i1, i2;
        double d[4], q[4], rot[3][3], temp1[4], temp2[4], 
         weight, xrot, xxyx, xxyy, xxyz, xyyx, xyyy, xyyz, xzyx, xzyy, 
         xzyz, yrot, zrot;
        double **c,**v;
        double ftemp[4][4];

        c = dmatrix(0,4,0,4);
        v = dmatrix(0,4,0,4);


        /*     build the upper triangle of the quadratic form matrix * */
        xxyx = 0.0e0;
        xxyy = 0.0e0;
        xxyz = 0.0e0;
        xyyx = 0.0e0;
        xyyy = 0.0e0;
        xyyz = 0.0e0;
        xzyx = 0.0e0;
        xzyy = 0.0e0;
        xzyz = 0.0e0;
        for( i = 1; i <= n1; i++ )
        {
                i1 = katoms[i][0];
                i2 = katoms[i][1];
                weight = atom[i1].atomwt;
                xxyx += weight*xyz1[i1][0]*xyz2[i2][0];
                xxyy += weight*xyz1[i1][1]*xyz2[i2][0];
                xxyz += weight*xyz1[i1][2]*xyz2[i2][0];
                xyyx += weight*xyz1[i1][0]*xyz2[i2][1];
                xyyy += weight*xyz1[i1][1]*xyz2[i2][1];
                xyyz += weight*xyz1[i1][2]*xyz2[i2][1];
                xzyx += weight*xyz1[i1][0]*xyz2[i2][2];
                xzyy += weight*xyz1[i1][1]*xyz2[i2][2];
                xzyz += weight*xyz1[i1][2]*xyz2[i2][2];
        }
        c[0][0] = xxyx + xyyy + xzyz;
        c[1][0] = xzyy - xyyz;
        c[1][1] = xxyx - xyyy - xzyz;
        c[2][0] = xxyz - xzyx;
        c[2][1] = xxyy + xyyx;
        c[2][2] = xyyy - xzyz - xxyx;
        c[3][0] = xyyx - xxyy;
        c[3][1] = xzyx + xxyz;
        c[3][2] = xyyz + xzyy;
        c[3][3] = xzyz - xxyx - xyyy;

        /*     diagonalize the quadratic form matrix
         * */
        jacobi( 4, 4, c, d, v, temp1, temp2 );

        /*     extract the desired quaternion
         * */
         for (i=0; i < 4; i++)
         {
            for(j=0; j < 4; j++)
               ftemp[j][i]=v[j][i];
         }
               
        q[0] = v[3][0];
        q[1] = v[3][1];
        q[2] = v[3][2];
        q[3] = v[3][3];

        /*     assemble rotation matrix that superimposes the molecules
         * */
        rot[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
        rot[0][1] = 2.0e0*(q[1]*q[2] - q[0]*q[3]);
        rot[0][2] = 2.0e0*(q[1]*q[3] + q[0]*q[2]);
        rot[1][0] = 2.0e0*(q[2]*q[1] + q[0]*q[3]);
        rot[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
        rot[1][2] = 2.0e0*(q[2]*q[3] - q[0]*q[1]);
        rot[2][0] = 2.0e0*(q[3]*q[1] - q[0]*q[2]);
        rot[2][1] = 2.0e0*(q[3]*q[2] + q[0]*q[1]);
        rot[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];

        /*     rotate second molecule to best fit with first molecule
         * */
        for( i = 1; i <= n1; i++ )
        {
                xrot = xyz2[i][0]*rot[0][0] + xyz2[i][1]*rot[1][0] + xyz2[i][2]*rot[2][0];
                yrot = xyz2[i][0]*rot[0][1] + xyz2[i][1]*rot[1][1] + xyz2[i][2]*rot[2][1];
                zrot = xyz2[i][0]*rot[0][2] + xyz2[i][1]*rot[1][2] + xyz2[i][2]*rot[2][2];
                xyz2[i][0] = xrot;
                xyz2[i][1] = yrot;
                xyz2[i][2] = zrot;
        }
        free_dmatrix(c,0,4,0,4);
        free_dmatrix(v,0,4,0,4);
        return;
}
// ===============================================
double distance(int ia,int ib)
{
    double dx,dy,dz,dist;
    dx = atom[ia].x - atom[ib].x;
    dy = atom[ia].y - atom[ib].y;
    dz = atom[ia].z - atom[ib].z;
    dist = sqrt(dx*dx + dy*dy + dz*dz);
    return(dist);
}
// ===============================================
void initialize_gmmx()
{
    int i, j, k;
    long int mask;

    mask = (1L << 0);
    gmmxring.nrings = 0;
    for (j=0; j < 4; j++)
    {
         gmmx_data.rng_size[j] = 0;
         gmmxring.nring_atoms[j] = 0;
         for (i=0; i <  30; i++)
         {
              gmmx_data.ring_atoms[i][j] = 0 ;
         }
         gmmx_data.clo_distance[j] = 0;
         gmmx_data.clo_ang[0][j] = 0;
         gmmx_data.clo_ang[1][j] = 0;
         gmmx_data.ring_resolution[j] = 0.0;
         gmmx_data.clo_bond[0][j] = 0 ;
         gmmx_data.clo_bond[1][j] = 0 ;
         gmmx_data.clo_bond[2][j] = 0 ;
         gmmx_data.clo_bond[3][j] = 0 ;
         gmmx_data.clo_bond[4][j] = 0 ;
         gmmx_data.clo_bond[5][j] = 0 ;
    }
            
    gmmx_data.method = 3;
    gmmx_data.iseed = 71277;
    gmmx_data.its = 990;
    gmmx_data.lnant = TRUE;
    gmmx_data.ecut = FALSE;
    gmmx_data.chig = TRUE;
    gmmx_data.bad15 = TRUE;
    gmmx_data.nopi = FALSE;
    gmmx_data.hybrid = TRUE;
     for (i=1; i <= natom; i++)
     {
         if (atom[i].flags & mask)
          {
           gmmx_data.nopi = TRUE;
           break;
          }
     }
    gmmx_data.hbond = TRUE;
    gmmx_data.heat = FALSE;
    
    gmmx_data.comp = 0;
    gmmx_data.qpmr = FALSE;
    gmmx_data.qdist = FALSE;
    gmmx_data.qang = FALSE;
    gmmx_data.qdihed = FALSE;
    gmmx_data.npmr = 0;
    gmmx_data.ndist = 0;
    gmmx_data.nang = 0;
    gmmx_data.ndihed = 0;

    gmmx_data.nrings = 0;
    gmmx_data.nbonds = 0;
    gmmx_data.ewindow = 3.5;
    gmmx_data.ewindow2 = 3.0;
    gmmx_data.boltz = 300.0;
    gmmx_data.ecutoff = 0.100;
    gmmx_data.bad15_cutoff = 3.00;
        
    gmmx_data.nsrms = 0;
    gmmx_data.nsrmsa = 0;
    gmmx_data.comp_method = 0;
    gmmx_data.comp_method2 = 1;
    gmmx_data.ermsa = 0.100;
    gmmx_data.crmsa = 0.250;
    gmmx_data.restart = FALSE;
    gmmx_data.include_file = FALSE;

    gmmx_data.max_search = 100000;
    
    gmmx_data.nrings = 0;
    for (i=0; i < 4; i++)
      gmmx_data.nring_bonds[i] = 0; 
    for (i=0; i < 30; i++)
    {
       gmmxring.nring_atoms[i] = 0;
       for (j=0; j < 30; j++)
         gmmxring.ring_atoms[i][j]=0;
       for (j=0; j < 4; j++)
       for (k=0; k < 3; k++)
        gmmx_data.ring_bond_data[j][i][k] = 0;
    }
    
}   
// ==========================================
#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

float ran1(int *idum)
{
    static long ix1,ix2,ix3;
    static float r[98];
    float temp;
    static int iff=0;
    int j;

    if (*idum < 0 || iff == 0)
    {
        iff = 1;
        ix1 = (IC1-(*idum)) %M1;
        ix1 = (IA1*ix1+IC1) %M1;
        ix2 = ix1 % M2;
        ix1 = (IA1*ix1+IC1) %M1;
        ix3 = ix1 % M3;
        for (j=1; j <= 97; j++)
        {
            ix1 = (IA1*ix1+IC1)%M1;
            ix2 = (IA2*ix2+IC2)%M2;
            r[j] = (ix1+ix2*RM2)*RM1;
        }
        *idum = 1;
    }
    ix1 = (IA1*ix1+IC1)%M1;
    ix2 = (IA2*ix2+IC2)%M2;
    ix3 = (IA3*ix3+IC3)%M3;
    j = 1 + ((97*ix3)/M3);
    if (j > 97 || j < 1) printf("Error in ran1\n");
    temp=r[j];
    r[j] = (ix1+ix2*RM2)*RM1;
    return (temp);
}
/* ================================================================= */
int check_bond(int ia1, int ia2, int *ibres)
{
    int i;
    int is_bond;
    int it1, it2, it;

    if (ia1 == 0 || ia2 == 0)
      return FALSE;
    is_bond = FALSE;
    for (i=0; i < MAXIAT; i++)
    {
        if (atom[ia1].iat[i] == ia2)
        {
            is_bond = TRUE;
            break;
        }
    }
    if (is_bond == FALSE)
      return FALSE;

    // need to check that it is not part of a ring
    // need to determine bond resolution
    it1 = atom[ia1].mmx_type;
    it2 = atom[ia2].mmx_type;
    if (it1 <= it2)
      it = it1*100 + it2;
    else
      it = it2*100 + it1;

    if (it == 203 || it == 322)
       *ibres = 30;
    else if (it == 606 || it == 615 || it == 634 ||
              it == 635 || it == 1515 || it == 1534 ||
              it == 3434 || it == 3435 || it == 3535 )
        *ibres = 30;
    else if (it == 309 || it == 306 )
        *ibres = 180;
    else if (it == 202)
        *ibres = 30;
    else if (it == 102)
        *ibres = 30;
    else
        *ibres = 120;
    return TRUE;       
}  
    
