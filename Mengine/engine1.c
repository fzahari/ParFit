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
#include "gmmx.h"
#include "pot.h"
#include "minim_values.h"
#include "minim_control.h"
#include "units.h"
#include "pistuf.h"
#include "control.h"
#include "pcmfile.h"
#include "engine1.h"
#include "user.h"

#define OPTIMIZE    1
#define SINGLE      2
#define SEARCH      3


#define NL gettoken(); if (pcmfile.head == 1000) goto L_30;

char Savename[80];

// ====================================================
void initialize_control()
{
    control.mode = 0;
    control.inftype = FTYPE_PCM;
    control.outftype = FTYPE_PCM;
    control.nconfout = -1;
    control.ffield = MMX;
    control.iprint = FALSE;
    minim_values.iprint = FALSE;
    control.jprint = FALSE;
    pistuf.jprint = FALSE;
    control.addconst = FALSE;
    control.optmet = 3;
    control.electro = 0;
    control.getcharge = 0;
    control.iuv = 0;
    control.pipl = 1;
    control.gbsa = FALSE;
    strcpy(control.infile,"");
    strcpy(control.outfile,"");
    strcpy(control.addparfile,"");
}
// =====================================
// parser commands
//      mode    Opt, Single, Search
//      infile  input filename - pcm, mm3, sdf
//      outfile output name != inputname, not changed
//      nconfout  - number of conformers to save - default is whatever you get
//      iprint overall printout
//      jprint  pi printout
//      addpar  added parameter file in prm format
//      FF force field
//      optmet - optimization method
//      electrostatics
//      getcharges
//      picontrol       uv , pipl
//      gbsa    integer
//      dielc   float

//      search uses standard gmmx input
//      smethod
//      seed
//      lnant
//      heat
//      dopi
//      chig
//      bad15
//      comptype
//      compmethod
//      hbond
//      stop nemin minconf dup max
//      ewindow
//      ecutoff
//      boltztemp
//      min_contact
//      crms
//      chk_hybrid
//      include
//      numrings, ringsize  standard ring information
//      numbonds standard bond parameters
//      second_cycle
//
//
// =======================================
void parse_control(char *filename)
{
    int sec_cycle,i,j,iz,ftype;
    FILE *infile;
    char dummy[80],dummy1[80],dummy2[80];

    infile = fopen(filename,"r");
    if (infile == NULL)
    {
         printf("Error opening control file: %s\n",filename);
         exit(0);
    }
    sec_cycle = FALSE;
    pcmfile.nocaps = FALSE;
    iz = 0;
//
L_30:   FetchRecord(infile,pcmfile.string);
        *(pcmfile.string + strlen(pcmfile.string) ) = '\0';
        if (feof(infile))
           goto L_DONE;
        pcmfile.head = 1;
//
L_40:
        NL;
        if (strcasecmp(pcmfile.token,"MODE") == 0)
        {
            NL;
            if (strcasecmp(pcmfile.token,"OPT") == 0)
               control.mode = OPTIMIZE;
           else if (strcasecmp(pcmfile.token,"SINGLE") == 0)
               control.mode = SINGLE;
           else if (strcasecmp(pcmfile.token,"SEARCH") == 0)
               control.mode = SEARCH;                    
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"INFILE") == 0 && sec_cycle == FALSE)
        {
            iz = sscanf(pcmfile.string,"%s %s %d",dummy,control.infile,&ftype);
            if (iz == 3)
               control.inftype = ftype;
            goto L_30;
        } else if (strcasecmp(pcmfile.token,"OUTFILE") == 0 && sec_cycle == FALSE)
        {
            sscanf(pcmfile.string,"%s %s %d",dummy,control.outfile,&ftype);
            if (iz == 3)
               control.outftype = ftype;
            goto L_30;
        } else if (strcasecmp(pcmfile.token,"NCONFOUT") == 0)
        {
            NL;
            control.nconfout = atoi(pcmfile.token);
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"IPRINT") == 0)
        {
            NL;  
            control.iprint = atoi(pcmfile.token);
            if (control.iprint == 0)
               minim_values.iprint = FALSE;
            else
               minim_values.iprint = TRUE;
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"JPRINT") == 0)
        {
            NL;
            control.jprint = atoi(pcmfile.token);
            if (control.jprint == 0)
               pistuf.jprint = FALSE;
            else
               pistuf.jprint = TRUE;
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"ADDPAR") == 0)
        {
            sscanf(pcmfile.string,"%s %s",dummy,control.addparfile);
            control.addconst = TRUE;
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"FORCEFIELD") == 0)
        {
            sscanf(pcmfile.string,"%s %s",dummy,pcmfile.token);
            if (strcasecmp(pcmfile.token,"MMX") == 0)
            {
                control.ffield = MMX;
                default_intype = MMX;
                default_outtype = MMX;
                strcpy(dummy,"mmxconst.prm");
                zero_data();
                read_datafiles(dummy);
            } else if (strcasecmp(pcmfile.token,"MM3") == 0)
            {
                control.ffield = MM3;
                default_intype = MM3;
                default_outtype = MM3;
                strcpy(dummy,"mm3.prm");
                zero_data();
                read_datafiles(dummy);
            } else if (strcasecmp(pcmfile.token,"MMFF94") == 0)
            {
                control.ffield = MMFF94;
                default_intype = MMFF94;
                default_outtype = MMFF94;
                strcpy(dummy,"mmff94.prm");
                zero_data();
                read_datafiles(dummy);
            }
	    goto L_40;
        } else if (strcasecmp(pcmfile.token,"OPTMET") == 0)
        {
            NL;
            control.optmet = atoi(pcmfile.token);
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"ELECTROSTATICS") == 0)
        {
            NL;
            control.electro = atoi(pcmfile.token);
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"GETCHARGES") == 0)
        {
            NL;
            control.getcharge = atoi(pcmfile.token);
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"PICONTROL") == 0)
        {
            NL;
            control.iuv = atoi(pcmfile.token);
            NL;
            control.pipl = atoi(pcmfile.token);
            goto L_40;
        } else if (strcasecmp(pcmfile.token,"DIELEC") == 0)
        {
            NL;
            units.dielec = atof(pcmfile.token);
			user.dielec = TRUE;
			NL;
        } else if (strcasecmp(pcmfile.token,"GBSA") == 0)
        {
            NL;
            i = atoi(pcmfile.token);
            if (i > 0)
                control.gbsa = TRUE;
			NL;
        } else if (strcasecmp(pcmfile.token,"METHOD") == 0)
        {
                NL;
                gmmx_data.method = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"SEED") == 0)
        {
                NL;
                gmmx_data.iseed = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"LNANT") == 0)
        {
                NL;
                gmmx_data.lnant = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"HEAT") == 0)
        {
                NL;
                gmmx_data.heat = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"DOPI") == 0)
        {
                NL;
                gmmx_data.nopi = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"CHIG") == 0)
        {
                NL;
                gmmx_data.chig = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"BAD") == 0)
        {
                NL;
                NL;
                gmmx_data.bad15 = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"COMP") == 0)
        {
                NL;
                gmmx_data.comp = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"HBOND") == 0)
        {
                NL;
                gmmx_data.hbond = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"COMPMET") == 0)
        {
                NL;
                if (sec_cycle == FALSE)
                  gmmx_data.comp_method = atoi(pcmfile.token);
                else
                  gmmx_data.comp_method2 = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"STOP") == 0)
        {
                NL;
                gmmx_data.kstop = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"MIN") == 0)
        {
                NL;
                gmmx_data.kmin = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"DUP") == 0)
        {
                NL;
                gmmx_data.kdup = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"MAX") == 0)
        {
                NL;
                gmmx_data.max_search = atoi(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"EWINDOW") == 0)
        {
                NL;
                if (sec_cycle == FALSE)
                  gmmx_data.ewindow = atof(pcmfile.token);
                else
                  gmmx_data.ewindow2 = atof(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"BOLTZMAN") == 0)
        {
                NL;
                gmmx_data.boltz = atof(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"MIN_CONTACT") == 0)
        {
                NL;
                gmmx_data.min_contact = atof(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"ECUTOFF") == 0)
        {
                NL;
//                if (sec_cycle == FALSE)
                  gmmx_data.ecutoff = atof(pcmfile.token);
//                else
//                  gmmx_data.ecutoff2 = atof(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"CRMS") == 0)
        {
                NL;
//                if (sec_cycle == FALSE)
                  gmmx_data.crmsa = atof(pcmfile.token);
//                else
//                  gmmx_data.crmsa2 = atof(pcmfile.token);
                goto L_40;
        } else if (strcasecmp(pcmfile.token,"CHKHYBRID") == 0)
        {
        } else if (strcasecmp(pcmfile.token,"INCLUDE") == 0)
        {
        } else if (strcasecmp(pcmfile.token,"NUMRINGS") == 0)
        {
                NL;
                gmmx_data.nrings = atoi(pcmfile.token);
                for (i=0; i < gmmx_data.nrings; i++)
                {
                    fgets(pcmfile.string,120,infile);
                    sscanf(pcmfile.string,"%s %d\n",dummy1, &gmmx_data.rng_size[i]);
                    for(j=0; j < gmmx_data.rng_size[i]; j++)
                       fscanf(infile,"%d ", &gmmx_data.ring_atoms[j][i]);
                    fgets(pcmfile.string,120,infile);
                    sscanf(pcmfile.string,"%s %d %s %d %d %d %d\n",dummy1, &gmmx_data.ring_resolution[i], dummy2, &gmmx_data.clo_bond[0][i],
                       &gmmx_data.clo_bond[1][i], &gmmx_data.clo_bond[2][i], &gmmx_data.clo_bond[3][i]);

                    fgets(pcmfile.string,120,infile);
                    sscanf(pcmfile.string,"%s %f %s %f %f\n",dummy1, &gmmx_data.clo_distance[i], dummy2, &gmmx_data.clo_ang[0][i],
                       &gmmx_data.clo_ang[1][i]);

                    fgets(pcmfile.string,120,infile);
                    sscanf(pcmfile.string,"%s %d\n",dummy1, &gmmx_data.nring_bonds[i]);
                    for (j=0; j < gmmx_data.nring_bonds[i]; j++)
                    {
                       fgets(pcmfile.string,120,infile);
                       sscanf(pcmfile.string,"%s %d %d %d\n", dummy1, &gmmx_data.ring_bond_data[i][j][0], &gmmx_data.ring_bond_data[i][j][1],
                          &gmmx_data.ring_bond_data[i][j][2]);
                    }
                }
                goto L_30;
        } else if (strcasecmp(pcmfile.token,"NUMBONDS") == 0)
        {
                NL;
                gmmx_data.nbonds = atoi(pcmfile.token);
                for (i=0; i < gmmx_data.nbonds; i++)
                {
                    fgets(pcmfile.string,120,infile);
                    sscanf(pcmfile.string,"%s %d %d %d %d %d\n", dummy1, &gmmx_data.bond_data[i][0], &gmmx_data.bond_data[i][1],
                      &gmmx_data.bond_data[i][2],&gmmx_data.bond_data[i][3],&gmmx_data.bond_data[i][4]);
//Hay debug
/*                  printf("MAIN - read in rot bond %2d", i);
                    printf("  %3d %3d %3d\n", gmmx_data.bond_data[i][0], gmmx_data.bond_data[i][1], gmmx_data.bond_data[i][2]);  */

                }
                goto L_30;
        } else if (strcasecmp(pcmfile.token,"SECOND_CYCLE") == 0)
        {
                sec_cycle = TRUE;
                FetchRecord(infile,pcmfile.string);
                sscanf(pcmfile.string,"%s %s\n",dummy1,gmmx_data.foname);
                FetchRecord(infile,pcmfile.string);
                sscanf(pcmfile.string,"%s %s\n",dummy1,gmmx_data.foname2);
                goto L_30;
        } else
           goto L_40;
// end of read
L_DONE:
    fclose(infile);
    return;
}
// ===============================================
// take infile name, strip off extension and generate infile.err
// open file
void get_errfile()
{
    int iz,i,j,ifound;

    strcpy(control.errorfile,"");
    iz = strlen(control.infile);
    ifound = FALSE;
    for (i=0; i < iz; i++)
    {
        if (control.infile[i] == '.')
        {
            for (j=0; j < i; j++)
            {
                control.errorfile[j] = control.infile[j];
            }
            control.errorfile[j+1] = '\0';
            ifound = TRUE;
            break;
        }
    }
    if (ifound == FALSE)
    {
        strcpy(control.errorfile,control.infile);
        strcat(control.errorfile,".err");
    } else
    {
        strcat(control.errorfile,".err");
    }
    errfile = fopen_path(Openbox.path,control.errorfile,"w");        
}         


