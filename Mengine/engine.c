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
#include "energies.h"
#include "pot.h"
#include "angles.h"
#include "attached.h"

#include "rings.h"
#include "torsions.h"
#include "nonbond.h"
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "pimatx.h"
#include "field.h"
#include "atom_k.h"
#include "cutoffs.h"
#include "gmmx.h"
#include "fix.h"
#include "minim_values.h"
#include "minim_control.h"
#include "files.h"
#include "solv.h"
#include "user.h"
#include "control.h"
#include "units.h"
#include "engine.h"

#define OPTIMIZE    1
#define SINGLE      2
#define SEARCH      3

       
/* ==================================================== */
int main(int argc, char *argv[])
{
    char filename[80];
    int  i, nheav,nhyd;
    
    initialize_pcmodel(); // set  basic defaults
    initialize();         // set structure defaults
    initialize_control(); // set control file defaults
    initialize_gmmx();    // set search defaults
    // check command line for filename
    if (argc > 1)
       strcpy(filename,argv[1]);
    else
    {
      strcpy(filename,"conpcm");
      //        printf("Enter control filename:");
      //  gets(filename);
    }

//  parse control file
    parse_control(filename);
		
// generate error or logfile name
    if (control.mode == SEARCH)
    {
        strcpy(Openbox.fname,control.infile);
        strcpy(gmmx_data.finame,control.infile);
        strcpy(gmmx_data.foname,control.outfile);
        minim_control.field = control.ffield;
        minim_control.method = control.optmet;
        if (control.addconst)
        {
            minim_control.added_const = TRUE;
            strcpy(minim_control.added_name,control.addparfile);
        }
        if (control.gbsa)
        {
            pot.use_solv = TRUE;
            solvent.type = 1;
        }
        pcmfin(1,0);
//        check_str();
//  decide on charges to use
        if (control.ffield == MM3)
		   minim_values.ndc = 0;
        if (control.getcharge)
		{
			if (control.getcharge == 1)
			{
				if (control.ffield == MM3)
				{
					minim_values.ndc = 4;
				}
			} else if (control.getcharge == 2)
			{
				// use gasteiger
				use_gast_chrg = TRUE;
			} else if (control.getcharge == 3)
			{
				// use external
				use_external_chrg = TRUE;
			} 
		}
//   
	find_rings();
        nheav = 0;
        nhyd = 0;
        for (i=1; i <= natom; i++)
        {
            if (atom[i].atomnum != 1 && atom[i].atomnum != 2)
               nheav++;
            else
               nhyd++;
        }
        if (gmmx_data.kstop == 0)
          gmmx_data.kstop = 5;
        if (gmmx_data.kmin == 0)
          gmmx_data.kmin = 5*nheav;
        if (gmmx_data.kdup == 0)
        {
           gmmx_data.kdup = 5*nheav/2;
           if (gmmx_data.kdup > 50)
              gmmx_data.kdup = 50;
        }
        
        run_gmmx();
    } else if (control.mode == OPTIMIZE || control.mode == SINGLE)
    {
        gmmx.run = FALSE;
        Openbox.ftype = control.inftype;
        Savebox.ftype = control.outftype;
        strcpy(Openbox.fname,control.infile);
        check_numfile(control.inftype,control.infile);
        for (i=1; i <= files.nfiles; i++)
        {
            initialize();
            if (control.inftype == FTYPE_PCM)
              pcmfin(i,0);
            else if (control.inftype == FTYPE_MM3)
              mm32mod(i,0);
            minim_control.field = control.ffield;
            minim_control.method = control.optmet;
            strcpy(Savebox.fname,control.outfile);
            if (control.addconst)
            {
                minim_control.added_const = TRUE;
                strcpy(minim_control.added_name,control.addparfile);
            }
            if (control.gbsa)
            {
                pot.use_solv = TRUE;
                solvent.type = 1;
            }
//  decide on charges to use
        if (control.ffield == MM3)
	       minim_values.ndc = 0;
	if (control.getcharge)
      	{
	     if (control.getcharge == 1)
	      {
       	      	if (control.ffield == MM3)
    		{
		     minim_values.ndc = 4;
		 }
	       } else if (control.getcharge == 2)
       	       {
		    // use gasteiger
		    use_gast_chrg = TRUE;
	     	} else if (control.getcharge == 3)
	      	{
		    // use external
		    use_external_chrg = TRUE;
		} 
	    }
            if (control.mode == OPTIMIZE)
               mmxsub(0);
            else
               mmxsub(1);
            if (i == 1)
               files.append = FALSE;
            else
               files.append = TRUE;
            pcmfout(1);
        }
    }
// quit
    fclose(pcmoutfile);
    exit(0);
    return TRUE;
}
// =======================================
void check_str()
{
    int i, good;
    good = TRUE;
    for (i=1; i <= natom; i++)
    {
        if (atom[i].iat[0] == 0)
        {
            printf("File: %s  atom: %d is not attached to anything\n",Openbox.fname,i);
             good = FALSE;
        }
    }
    if (good == FALSE)
        exit(0);
}
