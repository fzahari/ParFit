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

#include "energies.h"
#include "angles.h"
#include "rings.h"
#include "torsions.h"
#include "bonds_ff.h"
#include "derivs.h"
#include "hess.h"
#include "pot.h"
#include "gmmx.h"
#include "utility.h"
#include "field.h"
#include "optimize.h"
#include "pistuf.h"
#include "units.h"
#include "minim_control.h"
#include "dipmom.h"
#include "minim_values.h"
#include "solv.h"
#include "coordb.h"
#include "pcm7.h"

static int Missing_constants = FALSE;
    
//  ===================
void set_missing_constants(int result)
{
  Missing_constants = result;
}
int get_missing_constants()
{
  return (Missing_constants);
}
// =====================
void mmxsub(int ia)
{
    int nret;
    
    minim_control.type = ia;   // 0 to minimize, 1 for single point
    nret = setup_calculation();
    if (nret == FALSE)
    {
         end_calculation();
         return;
    }

    pcm7_min();
    end_calculation();   
}
// ================================================================        
void pcm7_min()
{
  int i;
  double etot;

      if(minim_values.iprint )
      {
         if (minim_control.field == MMX)
            fprintf(pcmoutfile,"\nForce Field: MMX\n");
         else if (minim_control.field == MM3)
            fprintf(pcmoutfile,"\nForce Field: MM3\n");
         else if (minim_control.field == MMFF94)
            fprintf(pcmoutfile,"\nForce Field: MMFF94\n");
         else if (minim_control.field == AMBER)
            fprintf(pcmoutfile,"\nForce Field: AMBER\n");
         else if (minim_control.field == OPLSAA)
            fprintf(pcmoutfile,"\nForce Field: Oplsaa\n");
         fprintf(pcmoutfile,"Structure: %s\n\n",Struct_Title);

         for (i=1; i <= natom; i++)
         {
             fprintf(pcmoutfile,"%4d  %-2s  %4d   %10f %10f %10f %ld\n",i,atom[i].name,atom[i].type,atom[i].x,atom[i].y,atom[i].z,atom[i].flags);
         }
                    
         fprintf(pcmoutfile,"\nInitial Energy : %f\n",energies.total);
         fprintf(pcmoutfile,"Str: %f    Bnd: %f    Tor: %f\n", energies.estr, energies.ebend, energies.etor);
         fprintf(pcmoutfile,"StrBnd: %f VDW: %f     QQ: %f\n",energies.estrbnd,energies.evdw+energies.e14+energies.ehbond, energies.eu);
         fprintf(pcmoutfile,"OOP: %f AA: %f     Strtor: %f\n",energies.eopb,energies.eangang, energies.estrtor);
         fprintf(pcmoutfile,"Urey: %f\n",energies.eurey);
         fprintf(pcmoutfile,"Dielectric Const: %f\n\n",units.dielec);
         if (pot.use_solv)fprintf(pcmoutfile,"Solv: %f\n\n",energies.esolv);
      }

// Hay debug
//         if (field.type == MMX || field.type == MM3) eheat();

      optimize_data.initial_energy = energies.total;
      optimize_data.initial_heat = get_heat_of_formation();
 
      if (minim_control.type == 1)  // single point calculation
         return;

       minimize();
      
      if (pot.use_picalc)
         pirite();
         
      if (minim_values.iprint == TRUE)
         etot = print_energy();
      else
         etot = energy();

//      if (field.type == MMX || field.type == MM3) eheat();

      optimize_data.final_energy = energies.total;
      optimize_data.final_heat = get_heat_of_formation();

//    compute dipole moment here
      dipolemom.xdipole = 0.0;
      dipolemom.ydipole = 0.0;
      dipolemom.zdipole = 0.0;
      if (pot.use_charge || pot.use_bufcharge || use_gast_chrg) charge_dipole();
      if (pot.use_dipole && !use_gast_chrg) dipole_dipole();
      dipolemom.total = sqrt(dipolemom.xdipole*dipolemom.xdipole +
                             dipolemom.ydipole*dipolemom.ydipole +
                             dipolemom.zdipole*dipolemom.zdipole);
      pcmfout(2);
    
      if( minim_values.iprint)
      {
         fprintf(pcmoutfile,"\n\nFinal Energy : %f\n",energies.total);
         fprintf(pcmoutfile,"Str: %f    Bnd: %f    Tor: %f\n", energies.estr, energies.ebend, energies.etor);
         fprintf(pcmoutfile,"StrBnd: %f VDW: %f     QQ: %f\n",energies.estrbnd,energies.evdw+energies.e14+energies.ehbond, energies.eu);
         fprintf(pcmoutfile,"OOP: %f AA: %f     Strtor: %f\n",energies.eopb,energies.eangang, energies.estrtor);
         fprintf(pcmoutfile,"Urey: %f\n",energies.eurey);
         fprintf(pcmoutfile,"Dielectric Const: %f\n\n",units.dielec);
         if (pot.use_solv)fprintf(pcmoutfile,"Solv: %f\n\n",energies.esolv);
         fprintf(pcmoutfile,"Dielectric: %f\n",units.dielec);
      }
}
// ==================================================================
int setup_calculation()
{
    int nRet = 0;
    double etot;
    if (minim_control.field == MM3)
    {
        // copy mm3 types into type file and retype
        pot.use_hbond = FALSE;
        hbondreset();
        hdel(1);
        set_field();
        generate_bonds();
        if (minim_control.added_const)
           get_added_const();
    } else if (minim_control.field == MMFF94)
    {
        // copy mm3 types into type file and retype
        pot.use_hbond = FALSE;
        pot.use_picalc = FALSE;
        hbondreset();
        pireset();
        hdel(1);
        set_field();
        generate_bonds();
        if (minim_control.added_const)
           get_added_const();
    } else
    {
        if (minim_control.added_const)
        {
           get_added_const();
        }
        type();
        set_field();
        generate_bonds();
    }

            
  // allocate memeory for derivatives
  get_memory();
 
  // setup bonds, angles, torsions, improper torsions and angles, allenes
  get_bonds();
  get_angles();
  get_torsions();
  get_coordbonds();

  if (high_coord.ncoord > 0)
      pot.use_highcoord = TRUE;

  attach();

  
  get_rings();
  // need allene

  // set active atoms
  set_active();

  // pi system
  orbital();
/*  assign local geometry potential function parameters  */
    Missing_constants = FALSE;
     errfile = fopen_path(pcwindir,"pcmod.err","w");
     if (pot.use_bond || pot.use_strbnd)  nRet = kbond();
     if (nRet == FALSE) // use end_calc to free memory
     {
         message_alert("Parameters missing. See pcmod.err for information","Error");
         optimize_data.param_avail = 1;
         fclose(errfile);
         energies.total = 9999.;
         return FALSE;
     }
     if (pot.use_angle || pot.use_strbnd || pot.use_angang) kangle();
     if (pot.use_angle || pot.use_opbend || pot.use_opbend_wilson) kopbend();
     if (pot.use_tors  || pot.use_strtor || pot.use_tortor) ktorsion();
     if (pot.use_strbnd) kstrbnd();
     if (pot.use_angang) kangang();
     if (pot.use_strtor) kstrtor();
      
      if (pot.use_coordb) kcoord_bonds();
      if (pot.use_buck || pot.use_hal || pot.use_lj) kvdw();
      if ((pot.use_charge || pot.use_bufcharge || pot.use_ewald) && !use_external_chrg && !use_gast_chrg) kcharge();
      if (pot.use_dipole&& !use_gast_chrg)  kdipole();
      if (pot.use_dipole&& !use_gast_chrg)  kcharge();
      if (pot.use_solv) ksolv();

      if (pot.use_picalc)
       {
           korbit();
           piseq(-1,0);
           piseq(0,0);
           if (pot.use_charge) piden();
       }
       if (Missing_constants == TRUE)
       {
           message_alert("Parameters missing. See pcmod.err for information","Error");
           fclose(errfile);
           optimize_data.param_avail = 1;
           energies.total = 9999.0;
           return (FALSE);
       } else
       {
           fclose(errfile);
           remove("pcmod.err");
       }    
       if (minim_values.iprint == TRUE)
       {
          etot = print_energy();          
       } else
          etot = energy();
       return TRUE;
}
// =================================================================
void end_calculation()
{

    if (pistuf.norb > 0)
        free_orbit();
        
    free_memory();
}
// ==================================================================
void gradient()
{
    int i, j;

      energies.estr = 0.0;
      energies.ebend = 0.0;
      energies.estrbnd = 0.0;
      energies.e14 = 0.0;
      energies.evdw = 0.0;
      energies.etor = 0.0;
      energies.eu = 0.0;
      energies.eopb = 0.0;
      energies.eangang = 0.0;
      energies.estrtor = 0.0;
      energies.ehbond = 0.0;
      energies.efix = 0.0;
      energies.eimprop = 0.0;
      energies.eimptors = 0.0;
      energies.eurey = 0.0;
      energies.esolv = 0.0;

      virial.virx = 0.0;
      virial.viry = 0.0;
      virial.virz = 0.0;
      
      for (i=1; i <= natom; i++)
           atom[i].energy = 0.0;
      
      for (i=1; i <= natom; i++)
      {
        for (j=0; j < 3; j++)
        {
            deriv.deb[i][j] = 0.0;
            deriv.dea[i][j] = 0.0;
            deriv.destb[i][j] = 0.0;
            deriv.deopb[i][j] = 0.0;
            deriv.detor[i][j] = 0.0;
            deriv.de14[i][j] = 0.0;
            deriv.devdw[i][j]= 0.0;
            deriv.deqq[i][j] = 0.0;
            deriv.deaa[i][j] = 0.0;
            deriv.destor[i][j] = 0.0;
            deriv.deimprop[i][j] = 0.0;
            deriv.dehb[i][j] = 0.0;
            deriv.desolv[i][j] = 0.0;
        }
      }
   
      if (pot.use_bond)ebond1();
      if (pot.use_angle)eangle1();
      if (pot.use_angang)eangang1();
      if (pot.use_opbend)eopbend1();
      if (pot.use_opbend_wilson) eopbend_wilson1();
      if (pot.use_tors)etorsion1();
      if (pot.use_strbnd)estrbnd1();
      if (pot.use_coordb)ecoord_bond1();
      if (pot.use_geom) egeom1();
     
      // mmff 
      if (pot.use_hal) ehal1();
      if (pot.use_bufcharge) ebufcharge1();
      // mmx
      if (pot.use_buck && field.type == MM3) ebuckmm31();

      if (pot.use_dipole && !use_gast_chrg)edipole1();
      if (pot.use_dipole && !use_gast_chrg)echrg_dipole1();
      if (pot.use_strtor) estrtor1();
      if (pot.use_solv) esolv1();
      if (pot.use_highcoord) ehigh_coord1();
      
      energies.total =  energies.estr + energies.ebend + energies.etor + energies.estrbnd + energies.e14+
                        energies.evdw + energies.eu + energies.ehbond + energies.eangang + energies.estrtor +
                        energies.eimprop + energies.eimptors + energies.eopb + energies.eurey + energies.esolv;
      for (i=1; i <= natom; i++)
      {
          for (j=0; j < 3; j++)
          {
              deriv.d1[i][j] = deriv.deb[i][j] + deriv.dea[i][j] + deriv.deaa[i][j] +
                               deriv.destb[i][j] + deriv.detor[i][j] + deriv.deopb[i][j] + deriv.dehb[i][j] +
                               deriv.destor[i][j] + deriv.deqq[i][j] + deriv.devdw[i][j] + deriv.de14[i][j] +
                               deriv.deimprop[i][j] + deriv.deub[i][j] + deriv.desolv[i][j];
         }
      }
}
// ==================================================================
double energy()
{
    double etot;
    int i;
//    struct timeb start,end;
    
      energies.total = 0.0;
      energies.estr = 0.0;
      energies.ebend = 0.0;
      energies.etor = 0.0;
      energies.estrbnd = 0.0;
      energies.evdw = 0.0;
      energies.e14 = 0.0;
      energies.ehbond = 0.0;
      energies.eu = 0.0;
      energies.eangang = 0.0;
      energies.estrtor = 0.0;
      energies.eimprop = 0.0;
      energies.eimptors = 0.0;
      energies.eopb = 0.0;
      energies.eurey = 0.0;
      energies.esolv = 0.0;

      for (i=1; i <= natom; i++)
           atom[i].energy = 0.0;

      if (pot.use_bond)ebond();
      if (pot.use_angle)eangle();
      if (pot.use_angang)eangang();
      if (pot.use_opbend)eopbend();
      if (pot.use_opbend_wilson) eopbend_wilson();
      if (pot.use_tors)etorsion();
      if (pot.use_strbnd)estrbnd();
      if (pot.use_geom) egeom();
      
      if (pot.use_coordb)ecoord_bond();
      // mmff 
      if (pot.use_hal) ehal();
      if (pot.use_bufcharge) ebufcharge();
      // mm3
      if (pot.use_buck && field.type == MM3) ebuckmm3();
            
      if (pot.use_dipole && !use_gast_chrg)edipole();
      if (pot.use_dipole && !use_gast_chrg)echrg_dipole();
      if (pot.use_strtor) estrtor();
      if (pot.use_solv) esolv();
      if (pot.use_highcoord) ehigh_coord();
                              
      energies.total =  energies.estr + energies.ebend + energies.etor + energies.estrbnd
             + energies.evdw + energies.e14 + energies.ehbond + energies.eu + energies.eangang +
               energies.estrtor + energies.eimprop + energies.eimptors + energies.eopb + energies.eurey + energies.esolv;
      etot = energies.total;
      return (etot);
}
/* =========================================================================== */ 
double print_energy()
{
    double etot;
    int i;
    
      energies.total = 0.0;
      energies.estr = 0.0;
      energies.ebend = 0.0;
      energies.etor = 0.0;
      energies.estrbnd = 0.0;
      energies.evdw = 0.0;
      energies.e14 = 0.0;
      energies.ehbond = 0.0;
      energies.eu = 0.0;
      energies.eangang = 0.0;
      energies.estrtor = 0.0;
      energies.eimprop = 0.0;
      energies.eimptors = 0.0;
      energies.eopb = 0.0;
      energies.eurey = 0.0;
      energies.esolv = 0.0;

      for (i=1; i <= natom; i++)
           atom[i].energy = 0.0;

      if (pot.use_bond)ebond();
      if (pot.use_angle)eangle();
      if (pot.use_angang)eangang();
      if (pot.use_opbend)eopbend();
      if (pot.use_opbend_wilson) eopbend_wilson();
      if (pot.use_tors)etorsion();
      if (pot.use_strbnd)estrbnd();
      if (pot.use_geom) egeom();
      
      if (pot.use_coordb)ecoord_bond();
      // mmff 
      if (pot.use_hal) ehal();
      if (pot.use_bufcharge) ebufcharge();
      
      // mmx && mm3
      if (pot.use_buck && field.type == MM3) ebuckmm3();      
      if (pot.use_charge || use_gast_chrg) echarge();

      if (pot.use_dipole && !use_gast_chrg)edipole();
      if (pot.use_dipole && !use_gast_chrg)echrg_dipole();
      if (pot.use_strtor) estrtor();
      if (pot.use_solv) esolv();
      if (pot.use_highcoord) ehigh_coord();
      energies.total =  energies.estr + energies.ebend + energies.etor + energies.estrbnd
             + energies.evdw + energies.e14 + energies.ehbond + energies.eu + energies.eangang +
               energies.estrtor + energies.eimprop + energies.eimptors + energies.eopb + energies.eurey + energies.esolv;
      etot = energies.total;
      
      return (etot);
}    

/* ================================================================== */

int iscoord_bond(int i, int j)
{
    int k;
    long int mask;

    mask = 1 << METCOORD_MASK ;

    for (k=0; k < MAXIAT; k++)
    {
        if (atom[i].iat[k] == j && atom[i].bo[k] == 9)
        {
            if (atom[i].type >= 300)
            {
                if (atom[j].type == 2 || atom[j].type == 4 || atom[j].type == 29 ||
                    atom[j].type == 30 || atom[j].type == 40 || atom[j].type == 48 )
                      return TRUE;
            }
            if (atom[j].type >= 300)
            {
                if (atom[i].type == 2 || atom[i].type == 4 || atom[i].type == 29 ||
                    atom[i].type == 30 || atom[i].type == 40 || atom[i].type == 48 )
                      return TRUE;
            }
        }
    }
// throughout Metal-lone pairs if coordinated bond
//  metal to donor atom
    if ( (atom[i].type >= 300) && atom[j].flags & mask)
      return TRUE;
    if ( (atom[j].type >= 300) && atom[i].flags & mask)
      return TRUE;
   
    return FALSE;
}

void get_coordbonds()
{
    int i, j, k, iatype, jatm;

    cbonds.nbnd = 0;
    pot.use_coordb = FALSE;
        
    for (i=1; i <= natom; i++)
    {
        if (atom[i].mmx_type >= 300)
        {
           for (j=0; j <= MAXIAT; j++)
           {
               if (atom[i].bo[j] == 9)  // coord bond
               {
                   iatype = atom[atom[i].iat[j]].mmx_type;
                   if (iatype == 2 || iatype == 4 || iatype == 29 ||
                       iatype == 30 || iatype == 40 || iatype == 48 )
                   {
                      cbonds.icbonds[cbonds.nbnd][0] = i;
                      cbonds.icbonds[cbonds.nbnd][1] = atom[i].iat[j];
                      cbonds.icbonds[cbonds.nbnd][2] = 0;
                      cbonds.nbnd++;
                   } else
                   {
                      jatm = atom[i].iat[j];
                      for (k = 0; k < MAXIAT; k++)
                      {
                          if (atom[atom[jatm].iat[k]].mmx_type == 20) // lp
                          {
                              cbonds.icbonds[cbonds.nbnd][0] = i;
                              cbonds.icbonds[cbonds.nbnd][1] = jatm;
                              cbonds.icbonds[cbonds.nbnd][2] = atom[jatm].iat[k];
                              cbonds.nbnd++;
                          }
                      }
                   }
                   if (cbonds.nbnd > MAXCBND)
                   {
                      fprintf(pcmoutfile,"Error - Too many coordinated bonds!\nProgram will now quit.\n");
                       exit(0);
                   }
               }
           }
        }
    }
    if (cbonds.nbnd > 0)
       pot.use_coordb = TRUE;
}
int isbond(int i, int j)
{
    int k;
    for (k=0; k < MAXIAT; k++)
    {
        if (atom[i].iat[k] == j && atom[i].bo[k] != 9)
            return TRUE;
    }
    return FALSE;
}

void get_bonds()
{
    int i, j;
    long int mask;

    bonds_ff.nbnd = 0;
    mask = 1L << 0;     // flag for pi atoms 
    
    for (i=1; i <= natom; i++)
    {
        for (j=i+1; j <= natom; j++)
        {
           if (isbond(i,j))
           {
              bonds_ff.i12[bonds_ff.nbnd][0] = i;
              bonds_ff.i12[bonds_ff.nbnd][1] = j;
              bonds_ff.nbnd++;
              if (bonds_ff.nbnd > MAXBND)
              {
                 fprintf(pcmoutfile,"Error - Too many bonds!\nProgram will now quit.\n");
                 exit(0);
               }
           }
        }
    }
}

