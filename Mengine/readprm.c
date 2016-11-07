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
#include "pot.h"
#include "field.h"
#include "atom_k.h"
#include "readprm.h"
#include "minim_control.h"
#include "minim_values.h"
#include "units.h"
#include "bondk.h"
#include "piatomk.h"
#include "angk.h"       
#include "crossterm.h"
#include "ooplane.h"
#include "tork.h"
#include "vdwk.h"               
#include "chargek.h"
#include "dipolek.h"
#include "eheat.h"        
#include "user.h"

void zero_data()
{
      int ii,j;

      bondk1.use_ring3 = FALSE;
      bondk1.use_ring4 = FALSE;
      bondk1.use_ring5 = FALSE;
      bondpk.use_pibond = FALSE;

      angk1.use_ang3 = FALSE;
      angk1.use_ang4 = FALSE;
      angk1.use_ang5 = FALSE;
      angf.use_angf = FALSE;

      torkn1.use_tor4 = FALSE;
      torkn1.use_tor5 = FALSE;

      strcpy(field.name,"");
      field.type = 0;
      field.bondunit = 0.0;
      field.bond_cubic = 0.0;
      field.bond_quartic = 0.0;
      field.angleunit = 0.0;
      field.angle_cubic = 0.0;
      field.angle_quartic = 0.0;
      field.angle_pentic = 0.0;
      field.angle_sextic = 0.0;
      field.str_bndunit = 0.0;
      field.ang_angunit = 0.0;
      field.str_torunit = 0.0;
      field.torsionunit = 0.0;
      strcpy(field.vdwtype,"");
      strcpy(field.radiusrule,"");
      strcpy(field.radiustype,"");
      strcpy(field.radiussize,"");
      strcpy(field.epsrule,"");
      field.vdw_14scale = 0.0;
      field.a_expterm = 0.0;
      field.b_expterm = 0.0;
      field.c_expterm = 0.0;
      field.chg_14scale = 0.0;
      field.dielectric = 0.0;
      field.pos_12scale = 0.0;
      field.pos_13scale = 0.0;

      units.bndunit = 1.0;
      units.cbnd = 0.0;
      units.qbnd = 0.0;
      units.angunit = 1.0 / (radian*radian);
      units.cang = 0.0;
      units.qang = 0.0;
      units.pang = 0.0;
      units.sang = 0.0;
      units.aaunit = 1.0 / (radian*radian);
      units.stbnunit = 1.0;
      units.ureyunit = 1.0;
      units.torsunit = 1.0;
      units.storunit = 1.0;
      units.v14scale = 1.0;
      units.pos_12scale = 1.0;
      units.pos_13scale = 1.0;
      units.aterm = 0.0;
      units.bterm = 0.0;
      units.cterm = 0.0;
//      units.dielec = 1.0;
      units.chgscale = 1.0;
//  atom type constants
      atom_k.natomtype = 0;
      charge_k.ncharge = 0;
      charge_k.nbndchrg = 0;
      charge_k.nbndchrgdel = 0;
      for (ii = 0; ii < MAXATOMTYPE; ii++)
      {
          atom_k.type[ii]= 0;
          atom_k.valency[ii]= 0;
          atom_k.tclass[ii] = 0;
          atom_k.tclass1[ii] = 0;
          atom_k.tclass2[ii] = 0;
          atom_k.number[ii] = 0;
          atom_k.ligands[ii] = 0;
          atom_k.weight[ii] = 0.00F;
          strcpy(atom_k.symbol[ii],"    ");
          strcpy(atom_k.description[ii],"                  ");
          charge_k.type[ii] = 0;
          charge_k.charge[ii] = 0.00F;
          charge_k.formchrg[ii] = 0.00F;
          charge_k.typechrg[ii] = 0.00F;
      }
//  bond constant data
      bondk1.nbnd = 0;
      vdwpr_k.nvdwpr = 0;
      dipole_k.ndipole = 0;
      electroneg.nelecti = 0;
      electroneg.nelectj = 0;       

      for( ii = 0; ii < MAXBONDCONST; ii++ )
      {
          strcpy(bondk1.kb[ii],"      ");
          bondk1.s[ii] = 0.;
          bondk1.t[ii] = 0.;
          strcpy(vdwpr_k.kv[ii],"      ");
          vdwpr_k.ia1[ii] = 0;
          vdwpr_k.ia2[ii] = 0;
          vdwpr_k.radius[ii] = 0.0;
          vdwpr_k.eps[ii] = 0.0;
          strcpy(dipole_k.kb[ii],"      ");
          dipole_k.bmom[ii] = 0.;
          charge_k.btype[ii] = 0;
          charge_k.bcharge[ii] = 0.0;
          charge_k.btypedel[ii] = 0;
          charge_k.bchargedel[ii] = 0.0;
      }
//  bond 3 constant
      bondk1.nbnd3 = 0;
      dipole_k.ndipole3 = 0;
      for( ii = 0; ii < MAXBOND3CONST; ii++ )
      {
          strcpy(bondk1.kb3[ii],"      ");
          bondk1.s3[ii] = 0.;
          bondk1.t3[ii] = 0.;
          strcpy(dipole_k.kb3[ii],"      ");
          dipole_k.bmom3[ii] = 0.;
      }
// bond 4 constant
      bondk1.nbnd4 = 0;
      dipole_k.ndipole4 = 0;
      for( ii = 0; ii < MAXBOND4CONST; ii++ )
      {
          strcpy(bondk1.kb4[ii],"      ");
          bondk1.s4[ii] = 0.;
          bondk1.t4[ii] = 0.;
          strcpy(dipole_k.kb4[ii],"      ");
          dipole_k.bmom4[ii] = 0.;
      }
// bond 5 constant
      bondk1.nbnd5 = 0;
      dipole_k.ndipole5 = 0;
      for( ii = 0; ii < MAXBOND5CONST; ii++ )
      {
          strcpy(bondk1.kb5[ii],"      ");
          bondk1.s5[ii] = 0.;
          bondk1.t5[ii] = 0.;
          strcpy(dipole_k.kb5[ii],"      ");
          dipole_k.bmom5[ii] = 0.;
      }
// deloc bond  constant
      bondk1.ndeloc = 0;
      for( ii = 0; ii < MAXBONDDELOC; ii++ )
      {
          strcpy(bondk1.kbdel[ii],"      ");
          bondk1.sdel[ii] = 0.;
          bondk1.tdel[ii] = 0.;
      }
// VDW constants
      vdw1.nvdw = 0;
      for( ii = 0; ii < MAXVDWCONST; ii++ )
      {
          vdw1.rad[ii] = 0.;
          vdw1.eps[ii] = 0.;
          vdw1.lpd[ii] = 0;
          vdw1.ihtyp[ii] = 0;
          vdw1.ihdon[ii] = 0;
          vdw1.alpha[ii] = 0.0;
          vdw1.n[ii] = 0.0;
          vdw1.a[ii] = 0.0;
          vdw1.g[ii] = 0.0;
          strcpy(vdw1.da[ii]," ");
      }
//  torsion
      torkn1.ntor = 0;
      for( ii = 0; ii < MAXTORCONST; ii++ )
      {
          strcpy(torkn1.kv[ii],"            ");
          torkn1.tv1[ii] = 0.;
          torkn1.tv2[ii] = 0.;
          torkn1.tv3[ii] = 0.;
          torkn1.tv4[ii] = 0.;
          torkn1.tv5[ii] = 0.;
          torkn1.tv6[ii] = 0.;
          torkn1.phase1[ii] = torkn1.phase2[ii] = torkn1.phase3[ii] = 0;
          torkn1.phase4[ii] = torkn1.phase5[ii] = torkn1.phase6[ii] = 0;
      }
//  torsion 4
      torkn1.ntor4 = 0;
      for( ii = 0; ii < MAXTOR4CONST; ii++ )
      {
          strcpy(torkn1.kv4[ii],"            ");
          torkn1.tv41[ii] = 0.;
          torkn1.tv42[ii] = 0.;
          torkn1.tv43[ii] = 0.;
          torkn1.phase41[ii] = torkn1.phase42[ii] = torkn1.phase43[ii] = 0;
      }
//  torsion 5
      torkn1.ntor5 = 0;
      for( ii = 0; ii < MAXTOR5CONST; ii++ )
      {
          strcpy(torkn1.kv5[ii],"            ");
          torkn1.tv51[ii] = 0.;
          torkn1.tv52[ii] = 0.;
          torkn1.tv53[ii] = 0.;
          torkn1.phase51[ii] = torkn1.phase52[ii] = torkn1.phase53[ii] = 0;
      }
//  torsion delocalized
      torkn1.ntordel = 0;
      for( ii = 0; ii < MAXTORDEL; ii++ )
      {
          strcpy(torkn1.kvdel[ii],"            ");
          torkn1.torindex[ii] = 0;
          torkn1.tvdel1[ii] = 0.;
          torkn1.tvdel2[ii] = 0.;
          torkn1.tvdel3[ii] = 0.;
          torkn1.phasedel1[ii] = torkn1.phasedel2[ii] = torkn1.phasedel3[ii] = 0;
      }
// fourier angles
      angf.nfang = 0;
      for (ii=0; ii < MAXANGCONST; ii++ )
      {
          strcpy(angf.kftype[ii],"         ");
          angf.fcon[ii] = 0.0;
          angf.fc0[ii] = 0.0;
          angf.fc1[ii] = 0.0;
          angf.fc2[ii] = 0.0;
          angf.fc3[ii] = 0.0;
          angf.fc4[ii] = 0.0;
          angf.fc5[ii] = 0.0;
          angf.fc6[ii] = 0.0;
      }
//  angle
      angk1.nang = 0;
      for( ii = 0; ii < MAXANGCONST; ii++ )
      {
          strcpy(angk1.ktype[ii],"         ");
          angk1.con[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.ang[ii][j] = 0.;
      }
//  angle 3
      angk1.nang3 = 0;
      for( ii = 0; ii < MAXANG3CONST; ii++ )
      {
          strcpy(angk1.ktype3[ii],"         ");
          angk1.con3[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.ang3[ii][j] = 0.;
      }
//  angle 4
      angk1.nang4 = 0;
      for( ii = 0; ii < MAXANG4CONST; ii++ )
      {
          strcpy(angk1.ktype4[ii],"         ");
          angk1.con4[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.ang4[ii][j] = 0.;
      }
//  angle 5
      angk1.nang5 = 0;
      for( ii = 0; ii < MAXANG5CONST; ii++ )
      {
          strcpy(angk1.ktype5[ii],"         ");
          angk1.con5[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.ang5[ii][j] = 0.;
      }
//  angle delocalized
      angk1.ndel = 0;
      for( ii = 0; ii < MAXANGDEL; ii++ )
      {
          strcpy(angk1.kdel[ii],"         ");
          angk1.condel[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.angdel[ii][j] = 0.;
      }
//  angle3 delocalized
      angk1.ndel3 = 0;
      for( ii = 0; ii < MAXANG3DEL; ii++ )
      {
          strcpy(angk1.kdel3[ii],"         ");
          angk1.condel3[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.angdel3[ii][j] = 0.;
      }
//  angle4 delocalized
      angk1.ndel4 = 0;
      for( ii = 0; ii < MAXANG4DEL; ii++ )
      {
          strcpy(angk1.kdel4[ii],"         ");
          angk1.condel4[ii] = 0.;
          for( j = 0; j < 3; j++ )
              angk1.angdel4[ii][j] = 0.;
      }
//  out of plane
      ooplane_k.nopbend = 0;
      for( ii = 0; ii < MAXOOP; ii++ )
      {
          strcpy(ooplane_k.iopb[ii],"            ");
          ooplane_k.copb[ii] = 0.;
      }
// cross terms
        crossterm_k.nangang = 0;
        crossterm_k.nstrbnd = 0;
        crossterm_k.nstrtor = 0;
        for (ii = 0; ii < MAXAA; ii++)
        {
            crossterm_k.ang_ang[ii] = 0;
            crossterm_k.aacon[ii][0] = 0.0;
            crossterm_k.aacon[ii][1] = 0.0;
            crossterm_k.aacon[ii][2] = 0.0;
        }
        for (ii = 0; ii < MAXSTBN; ii++)
        {
            strcpy(crossterm_k.stbn[ii],"         ");
            crossterm_k.stbnindex[ii] = 0;
            crossterm_k.stbncon[ii][0] = 0.0;
            crossterm_k.stbncon[ii][1] = 0.0;
            crossterm_k.stbncon[ii][2] = 0.0;
        }
        for (ii = 0; ii < MAXSTRTOR; ii++)
        {
            strcpy(crossterm_k.str_tor[ii],"      ");
            crossterm_k.str_torcon[ii] = 0.0;
        }
            
        
// heat parameters
        ehpara.neheat = 0;
        ehpara.nevhf = 0;
         for( ii = 0; ii < 220; ii++ )
        {
                ehpara.nn[ii] = 0;
                ehpara.ee[ii] = 0.;
                ehpara.estr[ii] = 0.;
                ehpara.et[ii] = 0.;
                strcpy( ehpara.cc[ii], "                " );
        }
        for( ii = 0; ii < 155; ii++ )
        {
                strcpy( ehpara.mm[ii], "                    " );
                strcpy( ehpara.kheat[ii], "            " );
                ehpara.ss[ii] = 0.;
                ehpara.ssless[ii] = 0.;
                ehpara.cent[ii] = 0;
        }
// pi heat parameters
       epiheat.npihf = 0;
        for( ii = 0; ii < 125; ii++ )
        {
                epiheat.nnp[ii] = 0;
                epiheat.eep[ii] = 0.;
                epiheat.aa[ii] = 0.;
                epiheat.bb[ii] = 0.;
                epiheat.ccc[ii] = 0.;
                strcpy( epiheat.ccp[ii], "                    " );
        }
// pi atom data
        piatomk.npiatom = 0;
        for (ii=0; ii < MAXPIATOM; ii++)
        {
            piatomk.kat[ii] = 0;
            piatomk.qp[ii] = 0;
            piatomk.q[ii] = 0.0;
            piatomk.ion[ii] = 0.0;
            piatomk.emz[ii] = 0.0;
            piatomk.zor[ii] = 0.0;
            piatomk.zz[ii] = 0.0;
            piatomk.w1[ii] = 0.0;
            piatomk.w2[ii] = 0.0;
            piatomk.w3[ii] = 0.0;
        }
// pi bond constants
        bondpk.npibond = 0;
        for (ii = 0; ii < MAXPIBONDCONST; ii++)
        {
                strcpy(bondpk.kb[ii],"      ");
                bondpk.bk[ii] = 0.0;
                bondpk.bl[ii] = 0.0;
                bondpk.bmom[ii] = 0.0;
                bondpk.sslop[ii] = 0.0;
                bondpk.tslop[ii] = 0.0;
                bondpk.tslop2[ii] = 0.0;
        }
        
// pi torsions
        torknp.npitor = 0;
        for (ii = 0; ii < MAXPITORCONST; ii++)
        {
                strcpy(torknp.kv[ii],"            ");
                torknp.ph1[ii] = 0;
                torknp.ph2[ii] = 0;
                torknp.ph3[ii] = 0;
                torknp.tv1[ii] = 0.0;
                torknp.tv2[ii] = 0.0;
                torknp.tv3[ii] = 0.0;
        }
// pi angles
        angkp1.npiang = 0;
        for (ii = 0; ii < MAXPIANGCONST; ii++)
        {
                strcpy(angkp1.kpa[ii],"         ");
                angkp1.pacon[ii] = 0.0;
                angkp1.panat[ii][0] = 0.0;
                angkp1.panat[ii][1] = 0.0;
                angkp1.panat[ii][2] = 0.0;
                angkp1.piobp[ii] = 0.0;
        }

}
/* =========================================== */
void get_added_const()
{
    int ia,ib,ic,id, ip1,ip2,ip3,iz,if3,if4,if5;
    float f1,f2,f3,f4;
    float v1[6];
    int se1[6];
    char pa[4],pb[4],pc[4],pd[4],pt[13];
    char line[151], iptemp[21], dummy[10], dummy1[10];
    FILE *wfile;

    wfile = fopen_path(minim_control.added_path,minim_control.added_name,"r");
    if (wfile == NULL)
    {
        message_alert("Error opening added constants file","Error");
        fprintf(pcmoutfile,"Error reading added constants file: %s\n",minim_control.added_name);
        return;
    }

    fprintf(pcmoutfile,"\nThe Following Added Parameters were read:\n\n");
    while ( FetchRecord(wfile,line) )
    {
        sscanf(line,"%s",iptemp);

        if (strcmp(iptemp,"bond") == 0)
        {
            sscanf(line,"%s %d %d %f %f",dummy,&ia,&ib,&f1,&f2);
            fprintf(pcmoutfile,"Bond: AtomTypes %d %d   Blen   %f    Bk %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE)
            {
              iz = check_bondconst(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                bondk1.s[bondk1.nbnd] = f1;
                bondk1.t[bondk1.nbnd] = f2;
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb[bondk1.nbnd],pt);
                bondk1.nbnd++;
              }
            } else
            {
                message_alert("Error reading bond atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"bond5") == 0)
        {
            iz = sscanf(line,"%s %d %d %f %f",dummy,&ia,&ib,&f1,&f2);
            fprintf(pcmoutfile,"Bond5: AtomTypes %d %d   Blen   %f    Bk %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && iz > 3)
            {
              iz = check_bond5const(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                bondk1.s5[bondk1.nbnd5] = f1;
                bondk1.t5[bondk1.nbnd5] = f2;
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb5[bondk1.nbnd5],pt);
                bondk1.nbnd5++;
              }
            } else
            {
                message_alert("Error reading bond atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"bond4") == 0)
        {
            iz = sscanf(line,"%s %d %d %f %f",dummy,&ia,&ib,&f1,&f2);
            fprintf(pcmoutfile,"Bond4: AtomTypes %d %d   Blen   %f    Bk %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && iz > 3)
            {
              iz = check_bond4const(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                bondk1.s4[bondk1.nbnd4] = f1;
                bondk1.t4[bondk1.nbnd4] = f2;
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb4[bondk1.nbnd4],pt);
                bondk1.nbnd4++;
              }
            } else
            {
                message_alert("Error reading bond atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"bond3") == 0)
        {
            sscanf(line,"%s %d %d %f %f",dummy,&ia,&ib,&f1,&f2);
            fprintf(pcmoutfile,"Bond3: AtomTypes %d %d   Blen   %f    Bk %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE)
            {
              iz = check_bond3const(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                bondk1.s3[bondk1.nbnd3] = f1;
                bondk1.t3[bondk1.nbnd3] = f2;
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb3[bondk1.nbnd3],pt);
                bondk1.nbnd3++;
              }
            } else
            {
                message_alert("Error reading bond atom types in added constants file","Error");
            }
                } else if (strcmp(iptemp,"bonddel") == 0)
                {
                } else if (strcmp(iptemp,"angdel") == 0)
                {
        } else if (strcmp(iptemp,"angle") == 0)
        {
            f1 =f2=f3=f4=0.0;
            sscanf(line,"%s %d %d %d %f %f %f %f",dummy,&ia,&ib,&ic,&f1,&f2,&f3,&f4);
            fprintf(pcmoutfile,"Angle: AtomTypes %d %d %d  Anat   %f    Acon %f \n",ia,ib,ic,f2,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_angleconst(ia,ib,ic,f1,f2,f3,f4);
              if (iz == FALSE)
              {
               angk1.con[angk1.nang] = f1;
               angk1.ang[angk1.nang][0] = f2; 
               angk1.ang[angk1.nang][1] = f3; 
               angk1.ang[angk1.nang][2] = f4; 
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
               strcpy(angk1.ktype[angk1.nang],pt);
               angk1.nang++;
              }
            } else
            {
                message_alert("Error reading Angle atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"angle5") == 0)
        {
            f1 =f2=f3=f4=0.0;
            sscanf(line,"%s %d %d %d %f %f %f %f",dummy,&ia,&ib,&ic,&f1,&f2,&f3,&f4);
            fprintf(pcmoutfile,"Angle5: AtomTypes %d %d %d  Anat   %f    Acon %f \n",ia,ib,ic,f2,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_angle5const(ia,ib,ic,f1,f2,f3,f4);
              if (iz == FALSE)
              {
               angk1.con5[angk1.nang5] = f1;
               angk1.ang5[angk1.nang5][0] = f2; 
               angk1.ang5[angk1.nang5][1] = f3; 
               angk1.ang5[angk1.nang5][2] = f4; 
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
               strcpy(angk1.ktype5[angk1.nang5],pt);
               angk1.nang5++;
              }
            } else
            {
                message_alert("Error reading Angle atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"angle4") == 0)
        {
            f1 =f2=f3=f4=0.0;
            sscanf(line,"%s %d %d %d %f %f %f %f",dummy,&ia,&ib,&ic,&f1,&f2,&f3,&f4);
            fprintf(pcmoutfile,"Angle4: AtomTypes %d %d %d  Anat   %f    Acon %f \n",ia,ib,ic,f2,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_angle4const(ia,ib,ic,f1,f2,f3,f4);
              if (iz == FALSE)
              {
               angk1.con4[angk1.nang4] = f1;
               angk1.ang4[angk1.nang4][0] = f2; 
               angk1.ang4[angk1.nang4][1] = f3; 
               angk1.ang4[angk1.nang4][2] = f4; 
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
               strcpy(angk1.ktype4[angk1.nang4],pt);
               angk1.nang4++;
              }
            } else
            {
                message_alert("Error reading Angle atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"angle3") == 0)
        {
            f1 =f2=f3=f4=0.0;
            sscanf(line,"%s %d %d %d %f %f %f %f",dummy,&ia,&ib,&ic,&f1,&f2,&f3,&f4);
            fprintf(pcmoutfile,"Angle3: AtomTypes %d %d %d  Anat   %f    Acon %f \n",ia,ib,ic,f2,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_angle3const(ia,ib,ic,f1,f2,f3,f4);
              if (iz == FALSE)
              {
               angk1.con3[angk1.nang3] = f1;
               angk1.ang3[angk1.nang3][0] = f2; 
               angk1.ang3[angk1.nang3][1] = f3; 
               angk1.ang3[angk1.nang3][2] = f4; 
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
               strcpy(angk1.ktype3[angk1.nang3],pt);
               angk1.nang3++;
              }
            } else
            {
                message_alert("Error reading Angle atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"torsion") == 0)
        {
            sscanf(line,"%s %d %d %d %d %f %d %f %d %f %d",dummy,&ia,&ib,&ic,&id,&f1,&ip1,&f2,
                &ip2,&f3,&ip3);
            fprintf(pcmoutfile,"Torsion: AtomTypes %d %d %d %d  V1   %f    V2 %f    V3 %f\n",ia,ib,ic,id,f1,f2,f3);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE && id < MAXATOMTYPE)
            {
              iz = check_torsionconst(ia,ib,ic,id,f1,f2,f3,ip1,ip2,ip3);
              if (iz == FALSE)
              {
                v1[0] = f1; v1[1] = f2; v1[2] = f3;
                v1[3] = v1[4] = v1[5] = 0.0;
                se1[0] = ip1; se1[1] = ip2; se1[2] = ip3;
                se1[3] = se1[4] = se1[5] = 0;
                torphase(6, v1, se1);
                torkn1.tv1[torkn1.ntor] = v1[0];
                torkn1.phase1[torkn1.ntor] = se1[0];
                torkn1.tv2[torkn1.ntor] = v1[1];
                torkn1.phase2[torkn1.ntor] = se1[1];
                torkn1.tv3[torkn1.ntor] = v1[2];
                torkn1.phase3[torkn1.ntor] = se1[2];
                torkn1.tv4[torkn1.ntor] = v1[3];
                torkn1.phase4[torkn1.ntor] = se1[3];
                torkn1.tv5[torkn1.ntor] = v1[4];
                torkn1.phase5[torkn1.ntor] = se1[4];
                torkn1.tv6[torkn1.ntor] = v1[5];
                torkn1.phase6[torkn1.ntor] = se1[5];
                 numeral(ia,pa,3);
                 numeral(ib,pb,3);
                 numeral(ic,pc,3);
                 numeral(id,pd,3);
                 strcpy(pt,pa);
                 strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 strcpy(torkn1.kv[torkn1.ntor],pt);
                 torkn1.ntor++; 
              }
            } else
            {
                message_alert("Error reading Torsion atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"torsion4") == 0)
        {
            sscanf(line,"%s %d %d %d %d %f %d %f %d %f %d",dummy,&ia,&ib,&ic,&id,&f1,&ip1,&f2,
                &ip2,&f3,&ip3);
            fprintf(pcmoutfile,"Torsion4: AtomTypes %d %d %d %d  V1   %f    V2 %f    V3 %f\n",ia,ib,ic,id,f1,f2,f3);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE && id < MAXATOMTYPE)
            {
              iz = check_torsion4const(ia,ib,ic,id,f1,f2,f3,ip1,ip2,ip3);
              if (iz == FALSE)
              {
                v1[0] = f1; v1[1] = f2; v1[2] = f3;
                v1[3] = v1[4] = v1[5] = 0.0;
                se1[0] = ip1; se1[1] = ip2; se1[2] = ip3;
                se1[3] = se1[4] = se1[5] = 0;
                torphase(6, v1, se1);
                torkn1.tv41[torkn1.ntor4] = v1[0];
                torkn1.phase41[torkn1.ntor4] = se1[0];
                torkn1.tv42[torkn1.ntor4] = v1[1];
                torkn1.phase42[torkn1.ntor4] = se1[1];
                torkn1.tv43[torkn1.ntor4] = v1[2];
                torkn1.phase43[torkn1.ntor4] = se1[2];
                 numeral(ia,pa,3);
                 numeral(ib,pb,3);
                 numeral(ic,pc,3);
                 numeral(id,pd,3);
                 strcpy(pt,pa);
                 strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 strcpy(torkn1.kv4[torkn1.ntor4],pt);
                 torkn1.ntor4++; 
              }
            } else
            {
                message_alert("Error reading Torsion atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"torsion5") == 0)
        {
            sscanf(line,"%s %d %d %d %d %f %d %f %d %f %d",dummy,&ia,&ib,&ic,&id,&f1,&ip1,&f2,
                &ip2,&f3,&ip3);
            fprintf(pcmoutfile,"Torsion5: AtomTypes %d %d %d %d  V1   %f    V2 %f    V3 %f\n",ia,ib,ic,id,f1,f2,f3);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE && id < MAXATOMTYPE)
            {
              iz = check_torsion5const(ia,ib,ic,id,f1,f2,f3,ip1,ip2,ip3);
              if (iz == FALSE)
              {
                v1[0] = f1; v1[1] = f2; v1[2] = f3;
                v1[3] = v1[4] = v1[5] = 0.0;
                se1[0] = ip1; se1[1] = ip2; se1[2] = ip3;
                se1[3] = se1[4] = se1[5] = 0;
                torphase(6, v1, se1);
                torkn1.tv51[torkn1.ntor5] = v1[0];
                torkn1.phase51[torkn1.ntor5] = se1[0];
                torkn1.tv52[torkn1.ntor5] = v1[1];
                torkn1.phase52[torkn1.ntor5] = se1[1];
                torkn1.tv53[torkn1.ntor5] = v1[2];
                torkn1.phase53[torkn1.ntor5] = se1[2];
                 numeral(ia,pa,3);
                 numeral(ib,pb,3);
                 numeral(ic,pc,3);
                 numeral(id,pd,3);
                 strcpy(pt,pa);
                 strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
                 strcpy(torkn1.kv5[torkn1.ntor5],pt);
                 torkn1.ntor5++; 
              }
            } else
            {
                message_alert("Error reading Torsion atom types in added constants file","Error");
            }
                } else if (strcmp(iptemp,"torsiondel") == 0)
                {
                } else if (strcmp(iptemp,"vdwmmff") == 0)
                {
                        sscanf( line, "%s %d %f %f %f %f %s", dummy, &ia, &f1, &f2, &f3, &f4, dummy1);
            fprintf(pcmoutfile,"VDWMMFF: AtomTypes %d   Alpha   %f    Atmp %f Ntmp %f  Gtmp %f Symbol %s\n",ia,f1,f2,f3,f4,dummy1);
                        vdw1.alpha[ia] = f1;
                        vdw1.n[ia] = f2;
                        vdw1.a[ia] = f3;
                        vdw1.g[ia] = f4;
                        strcpy(vdw1.da[ia],dummy1);
        } else if (strcmp(iptemp,"vdw") == 0)
        {
            sscanf(line,"%s %d %f %f %d %d %d",dummy,&ia,&f1,&f2,&if3,&if4,&if5);
            fprintf(pcmoutfile,"VDW: AtomTypes %d   Rad   %f    Eps %f \n",ia,f1,f2);
            vdw1.rad[ia] = f1;
            vdw1.eps[ia] = f2;
            vdw1.lpd[ia] = if3;
            vdw1.ihdon[ia] = if4;
            vdw1.ihtyp[ia] = if5;
        } else if (strcmp(iptemp,"vdwpr") == 0)
        {
            sscanf(line,"%s %d %d %f %f",dummy, &ia, &ib, &f1,&f2);
            fprintf(pcmoutfile,"VDWPR: AtomTypes %d %d  Rad   %f    Eps %f \n",ia,ib,f1,f2);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE)
            {
              iz = check_vdwpr(ia,ib,f1,f2);
              if (iz == FALSE)
              {
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                strcpy(pt,pa); strcat(pt,pb);
                strcpy(vdwpr_k.kv[vdwpr_k .nvdwpr],pt);
                vdwpr_k.ia1[vdwpr_k .nvdwpr] = ia;
                vdwpr_k.ia2[vdwpr_k .nvdwpr] = ib;
                vdwpr_k.radius[vdwpr_k .nvdwpr] = f1;
                vdwpr_k.eps[vdwpr_k .nvdwpr] = f2;
                vdwpr_k .nvdwpr++;
              }
            } else
            {
                message_alert("Error reading VDWPR atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"opbend") == 0)
        {
            sscanf( line, "%s %d %d %d %d %f", dummy, &ia, &ib, &ic,&id,&f1 );
            fprintf(pcmoutfile,"OPBEND: AtomTypes %d %d %d %d  Const  %f \n",ia,ib,ic,id,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_opbend(ia,ib,ic,id,f1);
              if (iz == FALSE)
              {
               numeral(ia,pa,3);
               numeral(ib,pb,3);
               numeral(ic,pc,3);
               numeral(id,pd,3);
               strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc); strcat(pt,pd);
               strcpy(ooplane_k.iopb[ooplane_k.nopbend],pt);
               ooplane_k.nopbend++;
              }
            } else
            {
                message_alert("Error reading Opbend atom types in added constants file","Error");
            }
        } else if (strcmp(iptemp,"strbnd") == 0)
        {
            id = 0;
            sscanf(line,"%s %d %d %d %f %f %f",dummy, &ia, &ib, &ic,
                &f1,&f2,&f3);
            fprintf(pcmoutfile,"STRBEND: AtomTypes %d %d Const  %f %f %f\n",ia,ib,f1,f2,f3);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE && ic < MAXATOMTYPE)
            {
              iz = check_strbnd(ia,ib,ic,f1,f2,f3);
              if (iz == FALSE)
              {
                numeral(ia,pa,3);
                numeral(ib,pb,3);
                numeral(ic,pc,3);
                strcpy(pt,pa); strcat(pt,pb); strcat(pt,pc);
                strcpy(crossterm_k.stbn[crossterm_k.nstrbnd],pt);
                crossterm_k.stbncon[crossterm_k.nstrbnd][0] = f1;
                crossterm_k.stbncon[crossterm_k.nstrbnd][1] = f2;
                crossterm_k.stbncon[crossterm_k.nstrbnd][2] = f3;
                crossterm_k.stbnindex[crossterm_k.nstrbnd] = 0;
                crossterm_k.nstrbnd++;
             }
            } else
            {
                message_alert("Error reading Strbend atom types in added constants file","Error");
            } 
        } else if (strcmp(iptemp,"charge") == 0)
        {
            sscanf( line,"%s %d %f", dummy, &ia, &f1);
            fprintf(pcmoutfile,"Charge  %d  %f\n",ia,f1);
            iz = check_charge(ia,f1);
            if (iz == FALSE)
            {
                charge_k.type[charge_k.ncharge] = ia;
                charge_k.charge[charge_k.ncharge] = f1;
                charge_k.ncharge++;
            }
        } else if (strcmp(iptemp,"mmffchrg") == 0)
        {
                        sscanf( line,"%s %d %f %f %f", dummy, &ia, &f1, &f2, &f3);
                        fprintf(pcmoutfile,"Mmffchrg %d %f %f %f\n ",ia,f1,f2,f3);
                        charge_k.type[ia] = ia;
                        charge_k.charge[ia] = f1;
                        charge_k.formchrg[ia] = f2;
                        charge_k.typechrg[ia] = f3;
    } else if (strcmp(iptemp,"bndchrg") == 0)
    {
                        sscanf(line,"%s %d %d %f",dummy,&ia,&ib,&f1);
                        fprintf(pcmoutfile,"Bndchrg %d %d  %f\n",ia,ib,f1);
                        iz = check_bndchrg(ia,ib,f1);
                        if (iz == FALSE)
                        {
                charge_k.btype[charge_k.nbndchrg] = ia*100+ib;
                                charge_k.bcharge[charge_k.nbndchrg] = f1;
                charge_k.nbndchrg++;
                        }
                } else if (strcmp(iptemp,"bndchrgdel") == 0)
                {
                        sscanf(line,"%s %d %d %f",dummy,&ia,&ib,&f1);
                        fprintf(pcmoutfile,"Bndchrgdel %d %d  %f\n",ia,ib,f1);
                        iz = check_bndchrgdel(ia,ib,f1);
                        if (iz == FALSE)
                        {
                                charge_k.btypedel[charge_k.nbndchrgdel] = ia*100+ib;
                                charge_k.bchargedel[charge_k.nbndchrgdel] = f1;
                                charge_k.nbndchrgdel++;
                        }
        } else if (strcmp(iptemp,"dipole") == 0)
        {
            sscanf(line,"%s %d %d %f",dummy,&ia,&ib,&f1);
            fprintf(pcmoutfile,"DIPOLE: AtomTypes %d %d  Bond Moment  %f \n",ia,ib,f1);
            if (ia < MAXATOMTYPE && ib < MAXATOMTYPE)
            {
              iz = check_dipoleconst(ia,ib,f1);
              if (iz == FALSE)
              {
                 dipole_k.bmom[dipole_k.ndipole] = f1;
                 numeral(ia,pa,3);
                 numeral(ib,pb,3);
                 strcpy(pt,pa); strcat(pt,pb);
                 strcpy(dipole_k.kb[dipole_k.ndipole],pt);
                 dipole_k.ndipole++;
              }
            } else
            {
                message_alert("Error reading Dipole atom types in added constants file","Error");
            }                        
        }
    }
    fclose(wfile);  
}
/* ========================================== */
int check_vdwpr(int ia, int ib,float f1, float f2)
{
    int i;
    char pa[4],pb[4],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < vdwpr_k .nvdwpr; i++)
    {
        if (strcmp(pt,vdwpr_k.kv[i]) == 0)
        {
            vdwpr_k.radius[i] = f1;
            vdwpr_k.eps[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_strbnd(int ia, int ib, int ic,float f1, float f2, float f3)
{
    int i;
    char pa[4],pb[4],pc[4],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < crossterm_k.nstrbnd; i++)
    {
        if (strcmp(pt,crossterm_k.stbn[i]) == 0)
        {
            crossterm_k.stbncon[i][0] = f1;
            crossterm_k.stbncon[i][1] = f2;
            crossterm_k.stbncon[i][2] = f3;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_opbend(int ia,int ib,int ic,int id, float ftemp)
{
    int i;
    char pa[4],pb[4],pc[4],pd[4],pt[13];
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    numeral(id,pd,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    strcat(pt,pd);
    for (i=0; i < ooplane_k.nopbend; i++)
    {
        if (strcmp(pt,ooplane_k.iopb[i]) == 0)
        {
            ooplane_k.copb[i] = ftemp;
            return TRUE;
        }
    }
    return FALSE;
}
// ===================================
int check_charge(int ia,float f1)
{
    int i;
    for (i=0; i < charge_k.ncharge; i++)
    {
        if (ia == charge_k.type[i])
        {
            charge_k.charge[i] = f1;
            return TRUE;
        }
    }
    return FALSE;
}
// ===================================
int check_bndchrg(int ia,int ib,float f1)
{
        int i, ktype;
        
        ktype = ia*100+ib;
    for (i=0; i < charge_k.nbndchrg; i++)
    {
        if (ktype == charge_k.btype[i])
        {        
                        charge_k.bcharge[i]  = f1;         
                        return TRUE;
        }
    }
    return FALSE;
}
// ===================================
int check_bndchrgdel(int ia,int ib,float f1)
{
        int i, ktype;
        
        ktype = ia*100+ib;
    for (i=0; i < charge_k.nbndchrgdel; i++)
    {
        if (ktype == charge_k.btypedel[i])
        {        
                        charge_k.bchargedel[i]  = f1;         
                        return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_dipoleconst(int ia,int ib,float bmom)
{
    int i;
    char pa[4],pb[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa); strcat(pt,pb);
    for (i=0; i < dipole_k.ndipole; i++)
    {
        if (strcmp(pt,dipole_k.kb[i]) == 0)
        {        
           dipole_k.bmom[i] = bmom;
           return TRUE;
        }
    }
    return FALSE;
}  
/* ========================================== */
int check_torsionconst(int ia,int ib,int ic,int id,float v1,float v2, float v3,
        int ip1, int ip2, int ip3)
{
    int i;
    char pa[4],pb[4],pc[4],pd[4],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    numeral(id,pd,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    strcat(pt,pd);
    for (i= 0; i < torkn1.ntor; i++)
    {
        if (strcmp(pt,torkn1.kv[i]) == 0)
        {
            torkn1.tv1[i] = v1;
            torkn1.tv2[i] = v2;
            torkn1.tv3[i] = v3;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_torsion4const(int ia,int ib,int ic,int id,float v1,float v2, float v3,
        int ip1, int ip2, int ip3)
{
    int i;
    char pa[4],pb[4],pc[4],pd[4],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    numeral(id,pd,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    strcat(pt,pd);
    for (i= 0; i < torkn1.ntor4; i++)
    {
        if (strcmp(pt,torkn1.kv4[i]) == 0)
        {
            torkn1.tv41[i] = v1;
            torkn1.tv42[i] = v2;
            torkn1.tv43[i] = v3;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */
int check_torsion5const(int ia,int ib,int ic,int id,float v1,float v2, float v3,
        int ip1, int ip2, int ip3)
{
    int i;
    char pa[4],pb[4],pc[4],pd[4],pt[13];

    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    numeral(id,pd,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    strcat(pt,pd);
    for (i= 0; i < torkn1.ntor5; i++)
    {
        if (strcmp(pt,torkn1.kv5[i]) == 0)
        {
            torkn1.tv51[i] = v1;
            torkn1.tv52[i] = v2;
            torkn1.tv53[i] = v3;
            return TRUE;
        }
    }
    return FALSE;
}
/* ========================================== */    
int check_angleconst(int ia,int ib,int ic,float f1, float f2,float f3,float f4)
{
    int i;
    char pa[4],pb[4],pc[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < angk1.nang; i++)
    {
        if (strcmp(pt,angk1.ktype[i]) == 0)
        {
            angk1.con[i] = f1;
            angk1.ang[i][0] = f2;
            angk1.ang[i][1] = f3;
            angk1.ang[i][2] = f4;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_angle5const(int ia,int ib,int ic,float f1, float f2,float f3,float f4)
{
    int i;
    char pa[4],pb[4],pc[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < angk1.nang5; i++)
    {
        if (strcmp(pt,angk1.ktype5[i]) == 0)
        {
            angk1.con5[i] = f1;
            angk1.ang5[i][0] = f2;
            angk1.ang5[i][1] = f3;
            angk1.ang5[i][2] = f4;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_angle4const(int ia,int ib,int ic,float f1, float f2,float f3,float f4)
{
    int i;
    char pa[4],pb[4],pc[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < angk1.nang4; i++)
    {
        if (strcmp(pt,angk1.ktype4[i]) == 0)
        {
            angk1.con4[i] = f1;
            angk1.ang4[i][0] = f2;
            angk1.ang4[i][1] = f3;
            angk1.ang4[i][2] = f4;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_angle3const(int ia,int ib,int ic,float f1, float f2,float f3,float f4)
{
    int i;
    char pa[4],pb[4],pc[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    numeral(ic,pc,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    strcat(pt,pc);
    for (i=0; i < angk1.nang3; i++)
    {
        if (strcmp(pt,angk1.ktype3[i]) == 0)
        {
            angk1.con3[i] = f1;
            angk1.ang3[i][0] = f2;
            angk1.ang3[i][1] = f3;
            angk1.ang3[i][2] = f4;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_bondconst(int ia,int ib,float f1, float f2)
{
    int i;
    char pa[4],pb[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < bondk1.nbnd; i++)
    {
        if (strcmp(pt,bondk1.kb[i]) == 0)
        {
            bondk1.s[i] = f1;
            bondk1.t[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_bond5const(int ia,int ib,float f1, float f2)
{
    int i;
    char pa[4],pb[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < bondk1.nbnd5; i++)
    {
        if (strcmp(pt,bondk1.kb5[i]) == 0)
        {
            bondk1.s5[i] = f1;
            bondk1.t5[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_bond4const(int ia,int ib,float f1, float f2)
{
    int i;
    char pa[4],pb[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < bondk1.nbnd4; i++)
    {
        if (strcmp(pt,bondk1.kb4[i]) == 0)
        {
            bondk1.s4[i] = f1;
            bondk1.t4[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}         
/* ========================================== */    
int check_bond3const(int ia,int ib,float f1, float f2)
{
    int i;
    char pa[4],pb[4],pt[13];
    
    numeral(ia,pa,3);
    numeral(ib,pb,3);
    strcpy(pt,pa);
    strcat(pt,pb);
    for (i=0; i < bondk1.nbnd3; i++)
    {
        if (strcmp(pt,bondk1.kb3[i]) == 0)
        {
            bondk1.s3[i] = f1;
            bondk1.t3[i] = f2;
            return TRUE;
        }
    }
    return FALSE;
}         
/* =============================================  */
void read_datafiles(char *datafile)
{
    char filename[256];
    char *ev = NULL;

       if ( (ev=getenv("ENGINE_DIR")) != NULL)
    {
           strcpy(filename,ev);
           strcat(filename,"/");
           strcat(filename,datafile);
	} else
          strcpy(filename,datafile);

    read_parameterfile(filename);

    return;
}
// =================================================
void read_parameterfile(char *string)
{
  char line[151], dumm[21], iptemp[21];
  char pa[4],pb[4],pt[7];
  char pc[4],pd[4],pang[10],ptor[13];
  int iz;
  int ihtype,lpde,ihdonor;
  float radius,epsilon,alpha,ntmp,atmp,gtmp;
  float v1[6];
  int   se1[6];
  int iser, inumber, iligand, iclass, iclass1, iclass2, ivalence;
  char symbol[3], descript[25];
  float wtt;
  float charge, fchrg,tchrg;
  int i, iatm1, iatm2, ii, ja, 
    jb, kk, kp1, kp2, kp3, kp4, kt1, kt2, kt3;
  FILE *datafile;
  

  datafile = fopen(string, "rt");
  
  if (datafile == NULL)
    {
      char message[80];
      sprintf(message,"Unable to open the data file %s.\nPlease check that this file was installed\n and is in the same directory as PCWIN",string);
      message_alert(message,"Read Parameter");
      printf("Unable to open the data file %s.\nPlease check that this file was installed\n and is in the same directory as PCWIN\n",string);
      fclose(datafile);
      return;
    }

  strcpy(descript,"");
  strcpy(symbol,"");
  while ( FetchRecord(datafile,line))
    {
      strcpy(iptemp,"");
      strcpy(dumm,"");
      sscanf(line,"%s",iptemp);
      /*  force field descriptors   */
            if (strcmp(iptemp,"forcefield") == 0)
            {
                sscanf(line,"%s %s",dumm, field.name);
                if (strcmp(field.name,"AMBER95") == 0)
                   field.type = AMBER;
                else if (strcmp(field.name,"MMX") == 0)
                   field.type = MMX;
                else if (strcmp(field.name,"MM2-1991") == 0)
                   field.type = MM2;
                else if (strcmp(field.name,"MM2-PLUS") == 0)
                   field.type = MM2;
                else if (strcmp(field.name,"MM3-1992") == 0 || strcmp(field.name,"MM3-1996") == 0 )
                   field.type = MM3;
                else if (strcmp(field.name,"CHARMM22-PROTEIN") == 0)
                   field.type = CHARMM;
                else if (strcmp(field.name,"MMFF94") == 0)
                   field.type = MMFF94;
                else if (strcmp(field.name,"OPLSAA") == 0)
                   field.type = OPLSAA;
                else
                   field.type = UNKNOWN;                   
            } else if (strcmp(iptemp,"bondunit") == 0)
            {
                sscanf(line,"%s %f",dumm, &field.bondunit);
                units.bndunit = field.bondunit;
            } else if (strcmp(iptemp,"bond-cubic") == 0)
            {
                sscanf(line,"%s %f",dumm, &field.bond_cubic);
                units.cbnd = field.bond_cubic;
            } else if (strcmp(iptemp,"bond-quartic") == 0)
            {
                sscanf(line,"%s %f",dumm, &field.bond_quartic);
                units.qbnd = field.bond_quartic;
            } else if (strcmp(iptemp,"angleunit") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.angleunit);
                  units.angunit = field.angleunit;
            } else if (strcmp(iptemp,"angle-cubic") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.angle_cubic);
                  units.cang = field.angle_cubic;
            } else if (strcmp(iptemp,"angle-quartic") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.angle_quartic);
                  units.qang = field.angle_quartic;
            } else if (strcmp(iptemp,"angle-pentic") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.angle_pentic);
                  units.pang = field.angle_pentic;
            } else if (strcmp(iptemp,"angle-sextic") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.angle_sextic);
                  units.sang = field.angle_sextic;
            } else if (strcmp(iptemp,"str-bndunit") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.str_bndunit);
                  units.stbnunit = field.str_bndunit;
            } else if (strcmp(iptemp,"ang-angunit") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.ang_angunit);
                  units.aaunit = field.ang_angunit;
            } else if (strcmp(iptemp,"torsionunit") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.torsionunit);
                  units.torsunit = field.torsionunit;
            } else if (strcmp(iptemp,"str-torunit") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.str_torunit);
                  units.storunit = field.str_torunit;
            } else if (strcmp(iptemp,"vdwtype") == 0)
            {
                  sscanf(line,"%s %s",dumm, field.vdwtype);
                  if (strcmp(field.vdwtype ,"LENNARD-JONES") == 0)
                    pot.use_lj = TRUE;
                  else if (strcmp(field.vdwtype ,"BUCKINGHAM") == 0)
                    pot.use_buck = TRUE;
                  else if (strcmp(field.vdwtype ,"BUFFERED-14-7") == 0)
                    pot.use_hal = TRUE;
                  else if (strcmp(field.vdwtype ,"BGAUSSIAN") == 0)
                    pot.use_gauss = TRUE;
            } else if (strcmp(iptemp,"radiustype") == 0)
            {
                  sscanf(line,"%s %s",dumm, field.radiustype);
            } else if (strcmp(iptemp,"radiussize") == 0)
            {
                  sscanf(line,"%s %s",dumm, field.radiussize);
            } else if (strcmp(iptemp,"radiusrule") == 0)
            {
                  sscanf(line,"%s %s",dumm, field.radiusrule);
            } else if (strcmp(iptemp,"epsilonrule") == 0)
            {
                  sscanf(line,"%s %s",dumm, field.epsrule);
            } else if (strcmp(iptemp,"a-expterm") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.a_expterm);
                  units.aterm = field.a_expterm;
            } else if (strcmp(iptemp,"b-expterm") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.b_expterm);
                  units.bterm = field.b_expterm;
            } else if (strcmp(iptemp,"c-expterm") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.c_expterm);
                  units.cterm = field.c_expterm;
            } else if (strcmp(iptemp,"vdw-14-scale") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.vdw_14scale);
                  units.v14scale = field.vdw_14scale;
            } else if (strcmp(iptemp,"pos-12scale") == 0)
            {
		  sscanf(line,"%s %f",dumm, &field.pos_12scale);
       		  units.pos_12scale = field.pos_12scale;
            } else if (strcmp(iptemp,"pos-13scale") == 0)
            {
      		  sscanf(line,"%s %f",dumm, &field.pos_13scale);
       		  units.pos_13scale = field.pos_13scale;
            } else if (strcmp(iptemp,"dielectric") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.dielectric);
                 if (user.dielec == FALSE)
                  units.dielec = field.dielectric;
            } else if (strcmp(iptemp,"chg-14-scale") == 0)
            {
                  sscanf(line,"%s %f",dumm, &field.chg_14scale);
                  units.chgscale = field.chg_14scale;
            } else if (strcmp(iptemp,"atom") == 0)    /* atom type descriptors */
            {
                if (atom_k.natomtype >= MAXATOMTYPE )
                {
                     message_alert("Error - Too many atom types","Read Parameter");
                     printf("Maximum number of atom types exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                }
                iser = -1;
		strcpy(descript,"");
		strcpy(symbol,"");
                sscanf(line, "%s %d %s %22c %d %f %d %d %d %d %d",dumm,&iser, symbol, descript, &inumber, &wtt,
                              &iligand, &ivalence, &iclass, &iclass1, &iclass2);
		if (iser != -1)
		  {
		    atom_k.type[iser] = iser;
		    strcpy(atom_k.symbol[iser], symbol);
		    descript[20] = '\0';
		    strcpy(atom_k.description[iser], descript);
		    //atom_k.description[iser][20] = '\0';
		    atom_k.number[iser]  = inumber;
		    atom_k.weight[iser]  = wtt;
		    atom_k.ligands[iser] = iligand;
		    atom_k.tclass[iser]   = iclass;
		    atom_k.valency[iser] = ivalence;
		    atom_k.tclass1[iser]  = iclass1;
		    atom_k.tclass2[iser]  = iclass2;
		    atom_k.natomtype++;
		  }
            } else if( strcmp(iptemp,"bond") == 0 )    /* bond constants  */
            {
                if (bondk1.nbnd >= MAXBONDCONST)
                {
                     message_alert("Error - Too many bond constants","Read Parameter");
                     printf("Maximum number of bond constants exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.s[bondk1.nbnd], &bondk1.t[bondk1.nbnd]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb[bondk1.nbnd],pt);
                bondk1.nbnd++;
            } else if( strcmp(iptemp,"bond3") == 0 )    /* bond constants  */
            {
                if (bondk1.nbnd3 >= MAXBOND3CONST)
                {
                     message_alert("Error - Too many bond3 constants","Read Parameter");
                     printf("Maximum number of bond3 constants exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.s3[bondk1.nbnd3], &bondk1.t3[bondk1.nbnd3]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb3[bondk1.nbnd3],pt);
//                bondk1.kb3[bondk1.nbnd3] = ii*100 + kk;
                bondk1.nbnd3++;
            } else if( strcmp(iptemp,"bond4") == 0 )    /* bond constants  */
            {
                if (bondk1.nbnd4 >= MAXBOND4CONST)
                {
                     message_alert("Error - Too many bond4 constants","Read Parameter");
                     printf("Maximum number of bond4 constants exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.s4[bondk1.nbnd4], &bondk1.t4[bondk1.nbnd4]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb4[bondk1.nbnd4],pt);
//                bondk1.kb4[bondk1.nbnd4] = ii*100 + kk;
                bondk1.nbnd4++;
            } else if( strcmp(iptemp,"bond5") == 0 )    /* bond constants  */
            {
                if (bondk1.nbnd5>= MAXBOND5CONST)
                {
                     message_alert("Error - Too many bond5 constants","Read Parameter");
                     printf("Maximum number of bond5 constants exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.s5[bondk1.nbnd5], &bondk1.t5[bondk1.nbnd5]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kb5[bondk1.nbnd5],pt);
                bondk1.nbnd5++;
            } else if( strcmp(iptemp,"bonddel") == 0 )    /* bond constants  */
            {
                if (bondk1.ndeloc>= MAXBONDDELOC)
                {
                     message_alert("Error - Too many delocalized bond constants","Read Parameter");
                     printf("Maximum number of delocalized bond constants exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                }
                        
                sscanf( line, "%s %d %d %f %f", dumm, &ii, &kk, 
                         &bondk1.sdel[bondk1.ndeloc], &bondk1.tdel[bondk1.ndeloc]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondk1.kbdel[bondk1.ndeloc],pt);
                bondk1.ndeloc++;
            } else if (strcmp(iptemp, "electroi") == 0)
            {
                if (electroneg.nelecti > 200)
                {
                     message_alert("Error - Too many electronegativity constants","Read Parameter");
                     printf("Maximum number of electronegativity constants exceeded. Others will be ignored\n");
                     return;
                }
                sscanf(line,"%s %d %d %d %f",dumm,&electroneg.itype[electroneg.nelecti],
                       &electroneg.ibond[electroneg.nelecti],
                       &electroneg.iattach[electroneg.nelecti],
                       &electroneg.icorr[electroneg.nelecti]);
                electroneg.nelecti++;
            } else if (strcmp(iptemp, "electroj") == 0)
            {
                if (electroneg.nelectj > 200)
                {
                     message_alert("Error - Too many electronegativity constants","Read Parameter");
                     printf("Maximum number of electronegativity constants exceeded. Others will be ignored\n");
                     return;
                }
                sscanf(line,"%s %d %d %d %f",dumm,&electroneg.jtype[electroneg.nelectj],
                       &electroneg.jbond[electroneg.nelectj],
                       &electroneg.jattach[electroneg.nelectj],
                       &electroneg.jcorr[electroneg.nelectj]);
                electroneg.nelectj++;
            } else if (strcmp(iptemp, "piatom") == 0)
            {
                if (piatomk.npiatom > MAXPIATOM)
                {
                     message_alert("Error - Too many piatom constants","Read Parameter");
                     printf("Maximum number of piatom constants exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                }
                sscanf(line,"%s %d %f %f %f %f %f %f %f %f %d", dumm, &piatomk.kat[piatomk.npiatom],
                  &piatomk.q[piatomk.npiatom],  &piatomk.ion[piatomk.npiatom],  &piatomk.emz[piatomk.npiatom],
                  &piatomk.zor[piatomk.npiatom],&piatomk.zz[piatomk.npiatom],  &piatomk.w1[piatomk.npiatom],
                  &piatomk.w2[piatomk.npiatom],  &piatomk.w3[piatomk.npiatom],&piatomk.qp[piatomk.npiatom]);
                piatomk.npiatom++;
            } else if (strcmp(iptemp, "pibond") == 0)
            {
                if (bondpk.npibond > MAXPIBONDCONST )
                {
                }
                bondpk.use_pibond = TRUE;
               sscanf(line,"%s %d %d %f %f %f %f %f %f",dumm, &ii, &kk,
                &bondpk.bk[bondpk.npibond], &bondpk.bl[bondpk.npibond], &bondpk.bmom[bondpk.npibond],
                &bondpk.sslop[bondpk.npibond], &bondpk.tslop[bondpk.npibond], &bondpk.tslop2[bondpk.npibond]);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa);
                strcat(pt,pb);
                strcpy(bondpk.kb[bondpk.npibond],pt);
                bondpk.npibond++;
            } else if (strcmp(iptemp,"piangle") == 0)
            {
                if (angkp1.npiang > MAXPIANGCONST)
                {
                }
                sscanf(line,"%s %d %d %d %f %f %f %f",dumm,&kt1,&kt2,&kt3,&angkp1.pacon[angkp1.npiang],
                    &angkp1.panat[angkp1.npiang][0],&angkp1.panat[angkp1.npiang][1],&angkp1.panat[angkp1.npiang][1]);
                numeral(kt1,pa,3);
                numeral(kt2,pb,3);
                numeral(kt3,pc,3);
                strcpy(pang,pa);
                strcat(pang,pb); strcat(pang,pc);
                strcpy(angkp1.kpa[angkp1.npiang],pang);
                angkp1.npiang++;
            } else if (strcmp(iptemp,"dipole") == 0)
            {
                if (dipole_k.ndipole > MAXBONDCONST)
                {
                    message_alert("Error - Too many bond dipoles","Read Parameter");
                    printf("Maximum number of bond dipole constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                sscanf( line,"%s %d %d %f",dumm, &ii, &kk, &radius);
                numeral(ii,pa,3);
                numeral(kk,pb,3);
                strcpy(pt,pa); strcat(pt,pb);
                strcpy(dipole_k.kb[dipole_k.ndipole],pt);
                dipole_k.bmom[dipole_k.ndipole] = radius;
                dipole_k.ndipole++;
            } else if (strcmp(iptemp, "charge") == 0)
            {
                if (charge_k.ncharge > MAXATOMTYPE)
                {
                    message_alert("Error - Too many atom charges","Read Parameter");
                    printf("Maximum number of atom charge constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                sscanf( line,"%s %d %f", dumm, &charge_k.type[charge_k.ncharge], &charge_k.charge[charge_k.ncharge]);
                charge_k.ncharge++;
            } else if (strcmp(iptemp, "mmffchrg") == 0)
            {
                if (charge_k.ncharge > MAXATOMTYPE)
                {
                    message_alert("Error - Too many atom charges","Read Parameter");
                    printf("Maximum number of atom charge constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                sscanf( line,"%s %d %f %f %f", dumm, &i, &charge, &fchrg, &tchrg);
                 charge_k.type[i] = i;
                 charge_k.charge[i] = charge;
                 charge_k.formchrg[i] = fchrg;
                 charge_k.typechrg[i] = tchrg;
            } else if (strcmp(iptemp,"bndchrgdel") == 0)
            {
                if (charge_k.nbndchrgdel > MAXBONDCONST)
                {
                    message_alert("Error - Too many bond charges","Read Parameter");
                    printf("Maximum number of bond charge constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                sscanf(line,"%s %d %d %f",dumm,&kp1,&kp2,&charge_k.bchargedel[charge_k.nbndchrgdel]);
                charge_k.btypedel[charge_k.nbndchrgdel] = kp1*100+kp2;
                charge_k.nbndchrgdel++;
            } else if (strcmp(iptemp,"bndchrg") == 0)
            {
                if (charge_k.nbndchrg > MAXBONDCONST)
                {
                    message_alert("Error - Too many bond charges","Read Parameter");
                    printf("Maximum number of bond charge constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                sscanf(line,"%s %d %d %f",dumm,&kp1,&kp2,&charge_k.bcharge[charge_k.nbndchrg]);
                charge_k.btype[charge_k.nbndchrg] = kp1*100+kp2;
                charge_k.nbndchrg++;
            }  else if (strcmp(iptemp,"vdwmmff") == 0)
            {
                if (vdw1.nvdw > MAXVDWCONST)
                {
                    message_alert("Error - Too many vdw constants","Read Parameter");
                    printf("Maximum number of vdw constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                sscanf( line, "%s %d %f %f %f %f %s", dumm, &i, &alpha, &ntmp, &atmp, &gtmp, symbol);
                vdw1.alpha[i] = alpha;
                vdw1.n[i] = ntmp;
                vdw1.a[i] = atmp;
                vdw1.g[i] = gtmp;
                strcpy(vdw1.da[i],symbol);
                vdw1.nvdw++;               
            }  else if( strcmp(iptemp,"vdw")  == 0 )
            {
                if (vdw1.nvdw > MAXVDWCONST)
                {
                    message_alert("Error - Too many vdw constants","Read Parameter");
                    printf("Maximum number of vdw constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                sscanf( line, "%s %d %f %f %d %d %d", dumm, &i, &radius, &epsilon, &lpde, &ihtype, &ihdonor);

                vdw1.rad[i] = radius;
                vdw1.eps[i] = epsilon;
                vdw1.lpd[i] = lpde;
                vdw1.ihdon[i] = ihdonor;
                vdw1.ihtyp[i] = ihtype;
                vdw1.nvdw++;
             }
        /*   torsion read */
             else if( strcmp(iptemp,"torsion4") == 0 )
             {
                 if (torkn1.ntor4 >= MAXTOR4CONST)
                 {
                      message_alert("Error - Too many torsion 4 constants","Read Parameter");
                      printf("Maximum number of torsion 4 constants exceeded. Others will be ignored\n");
                      fclose(datafile);
                      return;
                 }
                for(i=0; i < 6; i++)
                   v1[i] = se1[i] = 0.0;
                 iz = sscanf(line, "%s  %d %d %d %d %f %d %f %d %f %d %f %d %f %d %f %d", dumm, &kp1,
                         &kp2, &kp3, &kp4 ,
                         &v1[0], &se1[0],&v1[1], &se1[1],&v1[2], &se1[2],
                         &v1[3], &se1[3],&v1[4], &se1[4],&v1[5], &se1[5]);
                 torphase(6, v1, se1);
                 torkn1.tv41[torkn1.ntor4] = v1[0];
                 torkn1.tv42[torkn1.ntor4] = v1[1];
                 torkn1.tv43[torkn1.ntor4] = v1[2];
                 torkn1.phase41[torkn1.ntor4] = se1[0];
                 torkn1.phase42[torkn1.ntor4] = se1[1];
                 torkn1.phase43[torkn1.ntor4] = se1[2];

                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torkn1.kv4[torkn1.ntor4],ptor);
                 torkn1.ntor4++;
             }else if( strcmp(iptemp,"torsion5") == 0 )
             {
                 if (torkn1.ntor5 >= MAXTOR5CONST)
                 {
                      message_alert("Error - Too many torsion 5 constants","Read Parameter");
                      printf("Maximum number of torsion 5 constants exceeded. Others will be ignored\n");
                      fclose(datafile);
                      return;
                 }
                for(i=0; i < 6; i++)
                   v1[i] = se1[i] = 0.0;
                 iz = sscanf(line, "%s  %d %d %d %d %f %d %f %d %f %d %f %d %f %d %f %d", dumm, &kp1,
                         &kp2, &kp3, &kp4 ,
                         &v1[0], &se1[0],&v1[1], &se1[1],&v1[2], &se1[2],
                         &v1[3], &se1[3],&v1[4], &se1[4],&v1[5], &se1[5]);
                 torphase(6, v1, se1);
                 torkn1.tv51[torkn1.ntor5] = v1[0];
                 torkn1.tv52[torkn1.ntor5] = v1[1];
                 torkn1.tv53[torkn1.ntor5] = v1[2];
                 torkn1.phase51[torkn1.ntor5] = se1[0];
                 torkn1.phase52[torkn1.ntor5] = se1[1];
                 torkn1.phase53[torkn1.ntor5] = se1[2];
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torkn1.kv5[torkn1.ntor5],ptor);
                 torkn1.ntor5++;
             }else if( strcmp(iptemp,"torsiondel") == 0 )
             {
                 if (torkn1.ntordel >= MAXTORDEL)
                 {
                      message_alert("Error - Too many delocalized torsion  constants","Read Parameter");
                      printf("Maximum number of delocalized torsion  constants exceeded. Others will be ignored\n");
                      fclose(datafile);
                      return;
                 }
                for(i=0; i < 6; i++)
                   v1[i] = se1[i] = 0.0;
                 iz = sscanf(line, "%s  %d %d %d %d %d %f %d %f %d %f %d %f %d %f %d %f %d", dumm, &torkn1.torindex[torkn1.ntordel],
                         &kp1, &kp2, &kp3, &kp4 ,
                         &v1[0], &se1[0],&v1[1], &se1[1],&v1[2], &se1[2],
                         &v1[3], &se1[3],&v1[4], &se1[4],&v1[5], &se1[5]);
                 torphase(6, v1, se1);
                 torkn1.tvdel1[torkn1.ntordel] = v1[0];
                 torkn1.tvdel2[torkn1.ntordel] = v1[1];
                 torkn1.tvdel3[torkn1.ntordel] = v1[2];
                 torkn1.phasedel1[torkn1.ntordel] = se1[0];
                 torkn1.phasedel2[torkn1.ntordel] = se1[1];
                 torkn1.phasedel3[torkn1.ntordel] = se1[2];
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torkn1.kvdel[torkn1.ntordel],ptor);
                 torkn1.ntordel++;
             }else if( strcmp(iptemp, "torsion") == 0 )
             {
                if (torkn1.ntor >= MAXTORCONST)
                {
                    message_alert("Error - Too many torsion constants","Read Parameter");
                    printf("Maximum number of torsion constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                for(i=0; i < 6; i++)
                   v1[i] = se1[i] = 0.0;
                   
                sscanf( line, "%s %d %d %d %d %f %d %f %d %f %d %f %d %f %d %f %d", dumm,
                         &kp1, &kp2, &kp3, &kp4,
                         &v1[0], &se1[0],&v1[1], &se1[1],&v1[2], &se1[2],
                         &v1[3], &se1[3],&v1[4], &se1[4],&v1[5], &se1[5]);
                torphase(6, v1, se1);
                torkn1.tv1[torkn1.ntor] = v1[0];
                torkn1.phase1[torkn1.ntor] = se1[0];
                torkn1.tv2[torkn1.ntor] = v1[1];
                torkn1.phase2[torkn1.ntor] = se1[1];
                torkn1.tv3[torkn1.ntor] = v1[2];
                torkn1.phase3[torkn1.ntor] = se1[2];
                torkn1.tv4[torkn1.ntor] = v1[3];
                torkn1.phase4[torkn1.ntor] = se1[3];
                torkn1.tv5[torkn1.ntor] = v1[4];
                torkn1.phase5[torkn1.ntor] = se1[4];
                torkn1.tv6[torkn1.ntor] = v1[5];
                torkn1.phase6[torkn1.ntor] = se1[5];
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torkn1.kv[torkn1.ntor],ptor);

//                torkn1.kv[torkn1.ntor] = ((kp1*100 + kp2)*100 + kp3)*100 + kp4;
                torkn1.ntor++;
             }else if (strcmp(iptemp, "pitorsion") == 0)
             {
                 if (torknp.npitor > MAXPITORCONST)
                 {
                 }
                 sscanf( line,"%s %d %d %d %d %f %d %f %d %f %d", dumm, &kp1, &kp2, &kp3, &kp4,
                   &v1[0], &se1[0], &v1[1], &se1[1], &v1[2], &se1[2]);
                   
                 torphase(3, v1,se1);
                 torknp.tv1[torknp.npitor] = v1[0];
                 torknp.ph1[torknp.npitor] = se1[0];
                 torknp.tv2[torknp.npitor] = v1[1];
                 torknp.ph2[torknp.npitor] = se1[1];
                 torknp.tv3[torknp.npitor] = v1[2];
                 torknp.ph3[torknp.npitor] = se1[2];
                 numeral(kp1,pa,3);
                 numeral(kp2,pb,3);
                 numeral(kp3,pc,3);
                 numeral(kp4,pd,3);
                 strcpy(ptor,pa);
                 strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
                 strcpy(torknp.kv[torknp.npitor],ptor);
//                 torknp.kv[torknp.npitor] = ((kp1*100 + kp2)*100 + kp3)*100 + kp4;
                 torknp.npitor++; 
              } else if (strcmp(iptemp, "anglef") == 0)
              {
                  sscanf(line,"%s %d %d %d %f %f %f %f %f %f %f %f",dumm, &kt1,&kt2,&kt3, &angf.fcon[angf.nfang],
                    &angf.fc0[angf.nfang],&angf.fc1[angf.nfang],&angf.fc2[angf.nfang],&angf.fc3[angf.nfang],
                    &angf.fc4[angf.nfang],&angf.fc5[angf.nfang],&angf.fc6[angf.nfang]);
                  numeral(kt1,pa,3);
                  numeral(kt2,pb,3);
                  numeral(kt3,pc,3);
                  strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                  strcpy(angf.kftype[angf.nfang],pang);
                  angf.nfang++;
              } else if( strcmp(iptemp, "angle5") == 0 )
              {
                 if (angk1.nang5 >= MAXANG5CONST)
                 {
                     message_alert("Error - Too many angle 5 constants","Read Parameter");
                     printf("Maximum number of angle 5 constants exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                 }
                 sscanf( line, "%s %d %d %d %f %f %f %f", dumm, 
                         &kt1, &kt2, &kt3, &angk1.con5[angk1.nang5], &angk1.ang5[angk1.nang5][0], 
                         &angk1.ang5[angk1.nang5][1], &angk1.ang5[angk1.nang5][2] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.ktype5[angk1.nang5],pang);
                 angk1.nang5++;
             } else if( strcmp(iptemp, "angle4") == 0 )
             {
                 if (angk1.nang4 >= MAXANG4CONST)
                 {
                     message_alert("Error - Too many angle 4 constants","Read Parameter");
                     printf("Maximum number of angle 4 constants exceeded. Others will be ignored\n");
                     fclose(datafile);
                     return;
                 }
                 sscanf( line, "%s %d %d %d %f %f %f %f", dumm, 
                         &kt1, &kt2, &kt3, &angk1.con4[angk1.nang4], &angk1.ang4[angk1.nang4][0], 
                         &angk1.ang4[angk1.nang4][1], &angk1.ang4[angk1.nang4][2] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.ktype4[angk1.nang4],pang);
                 angk1.nang4++;
             } else if( strcmp(iptemp, "angle3") == 0 )
             {
                if (angk1.nang3 >= MAXANG3CONST )
                {
                    message_alert("Error - Too many angle 3 constants","Read Parameter");
                    printf("Maximum number of angle 3 constants exceeded. Others will be ignored\n");
                    fclose(datafile);
                    return;
                }
                sscanf( line, "%s  %d %d %d %f %f %f %f", dumm, 
                        &kt1, &kt2, &kt3, &angk1.con3[angk1.nang3], &angk1.ang3[angk1.nang3][0], 
                        &angk1.ang3[angk1.nang3][1], &angk1.ang3[angk1.nang3][2] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.ktype3[angk1.nang3],pang);
                angk1.nang3++;
           } else if( strcmp(iptemp, "angle") == 0 )
           {
               if (angk1.nang >= MAXANGCONST)
               {
                  message_alert("Error - Too many angle constants","Read Parameter");
                  printf("Maximum number of angle constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
               }
               sscanf( line, "%s  %d %d %d %f %f %f %f", dumm,
                        &kt1, &kt2, &kt3, &angk1.con[angk1.nang], &angk1.ang[angk1.nang][0], 
                        &angk1.ang[angk1.nang][1], &angk1.ang[angk1.nang][2] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.ktype[angk1.nang],pang);
               angk1.nang++;
           } else if( strcmp(iptemp, "angdel") == 0 )
           {
               if (angk1.ndel >= MAXANGDEL)
               {
                  message_alert("Error - Too many delocalized angle constants","Read Parameter");
                  printf("Maximum number of delocalized angle constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
               }
               sscanf( line, "%s  %d %d %d %f %f %f %f", dumm,
                        &kt1, &kt2, &kt3, &angk1.condel[angk1.ndel], &angk1.angdel[angk1.ndel][0], 
                        &angk1.angdel[angk1.ndel][1], &angk1.angdel[angk1.ndel][2] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.kdel[angk1.ndel],pang);
               angk1.ndel++;
           } else if( strcmp(iptemp, "ang3del") == 0 )
           {
               if (angk1.ndel3 >= MAXANG3DEL)
               {
                  message_alert("Error - Too many delocalized angle constants","Read Parameter");
                  printf("Maximum number of delocalized angle constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
               }
               sscanf( line, "%s  %d %d %d %f %f %f %f", dumm,
                        &kt1, &kt2, &kt3, &angk1.condel3[angk1.ndel3], &angk1.angdel3[angk1.ndel3][0], 
                        &angk1.angdel3[angk1.ndel3][1], &angk1.angdel3[angk1.ndel3][2] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.kdel3[angk1.ndel3],pang);
               angk1.ndel3++;
           } else if( strcmp(iptemp, "ang4del") == 0 )
           {
               if (angk1.ndel4 >= MAXANG4DEL)
               {
                  message_alert("Error - Too many delocalized angle constants","Read Parameter");
                  printf("Maximum number of delocalized angle constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
               }
               sscanf( line, "%s  %d %d %d %f %f %f %f", dumm,
                        &kt1, &kt2, &kt3, &angk1.condel4[angk1.ndel4], &angk1.angdel4[angk1.ndel4][0], 
                        &angk1.angdel4[angk1.ndel4][1], &angk1.angdel4[angk1.ndel4][2] );
                 numeral(kt1,pa,3);
                 numeral(kt2,pb,3);
                 numeral(kt3,pc,3);
                 strcpy(pang,pa); strcat(pang,pb); strcat(pang,pc);
                 strcpy(angk1.kdel4[angk1.ndel4],pang);
               angk1.ndel4++;
          } else if( strcmp(iptemp, "opbend") == 0 )
          {
               if (ooplane_k.nopbend >= MAXOOP)
               {
                  message_alert("Error - Too many out of plane bending constants","Read Parameter");
                  printf("Maximum number of out of plane bend constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
               }
               sscanf( line, "%s %d %d %d %d %f", dumm, &kp1, &kp2, &kp3,&kp4,
                         &ooplane_k.copb[ooplane_k.nopbend] );
               numeral(kp1,pa,3);
               numeral(kp2,pb,3);
               numeral(kp3,pc,3);
               numeral(kp4,pd,3);
               strcpy(ptor,pa); strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
               strcpy(ooplane_k.iopb[ooplane_k.nopbend],ptor);
               ooplane_k.nopbend++;
          } else if (strcmp(iptemp, "strbnd") == 0)
          {
              if (crossterm_k.nstrbnd > MAXSTBN)
              {
                  message_alert("Error - Too many str bend constants","Read Parameter");
                  printf("Maximum number of str bend constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
              }
              kp4 = 0;
              sscanf(line,"%s %d %d %d %f %f %f %d",dumm, &kp1, &kp2, &kp3,
                &crossterm_k.stbncon[crossterm_k.nstrbnd][0],
                &crossterm_k.stbncon[crossterm_k.nstrbnd][1],
                &crossterm_k.stbncon[crossterm_k.nstrbnd][2], &kp4);
               numeral(kp1,pa,3);
               numeral(kp2,pb,3);
               numeral(kp3,pc,3);
               strcpy(ptor,pa); strcat(ptor,pb); strcat(ptor,pc);
               strcpy(crossterm_k.stbn[crossterm_k.nstrbnd],ptor);
               crossterm_k.stbnindex[crossterm_k.nstrbnd] = kp4;
               crossterm_k.nstrbnd++;
          } else if (strcmp(iptemp, "angang") == 0)
          {
              if (crossterm_k.nangang > MAXAA)
              {
                  message_alert("Error - Too many ang ang constants","Read Parameter");
                  printf("Maximum number of ang ang constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
              }
              sscanf(line,"%s %d %f %f %f",dumm, &crossterm_k.ang_ang[crossterm_k.nangang],
               &crossterm_k.aacon[crossterm_k.nangang][0],
               &crossterm_k.aacon[crossterm_k.nangang][1],
               &crossterm_k.aacon[crossterm_k.nangang][2]);
              crossterm_k.nangang++;
          } else if (strcmp(iptemp, "strtors") == 0)
          {
              if (crossterm_k.nstrtor > MAXSTRTOR)
              {
                  message_alert("Error - Too many str tor constants","Read Parameter");
                  printf("Maximum number of str tor constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
              }
              sscanf(line,"%s %d %d %f",dumm, &kt1, &kt2, &crossterm_k.str_torcon[crossterm_k.nstrtor]);
              numeral(kt1,pa,3);
              numeral(kt2,pb,3);
              strcpy(pt,pa); strcat(pt,pb);
              strcpy(crossterm_k.str_tor[crossterm_k.nstrtor],pt);
              crossterm_k.nstrtor++;
          } else if (strcmp(iptemp, "vdwpr") == 0)
          {
              if (vdwpr_k .nvdwpr > MAXBONDCONST)
              {
                  message_alert("Error - Too many vdwpr constants","Read Parameter");
                  printf("Maximum number of vdwpr constants exceeded. Others will be ignored\n");
                  fclose(datafile);
                  return;
              }
              sscanf(line,"%s %d %d %f %f",dumm, &ii, &kk, &vdwpr_k.radius[vdwpr_k .nvdwpr],
              &vdwpr_k.eps[vdwpr_k .nvdwpr]);
              vdwpr_k.ia1[vdwpr_k .nvdwpr] = ii;
              vdwpr_k.ia2[vdwpr_k .nvdwpr] = kk;
              numeral(ii,pa,3);
              numeral(kk,pb,3);
              strcpy(pt,pa); strcat(pt,pb);
              strcpy(vdwpr_k.kv[vdwpr_k .nvdwpr],pt);
              vdwpr_k .nvdwpr++;
          } else if( strcmp(iptemp, "BDHF") == 0 )
          {
               if (ehpara.neheat >= 220)
               {
                   printf("Maximum number of BDHF constants exceeded. Others will be ignored\n");
                   fclose(datafile);
                   return;
               }
               sscanf( line, "%s %d %d %20c %f %f %f", dumm, &ja, &jb,
                         ehpara.cc[ehpara.neheat], &ehpara.ee[ehpara.neheat], &ehpara.estr[ehpara.neheat],
                         &ehpara.et[ehpara.neheat] );
               ehpara.nn[ehpara.neheat] = ja*100 + jb;
               ehpara.neheat++;
          } else if( strcmp(iptemp, "EVHF") == 0 )
          {
               if (ehpara.nevhf >= 155)
               {
                   printf("Maximum number of EVHF constants exceeded. Others will be ignored\n");
                   fclose(datafile);
                   return;
               }
               sscanf( line, "%s %d %d %d %d %d %20c %f %f", dumm, &ehpara.cent[ehpara.nevhf],
                         &kp1,&kp2,&kp3,&kp4, ehpara.mm[ehpara.nevhf],
                         &ehpara.ss[ehpara.nevhf], &ehpara.ssless[ehpara.nevhf]);
               numeral(kp1,pa,3);
               numeral(kp2,pb,3);
               numeral(kp3,pc,3);
               numeral(kp4,pd,3);
               strcpy(ptor,pa); strcat(ptor,pb); strcat(ptor,pc); strcat(ptor,pd);
               strcpy(ehpara.kheat[ehpara.nevhf],ptor);
               ehpara.nevhf++;
          } else if( strcmp(iptemp, "PIHF") == 0 )
          {
               if (epiheat.npihf >= 125)
               {
                   printf("Maximum number of PIHF constants exceeded. Others will be ignored\n");
                   fclose(datafile);
                   return;
               }
               sscanf( line, "%s %d %d %20c %f %f %f %f ", dumm, 
                        &iatm1, &iatm2, epiheat.ccp[epiheat.npihf], &epiheat.eep[epiheat.npihf], 
                        &epiheat.aa[epiheat.npihf], &epiheat.bb[epiheat.npihf], &epiheat.ccc[epiheat.npihf] );
               epiheat.nnp[epiheat.npihf] = iatm1*100 + iatm2;
               epiheat.npihf++;
          }
        }


        fclose(datafile);
/*
        printf(" field : %s\n",field.name);
        printf(" Atom Types: %d\n",atom_k.natomtype);
        printf(" Bonds: %d Bond3: %d Bond4: %d Bond5: %d\n",
         bondk1.nbnd,bondk1.nbnd3, bondk1.nbnd4, bondk1.nbnd5);
        printf(" Angle: %d Angle3: %d Angle4: %d Angle5: %d\n",
         angk1.nang,angk1.nang3, angk1.nang4, angk1.nang5);
        printf(" Torsion: %d  Torsion4: %d Torsion5: %d\n",
         torkn1.ntor, torkn1.ntor4, torkn1.ntor5);
        printf(" Vdw: %d OOP: %d Dipole: %d Charge: %d\n",
         vdw1.nvdw,ooplane_k.nopbend, dipole_k.ndipole, charge_k.ncharge);
        printf(" STBN: %d ANGANG: %d STRTOR: %d VDWPR: %d\n",
         crossterm_k.nstrbnd, crossterm_k.nangang, crossterm_k.nstrtor, vdwpr_k.nvdwpr);
*/
         
        return;

}


void torphase(int icount, float v1[6], int se[6])
{
    int i;
    float amp[6], phase[6];
    int fold[6];

    for (i=0; i < 6; i++)
    {
        fold[i] = 0;
        amp[i] = phase[i] = 0.0;
    }

    for (i=0; i < icount; i++)
    {
        amp[i] = v1[i];
        fold[i] = abs( (int)(2.0*se[i]) - (int) se[i]);
        if (se[i] > 0.0)
          phase[i] = 1.0;
        else
          phase[i] = -1.0;
        v1[i] = 0.0;
        se[i] = 0;
    }

    for (i=0; i < icount; i++)
    {
        if (fold[i] != 0 && fold[i] <= icount)
        {
            v1[fold[i]-1] = amp[i];
            se[fold[i]-1] = phase[i];
        }
    }
}

    