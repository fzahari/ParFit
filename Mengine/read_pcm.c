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
#include "fix.h"
#include "atom_k.h"
#include "optimize.h"
#include "pcmfile.h"
#include "minim_values.h"
#include "tsbond.h"
#include "files.h"
#include "pistuf.h"
#include "units.h"
#include "read_pcm.h"

static int oh2, oh3, oh6;
extern int mm3_mmx[];
extern int mmff_mmx[];

// ========================================================
void read_atomtypes(int field)
{
}
// ========================================================
void get_atomname(int type,char *name)
{
    int i;
    strcpy(name,"");
    for (i=1; i <= MAXATOMTYPE; i++)
    {
        if (type == atom_k.type[i])
        {
            strcpy(name,atom_k.symbol[i]);
            return;
        }
    }
}
// ========================================================
int mm3_mmxtype(int mm3type)
{
    int i;

    if (mm3type < 153)
       i = mm3_mmx[mm3type-1];
    else if (mm3type == 200)
       i = 60;
    else
       i = mm3type;
    return(i);
}
// ==================================================
int mmff_mmxtype(int mmfftype)
{
    int i;
    if (mmfftype < 100)
       i = mmff_mmx[mmfftype-1];
    else
       i = mmfftype;
    return(i);
}
// ==================================================
void pcmfout(int mode)
{
        FILE    *wfile;
        int icoord, imetflag, ispin,lps;
        int i,j, itype;

        /*  Mode = 1 normal packed write operation
         *  Mode = 2 fast no questions asked unpacked write of pcmod.bak
         *  Mode = 3 fast no questions ot file = wfile
         *  mode = 4 start of minim to pcmod.bak but ask questions about pi */

	wfile = NULL;
	itype = 0;
        if( mode == 1 )
        {
                /* assume append is set when we get a filename that exists */
                /* put up dialog box with current structure title, pi flags
                   and added constants question */
           if( files.append )
              wfile = fopen_path(Savebox.path,Savebox.fname,"a");
           else
              wfile = fopen_path(Savebox.path,Savebox.fname,"w");
        }
        else if( mode == 2 || mode == 4 )
        {
              wfile = fopen_path(pcwindir,"pcmod.bak","w");
        }
        else if( mode == 3 )
        {
           if( files.append )
              wfile = fopen_path(Savebox.path,Savebox.fname,"a");
           else
              wfile = fopen_path(Savebox.path,Savebox.fname,"w");
        }
        if (wfile == NULL)
        {
                message_alert("Error opening file in PCMFOUT. Cannot create file","PCM File");
                return;
        }
// remove lone pairs if not writing mmx types
        lps = FALSE;
        if (default_outtype != MMX)
        {
            for (i=1; i <= natom; i++)
            {
                if (atom[i].mmx_type == 20)
                {
                    lps = TRUE;
                    break;
                }
            }
        }
        if (lps == TRUE)
           hdel(1);
        
/* start of file writing  */
        fprintf(wfile,"{PCM %s\n",Struct_Title);
        fprintf(wfile,"NA %d\n",natom);
// optimization information
        fprintf(wfile,"OPT PARA %d CONV %d IE %8.3f FE %8.3f IH %8.3f FH %8.3f\n",optimize_data.param_avail,
           optimize_data.converge,optimize_data.initial_energy,optimize_data.final_energy,
           optimize_data.initial_heat,optimize_data.final_heat);    
/*  flags  */
        fprintf(wfile,"FL ");
        fprintf(wfile,"EINT%d ",minim_values.ndc);
        if( pistuf.iuv )
                fprintf(wfile,"UV%d ",pistuf.iuv);
        if( pistuf.nplane )
                fprintf(wfile,"PIPL%d ",pistuf.nplane);
        if( pistuf.jprint )
                fprintf(wfile,"PIPR%d ",pistuf.jprint);
        if( units.dielec != 1.5 )
                fprintf(wfile,"DIELC%f ",units.dielec);
        fprintf(wfile,"\n");
//  default output filetypes
        fprintf(wfile,"ATOMTYPES %d\n",default_outtype);
/*
        if( minim_values.nconst && mode > 1 )
                fprintf(wfile,"CO FILE %s\n",constf.cname);*/
/*  atoms and bonds lists */
        for( i = 1; i <= natom; i++ )
        {
            if (default_outtype == MMX)
            {
                itype = atom[i].mmx_type;
            } else if (default_outtype == MM3)
            {
                itype = atom[i].mm3_type;
            } else if (default_outtype == MMFF94)
            {
                itype = atom[i].mmff_type;
            }
            fprintf(wfile,"AT %d %d %8.4f %8.4f %8.4f",i,itype,atom[i].x,
                                atom[i].y, atom[i].z);
// bonds
                fprintf(wfile," B");
                for( j = 0; j < MAXIAT; j++ )
                {
                        if( atom[i].iat[j])
                        {
                                fprintf(wfile," %d %d",atom[i].iat[j],atom[i].bo[j]);
                        }
                }
                if( atom[i].flags & (1L << HBOND_MASK) )
                        fprintf(wfile," H ");
                if( atom[i].flags & (1L << PI_MASK) )
                        fprintf(wfile," P ");
                if( atom[i].mmx_type >= 300 )
                {
                        icoord = 0;
                        ispin = -1;
                        imetflag = 0;
                        if( atom[i].flags & (1L << SATMET_MASK) )
                                imetflag += 2;
                        if( atom[i].flags & (1L << GT18e_MASK) )
                                imetflag += 1;
                        if( imetflag == 2 )
                                icoord = 0;
                        if( imetflag == 1 )
                                icoord = 2;
                        if( imetflag == 3 )
                        {
                                icoord = 1;
                                imetflag = 0;
                                if( atom[i].flags & (1L << LOWSPIN_MASK) )
                                        imetflag += 2;
                                if( atom[i].flags & (1L << SQPLAN_MASK) )
                                        imetflag += 1;
                                if( imetflag == 2 )
                                        ispin = 0;
                                if( imetflag == 3 )
                                        ispin = 2;
                                if( imetflag == 1 )
                                        ispin = 1;
                        }
                        if( icoord )
                        {
                                fprintf(wfile," M %d",icoord);
                                if( ispin >= 0 )
                                        fprintf(wfile," %d",ispin);
                        }
                        fprintf(wfile," R %g",atom[i].radius);
                }
                if( atom[i].charge )
                {
                        fprintf(wfile," C %g",atom[i].charge);
                }
                fprintf(wfile,"\n");
        }
        if( fxtor.nfxtor == 1 || fxtor.nfxtor == 2 )
        {
                fprintf(wfile,"DD %d %d %d %d FROM %d TO %d BY %d\n",fxtor.iatom[0][0],  
                 fxtor.iatom[0][1], fxtor.iatom[0][2], fxtor.iatom[0][3], fxtor.start_ang[0],
                 fxtor.final_ang[0], fxtor.step[0] );
                if( fxtor.nfxtor == 2 )
                {
                     fprintf(wfile,"DD %d %d %d %d FROM %d TO %d BY %d\n",fxtor.iatom[1][0],  
                        fxtor.iatom[1][1], fxtor.iatom[1][2], fxtor.iatom[1][3], fxtor.start_ang[1],
                        fxtor.final_ang[1], fxtor.step[1] );
                }
        }
        if (fixdis.nfxstr)
        {
                for( i = 0; i < fixdis.nfxstr; i++ )
                {
                fprintf(wfile,"FIX DIS %d %d R %7.3f K %7.3f\n",fixdis.ifxstr[i][0], 
                                fixdis.ifxstr[i][1], fixdis.fxstrd[i], fixdis.fxstrc[i] );
                }
        }
        for (i=0; i < 15; i++)
        {
                if (ts_bondorder.fbnd[i] > 0.0F)
                {
                        fprintf(wfile,"FBND");
                        for (j=0; j<15; j++)
                        {
                                fprintf(wfile," %5.3f",ts_bondorder.fbnd[j]);
                        }
                        fprintf(wfile,"\n");
                        break;
                }
        }
        
        fprintf(wfile,"}\n");
        fclose(wfile);
} /* end of function */

/* ------------------------------- */

#define NL gettoken(); if (pcmfile.head == 1000) goto L_30;
#define ALPHABETIC      2
#define NUMERIC 1

void pcmfin(int nth,int isubred)
{
        unsigned int   substrismin;
        int             newatom, i, ibondcount, ibonded, ibondorder,
                        icoord, imetcount, ioffset, ispin, iprint,
                        isubmap[MAXSS], j, numbonds, niatom;
        int            newat, newtype, ifile, isMetal, biotype, itype;
        long int       mask;
        char           ianame[3];
        float          ddcount, xtmp, ytmp, ztmp;
        FILE           *pcminfile;

        imetcount = 0;
        isMetal = FALSE;
        ddcount = 0;
        substrismin = False;
        fxtor.nfxtor = 0;
        pcmfile.nocaps = False;
        numbonds = 0;
        ibondcount = 0;
	newatom = 0;
	newtype = 0;

        for (i = 0; i < MAXSS; i++) 
        {
                isubmap[i] = 0;
        }

       ioffset = 0;

        pcminfile = fopen_path(Openbox.path,Openbox.fname,"r");
         
        if (pcminfile == NULL)
        {
                message_alert("Error Opening PCM file","PCM Setup");
                fprintf(pcmoutfile,"PCMFIN error opening file %d %s\n",nth,Openbox.fname);
                return;
        }
        
        if (default_intype != MMX)
            read_atomtypes(default_intype);

        ifile = 0;
        if (nth != 1)
        {
                while ( FetchRecord(pcminfile,pcmfile.string) )
                {
                        if (strncasecmp(pcmfile.string,"{PCM",4) == 0)
                        {
                                ifile++;
                                if (ifile == nth)
                                        break;
                        }
                }
        }else
        {
                FetchRecord(pcminfile,pcmfile.string);
        }

/*  finished scaning to starting structure  now read file */
        if (isubred == 0)
        {
                sscanf(pcmfile.string, "{PCM %60c", Struct_Title);
                clean_string(Struct_Title);
        } else
        {
                if (natom < 1 )
                {
                        sscanf(pcmfile.string, "{PCM %60c", Struct_Title);
                        clean_string(Struct_Title);
                }   
        }
L_30:
        FetchRecord(pcminfile,pcmfile.string);
        *(pcmfile.string + strlen(pcmfile.string) ) = '\0';
        if (feof(pcminfile))
                goto L_170;
        pcmfile.head = 1;
        NL;
L_40:
        if (strcmp(pcmfile.token, "NA") == 0) 
        {
                gettoken();
                if (pcmfile.head == 1000)
                        goto L_30;
                niatom = atoi(pcmfile.token);
                if (niatom + ioffset > MAXATOM) 
                {
                        natom = ioffset;
                        fclose(pcminfile);
                        message_alert(" Maximum number of atoms exceeded in PCMFIN","PCM Setup");
                        fprintf(pcmoutfile," unable to continue reading, more than Max atoms %d\n",MAXATOM);
                        return;
                }
                niatom += ioffset;
                NL;
                goto L_30;
        } else if (strcmp(pcmfile.token,"ATOMTYPES") == 0)
        {
            gettoken();
            if (pcmfile.head == 1000)
               goto L_30;
            default_intype = atoi(pcmfile.token);
            read_atomtypes(default_intype);
            NL;
            goto L_30;
        }  else if (strcmp(pcmfile.token, "NT") == 0)  /* This is for Amino Acids */
        {
                NL;
                i = atoi(pcmfile.token);                /* i is atom number of N terminus */
                mask = 1L << NTERM;                         /* use bit 0 of substr 1 for N terminus*/
                atom[i+ioffset].flags |= mask;

                NL;
                if (strcmp(pcmfile.token, "CT") == 0)
                {
                        NL;
                        i = atoi(pcmfile.token);        /* i is atom number of C terminus  */
                        mask = 1L << CNTERM;                /* use bit 1 of substr 1 for C terminus */
                        atom[i+ioffset].flags |= mask;
                }
        }  else if (strcmp(pcmfile.token, "5PRIME") == 0)  /* This is for Nucleic Acids */
        {
                NL;
                i = atoi(pcmfile.token);                    /* i is atom number of 5prime terminus */
                mask = 1L << P5;                        
                atom[i+ioffset].flags |= mask;
        }else if (strcmp(pcmfile.token, "3PRIME") == 0)  /* This is for Nucleic Acids */
        {
                NL;
                i = atoi(pcmfile.token);        /* i is atom number of C terminus  */
                mask = 1L << P3;                /* use bit 1 of substr 1 for C terminus */
                atom[i+ioffset].flags |= mask;
        }  else if (strcmp(pcmfile.token, "FG") == 0)  /* This is for functional groups */
        {
                NL;
                i = atoi(pcmfile.token);                /* i is atom number of connect point */
                mask = 1L << DUMMY;                         /* use bit 0 of substr 1 for connect*/
                atom[i+ioffset].flags |= mask;

        }  else if (strcmp(pcmfile.token, "OT") == 0)   /*  This is for Sugars  */
        {
                NL;
                i = atoi(pcmfile.token);                /* i is atom number of O terminus */
                mask = 1L << NTERM;                                 /* use bit 0 of substr 1 for O terminus*/
                atom[i+ioffset].flags |= mask;

                NL;
                if (strcmp(pcmfile.token, "CT") == 0)
                {
                        NL;
                        i = atoi(pcmfile.token);        /* i is atom number of H terminus  */
                        mask = 1L << CNTERM;                         /* use bit 1 of substr 1 for H terminus */
                        atom[i+ioffset].flags |= mask;
                }

                NL;
                if (strcmp(pcmfile.token, "OTWO") == 0)
                {
                        NL;
                        i = atoi(pcmfile.token);        /* i is atom number of H terminus  */
                        oh2 = i+ioffset;
                }

                NL;
                if (strcmp(pcmfile.token, "OTHR") == 0)
                {
                        NL;
                        i = atoi(pcmfile.token);        /* i is atom number of H terminus  */
                        oh3 = i+ioffset;
                }

                NL;
                if (strcmp(pcmfile.token, "OSIX") == 0)
                {
                        NL;
                        i = atoi(pcmfile.token);        /* i is atom number of H terminus  */
                        oh6 = i+ioffset;
                }
        } else if (strcmp(pcmfile.token,"START") == 0) /* flag for polybuild */
        {
                NL;
                i = atoi(pcmfile.token);
                mask = 1L << DUMMY;
                atom[i+ioffset].flags |= mask;
                NL;   // END
                NL;
                i = atoi(pcmfile.token);
                mask = 1L << DUMMY;
                atom[i+ioffset].flags |= mask;
        } else if  (strcmp(pcmfile.token,"NMOL") == 0)
        {
            NL;
            NL;
            NL;
            NL;
            NL;
            NL;
            NL;
            goto L_30;
        } else if (strcmp(pcmfile.token, "OPT") == 0)
        {
            NL;  // para
            NL;
            optimize_data.param_avail = atoi(pcmfile.token);                   
            NL; // conv
            NL;
            optimize_data.converge = atoi(pcmfile.token);                   
            NL;  // ie
            NL;
            optimize_data.initial_energy = atof(pcmfile.token);                   
            NL;  // fe
            NL;
            optimize_data.final_energy = atof(pcmfile.token);                   
            NL;   // iH
            NL;
            optimize_data.initial_heat = atof(pcmfile.token);
            NL; // fH             
            NL;
            optimize_data.final_heat = atof(pcmfile.token);                   
            goto L_30;
        }  else if (strcmp(pcmfile.token, "SS") == 0)
        {
                goto L_30;
        }  else if (strcmp(pcmfile.token, "AT") == 0) 
        {
                NL;
                newat = atoi(pcmfile.token);   // atom serial number
                if (newat + ioffset > MAXATOM) 
                {
                        message_alert(" Maximum number of atoms exceeded in PCMFIN at atom read","PCM Setup");
                        fprintf(pcmoutfile,"atom read %d %d %d\n",newat,ioffset,natom);
                        natom = ioffset + newat - 1;
                        fclose(pcminfile);
                        exit(0);
                }
                newat += ioffset;
                NL;
                if (isdigit(*pcmfile.token)) 
                {
                        newtype = atoi(pcmfile.token);
                }       else   // character string for metal atom
                {
                        strcpy(ianame,pcmfile.token);
                        isMetal = TRUE;
                        i = strlen(ianame);
                        if ( i == 2) 
                            ianame[1] = tolower(ianame[1]);
                }
                NL;
                xtmp = atof(pcmfile.token);
                NL;
                ytmp = (atof(pcmfile.token));
                NL;
                ztmp = atof(pcmfile.token);
                NL;
                if (isMetal == FALSE)
                {
                    itype = newtype;
                    get_atomname(newtype,ianame);
                    if (default_intype == MM3)
                       itype = mm3_mmxtype(newtype);
                    else if (default_intype == MMFF94)
                       itype = mmff_mmxtype(newtype);
                       
                    newatom = make_atom(newtype, xtmp, ytmp, ztmp,ianame);
                    if (default_intype == MM3)
                        set_atomtype(newatom, itype,newtype,0,0,0);
                    else if (default_intype == MMFF94)
                        set_atomtype(newatom, itype,0,newtype,0,0);
                }else if (isMetal == TRUE)
                {
                    newatom = make_atom(0, xtmp, ytmp, ztmp,ianame);
                    isMetal = FALSE;
                }

L_90:
                if (strcmp(pcmfile.token, "B") == 0) 
                {
                        /* Reading bonds now */
                        NL;
                        j = 0;
        L_100:
                        j += 1;
                        if (pcmfile.state == NUMERIC) 
                        {
                                ibonded = atoi(pcmfile.token);
                                ibonded += ioffset;
                        } else 
                        {
                                goto L_90;
                        }
                        NL;
                        if (pcmfile.state == NUMERIC) 
                        {
                                ibondorder = atoi(pcmfile.token);
                        } else 
                        {
                                fprintf(stdout, "missing a bondorder, assuming 1\n");
                                goto L_40;
                        }
                        if (newatom < ibonded)
                                make_bond(newatom, ibonded, ibondorder);
                        if (ibondorder == 9)
                        {
                                atom[newatom].flags |= (1L << METCOORD_MASK);
                                atom[ibonded].flags |= (1L << METCOORD_MASK);
                        }       
                        NL;
                        goto L_100;
                }else if ( strcmp(pcmfile.token, "P")== 0)
                {
                        atom[newatom].flags |= (1L << PI_MASK);
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token, "D") == 0)   // biomolecule marking
                {
                        NL;
                        biotype = atoi(pcmfile.token);
                        atom[newatom].biotype =  biotype;                       
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token,"M")==0)
                {
                        NL;
                        icoord = atoi(pcmfile.token);
                        if (icoord == 0)
                        {
                                atom[newatom].flags |= (1L << SATMET_MASK);
                        } else if (icoord == 2)
                        {
                                atom[newatom].flags |= (1L << GT18e_MASK);
                        } else{
                                atom[newatom].flags |= (1L << GT18e_MASK);
                                atom[newatom].flags |= (1L << SATMET_MASK);
                        }
                        NL;
                        if (pcmfile.state != 1)
                                goto L_90;
                        ispin = atoi(pcmfile.token);
                        if (ispin == 0)
                        {
                                atom[newatom].flags |= (1L << LOWSPIN_MASK);
                        } else if (ispin == 1) 
                        {
                                atom[newatom].flags |= (1L << SQPLAN_MASK);
                        } else 
                        {
                                atom[newatom].flags |= (1L << LOWSPIN_MASK);
                                atom[newatom].flags |= (1L << SQPLAN_MASK);
                        }
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token,"S")== 0)
                {
                        NL;
                        if (pcmfile.state == NUMERIC)
                        {
                                j = atoi(pcmfile.token);
                                mask = (1L << j);
                                NL;
                                goto L_90;
                        } else
                        {
                                goto L_90;
                        }
                }else if (strcmp(pcmfile.token,"H")== 0)
                {
                        atom[newatom].flags |= (1L << HBOND_MASK);
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token,"C")==0)
                {
                        NL;
                        atom[newatom].charge = atof(pcmfile.token);
                        NL;
                        goto L_90;
                } else if (strcmp(pcmfile.token,"ML") == 0)
                {
                    NL;
                    j = atoi(pcmfile.token);
                    atom[newatom].molecule = j;
                    goto L_90;
                }else if (strcmp(pcmfile.token,"R")== 0)
                {
                        NL;
                        atom[newatom].radius = atof(pcmfile.token);
                        NL;
                        goto L_90;
                }               
        } else if (strncasecmp(pcmfile.token, "}",1) == 0) 
        {
                goto L_170;             /* end of structure */
        }  else if (strcmp(pcmfile.token, "FL") == 0) 
        {
                NL;
L_120:          if (strcmp(pcmfile.token, "PRINT") == 0) 
                {
                        iprint = 1;
                        NL;
                        if (pcmfile.state == ALPHABETIC)
                                goto L_120;
                        iprint = atoi(pcmfile.token);
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "DIELC") == 0) 
                {
                        NL;
                        if (pcmfile.state != 1)
                                goto L_120;
//                        units.dielec = atof(pcmfile.token);
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "UV") == 0) 
                {
                        NL;             
                        if (pcmfile.state != 1)
                                goto L_120;
                        pistuf.iuv = atoi(pcmfile.token);
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "PIPR") == 0) 
                {
                        NL;             
                        if (pcmfile.state != 1)
                                goto L_120;
                        pistuf.jprint = atoi(pcmfile.token);
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "PIPL") == 0) 
                {
                        NL;             
                        if (pcmfile.state != 1)
                                goto L_120;
                        pistuf.nplane = atoi(pcmfile.token);
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "EINT") == 0) 
                {
                        NL;             
                        if (pcmfile.state != 1)
                                goto L_120;
                        minim_values.ndc = atoi(pcmfile.token);
                        NL;
                        goto L_120;
                }
        }  else if (strcmp(pcmfile.token, "CO") == 0) 
        {
        } else if (strcmp(pcmfile.token, "FBND") == 0)
        {
                NL;
                ts_bondorder.fbnd[0] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[1] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[2] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[3] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[4] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[5] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[6] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[7] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[8] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[9] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[10] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[11] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[12] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[13] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[14] = atof(pcmfile.token);
        }  else if (strcmp(pcmfile.token, "FIX") == 0) 
        {
                NL;
                if (strcmp(pcmfile.token, "DIS") == 0) 
                {
                        NL;
                        fixdis.ifxstr[fixdis.nfxstr][0] = atoi(pcmfile.token);
                        NL;
                        fixdis.ifxstr[fixdis.nfxstr][1] = atoi(pcmfile.token);
                        NL;
L_1223:                 if (strcmp(pcmfile.token, "R") == 0) 
                        {
                                NL;
                                fixdis.fxstrd[fixdis.nfxstr] = atof(pcmfile.token);
                                NL;
                                goto L_1223;
                        } else if (strcmp(pcmfile.token, "K") == 0) 
                        {
                                NL;
                                fixdis.fxstrc[fixdis.nfxstr] = atof(pcmfile.token);
                                fixdis.nfxstr +=1;
                                NL;
                                goto L_1223;
                        }
                }
        }  else if (strcmp(pcmfile.token, "DD") == 0) 
        {
                if (isubred == 1)
                        goto L_30;

                NL;
                
                if (pcmfile.state == NUMERIC)
                {
                    fxtor.iatom[fxtor.nfxtor][0] = atoi(pcmfile.token);
                    NL;
                    fxtor.iatom[fxtor.nfxtor][1] = atoi(pcmfile.token);
                    NL;
                    fxtor.iatom[fxtor.nfxtor][2] = atoi(pcmfile.token);
                    NL;
                    fxtor.iatom[fxtor.nfxtor][3] = atoi(pcmfile.token);
                    NL;
L_150:              if (strcmp(pcmfile.token, "FROM") == 0)
                    {
                       NL;
                       fxtor.start_ang[fxtor.nfxtor] = atoi(pcmfile.token);
                       fxtor.curr_ang[fxtor.nfxtor] = fxtor.start_ang[fxtor.nfxtor];
                       NL;
                       goto L_150;
                    } else if (strcmp(pcmfile.token, "TO") == 0)
                    {
                       NL;
                       fxtor.final_ang[fxtor.nfxtor] = atoi(pcmfile.token);
                       NL;
                       goto L_150;
                    } else if (strcmp(pcmfile.token, "BY") == 0)
                    {
                       NL;
                       fxtor.step[fxtor.nfxtor] = atoi(pcmfile.token);
                       fxtor.nfxtor++;
                       NL;
                       goto L_150;                       
                    }
                }
        } 
        goto L_30;
L_170:
        fclose(pcminfile);
        for (i=1; i <= natom; i++)
        {
                if (atom[i].mmx_type == 5)
                {
                        flags.noh = True;
                        break;
                }
        }
        if (default_intype != MMX)
           read_atomtypes(MMX);
        return;
} 

void gettoken()
{
        int   foundtoken, i;

        pcmfile.state = 0;
        foundtoken = -1;
        strcpy(pcmfile.token, "");
        for (i = pcmfile.head - 1; i <= strlen(pcmfile.string); i++) {
                if ((*(pcmfile.string + i) != ':') &&
                    (*(pcmfile.string + i) != ' ') &&
                    (*(pcmfile.string + i) != '\0') &&
                    (*(pcmfile.string + i) != ';') &&
                    (*(pcmfile.string + i) != ',') &&
                    (*(pcmfile.string + i) != '\t')) {
                        if (foundtoken == -1) {
                                foundtoken = i;
                                if (isdigit(*(pcmfile.string + i))) {
                                        /* beginning to read a number */
                                        pcmfile.state = 1;
                                } else if ((*(pcmfile.string + i) == '.') ||
                                         (*(pcmfile.string + i) == '-')) {
                                        /* period or dash can't yet tell if
                                         * word or number */
                                        pcmfile.state = 3;
                                } else {
                                        /* must be a word */
                                        pcmfile.state = 2;
                                }
                        } else {
                                if (
                                    (
                                     (pcmfile.state == 2) &&
                                     (isdigit(*(pcmfile.string + i)))
                                     )
                                    ||
                                    (
                                     (pcmfile.state == 1) &&
                                     (isalpha(*(pcmfile.string + i)))
                                     )
                                    ) {
                                        foundtoken = i;
                                        goto L_30;
                                } else if (pcmfile.state == 3) {
                                        if (isdigit(*(pcmfile.string + i)) ||
                                            *(pcmfile.string + i) == '.') {
                                                pcmfile.state = 1;
                                        } else {
                                                pcmfile.state = 2;
                                        }
                                }
                        }
                        strncat(pcmfile.token, pcmfile.string + i, 1);
                } else if (foundtoken != -1) {
                        foundtoken = i;
                        goto L_30;
                }
        }
        if (foundtoken == -1) {
                pcmfile.head = 1000;
                return;
        }
L_30:
        pcmfile.head = foundtoken + 1;
        if (pcmfile.state == 2 && !pcmfile.nocaps)
        {
            for (i = 0; i < strlen(pcmfile.token); i++)
            {
                if (!isupper(*(pcmfile.token + i)))
                     *(pcmfile.token + i) = toupper(*(pcmfile.token + i));
            }
        }
        return;
} 
/* ============================================================ */
 FILE           *pcminfile;
/* ============================================================ */
void rdfile(int mode, int isub)
{
     if (mode == -1)
        fclose(pcminfile);
     else if (mode == 0)
     {
        pcminfile = fopen_path(Openbox.path,Openbox.fname, "r");
         
        if (pcminfile == NULL)
        {
                message_alert("Error Opening PCM file","Error");
                return;
        }
     } else
        rdpcm(mode,isub);
}
/* ============================================================ */
void rdpcm(int nth,int isubred)
{
        unsigned int   substrismin;
        int            i, ibondcount, ibonded, ibondorder,
                        imetcount, ioffset,
                        isubmap[MAXSS], j, numbonds, niatom;
        int            newat, newtype, isMetal;
        long int       mask;
        char           ianame[3];
        float          ddcount, xtmp, ytmp, ztmp;

        imetcount = 0;
        isMetal = FALSE;
        ddcount = 0;
        substrismin = False;
        pcmfile.nocaps = False;
        numbonds = 0;
        ibondcount = 0;
        for (i = 0; i < MAXSS; i++) 
        {
                isubmap[i] = 0;
        }

                ioffset = 0;


        FetchRecord(pcminfile,pcmfile.string);

/*  finished scaning to starting structure  now read file */
        if (isubred == 0)
        {
                sscanf(pcmfile.string, "{PCM %60c", Struct_Title);
                for(i=0; i < 60; i++)
                {
                    if (Struct_Title[i] == '\n')
                    {
                        Struct_Title[i] = '\0';
                        break;
                    }
                }
        } else
        {
                if (natom < 1 )
                {
                        sscanf(pcmfile.string, "{PCM %60c", Struct_Title);
                        for(i=0; i < 60; i++)
                        {
                           if (Struct_Title[i] == '\n')
                           {
                             Struct_Title[i] = '\0';
                             break;
                           }
                        }
                }     
        }
L_30:
        FetchRecord(pcminfile,pcmfile.string);
        *(pcmfile.string + strlen(pcmfile.string) ) = '\0';
        if (feof(pcminfile))
                goto L_170;
        pcmfile.head = 1;
        NL;
L_40:
        if (strcmp(pcmfile.token, "NA") == 0) 
        {
                gettoken();
                if (pcmfile.head == 1000)
                        goto L_30;
                niatom = atoi(pcmfile.token);
                if (niatom + ioffset > MAXATOM) 
                {
                        natom = ioffset;
                        fclose(pcminfile);
                        message_alert(" Maximum number of atoms exceeded in PCMFIN","Error");
                        fprintf(pcmoutfile," unable to continue reading, more than Max atoms %d\n",MAXATOM);
                        return;
                }
                niatom += ioffset;
                NL;
                goto L_30;
        } else if (strcmp(pcmfile.token,"ATOMTYPES") == 0)
        {
            gettoken();
            if (pcmfile.head == 1000)
               goto L_30;
            default_intype = atoi(pcmfile.token);
//            read_atomtypes(default_intype);
            NL;
            goto L_30;
        }  else if (strcmp(pcmfile.token, "NT") == 0)  /* This is for Amino Acids */
        {
                NL;
                i = atoi(pcmfile.token);                /* i is atom number of N terminus */
                mask = 1L << NTERM;                                 /* use bit 0 of substr 1 for N terminus*/
                atom[i+ioffset].flags |= mask;

                NL;
                if (strcmp(pcmfile.token, "CT") == 0)
                {
                        NL;
                        i = atoi(pcmfile.token);        /* i is atom number of C terminus  */
                        mask = 1L << CNTERM;                         /* use bit 1 of substr 1 for C terminus */
                        atom[i+ioffset].flags |= mask;
                }
        }  else if (strcmp(pcmfile.token, "OT") == 0)   /*  This is for Sugars  */
        {
                NL;
                i = atoi(pcmfile.token);                /* i is atom number of O terminus */
                mask = 1L << OTERM;                                 /* use bit 0 of substr 1 for O terminus*/
                atom[i+ioffset].flags |= mask;

                NL;
                if (strcmp(pcmfile.token, "CT") == 0)
                {
                        NL;
                        i = atoi(pcmfile.token);        /* i is atom number of H terminus  */
                        mask = 1L << COTERM;                         /* use bit 1 of substr 1 for H terminus */
                        atom[i+ioffset].flags |= mask;
                }
        } else if (strcmp(pcmfile.token,"START") == 0) /* flag for polybuild */
        {
                NL;
                i = atoi(pcmfile.token);
                mask = 1L << DUMMY;
                atom[i+ioffset].flags |= mask;
                NL;   // END
                NL;
                i = atoi(pcmfile.token);
                mask = 1L << DUMMY;
                atom[i+ioffset].flags |= mask;
                
        } else if (strcmp(pcmfile.token, "OPT") == 0)
        {
            NL;  // para
            NL;
            optimize_data.param_avail = atoi(pcmfile.token);                   
            NL; // conv
            NL;
            optimize_data.converge = atoi(pcmfile.token);                   
            NL;  // ie
            NL;
            optimize_data.initial_energy = atof(pcmfile.token);                   
            NL;  // fe
            NL;
            optimize_data.final_energy = atof(pcmfile.token);                   
            NL;   // iH
            NL;
            optimize_data.initial_heat = atof(pcmfile.token);
            NL; // fH             
            NL;
            optimize_data.final_heat = atof(pcmfile.token);                   
            goto L_30;
        }  else if (strcmp(pcmfile.token, "SS") == 0)
        {
                goto L_30;
        }  else if (strcmp(pcmfile.token, "AT") == 0) 
        {
                NL;
                newat = atoi(pcmfile.token);   // atom serial number
                if (newat + ioffset > MAXATOM) 
                {
                        message_alert(" Maximum number of atoms exceeded in PCMFIN at atom read","Error");
                        fprintf(pcmoutfile,"atom read %d %d %d\n",newat,ioffset,natom);
                        natom = ioffset + newat - 1;
                        fclose(pcminfile);
                        exit(0);
                }
                newat += ioffset;
                NL;
                if (isdigit(*pcmfile.token)) 
                {
                        newtype = atoi(pcmfile.token);
                }       else   // character string for metal atom
                {
                        strcpy(ianame,pcmfile.token);
                        isMetal = TRUE;
                }
                NL;
                xtmp = atof(pcmfile.token);
                NL;
                ytmp = (atof(pcmfile.token));
                NL;
                ztmp = atof(pcmfile.token);
                NL;
                atom[newat].x = xtmp;
                atom[newat].y = ytmp;
                atom[newat].z = ztmp;

L_90:
                if (strcmp(pcmfile.token, "B") == 0) 
                {
                        /* Reading bonds now */
                        NL;
                        j = 0;
        L_100:
                        j += 1;
                        if (pcmfile.state == NUMERIC) 
                        {
                                ibonded = atoi(pcmfile.token);
                                ibonded += ioffset;
                        } else 
                        {
                                goto L_90;
                        }
                        NL;
                        if (pcmfile.state == NUMERIC) 
                        {
                                ibondorder = atoi(pcmfile.token);
                        } else 
                        {
                                fprintf(stdout, "missing a bondorder, assuming 1\n");
                                goto L_40;
                        }
                        NL;
                        goto L_100;
                }else if ( strcmp(pcmfile.token, "P")== 0)
                {
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token, "D") == 0)   // biomolecule marking
                {
                        NL;
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token,"M")==0)
                {
                        NL;
                        NL;
                        if (pcmfile.state != 1)
                                goto L_90;
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token,"S")== 0)
                {
                        NL;
                        if (pcmfile.state == NUMERIC)
                        {
                                NL;
                                goto L_90;
                        } else
                        {
                                goto L_90;
                        }
                }else if (strcmp(pcmfile.token,"H")== 0)
                {
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token,"C")==0)
                {
                        NL;
                        NL;
                        goto L_90;
                }else if (strcmp(pcmfile.token,"R")== 0)
                {
                        NL;
                        NL;
                        goto L_90;
                }               
        } else if (strncmp(pcmfile.token, "}",1) == 0) 
        {
                goto L_170;             /* end of structure */
        }  else if (strcmp(pcmfile.token, "FL") == 0) 
        {
                NL;
L_120:  if (strcmp(pcmfile.token, "PRINT") == 0) 
                {
                        NL;
                        if (pcmfile.state == ALPHABETIC)
                                goto L_120;
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "DIELC") == 0) 
                {
                        NL;
                        if (pcmfile.state != 1)
                                goto L_120;
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "UV") == 0) 
                {
                        NL;             
                        if (pcmfile.state != 1)
                                goto L_120;
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "PIPR") == 0) 
                {
                        NL;             
                        if (pcmfile.state != 1)
                                goto L_120;
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "PIPL") == 0) 
                {
                        NL;             
                        if (pcmfile.state != 1)
                                goto L_120;
                        NL;
                        goto L_120;
                } else if (strcmp(pcmfile.token, "EINT") == 0) 
                {
                        NL;             
                        if (pcmfile.state != 1)
                                goto L_120;
                        NL;
                        goto L_120;
                }
        }  else if (strcmp(pcmfile.token, "CO") == 0) 
        {
        } else if (strcmp(pcmfile.token, "FBND") == 0)
        {
                NL;
                ts_bondorder.fbnd[0] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[1] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[2] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[3] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[4] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[5] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[6] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[7] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[8] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[9] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[10] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[11] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[12] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[13] = atof(pcmfile.token);
                NL;
                ts_bondorder.fbnd[14] = atof(pcmfile.token);
            goto L_30; // get another line
        }  else if (strcmp(pcmfile.token, "FIX") == 0) 
        {
                NL;
                if (strcmp(pcmfile.token, "DIS") == 0) 
                {
                        NL;
                        NL;
                        NL;
L_1223:         if (strcmp(pcmfile.token, "R") == 0) 
                        {
                                NL;
                                NL;
                                goto L_1223;
                        } else if (strcmp(pcmfile.token, "K") == 0) 
                        {
                                NL;
                                NL;
                                goto L_1223;
                        }
                } 
        }  else if (strcmp(pcmfile.token, "DD") == 0) 
        {
                if (isubred == 1)
                        goto L_30;

                NL;
                
                if (pcmfile.state == NUMERIC)
                {
                    NL;
                    NL;
                    NL;
                    NL;
L_150:              if (strcmp(pcmfile.token, "FROM") == 0)
                    {
                       NL;
                       NL;
                       goto L_150;
                    } else if (strcmp(pcmfile.token, "TO") == 0)
                    {
                       NL;
                       NL;
                       goto L_150;
                    } else if (strcmp(pcmfile.token, "BY") == 0)
                    {
                       NL;
                       NL;
                       goto L_150;                       
                    }
                }
        } 
        goto L_30;
L_170:
        return;
} 
