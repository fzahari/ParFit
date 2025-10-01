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
#include "utility.h"
#include "files.h"
#include "pcmfile.h"
#include "read_mm2.h"

#define NUMERIC 1

// ==============================================
double ReadMMFloat(char *inputline, int offset, int end)
{
    int i, count;
    char tmpname[15];
    float value;

    count = 0;
    for (i=offset; i < offset+end; i++)
    {
        tmpname[count] = inputline[i];
        count++;
    }
    sscanf(tmpname,"%f",&value);
    return (value);
}
// ==============================================
int ReadMMInt(char *inputline, int offset,int end)
{
    int i, value,count;
    char tmpname[15];

    count = 0;
    for (i=offset; i < offset+end; i++)
    {
        tmpname[count] = inputline[i];
        count++;
    }
    sscanf(tmpname,"%d",&value);
    return (value);
}

void write_mm3(void)
{
    int *logary;
    FILE *writefile;
    int i,j, method;
    int nsymm, nx, nrot, label, ndc, ncalc, hform, mvdw, ndrive;
    int iprint, nstr, init, nconst;
    int iset, jprint, ibeta, limit, iter, ndewar, lspi, nbpol, 
        lpchg, ipscale, icont;
    float tmax, del;
    int **attach, **iconn;
    int nattch, nconn, LPON;

/* check for metals and lone pairs */
    LPON = FALSE;
    for (i=1; i < natom; i++)
    {
        if (atom[i].mmx_type == 20)
           LPON = TRUE;
    }
    if (LPON == TRUE)
       hdel(1);
       
     attach = imatrix(0, natom, 0, 2);
     iconn = imatrix(0, 16, 0, (natom/2)+1);
     logary = ivector(0, natom+1);
      
/* setup logary lists  */
    method = 0; 
    for ( i=1; i < natom; i++)
    {
      logary[i] = FALSE;
      if( atom[i].flags & (1L << PI_MASK))
      {
          logary[i] = TRUE;
          method = 1;
      }
    }

    if( files.append )
      writefile = fopen_path(Savebox.path,Savebox.fname,"a");
    else
      writefile = fopen_path(Savebox.path,Savebox.fname,"w");

   if (writefile == NULL)
   {
       message_alert("Error open file in MOD2MM3. Null pointer","MM3 Setup");
       exit(0);
   }
/* generate mm2 attachments */
       generate_mm2lists(&nattch, &nconn, attach, iconn);

/* first line */
  iprint = 3;
  nstr = 0;
  init = 0;
  nconst = 0;
  tmax = 999.0;
  fprintf(writefile,"%-60s%d%4d %d  %d %d  %d%-5.0f\n",
          Struct_Title, method, natom, iprint, nstr, init, nconst, tmax);

/* second line if Method = 1 is line of logicals 
 * else is nconn line  */
   if (method == 1)
   {
      for( i = 1; i <= 60; i++ )
      {
         if (i > natom)
          fprintf(writefile,"F" );          
         else if (logary[i] == TRUE)
          fprintf(writefile,"T" );
         else if (logary[i] == FALSE)
          fprintf(writefile,"F" );
      }
      iset = 0;
      jprint = 0;
      ibeta = 0;
      limit = 5;
      iter = 8;
      ndewar = 0;
      lspi = 0;
      nbpol = 0;
      lpchg = 0;
      ipscale = 0;
      icont = 0;
      if (natom > 60)
      {
        for (i=61; i < natom; i++)
        {
                if (logary[i] == TRUE)
                   icont = 1;
        }
      }
      fprintf( writefile,"%2d%2d%2d%2d%2d%2d%2d%2d%1d%1d%1d%1d\n",
           1,iset, jprint, ibeta,
           limit, iter, ndewar, lspi, nbpol, lpchg, ipscale, icont);
      if( icont == 1)
      {
         for( i = 61; i <= 140; i++ )
         {
            if (i > natom)
              fprintf(writefile,"F" );             
            else if (logary[i] == TRUE)
              fprintf(writefile,"T");
            else if (logary[i] == FALSE)
              fprintf(writefile,"F" );
         }
         fprintf(writefile, "\n" );
      }
   }      
/*  nconn line  */
  nsymm = 0;
  del = 0.00008;
  nx = 0;
  nrot = 0;
  label = 0;
  ndc = 0;
  ncalc = 0;
  hform = 0;
  mvdw = 1;
  ndrive = 0;
  fprintf(writefile,"%1d%4d%5s%4.5f%8s%5d%5d%5d%5d%5d%5d%5d%5d%10d%5d\n",
          0, nconn, "", del, "", nattch, nsymm, nx, nrot, label, ndc, ncalc,
          hform, mvdw, ndrive);

/* iconn and nattch lines */
   for( i = 0; i < nconn; i++ )
   {
      for( j = 0; j < 16; j++ )
      {
         fprintf(writefile, "%5d", iconn[j][i] );
      }
      fprintf(writefile, "\n" );
   }
   if( nattch > 0 )
   {
      i = 0;
      for( j = 0; j < nattch; j++ )
      {
         i++;
         fprintf(writefile, "%5d%5d", attach[j][0], attach[j][1] );
         if ( i ==  8 )
         {
            fprintf(writefile, "\n" );
            i =0;
         }
      }
      if ( (nattch)%8 != 0)
        fprintf(writefile, "\n" );
   }
/* write atom information  */
/* could use routine here to reset atoms types to ensure that we have correct
   MM2 types   */

   for( i = 1; i <= natom; i++ )
   {
       fprintf(writefile,"%10.5f%10.5f%10.5f%5d\n", 
          atom[i].x, atom[i].y, atom[i].z, atom[i].mm3_type); 
   }
   fclose(writefile);

   free_imatrix( attach, 0, natom, 0, 2);
   free_imatrix(iconn ,0, 16, 0, (natom/2)+1);
   free_ivector(logary, 0, natom+1);

   if (LPON == TRUE)
      hadd();

}

void mm32mod(int isubred)
{
    char line[91], dummy1[100],atomname[3],typename[4];
    FILE *infile;
    int i,j, method, niatom, ijunk,itype,newatom,nnatch, iz, ic[16];
    int nsymm, nx, nrot, label, ndc, ncalc, mvdw, ndrive;
    int iprint, nstr, init, nconst, ibotptr, newtype;
    int iset, jprint, ibeta, limit, iter, ndewar, lspi, nbpol, 
        lpchg, ipscale, icont;
    float tmax, del, xtmp,ytmp,ztmp;
    int **iconn;
    int **attach;
    int nattch, nconn;


     if( !isubred )
     {
        flags.noh = True;
        ibotptr = 0;
     } else
        ibotptr = natom;

    infile = fopen_path(Openbox.path,Openbox.fname,"r");
   if (infile == NULL)
   {
       message_alert("Error open file in MM3 read. Null pointer","MM3 Setup");
       return;
   }

   FetchRecord(infile, line);
   sscanf(line,"%60c %d %4d %d  %d %d  %d %f",
          Struct_Title, &method, &niatom, &iprint, &nstr, &init, &nconst, &tmax);

        if (niatom < 0 || niatom > MAXATOM)
        {
            message_alert("Error reading number of atoms. File will not be read. Please check file type","File Read Error");
            fclose(infile);
            return;
        }

/* second line if Method = 1 is line of logicals else is nconn line  */
   if (method == 1)
   {
      FetchRecord(infile, line);
      for( i = 1; i <= 60; i++ )
      {
         if ( line[i-1] == 'T' || line[i-1] == 't' )
                atom[i+ibotptr].flags |= (1L << PI_MASK);
      }
      sscanf(line,"%60c %2d%2d%2d%2d%2d%2d%2d%2d%1d%1d%1d%1d",
           dummy1,&ijunk,&iset, &jprint, &ibeta,
           &limit, &iter, &ndewar, &lspi, &nbpol, &lpchg, &ipscale, &icont);
  }      
  FetchRecord(infile, line);
  sscanf(line,"%d %d  %f %d %d %d %d %d %d %d %d %d \n",
          &ijunk, &nconn,  &del, &nattch, &nsymm, &nx, &nrot, &label,
          &ndc, &ncalc, &mvdw, &ndrive);

     attach = imatrix(0, niatom, 0, 2);
     iconn = imatrix(0, 16, 0, nconn+1);

     for( i = 0; i < 2; i++ )
     {
         for( j = 0; j < niatom; j++ )
              attach[j][i] = 0;
     }
     for( i = 0; i < nconn+1; i++ )
     {
        for( j = 0; j < 16; j++ )
           iconn[j][i] = 0;
     }

/* iconn and nattch lines */
   for( i = 0; i < nconn; i++ )
   {
         FetchRecord(infile, line);
         iz = sscanf(line,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &ic[0], &ic[1],&ic[2],&ic[3],
            &ic[4],&ic[5],&ic[6],&ic[7], &ic[8], &ic[9],&ic[10],&ic[11], &ic[12], &ic[13],&ic[14],&ic[15] );
         for (j=0; j < iz; j++)
            iconn[j][i] = ic[j];
   }
   if( nattch > 0 )
   {
      ijunk = nattch/8;
      if (nattch%8 > 0)
        ijunk++;
      for( j = 0; j < ijunk; j++ )
      {
         if (j == ijunk-1 )
         {
             nnatch = nattch%8;
             if (nnatch == 0)
                 nnatch = 8;
         }else
         {
             nnatch = 8;
         }
         FetchRecord(infile, line);
         iz = sscanf(line,"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", &ic[0], &ic[1],&ic[2],&ic[3],
            &ic[4],&ic[5],&ic[6],&ic[7], &ic[8], &ic[9],&ic[10],&ic[11], &ic[12], &ic[13],&ic[14],&ic[15] );
         
         for( i = 0; i < iz/2 ; i++ )
         {
             attach[i+j*8][0] = ic[2*i];
             attach[i+j*8][1] = ic[2*i+1];
         }
      }
   }
/* write atom information  */
/* could use routine here to reset atoms types to ensure that we have correct
   MM2 types   */

   for( i = 1; i <= niatom; i++ )
   {     
       FetchRecord(infile, line);
//       sscanf(line,"%f %f %f",&xtmp,&ytmp,&ztmp);
       xtmp = ReadMMFloat(line, 0,10);
       ytmp = ReadMMFloat(line,11,10);
       ztmp = ReadMMFloat(line,21,10);
       atomname[0] = line[30];
       atomname[1] = line[31];
       if (atomname[1] == ' ')
          atomname[1] = '\0';
       else 
          atomname[2] = '\0';
       typename[0] = line[32];
       typename[1] = line[33];
       typename[2] = line[34];
       typename[3] = '\0';
       sscanf(typename,"%d",&itype);
       
       newtype = mm3_mmxtype(itype);
       newatom = make_atom(newtype,xtmp,ytmp,ztmp,atomname);
       set_atomdata(newatom,newtype,itype,0,0,0);
       if (newtype == 200) // lp ty
          strcpy(atom[newatom].name,"Du");
   }

   fclose(infile);
      
   /*     generate iat array */
   for( i = 0; i < nconn; i++ )
   {
       for( j = 0; j < 15; j++ )
       {
          if( !iconn[j + 1][i] )
               break;
             make_bond((iconn[j][i])+ibotptr,(iconn[j+1][i])+ibotptr,1);
       }
   }
   /*     add terminal bonds */
   for( i = 0; i < nattch; i++ )
   {
           make_bond((attach[i][0])+ibotptr,(attach[i][1])+ibotptr,1);
   }
   
   free_imatrix( attach, 0, niatom, 0, 2);
   free_imatrix(iconn ,0, 16, 0, nconn+1);

   set_atomtypes(MMX);
}
void generate_mm2lists(int *nattch,int *nconn, int **attach, int **iconn)
{
   int i, ii, j,k,l,nbnds;
   int *term, *usedb;

   term = ivector(0 , 2*natom);
   usedb = ivector(0, 2*natom);
   for (i=0; i < 2*natom; i++)
   {
       term[i] = FALSE;
       usedb[i] = FALSE;
   }

      nbnds = 0;
      *nconn = 0;
      *nattch = 0;

      for( ii = 0; ii < 2; ii++ )
      {
        for( i = 0; i < natom; i++ )
              attach[i][ii] = 0;
      }
      for( i = 0; i < (natom/2); i++ )
      {
         for( j = 0; j < 16; j++ )
           iconn[j][i] = 0;
      }
      
      for( i = 1; i <= natom; i++ )
      {
         for( j = 0; j < MAXIAT; j++ )
         {
           if( atom[i].iat[j] != 0 && atom[i].bo[j] != 9 )
           {
              if( atom[i].iat[j] >= i )
              {
                  usedb[nbnds] = FALSE;
                  term[nbnds] = FALSE;
                  if( atom[i].iat[1] == 0 || atom[atom[i].iat[j]].iat[1] == 0 )
                      term[nbnds] = TRUE;
                  nbnds = nbnds + 1;
                  if (nbnds > 2*natom)
                     message_alert("nbnds array exceeded in generate lists","MMX File Setup");
               }
           }
        }
      }
 /*  find first bond not in use and not terminal */
L_21:
      nbnds = 0;
      for( i = 1; i <= natom; i++ )
      {
         for( j = 0; j < MAXIAT; j++ )
         {
            if( atom[i].iat[j] != 0 && atom[i].bo[j] != 9 )
            {
               if( atom[i].iat[j] >= i )
               {
                  if( usedb[nbnds] == FALSE && term[nbnds] == FALSE )
                  {
                      iconn[0][*nconn] = i;
                      iconn[1][*nconn] = atom[i].iat[j];
                      usedb[nbnds] = TRUE;
                      nbnds = nbnds + 1;
                      goto L_51055;
                   }
                   nbnds = nbnds + 1;
               }
            }
        }
        /* no more bonds to add quit */
    }
    goto L_36;
L_51055:
    l = 1;
    nbnds = 0;
  /*   search bond list for bonds to add */
    for( j = 1; j <= natom; j++ )
    {
       for( k = 0; k < MAXIAT; k++ )
       {
          if( atom[j].iat[k] != 0 && atom[j].bo[k] != 9 )
          {
             if( atom[j].iat[k] >= j )
             {
               if( usedb[nbnds] == FALSE && term[nbnds] == FALSE )
               {
                  if( iconn[l][*nconn] == j )
                  {
                     if( l == 15 )
                     {
                        *nconn= *nconn+1;
                        iconn[0][*nconn] = j;
                        l = 0;
                     }
                     l = l + 1;
                     iconn[l][*nconn] = atom[j].iat[k];
                     usedb[nbnds] = TRUE;
                  }
              }
              nbnds = nbnds + 1;
           }
        }
     }
   }
   nbnds = 0;
   for( j = 1; j <= natom; j++ )
   {
      for( k = 0; k < MAXIAT; k++ )
      {
         if( atom[j].iat[k] != 0 && atom[j].bo[k] != 9 )
         {
            if( atom[j].iat[k] >= j )
            {
              if( usedb[nbnds] == FALSE && term[nbnds] == FALSE )
              {
                 if( iconn[l][*nconn] == atom[j].iat[k])
                 {
                    if( l == 15 )
                    {
                      *nconn = *nconn + 1;
                      iconn[0][*nconn] = atom[j].iat[k];
                      l = 0;
                    }
                    l = l + 1;
                    iconn[l][*nconn] = j;
                    *nconn = *nconn + 1;
                    usedb[nbnds] = TRUE;
                    nbnds = nbnds + 1;
                    goto L_21;
                 }
              }
              nbnds = nbnds + 1;
            }
          }
       }
    }
    *nconn = *nconn + 1;
    goto L_21;
L_36:

/*         generate mmp2-style attachment lists */
    *nattch = 0;
    for( i = 1; i <= natom; i++ )
    {
       if( atom[i].type != 0 )
       {
         if( atom[i].iat[1] == 0 )
         {
            attach[*nattch][0] = atom[i].iat[0];
            attach[*nattch][1] = i;
            *nattch = *nattch + 1;
         }
       }
    }

    free_ivector(term, 0, 2*natom);
    free_ivector(usedb, 0, 2*natom);
}

