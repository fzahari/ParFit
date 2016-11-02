#include "pcwin.h"
#include "pcmod.h"
#include "errno.h"
#include "files.h"
#include "read.h"

#define MAIN 0L
#define SUBSTRUCTURE 1L
#define MOPAC 0L
#define MNDO  1L

// ========================================
void clean_string(char *astring)
{
    int iz,i;
    iz = strlen(astring);
    for (i=0; i < iz; i++)
    {
        if (isprint(astring[i]) == 0)
          astring[i] = ' ';
    }
    for (i = iz -1; i >= 0; i--)
    {
        if (astring[i] != ' ')
        {
            astring[i+1] = '\0';
            break;
        }
    }
}
/*======================================================================*/
void * malloc_filename ( char *path , char *name )
{
  int ix, iz;
  char *filename;
  
  ix = 0;
  if (path != NULL)
     ix = strlen(path);
  iz = strlen(name);
  filename = malloc (ix + 1 + iz + 1) ;
  if ( ix  == 0 )
    {
      strcpy (filename,name);
    }
  else
    {
      sprintf (filename,"%s\\%s",path,name) ;
    }
  return ( (void *) filename ) ;
}
/*======================================================================*/
FILE * fopen_path ( char *path , char *name , char *mode )
{
  char *filename = malloc_filename(path,name) ;
  FILE *rval = fopen (filename,mode) ;
  if ( rval == NULL )
  {
     fprintf (pcmoutfile,"open of --%s-- failed\n",filename) ;
     fprintf(pcmoutfile,"Error %d %s\n",errno,strerror(errno));
     fclose(pcmoutfile);
     exit(0);
  }
  free (filename) ;
  return ( rval ) ;
}

/* ================================================== */
void    check_numfile(int ftype)
{
        int        igau, header, model;
        FILE       *infile;
        char      inputline[160];

        header = FALSE;
        model = FALSE;
        infile = fopen_path(Openbox.path,Openbox.fname,"r");
        if (infile == NULL)
        {
            message_alert("Error trying to check number of structures in file. Probably a bad filename.","Check Numfile");
            fprintf(pcmoutfile,"Error in check numfile:: %s %s\n",Openbox.path, Openbox.fname);
            files.nfiles = 0;
            return;
        }
        igau = 0;
        while ( fgets(inputline,159,infile) != NULL )
        {
            if (strncasecmp(inputline,"{PCM",4) == 0)
            {
                files.nfiles++;
            } else if (strncasecmp(inputline,"MMX ",4) == 0)
            {
                files.nfiles++;
            }else if (strncasecmp(inputline,"SDF ",4) == 0)
            {
                files.nfiles++;
            }
       }
        fclose(infile);
}
/* ------------------------------------*/
int get_bondorder(int ia,int ib)
{
    int i;
    for (i=0; i < MAXIAT; i++)
    {
        if (atom[ia].iat[i] == ib)
          return atom[ia].bo[i];
    }
    return 0;
}
/* ------------------------------------*/    
int is_bond(int ia, int ib)
{
    int i;
    for (i=0; i < MAXIAT; i++)
    {
        if (atom[ia].iat[i] == ib)
          return TRUE;
    }
    return FALSE;
}
/* ------------------------------------*/         
void mopaco(int istart,int ifinish)
{
     int i, ibi, ibj, it, j, jt, kk;
     float disij, dx, dy, dz, r, xii, yii, zii,rbdis[4];
     static float radii[110] =
        { 0.77, 0.77, 0.77, 0.77, 0.50, 0.66, 0.66, 0.70, 0.70, 0.70,
          0.64, 0.99, 1.14, 1.33, 1.04, 1.04, 1.04, 1.04, 1.17, 0.01,
          0.50, 0.77, 0.50, 0.50, 1.10, 0.88, 0.88, 0.50, 0.77, 0.77,
          1.22, 1.40, 1.20, 1.17, 1.37, 0.50, 0.70, 1.04, 1.17, 0.77,
          0.70, 0.66, 0.77, 0.70, 0.50, 0.66, 1.10, 0.77, 0.77, 0.77,
          0.77, 0.66, 0.70, 1.33, 0.66, 0.77, 0.77, 0.70, 0.01, 0.50,
          0.77, 0.77, 0.77, 0.70, 0.66, 0.70, 1.10, 0.66, 0.70, 0.50,
          0.77, 0.77, 0.77, 0.77, 0.50, 0.66, 0.66, 0.70, 0.70, 0.70,
          0.77, 0.77, 0.77, 0.77, 0.50, 0.66, 0.66, 0.70, 0.70, 0.70,
          0.77, 0.77, 0.77, 0.77, 0.50, 0.66, 0.66, 0.70, 0.70, 0.70,
          0.77, 0.77, 0.77, 0.77, 0.50, 0.66, 0.66, 0.70, 0.70, 0.70,};

        /*     SUBROUTINE SENDS BACK THE ACTUAL NUMBER OF ATOMS, AND
         *     DETERMINES THE BOND CONNECTIONS and bond orders TO ESTABLISH THE attached
         *     ATOM LIST. ATOM TYPES HAVE BEEN CHANGED TO MM2 TYPES.  THOSE THAT
         *     AREN'T IN MM2 ARE GIVEN #S = AT # + 60
         *     BOND ORDERS AND THE IATMS AND IBNDS LISTS ARE GENERATED
         *     NATOMS IS THE # OF ENTRIES IN THE Z MATRIX */

        rbdis[0] = 1.42;
        rbdis[1] = 1.25;
        rbdis[2] = 1.37;
        rbdis[3] = 1.25;

        /*  find all single bonds based on sum of covalent radii */
        for( i = istart; i <= (ifinish - 1L); i++ )
        {
                it = atom[i].mmx_type;
                for( j = i + 1L; j <= ifinish; j++ )
                {
                    if (!is_bond(i,j))
                    {
                        jt = atom[j].mmx_type;
                        dx = atom[i].x - atom[j].x;
                        dy = atom[i].y - atom[j].y;
                        dz = atom[i].z - atom[j].z;
                        r = sqrt( dx*dx + dy*dy + dz*dz );
                        if( it >= 80L && jt >= 80L )
                        {
                                if( r <= (atom[i].radius + atom[j].radius) )
                                        make_bond( i, j, 1 );
                        }else if( it >= 80L && jt < 80L )
                        {
                                if( r <= (atom[i].radius + radii[atom[j].mmx_type-1]) )
                                        make_bond( i, j, 1);

                        }else if( it < 80L && jt >= 80L )
                        {
                                if( r <= (radii[atom[i].mmx_type-1] + atom[j].radius + .1) )
                                        make_bond( i, j, 1 );
                        }else if( r <= (radii[atom[i].mmx_type-1] + radii[atom[j].mmx_type-1]+ .1) )
                        {
                                if( i < j )
                                        make_bond( i, j, 1 );
                        }
                    }
                 }
        }

        /*   find all double and triple bonds - search only types for
         *  carbon, oxygen, nitrogen, sulfur, phosphorous */
        for( i = istart; i <= (ifinish - 1); i++ )
        {
            it = atom[i].atomnum;
            if( it == 6 || it == 7 || it == 8 || it == 15 || it == 16 )
            {
                for( j = 0; j < MAXIAT; j++ )
                {
                     if( atom[i].iat[j] != 0 )
                     {
                        jt = atom[atom[i].iat[j]].atomnum;
                        if( jt == 6 || jt == 7 || jt == 8 || jt == 15 || jt == 16 )
                        {
                            xii = atom[i].x - atom[atom[i].iat[j]].x;
                            yii = atom[i].y - atom[atom[i].iat[j]].y;
                            zii = atom[i].z - atom[atom[i].iat[j]].z;
                            if( fabs( xii ) <= 2.3 )
                            {
                               if( fabs( yii ) <= 2.3 )
                               {
                                  if( fabs( zii ) <= 2.3 )
                                  {
                                     disij = sqrt( xii*xii + yii*yii + zii*zii );
                /*   c=o bond */
                                     if( ((it == 8 && jt == 6) || (it == 6 && jt == 8)) && disij <=  rbdis[3] + 0.05 )
                                     {
                                         if( i < atom[i].iat[j] ) 
                                            make_bond( i, atom[i].iat[j],1 );
                /*  alkyne */
                                     }else if( (it == 6 && jt == 6) &&  disij <= rbdis[1] )
                                     {
                                        if( i < atom[i].iat[j] )
                                            make_bond( i, atom[i].iat[j],2 );
                 /*  c=n */
                                     }else if( ((it == 7 && jt == 6) || (it == 6 && jt == 7)) 
                                                                  && disij <= rbdis[2] )
                                     {
                                        ibi = 0;
                                        for( kk = 0; kk < MAXIAT; kk++ )
                                        {
                                           if( atom[i].bo[kk] !=  9 )
                                              ibi = ibi + atom[i].bo[kk];
                                        }
                                        ibj = 0;
                                        for( kk = 0; kk < MAXIAT; kk++ )
                                        {
                                           if( atom[atom[i].iat[j]].bo[kk]  != 9 )
                                              ibj = ibj + atom[atom[i].iat[j]].bo[kk];
                                        }
                                        if( it == 7 )
                                        {
                                           if( ibi <= 2 && ibj <= 3 )
                                           {
                                             if( i < atom[i].iat[j] )
                                                 make_bond( i, atom[i].iat[j], 1 );
                                             break;
                                           }
                                         }else if( jt == 7 )
                                         {
                                            if( ibj <= 2 && ibi <= 3 )
                                            {
                                               if( i < atom[i].iat[j] )
                                                  make_bond( i, atom[i].iat[j], 1 );
                                               break;
                                            }
                                         }
     /*nitrile */                    }else if( ((it == 7 && jt == 6) || (it == 6 && jt == 7)) && disij <= 1.18 )
                                     {
                                          if( i < atom[i].iat[j] )
                                              make_bond( i, atom[i].iat[j], 2 );
     /* N=O */                       }else if ( ((it == 7 && jt == 8) || (it == 8 && jt == 7)) && disij <= 1.26)
                                     {
                                         ibi = 0;
                                         if (it == 7)
                                         {
                                             for (kk = 0; kk < MAXIAT; kk++)
                                             {
                                                 if (atom[i].bo[kk] != 9)
                                                   ibi += atom[i].bo[kk];
                                             }
                                         } else
                                         {
                                             for (kk = 0; kk < MAXIAT; kk++)
                                             {
                                                 if (atom[atom[i].iat[j]].bo[kk] != 9)
                                                   ibi += atom[atom[i].iat[j]].bo[kk];
                                             }
                                         }
                                         if (ibi <= 3)
                                         {
                                            if ( i < atom[i].iat[j])
                                                make_bond(i,atom[i].iat[j],1);
                                            break;
                                         }
     /* c=c bond */                  }else if( (it == 6 && jt == 6) && disij <= rbdis[0] )
                                     {
                                         ibi = 0;
                                         for( kk = 0; kk < MAXIAT; kk++ )
                                         {
                                            if( atom[i].bo[kk] != 9 )
                                                  ibi = ibi + atom[i].bo[kk];
                                         }
                                         ibj = 0;
                                         for( kk = 0; kk < MAXIAT; kk++ )
                                         {
                                            if( atom[atom[i].iat[j]].bo[kk] != 9 )
                                               ibj = ibj + atom[atom[i].iat[j]].bo[kk];
                                         }
                                         if( ibi <= 3 && ibj <= 3 )
                                         {
                                            if( i < atom[i].iat[j] )
                                               make_bond( i, atom[i].iat[j], 1 );
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
        type();
        return;
} 
/* =============================================  */
int FetchRecord(FILE *fp, char *buffer)
{
      int ch;
      char *ptr;
      
      ptr = buffer;
      do 
      {
         ch = getc(fp);
         if (ch >= ' ')
            *ptr++ = ch;
         else if (ch == '\n')
         {
            *ptr = 0;
            return TRUE;
         } else if (ch == '\r')
         {
            ch = getc(fp);
            if (ch != '\n')
                ungetc(ch,fp);
            *ptr = 0;
            return TRUE;
         } else if (ch == EOF)
         {
             *ptr = 0;
             return (ptr != buffer);
         } else
            *ptr++ = ch;
      } while (ptr < buffer+255);
      
      do 
      {
           ch = getc(fp);
       } while (ch !='\n' && ch != '\r' && ch != EOF);
       if (ch == '\r')
       {
           ch = getc(fp);
           if (ch != '\n')
              ungetc(ch,fp);
        }
        *ptr = 0;
        return TRUE;
} 



        
