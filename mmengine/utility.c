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
#include "utility.h"
void pcmfout(int);

void message_alert(char *, char *);

/*   Routines to allocate and free space  */
/*  ----------------------------------  */
float *vector( int nl, int nh)
{
        float *v;

        v = (float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
                if (!v)
        {
             message_alert("Error allocation fmemory. Please contact Serena Software. The program will now exit!","Utility");
             fprintf(pcmoutfile,"error allocating float vector\n");
             pcmfout(2);
             exit(0);
        }
        return v-nl;
}
/*  ----------------------------------  */
 int *ivector( int nl, int nh)
{
         int *v;

        v = ( int *)malloc((unsigned) (nh-nl+1)*sizeof( int));
        if (!v)
        {
             message_alert("Error allocation imemory. Please contact Serena Software. The program will now exit!","Utility");
             fprintf(pcmoutfile,"error allocating  int vector\n");
             pcmfout(2);
             exit(0);
        }
        return v-nl;
}
/*  ----------------------------------  */
 long int *ilvector( int nl, int nh)
{
        long int *v;

        v = (long int *)malloc((unsigned) (nh-nl+1)*sizeof(long int));
        if (!v)
        {
             message_alert("Error allocation imemory. Please contact Serena Software. The program will now exit!","Utility");
             fprintf(pcmoutfile,"error allocating  int vector\n");
             pcmfout(2);
             exit(0);
        }
        return v-nl;
}
/*  ----------------------------------  */
double *dvector( int nl, int nh)
{
        double *v;

        v = (double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
        if (!v)
        {
            char hatext[128];
            sprintf(hatext,"Error allocation dmemory. Size %d\n",nh);
            message_alert(hatext,"Utility");
            fprintf(pcmoutfile,"error allocating double vector\n");
             pcmfout(2);
            exit(0);
        }
        return v-nl;
}
/*  ----------------------------------  */

float **matrix( int nrl,  int nrh,  int ncl,  int nch)
{
         int i;
        float **m;

        m = (float **)malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
        if (!m)
        {
             message_alert("Error allocation farray. Please contact Serena Software. The program will now exit!","Utility");
             fprintf(pcmoutfile,"error allocating float array\n");
             pcmfout(2);
             exit(0);
        }
        m -= nrl;

        for (i = nrl; i <= nrh; i++)
        {
                m[i] = (float *)malloc((unsigned) (nch-ncl+1)*sizeof(float));
                if (!m[i])
                {
                        message_alert("Error allocation farray. Please contact Serena Software. The program will now exit!","Utility");
                        fprintf(pcmoutfile,"error allocating float array\n");
                        pcmfout(2);
                        exit(0);
                }
                m[i] -= ncl;
        }
        return m;
}
/*  ----------------------------------  */

double **dmatrix( int nrl,  int nrh,  int ncl,  int nch)
{
         int i;
        double **m;

        m = (double **)malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m)
        {
             message_alert("Error allocation darray. Please contact Serena Software. The program will now exit!","Utility");
             fprintf(pcmoutfile,"error allocating double array\n");
             pcmfout(2);
             exit(0);
        }
        m -= nrl;

        for (i = nrl; i <= nrh; i++)
        {
                m[i] = (double *)malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i])
                {
                        message_alert("Error allocation darray. Please contact Serena Software. The program will now exit!","Utility");
                        fprintf(pcmoutfile,"error allocating double array\n");
                        pcmfout(2);
                        exit(0);
                }
                m[i] -= ncl;
        }
        return m;
}
/*  ----------------------------------  */

 int **imatrix( int nrl,  int nrh,  int ncl,  int nch)
{
         int i, **m;

        m = ( int **)malloc((unsigned) (nrh-nrl+1)*sizeof( int*));
        if (!m)
        {
             message_alert("Error allocation iarray. Please contact Serena Software. The program will now exit!","Utility");
             fprintf(pcmoutfile,"error allocating array\n");
             pcmfout(2);
             exit(0);
        }
        m -= nrl;

        for (i = nrl; i <= nrh; i++)
        {
                m[i] = ( int *)malloc((unsigned) (nch-ncl+1)*sizeof( int));
                if (!m[i])
                {
                     message_alert("Error allocation farray. Please contact Serena Software. The program will now exit!","Utility");
                     fprintf(pcmoutfile,"error allocating array\n");
                     pcmfout(2);
                     exit(0);
                }
                m[i] -= ncl;
        }
        return m;
}
/*  ----------------------------------  */
void free_vector(float *v,  int nl,  int nh)
{
        free((char*) (v+nl));
}

void free_ivector( int *v,  int nl,  int nh)
{
        free((char*) (v+nl));
}
void free_ilvector( long int *v,  int nl,  int nh)
{
        free((char*) (v+nl));
}

void free_dvector(double *v,  int nl,  int nh)
{
        free((char*) (v+nl));
}

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
         int i;
        for (i=nrh; i >= nrl; i--)
        {
                free((char*) (m[i]+ncl));
        }
        free((char*) (m+nrl));
}

void free_imatrix( int **m, int nrl, int nrh, int ncl, int nch)
{
         int i;
        for (i=nrh; i >= nrl; i--)
        {
                free((char*) (m[i]+ncl));
        }
        free((char*) (m+nrl));
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
         int i;
        for (i=nrh; i >= nrl; i--)
        {
                free((char*) (m[i]+ncl));
        }
        free((char*) (m+nrl));
}
                        
