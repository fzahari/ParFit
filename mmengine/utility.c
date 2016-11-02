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
                        
