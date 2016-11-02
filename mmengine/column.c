#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "derivs.h"
#include "column.h"


void column(int nvar, int *h_init, int *h_stop,int *h_index,
             int *c_init,int *c_stop,int *c_index, int *c_value)
{
    int i,j,k,m;

    for (i=0; i < nvar; i++)
    {
        c_init[i] = 0;
        c_stop[i] = 0;
    }
// count number of elements in each column
    for (i=0; i < nvar; i++)
    {
        for (j=h_init[i]; j < h_stop[i]; j++)
        {
            k = h_index[j];
            c_stop[k]++;
        }
    }
// set each beginning marker past the last element for its column
    c_init[0] = c_stop[0] + 1;
    for (i=1; i <nvar; i++)
    {
        c_init[i] = c_init[i-1] + c_stop[i];
    }
// set column index scanning rows in reverse order
    for (i = nvar-1; i >= 0; i--)
    {
        for (j=h_init[i]; j < h_stop[i]; j++)
        {
            k = h_index[j];
            m = c_init[k]-1;
            c_init[k] = m;
            c_index[m] = i;
            c_value[m] = j;
        }
    }
    for(i=0; i < nvar; i++)
    {
       c_stop[i] = c_init[i] + c_stop[i]-1;
    }
}
            
        
