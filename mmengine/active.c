#include "pcwin.h"
#include "pcmod.h"
#include "active.h"

void set_active()
{
   int i;
   long int minflag;

   minflag = (1L << MIN_MASK);
   
   for (i=1; i <= natom; i++)
   {
       if (atom[i].flags & minflag)
          atom[i].use = FALSE;
       else
          atom[i].use = TRUE;
   }
} 
