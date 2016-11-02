#ifndef CROSSTERM_H
#define CROSSTERM_H

struct t_crossterm_k {
        int nangang, nstrbnd, nstrtor;
        int ang_ang[MAXAA], stbnindex[MAXSTBN] ;
        char str_tor[MAXSTRTOR][7], stbn[MAXSTBN][10];
        float aacon[MAXAA][3], stbncon[MAXSTBN][3], str_torcon[MAXSTRTOR];
        } crossterm_k;

#endif
