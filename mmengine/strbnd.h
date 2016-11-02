#ifndef STRBND_H
#define STRBND_H

struct t_strbnd {
        int nstrbnd, isb[MAXANG][3];
        float ksb1[MAXANG],ksb2[MAXANG];
        } strbnd;

#endif

void estrbnd(void);
void estrbnd1(void);
void estrbnd2(int);
void kstrbnd(void);
