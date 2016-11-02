#ifndef SOLVENT_H
#define SOLVENT_H

#define STILL 1
#define HCT   2

struct t_solvent {
    int type;
    double doffset, p1,p2,p3,p4,p5;
    double *shct,*asolv,*rsolv,*vsolv,*gpol,*rborn;
    } solvent;

#endif

void esolv(void);
void esolv1(void);
void esolv2(int);
void born1(void);
void born(void);
void ksolv(void);
int find_bond(int,int);
void message_alert(char *, char *);
