#ifndef ANGLES_H
#define ANGLES_H

#define HARMONIC  1
#define FOURIER   2

#define MAXMETALANG 2000

struct t_angles {
        int nang, i13[MAXANG][4], index[MAXANG];
        int nopb, angin[MAXANG],angtype[MAXANG];
        float acon[MAXANG], anat[MAXANG], copb[MAXANG];
        float fcon[MAXANG],c0[MAXANG],c1[MAXANG],c2[MAXANG],c3[MAXANG],c4[MAXANG],c5[MAXANG],c6[MAXANG];
        } angles;

struct t_high_coord {
        int ncoord, i13[MAXMETALANG][3];
        } high_coord;
#endif

void message_alert(char *, char *);
void get_angles(void);
int isangle(int, int);
