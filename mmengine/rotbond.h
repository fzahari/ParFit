#ifndef ROTBOND_H
#define ROTBOND_H

struct t_rotbond {
    int nbonds, bond[200][2], bond_res[200];
    int incl_amide, incl_alkene;
    } rotbond;

#endif

void find_rings(void);
void find_bonds(void);
void hdel(int);
void hadd(void);
int check_bond(int,int,int *);
void message_alert(char *, char *);
void matident(float (*)[4]);
void matrotsc(float (*)[4], int*, float*, float*);
void matrotat(float (*)[4], int*, float*);
void mattrans(float (*)[4], float*, float*, float*);
void matriply(float (*)[4], float (*)[4], float (*)[4]);
void matxform(float (*)[4], int);
void rotb(int, int, int, int, float);
void chain(int, int, int);
double dihdrl(int,int,int,int);
