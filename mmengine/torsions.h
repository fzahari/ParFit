#ifndef TORSIONS_H
#define TORSIONS_H

struct t_torsions {
        int ntor, **i14;
        float *v1,*v2,*v3,*v4,*v5,*v6;
        int   *ph1,*ph2,*ph3,*ph4,*ph5,*ph6;
        } torsions;

#endif

int isbond(int, int);
int isangle(int,int);
int istorsion(int, int);
void get_torsions(void);
int is_metal(char *);
int is_linear(int);
double dihdrl(int,int,int,int);
void etorsion(void);
void etorsion1(void);
void etorsion2(int);
int isbond(int,int);
int is_delocalbond(int,int);
void numeral(int,char *,int);
void angle(int,int,int,float *);
void four(char *,char *,char *,char *, char *);
int is_ring42(int,int);
int is_ring54(int, int, int, int);
void message_alert(char *, char *);
void  transomeg(float *, float *, float *, int *, int *, char *, char *, char *, char *,
         char *, char *, char *, char *, char *, char *, char *);
