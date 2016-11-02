#ifndef HESS_H
#define HESS_H

double hesscut;
    
struct t_hess {
    float **hessx, **hessy, **hessz;
    } hess; 
   
#endif 

void piseq(int, int);
void hessian(double *,int *, int *, int *,double *);
void ebond2(int);
void eangle2(int);
void estrbnd2(int);
void eangang2(int);
void eopbend2(int);
void etorsion2(int);
void estrtor2(int);
void ecoord_bond2(int);
void echarge2(int);
void edipole2(int);
void ehal2(int);
void ebufcharge2(int);
void eopbend_wilson2(int);      
void esolv2(int);
void ehigh_coord2(int);
void message_alert(char *,char *);


