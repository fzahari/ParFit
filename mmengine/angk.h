#ifndef ANGK_H
#define ANGK_H

struct t_angk1 {
        int use_ang3, use_ang4, use_ang5;
        int nang, nang3, nang4, nang5;
        int ndel, ndel3, ndel4;
        char  ktype[MAXANGCONST][10],ktype3[MAXANG3CONST][10],ktype4[MAXANG4CONST][10],ktype5[MAXANG5CONST][10];
        char  kdel[MAXANGDEL][10],kdel3[MAXANG3DEL][10],kdel4[MAXANG4DEL][10];
        float con[MAXANGCONST], ang[MAXANGCONST][3];
        float con3[MAXANG3CONST], ang3[MAXANG3CONST][3];        
        float con4[MAXANG4CONST], ang4[MAXANG4CONST][3];        
        float con5[MAXANG5CONST], ang5[MAXANG5CONST][3];
        float condel[MAXANGDEL],  angdel[MAXANGDEL][3];      
        float condel3[MAXANG3DEL], angdel3[MAXANG3DEL][3];      
        float condel4[MAXANG4DEL], angdel4[MAXANG4DEL][3];      
         } angk1;

struct t_angkp1 {
        int  npiang;
        char kpa[MAXPIANGCONST][10];
        float pacon[MAXPIANGCONST], panat[MAXPIANGCONST][3], piobp[MAXPIANGCONST];
        } angkp1;

struct t_angf {
        int use_angf, nfang;
        char kftype[MAXANGCONST][10];
        float fcon[MAXANGCONST],fc0[MAXANGCONST],fc1[MAXANGCONST],fc2[MAXANGCONST];
        float fc3[MAXANGCONST],fc4[MAXANGCONST],fc5[MAXANGCONST],fc6[MAXANGCONST];
        } angf;

#endif
