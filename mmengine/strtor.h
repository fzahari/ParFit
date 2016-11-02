#ifndef STRTOR_H
#define STRTOR_H

struct t_strtor {
    int nstrtor,ist[MAXTOR][2];
    float stcon[MAXTOR];
    }   strtor;

#endif

void estrtor(void);
void estrtor1(void);
void estrtor2(int);
void kstrtor(void);
