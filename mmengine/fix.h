// fix distantce data 
#ifndef FIX_H
#define FIX_H

#define MAXDIS 300

struct t_fixdis {
        int FixBond;
        int nfxstr, ifxstr[MAXDIS][2];
        float fxstrd[MAXDIS], fxstrc[MAXDIS];
        int fixatm;
        }       fixdis;

// fix angle data
struct t_fixangle {
        int nfxangle, ifxangle[10][3];
        float anat[10], acon[10];
        } fixangle;
        
// fix torsion data
struct t_fxtor {
        int nfxtor, iatom[10][4];
        int start_ang[10], curr_ang[10], final_ang[10];
        int step[10];
        } fxtor;
#endif


int find_bond(int,int);
int find_fixed_bond(int,int);
void egeom(void);
void egeom1(void);
void egeom2(int);
