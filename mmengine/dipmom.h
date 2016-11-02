#ifndef DIPMOM_H
#define DIPMOM_H

struct t_dipolemom {
        double total, xdipole, ydipole, zdipole;
       }  dipolemom;
struct t_dmomv {
        float xn,yn,zn,xp,yp,zp;
        }       dmomv;
struct t_quadmom {
        double xx,xy,xz,yx,yy,yz,zx,zy,zz;
       } quadmom;
 #endif

void charge_dipole(void);
void dipole_dipole(void);
void vibcharge_dipole(double **);
void vibdipole_dipole(double **);
