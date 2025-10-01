// Copyright (C) 2004-2016 Kevin E. Gilbert, dba Serena Software
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Kevin E. Gilbert
// Serena Software
// Box 3076
// Bloomington, IN 47402-3076
//
// gilbert@serenasoft.com

#include "pcwin.h"
#include "pcmod.h"
#include "bonds_ff.h"
#include "dipmom.h"

void charge_dipole(void)
{
    int i;
    double weight, xcenter, ycenter, zcenter;
    double xsum, ysum, zsum, debye;
    double cu;

    debye = 4.8033324;

    weight = 0.0;
    xcenter = 0.0;
    ycenter = 0.0;
    zcenter = 0.0;
    for (i=1; i <= natom; i++)
    {
        weight += atom[i].atomwt;
        xcenter += atom[i].x*atom[i].atomwt;
        ycenter += atom[i].y*atom[i].atomwt;
        zcenter += atom[i].z*atom[i].atomwt;
    }
    xcenter /= weight;
    ycenter /= weight;
    zcenter /= weight;

    xsum = 0.0;
    ysum = 0.0;
    zsum = 0.0;
    dmomv.xp = 0.0;
    dmomv.yp = 0.0;
    dmomv.xp = 0.0;
    dmomv.xn = 0.0;
    dmomv.yn = 0.0;
    dmomv.zn = 0.0;
    for(i=1; i <= natom; i++)
    {
        xsum += (atom[i].x-xcenter)*atom[i].charge;
        ysum += (atom[i].y-ycenter)*atom[i].charge;
        zsum += (atom[i].z-zcenter)*atom[i].charge;

        cu = atom[i].charge;
        if ( cu < 0.0)
        {
            dmomv.xn -= atom[i].x*cu;
            dmomv.yn -= atom[i].y*cu;
            dmomv.zn -= atom[i].z*cu;
        } else if (cu > 0.0)
        {
            dmomv.xp += atom[i].x*cu;
            dmomv.yp += atom[i].y*cu;
            dmomv.zp += atom[i].z*cu;
        }
        
    }
    dipolemom.xdipole = xsum*debye;
    dipolemom.ydipole = ysum*debye;
    dipolemom.zdipole = zsum*debye;
}

void dipole_dipole()
{
    int i,ia, ib;
    double xsum, ysum, zsum;
    double xi, yi, zi, fik, ri2;

    xsum = 0.0;
    ysum = 0.0;
    zsum = 0.0;
    for (i=0; i < bonds_ff.nbnd ; i++)
    {
        if (bonds_ff.bmom[i] != 0.0)
        {
            ia = bonds_ff.idpl[i][0];
            ib = bonds_ff.idpl[i][1];
            xi = atom[ib].x - atom[ia].x;
            yi = atom[ib].y - atom[ia].y;
            zi = atom[ib].z - atom[ia].z;
            ri2 = xi*xi + yi*yi + zi*zi;
            fik = bonds_ff.bmom[i]/sqrt(ri2);
            xsum += fik*xi;
            ysum += fik*yi;
            zsum += fik*zi;
        }
    }
    dipolemom.xdipole = xsum;
    dipolemom.ydipole = ysum;
    dipolemom.zdipole = zsum;

} 
/* ================================================= */
void vibcharge_dipole(double **acoord)
{
    int i;
    double weight, xcenter, ycenter, zcenter;
    double xsum, ysum, zsum, debye;
    double xy,xz,yx,yz,zx,zy;

    debye = 4.8033324;

    weight = 0.0;
    xcenter = 0.0;
    ycenter = 0.0;
    zcenter = 0.0;
    for (i=1; i <= natom; i++)
    {
        weight += atom[i].atomwt;
        xcenter += acoord[i][0]*atom[i].atomwt;
        ycenter += acoord[i][1]*atom[i].atomwt;
        zcenter += acoord[i][2]*atom[i].atomwt;
    }
    xcenter /= weight;
    ycenter /= weight;
    zcenter /= weight;

    xsum = 0.0;
    ysum = 0.0;
    zsum = 0.0;
    for(i=1; i <= natom; i++)
    {
        xsum += (acoord[i][0]-xcenter)*atom[i].charge;
        ysum += (acoord[i][1]-ycenter)*atom[i].charge;
        zsum += (acoord[i][2]-zcenter)*atom[i].charge;        
    }
    dipolemom.xdipole = xsum*debye;
    dipolemom.ydipole = ysum*debye;
    dipolemom.zdipole = zsum*debye;
// quadrapole moment
    xsum = 0.0;
    ysum = 0.0;
    zsum = 0.0;
    xy = xz = yx = yz = zx = zy = 0.0;
    for(i=1; i <= natom; i++)
    {
        xsum += (acoord[i][0]-xcenter)*(acoord[i][0]-xcenter)*atom[i].charge;
        ysum += (acoord[i][1]-ycenter)*(acoord[i][1]-ycenter)*atom[i].charge;
        zsum += (acoord[i][2]-zcenter)*(acoord[i][2]-zcenter)*atom[i].charge;        
        xy += (acoord[i][0]-xcenter)*(acoord[i][1]-ycenter)*atom[i].charge;
        xz += (acoord[i][0]-xcenter)*(acoord[i][2]-zcenter)*atom[i].charge;       
        yx += (acoord[i][1]-ycenter)*(acoord[i][0]-xcenter)*atom[i].charge;        
        yz += (acoord[i][1]-ycenter)*(acoord[i][2]-zcenter)*atom[i].charge;        
        zx += (acoord[i][2]-zcenter)*(acoord[i][0]-xcenter)*atom[i].charge;        
        zy += (acoord[i][2]-zcenter)*(acoord[i][1]-ycenter)*atom[i].charge;        
    }
    quadmom.xx = xsum;
    quadmom.xy = xy;
    quadmom.xz = xz;
    quadmom.yy = ysum;
    quadmom.yx = yx;
    quadmom.yz = yz;
    quadmom.zz = zsum;
    quadmom.zx = zx;
    quadmom.zy = zy;
    
}

void vibdipole_dipole(double **acoord)
{
    int i,ia, ib;
    double xsum, ysum, zsum,bmomj;
    double xi, yi, zi, fik, ri2;

    xsum = 0.0;
    ysum = 0.0;
    zsum = 0.0;
    for (i=0; i < bonds_ff.nbnd ; i++)
    {
         ia = bonds_ff.i12[i][0];
         ib = bonds_ff.i12[i][1];
         if (atom[ia].type == 1 && atom[ib].type == 5)
            bmomj = -0.25;
         else if (atom[ia].type == 1 && atom[ib].type == 36)
            bmomj = -0.25;
         else
            bmomj = bonds_ff.bmom[i];
         if (bmomj != 0.0)
         {
            xi = acoord[ib][0] - acoord[ia][0];
            yi = acoord[ib][1] - acoord[ia][1];
            zi = acoord[ib][2] - acoord[ia][2];
            ri2 = xi*xi + yi*yi + zi*zi;
            fik = bmomj/sqrt(ri2);
            xsum += fik*xi;
            ysum += fik*yi;
            zsum += fik*zi;
         }
    }
    dipolemom.xdipole = xsum;
    dipolemom.ydipole = ysum;
    dipolemom.zdipole = zsum;

} 
               
