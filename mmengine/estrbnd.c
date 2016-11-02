#include "pcwin.h"
#include "pcmod.h"

#include "energies.h"
#include "angles.h"
#include "bonds_ff.h"
#include "hess.h"
#include "derivs.h"
#include "units.h"
#include "minim_values.h"
#include "strbnd.h"
       
void estrbnd()
{
      int i,ij,j,k,ia,ib,ic;
      double e,dr,dt,angle,dot,cosine,dr1;
      double xia,yia,zia,xib,yib,zib;
      double xic,yic,zic;
      double rab,xab,yab,zab;
      double rcb,xcb,ycb,zcb;

      energies.estrbnd = 0;
      if (minim_values.iprint && strbnd.nstrbnd > 0)
      {
          fprintf(pcmoutfile,"\nStretch-Bend Terms\n");
          fprintf(pcmoutfile,"          At1     At2    At3      Angle     Anat   StbnCst1   StbnCst2  Estb\n");
      }
      
      for (i=0; i < strbnd.nstrbnd; i++)
      {
          ij = strbnd.isb[i][0];
          ia = angles.i13[ij][0];
          ib = angles.i13[ij][1];
          ic = angles.i13[ij][2];
          if (atom[ia].use || atom[ib].use || atom[ic].use)
          {
            xia = atom[ia].x;
            yia = atom[ia].y;
            zia = atom[ia].z;
            xib = atom[ib].x;
            yib = atom[ib].y;
            zib = atom[ib].z;
            xic = atom[ic].x;
            yic = atom[ic].y;
            zic = atom[ic].z;

            xab = xia - xib;
            yab = yia - yib;
            zab = zia - zib;
            rab = sqrt(xab*xab + yab*yab + zab*zab);
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
            if (rab*rcb != 0.0)
            {
               dot = xab*xcb + yab*ycb + zab*zcb;
               cosine = dot / (rab*rcb);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);
               dt = angle - angles.anat[ij];
               j = strbnd.isb[i][1];
               k = strbnd.isb[i][2];
               dr = 0.0;
               dr1 = 0.0;
               if (j > -1) dr = (rab - bonds_ff.bl[j])*strbnd.ksb1[i];
               if (k > -1) dr1 = (rcb - bonds_ff.bl[k])*strbnd.ksb2[i];
               e = units.stbnunit*dt*(dr+dr1);
               energies.estrbnd += e;
               atom[ia].energy += e;
               atom[ib].energy += e;
               atom[ic].energy += e;

               if (minim_values.iprint)
                 fprintf(pcmoutfile,"StrBnd: %2s(%-3d)-%2s(%-3d)-%2s(%-3d) : %-8.3f %-8.3f %-8.3f %-8.3f = %-8.4f\n",
                    atom[ia].name, ia,atom[ib].name, ib,atom[ic].name, ic,angle,angles.anat[ij],
                    strbnd.ksb1[i],strbnd.ksb2[i], e);
            }
          }
      }
}

void estrbnd1()
{
      int i,ij,j,k,ia,ib,ic;
      double e,dr,dr1,dt,angle,dot,cosine;
      double xia,yia,zia,xib,yib,zib;
      double xic,yic,zic;
      double rab,xab,yab,zab;
      double rcb,xcb,ycb,zcb;
      double xp,yp,zp,rp;
      double dedxia,dedyia,dedzia;
      double dedxib,dedyib,dedzib;
      double dedxic,dedyic,dedzic;
      double term,terma,termc, rab2, rcb2;
      double ksb1, ksb2;
      double dtemp;

      energies.estrbnd = 0;
      for (i=0; i <= natom; i++)
      {
          deriv.destb[i][0] = 0.0;
          deriv.destb[i][1] = 0.0;
          deriv.destb[i][2] = 0.0;
      }

      for (i=0; i < strbnd.nstrbnd; i++)
      {
          ij = strbnd.isb[i][0];
          ia = angles.i13[ij][0];
          ib = angles.i13[ij][1];
          ic = angles.i13[ij][2];
          if (atom[ia].use || atom[ib].use || atom[ic].use)
          {
            xia = atom[ia].x;
            yia = atom[ia].y;
            zia = atom[ia].z;
            xib = atom[ib].x;
            yib = atom[ib].y;
            zib = atom[ib].z;
            xic = atom[ic].x;
            yic = atom[ic].y;
            zic = atom[ic].z;

            xab = xia - xib;
            yab = yia - yib;
            zab = zia - zib;
            rab2 = xab*xab + yab*yab + zab*zab;
            rab = sqrt(xab*xab + yab*yab + zab*zab);
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb;
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
            xp = ycb*zab - zcb*yab;
            yp = zcb*xab - xcb*zab;
            zp = xcb*yab - ycb*xab;
            rp = sqrt(xp*xp + yp*yp + zp*zp);
            if (rp != 0.0)
            {
               dot = xab*xcb + yab*ycb + zab*zcb;
               cosine = dot / (rab*rcb);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);
               dt = angle - angles.anat[ij];

               dtemp = dot*dot/(rab2*rcb2);
               if (fabs(dtemp) >= 1.0)
                  dtemp = 0.9999*dtemp/fabs(dtemp);
               termc = sqrt(1.00 - dtemp);
               terma = rab*rcb;
               j = strbnd.isb[i][1];
               k = strbnd.isb[i][2];
               term = -units.stbnunit*radian; 
               if (j < 0)
               {
                  dr = 0.0;
                  ksb1 = 0.0;
               } else
               {
                  dr = strbnd.ksb1[i]*(rab - bonds_ff.bl[j]);
                  ksb1 = strbnd.ksb1[i];
               }
               if (k < 0)
               {
                  dr1 = 0.0;
                  ksb2 = 0.0;
               }  else
               {
                  dr1 = strbnd.ksb2[i]*(rcb - bonds_ff.bl[k]);
                  ksb2 = strbnd.ksb2[i];
               }
               e = units.stbnunit*dt*(dr+dr1);
               dedxia = (term *(xcb/terma - (xab*dot)/(rab2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb1*xab)/rab;
               dedyia = (term *(ycb/terma - (yab*dot)/(rab2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb1*yab)/rab;
               dedzia = (term *(zcb/terma - (zab*dot)/(rab2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb1*zab)/rab;

               dedxic = (term *(xab/terma - (xcb*dot)/(rcb2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb2*xcb)/rcb;
               dedyic = (term *(yab/terma - (ycb*dot)/(rcb2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb2*ycb)/rcb;
               dedzic = (term *(zab/terma - (zcb*dot)/(rcb2*terma))*(dr+dr1))/ termc + (units.stbnunit*dt*ksb2*zcb)/rcb;

               dedxib = ( (term /(rab*rcb))*( (xcb*dot)/rcb2 + (xab*dot)/rab2 + (2.0*xib-xia-xic))*(dr+dr1))/termc +
                             units.stbnunit*dt*( -(ksb1*xab)/rab - ksb2*xcb/rcb);
               dedyib = ( (term /(rab*rcb))*( (ycb*dot)/rcb2 + (yab*dot)/rab2 + (2.0*yib-yia-yic))*(dr+dr1))/termc +
                             units.stbnunit*dt*( -(ksb1*yab)/rab - ksb2*ycb/rcb);
               dedzib = ( (term /(rab*rcb))*( (zcb*dot)/rcb2 + (zab*dot)/rab2 + (2.0*zib-zia-zic))*(dr+dr1))/termc +
                             units.stbnunit*dt*( -(ksb1*zab)/rab - ksb2*zcb/rcb);
              
               energies.estrbnd += e;
               deriv.destb[ia][0] += dedxia;
               deriv.destb[ia][1] += dedyia;
               deriv.destb[ia][2] += dedzia;

               deriv.destb[ib][0] += dedxib;
               deriv.destb[ib][1] += dedyib;
               deriv.destb[ib][2] += dedzib;

               deriv.destb[ic][0] += dedxic;
               deriv.destb[ic][1] += dedyic;
               deriv.destb[ic][2] += dedzic;
               
               virial.virx += xab*dedxia + xcb*dedxic;
               virial.viry += yab*dedyia + ycb*dedyic;
               virial.virz += zab*dedzia + zcb*dedzic;
            }
          }
      }
}

void estrbnd2(int iatom)
{
    int i,j,k,ij;
    int ia,ib,ic;
    double dt,dr,angle,dot,cosine;
    double xia,yia,zia;
    double xib,yib,zib;
    double xic,yic,zic;
    double xab,yab,zab,rab;
    double xcb,ycb,zcb,rcb;
    double xp,yp,zp,rp,rp2;
    double terma,termc;
    double xrab,yrab,zrab,rab2;
    double xrcb,yrcb,zrcb,rcb2;
    double xabp,yabp,zabp;
    double xcbp,ycbp,zcbp;
    double ddtdxia,ddtdyia,ddtdzia;
    double ddtdxib,ddtdyib,ddtdzib;
    double ddtdxic,ddtdyic,ddtdzic;
    double ddrdxia,ddrdyia,ddrdzia;
    double ddrdxib,ddrdyib,ddrdzib;
    double ddrdxic,ddrdyic,ddrdzic;
    double dtxiaxia,dtxiayia,dtxiazia;
    double dtxibxib,dtxibyib,dtxibzib;
    double dtxicxic,dtxicyic,dtxiczic;
    double dtyiayia,dtyiazia,dtziazia;
    double dtyibyib,dtyibzib,dtzibzib;
    double dtyicyic,dtyiczic,dtziczic;
    double dtxibxia,dtxibyia,dtxibzia;
    double dtyibxia,dtyibyia,dtyibzia;
    double dtzibxia,dtzibyia,dtzibzia;
    double dtxibxic,dtxibyic,dtxibzic;
    double dtyibxic,dtyibyic,dtyibzic;
    double dtzibxic,dtzibyic,dtzibzic;
    double dtxiaxic,dtxiayic,dtxiazic;
    double dtyiaxic,dtyiayic,dtyiazic;
    double dtziaxic,dtziayic,dtziazic;
    double drxiaxia,drxiayia,drxiazia;
    double drxibxib,drxibyib,drxibzib;
    double drxicxic,drxicyic,drxiczic;
    double dryiayia,dryiazia,drziazia;
    double dryibyib,dryibzib,drzibzib;
    double dryicyic,dryiczic,drziczic;
    
      for (i=0; i < strbnd.nstrbnd; i++)
      {
          ij = strbnd.isb[i][0];
          ia = angles.i13[ij][0];
          ib = angles.i13[ij][1];
          ic = angles.i13[ij][2];
          if (iatom == ia || iatom == ib || iatom == ic)
          {
            xia = atom[ia].x;
            yia = atom[ia].y;
            zia = atom[ia].z;
            xib = atom[ib].x;
            yib = atom[ib].y;
            zib = atom[ib].z;
            xic = atom[ic].x;
            yic = atom[ic].y;
            zic = atom[ic].z;
            xab = xia - xib;
            yab = yia - yib;
            zab = zia - zib;
            rab = sqrt(xab*xab + yab*yab + zab*zab);
            xcb = xic - xib;
            ycb = yic - yib;
            zcb = zic - zib;
            rcb = sqrt(xcb*xcb + ycb*ycb + zcb*zcb);
            xp = ycb*zab - zcb*yab;
            yp = zcb*xab - xcb*zab;
            zp = xcb*yab - ycb*xab;
            rp = sqrt(xp*xp + yp*yp + zp*zp);
            if (rp != 0.0)
            {
              dot = xab*xcb + yab*ycb + zab*zcb;
              cosine = dot / (rab*rcb);
               if (cosine < -1.0)
                 cosine = -1.0;
               if (cosine > 1.0)
                 cosine = 1.0;
               angle = radian * acos(cosine);

               dt = angle - angles.anat[ij];
               terma = -radian / (rab*rab*rp);
               termc = radian / (rcb*rcb*rp);
               ddtdxia = terma * (yab*zp-zab*yp);
               ddtdyia = terma * (zab*xp-xab*zp);
               ddtdzia = terma * (xab*yp-yab*xp);
               ddtdxic = termc * (ycb*zp-zcb*yp);
               ddtdyic = termc * (zcb*xp-xcb*zp);
               ddtdzic = termc * (xcb*yp-ycb*xp);
               ddtdxib = -ddtdxia - ddtdxic;
               ddtdyib = -ddtdyia - ddtdyic;
               ddtdzib = -ddtdzia - ddtdzic;

               rab2 = 2.0 / (rab*rab);
               xrab = xab * rab2;
               yrab = yab * rab2;
               zrab = zab * rab2;
               rcb2 = 2.0 / (rcb*rcb);
               xrcb = xcb * rcb2;
               yrcb = ycb * rcb2;
               zrcb = zcb * rcb2;
               rp2 = 1.0 / (rp*rp);
               xabp = (yab*zp-zab*yp) * rp2;
               yabp = (zab*xp-xab*zp) * rp2;
               zabp = (xab*yp-yab*xp) * rp2;
               xcbp = (ycb*zp-zcb*yp) * rp2;
               ycbp = (zcb*xp-xcb*zp) * rp2;
               zcbp = (xcb*yp-ycb*xp) * rp2;

               dtxiaxia = terma*(xab*xcb-dot) + ddtdxia*(xcbp-xrab);
               dtxiayia = terma*(zp+yab*xcb) + ddtdxia*(ycbp-yrab);
               dtxiazia = terma*(zab*xcb-yp) + ddtdxia*(zcbp-zrab);
               dtyiayia = terma*(yab*ycb-dot) + ddtdyia*(ycbp-yrab);
               dtyiazia = terma*(xp+zab*ycb) + ddtdyia*(zcbp-zrab);
               dtziazia = terma*(zab*zcb-dot) + ddtdzia*(zcbp-zrab);
               dtxicxic = termc*(dot-xab*xcb) - ddtdxic*(xabp+xrcb);
               dtxicyic = termc*(zp-ycb*xab) - ddtdxic*(yabp+yrcb);
               dtxiczic = -termc*(yp+zcb*xab) - ddtdxic*(zabp+zrcb);
               dtyicyic = termc*(dot-yab*ycb) - ddtdyic*(yabp+yrcb);
               dtyiczic = termc*(xp-zcb*yab) - ddtdyic*(zabp+zrcb);
               dtziczic = termc*(dot-zab*zcb) - ddtdzic*(zabp+zrcb);
               dtxiaxic = terma*(yab*yab+zab*zab) - ddtdxia*xabp;
               dtxiayic = -terma*xab*yab - ddtdxia*yabp;
               dtxiazic = -terma*xab*zab - ddtdxia*zabp;
               dtyiaxic = -terma*xab*yab - ddtdyia*xabp;
               dtyiayic = terma*(xab*xab+zab*zab) - ddtdyia*yabp;
               dtyiazic = -terma*yab*zab - ddtdyia*zabp;
               dtziaxic = -terma*xab*zab - ddtdzia*xabp;
               dtziayic = -terma*yab*zab - ddtdzia*yabp;
               dtziazic = terma*(xab*xab+yab*yab) - ddtdzia*zabp;

               dtxibxia = -dtxiaxia - dtxiaxic;
               dtxibyia = -dtxiayia - dtyiaxic;
               dtxibzia = -dtxiazia - dtziaxic;
               dtyibxia = -dtxiayia - dtxiayic;
               dtyibyia = -dtyiayia - dtyiayic;
               dtyibzia = -dtyiazia - dtziayic;
               dtzibxia = -dtxiazia - dtxiazic;
               dtzibyia = -dtyiazia - dtyiazic;
               dtzibzia = -dtziazia - dtziazic;
               dtxibxic = -dtxicxic - dtxiaxic;
               dtxibyic = -dtxicyic - dtxiayic;
               dtxibzic = -dtxiczic - dtxiazic;
               dtyibxic = -dtxicyic - dtyiaxic;
               dtyibyic = -dtyicyic - dtyiayic;
               dtyibzic = -dtyiczic - dtyiazic;
               dtzibxic = -dtxiczic - dtziaxic;
               dtzibyic = -dtyiczic - dtziayic;
               dtzibzic = -dtziczic - dtziazic;
               dtxibxib = -dtxibxia - dtxibxic;
               dtxibyib = -dtxibyia - dtxibyic;
               dtxibzib = -dtxibzia - dtxibzic;
               dtyibyib = -dtyibyia - dtyibyic;
               dtyibzib = -dtyibzia - dtyibzic;
               dtzibzib = -dtzibzia - dtzibzic;

               j = strbnd.isb[i][1];
               k = strbnd.isb[i][2];
               dr = 0.0;
               terma = 0.0;
               termc = 0.0;
               if (j > -1)
               {
                  terma = units.stbnunit* strbnd.ksb1[i] ;
                  dr = terma * (rab-bonds_ff.bl[j]);
                  terma = terma / rab;
               }
               if (k > -1)
               {
                  termc = units.stbnunit* strbnd.ksb2[i];
                  dr = dr + termc*(rcb-bonds_ff.bl[k]);
                  termc = termc / rcb;
               }
               ddrdxia = terma * xab;
               ddrdyia = terma * yab;
               ddrdzia = terma * zab;
               ddrdxic = termc * xcb;
               ddrdyic = termc * ycb;
               ddrdzic = termc * zcb;
               ddrdxib = -ddrdxia - ddrdxic;
               ddrdyib = -ddrdyia - ddrdyic;
               ddrdzib = -ddrdzia - ddrdzic;

               xab = xab / rab;
               yab = yab / rab;
               zab = zab / rab;
               xcb = xcb / rcb;
               ycb = ycb / rcb;
               zcb = zcb / rcb;

               drxiaxia = terma * (1.0-xab*xab);
               drxiayia = -terma * xab*yab;
               drxiazia = -terma * xab*zab;
               dryiayia = terma * (1.0-yab*yab);
               dryiazia = -terma * yab*zab;
               drziazia = terma * (1.0-zab*zab);
               drxicxic = termc * (1.0-xcb*xcb);
               drxicyic = -termc * xcb*ycb;
               drxiczic = -termc * xcb*zcb;
               dryicyic = termc * (1.0-ycb*ycb);
               dryiczic = -termc * ycb*zcb;
               drziczic = termc * (1.0-zcb*zcb);
               drxibxib = drxiaxia + drxicxic;
               drxibyib = drxiayia + drxicyic;
               drxibzib = drxiazia + drxiczic;
               dryibyib = dryiayia + dryicyic;
               dryibzib = dryiazia + dryiczic;
               drzibzib = drziazia + drziczic;

            if (ia == iatom)
            {
               hess.hessx[ia][0] += dt*drxiaxia + dr*dtxiaxia + 2.0*ddtdxia*ddrdxia;
               hess.hessx[ia][1] += dt*drxiayia + dr*dtxiayia + ddtdxia*ddrdyia + ddtdyia*ddrdxia;
               hess.hessx[ia][2] += dt*drxiazia + dr*dtxiazia + ddtdxia*ddrdzia + ddtdzia*ddrdxia;
               hess.hessy[ia][0] += dt*drxiayia + dr*dtxiayia + ddtdyia*ddrdxia + ddtdxia*ddrdyia;
               hess.hessy[ia][1] += dt*dryiayia + dr*dtyiayia + 2.0*ddtdyia*ddrdyia;
               hess.hessy[ia][2] += dt*dryiazia + dr*dtyiazia + ddtdyia*ddrdzia + ddtdzia*ddrdyia;
               hess.hessz[ia][0] += dt*drxiazia + dr*dtxiazia + ddtdzia*ddrdxia + ddtdxia*ddrdzia;
               hess.hessz[ia][1] += dt*dryiazia + dr*dtyiazia + ddtdzia*ddrdyia + ddtdyia*ddrdzia;
               hess.hessz[ia][2] += dt*drziazia + dr*dtziazia + 2.0*ddtdzia*ddrdzia;
               hess.hessx[ib][0] += -dt*drxiaxia + dr*dtxibxia + ddtdxia*ddrdxib + ddtdxib*ddrdxia;
               hess.hessx[ib][1] += -dt*drxiayia + dr*dtxibyia + ddtdxia*ddrdyib + ddtdyib*ddrdxia;
               hess.hessx[ib][2] += -dt*drxiazia + dr*dtxibzia + ddtdxia*ddrdzib + ddtdzib*ddrdxia;
               hess.hessy[ib][0] += -dt*drxiayia + dr*dtyibxia + ddtdyia*ddrdxib + ddtdxib*ddrdyia;
               hess.hessy[ib][1] += -dt*dryiayia + dr*dtyibyia + ddtdyia*ddrdyib + ddtdyib*ddrdyia;
               hess.hessy[ib][2] += -dt*dryiazia + dr*dtyibzia + ddtdyia*ddrdzib + ddtdzib*ddrdyia;
               hess.hessz[ib][0] += -dt*drxiazia + dr*dtzibxia + ddtdzia*ddrdxib + ddtdxib*ddrdzia;
               hess.hessz[ib][1] += -dt*dryiazia + dr*dtzibyia + ddtdzia*ddrdyib + ddtdyib*ddrdzia;
               hess.hessz[ib][2] += -dt*drziazia + dr*dtzibzia + ddtdzia*ddrdzib + ddtdzib*ddrdzia;
               hess.hessx[ic][0] += dr*dtxiaxic + ddtdxia*ddrdxic + ddtdxic*ddrdxia;
               hess.hessx[ic][1] += dr*dtxiayic + ddtdxia*ddrdyic + ddtdyic*ddrdxia;
               hess.hessx[ic][2] += dr*dtxiazic + ddtdxia*ddrdzic + ddtdzic*ddrdxia;
               hess.hessy[ic][0] += dr*dtyiaxic + ddtdyia*ddrdxic + ddtdxic*ddrdyia;
               hess.hessy[ic][1] += dr*dtyiayic + ddtdyia*ddrdyic + ddtdyic*ddrdyia;
               hess.hessy[ic][2] += dr*dtyiazic + ddtdyia*ddrdzic + ddtdzic*ddrdyia;
               hess.hessz[ic][0] += dr*dtziaxic + ddtdzia*ddrdxic + ddtdxic*ddrdzia;
               hess.hessz[ic][1] += dr*dtziayic + ddtdzia*ddrdyic + ddtdyic*ddrdzia;
               hess.hessz[ic][2] += dr*dtziazic + ddtdzia*ddrdzic + ddtdzic*ddrdzia;
            } else if (ib == iatom)
            {
               hess.hessx[ib][0] += dt*drxibxib + dr*dtxibxib + 2.0*ddtdxib*ddrdxib;
               hess.hessx[ib][1] += dt*drxibyib + dr*dtxibyib + ddtdxib*ddrdyib + ddtdyib*ddrdxib;
               hess.hessx[ib][2] += dt*drxibzib + dr*dtxibzib + ddtdxib*ddrdzib + ddtdzib*ddrdxib;
               hess.hessy[ib][0] += dt*drxibyib + dr*dtxibyib + ddtdyib*ddrdxib + ddtdxib*ddrdyib;
               hess.hessy[ib][1] += dt*dryibyib + dr*dtyibyib + 2.0*ddtdyib*ddrdyib;
               hess.hessy[ib][2] += dt*dryibzib + dr*dtyibzib + ddtdyib*ddrdzib + ddtdzib*ddrdyib;
               hess.hessz[ib][0] += dt*drxibzib + dr*dtxibzib + ddtdzib*ddrdxib + ddtdxib*ddrdzib;
               hess.hessz[ib][1] += dt*dryibzib + dr*dtyibzib + ddtdzib*ddrdyib + ddtdyib*ddrdzib;
               hess.hessz[ib][2] += dt*drzibzib + dr*dtzibzib + 2.0*ddtdzib*ddrdzib;
               hess.hessx[ia][0] += -dt*drxiaxia + dr*dtxibxia + ddtdxib*ddrdxia + ddtdxia*ddrdxib;
               hess.hessx[ia][1] += -dt*drxiayia + dr*dtxibyia + ddtdxib*ddrdyia + ddtdyia*ddrdxib;
               hess.hessx[ia][2] += -dt*drxiazia + dr*dtxibzia + ddtdxib*ddrdzia + ddtdzia*ddrdxib;
               hess.hessy[ia][0] += -dt*drxiayia + dr*dtyibxia + ddtdyib*ddrdxia + ddtdxia*ddrdyib;
               hess.hessy[ia][1] += -dt*dryiayia + dr*dtyibyia + ddtdyib*ddrdyia + ddtdyia*ddrdyib;
               hess.hessy[ia][2] += -dt*dryiazia + dr*dtyibzia + ddtdyib*ddrdzia + ddtdzia*ddrdyib;
               hess.hessz[ia][0] += -dt*drxiazia + dr*dtzibxia + ddtdzib*ddrdxia + ddtdxia*ddrdzib;
               hess.hessz[ia][1] += -dt*dryiazia + dr*dtzibyia + ddtdzib*ddrdyia + ddtdyia*ddrdzib;
               hess.hessz[ia][2] += -dt*drziazia + dr*dtzibzia + ddtdzib*ddrdzia + ddtdzia*ddrdzib;
               hess.hessx[ic][0] += -dt*drxicxic + dr*dtxibxic + ddtdxib*ddrdxic + ddtdxic*ddrdxib;
               hess.hessx[ic][1] += -dt*drxicyic + dr*dtxibyic + ddtdxib*ddrdyic + ddtdyic*ddrdxib;
               hess.hessx[ic][2] += -dt*drxiczic + dr*dtxibzic + ddtdxib*ddrdzic + ddtdzic*ddrdxib;
               hess.hessy[ic][0] += -dt*drxicyic + dr*dtyibxic + ddtdyib*ddrdxic + ddtdxic*ddrdyib;
               hess.hessy[ic][1] += -dt*dryicyic + dr*dtyibyic + ddtdyib*ddrdyic + ddtdyic*ddrdyib;
               hess.hessy[ic][2] += -dt*dryiczic + dr*dtyibzic + ddtdyib*ddrdzic + ddtdzic*ddrdyib;
               hess.hessz[ic][0] += -dt*drxiczic + dr*dtzibxic + ddtdzib*ddrdxic + ddtdxic*ddrdzib;
               hess.hessz[ic][1] += -dt*dryiczic + dr*dtzibyic + ddtdzib*ddrdyic + ddtdyic*ddrdzib;
               hess.hessz[ic][2] += -dt*drziczic + dr*dtzibzic + ddtdzib*ddrdzic + ddtdzic*ddrdzib;
            }else if (ic == iatom)
            {
               hess.hessx[ic][0] += dt*drxicxic + dr*dtxicxic + 2.0*ddtdxic*ddrdxic;
               hess.hessx[ic][1] += dt*drxicyic + dr*dtxicyic + ddtdxic*ddrdyic + ddtdyic*ddrdxic;
               hess.hessx[ic][2] += dt*drxiczic + dr*dtxiczic + ddtdxic*ddrdzic + ddtdzic*ddrdxic;
               hess.hessy[ic][0] += dt*drxicyic + dr*dtxicyic + ddtdyic*ddrdxic + ddtdxic*ddrdyic;
               hess.hessy[ic][1] += dt*dryicyic + dr*dtyicyic + 2.0*ddtdyic*ddrdyic;
               hess.hessy[ic][2] += dt*dryiczic + dr*dtyiczic + ddtdyic*ddrdzic + ddtdzic*ddrdyic;
               hess.hessz[ic][0] += dt*drxiczic + dr*dtxiczic + ddtdzic*ddrdxic + ddtdxic*ddrdzic;
               hess.hessz[ic][1] += dt*dryiczic + dr*dtyiczic + ddtdzic*ddrdyic + ddtdyic*ddrdzic;
               hess.hessz[ic][2] += dt*drziczic + dr*dtziczic + 2.0*ddtdzic*ddrdzic;
               hess.hessx[ib][0] += -dt*drxicxic + dr*dtxibxic + ddtdxic*ddrdxib + ddtdxib*ddrdxic;
               hess.hessx[ib][1] += -dt*drxicyic + dr*dtxibyic + ddtdxic*ddrdyib + ddtdyib*ddrdxic;
               hess.hessx[ib][2] += -dt*drxiczic + dr*dtxibzic + ddtdxic*ddrdzib + ddtdzib*ddrdxic;
               hess.hessy[ib][0] += -dt*drxicyic + dr*dtyibxic + ddtdyic*ddrdxib + ddtdxib*ddrdyic;
               hess.hessy[ib][1] += -dt*dryicyic + dr*dtyibyic + ddtdyic*ddrdyib + ddtdyib*ddrdyic;
               hess.hessy[ib][2] += -dt*dryiczic + dr*dtyibzic + ddtdyic*ddrdzib + ddtdzib*ddrdyic;
               hess.hessz[ib][0] += -dt*drxiczic + dr*dtzibxic + ddtdzic*ddrdxib + ddtdxib*ddrdzic;
               hess.hessz[ib][1] += -dt*dryiczic + dr*dtzibyic + ddtdzic*ddrdyib + ddtdyib*ddrdzic;
               hess.hessz[ib][2] += -dt*drziczic + dr*dtzibzic + ddtdzic*ddrdzib + ddtdzib*ddrdzic;
               hess.hessx[ia][0] += dr*dtxiaxic + ddtdxic*ddrdxia + ddtdxia*ddrdxic;
               hess.hessx[ia][1] += dr*dtyiaxic + ddtdxic*ddrdyia + ddtdyia*ddrdxic;
               hess.hessx[ia][2] += dr*dtziaxic + ddtdxic*ddrdzia + ddtdzia*ddrdxic;
               hess.hessy[ia][0] += dr*dtxiayic + ddtdyic*ddrdxia + ddtdxia*ddrdyic;
               hess.hessy[ia][1] += dr*dtyiayic + ddtdyic*ddrdyia + ddtdyia*ddrdyic;
               hess.hessy[ia][2] += dr*dtziayic + ddtdyic*ddrdzia + ddtdzia*ddrdyic;
               hess.hessz[ia][0] += dr*dtxiazic + ddtdzic*ddrdxia + ddtdxia*ddrdzic;
               hess.hessz[ia][1] += dr*dtyiazic + ddtdzic*ddrdyia + ddtdyia*ddrdzic;
               hess.hessz[ia][2] += dr*dtziazic + ddtdzic*ddrdzia + ddtdzia*ddrdzic;
            }
            }
          }
      }
}
                  
