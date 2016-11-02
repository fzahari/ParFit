#include "pcwin.h"
#include "pcmod.h"
#include "field.h"
#include "atom_k.h"
#include "pcmwin1.h"

void set_atomtypes(int newtype)
{
    int i,badtype;
    char hatext[256];
    
    if (newtype == MMX)
    {
        for (i=1; i <= natom; i++)
        {
            atom[i].type = atom[i].mmx_type;
            atom[i].tclass = atom_k.tclass[atom[i].type];
        }   
    } else if (newtype == MM3)
    {
        for (i=1; i <= natom; i++)
        {
            atom[i].type = atom[i].mm3_type;
            atom[i].tclass = atom_k.tclass[atom[i].type];
        }
    } else if (newtype == MMFF94)
    {
        for (i=1; i <= natom; i++)
        {
            atom[i].type = atom[i].mmff_type;
            atom[i].tclass = atom_k.tclass[atom[i].type];
        }
    } else if (newtype == OPLSAA)
    {
        for (i=1; i <= natom; i++)
        {
            atom[i].type = atom[i].opls_type;
            atom[i].tclass = atom_k.tclass[atom[i].type];
        }
    } else if (newtype == AMBER)
    {
        for (i=1; i <= natom; i++)
        {
            atom[i].type = atom[i].amber_type;
            atom[i].tclass = atom_k.tclass[atom[i].type];
        }
    }
    badtype = FALSE;
    
    for (i=1; i <= natom; i++)
         if (atom[i].type == 0)
             badtype = TRUE;
     if (badtype == TRUE)
     {
           sprintf(hatext,"Error setting atom types in FField: %d \n",newtype);
           message_alert(hatext,"ERROR");
     }
}
void type()
{
    int i, j, jjbo, icount,ncarbon, nhyd;

    icount = 0;
    ncarbon = 0;
    nhyd = 0;
    for (i=1; i <= natom; i++)
    {
        if (atom[i].atomnum == 6)
        {
            jjbo = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                   jjbo += atom[i].bo[j];
                if (atom[atom[i].iat[j]].atomnum == 1)
                   nhyd++;
            }
            if ( jjbo < 3 && atom[i].mmx_type == 40 && icount > 0  )
            {
                quick_type();
                set_atomtypes(field.type);
                return;
            } else if (jjbo < 4 && atom[i].mmx_type != 40 && icount > 0)
            {
                quick_type();
                set_atomtypes(field.type);
                return;
            }else
            {
                icount++;
                ncarbon++;
            }
            if (icount > 15)
               break; 
        } else if (atom[i].atomnum == 8)
        {
            jjbo = 0;
            for (j=0; j < MAXIAT; j++)
            {
                if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
                   jjbo += atom[i].bo[j];
            }
            if (jjbo < 1 && icount > 0)
            {
                quick_type();
                set_atomtypes(field.type);
                return;
            } else if (jjbo == 1 && icount > 0 && atom[i].mmx_type == 6)
            {
                quick_type();
                set_atomtypes(field.type);
                return;                
            }else
                icount++;
            if (icount > 15)
               break; 
        }
        
    }
    if (ncarbon > 0 && nhyd == 0)  // pathological cases where there are no carbons, or no carbons that normally have hydrogens 
    {
        type_mmx();
        set_atomtypes(field.type);
        return;
    }
// normal typing routines   
    if (field.type == MMX)
    {
        type_mmx();
        set_atomtypes(MMX);
    } else if (field.type == MM3)
    {
        for (i=1; i <= natom; i++)
        {
            if (atom[i].mmx_type == 20)
            {
                hdel(1);
                break;
            }
        }
        type_mmx();
        set_atomtypes(MM3);
    } else if (field.type == MMFF94)
    {
        for (i=1; i <= natom; i++)
        {
            if (atom[i].mmx_type == 20)
            {
                hdel(1);
                break;
            }
        }
        type_mmx();
        set_atomtypes(MMFF94);
    } else if (field.type == AMBER)
    {
        for (i=1; i <= natom; i++)
        {
            if (atom[i].mmx_type == 20)
            {
                hdel(1);
                break;
            }
        }
        type_mmx();
        set_atomtypes(AMBER);
    } else if (field.type == OPLSAA)
    {
        for (i=1; i <= natom; i++)
        {
            if (atom[i].mmx_type == 20)
            {
                hdel(1);
                break;
            }
        }
        type_mmx();
        set_atomtypes(OPLSAA);
    }
    else
       message_alert("No typing rules for this force field are implemented!","ERROR");  
} 
/* ------------------------ */

double dihdrl(int i1,int i2,int i3,int i4)
{
        float           a1, a2, ang, b1, b2, c1, c2, siign,
                        x1, x3, x4, y1, y3, y4, z1, z3, z4;
        double          dihdrl_v,r1,r2;

        /****the arguments of this function are the atom indices of atoms a-b-c-d
         *     which determine dihedral angle dihdrl.
         *     dihdrl returns the degree value of the dihedral('torsional')
         *     angle defined by atoms i1,i2,i3, &i4.
         */
         
        dihdrl_v = 0.0;
        if (i1 <= 0 || i2 <= 0 || i3 <= 0 || i4 <= 0)
           return (dihdrl_v);
        if (i1 > natom || i2 > natom || i3 > natom || i4 > natom)
           return (dihdrl_v);        
           
        a1 = atom[i2].x;
        b1 = atom[i2].y;
        c1 = atom[i2].z;
        x1 = atom[i1].x - a1;
        y1 = atom[i1].y - b1;
        z1 = atom[i1].z - c1;
        x3 = atom[i3].x - a1;
        y3 = atom[i3].y - b1;
        z3 = atom[i3].z - c1;
        x4 = atom[i4].x - a1;
        y4 = atom[i4].y - b1;
        z4 = atom[i4].z - c1;
        a1 = y1 * z3 - y3 * z1;
        b1 = x3 * z1 - x1 * z3;
        c1 = x1 * y3 - x3 * y1;
        a2 = y4 * z3 - y3 * z4;
        b2 = x3 * z4 - x4 * z3;
        c2 = x4 * y3 - x3 * y4;
        r1 = sqrt(a1 * a1 + b1 * b1 + c1 * c1);
        r2 = sqrt(a2 * a2 + b2 * b2 + c2 * c2);
        if (r1 != 0.0 && r2 != 0.0)
           ang = (a1 * a2 + b1 * b2 + c1 * c2) /(r1*r2);
        else
           ang = 0.0;

        if (fabs(ang) > 1)
                ang = fabs(ang) / ang;
        ang = acos(ang);
        dihdrl_v = 57.29578F * ang;
        siign = x1 * a2 + y1 * b2 + z1 * c2;
        if (siign < 0.0F)
                dihdrl_v = -dihdrl_v;
        return (dihdrl_v);
}                               /* end of function */
/* --------------------------------- */
void avgleg()
{
        int i, iatyp, ib1, ib2, j;
        int hasSulf;
        float an, xx, xyzdis, yy, zz;

        /****  if there is no bond, average=1.0 *** */

        an = 0.0F;
        xyzdis = 0.0F;
        hasSulf = FALSE;
        for( i = 1; i <= natom; i++ )
        {
                if( atom[i].type && atom[i].type != 20 )
                {
                   for( j = 0; j < MAXIAT; j++ )
                   {
                      if( atom[i].iat[j] && atom[i].bo[j] != 9 )
                      {
        /*             count bonds only once - don't count bonds to hydrogen */
                          if( atom[i].iat[j] >= i )
                          {
                              ib1 = i;
                              iatyp = atom[i].type;
                              if (atom[i].atomnum == 16) hasSulf = TRUE;
                              if(  iatyp != 5 && iatyp != 21 && iatyp != 23 && iatyp != 24 && iatyp != 28 && iatyp != 36 && iatyp != 15 && iatyp != 38 )
                              {
                                   ib2 = atom[i].iat[j];
                                   iatyp = atom[atom[i].iat[j]].type;
                                   if( iatyp != 20 )
                                   {
                                      if( iatyp != 5 && iatyp != 21 && iatyp != 23 && iatyp != 24 && iatyp != 28 && iatyp != 36 && iatyp != 15 && iatyp != 38)
                                      {
                                            an += 1.0F;
                                            xx = atom[ib1].x - atom[ib2].x;
                                            yy = atom[ib1].y - atom[ib2].y;
                                            zz = atom[ib1].z - atom[ib2].z;
                                            xyzdis += sqrt( xx*xx + yy*yy + zz*zz );
                                       }
                                    }
                                }
                           }
                       }
                   }
                }
        }
        /* if an = 0 then there are no heavy atom bonds
         *  need to redo calculations based on c-h bonds */
        if( an != 0.0F )
        {
            if (hasSulf)
                flags.avleg = 1.73F/(xyzdis/an);        /* assume normal c-c distance   */
            else
                flags.avleg = 1.53F/(xyzdis/an);        /* assume normal c-c distance   */
        }
        else if( natom > 1 )
        {
                xyzdis = 0.0F;
                an = 0.0F;
                for( i = 1; i <= natom; i++ )
                {
                        ib1 = i;
                        for( j = 0; j < MAXIAT; j++ )
                        {
                                if( atom[i].iat[j] && atom[i].bo[j] != 9 )
                                {
                                        ib2 = atom[i].iat[j];
                                        if( ib2 > ib1 )
                                        {
                                                an += 1.0F;
                                                xx = atom[ib1].x - atom[ib2].x;
                                                yy = atom[ib1].y - atom[ib2].y;
                                                zz = atom[ib1].z - atom[ib2].z;
                                                xyzdis += sqrt( xx*xx + yy*yy + zz*zz );
                                        }
                                }
                        }
                }
                flags.avleg = 1.10F/(xyzdis/an); /* no c-c bonds assume c-h bond dist */
        }
        else
        {
                flags.avleg = 1.20F; /*  default when all else goes wrong */
        }
        return;
}
/* -------------------------- */
void center(int istruct)
{
        //  istruct = -2 do all
        //  istruct = 1 do main structrure 
        //  istruct = number do that substructure

        float centerx, centery, centerz;
        int i;

        centerx = 0.0F;
        centery = 0.0F;
        centerz = 0.0F;

        if (istruct == -2)
        {
                for (i=1; i <= natom; i++)
                {
              centerx += atom[i].x/natom;
          centery += atom[i].y/natom;
                  centerz += atom[i].z/natom;
                }
                for (i=1; i <= natom; i++)
                {
                        atom[i].x -= centerx;
                        atom[i].y -= centery;
                        atom[i].z -= centerz;
                }
        } else if (istruct == -1)      // main structure only
        {
                for (i=1; i <= natom; i++)
                {
                        if(atom[i].substr[0] == 0)
                        {
                      centerx += atom[i].x/natom;
                  centery += atom[i].y/natom;
                          centerz += atom[i].z/natom;
                        }
                }
                for (i=1; i <= natom; i++)
                {
                        if(atom[i].substr[0] == 0)
                        {
                                atom[i].x -= centerx;
                                atom[i].y -= centery;
                                atom[i].z -= centerz;
                        }
                }
        } else                  // substructure
        {
                for (i=1; i <= natom; i++)
                {
                        if (atom[i].substr[0] & (1L << istruct) )
                        {
                      centerx += atom[i].x/natom;
                  centery += atom[i].y/natom;
                          centerz += atom[i].z/natom;
                        }
                }
                for (i=1; i <= natom; i++)
                {
                        if (atom[i].substr[0] & (1L << istruct) )
                        {
                                atom[i].x -= centerx;
                                atom[i].y -= centery;
                                atom[i].z -= centerz;
                        }
                }
        }
}
/* ------------------------ */

double ooplane(int ia,int ib,int ic,int id)
{
      double angle,dot,cosine;
      double xia,yia,zia;
      double xib,yib,zib;
      double xic,yic,zic;
      double xid,yid,zid;
      double xab,yab,zab,rab2,rab;
      double xbc,ybc,zbc,rbc2,rbc;
      double xbd,ybd,zbd,rbd2,rbd;
      double sine;

        /****the arguments of this function are the atom indices of atoms a-b-c-d
         *     which determine dihedral angle dihdrl.
         *     dihdrl returns the degree value of the dihedral('torsional')
         *     angle defined by atoms i1,i2,i3, &i4.
         */
      angle = 0.0;
                 xia = atom[ia].x;
                 yia = atom[ia].y;
                 zia = atom[ia].z;
                 xib = atom[ib].x;
                 yib = atom[ib].y;
                 zib = atom[ib].z;
                 xic = atom[ic].x;
                 yic = atom[ic].y;
                 zic = atom[ic].z;
                 xid = atom[id].x;
                 yid = atom[id].y;
                 zid = atom[id].z;
                 
                 xab = xia - xib;
                 yab = yia - yib;
                 zab = zia - zib;
                 rab2 = (xab*xab+yab*yab+zab*zab);
                 rab = sqrt(rab2);
                 
                 xbc = xic - xib;
                 ybc = yic - yib;
                 zbc = zic - zib;
                 rbc2 = xbc*xbc+ybc*ybc+zbc*zbc;
                 rbc = sqrt(rbc2);
                 
                 if (rab2 != 0.0 && rbc2 != 0.0)
                 {
                    dot = xab*xbc + yab*ybc + zab*zbc;
                    cosine = dot /(rab*rbc);
                    if (cosine < -1.0)
                       cosine = -1.0;
                    if (cosine > 1.0)
                       cosine = 1.0;
                    sine = sqrt(1.00-cosine*cosine);

                    xbd = atom[id].x - atom[ib].x;
                    ybd = atom[id].y - atom[ib].y;
                    zbd = atom[id].z - atom[ib].z;
                    rbd2 = xbd*xbd+ybd*ybd+zbd*zbd;
                    rbd = sqrt(rbd2);
                    
                    cosine = ( xbd*(yab*zbc-zab*ybc) + ybd*(zab*xbc-xab*zbc) + zbd*(xab*ybc-yab*xbc))/((rab*rbc*rbd)*sine);

                    angle = radian * asin(cosine);
              }
        return (angle);
} 
/* ------------------------------------------------------------- */       
void quick_type()
{
int    i, iatik, iatkkk, iatomt, iatyp, ibt[MAXIAT], ibtk,
           ibut, idc, itia, j, ja, jji, k, k1, k2, katm, ki,
           kk, knt, l, lat, latm, lnt, k3, matm, jatm;

        /* *** loop over all atoms and determine type *** */
        for (i = 1; i <= natom; i++)
        {
            if (atom[i].mmx_type != 0)
            {
                iatyp = atom[i].mmx_type;
                if (atom[i].mmx_type == 22 || atom[i].mmx_type == 56 || atom[i].mmx_type == 57)
                        iatyp = 1;
                iatomt = iatyp;
                /* *** load all attached bond types into ibt() *** except for
                 * coordinated bonds */
                for (j = 0; j < MAXIAT; j++) 
                {
                        if (atom[i].bo[j] == 9 || atom[i].iat[j] == 0) 
                            ibt[j] = 0;
                        else 
                            ibt[j] = atom[i].bo[j];
                 }
                /* *** now branch on the atom type *** */
                if (iatyp < 5)
                        goto L_20;
                if (iatyp == 5 || iatyp == 21 || iatyp == 23 || iatyp == 24 || iatyp== 70)
                        goto L_100;
                if (iatyp == 6 || iatyp == 7)
                        goto L_140;
                if (((iatyp == 8 || iatyp == 9) || iatyp == 10) || iatyp == 37)
                        goto L_170;
                if ((((iatyp == 16 || iatyp == 17) || iatyp == 18) || iatyp == 15) || iatyp == 38)
                        goto L_230;
                if (iatyp == 34 || iatyp == 39)
                        goto L_250;
                if (iatyp == 25 || iatyp == 47)
                        goto L_270;
                if (iatyp == 44 || iatyp == 58)  // Aluminum
                        goto L_271;
                if (iatyp == 12 || iatyp == 74)  // chlorine
                        goto L_272;
                /* *** all other atoms have only one type *** */
                goto L_290;

                /* *** carbon *** */
L_20:
                iatomt = 1;
                for (k = 0; k < MAXIAT; k++)
                {
                    if (atom[i].iat[k] != 0 && atom[i].bo[k] >= 2 && atom[i].bo[k] != 9)
                    {
                        /* *** sp carbon *** */
                        if (atom[i].bo[k] == 3)
                        {
                            iatomt = 4;
                            goto L_290;
                        }
                        /* *** sp2 carbon *** */
                        ja = atom[i].iat[k];
                        iatomt = 2;
                        if (atom[ja].atomnum == 8 )
                        {
                                iatomt = 3;
                                /* check for ketenes and co2 */
                                kk = k + 1;
                                for (ki = kk; ki <= 10; ki++)
                                {
                                        if (atom[i].bo[ki] == 2)
                                        {
                                                iatomt = 4;
                                                goto L_290;
                                        }
                                }
                                goto L_290;
                        }
                        /* check for allene */
                        kk = k + 1;
                        for (ki = kk; ki <= 8; ki++)
                        {
                                if (ibt[ki] == 2)
                                {
                                        iatomt = 4;
                                        goto L_290;
                                }
                        }
                    }
                }
                /* *** 3 and 4 -mem ring carbons *** */
                if (iatomt != 1)
                        goto L_290;
                for (k = 0; k < MAXIAT; k++) 
                {
                        if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
                        {       
                                katm = atom[i].iat[k];
                                for (k1 = 0; k1 < MAXIAT; k1++) 
                                {
                                        if (atom[katm].iat[k1] != 0 && atom[katm].iat[k1] != i && atom[katm].bo[k1] != 9)
                                        {
                                                latm = atom[katm].iat[k1];
                                                for (k2 = 0; k2 < MAXIAT; k2++) 
                                                {
                                                        if (atom[latm].iat[k2] != 0 && atom[latm].iat[k2] != katm && atom[latm].bo[k2] != 9)
                                                        {
                                                                if (atom[latm].iat[k2] == i) 
                                                                {
                                                                        iatomt = 22;
                                                                        ibut = 1;
                                                                        goto L_290;
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
//  now check for four membered rings
                for (k = 0; k < MAXIAT; k++) 
                {
                        if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
                        {       
                                katm = atom[i].iat[k];
                                for (k1 = 0; k1 < MAXIAT; k1++) 
                                {
                                        if (atom[katm].iat[k1] != 0 && atom[katm].iat[k1] != i && atom[katm].bo[k1] != 9)
                                        {
                                                latm = atom[katm].iat[k1];
                                                for (k2 = 0; k2 < MAXIAT; k2++) 
                                                {
                                                        if (atom[latm].iat[k2] != 0 && atom[latm].iat[k2] != katm && atom[latm].iat[k2] != i  && atom[latm].bo[k2] != 9)
                                                        {
                                                                matm = atom[latm].iat[k2];
                                                                for (k3 = 0; k3 < MAXIAT; k3++)
                                                                {
                                                                        if (atom[matm].iat[k3] != 0 && atom[matm].iat[k3] != latm && atom[matm].bo[k3] != 9)
                                                                        {
                                                                                if (atom[matm].iat[k3] == i) 
                                                                                {
                                                                                        iatomt = 56;
                                                                                        goto L_290;
                                                                                }
                                                                        }
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
                goto L_290;

                /* *** hydrogen *** */
L_100:
                iatomt = 5;
                knt = 0;
                for (k=0;k < MAXIAT; k++)
                  if (atom[i].iat[k] != 0)
                     knt++;
                if (knt == 2)
                {
                   iatomt = 70;
                   goto L_290;
                }
                for (k = 0; k < MAXIAT; k++) 
                {
                   if (atom[i].iat[k] != 0) 
                   {
                       ja = atom[atom[i].iat[k]].mmx_type;
                       goto L_120;
                   }
                }
                goto L_290;
L_120:
                if (((ja == 8 || ja == 9) || ja == 37) || ja == 55) {
                        iatomt = 23;
                } else if (ja == 41 || ja == 46) {
                        iatomt = 24;
                } else if (ja == 6 || ja == 53) {
                        iatomt = 21;
                        for (j = 0; j < MAXIAT; j++) {
                                if (atom[atom[i].iat[k]].iat[j]
                                                        == 0)
                                        goto L_130;
                                if (atom[atom[i].iat[k]].iat[j]
                                                        == i)
                                        goto L_130;
                                if (atom[atom[atom[i].iat[k]].iat[j]].mmx_type == 3) {
                                        iatomt = 24;
                                } else if (atom[atom[atom[i].iat[k]].iat[j]].mmx_type == 2 ||
                                           atom[atom[atom[i].iat[k]].iat[j]].mmx_type == 40) {
                                        iatomt = 28;
                                }
                L_130:
                                ;
                        }
                }
                goto L_290;

                /* *** oxygen *** */
L_140:
                iatomt = 6;
                for (k = 0; k < MAXIAT; k++)
                {
                    if (atom[i].iat[k] != 0)
                    {
                        if (ibt[k] == 2)
                        {
                                iatomt = 7;
                                goto L_290;
                        } else if (ibt[k] == 3)
                        {
                            for (lat = 0; lat < MAXIAT; lat++)
                            {
                                 if (atom[i].iat[lat] != 0)
                                 {
                                     if (atom[atom[i].iat[lat]].mmx_type == 4)
                                     {
                                         iatomt = 46;
                                         goto L_290;
                                     }
                                 }
                            }
                        }
                    }
                }
                goto L_290;

                /* *** nitrogen *** */
L_170:
                iatomt = 8;
                ibtk = 0;
                /* determine maximum number of ligands attached to nitrogen */
                jji = 0;
                for (k = 0; k < MAXIAT; k++) 
                {
                        if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9)
                                jji += 1;
                }
                for (k = 0; k < MAXIAT; k++) 
                {
                        if (ibt[k] != 9) 
                        {
                             ibtk += ibt[k];
                        }
                        if (atom[i].iat[k] != 0) 
                        {
                                if (atom[atom[i].iat[k]].mmx_type== 20 ||
                                   (atom[atom[i].iat[k]].mmx_type >= 300) )
                                        goto L_190;
                        }
                }
                /* ammonium nitrogen - no lone pairs */
                if (ibtk == 4) 
                {
                        iatomt = 41;
                        goto L_290;
                }
                /* got a lone pair can not be ammonium */
L_190:
                for (k = 0; k < MAXIAT; k++) 
                {
                        /* azo, imine */
                        if (ibt[k] == 2) 
                        {
                                iatomt = 37;
                                goto L_290;
                                /* nitrile */
                        } else if (ibt[k] == 3) 
                        {
                                iatomt = 10;
                                goto L_290;
                        }
                }
                /* count the number of double bonds to an adjacent Sulfur or
                 * Selenium */
                idc = 0;
                for (k = 0; k < MAXIAT; k++) 
                {
                     iatik = atom[i].iat[k];
                     if (iatik != 0 && atom[i].bo[k] != 9) 
                     {
                         if (atom[iatik].atomnum == 16) 
                         {
                            for (kk = 0; kk < MAXIAT; kk++) 
                            {
                                 iatkkk = atom[iatik].iat[kk];
                                 if (iatkkk != 0 && atom[iatik].bo[kk] != 9) 
                                 {
                                     if (atom[iatik].bo[kk] > 1)
                                          idc += 1;
                                 }
                            }
                         }
                     }
                }
                if (idc == 1) 
                {
                    iatomt = 9;
                } else if (idc == 2) 
                {
                        iatomt = 8;
                }
                /* *** loop over all adjacent atoms to find enamines, amides */
                for (k = 0; k < MAXIAT; k++) 
                {
                   if (atom[i].iat[k] != 0 && atom[i].bo[k] != 9 )
                   {
                        itia = atom[atom[i].iat[k]].mmx_type;
                        if ( ((itia >= 300) || itia == 29 || itia == 30 || itia == 40  || itia == 26) && jji <= 3)
                        {
                            iatomt = 9;
                            goto L_290;
                        }
                        if (atom[atom[i].iat[k]].atomnum == 6) // adjacent carbon
                        {
                            jatm = atom[i].iat[k];
                            for (j=0; j < MAXIAT; j++)
                            {
                                if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] >= 2)
                                {
                                    iatomt = 9;
                                    goto L_290;
                                }
                            }
                        }
                        if (atom[atom[i].iat[k]].atomnum == 7) // adjacent nitrogen
                        {
                            jatm = atom[i].iat[k];
                            for (j=0; j < MAXIAT; j++)
                            {
                                if (atom[jatm].iat[j] != 0 && atom[jatm].bo[j] >= 2)
                                {
                                    iatomt = 9;
                                    goto L_290;
                                } else if (atom[jatm].iat[j] != 0 && atom[jatm].iat[j] != i)
                                {
                                    katm = atom[jatm].iat[j];
                                    if (atom[katm].atomnum == 6)
                                    {
                                        for (l=0; l < MAXIAT; l++)
                                        {
                                            if (atom[katm].iat[l] != 0 && atom[katm].iat[l] != jatm && atom[katm].bo[l] == 2)
                                            {
                                                iatomt = 9;
                                                goto L_290;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                   }
                }
                goto L_290;

                /* *** sulfur *** */
L_230:
                iatomt = 15;
                if (iatyp == 16)
                {
                        iatomt = 16;
                        goto L_290;
                }
                knt = 0;
                lnt = 0;
                for (k = 0; k < MAXIAT; k++) {
                        if (atom[i].iat[k] != 0) {
                                if (ibt[k] > 0 && atom[atom[i].iat[k]].mmx_type
                                                                 != 20) {
                                        lnt += 1;
                                }
                                if (ibt[k] > 1) {
                                        knt += 1;
                                        if (atom[atom[i].iat[k]].mmx_type
                                                          <= 4) {
                                                iatomt = 38;
                                                goto L_290;
                                        }
                                }
                        }
                }
                if (knt == 1) {
                        iatomt = 17;
                        if (lnt == 1)
                                iatomt = 38;
                }
                if (knt >= 2)
                        iatomt = 18;
                goto L_290;
                /* ** selenium */
L_250:
                iatomt = 34;
                for (k = 0; k < MAXIAT; k++) {
                        if (atom[i].iat[k] != 0) {
                                if (ibt[k] == 2) {
                                        iatomt = 39;
                                        goto L_290;
                                }
                        }
                }
                goto L_290;
                /* ** phosphorous */
L_270:
                iatomt = 25;
                knt = 0;
                for (j = 0; j < MAXIAT; j++) {
                        if (atom[i].iat[j] != 0)
                                knt += 1;
                }
                if (knt >= 5)
                        iatomt = 47;
                goto L_290;
                /*  aluminum  */
L_271:
               iatomt = 44;
               knt = 0;
                for (j = 0; j < MAXIAT; j++)
                {
                   if (atom[i].iat[j] != 0)
                      knt += 1;
                }
                if (knt == 4)
                   iatomt = 58;
                goto L_290;
           /*  Chlorine  */
L_272:
               iatomt = 12;
               knt = 0;
                for (j = 0; j < MAXIAT; j++)
                {
                   if (atom[i].iat[j] != 0)
                      knt += 1;
                }
                if (knt == 2)
                   iatomt = 74;
                goto L_290;
                /* *** typing is now done, load value and return *** */
L_290:
                if (field.type == MMX && iatomt < 300 )
                   set_atomdata(i,iatomt,0,0,0,0);
            }
        }
        return;
}
                 
