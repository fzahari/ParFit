#include "pcwin.h"
#include "pcmod.h"
#include "utility.h"
#include "energies.h"
#include "torsions.h"
#include "gmmx.h"
#include "pot.h"
#include "field.h"
#include "optimize.h"
#include "minim_values.h"
#include "dipmom.h"
#include "files.h"
#include "pcmfile.h"
#include "gmmx_run2.h"

static int *jjii;

int get_filename()
{
    int i;
    char filename[80];

    printf("Enter filename\n");
    gets(filename);
    strcpy(Openbox.fname,filename);
        strcpy(gmmx_data.foname,Openbox.fname);

        strcpy(gmmx_data.jobname, Openbox.fname);
        strcpy(Savebox.path,Openbox.path);
        for (i=0; i < strlen(gmmx_data.jobname); i++)
        {
            if ( gmmx_data.jobname[i] == '.')
            {
                gmmx_data.jobname[i] = '\0';
                break;
            }
        }
        strcpy(gmmx.tmpname, gmmx_data.jobname);
        strcat(gmmx.tmpname, "1.tmp");
    
        strcpy(gmmx.rmsname, gmmx_data.jobname);
        strcat(gmmx.rmsname, "1.rms");
    
        strcpy(gmmx.pkmname, gmmx_data.jobname);
        strcat(gmmx.pkmname, "2.pkm");
    
        strcpy(gmmx.rstname, gmmx_data.jobname);
        strcat(gmmx.rstname, "1.rst");
        
        strcpy(gmmx_data.foname2,gmmx_data.jobname);
        strcat(gmmx_data.foname2, "2.pcm");
    return (TRUE);
}

void gmmx2()
{
    int iftype,i, nrc, j,jj;
    double etot;

     gmmx.run = TRUE;
     gmmx.job_done = TRUE;
     gmmx_data.kmin = MAXFND;
     gmmx_data.kdup = MAXFND;
     gmmx_data.ecut = FALSE;
     gmmx_data.comp_method = gmmx_data.comp_method2;
     iftype = Openbox.ftype = FTYPE_PCM;
     files.nfiles = 0;
     check_numfile( iftype);

     gmmx_abort = FALSE;
/*  do calculations */
     strcpy(Openbox.fname,gmmx_data.foname);
     strcpy(Savebox.fname,gmmx.tmpname);
     initialize();
     pcmfin(1,0);
         
//  count number of attachments to each atom
    jjii = ivector(1,natom+1);
    for (i=1; i <= natom; i++)
    {
        jj = 0;
        for (j=0; j < MAXIAT; j++)
          if (atom[i].iat[j] != 0 && atom[i].bo[j] != 9)
            jj++;
        jjii[i] = jj;
    }

     if(gmmx_data.nopi)
       pimarkselection();  
     else
       pireset();    
                   
//  check chirality
     hster(0, jjii);
//  check epimers
     nrc = lepimer(0);
     nrc = bstereo(0);
//  check double bond stereochemistry

     files.batch = True;
     gmmx.nfound = 0;
     gmmx.neminim = 0;
     gmmx.eminim = 1000;
     gmmx.nsaved = 0;
     gmmx.nsflow = 1;
     gmmx.maxdupl = 0;
     
     rdfile(0);
    nrc = setup_calculation();
    if (nrc == FALSE)
    {
            message_alert("Error setting up GMMX Job","Error");
            return;
    }
     
     for (i = 1; i <= files.nfiles; i++)
     {
        rdfile(i);
        etot = energy();

// Hay debug
//      if (field.type == MMX || field.type == MM3)
//         eheat();

      optimize_data.initial_energy = energies.total;
      optimize_data.initial_heat = get_heat_of_formation();
      if (minim_values.iprint)
      {
         fprintf(pcmoutfile,"\nInitial Energy : %f\n",energies.total);
         fprintf(pcmoutfile,"Str: %f    Bnd: %f    Tor: %f\n", energies.estr, energies.ebend, energies.etor);
         fprintf(pcmoutfile,"StrBnd: %f VDW: %f     QQ: %f\n",energies.estrbnd,energies.evdw+energies.e14+energies.ehbond, energies.eu);
         fprintf(pcmoutfile,"OOP: %f AA: %f     Strtor: %f\n",energies.eopb,energies.eangang, energies.estrtor);
         fprintf(pcmoutfile,"Urey: %f\n",energies.eurey);
         if (pot.use_solv)fprintf(pcmoutfile,"Solv: %f\n\n",energies.esolv);
      }
/* write structure number information  */
        files.icurrent = i;
        generate_bonds();
         
        minimize();
        optimize_data.final_energy = energies.total;

// Hay debug
//      if (field.type == MMX || field.type == MM3)
//         eheat();

        optimize_data.final_heat = get_heat_of_formation();
        dipolemom.xdipole = 0.0;
        dipolemom.ydipole = 0.0;
        dipolemom.zdipole = 0.0;
        if (pot.use_charge) charge_dipole();
        if (pot.use_dipole) dipole_dipole();
        dipolemom.total = sqrt(dipolemom.xdipole*dipolemom.xdipole +
                             dipolemom.ydipole*dipolemom.ydipole +
                             dipolemom.zdipole*dipolemom.zdipole);
        
      if (minim_values.iprint)
      {
         fprintf(pcmoutfile,"\nFinal Energy : %f\n",energies.total);
         fprintf(pcmoutfile,"Str: %f    Bnd: %f    Tor: %f\n", energies.estr, energies.ebend, energies.etor);
         fprintf(pcmoutfile,"StrBnd: %f VDW: %f     QQ: %f\n",energies.estrbnd,energies.evdw+energies.e14+energies.ehbond, energies.eu);
         fprintf(pcmoutfile,"OOP: %f AA: %f     Strtor: %f\n",energies.eopb,energies.eangang, energies.estrtor);
         fprintf(pcmoutfile,"Urey: %f\n",energies.eurey);
         if (pot.use_solv)fprintf(pcmoutfile,"Solv: %f\n\n",energies.esolv);
      }
        gmmx.ecurrent = energies.total;
        gmmx.neminim++;
        nrc = check_id();
        if (nrc == TRUE)
        {
            if (gmmx.ecurrent < (gmmx.eminim+gmmx_data.ewindow2) )
            {
               sprintf(Struct_Title,"%i %9.3f %9.3f",gmmx.nfound + 1,gmmx.ecurrent,get_heat_of_formation());
               if (i == 1)
                  files.append = FALSE;
               else
                  files.append = TRUE;
                  
               pcmfout(1);
               
               gmmx.nfound++;
               gmmx.isitdupl = UNIQUE;
               if (gmmx.ecurrent < gmmx.eminim)
               {
                   gmmx.eminim = gmmx.ecurrent;
                   gmmx.nsflow = gmmx.nsaved;
               }
            }else
            {
                gmmx.isitdupl = UNIQUEHIGH;
            }
        }
        write_gmmx();
        nrc = check_stop();
        if (nrc == FALSE)
            break;

     }
     end_calculation();
     files.batch = False;
    // release all memory
    rdfile(-1);
    nrc = bstereo(-1);
    nrc = lepimer(-1);
    hster(-1, jjii);
    free_ivector(jjii, 1, natom+1);

   //  need to output results and do queries here
    if (gmmx.job_done == TRUE)
    {
        sort_file(gmmx.tmpname,gmmx_data.foname2);
        remove(gmmx.rmsname);
        remove(gmmx.tmpname);
    }
//modgraph    bkboltz();
    gmmx.run = FALSE;
    // convert to sdf
}
/* =============================================== */
void bkboltz()
{
    int i, j, k, ia1, ia2, iftype;
    float temp, blzk, sumprob, hfsumprob, elow, hflow,rt, pct, prob, hfprob;
    float dene, hfdene, eave, hfeave, dpmeave, dpmsqr, entropy;
    float entropy1, entropy2, hfpct;
    float cdis[50], cang[50], cdihed[50][36], cpmr[50], cdis16[50];
    double rdist, rdihed, rdist16;
    float rang, ang;
    int icount, idihed[50][36];
    char inputline[90],dummy1[20], dummy2[20];
    float energy, heat;
    
    FILE *wfile,*rfile;

    rfile = fopen_path(Openbox.path,gmmx_data.foname2,"r");

    icount = 0;
    rdist = 0.0;

    while ( fgets(inputline,90,rfile) != NULL)
    {
        if (strncasecmp(inputline,"{PCM",4) == 0)
        {
                sscanf(inputline,"%s %s %f %f",dummy1, dummy2, &energy, &heat);
                gmmx.econf[icount] = energy;
                gmmx.hfconf[icount] = heat;
                icount++;
        }
    }
    fclose(rfile);
    // reset nfound
    gmmx.nfound = icount;

    temp = gmmx_data.boltz;
    blzk = 1.988e-3;
    sumprob = 0.0;
    hfsumprob = 0.0;
    elow = 9999.0;
    hflow = 9999.0;
    rt = temp*blzk;

    if (gmmx_data.ndist > 0)
    {
        for (i=0; i < gmmx_data.ndist; i++)
        {
           cdis[i] = 0.0;
           cdis16[i] = 0.0;
         }
    }
    if (gmmx_data.nang > 0)
    {
        for (i=0; i < gmmx_data.nang; i++)
           cang[i] = 0.0;
    }
    if (gmmx_data.ndihed > 0)
    {
        for (i=0; i < gmmx_data.ndihed; i++)
         for (j=0; j < 36; j++)
          {
           cdihed[i][j] = 0.0;
           idihed[i][j] = 0;
          }
    }
    if (gmmx_data.npmr > 0)
    {
        for (i=0; i < gmmx_data.npmr; i++)
           cpmr[i] = 0.0;
    }
    
     wfile = fopen_path(Savebox.path,gmmx.pkmname,"w");

    // loop through econf and find lowest energy and hf
    for(i=0; i < gmmx.nfound; i++)
    {
        if (gmmx.econf[i] < elow)
            elow = gmmx.econf[i];
        if (gmmx.hfconf[i] < hflow)
            hflow = gmmx.hfconf[i];
    }

    for(i=0; i < gmmx.nfound; i++)
    {
        dene = gmmx.econf[i] - elow;
        prob = exp((-dene/rt));
        sumprob += prob;
        hfdene = gmmx.hfconf[i] - hflow;
        hfprob = exp((-hfdene/rt));
        hfsumprob += hfprob;
    }
    
    eave = 0.0;
    hfeave = 0.0;
    dpmeave = 0.0;
    dpmsqr = 0.0;
    entropy = 0.0;

    for (i=0; i < gmmx.nfound; i++)
    {
        dene = gmmx.econf[i] - elow;
        prob = exp(-dene/rt);
        pct = prob / sumprob;
        eave += gmmx.econf[i] * pct;
        
        hfdene = gmmx.hfconf[i] - hflow;
        hfprob = exp(-hfdene/rt);
        hfpct = hfprob / hfsumprob;
        hfeave += gmmx.hfconf[i] * hfpct;
        if( gmmx_data.heat)
        {
          dpmeave += gmmx.dpmconf[i] * hfpct;
          dpmsqr  += hfpct*gmmx.dpmconf[i]*gmmx.dpmconf[i];
          entropy += hfpct*log(hfpct);
        } else
        {
          dpmeave += gmmx.dpmconf[i] * pct;
          dpmsqr += pct*gmmx.dpmconf[i]*gmmx.dpmconf[i];
          entropy += pct*log(pct);
        }
    }

    fprintf(wfile," GMMX Conformation Search Output Results\n");
    fprintf(wfile," Temperature for Boltzmann Calculations: %6.2f\n",temp);
    fprintf(wfile," Average Energy: %8.3f Kcal/mol\n", eave);
    fprintf(wfile," Average Heat of Formation: %8.3f Kcal/mol\n", hfeave);
    fprintf(wfile," Average Dipole Moment: %6.2f\n",dpmeave);
    fprintf(wfile," Dipole Moment Squared: %6.2f\n",dpmsqr);

     iftype = Openbox.ftype = FTYPE_PCM;
     strcpy(Openbox.fname,gmmx_data.foname2);
     if (iftype == FTYPE_PCM)
     {
         rdfile(0);
     }

    for (i=0; i < gmmx.nfound; i++)
    {
        initialize();
        rdfile(i+1);
        dene = gmmx.econf[i] - elow;
        prob = exp((-dene/rt));
        pct = prob/sumprob *100.0;
        hfdene = gmmx.hfconf[i] - hflow;
        hfprob = exp((-hfdene/rt));
        hfpct = hfprob/hfsumprob * 100.0;

        if (gmmx_data.heat)
          pct = hfpct;

        fprintf(wfile,"Pop conf %5d = %7.2f  E = %8.3f DE = %5.3f HF = %8.3f DHF = %5.3f DPM = %5.3f\n",i+1, pct,
           gmmx.econf[i],dene,gmmx.hfconf[i], hfdene,gmmx.dpmconf[i]);
        
        if (gmmx_data.ndist > 0)
        {
            for(j=0; j < gmmx_data.ndist; j++)
            {
                rdist = distance(gmmx_data.dist_atoms[j][0], gmmx_data.dist_atoms[j][1]);
                rdist16 = 1/(rdist*rdist*rdist*rdist*rdist*rdist);
                cdis[j] += rdist*pct*0.01;
                cdis16[j] += rdist16*pct*0.01;
                fprintf(wfile," Distance ( %4d - %4d) = %6.2f Ang, 1/(r**6) = %7.5f\n",gmmx_data.dist_atoms[j][0], gmmx_data.dist_atoms[j][1],rdist,rdist16);
            }
        }

        if (gmmx_data.nang > 0)
        {
            for(j=0; j < gmmx_data.nang; j++)
            {
                angle(gmmx_data.ang_atoms[j][0], gmmx_data.ang_atoms[j][1], gmmx_data.ang_atoms[j][2], &rang);
                cang[j] += rang*pct*0.01;
                fprintf(wfile," Angle ( %4d - %4d - %4d) = %6.2f Deg\n",gmmx_data.ang_atoms[j][0], gmmx_data.ang_atoms[j][1],
                   gmmx_data.ang_atoms[j][2],rang);
            }
        }

        if (gmmx_data.ndihed > 0)
        {
            for(j=0; j < gmmx_data.ndihed; j++)
            {
                rdist = dihdrl(gmmx_data.dihed_atoms[j][0], gmmx_data.dihed_atoms[j][1],
                    gmmx_data.dihed_atoms[j][2], gmmx_data.dihed_atoms[j][3]);
                fprintf(wfile," Dihedral ( %4d - %4d - %4d - %4d) = %6.2f Deg\n",gmmx_data.dihed_atoms[j][0], gmmx_data.dihed_atoms[j][1],
                   gmmx_data.dihed_atoms[j][2],gmmx_data.dihed_atoms[j][3],rdist);
                ang = -180.0;
                for(k=0; k < 35; k++)
                 {
                  if (rdist >= ang && rdist < (ang + 10.0))
                   {
                      idihed[j][k] ++;
                      cdihed[j][k] += pct;
                   }
                   ang += 10.0;
                 }
                  if (rdist >= 170.0)
                   {
                      idihed[j][35] ++;
                      cdihed[j][35] += pct;
                   }
            }
        }

        if (gmmx_data.npmr > 0)
        {
            for(j=0; j < gmmx_data.npmr; j++)
            {
                ia1 = atom[gmmx_data.pmr_atoms[j][0]].iat[0];
                ia2 = atom[gmmx_data.pmr_atoms[j][1]].iat[0];
//                rdist = couple(gmmx_data.pmr_atoms[j][0], ia1, ia2,gmmx_data.pmr_atoms[j][1]);
                rdihed = dihdrl(gmmx_data.pmr_atoms[j][0], ia1, ia2,gmmx_data.pmr_atoms[j][1]);
                cpmr[j] += rdist*pct*0.01;
                fprintf(wfile," Coupling Constant ( %4d - %4d - %4d - %4d) = %5.2f Hz ( %6.2f deg)\n",gmmx_data.pmr_atoms[j][0],
                   ia1, ia2, gmmx_data.pmr_atoms[j][1],rdist, rdihed);
            }
        }
    }


    if (gmmx_data.ndist > 0)
    {
        fprintf(wfile,"\nThe Boltzmann Averaged distances are:\n");
        for (i=0; i < gmmx_data.ndist; i++)
           fprintf(wfile," Distance ( %4d - %4d ) = %6.2f Ang, 1/(r**6) = %7.5f\n",gmmx_data.dist_atoms[i][0], gmmx_data.dist_atoms[i][1],cdis[i],cdis16[i]);
    }

    if (gmmx_data.nang > 0)
    {
        fprintf(wfile,"\nThe Boltzmann Averaged Angles are:\n");
        for (i=0; i < gmmx_data.nang; i++)
           fprintf(wfile," Angle ( %4d - %4d - %4d ) = %6.2f Deg\n",gmmx_data.ang_atoms[i][0], gmmx_data.ang_atoms[i][1],
                gmmx_data.ang_atoms[i][2], cang[i]);
    }
           
    if (gmmx_data.ndihed > 0)
    {
        fprintf(wfile,"\nThe Boltzmann Averaged Dihedrals are:\n");
        for (i=0; i < gmmx_data.ndihed; i++)
         {
           fprintf(wfile," Angle ( %4d - %4d - %4d - %4d )\n",gmmx_data.dihed_atoms[i][0], gmmx_data.dihed_atoms[i][1],
                gmmx_data.dihed_atoms[i][2], gmmx_data.dihed_atoms[i][3]);
           ang = -180.0;
           for (j=0; j < 18; j++)
            {
                fprintf(wfile,"%6.1f - %6.1f %5d  %6.1f  %6.1f - %6.1f %5d  %6.1f\n",ang, ang + 10.0, idihed[i][j], cdihed[i][j],
                         ang + 180.0, ang + 190.0, idihed[i][j + 18], cdihed[i][j + 18]);
                ang = ang + 10.0;
             }
         }
    }

    if (gmmx_data.npmr > 0)
    {
        fprintf(wfile,"\nThe Boltzmann Averaged Coupling Constants are:\n");
        for (i=0; i < gmmx_data.npmr; i++)
           fprintf(wfile," Coupling Constant ( %4d - %4d ) = %5.2f Hz\n",gmmx_data.pmr_atoms[i][0], gmmx_data.pmr_atoms[i][1],cpmr[i]);
    }
    
    entropy1 = -1.9872*entropy;
    entropy2 = -rt*entropy;
    fprintf(wfile,"\nEntropy of mixing is: %6.2f cal/mol-K\n",entropy1);
    fprintf(wfile,"At %6.2f deg K,   TdS = %6.2f Kcal/mol\n",temp,entropy2);
    if (iftype == FTYPE_PCM)
        rdfile(-1);
    fclose(wfile);
}
// =============================================
/* ======================================================= */
void jacobi(int n,int np,double **a, double *d,double **v,double *b,double *z)
{
      int i,j,k,ip,iq,nrot,maxrot;
      double sm,tresh,s,c,t,theta,tau,h,g,p;

      maxrot = 100;
      nrot = 0;
      for (i=0; i < np; i++)
      {
          for (j=0; j < np; j++)
             v[j][i] = 0.0;
          v[i][i] = 1.0;
      }
      for (i=0; i < np; i++)
      {
          b[i] = a[i][i];
          d[i] = b[i];
          z[i] = 0.0;
      }
//
      for (i=0; i < maxrot; i++)
      {
          sm = 0.0;
          for (ip =0; ip < n-1; ip++)
          {
              for (iq=ip+1; iq < n; iq++)
              {
                  sm += fabs(a[iq][ip]);
              }
          }
          if (fabs(sm) < 0.001) goto L_10;
          if (i < 3)
             tresh = .2*sm/(n*n);
          else
             tresh = 0.0;
//
          for (ip =0; ip < n-1; ip++)
          {
              for (iq=ip+1; iq < n; iq++)
              {
                  g = 100.0 * fabs(a[iq][ip]);
                  if (i > 3 && fabs(d[ip])+g == fabs(d[ip])
                        && fabs(d[iq])+g == fabs(d[iq]))
                  {
                      a[iq][ip] = 0.0;
                  } else if (fabs(a[iq][ip]) > tresh)
                  {
                     h = d[iq] - d[ip];
                     if (fabs(h)+g == fabs(h))
                     {
                        t = a[iq][ip] / h;
                     }else
                     {
                        theta = 0.5*h / a[iq][ip];
                        t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
                        if (theta < 0.0)  t = -t;
                     }
                     c = 1.0 / sqrt(1.0+t*t);
                     s = t * c;
                     tau = s / (1.0+c);
                     h = t * a[iq][ip];
                     z[ip] -= h;
                     z[iq] += h;
                     d[ip] -= h;
                     d[iq] += h;
                     a[iq][ip] = 0.0;
                     for (j=0; j < ip-1; j++)
                     {
                        g = a[ip][j];
                        h = a[iq][j];
                        a[ip][j] = g - s*(h+g*tau);
                        a[iq][j] = h + s*(g-h*tau);
                     }
                     for (j=ip+1; j < iq-1; j++)
                     {
                        g = a[j][ip];
                        h = a[iq][j];
                        a[j][ip] = g - s*(h+g*tau);
                        a[iq][j] = h + s*(g-h*tau);
                     }
                     for (j=iq+1; j < n; j++)
                     {
                        g = a[j][ip];
                        h = a[j][iq];
                        a[j][ip] = g - s*(h+g*tau);
                        a[j][iq] = h + s*(g-h*tau);
                     }
                     for (j=0; j < n; j++)
                     {
                        g = v[ip][j];
                        h = v[iq][j];
                        v[ip][j] = g - s*(h+g*tau);
                        v[iq][j] = h + s*(g-h*tau);
                     }
                     nrot++;
                  }
              }
          }
          for (ip = 0; ip < n; ip++)
          {
            b[ip] = b[ip] + z[ip];
            d[ip] = b[ip];
            z[ip] = 0.0;
          }
      }
L_10:
      if (nrot == maxrot)
         message_alert("Error in Jacobi","Error");
// sort eigenvalues
      for (i=0; i < n-1; i++)
      {
         k = i;
         p = d[i];
         for (j=i+1; j < n; j++)
         {
            if (d[j] < p)
            {
               k = j;
               p = d[j];
            }
         }
         if (k != i)
         {
            d[k] = d[i];
            d[i] = p;
            for (j=0; j < n; j++)
            {
               p = v[i][j];
               v[i][j] = v[k][j];
               v[k][j] = p;
            }
         }
      }
}
// ============================================
void write_gmmx()
{
    char lpstring[30];
    float rate;

        switch (gmmx.isitdupl)
        {
           case UNIQUE:    
                strcpy(lpstring,"U: Kept");
                break;
           case UNIQUEHIGH:
                strcpy(lpstring,"U; Rej Energy");
                break;
           case IDENT:
                strcpy(lpstring,"Ident");
                break;
           case BADEPIMER:
                strcpy(lpstring,"R; Epimer");
                break;
           case BADCONDIST:
                strcpy(lpstring,"R; Const Dist");
                break;
           case BADCONTOR:
                strcpy(lpstring,"R; Const Tor");
                break;
           case BADDBOND:
                strcpy(lpstring,"R; DBond Isomer");
                break;
           case NUMISOMER:
                strcpy(lpstring,"R; Num Isomer");
                break;
           case ENANTIOMER:
                strcpy(lpstring,"R; Enantiomer");
                break;
           case BADHYBRID:
                strcpy(lpstring,"R; Hybrid");
                break;
           case BADENERGY:
                strcpy(lpstring,"R; High Energy");
                break;
           default:
                break;
        }
    if (gmmx.neminim > 0)
      rate = (float)gmmx.nfound/gmmx.neminim;
    else
      rate = 0.0;
    rate *= 100.;
//    if ( gmmx.neminim%5 == 0)
//    {
//        if (gmmx.smethod == BONDS)
//        {
//           printf("Bond - Nmin: %d Nfnd: %d Ecur: %f Emin: %f Rate: %f Result: %s  Dupl: %d   High: %d/%d\n",
//              gmmx.neminim,gmmx.nfound,gmmx.ecurrent,gmmx.eminim,rate,lpstring,gmmx.ndupl,gmmx.maxdupl,gmmx_data.kdup);
//        } else
//        {
//           printf("Cart - Nmin: %d Nfnd: %d Ecur: %f Emin: %f Rate: %f Result: %s  Dupl: %d   High: %d/%d\n",
//             gmmx.neminim,gmmx.nfound,gmmx.ecurrent,gmmx.eminim,rate,lpstring,gmmx.ndupl,gmmx.maxdupl,gmmx_data.kdup);
//        }
//    }
}

