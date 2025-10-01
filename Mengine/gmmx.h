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

#ifndef GMMX_H
#define GMMX_H

#define MAXFND         30000
#define UNIQUE         0
#define IDENT          1
#define UNIQUEHIGH     2
#define BADEPIMER      3
#define BADCONDIST     4
#define BADCONTOR      5
#define BADDBOND       6
#define NUMISOMER      7
#define ENANTIOMER     8
#define BADHYBRID      9
#define BADENERGY     10

#define BONDS          1
#define CART           2

#define RUN            90
#define STOP           99
#define PAUSE          98

int gmmx_abort;

struct t_gmmx_data {
    char jobname[256], finame[256], foname[256], foname2[256];
    char include_name[256];
    char gmmx_comment[64];
    int restart, include_file, include_ftype;
    int method, iseed, its;
    int lnant, heat, nopi, chig, bad15, hbond, ecut, hybrid;
    int kstop, kdup,kmin, max_search;
    int comp, ncompatoms, comp_atoms[50];
    int nsrms,nsrmsa,comp_method,comp_method2;
    int qpmr, qdist, qang, qdihed;
    int npmr, ndist, nang, ndihed;
    int pmr_atoms[50][2], dist_atoms[50][2];
    int ang_atoms[50][3], dihed_atoms[50][4];
    int nrings, rng_size[4], ring_atoms[30][4];
    int clo_bond[8][4], ring_resolution[4];
    float clo_ang[2][4], clo_distance[4];
    int nring_bonds[4], ring_bond_data[4][30][3];
    int nbonds, bond_data[50][5];
    float ewindow, boltz, min_contact, ecutoff, bad15_cutoff, ewindow2;
    float ermsa, crmsa;
    int ff;
    } gmmx_data;

struct t_gmmx {
    char tmpname[256], rmsname[256], pkmname[256], rstname[256];
    int run;
    int namvt, nrbond, smethod;
    int hybrida;
    float ebendold;
    int neminim, nfound, nsfmin, nsflow, ndupl,maxdupl, isitdupl;
    int nsaved, nlast, nlastk;
    int nsfound[MAXFND];
    double eminim, ecurrent, elast;
    float econf[MAXFND],hfconf[MAXFND],dpmconf[MAXFND];
    int job_done;
    } gmmx;

struct t_gmmxring{
    int nrings, nring_atoms[200],ring_atoms[200][30];
    int nbranch,branch_atoms[200],in_rings[200][10];
    int npairs,pair[200][4];
    int nfused,fused[200];
    }  gmmxring;

#endif
