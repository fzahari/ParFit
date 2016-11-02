#ifndef FIELD_H
#define FIELD_H

struct t_field {
  char    name[30];
  int     type;
  float   bondunit, bond_cubic, bond_quartic;
  float   angleunit, angle_cubic, angle_quartic, angle_pentic, angle_sextic;
  float   str_bndunit, ang_angunit, str_torunit, torsionunit;
  char    vdwtype[20], radiusrule[20], radiustype[20], radiussize[20], epsrule[20];
  float   vdw_14scale, a_expterm, b_expterm, c_expterm;
  float   chg_14scale, dielectric;
  float pos_12scale,pos_13scale;
        } field;
 #endif
void set_field(void);
void potoff(void);
