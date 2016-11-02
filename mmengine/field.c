#include "pcwin.h"
#include "pcmod.h"
#include "pot.h"
#include "field.h"
#include "minim_values.h"

void set_field()
{
      pot.use_bond = FALSE;
      pot.use_angle = FALSE;
      pot.use_strbnd = FALSE;
      pot.use_urey = FALSE;
      pot.use_angang = FALSE;
      pot.use_opbend = FALSE;
      pot.use_improp = FALSE;
      pot.use_imptor = FALSE;
      pot.use_tors = FALSE;
      pot.use_strtor = FALSE;
      pot.use_tortor = FALSE;
      pot.use_vdw = FALSE;
      pot.use_lj = FALSE;
      pot.use_buck = FALSE;
      pot.use_hal = FALSE;
      pot.use_gauss = FALSE;
      pot.use_charge = FALSE;
      pot.use_bufcharge = FALSE;
      pot.use_chrgdpl = FALSE;
      pot.use_dipole = FALSE;
      pot.use_polar = FALSE;
      pot.use_geom = FALSE;
      pot.use_extra = FALSE;
      pot.use_picalc = FALSE;
      pot.use_hbond = FALSE;
      pot.use_coordb = FALSE;
      pot.use_opbend_wilson = FALSE;
      pot.use_highcoord = FALSE;

    if (field.type == MMX )
    {
        pot.use_bond = TRUE;
        pot.use_angle = TRUE;
        pot.use_strbnd = TRUE;
        pot.use_opbend = TRUE;
        pot.use_opbend_wilson = FALSE;
        pot.use_tors = TRUE;
        pot.use_vdw = TRUE;
        pot.use_buck = TRUE;
        pot.use_lj = FALSE;
        pot.use_hal = FALSE;
        pot.use_gauss = FALSE;        
        pot.use_bufcharge= FALSE;
        if (minim_values.ndc == 4)
        {
             pot.use_charge = TRUE;
             pot.use_dipole = FALSE;
       }else
        {
             pot.use_charge = FALSE;
             pot.use_dipole = TRUE;
       }
        pot.use_urey = FALSE;
        
        pot.use_angang = FALSE;
        pot.use_improp = FALSE;
        pot.use_imptor = FALSE;
        pot.use_strtor = FALSE;
        pot.use_tortor = FALSE;
    } else if (field.type == MM3)
    {
        pot.use_bond = TRUE;
        pot.use_angle = TRUE;
        pot.use_strbnd = TRUE;
        pot.use_opbend = TRUE;
        pot.use_opbend_wilson = FALSE;
        pot.use_tors = TRUE;
        pot.use_vdw = TRUE;
        pot.use_buck = TRUE;
        pot.use_lj = FALSE;
        pot.use_hal = FALSE;
        pot.use_hbond = FALSE;
        
        pot.use_gauss = FALSE;
        if (minim_values.ndc == 4)
        {
             pot.use_charge = TRUE;
             pot.use_dipole = FALSE;
       }else
        {
             pot.use_charge = FALSE;
             pot.use_dipole = TRUE;
       }
        pot.use_bufcharge= FALSE;
        pot.use_urey = FALSE;
        pot.use_improp = FALSE;
        pot.use_imptor = FALSE;
        pot.use_angang = TRUE;
        pot.use_strtor = TRUE;
    } else if (field.type == MMFF94)
    {
        pot.use_bond = TRUE;
        pot.use_angle = TRUE;
        pot.use_strbnd = TRUE;
        pot.use_opbend = FALSE;
        pot.use_opbend_wilson = TRUE;
        pot.use_tors = TRUE;
        
        pot.use_vdw = FALSE;
        pot.use_buck = FALSE;
        pot.use_lj = FALSE;
        pot.use_hal = TRUE;
        
        pot.use_hbond = FALSE;
        
        pot.use_gauss = FALSE;
        pot.use_charge = FALSE;
        pot.use_bufcharge= TRUE;
        pot.use_dipole = FALSE;
        pot.use_urey = FALSE;
        pot.use_improp = FALSE;
        pot.use_imptor = FALSE;
        pot.use_angang = FALSE;
        pot.use_strtor = FALSE;
    }
}
void potoff()
{
      pot.use_bond = FALSE;
      pot.use_angle = FALSE;
      pot.use_strbnd = FALSE;
      pot.use_urey = FALSE;
      pot.use_angang = FALSE;
      pot.use_opbend = FALSE;
      pot.use_improp = FALSE;
      pot.use_imptor = FALSE;
      pot.use_tors = FALSE;
      pot.use_strtor = FALSE;
      pot.use_tortor = FALSE;
      pot.use_vdw = FALSE;
      pot.use_lj = FALSE;
      pot.use_buck = FALSE;
      pot.use_hal = FALSE;
      pot.use_gauss = FALSE;
      pot.use_charge = FALSE;
      pot.use_bufcharge = FALSE;
      pot.use_chrgdpl = FALSE;
      pot.use_dipole = FALSE;
      pot.use_polar = FALSE;
      pot.use_solv = FALSE;
      pot.use_geom = FALSE;
      pot.use_extra = FALSE;
      pot.use_picalc = FALSE;
      pot.use_hbond = FALSE;
      pot.use_coordb = FALSE;
      pot.use_opbend_wilson = FALSE;
}

