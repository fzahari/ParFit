#ifndef POT_H
#define POT_H

struct t_pot {
        int use_bond, use_angle, use_strbnd, use_urey, use_angang;
        int use_opbend, use_improp, use_imptor, use_tors, use_strtor;
        int use_tortor, use_vdw, use_lj, use_buck, use_hal, use_gauss;
        int use_charge, use_bufcharge, use_chrgdpl, use_dipole, use_polar, use_solv;
        int use_geom, use_extra, use_picalc, use_hbond, use_coordb, use_opbend_wilson;
        int use_bounds, use_image, use_replica, use_deform, use_ewald, use_highcoord;
        }  pot;

#endif

