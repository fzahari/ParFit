#ifndef ENERGIES_H
#define ENERGIES_H

struct  t_energies {
        double total, estr, ebend, estrbnd, e14, evdw, etor, eu, eopb,
               eangang, estrtor, ehbond, efix, eimprop, eimptors, eurey, esolv; } energies;

struct t_virial {
        double virx, viry, virz;
         } virial;
 #endif
