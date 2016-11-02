#ifndef ELEMENTTYPE_H
#define ELEMENTTYPE_H

struct ElementType { char symbol[3];
                     int atomnum;
                     float weight, covradius, vdwradius;
                     int s,p,d,f, type;
                   } Elements[];
 #endif
