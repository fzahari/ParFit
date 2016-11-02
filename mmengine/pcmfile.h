#ifndef PCMFILE_H
#define PCMFILE_H

struct t_pcmfile {
        char string[200];
        int head;
        char token[80];
        int state;
        unsigned int nocaps;
        }       pcmfile;

#endif
