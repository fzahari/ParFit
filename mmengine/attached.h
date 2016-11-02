#ifndef ATTACHED_H
#define ATTACHED_H

struct t_attached {
       int  *n13, *n14;
       int  **i13, **i14;  //**i13[12][MAXATOM ], **i14[144][MAXATOM ];
       } attached;

#endif
void attach(void);

