CC = gcc

# For FreeBSD uncomment the following
# CFLAGS = -g -Wpedantic -O -I/usr/X11R6/include -DFREEBSD

# For RedHat Linux; if you want warnings trys "-Wall" and 
# for debugging "-g".  Most of us (users) like it like this: 
CFLAGS = -g -O3 -ftree-vectorize -Wall
#CFLAGS = -g -O1 -Wall
#CFLAGS = -g -Wall
#CFLAGS = -g -Wall -Wstrict-prototypes -pg -O -I/usr/local/LessTif/Motif1.2/include
#CFLAGS = -g -Wall -Wstrict-prototypes -pg -O -I/usr/local/Motif2.0/include

# Red Hat with Open Group's Open Motif 
#XLIBDIRS= -L/usr/X11R6/lib

# EXTLIBS = -lm -lXi -lXext -lX11
#LDFLAGS =  ${LD_DIR1}libXm.so.1.0.2 ${LD_DIR1}libXt.so.6.0 ${LD_DIR1}libX11.so.6.1 -lm

# For Statically linked binary
LDFLAGS = -lm

# For normal shared library version
#LDFLAGS =  -lXm -lXt -lX11 -lm

RM = rm -f
CP = cp -f

# destination dir; set this for your setup
#BIN_DIR = /usr/local/bin

DST = engine

SRCS=  \
active.c  angles.c  attach.c  column.c  datastat.c  diagc.c \
dipmom.c  draw.c  eangang.c  eangle.c  ebond.c  ebuckmm3.c \
ebufchrg.c  echarge.c  echgdpl.c  ecoordb.c  edipole.c  egeom.c \
ehal.c  eheat.c  ehighcoord.c  eopbend.c  eobpw.c  engine.c \
engine1.c  esolv.c  estrbnd.c  estrtor.c  etorsion.c field.c \
fnd_bond.c  fnd_ring.c  get_mem.c  gmmx_run.c  gmx_run2.c \
hessian.c  initiali.c  l.c  k.c  kangang.c  kangle.c \
kbond.c  kcharge.c  kcoord.c  kdipole.c  kopend.c  ksolv.c \
kstrbnd.c  kstrtor.c  ktorsion.c  kvdw.c  mark.c  minimize.c \
numeral.c  orbit.c  pcm7.c  pcmwin1.c  pcmwin2.c  read.c \
read_pcm.c  readprm.c  read_mm2.c  rings.c  rotbond.c  search.c \
solve.c  tncg.c  torsions.c  type_mmx.c  u.c  utility.c 


OBJS= \
active.o  angles.o  attach.o  column.o  datastat.o  diagc.o \
dipmom.o  draw.o  eangang.o  eangle.o  ebond.o  ebuckmm3.o \
ebufchrg.o  echarge.o  echgdpl.o  ecoordb.o  edipole.o  egeom.o \
ehal.o  eheat.o  ehighcoord.o  eopbend.o  eobpw.o  engine.o \
engine1.o  esolv.o  estrbnd.o  estrtor.o  etorsion.o field.o \
fnd_bond.o  fnd_ring.o  get_mem.o  gmmx_run.o  gmx_run2.o \
hessian.o  initiali.o  l.o  k.o  kangang.o  kangle.o \
kbond.o  kcharge.o  kcoord.o  kdipole.o  kopend.o  ksolv.o \
kstrbnd.o  kstrtor.o  ktorsion.o  kvdw.o  mark.o  minimize.o \
numeral.o  orbit.o  pcm7.o  pcmwin1.o  pcmwin2.o  read.o \
read_pcm.o  readprm.o  read_mm2.o  rings.o  rotbond.o search.o \
solve.o  tncg.o  torsions.o  type_mmx.o  u.o  utility.o 

       
${DST}: $(OBJS)
	${CC} -g -o ${DST} ${OBJS} ${LDFLAGS}


install:
	@${CP} ${DST} ${BIN_DIR}
	@echo "Copying executable ${DST} to ${BIN_DIR}..."
	@echo "Don't forget to set the PCM_DIR environment variable."

clean:
	${RM} ${OBJS}
