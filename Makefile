
CC     = gcc
CFLAGS = -O -g3
BINDIR = ${HOME}/bin/${HOSTTYPE}

POTFITSRC = f1dim_r.c powell_lsq.c lubksb_r.c  mprove_r.c brent_r.c ludcmp_r.c linmin_r.c mnbrak_r.c nrutil_r.c force.c config.c param.c potential.c potfit.c 
POTFITHDR = potfit.h powell_lsq.h nrutil_r.h
POTFITOBJ = $(subst .c,.o,${POTFITSRC})

%.o: %.c ${POTFITHDR}
	${CC} ${CFLAGS} -c $(subst .o,.c,$@)

potfit: ${POTFITOBJ} ${POTFITHDR}
	${CC} ${CFLAGS} -o ${BINDIR}/$@ ${POTFITOBJ} -lm

clean:
	rm -f *.o *~ 

neightab: neightab.c
	${CC} ${CFLAGS} -o ${BINDIR}/neightab neightab.c -lm

pottrans: pottrans.c
	${CC} ${CFLAGS} -o ${BINDIR}/pottrans pottrans.c -lm
