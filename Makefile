
CC     = icc
CFLAGS = -O  -axK -openmp -ip -g -DEAM#(icc)
#CFLAGS = -g  -DEAM #gcc -O -g3
BINDIR = ${HOME}/bin/${HOSTTYPE}

POTFITSRC = f1dim_r.c powell_lsq.c lubksb_r.c  mprove_r.c brent_r.c ludcmp_r.c linmin_r.c mnbrak_r.c nrutil_r.c force.c config.c param.c potential.c potfit.c splines.c simann.c
POTFITHDR = potfit.h powell_lsq.h nrutil_r.h
POTFITOBJ = $(POTFITSRC:.c=.o)

.c.o: ${POTFITHDR}
	$(CC) $(CFLAGS) -c $<

#${POTFITOBJ}: $($@:.o=.c) ${POTFITHDR}
#	${CC} ${CFLAGS} -c $($@:.o=.c)

#potfit2: ${POTFITOBJ} ${POTFITHDR}
#	${CC} ${CFLAGS} -o ${BINDIR}/potfit_i ${POTFITOBJ} -lm

potfit: ${POTFITOBJ} ${POTFITHDR}
	${CC} ${CFLAGS} -o ${BINDIR}/$@ ${POTFITOBJ} -lm

clean:
	rm -f *.o *~ 

neightab: neightab.c
	${CC} ${CFLAGS} -o ${BINDIR}/neightab neightab.c -lm

pottrans: pottrans.c
	${CC} ${CFLAGS} -o ${BINDIR}/pottrans pottrans.c -lm

force2poscar: force2poscar.c
	${CC} ${CFLAGS} -o ${BINDIR}/force2poscar force2poscar.c
