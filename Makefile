
CC     = icc
CFLAGS = -O  -axK -openmp -ip -g #(icc)  -DSTRESS -DEAM -DPARABEL -DLIMIT
#CFLAGS = -g -DEAM -DSTRESS -DPARABEL  # -DEAM  -DLIMIT#gcc -O -g3
BINDIR = ${HOME}/bin/${HOSTTYPE}

POTFITSRC = f1dim_r.c powell_lsq.c lubksb_r.c  mprove_r.c brent_r.c ludcmp_r.c linmin_r.c mnbrak_r.c nrutil_r.c force.c config.c param.c potential.c potfit.c splines.c simann.c
POTFITHDR = potfit.h powell_lsq.h nrutil_r.h 
POTFITOBJ = $(POTFITSRC:.c=.o)

POTSCALESRC = param.c potential.c splines.c potscale.c nrutil_r.c 
POTSCALEHDR = potscale.h nrutil_r.h
POTSCALEOBJ = $(POTSCALESRC:.c=.o)

.c.o: ${POTFITHDR}
	$(CC) $(CFLAGS) -c $<

#${POTFITOBJ}: $($@:.o=.c) ${POTFITHDR}
#	${CC} ${CFLAGS} -c $($@:.o=.c)
potfit_par: ${POTFITOBJ} ${POTFITHDR}
	${CC} ${CFLAGS} -o ${BINDIR}/potfit_par ${POTFITOBJ} -lm

potfit2: ${POTFITOBJ} ${POTFITHDR}
	${CC} ${CFLAGS} -o ${BINDIR}/potfit_e ${POTFITOBJ} -lm

potfit: ${POTFITOBJ} ${POTFITHDR}
	${CC} ${CFLAGS} -o ${BINDIR}/$@ ${POTFITOBJ} -lm

#parab:  ${POTFITOBJ} ${POTFITHDR}
#	${CC} ${CFLAGS} -DPARABEL -o ${BINDIR}/$@ ${POTFITOBJ} -lm

clean:
	rm -f *.o *~ 

neightab: neightab.c
	${CC} ${CFLAGS} -o ${BINDIR}/neightab neightab.c -lm

pottrans: pottrans.c
	${CC} ${CFLAGS} -o ${BINDIR}/pottrans pottrans.c -lm

force2poscar: force2poscar.c
	${CC} ${CFLAGS} -o ${BINDIR}/force2poscar force2poscar.c

potscale: ${POTSCALEOBJ} ${POTSCALEHDR}
	${CC} ${CFLAGS} -o ${BINDIR}/$@ ${POTSCALEOBJ} -lm
