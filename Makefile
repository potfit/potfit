
CC     =gcc
CFLAGS = -g
BINDIR = ${HOME}/bin/${HOSTTYPE}

POTFITSRC = force.c config.c param.c potential.c potfit.c
POTFITHDR = potfit.h
POTFITOBJ := $(subst .c,.o,${POTFITSRC})

%.o: %.c ${POTFITHDR}
	${CC} ${CFLAGS} -c $(subst .o,.c,$@)

potfit: ${POTFITOBJ} ${POTFITHDR}
	${CC} ${CFLAGS} -o ${BINDIR}/$@ ${POTFITOBJ} -lm

clean:
	rm -f *.o *~ 
neightab:
	${CC} ${CFLAGS} -o ${BINDIR}/neightab neightab.c -lm
