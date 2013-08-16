CC     = gcc
CFLAGS = -g -O -DEAM -DNEWSCALE -DPOTSCALE
BINDIR = ${HOME}/bin/${HOSTTYPE}

AWKSCRIPTS = genpot.awk add_type.awk combine_eam.awk plotapot compnodes.awk osix2potfit

POTSCALESRC = param.c potential.c splines.c potscale.c
POTSCALEHDR = potscale.h
POTSCALEOBJ = $(POTSCALESRC:.c=.o)

.c.o: ${POTSCALEHDR}
	$(CC) $(CFLAGS) $(CINCLUDE) -c $<

.PHONY: force2imd vasp2force list_config potfit_setup makeapot

clean:
	rm -f *.o *~

pottrans: pottrans.c
	${CC} ${CFLAGS} -o ${BINDIR}/pottrans pottrans.c -lm

force2poscar: force2poscar.c
	${CC} ${CFLAGS} -o ${BINDIR}/force2poscar force2poscar.c

potscale: ${POTSCALEOBJ} ${POTSCALEHDR}
	${CC} ${CFLAGS} -o ${BINDIR}/$@ ${POTSCALEOBJ} -lm

vasp2force:
	cp vasp2force ${HOME}/bin/share

force2imd:
	cp force2imd ${HOME}/bin/share

list_config:
	cp list_config ${HOME}/bin/share

potfit_setup:
	cp potfit_setup ${HOME}/bin/share

makeapot:
	cp makeapot ${HOME}/bin/share

awkscripts:
	cd awk;	cp -pr ${AWKSCRIPTS} ${HOME}/bin/share;	cd -
