CC     = gcc
CFLAGS = -g -O -DEAM -DNEWSCALE -DPOTSCALE
#(icc) -DPOTSCALE -DPARABEL
BINDIR = ${HOME}/bin/${HOSTTYPE}

AWKSCRIPTS = genpot.awk add_type.awk combine_eam.awk plotapot compnodes osix2potfit

#POTFITSRC = f1dim_r.c powell_lsq.c lubksb_r.c  mprove_r.c brent_r.c ludcmp_r.c linmin_r.c mnbrak_r.c nrutil_r.c force.c config.c param.c potential.c potfit.c splines.c simann.c mpi_utils.c
#POTFITHDR = potfit.h powell_lsq.h nrutil_r.h
#POTFITOBJ = $(POTFITSRC:.c=.o)

POTSCALESRC = param.c potential.c splines.c potscale.c #nrutil_r.c
POTSCALEHDR = potscale.h #nrutil_r.h
POTSCALEOBJ = $(POTSCALESRC:.c=.o)

.c.o: ${POTSCALEHDR}
	$(CC) $(CFLAGS) $(CINCLUDE) -c $<

.PHONY: force2imd vasp2force vasp2energy list_config potfit_setup makeapot

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

vasp2force:
	cp vasp2force ${HOME}/bin/share

vasp2energy:
	cp vasp2energy ${HOME}/bin/share

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
