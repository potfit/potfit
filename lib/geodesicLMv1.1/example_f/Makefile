F = gfortran
FFLAG = -O2 -m64 -fPIC

LIB = -L../geodesicLM -lgeodesiclm

LAPACKLIBS = -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
             -Wl,--end-group -lpthread -Wl,--as-needed

LIB += ${LAPACKLIBS} 

OBJECTS = demo.o

demo: ${OBJECTS} 
	${F} ${FFLAG} -o demo ${OBJECTS} ${LIB} 


demo.o: demo.f90
	${F} ${FFLAG} -c demo.f90

.PHONY: clean
clean:
	rm *.o demo
