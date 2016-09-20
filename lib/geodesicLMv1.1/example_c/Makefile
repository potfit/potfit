FC = gfortran
CC = gcc
CFLAG = -O2 -g -Wall

LIB = -L../geodesicLM -lgeodesiclm -lm -lgfortran
LAPACKLIBS = -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
             -Wl,--end-group -lpthread -Wl,--as-needed
LIB +=${LAPACKLIBS} 

OBJECTS = demo.o

demo: ${OBJECTS} 
	${CC} ${CFLAG} -o demo ${OBJECTS} ${LIB} 

demo.o: demo.c
	${CC} ${CFLAG} -c demo.c

.PHONY: clean
clean:
	rm *.o demo
