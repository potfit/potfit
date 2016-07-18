F95=gfortran
F = gfortran
AR = ar
FFLAGS = -O2 -m64 -fPIC

.SUFFIXES: .f90 .f

.f90.o:
	$(F95) $(FFLAGS) -c $< -o $*.o

.f.o:
	$(F) $(FFLAGS) -c $< -o $*.o

.f90.a:
	$(F95) $(FFLAGS) -Iinclude -c $< -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o

.f.a:
	$(F) $(FFLAGS) -Iinclude -c $< -o $*.o
	$(AR) $(ARFLAGS) $@ $*.o

OBJS = accept.o dpmpar.o fdavv.o fdjac.o lambda.o geodesiclm.o updatejac.o converge.o dgqt.o destsv.o

libgeodesiclm.a: $(OBJS)
	$(AR) -r $@ $(OBJS)

.PHONY: clean
clean:
	rm -f *.a *.o
