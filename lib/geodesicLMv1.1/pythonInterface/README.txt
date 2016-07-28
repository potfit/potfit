README.txt for the python wrapper to the geodesic Levenberg-Marquardt routine.

To begin, compile the geodesic FORTRAN LM routine and archive the resulting files to libgeodesiclm.a  You should be able to do this with the included makefile in the geodesiclm package.  You will also need to have the BLAS and LAPACK libraries to link to.

To compile the python wrapper run:

f2py -c geodesiclm.pyf -Lpath_to_libgeodesiclm.a -lgeodesiclm  -Lpath_to_lapack -llapack -Lpath_to_blas -lblas

where path_to_libgeolevmar.a should be replaced with the path to the file libgeolevmar.a.  Similarly for path_to_lapack and path_to_blas.  For example, on my computer, I put all these files in C:\lib so that I run

python C:\Python27\Scripts\f2py.py -c geodesiclm.pyf -LC:\lib -lgeodesiclm  -LC:\lib -llapack -LC:\lib -lblas

This generates a file (_geodesiclm.pyd on windows or _geodesiclm.so on unix).  Place this file in the same folder as the included geodesiclm.py and you can call the python routine geodesiclm.

Sometimes, it is also necessary to link to the gfortran library:

f2py -c leastsq.pyf -Lpath_to_libgeodesiclm.a -lgeodesiclm  -Lpath_to_lapack -llapack -Lpath_to_blas -lblas -lgfortran

numpy often has its own library of lapack.  To see what is available type
f2py --help-link lapack
