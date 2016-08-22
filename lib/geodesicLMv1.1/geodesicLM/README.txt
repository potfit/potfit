This it the README for the geodesic Levenberg-Marquardt algorithm v1.0

Geodesic Levenberg-Marquardt is a variant of Levenberg-Marquardt that adds third order corrections to the proposed step from the directional second derivative (either by a finite difference estimate or an analytic evaluation).  The main routine is in the file geodesiclm.f90

The method makes use of the BLAS and LAPACK subroutines for matrix manipulation.  You should link to these libraries when compiling.

Although we have used this routine successfully in our own research, we do not guarantee that it is bug free.  If you encounter a bug (or an unexpected behavior) please let us know.  Send details about the bug to Mark Transtrum: mktranstrum@byu.edu

If you use this code, please acknowledge such by referencing one one of the following papers in any published work:
    
Transtrum M.K., Machta B.B., and Sethna J.P, Why are nonlinear fits to data so challenging?  Phys. Rev. Lett. 104, 060201 (2010)

Transtrum M.K., Machta B.B., and Sethna J.P., The geometry of nonlinear least squares with applications to sloppy model and optimization.  Phys. Rev. E. 80, 036701 (2011)

Disclaimer: No guarantee whatsoever is provided. No liability whatsoever is accepted for any loss or damage of any kind resulting from any defect or inaccuracy in this code. 
