! -*- f90 -*-            
! ****************************************
! Routine for calculating finite-difference second directional derivative

SUBROUTINE FDAvv(m,n,x,v,fvec,fjac, func,acc, jac_uptodate, h2)
  IMPLICIT NONE
  INTEGER m, n
  REAL (KIND=8) x(n), v(n), fvec(m), fjac(m,n), acc(m), xtmp(n), ftmp(m), h2
  LOGICAL jac_uptodate
  EXTERNAL func

  IF( jac_uptodate) THEN
     xtmp = x + h2*v
     CALL func(m,n,xtmp,ftmp)
     acc = (2.0d+0/h2)*( (ftmp - fvec)/h2 - MATMUL(fjac,v) )
  ELSE !if jacobian not up to date, do not use jacobian in F.D. (needs one more function call)
     xtmp = x + h2*v
     CALL func(m,n,xtmp,ftmp)
     xtmp = x - h2*v
     CALL func(m,n,xtmp,acc)
     acc = (ftmp - 2*fvec + acc)/(h2*h2)
  ENDIF
END SUBROUTINE FDAvv
