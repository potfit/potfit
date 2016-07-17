! -*- f90 -*-            
! ****************************************
! Routine for calculating finite-difference jacobian

SUBROUTINE FDJAC(m,n,x,fvec,fjac,func,eps,center_diff)
  IMPLICIT NONE
  INTEGER m,n, i
  REAL (KIND=8) x(n), dx(n), fvec(m), fjac(m,n), eps, epsmach, dpmpar
  LOGICAL center_diff
  REAL (KIND=8) h, temp1(m), temp2(m) 
  
  epsmach = dpmpar(1)
  IF(center_diff) THEN
     DO i = 1, n
        h = eps*ABS(x(i))
        IF (h < epsmach) h = eps
        dx(:) = 0.0D+00
        dx(i) = 0.5d+0*h
        CALL func(m,n,x+dx,temp1,0)
        CALL func(m,n,x-dx,temp2,0)
        fjac(:,i) = (temp1 - temp2)/h
     END DO
  ELSE
     DO i = 1, n
        h = eps*ABS(x(i))
        IF (h < epsmach) h = eps
        dx(:) = 0.0D+00
        dx(i) = h
        CALL func(m,n,x+dx,temp1,0)
        fjac(:,i) = (temp1 - fvec)/h
     END DO
  END IF
 
END SUBROUTINE FDJAC
