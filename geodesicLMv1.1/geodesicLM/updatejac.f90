! -*- f90 -*-
!****************************************
! Routine for rank-deficient jacobian update

SUBROUTINE UPDATEJAC(m,n,fjac, fvec, fvec_new, acc, v, a)
  IMPLICIT NONE
  INTEGER m, n, i, j
  REAL (KIND=8) fjac(m,n), fvec(m), fvec_new(m), acc(m)
  REAL (KIND=8) v(n), a(n), djac(m), v2(n), r1(m)

  r1 = fvec + 0.5*MATMUL(fjac,v) + 0.125d+0*acc
  djac = 2.0*(r1 - fvec - 0.5*MATMUL(fjac,v))/DOT_PRODUCT(v,v)
  DO i = 1,m
     DO j = 1,n
        fjac(i,j) = fjac(i,j) + djac(i)*0.5d+0*v(j)
     END DO
  END DO
  v2 = 0.5d+0*(v + a)
  djac = 0.5*(fvec_new - r1 - MATMUL(fjac,v2))/DOT_PRODUCT(v2,v2)
  DO i = 1,m
     DO j = 1,n
        fjac(i,j) = fjac(i,j) + djac(i)*v2(j)
     END DO
  END DO
END SUBROUTINE UPDATEJAC

  
