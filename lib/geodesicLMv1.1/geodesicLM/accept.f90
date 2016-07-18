! -*- f90 -*-
! file accept.f90


SUBROUTINE Acceptance(n,C, Cnew, Cbest, ibold, accepted, &
     & dtd, v, vold)
  IMPLICIT NONE
  INTEGER n,accepted, ibold
  REAL (KIND=8) C, Cnew, Cbest, beta
  REAL (KIND=8) dtd(n,n), v(n), vold(n)

  IF( Cnew .LE. C) THEN !! Accept all downhill steps
     accepted = MAX(accepted + 1, 1)
  ELSE
     !! Calculate beta
     IF (DOT_PRODUCT(vold,vold) .EQ. 0.0d+0) THEN
        beta = 1.0d+0
     ELSE
        beta = DOT_PRODUCT(v,MATMUL(dtd, vold))
        beta = beta/SQRT( DOT_PRODUCT(v,MATMUL(dtd,v)) * DOT_PRODUCT(vold,MATMUL(dtd,vold) ))
        beta = min(1.0d+0,1.0d+0-beta)
     END IF
     SELECT CASE (ibold)
     CASE(0) !! Only downhill steps 
        IF( Cnew .LE. C) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted - 1, -1)
        END IF
     CASE(1)
        IF(beta*Cnew .LE. Cbest) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted-1,-1)
        END IF
     CASE(2)
        IF(beta*beta*Cnew .LE. Cbest) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted-1,-1)
        END IF
     CASE(3)
        IF(beta*Cnew .LE. C) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted-1,-1)
        END IF
     CASE(4)
        IF(beta*beta*Cnew .LE. C) THEN
           accepted = MAX(accepted + 1, 1)
        ELSE
           accepted = MIN(accepted-1,-1)
        END IF
     END SELECT
  END IF
END SUBROUTINE Acceptance
