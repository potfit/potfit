! -*- f90 -*-
!****************************************
! Routines for updating lam

SUBROUTINE TrustRegion(n,m, fvec, fjac, dtd, delta, lam)
  !! Calls dgqt supplied by minpack to calculate the step and Lagrange multiplier
  IMPLICIT NONE
  INTEGER n, m, i, itmax, info
  REAL (KIND=8) fvec(m), fjac(m,n), dtd(n,n), delta, lam, v(n)
  REAL (KIND=8) z(n), wa1(n), wa2(n) !! Work arrays for dgqt
  REAL (KIND=8) rtol, atol, f
  REAL (KIND=8) jtilde(m,n), gradCtilde(n), g(n,n)

  !! Parameters for dgqt
  rtol = 1.0d-03
  atol = 1.0d-03 
  itmax = 10

  DO i = 1,n
     jtilde(:,i) = fjac(:,i)/SQRT(dtd(i,i)) !! This assumes that dtd is diagonal...
  END DO
  gradCtilde = MATMUL(fvec, jtilde)
  g = MATMUL(TRANSPOSE(jtilde), jtilde)
  CALL dgqt(n, g, n, gradCtilde, delta, rtol, atol, itmax, lam, f, v, info, i, z, wa1, wa2)
  !! Transform v back to non-dtd units
!!$  DO i = 1,n
!!$     v(i) = v(i)/SQRT(dtd(i,i))
!!$  END DO
  RETURN
END SUBROUTINE TrustRegion
!! traditional update methods

SUBROUTINE Updatelam_factor(lam, accepted, factoraccept, factorreject)
  !! Update lam based on accepted/rejected step
  IMPLICIT NONE
  INTEGER accepted
  REAL (KIND=8) lam, factoraccept, factorreject

  IF(accepted.GE.0) THEN
     lam = lam / factoraccept
  ELSE
     lam = lam * factorreject
  ENDIF
END SUBROUTINE Updatelam_factor

SUBROUTINE Updatelam_nelson(lam, accepted, factoraccept, factorreject, rho)
  !! Update method due to Nelson [ref]
  IMPLICIT NONE
  INTEGER accepted, i
  DOUBLE PRECISION lam, factoraccept, factorreject, nu, rho
  IF(accepted.GE.0) THEN
     lam = lam * MAX( 1.0d+0/factoraccept, 1.0d+0 - (2.0d+0*(rho - 0.5d+0))**3 )
  ELSE
     nu = factorreject
     DO i = 1,-1*accepted !! double nu for each rejection
        nu = nu*2.0d+0
     END DO
     lam = lam * nu
  ENDIF
END SUBROUTINE Updatelam_nelson


SUBROUTINE Updatelam_Umrigar(m,n,lam, accepted, v, vold, fvec, fjac, dtd, a_param,C,Cnew)
  !! Method due to Umrigar and Nightingale [unpublished]
  IMPLICIT NONE
  INTEGER n,m,accepted, info
  DOUBLE PRECISION lam,v(n),vold(n),fvec(m), fjac(m,n), g(n,n), dtd(n,n), C, Cnew
  DOUBLE PRECISION lamold, d1, d2, grad(n), factor, a_param
  DOUBLE PRECISION amemory, cos_on

  amemory = EXP(-1.0d+0/5.0d+0)
  cos_on = DOT_PRODUCT(v,MATMUL(dtd, vold))
  cos_on = cos_on/SQRT( DOT_PRODUCT(v, MATMUL(dtd,v))*DOT_PRODUCT(vold, MATMUL(dtd,vold)))
  IF( accepted.GE.0) THEN

     IF( Cnew .LE. C ) THEN
        IF( cos_on .GT. 0) THEN
           a_param = amemory*a_param + 1.0 - amemory
        ELSE
           a_param =  amemory*a_param + 0.5*(1.0 - amemory)
        END IF
     ELSE
        a_param =  amemory*a_param + 0.5*(1.0 - amemory)
     END IF

     factor = MIN( 100.0d+0, MAX(1.1d+0, 1.0d+0/(2.2D-16 + 1.0d+0-ABS(2.0d+0*a_param - 1.0d+0))**2))
     IF (Cnew .LE. C .AND. cos_on .GE. 0) THEN
        lam = lam/factor
     ELSEIF(Cnew .GT. C) THEN
        lam = lam*SQRT(factor)
     END IF

  ELSE
     a_param =  amemory*a_param
     factor = MIN( 100.0d+0, MAX(1.1d+0, 1.0d+0/(2.2D-16 + 1.0d+0-ABS(2.0d+0*a_param - 1.0d+0))**2))
     lamold = lam
     IF( cos_on .GT. 0 ) THEN
        lam = lam * SQRT(factor)
     ELSE
        lam = lam * factor
     END IF
     ! Check for a 10% change in drift
     ! Umrigar and Nightingal suggest a check that the the proposed change in lam actually produces a meaningful change in the step.
     ! But this code produces strange results in a few cases.  -MKT
!!$     IF(accepted .EQ. -1) THEN
!!$        d1 = SQRT(DOT_PRODUCT(v,v))
!!$        g = MATMUL(TRANSPOSE(fjac),fjac) + lam*dtd
!!$        CALL DPOTRF('U', n, g, n, info)
!!$        grad = MATMUL(fvec, fjac)
!!$        CALL DPOTRS('U', n, 1, g, n, grad, n, info)
!!$        d2 = SQRT(DOT_PRODUCT(grad, grad) )
!!$        IF( 10.0d+0*ABS(d2-d1) .LT. d2 ) lam = lam - 0.1*d2*(lamold - lam)/(d1 - d2)
!!$     END IF
  ENDIF
END SUBROUTINE Updatelam_Umrigar


!! Trust region update methods

SUBROUTINE Updatedelta_factor(delta, accepted, factoraccept, factorreject)
  !! Update lam based on accepted/rejected step
  IMPLICIT NONE
  INTEGER accepted
  REAL (KIND=8) delta, factoraccept, factorreject

  IF(accepted.GE.0) THEN
     delta = delta * factoraccept
  ELSE
     delta = delta / factorreject
  ENDIF
END SUBROUTINE Updatedelta_factor

SUBROUTINE Updatedelta_more(delta, lam, n, x, dtd, rho, C, Cnew, dirder, actred, av, avmax)
  IMPLICIT NONE
  INTEGER n
  DOUBLE PRECISION delta, lam, x(n), dtd(n,n), rho, C, Cnew, dirder, actred, av, avmax
  DOUBLE PRECISION pnorm, temp
  pnorm = SQRT(DOT_PRODUCT(x,MATMUL(dtd, x)))
  IF (rho .GT. 0.25d+0) THEN
     IF (lam .GT. 0.0d+0 .AND. rho .LT. 0.75d+0) THEN
        temp = 1.0d+0
     ELSE
        temp = 2.0d+0*pnorm/delta
     END IF
  ELSE
     IF ( actred .GE. 0.0d+0) THEN
        temp = 0.5d+0
     ELSE
        temp = 0.5d+0*dirder/(dirder + 0.5d+0*actred)
     END IF
     IF ( 0.01*Cnew .GE. C .OR. temp .LT. 0.1d+0) temp = 0.1d+0
  END IF
  !! We need to make sure that if acceleration is too big, we decrease teh step size
  IF (av .GT. avmax) THEN
     temp = MIN(temp,MAX(avmax/av,0.1d+0))
  END IF

  delta = temp*MIN(delta,10.0d+0*pnorm)
  lam = lam/temp
END SUBROUTINE Updatedelta_more


