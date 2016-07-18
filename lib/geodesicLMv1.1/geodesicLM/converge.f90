! -*- f90 -*-
! ****************************************
! Routine to Check for Convergence

SUBROUTINE convergence_check(m, n, converged, accepted, counter, C, Cnew, x, fvec, fjac, lam, xnew, &
     & nfev, maxfev, njev, maxjev, naev, maxaev, maxlam, minlam, artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,cos_alpha)
  IMPLICIT NONE
  INTEGER m,n, converged, accepted, counter, nfev, maxfev, njev, maxjev, naev, maxaev
  REAL (KIND=8) C, Cnew, x(n), fvec(m), fjac(m,n), xnew(n), grad(n), lam, rpar(m)
  REAL (KIND=8) maxlam, minlam, artol, Cgoal, gtol, xtol, xrtol, ftol, frtol, cos_alpha
  INTEGER i

!  The first few criteria should be checked every iteration, since
!  they depend on counts and the Jacobian but not the proposed step.

!  nfev
  IF(maxfev .GT. 0) THEN
     IF(nfev.GE.maxfev) THEN
        converged = -2
        counter = 0
        RETURN
     ENDIF
  ENDIF

!  njev
  IF(maxjev .GT. 0) THEN
     IF(njev.GE.maxjev) THEN
        converged = -3
        RETURN
     ENDIF
  ENDIF


!  naev
  IF(maxaev .GT. 0) THEN
     IF(naev.GE.maxaev) THEN
        converged = -4
        RETURN
     END IF
  ENDIF

!  maxlam
  IF(maxlam .GT. 0.0d+0) THEN
     IF(lam .GE. maxlam) THEN
        converged = -5
        RETURN
     END IF
  END IF

!  minlam
  IF(minlam .GT. 0.0d+0 .AND. lam .GT. 0.0d+0) THEN
     IF(lam .LE. minlam) THEN
        counter = counter + 1
        IF(counter .GE. 3) THEN
           converged = -6
           RETURN
        END IF
        RETURN
     END IF
  END IF

! artol -- angle between residual vector and tangent plane
  IF( artol .GT. 0.0d+0) THEN
     !! Only calculate the projection if artol > 0
     !! CALL projection(m,n,fvec, fjac, rpar,eps)
     !! cos_alpha = SQRT(DOT_PRODUCT(rpar,rpar)/DOT_PRODUCT(fvec,fvec))
     IF( cos_alpha .LE. artol) THEN
        converged = 1
        RETURN
     END IF
  END IF

! If gradient is small
  grad = -1.0d+0*MATMUL(fvec, fjac)
  IF(SQRT(DOT_PRODUCT(grad,grad)).LE.gtol) THEN
     converged = 3
     RETURN
  ENDIF

! If cost is sufficiently small
  IF (C .LT. Cgoal) THEN !! Check every iteration in order to catch a cost small on the first iteration
     converged = 2
     RETURN
  END IF


!  If step is not accepted, then don't check remaining criteria
  IF(accepted.LT.0) THEN
     counter = 0
     converged = 0
     RETURN
  ENDIF

! If step size is small
  if(SQRT(DOT_PRODUCT(x-xnew,x-xnew)).LT. xtol) THEN
     converged = 4
     RETURN
  ENDIF

! If each parameter is moving relatively small
  xrtolcheck: DO i = 1,n
     converged = 5
     IF(  ABS(x(i) - xnew(i)) .GT. xrtol*ABS(x(i)) .OR.  (xnew(i) .NE. xnew(i)) ) converged = 0 !! continue if big step or nan in xnew
     IF( converged .EQ. 0) EXIT xrtolcheck
  END DO xrtolcheck
  IF(converged .EQ. 5) RETURN

! If cost is not decreasing -- this can happen by accident, so we require that it occur three times in a row
  IF( (C - Cnew).LE. ftol .AND.(C-Cnew).GE.0.) THEN
     counter = counter + 1
     IF(counter .GE. 3) THEN
        converged = 6
        RETURN
     END IF
     RETURN
  ENDIF

! If cost is not decreasing relatively -- again can happen by accident so require three times in a row
  IF( (C - Cnew).LE.(frtol*C).AND.(C-Cnew).GE.0.) THEN
     counter = counter + 1
     IF(counter .GE. 3) THEN
        converged = 7
        RETURN
     ENDIF
     RETURN
  ENDIF

! If none of the above: continue
  counter = 0
  converged = 0
END SUBROUTINE convergence_check

