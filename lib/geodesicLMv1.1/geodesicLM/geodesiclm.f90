! -*- f90 -*-
! file leastsq.f90
! Main Geodesic-Bold-BroydenUpdate-Levenberg-Marquardt routine
! version 1.1

SUBROUTINE geodesiclm(func, jacobian, Avv, &
     & x, fvec, fjac, n, m, &
     & callback, info, &
     & analytic_jac, analytic_Avv, &
     & center_diff, h1, h2,&
     & dtd, damp_mode, &
     & niters, nfev, njev, naev, &
     & maxiter, maxfev, maxjev, maxaev, maxlam, minlam, &
     & artol, Cgoal, gtol, xtol, xrtol, ftol, frtol, &
     & converged, &
     & print_level, print_unit, &
     & imethod, iaccel, ibold, ibroyden, &
     & initialfactor, factoraccept, factorreject, avmax)

!*****************************************************************
! 
!    subroutine geodesicLM
!    
!    The purpose of geolevmar is to minimize the sum of the squares
!    of m nonlinear functions of n variables by a modification of
!    the Levenberg-Marquardt algorithm that utilizes the geodesic
!    acceleration step correction, bold acceptance criterion, and
!    a Broyden update of the jacobian matrix.  The method employs one
!    of several possible schemes for updating the Levenberg-Marquardt
!    parameter.  The user must provide a subroutine which calcualtes
!    the functions, and optionally the jacobian and a directional 
!    derivative of the functions.  The latter two will be estimated
!    by finite differences if not supplied.
!
!    If you use this code, please acknowledge such by referencing one
!    one of the following papers in any published work:
!    
!    Transtrum M.K., Machta B.B., and Sethna J.P, Why are nonlinear
!    fits to data so challenging?  Phys. Rev. Lett. 104, 060201 (2010)
!
!    Transtrum M.K., Machta B.B., and Sethna J.P., The geometry of
!    nonlinear least squares with applications to sloppy model and
!    optimization.  Phys. Rev. E. 80, 036701 (2011)
!
!
!    The subroutine statement is:
!
!    geodesicLM(func, jacobian, Avv, x, fvec, fjac, n, m, callback, info,
!              analytic_jac, analytic_Avv, center_diff, h1, h2,
!              dtd, damp_mode, niteres, nfev, njev, naev,
!              maxiters, maxfev, maxjev, maxaev, maxlam, minlam,
!              artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,
!              converged, print_level, print_unit,
!              imethod, iaccel, ibold, ibroyden,
!              initialfactor, factoraccept, factorreject, avmax)
!
!    where
!
!    func is a user supplied subroutine which calculates the functions and
!    should be written as follows:
!
!      subroutine func(m, n, x, fvec)
!      integer m, n
!      double precision x(n), fvec(m)
!      --------------------------------------------------------------------
!      calculates the function at x and returns their values in fvec
!      x, m, and n should be left unchanged
!      --------------------------------------------------------------------
!      end subroutine func
!
!    jacobian is a user supplied subroutine which calculates the jacobian of
!    of the functions if analytic_jac is .TRUE.
!    jacobian should be writen as follows     
!
!      subroutine jacobian(m, n, x, fjac)
!      integer m, n
!      double precision x(n), fjac(m,n)
!      --------------------------------------------------------------------
!      calculates the jacobian at x and returns their values in fjac
!      x, m, and n should be left unchanged
!      --------------------------------------------------------------------
!      end subroutine jacobian
!
!    Avv is a user supplied subroutine which calculates the directional
!    second derivative of the functions if analytic_Avv is .TRUE.
!    Avv should be writen as follows     
!
!      subroutine Avv(m, n, x, v, acc)
!      integer m, n
!      double precision x(n), v(n), acc(m)
!      --------------------------------------------------------------------
!      calculates the directional second derivative at x in the direction 
!      of v and returns the values in acc
!      x, v, m, and n should be left unchanged
!      --------------------------------------------------------------------
!      end subroutine Avv
!
!    x is an array of length n.  On input it contains an initial estimate of
!    the solution.  On exit, it contains the final estimate of the solution.
!
!    fvec is an output array of length m containing the funtion evaluation at
!    the final solution
!
!    fjac is an output array of dimension(m,n) containing the jacobian evaluation
!    the final solution.  The array MATMUL( TRANSPOSE(fjac), fjac) is an estimate
!    of the covariance matrix of parameters at the final solution
!
!    n an input integer set to the number of parameters
!
!    m an input integer set to the number of functions
!
!    callback a user supplied subroutine to be called after each iteration of the
!    algorithm.  
!    callback should be written as follows
!
!      subroutine callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info)
!      integer m, n, accepted, info
!      double precision x(n), v(n), a(n), fvec(m), fjac(m,n), acc(m), lam, dtd(n,n), fvec_new(m)
!      --------------------------------------------------------------------
!      m, n, x, v, a, fvec, fjac, acc, lam, dtd, fvec_new, accepted, should be left unchanged
!      On input, info = 0 and should be changed to a nonzero value if the user
!      wishes to terminate calculation
!      --------------------------------------------------------------------
!      end subroutine callback
!
!    info an output integer set to a nonzero value if the user terminated the routine
!    (see callback).
!
!    analytic_jac an input boolean set to .TRUE. if the subroutine jacobian calculates
!    the jacobian.  If .FALSE. then a finite difference estimate will be used.
!
!    analytic_Avv an input boolean set to .TRUE. if the subroutine Avv calculates
!    the directional second derivative.  If .FALSE. then a finite difference estimate
!    will be used.
!
!    center_diff an input boolean.  If finite differences are used to estimate the jacobian
!    then center differences will used if center_diff is .TRUE., otherwise, forward
!    differences will be used.  Note that center differences are more accurate by require
!    more function evaluations.
!
!    h1 an input double precision specifying the step size for the finite difference estimates
!    of the jacobian.
!
!    h2 an input double precision specifying the steps ize for the finite difference estiamtes
!    of the directional second derivative.
!
!    dtd a double precision array of dimension(n,n).  dtd is used as the damping matrix in the 
!    Levenberg-Marquardt routine.  It's exact treatment is specified by the damp_mode input.
!
!    damp_mode an input integer specifying the details of the LM damping as follows:
!      damp_mode = 0: dtd is set to the identity.
!      damp_mode = 1: dtd should be a positive definite, diagonal matrix whose entries are dynamically
!                updated based on the elements of the jacobian.
!
!    niters an output integer specifying the number of iterations of the algorithm.
!
!    nfev an output integer specifying the number of calls to func.  
!
!    njev an output integer specifying the number of calls to jacobian.
!
!    naev an output integer specifying the number of calls to Avv.
!
!    maxiter an input integer specifying the maximum number of routine iterations.
!
!    maxfev an input integer specifying the maximum number of function calls
!    if maxfev = 0, then there is no limit to the number of function calls.
!
!    maxjev an input integer specifying the maximum number of jacobian calls
!    if maxjev = 0, then there is no limit to the number of jacobian calls.
!
!    maxaev an input integer specifying the maximum number of Avv calls
!    if maxaev = 0, then there is no limit to the number of Avv calls.
!
!    maxlam an input double precision specifying the maximum allowed value of 
!    the damping term lambda. If this is negative, then there is no limit.
!
!    minlam an input double precision specifying the minimum allowed value of 
!    the damping term lambda. If lambda is smaller than this value for three consecutive steps
!    the routine terminates.  If this is negative, then there is no limit.
!
!    artol an input double precision.  The method will terminate when the cosine of the
!    angle between the residual vector and the range of the jacobian is less than artol.
!
!    Cgoal an input double precision.  The method will terminate when the cost (one half
!    the sum of squares of the function) falls below Cgoal.
!
!    gtol an input double precision.  The method will terminate when norm of Cost gradient 
!    falls below gtol.
!    
!    xtol an input double precision.  The method will terminate when parameters change by
!    less than xtol.
!
!    xrtol an input double precision.  The method will terminate if the relative change in
!    each of the parameters is less than xrtol.
!
!    ftol an input double precision.  The method will termiante if the Cost fails to decrease
!    by more than ftol for 3 consecutive iterations.
!
!    frtol an input double precision.  The method will terminate if the relative decrease in
!    Cost is less than frtol 3 consecutive iterations.
!
!    converged an output integer indicated the reason for termination:
!      converged = 1: artol 
!      converged = 2: Cgoal
!      converged = 3: gtol
!      converged = 4: xtol
!      converged = 5: xrtol
!      converged = 6: ftol
!      converged = 7: frtol
!      converged = -1: maxiters exeeded
!      converged = -2: maxfev exceeded
!      converged = -3: maxjev exceeded
!      converged = -4: maxaev exceeded
!      converged = -10: user requested termination in callback via info
!      converged = -11: Either the initial function evalaution or subsequent jacobian
!                       evaluations produced Nans.
!
!    print_level an input integer specifying the amount of details to be printed.
!    acceptable values range from 0 to 5, with larger number printing more details.
!
!    print_unit an input integer specifying the unit number details should be written to.
!
!    imethod an input integer specifying the method for updating the LM parameter
!      imethod = 0: adjusted by fixed factors after accepted/rejected steps
!      imethod = 1: adjusted as described in Nielson
!      imethod = 2: adjusted according to an unpublished method due to Cyrus Umrigar and Peter Nightingal
!      imethod = 10: step size Delta adjusted by fixed factors after accepted/rejected steps
!      imethod = 11: step size adjusted as described in More'
!
!    initialfactor an input double precision for specifying either the initial LM parameter
!    of the initial step size.
!
!    factoraccept an input double precision (larger than 1.0) specifying the factor by which
!    either the LM parameter or the step size will be adjusted after an accepted step if
!    imethod = 0 or 10
!
!    factorreject an input double precision (larger than 1.0) specifying the factor by which
!    either the LM parameter of the step size will be adjusted after a rejected step if
!    imethod = 0 or 10
!
!    avmax an input double precision specifying the maximum norm of the geodesic acceleration 
!    relative to the velocity vector.
!
!*****************************************************************

  IMPLICIT NONE
  !! Passed parameters
  EXTERNAL func, jacobian, Avv, callback

  REAL (KIND=8) x(n), fvec(m), fjac(m,n)
  INTEGER n, m
  LOGICAL analytic_jac, analytic_Avv, center_diff
  REAL (KIND=8) h1, h2
  REAL (KIND=8) dtd(n,n)
  INTEGER damp_mode, info
  INTEGER niters, nfev, njev, naev
  INTEGER maxiter, maxfev, maxaev, maxjev
  REAL (KIND=8) maxlam, minlam, artol, Cgoal, gtol, xtol, xrtol, ftol, frtol
  INTEGER converged
  INTEGER print_level, print_unit
  INTEGER iaccel, ibroyden, ibold, imethod
  REAL (KIND=8) avmax, initialfactor, factoraccept, factorreject

  !! Internal parameters

  REAL (KIND=8) acc(m), v(n), vold(n), a(n), lam, delta, cos_alpha, av
  REAL (KIND=8) fvec_new(m), fvec_best(m), C, Cnew, Cbest, Cold
  REAL (KIND=8) x_new(n), x_best(n)
  REAL (KIND=8) jtj(n,n), g(n,n)
  REAL (KIND=8) temp1, temp2, pred_red, dirder, actred, rho, a_param
  INTEGER i, j, istep, accepted, counter
  
  character(16) :: converged_info(-11:7)
  LOGICAL jac_uptodate, jac_force_update, valid_result

  ! strings for concluding print statement
  converged_info = '????????'
  converged_info(1) = 'artol reached'
  converged_info(2) = 'Cgoal reached'
  converged_info(3) = 'gtol reached'
  converged_info(4) = 'xtol reached'
  converged_info(5) = 'xrtol reached'
  converged_info(6) = 'ftol reached'
  converged_info(7) = 'frtol reached'
  converged_info(-1) = 'maxiters exeeded'
  converged_info(-2) = 'maxfev exceeded'
  converged_info(-3) = 'maxjev exceeded'
  converged_info(-4) = 'maxaev exceeded'
  converged_info(-10) = 'User Termination '
  converged_info(-11) = 'NaN Produced'

  IF(print_level .GE. 1) THEN
     WRITE(print_unit, *) "Optimizing with Geodesic-Levenberg-Marquardt algorithm, version 1.1"
     WRITE(print_unit, *) "Method Details:"
     WRITE(print_unit, *) "  Update method:   ", imethod
     WRITE(print_unit, *) "  acceleration:    ", iaccel
     WRITE(print_unit, *) "  Bold method:     ", ibold
     WRITE(print_unit, *) "  Broyden updates: ", ibroyden
     FLUSH(print_unit)
  ENDIF

  !! Initialize variables
  niters = 0
  nfev = 0
  naev = 0
  njev = 0
  converged = 0
  v(:) = 0.0d+0
  vold(:) = 0.0d+0
  a(:) = 0.0d+0
  cos_alpha = 1.0d+0
  av = 0.0d+0
  a_param = 0.5

  accepted = 0
  counter = 0
  CALL func(m,n,x,fvec)
  nfev = nfev + 1
  C = 0.5d+0*DOT_PRODUCT(fvec,fvec)
  IF(print_level .GE. 1) THEN
     !! times 2 to be compatible with potfit definition of Cost
     !!WRITE(print_unit, *) "  Initial Cost:    ", C
     WRITE(print_unit, *) "  Initial Cost:    ", 2.0*C
     FLUSH(print_unit)
  ENDIF
  valid_result = .TRUE.
  !! Check for nans in fvec
  checkfvec: DO i = 1,m     
     IF(fvec(i) /= fvec(i)) THEN
        valid_result = .FALSE.
        EXIT checkfvec
     END IF
  END DO checkfvec
  IF (.NOT. valid_result) THEN
     converged = -11
     maxiter = 0
  ENDIF
  Cbest = C
  fvec_best = fvec
  x_best = x
  IF(analytic_jac) THEN
     CALL jacobian(m,n,x,fjac)
     njev = njev + 1
  ELSE 
     CALL fdjac(m,n,x,fvec,fjac,func,h1,center_diff)
     IF (center_diff) THEN
        nfev = nfev + 2*n
     ELSE
        nfev = nfev + n
     ENDIF
  ENDIF
  jac_uptodate = .TRUE.
  jac_force_update = .FALSE.
  jtj = MATMUL(TRANSPOSE(fjac), fjac)

  !! Check fjac for nans
  valid_result = .TRUE.
  checkfjac_initial: DO i = 1,m     
     DO j = 1,n
        IF(fjac(i,j) /= fjac(i,j)) THEN
           valid_result = .FALSE.
           EXIT checkfjac_initial
        END IF
     END DO
  END DO checkfjac_initial
  IF( .NOT. valid_result) THEN
     converged = -11
     maxiter = 0
  ENDIF

  acc(:) = 0.0d+0
  a(:) = 0.0d+0

  !! Initialize scaling matrix
  IF(damp_mode.EQ.0) THEN
     dtd(:,:) = 0.0d+0
     DO i=1,n
        dtd(i,i) = 1.0d+0
     END DO
  ELSEIF(damp_mode.EQ.1) THEN
     DO i = 1,n
        dtd(i,i) = MAX(jtj(i,i),dtd(i,i))
     END DO
  ENDIF

  !! Initialize lambda
  IF(imethod .LT. 10) THEN
     lam = jtj(1,1)
     DO i = 2,n
        lam = MAX(jtj(i,i),lam)
     END DO
     lam = lam * initialfactor
  !! Initialize step bound if using trust region method
  ELSEIF(imethod .GE. 10) THEN
     delta = initialfactor*SQRT(DOT_PRODUCT(x,MATMUL(dtd,x)))
     lam = 1.0d+0
     IF(delta .EQ. 0.0d+0) delta = 100d+0
     IF( converged .EQ. 0) CALL TrustRegion(n,m,fvec,fjac,dtd,delta,lam) !! Do not call this if there were nans in either fvec or fjac
  ENDIF

  !! Main Loop
  main: DO istep=1, maxiter
     
     info = 0
     CALL callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info)
     IF( info .NE. 0) THEN
        converged = -10
        exit main
     ENDIF
     !! Update Functions

     !! Full or partial Jacobian Update?
     IF (accepted .GT. 0 .AND. ibroyden .LE. 0) jac_force_update = .TRUE.
     IF (accepted + ibroyden .LE. 0 .AND. .NOT. jac_uptodate) jac_force_update = .TRUE.  !Force jac update after too many failed attempts

     IF (accepted .GT. 0 .AND. ibroyden .GT. 0 .AND. .NOT. jac_force_update) THEN !! Rank deficient update of jacobian matrix
        CALL UPDATEJAC(m,n,fjac, fvec, fvec_new, acc, v, a)
        jac_uptodate = .FALSE.
     ENDIF

     IF( accepted .GT. 0) THEN !! Accepted step
        fvec = fvec_new
        x = x_new
        vold = v
        C = Cnew
        IF( C .LE. Cbest) THEN
           x_best = x
           Cbest = C
           fvec_best = fvec
        ENDIF
     ENDIF
     

     IF( jac_force_update ) THEN !! Full rank update of jacobian
        IF(analytic_jac) THEN
           CALL jacobian(m,n,x,fjac)
           njev = njev + 1
        ELSE
           CALL fdjac(m,n,x,fvec,fjac,func,h1,center_diff)
           IF (center_diff) THEN
              nfev = nfev + 2*n
           ELSE
              nfev = nfev + n
           ENDIF
        ENDIF
        jac_uptodate = .TRUE.
        jac_force_update = .FALSE.
     ENDIF

     !! Check fjac for nans
     valid_result = .TRUE.
     checkfjac: DO i = 1,m     
        DO j = 1,n
           IF(fjac(i,j) /= fjac(i,j)) THEN
              valid_result = .FALSE.
              EXIT checkfjac
           END IF
        END DO
     END DO checkfjac

     IF (valid_result) THEN !! If no nans in jacobian

        jtj = MATMUL(TRANSPOSE(fjac), fjac)
     
        !! Update Scaling/lam/TrustRegion
        IF(istep .GT. 1) THEN !! Only necessary after first step
           IF( damp_mode .EQ. 1) THEN
              DO i = 1,n
                 dtd(i,i) = MAX(jtj(i,i),dtd(i,i))
              END DO
           ENDIF
           !! Can add other lam-delta update methods 
           SELECT CASE( imethod )
           CASE( 0 )
              !! Update lam directly by fixed factors
              CALL Updatelam_factor(lam, accepted, factoraccept, factorreject)
           CASE( 1 )
              !! Update lam directly based on Gain Factor rho (see Nielson reference)
              CALL Updatelam_nelson(lam, accepted, factoraccept, factorreject, rho)
           CASE( 2 )
              !! Update lam directly using method of Umrigar and Nightingale [unpublished]
              CALL Updatelam_Umrigar(m,n,lam,accepted, v, vold, fvec, fjac, dtd, a_param, Cold, Cnew) 
           CASE(10)
              !! Update delta by fixed factors
              CALL UpdateDelta_factor(delta, accepted, factoraccept, factorreject)
              CALL TrustRegion(n,m,fvec, fjac, dtd, delta, lam)
           CASE(11)
              !! Update delta as described in More' reference
              CALL UpdateDelta_more(delta, lam, n, x, dtd, rho, C, Cnew, dirder, actred, av, avmax)
              CALL TrustRegion(n,m,fvec, fjac, dtd, delta, lam)
           END SELECT
        ENDIF

        !! Propose Step
        !! metric aray
        g = jtj + lam*dtd
        !! Cholesky decomposition
        CALL DPOTRF('U', n, g, n, info)
        !! CALL inv(n, g, info)
     ELSE !! If nans in jacobian
        converged = -11
        exit main
     ENDIF


     IF(info .EQ. 0) THEN  !! If matrix decomposition successful:
        !! v = -1.0d+0*MATMUL(g,MATMUL(fvec,fjac)) ! velocity
        v = -1.0d+0*MATMUL(fvec, fjac)
        CALL DPOTRS('U', n, 1, g, n, v, n, info)
        
        ! Calcualte the predicted reduction and the directional derivative -- useful for updating lam methods
        temp1 = 0.5d+0*DOT_PRODUCT(v,MATMUL(jtj, v))/C
        temp2 = 0.5d+0*lam*DOT_PRODUCT(v,MATMUL(dtd,v))/C
        pred_red = temp1 + 2.0d+0*temp2
        dirder = -1.0d+0*(temp1 + temp2)
        ! calculate cos_alpha -- cos of angle between step direction (in data space) and residual vector
        cos_alpha = ABS(DOT_PRODUCT(fvec, MATMUL(fjac, v)))
        cos_alpha = cos_alpha/SQRT( DOT_PRODUCT(fvec, fvec)*DOT_PRODUCT(MATMUL(fjac,v), MATMUL(fjac, v)) )
        IF ( imethod .LT. 10) delta = SQRT(DOT_PRODUCT(v, MATMUL(dtd, v)))  !! Update delta if not set directly
        !! update acceleration
        IF(iaccel .GT. 0) THEN
           IF( analytic_Avv ) THEN 
              CALL Avv(m,n,x,v,acc)
              naev = naev + 1
           ELSE
              CALL FDAvv(m,n,x,v,fvec, fjac, func, acc, jac_uptodate, h2)
              IF(jac_uptodate) THEN
                 nfev = nfev + 1
              ELSE 
                 nfev = nfev + 2 !! we don't use the jacobian if it is not up to date
              ENDIF
           ENDIF
           !! Check accel for nans
           valid_result = .TRUE.
           checkAccel: DO i = 1,m     
              IF(acc(i) /= acc(i)) THEN
                 valid_result = .FALSE.
                 EXIT checkAccel
              END IF
           END DO checkAccel
           IF (valid_result ) THEN
              a = -1.0d+0*MATMUL(acc, fjac)
              CALL DPOTRS('U', n, 1, g, n, a, n, info)
              !!a = -1.0d+0*MATMUL(g,MATMUL(acc,fjac))
           ELSE 
              a(:) = 0.0d+0 !! If nans in acc, we will ignore the acceleration term
           ENDIF
        ENDIF

        !! Evaluate at proposed step -- only necessary if av <= avmax
        av = SQRT(DOT_PRODUCT(a,MATMUL(dtd, a))/DOT_PRODUCT(v,MATMUL(dtd,v)))
        IF( av .LE. avmax) THEN
           x_new = x + v + 0.5d+0*a
           CALL func(m,n,x_new,fvec_new)
           nfev = nfev + 1
           Cnew = 0.5d+0*DOT_PRODUCT(fvec_new,fvec_new)
           Cold = C
           valid_result = .TRUE.
           !! Check for nans in fvec_new
           checkfvec_new: DO i = 1,m     
              IF(fvec_new(i) /= fvec_new(i)) THEN
                 valid_result = .FALSE.
                 EXIT checkfvec_new
              END IF
           END DO checkfvec_new
           IF (valid_result) THEN  !! If no nans, proceed as normal
              ! update rho and actred
              actred = 1.0d+0 - Cnew/C
              rho = 0.0d+0
              IF(pred_red .NE. 0.0d+0) rho = (1.0d+0 - Cnew/C)/pred_red
              !! Accept or Reject proposed step
              CALL Acceptance(n,C, Cnew, Cbest, ibold, accepted, dtd, v, vold)
           ELSE !! If nans in fvec_new, reject step
              actred = 0.0d+0
              rho = 0.0d+0
              accepted = MIN(accepted - 1, -1)
           ENDIF
        ELSE !! If acceleration too large, then reject
           accepted = MIN(accepted - 1, -1)
        ENDIF
     ELSE !! If matrix factorization fails, reject the proposed step
        accepted = MIN(accepted -1, -1)
     ENDIF

     !! Check Convergence
     IF (converged .EQ. 0) THEN
        CALL convergence_check(m, n, converged, accepted, counter, &
             & C, Cnew, x, fvec, fjac, lam, x_new, &
             & nfev, maxfev, njev, maxjev, naev, maxaev, maxlam, minlam, &
             & artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,cos_alpha)
        IF (converged .EQ. 1 .AND. .NOT. jac_uptodate) THEN  
           !! If converged by artol with an out of date jacobian, update the jacoban to confirm true convergence
           converged = 0
           jac_force_update = .TRUE.
        END IF
     ENDIF
     

     !! Print status
     IF (print_level .EQ. 2 .AND. accepted .GT. 0) THEN
        WRITE(print_unit, *) "  istep, nfev, njev, naev, accepted", istep, nfev, njev, naev, accepted
        !! times 2 to be compatible with potfit definition of Cost
        !!WRITE(print_unit, *) "  Cost, lam, delta", C, lam, delta
        WRITE(print_unit, *) "  Cost, lam, delta", 2.0*C, lam, delta
        WRITE(print_unit, *) "  av, cos alpha", av, cos_alpha
        FLUSH(print_unit)
     ELSEIF(print_level .EQ. 3) THEN
        WRITE(print_unit, *) "  istep, nfev, njev, naev, accepted", istep, nfev, njev, naev, accepted
        !! times 2 to be compatible with potfit definition of Cost
        !!WRITE(print_unit, *) "  Cost, lam, delta", C, lam, delta
        WRITE(print_unit, *) "  Cost, lam, delta", 2.0*C, lam, delta
        WRITE(print_unit, *) "  av, cos alpha", av, cos_alpha
        FLUSH(print_unit)
     ENDIF
     IF (print_level .EQ. 4 .AND. accepted .GT. 0) THEN
        WRITE(print_unit, *) "  istep, nfev, njev, naev, accepted", istep, nfev, njev, naev, accepted
        !! times 2 to be compatible with potfit definition of Cost
        !!WRITE(print_unit, *) "  Cost, lam, delta", C, lam, delta
        WRITE(print_unit, *) "  Cost, lam, delta", 2.0*C, lam, delta
        WRITE(print_unit, *) "  av, cos alpha", av, cos_alpha
        WRITE(print_unit, *) "  x = ", x
        WRITE(print_unit, *) "  v = ", v
        WRITE(print_unit, *) "  a = ", a
        FLUSH(print_unit)
     ELSEIF (print_level .EQ. 5) THEN
        WRITE(print_unit, *) "  istep, nfev, njev, naev, accepted", istep, nfev, njev, naev, accepted
        !! times 2 to be compatible with potfit definition of Cost
        !!WRITE(print_unit, *) "  Cost, lam, delta", C, lam, delta
        WRITE(print_unit, *) "  Cost, lam, delta", 2.0*C, lam, delta
        WRITE(print_unit, *) "  av, cos alpha", av, cos_alpha
        WRITE(print_unit, *) "  x = ", x
        WRITE(print_unit, *) "  v = ", v
        WRITE(print_unit, *) "  a = ", a
        FLUSH(print_unit)
     ENDIF

     ! If converged -- return
     IF(converged .NE. 0) THEN
        exit main
     ENDIF

     IF (accepted .GE. 0) jac_uptodate = .FALSE. !jacobian is now out of date
     
  END DO main
  ! end main loop
  
  ! If not converged
  IF(converged .EQ. 0) converged = -1
  niters = istep

  ! Return best fit found
  ! If the method converged, but final x is different from x_best -- what to do?
  x = x_best
  fvec = fvec_best

  IF(print_level .GE. 1) THEN
     WRITE(print_unit,*) "Optimization finished"
     WRITE(print_unit,*) "Results:"
     WRITE(print_unit,*) "  Converged:    ", converged_info(converged), converged
     !! times 2 to be compatible with potfit definition of Cost
     !!WRITE(print_unit,*) "  Final Cost: ", 0.5d+0*DOT_PRODUCT(fvec,fvec)
     !!WRITE(print_unit,*) "  Cost/DOF: ", 0.5d+0*DOT_PRODUCT(fvec,fvec)/(m-n)
     WRITE(print_unit,*) "  Final Cost: ", DOT_PRODUCT(fvec,fvec)
     WRITE(print_unit,*) "  Cost/DOF: ", DOT_PRODUCT(fvec,fvec)/(m-n)
     WRITE(print_unit,*) "  niters:     ", istep
     WRITE(print_unit,*) "  nfev:       ", nfev
     WRITE(print_unit,*) "  njev:       ", njev
     WRITE(print_unit,*) "  naev:       ", naev
     FLUSH(print_unit)
  ENDIF

END SUBROUTINE geodesiclm
     

     



