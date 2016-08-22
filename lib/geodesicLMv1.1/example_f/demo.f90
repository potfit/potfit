!==============================================================================
!
! an example to minimize Wood's function using geodesiclm 
! 
! Author: Mingjian Wen (wenxx151@umn.edu), University of Minnesota 
!
!==============================================================================

program demo
implicit none 


! variables
! input
integer n, m
double precision,allocatable:: x(:), fvec(:)
logical analytic_jac, analytic_Avv, center_diff
double precision h1, h2
double precision,allocatable:: dtd(:,:)
INTEGER damp_mode
INTEGER maxiter, maxfev, maxaev, maxjev
double precision maxlam, minlam, artol, Cgoal, gtol, xtol, xrtol, ftol, frtol
INTEGER print_level, print_unit
INTEGER iaccel, ibroyden, ibold, imethod
double precision avmax, initialfactor, factoraccept, factorreject

! output
double precision,allocatable:: fjac(:,:)
INTEGER info
INTEGER niters, nfev, njev, naev
INTEGER converged


!==============================================================================
! set parameters
!==============================================================================
n = 4
m = 6
allocate(x(n), fvec(m), dtd(n,n), fjac(m,n))


! initial guess of parameters
x = (/-3.0,-1.0,-3.0,-1.0/)


analytic_jac=.FALSE.
analytic_Avv=.FALSE.
center_diff=.FALSE.  

h1 = 1.0E-5 
h2 = 0.2
factoraccept = 5
factorreject = 2
maxlam = 1.E7
imethod = 0
initialfactor = 1
damp_mode = 0
maxiter = 500
maxfev = 0
Cgoal = 1E-5
maxjev = 0
iaccel = 0
avmax = 0.75
maxaev = 0
print_level =5
print_unit = 6
ibold = 0
ibroyden = 0
artol = 1.E-5
gtol = 1.5E-8
xtol = 1.E-10
xrtol = 1.5E-8
frtol = 1.5E-8


!==============================================================================
! call geodesiclm and write results
!==============================================================================
call geodesiclm(func, jacobian, Avv, &
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

! show result
write(*,"(a,4f10.5)", advance="yes") "The optimized values are:", x(1), x(2), x(3), x(4)


!==============================================================================
! functions 
!==============================================================================
! release memory
deallocate(x, fvec, dtd, fjac)

contains
! function to compute force 
! Wood's function 
subroutine func(m, n, x, fvec)
integer m, n
double precision x(n), fvec(m)
integer i
  fvec(1) = 10.0*(x(2) - x(1)**2)
  fvec(2) = 1.0 - x(1)
  fvec(3) = sqrt(90.0)*(x(4) - x(3)**2)
  fvec(4) = 1.0 - x(3)
  fvec(5) = sqrt(10.0)*(x(2) + x(4) - 2.0)
  fvec(6) = (x(2) - x(4))/sqrt(10.0)
  
  ! minus reference forces
  do i=1,m
    fvec(i) = fvec(i) - 0.0
  enddo

end subroutine func

! jacobian 
subroutine jacobian(m, n, x, fjac)
integer m, n
double precision x(n), fjac(m,n)
end subroutine jacobian

! Avv
subroutine Avv(m, n, x, v, acc)
integer m, n
double precision x(n), v(n), acc(m)
end subroutine Avv

! callback
subroutine callback(m,n,x,v,a,fvec,fjac,acc,lam,dtd,fvec_new,accepted,info)
integer m, n, accepted, info
double precision x(n), v(n), a(n), fvec(m), fjac(m,n), acc(m), lam, dtd(n,n), fvec_new(m)
end subroutine callback

end program demo
