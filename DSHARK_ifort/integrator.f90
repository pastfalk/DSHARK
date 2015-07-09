!> Calls numerical integration routines and returns the solution of the performed integration
!! \param n parameter needed for the integrand
!! \param h1 parameter needed for the integrand
!! \param h2 parameter needed for the integrand
!! \param h3 parameter needed for the integrand
!! \param case parameter needed for the integrand
!! \param sol solution of the numerical integration
subroutine integrator(n,h1,h2,h3,case,sol)
  use param_mod
  integer :: case, n
  complex :: h1
  real  :: h2
  integer :: h3

  real :: a
  complex :: sol, intgrnd
  integer :: lenaw
  real :: err, tiny
  real, allocatable :: aw(:)

  external intgrnd

  tiny = 1.0d-307
  lenaw=8000

  allocate(aw(0 : lenaw - 1))

  !lower boundary of the integration
  a=1.0

  !initialize double exponential quadrature method
  call intdeiini(lenaw, tiny, int_error, aw)

  !perform integration using double exponential quadrature
  !upper boundary of integration is infinity
  call intdei(intgrnd,n,h1,h2,h3,case, a, aw, sol, err)

end subroutine integrator
