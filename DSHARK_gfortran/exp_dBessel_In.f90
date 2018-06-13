!> Computes the first derivative of the modified Bessel function of the first kind using the relation dI_n/dz = 1/2*(I_(n-1) + I_(n+1))
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function exp_dBessel_In(n,z)
  implicit none
  integer :: n
  real :: z
  real :: exp_dBessel_In, exp_Bessel_In

  external exp_Bessel_In

  exp_dBessel_In=0.5*(exp_Bessel_In(n-1,z) + exp_Bessel_In(n+1,z))
  
end function exp_dBessel_In
