!> Computes the first derivative of the modified Bessel function of the first kind using the relation dI_n/dz = 1/2*(I_(n-1) + I_(n+1))
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function dBessel_In(n,z)

  implicit none
  integer :: n
  real :: z, dBessel_In, Bessel_In

  external Bessel_In

  dBessel_In=0.5*(Bessel_In(n-1,z) + Bessel_In(n+1,z))

end function dBessel_In
