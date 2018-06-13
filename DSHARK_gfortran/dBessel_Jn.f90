!> Computes the first derivative of the Bessel function of the first kind using the relation dJ_n/dz = 1/2*(J_(n-1) - J_(n+1))
!! \param n index of the Bessel function
!! \param z argument of the Bessel function 
function dBessel_Jn(n,z)
  implicit none
  integer :: n
  real :: z, dBessel_Jn

  dBessel_Jn=0.5*(BesJn(n-1,z) - BesJn(n+1,z))

end function dBessel_Jn
