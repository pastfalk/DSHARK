!> Provides functions which are integrated in subroutine integrator
!! \param z variable which is integrated over in subroutine integrator
!! \param n parameter which is determined in subroutine disp_det()
!! \param h1 parameter which is determined in subroutine disp_det()
!! \param h2 parameter which is determined in subroutine disp_det()
!! \param h3 parameter which is determined in subroutine disp_det()
!! \param case parameter which is determined in subroutine disp_det() and which decides over functional form of integrand
function intgrnd(z,n,h1,h2,h3,case)
  implicit none
  real :: z
  integer :: case, n
  complex :: h1
  real :: h2
  real :: h3
  complex :: intgrnd
  real :: dBessel_Jn
  complex :: Zk_func

  external dBessel_Jn
  external Zk_func
  
  !choose one of three different integrands
 
  if ( case.eq. 1) then
     !integrand appearing in epsilon_xx, epsilon_xz and epsilon_zz
     intgrnd = BesJn(n,h2*sqrt(z-1.0))*BesJn(n,h2*sqrt(z-1.0))/(z**h3) *Zk_func(h1/sqrt(z),h3-1.0)

  else if (case.eq.2) then
     !integrand appearing in epsilon_yy
     intgrnd = (z-1.0)*dBessel_Jn(n,h2*sqrt(z-1.0))*dBessel_Jn(n,h2*sqrt(z-1.0))/(z**h3) *Zk_func(h1/sqrt(z),h3-1.0)

  else if (case.eq.3) then
     !integrand appearing in epsilon_xy and epsilon_yz
     intgrnd = sqrt(z-1.0)*BesJn(n, h2*sqrt(z-1.0))*dBessel_Jn(n, h2*sqrt(z-1.0))/(z**h3) *Zk_func(h1/sqrt(z),h3-1.0)


  else if(case.eq.4) then
     !integrand for the Z function
     intgrnd=exp(-z*z) *(1.0/(z-h1) - 1.0/(z+h1))

  end if

end function intgrnd
