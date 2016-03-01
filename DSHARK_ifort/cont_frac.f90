subroutine cont_frac( zeta, Z_func )
!  use param_mod
  complex :: zeta
  complex :: Z_func
  real :: a


  a=0.5

  Z_func=zeta/(-zeta**2 + a/(1.0+(a/2.0)/(-zeta**2 +(a/4.0))))    



end subroutine cont_frac
