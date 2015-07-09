!> Computes the modified plasma dispersion function using the summation formula Eq. (32) in Summers et al. 1994
!! \param zeta argument of the modified plasma dispersion function
!! \param kap value of kappa parameter
function Zk_func(zeta,kap)

  use param_mod

  complex :: zeta, p, Zk_func
  real  :: fac_ratio
  integer :: r, j, kap
  
  p=(0.0,0.0)

  !summation formula for modified dispersion function is given by:
  !Z_kappa (zeta) = - (kappa-1/2)/(2*kappa^3/2) kappa!/(2*kappa)! * sum of (kappa+l)!/l! i^(kappa-l) *( 2/(zeta/sqrt(kappa) + i) over l = 0 to kappa

  !due to numerical reasons the appearing factorials are combined to:
  !kappa!/(2*kappa)! (kappa+l)!/l! = (l+1)*...*(kappa+l) / (kappa+1)*...*(kappa*kappa)

  do l=0,kap

     fac_ratio=1.0

     do r=1,kap
        fac_ratio=fac_ratio*1.0*(l+r)/(1.0*(kap+r))
     enddo

     p=p+fac_ratio*i**(kap-l) *(2.0/(zeta/sqrt(kap*1.0) + i))**(kap-l+1)

  end do

  Zk_func= -(kap-1.0/2.0)/(2.0*kap**(3.0/2.0)) * p


end function Zk_func
