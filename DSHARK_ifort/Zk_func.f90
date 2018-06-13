!> Computes the modified plasma dispersion function using the summation formula Eq.(32) in Summers et al. 1994, or - for non-integer kappas - uses Eq.(17) in Mace & Hellberg 1995
!! \param zeta argument of the modified plasma dispersion function
!! \param kap value of kappa parameter
function Zk_func(zeta,kappa_input)
  use param_mod
  implicit none
  complex :: zeta, p, Zk_func
  real  :: fac_ratio
  integer :: r, l
  real :: kappa_input
  integer :: kappa_int
  complex :: F21_val
  integer :: nfrac

  if(ceiling(kappa_input).eq.floor(kappa_input)) then

     kappa_int=nint(kappa_input)

     p=(0.0,0.0)

     !for integer kappa, summation formula of modified dispersion function is given by:
     !Z_kappa (zeta) = - (kappa-1/2)/(2*kappa^3/2) kappa!/(2*kappa)! * sum of (kappa+l)!/l! i^(kappa-l) *( 2/(zeta/sqrt(kappa) + i) over l = 0 to kappa
     !for numerical reasons, the appearing factorials are combined to:
     !kappa!/(2*kappa)! (kappa+l)!/l! = (l+1)*...*(kappa+l) / (kappa+1)*...*(kappa*kappa)

     do l=0,kappa_int

        fac_ratio=1.0

        do r=1,kappa_int
           fac_ratio=fac_ratio*1.0*(l+r)/(1.0*(kappa_int+r))
        enddo

        p=p+fac_ratio*i**(kappa_int-l) *(2.0/(zeta/sqrt(kappa_int*1.0) + i))**(kappa_int-l+1)

     end do

     Zk_func= -(kappa_int-1.0/2.0)/(2.0*kappa_int**(3.0/2.0)) * p

  else

     !for non-integer kappa, the modified dispersion function can be derived using the hypergeometric function 2F1

     if(kappa_input.lt.30.0) then
        nfrac=ceiling(kappa_input-1.0)
     else
        nfrac=ceiling(kappa_input-2.0)
     endif

     call F21(2.0_16*kappa_input+2.0_16,kappa_input+2.0_16,&
          & 0.5_16*(1.0+i*zeta/sqrt(kappa_input)),&
          & nfrac,F21_val)

     Zk_func=i*(kappa_input+0.5)*(kappa_input-0.5)/kappa_input**1.5 /(kappa_input+1.0) *F21_val

  endif


end function Zk_func
