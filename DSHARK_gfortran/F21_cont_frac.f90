!> Computes the hypergeometric function 2F1 from a continued fraction formula including nfrac terms
!! \param a argument of 2F1
!! \param b argument of 2F1
!! \param c argument of 2F1
!! \param z argument of 2F1
!! \param nfrac required number of terms to be included in continued fraction formula for an accurate evaluation of 2F1
function F21_cont_frac(a,b,c,z,nfrac)
  implicit none
  real(kind=16) :: a
  real(kind=16) :: b
  real(kind=16) :: c
  complex(kind=16) :: z
  integer :: nfrac,j, n
  complex :: F21_cont_frac
  complex(kind=16) ::  sol
  
  sol=1.0_16+((a+1.0_16*nfrac)/(1.0_16+1.0_16*nfrac))*((b+1.0_16*nfrac)/(c+1.0_16*nfrac))*z

  do j=nfrac,2,-1

     sol=1.0_16+( (a+1.0_16*j-1.0_16)/(1.0_16*j) *(b+1.0_16*j-1.0_16)/(c+1.0_16*j-1.0_16))*z -&
          & (a+1.0_16*j)/(1.0_16+1.0_16*j)*(b+1.0_16*j)/(c+1.0_16*j)*z/sol

  enddo

  sol=1.0_16+a*b/c*z/(1.0_16- 0.5_16*(a+1.0_16)*(b+1.0_16)/(c+1.0_16)*z/sol)

  F21_cont_frac=sol

end function F21_cont_frac
