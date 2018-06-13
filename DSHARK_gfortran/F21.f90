!> Computes the hypergeometric function 2F1(1,b;c;z) which is needed for evaluating the modified plasma dispersion function given by Eq.(17) in Mace & Hellberg 1995
!! \param b argument of 2F1
!! \param c argument of 2F1
!! \param z argument of 2F1
!! \param nfrac required number of terms to be included in continued fraction formula for an accurate evaluation of 2F1
!! \param F21_val estimated value of 2F1
subroutine F21(b,c,z,nfrac,F21_val)
  use param_mod
  implicit none
  real(kind=16) :: b
  real(kind=16) :: c
  complex(kind=16) :: z
  integer :: nfrac
  complex :: F21_val
  complex :: z1, z2
  real :: mag1, mag2
  real :: arg1, arg2
  complex :: F21_cont_frac

  external F21_cont_frac

  !'Pure' hypergeometric function for given arguments 1,b,c,z will not converge. Thus, have to use transformation formula given, e.g., in Abramowitz & Stegun 1964

  z1=z/(z-1.0)
  z2=1.0-1.0/z

  mag1=sqrt(real(z1)**2+aimag(z1)**2)
  mag2=sqrt(real(z2)**2+aimag(z2)**2)

  if(mag1.lt.mag2) then

     if(mag1.ge.1.0) then
        write(*,*) 'F21 does not converge'
        stop
     endif

     F21_val=F21_cont_frac(1.0_16,c-b,c,z/(z-1.0_16),nfrac)/(1.0_16-z)

  else

     arg1=abs(atan2(aimag(1.0-z),real(1.0-z)))
     arg2=abs(atan2(aimag(z),real(z)))
     
     if((arg1.lt.pi).and.(arg2.lt.pi)) then

        F21_val=((c-1)/(c-b-1))*F21_cont_frac(1.0_16,2.0_16-c,2.0_16-c+b,1.0_16-1.0_16/z,nfrac)/z +&
             &   (gamma(c)/(gamma(b)))*gamma(b-c+1) * z**(1.0_16-c) *(1.0_16-z)**(c-b-1.0_16)

     else
        write(*,*) 'F21 does not converge'
        stop
     endif


  endif


end subroutine F21
