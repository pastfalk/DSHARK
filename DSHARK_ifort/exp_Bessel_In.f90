!> Computes the exponentially scaled modified Bessel function of the first kind using the corresponding series expansion
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function exp_Bessel_In(n,z)
  use param_mod
  implicit none
  real :: expBes, old, error_old, error, z
  real :: exp_Bessel_In
  integer :: n, j
  real :: fac1, fac2
  integer :: m
  real :: gamma_frac
  integer :: zfrac

  external gamma_frac

  if(z.gt.3400.0) then
     write(*,*) 'Argument z in exp_Bessel_In too large! Try quadruple precision...'
     stop
  endif


  expBes=0.0
  error_old=1.0

  do j=0,10000

     old=expBes

     if(abs(n).eq.j) then
        expBes=expBes+2.0**abs(n)

     else if(abs(n).gt.j) then

        fac1=1.0

        do m=j+1,abs(n)
           fac1=fac1*(m*1.0/z)
        enddo
        expBes=expBes+fac1*2.0 **abs(n)

     else if(j.gt.abs(n)) then

        fac1=1.0

        do m=abs(n)+1,j
           fac1=fac1*(z/(1.0*m))
        enddo
        expBes=expBes+fac1*2.0 **abs(n)

     endif

     error=abs(1.0-old/expBes)

     if (error.lt.10.0**(-15)) then
        expBes=1.0/expBes
        exit
     else if(expBes.gt.10.0**300) then
        expBes=0.0
        exit
     endif

  enddo

  do m=1,10000 

     zfrac=2*m+abs(n)


     fac1=1.0

     do j=1,m

        fac1=fac1*(z*exp(-z/zfrac)/(2.0*j))

     enddo

     fac2=1.0

     do j=1,abs(n)

        fac2=fac2*(z*exp(-z/zfrac)/(2*(m+j)))
     enddo

     if(fac1*fac1*fac2.lt.10.0**(-300)) then
        cycle
     else
        old=expBes
        expBes=expBes+fac1*fac1*fac2
     endif

     error=abs(1.0-old/expBes)

     if (error.lt.10.0**(-15)) then
        exit
     endif

     error_old=error

  enddo


  exp_Bessel_In=expBes


end function exp_Bessel_In

