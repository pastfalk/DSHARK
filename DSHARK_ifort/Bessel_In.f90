!> Computes the modified Bessel function of the first kind
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function Bessel_In(n,z)

  real :: Bessel_In, old, error_old, error, z
  integer :: n, j
  real :: fac1, fac2
  integer :: m

  Bessel_In=0.0
  error_old=1.0

  if((z.lt.666.0)) then

     !determines modified Bessel function using the corresponding series expansion

     do m=0,1000  

        fac1=1.0
        
        do j=1,m
        
           fac1=fac1*(z/(2*j))

        enddo

        fac2=1.0

        do j=1,abs(n)

           fac2=fac2*(z/(2*(m+j)))

        enddo

        old=Bessel_In

        Bessel_In=Bessel_In+fac1*fac1*fac2


        error=abs(1.0-old/Bessel_In)

        if (error.lt.10.0**(-12)) then
           exit
        endif

        if ((m.ge.2).and.(abs(Bessel_In).le.10.0**(-50))) then
           Bessel_In=0.0
           exit
        endif

        error_old=error

     enddo


  else

     write(*,*) 'Warning: Bessel_In is getting really big now! You may change to kind=16 (quadruple precision).'

  endif


end function Bessel_In

