!> Computes the modified Bessel function of the first kind using relation I_n (z) = sum of z^j / j! * J_(n+j) (z) over j
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function Bessel_In(n,z)

  real :: Bessel_In, old, error_old, error, z
  integer :: n, j, fac, fac1, fac2
  integer :: m
  integer :: corr, corr2
  parameter(big=1.0e10)

  Bessel_In=0.0
  error_old=1.0


  if(z.lt.5.0) then

     !determine modified Bessel function using relation I_n (z) = sum of z^j / j! * J_(n+j) (z) over j
     !restrict sum to 33 iterations otherwise j! becomes to large
     !attention: this method gives inaccurate result when z becomes too large

     do j=0,33     

        !determine j!

        if(j.eq.0) then
           fac=1
        else
           fac=fac*j
        endif

        old=Bessel_In

        !compute j-th term of the sum and add it up
        Bessel_In = Bessel_In + z**j/fac * BesJn(abs(n)+j,z)


        !exit the iteration if new term gave relative contribution to the total sum of less than 10^-6
        !or if relative contribution increases after more than 3 iterations due to numerical issues
        error=abs(1.0-old/Bessel_In)

        if (error.lt.10.0**(-6)) then
           exit
        endif

        if ((j.ge.4).and.(error_old.lt.error)) then
           Bessel_In=old
           exit
        endif

        error_old=error

        !stop program execution if iteration did not converge after 33 steps
        if (j.eq.33) then
           write(*,*) 'Error: Bessel_In(n,z) does not converge'
           stop
        endif

     enddo


  else

     !determines modified Bessel function using the corresponding series expansion
     !this method is slower than the one given above, but gives accurate results also for high z

     corr=0
     corr2=0

     fac2=1.0
     fac1=1.0

     do m=1,abs(n) 
        fac2=fac2*m
        if(fac2.gt.big) then
           fac2=fac2/big
           corr2=corr2+1
        endif

     end do

     do j=0,1000  

        if(j.eq.0) then
           fac1=1.0
        else
           fac1=fac1*1.0*j
           fac2=fac2*(j+abs(n))*1.0
        endif

        if(fac1.gt.big) then
           fac1=fac1/big
           corr=corr+1
        endif

        if(fac2.gt.big) then
           fac2=fac2/big
           corr2=corr2+1
        endif

        old=Bessel_In

        Bessel_In = Bessel_In + 1.0/(big**corr)/fac1 * 1.0/(big**corr2)/fac2*&
             & (z/2.0)**(2*j+abs(n))


        error=abs(1.0-old/Bessel_In)

        if (error.lt.10.0**(-12)) then
           exit
        endif

        if ((j.ge.2).and.(abs(Bessel_In).le.10.0**(-50))) then
           Bessel_In=0.0
           exit
        endif


        if ((j.ge.4).and.(error_old.lt.error)) then
           Bessel_In=old
           exit
        endif

        error_old=error

     enddo


  endif


end function Bessel_In
