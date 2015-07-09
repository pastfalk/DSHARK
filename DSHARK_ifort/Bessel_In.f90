!> Computes the modified Bessel function of the first kind using relation I_n (z) = sum of z^j / j! * J_(n+j) (z) over j
!! \param n index of the modified Bessel function
!! \param z argument of the modified Bessel function
function Bessel_In(n,z)
  
  real :: Bessel_In, old, error_old, error, z
  integer :: n, j, fac

  Bessel_In=0.0
  error_old=1.0

  !determine modified Bessel function using relation I_n (z) = sum of z^j / j! * J_(n+j) (z) over j
  !restrict sum to 33 iterations otherwise j! becomes to large
  
  do j=0,33     

     !determine j!
  
     if(j.eq.0) then
        fac=1
     else
        fac=fac*j
     endif

     old=Bessel_In

     !compute j-th term of the sum and add it up
     Bessel_In = Bessel_In + z**j/fac * Bessel_Jn(abs(n)+j,z)
     

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

end function Bessel_In
