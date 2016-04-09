!> Rootfinding algorithm based on the Muller method presented, e.g., in Gerald & Wheatley (2003)
!! \param omega_start initial frequency guessas starting value for the iteration
!! \param k considered wavenumber
!! \param sol approximated root of the dispersion relation obtained from Muller iteration
subroutine muller(omega_start,k,sol)
  use param_mod
  implicit none

  real :: k
  complex :: omega_start, sol, disp_det
  complex, dimension(1:4) :: fx, omega
  complex :: a, b, c, d1, d2
  integer :: n, j

  external disp_det

  !Muller's method requires three starting points
  !one point is chosen by the user - take the other two points to be left and right of it 
  omega(1)=0.9* omega_start
  omega(2)=omega_start
  omega(3)=1.1*omega_start

  fx(1)=disp_det(omega(1),k)
  fx(2)=disp_det(omega(2),k)
  fx(3)=disp_det(omega(3),k)

  !restrict Muller iteration to 1000 steps

  do n=1,1000

     !determine coefficients from preceding three points

     a=((omega(2)-omega(3))*(fx(1)-fx(3))-(omega(1)-omega(3))*(fx(2)-fx(3)))/&
          & ((omega(1)-omega(3))*(omega(2)-omega(3))*(omega(1)-omega(2)))
     b=((omega(1)-omega(3))**2 * (fx(2)-fx(3)) - (omega(2) - omega(3))**2 *&
          & (fx(1)-fx(3)))/((omega(1)-omega(3))*(omega(2)-omega(3))*(omega(1)-omega(2)))
     c=fx(3)
    
     d1=b+sqrt(b**2 - 4.0*a*c)
     d2=b-sqrt(b**2 - 4.0*a*c)

     !compute new root from coefficients

     if  (abs(d1) .GE. abs(d2)) then
       
        omega(4)=omega(3)-2.0*c/d1

     else

        omega(4)=omega(3)-2.0*c/d2

     endif

    fx(4)=disp_det(omega(4),k)
    
    !measure the accuracy of iterated root and check exit-condition
    !user can choose between two different methods of measuring accuracy

    if(acc_measure==0) then

       !relative difference between two successive roots
       !more reliable than backward error, but more demanding (slower convergence of iteration)
       if ((abs(real(omega(4))/real(omega(3))-1.0) .lt. rf_error) .and. &
            (abs(aimag(omega(4))/aimag(omega(3))-1.0) .lt. rf_error)) exit

    else if(acc_measure==1) then

       !backward error, i.e., how close is function value to zero
       !can be problematic, when function very shallow around root
       if ((abs(real(fx(4))) .lt. rf_error).and.(abs(aimag(fx(4))) .lt. rf_error)) exit

    else

       write(*,*) 'Change accuracy measure to 0 or 1'
       write(*,*) '0 for relative error'
       write(*,*) '1 for backward error'
       stop

    endif

    !stop iteration if last step was ineffective
    if( (abs((real(fx(4))-real(fx(3)))/real(fx(4))) .lt. 10.0**(-12)).and. &
         ( abs((aimag(fx(4))-aimag(fx(3)))/aimag(fx(4))) .lt. 10.0**(-12))) then
       write(*,*) 'Last step in Muller iteration was ineffective'
       exit
    endif

    do j=1,3
       omega(j)=omega(j+1)
       fx(j)=fx(j+1)
    enddo

    if(n.eq.1000) then
       write(*,*) 'Error: Muller method did not converge'
       write(*,*) 'You may try lower accuracy'
       stop
    endif

  enddo

  !solution of root finding procedure
  sol=omega(4)


end subroutine muller
