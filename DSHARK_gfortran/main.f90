!> Scans through the requested wavenumber interval, computes corresponding frequencies and writes them to output file
program main
  use param_mod
  implicit none
  complex :: omega_start, increment
  real :: start, finish
  integer :: nk, ik
  real :: kstart, kend, dk
  real, allocatable, dimension (:) :: krange
  complex, allocatable, dimension (:) :: solution

  open(unit=7,status='unknown',file='omega.dat')

  call cpu_time(start)

  !read parameters from input.dat
  call read_data(omega_start, increment, kstart, kend, nk)
  
  allocate(krange(nk),solution(nk))

  dk=(kend-kstart)/(1.0*nk)

  do ik=1,nk
     krange(ik)=kstart+(ik-1)*dk
  enddo

  !scan through wavenumber interval
  do ik=1,nk

     write(*,*) ' '
     write(*,'(A7,I2,A10,F12.8)') '-------',ik,'------- k=', krange(ik)

     omega_start=real(omega_start)+i*abs(aimag(omega_start))
     

     !use Muller method to iterate root of dispersion relation
     call muller(omega_start,krange(ik),solution(ik))

     write(*,'(A9,E20.10,A9,E20.10)')  '   omega:', real(solution(ik)), '   gamma:',aimag(solution(ik))


     if ((ik.ge.3).and.(ik.lt.nk))  then

        !if three subsequent solutions omega(k) are found, use quadratic polynomial fit 
        !to guess next starting frequency for Muller iteration
        call polyfit(krange(ik-2:ik+1),solution(ik-2:ik),omega_start)

     else

        !for the first two solution omega(k) guess next starting frequency for Muller iteration
        !by raising the computed omega by an increment which is provided by the user
        omega_start=solution(ik)+increment

     end if

     write(7,'(F12.8,E20.10,E20.10)') krange(ik), real(solution(ik)), aimag(solution(ik))

  enddo

  call cpu_time(finish)
 
  write(*,*) 'Time elapsed: ', finish-start
 
  close(7)

  deallocate(krange,solution)
  deallocate(mu,q,kappa,dens)
  deallocate(beta_para,beta_perp,beta_ratio)


end program main
