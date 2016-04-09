!> Scans through the requested wavenumber interval, computes corresponding frequencies and writes them to output file
program main
  use param_mod
  implicit none

  complex :: omega_start, increment
  real :: start, finish
  integer :: ksteps, n
  real :: kstart, kend, dk
  real, allocatable, dimension (:) :: krange
  complex, allocatable, dimension (:) :: solution

  open(unit=7,status='unknown',file='omega.dat')

  call cpu_time(start)

  !read parameters from input.dat
  call read_data(omega_start, increment, kstart, kend, ksteps)
  
  allocate(krange(ksteps),solution(ksteps))

  dk=(kend-kstart)/(1.0*ksteps)
  do n=1,ksteps
     krange(n)=kstart+(n-1)*dk
  enddo

  !scan through wavenumber interval
  do n=1,ksteps

     write(*,*) ' '
     write(*,*) '-------',n,'-------'

     !use Muller method to iterate root of dispersion relation
     call muller(omega_start,krange(n),solution(n))

     write(*,*) krange(n), '  ', solution(n)

     if ((n .ge. 3).and.(n .lt. ksteps))  then

        !if three subsequent solutions omega(k) are found, use quadratic polynomial fit 
        !to guess next starting frequency for Muller iteration
        call polyfit(krange(n-2:n+1),solution(n-2:n),omega_start)

     else

        !for the first two solution omega(k) guess next starting frequency for Muller iteration
        !by raising the computed omega by an increment which is provided by the user
        omega_start=solution(n)+increment

     end if

     write(*,*) 'start: ',  omega_start

     write(7,103) krange(n), real(solution(n)), aimag(solution(n))
103  format(1x, 3(1x, 1e14.7))

  enddo

  call cpu_time(finish)
 
  write(*,*) 'Time elapsed: ', finish-start
 
  close(7)

  deallocate(krange,solution)
  deallocate(mu,q,kappa)
  deallocate(beta_para,beta_perp,beta_ratio)


end program main
