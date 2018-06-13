!> Reads the input parameters from input.dat
!! \param omega_start first initial frequency guess for the initial wavenumber
!! \param increment frequency increment suggested by the code user to determine the initial guesses for subsequent wavenumbers, before subroutine polyfit() takes over
!! \param kstart start value of the requested wavenumber interval
!! \param kend final value of the requested wavenumber interval
!! \param nk number of steps in the wavenumber interval provided
subroutine read_data(omega_start,increment,kstart,kend,nk)
  use param_mod
  implicit none
  real :: q_in, mu_in, dens_in
  real :: drift_in, beta_para_in, beta_perp_in
  real :: kappa_in
  real :: kstart, kend
  integer :: nk, n
  complex :: omega_start, increment
  real :: omega_r, omega_i, increment_r, increment_i

  !define namelists and read input data using namelists

  namelist /wavenumber/ &
     &  kstart, kend, nk

  namelist /initial_guess/ &
     &  omega_r, omega_i, increment_r, increment_i

  namelist /setup/ &
     &  Nspecies, theta, delta

  namelist /accuracy/ &
     &  rf_error, eps_error


  namelist /species/ &
     & q_in, mu_in, dens_in, drift_in,&
     & beta_para_in, beta_perp_in, kappa_in

  open(unit=17,status='old',file='input.dat')
  read(17,wavenumber)
  read(17,initial_guess)
  read(17,setup)
  read(17,accuracy)

  allocate(mu(Nspecies),drift(Nspecies),q(Nspecies))
  allocate(dens(Nspecies),kappa(Nspecies))
  allocate(beta_para(Nspecies),beta_perp(Nspecies),beta_ratio(Nspecies))

  do n=1,Nspecies
          
     read(17,species)

     q(n)=q_in
     mu(n)=mu_in
     dens(n)=dens_in
     drift(n)=drift_in
     beta_para(n)=beta_para_in
     beta_perp(n)=beta_perp_in
     kappa(n)=kappa_in
     
  enddo

  close(17)

  if(omega_r.eq.0.0) omega_r=10.0**(-12) !algorithm cannot handle exact zero
  if(omega_i.eq.0.0) omega_i=10.0**(-12) !algorithm cannot handle exact zero
  if(theta.eq.0.0) theta=0.0001          !algorithm cannot handle exact zero

  beta_ratio=beta_perp/beta_para

  omega_start=omega_r+i*omega_i
  increment=increment_r+i*increment_i
  theta=theta*2.0*pi/360.0

end subroutine read_data
