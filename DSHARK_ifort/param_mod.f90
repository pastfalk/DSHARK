!> Module containing all necessary parameters
module param_mod

implicit none

complex, parameter :: i=(0.0,1.0)
real, parameter :: pi= 4.*atan(1.)
integer :: Nspecies
real :: theta
integer :: NBessel
real, allocatable, dimension (:) :: mu
integer, allocatable, dimension (:) :: q
real, allocatable, dimension (:) :: beta_para
real, allocatable, dimension (:):: beta_perp
real, allocatable, dimension (:) :: beta_ratio
integer, allocatable, dimension (:) :: kappa
real :: delta 
real :: int_error
real :: rf_error
integer :: acc_measure


end module param_mod
