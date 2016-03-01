!> Computes the plasma dispersion function using J-pole expansion with J=8 (see for example Roennmark 1983)
!! \param zeta complex argument of the plasma dispersion function
function Z_func(zeta)
  use param_mod

  complex :: Z_func, zeta
  complex :: error_func



  if(abs(real(zeta)).lt.26.0) then

     call cerror(i*zeta,error_func)
     Z_func=i*sqrt(pi)*exp(-zeta*zeta)*(1.0+error_func)

  else
     call cont_frac(zeta,Z_func)

  endif



end function Z_func
