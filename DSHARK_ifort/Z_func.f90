!> Computes the plasma dispersion function using J-pole expansion with J=8 (see for example Roennmark 1983)
!! \param zeta complex argument of the plasma dispersion function
function Z_func(zeta)
  
  complex, dimension (8) :: bj, cj
  complex :: Z_func, zeta

  !plasma dispersion function can be approximated by J-pole expansion: Z(zeta) = sum over b_j/(zeta-c_j) from j = 1 to J
  !use coefficients bj and cj, e.g., given in Xie & Xiao 2014 for J=8

  Z_func=(0.0,0.0)

  bj(1)=(-0.01734012457471826,-0.04630639291680322)
  bj(2)=(-0.7399169923225014,0.8395179978099844)
  bj(3)=(5.840628642184073,0.9536009057643667)
  bj(4)=(-5.583371525286853,-11.20854319126599)
  bj(5)=(-0.01734012457471826,0.04630639291680322)
  bj(6)=(-0.7399169923225014,-0.8395179978099844)
  bj(7)=(5.840628642184073,-0.9536009057643667)
  bj(8)=(-5.583371525286853,11.20854319126599)

  cj(1)=(2.237687789201900,-1.625940856173727)
  cj(2)=(1.465234126106004,-1.789620129162444)
  cj(3)=(0.8392539817232638,-1.891995045765206)
  cj(4)=(0.2739362226285564,-1.941786875844713)
  cj(5)=(-2.237687789201900,-1.625940856173727)
  cj(6)=(-1.465234126106004,-1.789620129162444)
  cj(7)=(-0.8392539817232638,-1.891995045765206)
  cj(8)=(-0.2739362226285564,-1.941786875844713)

  do j=1,8

     if(zeta.eq.cj(j)) then
        write(*,*) 'Error: singularity in Z_func'
        stop
     else
        Z_func=Z_func + bj(j)/(zeta-cj(j))
     endif

  enddo


end function Z_func
