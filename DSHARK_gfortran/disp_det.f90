!> Computes the determinant of the dispersion tensor D
!! \param omega considered frequency
!! \param k considered wavenumber
function disp_det(omega,k)
  use param_mod

  complex :: epsilon_xx, epsilon_yy, epsilon_zz
  complex :: epsilon_xy, epsilon_xz, epsilon_yz

  complex :: intxx, intyy, intzz, intxy, intxz, intyz

  real :: lambda
  complex :: zeta

  complex :: omega, disp_det
  real :: k

  complex :: h1
  real :: h2
  integer :: h3

  integer :: m, n

  complex :: Z_func, dZ_func
  real :: Bessel_In, dBessel_In

  external Bessel_In
  external dBessel_In
  external Z_func
  external dZ_func

  epsilon_xx=(0.0,0.0);
  epsilon_yy=(0.0,0.0);
  epsilon_zz=(0.0,0.0);
  epsilon_xy=(0.0,0.0);
  epsilon_xz=(0.0,0.0);
  epsilon_yz=(0.0,0.0);

  !determine the dielectric tensor components for bi-kappa distributed particles (see Summers et al. (1994))
  !for very high kappa, consider dielectric tensor in the bi-Maxwellian limit (see, e.g., Brambilla (1998))

  !note that epsilon is renormalized: epsilon_ij -> epsilon_ij * (v_A ^2/c^2 * omega^2)  

  epsilon_xx=epsilon_xx+(delta*omega)**2 
  epsilon_yy=epsilon_yy+(delta*omega)**2 
  epsilon_zz=epsilon_zz+(delta*omega)**2 


  do m=1,Nspecies


     if (kappa(m).LT.50) then

        !bi-kappa scenario

        epsilon_xx=epsilon_xx+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2
        epsilon_yy=epsilon_yy+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2
        epsilon_zz=epsilon_zz+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2 *(k*sin(theta))**2 / (k*cos(theta))**2
        epsilon_zz=epsilon_zz+ 2.0*omega**2 *(2.0*kappa(m)-1.0)/(2.0*kappa(m)-3.0) *dens(m)**2 *q(m)**2 /&
             & beta_para(m) / (k *cos(theta))**2
        epsilon_xz=epsilon_xz- (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2 *(k*sin(theta))/(k*cos(theta))


        do n=-NBessel,NBessel

           h1=sqrt((2.0*kappa(m)+2.0)/(2.0*kappa(m)-3.0)) * (omega-n*q(m)*mu(m))/&
           & (sqrt(beta_para(m))*k*cos(theta)*sqrt(mu(m))) * sqrt(dens(m))
           h2=k*sin(theta)/sqrt(mu(m)*dens(m))/q(m) *sqrt((2.0*kappa(m)-3.0)*beta_perp(m)/2.0)
           h3=kappa(m)+2

           call integrator(n,h1,h2,h3,1,intxx)

           epsilon_xx=epsilon_xx+ 4.0*sqrt(2.0)*mu(m)**(3.0/2.0) *dens(m)**(5.0/2.0) *q(m)**4 * (kappa(m)-1.0/2.0)/&
                & (2.0*kappa(m)-3.0)**(3.0/2.0) * (kappa(m)+1)**(3.0/2.0) /(beta_perp(m)*&
                & (k*sin(theta))**2) /(sqrt(beta_para(m))*k*cos(theta))*&
                & n**2*(beta_ratio(m) * omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m)) *intxx

           call integrator(n,h1,h2,h3,2,intyy)

           epsilon_yy=epsilon_yy+ 2.0*sqrt(2.0) *sqrt(mu(m))*dens(m)**(3.0/2.0) * q(m)**2  *(kappa(m)-1.0/2.0)/&
                & sqrt(2.0*kappa(m)-3.0) * (kappa(m)+1)**(3.0/2.0) / (sqrt(beta_para(m)) *k*cos(theta)) *&
                & (beta_ratio(m) *omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m)) * intyy

           intzz=intxx

           epsilon_zz=epsilon_zz+ 4.0*sqrt(2.0)/sqrt( mu(m))*dens(m)**(5.0/2.0) * q(m)**2 * (kappa(m)-1.0/2.0)/& 
                & (2.0*kappa(m)-3.0)**(3.0/2.0) *(kappa(m)+1)**(3.0/2.0) /beta_perp(m) /&
                & sqrt(beta_para(m)) / (k*cos(theta))**3 * (beta_ratio(m) * omega - n*(beta_ratio(m) & 
                & -1.0) * mu(m)*q(m))*(omega-n*mu(m)*q(m))**2 * intzz

           call integrator(n,h1,h2,h3,3,intxy)

           epsilon_xy=epsilon_xy+4.0*i  * mu(m) *dens(m)**2 * q(m)**3  * 1.0/sqrt(beta_perp(m)) / &
                & sqrt(beta_para(m))/(k*cos(theta))/(k*sin(theta))*(kappa(m)-1.0/2.0)/(2.0*kappa(m)-3.0)*&
                & (kappa(m)+1)**(3.0/2.0) *n* (beta_ratio(m) *&
                & omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m)) * intxy

           intxz=intxx

           epsilon_xz=epsilon_xz+4.0*sqrt(2.0) * sqrt(mu(m)) *dens(m)**(5.0/2.0) * q(m)**3 *  (kappa(m)-1.0/2.0)/&
                & (2.0*kappa(m)-3.0)**(3.0/2.0)  *&
                & (kappa(m)+1)**(3.0/2.0) /beta_perp(m) /(k*sin(theta)) /sqrt(beta_para(m))/(k*cos(theta))**2 *&
                & n * (beta_ratio(m) * omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m))*(omega-n*mu(m)*q(m))*intxz

           intyz=intxy

           epsilon_yz=epsilon_yz- 4.0*i *dens(m)**2 *  q(m)**2 * (kappa(m)-1.0/2.0)/(2.0*kappa(m)-3.0) *&
                & (kappa(m)+1)**(3.0/2.0) /sqrt(beta_perp(m)) *1.0/((k*cos(theta))**2 *&
                & sqrt(beta_para(m)))* (beta_ratio(m) *&
                & omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m))*(omega-n*mu(m)*q(m))*intyz

        enddo

     else

        !bi-Maxwellian scenario
        
        epsilon_xx=epsilon_xx+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2
        epsilon_yy=epsilon_yy+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2

        lambda=(k*sin(theta))**2 * beta_perp(m) / (2.0 * q(m)**2 * mu(m) *dens(m))

        do n=-NBessel,NBessel

           zeta=(omega-n*mu(m)*q(m))/sqrt(beta_para(m))/(k*cos(theta))/sqrt(mu(m)) *sqrt(dens(m))


!           write(*,*) m, n, lambda, Bessel_In(n,lambda)


           epsilon_xx=epsilon_xx + sqrt(mu(m))*dens(m)**(3.0/2.0) *q(m)**2 /sqrt(beta_para(m))/(k*cos(theta))*exp(-lambda)*&
                & n**2 *Bessel_In(n,lambda) / lambda *(beta_ratio(m) * omega - (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta)

           epsilon_xy = epsilon_xy + i*sqrt(mu(m))*dens(m)**(3.0/2.0) *q(m)**2 /(k*cos(theta)) / sqrt(beta_para(m)) *n*&
                (dBessel_In(n,lambda)-Bessel_In(n,lambda))*exp(-lambda) *&
                (beta_ratio(m) *omega - (beta_ratio(m) -1.0)*n*mu(m)*q(m))* Z_func(zeta)
           
           epsilon_xz = epsilon_xz - mu(m)*dens(m)**2 *q(m)**3 / beta_perp(m)/ (k*sin(theta)) / (k*cos(theta))*&
                (beta_ratio(m) * omega-n*mu(m)*q(m)*(beta_ratio(m)-1.0))* n*&
                Bessel_In(n,lambda)*exp(-lambda)*dZ_func(zeta)           

           epsilon_yy = epsilon_yy + sqrt(mu(m)) *dens(m)**(3.0/2.0) * q(m)**2 / sqrt(beta_para(m)) / (k*cos(theta)) *&
                & (n**2 * Bessel_In(n,lambda) / lambda - 2.0*lambda*(dBessel_In(n,lambda)-Bessel_In(n,lambda)))*exp(-lambda) *&
                & (beta_ratio(m)*omega - (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta)

           epsilon_yz = epsilon_yz + i/2.0 *dens(m)* q(m) * ( k*sin(theta))/(k*cos(theta)) *&
                &  (beta_ratio(m)*omega-n*mu(m)*q(m)*(beta_ratio(m)-1.0))* (dBessel_In(n,lambda)-Bessel_In(n,lambda)) *&
                exp(-lambda) * dZ_func(zeta)

           epsilon_zz = epsilon_zz - dens(m)**2 * q(m)**2 *(omega-n*mu(m)*q(m)) / beta_perp(m) /&
                & (k*cos(theta))**2 *(omega*beta_ratio(m)-&
                (beta_ratio(m)-1.0)*n*mu(m)*q(m)) * Bessel_In(n,lambda)* exp(-lambda) * dZ_func(zeta)

        enddo

     endif

  enddo

  !the dispersion relation for a collisionless plasma can be generally written as:
  !0=det(kk - 1k^2 + omega^2/c^2 * epsilon)

  !to rewrite this in the form given below, use  k_x = k_para = k*cos(theta), k_y = 0, k_z = k_perp = k*sin(theta) and the symmetry relations of the bi-Maxwellian/bi-kappa dieletric tensor

  disp_det=(epsilon_xx-(k*cos(theta))**2)  *  (epsilon_yy-k**2)  *  (epsilon_zz-(k*sin(theta))**2 ) +&
       & 2.0*epsilon_xy*epsilon_yz*(epsilon_xz+k**2 * sin(theta)*cos(theta) ) -&
       & (epsilon_yy -k**2 )  *  (epsilon_xz + k**2 * sin(theta)*cos(theta) )**2 +&
       & (epsilon_xx-(k*cos(theta))**2 )* epsilon_yz**2 +&
       & (epsilon_zz - (k*sin(theta))**2  )*epsilon_xy**2


end function disp_det
