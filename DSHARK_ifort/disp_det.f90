!> Computes the determinant of the dispersion tensor D
!! \param omega considered frequency
!! \param k considered wavenumber
function disp_det(omega,k)
  use param_mod
  implicit none

  complex :: epsilon_xx, epsilon_yy, epsilon_zz
  complex :: epsilon_xy, epsilon_xz, epsilon_yz
  complex :: del_xx, del_yy, del_zz
  complex :: del_xy, del_xz, del_yz

  logical :: esc(6)

  complex :: intxxa, intyya, intxya
  complex :: intxxb, intyyb, intxyb

  real :: lambda
  complex :: zeta1, zeta2

  complex :: omega, disp_det
  real :: k

  complex :: h1a, h1b
  real :: h2
  integer :: h3

  integer :: m, n

  complex :: Z_func, dZ_func
  real :: exp_Bessel_In, exp_dBessel_In
  real :: expBes, expdBes

  external exp_Bessel_In
  external exp_dBessel_In

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

        esc=.true.

        !bi-kappa scenario

        epsilon_xx=epsilon_xx+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2
        epsilon_yy=epsilon_yy+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2
        epsilon_zz=epsilon_zz+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2 *(k*sin(theta))**2 / (k*cos(theta))**2
        epsilon_zz=epsilon_zz+ 2.0*omega**2 *(2.0*kappa(m)-1.0)/(2.0*kappa(m)-3.0) *dens(m)**2 *q(m)**2 /&
             & beta_para(m) / (k *cos(theta))**2
        epsilon_xz=epsilon_xz- (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2 *(k*sin(theta))/(k*cos(theta))

        n=0

        do while(esc(1).or.esc(2).or.esc(3).or.esc(4).or.esc(5).or.esc(6))


           if(n.eq.0) then

              h1a=sqrt((2.0*kappa(m)+2.0)/(2.0*kappa(m)-3.0)) * omega/&
                   & (sqrt(beta_para(m))*k*cos(theta)*sqrt(mu(m))) * sqrt(dens(m))
              h2=k*sin(theta)/sqrt(mu(m)*dens(m))/q(m) *sqrt((2.0*kappa(m)-3.0)*beta_perp(m)/2.0)
              h3=kappa(m)+2

              if(esc(1).or.esc(3).or.esc(5)) then
                 call integrator(n,h1a,h2,h3,1,intxxa)
              endif

              if(esc(2)) then
                 call integrator(n,h1a,h2,h3,2,intyya)
              endif

              if(esc(4).or.esc(6)) then
                 call integrator(n,h1a,h2,h3,3,intxya)
              endif

              epsilon_yy=epsilon_yy+ 2.0*sqrt(2.0) *sqrt(mu(m))*dens(m)**(3.0/2.0) * q(m)**2  *(kappa(m)-1.0/2.0)/&
                   & sqrt(2.0*kappa(m)-3.0) * (kappa(m)+1)**(3.0/2.0) / (sqrt(beta_para(m)) *k*cos(theta)) *&
                   & beta_ratio(m) *omega * intyya 

              epsilon_zz=epsilon_zz+ 4.0*sqrt(2.0)/sqrt( mu(m))*dens(m)**(5.0/2.0) * q(m)**2 * (kappa(m)-1.0/2.0)/& 
                   & (2.0*kappa(m)-3.0)**(3.0/2.0) *(kappa(m)+1)**(3.0/2.0) /beta_perp(m) /&
                   & sqrt(beta_para(m)) / (k*cos(theta))**3 *beta_ratio(m) * omega**3 * intxxa 

              epsilon_yz=epsilon_yz-4.0*i *dens(m)**2 *  q(m)**2 * (kappa(m)-1.0/2.0)/(2.0*kappa(m)-3.0) *&
                   & (kappa(m)+1)**(3.0/2.0) /sqrt(beta_perp(m)) *1.0/((k*cos(theta))**2 *&
                   & sqrt(beta_para(m)))* beta_ratio(m) *omega*omega *intxya 

           else

              h1a=sqrt((2.0*kappa(m)+2.0)/(2.0*kappa(m)-3.0)) * (omega-n*q(m)*mu(m))/&
                   & (sqrt(beta_para(m))*k*cos(theta)*sqrt(mu(m))) * sqrt(dens(m))
              h1b=sqrt((2.0*kappa(m)+2.0)/(2.0*kappa(m)-3.0)) * (omega+n*q(m)*mu(m))/&
                   & (sqrt(beta_para(m))*k*cos(theta)*sqrt(mu(m))) * sqrt(dens(m))
              h2=k*sin(theta)/sqrt(mu(m)*dens(m))/q(m) *sqrt((2.0*kappa(m)-3.0)*beta_perp(m)/2.0)
              h3=kappa(m)+2

              if(esc(1).or.esc(3).or.esc(5)) then
                 call integrator(n,h1a,h2,h3,1,intxxa)
                 call integrator(n,h1b,h2,h3,1,intxxb)
              endif

              if(esc(2)) then
                 call integrator(n,h1a,h2,h3,2,intyya)
                 call integrator(n,h1b,h2,h3,2,intyyb)
              endif

              if(esc(4).or.esc(6)) then
                 call integrator(n,h1a,h2,h3,3,intxya)
                 call integrator(n,h1b,h2,h3,3,intxyb)
              endif

              if(esc(1)) then

                 del_xx=4.0*sqrt(2.0)*mu(m)**(3.0/2.0) *dens(m)**(5.0/2.0) *q(m)**4 * (kappa(m)-1.0/2.0)/&
                      & (2.0*kappa(m)-3.0)**(3.0/2.0) * (kappa(m)+1)**(3.0/2.0) /(beta_perp(m)*&
                      & (k*sin(theta))**2) /(sqrt(beta_para(m))*k*cos(theta))*&
                      & n**2*(beta_ratio(m) * omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m)) *intxxa +&
                      
                      & 4.0*sqrt(2.0)*mu(m)**(3.0/2.0) *dens(m)**(5.0/2.0) *q(m)**4 * (kappa(m)-1.0/2.0)/&
                      & (2.0*kappa(m)-3.0)**(3.0/2.0) * (kappa(m)+1)**(3.0/2.0) /(beta_perp(m)*&
                      & (k*sin(theta))**2) /(sqrt(beta_para(m))*k*cos(theta))*&
                      & n**2*(beta_ratio(m) * omega + n*(beta_ratio(m)-1.0) * mu(m)*q(m)) *intxxb

                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_xx)/ real(epsilon_xx)).gt.eps_error).or. &
                      & (abs(aimag(del_xx)/aimag(epsilon_xx)).gt.eps_error)))) then

                    epsilon_xx=epsilon_xx+ del_xx

                 else
                    esc(1)=.false.
                 endif

              endif

              if(esc(2)) then

                 del_yy = 2.0*sqrt(2.0) *sqrt(mu(m))*dens(m)**(3.0/2.0) * q(m)**2  *(kappa(m)-1.0/2.0)/&
                      & sqrt(2.0*kappa(m)-3.0) * (kappa(m)+1)**(3.0/2.0) / (sqrt(beta_para(m)) *k*cos(theta)) *&
                      & (beta_ratio(m) *omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m)) * intyya +&
                      
                      & 2.0*sqrt(2.0) *sqrt(mu(m))*dens(m)**(3.0/2.0) * q(m)**2  *(kappa(m)-1.0/2.0)/&
                      & sqrt(2.0*kappa(m)-3.0) * (kappa(m)+1)**(3.0/2.0) / (sqrt(beta_para(m)) *k*cos(theta)) *&
                      & (beta_ratio(m) *omega + n*(beta_ratio(m)-1.0) * mu(m)*q(m)) * intyyb

                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_yy)/ real(epsilon_yy)).gt.eps_error).or. &
                      & (abs(aimag(del_yy)/aimag(epsilon_yy)).gt.eps_error)))) then

                    epsilon_yy=epsilon_yy+ del_yy

                 else
                    esc(2)=.false.
                 endif
              endif


              if(esc(3)) then

                 del_zz = 4.0*sqrt(2.0)/sqrt( mu(m))*dens(m)**(5.0/2.0) * q(m)**2 * (kappa(m)-1.0/2.0)/& 
                      & (2.0*kappa(m)-3.0)**(3.0/2.0) *(kappa(m)+1)**(3.0/2.0) /beta_perp(m) /&
                      & sqrt(beta_para(m)) / (k*cos(theta))**3 * (beta_ratio(m) * omega - n*(beta_ratio(m) & 
                      & -1.0) * mu(m)*q(m))*(omega-n*mu(m)*q(m))**2 * intxxa +&
                      
                      & 4.0*sqrt(2.0)/sqrt( mu(m))*dens(m)**(5.0/2.0) * q(m)**2 * (kappa(m)-1.0/2.0)/& 
                      & (2.0*kappa(m)-3.0)**(3.0/2.0) *(kappa(m)+1)**(3.0/2.0) /beta_perp(m) /&
                      & sqrt(beta_para(m)) / (k*cos(theta))**3 * (beta_ratio(m) * omega + n*(beta_ratio(m) & 
                      & -1.0) * mu(m)*q(m))*(omega+n*mu(m)*q(m))**2 * intxxb

                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_zz)/ real(epsilon_zz)).gt.eps_error).or. &
                      & (abs(aimag(del_zz)/aimag(epsilon_zz)).gt.eps_error)))) then

                    epsilon_zz=epsilon_zz+ del_zz

                 else
                    esc(3)=.false.
                 endif
              endif


              if(esc(4)) then

                 del_xy = 4.0*i  * mu(m) *dens(m)**2 * q(m)**3  * 1.0/sqrt(beta_perp(m)) / &
                      & sqrt(beta_para(m))/(k*cos(theta))/(k*sin(theta))*(kappa(m)-1.0/2.0)/(2.0*kappa(m)-3.0)*&
                      & (kappa(m)+1)**(3.0/2.0) *n* (beta_ratio(m) *&
                      & omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m)) * intxya +&
                      
                      & (-4.0)*i  * mu(m) *dens(m)**2 * q(m)**3  * 1.0/sqrt(beta_perp(m)) / &
                      & sqrt(beta_para(m))/(k*cos(theta))/(k*sin(theta))*(kappa(m)-1.0/2.0)/(2.0*kappa(m)-3.0)*&
                      & (kappa(m)+1)**(3.0/2.0) *n* (beta_ratio(m) *&
                      & omega + n*(beta_ratio(m)-1.0) * mu(m)*q(m)) * intxyb


                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_xy)/ real(epsilon_xy)).gt.eps_error).or. &
                      & (abs(aimag(del_xy)/aimag(epsilon_xy)).gt.eps_error)))) then

                    epsilon_xy=epsilon_xy+ del_xy

                 else
                    esc(4)=.false.
                 endif
              endif

              if(esc(5)) then

                 del_xz = 4.0*sqrt(2.0) * sqrt(mu(m)) *dens(m)**(5.0/2.0) * q(m)**3 *  (kappa(m)-1.0/2.0)/&
                      & (2.0*kappa(m)-3.0)**(3.0/2.0)  *&
                      & (kappa(m)+1)**(3.0/2.0) /beta_perp(m) /(k*sin(theta)) /sqrt(beta_para(m))/(k*cos(theta))**2 *&
                      & n * (beta_ratio(m) * omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m))*(omega-n*mu(m)*q(m))*intxxa +&
                      
                      & (-4.0)*sqrt(2.0) * sqrt(mu(m)) *dens(m)**(5.0/2.0) * q(m)**3 *  (kappa(m)-1.0/2.0)/&
                      & (2.0*kappa(m)-3.0)**(3.0/2.0)  *&
                      & (kappa(m)+1)**(3.0/2.0) /beta_perp(m) /(k*sin(theta)) /sqrt(beta_para(m))/(k*cos(theta))**2 *&
                      & n * (beta_ratio(m) * omega + n*(beta_ratio(m)-1.0) * mu(m)*q(m))*(omega+n*mu(m)*q(m))*intxxb         


                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_xz)/ real(epsilon_xz)).gt.eps_error).or. &
                      & (abs(aimag(del_xz)/aimag(epsilon_xz)).gt.eps_error)))) then

                    epsilon_xz=epsilon_xz+ del_xz

                 else
                    esc(5)=.false.
                 endif
              endif

              if(esc(6)) then

                 del_yz = (-4.0)*i *dens(m)**2 *  q(m)**2 * (kappa(m)-1.0/2.0)/(2.0*kappa(m)-3.0) *&
                      & (kappa(m)+1)**(3.0/2.0) /sqrt(beta_perp(m)) *1.0/((k*cos(theta))**2 *&
                      & sqrt(beta_para(m)))* (beta_ratio(m) *&
                      & omega - n*(beta_ratio(m)-1.0) * mu(m)*q(m))*(omega-n*mu(m)*q(m))*intxya +&
                      
                      & (-4.0)*i *dens(m)**2 *  q(m)**2 * (kappa(m)-1.0/2.0)/(2.0*kappa(m)-3.0) *&
                      & (kappa(m)+1)**(3.0/2.0) /sqrt(beta_perp(m)) *1.0/((k*cos(theta))**2 *&
                      & sqrt(beta_para(m)))* (beta_ratio(m) *&
                      & omega + n*(beta_ratio(m)-1.0) * mu(m)*q(m))*(omega+n*mu(m)*q(m))*intxyb

                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_yz)/ real(epsilon_yz)).gt.eps_error).or. &
                      & (abs(aimag(del_yz)/aimag(epsilon_yz)).gt.eps_error)))) then

                    epsilon_yz=epsilon_yz+ del_yz

                 else
                    esc(6)=.false.
                 endif
              endif

           endif

           n=n+1

        enddo

     else

        !bi-Maxwellian scenario

        esc=.true.

        epsilon_xx=epsilon_xx+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2
        epsilon_yy=epsilon_yy+ (beta_ratio(m)-1.0) *mu(m)*dens(m)*q(m)**2

        lambda=(k*sin(theta))**2 * beta_perp(m) / (2.0 * q(m)**2 * mu(m) *dens(m))

        n=0

        do while(esc(1).or.esc(2).or.esc(3).or.esc(4).or.esc(5).or.esc(6))

           expBes=exp_Bessel_In(n,lambda)
           expdBes=exp_dBessel_In(n,lambda)

           if(n.eq.0) then
              zeta1=omega/sqrt(beta_para(m))/(k*cos(theta))/sqrt(mu(m)) *sqrt(dens(m))

              epsilon_yy = epsilon_yy - sqrt(mu(m)) *dens(m)**1.5 * q(m)**2 / sqrt(beta_para(m)) / (k*cos(theta)) *&
                   & (2*lambda*(expdBes-expBes)) *(beta_ratio(m)*omega)*Z_func(zeta1)

              epsilon_yz = epsilon_yz + i/2.0 *dens(m)* q(m) * ( k*sin(theta))/(k*cos(theta)) *&
                   &  (beta_ratio(m)*omega)* (expdBes-expBes) *dZ_func(zeta1)


              epsilon_zz = epsilon_zz - dens(m)**2 * q(m)**2 *omega / beta_perp(m) /&
                   & (k*cos(theta))**2 * omega*beta_ratio(m) * expBes * dZ_func(zeta1) 

           else

              zeta1=(omega-n*mu(m)*q(m))/sqrt(beta_para(m))/(k*cos(theta))/sqrt(mu(m)) *sqrt(dens(m))
              zeta2=(omega+n*mu(m)*q(m))/sqrt(beta_para(m))/(k*cos(theta))/sqrt(mu(m)) *sqrt(dens(m))

              if(esc(1)) then

                 del_xx=sqrt(mu(m))*dens(m)**1.5 *q(m)**2 /sqrt(beta_para(m))/(k*cos(theta))*&
                      & n**2 *expBes / lambda *(beta_ratio(m) * omega - (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta1)+&
                      
                      & sqrt(mu(m))*dens(m)**1.5 *q(m)**2 /sqrt(beta_para(m))/(k*cos(theta))*&
                      & n**2 *expBes / lambda *(beta_ratio(m) * omega + (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta2)


                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_xx)/ real(epsilon_xx)).gt.eps_error).or. &
                      & (abs(aimag(del_xx)/aimag(epsilon_xx)).gt.eps_error)))) then

                    epsilon_xx=epsilon_xx+ del_xx

                 else
                    esc(1)=.false.
                 endif
              endif

              if(esc(2)) then

                 del_yy = sqrt(mu(m)) *dens(m)**1.5 * q(m)**2 / sqrt(beta_para(m)) / (k*cos(theta)) *&
                      & (n**2 * expBes / lambda - 2*lambda*(expdBes-expBes)) *&
                      & (beta_ratio(m)*omega - (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta1)+&
                      
                      sqrt(mu(m)) *dens(m)**1.5 * q(m)**2 / sqrt(beta_para(m)) / (k*cos(theta)) *&
                      & (n**2 * expBes / lambda - 2*lambda*(expdBes-expBes)) *&
                      & (beta_ratio(m)*omega + (beta_ratio(m)-1.0)*n*mu(m)*q(m))*Z_func(zeta2)

                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_yy)/ real(epsilon_yy)).gt.eps_error).or. &
                      & (abs(aimag(del_yy)/aimag(epsilon_yy)).gt.eps_error)))) then

                    epsilon_yy=epsilon_yy+ del_yy

                 else
                    esc(2)=.false.
                 endif
              endif

              if(esc(3)) then

                 del_zz = - dens(m)**2 * q(m)**2 *(omega-n*mu(m)*q(m)) / beta_perp(m) /&
                      & (k*cos(theta))**2 *(omega*beta_ratio(m)-&
                      & (beta_ratio(m)-1.0)*n*mu(m)*q(m)) * expBes * dZ_func(zeta1) -&
                      
                      & dens(m)**2 * q(m)**2 *(omega+n*mu(m)*q(m)) / beta_perp(m) /&
                      & (k*cos(theta))**2 *(omega*beta_ratio(m)+&
                      & (beta_ratio(m)-1.0)*n*mu(m)*q(m)) * expBes * dZ_func(zeta2)

                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_zz)/ real(epsilon_zz)).gt.eps_error).or. &
                      & (abs(aimag(del_zz)/aimag(epsilon_zz)).gt.eps_error)))) then

                    epsilon_zz=epsilon_zz+ del_zz

                 else
                    esc(3)=.false.
                 endif
              endif

              if(esc(4)) then

                 del_xy = i*sqrt(mu(m))*dens(m)**1.5 *q(m)**2 /(k*cos(theta)) / sqrt(beta_para(m)) *n*&
                      & (expdBes-expBes) *(beta_ratio(m) *omega - (beta_ratio(m) -1.0)*n*mu(m)*q(m))* Z_func(zeta1) +&
                      
                      & (-i)*sqrt(mu(m))*dens(m)**1.5 *q(m)**2 /(k*cos(theta)) / sqrt(beta_para(m)) *n*&
                      & (expdBes-expBes) *(beta_ratio(m) *omega + (beta_ratio(m) -1.0)*n*mu(m)*q(m))* Z_func(zeta2)


                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_xy)/ real(epsilon_xy)).gt.eps_error).or. &
                      & (abs(aimag(del_xy)/aimag(epsilon_xy)).gt.eps_error)))) then

                    epsilon_xy=epsilon_xy+ del_xy

                 else
                    esc(4)=.false.
                 endif
              endif

              if(esc(5)) then

                 del_xz = -mu(m)*dens(m)**2 *q(m)**3 / beta_perp(m)/ (k*sin(theta)) / (k*cos(theta))*&
                      & (beta_ratio(m) * omega-n*mu(m)*q(m)*(beta_ratio(m)-1.0))* n*expBes*dZ_func(zeta1)+&
                      
                      & mu(m)*dens(m)**2 *q(m)**3 / beta_perp(m)/ (k*sin(theta)) / (k*cos(theta))*&
                      & (beta_ratio(m) * omega+n*mu(m)*q(m)*(beta_ratio(m)-1.0))* n*expBes*dZ_func(zeta2)           


                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_xz)/ real(epsilon_xz)).gt.eps_error).or. &
                      & (abs(aimag(del_xz)/aimag(epsilon_xz)).gt.eps_error)))) then

                    epsilon_xz=epsilon_xz+ del_xz

                 else
                    esc(5)=.false.
                 endif
              endif

              if(esc(6)) then

                 del_yz = i/2.0 *dens(m)* q(m) * ( k*sin(theta))/(k*cos(theta)) *&
                      &  (beta_ratio(m)*omega-n*mu(m)*q(m)*(beta_ratio(m)-1.0))* (expdBes-expBes) *dZ_func(zeta1)+&
                      
                      & i/2.0 *dens(m)* q(m) * ( k*sin(theta))/(k*cos(theta)) *&
                      &  (beta_ratio(m)*omega+n*mu(m)*q(m)*(beta_ratio(m)-1.0))* (expdBes-expBes) *dZ_func(zeta2)

                 if((n.le.4).or.((n.gt.4).and.((abs(real(del_yz)/ real(epsilon_yz)).gt.eps_error).or. &
                      & (abs(aimag(del_yz)/aimag(epsilon_yz)).gt.eps_error)))) then

                    epsilon_yz=epsilon_yz+ del_yz

                 else
                    esc(6)=.false.
                 endif
              endif

           endif

           n=n+1

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
