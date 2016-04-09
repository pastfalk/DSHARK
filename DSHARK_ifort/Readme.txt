DSHARK - A dispersion solver for homogeneous plasmas with anisotropic kappa distributions 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a manual for the Fortran-90 code DSHARK containing a short description of the program, an explanation of the input parameters and some advices for an efficient usage of the code. Note that this code version is adapted to the GNU Fortran Compiler (gfortran).

General remarks
---------------
DSHARK is a dispersion relation solver which can derive the frequencies and growth rates of obliquely propagating waves in homogeneous, bi-kappa distributed plasmas with an arbitrary number of particle species. Kappa distributions are a generalization to the standard Maxwell-Boltzmann distributions, exhibiting long suprathermal tails. Kappa can adopt integer values 1 < kappa <= infinity, where the suprathermal tails are extending for a decreasing kappa parameter. For kappa -> infinity, the limiting case of Maxwell-Boltzmann distribution is recovered. Kappa distributions are omnipresent in space plasmas and can be directly observed in the solar wind.

In DSHARK, the dielectric tensors for both general bi-kappa distributions and bi-Maxwellians are implemented. The expressions for the dielectric tensor components for bi-kappa distributions are based on Summers et al. (1994). The components of the bi-Maxwellian dielectric tensor were derived from Brambilla (1998).

DSHARK is intended to enable a systematic study of parallel and oblique wave propagation in bi-kappa distributed plasmas. You are welcome to use and distribute this code, and to adapt DSHARK to your own purposes. If you publish some work, which is based on results obtained from DSHARK, please cite Astfalk et al. (2015) and make sure, that you mark any changes of the code with respect to the original source code.


The program structure
---------------------
The core of DSHARK is an iterative root finding algorithm enclosed by a loop of the considered wavenumber interval. Before the loop is started, all required parameters are read from an input file 'input.dat' by the routine 'read_parameters()' and all necessary data structures are initialized.

The input provides the iterative root finding algorithm with the initial frequency guesses for the first three wavenumbers of the considered interval. For all subsequent wavenumbers, the initial guesses are determined by the routine 'polyfit()' which uses quadratic polynomials to interpolate previous solutions.

Starting from the supplied initial guess and the given wavenumber, the routine 'muller()' iterates a complex root of the dispersion relation using the Muller method. Thus, for every iteration the determinant of the dispersion tensor has to be evaluated which requires the determination of the dielectric tensor components. This is done by the routine 'disp_det()'.

All necessary integrations are carried out by the routine 'integrator()' which uses a double exponential quadrature method developed by T. Ooura.

The evaluation of the plasma dispersion function and the modified plasma dispersion function is done by the separate functions 'Z_func()' and 'Zk_func'. After the loop successfully cycled through the wavenumber interval, all roots are printed to an output file.

(see Astfalk et al.(2015))


The input parameters
--------------------

&wavenumber

kstart - The lower border of the wavenumber interval the user intends to investigate.

kend   - The upper border of the wavenumber interval the user intends to investigate.

ksteps - Here, the user can choose the required wavenumber resolution. kstep determines at how many points DSHARK will evaluate the dispersion relation within the chosen wavenumber interval.


Note:
All wavenumbers are given in units of the inertial length, d, of the first particle species.



&initial_guess

omega_r     - The initial guess for the real frequency from which the Muller method, implemented in DSHARK, starts to find a root, omega(k), of the dispersion relation at the wavenumber k(1)=kstart.

omega_i     - The initial guess for the growth rate from which the Muller method, implemented in DSHARK, starts to find a root, omega(k), of the dispersion relation at the wavenumber k(1)=kstart.

increment_r - The frequency value by which the previously found root, omega(k), is incremented to provide the starting value for the next Muller iteration at the subsequent wavenumber, k(i+1)=k(i)+dk.

increment_i - The growth rate value by which the previously found root, omega(k), is incremented to provide the starting value for the next Muller iteration at the subsequent wavenumber, k(i+1)=k(i)+dk.


Note:
A proper initial guess is crucial for a successive root finding. If the guess lies too far away from the dispersion branch of interest, you may land on another branch. In general, DSHARK always converges to a certain root. It has to be figured out by the user, whether this is the root he was searching for.

The increments are only necessary for the initial guesses for the root finding at k(2) and k(3). For subsequent wavenumbers a quadratic polynomial approximation determines all following initial guesses. If dk is not too high, this works efficiently.

Both, frequencies and growth rates are always given in units of the gyro frequency of the first particle species.



&setup

Nspecies - The number of particle species the user wants to include.

theta 	 - The propagation angle of the waves, i.e. the angle between the wave vector k and the background magnetic field (which is aligned with the z-axis in the chosen coordinate system).

delta 	 - Ratio of gyro frequency and plasma frequency for the first particle species.


Note:

The parallel and perpendicular wavenumbers are given as k_para=k*cos(theta) and k_perp=k*sin(theta).

Delta gives a measure for the magnetization of the plasma. Low delta corresponds to weak, high delta corresponds to strong magnetization.



&accuracy

acc_measure - Determines the way, the accuracy of the Muller iterated roots is computed. Choose acc_measure=0 for the relative difference between two successive roots, whereas for acc_measure=1, DSHARK computes the backward error, i.e. how close the determinant of the dispersion tensor is to zero, for the current root.  The former is more reliable but also more demanding than the latter.

rf_error    - The 'root finding error' gives the exit-condition for the Muller iteration. It depends on the chosen accuracy measure. For acc_measure=0, an error of 1.0d-2 or 1.0d-3 produces good accuracy, whereas for acc_measure=1 an error of 1.0d-12 generally gives good results. But, of course, the choice depends on the accuracy requested by the user.

int_error   - The 'integration error' gives the exit condition for the numerical integrations. Generally, it is chosen to be 1.0d-15 which gives good accuracy and great performance.

eps_error   - The 'epsilon error' gives the exit condition for the sum over the Bessel index n. Once the relative contribution of the computed dielectric tensor components for a given n gets smaller than the given eps_error, the code exits the loop over the n in disp_det.f90. This ensures that the code 


Note:
If a solution seems fishy, play with these parameters and check whether the solution is numerically converged.

Choose the rf_error to be not too demanding, otherwise DSHARK may run into convergence problems.



&species

q_in         - Charge of the particles in units of the charge of the first particle species.

mu_in 	     - Mass of the particles in units of the mass of the first particle species.

dens_in	     - Density of the particles in units of the density of the first particle species.

beta_para_in - Beta parameter parallel to the background magnetic field.

beta_perp_in - Beta parameter perpendicular to the background magnetic field.

kappa_in     - Kappa parameter of the particles' distribution function.


Note:
If you need more than the two default particle species, just add additional parameter blocks below the two present &species blocks. The choice, which particle species is declared in the first &species block, is of major importance since the normalization of all output data depends on this choice. E.g., if you choose protons to be the first particle species, then all frequencies and growth rates will be given in units of the proton gyrofrequency and the wavenumbers will be in units of the proton inertial length.

The chosen kappa value significantly affects the performance of DSHARK. The larger the kappa parameter, the more demanding is the evaluation of the modified plasma dispersion function and, hence, the slower is the program execution. If kappa exceeds a certain limit, kappa = 50, the code will switch from  bi-kappa distributions to the bi-Maxwellian limit which can be solved much faster. If you are interested in the dispersion properties of bi-kappa distributions with kappa >= 50, you can avoid the switch to the bi-Maxwellian case by manually changing the default limit in disp_det.f90. 



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For more information about DSHARK, see Astfalk et al. (2015). If you have further questions concerning the usage of the code or if you like to discuss some general issues, feel free to write an email to patrick.astfalk@ipp.mpg.de.
