objects=param_mod.o main.o disp_det.o intgrnd.o Zk_func.o Z_func.o dZ_func.o  integrator.o intdeiini.o intdei.o dBessel_Jn.o muller.o polyfit.o read_data.o  Bessel_In.o  dBessel_In.o cerror.o cont_frac.o
f90comp = ifort
options=  -r8


dsolve: $(objects)
	$(f90comp) -o dsolve $(options) $(objects)

param_mod.mod: param_mod.o param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

param_mod.o: param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

main.o: param_mod.mod main.f90
	$(f90comp) -c $(options) main.f90

disp_det.o: param_mod.mod disp_det.f90
	$(f90comp) -c $(options) disp_det.f90

intgrnd.o: intgrnd.f90
	$(f90comp) -c $(options)  intgrnd.f90

integrator.o: param_mod.mod integrator.f90
	$(f90comp) -c $(options)  integrator.f90

intdeiini.o: intdeiini.f90
	$(f90comp) -c $(options)  intdeiini.f90

intdei.o: intdei.f90
	$(f90comp) -c $(options)  intdei.f90

dBessel_Jn.o: dBessel_Jn.f90
	$(f90comp) -c $(options)  dBessel_Jn.f90

Zk_func.o: param_mod.mod Zk_func.f90
	$(f90comp) -c $(options)  Zk_func.f90

Z_func.o: Z_func.f90
	$(f90comp) -c $(options)  Z_func.f90

dZ_func.o: dZ_func.f90
	$(f90comp) -c $(options)  dZ_func.f90

muller.o: param_mod.mod muller.f90
	$(f90comp) -c $(options)  muller.f90

polyfit.o: polyfit.f90
	$(f90comp) -c $(options)  polyfit.f90

read_data.o: param_mod.mod read_data.f90
	$(f90comp) -c $(options)  read_data.f90

Bessel_In.o: param_mod.mod  Bessel_In.f90
	$(f90comp) -c $(options)  Bessel_In.f90

dBessel_In.o: dBessel_In.f90
	$(f90comp) -c $(options)  dBessel_In.f90

cerror.o: param_mod.mod  cerror.f90
	$(f90comp) -c $(options)  cerror.f90

cont_frac.o: cont_frac.f90
	$(f90comp) -c $(options)  cont_frac.f90

clean:
	rm dsolve 
	rm param_mod.mod
	rm $(objects)
