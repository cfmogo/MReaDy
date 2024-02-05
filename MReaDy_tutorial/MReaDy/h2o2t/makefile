
#Comment number 1. Construction of the makefile

 PROGRAMS= suph2o2t_f31.o superficieh2o_f31.o superficieho2_f31.o superficieh2otrip_f31.o param4ct_f31.o
	
#Comment number 3. Variable FC will be the compiler to use
 FC=gfortran

#Variable FCFLAGS will be the opctions to pass to the compiler
 FCFLAGS= 

#Variable FLFLAGS will be the opctions to pass to the loader
 FLFLAGS= 

#Comment number 4. Create executables
  sh2o2t_f31:$(PROGRAMS)
	$(FC) -o sh2o2t_f31 $(PROGRAMS) $(FLFLAGS) 1>A1 2>A2;\
    
#	sh2o2t_f31;\

 param4ct_f31.mod :param4ct_f31.f90
	if [ -f param4ct_f31.mod ]; then rm param4ct_f31.mod; fi
	$(FC) -c $(FCFLAGS) param4ct_f31.f90

 suph2o2t_f31.o: suph2o2t_f31.f90 param4ct_f31.mod
	$(FC) -c $(FCFLAGS) suph2o2t_f31.f90


 superficieh2o_f31.o: superficieh2o_f31.f
	$(FC) -c $(FCFLAGS) superficieh2o_f31.f


 superficieho2_f31.o: superficieho2_f31.f
	$(FC) -c $(FCFLAGS) superficieho2_f31.f

 superficieho2trip_f31.o: superficieho2trip_f31.f
	$(FC) -c $(FCFLAGS) superficieho2trip_f31.f










