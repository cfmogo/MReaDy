
#Comment number 1. Construction of the makefile

PROGRAMS=	mready.o kdtree2_mr.o constants.o parameters.o  variables.o hhs.o \
                generate_HH.o generate_OOtrip.o generate_OOsing.o ohp.o oop.o \
          	potential_total.o read_input_file.o  \
          	verlet_alg.o write_geom.o  ho2.o generate_HO2.o \
		generate_H2O.o generate_H2Ot.o generate_O3.o generate_H3.o h2o.o h2ot.o \
		o3.o h3.o h4.o generate_H4.o generate_H2O2s.o suph2o2_10.o superficieh2o_10.o \
		superficieho2_10.o superficieh2otrip_10.o param4c_10.o generate_O4.o o4.o\
		h3o.o	generate_H3O.o ho3.o generate_HO3.o generate_H.o generate_O.o generate_Ar.o\
		desagrupa.o agrupa.o generate_HO.o o2d.o hht.o suph2o2t_f31.o \
		superficieh2o_f31.o superficieho2_f31.o superficieh2otrip_f31.o\
		param4ct_f31.o generate_H2O2t.o vohq.o o2q.o oh_pH.o switch.o\
		f40_Ar2.o f41_ArO.o f42_ArH.o f43_ArC.o 

#Comment number 3. Variable FC will be the compiler to use
#FC= gfortran  -fopenmp -finit-local-zero  -O0 # -fbounds-check  -finit-local-zero -O0 -g  -Wextra  -pedantic -fimplicit-none  -fbacktrace -Wall -std=f95
#Debug for ifort:
#FC= ifort -mcmodel=medium  -g -traceback -check all -check bounds -ftrapuv -debug all -gen-interfaces -warn interfaces -warn unused -qopenmp
#running ifort
#FC=ifort -mcmodel medium -shared-intel -real_size 64 -qopenmp #-O2 #-traceback -check all  #-parallel -O3# -qopt_report -qopt_report_phase All -c 

FC= gfortran  -g -finit-local-zero  -O0  -fopenmp  #-floop-nest-optimize # -fbounds-check  -finit-local-zero -O0 -g  -Wextra  -pedantic -fimplicit-none  -fbacktrace -Wall -std=f95



#Variable FLFLAGS will be the opctions to pass to the loader  
FLFLAGS= -llapack -lblas


#Comment number 4. Create executables


mready:$(PROGRAMS)
	$(FC) -o mready $(PROGRAMS) $(FLFLAGS) 1>A1 2>A2

constants.o :constants.f90 constants.mod
	$(FC) -c constants.f90

constants.mod :constants.f90 
	$(FC) -c constants.f90

parameters.o :parameters.f90
	$(FC) -c parameters.f90

parameters.mod :parameters.f90
	$(FC) -c parameters.f90

kdtree2_mr.o :kdtree2_mr.f90 kdtree2_module.mod
	$(FC) -c kdtree2_mr.f90


kdtree2_module.mod :kdtree2_mr.f90 
	$(FC) -c kdtree2_mr.f90

variables.o :variables.f90 parameters.o kdtree2_mr.o
	$(FC) -c variables.f90
	
variables.mod :variables.f90
	$(FC) -c variables.f90

f40_Ar2.o :Ar_Pots/f40_Ar2.f90 parameters.o
	$(FC) -c Ar_Pots/f40_Ar2.f90

f41_ArO.o :Ar_Pots/f41_ArO.f90 parameters.o
	$(FC) -c Ar_Pots/f41_ArO.f90

f42_ArH.o :Ar_Pots/f42_ArH.f90 parameters.o
	$(FC) -c Ar_Pots/f42_ArH.f90

f43_ArC.o :Ar_Pots/f43_ArC.f90 parameters.o
	$(FC) -c Ar_Pots/f43_ArC.f90

hhs.o: hhs.f parameters.f90 parameters.o
	$(FC) -c hhs.f
	
hht.o: hht.f90 parameters.f90  parameters.o
	$(FC) -c hht.f90

vohq.o: vohq.f90 parameters.f90 parameters.o
	$(FC) -c vohq.f90

ohp.o: ohp.f parameters.f90 parameters.o
	$(FC) -c ohp.f

oh_pH.o: oh_pH.f90 parameters.o
	$(FC) -c oh_pH.f90

o2d.o: o2d.f parameters.o
	$(FC) -c o2d.f

oop.o: oop.f parameters.o 
	$(FC) -c oop.f

o2q.o: o2q.f90 parameters.o
	$(FC) -c o2q.f90

generate_HH.o: generate_HH.f90 parameters.o constants.o variables.o 
	$(FC) -c generate_HH.f90

generate_HO.o: generate_HO.f90 parameters.o constants.o variables.o 
	$(FC) -c generate_HO.f90

generate_HO2.o: generate_HO2.f90 parameters.o constants.o variables.o 
	$(FC) -c generate_HO2.f90

generate_H2O.o: generate_H2O.f90 parameters.o constants.o variables.o
	$(FC) -c generate_H2O.f90

generate_H2Ot.o: generate_H2Ot.f90 parameters.o constants.o variables.o
	$(FC) -c generate_H2Ot.f90

generate_O3.o: generate_O3.f90 parameters.o constants.o variables.o
	$(FC) -c generate_O3.f90

generate_H3.o: generate_H3.f90 parameters.o constants.o variables.o
	$(FC) -c generate_H3.f90

generate_H4.o: generate_H4.f90 parameters.o constants.o variables.o
	$(FC) -c generate_H4.f90

generate_H2O2s.o: generate_H2O2s.f90 parameters.o constants.o variables.o
	$(FC) -c generate_H2O2s.f90

generate_H2O2t.o: generate_H2O2t.f90 parameters.o constants.o variables.o
	$(FC) -c generate_H2O2t.f90

generate_O4.o: generate_O4.f90 parameters.o constants.o variables.o
	$(FC) -c generate_O4.f90

generate_H3O.o: generate_H3O.f90 parameters.o constants.o variables.o
	$(FC) -c generate_H3O.f90

generate_HO3.o: generate_HO3.f90 parameters.o constants.o variables.o
	$(FC) -c generate_HO3.f90

generate_H.o: generate_H.f90 parameters.o constants.o variables.o
	$(FC) -c generate_H.f90

generate_O.o: generate_O.f90 parameters.o constants.o variables.o
	$(FC) -c generate_O.f90

generate_Ar.o: generate_Ar.f90 parameters.o constants.o variables.o
	$(FC) -c generate_Ar.f90

generate_OOtrip.o: generate_OOtrip.f90 parameters.o constants.o variables.o
	$(FC) -c generate_OOtrip.f90

generate_OOsing.o: generate_OOsing.f90 parameters.o constants.o variables.o
	$(FC) -c generate_OOsing.f90

potential_total.o: potential_total.f90 parameters.o constants.o variables.o kdtree2_mr.o
	$(FC) -c potential_total.f90

mready.o: mready.f90 parameters.o constants.o variables.o
	$(FC) -c mready.f90


ho2.o: ho2.f90 parameters.o constants.o variables.o
	$(FC) -c ho2.f90

h2o.o: h2o.f90 parameters.o constants.o variables.o
	$(FC) -c h2o.f90

h2ot.o: h2ot.f90 parameters.o constants.o variables.o
	$(FC) -c h2ot.f90

o3.o: o3.f90 parameters.o constants.o variables.o
	$(FC) -c o3.f90 

h3.o: h3.f90 parameters.o constants.o variables.o
	$(FC) -c h3.f90 

h4.o: h4.f90 parameters.o constants.o variables.o
	$(FC) -c h4.f90 

o4.o: o4.f90 parameters.o constants.o variables.o
	$(FC) -c o4.f90 

h3o.o: h3o.f90 parameters.o constants.o variables.o
	$(FC) -c h3o.f90

ho3.o: ho3.f90 parameters.o constants.o variables.o
	$(FC) -c ho3.f90

param4c_10.mod: param4c_10.f90
	$(FC) -c param4c_10.f90

param4c_10.o: param4c_10.f90
	$(FC) -c param4c_10.f90

suph2o2_10.o: h2o2/suph2o2_10.f90 param4c_10.o variables.o
	$(FC) -c h2o2/suph2o2_10.f90


superficieh2o_10.o: h2o2/superficieh2o_10.f
	$(FC) -c h2o2/superficieh2o_10.f


superficieho2_10.o: h2o2/superficieho2_10.f
	$(FC) -c h2o2/superficieho2_10.f


superficieh2otrip_10.o: h2o2/superficieh2otrip_10.f
	$(FC) -c h2o2/superficieh2otrip_10.f

read_input_file.o: read_input_file.f90 constants.o variables.o
	$(FC) -c read_input_file.f90

verlet_alg.o: verlet_alg.f90 parameters.o constants.o variables.o
	$(FC) -c verlet_alg.f90

 
write_geom.o: write_geom.f90 parameters.o constants.o variables.o
	$(FC) -c write_geom.f90


#Rot_analysis.o: Rot_analysis.f90 constants.o variables.o
#	$(FC) -c Rot_analysis.f90


desagrupa.o: desagrupa.f90 parameters.o constants.o variables.o
	$(FC) -c desagrupa.f90

agrupa.o:  agrupa.f90 parameters.o constants.o variables.o
	$(FC) -c agrupa.f90

param4ct_f31.mod: h2o2t/param4ct_f31.f90
	if [ -f param4ct_f31.mod ]; then rm param4ct_f31.mod; fi
	$(FC) -c h2o2t/param4ct_f31.f90

suph2o2t_f31.o: h2o2t/suph2o2t_f31.f90 param4ct_f31.mod
	$(FC) -c h2o2t/suph2o2t_f31.f90

superficieh2o_f31.o: h2o2t/superficieh2o_f31.f
	$(FC) -c h2o2t/superficieh2o_f31.f

superficieho2_f31.o: h2o2t/superficieho2_f31.f
	$(FC) -c h2o2t/superficieho2_f31.f

superficieh2otrip_f31.o: h2o2t/superficieh2otrip_f31.f
	$(FC) -c h2o2t/superficieh2otrip_f31.f

switch.o: switch.f90
	$(FC) -c switch.f90


clean: 
	rm -f nohup.out out/*.* fort.* 

kill: 
	rm -f *.o *.mod mready nohup.out out/*.* fort.* A1 A2
