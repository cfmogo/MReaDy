module chem_parameters
integer,parameter :: natoms_H=50
integer,parameter :: natoms_O=0
integer,parameter :: nmolec_HH=100
integer,parameter :: nmolec_OOtrip=100
integer,parameter :: nmolec_OOsing=0
integer,parameter :: nmolec_HO=0
integer,parameter :: nmolec_HO2=0
integer,parameter :: nmolec_H2O=0
integer,parameter :: nmolec_H2Ot=0
integer,parameter :: nmolec_O3=0
integer,parameter :: nmolec_H3=0
integer,parameter :: nmolec_H4=0
integer,parameter :: nmolec_H2O2s=0
integer,parameter :: nmolec_H2O2t=0
integer,parameter :: nmolec_O4=0
integer,parameter :: nmolec_H3O=0
integer,parameter :: nmolec_HO3=0
end module chem_parameters

module phys_parameters 
! Units nm K   
double precision, parameter :: temp=2000d0     ! Kelvin
integer         , parameter :: temp_distribution=1 ! 0 - uniform or 1 - Max. Bolt. dist.
double precision,parameter :: boxl=20d-9   ! Box size , normally d-9 (nano meters)
double precision,parameter :: out_bound=0.4d-9  !for molecule formation normally d-9 (nano meters)
double precision,parameter :: in_bound=0.3d-9   ! out dist. can not be smaller then in dist.
double precision,parameter :: rl_cut_off=1.5d-9   !for limit of long range distance interactions
end module phys_parameters

module sim_variables
double precision :: dt=0.05d-15        !Length of time-step, usually to 0.05d-15 (fs))
integer :: t_tot=100000    !Number of production steps 
integer ::tprint=10000         !steps after which to print geometry
integer :: tbackup=100                  !steps after which to print backupfiles
integer :: iseed(12)= (/5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5/)  !Seed for the random generator
integer :: bk_new=0               !Restart from intermediat point(0-no; 1,2,...n - yes)
integer :: which_step=0                !From which step?
double precision :: new_dt=0.05d-15        !Length of time-step, usually to 0.05d-15 (fs))
integer :: new_t_tot=100   !New Number of production steps 
integer :: new_tprint=1
integer :: new_tbackup=1
end module sim_variables


