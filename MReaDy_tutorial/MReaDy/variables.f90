module calc_parameters
use chem_parameters
integer, parameter :: npart=natoms_H+natoms_O+ &
           2*(nmolec_HH+nmolec_OOtrip+nmolec_OOsing+nmolec_HO)+ & 
     3*(nmolec_HO2+nmolec_H2O+nmolec_H2Ot+nmolec_O3+nmolec_H3)+ &
    4*(+nmolec_H4+nmolec_H2O2s+nmolec_H2O2t+nmolec_O4+nmolec_H3O+nmolec_HO3)


integer, parameter :: lgdiat=npart*(npart-1)/2
end module calc_parameters

module bk_variables
integer :: old_npart!,old_nunit_max
double precision ::old_temperature,old_boxl,old_out_bound, old_in_bound,old_rl_cut_off,old_dt


end module bk_variables


module variables
use calc_parameters
use kdtree2_module

implicit none

double precision    :: potential,massinput,ekin,r1te,r2te,r3te,potential_change
double precision    :: ekinr,ekinc
integer :: iget(12)



integer :: rl_part(0:npart)

double precision, dimension(npart)   :: axold,ayold,azold
double precision, dimension(npart)   :: ax,ay,az
double precision, dimension(npart)   :: vx,vy,vz
double precision, dimension(npart)   :: rx,ry,rz
double precision, dimension(npart)   :: rxold,ryold,rzold
double precision, dimension(npart)   :: dpotx,dpoty,dpotz 
double precision, dimension(npart)   :: dpotxd,dpotyd,dpotzd
double precision, dimension(npart)   :: mass,rmass
character*2     :: cname(0:npart)
integer         , dimension(npart,9) :: group,exgroup,group2




double precision :: potential_ant



integer, dimension(lgdiat,2) :: hh_pair,oo_pair,ho_pair,ArO_pair,Ar2_pair,ArH_pair
!integer, dimension :: hh_pair(lgdiatHH,2),oo_pair(lgdiatOO,2),ho_pair(lgdiatHO,2)
!integer, dimension :: ArO_pair(lgdiatArO,2),Ar2_pair(lgdiatArAr,2),ArH_pair(lgdiatArH,2)


double precision :: potentiald,potentiald_ant

integer :: tempo=-1
integer :: print_flag
double precision    ::total_real_time=0.0d0 ! total_real_time is the total time passed since first bk.

integer(kind=8) :: t_print_back_up,old_step_print
integer(kind=8) :: bk_old
character(len=12) ::    ini_bk_file_name,ex_file_name,geo_file_name,bk_new_string 
character(len=12) ::  bk_old_string  
integer :: mini ! for mininum configuration determination  



double precision, allocatable :: limits(:,:),fullxyz(:,:),pseudofullxyz(:,:)

double precision :: boxlim(1,2)


!for kdtree
integer, parameter :: nres=200
type(kdtree2), pointer :: tree
type(kdtree2_result),dimension(nres) :: results
integer,allocatable :: ref_fullxyz(:)


end module variables





