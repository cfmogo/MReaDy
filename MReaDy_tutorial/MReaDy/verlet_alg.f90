SUBROUTINE verlet_alg
use constants
use variables
use phys_parameters
use sim_variables
use omp_lib      
implicit none

double precision :: sumvsq,sumvx,sumvy,sumvz,dtsq
double precision :: dt2
integer :: i,j,k,t,icount 




icount=0
time:   do  t = 1, t_tot

tempo=t



print_flag=0 ! flag for printing times when there is a change in matrix group
total_real_time=total_real_time+dt



!#####according to tbackup,  write to backup file the positions
ekin=0.0d0
sumvsq = 0.0d0
sumvx  = 0.0d0
sumvy  = 0.0d0
sumvz  = 0.0d0
dtsq=(dt)**2
dt2=2.0d0*dt

 

rxold=rx
ryold=ry
rzold=rz

!$OMP  PARALLEL DO DEFAULT(SHARED) 
do i=1,npart
rx(i)=rx(i) + vx(i)*dt + 0.5d0*dtsq*ax(i)
ry(i)=ry(i) + vy(i)*dt + 0.5d0*dtsq*ay(i)
rz(i)=rz(i) + vz(i)*dt + 0.5d0*dtsq*az(i)
end do
!$OMP END PARALLEL DO 


axold=ax
ayold=ay
azold=az

potential=0.0d0
dpotx=0.0d0
dpoty=0.0d0
dpotz=0.0d0

call potential_total
   
!$OMP  PARALLEL do default(SHARED) reduction(+:ekin) 
do i=1,npart



ax(i) = - dpotx(i) * rmass(i)      
ay(i) = - dpoty(i) * rmass(i)
az(i) = - dpotz(i) * rmass(i)

vx(i)=vx(i) + 0.5d0 * dt * (axold(i) +  ax(i))
vy(i)=vy(i) + 0.5d0 * dt * (ayold(i) +  ay(i))
vz(i)=vz(i) + 0.5d0 * dt * (azold(i) +  az(i))

rx(i) = rx(i) - boxl*floor(rx(i)/boxl)
ry(i) = ry(i) - boxl*floor(ry(i)/boxl)
rz(i) = rz(i) - boxl*floor(rz(i)/boxl)


ekin=ekin+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)

end do
!$OMP END PARALLEL DO


group2=group

call desagrupa !######### where connections between molecules are undone#



call agrupa !########### where connections between molecules are created#


! Agrupa subroutine, after merging to units into one unit, will leave an empty (only zeros) line. 
! This forces to print out two times to time_change.csv below. 
! First when there is a real merging from agrupa script, and then, after a new step the matrix group 
! is updated in the desagrupa script, moving the new zero line to the end of the matrix. 
! The next part of the code,  is to get the agrupa matrix right, avoiding a second print. 

!**************************************************************
!update positions if there is an empty line group(i,8)=0
! and delete line





exgroup=0
j=0  !Number of lines from matrix exgroup               

prepartic0:     do i=1,npart
! Check if there is an empty line in group
if (group(i,2).eq.0) then     


!Run lines until emtpy line 
prepartic1:                     do k=1,j,1
!If there is a reference A BEFORE empty line calling a reference B after empty line
if(exgroup(k,8).gt.i) then
exgroup(k,8)=exgroup(k,8)-1   !update first reference A in exgroup       
end if
end do prepartic1

!Run lines after empty line
prepartic2:                     do k=i+1,npart
!If there is a reference A AFTER empty line calling a reference B after empty line 
if(group(k,8).gt.i) then
group(k,8)=group(k,8)-1         !update first reference A in group
end if

end do prepartic2

else
j=j+1
!if there is no empty line then copy line from group to exgroup 
exgroup(j,:)=group(i,:)


end if


end do prepartic0   




group=exgroup
!******************************************************************



parti1:do i=1,npart

!Check if there is creation or destruction of molecules
if (group2(i,2).ne.group(i,2))then
call time_change

icount=icount+1
exit
end if

end do parti1




!######## for printing the group matrix every nth change.###########
!if (mod(t,tprint)==0) then 
!
!call printgroup
!icount=0
!end if


!##############################################################################################################

!potential=0.0d0
!dpotx=0.0d0
!dpoty=0.0d0
!dpotz=0.0d0
!
!
!
!call potential_total




!-------------write to file the geometries of the particles---------------------


if(mod(t,int(tprint,8))==0)  call write_geom  

!--------according to tbackup,  write to backup file the positions------

if(mod(t,int(tbackup,8))==0) call write_bk()



end do time

print*,'acabei'    

END SUBROUTINE verlet_alg

subroutine system_mem_usage(valueRSS)
implicit none
!  use ifport !if on intel compiler
integer, intent(out) :: valueRSS

character(len=200):: filename=' '
character(len=80) :: line
character(len=8)  :: pid_char=' '
integer :: pid,getpid
logical :: ifxst

valueRSS=-1    ! return negative number if not found

!--- get process ID

pid=getpid()
write(pid_char,'(I8)') pid
filename='/proc/'//trim(adjustl(pid_char))//'/status'

!--- read system file

inquire (file=filename,exist=ifxst)
if (.not.ifxst) then
write (*,*) 'system file does not exist'
return
endif

open(unit=100, file=filename, action='read')
do
read (100,'(a)',end=120) line
if (line(1:6).eq.'VmRSS:') then
read (line(7:),*) valueRSS
exit
endif
enddo
120 continue
close(100)

return
end subroutine system_mem_usage

