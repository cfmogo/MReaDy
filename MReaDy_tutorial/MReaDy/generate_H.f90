subroutine generate_H(igen,xgat,ygat,zgat,vxgat,vygat,vzgat)
use variables
use constants
use phys_parameters
implicit none

!-----for positioning---------
double precision ,intent(out):: xgat,ygat,zgat
double precision      ::       rxij,ryij,rzij,delta_r_atoms,distmini
double precision       ::      phi,cosphi,sinphi,costeta,sinteta
integer       :: igen,i
!----for velocities      

double precision       ::      modvelocity,Ma,vel
double precision,intent(out)       ::        vxgat,vygat,vzgat
double precision       ::       xrand(3)

! H molecule parameters
Ma=massH11*uma

!      Determine a random point inside a square box of lengh boxl
1       call random_number(xrand)

        xgat=xrand(1)*boxl
        ygat=xrand(2)*boxl
        zgat=xrand(3)*boxl
      
!!#####################################
!      print*,'altered generate_H.f90'
!      xgat=boxl*0.3d0
!      ygat=boxl*0.5d0
!      zgat=boxl*0.5d0
!!#####################################
 
      ! To asure that points are not overlapped
                          
        do i=1,igen-1

            !determine distance between points      
            call mod_dist(xgat,ygat,zgat,rx(i),ry(i),rz(i),rxij,ryij,rzij,delta_r_atoms)                 

            !minimum distance between the new point and remaining atoms 
            distmini=4E-10

            !if deltar_r_atoms is less than distmini
            if(delta_r_atoms.lt.distmini) then
                goto 1 ! generates a new position      
            end if            
        enddo

!---------------------------------------------------------
!------------------------/for the velocity/---------------
!---------------------------------------------------------

        ekin=0.0d0

        select case (temp_distribution)

            case (0)                  

                modvelocity=sqrt(3*rkbolt*temp/(Ma))

            case (1)

                call velocity_dist_H(temp,vel)
                modvelocity = vel !for H

            case default

                stop 'error in selection of Kinetic Energy distribution'

        end select

        call random_number(xrand)

        phi=2.0d0*pi*xrand(1)      
        cosphi=cos(phi)                  
        sinphi=sin(phi)

        costeta=1.0d0-2.0d0*xrand(2)
        sinteta=sqrt(1.0d0-costeta**2)      
        
        vxgat=modvelocity*sinteta*cosphi
        vygat=modvelocity*sinteta*sinphi
        vzgat=modvelocity*costeta
      
!###############################
!      print*,'altered generate_H'       
!      vxgat=modvelocity
!      vygat=0.0D0
!      vzgat=0.0D0
!###############################

end subroutine generate_H 


subroutine velocity_dist_H(temperature,v)
use constants
implicit none
       
double precision, intent(in):: temperature
double precision, intent(out):: v
double precision ::value,pv,mass_tot
double precision ::xrand(3)

mass_tot=massH11*uma

rej:    do 
            call random_number(xrand)
            value = (xrand(1)/(1.0d0-xrand(1)))
            pv=value*exp(-value)*1.0d0/sqrt((1.0d0-xrand(1))**(3)*xrand(1))

            if (xrand(2)*1.8d0.le.pv) exit rej      
      
        end do rej

        v=sqrt(2.0d0*rkbolt*temperature*value/mass_tot)
    
end subroutine velocity_dist_H






