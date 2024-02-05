subroutine generate_H2O(iatomo,xg,yg,zg,vxg,vyg,vzg)
use variables
use constants
use phys_parameters
implicit none

!-----for positioning---------
double precision ,intent(out):: xg(3),yg(3),zg(3)
double precision      ::       rxij,ryij,rzij,delta_r_atoms,distmini
!double precision       ::      phi,cosphi,sinphi,costeta,sinteta,ang, cosang,sinang,sinang2      
double precision :: cosang,sinang,ang
integer       ::       i
!----for CM velocities      

double precision,intent(out)       ::        vxg(3),vyg(3),vzg(3)
integer,intent(in)        ::       iatomo
double precision       ::       xrand(3)

1      call random_number(xrand)
      
      xg(3)=xrand(1)*boxl
      yg(3)=xrand(2)*boxl 
      zg(3)=xrand(3)*boxl

!########################################
!      print*,'altered generate_H2O.f90'
!      xg(3)=0.5d0*boxl
!      yg(3)=0.5d0*boxl
!      zg(3)=0.5d0*boxl
!########################################
     ! To asure that points are not overlapped
   
      if (iatomo.gt.1) then
                           
            do i=1,iatomo-1
                !determine distance between points      
                call mod_dist(xg(3),yg(3),zg(3),rx(i),ry(i),rz(i),rxij,ryij,rzij,delta_r_atoms)                 

                !minimum distance between the new point and remaining atoms 
                distmini=6E-10

                !if deltar_r_atoms is less than distmini
                if(delta_r_atoms.lt.distmini) then
                    goto 1 ! generates a new position      
                end if            
            enddo
        end if
     
!     R1  H--H
!     R2  O--H
!     R3  O--H
 
!     the results should be:
!
!  r1, r2, r3
!   2.86083100000000            1.80979800000000              1.80979800000000     
!   1.513886678E-10            0.957703927E-10            0.957703927E-10
!   v 
! -0.370400689025461     

        cosang=(0.957703927E-10**2+0.957703927E-10**2-1.513886678E-10**2)/(2.0d0*0.957703927E-10*0.957703927E-10)
        ang=acos(cosang)
        sinang=sin(ang)
      
        xg(2)=xg(3)
        yg(2)=yg(3)+0.957703927E-10
        zg(2)=zg(3)

        xg(1)=xg(3)
        yg(1)=yg(3)+cosang*0.957703927E-10
        zg(1)=zg(3)+sinang*0.957703927E-10

        vxg(1)=0.0d0
        vyg(1)=0.0d0
        vzg(1)=0.0d0
        
        vxg(2)=0.0d0
        vyg(2)=0.0d0
        vzg(2)=0.0d0

        vxg(3)=0.0d0
        vyg(3)=0.0d0
        vzg(3)=0.0d0

end subroutine generate_H2O

