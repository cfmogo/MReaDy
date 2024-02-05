subroutine generate_O3(iatomo,xg,yg,zg,vxg,vyg,vzg)
use variables
use constants
use phys_parameters
implicit none

!-----for positioning---------
double precision ,intent(out):: xg(3),yg(3),zg(3)
double precision      ::       rxij,ryij,rzij,delta_r_atoms,distmini,r1eq,r2eq,r3eq
double precision       ::      ang, cosang,sinang      
integer       ::       i

!----for CM velocities      
double precision,intent(out)       ::        vxg(3),vyg(3),vzg(3)
integer,intent(in)        ::       iatomo
double precision       ::       xrand(3)

1   call random_number(xrand)

    xg(3)=xrand(1)*boxl
    yg(3)=xrand(2)*boxl
    zg(3)=xrand(3)*boxl

!       print*,'change in genO3'
!
!       xg(1)=0.2d0*boxl
!       yg(1)=0.5d0*boxl
!       zg(1)=0.5d0*boxl

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

    r1eq=2.3D0*5.2917721092d-11
    r2eq=2.3d0*5.2917721092d-11
    r3eq=sqrt(r1eq**2+r2eq**2-2.0d0*r1eq*r2eq*cos(117*acos(-1.0d0)/180.0d0))

    cosang=(r2eq**2+r3eq**2-r1eq**2)/(2.0d0*r2eq*r3eq)
    ang=acos(cosang)
    sinang=sin(ang)
      
    xg(2)=xg(3)
    yg(2)=yg(3)+r3eq
    zg(2)=zg(3)

    xg(1)=xg(3)
    yg(1)=yg(3)+cosang*r2eq
    zg(1)=zg(3)+sinang*r2eq

    vxg(1)=0.0d0
    vyg(1)=0.0d0
    vzg(1)=0.0d0
    
    vxg(2)=0.0d0
    vyg(2)=0.0d0
    vzg(2)=0.0d0

    vxg(3)=0.0d0
    vyg(3)=0.0d0
    vzg(3)=0.0d0

end subroutine generate_O3


