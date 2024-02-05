subroutine generate_HO2(iatomo,xg,yg,zg,vxg,vyg,vzg)
use variables
use constants
use phys_parameters
implicit none

!-----for positioning---------
double precision ,intent(out):: xg(3),yg(3),zg(3)
double precision      ::       rxij,ryij,rzij,delta_r_atoms,distmini
double precision       ::      ang, cosang,sinang
integer       ::       i

!----for CM velocities
double precision,intent(out)       ::        vxg(3),vyg(3),vzg(3)
integer,intent(in)        ::       iatomo
double precision       ::       xrand(3)

1   call random_number(xrand)

    xg(2)=xrand(1)*boxl
    yg(2)=xrand(2)*boxl
    zg(2)=xrand(3)*boxl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       print*,'change in genHO2'
!
!       xg(2)=0.5d0*boxl
!       yg(2)=0.5d0*boxl
!       zg(2)=0.5d0*boxl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! To asure that points are not overlapped
    if (iatomo.gt.1) then

        do i=1,iatomo-1

            !determine distance between points
            call mod_dist(xg(2),yg(2),zg(2),rx(i),ry(i),rz(i),rxij,ryij,rzij,delta_r_atoms)

            !minimum distance between the new point and remaining atoms
            distmini=6E-10

            !if deltar_r_atoms is less than distmini
            if(delta_r_atoms.lt.distmini) then
                goto 1 ! generates a new position
            end if
        enddo
    end if

    cosang=(1.336848472e-10**2+9.836576367e-11**2-1.866296924e-10**2)/(2.0d0*1.336848472e-10*9.836576367e-11)
    ang=acos(cosang)
    sinang=sin(ang)

    xg(1)=xg(2)
    yg(1)=yg(2)+cosang*9.836576367e-11
    zg(1)=zg(2)+sinang*9.836576367e-11

    xg(3)=xg(2)
    yg(3)=yg(2)+1.336848472e-10
    zg(3)=zg(2)

    vxg(1)=0.0d0
    vyg(1)=0.0d0
    vzg(1)=0.0d0

    vxg(2)=0.0d0
    vyg(2)=0.0d0
    vzg(2)=0.0d0

    vxg(3)=0.0d0
    vyg(3)=0.0d0
    vzg(3)=0.0d0

end subroutine generate_HO2


