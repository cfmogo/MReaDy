subroutine generate_HO3(iatomo,xg,yg,zg,vxg,vyg,vzg)
use variables
use constants
use phys_parameters
implicit none

!-----for positioning---------
double precision ,intent(out):: xg(4),yg(4),zg(4)
double precision      ::       rxij,ryij,rzij,delta_r_atoms,distmini!,rij,r1eq,r2eq,r3eq,Odist
!double precision       ::      phi,cosphi,sinphi,costeta,sinteta,ang, cosang,sinang,sinang2,um      
integer       ::       i

!----for CM velocities      
double precision,intent(out)       ::        vxg(4),vyg(4),vzg(4)
integer,intent(in)        ::       iatomo
double precision       ::       xrand(3)

1   call random_number(xrand)

    xg(4)=xrand(1)*boxl
    yg(4)=xrand(2)*boxl
    zg(4)=xrand(3)*boxl
    r1te=xg(4)
    r2te=yg(4)
    r3te=zg(4)

    ! To asure that points are not overlapped
      
    if (iatomo.gt.1) then
                       
        do i=1,iatomo-1

            !determine distance between points      
            call mod_dist(xg(4),yg(4),zg(4),rx(i),ry(i),rz(i),rxij,ryij,rzij,delta_r_atoms)                 
            
            !minimum distance between the new point and remaining atoms 
            distmini=2E-10
            
            !if deltar_r_atoms is less than distmini
            if(delta_r_atoms.lt.distmini) then
                goto 1 ! generates a new position      
            end if            
        enddo
    end if

    xg(1)=xg(4)  +0.22527835006D-09
    yg(1)=yg(4)  +0.35088331180D-09
    zg(1)=zg(4)  +0.13106365680D-09
    xg(2)=xg(4)  -0.33181413378D-10
    yg(2)=yg(4)  +0.10391949220D-09
    zg(2)=zg(4)  +0.64957876500D-10
    xg(3)=xg(4)  +0.13690605257D-10
    yg(3)=yg(4)  +0.21443492890D-09
    zg(3)=zg(4)  +0.22663448100D-10

    vxg(1)=0.0d0
    vyg(1)=0.0d0
    vzg(1)=0.0d0

    vxg(2)=0.0d0
    vyg(2)=0.0d0
    vzg(2)=0.0d0

    vxg(3)=0.0d0
    vyg(3)=0.0d0
    vzg(3)=0.0d0

    vxg(4)=0.0d0
    vyg(4)=0.0d0
    vzg(4)=0.0d0

end subroutine generate_HO3


