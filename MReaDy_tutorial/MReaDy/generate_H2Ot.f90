subroutine generate_H2Ot(iatomo,xg,yg,zg,vxg,vyg,vzg)
use variables
use constants
use phys_parameters
implicit none

!-----for positioning---------
double precision ,intent(out):: xg(3),yg(3),zg(3)
double precision      ::       rxij,ryij,rzij,delta_r_atoms,distmini !,rij,r1eq,r2eq,r3eq
!double precision       ::      phi,cosphi,sinphi,costeta,sinteta,ang, cosang,sinang,sinang2      
integer       ::       i

!----for CM velocities      
double precision,intent(out)       ::        vxg(3),vyg(3),vzg(3)
integer,intent(in)        ::       iatomo
double precision       ::       xrand(3)

1       call random_number(xrand)

        xg(1)=xrand(1)*boxl
        yg(1)=xrand(2)*boxl
        zg(1)=xrand(3)*boxl
     
        if (iatomo.gt.1) then
                           
            do i=1,iatomo-1

                !determine distance between points      
                call mod_dist(xg(1),yg(1),zg(1),rx(i),ry(i),rz(i),rxij,ryij,rzij,delta_r_atoms)                 

                !minimum distance between the new point and remaining atoms 
                distmini=6E-10
 
                !if deltar_r_atoms is less than distmini
                if(delta_r_atoms.lt.distmini) then
                    goto 1 ! generates a new position      
                end if            
            enddo
        end if
      
        xg(2)=xg(1)+0.0000000D+00
        yg(2)=yg(1)+0.3961125D-09
        zg(2)=zg(1)+0.0000000D+00
        xg(3)=xg(1)+0.0000000D+00
        yg(3)=yg(1)+0.2988780D-09
        zg(3)=zg(1)+0.0000000D+00
 
        vxg(1)=0.0d0
        vyg(1)=0.0d0
        vzg(1)=0.0d0
        
        vxg(2)=0.0d0
        vyg(2)=0.0d0
        vzg(2)=0.0d0

        vxg(3)=0.0d0
        vyg(3)=0.0d0
        vzg(3)=0.0d0

end subroutine generate_H2Ot


