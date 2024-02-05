program mready

use constants
use variables
use chem_parameters
use phys_parameters
use sim_variables
use calc_parameters
use omp_lib
implicit none
integer :: igen,i,iatomo
double precision :: xd(2),yd(2),zd(2),vxd(2),vyd(2),vzd(2)
double precision :: xt(3),yt(3),zt(3),vxt(3),vyt(3),vzt(3)
double precision :: xq(4),yq(4),zq(4),vxq(4),vyq(4),vzq(4)
double precision :: xqat,yqat,zqat,vxqat,vyqat,vzqat
potential_change=0.0d0
ekinr=0.0d0


!call getenv('OMP_NUM_THREADS',cnthre)
!read(cnthre,'(i)'),nthre
!call omp_set_num_threads(nthre)
!!write(*,*)"#Threads outside the parallel section:",omp_get_num_threads()
!

!!______________________________________________________________________
!CALL SYSTEM_CLOCK(tc1,nb_ticks_sec)
!ri1=real(tc1)
!
!n_t_s=real(nb_ticks_sec)
!!______________________________________________________________________



!!______________________________________________________________________
!rti1=real(tc1)
!CALL SYSTEM_CLOCK(tc1)
!rtf1=real(tc1)
!dti1=(rtf1-rti1)/n_t_s
!
!write(*,fmt='(A8,F12.6)') "read",dti1
!!______________________________________________________________________






potential_change=0.0d0


      !Determinacao de variaveis
      call read_input_file()


      if (bk_new.ge.1) go to 456
 !     print*,'after go to 456'



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

iatomo=0
cname=' '

      if (natoms_H.eq.0) then
      print*, 'atomos de H      =>           0'
      else
      print*, 'atomos de H      =>',natoms_H

H:      do igen= 1,natoms_H

      i=iatomo+igen

       call generate_H(i,xqat,yqat,zqat,vxqat,vyqat,vzqat)

      rx(i)=xqat
      ry(i)=yqat
      rz(i)=zqat


      vx(i)=vxqat
      vy(i)=vyqat
      vz(i)=vzqat


      cname(i)='Hy'

      mass(i)=massH11*uma
      
      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)

      end do H
      end if
      iatomo=iatomo+natoms_H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (natoms_O.eq.0) then
      print*, 'atomos de O      =>           0'
      else
      print*, 'atomos de O      =>',natoms_O

O:      do igen= 1,natoms_O
        i=iatomo+igen


       call generate_O(i,xqat,yqat,zqat,vxqat,vyqat,vzqat)

      rx(i)=xqat
      ry(i)=yqat
      rz(i)=zqat


      vx(i)=vxqat
      vy(i)=vyqat
      vz(i)=vzqat


      cname(i)='Ox'

      mass(i)=massO*uma
      
      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)     

      end do O
      end if

      iatomo=iatomo+natoms_O

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_HH.eq.0) then
      print*, 'moleculas de H2  =>           0'
      else
      print*, 'moleculas de H2 =>',nmolec_HH
HH:      do igen=1, nmolec_HH
      i=iatomo+igen*2-1

       call generate_HH(i,xd,yd,zd,vxd,vyd,vzd)
      rx(i)=xd(1)
      ry(i)=yd(1)
      rz(i)=zd(1)
      rx(i+1)=xd(2)
      ry(i+1)=yd(2)
      rz(i+1)=zd(2)

      vx(i)=vxd(1)
      vy(i)=vyd(1)
      vz(i)=vzd(1)
      vx(i+1)=vxd(2)
      vy(i+1)=vyd(2)
      vz(i+1)=vzd(2)
      cname(i)='Hy'
      cname(i+1)='Hy'
      mass(i)=massH11*uma
      mass(i+1)=massH11*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)


      end do HH
      end if

      iatomo=iatomo+nmolec_HH*2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_HO.eq.0) then
      print*, 'moleculas de HO  =>           0'
      else
      print*, 'moleculas de HO  =>',nmolec_HO

HO:      do igen= 1,nmolec_HO
      i=iatomo+igen*2-1

      call generate_HO(i,xd,yd,zd,vxd,vyd,vzd)
      rx(i)=xd(1)
      ry(i)=yd(1)
      rz(i)=zd(1)
      rx(i+1)=xd(2)
      ry(i+1)=yd(2)
      rz(i+1)=zd(2)

      vx(i)=vxd(1)
      vy(i)=vyd(1)
      vz(i)=vzd(1)
      vx(i+1)=vxd(2)
      vy(i+1)=vyd(2)
      vz(i+1)=vzd(2)
      cname(i)='Hy'
      cname(i+1)='Ox'
      mass(i)=massH11*uma
      mass(i+1)=massO*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)

      end do HO
      end if

      iatomo=iatomo+nmolec_HO*2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_OOsing.eq.0) then
      print*, 'moleculas de O2s =>           0'
      else
      print*, 'moleculas de O2s =>',nmolec_OOsing


OOs:do igen= 1,nmolec_OOsing

      i=iatomo+igen*2-1

      call generate_OOs(i,xd,yd,zd,vxd,vyd,vzd)

      rx(i)=xd(1)
      ry(i)=yd(1)
      rz(i)=zd(1)
      rx(i+1)=xd(2)
      ry(i+1)=yd(2)
      rz(i+1)=zd(2)

      vx(i)=vxd(1)
      vy(i)=vyd(1)
      vz(i)=vzd(1)
      vx(i+1)=vxd(2)
      vy(i+1)=vyd(2)
      vz(i+1)=vzd(2)
      cname(i)='Ox'
      cname(i+1)='Ox'
      mass(i)=massO*uma
      mass(i+1)=massO*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)

      end do OOs
      end if

      iatomo=iatomo+nmolec_OOsing*2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_OOtrip.eq.0) then
      print*, 'moleculas de O2t  =>          0'
      else
      print*, 'moleculas de O2t=>',nmolec_OOtrip

OOt:      do igen=1, nmolec_OOtrip
      i=iatomo + igen*2-1

      call generate_OOt(i,xd,yd,zd,vxd,vyd,vzd)

      rx(i)=xd(1)
      ry(i)=yd(1)
      rz(i)=zd(1)
      rx(i+1)=xd(2)
      ry(i+1)=yd(2)
      rz(i+1)=zd(2)

      vx(i)=vxd(1)
      vy(i)=vyd(1)
      vz(i)=vzd(1)
      vx(i+1)=vxd(2)
      vy(i+1)=vyd(2)
      vz(i+1)=vzd(2)
      cname(i)='Ox'
      cname(i+1)='Ox'
      mass(i)=massO*uma
      mass(i+1)=massO*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)

      end do OOt
      end if

      iatomo=iatomo + nmolec_OOtrip*2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (nmolec_H3.eq.0) then
      print*, 'moleculas de H3  =>           0'
      else
      print*, 'moleculas de H3 =>',nmolec_H3


H3:      do igen= 1,nmolec_H3
      i=iatomo+igen*3-2

       call generate_H3(i,xt,yt,zt,vxt,vyt,vzt)

      rx(i)=xt(1)
      ry(i)=yt(1)
      rz(i)=zt(1)

      rx(i+1)=xt(2)
      ry(i+1)=yt(2)
      rz(i+1)=zt(2)

      rx(i+2)=xt(3)
      ry(i+2)=yt(3)
      rz(i+2)=zt(3)

      vx(i)=vxt(1)
      vy(i)=vyt(1)
      vz(i)=vzt(1)

      vx(i+1)=vxt(2)
      vy(i+1)=vyt(2)
      vz(i+1)=vzt(2)

      vx(i+2)=vxt(3)
      vy(i+2)=vyt(3)
      vz(i+2)=vzt(3)



      cname(i)='Hy'
      cname(i+1)='Hy'
      cname(i+2)='Hy'

      mass(i)=massH11*uma
      mass(i+1)=massH11*uma
      mass(i+2)=massH11*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)



      end do H3
      end if

      iatomo=iatomo+nmolec_H3*3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_H2O.eq.0) then
      print*, 'moleculas de H2O  =>          0'
      else
      print*, 'moleculas de H2O =>',nmolec_H2O


H2O:      do igen= 1,nmolec_H2O
     i= iatomo+igen*3-2

       call generate_H2O(i,xt,yt,zt,vxt,vyt,vzt)

      rx(i)=xt(1)
      ry(i)=yt(1)
      rz(i)=zt(1)

      rx(i+1)=xt(2)
      ry(i+1)=yt(2)
      rz(i+1)=zt(2)

      rx(i+2)=xt(3)
      ry(i+2)=yt(3)
      rz(i+2)=zt(3)

      vx(i)=vxt(1)
      vy(i)=vyt(1)
      vz(i)=vzt(1)

      vx(i+1)=vxt(2)
      vy(i+1)=vyt(2)
      vz(i+1)=vzt(2)

      vx(i+2)=vxt(3)
      vy(i+2)=vyt(3)
      vz(i+2)=vzt(3)



      cname(i)='Hy'
      cname(i+1)='Hy'
      cname(i+2)='Ox'

      mass(i)=massH11*uma
      mass(i+1)=massH11*uma
      mass(i+2)=massO*uma


      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)

      end do H2O
      end if

      iatomo=iatomo+nmolec_H2O*3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_H2Ot.eq.0) then
      print*, 'moleculas de H2O trip  =>     0'
      else
      print*, 'moleculas de H2O trip =>',nmolec_H2Ot


H2Ot:      do igen= 1,nmolec_H2Ot
      i=iatomo+igen*3-2

       call generate_H2Ot(i,xt,yt,zt,vxt,vyt,vzt)

      rx(i)=xt(1)
      ry(i)=yt(1)
      rz(i)=zt(1)

      rx(i+1)=xt(2)
      ry(i+1)=yt(2)
      rz(i+1)=zt(2)

      rx(i+2)=xt(3)
      ry(i+2)=yt(3)
      rz(i+2)=zt(3)

      vx(i)=vxt(1)
      vy(i)=vyt(1)
      vz(i)=vzt(1)

      vx(i+1)=vxt(2)
      vy(i+1)=vyt(2)
      vz(i+1)=vzt(2)

      vx(i+2)=vxt(3)
      vy(i+2)=vyt(3)
      vz(i+2)=vzt(3)



      cname(i)='Hy'
      cname(i+1)='Hy'
      cname(i+2)='Ox'

       mass(i)=massH11*uma
       mass(i+1)=massH11*uma
       mass(i+2)=massO*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)

 
 
 
       end do H2Ot
      end if

      iatomo=iatomo+nmolec_H2Ot*3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_HO2.eq.0) then
      print*, 'moleculas de HO2 =>           0'
      else
      print*, 'moleculas de HO2 =>',nmolec_HO2

HO2:      do igen= 1,nmolec_HO2
     i=iatomo+igen*3-2

       call generate_HO2(i,xt,yt,zt,vxt,vyt,vzt)

      rx(i)=xt(1)
      ry(i)=yt(1)
      rz(i)=zt(1)

      rx(i+1)=xt(2)
      ry(i+1)=yt(2)
      rz(i+1)=zt(2)

      rx(i+2)=xt(3)
      ry(i+2)=yt(3)
      rz(i+2)=zt(3)

      vx(i)=vxt(1)
      vy(i)=vyt(1)
      vz(i)=vzt(1)

      vx(i+1)=vxt(2)
      vy(i+1)=vyt(2)
      vz(i+1)=vzt(2)

      vx(i+2)=vxt(3)
      vy(i+2)=vyt(3)
      vz(i+2)=vzt(3)


      cname(i)='Hy'
      cname(i+1)='Ox'
      cname(i+2)='Ox'

      mass(i)=massH11*uma
      mass(i+1)=massO*uma
      mass(i+2)=massO*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)


      end do HO2
      end if

      iatomo=iatomo+nmolec_HO2*3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_O3.eq.0) then
      print*, 'moleculas de O3  =>           0'
      else
      print*, 'moleculas de O3 =>',nmolec_O3


O3:      do igen= 1,nmolec_O3
      i=iatomo+igen*3-2

       call generate_O3(i,xt,yt,zt,vxt,vyt,vzt)

      rx(i)=xt(1)
      ry(i)=yt(1)
      rz(i)=zt(1)

      rx(i+1)=xt(2)
      ry(i+1)=yt(2)
      rz(i+1)=zt(2)

      rx(i+2)=xt(3)
      ry(i+2)=yt(3)
      rz(i+2)=zt(3)

      vx(i)=vxt(1)
      vy(i)=vyt(1)
      vz(i)=vzt(1)

      vx(i+1)=vxt(2)
      vy(i+1)=vyt(2)
      vz(i+1)=vzt(2)

      vx(i+2)=vxt(3)
      vy(i+2)=vyt(3)
      vz(i+2)=vzt(3)



      cname(i)='Ox'
      cname(i+1)='Ox'
      cname(i+2)='Ox'

      mass(i)=massO*uma
      mass(i+1)=massO*uma
      mass(i+2)=massO*uma


      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)




      end do O3
      end if

      iatomo=iatomo+nmolec_O3*3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_H4.eq.0) then
      print*, 'moleculas de H4  =>           0'
      else
      print*, 'moleculas de H4 =>',nmolec_H4


H4:      do igen= 1,nmolec_H4
      i=iatomo+igen*4-3

       call generate_H4(i,xq,yq,zq,vxq,vyq,vzq)

      rx(i)=xq(1)
      ry(i)=yq(1)
      rz(i)=zq(1)

      rx(i+1)=xq(2)
      ry(i+1)=yq(2)
      rz(i+1)=zq(2)

      rx(i+2)=xq(3)
      ry(i+2)=yq(3)
      rz(i+2)=zq(3)

      rx(i+3)=xq(4)
      ry(i+3)=yq(4)
      rz(i+3)=zq(4)


      vx(i)=vxq(1)
      vy(i)=vyq(1)
      vz(i)=vzq(1)

      vx(i+1)=vxq(2)
      vy(i+1)=vyq(2)
      vz(i+1)=vzq(2)

      vx(i+2)=vxq(3)
      vy(i+2)=vyq(3)
      vz(i+2)=vzq(3)

      vx(i+3)=vxq(4)
      vy(i+3)=vyq(4)
      vz(i+3)=vzq(4)

      cname(i)='Hy'
      cname(i+1)='Hy'
      cname(i+2)='Hy'
      cname(i+3)='Hy'

      mass(i)=massH11*uma
      mass(i+1)=massH11*uma
      mass(i+2)=massH11*uma
      mass(i+3)=massH11*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)
      ekinr=ekinr+0.5d0*mass(i+3)*(vx(i+3)**2+vy(i+3)**2+vz(i+3)**2)




      end do H4
      end if


      iatomo=iatomo+nmolec_H4*4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_H3O.eq.0) then
      print*, 'moleculas de H3O  =>          0'
      else
      print*, 'moleculas de H3O =>',nmolec_H3O


H3O:      do igen= 1,nmolec_H3O

       i=iatomo+igen*4-3


       call generate_H3O(i,xq,yq,zq,vxq,vyq,vzq)

      rx(i)=xq(1)
      ry(i)=yq(1)
      rz(i)=zq(1)

      rx(i+1)=xq(2)
      ry(i+1)=yq(2)
      rz(i+1)=zq(2)

      rx(i+2)=xq(3)
      ry(i+2)=yq(3)
      rz(i+2)=zq(3)

      rx(i+3)=xq(4)
      ry(i+3)=yq(4)
      rz(i+3)=zq(4)


      vx(i)=vxq(1)
      vy(i)=vyq(1)
      vz(i)=vzq(1)

      vx(i+1)=vxq(2)
      vy(i+1)=vyq(2)
      vz(i+1)=vzq(2)

      vx(i+2)=vxq(3)
      vy(i+2)=vyq(3)
      vz(i+2)=vzq(3)

      vx(i+3)=vxq(4)
      vy(i+3)=vyq(4)
      vz(i+3)=vzq(4)


      cname(i)='Hy'
      cname(i+1)='Hy'
      cname(i+2)='Hy'
      cname(i+3)='Ox'

       mass(i)=massH11*uma
       mass(i+1)=massH11*uma
       mass(i+2)=massH11*uma
       mass(i+3)=massO*uma


      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)
      ekinr=ekinr+0.5d0*mass(i+3)*(vx(i+3)**2+vy(i+3)**2+vz(i+3)**2)



 
 
       end do H3O
      end if

      iatomo=iatomo+nmolec_H3O*4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_H2O2s.eq.0) then
      print*, 'moleculas de H2O2  =>         0'
      else
      print*, 'moleculas de H2O2 =>',nmolec_H2O2s


H2O2s:      do igen= 1,nmolec_H2O2s
      i=iatomo+igen*4-3

       call generate_H2O2s(i,xq,yq,zq,vxq,vyq,vzq)

      rx(i)=xq(1)
      ry(i)=yq(1)
      rz(i)=zq(1)

      rx(i+1)=xq(2)
      ry(i+1)=yq(2)
      rz(i+1)=zq(2)

      rx(i+2)=xq(3)
      ry(i+2)=yq(3)
      rz(i+2)=zq(3)

      rx(i+3)=xq(4)
      ry(i+3)=yq(4)
      rz(i+3)=zq(4)


      vx(i)=vxq(1)
      vy(i)=vyq(1)
      vz(i)=vzq(1)

      vx(i+1)=vxq(2)
      vy(i+1)=vyq(2)
      vz(i+1)=vzq(2)

      vx(i+2)=vxq(3)
      vy(i+2)=vyq(3)
      vz(i+2)=vzq(3)

      vx(i+3)=vxq(4)
      vy(i+3)=vyq(4)
      vz(i+3)=vzq(4)


      cname(i)='Hy'
      cname(i+1)='Hy'
      cname(i+2)='Ox'
      cname(i+3)='Ox'

      mass(i)=massH11*uma
      mass(i+1)=massH11*uma
      mass(i+2)=massO*uma
      mass(i+3)=massO*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)
      ekinr=ekinr+0.5d0*mass(i+3)*(vx(i+3)**2+vy(i+3)**2+vz(i+3)**2)


      end do H2O2s
      end if

      iatomo=iatomo+nmolec_H2O2s*4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_H2O2t.eq.0) then
      print*, 'moleculas de H2O2t  =>        0'
      else
      print*, 'moleculas de H2O2t =>',nmolec_H2O2t


H2O2t:do igen= 1,nmolec_H2O2t
     i= iatomo+igen*4-3

       call generate_H2O2t(i,xq,yq,zq,vxq,vyq,vzq)

      rx(i)=xq(1)
      ry(i)=yq(1)
      rz(i)=zq(1)

      rx(i+1)=xq(2)
      ry(i+1)=yq(2)
      rz(i+1)=zq(2)

      rx(i+2)=xq(3)
      ry(i+2)=yq(3)
      rz(i+2)=zq(3)

      rx(i+3)=xq(4)
      ry(i+3)=yq(4)
      rz(i+3)=zq(4)


      vx(i)=vxq(1)
      vy(i)=vyq(1)
      vz(i)=vzq(1)

      vx(i+1)=vxq(2)
      vy(i+1)=vyq(2)
      vz(i+1)=vzq(2)

      vx(i+2)=vxq(3)
      vy(i+2)=vyq(3)
      vz(i+2)=vzq(3)

      vx(i+3)=vxq(4)
      vy(i+3)=vyq(4)
      vz(i+3)=vzq(4)


      cname(i)='Hy'
      cname(i+1)='Hy'
      cname(i+2)='Ox'
      cname(i+3)='Ox'

      mass(i)=massH11*uma
      mass(i+1)=massH11*uma
      mass(i+2)=massO*uma
      mass(i+3)=massO*uma

      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)
      ekinr=ekinr+0.5d0*mass(i+3)*(vx(i+3)**2+vy(i+3)**2+vz(i+3)**2)

      end do H2O2t
      end if

      iatomo=iatomo+nmolec_H2O2t*4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_HO3.eq.0) then
      print*, 'moleculas de HO3  =>          0'
      else
      print*, 'moleculas de HO3 =>',nmolec_HO3


HO3:      do igen= 1,nmolec_HO3
      i=iatomo+igen*4-3

       call generate_HO3(i,xq,yq,zq,vxq,vyq,vzq)

      rx(i)=xq(1)
      ry(i)=yq(1)
      rz(i)=zq(1)

      rx(i+1)=xq(2)
      ry(i+1)=yq(2)
      rz(i+1)=zq(2)

      rx(i+2)=xq(3)
      ry(i+2)=yq(3)
      rz(i+2)=zq(3)

      rx(i+3)=xq(4)
      ry(i+3)=yq(4)
      rz(i+3)=zq(4)


      vx(i)=vxq(1)
      vy(i)=vyq(1)
      vz(i)=vzq(1)

      vx(i+1)=vxq(2)
      vy(i+1)=vyq(2)
      vz(i+1)=vzq(2)

      vx(i+2)=vxq(3)
      vy(i+2)=vyq(3)
      vz(i+2)=vzq(3)

      vx(i+3)=vxq(4)
      vy(i+3)=vyq(4)
      vz(i+3)=vzq(4)


   cname(i)='Hy'
   cname(i+1)='Ox'
   cname(i+2)='Ox'
   cname(i+3)='Ox'

    mass(i)=massH11*uma
    mass(i+1)=massO*uma
    mass(i+2)=massO*uma
    mass(i+3)=massO*uma


    ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
    ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
    ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)
    ekinr=ekinr+0.5d0*mass(i+3)*(vx(i+3)**2+vy(i+3)**2+vz(i+3)**2)


    end do HO3
      end if

      iatomo=iatomo+nmolec_HO3*4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (nmolec_O4.eq.0) then
      print*, 'moleculas de O4  =>           0'
      else
      print*, 'moleculas de O4 =>',nmolec_O4


O4:      do igen= 1,nmolec_O4
      i=iatomo+igen*4-3

       call generate_O4(i,xq,yq,zq,vxq,vyq,vzq)

      rx(i)=xq(1)
      ry(i)=yq(1)
      rz(i)=zq(1)

      rx(i+1)=xq(2)
      ry(i+1)=yq(2)
      rz(i+1)=zq(2)

      rx(i+2)=xq(3)
      ry(i+2)=yq(3)
      rz(i+2)=zq(3)

      rx(i+3)=xq(4)
      ry(i+3)=yq(4)
      rz(i+3)=zq(4)


      vx(i)=vxq(1)
      vy(i)=vyq(1)
      vz(i)=vzq(1)

      vx(i+1)=vxq(2)
      vy(i+1)=vyq(2)
      vz(i+1)=vzq(2)

      vx(i+2)=vxq(3)
      vy(i+2)=vyq(3)
      vz(i+2)=vzq(3)

      vx(i+3)=vxq(4)
      vy(i+3)=vyq(4)
      vz(i+3)=vzq(4)


      cname(i)='Ox'
      cname(i+1)='Ox'
      cname(i+2)='Ox'
      cname(i+3)='Ox'

      mass(i)=massO*uma
      mass(i+1)=massO*uma
      mass(i+2)=massO*uma
      mass(i+3)=massO*uma


      ekinr=ekinr+0.5d0*mass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
      ekinr=ekinr+0.5d0*mass(i+1)*(vx(i+1)**2+vy(i+1)**2+vz(i+1)**2)
      ekinr=ekinr+0.5d0*mass(i+2)*(vx(i+2)**2+vy(i+2)**2+vz(i+2)**2)
      ekinr=ekinr+0.5d0*mass(i+3)*(vx(i+3)**2+vy(i+3)**2+vz(i+3)**2)






      end do O4
      end if

      iatomo=iatomo+nmolec_O4*4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
















!-------------------------------//-----------------------------------------------
      ! Determination and printing of the initial positions (new and old) and velocities
      !-------------------------------//-----------------------------------------------
      
!rti1=real(tc1)
!CALL SYSTEM_CLOCK(tc1)
!rtf1=real(tc1)
!dti1=(rtf1-rti1)/n_t_s
!
!write(*,fmt='(A8,F12.6)') "gera",dti1
!!______________________________________________________________________

      rmass=1/mass


!Saving constant parameters for backup
open(unit=20,file='out/ini.bk',access='stream')

write(unit=20)npart,temp,boxl,out_bound,in_bound,rl_cut_off,dt,cname,mass
close(20)








      call update_positions


       !---------------------------------------------------------------------------------



      !-------------Determination of the potential-----------------
      !---------------------------//-------------------------------


!            stop'problema'


456 continue



 


      call write_geom
      call time_change

!      call printgroup






      !-------------initializing molecules with verlet algorithm---------
      !-------------Velocities and directions have to be save -----------
      !-------------for later use after evolution of the system
      !---------------------------//-----------------------------------


      !-------------for moving molecules with verlet algorithm---------
      !---------------------------//-----------------------------------

!!______________________________________________________________________
!rti1=real(tc1)
!CALL SYSTEM_CLOCK(tc1)
!rtf1=real(tc1)
!dti1=(rtf1-rti1)/n_t_s
!
!!______________________________________________________________________

 
       call verlet_alg()

!!______________________________________________________________________
!rti1=rtf1
!CALL SYSTEM_CLOCK(tc1)
!rtf1=real(tc1)
!dti1=(rtf1-rti1)/n_t_s
!
!write(*,fmt='(A8,F12.6)') "Verlet",dti1
!!______________________________________________________________________

      
       
      
!!______________________________________________________________________
!CALL SYSTEM_CLOCK(tc1)
!rf1=real(tc1)
!dti1=(rf1-ri1)/n_t_s
!
!write(*,fmt='(A8,F12.6)') "total",dti1
!!______________________________________________________________________




      close(10)


 print*,'Program finished'

      
end


subroutine update_positions()
            use constants
            use variables
            use phys_parameters
            use sim_variables

            implicit none
            
            !Just saving the original rx, ry, rz positions
            ax=rx
            ay=ry
            az=rz

            rx=rx-vx*dt
            ry=ry-vy*dt
            rz=rz-vz*dt


           potential=0.0d0
           dpotx=0.0d0
           dpoty=0.0d0
           dpotz=0.0d0

           call potential_total


           ! retreiving the original r_xyz values
            rx=ax
            ry=ay            
            rz=az



            !Now the real acceleration
            ax = - dpotx * rmass       
            ay = - dpoty * rmass
            az = - dpotz * rmass

!            do i=1, size(dpotx)
!            print*,dpotx(i),dpoty(i),dpotz(i)
!            end do

!           print*,'rx'

            ! relocation of atoms inside the box
            rx = rx - boxl*floor(rx/boxl)
            ry = ry - boxl*floor(ry/boxl)
            rz = rz - boxl*floor(rz/boxl)

!            do i=1, size(rx)
!            print*,rx(i),ry(i),rz(i)
!            end do

end subroutine update_positions

!subroutine update_positions()
!            use constants
!            use variables
!            implicit none
!
!
!
!            rxold=rx-vx*dt
!            ryold=ry-vy*dt
!            rzold=rz-vz*dt
!
!            rxold = rxold - boxl*floor(rx/boxl)
!            rx = rx - boxl*floor(rx/boxl)
!            ryold = ryold - boxl*floor(ry/boxl)
!            ry = ry - boxl*floor(ry/boxl)
!            rzold = rzold - boxl*floor(rz/boxl)
!            rz = rz - boxl*floor(rz/boxl)
!
!end subroutine update_positions




















