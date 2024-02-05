subroutine potential_total!(i)
   use variables
   use phys_parameters
   use sim_variables
   use constants
   use omp_lib
   use kdtree2_module

   implicit none
!integer,intent(in) :: i
   integer ::ip, j, k, p, i

   double precision      :: rij, rxij, ryij, rzij, V, dV
   integer      :: gi3, gi4, gi5, gi6, gij, gkl
   double precision      :: r1, r2, r3, r4, r5, r6
   double precision      :: r1xij, r1yij, r1zij, r2xij, r2yij, r2zij, r3xij, r3yij, r3zij
   double precision      :: r4xij, r4yij, r4zij, r5xij, r5yij, r5zij, r6xij, r6yij, r6zij
   double precision      :: dvdr1, dvdr2, dvdr3, dvdr4, dvdr5, dvdr6
   double precision, dimension(12) :: dvdabccc
   double precision, dimension(6) ::r_4atoms, dvdr_4atoms
   integer :: ihh, ioo, iho

   integer :: nfound, soma

   double precision :: xyz(3, 1)

   integer :: npart_full

   integer :: refcount(3, npart)
   logical :: mask1(npart,17),mask2(3,npart)
  
   potential = 0.0d0

!construct array with group information for each atom
!$OMP PARALLEL DO

   do i = 1, npart
      rl_part(group(i, 3)) = i
      rl_part(group(i, 4)) = i
      rl_part(group(i, 5)) = i
      rl_part(group(i, 6)) = i
   end do
!$OMP END PARALLEL DO

   potential = 0.0d0
   dpotx = 0.0d0
   dpoty = 0.0d0
   dpotz = 0.0d0
   potentiald = 0.0d0

   dpotxd = 0.0d0
   dpotyd = 0.0d0
   dpotzd = 0.0d0
   rxij = 0.0d0
   ryij = 0.0d0
   rzij = 0.0d0
   rij = 0.0d0



!For getting the logical matrix for each of the present species

   do j = 1, 17
      do i = 1, npart
         mask1(i, j) = group(i, 2) /= j
      end do
   end do

!$OMP SECTIONS   private(gi3,gi4,rxij,ryij,rzij,rij,V,dv)  &
!$OMP REDUCTION(+:potential)

!$OMP SECTION

!!$OMP PARALLEL DO private(gi3,gi4,rxij,ryij,rzij,rij,V,dv) SHARED(rx,ry,rzi,dpotx,dpoty,dpotz) &
!!$OMP REDUCTION(+:potential)

   h2do: do i = 1, npart

      if (mask1(i, 1)) cycle h2do

      gi3 = group(i, 3)
      gi4 = group(i, 4)

      call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

      call pothhs_f1(rij, V, dV)

      dpotx(gi3) = dpotx(gi3) + dV*rxij/rij
      dpoty(gi3) = dpoty(gi3) + dV*ryij/rij
      dpotz(gi3) = dpotz(gi3) + dV*rzij/rij

      dpotx(gi4) = dpotx(gi4) - dV*rxij/rij
      dpoty(gi4) = dpoty(gi4) - dV*ryij/rij
      dpotz(gi4) = dpotz(gi4) - dV*rzij/rij

      potential = potential + V

   end do h2do
!!$OMP end parallel do

!$OMP SECTION

!!$OMP PARALLEL DO private(gi3,gi4,rxij,ryij,rzij,rij,V,dv) SHARED(rx,ry,rzi,dpotx,dpoty,dpotz) &
!!$OMP REDUCTION(+:potential)

   h2t: do i = 1, npart

      if (mask1(i, 14)) cycle h2t

      gi3 = group(i, 3)
      gi4 = group(i, 4)

      call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

      call pothht_f14(rij, V, dV)

      dpotx(gi3) = dpotx(gi3) + dV*rxij/rij
      dpoty(gi3) = dpoty(gi3) + dV*ryij/rij
      dpotz(gi3) = dpotz(gi3) + dV*rzij/rij

      dpotx(gi4) = dpotx(gi4) - dV*rxij/rij
      dpoty(gi4) = dpoty(gi4) - dV*ryij/rij
      dpotz(gi4) = dpotz(gi4) - dV*rzij/rij

      potential = potential + V

   end do h2t
!!$OMP end parallel do

!$OMP SECTION

!!$OMP PARALLEL DO private(gi3,gi4,rxij,ryij,rzij,rij,V,dv) SHARED(rx,ry,rzi,dpotx,dpoty,dpotz) &
!!$OMP REDUCTION(+:potential)

   o2pdo: do i = 1, npart
      if (mask1(i, 2)) cycle o2pdo

      gi3 = group(i, 3)
      gi4 = group(i, 4)

      call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

      call potoop_f2(rij, V, dV)

      dpotx(gi3) = dpotx(gi3) + dV*rxij/rij
      dpoty(gi3) = dpoty(gi3) + dV*ryij/rij
      dpotz(gi3) = dpotz(gi3) + dV*rzij/rij

      dpotx(gi4) = dpotx(gi4) - dV*rxij/rij
      dpoty(gi4) = dpoty(gi4) - dV*ryij/rij
      dpotz(gi4) = dpotz(gi4) - dV*rzij/rij

      potential = potential + V
!print*,'V',V/con
   end do o2pdo
!!$OMP end parallel do

!$OMP SECTION

!!$OMP PARALLEL DO private(gi3,gi4,rxij,ryij,rzij,rij,V,dv) SHARED(rx,ry,rzi,dpotx,dpoty,dpotz) &
!!$OMP REDUCTION(+:potential)
   o2ddo: do i = 1, npart
      if (mask1(i, 15)) cycle o2ddo

      gi3 = group(i, 3)
      gi4 = group(i, 4)

      call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

      call potood_f22(rij, V, dV)

      dpotx(gi3) = dpotx(gi3) + dV*rxij/rij
      dpoty(gi3) = dpoty(gi3) + dV*ryij/rij
      dpotz(gi3) = dpotz(gi3) + dV*rzij/rij

      dpotx(gi4) = dpotx(gi4) - dV*rxij/rij
      dpoty(gi4) = dpoty(gi4) - dV*ryij/rij
      dpotz(gi4) = dpotz(gi4) - dV*rzij/rij

      potential = potential + V

   end do o2ddo
!!$OMP end parallel do

!$OMP SECTION

!!$OMP PARALLEL DO private(gi3,gi4,rxij,ryij,rzij,rij,V,dv) SHARED(rx,ry,rz,dpotx,dpoty,dpotz) &
!!$OMP REDUCTION(+:potential)
   ohpdo: do i = 1, npart
      if (mask1(i, 3)) cycle ohpdo

      gi3 = group(i, 3)
      gi4 = group(i, 4)

      call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

      call potohp_f3(rij, V, dV)

      dpotx(gi3) = dpotx(gi3) + dV*rxij/rij
      dpoty(gi3) = dpoty(gi3) + dV*ryij/rij
      dpotz(gi3) = dpotz(gi3) + dV*rzij/rij

      dpotx(gi4) = dpotx(gi4) - dV*rxij/rij
      dpoty(gi4) = dpoty(gi4) - dV*ryij/rij
      dpotz(gi4) = dpotz(gi4) - dV*rzij/rij

      potential = potential + V

   end do ohpdo
!!$OMP end parallel do

!$OMP END SECTIONS

!!$OMP PARALLEL DO private(gi3,gi4,gi5,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!!$OMP REDUCTION(+:potential)
   ho2_2p: do i = 1, npart
      if (mask1(i, 4)) cycle ho2_2p

!The positions are predefined in the matrix "group" with H atoms first and O later

      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! O atom
      gi5 = group(i, 5)      ! O atom

      call vho2_f4_switch(i, V, dvdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      potential = potential + V

   end do ho2_2p
!!$OMP end parallel do

!for the h2o(X1A') surface

!$OMP PARALLEL DO private(gi3,gi4,gi5,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!$OMP REDUCTION(+:potential)
   h2os: do i = 1, npart
      if (mask1(i, 5)) cycle h2os

!The positions are predefined in the matrix "group" with H atoms first and O later
      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! H atom
      gi5 = group(i, 5)      ! O atom

      call vh2o_f5_switch(i, V, dvdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      potential = potential + V

   end do h2os
!$OMP end parallel do

!$OMP PARALLEL DO private(gi3,gi4,gi5,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!$OMP REDUCTION(+:potential)
   h2ot: do i = 1, npart
      if (mask1(i, 6)) cycle h2ot

!The positions are predefined in the matrix "group" with H atoms first and O later
      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! H atom
      gi5 = group(i, 5)      ! O atom

      call vh2ot_f6_switch(i, V, dvdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      potential = potential + V

   end do h2ot
!$OMP end parallel do

!$OMP PARALLEL DO private(gi3,gi4,gi5,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!$OMP REDUCTION(+:potential)
   o3_f7: do i = 1, npart
      if (mask1(i, 7)) cycle o3_f7

!The positions are predefined in the matrix "group" with H atoms first and O later
      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! H atom
      gi5 = group(i, 5)      ! O atom

      call VO3_f7_switch(i, V, dvdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      potential = potential + V

   end do o3_f7
!$OMP end parallel do

!for the H3 surface
!$OMP PARALLEL DO private(gi3,gi4,gi5,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!$OMP REDUCTION(+:potential)
   h3s: do i = 1, npart
      if (mask1(i, 8)) cycle h3s

!H3(1)
      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! H atom
      gi5 = group(i, 5)      ! H atom

      call bkmp2_f8_switch(i, V, dvdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      potential = potential + V

   end do h3s
!$OMP end parallel do



!for the H4 surface
!$OMP PARALLEL DO private(gi3,gi4,gi5,gi6,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!$OMP REDUCTION(+:potential)
   h4: do i = 1, npart

      if (mask1(i, 9)) cycle h4

!      if((tempo.eq.53).and.(tempo.le.54)) then
!      print*,'debuuuug'
!      print*,i,gi3,gi4,gi5,gi6
!      call pgroup(1,i)
!      print*,mask1(i,9)
!      end if

      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! H atom
      gi5 = group(i, 5)      ! H atom
      gi6 = group(i, 6)      ! H atom

      call h4bmkp_switch(i, V, dvdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      dpotx(gi6) = dpotx(gi6) + dvdabccc(10)
      dpoty(gi6) = dpoty(gi6) + dvdabccc(11)
      dpotz(gi6) = dpotz(gi6) + dvdabccc(12)

      potential = potential + V

   end do h4
!$OMP end parallel do

!for the H2O2s surface
!!$OMP PARALLEL DO private(gi3,gi4,gi5,gi6,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!!$OMP REDUCTION(+:potential)
   h2o2s: do i = 1, npart
      if (mask1(i, 10)) cycle h2o2s

      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! H atom
      gi5 = group(i, 5)      ! O atom
      gi6 = group(i, 6)      ! O atom

      call vh2o2_switch(i, V, dvdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      dpotx(gi6) = dpotx(gi6) + dvdabccc(10)
      dpoty(gi6) = dpoty(gi6) + dvdabccc(11)
      dpotz(gi6) = dpotz(gi6) + dvdabccc(12)

      potential = potential + V

   end do h2o2s
!!$OMP end parallel do

!for the O4 surface
!$OMP PARALLEL DO private(gi3,gi4,gi5,gi6,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!$OMP REDUCTION(+:potential)
   o4: do i = 1, npart
      if (mask1(i, 11)) cycle o4

      gi3 = group(i, 3)      ! O atom
      gi4 = group(i, 4)      ! O atom
      gi5 = group(i, 5)      ! O atom
      gi6 = group(i, 6)      ! O atom

      call vo4_f11_switch(i, V, dVdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      dpotx(gi6) = dpotx(gi6) + dvdabccc(10)
      dpoty(gi6) = dpoty(gi6) + dvdabccc(11)
      dpotz(gi6) = dpotz(gi6) + dvdabccc(12)

      potential = potential + V

   end do o4
!$OMP end parallel do

!for the H3O surface
!$OMP PARALLEL DO private(gi3,gi4,gi5,gi6,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!$OMP REDUCTION(+:potential)
   h3o: do i = 1, npart
      if (mask1(i, 12)) cycle h3o

      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! H atom
      gi5 = group(i, 5)      ! H atom
      gi6 = group(i, 6)      ! O atom

      call vh3o_switch(i, V, dvdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      dpotx(gi6) = dpotx(gi6) + dvdabccc(10)
      dpoty(gi6) = dpoty(gi6) + dvdabccc(11)
      dpotz(gi6) = dpotz(gi6) + dvdabccc(12)

      potential = potential + V

   end do h3o
!$OMP end parallel do

!for the HO3 surface
!$OMP PARALLEL DO private(gi3,gi4,gi5,gi6,dvdabccc,V) SHARED(dpotx,dpoty,dpotz) &
!$OMP REDUCTION(+:potential)
   ho3: do i = 1, npart
      if (mask1(i, 13)) cycle ho3

      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! O atom
      gi5 = group(i, 5)      ! O atom
      gi6 = group(i, 6)      ! O atom

      call VHO3_13_switch(i, V, dVdabccc)

      dpotx(gi3) = dpotx(gi3) + dvdabccc(1)
      dpoty(gi3) = dpoty(gi3) + dvdabccc(2)
      dpotz(gi3) = dpotz(gi3) + dvdabccc(3)

      dpotx(gi4) = dpotx(gi4) + dvdabccc(4)
      dpoty(gi4) = dpoty(gi4) + dvdabccc(5)
      dpotz(gi4) = dpotz(gi4) + dvdabccc(6)

      dpotx(gi5) = dpotx(gi5) + dvdabccc(7)
      dpoty(gi5) = dpoty(gi5) + dvdabccc(8)
      dpotz(gi5) = dpotz(gi5) + dvdabccc(9)

      dpotx(gi6) = dpotx(gi6) + dvdabccc(10)
      dpoty(gi6) = dpoty(gi6) + dvdabccc(11)
      dpotz(gi6) = dpotz(gi6) + dvdabccc(12)

      potential = potential + V

   end do ho3
!$OMP end parallel do

!for the H2O2t surface
   h2o2t: do i = 1, npart
      if (mask1(i, 16)) cycle h2o2t

      gi3 = group(i, 3)      ! H atom
      gi4 = group(i, 4)      ! H atom
      gi5 = group(i, 5)      ! O atom
      gi6 = group(i, 6)      ! O atom

      call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), r2xij, r2yij, r2zij, r2)
      call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi5), ry(gi5), rz(gi5), r3xij, r3yij, r3zij, r3)
      call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi6), ry(gi6), rz(gi6), r5xij, r5yij, r5zij, r5)
      call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), r6xij, r6yij, r6zij, r6)
      call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi6), ry(gi6), rz(gi6), r4xij, r4yij, r4zij, r4)
      call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi6), ry(gi6), rz(gi6), r1xij, r1yij, r1zij, r1)

      r_4atoms(1) = r1
      r_4atoms(2) = r2
      r_4atoms(3) = r3
      r_4atoms(4) = r4
      r_4atoms(5) = r5
      r_4atoms(6) = r6

      call h2o2surt_f31(r_4atoms, V, dvdr_4atoms)

      potential = potential + V

      dvdr1 = dvdr_4atoms(1)
      dvdr2 = dvdr_4atoms(2)
      dvdr3 = dvdr_4atoms(3)
      dvdr4 = dvdr_4atoms(4)
      dvdr5 = dvdr_4atoms(5)
      dvdr6 = dvdr_4atoms(6)

      dpotx(gi3) = dpotx(gi3) + dvdr2*r2xij/r2 + dvdr3*r3xij/r3 + dvdr5*r5xij/r5
      dpoty(gi3) = dpoty(gi3) + dvdr2*r2yij/r2 + dvdr3*r3yij/r3 + dvdr5*r5yij/r5
      dpotz(gi3) = dpotz(gi3) + dvdr2*r2zij/r2 + dvdr3*r3zij/r3 + dvdr5*r5zij/r5

      dpotx(gi4) = dpotx(gi4) - dvdr2*r2xij/r2 + dvdr4*r4xij/r4 + dvdr6*r6xij/r6
      dpoty(gi4) = dpoty(gi4) - dvdr2*r2yij/r2 + dvdr4*r4yij/r4 + dvdr6*r6yij/r6
      dpotz(gi4) = dpotz(gi4) - dvdr2*r2zij/r2 + dvdr4*r4zij/r4 + dvdr6*r6zij/r6

      dpotx(gi5) = dpotx(gi5) + dvdr1*r1xij/r1 - dvdr3*r3xij/r3 - dvdr6*r6xij/r6
      dpoty(gi5) = dpoty(gi5) + dvdr1*r1yij/r1 - dvdr3*r3yij/r3 - dvdr6*r6yij/r6
      dpotz(gi5) = dpotz(gi5) + dvdr1*r1zij/r1 - dvdr3*r3zij/r3 - dvdr6*r6zij/r6

      dpotx(gi6) = dpotx(gi6) - dvdr1*r1xij/r1 - dvdr4*r4xij/r4 - dvdr5*r5xij/r5
      dpoty(gi6) = dpoty(gi6) - dvdr1*r1yij/r1 - dvdr4*r4yij/r4 - dvdr5*r5yij/r5
      dpotz(gi6) = dpotz(gi6) - dvdr1*r1zij/r1 - dvdr4*r4zij/r4 - dvdr5*r5zij/r5

   end do h2o2t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!finding how many atoms that must be replicated on opposite side of the box


   refcount(1, :) = merge(1, 0, (rx .lt. rl_cut_off)) + merge(-1, 0, (rx .gt. (boxl - rl_cut_off)))
   refcount(2, :) = merge(1, 0, (ry .lt. rl_cut_off)) + merge(-1, 0, (ry .gt. (boxl - rl_cut_off)))
   refcount(3, :) = merge(1, 0, (rz .lt. rl_cut_off)) + merge(-1, 0, (rz .gt. (boxl - rl_cut_off)))

   mask2 = .false.
   mask2 = refcount .ne. 0

   npart_full = sum(2**(count(mask2, 1)))



!allocating memory for original and transposed atoms
   
   if (allocated(fullxyz)) then 
      deallocate(fullxyz)
   end if

   allocate (fullxyz(3, npart_full))

   if (allocated(ref_fullxyz))  then 
      deallocate(ref_fullxyz)
   end if

   allocate (ref_fullxyz(npart_full))

!This array will be used for saving information about beeing original or reflected atoms
!All the original atoms will keep their number
!Beyond npart there are reflections which will be set below

!$OMP PARALLEL DO
   do i = 1, npart
      ref_fullxyz(i) = i
   end do
!$OMP END PARALLEL DO

   fullxyz = 0.0d0
   fullxyz(1, 1:npart) = rx(1:npart)
   fullxyz(2, 1:npart) = ry(1:npart)
   fullxyz(3, 1:npart) = rz(1:npart)

!cycles for finding reflected atoms
   p = npart

   pontos: do ip = 1, npart

      x: do i = 0, 1
         if ((i .eq. 1) .and. (refcount(1, ip) .eq. 0)) cycle x
         y: do j = 0, 1
            if ((j .eq. 1) .and. (refcount(2, ip) .eq. 0)) cycle y
            z: do k = 0, 1
               if ((k .eq. 1) .and. (refcount(3, ip) .eq. 0)) cycle z

               if ((i + j + k) .eq. 0) cycle z
               p = p + 1
               fullxyz(1, p) = fullxyz(1, ip) + i*refcount(1, ip)*boxl
               fullxyz(2, p) = fullxyz(2, ip) + j*refcount(2, ip)*boxl
               fullxyz(3, p) = fullxyz(3, ip) + k*refcount(3, ip)*boxl
               ref_fullxyz(p) = ip

            end do z
         end do y
      end do x

   end do pontos

! Creating tree with atoms and refleted atoms

!   print*,'antes',associated(tree) 

   if (associated(tree).eqv. .true.) call kdtree2_destroy(tree)
!   print*,'antes',associated(tree) 
   tree => kdtree2_create(fullxyz, sort=.true., rearrange=.true.)
!   print*,'depois',associated(tree) 

   deallocate (fullxyz)

   soma = 0

   ihh = 0
   ioo = 0
   iho = 0

   r2 = rl_cut_off
   r2 = r2*r2

!loop for checking neighbours of each atom
   kdnear: do i = 1, npart

      results(:)%idx = 0
      results(:)%dis = 0.0d0

      xyz(1, 1) = rx(i)
      xyz(2, 1) = ry(i)
      xyz(3, 1) = rz(i)

      gij = i

      call kdtree2_r_nearest(tp=tree, qv=xyz(:, 1), r2=r2, nfound=nfound, nalloc=nres, results=results)

      soma = soma + nfound
      cnfound: do k = 2, nfound

         gkl = results(k)%idx
         ! If gkl is a reflected atom then it should be the original number
         gkl = ref_fullxyz(gkl)

         !If they are part of the same system then jumpo to next atom
         if (rl_part(gij) .eq. rl_part(gkl)) cycle

         if (cname(gij) == 'Hy' .and. cname(gkl) == 'Hy') then
            ihh = ihh + 1
            hh_pair(ihh, 1) = gij
            hh_pair(ihh, 2) = gkl
         end if
         if (cname(gij) == 'Ox' .and. cname(gkl) == 'Ox') then

            ioo = ioo + 1
            oo_pair(ioo, 1) = gij
            oo_pair(ioo, 2) = gkl

         end if
         if ((cname(gij) == 'Hy' .and. cname(gkl) == 'Ox') .or. &
             (cname(gij) == 'Ox' .and. cname(gkl) == 'Hy')) then

            iho = iho + 1
            ho_pair(iho, 1) = gij
            ho_pair(iho, 2) = gkl

         end if

      end do cnfound

   end do kdnear

   soma = soma - npart
!deallocate(ref_fullxyz)

   k = 0

!nthre=2
!call OMP_SET_NUM_THREADS(nthre)

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE(gij,gkl,rxij,ryij,rzij,rij,V,dV) REDUCTION(+:POTENTIALD)

   hhcalc: do i = 1, ihh

      if (hh_pair(i, 1) .eq. 0) cycle hhcalc

      gij = hh_pair(i, 1)
      gkl = hh_pair(i, 2)

      call mod_dist(rx(gij), ry(gij), rz(gij), rx(gkl), ry(gkl), rz(gkl), rxij, ryij, rzij, rij)

      call pothht_f14(rij, V, dV)

      dpotxd(gij) = dpotxd(gij) + (dV*rxij/rij)
      dpotyd(gij) = dpotyd(gij) + (dV*ryij/rij)
      dpotzd(gij) = dpotzd(gij) + (dV*rzij/rij)

      potentiald = potentiald + V/2.0d0

   end do hhcalc
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE(gij,gkl,rxij,ryij,rzij,rij,V,dV) REDUCTION(+:POTENTIALD)

   oocalc: do i = 1, ioo

      if (oo_pair(i, 1) .eq. 0) cycle oocalc

      gij = oo_pair(i, 1)
      gkl = oo_pair(i, 2)

      call mod_dist(rx(gij), ry(gij), rz(gij), rx(gkl), ry(gkl), rz(gkl), rxij, ryij, rzij, rij)

      call poto2q_f18(rij, V, dV)

      dpotxd(gij) = dpotxd(gij) + dV*rxij/rij
      dpotyd(gij) = dpotyd(gij) + dV*ryij/rij
      dpotzd(gij) = dpotzd(gij) + dV*rzij/rij

      potentiald = potentiald + V/2.0d0

   end do oocalc

!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED)  PRIVATE(gij,gkl,rxij,ryij,rzij,rij,V,dV) REDUCTION(+:POTENTIALD)
   hocalc: do i = 1, iho
      if (ho_pair(i, 1) .eq. 0) cycle

      gij = ho_pair(i, 1)
      gkl = ho_pair(i, 2)

      call mod_dist(rx(gij), ry(gij), rz(gij), rx(gkl), ry(gkl), rz(gkl), rxij, ryij, rzij, rij)

      call pothoq_f17(rij, V, dV)

      dpotxd(gij) = dpotxd(gij) + dV*rxij/rij
      dpotyd(gij) = dpotyd(gij) + dV*ryij/rij
      dpotzd(gij) = dpotzd(gij) + dV*rzij/rij

      potentiald = potentiald + V/2.0d0

   end do hocalc
!$OMP END PARALLEL DO

   potential = potential + potentiald

   dpotx = dpotx + dpotxd
   dpoty = dpoty + dpotyd
   dpotz = dpotz + dpotzd

   if (tempo.lt.1) then
!      call kdtree2_destroy(tree)
      deallocate (ref_fullxyz)
   end if

end subroutine potential_total

subroutine mod_dist_s(rx1, ry1, rz1, rx2, ry2, rz2, drx, dry, drz, mod_rdist)
   use variables

   implicit none
   double precision, intent(in) :: rx1, ry1, rz1, rx2, ry2, rz2
   double precision, intent(out):: drx, dry, drz
   double precision, intent(out):: mod_rdist

   drx = rx1 - rx2
   dry = ry1 - ry2
   drz = rz1 - rz2

!!for setting the minimum image convention
!drx = drx -anint(drx/boxl)*boxl
!dry = dry -anint(dry/boxl)*boxl
!drz = drz -anint(drz/boxl)*boxl

!for setting the cutoff radius
!Determine the vector between the ij particles
   mod_rdist = sqrt(drx**2 + dry**2 + drz**2)
!print*,''
!print*,rx2,ry2,rz2
!print*,drx,dry,drz

!print*,1,mod_rdist

end subroutine mod_dist_s

subroutine mod_dist1(rx1, ry1, rz1, rx2, ry2, rz2, drx, dry, drz, mod_rdist)
   use variables
   use phys_parameters
   implicit none
   double precision, intent(in) :: rx1, ry1, rz1, rx2, ry2, rz2
   double precision, intent(out):: drx, dry, drz
   double precision, intent(out):: mod_rdist

   drx = rx1 - rx2
   dry = ry1 - ry2
   drz = rz1 - rz2

!for setting the minimum image convention
   drx = drx - anint(drx/boxl)*boxl
   dry = dry - anint(dry/boxl)*boxl
   drz = drz - anint(drz/boxl)*boxl

!for setting the cutoff radius
!Determine the vector between the ij particles
   mod_rdist = sqrt(drx**2 + dry**2 + drz**2)
!print*,''
!print*,rx2,ry2,rz2
!print*,drx,dry,drz

!print*,1,mod_rdist

end subroutine mod_dist1

subroutine mod_dist2(rx1, ry1, rz1, rx2, ry2, rz2, drx, dry, drz, mod_rdist)
   use variables
   use phys_parameters

   implicit none
   double precision, intent(in) :: rx1, ry1, rz1, rx2, ry2, rz2
   double precision, intent(out):: drx, dry, drz
   double precision, intent(out):: mod_rdist

   drx = rx1 - rx2
   dry = ry1 - ry2
   drz = rz1 - rz2

!for setting the minimum image convention
   drx = drx - anint(drx/boxl)*boxl
   dry = dry - anint(dry/boxl)*boxl
   drz = drz - anint(drz/boxl)*boxl

!for setting the cutoff radius
!Determine the vector between the ij particles
   mod_rdist = sqrt(drx**2 + dry**2 + drz**2)
!print*,''
!print*,rx2,ry2,rz2
!print*,drx,dry,drz

!print*,2,mod_rdist

end subroutine mod_dist2

subroutine mod_dist3(rx1, ry1, rz1, rx2, ry2, rz2, drx, dry, drz, mod_rdist)
   use variables
   use phys_parameters

   implicit none
   double precision, intent(in) :: rx1, ry1, rz1, rx2, ry2, rz2
   double precision, intent(out):: drx, dry, drz
   double precision, intent(out):: mod_rdist

   drx = rx1 - rx2
   dry = ry1 - ry2
   drz = rz1 - rz2

!for setting the minimum image convention
   drx = drx - anint(drx/boxl)*boxl
   dry = dry - anint(dry/boxl)*boxl
   drz = drz - anint(drz/boxl)*boxl

!for setting the cutoff radius
!Determine the vector between the ij particles
   mod_rdist = sqrt(drx**2 + dry**2 + drz**2)
!print*,''
!print*,rx2,ry2,rz2
!print*,drx,dry,drz

!print*,3,mod_rdist

end subroutine mod_dist3
subroutine mod_dist4(rx1, ry1, rz1, rx2, ry2, rz2, drx, dry, drz, mod_rdist)
   use variables
   use phys_parameters

   implicit none
   double precision, intent(in) :: rx1, ry1, rz1, rx2, ry2, rz2
   double precision, intent(out):: drx, dry, drz
   double precision, intent(out):: mod_rdist

   drx = rx1 - rx2
   dry = ry1 - ry2
   drz = rz1 - rz2

!for setting the minimum image convention
   drx = drx - anint(drx/boxl)*boxl
   dry = dry - anint(dry/boxl)*boxl
   drz = drz - anint(drz/boxl)*boxl

!for setting the cutoff radius
!Determine the vector between the ij particles
   mod_rdist = sqrt(drx**2 + dry**2 + drz**2)
!print*,''
!print*,rx2,ry2,rz2
!print*,drx,dry,drz

!print*,4,mod_rdist

end subroutine mod_dist4

subroutine mod_dist(rx1, ry1, rz1, rx2, ry2, rz2, drx, dry, drz, mod_rdist)
   use variables
   use phys_parameters

   implicit none
   double precision, intent(in) :: rx1, ry1, rz1, rx2, ry2, rz2
   double precision, intent(out):: drx, dry, drz
   double precision, intent(out):: mod_rdist

   drx = rx1 - rx2
   dry = ry1 - ry2
   drz = rz1 - rz2

!for setting the minimum image convention
   drx = drx - anint(drx/boxl)*boxl
   dry = dry - anint(dry/boxl)*boxl
   drz = drz - anint(drz/boxl)*boxl

!for setting the cutoff radius
!Determine the vector between the ij particles
   mod_rdist = sqrt(drx**2 + dry**2 + drz**2)
!print*,''
!print*,rx2,ry2,rz2
!print*,drx,dry,drz

!print*,mod_rdist

end subroutine mod_dist

