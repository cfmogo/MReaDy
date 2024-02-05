subroutine h4bmkp_switch(i, V, dv)
   use phys_parameters
   use variables
   use constants

   double precision, dimension(12) :: dV1, dV
   integer      :: gi3, gi4, gi5, gi6
   integer      :: i, j, k, l
! double precision, dimension(3) :: acc,bcc,ccc,dcc
   double precision :: dist
   double precision :: rxij, ryij, rzij, rij
   double precision :: r1xij, r1yij, r1zij
   double precision :: r2xij, r2yij, r2zij
   double precision :: r3xij, r3yij, r3zij
   double precision :: r4xij, r4yij, r4zij
   double precision :: r5xij, r5yij, r5zij
   double precision :: r6xij, r6yij, r6zij
   double precision :: dv2dr1, dv2dr2, dv2dr3, dv2dr4, dv2dr5, dv2dr6
   double precision :: r1, r2, r3, r4, r5, r6

   double precision, dimension(6) :: rin!,fn
   double precision :: fs!gama,rzero,fs
   double precision, dimension(3) :: rH3, dvdrH3!dvdrV2
   double precision :: V2, V1, V, dV2, dV3, dvdiat, V4, vdiat, v3
   double precision :: f3, f4
!double precision :: dr1,dr2,dr3,dr4
   double precision :: dfsdr1, dfsdr2, dfsdr3, dfsdr4, dfsdr5, dfsdr6
   double precision, dimension(4, 3) :: vect_4atoms, dvdcc
!double precision, dimension(6) ::r_4atoms,dvdr_4atoms
!double precision :: Derr

   fs = 0.0d0
   V = 0.0d0
   V1 = 0.0d0
   V2 = 0.0d0
   V3 = 0.0d0
   dv1 = 0.0d0
   dv2 = 0.0d0
   dv3 = 0.0d0
   dv4 = 0.0d0
   dv = 0.0d0
   V4 = 0.0d0
   vect_4atoms = 0.0d0
   dvdcc = 0.0d0

   dfsdr1 = 0.0d0
   dfsdr2 = 0.0d0
   dfsdr3 = 0.0d0
   dfsdr4 = 0.0d0
   dfsdr5 = 0.0d0
   dfsdr6 = 0.0d0

   r1xij = 0.0d0
   r1yij = 0.0d0
   r1zij = 0.0d0

   r2xij = 0.0d0
   r2yij = 0.0d0
   r2zij = 0.0d0

   r3xij = 0.0d0
   r3yij = 0.0d0
   r3zij = 0.0d0

   r4xij = 0.0d0
   r4yij = 0.0d0
   r4zij = 0.0d0

   r5xij = 0.0d0
   r5yij = 0.0d0
   r5zij = 0.0d0

   r6xij = 0.0d0
   r6yij = 0.0d0
   r6zij = 0.0d0

   dv2dr1 = 0.0d0
   dv2dr2 = 0.0d0
   dv2dr3 = 0.0d0
   dv2dr4 = 0.0d0
   dv2dr5 = 0.0d0
   dv2dr6 = 0.0d0

!Geting the coordenates
!print*,'gi3,gi4,gi5,gi6'
!print*,gi3,gi4,gi5,gi6
!print*,'coordinout'
!print*,rx(gi3),ry(gi3),rz(gi3)
!print*,rx(gi4),ry(gi4),rz(gi4)
!print*,rx(gi5),ry(gi5),rz(gi5)
!print*,rx(gi6),ry(gi6),rz(gi6)



!The positions are predefined in the matrix "group"
   gi3 = group(i, 3)      ! H atom
   gi4 = group(i, 4)      ! H atom
   gi5 = group(i, 5)      ! H atom
   gi6 = group(i, 6)      ! H atom

   vect_4atoms(1, 1) = rx(gi3)
   vect_4atoms(1, 2) = ry(gi3)
   vect_4atoms(1, 3) = rz(gi3)

   vect_4atoms(2, 1) = rx(gi4)
   vect_4atoms(2, 2) = ry(gi4)
   vect_4atoms(2, 3) = rz(gi4)

   vect_4atoms(3, 1) = rx(gi5)
   vect_4atoms(3, 2) = ry(gi5)
   vect_4atoms(3, 3) = rz(gi5)

   vect_4atoms(4, 1) = rx(gi6)
   vect_4atoms(4, 2) = ry(gi6)
   vect_4atoms(4, 3) = rz(gi6)

   do k = 1, 3
   do l = 1, 3
   do j = l + 1, 4

      dist = vect_4atoms(l, k) - vect_4atoms(j, k)
      if (dist .gt. boxl/2.0d0) then
         vect_4atoms(j, k) = vect_4atoms(j, k) + boxl
      elseif (dist .lt. -(boxl/2.0d0)) then
         vect_4atoms(l, k) = vect_4atoms(l, k) + boxl
      end if
   end do
   end do
   end do

!       The quadriatomic surface

   call h4bmkp_cc(vect_4atoms, V1, dVdcc, 1)

   dv1(1) = dv1(1) + dvdcc(1, 1)
   dv1(2) = dv1(2) + dvdcc(1, 2)
   dv1(3) = dv1(3) + dvdcc(1, 3)

   dv1(4) = dv1(4) + dvdcc(2, 1)
   dv1(5) = dv1(5) + dvdcc(2, 2)
   dv1(6) = dv1(6) + dvdcc(2, 3)

   dv1(7) = dv1(7) + dvdcc(3, 1)
   dv1(8) = dv1(8) + dvdcc(3, 2)
   dv1(9) = dv1(9) + dvdcc(3, 3)

   dv1(10) = dv1(10) + dvdcc(4, 1)
   dv1(11) = dv1(11) + dvdcc(4, 2)
   dv1(12) = dv1(12) + dvdcc(4, 3)

! Checking inner distance

   call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), &
                 rz(gi4), rxij, ryij, rzij, rij)

   rin(2) = rij
   r2 = rij
   r2xij = rxij
   r2yij = ryij
   r2zij = rzij

   call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi5), ry(gi5), &
                 rz(gi5), rxij, ryij, rzij, rij)

   rin(3) = rij
   r3 = rij
   r3xij = rxij
   r3yij = ryij
   r3zij = rzij

   call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi6), ry(gi6), &
                 rz(gi6), rxij, ryij, rzij, rij)

   rin(5) = rij
   r5 = rij
   r5xij = rxij
   r5yij = ryij
   r5zij = rzij

   call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), &
                 rz(gi5), rxij, ryij, rzij, rij)

   rin(6) = rij
   r6 = rij
   r6xij = rxij
   r6yij = ryij
   r6zij = rzij

   call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi6), ry(gi6), &
                 rz(gi6), rxij, ryij, rzij, rij)

   rin(4) = rij
   r4 = rij
   r4xij = rxij
   r4yij = ryij
   r4zij = rzij

   call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi6), ry(gi6), &
                 rz(gi6), rxij, ryij, rzij, rij)

   rin(1) = rij
   r1 = rij
   r1xij = rxij
   r1yij = ryij
   r1zij = rzij

!'H4 -> H(4) + H3'

   if ((r1 .ge. in_bound) .and. (r4 .ge. in_bound) .and. &
       (r5 .ge. in_bound)) then
      if (group(i, 9) .eq. 0) group(i, 9) = 11
      if (group(i, 9) .ne. 11) then
         print *, 'exchange in molecular complex 1 ', tempo
         group(i, 9) = 11
      end if

!            print*,'tempo',tempo

! H3 potential
      rH3(1) = rin(2)
      rH3(2) = rin(3)
      rH3(3) = rin(6)

      call bkmp2_f8(rH3, V2, dvdrH3, 1)

      dv2dr2 = dvdrH3(1)
      dv2dr3 = dvdrH3(2)
      dv2dr6 = dvdrH3(3)

! Switch function

      fs = f3(rin(1), rin(4), rin(5))
      call df3(rin(1), rin(4), rin(5), dfsdr1, dfsdr4, dfsdr5)

! for the diatomic potencials

      call pothht_f14(rin(5), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r5xij/r5
      dv(2) = dv(2) + fs*dVdiat*r5yij/r5
      dv(3) = dv(3) + fs*dVdiat*r5zij/r5

      dv(10) = dv(10) - fs*dVdiat*r5xij/r5
      dv(11) = dv(11) - fs*dVdiat*r5yij/r5
      dv(12) = dv(12) - fs*dVdiat*r5zij/r5

      V4 = V4 + Vdiat

      call pothht_f14(rin(4), Vdiat, dVdiat)

      dv(4) = dv(4) + fs*dVdiat*r4xij/r4
      dv(5) = dv(5) + fs*dVdiat*r4yij/r4
      dv(6) = dv(6) + fs*dVdiat*r4zij/r4

      dv(10) = dv(10) - fs*dVdiat*r4xij/r4
      dv(11) = dv(11) - fs*dVdiat*r4yij/r4
      dv(12) = dv(12) - fs*dVdiat*r4zij/r4

      V4 = V4 + Vdiat

      call pothht_f14(rin(1), Vdiat, dVdiat)

      dv(7) = dv(7) + fs*dVdiat*r1xij/r1
      dv(8) = dv(8) + fs*dVdiat*r1yij/r1
      dv(9) = dv(9) + fs*dVdiat*r1zij/r1

      dv(10) = dv(10) - fs*dVdiat*r1xij/r1
      dv(11) = dv(11) - fs*dVdiat*r1yij/r1
      dv(12) = dv(12) - fs*dVdiat*r1zij/r1

      V4 = V4 + Vdiat

      V = V1*(1 - fs) + (V2 + V4)*fs

      dv = dv + dv1*(1 - fs)

      dv(1) = dv(1) + (v4 + v2 - v1)*(dfsdr5*r5xij/r5)
      dv(2) = dv(2) + (v4 + v2 - v1)*(dfsdr5*r5yij/r5)
      dv(3) = dv(3) + (v4 + v2 - v1)*(dfsdr5*r5zij/r5)

      dv(4) = dv(4) + (v4 + v2 - v1)*(dfsdr4*r4xij/r4)
      dv(5) = dv(5) + (v4 + v2 - v1)*(dfsdr4*r4yij/r4)
      dv(6) = dv(6) + (v4 + v2 - v1)*(dfsdr4*r4zij/r4)

      dv(7) = dv(7) + (v4 + v2 - v1)*(dfsdr1*r1xij/r1)
      dv(8) = dv(8) + (v4 + v2 - v1)*(dfsdr1*r1yij/r1)
      dv(9) = dv(9) + (v4 + v2 - v1)*(dfsdr1*r1zij/r1)

      dv(10) = dv(10) + (v4 + v2 - v1)*(-dfsdr1*r1xij/r1 - dfsdr4*r4xij/r4 - dfsdr5*r5xij/r5)
      dv(11) = dv(11) + (v4 + v2 - v1)*(-dfsdr1*r1yij/r1 - dfsdr4*r4yij/r4 - dfsdr5*r5yij/r5)
      dv(12) = dv(12) + (v4 + v2 - v1)*(-dfsdr1*r1zij/r1 - dfsdr4*r4zij/r4 - dfsdr5*r5zij/r5)

      dv(1) = dv(1) + fs*(dv2dr2*r2xij/r2 + dv2dr3*r3xij/r3)
      dv(2) = dv(2) + fs*(dv2dr2*r2yij/r2 + dv2dr3*r3yij/r3)
      dv(3) = dv(3) + fs*(dv2dr2*r2zij/r2 + dv2dr3*r3zij/r3)

      dv(4) = dv(4) + fs*(-dv2dr2*r2xij/r2 + dv2dr6*r6xij/r6)
      dv(5) = dv(5) + fs*(-dv2dr2*r2yij/r2 + dv2dr6*r6yij/r6)
      dv(6) = dv(6) + fs*(-dv2dr2*r2zij/r2 + dv2dr6*r6zij/r6)

      dv(7) = dv(7) + fs*(-dv2dr3*r3xij/r3 - dv2dr6*r6xij/r6)
      dv(8) = dv(8) + fs*(-dv2dr3*r3yij/r3 - dv2dr6*r6yij/r6)
      dv(9) = dv(9) + fs*(-dv2dr3*r3zij/r3 - dv2dr6*r6zij/r6)

!       print*,'dv',dv

!H4 -> H(1) + H3
   else if ((r2 .ge. in_bound) .and. (r3 .ge. in_bound) .and. &
            (r5 .ge. in_bound)) then

      if (group(i, 9) .eq. 0) group(i, 9) = 12

      if (group(i, 9) .ne. 12) then
         print *, 'Change in group'
         group(i, 9) = 12
      end if

!        print*,'here we are'
!H4 -> H(1) + H3

!print*,"mid H(1) + H3",tempo
      print *, 'Exchange in channels. Time:', t

! H3 potential
      rH3(1) = rin(6)
      rH3(2) = rin(1)
      rH3(3) = rin(4)

      call bkmp2_f8(rH3, V2, dvdrH3, 1)

      dv2dr6 = dvdrH3(1)
      dv2dr1 = dvdrH3(2)
      dv2dr4 = dvdrH3(3)

! Switch function

      fs = f3(rin(2), rin(3), rin(5))
      call df3(rin(2), rin(3), rin(5), dfsdr2, dfsdr3, dfsdr5)

! for the diatomic potencials

      call pothht_f14(rin(5), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r5xij/r5
      dv(2) = dv(2) + fs*dVdiat*r5yij/r5
      dv(3) = dv(3) + fs*dVdiat*r5zij/r5

      dv(10) = dv(10) - fs*dVdiat*r5xij/r5
      dv(11) = dv(11) - fs*dVdiat*r5yij/r5
      dv(12) = dv(12) - fs*dVdiat*r5zij/r5

      V4 = V4 + Vdiat

      call pothht_f14(rin(2), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r2xij/r2
      dv(2) = dv(2) + fs*dVdiat*r2yij/r2
      dv(3) = dv(3) + fs*dVdiat*r2zij/r2

      dv(4) = dv(4) - fs*dVdiat*r2xij/r2
      dv(5) = dv(5) - fs*dVdiat*r2yij/r2
      dv(6) = dv(6) - fs*dVdiat*r2zij/r2

      V4 = V4 + Vdiat

      call pothht_f14(rin(3), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r3xij/r3
      dv(2) = dv(2) + fs*dVdiat*r3yij/r3
      dv(3) = dv(3) + fs*dVdiat*r3zij/r3

      dv(7) = dv(7) - fs*dVdiat*r3xij/r3
      dv(8) = dv(8) - fs*dVdiat*r3yij/r3
      dv(9) = dv(9) - fs*dVdiat*r3zij/r3

      V4 = V4 + Vdiat

      V = V1*(1 - fs) + (V2 + V4)*fs

      dv = dv + dv1*(1 - fs)

      dv(1) = dv(1) + (v4 + v2 - v1)*(dfsdr2*r2xij/r2 + dfsdr3*r3xij/r3 + dfsdr5*r5xij/r5)
      dv(2) = dv(2) + (v4 + v2 - v1)*(dfsdr2*r2yij/r2 + dfsdr3*r3yij/r3 + dfsdr5*r5yij/r5)
      dv(3) = dv(3) + (v4 + v2 - v1)*(dfsdr2*r2zij/r2 + dfsdr3*r3zij/r3 + dfsdr5*r5zij/r5)

      dv(4) = dv(4) + (v4 + v2 - v1)*(-dfsdr2*r2xij/r2)
      dv(5) = dv(5) + (v4 + v2 - v1)*(-dfsdr2*r2yij/r2)
      dv(6) = dv(6) + (v4 + v2 - v1)*(-dfsdr2*r2zij/r2)

      dv(7) = dv(7) + (v4 + v2 - v1)*(-dfsdr3*r3xij/r3)
      dv(8) = dv(8) + (v4 + v2 - v1)*(-dfsdr3*r3yij/r3)
      dv(9) = dv(9) + (v4 + v2 - v1)*(-dfsdr3*r3zij/r3)

      dv(10) = dv(10) + (v4 + v2 - v1)*(-dfsdr5*r5xij/r5)
      dv(11) = dv(11) + (v4 + v2 - v1)*(-dfsdr5*r5yij/r5)
      dv(12) = dv(12) + (v4 + v2 - v1)*(-dfsdr5*r5zij/r5)

      dv(4) = dv(4) + fs*(dv2dr4*r4xij/r4 + dv2dr6*r6xij/r6)
      dv(5) = dv(5) + fs*(dv2dr4*r4yij/r4 + dv2dr6*r6yij/r6)
      dv(6) = dv(6) + fs*(dv2dr4*r4zij/r4 + dv2dr6*r6zij/r6)

      dv(7) = dv(7) + fs*(dv2dr1*r1xij/r1 - dv2dr6*r6xij/r6)
      dv(8) = dv(8) + fs*(dv2dr1*r1yij/r1 - dv2dr6*r6yij/r6)
      dv(9) = dv(9) + fs*(dv2dr1*r1zij/r1 - dv2dr6*r6zij/r6)

      dv(10) = dv(10) + fs*(-dv2dr1*r1xij/r1 - dv2dr4*r4xij/r4)
      dv(11) = dv(11) + fs*(-dv2dr1*r1yij/r1 - dv2dr4*r4yij/r4)
      dv(12) = dv(12) + fs*(-dv2dr1*r1zij/r1 - dv2dr4*r4zij/r4)

!'H4 -> H(2) + H3
   else if ((r2 .ge. in_bound) .and. (r4 .ge. in_bound) .and. &
            (r6 .ge. in_bound)) then

      if (group(i, 9) .eq. 0) group(i, 9) = 13

      if (group(i, 9) .ne. 13) then
         print *, 'Change in group'
         group(i, 9) = 13
      end if

!        print*,'here we are'
!H4 -> H(2) + H3

      print *, "mid H(2) + H3", tempo
!print*,'Exchange in channels. Time:',t

! H3 potential
      rH3(1) = rin(3)
      rH3(2) = rin(5)
      rH3(3) = rin(1)

      call bkmp2_f8(rH3, V2, dvdrH3, 1)

      dv2dr3 = dvdrH3(1)
      dv2dr5 = dvdrH3(2)
      dv2dr1 = dvdrH3(3)

! Switch function

      fs = f3(rin(2), rin(4), rin(6))
      call df3(rin(2), rin(4), rin(6), dfsdr2, dfsdr4, dfsdr6)

! for the diatomic potencials

      call pothht_f14(rin(4), Vdiat, dVdiat)

      dv(4) = dv(4) + fs*dVdiat*r4xij/r4
      dv(5) = dv(5) + fs*dVdiat*r4yij/r4
      dv(6) = dv(6) + fs*dVdiat*r4zij/r4

      dv(10) = dv(10) - fs*dVdiat*r4xij/r4
      dv(11) = dv(11) - fs*dVdiat*r4yij/r4
      dv(12) = dv(12) - fs*dVdiat*r4zij/r4

      V4 = V4 + Vdiat

      call pothht_f14(rin(2), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r2xij/r2
      dv(2) = dv(2) + fs*dVdiat*r2yij/r2
      dv(3) = dv(3) + fs*dVdiat*r2zij/r2

      dv(4) = dv(4) - fs*dVdiat*r2xij/r2
      dv(5) = dv(5) - fs*dVdiat*r2yij/r2
      dv(6) = dv(6) - fs*dVdiat*r2zij/r2

      V4 = V4 + Vdiat

      call pothht_f14(rin(6), Vdiat, dVdiat)

      dv(4) = dv(4) + fs*dVdiat*r6xij/r6
      dv(5) = dv(5) + fs*dVdiat*r6yij/r6
      dv(6) = dv(6) + fs*dVdiat*r6zij/r6

      dv(7) = dv(7) - fs*dVdiat*r6xij/r6
      dv(8) = dv(8) - fs*dVdiat*r6yij/r6
      dv(9) = dv(9) - fs*dVdiat*r6zij/r6

      V4 = V4 + Vdiat

      V = V1*(1 - fs) + (V2 + V4)*fs

      dv = dv + dv1*(1 - fs)

      dv(1) = dv(1) + (v4 + v2 - v1)*(dfsdr2*r2xij/r2)
      dv(2) = dv(2) + (v4 + v2 - v1)*(dfsdr2*r2yij/r2)
      dv(3) = dv(3) + (v4 + v2 - v1)*(dfsdr2*r2zij/r2)

      dv(4) = dv(4) + (v4 + v2 - v1)*(-dfsdr2*r2xij/r2 + dfsdr4*r4xij/r4 + dfsdr6*r6xij/r6)
      dv(5) = dv(5) + (v4 + v2 - v1)*(-dfsdr2*r2yij/r2 + dfsdr4*r4yij/r4 + dfsdr6*r6yij/r6)
      dv(6) = dv(6) + (v4 + v2 - v1)*(-dfsdr2*r2zij/r2 + dfsdr4*r4zij/r4 + dfsdr6*r6zij/r6)

      dv(7) = dv(7) + (v4 + v2 - v1)*(-dfsdr6*r6xij/r6)
      dv(8) = dv(8) + (v4 + v2 - v1)*(-dfsdr6*r6yij/r6)
      dv(9) = dv(9) + (v4 + v2 - v1)*(-dfsdr6*r6zij/r6)

      dv(10) = dv(10) + (v4 + v2 - v1)*(-dfsdr4*r4xij/r4)
      dv(11) = dv(11) + (v4 + v2 - v1)*(-dfsdr4*r4yij/r4)
      dv(12) = dv(12) + (v4 + v2 - v1)*(-dfsdr4*r4zij/r4)

      dv(1) = dv(1) + fs*(dv2dr3*r3xij/r3 + dv2dr5*r5xij/r5)
      dv(2) = dv(2) + fs*(dv2dr3*r3yij/r3 + dv2dr5*r5yij/r5)
      dv(3) = dv(3) + fs*(dv2dr3*r3zij/r3 + dv2dr5*r5zij/r5)

      dv(7) = dv(7) + fs*(dv2dr1*r1xij/r1 - dv2dr3*r3xij/r3)
      dv(8) = dv(8) + fs*(dv2dr1*r1yij/r1 - dv2dr3*r3yij/r3)
      dv(9) = dv(9) + fs*(dv2dr1*r1zij/r1 - dv2dr3*r3zij/r3)

      dv(10) = dv(10) + fs*(-dv2dr1*r1xij/r1 - dv2dr5*r5xij/r5)
      dv(11) = dv(11) + fs*(-dv2dr1*r1yij/r1 - dv2dr5*r5yij/r5)
      dv(12) = dv(12) + fs*(-dv2dr1*r1zij/r1 - dv2dr5*r5zij/r5)

!'H -> H(3) + H3
   else if ((r1 .ge. in_bound) .and. (r3 .ge. in_bound) .and. &
            (r6 .ge. in_bound)) then

      if (group(i, 9) .eq. 0) group(i, 9) = 14

      if (group(i, 9) .ne. 14) then
         print *, 'Change in group'
         group(i, 9) = 14
      end if
!        print*,'here we are'
!H4 -> H(3) + H3

!             print*,"mid H(2) + H3",tempo
!print*,'Exchange in channels. Time:',t

! H3 potential
      rH3(1) = rin(2)
      rH3(2) = rin(4)
      rH3(3) = rin(5)

      call bkmp2_f8(rH3, V2, dvdrH3, 1)

      dv2dr2 = dvdrH3(1)
      dv2dr4 = dvdrH3(2)
      dv2dr5 = dvdrH3(3)

! Switch function

      fs = f3(rin(1), rin(3), rin(6))
      call df3(rin(1), rin(3), rin(6), dfsdr1, dfsdr3, dfsdr6)

! for the diatomic potencials

      call pothht_f14(rin(1), Vdiat, dVdiat)

      dv(7) = dv(7) + fs*dVdiat*r1xij/r1
      dv(8) = dv(8) + fs*dVdiat*r1yij/r1
      dv(9) = dv(9) + fs*dVdiat*r1zij/r1

      dv(10) = dv(10) - fs*dVdiat*r1xij/r1
      dv(11) = dv(11) - fs*dVdiat*r1yij/r1
      dv(12) = dv(12) - fs*dVdiat*r1zij/r1

      V4 = V4 + Vdiat

      call pothht_f14(rin(3), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r3xij/r3
      dv(2) = dv(2) + fs*dVdiat*r3yij/r3
      dv(3) = dv(3) + fs*dVdiat*r3zij/r3

      dv(7) = dv(7) - fs*dVdiat*r3xij/r3
      dv(8) = dv(8) - fs*dVdiat*r3yij/r3
      dv(9) = dv(9) - fs*dVdiat*r3zij/r3

      V4 = V4 + Vdiat

      call pothht_f14(rin(6), Vdiat, dVdiat)

      dv(4) = dv(4) + fs*dVdiat*r6xij/r6
      dv(5) = dv(5) + fs*dVdiat*r6yij/r6
      dv(6) = dv(6) + fs*dVdiat*r6zij/r6

      dv(7) = dv(7) - fs*dVdiat*r6xij/r6
      dv(8) = dv(8) - fs*dVdiat*r6yij/r6
      dv(9) = dv(9) - fs*dVdiat*r6zij/r6

      V4 = V4 + Vdiat

      V = V1*(1 - fs) + (V2 + V4)*fs

      dv = dv + dv1*(1 - fs)

      dv(1) = dv(1) + (v4 + v2 - v1)*(dfsdr3*r3xij/r3)
      dv(2) = dv(2) + (v4 + v2 - v1)*(dfsdr3*r3yij/r3)
      dv(3) = dv(3) + (v4 + v2 - v1)*(dfsdr3*r3zij/r3)

      dv(4) = dv(4) + (v4 + v2 - v1)*(dfsdr6*r6xij/r6)
      dv(5) = dv(5) + (v4 + v2 - v1)*(dfsdr6*r6yij/r6)
      dv(6) = dv(6) + (v4 + v2 - v1)*(dfsdr6*r6zij/r6)

      dv(7) = dv(7) + (v4 + v2 - v1)*(dfsdr1*r1xij/r1 - dfsdr3*r3xij/r3 - dfsdr6*r6xij/r6)
      dv(8) = dv(8) + (v4 + v2 - v1)*(dfsdr1*r1yij/r1 - dfsdr3*r3yij/r3 - dfsdr6*r6yij/r6)
      dv(9) = dv(9) + (v4 + v2 - v1)*(dfsdr1*r1zij/r1 - dfsdr3*r3zij/r3 - dfsdr6*r6zij/r6)

      dv(10) = dv(10) + (v4 + v2 - v1)*(-dfsdr1*r1xij/r1)
      dv(11) = dv(11) + (v4 + v2 - v1)*(-dfsdr1*r1yij/r1)
      dv(12) = dv(12) + (v4 + v2 - v1)*(-dfsdr1*r1zij/r1)

      dv(1) = dv(1) + fs*(dv2dr2*r2xij/r2 + dv2dr5*r5xij/r5)
      dv(2) = dv(2) + fs*(dv2dr2*r2yij/r2 + dv2dr5*r5yij/r5)
      dv(3) = dv(3) + fs*(dv2dr2*r2zij/r2 + dv2dr5*r5zij/r5)

      dv(4) = dv(4) + fs*(-dv2dr2*r2xij/r2 + dv2dr4*r4xij/r4)
      dv(5) = dv(5) + fs*(-dv2dr2*r2yij/r2 + dv2dr4*r4yij/r4)
      dv(6) = dv(6) + fs*(-dv2dr2*r2zij/r2 + dv2dr4*r4zij/r4)

      dv(10) = dv(10) + fs*(-dv2dr4*r4xij/r4 - dv2dr5*r5xij/r5)
      dv(11) = dv(11) + fs*(-dv2dr4*r4yij/r4 - dv2dr5*r5yij/r5)
      dv(12) = dv(12) + fs*(-dv2dr4*r4zij/r4 - dv2dr5*r5zij/r5)

   else if ((r1 .ge. in_bound) .and. (r2 .ge. in_bound) .and. &
            (r5 .ge. in_bound) .and. (r6 .ge. in_bound)) then
!              print*,'exit'

! 'H4 -> H2(1-3) + H2(2-4)'
      if (group(i, 9) .eq. 0) group(i, 9) = 15

      if (group(i, 9) .ne. 15) then
         print *, 'Change in group'
!        stop 'here'
         group(i, 9) = 15
      end if
! H2 function
      call pothhs_f1(r3, V2, dV2)
! H2 function
      call pothhs_f1(r4, V3, dV3)

! Switch function
      fs = f4(rin(1), rin(2), rin(5), rin(6))
      call df4(rin(1), rin(2), rin(5), rin(6), dfsdr1, dfsdr2, dfsdr5, dfsdr6)

! for the diatomic potencials

      call pothht_f14(rin(5), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r5xij/r5
      dv(2) = dv(2) + fs*dVdiat*r5yij/r5
      dv(3) = dv(3) + fs*dVdiat*r5zij/r5

      dv(10) = dv(10) - fs*dVdiat*r5xij/r5
      dv(11) = dv(11) - fs*dVdiat*r5yij/r5
      dv(12) = dv(12) - fs*dVdiat*r5zij/r5

      V4 = V4 + Vdiat

      call pothht_f14(rin(2), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r2xij/r2
      dv(2) = dv(2) + fs*dVdiat*r2yij/r2
      dv(3) = dv(3) + fs*dVdiat*r2zij/r2

      dv(4) = dv(4) - fs*dVdiat*r2xij/r2
      dv(5) = dv(5) - fs*dVdiat*r2yij/r2
      dv(6) = dv(6) - fs*dVdiat*r2zij/r2

      V4 = V4 + Vdiat

      call pothht_f14(rin(6), Vdiat, dVdiat)

      dv(4) = dv(4) + fs*dVdiat*r6xij/r6
      dv(5) = dv(5) + fs*dVdiat*r6yij/r6
      dv(6) = dv(6) + fs*dVdiat*r6zij/r6

      dv(7) = dv(7) - fs*dVdiat*r6xij/r6
      dv(8) = dv(8) - fs*dVdiat*r6yij/r6
      dv(9) = dv(9) - fs*dVdiat*r6zij/r6

      V4 = V4 + Vdiat

      call pothht_f14(rin(1), Vdiat, dVdiat)

      dv(7) = dv(7) + fs*dVdiat*r1xij/r1
      dv(8) = dv(8) + fs*dVdiat*r1yij/r1
      dv(9) = dv(9) + fs*dVdiat*r1zij/r1

      dv(10) = dv(10) - fs*dVdiat*r1xij/r1
      dv(11) = dv(11) - fs*dVdiat*r1yij/r1
      dv(12) = dv(12) - fs*dVdiat*r1zij/r1

      V4 = V4 + Vdiat

      V = V1*(1 - fs) + (V2 + V3 + V4)*fs

      dv = dv + dv1*(1 - fs)

      dv(1) = dv(1) + (v2 + v3 + v4 - v1)*(dfsdr2*r2xij/r2 + dfsdr5*r5xij/r5)
      dv(2) = dv(2) + (v2 + v3 + v4 - v1)*(dfsdr2*r2yij/r2 + dfsdr5*r5yij/r5)
      dv(3) = dv(3) + (v2 + v3 + v4 - v1)*(dfsdr2*r2zij/r2 + dfsdr5*r5zij/r5)

      dv(4) = dv(4) + (v2 + v3 + v4 - v1)*(-dfsdr2*r2xij/r2 + dfsdr6*r6xij/r6)
      dv(5) = dv(5) + (v2 + v3 + v4 - v1)*(-dfsdr2*r2yij/r2 + dfsdr6*r6yij/r6)
      dv(6) = dv(6) + (v2 + v3 + v4 - v1)*(-dfsdr2*r2zij/r2 + dfsdr6*r6zij/r6)

      dv(7) = dv(7) + (v2 + v3 + v4 - v1)*(dfsdr1*r1xij/r1 - dfsdr6*r6xij/r6)
      dv(8) = dv(8) + (v2 + v3 + v4 - v1)*(dfsdr1*r1yij/r1 - dfsdr6*r6yij/r6)
      dv(9) = dv(9) + (v2 + v3 + v4 - v1)*(dfsdr1*r1zij/r1 - dfsdr6*r6zij/r6)

      dv(10) = dv(10) + (v2 + v3 + v4 - v1)*(-dfsdr1*r1xij/r1 - dfsdr5*r5xij/r5)
      dv(11) = dv(11) + (v2 + v3 + v4 - v1)*(-dfsdr1*r1yij/r1 - dfsdr5*r5yij/r5)
      dv(12) = dv(12) + (v2 + v3 + v4 - v1)*(-dfsdr1*r1zij/r1 - dfsdr5*r5zij/r5)

      dv(1) = dv(1) + fs*(dv2*r3xij/r3)
      dv(2) = dv(2) + fs*(dv2*r3yij/r3)
      dv(3) = dv(3) + fs*(dv2*r3zij/r3)

      dv(4) = dv(4) + fs*(dv3*r4xij/r4)
      dv(5) = dv(5) + fs*(dv3*r4yij/r4)
      dv(6) = dv(6) + fs*(dv3*r4zij/r4)

      dv(7) = dv(7) + fs*(-dv2*r3xij/r3)
      dv(8) = dv(8) + fs*(-dv2*r3yij/r3)
      dv(9) = dv(9) + fs*(-dv2*r3zij/r3)

      dv(10) = dv(10) + fs*(-dv3*r4xij/r4)
      dv(11) = dv(11) + fs*(-dv3*r4yij/r4)
      dv(12) = dv(12) + fs*(-dv3*r4zij/r4)

   else if ((r1 .ge. in_bound) .and. (r2 .ge. in_bound) .and. &
            (r3 .ge. in_bound) .and. (r4 .ge. in_bound)) then
!              print*,'exit5'

! 'H4 -> H2(2-3) + H2(1-4)'
      if (group(i, 9) .eq. 0) group(i, 9) = 16

      if (group(i, 9) .ne. 16) then
         print *, 'Change in group'
         group(i, 9) = 16
      end if

! H2 function
      call pothhs_f1(r6, V2, dV2)
! OH function
      call pothhs_f1(r5, V3, dV3)

! Switch function
      fs = f4(rin(1), rin(2), rin(3), rin(4))
      call df4(rin(1), rin(2), rin(3), rin(4), dfsdr1, dfsdr2, dfsdr3, dfsdr4)

! for the diatomic potencials

      call pothht_f14(rin(2), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r2xij/r2
      dv(2) = dv(2) + fs*dVdiat*r2yij/r2
      dv(3) = dv(3) + fs*dVdiat*r2zij/r2

      dv(4) = dv(4) - fs*dVdiat*r2xij/r2
      dv(5) = dv(5) - fs*dVdiat*r2yij/r2
      dv(6) = dv(6) - fs*dVdiat*r2zij/r2

      V4 = V4 + Vdiat

      call pothht_f14(rin(3), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r3xij/r3
      dv(2) = dv(2) + fs*dVdiat*r3yij/r3
      dv(3) = dv(3) + fs*dVdiat*r3zij/r3

      dv(7) = dv(7) - fs*dVdiat*r3xij/r3
      dv(8) = dv(8) - fs*dVdiat*r3yij/r3
      dv(9) = dv(9) - fs*dVdiat*r3zij/r3

      V4 = V4 + Vdiat

      call pothht_f14(rin(4), Vdiat, dVdiat)

      dv(4) = dv(4) + fs*dVdiat*r4xij/r4
      dv(5) = dv(5) + fs*dVdiat*r4yij/r4
      dv(6) = dv(6) + fs*dVdiat*r4zij/r4

      dv(10) = dv(10) - fs*dVdiat*r4xij/r4
      dv(11) = dv(11) - fs*dVdiat*r4yij/r4
      dv(12) = dv(12) - fs*dVdiat*r4zij/r4

      V4 = V4 + Vdiat

      call pothht_f14(rin(1), Vdiat, dVdiat)

      dv(7) = dv(7) + fs*dVdiat*r1xij/r1
      dv(8) = dv(8) + fs*dVdiat*r1yij/r1
      dv(9) = dv(9) + fs*dVdiat*r1zij/r1

      dv(10) = dv(10) - fs*dVdiat*r1xij/r1
      dv(11) = dv(11) - fs*dVdiat*r1yij/r1
      dv(12) = dv(12) - fs*dVdiat*r1zij/r1

      V4 = V4 + Vdiat

      V = V1*(1 - fs) + (V2 + V3 + V4)*fs

      dv = dv + dv1*(1 - fs)

      dv(1) = dv(1) + (v2 + v3 + v4 - v1)*(dfsdr2*r2xij/r2 + dfsdr3*r3xij/r3)
      dv(2) = dv(2) + (v2 + v3 + v4 - v1)*(dfsdr2*r2yij/r2 + dfsdr3*r3yij/r3)
      dv(3) = dv(3) + (v2 + v3 + v4 - v1)*(dfsdr2*r2zij/r2 + dfsdr3*r3zij/r3)

      dv(4) = dv(4) + (v2 + v3 + v4 - v1)*(-dfsdr2*r2xij/r2 + dfsdr4*r4xij/r4)
      dv(5) = dv(5) + (v2 + v3 + v4 - v1)*(-dfsdr2*r2yij/r2 + dfsdr4*r4yij/r4)
      dv(6) = dv(6) + (v2 + v3 + v4 - v1)*(-dfsdr2*r2zij/r2 + dfsdr4*r4zij/r4)

      dv(7) = dv(7) + (v2 + v3 + v4 - v1)*(dfsdr1*r1xij/r1 - dfsdr3*r3xij/r3)
      dv(8) = dv(8) + (v2 + v3 + v4 - v1)*(dfsdr1*r1yij/r1 - dfsdr3*r3yij/r3)
      dv(9) = dv(9) + (v2 + v3 + v4 - v1)*(dfsdr1*r1zij/r1 - dfsdr3*r3zij/r3)

      dv(10) = dv(10) + (v2 + v3 + v4 - v1)*(-dfsdr1*r1xij/r1 - dfsdr4*r4xij/r4)
      dv(11) = dv(11) + (v2 + v3 + v4 - v1)*(-dfsdr1*r1yij/r1 - dfsdr4*r4yij/r4)
      dv(12) = dv(12) + (v2 + v3 + v4 - v1)*(-dfsdr1*r1zij/r1 - dfsdr4*r4zij/r4)

      dv(1) = dv(1) + fs*(dv3*r5xij/r5)
      dv(2) = dv(2) + fs*(dv3*r5yij/r5)
      dv(3) = dv(3) + fs*(dv3*r5zij/r5)

      dv(4) = dv(4) + fs*(dv2*r6xij/r6)
      dv(5) = dv(5) + fs*(dv2*r6yij/r6)
      dv(6) = dv(6) + fs*(dv2*r6zij/r6)

      dv(7) = dv(7) + fs*(-dv2*r6xij/r6)
      dv(8) = dv(8) + fs*(-dv2*r6yij/r6)
      dv(9) = dv(9) + fs*(-dv2*r6zij/r6)

      dv(10) = dv(10) + fs*(-dv3*r5xij/r5)
      dv(11) = dv(11) + fs*(-dv3*r5yij/r5)
      dv(12) = dv(12) + fs*(-dv3*r5zij/r5)

   else if ((r3 .ge. in_bound) .and. (r4 .ge. in_bound) .and. &
            (r5 .ge. in_bound) .and. (r6 .ge. in_bound)) then
!             print*,'exit4'

      if (group(i, 9) .eq. 0) group(i, 9) = 17

      if (group(i, 9) .ne. 17) then
         print *, 'Change in group'
         group(i, 9) = 17
      end if
! H2 function
      call pothhs_f1(r2, V2, dV2)
! H2 function
      call pothhs_f1(r1, V3, dV3)

! Switch function
      fs = f4(rin(3), rin(4), rin(5), rin(6))
      call df4(rin(3), rin(4), rin(5), rin(6), dfsdr3, dfsdr4, dfsdr5, dfsdr6)

! for the diatomic potencials

      call pothht_f14(rin(3), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r3xij/r3
      dv(2) = dv(2) + fs*dVdiat*r3yij/r3
      dv(3) = dv(3) + fs*dVdiat*r3zij/r3

      dv(7) = dv(7) - fs*dVdiat*r3xij/r3
      dv(8) = dv(8) - fs*dVdiat*r3yij/r3
      dv(9) = dv(9) - fs*dVdiat*r3zij/r3
      V4 = V4 + Vdiat

      call pothht_f14(rin(5), Vdiat, dVdiat)

      dv(1) = dv(1) + fs*dVdiat*r5xij/r5
      dv(2) = dv(2) + fs*dVdiat*r5yij/r5
      dv(3) = dv(3) + fs*dVdiat*r5zij/r5

      dv(10) = dv(10) - fs*dVdiat*r5xij/r5
      dv(11) = dv(11) - fs*dVdiat*r5yij/r5
      dv(12) = dv(12) - fs*dVdiat*r5zij/r5

      V4 = V4 + Vdiat

      call pothht_f14(rin(6), Vdiat, dVdiat)

      dv(4) = dv(4) + fs*dVdiat*r6xij/r6
      dv(5) = dv(5) + fs*dVdiat*r6yij/r6
      dv(6) = dv(6) + fs*dVdiat*r6zij/r6

      dv(7) = dv(7) - fs*dVdiat*r6xij/r6
      dv(8) = dv(8) - fs*dVdiat*r6yij/r6
      dv(9) = dv(9) - fs*dVdiat*r6zij/r6

      V4 = V4 + Vdiat

      call pothht_f14(rin(4), Vdiat, dVdiat)

      dv(4) = dv(4) + fs*dVdiat*r4xij/r4
      dv(5) = dv(5) + fs*dVdiat*r4yij/r4
      dv(6) = dv(6) + fs*dVdiat*r4zij/r4

      dv(10) = dv(10) - fs*dVdiat*r4xij/r4
      dv(11) = dv(11) - fs*dVdiat*r4yij/r4
      dv(12) = dv(12) - fs*dVdiat*r4zij/r4

      V4 = V4 + Vdiat

      V = V1*(1 - fs) + (V2 + V3 + V4)*fs
      dv = dv + dv1*(1 - fs)
!       print*,'dv depois 0',dv
      dv(1) = dv(1) + (v2 + v3 + v4 - v1)*(dfsdr3*r3xij/r3 + dfsdr5*r5xij/r5)
!        print*
!        print*, v2,v3,v4,v1,  dfsdr3 ,  r3xij , r3 , dfsdr5 , r5xij , r5

      dv(2) = dv(2) + (v2 + v3 + v4 - v1)*(dfsdr3*r3yij/r3 + dfsdr5*r5yij/r5)
      dv(3) = dv(3) + (v2 + v3 + v4 - v1)*(dfsdr3*r3zij/r3 + dfsdr5*r5zij/r5)

      dv(4) = dv(4) + (v2 + v3 + v4 - v1)*(dfsdr4*r4xij/r4 + dfsdr6*r6xij/r6)
      dv(5) = dv(5) + (v2 + v3 + v4 - v1)*(dfsdr4*r4yij/r4 + dfsdr6*r6yij/r6)
      dv(6) = dv(6) + (v2 + v3 + v4 - v1)*(dfsdr4*r4zij/r4 + dfsdr6*r6zij/r6)

      dv(7) = dv(7) + (v2 + v3 + v4 - v1)*(-dfsdr3*r3xij/r3 - dfsdr6*r6xij/r6)
      dv(8) = dv(8) + (v2 + v3 + v4 - v1)*(-dfsdr3*r3yij/r3 - dfsdr6*r6yij/r6)
      dv(9) = dv(9) + (v2 + v3 + v4 - v1)*(-dfsdr3*r3zij/r3 - dfsdr6*r6zij/r6)

      dv(10) = dv(10) + (v2 + v3 + v4 - v1)*(-dfsdr4*r4xij/r4 - dfsdr5*r5xij/r5)
      dv(11) = dv(11) + (v2 + v3 + v4 - v1)*(-dfsdr4*r4yij/r4 - dfsdr5*r5yij/r5)
      dv(12) = dv(12) + (v2 + v3 + v4 - v1)*(-dfsdr4*r4zij/r4 - dfsdr5*r5zij/r5)
!            print*,'dv depois  1',dv

      dv(1) = dv(1) + fs*(dv2*r2xij/r2)
      dv(2) = dv(2) + fs*(dv2*r2yij/r2)
      dv(3) = dv(3) + fs*(dv2*r2zij/r2)

      dv(4) = dv(4) + fs*(-dv2*r2xij/r2)
      dv(5) = dv(5) + fs*(-dv2*r2yij/r2)
      dv(6) = dv(6) + fs*(-dv2*r2zij/r2)

      dv(7) = dv(7) + fs*(dv3*r1xij/r1)
      dv(8) = dv(8) + fs*(dv3*r1yij/r1)
      dv(9) = dv(9) + fs*(dv3*r1zij/r1)

      dv(10) = dv(10) + fs*(-dv3*r1xij/r1)
      dv(11) = dv(11) + fs*(-dv3*r1yij/r1)
      dv(12) = dv(12) + fs*(-dv3*r1zij/r1)

!            print*,'dv depois  2',dv

   else!  end if

      v = v1

      dv = dv1

      group(i, 9) = 0

   end if

end subroutine

!  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!  _/                                                                _/
!  _/        surface 400s4369z H4 surface evaluation routines        _/
!  _/        (4004369.26), also known as the BMKP H4 surface         _/
!  _/        ---Arnold I. Boothroyd:  Version of 24 May 2002         _/
!  _/        ---Modified from original version of 6 July 2001        _/
!  _/            (for better extrapolation to short distances)       _/
!  _/                                                                _/
!  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!
!
!  Any QUESTIONS/PROBLEMS/COMMENTS concerning this programme can be
!  addressed to: boothroy@cita.utoronto.ca or pgmartin@cita.utoronto.ca
!
!
!**********************************************************************
!
!  Note that a short program to invoke and test the H4 surface routines
!  in this file are contained in the separate file  h4bmkp_test.f
!
!**********************************************************************
!--UPDATED 24 May 2002 (from the original version of 6 July 2001):
!   * Improved extrapolation to short distances, and interface provided
!  to allow user to choose between this improved extrapolation and the
!  original version of the surface --- also allows the user to choose
!  whether the program should stop on geometry-error or just return a
!  large value as an error flag.
!   * Fixed the asymptotic effect of epsilon in the London equation by
!  adding 1 microhartree, so that the energy of separated H + H + H + H
!  is now exactly zero (not -1 microhartree) --- this tiny increase in
!  the returned energy value is negligible wherever the use of the BMKP
!  H4 surface is appropriate (recall that, in the BMKP H4 surface, the
!  out H_2 + H_2 van der Waals well is cut off exponentially).  Note
!  that this change has no effect on the derivatives of the surface.
!**********************************************************************
!
!  ALL REAL VARIABLES ARE DOUBLE PRECISION (double precision)
!
! The H4 surface is continuous, as are its first and second derivatives
!  (except at internuclear distance r_i = 0.045156498889612436, where
!  the first derivative is continous but the second derivative is not).
!  If any distance r_i is too small (by default, less than 0.001 bohr),
!  then the surface is not evaluated: the derivatives are set to zero
!  and a large value (by default 999999.0) is returned instead.
!
!  Note that separated H + H + H + H is considered to have zero energy.
!
!----------------------------------------------------------------------
!
! NOTE that, as one or more of the 6 distances ri is reduced below
!  about 0.4 bohr, the original fitted surface can in some cases "turn
!  over" (and become smaller) after reaching a maximum value of about
!  0.9 hartree (this is a spurious effect due to extrapolation below
!  the shortest distances to which the surface was fitted); the fitted
!  surface can even become negative as the distance(s) are reduced
!  below about 0.3 bohr.  The NEW default (error level = 1: see below)
!  is to switch over to the Johnson-corrected London equation as the
!  shortest distance is reduced from 0.5 to 0.4 bohrs; but the user
!  can reset the error level to choose to handle short distances in
!  different ways, by calling the subroutine h4bmkp_set_err:
!
!    subroutine h4bmkp_set_err( level )
!    ==================================
!  Set the behavior on error, according to the (integer) input level:
!
!    level = 3 : if any distance r_i < 0.45 bohrs, it is considered to
!                be an error (and the behavior is the same as if an
!                invalid geometry was input to the h4bmkp_r interface):
!                derivatives are set to zero, and to flag this error a
!                value of Vfit = 99.0 is returned
!
!    level = 2 : if any distance r_i < 0.301 (or if it is an invalid
!                geometry, for h4bmkp_r), then return Vfit = 99.0
!          *** NOT RECOMMENDED: this case is the same as the behavior
!                of the BMKP H4 surface version prior to 24 May 2002.
!
!==> level = 1 (default) : if any r_i < 0.001 (or invalid geometry,
!                for h4bmkp_r), then return Vfit = 9999999.0;
!          *** AND ALSO FIX THE EXTRAPOLATION: switch over from fitted
!                surface to (Johnson-corrected) London surface as the
!                smallest of the distances r_i drops from 0.5 to 0.4
!
!    level = 0 : if r_i < 0.001 (or invalid), return Vfit = 9999999.0
!                but use fitted surface to extrapolate short distances
!                (NOT RECOMMENDED: spurious effects at short distances)
!
!    level = -1, -2, or -3 : the same limits as level_err = 1, 2, or 3,
!                respectively; but an error causes the program to STOP
!
!  Note that the current values of the error level, the minimum allowed
!  distance, and the error-return flag value can be checked via:
!
!    subroutine h4bmkp_ask_err( level, r_err, err_flag )
!    ===================================================
!
!----------------------------------------------------------------------
!
! NOTE that more than one interface is supplied here to the analytic
!  H4 BMKP surface; these subroutines are listed briefly below, and
!  described in more detail further below and in the comments at the
!  beginning of each subroutine (note that for the variables below
!  ideriv, ncalls, and iprint are integers, with all others being
!  double precision real variables):
!
!    subroutine h4bmkp( cc, Vfit, r, dVdcc, dVdr, Derr, ideriv )
!    ===========================================================
!  ---[general interface]: input cc(4,3), ideriv
!  ---output Vfit, r(6), dVdcc(4,3), dVdr(6), Derr
!
!    subroutine h4eval(cc,xmas,ideriv,ncalls,iprint,Vfit,dVdr,surfidp)
!    =================================================================
!  ---[old Keogh interface]: input cc(4,3),xmas(4),ideriv,ncalls,iprint
!  ---output Vfit,dVdr(6),surfidp
!
!    subroutine h4bmkp_cc( cc, Vfit, dVdcc, ideriv )
!    ===============================================
!  ---[minimal cartesian coordinate interface]: input cc(4,3), ideriv
!  ---output Vfit, dVdcc(4,3)
!
!    subroutine h4bmkp_r( r, Vfit, dVdr, Derr, ideriv )
!    ==================================================
!  ---[minimal internuclear-distance interface]: input r(6), ideriv
!  ---output Vfit, dVdr(6), Derr  <<<THIS INTERFACE NOT RECOMMENDED>>>
!
!----------------------------------------------------------------------
!
! NOTE that exactly-similar interfaces are also provided that return
!  the Johnson-corrected H4 London surface and its derivatives (with
!  the cusps rounded off by one microhartree), namely:
!
!  VlondJ_h4bmkp  VlondJ_h4eval  VlondJ_cc_h4bmkp  VlondJ_r_h4bmkp
!
!----------------------------------------------------------------------
!
! NOTE that the Schwenke H-H interaction energy may be obtained from:
!
!    subroutine vH2opt95_h4bmkp( rHH, E, ideriv )
!    ============================================
!  ---input rHH, ideriv
!  ---output E(3)
!    where, in the three-element array E:
!    E(1) returns the Schwenke H-H energy for H2-molecule size rHH;
!    if ideriv > 0, E(2) returns first derivative with respect to rHH;
!    E(3) is always set to zero.
!
!----------------------------------------------------------------------
!
!  NOTE on the accuracy of the derivatives of the BMKP H4 surface:
!  ***************************************************************
!
!  The analytic derivatives dVdcc() have been checked against numerical
!  derivatives (for 74903 of the fitted geometries --- all except the
!  very-compact non-ab-initio conformations), and agree to within the
!  accuracy of the numerical derivatives (typically, of order 1.E-10
!  hartrees/bohr, and never worse than about 1.E-5 hartrees/bohr, even
!  for H3 + H at the conical intersection).  Note that each geometry
!  was both tested 'as is' and also randomly rotated in coordinate
!  space, and reduced-size geometries were tested to confirm that the
!  derivatives were O.K. in the short-distance extrapolation region.
!
!  The analytic derivatives dVdr() were checked by using dr/dcc to
!  transform them into derivatives with respect to cartesian coords,
!  and comparing them to dVdcc().  Agreement is typically to within
!  1.E-9 hartrees/bohr, the worst case found among the above tested
!  geometries being 3.9E-5 hartrees/bohr.  NOTE that for some planar
!  and linear geometries, where the six interatomic distances r() are
!  not independent, it may be necessary to obtain dVdr derivatives
!  from from slightly shifted versions of the input coordinates, but
!  this does not significantly degrade the agreement in most cases
!  (the worst being the above-mentioned case of 3.9E-5 hartrees/bohr);
!  recall, however, that for such cases there is no unique correct
!  dVdr() vector, due to the non-independence of the distances r().
!  ---NOTE: between 16 and 184 additional surface evaluations may be
!  required when the shifted-geometry versions of dVdr() must be used.
!
!**********************************************************************
!
!
!               MORE DETAILED DESCRIPTION OF INTERFACES:
!               ****************************************
!
!
!     subroutine h4bmkp( cc, Vfit, r, dVdcc, dVdr, Derr, ideriv )
!     ===========================================================
!
!  H4 BMKP surface (A.I.Boothroyd, P.G.Martin): version of 22 May 2002
!
!  Input variables:
!  ================
!
!    cc(4,3) : Cartesian coordinates in bohrs: cc(i,j) = X_j for atom_i
!
!    ideriv : IF ideriv = 0 (or less), then compute only Vfit and r ;
!             IF ideriv = 1, then compute only Vfit, r, and dVdcc ;
!             IF ideriv = 2, then compute Vfit, r, dVdcc, dVdr, Derr
!               (NOTE that computing dVdr can increase the CPU-time
!               by up to two orders of magnitude for some planar
!               geometries, where a number of shifted-geometry cases
!               may be needed to get good dVdr values --- see below) ;
!             IF ideriv > 2, then do the same as for ideriv = 2, BUT:
!               if it is necessary to "fix up" the dVdr derivatives by
!               using slightly shifted geometries, then make a final
!               call with the unshifted geometry to reset the internal
!               common block values --- this is useful ONLY if you
!               access these H4 BMKP common blocks or internal H4 BMKP
!               subroutines (those called "<something>_h4bmkp(...)")
!               from your own routines outside this file.
!
!             If derivatives are not computed, they are set to zero.
!
!  Output variables:
!  =================
!
!    Vfit : the H4 surface energy in hartrees; NOTE that if any of the
!             internuclear distances is less than the specified value
!             of rmin_err (default 0.001: see subroutine h4bmkp_set_err
!             above), then it is an error, and Vfit is seet to the
!             value err_return, i.e., to 9999999.0 by default (or to
!             99.0, if the user has set the error level to 2 or -2).
!
!    r(6) : the 6 internuclear distances rAB, rAC, rAD, rBC, rBD, rCD
!
!    dVdcc(4,3) : the derivatives of Vfit in hartrees/bohr with respect
!                   to the 12 cartesian coordinates (see ideriv)
!
!    dVdr(6) : the derivatives of Vfit in hartrees/bohr with respect
!                to the 6 internuclear distances r() (see ideriv).
!              NOTE that the 6 distances r() are NOT independent for
!                planar or linear geometries, and two sets of dVdr()
!                that differ by large amounts (of order unity, or even
!                by orders of magnitude) may both agree equally well
!                (of order 1.E-5 or better) with the cartesian coord
!                derivatives dVdcc().
!
!    Derr : this gives the largest absolute error (in hartrees/bohr) in
!             the dVdr() values, as estimated by transforming them into
!             cartesian coordinate derivatives and comparing the values
!             with dVdcc().  Except for some planar geometries (most
!             often, where three or four of the H-atoms are in a line),
!             this difference Derr should be negligibly small.
!           ---For geometry-errors, Derr is set to err_return (999999.)
!
!----------------------------------------------------------------------
!
!
!     subroutine h4eval(cc,xmas,ideriv,ncalls,iprint,Vfit,dVdr,surfidp)
!     =================================================================
!
!  Bill Keogh's old H4 interface, now returning the BMKP H4 surface.
!
!  Input:  cc(4,3) : cartesian coords (bohrs): cc(i,j) = X_j for atom_i
!          xmas(4) : is ignored
!          ideriv : flag telling whether to compute the derivatives:
!                    ideriv = 0 : no derivatives are returned
!                    ideriv > 0 : derivatives dVdr() are returned, with
!                               respect to r12, r13, r14, r23, r24, r34
!          ncalls : is ignored
!          iprint : is ignored
!
!  Output: Vfit : H4 BMKP surface energy in hartrees
!          dVdr(6) : derivatives in hartrees/bohr (if ideriv > 0)
!          surfidp : numeric surface ID: returns a value of 4004369.26
!
!----------------------------------------------------------------------
!
!
!     subroutine h4bmkp_cc( cc, Vfit, dVdcc, ideriv )
!     ===============================================
!
!  Minimal cartesian coordinate interface to the H4 BMKP surface.
!
!  Input: cc(4,3) : Cartesian coords (bohrs): cc(i,j) = X_j for atom_i
!         ideriv : if ideriv > 0, then derivatives dVdcc() are returned
!
!  Output: Vfit : H4 BMKP surface energy in hartrees
!          dVdcc(4,3) : derivatives w.r.t. cartesian coords, in
!                         hartrees/bohr (only calculated if ideriv > 0)
!
!----------------------------------------------------------------------
!
!
!     subroutine h4bmkp_r( r, Vfit, dVdr, Derr, ideriv )
!     ==================================================
!
!  Minimal internuclear-distance interface to the H4 BMKP surface.
!
!  Input: r(6) : internuclear distances (bohrs) r12,r13,r14,r23,r24,r34
!         ideriv : if ideriv > 0, then derivatives dVdr() are returned
!
!  Output: Vfit : H4 BMKP surface energy in hartrees; if distances r()
!                   are invalid (do not correspond to any conformation
!                   of 4 H-atoms, or any distance is < rmin_err),
!                   then a value of Vfit = err_return is returned.
!          dVdr(6) : derivatives w.r.t. distances, in hartrees/bohr
!                         (only calculated if ideriv > 0)
!          Derr : returns err_return for invalid set of distances r();
!                   otherwise, gives the largest absolute dVdr() error
!
!
!**********************************************************************
!
subroutine h4bmkp_set_err(level)
!     ==================================
!----------------------------------------------------------------------
!  Set the error level (also controls short-distance extrapolation).
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
!__A large positive or negative value yields the default extrapolation.
!
   if (level .ge. 9) then
      level_err = 1
   else if (level .le. -9) then
      level_err = -1
   else
      level_err = max(-3, min(3, level))
   end if
!
!__Set the shortest allowed distance and the return error flag value.
!
   if (abs(level_err) .le. 1) then
      err_return = 9999999.d0
      rmin_err = 0.001d0
   else
      err_return = 99.d0
      if (abs(level_err) .le. 2) then
         rmin_err = 0.301d0
      else
         rmin_err = 0.45d0
      end if
   end if
!
   return
end
!
!**********************************************************************
!
subroutine h4bmkp_ask_err(level, r_err, err_flag)
!     ===================================================
!----------------------------------------------------------------------
!  Return the error level, error-distance, and error-return flag value.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
   level = level_err
   r_err = rmin_err
   err_flag = err_return
!
   return
end
!
!**********************************************************************
!
subroutine h4bmkp(ccin, Vfit, r, dVdcc, dVdr, Derr, ideriv)
!     ===========================================================
!----------------------------------------------------------------------
!  H4 BMKP surface (A.I.Boothroyd, P.G.Martin): version of 22 May 2002
!----------------------------------------------------------------------
!
!  Input variables:
!  ================
!
!    cc(4,3) : Cartesian coordinates in bohrs: cc(i,j) = X_j for atom_i
!
!    ideriv : IF ideriv = 0 (or less), then compute only Vfit and r ;
!             IF ideriv = 1, then compute only Vfit, r, and dVdcc ;
!             IF ideriv = 2, then compute Vfit, r, dVdcc, dVdr, Derr
!               (NOTE that computing dVdr can increase the CPU-time
!               by up to two orders of magnitude for some planar
!               geometries, where a number of shifted-geometry cases
!               may be needed to get good dVdr values --- see below) ;
!             IF ideriv > 2, then do the same as for ideriv = 2, BUT:
!               if it is necessary to "fix up" the dVdr derivatives by
!               using slightly shifted geometries, then make a final
!               call with the unshifted geometry to reset the internal
!               common block values --- this is useful ONLY if you
!               access these H4 BMKP common blocks or internal H4 BMKP
!               subroutines (those called "<something>_h4bmkp(...)")
!               from your own routines outside this file.
!
!             If derivatives are not computed, they are set to zero.
!
!  Output variables:
!  =================
!
!    Vfit : the H4 surface energy in hartrees; NOTE that if any of the
!             internuclear distances is less than the specified value
!             of rmin_err (default 0.001: see subroutine h4bmkp_set_err
!             above), then it is an error, and Vfit is seet to the
!             value err_return, i.e., to 9999999.0 by default (or to
!             99.0, if the user has set the error level to 2 or -2).
!
!    r(6) : the 6 internuclear distances rAB, rAC, rAD, rBC, rBD, rCD
!
!    dVdcc(4,3) : the derivatives of Vfit in hartrees/bohr with respect
!                   to the 12 cartesian coordinates (see ideriv)
!
!    dVdr(6) : the derivatives of Vfit in hartrees/bohr with respect
!                to the 6 internuclear distances r() (see ideriv).
!              NOTE that the 6 distances r() are NOT independent for
!                planar or linear geometries, and two sets of dVdr()
!                that differ by large amounts (of order unity, or even
!                by orders of magnitude) may both agree equally well
!                (of order 1.E-5 or better) with the cartesian coord
!                derivatives dVdcc().
!
!    Derr : this gives the largest absolute error (in hartrees/bohr) in
!             the dVdr() values, as estimated by transforming them into
!             cartesian coordinate derivatives and comparing the values
!             with dVdcc().  Except for some planar geometries (most
!             often, where three or four of the H-atoms are in a line),
!             this difference Derr should be negligibly small.
!           ---For geometry-errors, Derr is set to err_return (999999.)
!
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension ccin(4, 3), cc(4, 3), r(6), dVdcc(4, 3), dVdr(6)
!
!__COMMON /c_geometry_h4bmkp/ : geometry description variables:
!   geo(3,6) : for each of the three ways (ig=1,2,3) of dividing the four
!      H-atoms into two pairs (Schwenke's decomposition): geo(ig,1) = r_a ,
!      geo(ig,2) = r_b ,  geo(ig,3) = R ,  geo(ig,4) = cos(theta1) ,
!      geo(ig,5) = cos(theta2) ,  geo(ig,6) = cos(phi) ;  NOTE that some
!      additional necessary values are stowed in common /c_geo_sin/ .
!   dridcc(6,12) : derivatives of the six distances r() with respect to the
!      twelve cartesian coordinates, in same order as stored internally in the
!      matrix cc(), i.e., in the order:  cc(1,1), cc(2,1), cc(3,1), cc(4,1),
!      cc(1,2), cc(2,2), cc(3,2), cc(4,2), cc(1,3), cc(2,3), cc(3,3), cc(4,3)
!   dgeodcc(3,6,12) : derivatives of Schwenke's decomposition values, in
!      geo(3,6), with respect to the twelve cartesian coordinates cc().
!   fSA(3) : relative "goodness" values for the three decompositions in geo().
!   sumVp(3) : Schwenke-type selector functions (using pair-wise anti-bonding
!      potentials), for each of the three decompositions in geo().
!   dfSAdcc(3,12) : derivatives of the relative "goodness" values for the
!      three decompositions in geo(), with respect to the cartesian coords.
!   dfSAdri(3,6) : derivatives of the relative "goodness" values for the
!      three decompositions in geo(), with respect to the distances r().
!   dsumVpdri(3,6) : derivatives of the three Schwenke-type selector functions
!      in sumVp(), with respect to the distances r().
!   Yy(3,6) : the Schwenke-type "curly-Y" angular functions, composed from Ylm.
!   dYydgeo(3,6,6) : derivatives of the "curly-Y" angular functions, with
!      respect to Schwenke's decomposition variables geo() of the same way ig:
!      dYydgeo(ig,i,j) = d Y(ig,i) / d geo(ig,j).
!   dgeodri(3,6,6) :  derivatives of Schwenke's decomposition values, in
!      geo(3,6), with respect to the distances r().
!   jco_ri(6,3) : the values jco_ri(*,ig) are the distance re-orderings needed
!      for input to Schwenke/Keogh routine that computes dgeodri(ig,*,*).
!   jco_geo(6) : re-orderings of geo(ig,*) variables for the equivalent ones
!      used in the Schwenke/Keogh routine that computes dgeodri(ig,*,*).
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!                              ! max allowed diff dVdr() vs. dVdcc()
   parameter(tol_Derr=1.d-6)
!
!__Get the H4 surface energy and any required derivatives.
!
   Derr = 0.d0
!

   cc = ccin/5.2917721092d-11

   call surftry_h4bmkp(cc, Vfit, r, dVdcc, dVdr, ideriv)
!                                                ! geom error?
   if (Vfit .eq. err_return) then
      if (level_err .lt. 0) stop ' STOP -- H4 distance error. '
      Derr = err_return
      return
   end if
!
!__If derivatives with respect to internuclear distances r() were required:
!
   if (ideriv .gt. 1) then
!
!__check them against derivatives with respect to cartesian coordiantes,
!
      k = 0
      do jj = 1, 3
      do ii = 1, 4
         k = k + 1
         dcc = 0.d0
         do j = 1, 6
            dcc = dcc + dVdr(j)*dridcc(j, k)
         end do
         Derr = max(Derr, abs(dVdcc(ii, jj) - dcc))
      end do
      end do
!
!__and if there is a significant mismatch, attempt to do better (by getting
!  derivatives at slightly shifted versions of the input conformation).
!
      if (Derr .gt. tol_Derr) &
         call fixdVdr_h4bmkp(cc, r, dVdcc, dVdr, Derr, ideriv)
!
   end if
!
   Vfit = Vfit*4.3597482D-18
   dVdr = dVdr*4.3597482D-18/5.2917721092d-11
   dvdcc = dvdcc*4.3597482D-18/5.2917721092d-11
!      print*,'cc'
!      print*,cc

!      print*,'Vfit=',Vfit

!      print*,'r'
!      print*,r

!      print*,'dVdcc'
!      print*,dVdcc

!      print*,'dVdr'
!      print*,dVdr

!      print*,'Derr'
!      print*,Derr

   return
end
!
!**********************************************************************
!

!**********************************************************************
!
subroutine h4bmkp_cc(cc, Vfit, dVdcc, ideriv)
!     ===============================================
!----------------------------------------------------------------------
!  Minimal cartesian coordinate interface to the H4 BMKP surface.
!----------------------------------------------------------------------
!
!  Input: cc(4,3) : Cartesian coords (bohrs): cc(i,j) = X_j for atom_i
!         ideriv : if ideriv > 0, then derivatives dVdcc() are returned
!
!  Output: Vfit : H4 BMKP surface energy in hartrees
!          dVdcc(4,3) : derivatives w.r.t. cartesian coords, in
!                         hartrees/bohr (only calculated if ideriv > 0)
!
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension cc(4, 3), dVdcc(4, 3), r(6), dVdr(6)
!
   call h4bmkp(cc, Vfit, r, dVdcc, dVdr, Derr, min(ideriv, 1))
!
   return
end
!
!**********************************************************************
!
subroutine h4bmkp_r(rin, Vfit, dVdr, Derr, ideriv)
!     ==================================================
!----------------------------------------------------------------------
!  Minimal internuclear-distance interface to the H4 BMKP surface.
!----------------------------------------------------------------------
!
!  Input: r(6) : internuclear distances (bohrs) r12,r13,r14,r23,r24,r34
!         ideriv : if ideriv > 0, then derivatives dVdr() are returned
!
!  Output: Vfit : H4 BMKP surface energy in hartrees; if distances r()
!                   are invalid (do not correspond to any conformation
!                   of 4 H-atoms, or any distance is < rmin_err),
!                   then a value of Vfit = err_return is returned.
!          dVdr(6) : derivatives w.r.t. distances, in hartrees/bohr
!                         (only calculated if ideriv > 0)
!          Derr : returns err_return for invalid set of distances r();
!                   otherwise, gives the largest absolute dVdr() error
!
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(6), dVdr(6), rin(6)
!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
   dimension cc(4, 3), dVdcc(4, 3), ri(6), dVdri(6), iri(6)
!
!__Get cartesian coordinates from the six distances, checking that the
!  distances correspond to a valid geometry:
!
   r = rin!/5.2917721092d-11
!print*,r
   call getccofr_h4bmkp(r, cc, ivalid, errdr, iri)
!

   if (ivalid .lt. 0) then
      if (level_err .lt. 0) stop ' STOP -- H4 distance error. '
      Vfit = err_return
      Derr = err_return
      return
   end if
!
!      print*,'ideriv=',ideriv
   if (ideriv .le. 0) then
      ider_use = 0
   else
      ider_use = max(ideriv, 2)
   end if
!
!__Get the H4 BMKP energy (and derivatives, if necessary):
!
!      print*,'cc=',cc

   call h4bmkp(cc, Vfit, ri, dVdcc, dVdri, Derr, ider_use)
!      PRINT*,'ri',ri
!
!__Reorder the derivatives, to correspond to the input distances:
!
!      print*,'Vfit=',Vfit
!      Vfit=Vfit*4.3597482D-18

   do i = 1, 6
      dVdr(iri(i)) = dVdri(i)!*4.3597482D-18/5.2917721092d-11
   end do
!
!      print*,Dvdr
   return
end
!
!**********************************************************************
!
subroutine VlondJ_h4bmkp(cc, Vlond, r, dVdcc, dVdr, Derr, ideriv)
!     ===========================================================
!----------------------------------------------------------------------
!  Johnson-corrected H4 London equation (1 microhartree cusp-rounding).
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension cc(4, 3), r(6), dVdcc(4, 3), dVdr(6)
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
   common/c_rcl_h4bmkp/rlo_cl, rhi_cl, rmid_cl, dinv_cl, oinv_cl
   save/c_rcl_h4bmkp/
!
   common/c_init_h4bmkp/init_h4bmkp, iread_h4bmkp
   save/c_init_h4bmkp/
!
!__There will never be a problem with the London-equation derivatives:
!
   Derr = 0.d0
!
!__On first call to this routine, do all necessary initializations:
!
   if (init_h4bmkp .eq. 0) call parinit_h4bmkp
!
!__This is a fudge to set fnlclose to zero and get the London equation:
!
   levtmp = level_err
   rtmp = rlo_cl
   if (levtmp .ge. 0) then
      level_err = 1
   else
      level_err = -1
   end if
   rlo_cl = 1.d36
!
!__Get the London energy (and its derivatives, if necessary) by this fudge:
!
   call surftry_h4bmkp(cc, Vlond, r, dVdcc, dVdr, ideriv)
!
!__Restore the fudged values:
!
   level_err = levtmp
   rlo_cl = rtmp
!
!__If level_err < 0, then one should stop on error:
!
   if (Vlond .eq. err_return .and. level_err .lt. 0) &
      stop ' STOP -- H4 distance error. '
!
   return
end
!
!**********************************************************************
!

!**********************************************************************
!
subroutine VlondJ_cc_h4bmkp(cc, Vlond, dVdcc, ideriv)
!     =======================================================
!----------------------------------------------------------------------
!  Minimal cartesian coordinate Johnson-corrected H4 London equation.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension cc(4, 3), dVdcc(4, 3), r(6), dVdr(6)
!
   call VlondJ_h4bmkp(cc, Vlond, r, dVdcc, dVdr, Derr, min(ideriv, 1))
!
   return
end
!
!**********************************************************************
!
subroutine VlondJ_r_h4bmkp(r, Vlond, dVdr, Derr, ideriv)
!     ==========================================================
!----------------------------------------------------------------------
!  Minimal internuclear-distance Johnson-corrected H4 London equation.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(6), dVdr(6)
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   common/c_fnlclose_h4bmkp/fnlclose, dfnlclosedri(6)
   save/c_fnlclose_h4bmkp/
!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
   common/c_rcl_h4bmkp/rlo_cl, rhi_cl, rmid_cl, dinv_cl, oinv_cl
   save/c_rcl_h4bmkp/
!
   common/c_init_h4bmkp/init_h4bmkp, iread_h4bmkp
   save/c_init_h4bmkp/
!
   dimension cc(4, 3), dVlondfSA(3), iri(6)
!
!__On first call to this routine, do all necessary initializations:
!
   if (init_h4bmkp .eq. 0) call parinit_h4bmkp
!
!__Check that the six distances correspond to a valid geometry:
!
   call getccofr_h4bmkp(r, cc, ivalid, errdr, iri)
!
   if (ivalid .lt. 0) then
      if (level_err .lt. 0) stop ' STOP -- H4 distance error. '
      Vfit = err_return
      Derr = err_return
      return
   end if
!
!__There will never be a problem with the London-equation derivatives:
!
   Derr = 0.d0
!
!__Set fnlclose and fSA() to zero:
!
   fnlclose = 0.d0
   do i = 1, 6
      dfnlclosedri(i) = 0.d0
   end do
   fSA(1) = 0.d0
   fSA(2) = 0.d0
   fSA(3) = 0.d0
!
!__Get the London energy (and its derivatives, if necessary):
!
   call comlon_h4bmkp(r, ideriv, Vlond, dVdr, dVlondfSA)
!
   return
end
!
!**********************************************************************
!
subroutine fixdVdr_h4bmkp(cc, r, dVdcc, dVdr, Derr, ideriv)
!----------------------------------------------------------------------
!  Try to fix up erroneous derivatives dVdr() via coordinate shifts.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension cc(4, 3), r(6), dVdcc(4, 3), dVdr(6)
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   common/c_rcalc_indices_h4bmkp/ipart1(6), ipart2(6)
   save/c_rcalc_indices_h4bmkp/
!
   common/c_reverse_indices_h4bmkp/ir_ij(4, 4)
   save/c_reverse_indices_h4bmkp/
!
   common/c_cross_h4bmkp/pt(2, 3), ptm(2)
   save/c_cross_h4bmkp/
!
   parameter(v_big=1.d+36, v_bigger=1.01d+36)
!
   parameter(tol_Derr=1.d-6, tol_pt=1.d-7, dsh_lo=1.d-7, &
             dsh_himin=100.d0*dsh_lo)
!
   dimension rvp(2, 3), RRv(3), dpdrvp(2, 3), dpdRRv(3)
!
   dimension dridcc_0(6, 12), dVdcc_0(12), dcci(3), ccsh(4, 3), &
      rsh(6), dVdrsh(6, 9, 3), dVdccsh(12), dersh(9, 3)
!
   dimension ia(4)
!
   common/c_fixdVdr_h4bmkp/i1sh(5:8), i2sh(5:8)
   save/c_fixdVdr_h4bmkp/
!-debug[
   common/c_debug_h4bmkp/i_debug_h4bmkp, kloop, derr_0
   save/c_debug_h4bmkp/
!
!      data i_debug_h4bmkp/0/
!-debug]
!
!
!__First, get a direction dcci() along which to shift one or two of the atoms.
!  For a planar or near-planar geometry, this should be perpendicular to the
!  plane; for a linear geometry, it should be perpendicular to the line.
!
!                  ! find endpoints of the longest of the 6 distances r()
   rmax = r(1)
   klong = 1
   do k = 2, 6
      if (r(k) .gt. r(klong)) klong = k
   end do
   ia(1) = ipart1(klong)
   ia(4) = ipart2(klong)
!                     ! pair each endpoint with one of the remaining atoms
   knot = 7 - klong
   ia(2) = ipart1(knot)
   ia(3) = ipart2(knot)
   k12 = ir_ij(ia(1), ia(2))
   k13 = ir_ij(ia(1), ia(3))
   k24 = ir_ij(ia(2), ia(4))
   k34 = ir_ij(ia(3), ia(4))
!                                               ! (the closer one)
   if (r(k13) .le. r(k12) .and. r(k24) .le. r(k34)) then
      ia(3) = ipart1(knot)
      ia(2) = ipart2(knot)
      k12 = ir_ij(ia(1), ia(2))
      k13 = ir_ij(ia(1), ia(3))
      k24 = ir_ij(ia(2), ia(4))
      k34 = ir_ij(ia(3), ia(4))
   end if
!                  ! get cross products of pairs vs. longest distance
   do i = 1, 3
      RRv(i) = cc(ia(4), i) - cc(ia(1), i)
      rvp(1, i) = cc(ia(2), i) - cc(ia(1), i)
      rvp(2, i) = cc(ia(4), i) - cc(ia(3), i)
   end do
!
   call cross_h4bmkp(rvp, RRv, cosphi, dpdrvp, dpdRRv, 0, iden)
!
!                               ! choose the longer cross product vector
   if (ptm(1) .ge. ptm(2)) then
      ipt = 1
   else
      ipt = 2
   end if
!                              ! and normalize it, if it is non-zero
   if (ptm(ipt) .ge. tol_pt) then
   do i = 1, 3
      dcci(i) = pt(ipt, i)/ptm(ipt)
   end do
!               ! or, if it is zero, get a perpendicular to longest distance
   else
   if (abs(RRv(1)) .ge. max(abs(RRv(2)), abs(RRv(3)))) then
      j1 = 1
      if (abs(RRv(2)) .ge. abs(RRv(3))) then
         j2 = 2
         j3 = 3
      else
         j2 = 3
         j3 = 2
      end if
   else if (abs(RRv(2)) .ge. abs(RRv(3))) then
      j1 = 2
      if (abs(RRv(1)) .ge. abs(RRv(3))) then
         j2 = 1
         j3 = 3
      else
         j2 = 3
         j3 = 1
      end if
   else
      j1 = 3
      if (abs(RRv(1)) .ge. abs(RRv(2))) then
         j2 = 1
         j3 = 2
      else
         j2 = 2
         j3 = 1
      end if
   end if
   d = sqrt(RRv(j1)**2 + RRv(j2)**2)
   if (d .gt. 0.d0) then
      dcci(j1) = -RRv(j2)/d
      dcci(j2) = RRv(j1)/d
   else
      dcci(j1) = -1.d0
      dcci(j2) = 1.d0
   end if
   dcci(j3) = 0.d0
   end if
!-debug[
   derr_0 = Derr
   if (i_debug_h4bmkp .gt. 0) then
      write (66, '("+++++++ dcci",3f15.10,"  initial d,D:")') &
         dcci(1), dcci(2), dcci(3)
      write (66, '("   d",1p,e25.17," D",6e25.17)') &
         Derr, (dVdr(i), i=1, 6)
   end if
!-debug]
!
!__Store the original derivatives of E and of r() w.r.t. cartesian coords.
!
   j = 0
   do jj = 1, 3
   do ii = 1, 4
      j = j + 1
      dVdcc_0(j) = dVdcc(ii, jj)
      ccsh(ii, jj) = cc(ii, jj)
   end do
   end do
!
   do i = 1, 6
   do j = 1, 12
      dridcc_0(i, j) = dridcc(i, j)
   end do
   end do
!
!__Try shifting one or two of the atoms and getting "shifted" derivatives.
!
   dsh_at = dsh_lo
   dershmin_1 = v_big
   dershmin_2 = v_bigger
   dsh_hi = max(dsh_himin, min(1.d0, 0.001d0*r(klong), 0.1d0 &
                               *min(r(knot), r(k12), r(k13), r(k24), r(k34))))
   nok_1 = 0
   nok_2 = 0
   kloop = 0
   kloop_hi = 30
!
!__Try shift sizes starting from 1.E-7 bohrs, but no higher than a maximum
!  size of  min{ 0.001 * r(longest) , 0.1 * r(shortest) , 1.0 bohr }  ,
!  increasing the shift size by a factor of 2 each time; also stop if dVdr
!  error is no longer decreasing and has been less than 1.E-6 hartree/bohr
!  for two separate shift sizes.
!
   do while (dsh_at .le. dsh_hi .and. (nok_1 .eq. 0 .or. &
                                       nok_2 .eq. 0 .or. dershmin_2 .lt. dershmin_1) .and. &
             (kloop .lt. kloop_hi .or. dershmin_1 .gt. 1.d-3))
!
      kloop = kloop + 1
!                        ! if just-computed dVdr() vector is better than
!                        ! stored dVdr(), then store it instead
!
      if (kloop .eq. 2 .or. dershmin_2 .lt. dershmin_1) then
         dsh_1 = dsh_2
         dershmin_1 = dershmin_2
         nok_1 = nok_2
         do k = 1, 9
            dersh(k, 1) = dersh(k, 2)
            do j = 1, 6
               dVdrsh(j, k, 1) = dVdrsh(j, k, 2)
            end do
         end do
         if (dershmin_1 .lt. 1.d-5) then
            kloop_hi = kloop + 5
         else if (dershmin_1 .lt. 1.d-4) then
            kloop_hi = kloop + 7
         else if (dershmin_1 .lt. 1.d-3) then
            kloop_hi = kloop + 10
         else
            kloop_hi = 30
         end if
      end if
!                  ! with the new shift size  dsh_at,
      nok_2 = 0
      dsh_2 = dsh_at
!                  ! get derivadives dVdr() at 8 different shift types:
      do k = 1, 8
!                              ! get coord shifts of current type
         if (k .le. 4) then
            ccsh(k, 1) = cc(k, 1) + dcci(1)*dsh_at
            ccsh(k, 2) = cc(k, 2) + dcci(2)*dsh_at
            ccsh(k, 3) = cc(k, 3) + dcci(3)*dsh_at
         else
            i1 = ia(i1sh(k))
            i2 = ia(i2sh(k))
            ccsh(i1, 1) = cc(i1, 1) + dcci(1)*dsh_at
            ccsh(i1, 2) = cc(i1, 2) + dcci(2)*dsh_at
            ccsh(i1, 3) = cc(i1, 3) + dcci(3)*dsh_at
            ccsh(i2, 1) = cc(i2, 1) + dcci(1)*dsh_at
            ccsh(i2, 2) = cc(i2, 2) + dcci(2)*dsh_at
            ccsh(i2, 3) = cc(i2, 3) + dcci(3)*dsh_at
         end if
!                                                ! get dVdr
         call surftry_h4bmkp(ccsh, Esh, rsh, dVdccsh, &
                             dVdrsh(1, k, 2), 2)
!                            ! and its error relative to unshifted dVdcc
         dersh(k, 2) = 0.d0
         do j = 1, 12
            dcc = 0.d0
            do i = 1, 6
               dcc = dcc + dVdrsh(i, k, 2)*dridcc_0(i, j)
            end do
            dersh(k, 2) = max(dersh(k, 2), abs(dVdcc_0(j) - dcc))
         end do
         if (k .eq. 1) then
            dershmin_2 = dersh(1, 2)
         else
            dershmin_2 = min(dershmin_2, dersh(k, 2))
         end if
!                                                ! in tolerance?
         if (dersh(k, 2) .le. tol_Derr) nok_2 = nok_2 + 1
!                                                ! reset coords
         if (k .le. 4) then
            ccsh(k, 1) = cc(k, 1)
            ccsh(k, 2) = cc(k, 2)
            ccsh(k, 3) = cc(k, 3)
         else
            ccsh(i1, 1) = cc(i1, 1)
            ccsh(i1, 2) = cc(i1, 2)
            ccsh(i1, 3) = cc(i1, 3)
            ccsh(i2, 1) = cc(i2, 1)
            ccsh(i2, 2) = cc(i2, 2)
            ccsh(i2, 3) = cc(i2, 3)
         end if
!
      end do
!                  ! get a 9th dVdr() by averaging best of the above 8
      do i = 1, 6
         dVdrsh(i, 9, 2) = 0.d0
      end do
      nuse = 0
      dlim = 10.d0*dershmin_2
      do k = 1, 8
      if (dersh(k, 2) .le. dlim) then
         nuse = nuse + 1
         do i = 1, 6
            dVdrsh(i, 9, 2) = dVdrsh(i, 9, 2) + dVdrsh(i, k, 2)
         end do
      end if
      end do
      do i = 1, 6
         dVdrsh(i, 9, 2) = dVdrsh(i, 9, 2)/nuse
      end do
      dersh(9, 2) = 0.d0
      do j = 1, 12
         dcc = 0.d0
         do i = 1, 6
            dcc = dcc + dVdrsh(i, 9, 2)*dridcc_0(i, j)
         end do
         dersh(9, 2) = max(dersh(9, 2), abs(dVdcc_0(j) - dcc))
      end do
!                                                ! in tolerance?
      if (dersh(9, 2) .le. tol_Derr) nok_2 = nok_2 + 1
!
!-debug[
      if (i_debug_h4bmkp .gt. 0) then
         write (66, '("+++++++ kloop",i4," dx,dmin",1p,2e25.17)') &
            kloop, dsh_2, dershmin_2
         do k = 1, 9
            write (66, '(i2," d",1p,e25.17," D",6e25.17)') &
               k, dersh(k, 2), (dVdrsh(i, k, 2), i=1, 6)
         end do
      end if
!-debug]
!                        ! double size of coord shift, for the next loop
      dsh_at = dsh_at*2.d0
!
   end do
!
!__Try extrapolation to zero shift (this generally does not work very well),
!  and store the best set of derivatives dVdr() to be returned
!
   f_2 = -dsh_1/(dsh_2 - dsh_1)
   f_1 = 1.d0 - f_2
!
   do k = 1, 9
!                  ! extrapolate this set of dVdr() to zero coord shift
      do i = 1, 6
         dVdrsh(i, k, 3) = f_1*dVdrsh(i, k, 1) + f_2*dVdrsh(i, k, 2)
      end do
!                        ! and get its error relative to dVdcc()
      dersh(k, 3) = 0.d0
      do j = 1, 12
         dcc = 0.d0
         do i = 1, 6
            dcc = dcc + dVdrsh(i, k, 3)*dridcc_0(i, j)
         end do
         dersh(k, 3) = max(dersh(k, 3), abs(dVdcc_0(j) - dcc))
      end do
!                  ! check for improved errors among current shift type
      do j = 1, 3
      if (dersh(k, j) .lt. Derr) then
         Derr = dersh(k, j)
         do i = 1, 6
            dVdr(i) = dVdrsh(i, k, j)
         end do
      end if
      end do
!
   end do
!-debug[
   if (i_debug_h4bmkp .gt. 0) then
      write (66, '("******* k_mat",i2," dx,dmin",1p,2e25.17)') &
         1, dsh_1, dershmin_1
      do k = 1, 9
         write (66, '(i2," d",1p,e25.17," D",6e25.17)') &
            k, dersh(k, 1), (dVdrsh(i, k, 1), i=1, 6)
      end do
      write (66, '("******* k_mat",i2," dx,dmin",1p,2e25.17)') &
         2, dsh_2, dershmin_2
      do k = 1, 9
         write (66, '(i2," d",1p,e25.17," D",6e25.17)') &
            k, dersh(k, 2), (dVdrsh(i, k, 2), i=1, 6)
      end do
      write (66, '("******* k_mat",i2," f_1,f_2",1p,2e25.17)') &
         3, f_1, f_2
      do k = 1, 9
         write (66, '(i2," d",1p,e25.17," D",6e25.17)') &
            k, dersh(k, 3), (dVdrsh(i, k, 3), i=1, 6)
      end do
      write (66, '("******* BEST returned:")')
      write (66, '("   d",1p,e25.17," D",6e25.17)') &
         Derr, (dVdr(i), i=1, 6)
      write (66, '("*******")')
   end if
!-debug]
!
!__If ideriv = 3 or more, then make a final call to the surface routine with
!  the original (unshifted) coordinates (discarding resulting dVdr values),
!  to reset internal common blocks to the values appropriate to these
!  (unshifted) coordinates.  This is useful ONLY if you later access any of
!  these common blocks (or internal H4 BMKP subroutines) from OUTSIDE of the
!  H4 BMKP routines in this file.
!
   if (ideriv .gt. 2) then
      call surftry_h4bmkp(cc, Esh, rsh, dVdccsh, dVdrsh(1, 1, 1), 2)
   end if
!
   return
end
!
!**********************************************************************
!
subroutine surftry_h4bmkp(cc, Vfit, r, dVdcc, dVdr, ideriv)
!----------------------------------------------------------------------
!  Compute the H4 BMKP surface for the given cartesian coordinates.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)
!
   dimension cc(4, 3), r(6), dVdcc(4, 3), dVdr(6)
!
   parameter(nterms=8, maxl=7*224 + 1350)
   common/c_coef_h4bmkp/coef(maxl), vsca(nterms), beta(nterms), &
      pow(nterms), iv(nterms), imin(nterms), imax(nterms), &
      kmin(nterms), kmax(nterms), irho(nterms), mlin, njac, &
      ikmin_lo, ikmax_hi
   save/c_coef_h4bmkp/
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   common/c_init_h4bmkp/init_h4bmkp, iread_h4bmkp
   save/c_init_h4bmkp/
!
!-old;      dimension dVlondri(6), dVlondfSA(3), dVnsbdri(6),
!-old;     $     dVlindri(3,6), dVlindgeo(3,6), dVmbedri(6), Vlinear(3)
!
   dimension dVlondri(6), dVlondfSA(3), dVnsbdri(6), &
      dVlindri(3, 6), dVlindgeo(3, 6), dVmbedri(6)
!
   common/c_Vparts_h4bmkp/Vlinear(3), Vlinmbe
   save/c_Vparts_h4bmkp/
!
   common/c_fnlclose_h4bmkp/fnlclose, dfnlclosedri(6)
   save/c_fnlclose_h4bmkp/
!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
!__On first call to this routine, do all necessary initializations
!
   if (init_h4bmkp .eq. 0) call parinit_h4bmkp
!
!                        ! always initialize derivatives to zero
   do i = 1, 4
      dVdcc(i, 1) = 0.d0
      dVdcc(i, 2) = 0.d0
      dVdcc(i, 3) = 0.d0
   end do
   do i = 1, 6
      dVdr(i) = 0.d0
   end do
!
!__Get the ra,rb,R geometry; check for error (any distance < rmin_err):
!
   call geomet_h4bmkp(r, cc, ideriv, ierr)
!                                    ! if distance too small, return
   if (ierr .gt. 0) then
      if (level_err .lt. 0) stop ' STOP -- H4 distance error. '
      Vfit = err_return

      return
   end if
!                              ! if we need more than the London eqn:
   if (fnlclose .gt. 0.d0) then
!
!__Get the geometry fractions (goodness factors) fSA:
!
      call vcorrec_h4bmkp(ideriv)
!
!__Calculate sumVp (used in Vcorr):
!
      call vpcalc_h4bmkp(r, ideriv)
!
!__Work out the Yy values (from spherical harmonics) for the 3 geometries:
!
      call getYy_h4bmkp(ideriv)
!
   end if
!
!__Calculate VlonH4, using composite London surface, possibly using fSA
!
   call comlon_h4bmkp(r, ideriv, VlonH4, dVlondri, dVlondfSA)
!
!                              ! if we need only the London eqn:
   if (fnlclose .eq. 0.d0) then
!
      Vnonlondon = 0.d0
      Vnsb = 0.d0
      do ig = 1, 3
         fSA(ig) = 0.d0
         Vlinear(ig) = 0.d0
      end do
      Vlinmbe = 0.d0
!                        ! else, if we need more than the London eqn:
   else
!
!__Calculate Vnsb, the H3 contributions (from the BKMP2 H3 terms)
!
      call Vnsb_h4bmkp(r, ideriv, Vnsb, dVnsbdri)
!
!__Get the linear-parameter correction terms, including any MBE term:
!
      call Vlinear_h4bmkp(r, ideriv, Vlinear, Vlinmbe, &
                          dVlindri, dVlindgeo, dVmbedri)
!
!__Add up all the pieces of the potential:
!
      Vnonlondon = Vnsb + fSA(1)*Vlinear(1) &
                   + fSA(2)*Vlinear(2) + fSA(3)*Vlinear(3) + Vlinmbe
!
   end if
!
   Vfit = VlonH4 + Vnonlondon*fnlclose
!
!__Calculate the derivatives, if they are required:
!
   if (ideriv .gt. 0) then
!                              ! usual case: all terms contribute:
      if (fnlclose .gt. 0.d0) then
!
!__Calculate the derivatives with respect to the Cartesian coordinates
!
         j = 0
         do jj = 1, 3
         do ii = 1, 4
            j = j + 1
            dcc = 0.d0
            do ig = 1, 3
            if (fSA(ig) .gt. 0.d0) then
               dcc = dcc + (Vlinear(ig)*fnlclose &
                            + dVlondfSA(ig))*dfSAdcc(ig, j)
               do k = 1, 6
                  dcc = dcc &
                        + (dVlindgeo(ig, k)*dgeodcc(ig, k, j) &
                           + dVlindri(ig, k)*dridcc(k, j)) &
                        *fSA(ig)*fnlclose
               end do
            end if
            end do
            do k = 1, 6
               dcc = dcc + (dVlondri(k) &
                            + (dVnsbdri(k) + dVmbedri(k))*fnlclose &
                            + Vnonlondon*dfnlclosedri(k)) &
                     *dridcc(k, j)
            end do
            dVdcc(ii, jj) = dcc
         end do
         end do
!
!__Calculate the derivatives with respect to the 6 interatomic distances
!
         if (ideriv .gt. 1) then
!
            do j = 1, 6
               dcc = 0.d0
               do ig = 1, 3
               if (fSA(ig) .gt. 0.d0) then
                  dcc = dcc + (Vlinear(ig)*fnlclose &
                               + dVlondfSA(ig))*dfSAdri(ig, j) &
                        + dVlindri(ig, j)*fSA(ig)*fnlclose
                  do k = 1, 6
                     dcc = dcc + dVlindgeo(ig, k) &
                           *dgeodri(ig, k, j)*fSA(ig)*fnlclose
                  end do
               end if
               end do
               dVdr(j) = dcc + dVlondri(j) &
                         + (dVnsbdri(j) + dVmbedri(j))*fnlclose &
                         + Vnonlondon*dfnlclosedri(j)
            end do
!
         end if
!
!__If abs(level_err) = 1 and the smallest distance is 0.4 bohrs or less,
!  then fnlclose = 0.0: only the Johnson-corrected London equation is used:
!
      else
!                                    ! get cartesian derivatives
         j = 0
         do jj = 1, 3
         do ii = 1, 4
            j = j + 1
            dcc = 0.d0
            do k = 1, 6
               dcc = dcc + dVlondri(k)*dridcc(k, j)
            end do
            dVdcc(ii, jj) = dcc
         end do
         end do
!                                    ! get distance derivatives
         if (ideriv .gt. 1) then
         do j = 1, 6
            dVdr(j) = dVlondri(j)
         end do
         end if
!
      end if
!
   end if
!
!  The following is necessary to avoid an IRIX-mips compiler bug; without it,
!  derivatives return values of "not-a-number" for optimization -O0 or -O1.
!
   vr_tmp = dVdr(1)
   vc_tmp = dVdcc(1, 1)
!

!      Vfit = Vfit*4.3597482D-18

   return
end
!
!**********************************************************************
!
block data coef_h4bmkp
!----------------------------------------------------------------------
!  Initial values of common-block variables, where required.
!----------------------------------------------------------------------
!
implicit double precision(a - h, o - z)
implicit integer(i - n)
!                                          ! initialization flags
common/c_init_h4bmkp/init_h4bmkp, iread_h4bmkp
save/c_init_h4bmkp/
!                                    ! main surface coefficients:
common/c_debug_h4bmkp/i_debug_h4bmkp, kloop, derr_0
save/c_debug_h4bmkp/
!
data i_debug_h4bmkp/0/
parameter(n_p=788)
!

common/c_fixdVdr_h4bmkp/i1sh(5:8), i2sh(5:8)
save/c_fixdVdr_h4bmkp/
!

data i1sh/1, 1, 2, 2/, i2sh/3, 4, 3, 4/

parameter(nterms=8, maxl=7*224 + 1350)
common/c_coef_h4bmkp/coef(maxl), vsca(nterms), beta(nterms), &
   pow(nterms), iv(nterms), imin(nterms), imax(nterms), &
   kmin(nterms), kmax(nterms), irho(nterms), mlin, njac, &
   ikmin_lo, ikmax_hi
save/c_coef_h4bmkp/
!                                          ! mbe beta and indices:
parameter(max_hi=12, maxsum_hi=13, &
          kind3hi=max_hi*max_hi, kind4hi=kind3hi*max_hi, &
          kind5hi=3000, kind6hi=3000)
common/c_mbe_index_h4bmkp/beta_p, mbe_perm(6, 24), nterms_use, &
   indmin, indmax, maxsum, k1_lo, k1_hi, k2_lo(max_hi), &
   k2_hi(max_hi), k3_lo(kind3hi), k3_hi(kind3hi), &
   k4_lo(kind4hi), k4_hi(kind4hi), k5_lo(kind5hi), &
   k5_hi(kind5hi), k6_lo(kind6hi), k6_hi(kind6hi)
save/c_mbe_index_h4bmkp/
!                                                 ! beta_J for fSA
common/c_vcorrec_pars_h4bmkp/beta_j
save/c_vcorrec_pars_h4bmkp/
!                                                 ! London factors
parameter(widlondinv=1.d0/3.75d0)
!
common/c_lon_h4bmkp/eps_lond, sflo_lond, sfhi_lond, sfmid_lond, &
   dhsf_lond, dhsfinv_lond, rdelta_lond, iswtyp_lond
save/c_lon_h4bmkp/
!                                                  ! Vnsb switch
common/c_h3nsb_switch_h4bmkp/s0_nsb_1, st_nsb_1, sh_nsb_1, &
   s0_nsb_2, st_nsb_2, sh_nsb_2
save/c_h3nsb_switch_h4bmkp/
!                                                ! for sumVp
common/c_vpcalc_par_h4bmkp/ap_vpcalc_2, ap_vpcalc_3
save/c_vpcalc_par_h4bmkp/
!                              ! coefficient reordering indices
parameter(n_ord=27)
common/c_order_coef_h4bmkp/iord(n_ord)
save/c_order_coef_h4bmkp/
!                                             ! get r() from cc(),
common/c_rcalc_indices_h4bmkp/ipart1(6), ipart2(6)
save/c_rcalc_indices_h4bmkp/
!                                            ! and reverse indices
common/c_reverse_indices_h4bmkp/ir_ij(4, 4)
save/c_reverse_indices_h4bmkp/
!                                                      ! geom,
common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
   fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
   Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
save/c_geometry_h4bmkp/
!                                                ! and more geom
common/c_geo_sin_h4bmkp/geo_sin(3, 0:5)
save/c_geo_sin_h4bmkp/
!                                                     ! derivs
common/c_geomet_local_h4bmkp/drvdcc(6, 3, 12), dcomdcc(6, 3, 12)
save/c_geomet_local_h4bmkp/
!                               ! some default selector function indices
!
common/c_Achoice_h4bmkp/isqtet, igen1, icomp, inonlinh3, igen2
save/c_Achoice_h4bmkp/
!
common/c_chkh4geom_h4bmkp/ii2(6), ii3(6), ii4(6), ii5(6)
save/c_chkh4geom_h4bmkp/
!
common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
save/c_levelerr_h4bmkp/
!                                                     ! fnlclose
common/c_rcl_h4bmkp/rlo_cl, rhi_cl, rmid_cl, dinv_cl, oinv_cl
save/c_rcl_h4bmkp/
!                                                     ! for VH2
common/c_vh2_rlo_h4bmkp/rsw_h2, rlo_h2, aa_h2, bb_h2, cc_h2
save/c_vh2_rlo_h4bmkp/
!                                    ! for old interface
common/c_oldargs_h4bmkp/surfid
save/c_oldargs_h4bmkp/
!
!__At begining, have not yet initialized; have not read in a surface:
!
data init_h4bmkp/0/, iread_h4bmkp/0/
!
!__Linear coefficients for the H4 BMKP surface:
!
data(coef(i), i=1, 30)/ &
   11.41656861561020d0, 2.461981010111942d0, &
   -1.467350524873491d0, &
   0.d0, 0.4320536169168858d0, -.2425867014300068d0, &
   1.369331303380642d0, -1.433069016328931d0, &
   0.3878029623800905d0, &
   0.2243777140165565d0, 0.d0, 0.d0, &
   0.9565185801891209d0, 0.d0, 0.1679487106395962d0, &
   -.8043324172825887d-1, -.5553208620290703d-1, &
   0.4021797399144415d-1, &
   0.d0, 0.d0, 0.1407145922487997d0, &
   0.d0, -.5119002797896408d-1, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.2732447662063222d-1, 0.d0/
!
data(coef(i), i=31, 60)/ &
   0.2569544783161231d-2, 0.d0, -.3351870046177481d-1, &
   0.d0, 0.d0, -.2326424596451241d-2, &
   0.d0, 0.d0, 0.d0, &
   -.5881369223191887d-2, 0.d0, -.4674389804646140d-1, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.2235313229563889d-1, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   -.4422143271573974d-2, 0.2313604270707399d-2, 0.d0, &
   -.2588537141170540d-1, 0.1826575181931046d-1, 0.d0, &
   0.d0, -.1007716530727486d-1, -.3307300279251005d-2/
!
data(coef(i), i=61, 90)/ &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, -.4247747992767682d-2, 0.d0, &
   0.2044170535757541d-2, 0.d0, 0.d0, &
   0.1522380389286762d0, 0.d0, 0.2193168336511549d-1, &
   -.3864155496673826d0, 0.2952821198469801d-1, 0.d0, &
   -.5122925515081258d-1, -.7933584692639094d-1, 0.d0, &
   -.9640444763980278d-1, 0.d0, 0.d0, &
   0.d0, 0.d0, -.1766881903778219d-1, &
   -.7235812426737219d-1, 0.1185383058919827d0, &
   -.1391655300384850d-1/
!
data(coef(i), i=91, 120)/ &
   0.d0, 0.6580256079215258d-2, -.6442806985257051d-2, &
   0.d0, 0.d0, 0.1061238256124183d-1, &
   -.7080023799253410d-2, -.9177354648147840d-2, 0.d0, &
   -.1189054490030095d0, 0.8131113974522991d-1, &
   -.6028550495422950d-1, &
   -.1340099467289728d-1, 0.1188970731299674d-1, 0.d0, &
   0.1216426793454426d0, -.8154964357971532d-1, &
   0.6642703603111079d-1, &
   0.1318918383813553d-1, -.1392311809243365d-1, &
   0.2084072455893501d-2, &
   -.2821790556887770d-1, 0.1836263349587651d-1, &
   -.1461994658707290d-1, &
   -.2933600332178104d-2, 0.3279655300405622d-2, &
   -.7536040194164841d-3, &
   0.9662904089135244d-3, 0.d0, -.2353041625429331d-2/
!
data(coef(i), i=121, 150)/ &
   0.d0, 0.5611934951252549d-3, 0.9564226935654422d-3, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, -.6829845855150903d-3, &
   0.d0, 0.d0, 0.1903275078627081d-3, &
   -.3308243596276508d-4, 0.d0, 0.7794519273561473d-4, &
   0.d0, 0.d0, 0.2101122114040021d-2, &
   0.d0, -.4040602236056038d-3, -.2587008024119695d-3, &
   -.6364645397740569d-3, 0.d0, -.9353055974199075d-3, &
   0.d0, 0.d0, 0.4711266249527347d-3, &
   0.d0, 0.1569194222582758d-3, 0.d0/
!
data(coef(i), i=151, 180)/ &
   0.d0, 0.d0, -.7820219543883007d-4, &
   0.d0, 0.d0, 0.2074127775895557d-3, &
   0.d0, 0.d0, 0.3488603138430314d-3, &
   0.d0, 0.d0, -.4895717143534099d-3, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   -.1371570251840013d-4, 0.7871566414221949d-4, &
   -.8492275975618101d-4, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, -.1708344657789834d-2/
!
data(coef(i), i=181, 210)/ &
   0.d0, 0.3323527691292804d-3, -.7333214741878197d-4, &
   0.d0, 0.6460682421301904d-3, 0.d0, &
   -.2277390762352338d-3, 0.d0, 0.9128512806546217d-3, &
   -.6958228400160269d-3, 0.2660140098402603d-3, 0.d0, &
   0.5383770268573472d-3, -.5860104570730753d-3, 0.d0, &
   0.d0, 0.d0, -.7050220161376720d-4, &
   0.d0, -.1849834003028463d0, 0.d0, &
   0.4620812091093298d0, -.1757213672839463d0, 0.d0, &
   0.d0, -.4057065045231522d0, 0.5750547530386685d0, &
   -.8603309454113299d-1, 0.d0, 0.d0/
!
data(coef(i), i=211, 240)/ &
   0.2344020773325997d0, 0.d0, -.8795349352684374d-1, &
   0.1264880308849537d-1, 0.d0, 0.2668482770587186d-2, &
   0.d0, 0.d0, -.1100610087674357d0, &
   0.d0, 0.4198142714188278d-1, -.5309915669636895d-1, &
   0.d0, 0.8721717815328730d-2, 0.4286969491463773d-1, &
   0.d0, -.8749194176570741d-2, 0.2024906272110951d-1, &
   0.9083193758180202d-2, -.7903699985451101d-2, 0.d0, &
   0.d0, 0.d0, -.2249087223226226d-2, &
   0.d0, -.7488035985516799d-1, -.7884539373833625d-1, &
   0.4784488215937220d-1, 0.d0, 0.7757934134790705d-2/
!
data(coef(i), i=241, 270)/ &
   0.5405421234966801d-1, 0.8802835957276156d-1, 0.d0, &
   -.4030577591563147d-1, 0.d0, 0.d0, &
   -.2984565152037004d-1, -.1727219033185700d-1, &
   0.1182144991915245d-1, &
   0.7938201173244584d-2, -.1949620877956762d-2, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, -.3235243307954152d-3, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0/
!
data(coef(i), i=271, 300)/ &
   -.7620268765924905d-1, 0.2615432676679663d0, 0.d0, &
   -.3733918059753137d0, 0.d0, 0.d0, &
   0.2108908100274858d0, -.3908497880553245d-1, &
   -.4135382830384930d-1, &
   0.1522872914884577d0, -.2855267631463690d0, &
   0.1187218798350170d-1, &
   0.d0, 0.1720929635828025d0, 0.d0, &
   -.4762896947490703d-1, 0.d0, 0.1296479727681136d-1, &
   0.d0, 0.3875683711725621d-1, 0.d0, &
   0.d0, -.2728020475068825d-1, -.3568966505779430d-2, &
   0.2391400007431497d-2, 0.d0, 0.d0, &
   0.7092715349016440d-1, -.2171083138201483d-1, 0.d0/
!
data(coef(i), i=301, 330)/ &
   0.d0, 0.3703995087126884d-2, -.2434327063040666d-2, &
   -.9427841028470206d-1, 0.d0, 0.d0, &
   0.1294971819195222d-2, 0.9682682112195861d-3, 0.d0, &
   0.1631275239132733d-1, 0.3636984154327500d-2, &
   -.1340984765576260d-2, &
   0.d0, -.5956593779898524d-3, 0.2314568357119641d-3, &
   0.d0, 0.d0, -.1992197266813824d-2, &
   0.d0, 0.d0, -.1694357150309008d-3, &
   0.d0, 0.7029485453338565d-3, 0.d0, &
   0.d0, 0.d0, 0.1070975861800831d-3, &
   -.2113079739386820d-3, 0.d0, 0.d0/
!
data(coef(i), i=331, 360)/ &
   -.2539283910284728d-4, 0.d0, 0.d0, &
   0.d0, 0.1099815327496916d-2, -.1765748911285442d-2, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.1998984184634926d-4, 0.d0, 0.d0, &
   0.d0, -.6705593097363792d-4, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   -.7370107278577147d-4, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0/
!
data(coef(i), i=361, 390)/ &
   0.d0, 0.d0, -.1226471127884504d-3, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.3602375326308723d-4, &
   -.3519096757861601d-2, -.1411948724683564d-2, 0.d0, &
   -.1573773984487967d-1, 0.1185480315550147d-1, 0.d0, &
   0.5103394787098075d-2, -.3614550200861926d-2, &
   -.1770324228084872d-3, &
   0.d0, -.1992171779083550d-2, 0.d0, &
   0.1270692393828489d-1, -.4921205014438497d-2, &
   0.2064885831382025d-3, &
   -.2060673869866738d-2, 0.1474585053339496d-2, 0.d0, &
   0.d0, 0.9978301872727278d-3, -.3964660536460676d-4/
!
data(coef(i), i=391, 420)/ &
   -.1773396097672859d-2, 0.d0, 0.d0, &
   0.1172919080826781d-3, -.9737749811915252d-4, 0.d0, &
   -45.26288823530060d0, -8.088211185381054d0, &
   3.668642918880288d0, &
   0.4786804751301441d0, -1.707321870517821d0, &
   1.088841015252361d0, &
   -3.303208576830779d0, 3.936877551205625d0, 0.d0, &
   0.4179110229268445d0, 0.d0, -.3959245897503227d0, &
   3.770841516697761d0, 0.d0, -1.114631944449590d0, &
   -.8088285725825271d-1, 0.1566691026629507d0, &
   -.1855320817030686d-1, &
   0.d0, -.6410649862995392d-1, -.3266963544501503d0, &
   0.4155996356413736d-1, 0.d0, 0.d0/
!
data(coef(i), i=421, 450)/ &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, -.5501384523353580d-2, 0.1060007442378859d0, &
   -.1608713465509136d-1, 0.1437949646266004d-1, 0.d0, &
   0.d0, 0.3394636523584611d-1, -.1134732985014914d0, &
   0.d0, 0.d0, 0.1277036538813291d0, &
   0.d0, 0.d0, 0.d0, &
   0.8296264006392355d-2, -.3992780588062202d-1, &
   -.2201005099651295d-1, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.1166889972496094d-1, -.5709122448787395d-2/
!
data(coef(i), i=451, 480)/ &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.5276267894590526d-1, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.5828792976399929d-2, &
   0.d0, 0.d0, 0.d0, &
   -.6901495268242744d0, 0.d0, 0.d0, &
   1.463682894109018d0, 0.d0, 0.d0, &
   0.6284874210232054d0, 0.7177356034107231d-1, &
   0.1114026925935851d0, &
   0.7639309350060595d0, 0.d0, 0.d0/
!
data(coef(i), i=481, 510)/ &
   0.d0, 0.d0, 0.d0, &
   -.3105415152600133d0, -.1846415240159753d0, &
   -.1857352512650449d-1, &
   0.d0, 0.d0, 0.d0, &
   -.1474876842579986d0, -.3788771252443011d-1, &
   0.1371285727211139d-1, &
   0.d0, 0.6625356371830053d-1, -.1026112057095907d-1, &
   0.d0, -.1105251492579490d-5, 0.d0, &
   0.d0, -.6418485229163887d-6, 0.1132970198243758d-5, &
   0.3088336285550473d-5, 0.6006136636933950d-6, &
   0.3965148926799415d-6, &
   0.d0, 0.d0, -.4447322393731460d-6, &
   -.1699716222759334d-5, 0.d0, 0.d0/
!
data(coef(i), i=511, 540)/ &
   0.4528105444195479d-7, 0.d0, 0.4873262657912689d-7, &
   0.d0, -.1069438713799491d-5, 0.1111275907252112d-5, &
   0.d0, 0.6396482750065100d-7, 0.d0, &
   0.2620093358105838d-5, 0.d0, -.6800012598833259d-6, &
   0.d0, 0.d0, 0.d0, &
   -.1859877217865141d-5, 0.5267767186786455d-6, 0.d0, &
   -.4490791393108958d-7, 0.d0, 0.d0, &
   -.1808544844425166d-5, 0.d0, 0.1503227604390245d-5, &
   -.8450361230700443d-7, 0.d0, 0.d0, &
   0.1594595406051139d-5, 0.d0, -.1658370074628046d-5/
!
data(coef(i), i=541, 570)/ &
   0.d0, 0.1504941865060415d-6, -.1372299792933966d-7, &
   -.8036673503006943d-6, 0.4019859576399040d-6, 0.d0, &
   -.4451491359093753d-7, 0.d0, 0.d0, &
   0.4390555540312022d-5, -.1988249891489124d-5, 0.d0, &
   -.2560677863118263d-6, 0.9815971961161848d-6, &
   -.8183070027992053d-6, &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.1870275113972953d-6, &
   0.d0, 0.d0, 0.d0, &
   0.1899881667341878d-7, -.4642199767652417d-7, 0.d0, &
   -.3818998927536044d-5, 0.d0, 0.d0/
!
data(coef(i), i=571, 600)/ &
   0.4643005856103144d-5, 0.d0, 0.d0, &
   -.6998860287758503d-6, 0.d0, 0.d0, &
   0.d0, 0.1159361237893339d-5, 0.d0, &
   -.3532821151131939d-5, 0.d0, -.1119561115279527d-6, &
   0.6730732445924999d-6, 0.d0, 0.d0, &
   0.1743906906514087d-5, -.1011848690997574d-5, &
   0.8010838221017618d-7, &
   0.1751392240339859d-6, 0.1887909194821102d-6, &
   0.1006090339987867d-7, &
   -.1179071341543131d-6, 0.d0, 0.d0, &
   -.2577636134332770d0, 0.d0, 0.d0, &
   0.d0, 0.2950066706001413d-2, -.8436249003701826d-3/
!
data(coef(i), i=601, 630)/ &
   0.9650392433268715d-1, 0.1180865411962495d-1, &
   -.1890768945728852d-1, &
   0.d0, 0.d0, 0.1396056667617217d-3, &
   -.4784573309839698d-2, -.2307816824880044d-2, &
   0.2236248930560251d-2, &
   0.2741664444144839d-3, -.1050250271971434d-3, 0.d0, &
   0.d0, -.1319897018744042d-2, -.4284037762263289d-2, &
   0.9989080207036391d-3, 0.d0, -.1884447664581848d-3, &
   0.1495754666678066d-3, 0.d0, 0.4216827533402400d-2, &
   -.5873115173081413d-3, 0.d0, 0.d0, &
   0.d0, 0.d0, -.4041229946271942d-3, &
   0.6443324536316071d-4, -.8167457890020569d-5, &
   0.1056967699024861d-4/
!
data(coef(i), i=631, 660)/ &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.2803923282671573d-3, &
   0.1255823209120723d-3, 0.d0, 0.d0, &
   0.3631834621471845d-4, -.1237743740137913d-3, &
   -.7679069449277003d-4, &
   -.5544353440897096d-4, 0.d0, 0.d0, &
   0.d0, 0.2226610345487547d-4, 0.d0, &
   0.1102858209688644d-2, -.9707421030189208d-4, 0.d0, &
   0.d0, 0.d0, -.1788348599144874d-4, &
   -.8680038386835745d-4, 0.d0, 0.7502974676591654d-4, &
   0.d0, 0.d0, 0.d0/
!
data(coef(i), i=661, 690)/ &
   0.d0, 0.d0, 0.d0, &
   0.d0, 0.d0, -.1140689346615696d-5, &
   0.1338360967887129d-1, -.4734314575965839d-2, 0.d0, &
   0.1164975543347115d-1, 0.2368768109003181d-2, &
   0.2801552358660103d-3, &
   0.7028567446182806d-3, -.3124827789848925d-3, &
   -.3151243066290188d-3, &
   0.d0, 0.5887365420884373d-3, 0.d0, &
   -.8018026247461792d-2, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.8600408000548635d-4, &
   -.2211835739470488d-3, 0.d0, 0.3075426445130659d-4, &
   0.8401140913442456d-3, 0.d0, -.4243877152533952d-4/
!
data(coef(i), i=691, 720)/ &
   0.d0, 0.d0, 0.d0, &
   0.d0, 4.570473204155594d0, -.2798626320494389d0, &
   -6.940100755011418d0, 0.d0, 64.98997989641670d0, &
   1.715566496848222d0, -7.759401344839803d0, &
   -7.307106314539124d0, &
   -11.18609692968977d0, -92.20698555542199d0, &
   29.72626438699660d0, &
   11.98208769527558d0, -39.24233472335177d0, &
   137.9163648926951d0, &
   65.69977054586563d0, 0.d0, 214.0380086983486d0, &
   55.46915239130455d0, -166.2744780226556d0, 0.d0, &
   40.81904694776546d0, 0.d0, -62.76548193762972d0, &
   -37.71415377213584d0, -167.4967713102476d0, &
   -345.0583498103054d0/
!
data(coef(i), i=721, 750)/ &
   -192.0902057087892d0, -79.86702739644112d0, 0.d0, &
   0.d0, 413.0483295509064d0, 0.d0, &
   39.64091918839881d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   96.66100850781979d0, -169.3760884098744d0, 0.d0, &
   573.3282420428860d0, -3.671947638444790d0, &
   25.53422541109841d0, &
   -19.08478197765944d0, 14.20835333585394d0, 0.d0, &
   283.3923273325448d0, -217.6861835311140d0, &
   -62.28512734587679d0, &
   -62.71380848477281d0, -21.74157883889574d0, &
   86.04408325527467d0, &
   -323.5339117089853d0, 260.0626158828048d0, &
   -103.0644269330105d0/
!
data(coef(i), i=751, 780)/ &
   -271.2087580376460d0, 0.d0, 0.d0, &
   0.d0, 0.d0, 0.d0, &
   69.07644925729612d0, -52.23384556701277d0, 0.d0, &
   228.6493095830203d0, 124.2085335019347d0, &
   -85.21025019384091d0, &
   496.0374997733861d0, -684.4131378969931d0, &
   415.4879765735129d0, &
   890.6952128611791d0, -387.9510910233527d0, &
   230.2546692076336d0, &
   -101.1696229503429d0, 0.d0, -118.3505341147770d0, &
   0.d0, -734.6888505805208d0, 0.d0, &
   -387.8540736616000d0, 0.d0, 0.d0, &
   -249.9596065837551d0, 0.d0, 0.d0/
!
data(coef(i), i=781, n_p)/ &
   -99.09862868012513d0, 0.d0, 0.d0, &
   0.d0, -155.3886754928112d0, 563.8830985327528d0, &
   -94.35505172463631d0, 189.7518257337113d0/
!
!__Some non-linear parameters for the H4 BMKP surface:
!                                          ! rhoshifts must be 0.0
data vsca/nterms*0.d0/
!                                          ! selector betas
data beta/5*0.006d0, 0.003d0, 0.006d0, 0.0d0/
!                                          ! rho-powers must be 3.
data pow/nterms*3.d0/
!
!__Linear-term type-flags for the H4 BMKP surface:
!
data iv/2, 6*3, 100/
!
!__Lower and upper limits of linear-term power-summation ranges:
!
data imin/nterms*0/
!
data imax/7*2, 3/
!
data kmin/nterms*0/
!
data kmax/7*2, 8/
!
!__Cut-off function type-flags must all be unity:
!
data irho/nterms*1/
!
!__There are 788 linear parameters, and 7 Jacobi-type terms (that use powers of
!  r_a, r_b, R, and "curly-Y" functions of the angles theta1, theta2, phi):
!
data mlin/n_p/, njac/7/
!
!__MBE beta_p and indices for mbe-term; these indices were generated using the
!  subroutine "setmbe_h4bmkp()", then slightly modified for added efficiency:
!
data beta_p/0.8d0/
!
data mbe_perm/ &
   1, 2, 3, 4, 5, 6, 1, 3, 2, 5, 4, 6, 2, 1, 3, 4, 6, 5, 2, 3, 1, 6, 4, 5, &
   3, 1, 2, 5, 6, 4, 3, 2, 1, 6, 5, 4, 1, 4, 5, 2, 3, 6, 1, 5, 4, 3, 2, 6, &
   4, 1, 5, 2, 6, 3, 4, 5, 1, 6, 2, 3, 5, 1, 4, 3, 6, 2, 5, 4, 1, 6, 3, 2, &
   2, 4, 6, 1, 3, 5, 2, 6, 4, 3, 1, 5, 4, 2, 6, 1, 5, 3, 4, 6, 2, 5, 1, 3, &
   6, 2, 4, 3, 5, 1, 6, 4, 2, 5, 3, 1, 3, 5, 6, 1, 2, 4, 3, 6, 5, 2, 1, 4, &
   5, 3, 6, 1, 4, 2, 5, 6, 3, 4, 1, 2, 6, 3, 5, 2, 4, 1, 6, 5, 3, 4, 2, 1/
!
data nterms_use/95/, indmin/0/, indmax/3/, maxsum/8/
!
data k1_lo/1/, k1_hi/3/
!
data k2_lo/max_hi*1/
!
data(k2_hi(i), i=1, 12)/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/
!
data k3_lo/kind3hi*0/
!
data(k3_hi(i), i=1, 6)/1, 1, 2, 1, 2, 2/
!
data k4_lo/kind4hi*0/
!
data(k4_hi(i), i=1, 15)/0, 1, 1, 1, 1, 2, 2, 1, 1, 2, &
   2, 1, 1, 1, 0/
!
data(k5_lo(i), i=1, 32)/1, 0, 0, 0, 0, 0, 0, 1, 1, 0, &
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, &
   0, 0/
!
data(k5_hi(i), i=1, 32)/1, 0, 1, 1, 0, 0, 1, 2, 1, 2, &
   2, 1, 0, 1, 0, 1, 0, 0, 1, 2, 1, 0, 2, 1, 0, 0, 0, 2, 1, 1, &
   0, 0/
!
data(k6_lo(i), i=1, 52)/0, 0, 0, 0, 1, 0, 1, 0, 0, 0, &
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, &
   0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
   0, 0/
!
data(k6_hi(i), i=1, 52)/1, 0, 0, 1, 2, 2, 1, 2, 2, 2, &
   1, 2, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 3, 2, 3, 2, &
   1, 3, 2, 1, 2, 1, 1, 2, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, &
   0, 0/
!
!__The non-linear parameter beta_J for the goodness-of-choice factors fSA:
!
data beta_j/10.d0/
!
!__The parameters for the cusp-rounding part of the London equation:
!
data eps_lond/0.775d0/, sflo_lond/0.d0/, &
   sfhi_lond/7.5d0/, sfmid_lond/3.75d0/, &
   dhsf_lond/3.75d0/, dhsfinv_lond/widlondinv/, &
   rdelta_lond/3.d0/, iswtyp_lond/4/
!
!__The parameters for the weight/switch among H3 choices for Vnsb:
!
data s0_nsb_1/10.04d0/, st_nsb_1/0.312d0/, &
   sh_nsb_1/1.61d0/, s0_nsb_2/4.03d0/, &
   st_nsb_2/0.53d0/, sh_nsb_2/1.079d0/
!
!__The parameters for the pairwise non-bonding potential for sumVp (for Vcorr):
!
data ap_vpcalc_2/0.2593d0/, ap_vpcalc_3/1.6813d0/
!
!__The coefficients required for each q1q2mu subset, in the order in which
!  the actual terms are computed (takes care of symmetry for off-diagonal):
!
data iord/1, 2, 4, 2, 3, 5, 4, 5, 6, &
   7, 8, 10, 8, 9, 11, 10, 11, 12, &
   13, 14, 16, 14, 15, 17, 16, 17, 18/
!
!__The indices of the Cartesian coordinate array cc(4,3) used to get r(i):
!
data ipart1/1, 1, 1, 2, 2, 3/
data ipart2/2, 3, 4, 3, 4, 4/
!
!__The ir_ij(i,j) = the index k of the distance r(k) between atoms i and j
!
data ir_ij/0, 1, 2, 3, 1, 0, 4, 5, 2, 4, 0, 6, 3, 5, 6, 0/
!
!__Initialize derivatives to zero, to take care of those that are zero by
!  definition and thus are not explicitly set:
!
data geo/18*0.d0/, dridcc/72*0.d0/, dgeodcc/216*0.d0/, &
   fSA/3*0.d0/, sumVp/3*0.d0/, dfSAdcc/36*0.d0/, &
   dfSAdri/18*0.d0/, dsumVpdri/18*0.d0/, Yy/18*0.d0/, &
   dYydgeo/108*0.d0/, dgeodri/108*0.d0/, &
   jco_ri/1, 2, 3, 4, 5, 6, 2, 1, 3, 4, 6, 5, 3, 1, 2, 5, 6, 4/, &
   jco_geo/3, 1, 2, 4, 5, 6/
!
data geo_sin/18*0.d0/
!
!__Make sure that derivatives drvdcc are zero where they should be:
!
data drvdcc/216*0.d0/, dcomdcc/216*0.d0/
!
!__These indices result in the default set of selector functions A():
!
data isqtet/3/, igen1/5/, icomp/0/, inonlinh3/6/, &
   igen2/8/
!
!__These arrays give distance indices for use in getccofr_h4bmkp:
!
data ii2/2, 1, 1, 1, 1, 2/
data ii3/3, 3, 2, 5, 4, 4/
data ii4/4, 4, 5, 2, 3, 3/
data ii5/5, 6, 6, 6, 6, 5/
!
!__Data for error-checking:
!
data err_return/999999.d0/, rmin_err/0.001d0/, &
   level_err/1/
!
!  This is for fnlclose (small-distance London-switchover factor for
!  level_err = 1); all but rlo_cl and rhi_cl will be recomputed:
!
data rlo_cl/0.4d0/, rhi_cl/0.5d0/, rmid_cl/0.45d0/, &
   dinv_cl/20.d0/, oinv_cl/10.d0/
!
!  This is for small-distance VH2; all but rsw_h2 will be recomputed:
!
data rsw_h2/0.1193892052031558d0/, &
   rlo_h2/0.04515649888961240d0/, &
   aa_h2/0.42018214292507589d0/, &
   bb_h2/0.052365088838362515d0/, &
   cc_h2/-0.0011823120379919859d0/
!
!__This is used for the old interface to the H4 surface:
!
data surfid/4004369.26d0/
!
end
!
!**********************************************************************
!
subroutine parinit_h4bmkp
!----------------------------------------------------------------------
!  Perform initializations.
!----------------------------------------------------------------------
!
!  History of Aselect for each of the 8 terms in Vlinear:
!    (note that old iprnt(18) = 4 unless a different surface is read in):
!
!    iprnt(18) =  [OldName]   1          2          3          4: BMKP H4
!                  ------     ------     ------     ------     ----------
!        iV(1) --> Acorr  --> Acorr  --> Acorr  --> Acorr  --> Acorr
!        iV(2) --> Atet   --> Akite  --> Akite  --> Akite  --> Akite
!        iV(3) --> Asq    --> Atetsq --> Atetsq --> AnlH3H --> Atetsq
!        iV(4) --> AH3H   --> AH3H   --> AH3H   --> AH3H   --> AH3H
!        iV(5) --> Agen=1 --> Agen=1 --> AnlH3H --> Agen=1 --> Agen=1.0
!        iV(6) --> Acomp  --> Acomp  --> Acomp  --> Acomp  --> AnonlinH3H
!        iV(7) --> AH2H2  --> AH2H2  --> AH2H2  --> AH2H2  --> AH2H2
!        iV(8) --> Aeight --> AnlH3H --> 1.0    --> 1.0    --> 1.0
!
!  Allowed values and meanings for iV(i):  ( iV(i) = 0 means term not used )
!
!        iV(1) > 0 (=2): use Vcorr: BMKP H4
!        iV(2) > 0 (=3): use Vtet --> Vkite: BMKP H4
!        iV(3) = 1: 5-term Asq: NOT IMPLEMENTED in h4bmkp.f;
!              = 2: 10-term Vsq: NOT IMPLEMENTED in h4bmkp.f;
!              > 2 (=3): use Vsq --> Vtetsq: BMKP H4
!        iV(4) > 0 (=3): use Vh3h: BMKP H4
!        iV(5) > 0 (=3): use Vgen: BMKP H4
!        iV(6) = 94: use V_j94: NOT IMPLEMENTED in h4bmkp.f;
!              > 0 and not 94 (=3): use Vcomp --> VnonlinH3H: BMKP H4
!        iV(7) > 0 (=3): use Vh2h2: BMKP H4
!        iV(8) = 1: one-parameter smoothed step: NOT IMPLEMENTED in h4bmkp.f;
!              = 2 or 3: use Veight --> Vgen;
!              = 94: use V_j94: NOT IMPLEMENTED in h4bmkp.f;
!              > 3 and not 94: use Aguado-type many-body expansion: BMKP H4:
!                              (note that only iprnt(17) = 0 is allowed now);
!           --- some MBE iV(8) values have following effects:
!        iV(8) = 865: indmin = 0, indmax = 12, maxsum = 12 (865 terms)
!        iV(8) = 124: indmin = 0, indmax = 8,  maxsum = 12 (851 terms)
!        iV(8) = 95 : indmin = 0, indmax = 4,  maxsum = 13 (542 terms)
!        iV(8) = 100: indmin = 0, indmax = 3,  maxsum = 8  (95 terms): BMKP H4
!        iV(8) = 140: indmin = 0, indmax = 3,  maxsum = 9  (130 terms)
!        iV(8) = 200: indmin = 0, indmax = 4,  maxsum = 9  (190 terms)
!           --- otherwise: indmin = imin(8), indmax = imax(8), maxsum = kmax(8)
!           --- NOTE that for the specific cases above, the values of imin(8),
!                    imax(8), and kmax(8) are ignored
!
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   parameter(pi=3.141592653589793d0, &
             z3=3.d0, z4=4.d0, z5=5.d0, z24=24.d0, z96=96.d0, &
             zfm0=z5/(z4*pi), zfm1=z5/(z24*pi), zfm2=z5/(z96*pi), &
             zfl0m0=1.d0/(4.d0*pi), z2=2.d0)
!
   common/c_curlyY_par_h4bmkp/z4pi, z4pi_rt2, rt_zfl0m0, &
      rt_zfm0, rt_zfm1_3, rt_zfm2_3
   save/c_curlyY_par_h4bmkp/
!
   parameter(nterms=8, maxl=7*224 + 1350)
   common/c_coef_h4bmkp/coef(maxl), vsca(nterms), beta(nterms), &
      pow(nterms), iv(nterms), imin(nterms), imax(nterms), &
      kmin(nterms), kmax(nterms), irho(nterms), mlin, njac, &
      ikmin_lo, ikmax_hi
   save/c_coef_h4bmkp/
!
   common/c_lon_h4bmkp/eps_lond, sflo_lond, sfhi_lond, sfmid_lond, &
      dhsf_lond, dhsfinv_lond, rdelta_lond, iswtyp_lond
   save/c_lon_h4bmkp/
!
   parameter(max_hi=12, &
             kind3hi=max_hi*max_hi, kind4hi=kind3hi*max_hi, &
             kind5hi=3000, kind6hi=3000)
   common/c_mbe_index_h4bmkp/beta_p, mbe_perm(6, 24), nterms_use, &
      indmin, indmax, maxsum, k1_lo, k1_hi, k2_lo(max_hi), &
      k2_hi(max_hi), k3_lo(kind3hi), k3_hi(kind3hi), &
      k4_lo(kind4hi), k4_hi(kind4hi), k5_lo(kind5hi), &
      k5_hi(kind5hi), k6_lo(kind6hi), k6_hi(kind6hi)
   save/c_mbe_index_h4bmkp/
!
   parameter(n_ord=27)
   common/c_order_coef_h4bmkp/iord(n_ord)
   save/c_order_coef_h4bmkp/
!
   parameter(n_ord_max=64, nc_ord=n_ord_max*4*nterms)
   common/c_ordcoef_h4bmkp/ordcoef(nc_ord), ncol(nterms), &
      nLL(nterms), nLL4(nterms), nkeep(nterms), ieqp(nterms)
   save/c_ordcoef_h4bmkp/
!
   common/c_Achoice_h4bmkp/isqtet, igen1, icomp, inonlinh3, igen2
   save/c_Achoice_h4bmkp/
!
   common/c_rcl_h4bmkp/rlo_cl, rhi_cl, rmid_cl, dinv_cl, oinv_cl
   save/c_rcl_h4bmkp/
!
   common/c_vh2_rlo_h4bmkp/rsw_h2, rlo_h2, aa_h2, bb_h2, cc_h2
   save/c_vh2_rlo_h4bmkp/
!
   common/c_init_h4bmkp/init_h4bmkp, iread_h4bmkp
   save/c_init_h4bmkp/
!
   parameter(nqqm_orig=18)
!
   dimension E_h2(3)
!
!__Get small-distance VH2 parameters (these should be essentially
!  the same as the values stored in the block data; the value of
!  rsw_h2 stored there is such that the second derivative will
!  remain continuous at rsw_h2, although not at rlo_h2):
!
   rtmp = rsw_h2
   rsw_h2 = rlo_h2 + (rtmp - rlo_h2)*0.01d0
   call vH2opt95_h4bmkp(rtmp, E_h2, 1)
   rsw_h2 = rtmp
!
   f_h2 = (-1.d0/rsw_h2**2 - E_h2(2)) &
          /(E_h2(1) - 1.d0/rsw_h2)
!
   aa_h2 = 3.d0/rsw_h2**4 - f_h2/rsw_h2**3
   bb_h2 = 2.d0*f_h2/rsw_h2**2 - 4.d0/rsw_h2**3
   cc_h2 = 1.d0/rsw_h2**2 - f_h2/rsw_h2
!
   rlo_h2 = (-bb_h2 - sqrt(bb_h2**2 - 4.d0*aa_h2*cc_h2)) &
            /(2.d0*aa_h2)
!
   cc_h2 = (E_h2(1) - 1.d0/rsw_h2) &
           /(1.d0/rsw_h2**3 - 2.d0/(rlo_h2*rsw_h2**2) &
             + 1.d0/(rsw_h2*rlo_h2**2))
   bb_h2 = -2.d0*cc_h2/rlo_h2
   aa_h2 = 1.d0 + cc_h2/rlo_h2**2
!
!__Get small-distance London-switchover parameters (these should be
!  the same as the values stored in the block data):
!
   rmid_cl = (rlo_cl + rhi_cl)/2.d0
   dinv_cl = 2.d0/(rhi_cl - rlo_cl)
   oinv_cl = dinv_cl/2.d0
!
!__Get constants for use in the curlyY computations.
!
   z4pi = z4*pi
   z4pi_rt2 = z4pi*sqrt(z2)
   rt_zfl0m0 = sqrt(zfl0m0)
   rt_zfm0 = sqrt(zfm0)
   rt_zfm1_3 = sqrt(zfm1)*z3
   rt_zfm2_3 = sqrt(zfm2)*z3
!
!__Check whether there is an MBE-type term (get number of Jacobi-terms).
!
   if (iV(8) .le. 3) then
      njac = nterms
   else
      njac = nterms - 1
   end if
!
!__Check whether the variable selector function indices are needed
!
   if (isqtet .gt. 0) then
      if (iV(isqtet) .le. 0) isqtet = 0
   end if
   if (igen1 .gt. 0) then
      if (iV(igen1) .le. 0) igen1 = 0
   end if
   if (icomp .gt. 0) then
      if (iV(icomp) .le. 0) icomp = 0
   end if
   if (inonlinh3 .gt. 0) then
      if (iV(inonlinh3) .le. 0) inonlinh3 = 0
   end if
   if (igen2 .gt. 0) then
      if (iV(igen2) .le. 0) igen2 = 0
   end if
!
!__Get coefficients in the order in which the first 4 sets of q1q2mu terms are
!  computed (note that there are some duplications of the coefficient values):
!
   nsord = 0
   nsc = 0
!
   ikmax_hi = -9
   ikmin_lo = 9
   nshift = 0
   mshift = 0
   nprev = 0
   n_terms = 0
!                  ! for each of the Jacobi-type terms
   do n = 1, njac
!                          ! IF term does not exist, has no coefficients
      if (iV(n) .le. 0) then
!
         ncol(n) = 0
         nLL4(n) = 0
         ieqp(n) = 0
!                        ! OTHERWISE, if term does exist:
      else
!                  ! check whether index ranges are same as previous term
         ieqp(n) = 1
         if (nprev .gt. 0) then
            if (imin(n) .ne. imin(nprev) .or. &
                imax(n) .ne. imax(nprev) .or. &
                kmin(n) .ne. kmin(nprev) .or. &
                kmax(n) .ne. kmax(nprev)) ieqp(n) = 0
         end if
!                                              ! index extremes
         ikmin_lo = min(ikmin_lo, imin(n), kmin(n))
         if (ikmin_lo .lt. -8) stop &
            ' STOP -- Error: parinit_h4bmkp: imin,kmin < -8 '
         ikmax_hi = max(ikmax_hi, imax(n), kmax(n))
         if (ikmax_hi .gt. 9) stop &
            ' STOP -- Error: parinit_h4bmkp: imax,kmax > 9 '
!
!                              ! numbers of values in the summations
         nk = kmax(n) - kmin(n) + 1
         ni = imax(n) - imin(n) + 1
         nisq = ni*ni
         nLL(n) = nk*nisq
         nLL4(n) = 4*nLL(n)
         nqqm = (nk*ni*(ni + 1))/2
         nkeep(n) = nqqm*4
         ncol(n) = nkeep(n) + nLL(n)
         n_terms = n_terms + ncol(n)
         if (n_terms .gt. maxl) stop ' STOP -- nterms > maxl '
         ikeep = nshift
         Lshift = 0
!                        ! get version of coefficients in order needed
         do L = 1, 4
            kshift = Lshift
            do k = 1, nk
            do i = 1, ni
            do j = 1, i
               ij = (i - 1)*ni + j + kshift
               ikeep = ikeep + 1
               ordcoef(ij + mshift) = coef(ikeep)
               if (i .ne. j) then
                  ji = (j - 1)*ni + i + kshift
                  ordcoef(ji + mshift) = coef(ikeep)
               end if
            end do
            end do
            kshift = kshift + nisq
            end do
            Lshift = Lshift + nLL(n)
         end do
         mshift = mshift + nLL4(n)
         nshift = nshift + ncol(n)
!                                  ! and check them, for BMKP H4 case
         if (iread_h4bmkp .eq. 0) then
         do k = 1, 4
         do j = 1, n_ord
         if (ordcoef(j + nsord) .ne. &
             coef(iord(j) + nsc)) then
            write (6, *) 'n,nsord,nsc,k,j,iord(j)', &
               n, nsord, nsc, k, j, iord(j)
            stop ' STOP -- Error: ordcoef mismatch. '
         end if
         end do
         nsord = nsord + n_ord
         nsc = nsc + nqqm_orig
         end do
         nsc = nsc + n_ord
         end if
!                  ! this is now the latest existing Jacobi-type term
         nprev = n
!
      end if
!
   end do
!
!__Check that the long-distance cut-off term has the right form:
!
   do i = 1, nterms
      if (vsca(i) .ne. 0.d0) stop &
         ' STOP -- Error: parinit_h4bmkp: rhoshift non-zero. '
      if (abs(pow(i) - 3.d0) .gt. 1.d-6) stop &
         ' STOP -- Error: parinit_h4bmkp: pow (of rho) is not 3. '
   end do
!
!__Check for non-allowed Jacobi term types.
!
   if (iV(3) .eq. 1 .or. iV(3) .eq. 2) stop &
      ' STOP -- Error: parinit_h4bmkp: bad iV(3) = 1 or 2 '
   if (iV(6) .eq. 94) stop &
      ' STOP -- Error: parinit_h4bmkp: bad iV(6) = 94 '
   if (iV(8) .eq. 1 .or. iV(8) .eq. 94) stop &
      ' STOP -- Error: parinit_h4bmkp: bad iV(8) = 1 or 94 '
!
!__Check for non-allowed London cusp-rounding types.
!
   if (iswtyp_lond .ne. 0 .and. iswtyp_lond .ne. 4) stop &
      ' STOP -- Error: parinit_h4bmkp: iswtyp_lond not 0 or 4 '
!
!__Set the midpoint, half-width, and inverse half-width of London cusp-rounding
!
   if (iswtyp_lond .gt. 0) then
      if (sflo_lond .ge. sfhi_lond - 0.1d0) stop &
         ' STOP -- Error: parinit_h4bmkp: Lond:sfhi - sflo < 0.1 '
      sfmid_lond = (sflo_lond + sfhi_lond)*0.5d0
      dhsf_lond = (sfhi_lond - sflo_lond)*0.5d0
      dhsfinv_lond = 1.d0/dhsf_lond
   end if
!
!__If a surface has been read in (non-default, i.e., not BMKP H4 surface),
!  and it has an MBE-type term, then get the necessary MBE-index limits:
!
   if (iread_h4bmkp .ne. 0 .and. njac .lt. nterms) then
!                                                ! main limits
      indmin = 0
      if (iV(8) .eq. 95) then
         indmax = 4
         maxsum = 13
      else if (iV(8) .eq. 124) then
         indmax = 8
         maxsum = 12
      else if (iV(8) .eq. 865) then
         indmax = 12
         maxsum = 12
      else if (iV(8) .eq. 100) then
         indmax = 3
         maxsum = 8
      else if (iV(8) .eq. 140) then
         indmax = 3
         maxsum = 9
      else if (iV(8) .eq. 200) then
         indmax = 4
         maxsum = 9
      else
         indmin = imin(nterms)
         indmax = imax(nterms)
         maxsum = kmax(nterms)
      end if
!                                    ! the rest of the limits
      call setmbe_h4bmkp(nterms_mbe)
!                                          ! too many coeffs?
      if (n_terms + nterms_mbe .gt. maxl) stop &
         ' STOP -- Error: parinit_h4bmkp: nterms > maxl '
!
   end if
!
!__The initialization has been performed; should only be repeated if a surface
!  is read in.
!
   init_h4bmkp = 1
   iread_h4bmkp = max(iread_h4bmkp, 0)
!
   return
end
!
!**********************************************************************
!
subroutine setmbe_h4bmkp(nterms_mbe)
!----------------------------------------------------------------------
!  Get index-limits for summations required for MBE-type term.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   parameter(max_hi=12, maxsum_hi=13, &
             kind3hi=max_hi*max_hi, kind4hi=kind3hi*max_hi, &
             kind5hi=3000, kind6hi=3000)
   common/c_mbe_index_h4bmkp/beta_p, mbe_perm(6, 24), nterms_use, &
      indmin, indmax, maxsum, k1_lo, k1_hi, k2_lo(max_hi), &
      k2_hi(max_hi), k3_lo(kind3hi), k3_hi(kind3hi), &
      k4_lo(kind4hi), k4_hi(kind4hi), k5_lo(kind5hi), &
      k5_hi(kind5hi), k6_lo(kind6hi), k6_hi(kind6hi)
   save/c_mbe_index_h4bmkp/
!
   dimension izero(0:max_hi), in0low(0:6)
!                                    ! check for bad index ranges
   if (indmin .lt. 0) stop &
      ' STOP -- Error: setmbe_h4bmkp: MBE: indmin=imin(8) < 0 '
   if (indmax .gt. max_hi) stop &
      ' STOP -- Error: setmbe_h4bmkp: MBE: indmax=imax(8) > 12 '
   if (maxsum .gt. maxsum_hi) stop &
      ' STOP -- Error: setmbe_h4bmkp: MBE: maxsum=kmax(8) > 13 '
   if (indmax .gt. maxsum) stop &
      ' STOP -- Error: setmbe_h4bmkp: MBE: indmax > maxsum '
!
!__Some initializations needed for finding MBE-ranges:
!
   izero(0) = 1
   do i = 1, max_hi
      izero(i) = 0
   end do
!
   n0max = 3
   isum_5 = 1
   isum_maxsum_5 = isum_5 - maxsum
!
   do i = 0, 6
   if (i .lt. n0max) then
      in0low(i) = max(indmin, 0)
   else
      in0low(i) = max(indmin, 1)
   end if
   end do
!
   k = 0
   nterms_use = 0
   kind3 = 0
   kind4 = 0
   kind5 = 0
   kind6 = 0
!
   k1_lo = in0low(6)
   k1_hi = min(indmax, maxsum + n0max - 5)
!
!__Go through the MBE-summations, to get index-limits for allowed terms in sum:
!
   do i1 = k1_lo, k1_hi
!
      n1 = izero(i1)
      k1 = maxsum - i1
      k2_lo(i1) = in0low(n1)
      k2_hi(i1) = min(i1, k1, indmax)
!
      do i2 = k2_lo(i1), k2_hi(i1)
!
         n2 = n1 + izero(i2)
         k2 = k1 - i2
         i3_lo = in0low(n2)
         if (i1 + i2 .eq. 0) i3_lo = in0low(6)
         kind3 = kind3 + 1
         if (kind3 .gt. kind3hi) stop &
            ' STOP -- Error: setmbe_h4bmkp: kind3 too big '
         k3_lo(kind3) = i3_lo
         k3_hi(kind3) = min(i2, k2, indmax)
!
         do i3 = i3_lo, k3_hi(kind3)
!
            n3 = n2 + izero(i3)
            k3 = k2 - i3
            kind4 = kind4 + 1
            if (kind4 .gt. kind4hi) stop &
               ' STOP -- Error: setmbe_h4bmkp: kind4 too big '
            k4_lo(kind4) = in0low(n3)
            k4_hi(kind4) = min(i2, k3, indmax)
!
            do i4 = k4_lo(kind4), k4_hi(kind4)
!
               n4 = n3 + izero(i4)
               k4 = k3 - i4
               i5_lo = in0low(n4)
               if (i1 + i4 .eq. 0) i5_lo = in0low(6)
               i5_lo = max(i5_lo, isum_maxsum_5 + k4)
               if (i4 .eq. i2) then
                  i5_hi = min(i3, k4, indmax)
               else
                  i5_hi = min(i2, k4, indmax)
               end if
               if (i3 .eq. i2) i5_hi = min(i5_hi, i4)
               if (i4 .gt. i3) i5_hi = min(i5_hi, i2 - 1)
               kind5 = kind5 + 1
               if (kind5 .gt. kind5hi) stop &
                  ' STOP -- Error: setmbe_h4bmkp: kind5 too big '
               k5_lo(kind5) = i5_lo
               k5_hi(kind5) = i5_hi
!
               do i5 = i5_lo, i5_hi
!
                  n5 = n4 + izero(i5)
                  k5 = k4 - i5
                  i6_lo = in0low(n5)
                  if (i2 + i4 .eq. 0 .or. i3 + i5 .eq. 0) &
                     i6_lo = max0(i6_lo, 1)
                  i6_lo = max(i6_lo, isum_maxsum_5 + k4, &
                              isum_maxsum_5 + k5 + i1, &
                              isum_maxsum_5 + k5 + i2, &
                              isum_maxsum_5 + k5 + i3, &
                              isum_maxsum_5 + k5 + i4)
                  if (i4 .le. i3) then
                     i6_hi = min(i1, k5, indmax)
                  else
                     i6_hi = min(i1 - 1, k5, indmax)
                  end if
                  if (i1 .eq. i4) then
                     i6_hi = min0(i6_hi, i3, i5)
                  else if (i1 .eq. i2) then
                     i6_hi = min0(i6_hi, i5)
                  end if
                  kind6 = kind6 + 1
                  if (kind6 .gt. kind6hi) stop &
                     'STOP -- Error: setmbe_h4bmkp: kind6 too big'
                  k6_lo(kind6) = i6_lo
                  k6_hi(kind6) = i6_hi
!
                  if (i6_hi .lt. i6_lo) then
!
                     if (i5 .eq. k5_lo(kind5)) then
                        k5_lo(kind5) = k5_lo(kind5) + 1
                        kind6 = kind6 - 1
                     else if (i5 .eq. k5_hi(kind5)) then
                        k5_hi(kind5) = k5_hi(kind5) - 1
                        kind6 = kind6 - 1
                     end if
!
                  else
!
                     nterms_use = nterms_use + i6_hi - i6_lo + 1
!
                  end if
!
               end do
!
               if (i5_hi .lt. i5_lo) then
               if (i4 .eq. k4_lo(kind4)) then
                  k4_lo(kind4) = k4_lo(kind4) + 1
                  kind5 = kind5 - 1
               else if (i4 .eq. k4_hi(kind4)) then
                  k4_hi(kind4) = k4_hi(kind4) - 1
                  kind5 = kind5 - 1
               end if
               end if
!
            end do
!
            if (k4_hi(kind4) .lt. k4_lo(kind4)) then
            if (i3 .eq. k3_lo(kind3)) then
               k3_lo(kind3) = k3_lo(kind3) + 1
               kind4 = kind4 - 1
            else if (i3 .eq. k3_hi(kind3)) then
               k3_hi(kind3) = k3_hi(kind3) - 1
               kind4 = kind4 - 1
            end if
            end if
!
         end do
!
         if (k3_hi(kind3) .lt. k3_lo(kind3)) then
         if (i2 .eq. k2_lo(i1)) then
            k2_lo(i1) = k2_lo(i1) + 1
            kind3 = kind3 - 1
         else if (i2 .eq. k2_hi(i1)) then
            k2_hi(i1) = k2_hi(i1) - 1
            kind3 = kind3 - 1
         end if
         end if
!
      end do
!
   end do
!
!__Return the total number of terms in the MBE-summation, i.e., the number of
!  linear coefficients required for the MBE term.
!
   nterms_mbe = nterms_use
!
   return
end
!
!**********************************************************************
!
subroutine jcoord_h4bmkp(RIJ, Q, DERR)
!----------------------------------------------------------------------
!  Schwenke/Keogh routine that computes values for dgeodri(ig,*,*).
!----------------------------------------------------------------------
!     in:  Rij(6)
!     out: Q(6),dErr(6,6)
!----------------------------------------------------------------------
! mar30/93 ... modified to be used in Keogh's H4 programme
!              this version assumes ntraj=1
!              riji,rij2 now calculated instead of passed in
!  Schwenke's routine
!  1.d-39 caused underflow problems ... changed to epsilo by W.J.Keogh
!  calculate the jacobi coordinates from the 6 internuclear distances
!  also calculate the derivatives
!  this module call no others
! RIJ(I) ARE THE INTERPAIR DISTANCES, Rij(1)=r12 Rij(2)=r13 Rij(3)=r14
!                                     Rij(4)=r23 Rij(5)=r24 Rij(6)=r34
! Q(I) ARE THE JACOBI COORDINATES, Q(1)=R Q(2)=r1 Q(3)=r2
!                                  Q(4)=cosT1 Q(5)=cosT2 Q(6)=cosPhi
! derr(I,J) IS THE DERIVATIVE OF THE I'TH JACOBI COORDINATE WRT THE
!           J'TH INTERPAIR DISTANCE.
! RIJI ARE THE INVERSES OF THE INTERPAIR DISTANCES
! RIJ2 ARE THE SQUARES OF THE INTERPAIR DISTANCES
!
!   CONVERT FROM INTERPAIR DISTANCES TO JACOBI COORDINATES.
!   MOLECULE 1 IS ATOMS 1 AND 2, AND THE MOLECULES ARE HOMONUCLEAR.
!   ALSO CALCULATE THE DERIVATIVES WRT THE INTERPAIR DISTANCES.
!
!   NOTE THAT FOR CERTAIN GEOMENTRIES (EG LINEAR), THIS ROUTINE WILL
!   FAIL AND DIVIDE BY ZERO.  (Division by zero is now avoided, but
!   erroneous values may still be returned for some planar and linear
!   geometries.)
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)
!
!-not;      logical pf
   PARAMETER(zero=0.d0, half=0.5d0, one=1.d0)
!
   dimension X(3, 4), XD(3, 4, 6), DER(6), DD(9, 6), &
      RIJ(6), Q(6), DERR(6, 6), RIJI(6), &
      RIJ2(6)
!
   parameter(epsilo=1.d-30)
!
   do i = 1, 6
      riji(i) = 1.d0/rij(i)
      rij2(i) = rij(i)*rij(i)
   end do
!-not;      pf = .false.
!-not;      if(pf)then
!-not;         write(6,*) ' entering subr.jacoord,  eps=',epsilo
!-not;         write(6,6000) rij
!-not;      end if
!-not; 6000 format('Rij = ',6(1x,g12.6))
   DO I = 1, 4
   DO J = 1, 3
      X(J, I) = zero
   end do
   end do
   DO I = 1, 6
   DO J = 1, 4
   DO K = 1, 3
      XD(K, J, I) = zero
   end do
   end do
   end do
   X(3, 2) = RIJ(1)
   XD(3, 2, 1) = one
   CC = -(RIJ2(4) - RIJ2(1) - RIJ2(2))*half &
        *(RIJI(1)*RIJI(2))
   DER(1) = (RIJI(2)) - CC*RIJI(1)
   DER(2) = (RIJI(1)) - CC*RIJI(2)
   DER(3) = zero
   DER(4) = -RIJ(4)*(RIJI(1)*RIJI(2))
   DER(5) = zero
   DER(6) = zero
   CCS = one - CC*CC
   X(3, 3) = RIJ(2)*CC
   DO I = 1, 6
      XD(3, 3, I) = RIJ(2)*DER(I)
   end do
!-not;      if(pf) write(6,*) 'jcoord: loop 1'
   XD(3, 3, 2) = XD(3, 3, 2) + CC
   CCS = SQRT(ABS(CCS))
   CC = RIJ(2)*CC/(CCS + epsilo)
   X(2, 3) = RIJ(2)*CCS
   DO I = 1, 6
      XD(2, 3, I) = -CC*DER(I)
   end do
!-not;      if(pf) write(6,*) 'jcoord: loop 2'
   XD(2, 3, 2) = XD(2, 3, 2) + CCS
   X(3, 4) = (RIJ2(3) + RIJ2(1) - RIJ2(5))*half*RIJI(1)
   XD(3, 4, 1) = one - X(3, 4)*RIJI(1)
   XD(3, 4, 3) = RIJ(3)*RIJI(1)
   XD(3, 4, 5) = -RIJ(5)*RIJI(1)
   DER(1) = one/(X(2, 3) + epsilo)
   X(2, 4) = (RIJ2(3) + X(2, 3)**2 + X(3, 3)**2 - 2E0*X(3, 4)* &
              X(3, 3) - RIJ2(6))*half*DER(1)
   DO I = 1, 6
      XD(2, 4, I) = (X(2, 3)*XD(2, 3, I) + X(3, 3)*(XD(3, 3, I) &
                                                    - XD(3, 4, I)) - X(3, 4)*XD(3, 3, I) - X(2, 4) &
                     *XD(2, 3, I))*DER(1)
   end do
!-not;      if(pf) write(6,*) 'jcoord: loop 3'
   XD(2, 4, 3) = XD(2, 4, 3) + RIJ(3)*DER(1)
   XD(2, 4, 6) = XD(2, 4, 6) - RIJ(6)*DER(1)
   X(1, 4) = (RIJ2(3) - X(2, 4)**2 - X(3, 4)**2)
   CC = one/(SQRT(ABS(X(1, 4))) + epsilo)
   DO I = 1, 6
      XD(1, 4, I) = -(X(2, 4)*XD(2, 4, I) + X(3, 4)*XD(3, 4, I))*CC
   end do
!-not;      if(pf) write(6,*) 'jcoord: loop 4 '
   XD(1, 4, 3) = XD(1, 4, 3) + RIJ(3)*CC
   X(1, 4) = one/(CC + epsilo)
   RX = half*(X(1, 3) + X(1, 4))
   RY = half*(X(2, 3) + X(2, 4))
   RZ = half*(X(3, 3) + X(3, 4) - RIJ(1))
   DO I = 1, 6
      DD(1, I) = half*(XD(1, 3, I) + XD(1, 4, I))
      DD(2, I) = half*(XD(2, 3, I) + XD(2, 4, I))
      DD(3, I) = half*(XD(3, 3, I) + XD(3, 4, I))
   end do
!-not;      if(pf) write(6,*) 'jcoord: Q(1) and derr(1,i) '
   DD(3, 1) = DD(3, 1) - half
   Q(1) = SQRT(ABS(RX*RX + RY*RY + RZ*RZ))
   CC = one/(Q(1) + epsilo)
   DO I = 1, 6
      DERR(1, I) = (DD(1, I)*RX + DD(2, I)*RY &
                    + DD(3, I)*RZ)*CC
   end do
!-not;      if(pf)         write(6,*) 'jcoord: Q(2) and derr(2,i) '
   R1Z = RIJ(1)
   Q(2) = RIJ(1)
   DO I = 1, 6
      DERR(2, I) = zero
   end do
!-not;      if(pf)         write(6,*) 'jcoord: loop 6'
   DERR(2, 1) = one
   R2X = X(1, 4) - X(1, 3)
   R2Y = X(2, 4) - X(2, 3)
   R2Z = X(3, 4) - X(3, 3)
   DO I = 1, 6
      DD(4, I) = XD(1, 4, I) - XD(1, 3, I)
      DD(5, I) = XD(2, 4, I) - XD(2, 3, I)
      DD(6, I) = XD(3, 4, I) - XD(3, 3, I)
   end do
!-not;      if(pf)         write(6,*) 'jcoord: Q(3) and derr(3,i) '
   Q(3) = RIJ(6)
   CCS = one/(Q(3) + epsilo)
   DO I = 1, 6
      DERR(3, I) = (DD(4, I)*R2X + DD(5, I)*R2Y &
                    + DD(6, I)*R2Z)*CCS
   end do
!-not;      if(pf)         write(6,*) 'jcoord: Q(4) and derr(4,i) '
   Q(4) = RZ*CC
   DO I = 1, 6
      DERR(4, I) = (DD(3, I) - RZ*DERR(1, I)*CC)*CC
   end do
!-not;      if(pf)         write(6,*) 'jcoord: Q(5) and derr(5,i) '
   CCS = Q(1)*Q(3)
   CCS = one/(CCS + epsilo)
   Q(5) = (R2X*RX + R2Y*RY + R2Z*RZ)*CCS
   DO I = 1, 6
      DERR(5, I) = (DD(4, I)*RX + DD(1, I)*R2X + DD(5, I)*RY &
                    + DD(2, I)*R2Y + DD(6, I)*RZ + DD(3, I)*R2Z &
                    - (Q(1)*DERR(3, I) + DERR(1, I)*Q(3))*Q(5)) &
                   *CCS
   end do
!-not;      if(pf)         write(6,*) 'jcoord: loop 9'
   CM1X = -R1Z*RY
   CM1Y = R1Z*RX
   CM1Z = zero
   CM2X = R2Y*RZ - R2Z*RY
   CM2Y = R2Z*RX - R2X*RZ
   CM2Z = R2X*RY - R2Y*RX
   DO I = 1, 6
      DD(7, I) = R2Y*DD(3, I) + DD(5, I)*RZ - R2Z*DD(2, I) &
                 - DD(6, I)*RY
      DD(8, I) = R2Z*DD(1, I) + DD(6, I)*RX - R2X*DD(3, I) &
                 - DD(4, I)*RZ
      DD(9, I) = R2X*DD(2, I) + DD(4, I)*RY - R2Y*DD(1, I) &
                 - DD(5, I)*RX
      DD(4, I) = -R1Z*DD(2, I)
      DD(5, I) = R1Z*DD(1, I)
      DD(6, I) = zero
   end do
!-not;      if(pf)         write(6,*) 'jcoord: loop 10'
   DD(4, 1) = DD(4, 1) - RY
   DD(5, 1) = DD(5, 1) + RX
   R2X = SQRT(CM1X**2 + CM1Y**2 + CM1Z**2)
   R2Y = SQRT(CM2X**2 + CM2Y**2 + CM2Z**2)
   CC = one/(R2X + epsilo)
   CCS = one/(R2Y + epsilo)
   DO I = 1, 6
      DD(1, I) = (CM1X*DD(4, I) + CM1Y*DD(5, I) &
                  + CM1Z*DD(6, I))*CC
      DD(2, I) = (CM2X*DD(7, I) + CM2Y*DD(8, I) &
                  + CM2Z*DD(9, I))*CCS
   end do
!-not;      if(pf)         write(6,*) 'jcoord: Q(6) and derr(6,i) '
   RZ = CC*CCS
   Q(6) = (CM1X*CM2X + CM1Y*CM2Y + CM1Z*CM2Z)*RZ
   DO I = 1, 6
      DERR(6, I) = (CM1X*DD(7, I) + DD(4, I)*CM2X &
                    + CM1Y*DD(8, I) + DD(5, I)*CM2Y &
                    + CM1Z*DD(9, I) + DD(6, I)*CM2Z &
                    - (R2X*DD(2, I) + DD(1, I)*R2Y)*Q(6))*RZ
   end do
!-not;      if(pf)       write(6,*) ' leaving  jacoord'
   return
end
!
!**********************************************************************
!
subroutine geomet_h4bmkp(r, cc, ideriv, ierr)
!----------------------------------------------------------------------
!  Calculate values and derivatives for Jacobi coordinates in geo(3,6).
!----------------------------------------------------------------------
!  First, calculates the 6 internuclear distances r(i) in order, from
!  the cartesian coordinates cc(4,3); then, calculates ra, rb, R,
!  cos(theta1), cos(theta2) and cos(phi) for each of the 3 ways of
!  getting H2+H2 from H4 (get derivatives too, if required):
!    geo(1,i),i=1,6 contains the values for AB+CD
!    geo(2,i)                               AC+BD
!    geo(3,i)                               AD+BC
!  For each of the three pairings (ig=1,2,3) of H4:
!    geo(ig,1)=ra            geo(ig,2)=rb            geo(ig,3)=R
!    geo(ig,4)=cos(theta1)   geo(ig,5)=cos(theta2)   geo(ig,6)=cos(phi)
!  r(1) = r12 = rAB,   r(2) = r13 = rAC,   r(3) = r14 = rAD,
!  r(4) = r23 = rBC,   r(5) = r24 = rBD,   r(6) = r34 = rCD
!  rv(k,i) contains the six internuclear vectors
!  rvp(k,i) contains the 2 rv vectors to be passed to DOT and CROSS
!  com(   ) contains the coords of the center-of-mass of each atom pair
!  RR(i)    contains the vector joining appropriate centres of mass
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension cc(4, 3), r(6)
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   common/c_rcalc_indices_h4bmkp/ipart1(6), ipart2(6)
   save/c_rcalc_indices_h4bmkp/
!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
   parameter(epsilon=1.d-06, zh=0.5d0)
!
!-old;      parameter ( rmin_err = 0.301d0, tolcos = 1.d0 - 5.d-13 )
!-old;      parameter ( rmin_err = 0.301d0, tolcos = 1.d0 - 5.d-15 )
!-prev-no-level;      parameter ( rmin_err = 0.301d0, tolsin = 1.d-7 )
   parameter(tolsin=1.d-7)
!
   dimension RRv(3), rvp(2, 3), rv(6, 3), com(6, 3), dpdrvp(2, 3), &
      dpdRRv(3), rjac(6), qjac(6), dqdrjac(6, 6)
!
   common/c_geo_sin_h4bmkp/geo_sin(3, 0:5)
   save/c_geo_sin_h4bmkp/
!
   common/c_cross_h4bmkp/pt(2, 3), ptm(2)
   save/c_cross_h4bmkp/
!
   common/c_geomet_local_h4bmkp/drvdcc(6, 3, 12), dcomdcc(6, 3, 12)
   save/c_geomet_local_h4bmkp/
!
   common/c_fnlclose_h4bmkp/fnlclose, dfnlclosedri(6)
   save/c_fnlclose_h4bmkp/
!
!__Calculate six internuclear vectors and distances based on cartesian coords:
!

!   print*,'ccc',cc
!   print*,'rrr',r
   do k = 1, 6
      i1 = ipart1(k)
      i2 = ipart2(k)
      rv(k, 1) = cc(i2, 1) - cc(i1, 1)
      rv(k, 2) = cc(i2, 2) - cc(i1, 2)
      rv(k, 3) = cc(i2, 3) - cc(i1, 3)
      r(k) = sqrt(rv(k, 1)**2 + rv(k, 2)**2 + rv(k, 3)**2)
   end do
!
!__Get small-distance correction factor, and smallest-distance index ierr:
!
   call fnlclose_h4bmkp(r, ideriv, ierr)
!
!__Check for error (smallest distance too short):
!
   if (r(ierr) .lt. rmin_err) return
!
   ierr = 0
!
!__Calculate cartesian coordinates of center-of-mass of each atom pair:
!
   do i = 1, 3
      com(1, i) = zh*rv(1, i)
      com(2, i) = zh*rv(2, i)
      com(3, i) = zh*rv(3, i)
      com(4, i) = rv(1, i) + zh*rv(4, i)
      com(5, i) = rv(1, i) + zh*rv(5, i)
      com(6, i) = rv(2, i) + zh*rv(6, i)
   end do
!
!__Calculate derivatives of r(k), rv(k,j), com(i,j) with respect to cc(m,n):
!
   if (ideriv .gt. 0) then
!
      do k = 1, 6
!
         i1 = ipart1(k)
         i2 = ipart2(k)
!
         if (r(k) .gt. 0.d0) then
            dridcc(k, i2) = rv(k, 1)/r(k)
            dridcc(k, i2 + 4) = rv(k, 2)/r(k)
            dridcc(k, i2 + 8) = rv(k, 3)/r(k)
         else
            dridcc(k, i2) = 1.d0
            dridcc(k, i2 + 4) = 1.d0
            dridcc(k, i2 + 8) = 1.d0
         end if
         dridcc(k, i1) = -dridcc(k, i2)
         dridcc(k, i1 + 4) = -dridcc(k, i2 + 4)
         dridcc(k, i1 + 8) = -dridcc(k, i2 + 8)
!
         drvdcc(k, 1, i1) = -1.d0
         drvdcc(k, 1, i2) = 1.d0
         drvdcc(k, 2, i1 + 4) = -1.d0
         drvdcc(k, 2, i2 + 4) = 1.d0
         drvdcc(k, 3, i1 + 8) = -1.d0
         drvdcc(k, 3, i2 + 8) = 1.d0
!
      end do
!
      do i = 1, 3
      do j = 1, 12
         dcomdcc(1, i, j) = zh*drvdcc(1, i, j)
         dcomdcc(2, i, j) = zh*drvdcc(2, i, j)
         dcomdcc(3, i, j) = zh*drvdcc(3, i, j)
         dcomdcc(4, i, j) = drvdcc(1, i, j) + zh*drvdcc(4, i, j)
         dcomdcc(5, i, j) = drvdcc(1, i, j) + zh*drvdcc(5, i, j)
         dcomdcc(6, i, j) = drvdcc(2, i, j) + zh*drvdcc(6, i, j)
      end do
      end do
!
   end if
!
!__If only the London equation is being used, we don't need anything else.
!
   if (fnlclose .eq. 0.d0) return
!
!__Calculate ra, rb, R, cos(theta1), cos(theta2), cos(phi) for each of the
!  three ways ig of dividing the H4 into H2+H2 (i.e., AB+CD, AC+BD, AD+BC).
!
!  Note that, if R < 1.E-6 bohrs, then this way of decomposing into ra, rb, R
!  gets essentially zero weight, and there is no point in calculating angles.
!
   do ig = 1, 3
!                  ! get the decomposition indices: ra = r(l1), rb = r(l2)
      l1 = ig
      l2 = 7 - ig
!                                    ! get vectors of R, ra, rb
      do i = 1, 3
         RRv(i) = com(l2, i) - com(l1, i)
         rvp(1, i) = rv(l1, i)
         rvp(2, i) = rv(l2, i)
      end do
!                                          ! and their magnitudes
      r1 = r(l1)
      r2 = r(l2)
      R3 = sqrt(RRv(1)**2 + RRv(2)**2 + RRv(3)**2)
!                                          ! store ra, rb, R
      geo(ig, 1) = r1
      geo(ig, 2) = r2
      geo(ig, 3) = R3
!                              ! IF R(ig) < 1.E-6, set all else to zero:
      if (R3 .le. epsilon) then
!
         geo(ig, 4) = 0.d0
         geo(ig, 5) = 0.d0
         geo(ig, 6) = 0.d0
         do i = 0, 5
            geo_sin(ig, i) = 0.d0
         end do
!
         if (ideriv .gt. 0) then
!
            do i = 1, 6
            do j = 1, 12
               dgeodcc(ig, i, j) = 0.d0
            end do
            do j = 1, 6
               dgeodri(ig, i, j) = 0.d0
            end do
            end do
!
         end if
!                  ! OTHERWISE, i.e., for R(ig) > 1.E-6: calculate stuff:
      else
!                                            ! get theta cosines
         costh1 = (rv(l1, 1)*RRv(1) + rv(l1, 2)*RRv(2) &
                   + rv(l1, 3)*RRv(3))/(R3*r1)
         costh2 = (rv(l2, 1)*RRv(1) + rv(l2, 2)*RRv(2) &
                   + rv(l2, 3)*RRv(3))/(R3*r2)
!
!__Get the sines of theta1 and theta2 from the cross products; then:
!  if  min{ sin(theta1) , sin(theta2) }  >  1.E-7  , then get cos(phi) value:
!
         call cross_h4bmkp(rvp, RRv, cosphi, dpdrvp, dpdRRv, ideriv, iden)
!
!__Check whether sines are below threshold value of 1.E-7, and also store the
!  squares of sines and cosines; recalculate sines or cosines more accurately:
!
         isin1 = 1
         if (abs(costh1) .lt. 0.7d0) then
            geo_sin(ig, 0) = costh1**2
            geo_sin(ig, 2) = 1.d0 - geo_sin(ig, 0)
            geo_sin(ig, 4) = sqrt(geo_sin(ig, 2))
         else
            geo_sin(ig, 4) = ptm(1)/(r1*R3)
            geo_sin(ig, 2) = geo_sin(ig, 4)**2
            if (geo_sin(ig, 4) .lt. tolsin) isin1 = 0
            geo_sin(ig, 0) = 1.d0 - geo_sin(ig, 2)
            if (costh1 .ge. 0.d0) then
               costh1 = sqrt(geo_sin(ig, 0))
            else
               costh1 = -sqrt(geo_sin(ig, 0))
            end if
         end if
!
         isin2 = 1
         if (abs(costh2) .lt. 0.7d0) then
            geo_sin(ig, 1) = costh2**2
            geo_sin(ig, 3) = 1.d0 - geo_sin(ig, 1)
            geo_sin(ig, 5) = sqrt(geo_sin(ig, 3))
         else
            geo_sin(ig, 5) = ptm(2)/(r2*R3)
            geo_sin(ig, 3) = geo_sin(ig, 5)**2
            if (geo_sin(ig, 5) .lt. tolsin) isin2 = 0
            geo_sin(ig, 1) = 1.d0 - geo_sin(ig, 5)**2
            if (costh2 .ge. 0.d0) then
               !print*,'in h4s',geo_sin(ig,1)
               costh2 = sqrt(geo_sin(ig, 1))
            else
               costh2 = -sqrt(geo_sin(ig, 1))
            end if
         end if
!
         geo(ig, 4) = costh1
         geo(ig, 5) = costh2
!
!__If one of the sines of theta1,2 is < 1.E-7, or if cos(phi) could not be
!  calculated (e.g., due to short distance R and small sines), set cos(phi) = 0
!
         if (isin1 .eq. 0 .or. isin2 .eq. 0 .or. iden .eq. 0) then
!
            cosphi = 0.d0
            iden = 0
!
         end if
!
         geo(ig, 6) = cosphi
!
!__Calculate derivatives of the Jacobi coordinates geo(ig,*), if necessary:
!
         if (ideriv .gt. 0) then
!                              ! derivatives w.r.t cartesian coords:
            do j = 1, 12
!
               dRRv1 = dcomdcc(l2, 1, j) - dcomdcc(l1, 1, j)
               dRRv2 = dcomdcc(l2, 2, j) - dcomdcc(l1, 2, j)
               dRRv3 = dcomdcc(l2, 3, j) - dcomdcc(l1, 3, j)
!
               dgeodcc(ig, 1, j) = dridcc(l1, j)
!
               dgeodcc(ig, 2, j) = dridcc(l2, j)
!
               dgeodcc(ig, 3, j) = (RRv(1)*dRRv1 &
                                    + RRv(2)*dRRv2 + RRv(3)*dRRv3)/R3
!
               if (isin1 .eq. 0) then
                  dgeodcc(ig, 4, j) = 0.d0
               else
                  dgeodcc(ig, 4, j) = (drvdcc(l1, 1, j)*RRv(1) &
                                       + rv(l1, 1)*dRRv1 + drvdcc(l1, 2, j)*RRv(2) &
                                       + rv(l1, 2)*dRRv2 + drvdcc(l1, 3, j)*RRv(3) &
                                       + rv(l1, 3)*dRRv3)/(R3*r1) &
                                      - costh1*(dgeodcc(ig, 3, j)/R3 &
                                                + dridcc(l1, j)/r1)
               end if
!
               if (isin2 .eq. 0) then
                  dgeodcc(ig, 5, j) = 0.d0
               else
                  dgeodcc(ig, 5, j) = (drvdcc(l2, 1, j)*RRv(1) &
                                       + rv(l2, 1)*dRRv1 + drvdcc(l2, 2, j)*RRv(2) &
                                       + rv(l2, 2)*dRRv2 + drvdcc(l2, 3, j)*RRv(3) &
                                       + rv(l2, 3)*dRRv3)/(R3*r2) &
                                      - costh2*(dgeodcc(ig, 3, j)/R3 &
                                                + dridcc(l2, j)/r2)
               end if
!
               if (iden .eq. 0) then
                  dgeodcc(ig, 6, j) = 0.d0
               else
                  dgeodcc(ig, 6, j) = dpdrvp(1, 1)*drvdcc(l1, 1, j) &
                                      + dpdrvp(1, 2)*drvdcc(l1, 2, j) &
                                      + dpdrvp(1, 3)*drvdcc(l1, 3, j) &
                                      + dpdrvp(2, 1)*drvdcc(l2, 1, j) &
                                      + dpdrvp(2, 2)*drvdcc(l2, 2, j) &
                                      + dpdrvp(2, 3)*drvdcc(l2, 3, j) &
                                      + dpdRRv(1)*dRRv1 &
                                      + dpdRRv(2)*dRRv2 &
                                      + dpdRRv(3)*dRRv3
               end if
!
            end do
! derivatives w.r.t distances r():
            if (ideriv .gt. 1) then
!                                    ! get re-ordered distances
               do j = 1, 6
                  rjac(j) = r(jco_ri(j, ig))
               end do
!                                          ! get derivatives,
               call jcoord_h4bmkp(rjac, qjac, dqdrjac)
!                                          ! and re-order them
               do i = 1, 6
               do j = 1, 6
                  dgeodri(ig, jco_geo(i), jco_ri(j, ig)) = &
                     dqdrjac(i, j)
               end do
               end do
!                                    ! if arbitrary cos(phi):
               if (iden .eq. 0) then
               do j = 1, 6
                  dgeodri(ig, 6, j) = 0.d0
               end do
               end if
!
            end if
!
         end if
!
      end if
!
   end do
!
   return
end
!
!**********************************************************************
!
subroutine fnlclose_h4bmkp(r, ideriv, i_small)
!----------------------------------------------------------------------
!  Calculate close-distance switch factor, if necessary.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(6)
!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
   common/c_fnlclose_h4bmkp/fnlclose, dfnlclosedri(6)
   save/c_fnlclose_h4bmkp/
!
   common/c_rcl_h4bmkp/rlo_cl, rhi_cl, rmid_cl, dinv_cl, oinv_cl
   save/c_rcl_h4bmkp/
!
   dimension fnl_i(6), dfnldr_i(6)
!
   fnlclose = 1.d0
   i_small = 1
!
   do k = 1, 6
      if (r(k) .lt. r(i_small)) i_small = k
      dfnlclosedri(k) = 0.d0
   end do
!
!__For level_err = 1, get small-distance correction factor:
!
   if (level_err .ne. 1 .and. level_err .ne. -1) return
!
   if (r(i_small) .le. rlo_cl) then
!
      fnlclose = 0.d0
!
   else if (r(i_small) .lt. rhi_cl) then
!
      k = 1
      do while (fnlclose .gt. 0.d0 .and. k .le. 6)
      if (r(k) .ge. rhi_cl) then
         fnl_i(k) = 1.d0
      else if (r(k) .le. rlo_cl) then
         fnl_i(k) = 0.d0
      else if (r(k) .le. rmid_cl) then
         o_mxd = 0.5d0 - (r(k) - rmid_cl)*oinv_cl
         o_pxd = 1.d0 + (r(k) - rmid_cl)*dinv_cl
         fnl_i(k) = o_pxd**3*o_mxd
         if (ideriv .gt. 0) dfnldr_i(k) = oinv_cl &
                                          *o_pxd**2*(6.d0*o_mxd - o_pxd)
      else
         o_mxd = 1.d0 - (r(k) - rmid_cl)*dinv_cl
         o_pxd = 0.5d0 + (r(k) - rmid_cl)*oinv_cl
         fnl_i(k) = 1.d0 - o_mxd**3*o_pxd
         if (ideriv .gt. 0) dfnldr_i(k) = oinv_cl &
                                          *o_mxd**2*(6.d0*o_pxd - o_mxd)
      end if
      fnlclose = fnlclose*fnl_i(k)
      k = k + 1
      end do
!
      if (ideriv .gt. 0 .and. fnlclose .gt. 0.d0) then
!
         do k = 1, 6
         if (r(k) .ge. rhi_cl) then
            dfnlclosedri(k) = 0.d0
         else
            dfnlclosedri(k) = fnlclose*dfnldr_i(k)/fnl_i(k)
         end if
         end do
!
      end if
!
   end if
!
   return
end
!
!**********************************************************************
!
subroutine vcorrec_h4bmkp(ideriv)
!----------------------------------------------------------------------
!  Calculate (normalized) geometry fractions (goodness factors) fSA.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   common/c_vcorrec_pars_h4bmkp/beta_j
   save/c_vcorrec_pars_h4bmkp/
!
   parameter(epsilon=1.d-06, tolfsa=1.d-18)
!
   dimension sa(3), eta(3), etarg(3), dsadcc(3, 12), dsadri(3, 6)
!
!__Calculate unnormalized geometry fractions (goodness factors) sa:
!
   do ig = 1, 3
!                                   ! if R(ig) < 1.E-6 : sa(ig) = 0.0
      if (geo(ig, 3) .le. epsilon) then
!
         sa(ig) = 0.d0
!                        ! else: get value sa(ig) = exp( -beta_J * eta )
      else
!
         etarg(ig) = 1.d0/geo(ig, 1)**2 + 1.d0/geo(ig, 2)**2
         eta(ig) = 1.d0/(geo(ig, 3)*sqrt(etarg(ig)))
         sa(ig) = exp(-beta_j*eta(ig))
!
      end if
!
   end do
!
!__Set any negligibly small sa(ig) values (relative to other sa's) to zero.
!
   sasumtol = (sa(1) + sa(2) + sa(3))*tolfsa
!
   do ig = 1, 3
      if (sa(ig) .lt. sasumtol) sa(ig) = 0.d0
   end do
!
!__Calculate fraction (normalized) sa values fSA:
!
   sasum = sa(1) + sa(2) + sa(3)
!
   do ig = 1, 3
      fSA(ig) = sa(ig)/sasum
   end do
!
!__Calculate derivatives of fSA, if required:
!
   if (ideriv .gt. 0) then
!                              ! need all dsadcc,ri to get dfSAdcc,ri
      do ig = 1, 3
!
         if (sa(ig) .eq. 0.d0) then
!
            do j = 1, 12
               dfSAdcc(ig, j) = 0.d0
               dsadcc(ig, j) = 0.d0
            end do
!
            do j = 1, 6
               dfSAdri(ig, j) = 0.d0
               dsadri(ig, j) = 0.d0
            end do
!
         else
!
            dsadgeo_1 = -beta_j*sa(ig)*eta(ig) &
                        /(etarg(ig)*geo(ig, 1)**3)
            dsadgeo_2 = -beta_j*sa(ig)*eta(ig) &
                        /(etarg(ig)*geo(ig, 2)**3)
            dsadgeo_3 = beta_j*sa(ig)*eta(ig)/geo(ig, 3)
!
            do j = 1, 12
               dsadcc(ig, j) = dsadgeo_1*dgeodcc(ig, 1, j) &
                               + dsadgeo_2*dgeodcc(ig, 2, j) &
                               + dsadgeo_3*dgeodcc(ig, 3, j)
            end do
!
            if (ideriv .gt. 1) then
            do j = 1, 6
               dsadri(ig, j) = dsadgeo_1*dgeodri(ig, 1, j) &
                               + dsadgeo_2*dgeodri(ig, 2, j) &
                               + dsadgeo_3*dgeodri(ig, 3, j)
            end do
            end if
!
         end if
!
      end do
!
      do ig = 1, 3
!
         if (sa(ig) .gt. 0.d0) then
!                                          ! get dfSAdcc
            do j = 1, 12
               dfSAdcc(ig, j) = (dsadcc(ig, j) - fSA(ig) &
                                 *(dsadcc(1, j) + dsadcc(2, j) + dsadcc(3, j))) &
                                /sasum
            end do
!                                          ! get dfSAdri
            if (ideriv .gt. 1) then
            do j = 1, 6
               dfSAdri(ig, j) = (dsadri(ig, j) - fSA(ig) &
                                 *(dsadri(1, j) + dsadri(2, j) &
                                   + dsadri(3, j)))/sasum
            end do
            end if
!
         end if
!
      end do
!
   end if
!
   return
end
!
!**********************************************************************
!
subroutine cross_h4bmkp(rvp, RRv, cosphi, dpdrvp, dpdRRv, ideriv, iden)
!----------------------------------------------------------------------
!  Cross products of ra and rb with R; dot product of these -> cos(phi)
!----------------------------------------------------------------------
!  Calculate factors for sines of theta1 and theta2, then calculate the
!  cosine of the angle phi [if sin(theta1,2) > 0], and its derivatives.
!   rvp contains the two appropriate internuclear vectors
!   RRv contains the appropriate vector joining the c-o-m's
!   pt  contains the two cross products
!   ptm contains the magnitudes of the 2 cross-product vectors
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension rvp(2, 3), RRv(3), dpdrvp(2, 3), dpdRRv(3)
!
   common/c_cross_h4bmkp/pt(2, 3), ptm(2)
   save/c_cross_h4bmkp/
!
!-old;      parameter ( toldenom = 1.d-14 )
   parameter(toldenom=1.d-15)
!
!__Calculate magnitudes of both resulting cross product vectors:
!
   pt(1, 1) = rvp(1, 2)*RRv(3) - RRv(2)*rvp(1, 3)
   pt(1, 2) = RRv(1)*rvp(1, 3) - rvp(1, 1)*RRv(3)
   pt(1, 3) = rvp(1, 1)*RRv(2) - RRv(1)*rvp(1, 2)
!
   ptm1sq = pt(1, 1)**2 + pt(1, 2)**2 + pt(1, 3)**2
   ptm(1) = sqrt(ptm1sq)
!
   pt(2, 1) = rvp(2, 2)*RRv(3) - RRv(2)*rvp(2, 3)
   pt(2, 2) = RRv(1)*rvp(2, 3) - rvp(2, 1)*RRv(3)
   pt(2, 3) = rvp(2, 1)*RRv(2) - RRv(1)*rvp(2, 2)
!
   ptm2sq = pt(2, 1)**2 + pt(2, 2)**2 + pt(2, 3)**2
   ptm(2) = sqrt(ptm2sq)
!
!__Calculate cos(phi) only if both terms in denominator are > zero:
!
   denom = ptm(1)*ptm(2)
!
   if (min(ptm(1), ptm(2), denom) .gt. toldenom) then
!
      cosphi = (pt(1, 1)*pt(2, 1) + pt(1, 2)*pt(2, 2) &
                + pt(1, 3)*pt(2, 3))/denom
      iden = 1
!                              ! derivatives
      if (ideriv .gt. 0) then
!
         dpdrvp(1, 1) = &
            (RRv(2)*pt(2, 3) - RRv(3)*pt(2, 2))/denom &
            - cosphi*(RRv(2)*pt(1, 3) - RRv(3)*pt(1, 2)) &
            /ptm1sq
         dpdrvp(1, 2) = &
            (RRv(3)*pt(2, 1) - RRv(1)*pt(2, 3))/denom &
            - cosphi*(RRv(3)*pt(1, 1) - RRv(1)*pt(1, 3)) &
            /ptm1sq
         dpdrvp(1, 3) = &
            (RRv(1)*pt(2, 2) - RRv(2)*pt(2, 1))/denom &
            - cosphi*(RRv(1)*pt(1, 2) - RRv(2)*pt(1, 1)) &
            /ptm1sq
!
         dpdrvp(2, 1) = &
            (RRv(2)*pt(1, 3) - RRv(3)*pt(1, 2))/denom &
            - cosphi*(RRv(2)*pt(2, 3) - RRv(3)*pt(2, 2)) &
            /ptm2sq
         dpdrvp(2, 2) = &
            (RRv(3)*pt(1, 1) - RRv(1)*pt(1, 3))/denom &
            - cosphi*(RRv(3)*pt(2, 1) - RRv(1)*pt(2, 3)) &
            /ptm2sq
         dpdrvp(2, 3) = &
            (RRv(1)*pt(1, 2) - RRv(2)*pt(1, 1))/denom &
            - cosphi*(RRv(1)*pt(2, 2) - RRv(2)*pt(2, 1)) &
            /ptm2sq
!
         dpdRRv(1) = (rvp(1, 3)*pt(2, 2) - rvp(1, 2)*pt(2, 3) &
                      + rvp(2, 3)*pt(1, 2) - rvp(2, 2)*pt(1, 3))/denom &
                     - cosphi*( &
                     (rvp(1, 3)*pt(1, 2) - rvp(1, 2)*pt(1, 3))/ptm1sq &
                     + (rvp(2, 3)*pt(2, 2) - rvp(2, 2)*pt(2, 3))/ptm2sq &
                     )
         dpdRRv(2) = (rvp(1, 1)*pt(2, 3) - rvp(1, 3)*pt(2, 1) &
                      + rvp(2, 1)*pt(1, 3) - rvp(2, 3)*pt(1, 1))/denom &
                     - cosphi*( &
                     (rvp(1, 1)*pt(1, 3) - rvp(1, 3)*pt(1, 1))/ptm1sq &
                     + (rvp(2, 1)*pt(2, 3) - rvp(2, 3)*pt(2, 1))/ptm2sq &
                     )
         dpdRRv(3) = (rvp(1, 2)*pt(2, 1) - rvp(1, 1)*pt(2, 2) &
                      + rvp(2, 2)*pt(1, 1) - rvp(2, 1)*pt(1, 2))/denom &
                     - cosphi*( &
                     (rvp(1, 2)*pt(1, 1) - rvp(1, 1)*pt(1, 2))/ptm1sq &
                     + (rvp(2, 2)*pt(2, 1) - rvp(2, 1)*pt(2, 2))/ptm2sq &
                     )
!
      end if
!
   else
!                        ! (set cosphi to arbitrary value for zero denominator)
      cosphi = 0.d0
      iden = 0
!
   end if
!
   return
end
!
!**********************************************************************
!
subroutine Vnsb_h4bmkp(r, ideriv, Vnsb, dVnsbdri)
!----------------------------------------------------------------------
!  Calculate the Vnsb term, the contribution from the BKMP2 H3 surface.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)
!
   dimension r(6), dVnsbdri(6)
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   dimension dVabc(6), dVabd(6), dVacd(6), dVbcd(6)
!
!__Get contributions (including switch/weight factors) from each way of
!  subdividing H4 into H3 + H, and sum them up:
!
   r1 = r(1)
   r2 = r(2)
   r3 = r(3)
   r4 = r(4)
   r5 = r(5)
   r6 = r(6)
   call Vxyzcal_h4bmkp(r1, r4, r2, r3, r5, r6, ideriv, Vabc, dVabc)
   call Vxyzcal_h4bmkp(r1, r5, r3, r2, r4, r6, ideriv, Vabd, dVabd)
   call Vxyzcal_h4bmkp(r2, r6, r3, r1, r4, r5, ideriv, Vacd, dVacd)
   call Vxyzcal_h4bmkp(r4, r6, r5, r1, r2, r3, ideriv, Vbcd, dVbcd)
!
   Vnsb = Vabc + Vabd + Vacd + Vbcd
!
!__Compute the derivatives, if necessary (note re-ordering of distances):
!
   if (ideriv .gt. 0) then
!
      dVnsbdri(1) = dVabc(1) + dVabd(1) + dVacd(4) + dVbcd(4)
      dVnsbdri(2) = dVabc(3) + dVabd(4) + dVacd(1) + dVbcd(5)
      dVnsbdri(3) = dVabc(4) + dVabd(3) + dVacd(3) + dVbcd(6)
      dVnsbdri(4) = dVabc(2) + dVabd(5) + dVacd(5) + dVbcd(1)
      dVnsbdri(5) = dVabc(5) + dVabd(2) + dVacd(6) + dVbcd(3)
      dVnsbdri(6) = dVabc(6) + dVabd(6) + dVacd(2) + dVbcd(2)
!
   end if
!
   return
end
!
!**********************************************************************
!
subroutine vpcalc_h4bmkp(r, ideriv)
!----------------------------------------------------------------------
!  Caluclate Schwenke-type selector functions sumVp(ig).
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

   dimension r(6)
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
!__These two nonlinear parameters are set in block data coef_h4bmkp
!
   common/c_vpcalc_par_h4bmkp/ap_vpcalc_2, ap_vpcalc_3
   save/c_vpcalc_par_h4bmkp/
!
   dimension vp(6)
!                  ! the 6 pairwise non-bonding potential values Vp:
   do j = 1, 6
!
      vp(j) = exp(-ap_vpcalc_2*r(j)**ap_vpcalc_3)
!
   end do
!                                    ! appropriate sums of these
   sumvp(1) = vp(2) + vp(5) + vp(3) + vp(4)
   sumvp(2) = vp(1) + vp(6) + vp(3) + vp(4)
   sumvp(3) = vp(1) + vp(6) + vp(2) + vp(5)
!                                    ! and derivatives, if needed
   if (ideriv .gt. 0) then
!
      do j = 1, 6
!
         dvpdri = -ap_vpcalc_2*ap_vpcalc_3 &
                  *r(j)**(ap_vpcalc_3 - 1.d0)*vp(j)
!
         if (j .ne. 1 .and. j .ne. 6) dsumVpdri(1, j) = dvpdri
         if (j .ne. 2 .and. j .ne. 5) dsumVpdri(2, j) = dvpdri
         if (j .ne. 3 .and. j .ne. 4) dsumVpdri(3, j) = dvpdri
!
      end do
!
   end if
!
   return
end
!
!**********************************************************************
!
subroutine comlon_h4bmkp(r, ideriv, vlon, dVlondri, dVlondfSA)
!----------------------------------------------------------------------
!  The H4 London equation, with geometry-dependent cusp-rounding.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(6), dVlondri(6), dVlondfSA(3)
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
!__These now have data values in block data coef_h4bmkp:
!
   common/c_lon_h4bmkp/eps_lond, sflo_lond, sfhi_lond, sfmid_lond, &
      dhsf_lond, dhsfinv_lond, rdelta_lond, iswtyp_lond
   save/c_lon_h4bmkp/
!
   common/c_fnlclose_h4bmkp/fnlclose, dfnlclosedri(6)
   save/c_fnlclose_h4bmkp/
!
   dimension vQ(6), vJ(6), dQdri(6), dJdri(6), Etrip(3), Esing(3)
!
   dimension rpd(6), deps2dri(6), deps2dfSA(3), dsumJhalfdri(6)
!
   parameter(f_1o4=1.d0/4.d0, eps_h3lond=1.d-6)
!
   parameter(half=0.5d0)
!
!__Calculate cusp-rounding values as a function of the current geometry.
!  Note that, initially, deps2...() are actually derivatives of eps;
!  only later are they revised to become derivatives of eps2 (=eps1**2).
!
   do i = 1, 6
      rpd(i) = r(i) + rdelta_lond
      deps2dri(i) = 0.d0
   end do
   deps2dfSA(1) = 0.d0
   deps2dfSA(2) = 0.d0
   deps2dfSA(3) = 0.d0
!
   if (iswtyp_lond .eq. 0 .or. fnlclose .eq. 0.d0) then
!
      eps = eps_h3lond
!                        ! usual case: iswtyp_lond = 4
   else
!
      sfterm1 = rpd(2)*rpd(3)*rpd(4)*rpd(5)
      sfterm2 = rpd(1)*rpd(3)*rpd(4)*rpd(6)
      sfterm3 = rpd(1)*rpd(2)*rpd(5)*rpd(6)
!
      s_f = (fSA(1)*sfterm1 + fSA(2)*sfterm2 &
             + fSA(3)*sfterm3)**f_1o4
!
!-old;         s_f = ( fSA(1) * ( r(2) + rdelta_lond )
!-old;     $        * ( r(3) + rdelta_lond ) * ( r(4) + rdelta_lond )
!-old;     $        * ( r(5) + rdelta_lond )
!-old;     $        + fSA(2) * ( r(1) + rdelta_lond )
!-old;     $        * ( r(3) + rdelta_lond ) * ( r(4) + rdelta_lond )
!-old;     $        * ( r(6) + rdelta_lond )
!-old;     $        + fSA(3) * ( r(1) + rdelta_lond )
!-old;     $        * ( r(2) + rdelta_lond ) * ( r(5) + rdelta_lond )
!-old;     $        * ( r(6) + rdelta_lond ) )**f_1o4
!
      if (s_f .le. sflo_lond) then
         eps = eps_lond
      else if (s_f .ge. sfhi_lond) then
         eps = eps_h3lond
      else
         o_mxd = 1.d0 - (s_f - sfmid_lond)*dhsfinv_lond
         o_pxd = 1.d0 + (s_f - sfmid_lond)*dhsfinv_lond
         if (s_f .le. sfmid_lond) then
            sw_use = o_pxd**3*o_mxd*0.5d0
            o_msw = 1.d0 - sw_use
            if (ideriv .gt. 0) derfac = 0.5d0*dhsfinv_lond &
                                        *o_pxd**2*(3.d0*o_mxd - o_pxd) &
                                        *(eps_h3lond - eps_lond)
         else
            o_msw = o_mxd**3*o_pxd*0.5d0
            sw_use = 1.d0 - o_msw
            if (ideriv .gt. 0) derfac = 0.5d0*dhsfinv_lond &
                                        *o_mxd**2*(o_mxd - 3.d0*o_pxd) &
                                        *(eps_lond - eps_h3lond)
         end if
         eps = o_msw*eps_lond + sw_use*eps_h3lond
         if (ideriv .gt. 0) then
            derfac = derfac*f_1o4/s_f**3
            deps2dri(1) = (fSA(2)*rpd(3)*rpd(4) &
                           + fSA(3)*rpd(2)*rpd(5))*rpd(6)*derfac
            deps2dri(2) = (fSA(1)*rpd(3)*rpd(4) &
                           + fSA(3)*rpd(1)*rpd(6))*rpd(5)*derfac
            deps2dri(3) = (fSA(1)*rpd(2)*rpd(5) &
                           + fSA(2)*rpd(1)*rpd(6))*rpd(4)*derfac
            deps2dri(4) = (fSA(1)*rpd(2)*rpd(5) &
                           + fSA(2)*rpd(1)*rpd(6))*rpd(3)*derfac
            deps2dri(5) = (fSA(1)*rpd(3)*rpd(4) &
                           + fSA(3)*rpd(1)*rpd(6))*rpd(2)*derfac
            deps2dri(6) = (fSA(2)*rpd(3)*rpd(4) &
                           + fSA(3)*rpd(2)*rpd(5))*rpd(1)*derfac
            deps2dfSA(1) = sfterm1*derfac
            deps2dfSA(2) = sfterm2*derfac
            deps2dfSA(3) = sfterm3*derfac
         end if
      end if
!
   end if
!
   eps1 = eps_h3lond + fnlclose*(eps - eps_h3lond)
   eps2 = eps1**2
!
!__Get derivatives of eps1**2:
!
   if (ideriv .gt. 0 .and. fnlclose .gt. 0.d0) then
   do i = 1, 6
      deps2dri(i) = 2.d0*eps1*(fnlclose*deps2dri(i) &
                               + (eps - eps_h3lond)*dfnlclosedri(i))
   end do
   deps2dfSA(1) = 2.d0*eps1*fnlclose*deps2dfSA(1)
   deps2dfSA(2) = 2.d0*eps1*fnlclose*deps2dfSA(2)
   deps2dfSA(3) = 2.d0*eps1*fnlclose*deps2dfSA(3)
   end if
!
!__Calculate the singlet and triplet energy for each of 6 distances:
!
   do i = 1, 6
!
!__Use Schwenke's H2 potential for the singlet state:
!
      call vH2opt95_h4bmkp(r(i), Esing, ideriv)
!
!__Use our triplet equation with the johnson correction:
!
      call triplet95_h4bmkp(r(i), Etrip, ideriv)
!
      vQ(i) = half*(Esing(1) + Etrip(1))
      vJ(i) = half*(Esing(1) - Etrip(1))
!
      if (ideriv .gt. 0) then
         dQdri(i) = half*(Esing(2) + Etrip(2))
         dJdri(i) = half*(Esing(2) - Etrip(2))
      end if
!
   end do
!
!__The H4 london equation:
!
   sumQ = vQ(1) + vQ(2) + vQ(3) + vQ(4) + vQ(5) + vQ(6)
!
   termJ1 = vJ(1) + vJ(6) - vJ(2) - vJ(5)
   termJ2 = vJ(2) + vJ(5) - vJ(3) - vJ(4)
   termJ3 = vJ(3) + vJ(4) - vJ(1) - vJ(6)
!
   sumJhalf = half*(termJ1**2 + termJ2**2 + termJ3**2)
!
   termJsub = sqrt(eps2 + sumJhalf)
!
   Vlon = sumQ - termJsub + eps_h3lond
!
!__The derivatives, if necessary:
!
   if (ideriv .gt. 0) then
!
      dsumJhalfdri(1) = (termJ1 - termJ3)*dJdri(1)
      dsumJhalfdri(2) = (termJ2 - termJ1)*dJdri(2)
      dsumJhalfdri(3) = (termJ3 - termJ2)*dJdri(3)
      dsumJhalfdri(4) = (termJ3 - termJ2)*dJdri(4)
      dsumJhalfdri(5) = (termJ2 - termJ1)*dJdri(5)
      dsumJhalfdri(6) = (termJ1 - termJ3)*dJdri(6)
!                                             ! Note: termJsub > 0
      derfac = -0.5d0/termJsub
!
      do i = 1, 6
         dVlondri(i) = dQdri(i) + (deps2dri(i) &
                                   + dsumJhalfdri(i))*derfac
      end do
!
      dVlondfSA(1) = deps2dfSA(1)*derfac
      dVlondfSA(2) = deps2dfSA(2)*derfac
      dVlondfSA(3) = deps2dfSA(3)*derfac
!
   end if
!
   return
end
!
!**********************************************************************
!
subroutine getYy_h4bmkp(ideriv)
!----------------------------------------------------------------------
!  Calculate the "curly-Y" functions YY(3,6) of theta1, theta2, phi.
!----------------------------------------------------------------------
!
!  For Yy(ig,k), the value of k indexes the set of indices q1,q2,mu:
!
!   k=1 ... 000 column    k=2 ... 220 column    k=3 ... 221 column
!   k=4 ... 222 column    k=5 ... 200 column
!   k=6 ... 020 column {must reverse the ij indices in the ijk sum}
!
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
   common/c_curlyY_par_h4bmkp/z4pi, z4pi_rt2, rt_zfl0m0, &
      rt_zfm0, rt_zfm1_3, rt_zfm2_3
   save/c_curlyY_par_h4bmkp/
!
   common/c_geo_sin_h4bmkp/geo_sin(3, 0:5)
   save/c_geo_sin_h4bmkp/
!
   parameter(zh=0.5d0, z3h=1.5d0, twoz3h=3.d0, tol=1.d-7)
!
   do ig = 1, 3
!
      if (fSA(ig) .gt. 0.d0) then
!
         ct1 = geo(ig, 4)
         ct2 = geo(ig, 5)
!-old;            ct1sq = ct1**2
!-old;            ct2sq = ct2**2
         ct1sq = geo_sin(ig, 0)
         ct2sq = geo_sin(ig, 1)
         cphi = geo(ig, 6)
!-old;            st1 = sqrt( max( 1.d0 - ct1**2 , 0.d0 ) )
!-old;            st2 = sqrt( max( 1.d0 - ct2**2 , 0.d0 ) )
         st1 = geo_sin(ig, 4)
         st2 = geo_sin(ig, 5)
         c2phi = 2.d0*cphi**2 - 1.d0
!
!__Get individual partial-spherical-harmonics of cos(theta1) and cos(theta2):
!
!                                                ! Y00
         y00 = rt_zfl0m0
!                                                ! Y20
         y20_1 = rt_zfm0*(z3h*ct1sq - zh)
         y20_2 = rt_zfm0*(z3h*ct2sq - zh)
!                                                ! Y21
         y21_1 = rt_zfm1_3*ct1*st1
         y21_2 = rt_zfm1_3*ct2*st2
!                                                ! Y22
!-old;            y22_1 = rt_zfm2_3 * ( 1.d0 - ct1**2 )
!-old;            y22_2 = rt_zfm2_3 * ( 1.d0 - ct2**2 )
         y22_1 = rt_zfm2_3*geo_sin(ig, 2)
         y22_2 = rt_zfm2_3*geo_sin(ig, 3)
!
!__Now work out curlyY values for each case :
!                                                ! Yy000
         Yy(ig, 1) = z4pi*y00**2
!                                                ! Yy220
         Yy(ig, 2) = z4pi*y20_1*y20_2
!                                                ! Yy221
         Yy(ig, 3) = z4pi_rt2*y21_1*y21_2*cphi
!                                                ! Yy222
         Yy(ig, 4) = z4pi_rt2*y22_1*y22_2*c2phi
!                                                ! Yy200
         Yy(ig, 5) = z4pi*y20_1*y00
!                                                ! Yy020
         Yy(ig, 6) = z4pi*y00*y20_2
!
!__Derivatives.  Note that dYydgeo(ig,3,4) and dYydgeo(ig,3,5) approach
!  negative infinity as cos(theta) approaches 1 or -1, but the derivative with
!  respect to theta remains finite.  Note however that the value of cos(phi)
!  for these cases is undefined; if set to zero, these derivatives then also
!  go to zero, and one need not worry about them.  Even when they are large,
!  they are multiplied by a small value, and never result in erroneous values.
!
         if (ideriv .gt. 0) then
!                                                      ! /dgeo
            dYydgeo(ig, 2, 4) = z4pi*rt_zfm0*twoz3h*ct1*y20_2
!
            dYydgeo(ig, 2, 5) = z4pi*y20_1*rt_zfm0*twoz3h*ct2
!
            if (st1 .gt. tol) then
               dYydgeo(ig, 3, 4) = z4pi_rt2*rt_zfm1_3*(st1 &
                                                       - ct1sq/st1)*y21_2*cphi
            else
               dYydgeo(ig, 3, 4) = 0.d0
            end if
!
            if (st2 .gt. tol) then
               dYydgeo(ig, 3, 5) = z4pi_rt2*y21_1*rt_zfm1_3 &
                                   *(st2 - ct2sq/st2)*cphi
            else
               dYydgeo(ig, 3, 5) = 0.d0
            end if
!
            dYydgeo(ig, 3, 6) = z4pi_rt2*y21_1*y21_2
!
            dYydgeo(ig, 4, 4) = -2.d0*z4pi_rt2 &
                                *rt_zfm2_3*ct1*y22_2*c2phi
!
            dYydgeo(ig, 4, 5) = -2.d0*z4pi_rt2 &
                                *y22_1*rt_zfm2_3*ct2*c2phi
!
            dYydgeo(ig, 4, 6) = z4pi_rt2*y22_1*y22_2*4.d0*cphi
!
            dYydgeo(ig, 5, 4) = z4pi*rt_zfm0*twoz3h*ct1*y00
!
            dYydgeo(ig, 6, 5) = z4pi*y00*rt_zfm0*twoz3h*ct2
!
         end if
!
      end if
!
   end do
!
   return
end
!
!**********************************************************************
!
subroutine Vlinear_h4bmkp(r, ideriv, Vlinear, Vlinmbe, &
                          dVlindri, dVlindgeo, dVmbedri)
!----------------------------------------------------------------------
!  Get the H4 linear correction terms and their derivatives.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(6), Vlinear(3), dVlindri(3, 6), dVlindgeo(3, 6), &
      dVmbedri(6)
!
   parameter(nterms=8, maxl=7*224 + 1350)
   common/c_coef_h4bmkp/coef(maxl), vsca(nterms), beta(nterms), &
      pow(nterms), iv(nterms), imin(nterms), imax(nterms), &
      kmin(nterms), kmax(nterms), irho(nterms), mlin, njac, &
      ikmin_lo, ikmax_hi
   save/c_coef_h4bmkp/
!
   parameter(n_ord_max=64, nc_ord=n_ord_max*4*nterms)
   common/c_ordcoef_h4bmkp/ordcoef(nc_ord), ncol(nterms), &
      nLL(nterms), nLL4(nterms), nkeep(nterms), ieqp(nterms)
   save/c_ordcoef_h4bmkp/
!
   dimension temp(n_ord_max), tcor(n_ord_max), &
      teda(n_ord_max), tcda(n_ord_max), &
      tedb(n_ord_max), tcdb(n_ord_max), &
      tedR(n_ord_max), tcdR(n_ord_max)
!
   parameter(nh_lo=-9, nh_hi=9)
   dimension r_a(nh_lo:nh_hi), r_b(nh_lo:nh_hi), r_R(nh_lo:nh_hi)
!
   dimension et(nterms), detdri(nterms, 6), &
      vs(6), vsa(6), vsb(6), vsR(6)
!
   common/c_geometry_h4bmkp/geo(3, 6), dridcc(6, 12), dgeodcc(3, 6, 12), &
      fSA(3), sumVp(3), dfSAdcc(3, 12), dfSAdri(3, 6), dsumVpdri(3, 6), &
      Yy(3, 6), dYydgeo(3, 6, 6), dgeodri(3, 6, 6), jco_ri(6, 3), jco_geo(6)
   save/c_geometry_h4bmkp/
!
!-old;      dimension Aselect(nterms), dAdri(nterms,6)
!
   dimension dAdri(nterms, 6)
!
   common/c_Aselect_h4bmkp/Aselect(nterms)
   save/c_Aselect_h4bmkp/
!
   parameter(R3toler=1.d-6)
!
!__Initialize to zero the MBE term Vmbe and linear-term-sums Vlinear(ig),ig=1,3
!
   Vmbe = 0.d0
   Vlinear(1) = 0.d0
   Vlinear(2) = 0.d0
   Vlinear(3) = 0.d0
!
!__Get rbar, rho, and selector functions A (except Acorr) and derivatives:
!
   call Asqsub_h4bmkp(r, ideriv, rbar, rho, Aselect, dAdri)
!
!__Calculate the exponential cut-off term for each correction term:
!
!                              ! calculate value for first term
   et(1) = exp(-beta(1)*rho**3)
!                              ! including derivatives, if necessary
   if (ideriv .gt. 0) then
!                                    ! note drho/dr(i) = r(i)/rho
      dds = -3.d0*beta(1)*rho*et(1)
      do i = 1, 6
         detdri(1, i) = dds*r(i)
         dVlindri(1, i) = 0.d0
         dVlindri(2, i) = 0.d0
         dVlindri(3, i) = 0.d0
         dVlindgeo(1, i) = 0.d0
         dVlindgeo(2, i) = 0.d0
         dVlindgeo(3, i) = 0.d0
      end do
!
   end if
!                  ! calculate values of exponential cut-off for other terms
   do kt = 2, 8
!                              ! case where beta(i) = 0.0
      if (beta(kt) .le. 0.d0) then
         et(kt) = 1.d0
         do i = 1, 6
            detdri(kt, i) = 0.d0
         end do
!                                    ! case where beta(i) = beta(1)
      else if (beta(kt) .eq. beta(1)) then
         et(kt) = et(1)
         do i = 1, 6
            detdri(kt, i) = detdri(1, i)
         end do
!                                    ! new beta(i) > 0.0
      else
         et(kt) = exp(-beta(kt)*rho**3)
         if (ideriv .gt. 0) then
            dds = -3.d0*beta(kt)*rho*et(kt)
            do i = 1, 6
               detdri(kt, i) = dds*r(i)
            end do
         end if
      end if
!
   end do
!
!__For each way ig=1,2,3 of dividing H4 into r_a, r_b, R: calculate the Vcorr
!  correction term, and V_2 through V_Njac (where Njac = 7 for BMKP H4, or 8):
!
!                                         ! lowest needed power iklo
   if (ikmin_lo .eq. 0 .and. ideriv .gt. 0) then
      iklo = 0
      r_a(-1) = 1.d0
      r_b(-1) = 1.d0
      r_B(-1) = 1.d0
   else
      iklo = ikmin_lo - min(ideriv, 1)
   end if
!                  ! (r**0 = 1.0)
   r_a(0) = 1.d0
   r_b(0) = 1.d0
   r_R(0) = 1.d0
!                   ! For each way ig=1,2,3 of dividing H4 into r_a, r_b, R:
   do ig = 1, 3
!                              ! If this way has non-zero weight fSA:
      if (fsa(ig) .gt. 0.d0) then
!                               ! calculate needed powers of r_a, r_b, R
         r_a(1) = geo(ig, 1)
         r_b(1) = geo(ig, 2)
         r_R(1) = geo(ig, 3)
         jm1 = 1
         do j = 2, ikmax_hi
            r_a(j) = r_a(jm1)*r_a(1)
            r_b(j) = r_b(jm1)*r_b(1)
            r_R(j) = r_R(jm1)*r_R(1)
            jm1 = j
         end do
         if (iklo .lt. 0) then
            jp1 = 0
            R_3 = max(r_R(1), R3toler)
            do j = -1, iklo, -1
               r_a(j) = r_a(jp1)/r_a(1)
               r_b(j) = r_b(jp1)/r_b(1)
               r_R(j) = r_R(jp1)/R_3
               jp1 = j
            end do
         end if
!                  ! no coefs have been used yet (re-ordered or not)
         nd0 = 0
         ndone = 0
!                              ! if first linear term Vcorr exists:
         if (iV(1) .gt. 0) then
!                              ! get factors r_a**i * r_b**j * R**k
            j = 0
            do kk = kmin(1), kmax(1)
            do ii = imin(1), imax(1)
            do jj = imin(1), imax(1)
               j = j + 1
               tcor(j) = r_a(ii)*r_b(jj)*r_R(kk)
               temp(j) = r_a(jj)*r_b(ii)*r_R(kk)
               if (ideriv .gt. 0) then
                  tcda(j) = ii*r_a(ii - 1)*r_b(jj)*r_R(kk)
                  tcdb(j) = jj*r_a(ii)*r_b(jj - 1)*r_R(kk)
                  tcdR(j) = kk*r_a(ii)*r_b(jj)*r_R(kk - 1)
                  teda(j) = jj*r_a(jj - 1)*r_b(ii)*r_R(kk)
                  tedb(j) = ii*r_a(jj)*r_b(ii - 1)*r_R(kk)
                  tedR(j) = kk*r_a(jj)*r_b(ii)*r_R(kk - 1)
               end if
            end do
            end do
            end do
!                          ! get starting positions for re-ordered coefs
!                          !  for q1q2mu = 000 (nd0), 220, 221, 222
            nd1 = nd0 + nLL(1)
            nd2 = nd1 + nLL(1)
            nd3 = nd2 + nLL(1)
!                          ! get starting position for non-re-ordered
!                          ! coefs (for both of q1q2mu = 020, 200)
            nd4 = nkeep(1)
!                        ! initialize sums to zero
            do i = 1, 6
               vs(i) = 0.d0
            end do
!                        ! sum over i,j,k of  r_a**i * r_b**j * R**k
!                        !                      * Coef_term1(q1q2mu,ijk)
            do j = 1, nLL(1)
               vs(1) = vs(1) + tcor(j)*ordcoef(j)
               vs(2) = vs(2) + tcor(j)*ordcoef(j + nd1)
               vs(3) = vs(3) + tcor(j)*ordcoef(j + nd2)
               vs(4) = vs(4) + tcor(j)*ordcoef(j + nd3)
               vs(5) = vs(5) + tcor(j)*coef(j + nd4)
               vs(6) = vs(6) + temp(j)*coef(j + nd4)
            end do
!                                          ! combine with curly-Y,
            v_s = vs(1)*YY(ig, 1) + vs(2)*YY(ig, 2) &
                  + vs(3)*YY(ig, 3) + vs(4)*YY(ig, 4) &
                  + vs(5)*YY(ig, 5) + vs(6)*YY(ig, 6)
!                                            ! sumVp, and cut-off
            Vlinear(ig) = sumVp(ig)*v_s*et(1)
!                                           ! derivatives, if needed
            if (ideriv .gt. 0) then
!
               do i = 1, 6
                  vsa(i) = 0.d0
                  vsb(i) = 0.d0
                  vsR(i) = 0.d0
               end do
!                                    ! accumulate derivative sums
               do j = 1, nLL(1)
                  vsa(1) = vsa(1) + tcda(j)*ordcoef(j)
                  vsa(2) = vsa(2) + tcda(j)*ordcoef(j + nd1)
                  vsa(3) = vsa(3) + tcda(j)*ordcoef(j + nd2)
                  vsa(4) = vsa(4) + tcda(j)*ordcoef(j + nd3)
                  vsa(5) = vsa(5) + tcda(j)*coef(j + nd4)
                  vsa(6) = vsa(6) + teda(j)*coef(j + nd4)
                  vsb(1) = vsb(1) + tcdb(j)*ordcoef(j)
                  vsb(2) = vsb(2) + tcdb(j)*ordcoef(j + nd1)
                  vsb(3) = vsb(3) + tcdb(j)*ordcoef(j + nd2)
                  vsb(4) = vsb(4) + tcdb(j)*ordcoef(j + nd3)
                  vsb(5) = vsb(5) + tcdb(j)*coef(j + nd4)
                  vsb(6) = vsb(6) + tedb(j)*coef(j + nd4)
                  vsR(1) = vsR(1) + tcdR(j)*ordcoef(j)
                  vsR(2) = vsR(2) + tcdR(j)*ordcoef(j + nd1)
                  vsR(3) = vsR(3) + tcdR(j)*ordcoef(j + nd2)
                  vsR(4) = vsR(4) + tcdR(j)*ordcoef(j + nd3)
                  vsR(5) = vsR(5) + tcdR(j)*coef(j + nd4)
                  vsR(6) = vsR(6) + tedR(j)*coef(j + nd4)
               end do
!                        ! derivatives w.r.t. r(*) at fixed geo(ig,*)
               do j = 1, 6
                  dVlindri(ig, j) = (sumVp(ig)*detdri(1, j) &
                                     + dsumVpdri(ig, j)*et(1))*v_s
               end do
!                                           ! and geo() at fixed r()
               dVlindgeo(ig, 1) = sumVp(ig)*et(1) &
                                  *(vsa(1)*YY(ig, 1) + vsa(2)*YY(ig, 2) &
                                    + vsa(3)*YY(ig, 3) + vsa(4)*YY(ig, 4) &
                                    + vsa(5)*YY(ig, 5) + vsa(6)*YY(ig, 6))
!
               dVlindgeo(ig, 2) = sumVp(ig)*et(1) &
                                  *(vsb(1)*YY(ig, 1) + vsb(2)*YY(ig, 2) &
                                    + vsb(3)*YY(ig, 3) + vsb(4)*YY(ig, 4) &
                                    + vsb(5)*YY(ig, 5) + vsb(6)*YY(ig, 6))
!
               dVlindgeo(ig, 3) = sumVp(ig)*et(1) &
                                  *(vsR(1)*YY(ig, 1) + vsR(2)*YY(ig, 2) &
                                    + vsR(3)*YY(ig, 3) + vsR(4)*YY(ig, 4) &
                                    + vsR(5)*YY(ig, 5) + vsR(6)*YY(ig, 6))
!
               dVlindgeo(ig, 4) = (vs(2)*dYydgeo(ig, 2, 4) &
                                   + vs(3)*dYydgeo(ig, 3, 4) &
                                   + vs(4)*dYydgeo(ig, 4, 4) &
                                   + vs(5)*dYydgeo(ig, 5, 4))*sumVp(ig)*et(1)
!
               dVlindgeo(ig, 5) = (vs(2)*dYydgeo(ig, 2, 5) &
                                   + vs(3)*dYydgeo(ig, 3, 5) &
                                   + vs(4)*dYydgeo(ig, 4, 5) &
                                   + vs(6)*dYydgeo(ig, 6, 5))*sumVp(ig)*et(1)
!
               dVlindgeo(ig, 6) = (vs(3)*dYydgeo(ig, 3, 6) &
                                   + vs(4)*dYydgeo(ig, 4, 6))*sumVp(ig)*et(1)
!
            end if
!                        ! have now used the coefs for first linear term
            ndone = ncol(1)
            nd0 = nLL4(1)
!
         end if
!                        ! for all other linear Jacobi terms:
         do k = 2, njac
!                              ! if this term exists, compute it:
            if (iV(k) .gt. 0) then
!                                   ! if index ranges differ from prev
               if (ieqp(k) .eq. 0) then
!                                    ! then recompute power factors
                  j = 0
                  do kk = kmin(k), kmax(k)
                  do ii = imin(k), imax(k)
                  do jj = imin(k), imax(k)
                     j = j + 1
                     tcor(j) = r_a(ii)*r_b(jj)*r_R(kk)
                     temp(j) = r_a(jj)*r_b(ii)*r_R(kk)
                     if (ideriv .gt. 0) then
                        tcda(j) = ii*r_a(ii - 1)*r_b(jj) &
                                  *r_R(kk)
                        tcdb(j) = jj*r_a(ii)*r_b(jj - 1) &
                                  *r_R(kk)
                        tcdR(j) = kk*r_a(ii)*r_b(jj) &
                                  *r_R(kk - 1)
                        teda(j) = jj*r_a(jj - 1)*r_b(ii) &
                                  *r_R(kk)
                        tedb(j) = ii*r_a(jj)*r_b(ii - 1) &
                                  *r_R(kk)
                        tedR(j) = kk*r_a(jj)*r_b(ii) &
                                  *r_R(kk - 1)
                     end if
                  end do
                  end do
                  end do
!
               end if
!                              ! get starting positions for coefs
               nd1 = nd0 + nLL(k)
               nd2 = nd1 + nLL(k)
               nd3 = nd2 + nLL(k)
               nd4 = ndone + nkeep(k)
!
               do i = 1, 6
                  vs(i) = 0.d0
               end do
!                              ! sums of r_a**i * r_b**j * R**k * Coef
               do j = 1, nLL(k)
                  vs(1) = vs(1) + tcor(j)*ordcoef(j + nd0)
                  vs(2) = vs(2) + tcor(j)*ordcoef(j + nd1)
                  vs(3) = vs(3) + tcor(j)*ordcoef(j + nd2)
                  vs(4) = vs(4) + tcor(j)*ordcoef(j + nd3)
                  vs(5) = vs(5) + tcor(j)*coef(j + nd4)
                  vs(6) = vs(6) + temp(j)*coef(j + nd4)
               end do
!                                                ! curly-Y
               v_s = vs(1)*YY(ig, 1) + vs(2)*YY(ig, 2) &
                     + vs(3)*YY(ig, 3) + vs(4)*YY(ig, 4) &
                     + vs(5)*YY(ig, 5) + vs(6)*YY(ig, 6)
!                                               ! add this term to
!                                                      ! Vlinear
               Vlinear(ig) = Vlinear(ig) + v_s*Aselect(k)*et(k)
!
               if (ideriv .gt. 0) then
!
                  do i = 1, 6
                     vsa(i) = 0.d0
                     vsb(i) = 0.d0
                     vsR(i) = 0.d0
                  end do
!                                    ! accumulate derivative sums
                  do j = 1, nLL(k)
                     vsa(1) = vsa(1) + tcda(j)*ordcoef(j + nd0)
                     vsa(2) = vsa(2) + tcda(j)*ordcoef(j + nd1)
                     vsa(3) = vsa(3) + tcda(j)*ordcoef(j + nd2)
                     vsa(4) = vsa(4) + tcda(j)*ordcoef(j + nd3)
                     vsa(5) = vsa(5) + tcda(j)*coef(j + nd4)
                     vsa(6) = vsa(6) + teda(j)*coef(j + nd4)
                     vsb(1) = vsb(1) + tcdb(j)*ordcoef(j + nd0)
                     vsb(2) = vsb(2) + tcdb(j)*ordcoef(j + nd1)
                     vsb(3) = vsb(3) + tcdb(j)*ordcoef(j + nd2)
                     vsb(4) = vsb(4) + tcdb(j)*ordcoef(j + nd3)
                     vsb(5) = vsb(5) + tcdb(j)*coef(j + nd4)
                     vsb(6) = vsb(6) + tedb(j)*coef(j + nd4)
                     vsR(1) = vsR(1) + tcdR(j)*ordcoef(j + nd0)
                     vsR(2) = vsR(2) + tcdR(j)*ordcoef(j + nd1)
                     vsR(3) = vsR(3) + tcdR(j)*ordcoef(j + nd2)
                     vsR(4) = vsR(4) + tcdR(j)*ordcoef(j + nd3)
                     vsR(5) = vsR(5) + tcdR(j)*coef(j + nd4)
                     vsR(6) = vsR(6) + tedR(j)*coef(j + nd4)
                  end do
!                              ! add derivatives of this term, for r()
                  do j = 1, 6
                     dVlindri(ig, j) = dVlindri(ig, j) &
                                       + (Aselect(k)*detdri(k, j) &
                                          + dAdri(k, j)*et(k))*v_s
                  end do
!                                          ! and for geo(ig,*)
                  dVlindgeo(ig, 1) = dVlindgeo(ig, 1) &
                                     + (vsa(1)*YY(ig, 1) + vsa(2)*YY(ig, 2) &
                                        + vsa(3)*YY(ig, 3) + vsa(4)*YY(ig, 4) &
                                        + vsa(5)*YY(ig, 5) + vsa(6)*YY(ig, 6)) &
                                     *Aselect(k)*et(k)
!
                  dVlindgeo(ig, 2) = dVlindgeo(ig, 2) &
                                     + (vsb(1)*YY(ig, 1) + vsb(2)*YY(ig, 2) &
                                        + vsb(3)*YY(ig, 3) + vsb(4)*YY(ig, 4) &
                                        + vsb(5)*YY(ig, 5) + vsb(6)*YY(ig, 6)) &
                                     *Aselect(k)*et(k)
!
                  dVlindgeo(ig, 3) = dVlindgeo(ig, 3) &
                                     + (vsR(1)*YY(ig, 1) + vsR(2)*YY(ig, 2) &
                                        + vsR(3)*YY(ig, 3) + vsR(4)*YY(ig, 4) &
                                        + vsR(5)*YY(ig, 5) + vsR(6)*YY(ig, 6)) &
                                     *Aselect(k)*et(k)
!
                  dVlindgeo(ig, 4) = dVlindgeo(ig, 4) &
                                     + (vs(2)*dYydgeo(ig, 2, 4) &
                                        + vs(3)*dYydgeo(ig, 3, 4) &
                                        + vs(4)*dYydgeo(ig, 4, 4) &
                                        + vs(5)*dYydgeo(ig, 5, 4)) &
                                     *Aselect(k)*et(k)
!
                  dVlindgeo(ig, 5) = dVlindgeo(ig, 5) &
                                     + (vs(2)*dYydgeo(ig, 2, 5) &
                                        + vs(3)*dYydgeo(ig, 3, 5) &
                                        + vs(4)*dYydgeo(ig, 4, 5) &
                                        + vs(6)*dYydgeo(ig, 6, 5)) &
                                     *Aselect(k)*et(k)
!
                  dVlindgeo(ig, 6) = dVlindgeo(ig, 6) &
                                     + (vs(3)*dYydgeo(ig, 3, 6) &
                                        + vs(4)*dYydgeo(ig, 4, 6)) &
                                     *Aselect(k)*et(k)
!
               end if
!                                ! have now used coefs for this term
               ndone = ndone + ncol(k)
               nd0 = nd0 + nLL4(k)
!
            end if
!
         end do
!
      end if
!
   end do
!
!__If needed, get the many body expansion (MBE) term Vmbe (and its derivatives)
!
   if (njac .lt. nterms .and. iV(nterms) .gt. 0) then
!
      call mbe_h4bmkp(r, ideriv, ndone, Vmbe, dVmbedri)
!
      Vlinmbe = Vmbe*Aselect(nterms)*et(nterms)
!
      if (ideriv .gt. 0) then
      do i = 1, 6
         dVmbedri(i) = dVmbedri(i)*Aselect(nterms)*et(nterms) &
                       + (dAdri(nterms, i)*et(nterms) &
                          + Aselect(nterms)*detdri(nterms, i))*Vmbe
      end do
      end if
!
   end if
!
   return
end
!
!**********************************************************************
!
subroutine Asqsub_h4bmkp(r, ideriv, rbar, rho, Aselect, dAdri)
!----------------------------------------------------------------------
!  Calculate rbar, rho, and selector functions A and their derivatives.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   parameter(nterms=8, maxl=7*224 + 1350)
   common/c_coef_h4bmkp/coef(maxl), vsca(nterms), beta(nterms), &
      pow(nterms), iv(nterms), imin(nterms), imax(nterms), &
      kmin(nterms), kmax(nterms), irho(nterms), mlin, njac, &
      ikmin_lo, ikmax_hi
   save/c_coef_h4bmkp/
!
   dimension r(6), Aselect(nterms), dAdri(nterms, 6)
!
   common/c_Achoice_h4bmkp/isqtet, igen1, icomp, inonlinh3, igen2
   save/c_Achoice_h4bmkp/
!
   parameter(fside=0.87867965644035752d0, fs2=2.d0*fside, &
             fdiag=1.2426406871192852d0, fd2=2.d0*fdiag)
!
   parameter(rr=1.15d0, rt=1.25d0, delta=0.1d0)
!
   parameter(f_1_3=1.d0/3.d0, f_2_3=2.d0/3.d0, &
             f13fs2=f_1_3*fs2, twomf13fs2=2.d0 - f13fs2, &
             f13fd2=f_1_3*fd2, twomf13fd2=2.d0 - f13fd2, &
             f_8_3=8.d0/3.d0)
!
!__Calculate average distance rbar, and rho (proportional to rms distance)
!
   r1 = r(1)
   r2 = r(2)
   r3 = r(3)
   r4 = r(4)
   r5 = r(5)
   r6 = r(6)
   rbar = (r1 + r2 + r3 + r4 + r5 + r6)/6.d0
   rbar2 = rbar*rbar
   rbar4 = rbar2*rbar2
   rbar8 = rbar4*rbar4
   rbar16 = rbar8*rbar8
   rho = sqrt(r1**2 + r2**2 + r3**2 + r4**2 + r5**2 + r6**2)
!
!__Calculate Atetsq, if needed (yes for BMKP H4):
!
   if (isqtet .gt. 0) then
!
!____calculate Atet (Atet=0 for tetrahedra):
!
      Atet = (r1 - rbar)**2 + (r2 - rbar)**2 + (r3 - rbar)**2 &
             + (r4 - rbar)**2 + (r5 - rbar)**2 + (r6 - rbar)**2
      Atetden = delta*rbar2 + Atet
      Atetx = rbar4*delta/Atetden
!
!____calculate Asq: Asq is a measure of 'non-squareness': (Asq=0 for squares)
!
      v16s = r1 + r6 - fs2*rbar
      v16d = r1 + r6 - fd2*rbar
      v25s = r2 + r5 - fs2*rbar
      v25d = r2 + r5 - fd2*rbar
      v34s = r3 + r4 - fs2*rbar
      v34d = r3 + r4 - fd2*rbar
!
      Asq0 = ((v16s*v16d)**2 + (v25s*v25d)**2 &
              + (v34s*v34d)**2)/rbar2
!
      Asq2 = Asq0 + (r1 - r6)**2 + (r2 - r5)**2 + (r3 - r4)**2
!
      Asqden = delta*rbar2 + Asq2
      Asq2x = rbar4*delta/Asqden
!
      Aselect(isqtet) = Atetx + Asq2x
!
      if (ideriv .gt. 0) then
!
         ddr1 = (twomf13fs2*v16s*v16d**2 &
                 + twomf13fd2*v16s**2*v16d &
                 - f13fs2*(v25s*v25d**2 + v34s*v34d**2) &
                 - f13fd2*(v25s**2*v25d + v34s**2*v34d)) &
                /rbar2
!
         ddr2 = (twomf13fs2*v25s*v25d**2 &
                 + twomf13fd2*v25s**2*v25d &
                 - f13fs2*(v16s*v16d**2 + v34s*v34d**2) &
                 - f13fd2*(v16s**2*v16d + v34s**2*v34d)) &
                /rbar2
!
         ddr3 = (twomf13fs2*v34s*v34d**2 &
                 + twomf13fd2*v34s**2*v34d &
                 - f13fs2*(v25s*v25d**2 + v16s*v16d**2) &
                 - f13fd2*(v25s**2*v25d + v16s**2*v16d)) &
                /rbar2
!
         dds = -Asq2x/Asqden
!
         dAdri(isqtet, 4) = (ddr3 + 2.d0*(r4 - r3))*dds
!
         dAdri(isqtet, 5) = (ddr2 + 2.d0*(r5 - r2))*dds
!
         dAdri(isqtet, 6) = (ddr1 + 2.d0*(r6 - r1))*dds
!
         dAdri(isqtet, 1) = (ddr1 + 2.d0*(r1 - r6))*dds
!
         dAdri(isqtet, 2) = (ddr2 + 2.d0*(r2 - r5))*dds
!
         dAdri(isqtet, 3) = (ddr3 + 2.d0*(r3 - r4))*dds
!
         ddd = (f_2_3/rbar + f_1_3*(Asq0/rbar &
                                    - delta*rbar)/Asqden)*Asq2x &
               + (f_2_3/rbar - delta*f_1_3*rbar/Atetden) &
               *Atetx
!
         dds = 2.d0*Atetx/Atetden
!
         do j = 1, 6
            dAdri(isqtet, j) = dAdri(isqtet, j) + ddd &
                               - (r(j) - rbar)*dds
         end do
!
      end if
!
   end if
!
!__Calculate Ah3h, if needed (yes for BMKP H4):
!
   if (inonlinh3 .gt. 0 .or. iV(4) .gt. 0) then
!
      w1 = (r1 + r2 + r4) - (r3 + r5 + r6)
      w2 = (r1 + r3 + r5) - (r2 + r4 + r6)
      w3 = (r2 + r3 + r6) - (r1 + r4 + r5)
      w4 = (r4 + r5 + r6) - (r1 + r2 + r3)
      Ah3h = w1*w1 + w2*w2 + w3*w3 + w4*w4
!
      Aselect(4) = Ah3h
!
      if (ideriv .gt. 0) then
!
         dAdri(4, 1) = 2.d0*(w1 + w2 - w3 - w4)
         dAdri(4, 2) = 2.d0*(w1 - w2 + w3 - w4)
         dAdri(4, 3) = 2.d0*(w2 - w1 + w3 - w4)
         dAdri(4, 4) = 2.d0*(w1 - w2 - w3 + w4)
         dAdri(4, 5) = 2.d0*(w2 - w1 - w3 + w4)
         dAdri(4, 6) = 2.d0*(w3 + w4 - w1 - w2)
!
      end if
!
   end if
!
!__Calculate Anonlinh3h, if needed (yes for BMKP H4):
!
   if (inonlinh3 .gt. 0) then
!
      wn1a = r1 + r2 - r4
      wn1b = r1 + r4 - r2
      wn1c = r2 + r4 - r1
      wn1 = wn1a*wn1b*wn1c
!
      wn2a = r1 + r3 - r5
      wn2b = r1 + r5 - r3
      wn2c = r3 + r5 - r1
      wn2 = wn2a*wn2b*wn2c
!
      wn3a = r2 + r3 - r6
      wn3b = r2 + r6 - r3
      wn3c = r3 + r6 - r2
      wn3 = wn3a*wn3b*wn3c
!
      wn4a = r4 + r5 - r6
      wn4b = r4 + r6 - r5
      wn4c = r5 + r6 - r4
      wn4 = wn4a*wn4b*wn4c
!
      Anonlinh3h = (wn1**2*w1**12 + wn2**2*w2**12 &
                    + wn3**2*w3**12 + wn4**2*w4**12)/rbar16
!
      Aselect(inonlinh3) = Anonlinh3h
!
      if (ideriv .gt. 0) then
!
         wn1_w1_11 = wn1*w1**11
         wn2_w2_11 = wn2*w2**11
         wn3_w3_11 = wn3*w3**11
         wn4_w4_11 = wn4*w4**11
!
         dAdri(inonlinh3, 1) = (2.d0*(wn1b*wn1c + wn1a*wn1c &
                                      - wn1a*wn1b)*w1 + 12.d0*wn1)*wn1_w1_11 &
                               + (2.d0*(wn2b*wn2c + wn2a*wn2c &
                                        - wn2a*wn2b)*w2 + 12.d0*wn2)*wn2_w2_11 &
                               - 12.d0*(wn3*wn3_w3_11 + wn4*wn4_w4_11)
!
         dAdri(inonlinh3, 2) = (2.d0*(wn1b*wn1c - wn1a*wn1c &
                                      + wn1a*wn1b)*w1 + 12.d0*wn1)*wn1_w1_11 &
                               + (2.d0*(wn3b*wn3c + wn3a*wn3c &
                                        - wn3a*wn3b)*w3 + 12.d0*wn3)*wn3_w3_11 &
                               - 12.d0*(wn2*wn2_w2_11 + wn4*wn4_w4_11)
!
         dAdri(inonlinh3, 3) = (2.d0*(wn2b*wn2c - wn2a*wn2c &
                                      + wn2a*wn2b)*w2 + 12.d0*wn2)*wn2_w2_11 &
                               + (2.d0*(wn3b*wn3c - wn3a*wn3c &
                                        + wn3a*wn3b)*w3 + 12.d0*wn3)*wn3_w3_11 &
                               - 12.d0*(wn1*wn1_w1_11 + wn4*wn4_w4_11)
!
         dAdri(inonlinh3, 4) = (2.d0*(wn1a*wn1b + wn1a*wn1c &
                                      - wn1b*wn1c)*w1 + 12.d0*wn1)*wn1_w1_11 &
                               + (2.d0*(wn4b*wn4c + wn4a*wn4c &
                                        - wn4a*wn4b)*w4 + 12.d0*wn4)*wn4_w4_11 &
                               - 12.d0*(wn2*wn2_w2_11 + wn3*wn3_w3_11)
!
         dAdri(inonlinh3, 5) = (2.d0*(wn2a*wn2b + wn2a*wn2c &
                                      - wn2b*wn2c)*w2 + 12.d0*wn2)*wn2_w2_11 &
                               + (2.d0*(wn4b*wn4c - wn4a*wn4c &
                                        + wn4a*wn4b)*w4 + 12.d0*wn4)*wn4_w4_11 &
                               - 12.d0*(wn1*wn1_w1_11 + wn3*wn3_w3_11)
!
         dAdri(inonlinh3, 6) = (2.d0*(wn3a*wn3b + wn3a*wn3c &
                                      - wn3b*wn3c)*w3 + 12.d0*wn3)*wn3_w3_11 &
                               + (2.d0*(wn4a*wn4b + wn4a*wn4c &
                                        - wn4b*wn4c)*w4 + 12.d0*wn4)*wn4_w4_11 &
                               - 12.d0*(wn1*wn1_w1_11 + wn2*wn2_w2_11)
!
         ddd = f_8_3*Anonlinh3h/rbar
!
         do j = 1, 6
            dAdri(inonlinh3, j) = dAdri(inonlinh3, j)/rbar16 - ddd
         end do
!
      end if
!
   end if
!
!__Calculate Ah2h2, if needed (yes for BMKP H4):
!
   if (iV(7) .gt. 0) then
!
      r1pr6 = r1 + r6
      r2pr5 = r2 + r5
      r3pr4 = r3 + r4
      r2345 = r2 + r3 + r4 + r5
      r1346 = r1 + r3 + r4 + r6
      r1256 = r1 + r2 + r5 + r6
!
      Ah2h2 = r1pr6**2 - r2345**2 + r2pr5**2 - r1346**2 &
              + r3pr4**2 - r1256**2
!
      Aselect(7) = Ah2h2
!
      if (ideriv .gt. 0) then
!
         dAdri(7, 1) = 2.d0*(r1pr6 - r1346 - r1256)
         dAdri(7, 2) = 2.d0*(r2pr5 - r2345 - r1256)
         dAdri(7, 3) = 2.d0*(r3pr4 - r2345 - r1346)
         dAdri(7, 4) = 2.d0*(r3pr4 - r2345 - r1346)
         dAdri(7, 5) = 2.d0*(r2pr5 - r2345 - r1256)
         dAdri(7, 6) = 2.d0*(r1pr6 - r1346 - r1256)
!
      end if
!
   end if
!
!__Calculate Agen, if needed (yes for BMKP H4):
!
   if (igen1 .gt. 0) then
!
      Agen = 1.d0
      Aselect(igen1) = Agen
!
      if (ideriv .gt. 0) then
      do j = 1, 6
         dAdri(igen1, j) = 0.d0
      end do
      end if
!
   end if
!
!__Calculate Akite, if needed (yes for BMKP H4):
!
   if (iV(2) .gt. 0) then
!
      ta1d = (r1 - r5)**2 + rbar2
      ta2d = (r2 - r6)**2 + rbar2
      rat43sq = (r4/r3)**2
      ta = rat43sq*(rbar2/ta1d + rbar2/ta2d)
!
      tb1d = (r1 - r2)**2 + rbar2
      tb2d = (r5 - r6)**2 + rbar2
      rat34sq = (r3/r4)**2
      tb = rat34sq*(rbar2/tb1d + rbar2/tb2d)
! --> 61:
      tc1d = (r2 - r4)**2 + rbar2
      tc2d = (r3 - r5)**2 + rbar2
!-err;         rat16sq = (r1/r6)**2
!-err;         tc = rat16sq * ( rbar2 / tc1d + rbar2 / tc2d )
      rat61sq = (r6/r1)**2
      tc = rat61sq*(rbar2/tc1d + rbar2/tc2d)
! --> 16:
      td1d = (r2 - r3)**2 + rbar2
      td2d = (r4 - r5)**2 + rbar2
!-err;         rat61sq = (r6/r1)**2
!-err;         td = rat61sq * ( rbar2 / td1d + rbar2 / td2d )
      rat16sq = (r1/r6)**2
      td = rat16sq*(rbar2/td1d + rbar2/td2d)
! --> 52:
      te1d = (r1 - r4)**2 + rbar2
      te2d = (r3 - r6)**2 + rbar2
!-err;         rat25sq = (r2/r5)**2
!-err;         te = rat25sq * ( rbar2 / te1d + rbar2 / te2d )
      rat52sq = (r5/r2)**2
      te = rat52sq*(rbar2/te1d + rbar2/te2d)
! --> 25:
      tf1d = (r1 - r3)**2 + rbar2
      tf2d = (r4 - r6)**2 + rbar2
!-err;         rat52sq = (r5/r2)**2
!-err;         tf = rat52sq * ( rbar2 / tf1d + rbar2 / tf2d )
      rat25sq = (r2/r5)**2
      tf = rat25sq*(rbar2/tf1d + rbar2/tf2d)
!
      Akite = ta + tb + tc + td + te + tf
!
      Aselect(2) = Akite
!
      if (ideriv .gt. 0) then
!
!-err;            ddd = f_1_3 * Akite / rbar - rbar2 * f_1_3 * rbar
!-err;     $           * ( rat43sq / ta1d**2 + rat43sq / ta2d**2
!-err;     $           + rat34sq / tb1d**2 + rat34sq / tb2d**2
!-err;     $           + rat16sq / tc1d**2 + rat16sq / tc2d**2
!-err;     $           + rat61sq / td1d**2 + rat61sq / td2d**2
!-err;     $           + rat25sq / te1d**2 + rat25sq / te2d**2
!-err;     $           + rat52sq / tf1d**2 + rat52sq / tf2d**2 )
!-err;c
!-err;            dAdri(2,1) = ddd + 2.d0 * ( ( tc - td ) / r1
!-err;     $           - ( rat43sq * (r1-r5) / ta1d**2
!-err;     $           + rat34sq * (r1-r2) / tb1d**2
!-err;     $           + rat25sq * (r1-r4) / te1d**2
!-err;     $           + rat52sq * (r1-r3) / tf1d**2 ) * rbar2 )
!-err;c
!-err;            dAdri(2,2) = ddd + 2.d0 * ( ( te - tf ) / r2
!-err;     $           - ( rat43sq * (r2-r6) / ta2d**2
!-err;     $           - rat34sq * (r1-r2) / tb1d**2
!-err;     $           + rat16sq * (r2-r4) / tc1d**2
!-err;     $           + rat61sq * (r2-r3) / td1d**2 ) * rbar2 )
!-err;c
!-err;            dAdri(2,3) = ddd + 2.d0 * ( ( tb - ta ) / r3
!-err;     $           - ( rat16sq * (r3-r5) / tc2d**2
!-err;     $           - rat61sq * (r2-r3) / td1d**2
!-err;     $           + rat25sq * (r3-r6) / te2d**2
!-err;     $           - rat52sq * (r1-r3) / tf1d**2 ) * rbar2 )
!-err;c
!-err;            dAdri(2,4) = ddd + 2.d0 * ( ( ta - tb ) / r4
!-err;     $           - ( rat61sq * (r4-r5) / td2d**2
!-err;     $           - rat16sq * (r2-r4) / tc1d**2
!-err;     $           - rat25sq * (r1-r4) / te1d**2
!-err;     $           + rat52sq * (r4-r6) / tf2d**2 ) * rbar2 )
!-err;c
!-err;            dAdri(2,5) = ddd + 2.d0 * ( ( tf - te ) / r5
!-err;     $           - ( rat34sq * (r5-r6) / tb2d**2
!-err;     $           - rat43sq * (r1-r5) / ta1d**2
!-err;     $           - rat16sq * (r3-r5) / tc2d**2
!-err;     $           - rat61sq * (r4-r5) / td2d**2 ) * rbar2 )
!-err;c
!-err;            dAdri(2,6) = ddd + 2.d0 * ( ( td - tc ) / r6
!-err;     $           - ( rat43sq * (r6-r2) / ta2d**2
!-err;     $           - rat34sq * (r5-r6) / tb2d**2
!-err;     $           - rat25sq * (r3-r6) / te2d**2
!-err;     $           - rat52sq * (r4-r6) / tf2d**2 ) * rbar2 )
!
         ddd = f_1_3*Akite/rbar - rbar2*f_1_3*rbar &
               *(rat43sq/ta1d**2 + rat43sq/ta2d**2 &
                 + rat34sq/tb1d**2 + rat34sq/tb2d**2 &
                 + rat61sq/tc1d**2 + rat61sq/tc2d**2 &
                 + rat16sq/td1d**2 + rat16sq/td2d**2 &
                 + rat52sq/te1d**2 + rat52sq/te2d**2 &
                 + rat25sq/tf1d**2 + rat25sq/tf2d**2)
!
         dAdri(2, 1) = ddd + 2.d0*((td - tc)/r1 &
                                   - (rat43sq*(r1 - r5)/ta1d**2 &
                                      + rat34sq*(r1 - r2)/tb1d**2 &
                                      + rat52sq*(r1 - r4)/te1d**2 &
                                      + rat25sq*(r1 - r3)/tf1d**2)*rbar2)
!
         dAdri(2, 2) = ddd + 2.d0*((tf - te)/r2 &
                                   - (rat43sq*(r2 - r6)/ta2d**2 &
                                      - rat34sq*(r1 - r2)/tb1d**2 &
                                      + rat61sq*(r2 - r4)/tc1d**2 &
                                      + rat16sq*(r2 - r3)/td1d**2)*rbar2)
!
         dAdri(2, 3) = ddd + 2.d0*((tb - ta)/r3 &
                                   - (rat61sq*(r3 - r5)/tc2d**2 &
                                      - rat16sq*(r2 - r3)/td1d**2 &
                                      + rat52sq*(r3 - r6)/te2d**2 &
                                      - rat25sq*(r1 - r3)/tf1d**2)*rbar2)
!
         dAdri(2, 4) = ddd + 2.d0*((ta - tb)/r4 &
                                   - (rat16sq*(r4 - r5)/td2d**2 &
                                      - rat61sq*(r2 - r4)/tc1d**2 &
                                      - rat52sq*(r1 - r4)/te1d**2 &
                                      + rat25sq*(r4 - r6)/tf2d**2)*rbar2)
!
         dAdri(2, 5) = ddd + 2.d0*((te - tf)/r5 &
                                   - (rat34sq*(r5 - r6)/tb2d**2 &
                                      - rat43sq*(r1 - r5)/ta1d**2 &
                                      - rat61sq*(r3 - r5)/tc2d**2 &
                                      - rat16sq*(r4 - r5)/td2d**2)*rbar2)
!
         dAdri(2, 6) = ddd + 2.d0*((tc - td)/r6 &
                                   - (rat43sq*(r6 - r2)/ta2d**2 &
                                      - rat34sq*(r5 - r6)/tb2d**2 &
                                      - rat52sq*(r3 - r6)/te2d**2 &
                                      - rat25sq*(r4 - r6)/tf2d**2)*rbar2)
!
      end if
!
   end if
!
!__Calculate Acomp, if needed (NOT needed for BMKP H4 surface):
!
   if (icomp .gt. 0) then
!
      Acomp = 0.d0
      do i = 1, 6
         dAdri(icomp, i) = 0.d0
         if (r(i) .lt. rr) then
            ri = r(i)
            Acomp = Acomp + (rr - ri)**3/(rt - ri)
            if (ideriv .gt. 0) &
               dAdri(icomp, i) = (rr - ri)**3/(rt - ri)**2 &
                                 - 3.d0*(rr - ri)**2/(rt - ri)
         end if
      end do
!
      Aselect(icomp) = Acomp
!
   end if
!
!__Calclulate Aeight, if needed (yes for BMKP H4, since iflag_Aterms = 4)
!
   if (igen2 .gt. 0) then
!
      Aselect(igen2) = 1.d0
!
      if (ideriv .gt. 0) then
      do j = 1, 6
         dAdri(igen2, j) = 0.d0
      end do
      end if
!
   end if
!
   return
end
!
!**********************************************************************
!
subroutine mbe_h4bmkp(r, ideriv, ndone, Vmbe, dVmbedri)
!----------------------------------------------------------------------
!  Calculate Aguado-type many-body-expansion term Vmbe and derivatives.
!----------------------------------------------------------------------
!
!  Permutations of distances resulting from the permutations of atoms:
!
!  H-permutation  Rab  Rac  Rad  Rbc  Rbd  Rcd
!  ------------- --1- --2- --3- --4- --5- --6-
!   1      abcd: r(1) r(2) r(3) r(4) r(5) r(6)
!   2      abdc:   1    3    2    5    4    6
!   3      acbd:   2    1    3    4    6    5
!   4      acdb:   2    3    1    6    4    5
!   5      adbc:   3    1    2    5    6    4
!   6      adcb:   3    2    1    6    5    4
!   7      bacd:   1    4    5    2    3    6
!   8      badc:   1    5    4    3    2    6
!   9      bcad:   4    1    5    2    6    3
!  10      bcda:   4    5    1    6    2    3
!  11      bdac:   5    1    4    3    6    2
!  12      bdca:   5    4    1    6    3    2
!  13      cabd:   2    4    6    1    3    5
!  14      cadb:   2    6    4    3    1    5
!  15      cbad:   4    2    6    1    5    3
!  16      cbda:   4    6    2    5    1    3
!  17      cdab:   6    2    4    3    5    1
!  18      cdba:   6    4    2    5    3    1
!  19      dabc:   3    5    6    1    2    4
!  20      dacb:   3    6    5    2    1    4
!  21      dbac:   5    3    6    1    4    2
!  22      dbca:   5    6    3    4    1    2
!  23      dcab:   6    3    5    2    4    1
!  24      dcba:   6    5    3    4    2    1
!
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(6), dVmbedri(6)
!
!  Note: values of beta_p and indices stored in block data coef_h4bmkp
!
   parameter(max_hi=12, &
             kind3hi=max_hi*max_hi, kind4hi=kind3hi*max_hi, &
             kind5hi=3000, kind6hi=3000)
   common/c_mbe_index_h4bmkp/beta_p, mbe_perm(6, 24), nterms_use, &
      indmin, indmax, maxsum, k1_lo, k1_hi, k2_lo(max_hi), &
      k2_hi(max_hi), k3_lo(kind3hi), k3_hi(kind3hi), &
      k4_lo(kind4hi), k4_hi(kind4hi), k5_lo(kind5hi), &
      k5_hi(kind5hi), k6_lo(kind6hi), k6_hi(kind6hi)
   save/c_mbe_index_h4bmkp/
!
   parameter(nterms=8, maxl=7*224 + 1350)
   common/c_coef_h4bmkp/coef(maxl), vsca(nterms), beta(nterms), &
      pow(nterms), iv(nterms), imin(nterms), imax(nterms), &
      kmin(nterms), kmax(nterms), irho(nterms), mlin, njac, &
      ikmin_lo, ikmax_hi
   save/c_coef_h4bmkp/
!
   dimension pi(6), ppow(0:max_hi, 6)
!
   dimension val_1(24), val_2(24), val_3(24), val_4(24), val_5(24), &
      val_6(24), dlnvdri(24, 6), dlnpdri(6)
!
!__Get necessary powers of [ r_i * exp( -beta_p * r_i) ] for each distance r_i.
!
   do i = 1, 6
      pi(i) = r(i)*exp(-beta_p*r(i))
      ppow(0, i) = 1.d0
      ppow(1, i) = pi(i)
      do j = 2, indmax
         ppow(j, i) = ppow(j - 1, i)*pi(i)
      end do
      dVmbedri(i) = 0.d0
      dlnpdri(i) = 1.d0/r(i) - beta_p
   end do
!
!__Accumulate MBE term Vmbe (and derivatives, if necessary):
!
   Vmbe = 0.d0
!                  ! current position in coefficent array
   kcoef = ndone
!                  ! positions in sum-index-limits arrays
   kind3 = 0
   kind4 = 0
   kind5 = 0
   kind6 = 0
!                        ! first index (power) i1:
   do i1 = k1_lo, k1_hi
!                        ! begin factors for each of the 24 permutations
      do i = 1, 24
         val_1(i) = ppow(i1, mbe_perm(1, i))
      end do
      if (ideriv .gt. 0) then
      if (i1 .eq. 0) then
      do i = 1, 24
         dlnvdri(i, mbe_perm(1, i)) = 0.d0
      end do
      else
      do i = 1, 24
         dlnvdri(i, mbe_perm(1, i)) = &
            i1*dlnpdri(mbe_perm(1, i))
      end do
      end if
      end if
!                              ! second index (power) i2:
      do i2 = k2_lo(i1), k2_hi(i1)
!                              ! continue factors for each permutation
         do i = 1, 24
            val_2(i) = val_1(i)*ppow(i2, mbe_perm(2, i))
         end do
         if (ideriv .gt. 0) then
         if (i2 .eq. 0) then
         do i = 1, 24
            dlnvdri(i, mbe_perm(2, i)) = 0.d0
         end do
         else
         do i = 1, 24
            dlnvdri(i, mbe_perm(2, i)) = &
               i2*dlnpdri(mbe_perm(2, i))
         end do
         end if
         end if
!                        ! index for third-index-limits arrays
         kind3 = kind3 + 1
!                                    ! third index (power) i3:
         do i3 = k3_lo(kind3), k3_hi(kind3)
!
            do i = 1, 24
               val_3(i) = val_2(i)*ppow(i3, mbe_perm(3, i))
            end do
            if (ideriv .gt. 0) then
            if (i3 .eq. 0) then
            do i = 1, 24
               dlnvdri(i, mbe_perm(3, i)) = 0.d0
            end do
            else
            do i = 1, 24
               dlnvdri(i, mbe_perm(3, i)) = &
                  i3*dlnpdri(mbe_perm(3, i))
            end do
            end if
            end if
!
            kind4 = kind4 + 1
!                                      ! fourth index (power) i4:
            do i4 = k4_lo(kind4), k4_hi(kind4)
!
               do i = 1, 24
                  val_4(i) = val_3(i)*ppow(i4, mbe_perm(4, i))
               end do
               if (ideriv .gt. 0) then
               if (i4 .eq. 0) then
               do i = 1, 24
                  dlnvdri(i, mbe_perm(4, i)) = 0.d0
               end do
               else
               do i = 1, 24
                  dlnvdri(i, mbe_perm(4, i)) = &
                     i4*dlnpdri(mbe_perm(4, i))
               end do
               end if
               end if
!
               kind5 = kind5 + 1
!                                          ! fifth index (power) i5:
               do i5 = k5_lo(kind5), k5_hi(kind5)
!
                  do i = 1, 24
                     val_5(i) = val_4(i)*ppow(i5, mbe_perm(5, i))
                  end do
                  if (ideriv .gt. 0) then
                  if (i5 .eq. 0) then
                  do i = 1, 24
                     dlnvdri(i, mbe_perm(5, i)) = 0.d0
                  end do
                  else
                  do i = 1, 24
                     dlnvdri(i, mbe_perm(5, i)) = &
                        i5*dlnpdri(mbe_perm(5, i))
                  end do
                  end if
                  end if
!
                  kind6 = kind6 + 1
!                                          ! sixth index (power) i6:
                  do i6 = k6_lo(kind6), k6_hi(kind6)
!                                          ! index for coefficient
                     kcoef = kcoef + 1
!                                            ! if coef is non-zero
                     if (coef(kcoef) .ne. 0.d0) then
!                                             ! accumulate factors
                        value = 0.d0
                        do i = 1, 24
                           val_6(i) = &
                              val_5(i)*ppow(i6, mbe_perm(6, i))
                           value = value + val_6(i)
                        end do
                        Vmbe = Vmbe + value*coef(kcoef)
!                                                ! and derivatives
                        if (ideriv .gt. 0) then
!
                           if (i6 .eq. 0) then
                           do i = 1, 24
                              dlnvdri(i, mbe_perm(6, i)) = 0.d0
                           end do
                           else
                           do i = 1, 24
                              dlnvdri(i, mbe_perm(6, i)) = &
                                 i6*dlnpdri(mbe_perm(6, i))
                           end do
                           end if
!
                           do j = 1, 6
                              value = 0.d0
                              do i = 1, 24
                                 value = value &
                                         + val_6(i)*dlnvdri(i, j)
                              end do
                              dVmbedri(j) = dVmbedri(j) &
                                            + value*coef(kcoef)
                           end do
!
                        end if
!
                     end if
!
                  end do
!
               end do
!
            end do
!
         end do
!
      end do
!
   end do
!
   return
!
end
!
!**********************************************************************
!
subroutine Vxyzcal_h4bmkp(rt1, rt2, rt3, rh1, rh2, rh3, ideriv, Vxyz, dV)
!----------------------------------------------------------------------
!  Calculate the 3-body term for given subdivision H3+H (from BKMP2 H3)
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension dV(6)
!
!__These switch values are now stored in block data coef_h4bmkp
!
   common/c_h3nsb_switch_h4bmkp/s0_nsb_1, st_nsb_1, sh_nsb_1, &
      s0_nsb_2, st_nsb_2, sh_nsb_2
   save/c_h3nsb_switch_h4bmkp/
!
   dimension dVasbe(3), dCalasbe(3), rt(3)
!
!__Get H3 weight/switch factors S_1 and S_2
!
   rotsq = 1.d0/rt1**2 + 1.d0/rt2**2 + 1.d0/rt3**2
   rohsq = 1.d0/rh1**2 + 1.d0/rh2**2 + 1.d0/rh3**2
!
   rot = sqrt(rotsq)
   roh = sqrt(rohsq)
!
   es1 = s0_nsb_1*exp(st_nsb_1/rot - sh_nsb_1/roh)
   es2 = s0_nsb_2*exp(st_nsb_2/rot - sh_nsb_2/roh)
!
   S_1 = 1.d0/(1.d0 + es1)
   S_2 = 1.d0/(1.d0 + es2)
!
   rt(1) = rt1
   rt(2) = rt2
   rt(3) = rt3
!
!__Get the H3 terms
!
   call h3bend_h4bmkp(rt, ideriv, Vasbe, dVasbe, Calasbe, dCalasbe)
!
!__Combine the H3 terms with the switch/weight factors
!
   Vxyz = S_1*Vasbe + S_2*Calasbe
!
!__Get the derivatives, if necessary
!
   if (ideriv .gt. 0) then
!
      f1 = Vasbe*S_1**2*es1
      f2 = Calasbe*S_2**2*es2
!
      derfac = (f1*st_nsb_1 + f2*st_nsb_2)/(rot*rotsq)
!
      dV(1) = S_1*dVasbe(1) + S_2*dCalasbe(1) - derfac/rt1**3
      dV(2) = S_1*dVasbe(2) + S_2*dCalasbe(2) - derfac/rt2**3
      dV(3) = S_1*dVasbe(3) + S_2*dCalasbe(3) - derfac/rt3**3
!
      derfac = (f1*sh_nsb_1 + f2*sh_nsb_2)/(roh*rohsq)
!
      dV(4) = derfac/rh1**3
      dV(5) = derfac/rh2**3
      dV(6) = derfac/rh3**3
!
   end if
!
   return
end
!
!**********************************************************************
!
!  This subset of the BKMP2 H3 surface computes only non-London terms,
!  for use in the BMKP analytic H4 surface.  The H2 singlet and triplet
!  routines below are called by the H4 surface, for its London equation.
!
!  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!  _/          surface950621 H3 surface evaluation routines          _/
!  _/              (also known as the bkmp2 H3 surface)              _/
!  _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
!
!
subroutine h3bend_h4bmkp(r, ideriv, Vasbe, dVasbe, Calasbe, dCalasbe)
!     ----------------------------------------------------------------
!
! updated 10june 1996 from bkmp1 (surf 706) to bkmp2
!----------------------------------------------------------------------
! Calculate total H3 potential from all of its parts
! if (ideriv.gt.0) also calculate the dV/dr derivatives
! all distances are in bohrs and all energies are in hartrees
!
! For a discussion of this surface, see:
!    A.I.Boothroyd, W.J.Keogh, P.G.Martin, M.R.Peterson
!    Journal of Chemical Physics 95 pp 4343-4359 (Sept15/1991)
!    and JCP 104 pp 7139-7152 (May8/1996)
!
! Note: this file contains parameters for a surface refitted
!       on june21/95 to a set of several thousand ab initio
!       points.  the 'surface706' surface parameters have been
!       removed.  The routine names have been modified
!       slightly (usually a '95' appended).
!
! Note: the surface parameter values as published in 1991 lead to
!       an anomolously deep van der Waals well for a very compact
!       H2 molecule (say r=0.8).  After that paper was submitted,
!       this problem was fixed and the corrected Cbend coefficients
!       are used in this version of the surface (version 706).
!
! any QUESTIONS/PROBLEMS/COMMENTS concerning this programme can be
! addressed to :  boothroy@cita.utoronto.ca  or  pgmartin@cita.utoronto.ca
!
! version:
! apr12/95 ... parameters for surface850308 added
! jul27/91 ... surf706d.out Cbend values put in
!----------------------------------------------------------------------
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

   dimension r(3), dVasbe(3), dCalasbe(3)
   dimension dVlon(3), dVas(3), dVbnda(3), dVbndb(3), &
      dCal(3), dCas(3), dCbnda(3), dCbndb(3)
   dimension dT(3), T(3), dAcalc(3)
!
!  zero everything to avoid any 'funny' values:
!
   Vlon = 0.d0
   Vas = 0.d0
   Vbnda = 0.d0
   Vbndb = 0.d0
   Cal = 0.d0
   Cas = 0.d0
   cbnda = 0.d0
   cbndb = 0.d0
!
   do i = 1, 3
      dVlon(i) = 0.d0
      dVas(i) = 0.d0
      dVbnda(i) = 0.d0
      dVbndb(i) = 0.d0
      dCal(i) = 0.d0
      dCas(i) = 0.d0
      dCbnda(i) = 0.d0
      dCbndb(i) = 0.d0
      dVasbe(i) = 0.d0
      dCalasbe(i) = 0.d0
      dAcalc(i) = 0.d0
   end do
!
!-noneed;      call h3lond95_h4bmkp( r, Vlon, dVlon )
!
   call vascal95_h4bmkp(r, Vas, dVas, ideriv, Acalc, dAcalc)
!
!    now do any corrections required for compact geometries:
!
   call compac95_h4bmkp(r, icompc, T, dT)
!
!    oct.3/90 compact routines only called for compact geometries:
!
   if (icompc .gt. 0) then
      call csym95_h4bmkp(r, cal, dCal)
      call casym95_h4bmkp(r, Cas, dCas, ideriv, T, dT, Acalc, dAcalc)
   end if
!
   call vbcb95_h4bmkp(r, icompc, T, dT, ideriv, &
                      Vbnda, Vbndb, dVbnda, dVbndb, &
                      Cbnda, Cbndb, dCbnda, dCbndb)
!
!    add up the various parts of the potential:
!
!
   Vasbe = Vas + Vbnda + Vbndb
   Calasbe = Cal + Cas + Cbnda + Cbndb
!
!    add up the various parts of the derivative:
!
   if (ideriv .gt. 0) then
      dVasbe(1) = dVas(1) + dVbnda(1) + dVbndb(1)
      dVasbe(2) = dVas(2) + dVbnda(2) + dVbndb(2)
      dVasbe(3) = dVas(3) + dVbnda(3) + dVbndb(3)
      dCalasbe(1) = dCal(1) + dCas(1) + dCbnda(1) + dCbndb(1)
      dCalasbe(2) = dCal(2) + dCas(2) + dCbnda(2) + dCbndb(2)
      dCalasbe(3) = dCal(3) + dCas(3) + dCbnda(3) + dCbndb(3)
   end if
!
   return
end
!
!**********************************************************************
!
subroutine triplet95_h4bmkp(r, E3, ideriv)
!--------------------------------------------------------------------
! apr12/95 surface950308 values added
! oct04/90 surface626 values added
!  H2 triplet curve and derivatives:
!  calculates triplet potential and first derivative
!  uses truhlar horowitz equation with our extension
!  uses the Johnson correction at short distances (r < rr)
!     if r .ge. rr         use modified t/h triplet equation
!     if r .le. rl         use the Johnson correction
!     in between           use the transition equation
!--------------------------------------------------------------------
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

   dimension E3(3), E1(3)
   parameter(rl=0.95d0, rr=1.15d0)
   parameter(z1=1.d0, z2=2.d0)
! triplet values jun21
   parameter( &
      a1=-0.0298546962d0, a2=-23.9604445036d0, a3=-42.5185569474d0, &
      a4=2.0382390988d0, a5=-11.5214861455d0, a6=1.5309487826d0, &
      C1=-0.4106358351531854d0, C2=-0.0770355790707090d0, &
      C3=0.4303193846943223d0)
!
   E3(3) = 0.d0
   E3(2) = 0.d0
!
   if (r .ge. rr) then
!       modified truhlar/horowitz triplet equation:
      exdr = exp(-a4*r)
      rsq = r*r
      ra6 = r**(-a6)
      E3(1) = a1*(a2 + r + a3*rsq + a5*ra6)*exdr
!       first derivative of triplet curve:
      if (ideriv .gt. 0) then
         ra61 = r**(-a6 - z1)
         E3(2) = a1*exdr*(z1 - a2*a4 + (z2*a3 - a4)*r &
                          - a3*a4*rsq - a5*a6*ra61 - a4*a5*ra6)
      end if
   else
      dr = r - rl
      call vh2opt95_h4bmkp(r, e1, ideriv)
      if (r .le. rl) then
!          Johnson triplet equation:
         E3(1) = E1(1) + c2*dr + c3
         if (ideriv .gt. 0) E3(2) = E1(2) + c2
      else
!          Transition equation:
         E3(1) = E1(1) + c1*dr*dr*dr + c2*dr + c3
         if (ideriv .gt. 0) E3(2) = E1(2) + 3.d0*c1*dr*dr + c2
      end if
   end if
!
   return
end
!
!**********************************************************************
!
subroutine vH2opt95_h4bmkp(r, E, ideriv)
!-----------------------------------------------------------------
! jul09/90 ... super duper speedy version
! self-contained version of schwenke's H2 potential
! all distances in bohrs and all energies in Hartrees
! (1st deriv added on May 2 1989; 2nd deriv on May 28 1989)
!-----------------------------------------------------------------
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   common/c_vh2_rlo_h4bmkp/rsw_h2, rlo_h2, aa_h2, bb_h2, cc_h2
   save/c_vh2_rlo_h4bmkp/
!
   dimension E(3)
!
   parameter( &
      a0=0.03537359271649620d0, a1=2.013977588700072d0, &
      a2=-2.827452449964767d0, a3=2.713257715593500d0, &
      a4=-2.792039234205731d0, a5=2.166542078766724d0, &
      a6=-1.272679684173909d0, a7=0.5630423099212294d0, &
      a8=-0.1879397372273814d0, a9=0.04719891893374140d0, &
      a10=-0.008851622656489644d0, a11=0.001224998776243630d0, &
      a12=-1.227820520228028d-04, a13=8.638783190083473d-06, &
      a14=-4.036967926499151d-07, a15=1.123286608335365d-08, &
      a16=-1.406619156782167d-10)
!
   parameter(r0=3.5284882d0, DD=0.160979391d0, &
             c6=6.499027d0, c8=124.3991d0, c10=3285.828d0)
!
   E(2) = 0.d0
   E(3) = 0.d0
   r2 = r*r
   r3 = r2*r
   r4 = r3*r
!
! Small distance cases (r < rsw_h2): switch from VH2_schwenke to 1/r:
!
!   NOTE: Vh2 and the 1st derivative are continuous everywhere;
!     at r = rsw_h2 = 0.1193892052031558, 2nd derivative is continuous;
!     at r = rlo_h2 = 0.045156498889612436, 2nd derivative is badly
!                                   discontinuous (by a factor of two).
!
   if (r .le. rlo_h2) then
      e(1) = 1.d0/r
      if (ideriv .gt. 0) then
         e(2) = -1.d0/r2
!-noneed;            if ( ideriv .ge. 2 ) e(3) = 2.d0 / r3
      end if
      return
   else if (r .lt. rsw_h2) then
      e(1) = aa_h2/r + bb_h2/r2 + cc_h2/r3
      if (ideriv .gt. 0) then
         e(2) = -aa_h2/r2 - 2.d0*bb_h2/r3 - 3.d0*cc_h2/r4
!-noneed;            if ( ideriv .ge. 2 ) e(3) = 2.d0 * aa_h2 / r3
!-noneed;     $           + 6.d0 * bb_h2 / r4 + 12.d0 * cc_h2 / ( r4 * r )
      end if
      return
   end if
!
!  Standard Schwenke curve:
!
   r5 = r4*r
   r6 = r5*r
   r7 = r6*r
   r8 = r7*r
   r9 = r8*r
   r10 = r9*r
   r11 = r10*r
   r12 = r11*r
   r13 = r12*r
   r14 = r13*r
   r15 = r14*r
   r02 = r0*r0
   r04 = r02*r02
   r06 = r04*r02
   rr2 = r2 + r02
   rr4 = r4 + r04
   rr6 = r6 + r06
   rr25 = rr2*rr2*rr2*rr2*rr2
!     general term:  a(i)*r(i-1),  i=0,16
   alphaR = a0/r + a1 &
            + a2*r + a3*r2 + a4*r3 + a5*r4 + a6*r5 &
            + a7*r6 + a8*r7 + a9*r8 + a10*r9 + a11*r10 &
            + a12*r11 + a13*r12 + a14*r13 + a15*r14 + a16*r15
   exalph = exp(alphaR)
   vsr = DD*(exalph - 1.d0)*(exalph - 1.d0) - DD
   vlr = -C6/rr6 - C8/(rr4*rr4) - C10/rr25
   E(1) = vsr + vlr
!
!  calculate first derivative if required:
!
   if (ideriv .gt. 0) then
!
      rr26 = rr25*rr2
      rr43 = rr4*rr4*rr4
!        general term:  (i-1)*a(i)*r**(i-2)  ,  i=0,16
      dalphR = -a0/r2 + a2 + 2.d0*a3*r &
               + 3.d0*a4*r2 + 4.d0*a5*r3 + 5.d0*a6*r4 &
               + 6.d0*a7*r5 + 7.d0*a8*r6 + 8.d0*a9*r7 &
               + 9.d0*a10*r8 + 10.d0*a11*r9 + 11.d0*a12*r10 &
               + 12.d0*a13*r11 + 13.d0*a14*r12 + 14.d0*a15*r13 &
               + 15.d0*a16*r14
      dvsr = 2.d0*DD*(exalph - 1.d0)*exalph*dalphR
      dvlr = 6.d0*C6*r5/(rr6*rr6) &
             + 8.d0*C8*r3/rr43 + 10.d0*C10*r/rr26
      E(2) = dvsr + dvlr
!
!-noneed;C  calculate second derivative if required:
!-noneed;c
!-noneed;         if(ideriv.ge.2)then
!-noneed;            r10  = r6*r4
!-noneed;            rr27 = rr26*rr2
!-noneed;            rr44 = rr43*rr4
!-noneed;            rr62 = rr6*rr6
!-noneed;            rr63 = rr62*rr6
!-noneed;c     general term: (i-1)*(i-2)*a(i)*r**(i-3),    i=0,16
!-noneed;            ddalph = 2.d0*a0/r3 +  2.d0*a3      +  6.d0*a4 *r
!-noneed;     .           + 12.d0*a5 *r2    + 20.d0*a6 *r3  + 30.d0*a7 *r4
!-noneed;     .           + 42.d0*a8 *r5    + 56.d0*a9 *r6  + 72.d0*a10*r7
!-noneed;     .           + 90.d0*a11*r8    +110.d0*a12*r9  +132.d0*a13*r10
!-noneed;     .           +156.d0*a14*r11   +182.d0*a15*r12 +210.d0*a16*r13
!-noneed;            ddvsr = 2.d0 *DD *exalph
!-noneed;     .           * ( (2.d0*exalph-1.d0)*dalphR*dalphR
!-noneed;     .           + (exalph-1.d0)*ddalph )
!-noneed;            ddvlr =- 72.d0* C6 *r10 / rr63 -96.d0 *C8 *r6 / rr44
!-noneed;     .           -120.d0*C10 *r2  / rr27 +30.d0 *C6 *r4 / rr62
!-noneed;     .           + 24.d0* C8 *r2  / rr43 +10.d0 *C10    / rr26
!-noneed;            E(3) = ddvsr + ddvlr
!-noneed;         endif
!
   end if
!
   return
end
!
!**********************************************************************
!
!-noneed;      subroutine h3lond95_h4bmkp(r,Vlon,dVlon)
!-noneed;c---------------------------------------------------------------------
!-noneed;c version of may 12/90 ... derivatives corrected (0.5 changed to 0.25)
!-noneed;c calculates the h3 london terms and derivatives
!-noneed;c modified oct 7/89 to include eps**2 term which rounds off the
!-noneed;c cusp in the h3 potential which occurs at equilateral triangle
!-noneed;c configurations
!-noneed;c---------------------------------------------------------------------
!-noneed;       implicit double precision (a-h,o-z))
!-noneed;      double precision Q(3),J(3),Jt
!-noneed;      dimension r(3),e1(3),e3(3),esing(3),etrip(3),dVlon(3)
!-noneed;      dimension de1(3),de3(3)
!-noneed;c     dimension dJ1(3),dJ2(3),dJ3(3)
!-noneed;      parameter( half=0.5d0, two=2.d0 , eps2=1.d-12 )
!-noneed;      ipr = 0
!-noneed;      do i=1,3
!-noneed;         call vh2opt95_h4bmkp(r(i),esing,2)
!-noneed;        e1(i) = esing(1)
!-noneed;       de1(i) = esing(2)
!-noneed;       call triplet95_h4bmkp(r(i),etrip,2)
!-noneed;       if(ipr.gt.0) write(7,7400) i,esing,etrip
!-noneed;        e3(i) = etrip(1)
!-noneed;       de3(i) = etrip(2)
!-noneed;         q(i)  = half*(E1(i) + E3(i))
!-noneed;         j(i)  = half*(E1(i) - E3(i))
!-noneed;      enddo
!-noneed;      sumQ  =   q(1) + q(2) + q(3)
!-noneed;      sumJ  =   abs( j(2)-j(1) )**2
!-noneed;     .        + abs( j(3)-j(2) )**2
!-noneed;     .        + abs( j(3)-j(1) )**2
!-noneed;      Jt     = half*sumJ + eps2
!-noneed;      rootJt = sqrt(Jt)
!-noneed;      Vlon   = sumQ - rootJt
!-noneed;      if(ipr.gt.0) then
!-noneed;       write(7,7410) sumq,sumj
!-noneed;         write(7,7420) vlon,rootJt
!-noneed;      endif
!-noneed;c  calculate the derivatives with respect to r(i):
!-noneed;      dVlon(1) = half*(dE1(1)+de3(1))
!-noneed;     .         - 0.25d0*(two*j(1)-j(2)-j(3))*(dE1(1)-de3(1))/rootJt
!-noneed;      dVlon(2) = half*(dE1(2)+de3(2))
!-noneed;     .         - 0.25d0*(two*j(2)-j(3)-j(1))*(dE1(2)-de3(2))/rootJt
!-noneed;      dVlon(3) = half*(dE1(3)+de3(3))
!-noneed;     .         - 0.25d0*(two*j(3)-j(1)-j(2))*(dE1(3)-de3(3))/rootJt
!-noneed;      if(ipr.gt.0) then
!-noneed;         write(7,7000) r,e1,e3,vlon
!-noneed;       write(7,7100) q,j
!-noneed;       write(7,7200) dVlon
!-noneed;      endif
!-noneed; 7000 format('             r = ',3(1x,f12.6),/,
!-noneed;     .       '            E1 = ',3(1x,f12.8),/,
!-noneed;     .       '            E3 = ',3(1x,f12.8),/,
!-noneed;     .       '          vlon = ',1x,f12.8)
!-noneed; 7100 format(13x,'q = ',3(1x,e12.6),/,13x,'j = ',3(1x,e12.6))
!-noneed; 7200 format('         dVlon = ',3(1x,g12.6))
!-noneed; 7400 format('from subr.london: ',/,
!-noneed;     .       '  using r',i1,':   Esinglet=',3(1x,f12.8),/,
!-noneed;     .       '         ',1x,'    Etriplet=',3(1x,f12.8))
!-noneed; 7410 format('         sumQ = ',g12.6,'        sumJ = ',g12.6)
!-noneed; 7420 format('         Vlon = ',f12.8,'      rootJt = ',g12.6)
!-noneed;      return
!-noneed;      end
!
!**********************************************************************
!
subroutine Vascal95_h4bmkp(rpass, Vas, dVas, ideriv, A, dA)
!------------------------------------------------------------------
! version of apr12/95  950308 values
! version of oct11/90  fit632c.out values
! version of oct5/90   surf626 values
!  calculate the asymmetric correction term and its derivatives
!  see equations [14] to [16] of truhlar/horowitz 1978 paper
!------------------------------------------------------------------
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

   dimension dVas(3), dA(3), dS(3), rpass(3)
! vasym values jun21
   parameter( &
      aa1=0.3788951192d-02, aa2=0.1478100901d-02, &
      aa3=-.1848513849d-03, aa4=0.9230803609d-05, &
      aa5=-.1293180255d-06, aa6=0.5237179303d+00, &
      aa7=-.1112326215d-02)
!
   r1 = rpass(1)
   r2 = rpass(2)
   r3 = rpass(3)
   R = r1 + r2 + r3
   Rsq = R*R
   Rcu = Rsq*R
!
!   calculate the Vas term first (eq.14 of truhlar/horowitz)
!
   A = (r1 - r2)*(r2 - r3)*(r3 - r1)
   if (A .lt. 0.d0) then
      isgn_A = -1
      A = -A
   else
      isgn_A = 1
   end if
!
   A2 = A*A
   A3 = A2*A
   A4 = A3*A
   A5 = A4*A
   exp1 = exp(-aa1*Rcu)
   exp6 = exp(-aa6*R)
   S = aa2*A2 + aa3*A3 + aa4*A4 + aa5*A5
   Vas = S*exp1 + aa7*A2*exp6/R
!
   if (ideriv .gt. 0) then
   if (isgn_A .gt. 0) then
      dA(1) = (r2 + r3 - 2.d0*r1)*(r2 - r3)
      dA(2) = (r3 + r1 - 2.d0*r2)*(r3 - r1)
      dA(3) = (r1 + r2 - 2.d0*r3)*(r1 - r2)
   else
      dA(1) = (2.d0*r1 - r2 - r3)*(r2 - r3)
      dA(2) = (2.d0*r2 - r3 - r1)*(r3 - r1)
      dA(3) = (2.d0*r3 - r1 - r2)*(r1 - r2)
   end if
   dS(1) = (2.d0*aa2*A + 3.d0*aa3*A2 &
            + 4.d0*aa4*A3 + 5.d0*aa5*A4)*dA(1)
   dVas(1) = -3.d0*aa1*Rsq*S*exp1 + dS(1)*exp1 &
             - aa7*A2*exp6/Rsq + 2.d0*aa7*A*dA(1)*exp6/R &
             - aa6*aa7*A2*exp6/R
   dS(2) = (2.d0*aa2*A + 3.d0*aa3*A2 &
            + 4.d0*aa4*A3 + 5.d0*aa5*A4)*dA(2)
   dVas(2) = -3.d0*aa1*Rsq*S*exp1 + dS(2)*exp1 &
             - aa7*A2*exp6/Rsq + 2.d0*aa7*A*dA(2)*exp6/R &
             - aa6*aa7*A2*exp6/R
   dS(3) = (2.d0*aa2*A + 3.d0*aa3*A2 &
            + 4.d0*aa4*A3 + 5.d0*aa5*A4)*dA(3)
   dVas(3) = -3.d0*aa1*Rsq*S*exp1 + dS(3)*exp1 &
             - aa7*A2*exp6/Rsq + 2.d0*aa7*A*dA(3)*exp6/R &
             - aa6*aa7*A2*exp6/R
   end if
   return
end
!
!**********************************************************************
!
subroutine compac95_h4bmkp(r, icompc, T, dT)
!---------------------------------------------------------------
!  version of may 14/90 ... calculates T and dT values also
!  decide whether or not this particular geometry is compact,
!  that is, are any of the three distances smaller than the
!  rr value from the johnson correction.
!-----------------------------------------------------------------
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

   dimension r(3), T(3), dT(3)
   parameter(rr=1.15d0, rp=1.25d0)
!
! calculate the T(i) values, for any compact geometries:
!
   icompc = 0
!
   do i = 1, 3
   if (r(i) .lt. rr) then
      icompc = icompc + 1
      top = rr - r(i)
      bot = rp - r(i)
      top2 = top*top
      top3 = top2*top
      bot2 = bot*bot
      T(i) = top3/bot
      dT(i) = (top/bot - 3.d0)*top2/bot
   else
      T(i) = 0.d0
      dT(i) = 0.d0
   end if
   end do
!
   return
end
!
!**********************************************************************
!
subroutine casym95_h4bmkp(r, Cas, dCas, ideriv, T, dT, A, dA)
!---------------------------------------------------------------
! apr12/95 ... surface950308 parameters added
! version of sept14/90
!  the compact asymmetric correction term and derivatives
!---------------------------------------------------------------
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(3), T(3)
   dimension dPr(3), dT(3), dsumT(3), dterm1(3), determ(3), dseries(3), &
      dCas(3), dA(3)
!
! Casym values jun21
   parameter( &
      u1=0.2210243144d+00, u2=0.4367417579d+00, u3=0.6994985432d-02, &
      u4=0.1491096501d+01, u5=0.1602896673d+01, u6=-.2821747323d+01, &
      u7=0.4948310833d+00, u8=-.3540394679d-01, u9=-.3305809954d+01, &
      u10=0.3644382172d+01, u11=-.9997570970d+00, u12=0.7989919534d-01, &
      u13=-.1075807322d-02)
!
   cas = 0.d0
   dCas(1) = 0.d0
   dCas(2) = 0.d0
   dCas(3) = 0.d0
   A2 = A*A
   sumT = T(1) + T(2) + T(3)
   Sr = r(1) + r(2) + r(3)
   Pr = r(1)*r(2)*r(3)
   Sr2 = Sr*Sr
   Sr3 = Sr2*Sr
   Pr2 = Pr*Pr
   Pr3 = Pr2*Pr
!    write out the series explicitly:
   series = 1.d0 + u4/Pr2 + u5/Pr + u6 + u7*Pr + u8*Pr2 &
            + A*(u9/Pr2 + u10/Pr + u11 + u12*Pr + u13*Pr2)
   term1 = u1/Pr**u2
   eterm = exp(-u3*Sr3)
   Cas = sumt*a2*term1*series*eterm
   if (ideriv .gt. 0) then
      dPr(1) = r(2)*r(3)
      dPr(2) = r(3)*r(1)
      dPr(3) = r(1)*r(2)
      do i = 1, 3
         dsumt(i) = dt(i)
         dterm1(i) = -1.d0*u1*u2*Pr**(-u2 - 1.d0)*dPr(i)
         determ(i) = eterm*(-3.d0*u3*Sr2)
         dseries(i) = &
            dPr(i)*(-2.d0*u4/Pr3 - u5/Pr2 + u7 + 2.d0*u8*Pr) &
            + dA(i)*(u9/Pr2 + u10/Pr + u11 + u12*Pr + u13*Pr2) &
            + A*dPr(i)*(-2.d0*u9/Pr3 - u10/Pr2 + u12 + 2.d0*u13*Pr)
         dCas(i) = &
            dsumt(i)*a2*term1*series*eterm &
            + 2.d0*a*dA(i)*sumt*term1*series*eterm &
            + dterm1(i)*sumt*a2*series*eterm &
            + dseries(i)*sumt*a2*term1*eterm &
            + determ(i)*sumt*a2*term1*series
      end do
   end if
   return
end
!
!**********************************************************************
!
subroutine csym95_h4bmkp(r, Cal, dCal)
!---------------------------------------------------------------
! version of apr12/95  950308 values
! version of oct12/90  surf636 values
!  calculate the 'compact all' correction term and derivatives
!  a correction term (added sept 11/89), which adds a small
!  correction to all compact geometries
!---------------------------------------------------------------
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(3), dCal(3), G(3)
   dimension dG(3), sumv(3)
!
   parameter(rr=1.15d0, rp=1.25d0)
!
! Csym values jun21
   parameter( &
      v1=-.2071708868d+00, v2=-.5672350377d+00, v3=0.9058780367d-02)
!
   Sr = r(1) + r(2) + r(3)
   Sr2 = Sr*Sr
   Sr3 = Sr*Sr2
   exp3 = exp(-v3*Sr3)
   dexp3 = -3.d0*v3*Sr2*exp3
   do i = 1, 3
      ri = r(i)
      if (ri .lt. rr) then
         rrri = rr - ri
         rrri2 = rrri*rrri
         rrri3 = rrri*rrri2
         rpri = rp - ri
         rpri2 = rpri*rpri
         sumv(i) = v1 + v1*v2*ri
         G(i) = (rrri3/rpri)*sumv(i)
         dG(i) = (rrri3/rpri2)*sumv(i) &
                 - 3.d0*(rrri2/rpri)*sumv(i) &
                 + (rrri3/rpri)*v1*v2
      else
         G(i) = 0.d0
         dG(i) = 0.d0
      end if
   end do
   sumG = G(1) + G(2) + G(3)
   cal = sumG*exp3
   dCal(1) = dG(1)*exp3 + sumG*dexp3
   dCal(2) = dG(2)*exp3 + sumG*dexp3
   dCal(3) = dG(3)*exp3 + sumG*dexp3
   return
end
!
!**********************************************************************
!
subroutine vbcb95_h4bmkp(rpass, icompc, T, dT, ideriv, &
                         Vbnda, Vbndb, dVbnda, dVbndb, &
                         Cbnda, Cbndb, dCbnda, dCbndb)
!---------------------------------------------------------------
! apr12/95 surface950308 values added
! jul24/91 fit705 cbend values added
! feb27/91 modified to match equation in h3 paper more closely
! nov 4/90 Vbend coefficients now a,g  Cbend coeff's still c,d
! in this version, the derivatives are always calculated
! subroutines vbend, cbend and B1ab all combined into this one
! module in order to improve efficiency by not having to pass
! around B1a,B1b,B2,B3a,B3b functions and derivatives
!  B1a = 1 - sum of [ P1(cos(theta(i))) ]
!  B1b = 1 - sum of [ P3(cos(theta(i))) ]
!--------------------------------------------------------------
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
!     arrays required for B1,B2,B3 calculations:
   dimension rpass(3), dB1a(3), dB1b(3), dB2(3), dB3(6)
!     dimension th(3)
!     arrays required for Vbend calculations:
   dimension dVbnda(3), dVbndb(3)
   dimension dVb1(3), dvb2(3), dvb3(3), dvb4(3), dvb5(3)
!     dimension dVbnd(3)
   dimension Vbnd(2)
   dimension vb(2, 25)
!     arrays required for Cbend calculations:
   dimension Cbnds(2), dCbnds(2, 3)
   dimension dCb1(3), dCb2(3), dCb3(3), dCb4(3), dCb5(3), dCb6(3), &
      dCb7(3), dCb8(3)
   dimension T(3), dP(3)
   dimension dT(3), dCbnda(3), dCbndb(3)
   dimension cb(2, 25)
   parameter(z58=0.625d0, z38=0.375d0)
!     parameters required for Vbend calculations:
   parameter(beta1=0.52d0, beta2=0.052d0, beta3=0.79d0)
! vbenda terms jun21         expterms= 0.52d0, 0.052d0, 0.79d0
   parameter( &
      a11=-.1838073394d+03, a12=0.1334593242d+02, a13=-.2358129537d+00, &
      a21=-.4668193478d+01, a22=0.7197506670d+01, a23=0.2162004275d+02, &
      a24=0.2106294028d+02, a31=0.4242962586d+01, a32=0.4453505045d+01, &
      a41=-.1456918088d+00, a42=-.1692657366d-01, a43=0.1279520698d+01, &
      a44=-.4898940075d+00, a51=0.1742295219d+03, a52=0.3142175348d+02, &
      a53=0.5152903406d+01)
! vbendb terms jun21
   parameter( &
      g11=-.4765732725d+02, g12=0.3648933563d+01, g13=-.7141145244d-01, &
      g21=0.1002349176d-01, g22=0.9989856329d-02, g23=-.4161953634d-02, &
      g24=0.9075807910d-03, g31=-.2693628729d+00, g32=-.1399065763d-01, &
      g41=-.1417634346d-01, g42=-.4870024792d-03, g43=0.1312231847d+00, &
      g44=-.4409850519d-01, g51=0.5382970863d+02, g52=0.4587102824d+01, &
      g53=0.1768550515d+01)
!
! cbenda terms jun21/95
   parameter( &
      c11=0.1860299931d+04, c12=-.6134458037d+03, c13=0.7337207161d+02, &
      c14=-.2676717625d+04, c15=0.1344099415d+04, c21=0.1538913137d+03, &
      c22=0.4348007369d+02, c23=0.1719720677d+03, c24=0.2115963042d+03, &
      c31=-.7026089414d+02, c32=-.1300938992d+03, c41=0.1310273564d+01, &
      c42=-.6175149574d+00, c43=-.2679089358d+02, c44=0.5577477171d+01, &
      c51=-.3543353539d+04, c52=-.3740709591d+03, c53=0.7979303144d+02, &
      c61=-.1104230585d+04, c62=0.4603572025d+04, c63=-.5593496634d+04, &
      c71=-.1069406434d+02, c72=0.1021807153d+01, c73=0.6669828341d-01, &
      c74=0.4168542348d+02, c75=0.1751608567d+02, c81=0.9486883238d+02, &
      c82=-.1519334221d+02, c83=0.4024697252d+04, c84=-.2225159395d+02)
!           (last 4 parameters renumbered)
! cbendb terms jun21/95
   parameter( &
      d11=0.4203543357d+03, d12=-.4922474096d+02, d13=0.3362942544d+00, &
      d14=-.3827423082d+03, d15=0.1746726001d+03, d21=0.1699995737d-01, &
      d22=0.1513036778d-01, d23=0.2659119354d-01, d24=-.5760387483d-02, &
      d31=0.1020622621d+02, d32=0.1050536271d-01, d41=0.6836172780d+00, &
      d42=-.1627858240d+00, d43=-.6925485045d+01, d44=0.1632567385d+01, &
      d51=0.1083595009d+04, d52=0.4641431791d+01, d53=-.8233144461d+00, &
      d61=-.6157225942d+02, d62=0.3094361471d+03, d63=-.3299631143d+03, &
      d71=0.8866227120d+01, d72=-.1382126854d+01, d73=0.7620770145d-01, &
      d74=-.5145757859d+02, d75=0.2046097265d+01, d81=0.2540775558d+01, &
      d82=-.4889246569d+00, d83=-.1127439280d+04, d84=-.2269932295d+01)
!           (last 4 parameters renumbered)
!
!   feb27/91 new c51 = c51+c83;   new d51 = d51+d83
!
   cx1 = c51 + c83
   dx1 = d51 + d83
   r1 = rpass(1)
   r2 = rpass(2)
   r3 = rpass(3)
   t1 = r1*r1 - r2*r2 - r3*r3
   t2 = r2*r2 - r3*r3 - r1*r1
   t3 = r3*r3 - r1*r1 - r2*r2
!  calculate the cosines of the three internal angles:
   c1 = t1/(-2.d0*r2*r3)
   c2 = t2/(-2.d0*r3*r1)
   c3 = t3/(-2.d0*r1*r2)
   sum = c1 + c2 + c3
   B1a = 1.d0 - sum
   c1cube = c1*c1*c1
   c2cube = c2*c2*c2
   c3cube = c3*c3*c3
   cos3t1 = 4.d0*c1cube - 3.d0*c1
   cos3t2 = 4.d0*c2cube - 3.d0*c2
   cos3t3 = 4.d0*c3cube - 3.d0*c3
   sumb = cos3t1 + cos3t2 + cos3t3
   B1b = 1.d0 - (z58*sumb + z38*sum)
!  calculate derivatives if desired
   if (ideriv .gt. 0) then
      dc1dr1 = -r1/(r2*r3)
      dc2dr2 = -r2/(r1*r3)
      dc3dr3 = -r3/(r1*r2)
      dc1dr2 = (t1/(r2*r2) + 2.d0)/(2.d0*r3)
      dc1dr3 = (t1/(r3*r3) + 2.d0)/(2.d0*r2)
      dc2dr1 = (t2/(r1*r1) + 2.d0)/(2.d0*r3)
      dc2dr3 = (t2/(r3*r3) + 2.d0)/(2.d0*r1)
      dc3dr1 = (t3/(r1*r1) + 2.d0)/(2.d0*r2)
      dc3dr2 = (t3/(r2*r2) + 2.d0)/(2.d0*r1)
      db1a(1) = -1.d0*(dc1dr1 + dc2dr1 + dc3dr1)
      db1a(2) = -1.d0*(dc1dr2 + dc2dr2 + dc3dr2)
      db1a(3) = -1.d0*(dc1dr3 + dc2dr3 + dc3dr3)
      d1 = 12.d0*c1*c1 - 3.d0
      d2 = 12.d0*c2*c2 - 3.d0
      d3 = 12.d0*c3*c3 - 3.d0
      db1b(1) = -z58*(d1*dc1dr1 + d2*dc2dr1 + d3*dc3dr1) &
                - z38*(dc1dr1 + dc2dr1 + dc3dr1)
      db1b(2) = -z58*(d1*dc1dr2 + d2*dc2dr2 + d3*dc3dr2) &
                - z38*(dc1dr2 + dc2dr2 + dc3dr2)
      db1b(3) = -z58*(d1*dc1dr3 + d2*dc2dr3 + d3*dc3dr3) &
                - z38*(dc1dr3 + dc2dr3 + dc3dr3)
   end if
!
!  calculate the quantities used by both Vbend and Cbend:
!
   R = r1 + r2 + r3
   Rsq = R*R
   B2 = 1.d0/r1 + 1.d0/r2 + 1.d0/r3
   B3 = (r2 - r1)*(r2 - r1) + (r3 - r2)*(r3 - r2) + (r1 - r3)*(r1 - r3)
   eps2 = 1.d-12
   B3b = sqrt(B3 + eps2)
   exp1 = exp(-beta1*R)
   exp2 = exp(-beta2*Rsq)
   exp7 = exp(-beta3*R)
   dexp1 = -beta1*exp1
   dexp2 = -2.d0*beta2*R*exp2
   dexp7 = -beta3*exp7
!   do the vbnda calculations:
   B1 = B1a
   B12 = B1*B1
   B13 = B12*B1
   B14 = B13*B1
   B15 = B14*B1
   asum = a11 + a12*R + a13*Rsq
   bsum = a21*B12 + a22*B13 + a23*B14 + a24*B15
   csum = a31*B1*exp1 + a32*B12*exp2
   dsum1 = a41*exp1 + a42*exp2
   dsum2 = a43*exp1 + a44*exp2
   fsum = a51 + a52*R + a53*Rsq
   vb(1, 1) = B1*asum*exp1
   vb(1, 2) = bsum*exp2
   vb(1, 3) = B2*csum
   vb(1, 4) = B1*B3*dsum1 + B1*B3b*dsum2
   vb(1, 5) = B1*fsum*exp7
   Vbnd(1) = vb(1, 1) + vb(1, 2) + vb(1, 3) + vb(1, 4) + vb(1, 5)
!
   if (ideriv .gt. 0) then
!
      dB2(1) = -1.d0/(r1*r1)
      dB2(2) = -1.d0/(r2*r2)
      dB2(3) = -1.d0/(r3*r3)
      dB3(1) = 4.d0*r1 - 2.d0*r2 - 2.d0*r3
      dB3(2) = 4.d0*r2 - 2.d0*r3 - 2.d0*r1
      dB3(3) = 4.d0*r3 - 2.d0*r1 - 2.d0*r2
      dB3(4) = 0.5d0*dB3(1)/B3b
      dB3(5) = 0.5d0*dB3(2)/B3b
      dB3(6) = 0.5d0*dB3(3)/B3b
      dasum = a12 + 2.d0*a13*R
      dbsum = 2.d0*a21*B1 + 3.d0*a22*B12 &
              + 4.d0*a23*B13 + 5.d0*a24*B14
      ddsum1 = a41*dexp1 + a42*dexp2
      ddsum2 = a43*dexp1 + a44*dexp2
      dfsum = a52 + 2.d0*a53*R
      dVb1(1) = dB1a(1)*asum*exp1 + B1*dasum*exp1 + B1*asum*dexp1
      dVb1(2) = dB1a(2)*asum*exp1 + B1*dasum*exp1 + B1*asum*dexp1
      dVb1(3) = dB1a(3)*asum*exp1 + B1*dasum*exp1 + B1*asum*dexp1
!
      dVb2(1) = dbsum*dB1a(1)*exp2 + bsum*dexp2
      dVb2(2) = dbsum*dB1a(2)*exp2 + bsum*dexp2
      dVb2(3) = dbsum*dB1a(3)*exp2 + bsum*dexp2
!
! calculate the Vb3 derivatives:
      dVb3(1) = dB2(1)*csum + B2*(a31*dB1a(1)*exp1 + a31*B1*dexp1 &
                                  + 2.d0*a32*B1*dB1a(1)*exp2 + a32*B12*dexp2)
      dVb3(2) = dB2(2)*csum + B2*(a31*dB1a(2)*exp1 + a31*B1*dexp1 &
                                  + 2.d0*a32*B1*dB1a(2)*exp2 + a32*B12*dexp2)
      dVb3(3) = dB2(3)*csum + B2*(a31*dB1a(3)*exp1 + a31*B1*dexp1 &
                                  + 2.d0*a32*B1*dB1a(3)*exp2 + a32*B12*dexp2)
!
! calculate the Vb4 derivatives (may 27/90):
      dVb4(1) = dB1a(1)*B3*dsum1 + B1*dB3(1)*dsum1 + B1*B3*ddsum1 &
                + dB1a(1)*B3b*dsum2 + B1*dB3(4)*dsum2 + B1*B3b*ddsum2
      dVb4(2) = dB1a(2)*B3*dsum1 + B1*dB3(2)*dsum1 + B1*B3*ddsum1 &
                + dB1a(2)*B3b*dsum2 + B1*dB3(5)*dsum2 + B1*B3b*ddsum2
      dVb4(3) = dB1a(3)*B3*dsum1 + B1*dB3(3)*dsum1 + B1*B3*ddsum1 &
                + dB1a(3)*B3b*dsum2 + B1*dB3(6)*dsum2 + B1*B3b*ddsum2
!
! calculate the Vb5 derivatives:
      dVb5(1) = dB1a(1)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
      dVb5(2) = dB1a(2)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
      dVb5(3) = dB1a(3)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
!
! calculate the overall derivatives:
      dVbnda(1) = dvb1(1) + dvb2(1) + dvb3(1) + dvb4(1) + dvb5(1)
      dVbnda(2) = dvb1(2) + dvb2(2) + dvb3(2) + dvb4(2) + dvb5(2)
      dVbnda(3) = dvb1(3) + dvb2(3) + dvb3(3) + dvb4(3) + dvb5(3)
!
   end if
!
   Vbnda = Vbnd(1)
!
!---- now do the vbendb calculations --------
   B1 = B1b
   B12 = B1*B1
   B13 = B12*B1
   B14 = B13*B1
   B15 = B14*B1
   asum = g11 + g12*R + g13*Rsq
   bsum = g21*B12 + g22*B13 + g23*B14 + g24*B15
   csum = g31*B1*exp1 + g32*B12*exp2
   dsum1 = g41*exp1 + g42*exp2
   dsum2 = g43*exp1 + g44*exp2
   fsum = g51 + g52*R + g53*Rsq
   vb(2, 1) = B1*asum*exp1
   vb(2, 2) = bsum*exp2
   vb(2, 3) = B2*csum
   vb(2, 4) = B1*B3*dsum1 + B1*B3b*dsum2
   vb(2, 5) = B1*fsum*exp7
   Vbnd(2) = vb(2, 1) + vb(2, 2) + vb(2, 3) + vb(2, 4) + vb(2, 5)
!
   if (ideriv .gt. 0) then
!
      dasum = g12 + 2.d0*g13*R
      dbsum = 2.d0*g21*B1 + 3.d0*g22*B12 &
              + 4.d0*g23*B13 + 5.d0*g24*B14
      ddsum1 = g41*dexp1 + g42*dexp2
      ddsum2 = g43*dexp1 + g44*dexp2
      dfsum = g52 + 2.d0*g53*R
      dVb1(1) = dB1b(1)*asum*exp1 + B1*dasum*exp1 + B1*asum*dexp1
      dVb1(2) = dB1b(2)*asum*exp1 + B1*dasum*exp1 + B1*asum*dexp1
      dVb1(3) = dB1b(3)*asum*exp1 + B1*dasum*exp1 + B1*asum*dexp1
!
      dVb2(1) = dbsum*dB1b(1)*exp2 + bsum*dexp2
      dVb2(2) = dbsum*dB1b(2)*exp2 + bsum*dexp2
      dVb2(3) = dbsum*dB1b(3)*exp2 + bsum*dexp2
!
! calculate the Vb3 derivatives:
      dVb3(1) = dB2(1)*csum + B2*(g31*dB1b(1)*exp1 + g31*B1*dexp1 &
                                  + 2.d0*g32*B1*dB1b(1)*exp2 + g32*B12*dexp2)
      dVb3(2) = dB2(2)*csum + B2*(g31*dB1b(2)*exp1 + g31*B1*dexp1 &
                                  + 2.d0*g32*B1*dB1b(2)*exp2 + g32*B12*dexp2)
      dVb3(3) = dB2(3)*csum + B2*(g31*dB1b(3)*exp1 + g31*B1*dexp1 &
                                  + 2.d0*g32*B1*dB1b(3)*exp2 + g32*B12*dexp2)
!
! calculate the Vb4 derivatives (may 27/90):
      dVb4(1) = dB1b(1)*B3*dsum1 + B1*dB3(1)*dsum1 + B1*B3*ddsum1 &
                + dB1b(1)*B3b*dsum2 + B1*dB3(4)*dsum2 + B1*B3b*ddsum2
      dVb4(2) = dB1b(2)*B3*dsum1 + B1*dB3(2)*dsum1 + B1*B3*ddsum1 &
                + dB1b(2)*B3b*dsum2 + B1*dB3(5)*dsum2 + B1*B3b*ddsum2
      dVb4(3) = dB1b(3)*B3*dsum1 + B1*dB3(3)*dsum1 + B1*B3*ddsum1 &
                + dB1b(3)*B3b*dsum2 + B1*dB3(6)*dsum2 + B1*B3b*ddsum2
!
! calculate the Vb5 derivatives:
      dVb5(1) = dB1b(1)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
      dVb5(2) = dB1b(2)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
      dVb5(3) = dB1b(3)*fsum*exp7 + B1*dfsum*exp7 + B1*fsum*dexp7
!
! calculate the overall derivatives:
      dVbndb(1) = dvb1(1) + dvb2(1) + dvb3(1) + dvb4(1) + dvb5(1)
      dVbndb(2) = dvb1(2) + dvb2(2) + dvb3(2) + dvb4(2) + dvb5(2)
      dVbndb(3) = dvb1(3) + dvb2(3) + dvb3(3) + dvb4(3) + dvb5(3)
!
   end if
!
   Vbndb = Vbnd(2)
!
!  now calculate the Cbend correction terms:
!  Cbnda uses the B1a formula and Cbndb uses the B1b formula
!
   if (icompc .eq. 0) return
!
   sumT = T(1) + T(2) + T(3)
   Rcu = Rsq*R
   P = r1*r2*r3
   Psq = P*P
   Pcu = Psq*P
   dP(1) = r2*r3
   dP(2) = r3*r1
   dP(3) = r1*r2
   Cbnds(1) = 0.d0
   Cbnds(2) = 0.d0
   dCbnda(1) = 0.d0
   dCbnda(2) = 0.d0
   dCbnda(3) = 0.d0
   dCbndb(1) = 0.d0
   dCbndb(2) = 0.d0
   dCbndb(3) = 0.d0
!  calculate exponentials and derivatives (common to cbnda and cbndb):
   exp7 = exp(-beta2*Rcu)
   dexp7 = -3.d0*beta2*Rsq*exp7
! ----- calculate the Cbnda correction term ----------------
   B1 = B1a
   B12 = B1*B1
   B13 = B12*B1
   B14 = B13*B1
   B15 = B14*B1
   asum = c11 + c12*R + c13*Rsq + c14/R + c15/Rsq
   bsum = c21*B12 + c22*B13 + c23*B14 + c24*B15
   csum = c31*B1*exp1 + c32*B12*exp2
   dsum1 = c41*exp1 + c42*exp2
   dsum2 = c43*exp1 + c44*exp2
   fsum = cx1 + c52*R + c53*Rsq
   gsum = c61 + c62/R + c63/Rsq
   aasum = c71 + c72*P + c73*Psq + c74/P + c75/Psq
   ffsum = c81*P + c82*Psq + c84/Psq
   dasum = c12 + 2.d0*c13*R - c14/Rsq - 2.d0*c15/Rcu
   dbsum = 2.d0*c21*B1 + 3.d0*c22*B12 + 4.d0*c23*B13 + 5.d0*c24*B14
   ddsum1 = c41*dexp1 + c42*dexp2
   ddsum2 = c43*dexp1 + c44*dexp2
   dfsum = c52 + 2.d0*c53*R
   dgsum = -c62/Rsq - 2.d0*c63/Rcu
   daasum = c72 + 2.d0*c73*P - c74/Psq - 2.d0*c75/Pcu
   dffsum = c81 + 2.d0*c82*P - 2.d0*c84/Pcu
   cb(1, 1) = B1*asum*exp1/P
   cb(1, 2) = bsum*exp2
   cb(1, 3) = B2*csum
   cb(1, 4) = B1*B3*dsum1 + B1*B3b*dsum2
   cb(1, 5) = B1*fsum*exp7/P
   cb(1, 6) = B1*gsum*exp7
   cb(1, 7) = B1*aasum*exp2
   cb(1, 8) = B1*ffsum*exp7
!
   Cbnds(1) = cb(1, 1) + cb(1, 2) + cb(1, 3) + cb(1, 4) + cb(1, 5) &
              + cb(1, 6) + cb(1, 7) + cb(1, 8)
!
!       calculate the derivatives:
!
   if (ideriv .gt. 0) then
!
      dCb1(1) = dB1a(1)*asum*exp1/P + B1*dasum*exp1/P &
                + B1*asum*dexp1/P - B1*asum*exp1*dP(1)/Psq
      dCb1(2) = dB1a(2)*asum*exp1/P + B1*dasum*exp1/P &
                + B1*asum*dexp1/P - B1*asum*exp1*dP(2)/Psq
      dCb1(3) = dB1a(3)*asum*exp1/P + B1*dasum*exp1/P &
                + B1*asum*dexp1/P - B1*asum*exp1*dP(3)/Psq
!
      dCb2(1) = dbsum*dB1a(1)*exp2 + bsum*dexp2
      dCb2(2) = dbsum*dB1a(2)*exp2 + bsum*dexp2
      dCb2(3) = dbsum*dB1a(3)*exp2 + bsum*dexp2
!
      dCb3(1) = dB2(1)*csum + B2*(c31*dB1a(1)*exp1 + c31*B1*dexp1 &
                                  + c32*2.d0*B1*dB1a(1)*exp2 + c32*B12*dexp2)
      dCb3(2) = dB2(2)*csum + B2*(c31*dB1a(2)*exp1 + c31*B1*dexp1 &
                                  + c32*2.d0*B1*dB1a(2)*exp2 + c32*B12*dexp2)
      dCb3(3) = dB2(3)*csum + B2*(c31*dB1a(3)*exp1 + c31*B1*dexp1 &
                                  + c32*2.d0*B1*dB1a(3)*exp2 + c32*B12*dexp2)
!          dCb4 equations may 27/90:
      dCb4(1) = dB1a(1)*B3*dsum1 + B1*dB3(1)*dsum1 + B1*B3*ddsum1 &
                + dB1a(1)*B3b*dsum2 + B1*dB3(4)*dsum2 + B1*B3b*ddsum2
      dCb4(2) = dB1a(2)*B3*dsum1 + B1*dB3(2)*dsum1 + B1*B3*ddsum1 &
                + dB1a(2)*B3b*dsum2 + B1*dB3(5)*dsum2 + B1*B3b*ddsum2
      dCb4(3) = dB1a(3)*B3*dsum1 + B1*dB3(3)*dsum1 + B1*B3*ddsum1 &
                + dB1a(3)*B3b*dsum2 + B1*dB3(6)*dsum2 + B1*B3b*ddsum2
!          dCb5 equations may 27/90: (corrected jun01)
      dCb5(1) = dB1a(1)*fsum*exp7/P + B1*dfsum*exp7/P &
                + B1*fsum*dexp7/P - B1*fsum*exp7*dP(1)/Psq
      dCb5(2) = dB1a(2)*fsum*exp7/P + B1*dfsum*exp7/P &
                + B1*fsum*dexp7/P - B1*Fsum*exp7*dP(2)/Psq
      dCb5(3) = dB1a(3)*fsum*exp7/P + B1*dfsum*exp7/P &
                + B1*fsum*dexp7/P - B1*Fsum*exp7*dP(3)/Psq
!
      dCb6(1) = dB1a(1)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
      dCb6(2) = dB1a(2)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
      dCb6(3) = dB1a(3)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
!
      dCb7(1) = dB1a(1)*aasum*exp2 + B1*daasum*dP(1)*exp2 &
                + B1*aasum*dexp2
      dCb7(2) = dB1a(2)*aasum*exp2 + B1*daasum*dP(2)*exp2 &
                + B1*aasum*dexp2
      dCb7(3) = dB1a(3)*aasum*exp2 + B1*daasum*dP(3)*exp2 &
                + B1*aasum*dexp2
      dCb8(1) = dB1a(1)*ffsum*exp7 + B1*dffsum*dP(1)*exp7 &
                + B1*ffsum*dexp7
      dCb8(2) = dB1a(2)*ffsum*exp7 + B1*dffsum*dP(2)*exp7 &
                + B1*ffsum*dexp7
      dCb8(3) = dB1a(3)*ffsum*exp7 + B1*dffsum*dP(3)*exp7 &
                + B1*ffsum*dexp7
      dCbnds(1, 1) = dCb1(1) + dCb2(1) + dCb3(1) + dCb4(1) + dCb5(1) &
                     + dCb6(1) + dCb7(1) + dCb8(1)
      dCbnds(1, 2) = dCb1(2) + dCb2(2) + dCb3(2) + dCb4(2) + dCb5(2) &
                     + dCb6(2) + dCb7(2) + dCb8(2)
      dCbnds(1, 3) = dCb1(3) + dCb2(3) + dCb3(3) + dCb4(3) + dCb5(3) &
                     + dCb6(3) + dCb7(3) + dCb8(3)
!
   end if
!
! ----- calculate the Cbndb correction term ----------------
   B1 = B1b
   B12 = B1*B1
   B13 = B12*B1
   B14 = B13*B1
   B15 = B14*B1
   asum = d11 + d12*R + d13*Rsq + d14/R + d15/Rsq
   bsum = d21*B12 + d22*B13 + d23*B14 + d24*B15
   csum = d31*B1*exp1 + d32*B12*exp2
   dsum1 = d41*exp1 + d42*exp2
   dsum2 = d43*exp1 + d44*exp2
   fsum = dx1 + d52*R + d53*Rsq
   gsum = d61 + d62/R + d63/Rsq
   aasum = d71 + d72*P + d73*Psq + d74/P + d75/Psq
   ffsum = d81*P + d82*Psq + d84/Psq
   dasum = d12 + 2.d0*d13*R - d14/Rsq - 2.d0*d15/Rcu
   dbsum = 2.d0*d21*B1 + 3.d0*d22*B12 + 4.d0*d23*B13 + 5.d0*d24*B14
   ddsum1 = d41*dexp1 + d42*dexp2
   ddsum2 = d43*dexp1 + d44*dexp2
   dfsum = d52 + 2.d0*d53*R
   dgsum = -d62/Rsq - 2.d0*d63/Rcu
   daasum = d72 + 2.d0*d73*P - d74/Psq - 2.d0*d75/Pcu
   dffsum = d81 + 2.d0*d82*P - 2.d0*d84/Pcu
   cb(2, 1) = B1*asum*exp1/P
   cb(2, 2) = bsum*exp2
   cb(2, 3) = B2*csum
   cb(2, 4) = B1*B3*dsum1 + B1*B3b*dsum2
   cb(2, 5) = B1*fsum*exp7/P
   cb(2, 6) = B1*gsum*exp7
   cb(2, 7) = B1*aasum*exp2
   cb(2, 8) = B1*ffsum*exp7
!
   Cbnds(2) = cb(2, 1) + cb(2, 2) + cb(2, 3) + cb(2, 4) + cb(2, 5) &
              + cb(2, 6) + cb(2, 7) + cb(2, 8)
!
!       calculate the derivatives:
!
   if (ideriv .gt. 0) then
!
      dCb1(1) = dB1b(1)*asum*exp1/P + B1*dasum*exp1/P &
                + B1*asum*dexp1/P - B1*asum*exp1*dP(1)/Psq
      dCb1(2) = dB1b(2)*asum*exp1/P + B1*dasum*exp1/P &
                + B1*asum*dexp1/P - B1*asum*exp1*dP(2)/Psq
      dCb1(3) = dB1b(3)*asum*exp1/P + B1*dasum*exp1/P &
                + B1*asum*dexp1/P - B1*asum*exp1*dP(3)/Psq
!
      dCb2(1) = dbsum*dB1b(1)*exp2 + bsum*dexp2
      dCb2(2) = dbsum*dB1b(2)*exp2 + bsum*dexp2
      dCb2(3) = dbsum*dB1b(3)*exp2 + bsum*dexp2
!
      dCb3(1) = dB2(1)*csum + B2*(d31*dB1b(1)*exp1 + d31*B1*dexp1 &
                                  + d32*2.d0*B1*dB1b(1)*exp2 + d32*B12*dexp2)
      dCb3(2) = dB2(2)*csum + B2*(d31*dB1b(2)*exp1 + d31*B1*dexp1 &
                                  + d32*2.d0*B1*dB1b(2)*exp2 + d32*B12*dexp2)
      dCb3(3) = dB2(3)*csum + B2*(d31*dB1b(3)*exp1 + d31*B1*dexp1 &
                                  + d32*2.d0*B1*dB1b(3)*exp2 + d32*B12*dexp2)
!          dCb4 equations may 27/90:
      dCb4(1) = dB1b(1)*B3*dsum1 + B1*dB3(1)*dsum1 + B1*B3*ddsum1 &
                + dB1b(1)*B3b*dsum2 + B1*dB3(4)*dsum2 + B1*B3b*ddsum2
      dCb4(2) = dB1b(2)*B3*dsum1 + B1*dB3(2)*dsum1 + B1*B3*ddsum1 &
                + dB1b(2)*B3b*dsum2 + B1*dB3(5)*dsum2 + B1*B3b*ddsum2
      dCb4(3) = dB1b(3)*B3*dsum1 + B1*dB3(3)*dsum1 + B1*B3*ddsum1 &
                + dB1b(3)*B3b*dsum2 + B1*dB3(6)*dsum2 + B1*B3b*ddsum2
!          dCb5 equations may 27/90:
      dCb5(1) = dB1b(1)*fsum*exp7/P + B1*dfsum*exp7/P &
                + B1*fsum*dexp7/P - B1*fsum*exp7*dP(1)/Psq
      dCb5(2) = dB1b(2)*fsum*exp7/P + B1*dfsum*exp7/P &
                + B1*fsum*dexp7/P - B1*fsum*exp7*dP(2)/Psq
      dCb5(3) = dB1b(3)*fsum*exp7/P + B1*dfsum*exp7/P &
                + B1*fsum*dexp7/P - B1*fsum*exp7*dP(3)/Psq
!
      dCb6(1) = dB1b(1)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
      dCb6(2) = dB1b(2)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
      dCb6(3) = dB1b(3)*gsum*exp7 + B1*dgsum*exp7 + B1*gsum*dexp7
!
      dCb7(1) = dB1b(1)*aasum*exp2 + B1*daasum*dP(1)*exp2 &
                + B1*aasum*dexp2
      dCb7(2) = dB1b(2)*aasum*exp2 + B1*daasum*dP(2)*exp2 &
                + B1*aasum*dexp2
      dCb7(3) = dB1b(3)*aasum*exp2 + B1*daasum*dP(3)*exp2 &
                + B1*aasum*dexp2
      dCb8(1) = dB1b(1)*ffsum*exp7 + B1*dffsum*dP(1)*exp7 &
                + B1*ffsum*dexp7
      dCb8(2) = dB1b(2)*ffsum*exp7 + B1*dffsum*dP(2)*exp7 &
                + B1*ffsum*dexp7
      dCb8(3) = dB1b(3)*ffsum*exp7 + B1*dffsum*dP(3)*exp7 &
                + B1*ffsum*dexp7
      dCbnds(2, 1) = dCb1(1) + dCb2(1) + dCb3(1) + dCb4(1) + dCb5(1) &
                     + dCb6(1) + dCb7(1) + dCb8(1)
      dCbnds(2, 2) = dCb1(2) + dCb2(2) + dCb3(2) + dCb4(2) + dCb5(2) &
                     + dCb6(2) + dCb7(2) + dCb8(2)
      dCbnds(2, 3) = dCb1(3) + dCb2(3) + dCb3(3) + dCb4(3) + dCb5(3) &
                     + dCb6(3) + dCb7(3) + dCb8(3)
!
!  calculate the total derivative from the pieces:
      dCbnda(1) = dT(1)*Cbnds(1) + sumT*dCbnds(1, 1)
      dCbndb(1) = dT(1)*Cbnds(2) + sumT*dCbnds(2, 1)
      dCbnda(2) = dT(2)*Cbnds(1) + sumT*dCbnds(1, 2)
      dCbndb(2) = dT(2)*Cbnds(2) + sumT*dCbnds(2, 2)
      dCbnda(3) = dT(3)*Cbnds(1) + sumT*dCbnds(1, 3)
      dCbndb(3) = dT(3)*Cbnds(2) + sumT*dCbnds(2, 3)
!
   end if
!
   Cbnda = sumT*Cbnds(1)
   Cbndb = sumT*Cbnds(2)
!
   return
!
end
!
!**********************************************************************
!
subroutine getccofr_h4bmkp(r, cc, ivalid, errdr, iri)
!----------------------------------------------------------------------
!  Get cartesian coords for a of distances r(); check for error.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   dimension r(6), cc(4, 3), iri(6)
!
   parameter(eps_err=1.d-5, eps_chk=1.d-14)
!
   common/c_levelerr_h4bmkp/err_return, rmin_err, level_err
   save/c_levelerr_h4bmkp/
!
   common/c_chkh4geom_h4bmkp/ii2(6), ii3(6), ii4(6), ii5(6)
   save/c_chkh4geom_h4bmkp/
!

   do i = 1, 4
      cc(i, 1) = 0.d0
      cc(i, 2) = 0.d0
      cc(i, 3) = 0.d0
   end do
!                              ! find longest distance i1
   i1 = 1
   do j = 1, 6
   if (r(j) .lt. rmin_err) then
      errdr = err_return
      ivalid = -j
      return
   end if
   if (r(j) .gt. r(i1)) i1 = j
   end do
!                  ! complement i6, triangles i2,i4,i1 and i3,i5,i1
   i2 = ii2(i1)
   i3 = ii3(i1)
   i4 = ii4(i1)
   i5 = ii5(i1)
   i6 = 7 - i1
!
   iri(1) = i1
   iri(2) = i2
   iri(3) = i3
   iri(4) = i4
   iri(5) = i5
   iri(6) = i6
!
!-noneed;      if ( r(i1) .lt. 0.01d0 ) then
!-noneed;         ivalid = -1
!-noneed;         errdr = 99.d0
!-noneed;         return
!-noneed;      endif
!
!__Check the H3 sub-conformations for linearity or for error.
!
   call chkh3_h4bmkp(r(i1), r(i2), r(i4), ivalid1, errdr1)
   call chkh3_h4bmkp(r(i1), r(i3), r(i5), ivalid2, errdr2)
   call chkh3_h4bmkp(r(i2), r(i3), r(i6), ivalid3, errdr3)
   call chkh3_h4bmkp(r(i4), r(i5), r(i6), ivalid4, errdr4)
!
   ivalid = min(ivalid1, ivalid2, ivalid3, ivalid4)
   errdr = max(errdr1, errdr2, errdr3, errdr4)
!                                          ! bad H3 sub-geom?
   if (ivalid .lt. 0) then
      ivalid = -7
      return
   end if
!
!__Get the H4 conformation cartesian coordinates, checking for error.
!
   if (errdr1 .ge. 0.d0) then
      xa = r(i2)
      ya = 0.d0
   else
      xa = (r(i1)**2 + r(i2)**2 - r(i4)**2)*0.5d0/r(i1)
      ya = sqrt(max(r(i2)**2 - xa**2, 0.d0))
   end if
!
   if (errdr2 .ge. 0.d0) then
      xb = r(i3)
      yb = 0.d0
   else
      xb = (r(i1)**2 + r(i3)**2 - r(i5)**2)*0.5d0/r(i1)
      yb = sqrt(max(r(i3)**2 - xb**2, 0.d0))
   end if
!                                     ! min and max r(i6) values
   r6_lo = sqrt((xb - xa)**2 + (yb - ya)**2)
   r6_hi = sqrt((xb - xa)**2 + (yb + ya)**2)
!
   errdr = max(r(i6) - r6_hi, r6_lo - r(i6), errdr)
!                                                ! bad r(i6) ?
   if (errdr .gt. eps_err) then
      ivalid = -8
      return
   else if (errdr .gt. eps_chk) then
      ivalid = 0
   end if
!                  ! finish getting H4 cartesian coordinates
   cc(2, 1) = r(i1)
   cc(3, 1) = xa
   cc(4, 1) = xb
   cc(3, 2) = ya
!
   if (min(ya, yb) .eq. 0.d0) then
      cc(4, 2) = yb
   else
      cc(4, 2) = (ya**2 + yb**2 + (xb - xa)**2 - r(i6)**2)*0.5d0 &
                 /ya
      cc(4, 3) = sqrt(max(0.d0, yb**2 - cc(4, 2)**2))
   end if
!
   return
end
!
!**********************************************************************
!
subroutine chkh3_h4bmkp(r1, r2, r3, ivalid, errdr)
!----------------------------------------------------------------------
!  Check for validity of input H3 geometry.
!----------------------------------------------------------------------
!
   implicit double precision(a - h, o - z)
   implicit integer(i - n)

!
   parameter(eps_err=1.d-5, eps_chk=1.d-14)
!
   rhi = max(r1, r2)
   if (r3 .ge. rhi) then
      rmid = rhi
      rhi = r3
   else
      rmid = max(r3, min(r1, r2))
   end if
   rlo = min(r1, r2, r3)
   errdr = rhi - rmid - rlo
   if (errdr .gt. eps_err) then
      ivalid = -1
   else if (errdr .gt. eps_chk) then
      ivalid = 0
   else
      ivalid = 1
   end if
!
   return
end
!
