subroutine read_input_file()

   use constants
   use variables
   use phys_parameters
   use bk_variables
   use sim_variables
   use chem_parameters
   implicit none

   double precision :: xrand(3)

   integer      :: io
   integer      :: i, j

   integer(kind=8)::old_step
   integer :: g_elements, at_elements

!character(80) :: rubbish,answer
!for checking optimization
!allocate(times(0:t_tot,17))

   print *, 'reading file input_file'

   if (out_bound .lt. in_bound) stop 'outer dist. can not be smaller then inner dist. in input file'

   call random_seed(put=iseed)

   write (bk_new_string, '(i12)') bk_new
   bk_new_string = adjustl(bk_new_string)

   write (bk_old_string, '(i12)') bk_new - 1

   bk_old_string = adjustl(bk_old_string)

   open (unit=41, file='out/geo'//trim(bk_new_string)//'.bk', status='REPLACE')
   close (41)

   if (bk_new .ge. 1) then

      dt = new_dt
      t_tot = new_t_tot
      tprint = new_tprint
      tbackup = new_tbackup

      open (unit=20, file='out/ini.bk', access="stream")

      read (unit=20) old_npart, old_temperature, old_boxl, &
         old_out_bound, old_in_bound, old_rl_cut_off, old_dt, cname, mass
      rmass = 1.0d0/mass

      write (*, *) 'Before restarting check initial configuration!'

      write (unit=*, fmt='(A20,I15)') ' which_step=', which_step
      write (unit=*, fmt='(A20,I15)') ' t_tot=', t_tot
      write (unit=*, fmt='(A20,I15)') 'old_npart=', old_npart
      write (unit=*, fmt='(A20,I15)') 'old_temp(K)=', int(old_temperature)
      write (unit=*, fmt='(A20,F15.3)') 'old_boxl(nm)=', old_boxl/rnm
      write (unit=*, fmt='(A20,F15.3)') 'old_out_bound(nm)=', old_out_bound/rnm
      write (unit=*, fmt='(A20,F15.3)') 'old_in_bound(nm)=', old_in_bound/rnm
      write (unit=*, fmt='(A20,F15.3)') 'old_rl_cut_off(nm)=', old_rl_cut_off/rnm
      write (unit=*, fmt='(A20,F15.3)') 'old_dt(fs)=', old_dt/fentosec
!write(unit=*,fmt='(A20,I15)')'old_nunit_max=',old_nunit_max

      close (20)

      group = 0
      group2 = 0
      exgroup = 0

      open (unit=21, file='out/geo'//trim(bk_old_string)//'.bk', access='stream')

      io = 0
      bk_loop: do
         read (unit=21, IOSTAT=io) old_step

         print *, 'old_step', old_step

         IF (io > 0) STOP 'Check input.  Something was wrong in reading old_step'

         IF (io < 0) THEN
            print *, 'Backup file finished  and new step was not found!!'
            print *, 'Restart and continue from last step?', old_step
            stop
         END IF

         print *, 'reading step', old_step, 'from backup'
         print *, 'OK'
         if (old_step .ge. (which_step + 2)) then
            print *, 'old_step', old_step, 'which_step', which_step
            exit bk_loop
            print *, 'exit 1'
         else
            print *, 'reading'
         end if
         read (unit=21, IOSTAT=io) total_real_time
         if (io .ne. 0) stop '0 reading backup file problem'
         print *, 'trt', total_real_time
         read (unit=21, IOSTAT=io) group
         if (io .ne. 0) stop '1 reading backup file problem'
!read(unit=21,IOSTAT=io)cname
!if (io.ne.0) stop '2 reading backup file problem'
         read (unit=21, IOSTAT=io) rx, ry, rz
         if (io .ne. 0) stop '4 reading backup file problem'

         read (unit=21, IOSTAT=io) vx, vy, vz
         if (io .ne. 0) stop '5 reading backup file problem'

         read (unit=21, IOSTAT=io) ax, ay, az
         if (io .ne. 0) stop '4 reading backup file problem'

         read (unit=21, IOSTAT=io) ekin
         if (io .ne. 0) stop '5 reading backup file problem'

         read (unit=21, IOSTAT=io) potential
         if (io .ne. 0) stop '7 reading backup file problem'

         read (unit=21, IOSTAT=io) iget
         if (io .ne. 0) stop '10 reading backup file problem'

!Saving old_step for printing porpuses, used in write_geom.f90
         old_step_print = old_step

         if (old_step .eq. which_step) exit bk_loop

      end do bk_loop

      close (21)

      close (8)

      call random_seed(put=iget)

      call potential_total

!write (*,*) "##########from recalculate#################################"
!call printgroup
!write (*,*) "##########from recalculate#################################"

      go to 123

   end if

   group = 0
   group2 = 0
   exgroup = 0

   g_elements = 0  !counting molecules
   at_elements = 0 !counting atoms

!Number of particles per molecule
!Type of Particle
! and loop for filling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!for H
   group(g_elements + 1:g_elements + natoms_H, 1) = 1
   group(g_elements + 1:g_elements + natoms_H, 2) = -1

   do j = 1, natoms_H
      i = g_elements + j
      group(i, 3) = at_elements + 1
      at_elements = at_elements + 1

   end do

   g_elements = g_elements + natoms_H

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!for O

   group(g_elements + 1:g_elements + natoms_O, 1) = 1
   group(g_elements + 1:g_elements + natoms_O, 2) = -2

   do j = 1, natoms_O
      i = g_elements + j
      group(i, 3) = at_elements + 1

      call random_number(xrand(1))
!print*,'xrand(1)',xrand(1)
      if (xrand(1) .le. 1.0d0/2.0d0) then
         group(i, 7) = 3
      else
         group(i, 7) = 3
      end if

      at_elements = at_elements + 1
   end do

   g_elements = g_elements + natoms_O

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for HH
   group(g_elements + 1:g_elements + nmolec_HH, 1) = 2
   group(g_elements + 1:g_elements + nmolec_HH, 2) = 1
   do j = 1, nmolec_HH
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      at_elements = at_elements + 2
   end do

   g_elements = g_elements + nmolec_HH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for HO
   group(g_elements + 1:g_elements + nmolec_HO, 1) = 2
   group(g_elements + 1:g_elements + nmolec_HO, 2) = 3

   do j = 1, nmolec_HO
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      at_elements = at_elements + 2
   end do

   g_elements = g_elements + nmolec_HO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for O2sing

   group(g_elements + 1:g_elements + nmolec_OOsing, 1) = 2
   group(g_elements + 1:g_elements + nmolec_OOsing, 2) = 15

   do j = 1, nmolec_OOsing
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      at_elements = at_elements + 2
   end do

   g_elements = g_elements + nmolec_OOsing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for OOt
   group(g_elements + 1:g_elements + nmolec_OOtrip, 1) = 2
   group(g_elements + 1:g_elements + nmolec_OOtrip, 2) = 2

   do j = 1, nmolec_OOtrip
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      at_elements = at_elements + 2
   end do

   g_elements = g_elements + nmolec_OOtrip

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for H3
   group(g_elements + 1:g_elements + nmolec_H3, 1) = 3
   group(g_elements + 1:g_elements + nmolec_H3, 2) = 8
   do j = 1, nmolec_H3
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      at_elements = at_elements + 3
   end do

   group(:, 6) = 0

   g_elements = g_elements + nmolec_H3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for H2O
   group(g_elements + 1:g_elements + nmolec_H2O, 1) = 3
   group(g_elements + 1:g_elements + nmolec_H2O, 2) = 5

   do j = 1, nmolec_H2O
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      at_elements = at_elements + 3
   end do

   g_elements = g_elements + nmolec_H2O

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for H2Ot
   group(g_elements + 1:g_elements + nmolec_H2Ot, 1) = 3
   group(g_elements + 1:g_elements + nmolec_H2Ot, 2) = 6

   do j = 1, nmolec_H2Ot
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      at_elements = at_elements + 3
   end do

   g_elements = g_elements + nmolec_H2Ot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for HO2

   group(g_elements + 1:g_elements + nmolec_HO2, 1) = 3
   group(g_elements + 1:g_elements + nmolec_HO2, 2) = 4

   do j = 1, nmolec_HO2
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      at_elements = at_elements + 3
   end do

   g_elements = g_elements + nmolec_HO2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for O3
   group(g_elements + 1:g_elements + nmolec_O3, 1) = 3
   group(g_elements + 1:g_elements + nmolec_O3, 2) = 7
   do j = 1, nmolec_O3
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      at_elements = at_elements + 3
   end do

   g_elements = g_elements + nmolec_O3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for H4

   group(g_elements + 1:g_elements + nmolec_H4, 1) = 4
   group(g_elements + 1:g_elements + nmolec_H4, 2) = 9

   do j = 1, nmolec_H4
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      group(i, 6) = at_elements + 4
      at_elements = at_elements + 4
   end do

   g_elements = g_elements + nmolec_H4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for H3O

   group(g_elements + 1:g_elements + nmolec_H3O, 1) = 4
   group(g_elements + 1:g_elements + nmolec_H3O, 2) = 12
   do j = 1, nmolec_H3O
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      group(i, 6) = at_elements + 4
      at_elements = at_elements + 4
   end do

   g_elements = g_elements + nmolec_H3O

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for H2O2

   group(g_elements + 1:g_elements + nmolec_H2O2s, 1) = 4
   group(g_elements + 1:g_elements + nmolec_H2O2s, 2) = 10
   do j = 1, nmolec_H2O2s
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      group(i, 6) = at_elements + 4
      at_elements = at_elements + 4
   end do

   g_elements = g_elements + nmolec_H2O2s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for H2O2t

   group(g_elements + 1:g_elements + nmolec_H2O2t, 1) = 4
   group(g_elements + 1:g_elements + nmolec_H2O2t, 2) = 16

   do j = 1, nmolec_H2O2t
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      group(i, 6) = at_elements + 4
      at_elements = at_elements + 4
   end do

   g_elements = g_elements + nmolec_H2O2t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for HO3
   group(g_elements + 1:g_elements + nmolec_HO3, 1) = 4
   group(g_elements + 1:g_elements + nmolec_HO3, 2) = 13

   do j = 1, nmolec_HO3
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      group(i, 6) = at_elements + 4
      at_elements = at_elements + 4
   end do

   g_elements = g_elements + nmolec_HO3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!for O4

   group(g_elements + 1:g_elements + nmolec_O4, 1) = 4
   group(g_elements + 1:g_elements + nmolec_O4, 2) = 11

   do j = 1, nmolec_O4
      i = g_elements + j
      group(i, 3) = at_elements + 1
      group(i, 4) = at_elements + 2
      group(i, 5) = at_elements + 3
      group(i, 6) = at_elements + 4
      at_elements = at_elements + 4
   end do

   g_elements = g_elements + nmolec_O4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   group(:, 8) = 0

   write (*, *) "Your choice of parameters:"
   write (*, *) "Number of H atoms", natoms_H
   write (*, *) "Number of O atoms", natoms_O
   write (*, *) "Number of H2 molecules", nmolec_HH
   write (*, *) "Number of O2 molecules", nmolec_OOtrip
   write (*, *) "Number of O2 molecules", nmolec_OOsing
   write (*, *) "Number of HO molecules", nmolec_HO
   write (*, *) "Number of HO2 molecules", nmolec_HO2
   write (*, *) "Number of H2O molecules", nmolec_H2O
   write (*, *) "Number of H2O molecules", nmolec_H2Ot
   write (*, *) "Number of O3 molecules", nmolec_O3
   write (*, *) "Number of H3 molecules", nmolec_H3
   write (*, *) "Number of H4 molecules", nmolec_H4
   write (*, *) "Number of H2O2s molecules", nmolec_H2O2s
   write (*, *) "Number of H2O2t molecules", nmolec_H2O2t
   write (*, *) "Number of O4 molecules", nmolec_O4
   write (*, *) "Number of H3O molecules", nmolec_H3O
   write (*, *) "Number of HO3 molecules", nmolec_HO3

   write (*, *) "total number of particles", npart
   write (*, *) "Temperature distribution", temp_distribution
   write (*, *) "Temperature [K]", temp
   write (*, *) "Box size [nm]:", boxl/rnm
   write (*, *) "Cut distance [nm]:", out_bound/rnm
   write (*, *) "Mix distance [nm]:", in_bound/rnm
   write (*, *) "Long range cut off [nm]:", rl_cut_off/rnm
   write (*, *) "Integration step [fs]:", dt/fentosec
   write (*, *) "Total integration steps ", t_tot
   write (*, *) "Total time [fs]", dt*t_tot/fentosec
   write (*, *) "Printing step", tprint
   write (*, *) "Random iniciator", iseed
! para backup
123 continue

!      open(unit=41,file='out/geo'//trim(bk_new_string)//'.bk',form='unformatted')
   open (unit=41, file='out/geo'//trim(bk_new_string)//'.bk', access='stream')

   close (41)

!For constant parameters
   open (unit=20, file='out/ini'//trim(bk_new_string)//'.bk', access='stream')

   write (unit=20) npart, temp, boxl, out_bound, in_bound, rl_cut_off, dt!,nunit_max
   close (20)

!if (bk_new.ge.1) then
!write(*,*)'new'
!do i = 1, npart
!
!write (*,fmt='(A3,3F18.12,I8)')cname(i),(rx(i))/rnm,(ry(i))/rnm,(rz(i))/rnm,i
!
!end do
!end if

END subroutine read_input_file
