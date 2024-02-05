subroutine desagrupa
   use variables
   use constants
   use phys_parameters
   implicit none

   integer      :: i, j, m, switchout, it
   double precision      :: rxij, ryij, rzij, rij, r1, r2, r3, r4, r5, r6
   double precision      :: potential_1, potential_2
   integer      :: gi3, gi4, gi5, gi6, mgj8, mgi8
   character*50 :: texto !######### for printing with pgroup ##########


   interface
      subroutine pgroup(tipo, i, event_str)

         integer, intent(in):: i, tipo
         character*50, optional :: event_str
      end subroutine pgroup
   end interface

! to check if connection of molecules due to a nonreacting collision is still valid

   parti_a: do i = 1, npart - 1
   if (group(i, 8) .ne. 0) then
      mgi8 = group(i, 8)
      call out_bound8(i, mgi8, switchout)
      if (switchout == 1) then
         call print_line
         texto = 'Quebra de complexo não ligante'
         call pgroup(1, i)
         call pgroup(1, mgi8)
         call print_line

         group(mgi8, 8) = 0
         group(i, 8) = 0

!for checking OH states
         if (group(i, 2) .eq. 3) call write_diat_posi_velo(3, group(i, 3), group(i, 4))
         if (group(mgi8, 2) .eq. 3) call write_diat_posi_velo(3, group(mgi8, 3), group(mgi8, 4))

      end if

   end if

   end do parti_a

   exgroup = 0
   m = 0
   j = 1

   partic: do i = 1, npart

      select case (group(i, 2))

      case (0)

         continue

      case (-2, -1)

         exgroup(j, :) = group(i, :)
         j = j + 1

         cycle partic

      case (1) ! H2

         exgroup(j, :) = group(i, :)

         j = j + 1

         gi3 = group(i, 3)
         gi4 = group(i, 4)

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         if (rij .ge. out_bound) then

            write (unit=10, fmt='(a)') 'H2 -> H(1) + H(2)'
            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2 -> H(1) + H(2)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if
         cycle partic

      case (14)

         exgroup(j, :) = group(i, :)

         j = j + 1

         gi3 = group(i, 3)
         gi4 = group(i, 4)

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         if (rij .ge. out_bound) then

            write (unit=10, fmt='(a)') 'H2t -> H(1) + H(2)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2t -> H(1) + H(2)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if
         cycle partic

      case (2) ! O2 (3Σ)

         exgroup(j, :) = group(i, :)

         j = j + 1

         gi3 = group(i, 3)
         gi4 = group(i, 4)

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         if (rij .ge. out_bound) then
            m = m + 1

            print *, 'corte desagrupa'

            write (unit=10, fmt='(a)') 'O2 -> O(1) + O(2)'

            j = j - 1

            call print_line
            texto = 'O2 -> O(1) + O(2)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -2
            exgroup(j, 3) = gi3
            exgroup(j, 4:6) = 0
            exgroup(j, 7) = 3
            exgroup(j, 8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if
         cycle partic
      case (15) ! O2(1Δ)

         exgroup(j, :) = group(i, :)

         j = j + 1

         gi3 = group(i, 3)
         gi4 = group(i, 4)

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         if (rij .ge. out_bound) then
            m = m + 1

            print *, 'corte desagrupa'

            write (unit=10, fmt='(a)') 'O2(1Delta) -> O(1) + O(2)'

            j = j - 1

            call print_line
            texto = 'O2(1Delta) -> O(1) + O(2)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -2
            exgroup(j, 3) = gi3
            exgroup(j, 4:6) = 0
            exgroup(j, 7) = 3
            exgroup(j, 8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         cycle partic

      case (3) !  OH(2PI)

         exgroup(j, :) = group(i, :)

         j = j + 1

         gi3 = group(i, 3)
         gi4 = group(i, 4)

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         if (rij .ge. out_bound) then

            m = m + 1
            write (unit=10, fmt='(a)') 'OH(2PI) -> H(1) + O(2)'
            j = j - 1

            !for checking OH states
            call write_diat_posi_velo(4,gi3,gi4)

            call print_line
            texto = 'OH(2PI) -> H(1) + O(2)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if
         cycle partic

      case (17) !  OH(4)
         exgroup(j, :) = group(i, :)

         j = j + 1

         gi3 = group(i, 3)
         gi4 = group(i, 4)

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         if (rij .ge. out_bound) then

            write (unit=10, fmt='(a)') 'OH(4) -> H(1) + O(2)'

            !for checking OH states
            call write_diat_posi_velo(4,gi3,gi4)

            m = m + 1
            j = j - 1

            call print_line
            texto = 'OH(4) -> H(1) + O(2)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if
         cycle partic

      case (4)      !for the ho2 surface

         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! O atom
         gi5 = group(i, 5)      ! O atom

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)
         r1 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)
         r2 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi3), ry(gi3), rz(gi3), rxij, ryij, rzij, rij)
         r3 = rij

         exgroup(j, :) = group(i, :)
         j = j + 1

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'HO2 -> O(2) + OH'

            m = m + 1
            j = j - 1

            call print_line
            texto =  'HO2 -> O(2) + OH'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi3, gi5)

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r1 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'HO2 -> H(1) + O2'
            m = m + 1
            j = j - 1

            call print_line
            texto = 'HO2 -> H(1) + O2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'HO2 -> O(3) + OH'

            m = m + 1
            j = j - 1

            call print_line
            texto =   'HO2 -> O(3) + OH'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi3, gi4)

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .ge. out_bound)) stop 'wrong geo'

         cycle partic

      case (5)      !for the h2o surfaces

         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! H atom
         gi5 = group(i, 5)      ! O atom

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r1 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi3), ry(gi3), rz(gi3), rxij, ryij, rzij, rij)

         r3 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .lt. out_bound)) then

            m = m + 1
            write (unit=10, fmt='(a)') 'H2O -> H(2) + OH'

            j = j - 1

            call print_line
            texto = 'H2O -> H(2) + OH'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi3, gi5)

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r1 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'H2O -> H(a) + OH'
            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O -> H(a) + OH'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 3
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi4, gi5)

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'H2O -> O(b) + H2'
            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O -> O(b) + H2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 1
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .ge. out_bound)) stop 'wrong geo'
         cycle partic

      case (6)      !for the h2otrip surfaces

         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! H atom
         gi5 = group(i, 5)      ! O atom

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r1 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi3), ry(gi3), rz(gi3), rxij, ryij, rzij, rij)

         r3 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'H2O(t) -> H(a) + OH'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O(t) -> H(a) + OH'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi3, gi5)

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r1 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'H2O(t) -> H(b) + OH'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O(t) -> H(b) + OH'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 3
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

!for checking OH states
            call write_diat_posi_velo(1, gi4, gi5)

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'H2O(t) -> O(3) + H2'

            m = m + 1
            j = j - 1

!            poti1 = vh2ot_f6(r2, r1, r3)/con

!            call pothhs_f1(r2, V, dV)
!            poti2 = V/con



            call print_line
            texto = 'H2O(t) -> O(3) + H2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .ge. out_bound)) stop 'wrong geo'
         cycle partic

      case (7)      !for the O3

         gi3 = group(i, 3)      ! O atom
         gi4 = group(i, 4)      ! O atom
         gi5 = group(i, 5)      ! O atom

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r1 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi3), ry(gi3), rz(gi3), rxij, ryij, rzij, rij)

         r3 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'O3 -> O(2) + O2'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O3 -> O(2) + O2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 2
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r1 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'O3 -> O(1) + O2'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O3 -> O(1) + O2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -2
            exgroup(j, 3) = gi3
            exgroup(j, 4:6) = 0
            exgroup(j, 7) = 3
            exgroup(j, 8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .lt. out_bound)) then

            m = m + 1
            write (unit=10, fmt='(a)') 'O3 -> O(3) + O2'

            j = j - 1

            print *, 'tres'

            group(exgroup(j, 8), 8) = 0

            call print_line
            texto = 'O3 -> O(3) + O2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 2
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = 0
            exgroup(j, 6) = 0
            exgroup(j, 7) = 0
            exgroup(j, 8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4) = 0
            exgroup(j + 1, 5) = 0
            exgroup(j + 1, 6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .ge. out_bound)) stop 'wrong geo'

         cycle partic

      case (8)      !for the H3

         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! H atom
         gi5 = group(i, 5)      ! H atom

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r1 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi3), ry(gi3), rz(gi3), rxij, ryij, rzij, rij)

         r3 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .lt. out_bound)) then

            m = m + 1
            write (unit=10, fmt='(a)') 'H3 -> H(2) + H2'

            j = j - 1

            call print_line
            texto = 'H3 -> H(2) + H2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!  update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r1 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'H3 -> H(1) + H2'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3 -> H(1) + H2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .lt. out_bound)) then

            write (unit=10, fmt='(a)') 'H3 -> H(3) + H2'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3 -> H(3) + H2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r1 .ge. out_bound) .and. (r2 .ge. out_bound)) stop 'wrong geo'
         cycle partic

      case (9)      !for the H4 surface

!The positions are predefined in the matrix "group" with H atoms first and O later
         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! H atom
         gi5 = group(i, 5)      ! H atom
         gi6 = group(i, 6)      ! H atom

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r3 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r5 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r6 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r4 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r1 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'H4 -> H(1) + H3'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H4 -> H(1) + H3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 3
            exgroup(j + 1, 2) = 8
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5) = gi6
            exgroup(j + 1, 6:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'H4 -> H(2) + H3'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H4 -> H(2) + H3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 8
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'H4 -> H(3) + H3'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H4 -> H(3) + H3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 8
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
!                  stop

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'H4 -> H(4) + H3'

            m = m + 1
            j = j - 1

            call print_line
            call pgroup(3, j, texto)
            texto = 'H4 -> H(4) + H3'
            exgroup(j, 1) = 3
            exgroup(j, 2) = 8
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi5
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi6
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
!                  stop

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'H4 -> H2(1-3) + H2(2-4)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H4 -> H2(1-3) + H2(2-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r4 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'H4 -> H2(1-4) + H2(2-3)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H4 -> H2(1-4) + H2(2-3)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi6
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'H4 -> H2(1-4) + H2(2-3)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H4 -> H2(1-4) + H2(2-3)'
            call pgroup(3, j, texto)

            exgroup(j, 7) = 0
            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 1
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         cycle partic

      case (10)      !for the H2O2 (1A) surface

!The positions are predefined in the matrix "group" with H atoms first and O later
         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! H atom
         gi5 = group(i, 5)      ! O atom
         gi6 = group(i, 6)      ! O atom

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r3 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r5 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r6 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r4 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r1 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(1A) -> H + HO2'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(1A) -> H + HO2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 3
            exgroup(j + 1, 2) = 4
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5) = gi6
            exgroup(j + 1, 6:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(1A) -> H + HO2'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(1A) -> H + HO2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 4
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(1A) -> O + H2O'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(1A) -> O + H2O'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 5
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 1
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(1A) -> O + H2O'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(1A) -> O + H2O'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 5
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi5
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi6
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 1
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(1A) -> OH + OH'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(1A) -> OH + OH'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 3
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

!for checking OH states
            call write_diat_posi_velo(1, gi4, gi6)
            call write_diat_posi_velo(1, gi3, gi5)

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r4 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(1A) -> OH + OH'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(1A) -> OH + OH'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi6
            exgroup(j, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi3, gi6)

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 3
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi4, gi5)

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if
         if ((r3 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(1A) -> H2 + O2'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(1A) -> H2 + O2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 15! forOD
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         cycle partic

      case (16)      !for the H2O2 (3A) surface

!The positions are predefined in the matrix "group" with H atoms first and O later
         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! H atom
         gi5 = group(i, 5)      ! O atom
         gi6 = group(i, 6)      ! O atom

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r3 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r5 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r6 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r4 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r1 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(3A) -> H(2S) + HO2(2A")'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(3A) -> H(2S) + HO2(2A")'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 3
            exgroup(j + 1, 2) = 4
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5) = gi6
            exgroup(j + 1, 6:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(3A) -> H(2S) + HO2(2A")'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(3A) -> H(2S) + HO2(2A")'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 4
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(3A) -> O(3P) + H2O(1A`)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(3A) -> O(3P) + H2O(1A`)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 5
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(3A) -> O(3P) + H2O(1A`)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(3A) -> O(3P) + H2O(1A`)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 5
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi5
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi6
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(3A) -> OH(2Pi) + OH(2Pi)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(3A) -> OH(2Pi) + OH(2Pi)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi3, gi5)

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 3
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi4, gi6)

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r4 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(3A) -> OH(2Pi) + OH(2Pi)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(3A) -> OH(2Pi) + OH(2Pi)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi6
            exgroup(j, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi3, gi6)

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 3
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

!for checking OH states
            call write_diat_posi_velo(1, gi4, gi5)

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(A,I10)') 'tempo=', (tempo + old_step_print)
            write (unit=10, fmt='(a)') 'H2O2(3A) -> H2(1Sigma) + O2(3Sigma)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H2O2(3A) -> H2(1Sigma) + O2(3Sigma)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if
         cycle partic

      case (11)      !for the O4 surface

!The positions are predefined in the matrix "group"
         gi3 = group(i, 3)      ! O atom
         gi4 = group(i, 4)      ! O atom
         gi5 = group(i, 5)      ! O atom
         gi6 = group(i, 6)      ! O atom

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r1 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r4 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r5 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r6 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r3 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r1 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'O4 -> O(1) + O3'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O4 -> O(1) + O3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -2
            exgroup(j, 3) = gi3
            exgroup(j, 4:6) = 0
            exgroup(j, 7) = 3
            exgroup(j, 8) = 0

            exgroup(j + 1, 1) = 3
            exgroup(j + 1, 2) = 7
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5) = gi6
            exgroup(j + 1, 6:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'O4 -> O(2) + O3'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O4 -> O(2) + O3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 7
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'O4 -> O(3) + O3'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O4 -> O(3) + O3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 7
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'O4 -> O(4) + O3'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O4 -> O(4) + O3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 7
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi5
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi6
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'O4 -> O2(1-3) + O2(2-4)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O4 -> O2(1-3) + O2(2-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 2
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r4 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'O4 -> O2(1-4) + O2(2-3)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O4 -> O2(1-4) + O2(2-3)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 2
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi6
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            write (unit=10, fmt='(a)') 'O4 -> O2(1-2) + O2(3-4)'

            m = m + 1
            j = j - 1

            call print_line
            texto = 'O4 -> O2(1-2) + O2(3-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 2
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

!update nonreactive connections
            call check_connections(i)

         end if

         cycle partic

      case (12)      !for the H3O surface

!The positions are predefined in the matrix "group" with H atoms first and O later
         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! H atom
         gi5 = group(i, 5)      ! H atom
         gi6 = group(i, 6)      ! O atom

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)

         r2 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r3 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r5 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)

         r6 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r4 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r1 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3O -> H(1) + H2O'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 3
            exgroup(j + 1, 2) = 5
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5) = gi6
            exgroup(j + 1, 6:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'H3O -> H(1) + H2O'

            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3O -> H(2) + H2O'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 5
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'H3O -> H(2) + H2O'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3O -> H(3) + H2O'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 5
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -1
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

            write (unit=10, fmt='(a)') 'H3O -> H(3) + H2O'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3O -> O + H3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 8
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi5
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi6
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

            write (unit=10, fmt='(a)') 'H3O -> O + H3'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3O -> H2(1-3) + HO(2-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 3
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

!for checking HO states
            call write_diat_posi_velo(1, gi4, gi6)

            j = j + 2
            write (unit=10, fmt='(a)') 'H3O -> H2(1-3) + HO(2-4)'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r4 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3O -> H2(2-3) + HO(1-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi6
            exgroup(j, 5:8) = 0

!for checking HO states
            call write_diat_posi_velo(1, gi3, gi6)

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 1
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'H3O -> H2(2-3) + HO(1-4)'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'H3O -> H2(1-2) + HO(3-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 1
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 3
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

!for checking HO states
            call write_diat_posi_velo(1, gi5, gi6)

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'H3O -> H2(1-2) + HO(3-4)'

!update nonreactive connections
            call check_connections(i)

         end if

         cycle partic
      case (13)      !for the HO3 surface

!The positions are predefined in the matrix "group" with H atoms first and O later
         gi3 = group(i, 3)      ! H atom
         gi4 = group(i, 4)      ! O atom
         gi5 = group(i, 5)      ! O atom
         gi6 = group(i, 6)      ! O atom

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi4), ry(gi4), rz(gi4), rxij, ryij, rzij, rij)
         r2 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)
         r3 = rij

         call mod_dist(rx(gi3), ry(gi3), rz(gi3), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)
         r5 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi5), ry(gi5), rz(gi5), rxij, ryij, rzij, rij)
         r6 = rij

         call mod_dist(rx(gi4), ry(gi4), rz(gi4), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)
         r4 = rij

         call mod_dist(rx(gi5), ry(gi5), rz(gi5), rx(gi6), ry(gi6), rz(gi6), rxij, ryij, rzij, rij)

         r1 = rij

         exgroup(j, :) = group(i, :)

         j = j + 1

         if ((r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            m = m + 1

            j = j - 1

            call print_line
            texto = 'HO3 -> H(1) + O3'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 1
            exgroup(j, 2) = -1
            exgroup(j, 3) = gi3
            exgroup(j, 4:8) = 0

            exgroup(j + 1, 1) = 3
            exgroup(j + 1, 2) = 7
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5) = gi6
            exgroup(j + 1, 6:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'HO3 -> H(1) + O3'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r2 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'HO3 -> O(2) + HO2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 4
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'HO3 -> O(2) + HO2'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'HO3 -> O(3) + HO2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 4
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi6
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

            write (unit=10, fmt='(a)') 'HO3 -> O(3) + HO2'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'HO3 -> O(4) + HO2'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 3
            exgroup(j, 2) = 4
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5) = gi5
            exgroup(j, 6:8) = 0

            exgroup(j + 1, 1) = 1
            exgroup(j + 1, 2) = -2
            exgroup(j + 1, 3) = gi6
            exgroup(j + 1, 4:6) = 0
            exgroup(j + 1, 7) = 3
            exgroup(j + 1, 8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2

            write (unit=10, fmt='(a)') 'HO3 -> O(4) + HO2'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'HO3 -> HO(1-3) + O2(2-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi5
            exgroup(j, 5:8) = 0

!for checking HO states
            call write_diat_posi_velo(1, gi3, gi5)

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'HO3 -> HO(1-3) + O2(2-4)'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r1 .ge. out_bound) .and. (r2 .ge. out_bound) .and. (r3 .ge. out_bound) .and. (r4 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'HO3 -> O2(2-3) + HO(1-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi6
            exgroup(j, 5:8) = 0

!for checking HO states
            call write_diat_posi_velo(1, gi3, gi6)

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi4
            exgroup(j + 1, 4) = gi5
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'HO3 -> O2(2-3) + HO(1-4)'

!update nonreactive connections
            call check_connections(i)

         end if

         if ((r3 .ge. out_bound) .and. (r4 .ge. out_bound) .and. (r5 .ge. out_bound) .and. (r6 .ge. out_bound)) then

            m = m + 1
            j = j - 1

            call print_line
            texto = 'HO3 -> HO(1-2) + O2(3-4)'
            call pgroup(3, j, texto)

            exgroup(j, 1) = 2
            exgroup(j, 2) = 3
            exgroup(j, 3) = gi3
            exgroup(j, 4) = gi4
            exgroup(j, 5:8) = 0

!for checking HO states
            call write_diat_posi_velo(1, gi3, gi4)

            exgroup(j + 1, 1) = 2
            exgroup(j + 1, 2) = 2
            exgroup(j + 1, 3) = gi5
            exgroup(j + 1, 4) = gi6
            exgroup(j + 1, 5:8) = 0

            call print_points
            call pgroup(3, j)
            call pgroup(3, j + 1)
            call print_line

            j = j + 2
            write (unit=10, fmt='(a)') 'HO3 -> HO(1-2) + O2(3-4)'

!update nonreactive connections
            call check_connections(i)

         end if

         cycle partic

      end select

   end do partic

   if (m .ne. 0) then
!    print*,'m=',m
      call potential_total
      potential_1 = potential
      write (unit=10, fmt='(a,i11,/,a,D20.11)') 'Step=', (tempo + old_step_print), "V=", potential/con

      group = exgroup

!    print*,'desagrupa depois'

      call potential_total
      potential_2 = potential
      potential_change = potential_change + potential_2 - potential_1
      write (unit=10, fmt='(a,D20.11)') "V=", potential/con
      write (unit=10, fmt='(a,D20.11)') "V_change=", potential_change/con

   end if

   check: do it = 1, npart - 1

      if (group(it, 8) /= 0) then
         mgj8 = group(it, 8)

         if (group(mgj8, 8) == 0) then

            group(it, 8) = 0

         end if
      end if

   end do check

end subroutine desagrupa

subroutine out_bound8(iin, kin, swout)
   use variables
   use constants
   use phys_parameters
   implicit none
   integer, intent(in) :: iin, kin
   double precision      :: rxij, ryij, rzij, rij
   integer, intent(out) ::swout
   integer      :: switch, i, j, gi, gk

   swout = 1

   test: do i = 3, 2 + group(iin, 1)
   do j = 3, 2 + group(kin, 1)

      gi = group(iin, i)
      gk = group(kin, j)

      call mod_dist(rx(gi), ry(gi), rz(gi), rx(gk), ry(gk), rz(gk), rxij, ryij, rzij, rij)
!Verifica distancias em relacao a distancia de aglumeracao

      if (rij .gt. (out_bound)) then

         switch = 1
      else
         switch = 0

      end if

      swout = swout*switch

   end do

   end do test

end subroutine out_bound8

subroutine check_connections(il)
!subroutine goes through the connections column and
!updates values the values of the moleculesi/atoms lines.
!-exgroup(i,8) if it is already copyed to exgroup (im<il)
!-group(i,8) there is no copy yet (im>il).
   use variables
   integer :: il, im

   do im = 1, npart
   if (im .lt. il) then

      if (exgroup(im, 8) .gt. il) exgroup(im, 8) = exgroup(im, 8) + 1

   else

      if (group(im, 8) .gt. il) group(im, 8) = group(im, 8) + 1

   end if
   end do

end subroutine check_connections

