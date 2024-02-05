subroutine agrupa
   use variables
   use phys_parameters
   use sim_variables
   use constants

   implicit none
   integer      :: i, k, iatom, katom, mgi8, mgk8, j

   double precision      :: xrand(3)
   integer      :: iswitch
   character*50 :: texto

!kdtree variables
   double precision :: r2, rij
   double precision :: xyz(3, 1)

   integer :: nfound

   interface

      subroutine pgroup(tipo, i, event_str)

         integer, intent(in):: i, tipo
         character*50, optional :: event_str

      end subroutine pgroup

   end interface

   exgroup = 0

   r2 = out_bound
   r2 = r2*r2

   atomcycle: do iatom = 1, npart

      xyz(1, 1) = rx(iatom)
      xyz(2, 1) = ry(iatom)
      xyz(3, 1) = rz(iatom)

!calling the kdtree
      call kdtree2_r_nearest(tp=tree, qv=xyz(:, 1), r2=r2, nfound=nfound, nalloc=nres, results=results)

!Up to here we are working with the sequence of atoms,
!but in the following we have to work with the index in the group matrix.
!cycle for the found neighbours

      filter: do j = 2, nfound

! Getting the neighbours reference
         katom = results(j)%idx

! If katom is a reflected atom then it should be the original number
         katom = ref_fullxyz(katom)

         if (katom .eq. 0) then
            print *, 'something wrong with', i
            stop
         end if

         i = rl_part(iatom)
         k = rl_part(katom)

         if (i .eq. k) cycle filter ! filter

         mgi8 = group(i, 8)
         mgk8 = group(k, 8)

         rij = sqrt(results(j)%dis)

! If both molecules are part of the same system then jump
         if ((k .eq. mgi8) .and. (i .eq. mgk8)) cycle filter
!if (mgi8.eq.mgk8) cycle filter ! !filter

!Default: Normal collision
         iswitch = 1

!Se houver um atomo ja muito dentro fica colisao nao reactiva
         if (rij .lt. (out_bound*0.90d0)) then
            iswitch = -1
         end if

!Check if distance is smaller than mix group reacting distance

! if the iswitch equals -1 then the particles are to close to each other

         if (iswitch == -1) then

            if (mgi8 .ne. 0) then
               if (group(mgi8, 2) .eq. 4) cycle filter
            end if

            if (mgk8 .ne. 0) then
               if (group(mgk8, 2) .eq. 4) cycle filter
            end if

            group(i, 8) = 0
            group(k, 8) = 0

            if (mgi8 .ne. 0) group(mgi8, 8) = 0
            if (mgk8 .ne. 0) group(mgk8, 8) = 0

            call unit_switch(i, k)
            cycle filter

         end if

! if the iswitch equals zero then the particles are not colliding
         if (iswitch == 0) cycle filter ! filter

!########### DEFAULT->  IF ISWITCH==1 the particles are colliding

!If there are any molecules associated with i and k, this connection will be lost
         if ((mgi8 + mgk8) /= 0) then

!#### if there is already information of a non reactive ######################3
!#### collision with HO2 then jump particle ######################
            if (mgi8 .ne. 0) then
               if (group(mgi8, 2) .eq. 4) cycle filter
            end if

            if (mgk8 .ne. 0) then
               if (group(mgk8, 2) .eq. 4) cycle filter ! filter
            end if

            group(i, 8) = 0
            group(k, 8) = 0

            if (mgi8 .ne. 0) group(mgi8, 8) = 0
            if (mgk8 .ne. 0) group(mgk8, 8) = 0

         end if

! Begins the cases for collisions
         select case (group(i, 2))
         case (-3) ! Argon

               texto = 'N-REACT ->  Ar +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

            cycle filter

         case (-2) !Ox
            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'N-REACT -> O + Ar'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) !Ox + Ox

               if ((group(i, 7) + group(k, 7)) .eq. 6) then !if both Ox are O(3p)

                  call random_number(xrand)

                  if (xrand(1) .le. (3.0d0/81.0d0)) then
                     texto = 'REACT -> O(trip) + O(trip) '

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 2
                     group(i, 2) = 2
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(k, 3)
                     group(i, 5:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                  else

                     texto = 'N-REACT -> O(trip) + O(trip) '

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if

               else ! O(1D) + O(1D) or O(1D) + O(3P)

                  texto = 'N-REACT -> O(1D) + O(1D) '

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if
!!!
!!!
!!!  !      if ((xrand(1).gt.(3.0d0/81.0d0)).and.(xrand(1).lt.(8.0d0/81.0d0))) then
!!!  !
!!!  !
!!!  !
!!!  !
!!!  !
!!!  !      group(i,1)=2
!!!  !      group(i,2)=15
!!!  !      group(i,3)=group(i,3)
!!!  !      group(i,4)=group(k,3)
!!!  !      group(i,5:8)=0
!!!  !
!!!  !      group(k,:)=0
!!!  !
!!!  !
!!!  !
!!!  !      print*,'agrupa depois'
!!!  !
!!!  !      print*,'unidades'
!!!  !      print*,'i',i,'k',k
!!!  !      call printgroup
!!!  !
!!!  !
!!!  !      end if
!!!
!!!
!!!  !      cycle filter
!!!
!!!  !      end if
!!!
!!!
            case (-1) !Ox + Hy

               if (group(i, 7) .eq. 3) then

                  call random_number(xrand)

                  if (xrand(1) .le. (2.0d0/9.0d0)) then

                     write (unit=10, fmt='(a)') 'Ox + Hy -> OH'

                     texto = 'REACT -> Ox + Hy'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 2
                     group(i, 2) = 3
                     group(i, 4) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 5:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     !      for checking HO states
                     call write_diat_posi_velo(5, group(i, 3), group(i, 4))

                  else if ((xrand(1) .gt. (2.0d0/9.0d0)) .and. (xrand(1) .le. (6.0d0/9.0d0))) then

                     write (unit=10, fmt='(a)') 'Ox + Hy -> OH(4)'

                     texto = 'REACT -> Ox + Hy -> OH(4)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 2
                     group(i, 2) = 17
                     group(i, 4) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 5:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     !      for checking HO states
                     call write_diat_posi_velo(5, group(i, 3), group(i, 4))

                  else

                     texto = 'REACT -> Ox + Hy -> OH(4)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                  end if

               else

                  texto = 'N-REACT -> Ox + Hy'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

               end if

               cycle filter

            case (0) !Ox + nothing

               cycle filter

            case (1) !Ox + H2

               select case (group(i, 7))

               case (3) !Ox(3p) + H2

                  call random_number(xrand)

                  ! value altered for compensation the missing surface
                  !if (xrand(1).le.(1.0d0/3.0d0)) then            ! from the surface H2O it is possible to see the different channels
                  if (xrand(1) .le. (3.0d0/6.0d0)) then !It should be 1/3. There is another surface to be implemented. Here the probabilities are changed

                     texto = 'REACT -> O(trip) + H2'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 3
                     group(i, 2) = 6
                     group(i, 5) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 6:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> O(trip) + H2'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if

               case (1) !Ox(1d) + H2

                  call random_number(xrand)

                  if (xrand(1) .le. (1.0d0/5.0d0)) then

                     texto = 'REACT -> Ox(1d) + H2'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 3
                     group(i, 2) = 5
                     group(i, 5) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 6:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> Ox(1d) + H2'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if

               end select

            case (2) !Ox + O2

               if (group(i, 7) .eq. 3) then

                  call random_number(xrand)
                  if (xrand(1) .le. (1.0d0/27.0d0)) then

                     texto = 'REACT -> Ox + O2'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 3
                     group(i, 2) = 7
                     group(i, 5) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 6:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> Ox + O2'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if
               end if

               if (group(i, 7) .eq. 1) then

                  texto = 'N-REACT -> Ox + O2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (3) ! Ox + HO

               if (group(i, 7) .eq. 3) then

                  call random_number(xrand)
                  if (xrand(1) .le. (1.0d0/18.0d0)) then

                     !for checking HO states
                     call write_diat_posi_velo(0, group(k, 3), group(k, 4))

                     texto = 'REACT -> O(trip) + HO'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 3
                     group(i, 2) = 4
                     group(i, 5) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 6:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

!for checking HO states
                     call write_diat_posi_velo(2, group(k, 3), group(k, 4))

                     texto = 'N-REACT -> O(trip) + HO'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

!call printgroup

                     cycle filter

                  end if
               end if

               if (group(i, 7) .eq. 1) then

!for checking HO states
                  call write_diat_posi_velo(2, group(k, 3), group(k, 4))

                  texto = 'N-REACT -> O(sing) + HO'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

            case (4) ! Ox + HO2

               if (group(i, 7) .eq. 3) then
                  call random_number(xrand)
                  if (xrand(1) .le. (1.0d0/9.0d0)) then

                     texto = 'REACT ->  Ox + HO2'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 13
                     group(i, 6) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 5) = group(k, 5)
                     group(i, 7:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT ->  Ox + HO2'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if
               end if

               if (group(i, 7) .eq. 1) then

                  texto = 'N-REACT ->  Ox + HO2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

               cycle filter

            case (5) ! Ox +      H2O(1A')
               if (group(i, 7) .eq. 1) then  ! Ox(1D) +      H2O(~X1A')

                  call random_number(xrand)
                  if (xrand(1) .le. 0.30d0) then   !original was 0.20d0, but reactivity was to low. Other surfaces must be included later.
                     texto = 'REACT -> Ox + H2O(1A")'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 10              ! H2O2(1A)
                     group(i, 6) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 5) = group(k, 5)
                     group(i, 7:8) = 0
                     group(i, 9) = 302
                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> Ox + H2O(1A")'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if
               end if

               if (group(i, 7) .eq. 3) then      ! Ox(3P) +      H2O(~X1A')

                  texto = 'N-REACT -> Ox + H2O(1A")'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (6) ! Ox + H2O(3A'')

               if (group(i, 7) .eq. 1) then ! Ox(1D) + H2O(3A”)

                  texto = 'N-REACT -> Ox + H2O(3A")'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

               if (group(i, 7) .eq. 3) then  ! Ox(3P) + H2O(3A”)

                  call random_number(xrand)
                  if (xrand(1) .le. (3.0d0/54.0d0)) then  !Should be 1/27

                     texto = 'REACT -> Ox(3P) + H2O(3A")'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 10                  ! H2O2(1A)
                     group(i, 6) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 5) = group(k, 5)
                     group(i, 7:8) = 0
                     group(i, 9) = 301
                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> Ox(3P) + H2O(3A")'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if
               end if

            case (7) ! Ox + O3

               if (group(i, 7) .eq. 3) then
                  call random_number(xrand)
                  if (xrand(1) .le. (1.0d0/3.0d0)) then

                     texto = 'REACT ->  Ox + O3'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 11
                     group(i, 6) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 5) = group(k, 5)
                     group(i, 7:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT ->  Ox + O3'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter
                  end if

               end if
               if (group(i, 7) .eq. 1) then

                  texto = 'N-REACT ->  Ox + O3'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

               cycle filter

            case (8)  ! Ox + H3

               if (group(i, 7) .eq. 3) then

                  call random_number(xrand)
                  if (xrand(1) .le. (1.0d0/9.0d0)) then

                     texto = 'REACT ->  Ox + H3'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 12
                     group(i, 6) = group(i, 3)
                     group(i, 3) = group(k, 3)
                     group(i, 4) = group(k, 4)
                     group(i, 5) = group(k, 5)
                     group(i, 7:8) = 0
                     group(i, 9) = 1

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  end if
               end if

               if (group(i, 7) .eq. 1) then

                  texto = 'N-REACT ->  Ox + H3'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (15) !Ox + O2(1D)

               if (group(i, 7) .eq. 1) then      !Ox(1D) + O2(1D)

                  print *, 'N-REACT'
                  call unit_switch(i, k)
                  print *, 'unit_switch depois'
                  print *, 'unidades'
                  print *, 'i', i, 'k', k
                  !call printgroup

                  cycle filter
               end if

               if (group(i, 7) .eq. 3) then      !Ox(3P) + O2(1D)

                  texto = 'N-REACT -> Ox(3P) + O2(1D)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if
            case (9:13, 16:max_cases)

            case default

               texto = 'N-REACT -> O + Over 4 atoms'
               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter ! filter
            end select

         case (-1)      ! H
            select case (group(k, 2))

            case (-3) !  H + Ar

               texto = 'N-REACT -> H + Ar'
               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) !H + O
               if (group(k, 7) .eq. 3) then

                  call random_number(xrand)

                  if (xrand(1) .le. (2.0d0/9.0d0)) then

!for checking HO states
                     call write_diat_posi_velo(5, group(i, 3), group(k, 3))

                     texto = 'REACT -> H + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 2
                     group(i, 2) = 3
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(k, 3)
                     group(i, 5:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                  end if

                  if ((xrand(1) .gt. (2.0d0/9.0d0)) .and. (xrand(1) .le. (6.0d0/9.0d0))) then

                     texto = 'REACT ->  H + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 2
                     group(i, 2) = 17
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(k, 3)
                     group(i, 5:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                  end if

                  if (xrand(1) .ge. (6.0d0/9.0d0)) then

                     texto = 'N-REACT -> H + O(3p)'
                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter ! filter

                  end if

                  cycle filter

               end if

               if (group(k, 7) .eq. 1) then

                     texto = 'N-REACT -> H + O(1d)'
                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                 cycle filter

               end if

            case (-1) ! H + H

               call random_number(xrand)
               if (xrand(1) .le. 0.25d0) then

                  texto = 'REACT -> H + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)
                      
                  group(i, 1) = 2
                  group(i, 2) = 1
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(k, 3)
                  group(i, 5:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'REACT -> H + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)
 
                  group(i, 1) = 2
                  group(i, 2) = 14
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(k, 3)
                  group(i, 5:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               end if

            case (0)
               cycle filter
            case (1) ! H + H2

               texto = 'REACT -> H + H2'

               call pgroup(1, i, texto)
               call pgroup(1, k)
 
               group(i, 1) = 3
               group(i, 2) = 8
               group(i, 3) = group(i, 3)
               group(i, 4) = group(k, 3)
               group(i, 5) = group(k, 4)
               group(i, 6:8) = 0

               group(k, :) = 0

               call print_points
               call pgroup(1, i)
               call print_line

               cycle filter

            case (2) ! H + O2

               call random_number(xrand)

               if (xrand(1) .le. 1.0d0/3.0d0) then
                  texto = 'REACT -> H + O2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 3
                  group(i, 2) = 4
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(k, 3)
                  group(i, 5) = group(k, 4)
                  group(i, 6:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter
               else

                  texto = 'N-REACT -> H + O2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

!print*,'unit_switch depois'
!print*,'unidades'
!print*,'i',i,'k',k
!call printgroup

                  cycle filter ! filter

               end if

            case (3)      ! H + HO

               call random_number(xrand)
               if (xrand(1) .le. 0.125d0) then

!for checking HO states
                  call write_diat_posi_velo(0, group(k, 3), group(k, 4))

                  texto = 'REACT -> H + HO'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 3
                  group(i, 2) = 5
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(k, 3)
                  group(i, 5) = group(k, 4)
                  group(i, 6:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               end if

               if ((xrand(1) .gt. 0.125d0) .and. (xrand(1) .le. 0.5d0)) then  ! Should be 3/8. Percentage is altered for counting with the extra surfaces missing
!!!
!for checking HO states
                  call write_diat_posi_velo(0, group(k, 3), group(k, 4))

                  texto = 'REACT -> H + HO'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 3
                  group(i, 2) = 6
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(k, 3)
                  group(i, 5) = group(k, 4)
                  group(i, 6:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

!for checking HO states
                  call write_diat_posi_velo(2, group(k, 3), group(k, 4))

                  texto = 'N-REACT -> H + OH'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

!print*,'unit_switch depois'
!print*,'unidades'
!print*,'i',i,'k',k
!call printgroup

                  cycle filter ! filter

               end if

            case (4) !H + HO2

               call random_number(xrand)
               if (xrand(1) .lt. 0.375d0) then !  it shout be 0,25. It was changed to compensate the low reactivity of H2O2(1A)

                  texto = 'REACT -> H + HO2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 10
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(k, 3)
                  group(i, 5) = group(k, 4)
                  group(i, 6) = group(k, 5)
                  group(i, 7:8) = 0
                  group(i, 9) = 310
                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> H + HO2'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

!print*,'unit_switch depois'
!print*,'unidades'
!print*,'i',i,'k',k
!call printgroup

                  cycle filter ! filter

!                        group(i,1)=4
!                        group(i,2)=16         ! (H + HO2(2A''))* !! As to be correct to H2O2t case(16)
!                        group(i,3)=group(i,3)
!                        group(i,4)=group(k,3)
!                        group(i,5)=group(k,4)
!                        group(i,6)=group(k,5)
!                        group(i,7:8)=0
!
!                        group(k,:)=0
!
!
!
!                        print*,'agrupa depois'
!
!                        print*,'unidades'
!                        print*,'i',i,'k',k
!                        !call printgroup
!
!                        cycle filter
               end if

            case (5) ! H + H2O(1A')

               texto = 'REACT -> H + H2O(1A)'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               group(i, 1) = 4
               group(i, 2) = 12            !H3O
               group(i, 3) = group(i, 3)
               group(i, 4) = group(k, 3)
               group(i, 5) = group(k, 4)
               group(i, 6) = group(k, 5)
               group(i, 7:8) = 0
               group(i, 9) = 2

               group(k, :) = 0

               call print_points
               call pgroup(1, i)
               call print_line

               cycle filter

            case (6) ! H + H2O(3A'')

               call random_number(xrand)
               if (xrand(1) .le. (1.0d0/3.0d0)) then

                  texto = 'REACT -> H + H2O(1A")'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 12
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(k, 3)
                  group(i, 5) = group(k, 4)
                  group(i, 6) = group(k, 5)
                  group(i, 7:8) = 0
                  group(i, 9) = 3

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT ->  Ox + H2O(3A")'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

            case (7) ! H + O3
               texto = 'REACT ->  Ox + O3'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               group(i, 1) = 2
               group(i, 2) = 13
               group(i, 3) = group(i, 3)
               group(i, 4) = group(k, 3)
               group(i, 5) = group(k, 4)
               group(i, 6) = group(k, 5)
               group(i, 7:8) = 0

               group(k, :) = 0

               call print_points
               call pgroup(1, i)
               call print_line

               cycle filter

            case (8) ! H + H3

               call random_number(xrand)
               if (xrand(1) .lt. 0.25d0) then

                  texto = 'REACT ->  H + H3'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 9
                  group(i, 6) = group(i, 3)
                  group(i, 3) = group(k, 3)
                  group(i, 4) = group(k, 4)
                  group(i, 5) = group(k, 5)
                  group(i, 7:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT ->  H + H3'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (15) ! H + O2(1D)

               texto = 'N-REACT ->  H + O2(1D)'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter
            case (9:13, 16:max_cases)
            case default
               texto = 'N-REACT -> H + Over 4 atoms'
               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

!print*,'unit_switch depois'
!print*,'unidades'
!print*,'i',i,'k',k
!call printgroup

               cycle filter ! filter
            end select
         case (0)
            cycle filter

         case (1) ! H2
            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'N-REACT ->  Argon +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) ! H2 + Ox

               select case (group(k, 7))

               case (3) ! H2 + Ox(3p)

                  call random_number(xrand)

                  ! for compentation of the missing surface
                  !                       if (xrand(1).le.(1.0d0/3.0d0)) then            ! from the surface H2O it is possible to see the different channels

                  if (xrand(1) .le. (3.0d0/6.0d0)) then        !It should be 1/3. There is another surface to be implemented. Here the probabilities are changed  !

                     texto = 'REACT -> H2 + Ox(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 3
                     group(i, 2) = 6
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(k, 3)
                     group(i, 6:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> H2 + Ox(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if

               case (1) ! H2 + Ox(1d)

                  call random_number(xrand)

                  if (xrand(1) .le. (1.0d0/5.0d0)) then

                     texto = 'REACT -> H2 + Ox(1d)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 3
                     group(i, 2) = 5
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(k, 3)
                     group(i, 6:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> H2 + Ox(1d)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if

               end select

            case (-1)      !H2 + H

               texto = 'REACT -> H2 + H'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               group(i, 1) = 3
               group(i, 2) = 8
               group(i, 3) = group(i, 3)
               group(i, 4) = group(i, 4)
               group(i, 5) = group(k, 3)
               group(i, 6:8) = 0

               group(k, :) = 0

               call print_points
               call pgroup(1, i)
               call print_line

               cycle filter

            case (0)

               cycle filter
            case (1) ! H2 + H2

               texto = 'REACT -> H2 + H2'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               group(i, 1) = 4
               group(i, 2) = 9
               group(i, 3) = group(i, 3)
               group(i, 4) = group(i, 4)
               group(i, 5) = group(k, 3)
               group(i, 6) = group(k, 4)
               group(i, 7:8) = 0

               group(k, :) = 0

               call print_points
               call pgroup(1, i)
               call print_line

               cycle filter
            case (2) ! H2 + O2(3Sigma)

               texto = 'N-REACT ->  Ox + H3'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

               !group(i,1)=4
               !group(i,2)=16     !      (O2 (3Σ) +H2)*  has to be corrected to H2O2t (case(16))
               !group(i,3)=group(i,3)
               !group(i,4)=group(i,4)
               !group(i,5)=group(k,3)
               !group(i,6)=group(k,4)
               !group(i,7:8)=0
               !
               !group(k,:)=0
               !
               !print*,'agrupa depois'
               !
               !print*,'unidades'
               !print*,'i',i,'k',k
               !
               !cycle filter

            case (15) ! H2 + O2(1Delta)

               call random_number(xrand)

               if (xrand(1) .le. (0.75d0)) then !it should be 1/2

                  texto = 'REACT -> H2 + O2(1Delta)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 10
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(i, 4)
                  group(i, 5) = group(k, 3)
                  group(i, 6) = group(k, 4)
                  group(i, 7:8) = 0
                  group(i, 9) = 340

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> H2 + O2(1Delta)'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (3) !H2 + OH

               call random_number(xrand)

               if (xrand(1) .le. (0.5d0)) then

                  texto = 'REACT ->  H2 + OH'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

!for checking HO states
                  call write_diat_posi_velo(0, group(k, 3), group(k, 4))

                  group(i, 1) = 4
                  group(i, 2) = 12
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(i, 4)
                  group(i, 5) = group(k, 3)
                  group(i, 6) = group(k, 4)
                  group(i, 7:8) = 0
                  group(i, 9) = 4

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> H2 + OH'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (4:13, 16:max_cases)

               texto = 'N-REACT ->  H2 +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter
            end select

         case (2)
            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'N-REACT ->  O2(3Sigma) + Ar'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2)      !O2(3Sigma) + O

               if (group(k, 7) .eq. 3) then

                  call random_number(xrand)
                  if (xrand(1) .le. (1.0d0/27.0d0)) then

                     texto = 'REACT -> O2(3Sigma) + O'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 3
                     group(i, 2) = 7
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(k, 3)
                     group(i, 6:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> O2(3Sigma) + O(trip)'
                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

!print*,'unit_switch depois'
!print*,'unidades'
!print*,'i',i,'k',k
!call printgroup

                     cycle filter ! filter

                  end if
               end if

               if (group(k, 7) .eq. 1) then

                  texto = 'N-REACT -> O2(1Delta) + O(trip)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

            case (-1) !O2(3Sigma) + H

               call random_number(xrand)

               if (xrand(1) .le. (1.0d0/3.0d0)) then

                  texto = 'REACT -> O2(3Sigma) + H'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 3
                  group(i, 2) = 4
                  group(i, 5) = group(i, 4)
                  group(i, 4) = group(i, 3)
                  group(i, 3) = group(k, 3)
                  group(i, 6:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> O2(3Sigma) + H'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

!print*,'unit_switch depois'
!print*,'unidades'
!print*,'i',i,'k',k
!call printgroup

                  cycle filter ! filter

               end if

            case (0)

               cycle filter
            case (1) !O2(3Sigma) + H2

               texto = 'N-REACT -> O2(3Sigma) + H2'
               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

               !
               !                        group(i,1)=4
               !                        group(i,2)=16 !    (O2 (3Σ) +H2)*      As to be corrected to H2O2t case(16)
               !                        group(i,5)=group(i,3)
               !                        group(i,6)=group(i,4)
               !                        group(i,3)=group(k,3)
               !                        group(i,4)=group(k,4)
               !                        group(i,7:8)=0
               !
               !                        group(k,:)=0
               !
               !                        print*,'agrupa depois'
               !
               !                        print*,'unidades'
               !                        print*,'i',i,'k',k
               !                        !call printgroup
               !
               !                        cycle filter
            case (2) !O2(3Sigma) + O2(3Sigma)

               call random_number(xrand)

               if (xrand(1) .le. (1.0d0/3.0d0)) then

                  texto = 'REACT -> O2(3Sigma) + O2(3Sigma)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 11
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(i, 4)
                  group(i, 5) = group(k, 3)
                  group(i, 6) = group(k, 4)
                  group(i, 7:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> O2(3Sigma) + O2(3Sigma)'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter ! filter

               end if

            case (3) !O2(3Sigma) + HO

               call random_number(xrand)
               if (xrand(1) .le. (1.0d0/6.0d0)) then

                  texto = 'REACT -> O2(3Sigma) + HO'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

!for checking HO states
                  call write_diat_posi_velo(0, group(k, 3), group(k, 4))

                  group(i, 1) = 4
                  group(i, 2) = 13
                  group(i, 5) = group(i, 3)
                  group(i, 6) = group(i, 4)
                  group(i, 3) = group(k, 3)
                  group(i, 4) = group(k, 4)
                  group(i, 7:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

!for checking HO states
                  call write_diat_posi_velo(2, group(k, 3), group(k, 4))

                  texto = 'N-REACT -> O2(3Sigma) + OH'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  !print*,'unit_switch depois'
                  !print*,'unidades'
                  !print*,'i',i,'k',k
                  !call printgroup

                  cycle filter ! filter
               end if

            case (15) !O2(3Sigma) + O2(1Delta)

               texto = 'N-REACT -> O2(3Sigma) + O2(1Delta)'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (4:13, 16:max_cases)

               texto = 'N-REACT -> O2(3Sigma) +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

!print*,'unit_switch depois'
!print*,'unidades'
!print*,'i',i,'k',k
!call printgroup

               cycle filter ! filter
            end select

         case (15)
            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'REACT -> O2(1D) + Ar'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               !call printgroup

               cycle filter

            case (-2)      !O2(1D) + O

               if (group(k, 7) .eq. 1) then      !O2(1D) + Ox(sing)

                  texto = 'N-REACT ->  O2(1D) + Ox'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

               if (group(k, 7) .eq. 3) then      !O2(1D) + Ox(trip)

                  texto = 'N-REACT -> O2(1D) + Ox(trip)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

               if (group(k, 7) .eq. 0) then
                  print *, k
                  stop 'no information about sing. or trip. state'
               end if

               cycle filter

            case (-1) !O2(1D) + H

               texto = 'N-REACT -> O2(1D) + H'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (0)

               cycle filter
            case (1) !O2(1D) + H2

               call random_number(xrand)

               if (xrand(1) .le. (3.0d0/4.0d0)) then

                  texto = 'REACT -> O2(1D) + H2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 10
                  group(i, 5) = group(i, 3)
                  group(i, 6) = group(i, 4)
                  group(i, 3) = group(k, 3)
                  group(i, 4) = group(k, 4)
                  group(i, 7:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> O2(1D) + H2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

            case (2) !O2(1Delta) + O2(3Sigma)

               texto = 'N-REACT ->  O2(1Delta) + O2(3Sigma)'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (3) !O2(1D) + HO

!for checking HO states
               call write_diat_posi_velo(2, group(k, 3), group(k, 4))

               texto = 'N-REACT -> O2(1D) + HO'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (15) !O2(1D) + O2(1D)

               texto = 'N-REACT -> O2(1D) + O2(1D)'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (4:13, 16:max_cases)

               texto = 'N-REACT -> O2(1D) +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            end select

         case (3)      ! HO
            select case (group(k, 2))

            case (-3) ! Argon

               !for checking HO states
               call write_diat_posi_velo(2, group(i, 3), group(i, 4))

               texto = 'N-REACT ->  HO + Ar'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) ! HO + O

               if (group(k, 7) .eq. 3) then

                  call random_number(xrand)
                  if (xrand(1) .le. (1.0d0/18.0d0)) then

!for checking HO states
                     call write_diat_posi_velo(0, group(i, 3), group(i, 4))

                     texto = 'REACT -> HO + O(trip)'
                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 3
                     group(i, 2) = 4      !HO2(2A'')
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(k, 3)
                     group(i, 6) = 0
                     group(i, 7:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

!for checking HO states
                     call write_diat_posi_velo(2, group(i, 3), group(i, 4))

                     texto = 'N-REACT -> HO + O(trip)'
                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     !print*,'unit_switch depois'
                     !print*,'unidades'
                     !print*,'i',i,'k',k
                     !call printgroup

                     cycle filter ! filter

                  end if
               end if

               if (group(k, 7) .eq. 1) then

!for checking HO states
                  call write_diat_posi_velo(2, group(i, 3), group(i, 4))

                  texto = 'N-REACT -> HO + O(sing)'
                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  !print*,'unit_switch depois'
                  !print*,'unidades'
                  !print*,'i',i,'k',k
                  !call printgroup

                  cycle filter ! filter

               end if

            case (-1) !HO + H

               call random_number(xrand)
               if (xrand(1) .le. (1.0d0/8.0d0)) then

                  !for checking HO states
                  call write_diat_posi_velo(0, group(i, 3), group(i, 4))

                  texto = 'REACT -> HO + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 3
                  group(i, 2) = 5
                  group(i, 5) = group(i, 4)
                  group(i, 4) = group(i, 3)
                  group(i, 3) = group(k, 3)
                  group(i, 6:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else if ((xrand(1) .gt. (1.0d0/8.0d0)) .and. (xrand(1) .le. (4.0d0/8.0d0))) then

                  !for checking HO states
                  call write_diat_posi_velo(2, group(i, 3), group(i, 4))

                  texto = 'REACT -> HO + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 3
                  group(i, 2) = 6
                  group(i, 5) = group(i, 4)
                  group(i, 4) = group(i, 3)
                  group(i, 3) = group(k, 3)
                  group(i, 6:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  !for checking HO states
                  call write_diat_posi_velo(2, group(i, 3), group(i, 4))

                  texto = 'N-REACT ->  HO + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (0)

               cycle filter
            case (1)       !HO + H2

               call random_number(xrand)

               if (xrand(1) .le. (1.0d0/2.0d0)) then

!for checking HO states
                  call write_diat_posi_velo(0, group(i, 3), group(i, 4))

                  texto = 'REACT -> HO + H2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 12
                  group(i, 5) = group(i, 3)
                  group(i, 6) = group(i, 4)
                  group(i, 3) = group(k, 3)
                  group(i, 4) = group(k, 4)
                  group(i, 7:8) = 0
                  group(i, 9) = 5

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

!for checking HO states
                  call write_diat_posi_velo(2, group(i, 3), group(i, 4))

                  texto = 'N-REACT ->  HO + H2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

            case (2)      !HO + O2

               call random_number(xrand)

               if (xrand(1) .le. (1.0d0/6.0d0)) then

                  !for checking HO states
                  call write_diat_posi_velo(0, group(i, 3), group(i, 4))

                  texto = 'REACT ->  HO + O2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 13
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(i, 4)
                  group(i, 5) = group(k, 3)
                  group(i, 6) = group(k, 4)
                  group(i, 7:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  !for checking HO states
                  call write_diat_posi_velo(2, group(i, 3), group(i, 4))

                  texto = 'N-REACT ->  HO + O2'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

            case (3) !HO + HO

               call random_number(xrand)
               if (xrand(1) .lt. 0.09375d0) then !it should be 1/16

                  !for checking HO states
                  call write_diat_posi_velo(0, group(i, 3), group(i, 4))

                  texto = 'REACT ->  HO + HO'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 10
                  group(i, 3) = group(i, 3)
                  group(i, 5) = group(i, 4)
                  group(i, 4) = group(k, 3)
                  group(i, 6) = group(k, 4)
                  group(i, 7:8) = 0
                  group(i, 9) = 320

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  !for checking HO states
                  call write_diat_posi_velo(2, group(i, 3), group(i, 4))

                  texto = 'N-REACT ->  HO + HO'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (15)      !HO + O2(1Delta)

               !for checking HO states
               call write_diat_posi_velo(2, group(i, 3), group(i, 4))

               texto = 'N-REACT -> HO + O2(1D)'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (4:13, 16:max_cases)

!for checking HO states
               call write_diat_posi_velo(2, group(i, 3), group(i, 4))

               texto = 'N-REACT -> HO +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            end select

         case (4)      !HO2

            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'N-REACT ->  HO2(2A") + Ar'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) !HO2(2A") + O

               if (group(k, 7) .eq. 3) then

                  call random_number(xrand)
                  if (xrand(1) .le. (1.0d0/9.0d0)) then

                     texto = 'REACT ->  HO2(2A") + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 13            !HO3(2A)
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(i, 5)
                     group(i, 6) = group(k, 3)
                     group(i, 7:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter
                  else

                     texto = 'N-REACT ->  HO2(2A") + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if
               end if

               if (group(k, 7) .eq. 1) then

                  texto = 'N-REACT ->  HO2(2A") + O(1d)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (-1) !HO2 + H

               call random_number(xrand)
               if (xrand(1) .lt. 0.375d0) then

                  texto = 'REACT -> HO2 + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 10
                  group(i, 6) = group(i, 5)
                  group(i, 5) = group(i, 4)
                  group(i, 4) = group(i, 3)
                  group(i, 3) = group(k, 3)
                  group(i, 7:8) = 0
                  group(i, 9) = 310

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> HO2 + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter ! filter
               end if

            case (0)
               cycle filter

            case (1:max_cases)

               texto = 'N-REACT ->  HO2 +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

!      for checking HO states
               if (group(k, 2) .eq. 3) call write_diat_posi_velo(2, group(k, 3), group(k, 4))

!call printgroup
               cycle filter ! filter
            end select

         case (5) ! H2O(~X1A')

            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'N-REACT ->  H2O(~X1A) + Ar'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) !H2O(~X1A') + O

               if (group(k, 7) .eq. 3) then       !H2O(~X1A') + O (3P)

                  texto = 'N-REACT ->  H2O(~X1A) + O(3p)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

!                        group(i,1)=4
!                        group(i,2)=16   !  (H2O(~X1A') +O(3P))*  has to be corrected to H2O2t case(16)
!                        group(i,3)=group(i,3)
!                        group(i,4)=group(i,4)
!                        group(i,5)=group(i,5)
!                        group(i,6)=group(k,3)
!                        group(k,7:8)=0
!
!                        group(k,:)=0
!
!                        print*,'agrupa depois'
!
!                        print*,'unidades'
!                        print*,'i',i,'k',k
!                        !call printgroup
!
!                        cycle filter

               end if

               if (group(k, 7) .eq. 1) then      !H2O(~X1A') + O (1D) !original was 0.20d0, but reactivity was to low. Other surfaces must be included later.

                  call random_number(xrand)
                  if (xrand(1) .lt. 0.3d0) then

                     texto = 'REACT ->  H2O(~X1A) + O(1d)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 10                  !H2O2(1A)
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(i, 5)
                     group(i, 6) = group(k, 3)
                     group(i, 7:8) = 0
                     group(i, 9) = 302
                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT ->  H2O(~X1A) + O(1d)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if

                  cycle filter
               end if

            case (-1) !H2O(~X1A') + H

               !                        call random_number(xrand)
               !                        if (xrand(1).le.(1.0d0/3.0d0)) then

               texto = 'REACT -> H2O(~X1A) + H'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               group(i, 1) = 4
               group(i, 2) = 12            !H3O
               group(i, 6) = group(i, 5)
               group(i, 5) = group(i, 4)
               group(i, 4) = group(i, 3)
               group(i, 3) = group(k, 3)
               group(i, 7:8) = 0
               group(i, 9) = 2

               group(k, :) = 0

               call print_points
               call pgroup(1, i)
               call print_line
!call printgroup

               cycle filter
               !                        else
               !
               !
               !                        print*,'N-REACT'
               !                        call unit_switch(i,k)
               !                        print*,'unit_switch depois'
               !                        print*,'unidades'
               !                        print*,'i',i,'k',k
               !                        !call printgroup

               !                        cycle filter
               !
               !                        end if

            case (0)
               cycle filter

            case (1:max_cases)

               texto = 'N-REACT -> HO2 +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

!      for checking HO states
               if (group(k, 2) .eq. 3) call write_diat_posi_velo(2, group(k, 3), group(k, 4))

               !call printgroup
               cycle filter
            end select
         case (6)      ! H2O(3A)

            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'REACT -> HO2 + H'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) ! H2O(3A) + O
               if (group(k, 7) .eq. 3) then      ! H2O(3A) + O(3P)

                  call random_number(xrand)

                  if (xrand(1) .le. (3.0d0/54.0d0)) then  !Should be 1/27  !

                     texto = 'REACT -> HO2 + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 10                  !H2O2(1A)
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(i, 5)
                     group(i, 6) = group(k, 3)
                     group(i, 7:8) = 0
                     group(i, 9) = 301

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> HO2 + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if
               end if

               if (group(k, 7) .eq. 1) then

                  texto = 'N-REACT -> HO2 + O(1d)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

            case (-1) ! H2O(3A) + H

               call random_number(xrand)
               if (xrand(1) .le. (1.0d0/3.0d0)) then

                  texto = 'REACT -> H2O(3A) + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 12            !H3O
                  group(i, 6) = group(i, 5)
                  group(i, 5) = group(i, 4)
                  group(i, 4) = group(i, 3)
                  group(i, 3) = group(k, 3)
                  group(i, 7:8) = 0
                  group(i, 9) = 3

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> H2O(3A) + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter
               end if

            case (0)
               cycle filter

            case (1:max_cases)

               texto = 'N-REACT -> H2O(3A) +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

!      for checking HO states
               if (group(k, 2) .eq. 3) call write_diat_posi_velo(2, group(k, 3), group(k, 4))
               !call printgroup
               cycle filter
            end select

         case (7)      ! O3
            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'N-REACT -> O3 + Ar'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) ! O3 + O(3P)

               if (group(k, 7) .eq. 3) then
                  call random_number(xrand)

                  if (xrand(1) .le. (1.0d0/3.0d0)) then

                     texto = 'REACT -> O3 + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 11            !O4
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(i, 5)
                     group(i, 6) = group(k, 3)
                     group(i, 7:8) = 0

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> O3 + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if
               end if

               if (group(k, 7) .eq. 1) then

                  texto = 'N-REACT -> O3 + O(1d)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (-1) ! O3 + H

               texto = 'REACT -> O3 + H'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               group(i, 1) = 4
               group(i, 2) = 13      !HO3(2A)
               group(i, 6) = group(i, 5)
               group(i, 5) = group(i, 4)
               group(i, 4) = group(i, 3)
               group(i, 3) = group(k, 3)
               group(i, 7:8) = 0

               group(k, :) = 0

               call print_points
               call pgroup(1, i)
               call print_line

               cycle filter

            case (0)
               cycle filter

            case (1:max_cases)

               texto = 'N-REACT -> O3 +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

!      for checking HO states
               if (group(k, 2) .eq. 3) call write_diat_posi_velo(2, group(k, 3), group(k, 4))
               !call printgroup
               cycle filter
            end select

         case (8) !H3
            select case (group(k, 2))

            case (-3) ! Argon

               texto = 'N-REACT -> H3 + Ar'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

               cycle filter

            case (-2) !H3 + O

               if (group(k, 7) .eq. 3) then

                  call random_number(xrand)
                  !                        if (xrand(1).le.(1.0d0/9.0d0)) then
                  if (xrand(1) .le. (1.0d0)) then

                     texto = 'REACT -> H3 + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     group(i, 1) = 4
                     group(i, 2) = 12
                     group(i, 3) = group(i, 3)
                     group(i, 4) = group(i, 4)
                     group(i, 5) = group(i, 5)
                     group(i, 6) = group(k, 3)
                     group(i, 7:8) = 0
                     group(i, 9) = 1

                     group(k, :) = 0

                     call print_points
                     call pgroup(1, i)
                     call print_line

                     cycle filter

                  else

                     texto = 'N-REACT -> H3 + O(3p)'

                     call pgroup(1, i, texto)
                     call pgroup(1, k)

                     call unit_switch(i, k)

                     call print_points
                     call pgroup(1, i)
                     call pgroup(1, k)
                     call print_line

                     cycle filter

                  end if
               end if

               if (group(k, 7) .eq. 1) then

                  texto = 'N-REACT -> H3 + O(1d)'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if

            case (-1) ! H3 + H

               call random_number(xrand)
               if (xrand(1) .lt. 0.25d0) then

                  texto = 'REACT -> H3 + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  group(i, 1) = 4
                  group(i, 2) = 9
                  group(i, 3) = group(i, 3)
                  group(i, 4) = group(i, 4)
                  group(i, 5) = group(i, 5)
                  group(i, 6) = group(k, 3)
                  group(i, 7:8) = 0

                  group(k, :) = 0

                  call print_points
                  call pgroup(1, i)
                  call print_line

                  cycle filter

               else

                  texto = 'N-REACT -> H3 + H'

                  call pgroup(1, i, texto)
                  call pgroup(1, k)

                  call unit_switch(i, k)

                  call print_points
                  call pgroup(1, i)
                  call pgroup(1, k)
                  call print_line

                  cycle filter

               end if
            case (0)
               cycle filter

            case (1:max_cases)

               texto = 'N-REACT -> H3 +'

               call pgroup(1, i, texto)
               call pgroup(1, k)

               call unit_switch(i, k)

               call print_points
               call pgroup(1, i)
               call pgroup(1, k)
               call print_line

!      for checking HO states
               if (group(k, 2) .eq. 3) call write_diat_posi_velo(2, group(k, 3), group(k, 4))
!call printgroup
               cycle filter

            end select

         case (9:13, 16:max_cases)

            texto = 'N-REACT -> Over 4 atoms'

            call pgroup(1, i, texto)
            call pgroup(1, k)

            call unit_switch(i, k)

            call print_points
            call pgroup(1, i)
            call pgroup(1, k)
            call print_line

            cycle filter
         end select

!$OMP PARALLEL DO
         do k = 1, npart
            rl_part(group(k, 3)) = k
            rl_part(group(k, 4)) = k
            rl_part(group(k, 5)) = k
            rl_part(group(k, 6)) = k
         end do
!$OMP END PARALLEL DO

      end do filter

   end do atomcycle

   call kdtree2_destroy(tree)
   deallocate (ref_fullxyz)

!!!        end do filter
!!!    end do search
!!!   end do filter

   exgroup = 0

end subroutine agrupa

subroutine unit_switch(iin, kin)
   use variables
   use phys_parameters
   use sim_variables
   implicit none

   integer, intent(in) :: iin, kin

   group(iin, 8) = kin
   group(kin, 8) = iin

end subroutine unit_switch

subroutine rswitch(iin, kin, switchout)
   use variables
   use phys_parameters
   use sim_variables
   use constants
   implicit none
   integer, intent(in) :: iin, kin
   double precision      :: rxij, ryij, rzij, rij
   integer, intent(out) ::switchout
   integer      :: ni, nj, gi, gk

   switchout = 0
!        print*,tempo
!        !call printgroup
!            print*,'iin,kin',iin,kin

   test: do ni = 3, 2 + group(iin, 1)

!           print*,'ni,2+group(iin,1)',ni,2+group(iin,1)

      do nj = 3, 2 + group(kin, 1)
!            print*,'ni,2+group(iin,1)',ni,2+group(iin,1)

         gi = group(iin, ni)
         gk = group(kin, nj)
!            print*,gi,gk

         call mod_dist(rx(gi), ry(gi), rz(gi), rx(gk), ry(gk), rz(gk), rxij, ryij, rzij, rij)

!Verifica distancias em relacao a distancia de aglumeracao

         if (rij .lt. out_bound) then
            switchout = 1
!Se houver um atomo ja muito dentro fica colisao nao reactiva
            if (rij .lt. (out_bound*0.90d0)) then
               switchout = -1
               exit test
            end if
         end if

      end do
   end do test
end subroutine rswitch

