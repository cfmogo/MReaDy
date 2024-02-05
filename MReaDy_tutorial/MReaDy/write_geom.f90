SUBROUTINE write_geom
   use constants
   use variables
   use phys_parameters
   use sim_variables
   implicit none
   
   integer :: i
  

!   reactivas
   if ((tempo) .eq. -1) then
      open (unit=30, file='out/geo'//trim(bk_new_string)//'.xyz', status='replace')
      tempo = 0
   else
      open (unit=30, file='out/geo'//trim(bk_new_string)//'.xyz', position="APPEND")
   end if
   write (30, *) npart

   write (unit=30, fmt='(A5,I10,F11.2,A1,3(A4,E16.8))') "Step", tempo + old_step_print, &
      total_real_time/fentosec, "fs ", "V=", potential/con, "T=", ekin/con, &
      "U=", (potential + ekin)/con

   writting: do i = 1, npart

      write (unit=30, fmt='(A3,4F20.12,I8)') cname(i), (rx(i))/rnm, (ry(i))/rnm, &
         (rz(i))/rnm, sqrt((vx(i))**2 + (vy(i))**2 + (vz(i))**2), i

   end do writting

   close (30)

END SUBROUTINE write_geom

subroutine print_line
   write (*, *) '-------------------------------------------------------------------------------'
end subroutine print_line

subroutine print_points
   write (*, *) '...............................................................................'
end subroutine print_points

subroutine printgroup
   use variables
   use constants
   implicit none

   integer      ::      i, i2, i3, i4, i5, i6, i7, i8
   character*24      ::str, uname(-2:max_cases)
   data uname/'O', 'H', ' ', 'H2s', 'O2(3S)', 'OH(2Pi)', 'HO2', 'H2O(1A`)', &
      'H2O(3A``)', 'O3', 'H3', 'H4', 'H2O2(1A)', 'O4', 'H3O', 'HO3', 'H2t', &
      'O2(1D)', 'H2O2(3A)', 'HO(4)'/

   character*3      :: str3, str4, str5, str6, str7, str8
   write (*, *) 'Step', tempo + old_step_print
   write (*, fmt='(A15,F12.2,A4)') 'Total time ', total_real_time/fentosec, 'fs'
   write (*, *) '-------------------------------------------------------------------------------'

   writting: do i = 1, npart
      if (group(i, 2) == 0) cycle

      i2 = group(i, 2)
      i3 = group(i, 3)
      i4 = group(i, 4)
      i5 = group(i, 5)
      i6 = group(i, 6)
      i7 = group(i, 7)
      i8 = group(i, 8)
      str = ' '//(uname(i2))
      str3 = ' '//cname(i3)
      str4 = ' '//cname(i4)
      str5 = ' '//cname(i5)
      str6 = ' '//cname(i6)
      str7 = ' '//cname(i7)
      str8 = ' '//cname(i8)
      write (*, fmt='(A1,I4,I2,I6,A25,4(I5,A3),I5,I5,A1)') '|', i, group(i, 1:2), &
         str, i3, str3, i4, str4, i5, str5, i6, str6, group(i, 7:8), '|'
   end do writting
   write (*, *) '-------------------------------------------------------------------------------'
   write (*, *) ''
end subroutine printgroup

!********** for printing part of the group/group2/exgroup matrix which changed *************************

subroutine pgroup(tipo, i, event_str)
   use variables
   use phys_parameters
   use sim_variables
   use constants
   implicit none
   integer, intent(in):: i, tipo
   character*50, optional :: event_str
   integer      ::      i2, i3, i4, i5, i6, i7, i8
   character*24      ::str, uname(-2:max_cases)
   data uname/'O', 'H', ' ', 'H2s', 'O2(3S)', 'OH(2Pi)', 'HO2', 'H2O(1A`)', &
      'H2O(3A``)', 'O3', 'H3', 'H4', 'H2O2(1A)', 'O4', 'H3O', 'HO3', 'H2t', &
      'O2(1D)', 'H2O2(3A)', 'HO(4)'/

   character*3      :: str3, str4, str5, str6, str7, str8

      if (print_flag .eq. 0) then
         write (*, fmt='(A6,I15,A15,F12.2,A4)') 'Step', tempo + old_step_print, 'Total time ', total_real_time/fentosec, 'fs'
         print_flag = 1
      end if


   if (present(event_str)) write (*, *) event_str

   write (*, *) '-------------------------------------------------------------------------------'

   if (tipo .eq. 1) then

      i2 = group(i, 2)
      i3 = group(i, 3)
      i4 = group(i, 4)
      i5 = group(i, 5)
      i6 = group(i, 6)
      i7 = group(i, 7)
      i8 = group(i, 8)
      str = ' '//(uname(i2))
      str3 = ' '//cname(i3)
      str4 = ' '//cname(i4)
      str5 = ' '//cname(i5)
      str6 = ' '//cname(i6)
      str7 = ' '//cname(i7)
      str8 = ' '//cname(i8)

      write (*, fmt='(A1,I4,I2,I6,A25,4(I5,A3),I5,I5,A1)') '|', i, group(i, 1:2), &
         str, i3, str3, i4, str4, i5, str5, i6, str6, group(i, 7:8), '|'

   else if (tipo .eq. 2) then

      i2 = group2(i, 2)
      i3 = group2(i, 3)
      i4 = group2(i, 4)
      i5 = group2(i, 5)
      i6 = group2(i, 6)
      i7 = group2(i, 7)
      i8 = group2(i, 8)
      str = ' '//(uname(i2))
      str3 = ' '//cname(i3)
      str4 = ' '//cname(i4)
      str5 = ' '//cname(i5)
      str6 = ' '//cname(i6)
      str7 = ' '//cname(i7)
      str8 = ' '//cname(i8)

      write (*, fmt='(A1,I4,I2,I6,A25,4(I5,A3),I5,I5,A1)') '|', i, group2(i, 1:2), &
         str, i3, str3, i4, str4, i5, str5, i6, str6, group2(i, 7:8), '|'

   else if (tipo .eq. 3) then

      i2 = exgroup(i, 2)
      i3 = exgroup(i, 3)
      i4 = exgroup(i, 4)
      i5 = exgroup(i, 5)
      i6 = exgroup(i, 6)
      i7 = exgroup(i, 7)
      i8 = exgroup(i, 8)
      str = ' '//(uname(i2))
      str3 = ' '//cname(i3)
      str4 = ' '//cname(i4)
      str5 = ' '//cname(i5)
      str6 = ' '//cname(i6)
      str7 = ' '//cname(i7)
      str8 = ' '//cname(i8)

      write (*, fmt='(A1,I4,I2,I6,A25,4(I5,A3),I5,I5,A1)') '|', i, exgroup(i, 1:2), &
         str, i3, str3, i4, str4, i5, str5, i6, str6, exgroup(i, 7:8), '|'
   end if

end subroutine pgroup

subroutine write_bk()
   use constants
   use variables
   use phys_parameters
   use sim_variables
   implicit none
  

   open (unit=41, file='out/geo'//trim(bk_new_string)//'.bk', position="APPEND", access="stream")
   write (unit=41) tempo + old_step_print
   write (unit=41) total_real_time
   write (unit=41) group

   write (unit=41) rx, ry, rz
   write (unit=41) vx, vy, vz
   write (unit=41) ax, ay, az

   write (unit=41) ekin
!!
   write (unit=41) potential

   call random_seed(get=iget)
   write (unit=41) iget

   close (41)

end subroutine write_bk

subroutine time_change ! for follow the reaction changes through time
   use constants
   use variables
   use phys_parameters
   use sim_variables
   implicit none
   integer      ::      i, countcases(-2:max_cases)
   character*24      ::uname(-2:max_cases)
   character*30  :: formato

   data uname/'O', 'H', 'Empty', 'H2s', 'O2(3S)', 'OH(2Pi)', 'HO2', 'H2O(1A`)', &
      'H2O(3A``)', 'O3', 'H3', 'H4', 'H2O2(1A)', 'O4', 'H3O', 'HO3', 'H2t', &
      'O2(1D)', 'H2O2(3A)', 'HO(4)'/
   countcases = 0
!$OMP parallel do shared(group) reduction(+:countcases)
   do i = 1, npart
!                  print*,'i=i',i
!                  print*,group(i,2)
!                  print*,countcases(group(i,2))
      countcases(group(i, 2)) = countcases(group(i, 2)) + 1
!                  print*,'i=',i
!                  print*,group(i,2)
!                  print*,countcases(group(i,2))
      countcases(0) = 0
   end do
!$OMP end parallel do
   print *, 'tempo', tempo
   if ((tempo) .eq. 0) then
      open (unit=51, file='out/time_changes'//trim(bk_new_string)//'.cvs', status='replace')

      write (51, '(2(a13),20(a9,a1))') 'Steps;', 'Time/fs;', (adjustr(trim(uname(i))), ';', i=-2, max_cases)
!            write(51,*)old_step_print*dt/fentosec,';', (countcases(i),';',i=-2,max_cases)

      write (51, '(i13,a1,f12.2,a1,20(i9,a1))') old_step_print, ';', total_real_time/fentosec, ';', &
         (countcases(i), ';', i=-2, max_cases)
   else
      open (unit=51, file='out/time_changes'//trim(bk_new_string)//'.cvs', position="APPEND")
      write (formato, '(a17,i5,a8)') '(i13,a1,f12.2,a1,', max_cases + 3, '(i9,a1))'

      write (51, formato) tempo + old_step_print, ';', total_real_time/fentosec, ';', &
         (countcases(i), ';', i=-2, max_cases)
   end if
   close (51)
end subroutine


subroutine write_diat_posi_velo(state, atom1, atom2)

   use variables
   use constants
   integer, intent(in)   ::  atom1, atom2, state

   open (unit=80, file='out/diat_state_'//trim(bk_new_string)//'.dat', position="APPEND")
   write (80, *) '-----------------------------------------------------------'
   write (80, fmt='(a53,i3)') 'Dest(0) Creat(1)/Coll(3)/Desintegrate(4)/Integ(5)->', state
   write (80, *) tempo + old_step_print, total_real_time/fentosec
   write (80, fmt='(a6,i4,a6,i4)') 'atoms ', atom1, ' and ', atom2
   write (80, *) 'position and velocity for H atom'
!      write(80,*) 'atom 1 type= ',cname(atom1)
   write (80, fmt='(3D20.12)') rx(atom1), ry(atom1), rz(atom1)
   write (80, fmt='(3D20.12)') vx(atom1), vy(atom1), vz(atom1)
!      write(80,*) 'atom 2 type= ',cname(atom2)
   write (80, *) 'position and velocity for O atom'
   write (80, fmt='(3D20.12)') rx(atom2), ry(atom2), rz(atom2)
   write (80, fmt='(3D20.12)') vx(atom2), vy(atom2), vz(atom2)

   close (80)

end subroutine

