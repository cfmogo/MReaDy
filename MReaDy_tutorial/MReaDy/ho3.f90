       subroutine VHO3_13_switch(i,V,dv)
       use variables
       use constants
       double precision, dimension(12) :: dV1,dV
       integer      :: gi3,gi4,gi5,gi6
       integer      :: i!,j,k,l
!       double precision, dimension(3) :: acc,bcc,ccc,dcc
!       double precision :: dist
       double precision :: rxij,ryij,rzij,rij
       double precision :: r1xij,r1yij,r1zij
       double precision :: r2xij,r2yij,r2zij
       double precision :: r3xij,r3yij,r3zij
       double precision :: r4xij,r4yij,r4zij
       double precision :: r5xij,r5yij,r5zij
       double precision :: r6xij,r6yij,r6zij
       double precision :: dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6
       double precision :: r1,r2,r3,r4,r5,r6

       double precision, dimension(6) :: rin!,fn
       double precision :: fs! gama,rzero,fs
!       double precision, dimension(3) :: rH3,dvdrV2,dvdrH3
       double precision :: V2,V1,V,dV2,dV3,dvdiat,V4 ,vdiat ,v3
       double precision :: f3,f4
!       double precision :: dr1,dr2,dr3,dr4
       double precision :: dfsdr1,dfsdr2,dfsdr3,dfsdr4,dfsdr5,dfsdr6
       double precision, dimension(4,3) :: vect_4atoms ,dvdcc
       double precision, dimension(6) ::r_4atoms,dvdr_4atoms
       double precision :: VHO3_13,vo3_f7,vho2_f4
       double precision :: dv1dr1,dv1dr2,dv1dr3,dv1dr4,dv1dr5,dv1dr6



      fs=0.0d0
       V=0.0d0
      V1=0.0d0
      V2=0.0d0
      V3=0.0d0
      dv1=0.0d0
      dv2=0.0d0
      dv3=0.0d0
      dv4=0.0d0
      dv=0.0d0
      V4=0.0d0
      vect_4atoms=0.0d0
      dvdcc=0.0d0

            dfsdr1=0.0d0
            dfsdr2=0.0d0
            dfsdr3=0.0d0
            dfsdr4=0.0d0
            dfsdr5=0.0d0
            dfsdr6=0.0d0

            r1xij=0.0d0
            r1yij=0.0d0
            r1zij=0.0d0

            r2xij=0.0d0
            r2yij=0.0d0
            r2zij=0.0d0

            r3xij=0.0d0
            r3yij=0.0d0
            r3zij=0.0d0

            r4xij=0.0d0
            r4yij=0.0d0
            r4zij=0.0d0

            r5xij=0.0d0
            r5yij=0.0d0
            r5zij=0.0d0

            r6xij=0.0d0
            r6yij=0.0d0
            r6zij=0.0d0

            dv2dr1=0.0d0
            dv2dr2=0.0d0
            dv2dr3=0.0d0
            dv2dr4=0.0d0
            dv2dr5=0.0d0
            dv2dr6=0.0d0

            dv1dr1=0.0d0
            dv1dr2=0.0d0
            dv1dr3=0.0d0
            dv1dr4=0.0d0
            dv1dr5=0.0d0
            dv1dr6=0.0d0




            !Geting the coordenates

            !The positions are predefined in the matrix "group"
            gi3=group(i,3)      ! O atom
            gi4=group(i,4)      ! O atom
            gi5=group(i,5)      ! O atom
            gi6=group(i,6)      ! O atom



!      Print*,rx(gi3),ry(gi3),rz(gi3)
!      Print*,rx(gi4),ry(gi4),rz(gi4)
!      Print*,rx(gi5),ry(gi5),rz(gi5)
!      Print*,rx(gi6),ry(gi6),rz(gi6)



















        ! Checking inner distance



               call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi4),ry(gi4), &
        rz(gi4),rxij,ryij,rzij,rij)

            rin(2)=rij
            r2=rij
            r2xij=rxij
            r2yij=ryij
            r2zij=rzij

               call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi5),ry(gi5), &
        rz(gi5),rxij,ryij,rzij,rij)

            rin(3)=rij
            r3=rij
            r3xij=rxij
            r3yij=ryij
            r3zij=rzij


              call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi6),ry(gi6), &
        rz(gi6),rxij,ryij,rzij,rij)

            rin(5)=rij
            r5=rij
            r5xij=rxij
            r5yij=ryij
            r5zij=rzij


               call mod_dist(rx(gi4),ry(gi4),rz(gi4),rx(gi5),ry(gi5), &
        rz(gi5),rxij,ryij,rzij,rij)

            rin(6)=rij
            r6=rij
            r6xij=rxij
            r6yij=ryij
            r6zij=rzij



               call mod_dist(rx(gi4),ry(gi4),rz(gi4),rx(gi6),ry(gi6), &
        rz(gi6),rxij,ryij,rzij,rij)

            rin(4)=rij
            r4=rij
            r4xij=rxij
            r4yij=ryij
            r4zij=rzij

               call mod_dist(rx(gi5),ry(gi5),rz(gi5),rx(gi6),ry(gi6), &
        rz(gi6),rxij,ryij,rzij,rij)

            rin(1)=rij
            r1=rij
            r1xij=rxij
            r1yij=ryij
            r1zij=rzij

            r_4atoms(1)=r1
            r_4atoms(2)=r2
            r_4atoms(3)=r3
            r_4atoms(4)=r4
            r_4atoms(5)=r5
            r_4atoms(6)=r6


      V1=VHO3_13(r_4atoms(2),r_4atoms(3),r_4atoms(5),r_4atoms(6),r_4atoms(4),r_4atoms(1))

      call DVHO3_13(r_4atoms(2),r_4atoms(3),r_4atoms(5),r_4atoms(6),r_4atoms(4),r_4atoms(1),dvdr_4atoms)


            dv1dr2=dvdr_4atoms(1)
            dv1dr3=dvdr_4atoms(2)
            dv1dr5=dvdr_4atoms(3)
            dv1dr6=dvdr_4atoms(4)
            dv1dr4=dvdr_4atoms(5)
            dv1dr1=dvdr_4atoms(6)


!             print*,'V1',V1/con
!             print*,dvdr_4atoms


      dv1(1)  =  dv1(1)  + dv1dr2 * r2xij / r2 + dv1dr3 * r3xij / r3 + dv1dr5 * r5xij / r5
      dv1(2)  =  dv1(2)  + dv1dr2 * r2yij / r2 + dv1dr3 * r3yij / r3 + dv1dr5 * r5yij / r5
      dv1(3)  =  dv1(3)  + dv1dr2 * r2zij / r2 + dv1dr3 * r3zij / r3 + dv1dr5 * r5zij / r5

      dv1(4)  =  dv1(4)  - dv1dr2 * r2xij / r2 + dv1dr4 * r4xij / r4 + dv1dr6 * r6xij / r6
      dv1(5)  =  dv1(5)  - dv1dr2 * r2yij / r2 + dv1dr4 * r4yij / r4 + dv1dr6 * r6yij / r6
      dv1(6)  =  dv1(6)  - dv1dr2 * r2zij / r2 + dv1dr4 * r4zij / r4 + dv1dr6 * r6zij / r6

      dv1(7)  =  dv1(7)  + dv1dr1 * r1xij / r1 - dv1dr3 * r3xij / r3 - dv1dr6 * r6xij / r6
      dv1(8)  =  dv1(8)  + dv1dr1 * r1yij / r1 - dv1dr3 * r3yij / r3 - dv1dr6 * r6yij / r6
      dv1(9)  =  dv1(9)  + dv1dr1 * r1zij / r1 - dv1dr3 * r3zij / r3 - dv1dr6 * r6zij / r6

      dv1(10) =  dv1(10) - dv1dr1 * r1xij / r1 - dv1dr4 * r4xij / r4 - dv1dr5 * r5xij / r5
      dv1(11) =  dv1(11) - dv1dr1 * r1yij / r1 - dv1dr4 * r4yij / r4 - dv1dr5 * r5yij / r5
      dv1(12) =  dv1(12) - dv1dr1 * r1zij / r1 - dv1dr4 * r4zij / r4 - dv1dr5 * r5zij / r5


!          print*,'V1A',V1/con
!          stop


               !'HO3 -> O(4) + HO2'

            if ((r1.ge.in_bound).and.(r4.ge.in_bound).and. &
        (r5.ge.in_bound)) then

            if(group(i,9).eq.0) group(i,9)=11
            if(group(i,9).ne.11) then
            print*,'exchange in molecular complex 1 ',tempo
            group(i,9)=11
            end if


        ! O3 potential

       V2=vho2_f4(r6,r2,r3)
       call dervho2_f4(r6,r2,r3,dv2dr6,dv2dr2,dv2dr3)



		! Switch function


        fs=f3(rin(1),rin(4),rin(5))
        call df3(rin(1),rin(4),rin(5),dfsdr1,dfsdr4,dfsdr5)



        ! for the diatomic potencials

            call pothoq_f17(rin(5),Vdiat,dVdiat)

!            print*,'rin(5),Vdiat'
!            print*,rin(5),Vdiat



            dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
            dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
            dv(3) = dv(3)  + fs * dVdiat*r5zij/r5

            dv(10)= dv(10) - fs * dVdiat*r5xij/r5
            dv(11)= dv(11) - fs * dVdiat*r5yij/r5
            dv(12)= dv(12) - fs * dVdiat*r5zij/r5


            V4 = V4 + Vdiat


            call poto2q_f18(rin(4),Vdiat,dVdiat)

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat



            dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
            dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
            dv(6) = dv(6)  + fs * dVdiat*r4zij/r4

            dv(10)= dv(10) - fs * dVdiat*r4xij/r4
            dv(11)= dv(11) - fs * dVdiat*r4yij/r4
            dv(12)= dv(12) - fs * dVdiat*r4zij/r4


            V4 = V4 + Vdiat



            call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(1),Vdiat'
!            print*,rin(1),Vdiat



            dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
            dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
            dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

            dv(10)= dv(10) - fs * dVdiat*r1xij/r1
            dv(11)= dv(11) - fs * dVdiat*r1yij/r1
            dv(12)= dv(12) - fs * dVdiat*r1zij/r1


            V4 = V4 + Vdiat

        V=V1*(1-fs)+(V2+V4)*fs
!        print*,rin

       dv=dv + dv1*(1-fs)




       dv(1)  = dv(1) + (v4 + v2 - v1) * ( dfsdr5 * r5xij / r5)
       dv(2)  = dv(2) + (v4 + v2 - v1) * ( dfsdr5 * r5yij / r5)
       dv(3)  = dv(3) + (v4 + v2 - v1) * ( dfsdr5 * r5zij / r5)

       dv(4)  = dv(4) + (v4 + v2 - v1) * ( dfsdr4 * r4xij / r4)
       dv(5)  = dv(5) + (v4 + v2 - v1) * ( dfsdr4 * r4yij / r4)
       dv(6)  = dv(6) + (v4 + v2 - v1) * ( dfsdr4 * r4zij / r4)

       dv(7)  = dv(7) + (v4 + v2 - v1) * ( dfsdr1 * r1xij / r1)
       dv(8)  = dv(8) + (v4 + v2 - v1) * ( dfsdr1 * r1yij / r1)
       dv(9)  = dv(9) + (v4 + v2 - v1) * ( dfsdr1 * r1zij / r1)

       dv(10) = dv(10) + (v4 +v2 - v1) * (-dfsdr1 * r1xij / r1 - dfsdr4 * r4xij / r4 - dfsdr5 * r5xij / r5)
       dv(11) = dv(11) + (v4 +v2 - v1) * (-dfsdr1 * r1yij / r1 - dfsdr4 * r4yij / r4 - dfsdr5 * r5yij / r5)
       dv(12) = dv(12) + (v4 +v2 - v1) * (-dfsdr1 * r1zij / r1 - dfsdr4 * r4zij / r4 - dfsdr5 * r5zij / r5)


       dv(1)  = dv(1)  + fs * ( dv2dr2 * r2xij / r2 + dv2dr3 * r3xij /r3)
       dv(2)  = dv(2)  + fs * ( dv2dr2 * r2yij / r2 + dv2dr3 * r3yij /r3)
       dv(3)  = dv(3)  + fs * ( dv2dr2 * r2zij / r2 + dv2dr3 * r3zij /r3)

       dv(4)  = dv(4)  + fs * (-dv2dr2 * r2xij / r2 + dv2dr6 * r6xij /r6)
       dv(5)  = dv(5)  + fs * (-dv2dr2 * r2yij / r2 + dv2dr6 * r6yij /r6)
       dv(6)  = dv(6)  + fs * (-dv2dr2 * r2zij / r2 + dv2dr6 * r6zij /r6)

       dv(7)  = dv(7)  + fs * (-dv2dr3 * r3xij / r3 - dv2dr6 * r6xij /r6)
       dv(8)  = dv(8)  + fs * (-dv2dr3 * r3yij / r3 - dv2dr6 * r6yij /r6)
       dv(9)  = dv(9)  + fs * (-dv2dr3 * r3zij / r3 - dv2dr6 * r6zij /r6)


!       print*,'dv',dv


            !HO3 -> H(1) + O3
         else if ((r2.ge.in_bound).and.(r3.ge.in_bound).and. &
        (r5.ge.in_bound)) then

        if(group(i,9).eq.0)  group(i,9)=12

        if(group(i,9).ne.12) then
        print*,'Change in group'
        group(i,9)=12
        end if


         print*,"mid H(1) + O3 tempo=",tempo
!        print*,'Exchange in channels. Time:',t



        ! O3 potential

       V2=vo3_f7(r6,r4,r1)

       call dero3_f7(r6,r4,r1,dv2dr6,dv2dr4,dv2dr1)




        ! Switch function


        fs=f3(rin(2),rin(3),rin(5))
        call df3(rin(2),rin(3),rin(5),dfsdr2,dfsdr3,dfsdr5)






        ! for the diatomic potencials

            call pothoq_f17(rin(5),Vdiat,dVdiat)

!            print*,'rin(5),Vdiat'
!            print*,rin(5),Vdiat



            dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
            dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
            dv(3) = dv(3)  + fs * dVdiat*r5zij/r5

            dv(10)= dv(10) - fs * dVdiat*r5xij/r5
            dv(11)= dv(11) - fs * dVdiat*r5yij/r5
            dv(12)= dv(12) - fs * dVdiat*r5zij/r5


            V4 = V4 + Vdiat


            call pothoq_f17(rin(2),Vdiat,dVdiat)

!            print*,'rin(2),Vdiat'
!            print*,rin(2),Vdiat

            dv(1) = dv(1) + fs * dVdiat*r2xij/r2
            dv(2) = dv(2) + fs * dVdiat*r2yij/r2
            dv(3) = dv(3) + fs * dVdiat*r2zij/r2

            dv(4) = dv(4) - fs * dVdiat*r2xij/r2
            dv(5) = dv(5) - fs * dVdiat*r2yij/r2
            dv(6) = dv(6) - fs * dVdiat*r2zij/r2


            V4 = V4 + Vdiat


            call pothoq_f17(rin(3),Vdiat,dVdiat)

!            print*,'rin(3),Vdiat'
!            print*,rin(3),Vdiat



            dv(1) = dv(1) + fs * dVdiat*r3xij/r3
            dv(2) = dv(2) + fs * dVdiat*r3yij/r3
            dv(3) = dv(3) + fs * dVdiat*r3zij/r3

            dv(7) = dv(7) - fs * dVdiat*r3xij/r3
            dv(8) = dv(8) - fs * dVdiat*r3yij/r3
            dv(9) = dv(9) - fs * dVdiat*r3zij/r3


            V4 = V4 + Vdiat

        V=V1*(1-fs)+(V2+V4)*fs
!        print*,rin

      dv=dv + dv1*(1-fs)


      dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr2 *  r2xij / r2 + dfsdr3 * r3xij / r3 + dfsdr5 * r5xij / r5)
      dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr2 *  r2yij / r2 + dfsdr3 * r3yij / r3 + dfsdr5 * r5yij / r5)
      dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr2 *  r2zij / r2 + dfsdr3 * r3zij / r3 + dfsdr5 * r5zij / r5)

      dv(4) = dv(4) + (v4 + v2 - v1) * (-dfsdr2 *  r2xij / r2)
      dv(5) = dv(5) + (v4 + v2 - v1) * (-dfsdr2 *  r2yij / r2)
      dv(6) = dv(6) + (v4 + v2 - v1) * (-dfsdr2 *  r2zij / r2)

      dv(7) = dv(7) + (v4 + v2 - v1) * (-dfsdr3 *  r3xij / r3)
      dv(8) = dv(8) + (v4 + v2 - v1) * (-dfsdr3 *  r3yij / r3)
      dv(9) = dv(9) + (v4 + v2 - v1) * (-dfsdr3 *  r3zij / r3)

      dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr5 *  r5xij / r5)
      dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr5 *  r5yij / r5)
      dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr5 *  r5zij / r5)


      dv(4) = dv(4) + fs * ( dv2dr4 * r4xij / r4 + dv2dr6 * r6xij / r6)
      dv(5) = dv(5) + fs * ( dv2dr4 * r4yij / r4 + dv2dr6 * r6yij / r6)
      dv(6) = dv(6) + fs * ( dv2dr4 * r4zij / r4 + dv2dr6 * r6zij / r6)

      dv(7) = dv(7)  + fs * ( dv2dr1 * r1xij / r1 - dv2dr6 * r6xij / r6)
      dv(8) = dv(8)  + fs * ( dv2dr1 * r1yij / r1 - dv2dr6 * r6yij / r6)
      dv(9) = dv(9)  + fs * ( dv2dr1 * r1zij / r1 - dv2dr6 * r6zij / r6)

      dv(10)= dv(10)  + fs * (-dv2dr1 * r1xij / r1 - dv2dr4 * r4xij / r4)
      dv(11)= dv(11)  + fs * (-dv2dr1 * r1yij / r1 - dv2dr4 * r4yij / r4)
      dv(12)= dv(12)  + fs * (-dv2dr1 * r1zij / r1 - dv2dr4 * r4zij / r4)




            !'HO3 -> O(2) + HO2
        else if ((r2.ge.in_bound).and.(r4.ge.in_bound).and. &
        (r6.ge.in_bound)) then

        if(group(i,9).eq.0)  group(i,9)=13

        if(group(i,9).ne.13) then
        print*,'Change in group'
        group(i,9)=13
        end if


!        print*,'here we are'
            !HO3 -> O(2) + HO2

!             print*,'mid O(2) + HO2',tempo
        !print*,'Exchange in channels. Time:',t

        ! HO2 potential




       V2=vho2_f4(r4,r2,r5)
       call dervho2_f4(r4,r2,r5,dv2dr4,dv2dr2,dv2dr5)










        ! Switch function


        fs=f3(rin(1),rin(3),rin(6))
        call df3(rin(1),rin(3),rin(6),dfsdr1,dfsdr3,dfsdr6)


!        print*,'fs=',fs,dfsdr2,dfsdr4,dfsdr6





        ! for the diatomic potencials

            call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat



            dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
            dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
            dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

            dv(10)= dv(10) - fs * dVdiat*r1xij/r1
            dv(11)= dv(11) - fs * dVdiat*r1yij/r1
            dv(12)= dv(12) - fs * dVdiat*r1zij/r1


            V4 = V4 + Vdiat


            call pothoq_f17(rin(3),Vdiat,dVdiat)


!            print*,'rin(2),Vdiat'
!            print*,rin(2),Vdiat

            dv(1) = dv(1) + fs * dVdiat*r3xij/r3
            dv(2) = dv(2) + fs * dVdiat*r3yij/r3
            dv(3) = dv(3) + fs * dVdiat*r3zij/r3

            dv(7) = dv(7) - fs * dVdiat*r3xij/r3
            dv(8) = dv(8) - fs * dVdiat*r3yij/r3
            dv(9) = dv(9) - fs * dVdiat*r3zij/r3


            V4 = V4 + Vdiat


            call poto2q_f18(rin(6),Vdiat,dVdiat)

!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat



            dv(4) = dv(4) + fs * dVdiat*r6xij/r6
            dv(5) = dv(5) + fs * dVdiat*r6yij/r6
            dv(6) = dv(6) + fs * dVdiat*r6zij/r6

            dv(7) = dv(7) - fs * dVdiat*r6xij/r6
            dv(8) = dv(8) - fs * dVdiat*r6yij/r6
            dv(9) = dv(9) - fs * dVdiat*r6zij/r6


            V4 = V4 + Vdiat

        V=V1*(1-fs)+(V2+V4)*fs
!        print*,rin

      dv=dv + dv1*(1-fs)


      dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr2 *  r2xij / r2)
      dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr2 *  r2yij / r2)
      dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr2 *  r2zij / r2)

      dv(4) = dv(4) + (v4 + v2 - v1) * (-dfsdr2 *  r2xij / r2 + dfsdr4 *  r4xij / r4 + dfsdr6 *  r6xij / r6)
      dv(5) = dv(5) + (v4 + v2 - v1) * (-dfsdr2 *  r2yij / r2 + dfsdr4 *  r4yij / r4 + dfsdr6 *  r6yij / r6)
      dv(6) = dv(6) + (v4 + v2 - v1) * (-dfsdr2 *  r2zij / r2 + dfsdr4 *  r4zij / r4 + dfsdr6 *  r6zij / r6)

      dv(7) = dv(7) + (v4 + v2 - v1) * (-dfsdr6 *  r6xij / r6)
      dv(8) = dv(8) + (v4 + v2 - v1) * (-dfsdr6 *  r6yij / r6)
      dv(9) = dv(9) + (v4 + v2 - v1) * (-dfsdr6 *  r6zij / r6)

      dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr4 *  r4xij / r4)
      dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr4 *  r4yij / r4)
      dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr4 *  r4zij / r4)


      dv(1) = dv(1) + fs * ( dv2dr2 * r2xij / r2 + dv2dr5 * r5xij / r5)
      dv(2) = dv(2) + fs * ( dv2dr2 * r2yij / r2 + dv2dr5 * r5yij / r5)
      dv(3) = dv(3) + fs * ( dv2dr2 * r2zij / r2 + dv2dr5 * r5zij / r5)

      dv(4)  = dv(4)  + fs * (-dv2dr2 * r2xij / r2 + dv2dr4 * r4xij / r4)
      dv(5)  = dv(5)  + fs * (-dv2dr2 * r2yij / r2 + dv2dr4 * r4yij / r4)
      dv(6)  = dv(6)  + fs * (-dv2dr2 * r2zij / r2 + dv2dr4 * r4zij / r4)

      dv(10)  = dv(10)  + fs * (-dv2dr4 * r4xij / r4 - dv2dr5 * r5xij / r5)
      dv(11)  = dv(11)  + fs * (-dv2dr4 * r4yij / r4 - dv2dr5 * r5yij / r5)
      dv(12)  = dv(12)  + fs * (-dv2dr4 * r4zij / r4 - dv2dr5 * r5zij / r5)

            !'HO3 -> O(3) + HO2
        else if ((r1.ge.in_bound).and.(r3.ge.in_bound).and. &
        (r6.ge.in_bound)) then

        if(group(i,9).eq.0)  group(i,9)=14

        if(group(i,9).ne.14) then
        print*,'Change in group'
        group(i,9)=14
        end if
           !HO3 -> O(3) + HO2




!             print*,"mid O(3) + HO2",tempo
        !print*,'Exchange in channels. Time:',t

        ! H3 potential

       V2=vho2_f4(r4,r2,r5)
       call dervho2_f4(r4,r2,r5,dv2dr4,dv2dr2,dv2dr5)


        ! Switch function


        fs=f3(rin(1),rin(3),rin(6))
        call df3(rin(1),rin(3),rin(6),dfsdr1,dfsdr3,dfsdr6)


!        print*,'fs=',fs,dfsdr1,dfsdr3,dfsdr6





        ! for the diatomic potencials

            call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(1),Vdiat'
!            print*,rin(1),Vdiat



            dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
            dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
            dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

            dv(10)= dv(10) - fs * dVdiat*r1xij/r1
            dv(11)= dv(11) - fs * dVdiat*r1yij/r1
            dv(12)= dv(12) - fs * dVdiat*r1zij/r1


            V4 = V4 + Vdiat


            call pothoq_f17(rin(3),Vdiat,dVdiat)


!            print*,'rin(3),Vdiat'
!            print*,rin(3),Vdiat

            dv(1) = dv(1) + fs * dVdiat*r3xij/r3
            dv(2) = dv(2) + fs * dVdiat*r3yij/r3
            dv(3) = dv(3) + fs * dVdiat*r3zij/r3

            dv(7) = dv(7) - fs * dVdiat*r3xij/r3
            dv(8) = dv(8) - fs * dVdiat*r3yij/r3
            dv(9) = dv(9) - fs * dVdiat*r3zij/r3


            V4 = V4 + Vdiat


            call poto2q_f18(rin(6),Vdiat,dVdiat)

!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat



            dv(4) = dv(4) + fs * dVdiat*r6xij/r6
            dv(5) = dv(5) + fs * dVdiat*r6yij/r6
            dv(6) = dv(6) + fs * dVdiat*r6zij/r6

            dv(7) = dv(7) - fs * dVdiat*r6xij/r6
            dv(8) = dv(8) - fs * dVdiat*r6yij/r6
            dv(9) = dv(9) - fs * dVdiat*r6zij/r6


            V4 = V4 + Vdiat

        V=V1*(1-fs)+(V2+V4)*fs
!        print*,rin
!        print*,'V=',V,' V1=',V1,' V2=',V2

      dv=dv + dv1*(1-fs)


      dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr3 *  r3xij / r3)
      dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr3 *  r3yij / r3)
      dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr3 *  r3zij / r3)

      dv(4) = dv(4) + (v4 + v2 - v1) * ( dfsdr6 *  r6xij / r6)
      dv(5) = dv(5) + (v4 + v2 - v1) * ( dfsdr6 *  r6yij / r6)
      dv(6) = dv(6) + (v4 + v2 - v1) * ( dfsdr6 *  r6zij / r6)

      dv(7) = dv(7) + (v4 + v2 - v1) * ( dfsdr1 *  r1xij / r1 - dfsdr3 *  r3xij / r3 - dfsdr6 *  r6xij / r6)
      dv(8) = dv(8) + (v4 + v2 - v1) * ( dfsdr1 *  r1yij / r1 - dfsdr3 *  r3yij / r3 - dfsdr6 *  r6yij / r6)
      dv(9) = dv(9) + (v4 + v2 - v1) * ( dfsdr1 *  r1zij / r1 - dfsdr3 *  r3zij / r3 - dfsdr6 *  r6zij / r6)

      dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr1 *  r1xij / r1)
      dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr1 *  r1yij / r1)
      dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr1 *  r1zij / r1)


      dv(1) = dv(1) + fs * ( dv2dr2 * r2xij / r2 + dv2dr5 * r5xij / r5)
      dv(2) = dv(2) + fs * ( dv2dr2 * r2yij / r2 + dv2dr5 * r5yij / r5)
      dv(3) = dv(3) + fs * ( dv2dr2 * r2zij / r2 + dv2dr5 * r5zij / r5)

      dv(4)  = dv(4)  + fs * (-dv2dr2 * r2xij / r2 + dv2dr4 * r4xij / r4)
      dv(5)  = dv(5)  + fs * (-dv2dr2 * r2yij / r2 + dv2dr4 * r4yij / r4)
      dv(6)  = dv(6)  + fs * (-dv2dr2 * r2zij / r2 + dv2dr4 * r4zij / r4)

      dv(10)  = dv(10)  + fs * (-dv2dr4 * r4xij / r4 - dv2dr5 * r5xij / r5)
      dv(11)  = dv(11)  + fs * (-dv2dr4 * r4yij / r4 - dv2dr5 * r5yij / r5)
      dv(12)  = dv(12)  + fs * (-dv2dr4 * r4zij / r4 - dv2dr5 * r5zij / r5)


         else if ((r1.ge.in_bound).and.(r2.ge.in_bound).and. &
        (r5.ge.in_bound).and.(r6.ge.in_bound)) then
!              print*,'exit'



             ! 'HO3 -> HO(1-3) + O2(2-4)'
        if(group(i,9).eq.0)  group(i,9)=15

        if(group(i,9).ne.15) then
        print*,'Change in group'
!        stop 'here'
        group(i,9)=15
        end if





         	! HO function
        	call potohp_f3(r3,V2,dV2)
            ! O2 function
            call potoop_f2(r4,V3,dV3)



       ! Switch function
        fs=f4(rin(1),rin(2),rin(5),rin(6))
        call df4(rin(1),rin(2),rin(5),rin(6),dfsdr1,dfsdr2,dfsdr5,dfsdr6)



        ! for the diatomic potencials

            call pothoq_f17(rin(5),Vdiat,dVdiat)

!            print*,'rin(5),Vdiat'
!            print*,rin(5),Vdiat



            dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
            dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
            dv(3) = dv(3)  + fs * dVdiat*r5zij/r5

            dv(10)= dv(10) - fs * dVdiat*r5xij/r5
            dv(11)= dv(11) - fs * dVdiat*r5yij/r5
            dv(12)= dv(12) - fs * dVdiat*r5zij/r5


            V4 = V4 + Vdiat



            call pothoq_f17(rin(2),Vdiat,dVdiat)


!            print*,'rin(2),Vdiat'
!            print*,rin(2),Vdiat

            dv(1) = dv(1) + fs * dVdiat*r2xij/r2
            dv(2) = dv(2) + fs * dVdiat*r2yij/r2
            dv(3) = dv(3) + fs * dVdiat*r2zij/r2

            dv(4) = dv(4) - fs * dVdiat*r2xij/r2
            dv(5) = dv(5) - fs * dVdiat*r2yij/r2
            dv(6) = dv(6) - fs * dVdiat*r2zij/r2


            V4 = V4 + Vdiat



            call poto2q_f18(rin(4),Vdiat,dVdiat)

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat



            dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
            dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
            dv(6) = dv(6)  + fs * dVdiat*r4zij/r4

            dv(10)= dv(10) - fs * dVdiat*r4xij/r4
            dv(11)= dv(11) - fs * dVdiat*r4yij/r4
            dv(12)= dv(12) - fs * dVdiat*r4zij/r4


            V4 = V4 + Vdiat



            call poto2q_f18(rin(6),Vdiat,dVdiat)


!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat

            dv(4) = dv(4) + fs * dVdiat*r6xij/r6
            dv(5) = dv(5) + fs * dVdiat*r6yij/r6
            dv(6) = dv(6) + fs * dVdiat*r6zij/r6

            dv(7) = dv(7) - fs * dVdiat*r6xij/r6
            dv(8) = dv(8) - fs * dVdiat*r6yij/r6
            dv(9) = dv(9) - fs * dVdiat*r6zij/r6


            V4 = V4 + Vdiat



            call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(1),Vdiat'
!            print*,rin(1),Vdiat



            dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
            dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
            dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

            dv(10)= dv(10) - fs * dVdiat*r1xij/r1
            dv(11)= dv(11) - fs * dVdiat*r1yij/r1
            dv(12)= dv(12) - fs * dVdiat*r1zij/r1


            V4 = V4 + Vdiat





        V=V1*(1-fs)+(V2+V3+V4)*fs

       dv=dv + dv1*(1-fs)

       dv(1)  =   dv(1)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2xij / r2 + dfsdr5 * r5xij / r5)
       dv(2)  =   dv(2)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2yij / r2 + dfsdr5 * r5yij / r5)
       dv(3)  =   dv(3)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2zij / r2 + dfsdr5 * r5zij / r5)

       dv(4)  =   dv(4)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2xij / r2 + dfsdr6 * r6xij / r6)
       dv(5)  =   dv(5)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2yij / r2 + dfsdr6 * r6yij / r6)
       dv(6)  =   dv(6)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2zij / r2 + dfsdr6 * r6zij / r6)

       dv(7)  =   dv(7)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1xij / r1 - dfsdr6 * r6xij / r6)
       dv(8)  =   dv(8)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1yij / r1 - dfsdr6 * r6yij / r6)
       dv(9)  =   dv(9)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1zij / r1 - dfsdr6 * r6zij / r6)

       dv(10) =   dv(10) + (v2+v3+v4-v1) * (-dfsdr1 *  r1xij / r1 - dfsdr5 * r5xij / r5)
       dv(11) =   dv(11) + (v2+v3+v4-v1) * (-dfsdr1 *  r1yij / r1 - dfsdr5 * r5yij / r5)
       dv(12) =   dv(12) + (v2+v3+v4-v1) * (-dfsdr1 *  r1zij / r1 - dfsdr5 * r5zij / r5)


       dv(1)  =   dv(1)  + fs * ( dv2 * r3xij / r3)
       dv(2)  =   dv(2)  + fs * ( dv2 * r3yij / r3)
       dv(3)  =   dv(3)  + fs * ( dv2 * r3zij / r3)

       dv(4)  =   dv(4)  + fs * ( dv3 * r4xij / r4)
       dv(5)  =   dv(5)  + fs * ( dv3 * r4yij / r4)
       dv(6)  =   dv(6)  + fs * ( dv3 * r4zij / r4)

       dv(7)  =   dv(7)  + fs * (-dv2 * r3xij / r3)
       dv(8)  =   dv(8)  + fs * (-dv2 * r3yij / r3)
       dv(9)  =   dv(9)  + fs * (-dv2 * r3zij / r3)

       dv(10) =   dv(10) + fs * (-dv3 * r4xij / r4)
       dv(11) =   dv(11) + fs * (-dv3 * r4yij / r4)
       dv(12) =   dv(12) + fs * (-dv3 * r4zij / r4)









         else if ((r1.ge.in_bound).and.(r2.ge.in_bound).and. &
        (r3.ge.in_bound).and.(r4.ge.in_bound)) then
!              print*,'exit5'



                    ! 'HO3 -> O2(2-3) + HO(1-4)'
        if(group(i,9).eq.0)  group(i,9)=16

        if(group(i,9).ne.16) then
        print*,'Change in group'
        group(i,9)=16
        end if



        	! O2 function
        	call potoop_f2(r6,V2,dV2)
            ! O2 function
            call potohp_f3(r5,V3,dV3)



       ! Switch function
        fs=f4(rin(1),rin(2),rin(3),rin(4))
        call df4(rin(1),rin(2),rin(3),rin(4),dfsdr1,dfsdr2,dfsdr3,dfsdr4)



        ! for the diatomic potencials

            call pothoq_f17(rin(2),Vdiat,dVdiat)


!            print*,'rin(2),Vdiat'
!            print*,rin(2),Vdiat

            dv(1) = dv(1) + fs * dVdiat*r2xij/r2
            dv(2) = dv(2) + fs * dVdiat*r2yij/r2
            dv(3) = dv(3) + fs * dVdiat*r2zij/r2

            dv(4) = dv(4) - fs * dVdiat*r2xij/r2
            dv(5) = dv(5) - fs * dVdiat*r2yij/r2
            dv(6) = dv(6) - fs * dVdiat*r2zij/r2


            V4 = V4 + Vdiat


            call pothoq_f17(rin(3),Vdiat,dVdiat)


!            print*,'rin(3),Vdiat'
!            print*,rin(3),Vdiat

            dv(1) = dv(1) + fs * dVdiat*r3xij/r3
            dv(2) = dv(2) + fs * dVdiat*r3yij/r3
            dv(3) = dv(3) + fs * dVdiat*r3zij/r3

            dv(7) = dv(7) - fs * dVdiat*r3xij/r3
            dv(8) = dv(8) - fs * dVdiat*r3yij/r3
            dv(9) = dv(9) - fs * dVdiat*r3zij/r3


            V4 = V4 + Vdiat



            call poto2q_f18(rin(4),Vdiat,dVdiat)

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat



            dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
            dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
            dv(6) = dv(6)  + fs * dVdiat*r4zij/r4

            dv(10)= dv(10) - fs * dVdiat*r4xij/r4
            dv(11)= dv(11) - fs * dVdiat*r4yij/r4
            dv(12)= dv(12) - fs * dVdiat*r4zij/r4


            V4 = V4 + Vdiat



            call poto2q_f18(rin(6),Vdiat,dVdiat)


!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat

            dv(4) = dv(4) + fs * dVdiat*r6xij/r6
            dv(5) = dv(5) + fs * dVdiat*r6yij/r6
            dv(6) = dv(6) + fs * dVdiat*r6zij/r6

            dv(7) = dv(7) - fs * dVdiat*r6xij/r6
            dv(8) = dv(8) - fs * dVdiat*r6yij/r6
            dv(9) = dv(9) - fs * dVdiat*r6zij/r6


            V4 = V4 + Vdiat






            call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(1),Vdiat'
!            print*,rin(1),Vdiat



            dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
            dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
            dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

            dv(10)= dv(10) - fs * dVdiat*r1xij/r1
            dv(11)= dv(11) - fs * dVdiat*r1yij/r1
            dv(12)= dv(12) - fs * dVdiat*r1zij/r1


            V4 = V4 + Vdiat





        V=V1*(1-fs)+(V2+V3+V4)*fs

       dv=dv + dv1*(1-fs)

       dv(1)  =   dv(1)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2xij / r2 + dfsdr3 * r3xij / r3)
       dv(2)  =   dv(2)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2yij / r2 + dfsdr3 * r3yij / r3)
       dv(3)  =   dv(3)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2zij / r2 + dfsdr3 * r3zij / r3)

       dv(4)  =   dv(4)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2xij / r2 + dfsdr4 * r4xij / r4)
       dv(5)  =   dv(5)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2yij / r2 + dfsdr4 * r4yij / r4)
       dv(6)  =   dv(6)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2zij / r2 + dfsdr4 * r4zij / r4)

       dv(7)  =   dv(7)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1xij / r1 - dfsdr3 * r3xij / r3)
       dv(8)  =   dv(8)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1yij / r1 - dfsdr3 * r3yij / r3)
       dv(9)  =   dv(9)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1zij / r1 - dfsdr3 * r3zij / r3)

       dv(10) =   dv(10) + (v2+v3+v4-v1) * (-dfsdr1 *  r1xij / r1 - dfsdr4 * r4xij / r4)
       dv(11) =   dv(11) + (v2+v3+v4-v1) * (-dfsdr1 *  r1yij / r1 - dfsdr4 * r4yij / r4)
       dv(12) =   dv(12) + (v2+v3+v4-v1) * (-dfsdr1 *  r1zij / r1 - dfsdr4 * r4zij / r4)


       dv(1)  =   dv(1)  + fs * ( dv3 * r5xij / r5)
       dv(2)  =   dv(2)  + fs * ( dv3 * r5yij / r5)
       dv(3)  =   dv(3)  + fs * ( dv3 * r5zij / r5)

       dv(4)  =   dv(4)  + fs * ( dv2 * r6xij / r6)
       dv(5)  =   dv(5)  + fs * ( dv2 * r6yij / r6)
       dv(6)  =   dv(6)  + fs * ( dv2 * r6zij / r6)

       dv(7)  =   dv(7)  + fs * (-dv2 * r6xij / r6)
       dv(8)  =   dv(8)  + fs * (-dv2 * r6yij / r6)
       dv(9)  =   dv(9)  + fs * (-dv2 * r6zij / r6)

       dv(10) =   dv(10) + fs * (-dv3 * r5xij / r5)
       dv(11) =   dv(11) + fs * (-dv3 * r5yij / r5)
       dv(12) =   dv(12) + fs * (-dv3 * r5zij / r5)


         else if ((r3.ge.in_bound).and.(r4.ge.in_bound).and. &
        (r5.ge.in_bound).and.(r6.ge.in_bound)) then
!             print*,'exit4'



!             print*,'HO3 -> HO(1-2) + O2(3-4)'
!


        if(group(i,9).eq.0)  group(i,9)=17

        if(group(i,9).ne.17) then
        print*,'Change in group'
        group(i,9)=17
        end if
         	! HO function
        	call potohp_f3(r2,V2,dV2)
            ! O2 function
            call potoop_f2(r1,V3,dV3)



       ! Switch function
        fs=f4(rin(3),rin(4),rin(5),rin(6))
        call df4(rin(3),rin(4),rin(5),rin(6),dfsdr3,dfsdr4,dfsdr5,dfsdr6)


        ! for the diatomic potencials

            call pothoq_f17(rin(3),Vdiat,dVdiat)

            V4 = V4 + Vdiat


!            print*,'rin(3),Vdiat'
!            print*,rin(3),Vdiat/con

            dv(1) = dv(1) + fs * dVdiat*r3xij/r3
            dv(2) = dv(2) + fs * dVdiat*r3yij/r3
            dv(3) = dv(3) + fs * dVdiat*r3zij/r3

            dv(7) = dv(7) - fs * dVdiat*r3xij/r3
            dv(8) = dv(8) - fs * dVdiat*r3yij/r3
            dv(9) = dv(9) - fs * dVdiat*r3zij/r3


            call pothoq_f17(rin(5),Vdiat,dVdiat)

            V4 = V4 + Vdiat

!            print*,'rin(5),Vdiat'
!            print*,rin(5),Vdiat/con



            dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
            dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
            dv(3) = dv(3)  + fs * dVdiat*r5zij/r5

            dv(10)= dv(10) - fs * dVdiat*r5xij/r5
            dv(11)= dv(11) - fs * dVdiat*r5yij/r5
            dv(12)= dv(12) - fs * dVdiat*r5zij/r5







            call poto2q_f18(rin(6),Vdiat,dVdiat)

            V4 = V4 + Vdiat

!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat/con

            dv(4) = dv(4) + fs * dVdiat*r6xij/r6
            dv(5) = dv(5) + fs * dVdiat*r6yij/r6
            dv(6) = dv(6) + fs * dVdiat*r6zij/r6

            dv(7) = dv(7) - fs * dVdiat*r6xij/r6
            dv(8) = dv(8) - fs * dVdiat*r6yij/r6
            dv(9) = dv(9) - fs * dVdiat*r6zij/r6




            call poto2q_f18(rin(4),Vdiat,dVdiat)

            V4 = V4 + Vdiat

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat/con



            dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
            dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
            dv(6) = dv(6)  + fs * dVdiat*r4zij/r4

            dv(10)= dv(10) - fs * dVdiat*r4xij/r4
            dv(11)= dv(11) - fs * dVdiat*r4yij/r4
            dv(12)= dv(12) - fs * dVdiat*r4zij/r4



!            print*,'V2',V2/con
!            print*,'V3',V3/con
!            print*,'V4',V4/con



        V=V1*(1-fs)+(V2+V3+V4)*fs
!       print*,'V1,(V2+V3+V4),fs'
!       print*,V1/con,(V2+V3+V4)/con,fs

        dv=dv + dv1*(1-fs)
!       print*,'dv depois 0',dv
       dv(1)  =   dv(1)  + (v2+v3+v4-v1) * ( dfsdr3 *  r3xij / r3 + dfsdr5 * r5xij / r5)
!        print*
!        print*, v2,v3,v4,v1,  dfsdr3 ,  r3xij , r3 , dfsdr5 , r5xij , r5

       dv(2)  =   dv(2)  + (v2+v3+v4-v1) * ( dfsdr3 *  r3yij / r3 + dfsdr5 * r5yij / r5)
       dv(3)  =   dv(3)  + (v2+v3+v4-v1) * ( dfsdr3 *  r3zij / r3 + dfsdr5 * r5zij / r5)

       dv(4)  =   dv(4)  + (v2+v3+v4-v1) * ( dfsdr4 * r4xij / r4 + dfsdr6 * r6xij / r6)
       dv(5)  =   dv(5)  + (v2+v3+v4-v1) * ( dfsdr4 * r4yij / r4 + dfsdr6 * r6yij / r6)
       dv(6)  =   dv(6)  + (v2+v3+v4-v1) * ( dfsdr4 * r4zij / r4 + dfsdr6 * r6zij / r6)

       dv(7)  =   dv(7)  + (v2+v3+v4-v1) * (-dfsdr3 * r3xij / r3 - dfsdr6 * r6xij / r6)
       dv(8)  =   dv(8)  + (v2+v3+v4-v1) * (-dfsdr3 * r3yij / r3 - dfsdr6 * r6yij / r6)
       dv(9)  =   dv(9)  + (v2+v3+v4-v1) * (-dfsdr3 * r3zij / r3 - dfsdr6 * r6zij / r6)

       dv(10) =   dv(10) + (v2+v3+v4-v1) * (-dfsdr4 *  r4xij / r4 - dfsdr5 * r5xij / r5)
       dv(11) =   dv(11) + (v2+v3+v4-v1) * (-dfsdr4 *  r4yij / r4 - dfsdr5 * r5yij / r5)
       dv(12) =   dv(12) + (v2+v3+v4-v1) * (-dfsdr4 *  r4zij / r4 - dfsdr5 * r5zij / r5)
!            print*,'dv depois  1',dv



       dv(1)  =   dv(1)  + fs * ( dv2 * r2xij / r2)
       dv(2)  =   dv(2)  + fs * ( dv2 * r2yij / r2)
       dv(3)  =   dv(3)  + fs * ( dv2 * r2zij / r2)

       dv(4)  =   dv(4)  + fs * (-dv2 * r2xij / r2)
       dv(5)  =   dv(5)  + fs * (-dv2 * r2yij / r2)
       dv(6)  =   dv(6)  + fs * (-dv2 * r2zij / r2)

       dv(7)  =   dv(7)  + fs * ( dv3 * r1xij / r1)
       dv(8)  =   dv(8)  + fs * ( dv3 * r1yij / r1)
       dv(9)  =   dv(9)  + fs * ( dv3 * r1zij / r1)

       dv(10) =   dv(10) + fs * (-dv3 * r1xij / r1)
       dv(11) =   dv(11) + fs * (-dv3 * r1yij / r1)
       dv(12) =   dv(12) + fs * (-dv3 * r1zij / r1)

!            print*,'dv depois  2',dv



         else!  end if

         v=v1



          dv  = dv1

          group(i,9)=0

      end if

                  end subroutine




!

!
! *********************** All units are in a.u. ******************************
! *      FUNCTION VHO3_13(RAB,RAC,RAD,RBC,RBD,RCD)                               *
! *      subroutine DVHO3_13(RAB,RAC,RAD,RBC,RBD,RCD,DVR(1---6))                 *
! *The label A stands for H atom and the others B,C,D refer to three O atoms.*
! ****************************************************************************
! =                 March 1996, Coimbra                                      =
! ============================================================================
!  TOTAL PES OF HO3(a.u.): function VHO3_13(Ri)
       FUNCTION VHO3_13(RABin,RACin,RADin,RBCin,RBDin,RCDin)
       IMPLICIT REAL*8(A-H,O-Z)
       RAB=RABin/5.2917721092d-11
       RAC=RACin/5.2917721092d-11
       RAD=RADin/5.2917721092d-11
       RBC=RBCin/5.2917721092d-11
       RBD=RBDin/5.2917721092d-11
       RCD=RCDin/5.2917721092d-11


           PRECEPS=1.0D-20

       IF (RAB.LT.PRECEPS) RAB=PRECEPS
       IF (RAC.LT.PRECEPS) RAC=PRECEPS
       IF (RAD.LT.PRECEPS) RAD=PRECEPS
       IF (RBC.LT.PRECEPS) RBC=PRECEPS
       IF (RBD.LT.PRECEPS) RBD=PRECEPS
       IF (RCD.LT.PRECEPS) RCD=PRECEPS



       TWO=VOH_13(RAB)+VOH_13(RAC)+VOH_13(RAD)+VOO_13(RBC)+ &
       VOO_13(RBD)+VOO_13(RCD)
       THREE=VO3_13(RBC,RBD,RCD)+VHO2_13(RBC,RAB,RAC)+ &
       VHO2_13(RBD,RAB,RAD)+VHO2_13(RCD,RAC,RAD)
       FOURE=ELECHO3_13(RAB,RAC,RAD,RBC,RBD,RCD)
       FOURHF=FVHO3_13(RAB,RAC,RAD,RBC,RBD,RCD)
       VHO3_13=TWO+THREE+FOURE+FOURHF
       VHO3_13=VHO3_13*4.3597482D-18
!W
!        PRINT*,'VHO3_13=',VHO3
!        PRINT*,'TWO=',TWO
!        PRINT*,'THREE=',THREE
!        PRINT*,'FOURE=',FOURE
!        PRINT*,'FOURHF=',FOURHF
!W
       RETURN
       END

!        function derho2_13(r1,r2,r3,i)
!        implicit real*8(a-h,o-z)
!
!        r4=20.d0
!        if(i.eq.1)then
!        derho2_13=dvho4(r1,r4,r4,r2,r4,r4,r3,r4,r4,r4,1)
!        else if(i.eq.2)then
!        derho2_13=dvho4(r1,r4,r4,r2,r4,r4,r3,r4,r4,r4,4)
!        else if(i.eq.3)then
!        derho2_13=dvho4(r1,r4,r4,r2,r4,r4,r3,r4,r4,r4,7)
!        else
!        endif
!        end

!        function ho2_13(r1,r2,r3)
!        implicit real*8(a-h,o-z)
!        r4=20.d0
!        ho2_13=vho4(r1,r4,r4,r2,r4,r4,r3,r4,r4,r4)
!        end





      SUBROUTINE DVHO3_13(RABin,RACin,RADin,RBCin,RBDin,RCDin,DVR)
!  THE FIRST PARTIAL DERIVATIVES OF TOTAL PES OF HO3(a.u.)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DVR(6),WORK1(6),WORK2(6)
      con=4.3597482D-18/5.2917721092d-11
      RAB=RABin/5.2917721092d-11
      RAC=RACin/5.2917721092d-11
      RAD=RADin/5.2917721092d-11
      RBC=RBCin/5.2917721092d-11
      RBD=RBDin/5.2917721092d-11
      RCD=RCDin/5.2917721092d-11


      DVR(1)=DVOH_13(RAB)
      DVR(2)=DVOH_13(RAC)
      DVR(3)=DVOH_13(RAD)
      DVR(4)=DVOO_13(RBC)
      DVR(5)=DVOO_13(RBD)
      DVR(6)=DVOO_13(RCD)
      DVR(4)=DVR(4)+DVO3_13(RBC,RBD,RCD,1)
      DVR(5)=DVR(5)+DVO3_13(RBC,RBD,RCD,2)
      DVR(6)=DVR(6)+DVO3_13(RBC,RBD,RCD,3)
      CALL DERIVHO2_13(RBC,RAB,RAC,DVRBC,DVRAB,DVRAC)
      DVR(4)=DVR(4)+DVRBC
      DVR(1)=DVR(1)+DVRAB
      DVR(2)=DVR(2)+DVRAC
      CALL DERIVHO2_13(RBD,RAB,RAD,DVRBD,DVRAB,DVRAD)
      DVR(5)=DVR(5)+DVRBD
      DVR(1)=DVR(1)+DVRAB
      DVR(3)=DVR(3)+DVRAD
        CALL DERIVHO2_13(RCD,RAC,RAD,DVRCD,DVRAC,DVRAD)
      DVR(6)=DVR(6)+DVRCD
      DVR(2)=DVR(2)+DVRAC
      DVR(3)=DVR(3)+DVRAD
      CALL DELECHO3_13(RAB,RAC,RAD,RBC,RBD,RCD,WORK1)
           CALL DFVHO3_13(RAB,RAC,RAD,RBC,RBD,RCD,WORK2)
      DO 10 I=1,6
      DVR(I)=DVR(I)+WORK1(I)+WORK2(I)

 10      CONTINUE
      DVR=DVR*con
      RETURN
      END

      FUNCTION VHO2_13(ROO,RHO1,RHO2)
!=====================================================================
!  This is the latest DMBE potential energy surface for HO2
!                           March 29, 1995
!=====================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      DIMENSION X(3)
      COMMON/COEFF_13/C(52)
      COMMON/THRBOD_13/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3
      COMMON/REFGEO_13/R10,R20,R30
!      COMMON/TEST_13/ITEST
!     ****************************************************************
      R1=ROO
      R2=RHO1
      R3=RHO2
      Q1=1.0D0/DSQRT(3.0D0)*(R1+R2+R3)
      Q2=1.0D0/DSQRT(2.0D0)*(R2-R3)
      Q3=1.0D0/DSQRT(6.0D0)*(2.0D0*R1-R2-R3)
      VHO2_13=THREBQ_13(Q1,Q2,Q3)+ &
        EXDIS_13(R1,R2,R3)+ELECT_13(R1,R2,R3)
      RETURN
      END


      FUNCTION THREBQ_13(Q1,Q2,Q3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/COEFF_13/C(52)
      COMMON/THRBOD_13/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3
      COMMON/REFGEO_13/R10,R20,R30
!     ****************************************************************
      Q12=Q1*Q1
      Q13=Q12*Q1
      Q14=Q13*Q1
      Q15=Q14*Q1
      Q16=Q15*Q1
      Q22=Q2*Q2
      Q32=Q3*Q3
      TQ1=Q22+Q32
      TQ2=Q32-3.0D0*Q22
      TQ3=Q22-Q32
      TQ12=TQ1*TQ1
      TQ13=TQ12*TQ1
      TQ22=TQ2*TQ2
      S1=R1-R10
      S2=R2-R20
      S3=R3-R30
      POLQ=C(1)*Q1+C(2)*Q12+C(3)*TQ1+C(4)*Q13+C(5)*Q1*TQ1+ &
      C(6)*Q3*TQ2+C(7)*Q14+C(8)*Q12*TQ1+C(9)*TQ1**2+C(10)*Q1*Q3*TQ2+ &
      C(11)*Q3+C(12)*Q1*Q3+C(13)*TQ3+C(14)*Q12*Q3+C(15)*Q1*TQ3+ &
      C(16)*Q3*TQ1+C(17)*Q13*Q3+C(18)*Q12*TQ3+C(19)*Q1*Q3*TQ1+ &
      C(20)*Q32*TQ2+C(21)*TQ1*TQ3+C(22)+C(23)*Q15+C(24)*Q13*TQ1+ &
      C(25)*Q1*TQ12+C(26)*Q12*Q3*TQ2+C(27)*Q3*TQ1*TQ2+C(28)*Q14*Q3+ &
      C(29)*Q13*TQ3+C(30)*Q12*Q3*TQ1+C(31)*Q1*Q32*TQ2+C(32)*Q1*TQ1*TQ3+ &
      C(33)*Q3*TQ12+C(34)*Q3*TQ2*TQ3+C(35)*Q16+C(36)*Q14*TQ1+ &
      C(37)*Q12*TQ12+C(38)*Q13*Q3*TQ2+C(39)*Q1*Q3*TQ1*TQ2+C(40)*TQ13+ &
      C(41)*Q32*TQ22+C(42)*Q15*Q3+C(43)*Q14*TQ3+C(44)*Q13*Q3*TQ1+ &
      C(45)*Q12*Q32*TQ2+C(46)*Q12*TQ1*TQ3+C(47)*Q1*Q3*TQ12+ &
      C(48)*Q1*Q3*TQ2*TQ3+C(49)*Q32*TQ1*TQ2+C(50)*TQ12*TQ3
      DECAY1=1.0D0-DTANH(C(51)*S1)
      DECAY2=1.0D0-DTANH(C(52)*S2)
      DECAY3=1.0D0-DTANH(C(52)*S3)
      THREBQ_13=POLQ*DECAY1*DECAY2*DECAY3
      RETURN
      END


      FUNCTION VOH_13(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ****************************************************************
      VOH_13=EHFOH_13(R)+DISOH_13(R)
      RETURN
      END

      FUNCTION EHFOH_13(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DATA D,ASV/0.13825385D0,2.6564788D0,1.7450528D0,0.71014391D0, &
                 2.5453276D0/
!     ****************************************************************
      X=R-RMOH
      R2=X*X
      R3=R2*X
      EHFOH_13=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*DEXP(-ASV(4)*X)
      RETURN
      END


      FUNCTION DISOH_13(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_13/COO(10),COH(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DISOH_13=DISP_13(R,COH(6),COH(8),COH(10),R0OH,RMOH)
      RETURN
      END


      FUNCTION VOO_13(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ****************************************************************
      VOO_13=EHFOO_13(R)+DISOO_13(R)
      RETURN
      END


      FUNCTION EHFOO_13(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DATA D,ASV/0.14291202D0,3.6445906D0,3.9281238D0,2.0986689D0, &
                 3.3522498D0/
!     ****************************************************************
      X=R-RMOO
      R2=X*X
      R3=R2*X
      EHFOO_13=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*DEXP(-ASV(4)*X)
      RETURN
      END


      FUNCTION DISOO_13(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_13/COO(10),COH(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DISOO_13=DISP_13(R,COO(6),COO(8),COO(10),R0OO,RMOO)
      RETURN
      END


      FUNCTION CEF_13(CAS,RK01,RK11,RK02,RK12,RE1,RE2,R1,R2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CEF_13=0.5D0*CAS*((1.0D0-RK01*DEXP(-RK11*(R1-RE1)))* &
       DTANH(RK12*R2)+(1.0D0-RK02*DEXP(-RK12*(R2-RE2)))*DTANH(RK11*R1))
      RETURN
      END


      FUNCTION EXDIS_13(R1,R2,R3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_13/COO(10),COH(10)
      COMMON/RKVAL_13/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DISCO_13/CEFOO(10),CEFOH2(10),CEFOH3(10), &
       CEDOO(10),CEDOH2(10)
      COMMON/DISCO2_13/CEDOH3(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DO 10 IN=6,10,2
      CEFOO(IN)=CEF_13(COO(IN),RK0OH(IN),RK1OH(IN),RK0OH(IN),RK1OH(IN), &
         RMOH,RMOH,R2,R3)
      CEDOO(IN)=CEFOO(IN)-COO(IN)
      CEFOH2(IN)=CEF_13(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN), &
         RMOO,RMOH,R1,R3)
      CEDOH2(IN)=CEFOH2(IN)-COH(IN)
      CEFOH3(IN)=CEF_13(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN), &
         RMOO,RMOH,R1,R2)
      CEDOH3(IN)=CEFOH3(IN)-COH(IN)
   10 CONTINUE
      EXDIS_13=DISP_13(R1,CEDOO(6),CEDOO(8),CEDOO(10),R0OO,RMOO) &
           +DISP_13(R2,CEDOH2(6),CEDOH2(8),CEDOH2(10),R0OH,RMOH) &
           +DISP_13(R3,CEDOH3(6),CEDOH3(8),CEDOH3(10),R0OH,RMOH)
      RETURN
      END


      FUNCTION ELECT_13(R1,R2,R3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/POLAR_13/C4,C5
      COMMON/RKVAL_13/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
      COMMON/DAMPC_13/ADAMP(10),BDAMP(10)
      COMMON/WELECT_13/C4OHR2,C5OHR2,C4OHR3,C5OHR3,C4OO,C5OO,TERM4,TERM5
!     ****************************************************************
      EPSD=1.0D-30
      IF (R1.LT.EPSD) R1=EPSD
        IF (R2.LT.EPSD) R2=EPSD
        IF (R3.LT.EPSD) R3=EPSD

      C42=C4
      C43=C4
      C52=C5
      C53=C5
      R23=R2**3
      R24=R23*R2
      R33=R3**3
      R34=R33*R3
      R14=R1**4
      R15=R14*R1
      R25=R24*R2
      R35=R34*R3
      RMQ=RMOH**4
      RMQ5=0.50D0/RMQ
      RMR3=RMQ5*R34
      RMR2=RMQ5*R24
      RMR33=RMQ5*R33
      RMR23=RMQ5*R23
      TAO=DTANH(RK1OO(4)*R1)
      TAH2=DTANH(RK1OH(4)*R2)
      TAH3=DTANH(RK1OH(4)*R3)
      EX3=DEXP(-RK1OH(4)*(R3-RMOH))
      EX2=DEXP(-RK1OH(4)*(R2-RMOH))
      R3E3=RMR3*EX3
      R2E2=RMR2*EX2
      CRE43=C4*R3E3
      CRE42=C4*R2E2
      CRE53=C5*R3E3
      CRE52=C5*R2E2
      C4OHR2=CRE43*TAO
      C5OHR2=CRE53*TAO
      C4OHR3=CRE42*TAO
      C5OHR3=CRE52*TAO
      C4OO=CRE43*TAH2+CRE42*TAH3
      C5OO=CRE53*TAH2+CRE52*TAH3
      RROH2=2.0D0*R2/(RMOH+2.5D0*R0OH)
      RROH3=2.0D0*R3/(RMOH+2.5D0*R0OH)
      RROO=2.0D0*R1/(RMOO+2.5D0*R0OO)
      TERM4=C4OO/R14*(1.0D0-DEXP(-ADAMP(4)*RROO-BDAMP(4)*RROO**2))**4+ &
       C4OHR2/R24*(1.0D0-DEXP(-ADAMP(4)*RROH2-BDAMP(4)*RROH2**2))**4+ &
       C4OHR3/R34*(1.0D0-DEXP(-ADAMP(4)*RROH3-BDAMP(4)*RROH3**2))**4
      TERM5=C5OO/R15*(1.0D0-DEXP(-ADAMP(5)*RROO-BDAMP(5)*RROO**2))**5+ &
       C5OHR2/R25*(1.0D0-DEXP(-ADAMP(5)*RROH2-BDAMP(5)*RROH2**2))**5+ &
       C5OHR3/R35*(1.0D0-DEXP(-ADAMP(5)*RROH3-BDAMP(5)*RROH3**2))**5
      ELECT_13=TERM4+TERM5
      RETURN
      END


      FUNCTION DISP_13(R,C6,C8,C10,R0,RM)
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DAMPC_13/ADAMP(10),BDAMP(10)
!     ****************************************************************
      IF (R.EQ.0.0D0) THEN
      DISP_13=0.0D0
      ELSE
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      D6=(1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6
      D8=(1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8
      D10=(1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10
      DISP_13=-C6/R6*D6-C8/R8*D8-C10/R10*D10
        ENDIF
        RETURN
      END

      function qderho2_13(r1,r2,r3,i)
        implicit real*8(a-h,o-z)
        call DERIVHO2_13(r1,r2,r3,d1,d2,d3)
        if(i.eq.1)then
        qderho2_13=d1

        else if(i.eq.2)then
        qderho2_13=d2

        else if(i.eq.3)then
        qderho2_13=d3

        else
        endif
        return
        end


      SUBROUTINE DERIVHO2_13(ROO,RHO1,RHO2,DD1,DD2,DD3)
!     ******************************************************************
!
!     TO COMPUTE THE HO2 SURFACE AND ITS DERIVATIVES
!
!     ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/COEFF_13/C(52)
      COMMON/THRBOD_13/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3
      COMMON/REFGEO_13/R10,R20,R30
      common/novo/NV,NJ,IREACT,IDEX,IE
!     ****************************************************************
      R1=ROO
      R2=RHO1
      R3=RHO2
      Q1=1.0D0/DSQRT(3.0D0)*(R1+R2+R3)
      Q2=1.0D0/DSQRT(2.0D0)*(R2-R3)
      Q3=1.0D0/DSQRT(6.0D0)*(2.0D0*R1-R2-R3)
      F=THREBQ_13(Q1,Q2,Q3)
      S1=R1-R10
      S2=R2-R20
      S3=R3-R30
      Q12=Q1*Q1
      Q13=Q12*Q1
      Q14=Q13*Q1
      Q15=Q14*Q1
      Q16=Q15*Q1
      Q22=Q2*Q2
      Q32=Q3*Q3
      TQ1=Q22+Q32
      TQ2=Q32-3.0D0*Q22
      TQ3=Q22-Q32
      TQ12=TQ1*TQ1
      TQ13=TQ12*TQ1
      TQ22=TQ2*TQ2
      DQ1R1=1.0D0/DSQRT(3.0D0)
      DQ1R2=DQ1R1
      DQ1R3=DQ1R1
      DQ2R1=0.0D0
      DQ2R2=1.0D0/DSQRT(2.0D0)
      DQ2R3=-DQ2R2
      DQ3R1=2.0D0/DSQRT(6.0D0)
      DQ3R2=-0.5D0*DQ3R1
      DQ3R3=DQ3R2
      DTQ2=2.0D0*Q2
      DTQ3=2.0D0*Q3
      DPOQ1=C(1)+2.0D0*C(2)*Q1+3.0D0*C(4)*Q12+C(5)*TQ1+4.0D0*C(7)*Q13+ &
       2.0D0*C(8)*Q1*TQ1+C(10)*Q3*TQ2+C(12)*Q3+C(14)*Q1*DTQ3+C(15)*TQ3+ &
        3.0D0*C(17)*Q12*Q3+2.0D0*C(18)*Q1*TQ3+C(19)*Q3*TQ1+5.0D0* &
       C(23)*Q14+3.0D0*C(24)*Q12*TQ1+C(25)*TQ12+C(26)*Q1*DTQ3*TQ2+4.0D0 &
        *C(28)*Q13*Q3+3.0D0*C(29)*Q12*TQ3+C(30)*Q1*DTQ3*TQ1+C(31)*Q32* &
        TQ2+C(32)*TQ1*TQ3+6.0D0*C(35)*Q15+4.0D0*C(36)*Q13*TQ1+2.0D0* &
        C(37)*Q1*TQ12+3.0D0*C(38)*Q12*Q3*TQ2+C(39)*Q3*TQ1*TQ2+5.0D0* &
        C(42)*Q14*Q3+4.0D0*C(43)*Q13*TQ3+3.0D0*C(44)*Q12*Q3*TQ1+2.0D0* &
        C(45)*Q1*Q32*TQ2+2.0D0*C(46)*Q1*TQ1*TQ3+C(47)*Q3*TQ12+C(48)*Q3* &
        TQ2*TQ3
      DPOQ2=C(3)*DTQ2+C(5)*Q1*DTQ2-C(6)*Q3*6.0D0*Q2+C(8)*Q12*DTQ2 &
         +C(9)*2.0D0*TQ1*DTQ2-C(10)*Q1*Q3*6.0D0*Q2+C(13)*DTQ2+ &
         C(15)*Q1*DTQ2+C(16)*Q3*DTQ2+C(18)*Q12*DTQ2+ &
         C(19)*Q1*Q3*DTQ2-C(20)*Q32*6.0D0*Q2+C(21)*DTQ2*TQ3+ &
         C(21)*TQ1*DTQ2+C(24)*Q13*DTQ2+C(25)*Q1*2.0D0*TQ1*DTQ2- &
         C(26)*Q12*Q3*6.0D0*Q2+C(27)*Q3*DTQ2*TQ2-C(27)*Q3*TQ1*6.0D0* &
         Q2+C(29)*Q13*DTQ2+C(30)*Q12*Q3*DTQ2-C(31)*Q1*Q32*6.0D0*Q2 &
         +C(32)*Q1*DTQ2*(TQ1+TQ3)+C(33)*DTQ3*TQ1*DTQ2+C(34)*Q3*DTQ2* &
         (TQ2-3.0D0*TQ3)+C(36)*Q14*DTQ2+4.0D0*C(37)*Q12*Q2*TQ1- &
         C(38)*Q13*Q3*6.0D0*Q2+C(39)*Q1*Q3*DTQ2*(TQ2-3.0D0*TQ1)+ &
         C(40)*3.0D0*TQ12*DTQ2-C(41)*Q32*TQ2*6.0D0*DTQ2+ &
         C(43)*Q14*DTQ2+C(44)*Q13*Q3*DTQ2-C(45)*Q12*Q32*6.0D0*Q2+ &
         C(46)*Q12*DTQ2*(TQ3+TQ1)+C(47)*Q1*DTQ3*TQ1*DTQ2+ &
        C(48)*Q1*Q3*DTQ2*(TQ2-3.0D0*TQ3)+C(49)*Q32*DTQ2*(TQ2-3.0D0*TQ1) &
         +C(50)*DTQ2*TQ1*(2.0D0*TQ3+TQ1)
      DPOQ3=C(3)*DTQ3+C(5)*Q1*DTQ3+C(6)*TQ2+C(6)*Q3*DTQ3+C(8)*Q12*DTQ3+ &
         C(9)*TQ1*2.0D0*DTQ3+C(10)*Q1*(TQ2+Q3*DTQ3)+C(11)+ &
         C(12)*Q1-C(13)*DTQ3+C(14)*Q12-C(15)*Q1*DTQ3+C(16)*(TQ1+ &
         Q3*DTQ3)+C(17)*Q13-C(18)*Q12*DTQ3+C(19)*Q1*(TQ1+Q3*DTQ3) &
         +C(20)*DTQ3*(TQ2+Q32)+C(21)*DTQ3*(TQ3-TQ1) &
        +C(24)*Q13*DTQ3+C(25)*Q1*2.0D0*TQ1*DTQ3+C(26)*Q12*(TQ2+Q3*DTQ3) &
         +C(27)*(TQ1*TQ2+Q3*DTQ3*(TQ2+TQ1))+C(28)*Q14-C(29)*Q13*DTQ3 &
         +C(30)*Q12*(TQ1+Q3*DTQ3)+C(31)*Q1*DTQ3*(TQ2+Q32)+C(32)*Q1* &
         DTQ3*(TQ3-TQ1)+C(33)*TQ1*(TQ1+DTQ3*DTQ3)+C(34)*TQ2*TQ3+ &
        C(34)*Q3*DTQ3*(TQ3-TQ2)+C(36)*Q14*DTQ3+C(37)*Q12*2.0D0*TQ1*DTQ3 &
         +C(38)*Q13*(TQ2+Q3*DTQ3)+C(39)*Q1*TQ1*TQ2+C(39)*Q1*Q3*DTQ3* &
        (TQ1+TQ2)+C(40)*3.0D0*TQ12*DTQ3+C(41)*DTQ3*TQ2*(TQ2+2.0D0*Q32)+ &
        C(42)*Q15-C(43)*Q14*DTQ3+C(44)*Q13*(TQ1+Q3*DTQ3)+C(45)*Q12*DTQ3 &
         *(TQ2+Q32)+C(46)*Q12*DTQ3*(TQ3-TQ1)+C(47)*Q1*(TQ12+TQ1*4.0D0* &
         Q32)+C(48)*Q1*TQ2*TQ3+C(48)*Q1*Q3*DTQ3*(TQ3-TQ2)+C(49)*DTQ3* &
        TQ1*TQ2+C(49)*Q32*DTQ3*(TQ1+TQ2)+C(50)*DTQ3*TQ1*(2.0D0*TQ3-TQ1)
      CALL  DEREXDIS_13(R1,R2,R3,DER1,DER2,DER3)
      CALL  DERELECT_13(R1,R2,R3,DERI1,DERI2,DERI3)
!      print *,'  r ', r1,r2,r3
!      print *,' c(51) ',c(51), '  s1 ',s1
      DD1=((DPOQ1*DQ1R1+DPOQ2*DQ2R1+DPOQ3*DQ3R1)*DECAY1-POLQ*C(51)/ &
         DCOSH(C(51)*S1)**2)*DECAY2*DECAY3+DER1+DERI1
      DD2=((DPOQ1*DQ1R2+DPOQ2*DQ2R2+DPOQ3*DQ3R2)*DECAY2-POLQ*C(52)/ &
         DCOSH(C(52)*S2)**2)*DECAY1*DECAY3+DER2+DERI2
      DD3=((DPOQ1*DQ1R3+DPOQ2*DQ2R3+DPOQ3*DQ3R3)*DECAY3-POLQ*C(52)/ &
         DCOSH(C(52)*S3)**2)*DECAY1*DECAY2+DER3+DERI3
      RETURN
      END


      SUBROUTINE DEREXDIS_13(R1,R2,R3,DER1,DER2,DER3)
!     ***************************************************************
!     TO COMPUTE THE DERIVATIVES OF THE EXCHANGE - DISPERSION TERM
!     ***************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CEFOO(10),CEFOH2(10),CEFOH3(10),CEDOO(10),CEDOH2(10)
      DIMENSION CEDOH3(10)
      COMMON/DISPC_13/COO(10),COH(10)
      COMMON/RKVAL_13/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
      COMMON/DAMPC_13/ADAMP(10),BDAMP(10)
!     ***************************************************************
      T1=R1-RMOO
      T5=R2-RMOH
      T6=R3-RMOH
      T7=0.5D0*COO(6)
      T8=0.5D0*COO(8)
      T9=0.5D0*COO(10)
      T10=0.5D0*COH(6)
      T11=0.5D0*COH(8)
      T12=0.5D0*COH(10)
      T13=RK1OH(6)*RK0OH(6)
      T14=RK1OH(8)*RK0OH(8)
      T15=RK1OH(10)*RK0OH(10)
      T16=RK1OO(6)*RK0OO(6)
      T18=RK1OO(8)*RK0OO(8)
      T19=RK1OO(10)*RK0OO(10)
      T20=DEXP(-RK1OH(6)*T5)
      T21=DEXP(-RK1OH(6)*T6)
      T22=DEXP(-RK1OH(8)*T5)
      T23=DEXP(-RK1OH(8)*T6)
      T24=DEXP(-RK1OH(10)*T5)
      T25=DEXP(-RK1OH(10)*T6)
      T26=DEXP(-RK1OO(6)*T1)
      T27=DEXP(-RK1OO(8)*T1)
      T28=DEXP(-RK1OO(10)*T1)
      T29=DCOSH(RK1OH(6)*R2)
      T30=DCOSH(RK1OH(6)*R3)
      T31=DCOSH(RK1OH(8)*R2)
      T32=DCOSH(RK1OH(8)*R3)
      T33=DCOSH(RK1OH(10)*R2)
      T34=DCOSH(RK1OH(10)*R3)
      T35=DCOSH(RK1OO(6)*R1)
      T36=DCOSH(RK1OO(8)*R1)
      T37=DCOSH(RK1OO(10)*R1)
      T38=DTANH(RK1OH(6)*R2)
      T39=DTANH(RK1OH(6)*R3)
      T40=DTANH(RK1OH(8)*R2)
      T41=DTANH(RK1OH(8)*R3)
      T42=DTANH(RK1OH(10)*R2)
      T43=DTANH(RK1OH(10)*R3)
      T44=DTANH(RK1OO(6)*R1)
      T45=DTANH(RK1OO(8)*R1)
      T46=DTANH(RK1OO(10)*R1)
      DO 10 IN=6,10,2
      CEFOO(IN)=CEF_13(COO(IN),RK0OH(IN),RK1OH(IN),RK0OH(IN),RK1OH(IN), &
             RMOH,RMOH,R2,R3)
      CEDOO(IN)=CEFOO(IN)-COO(IN)
      CEFOH2(IN)=CEF_13(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN), &
             RMOO,RMOH,R1,R3)
      CEDOH2(IN)=CEFOH2(IN)-COH(IN)
      CEFOH3(IN)=CEF_13(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN), &
             RMOO,RMOH,R1,R2)
      CEDOH3(IN)=CEFOH3(IN)-COH(IN)
   10 CONTINUE
      RR1=2.0D0*R1/(RMOO+2.5D0*R0OO)
      T6R1=((1.0D0-EXP(-ADAMP(6)*RR1-BDAMP(6)*RR1**2))/R1)**6
      T8R1=((1.0D0-EXP(-ADAMP(8)*RR1-BDAMP(8)*RR1**2))/R1)**8
      T10R1=((1.0D0-EXP(-ADAMP(10)*RR1-BDAMP(10)*RR1**2))/R1)**10
      RR2=2.0D0*R2/(RMOH+2.5D0*R0OH)
      T6R2=((1.0D0-EXP(-ADAMP(6)*RR2-BDAMP(6)*RR2**2))/R2)**6
      T8R2=((1.0D0-EXP(-ADAMP(8)*RR2-BDAMP(8)*RR2**2))/R2)**8
      T10R2=((1.0D0-EXP(-ADAMP(10)*RR2-BDAMP(10)*RR2**2))/R2)**10
      RR3=2.0D0*R3/(RMOH+2.5D0*R0OH)
      T6R3=((1.0D0-EXP(-ADAMP(6)*RR3-BDAMP(6)*RR3**2))/R3)**6
      T8R3=((1.0D0-EXP(-ADAMP(8)*RR3-BDAMP(8)*RR3**2))/R3)**8
      T10R3=((1.0D0-EXP(-ADAMP(10)*RR3-BDAMP(10)*RR3**2))/R3)**10
      DC61R2=T7*(T13 &
            *T20*T39+ &
            (1.0D0-RK0OH(6)*T21) &
            *RK1OH(6)/T29**2)
      DC61R3=T7*(T13 &
            *T21*T38+ &
            (1.0D0-RK0OH(6)*T20) &
            *RK1OH(6)/T30**2)
      DC81R2=T8*(T14 &
            *T22*T41+ &
            (1.0D0-RK0OH(8)*T23) &
            *RK1OH(8)/T31**2)
      DC81R3=T8*(T14 &
            *T23*T40+ &
            (1.0D0-RK0OH(8)*T22) &
            *RK1OH(8)/T32**2)
      D101R2=T9*(T15 &
            *T24*T43+ &
            (1.0D0-RK0OH(10)*T25) &
            *RK1OH(10)/T33**2)
      D101R3=T9*(T15 &
            *T25*T42+ &
            (1.0D0-RK0OH(10)*T24) &
            *RK1OH(10)/T34**2)
      DC62R1=T10*(T16 &
            *T26*T39+ &
            (1.0D0-RK0OH(6)*T21) &
            *RK1OO(6)/T35**2)
      DC62R3=T10*(T13 &
            *T21*T44+ &
            (1.0D0-RK0OO(6)*T26) &
            *RK1OH(6)/T30**2)
      DC82R1=T11*(T18 &
            *T27*T41+ &
            (1.0D0-RK0OH(8)*T23) &
            *RK1OO(8)/T36**2)
      DC82R3=T11*(T14 &
            *T23*T45+ &
            (1.0D0-RK0OO(8)*T27) &
            *RK1OH(8)/T32**2)
      D102R1=T12*(T19 &
            *T28*T43+ &
            (1.0D0-RK0OH(10)*T25) &
            *RK1OO(10)/T37**2)
      D102R3=T12*(T15 &
            *T25*T46+ &
            (1.0D0-RK0OO(10)*T28) &
            *RK1OH(10)/T34**2)
      DC63R1=T10*(T16 &
            *T26*T38+ &
            (1.0D0-RK0OH(6)*T20) &
            *RK1OO(6)/T35**2)
      DC63R2=T10*(T13 &
            *T20*T44+ &
            (1.0D0-RK0OO(6)*T26) &
            *RK1OH(6)/T29**2)
      DC83R1=T11*(T18 &
            *T27*T40+ &
            (1.0D0-RK0OH(8)*T22) &
            *RK1OO(8)/T36**2)
      DC83R2=T11*(T14 &
            *T22*T45+ &
            (1.0D0-RK0OO(8)*T27) &
            *RK1OH(8)/T31**2)
      D103R1=T12*(T19 &
            *T28*T42+ &
            (1.0D0-RK0OH(10)*T24) &
            *RK1OO(10)/T37**2)
      D103R2=T12*(T15 &
            *T24*T46+ &
            (1.0D0-RK0OO(10)*T28) &
            *RK1OH(10)/T33**2)
      DER1=DDISP_13(R1,CEDOO(6),CEDOO(8),CEDOO(10),R0OO,RMOO)-DC62R1* &
               T6R2-DC63R1*T6R3-DC82R1*T8R2-DC83R1*T8R3-D102R1*T10R2- &
           D103R1*T10R3
      DER2=DDISP_13(R2,CEDOH2(6),CEDOH2(8),CEDOH2(10),R0OH,RMOH)-DC61R2* &
               T6R1-DC63R2*T6R3-DC81R2*T8R1-DC83R2*T8R3-D101R2*T10R1- &
           D103R2*T10R3
      DER3=DDISP_13(R3,CEDOH3(6),CEDOH3(8),CEDOH3(10),R0OH,RMOH)-DC61R3* &
               T6R1-DC62R3*T6R2-DC81R3*T8R1-DC82R3*T8R2-D101R3*T10R1- &
           D102R3*T10R2
      RETURN
      END

      SUBROUTINE DERELECT_13(R1,R2,R3,DER1,DER2,DER3)
!     ****************************************************************
!     TO COMPUTE THE DERIVATIVES OF THE ELECTROSTATIC TERM
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/POLAR_13/C4,C5
      COMMON/RKVAL_13/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
      COMMON/DAMPC_13/ADAMP(10),BDAMP(10)
!     ****************************************************************
      C42=C4
      C43=C4
      C52=C5
      C53=C5
      RR1=2.0D0*R1/(RMOO+2.5D0*R0OO)
      D4R1=1.0D0-DEXP(-ADAMP(4)*RR1-BDAMP(4)*RR1*RR1)
      D5R1=1.0D0-DEXP(-ADAMP(5)*RR1-BDAMP(5)*RR1*RR1)
      T4R1=(D4R1/R1)**4
      T5R1=(D5R1/R1)**5
      RR2=2.0D0*R2/(RMOH+2.5D0*R0OH)
      D4R2=1.0D0-DEXP(-ADAMP(4)*RR2-BDAMP(4)*RR2*RR2)
      D5R2=1.0D0-DEXP(-ADAMP(5)*RR2-BDAMP(5)*RR2*RR2)
      T4R2=(D4R2/R2)**4
      T5R2=(D5R2/R2)**5
      RR3=2.0D0*R3/(RMOH+2.5D0*R0OH)
      D4R3=1.0D0-DEXP(-ADAMP(4)*RR3-BDAMP(4)*RR3*RR3)
      D5R3=1.0D0-DEXP(-ADAMP(5)*RR3-BDAMP(5)*RR3*RR3)
      T4R3=(D4R3/R3)**4
      T5R3=(D5R3/R3)**5
      R23=R2**3
      R24=R23*R2
      R33=R3**3
      R34=R33*R3
      R14=R1**4
      R15=R14*R1
      R25=R24*R2
      R35=R34*R3
      RMQ=RMOH**4
      RMQ5=0.50D0/RMQ
      RMR3=RMQ5*R34
      RMR2=RMQ5*R24
      RMR33=RMQ5*R33
      RMR23=RMQ5*R23
      TAO=DTANH(RK1OO(4)*R1)
      TAH2=DTANH(RK1OH(4)*R2)
      TAH3=DTANH(RK1OH(4)*R3)
      EX3=DEXP(-RK1OH(4)*(R3-RMOH))
      EX2=DEXP(-RK1OH(4)*(R2-RMOH))
      COR=1.0D0/DCOSH(RK1OO(4)*R1)
      COR2=1.0D0/DCOSH(RK1OH(4)*R2)
      COR3=1.0D0/DCOSH(RK1OH(4)*R3)
      R3E3=RMR3*EX3
      R2E2=RMR2*EX2
      RKOCC=RK1OO(4)*COR*COR
      RKHC3=RK1OH(4)*COR3*COR3
      RKHC2=RK1OH(4)*COR2*COR2
      CRE43=C4*R3E3
      CRE42=C4*R2E2
      CRE53=C5*R3E3
      CRE52=C5*R2E2
      DIVR3=4.0D0/R3-RK1OH(4)
      DIVR2=4.0D0/R2-RK1OH(4)
      DIF3=4.0D0-R3*RK1OH(4)
      DIF2=4.0D0-R2*RK1OH(4)
      C4OHR2=CRE43*TAO
      C5OHR2=CRE53*TAO
      C4OHR3=CRE42*TAO
      C5OHR3=CRE52*TAO
      C4OO=CRE43*TAH2+CRE42*TAH3
      C5OO=CRE53*TAH2+CRE52*TAH3
      DC42R3=C4OHR2*DIVR3
      DC52R3=C5OHR2*DIVR3
      DC43R2=C4OHR3*DIVR2
      DC53R2=C5OHR3*DIVR2
      DC42R1=CRE43*RKOCC
      DC52R1=CRE53*RKOCC
      DC43R1=CRE42*RKOCC
      DC53R1=CRE52*RKOCC
      DC42R2=0.0D0
      DC52R2=0.0D0
      DC43R3=0.0D0
      DC53R3=0.0D0
      DC41R3=C4*RMR33*EX3*TAH2*DIF3+C4*RMR2*EX2*RKHC3
      DC41R2=C4*RMR23*EX2*TAH3*DIF2+C4*RMR3*EX3*RKHC2
      DC51R3=C5*RMR33*EX3*TAH2*DIF3+C5*RMR2*EX2*RKHC3
      DC51R2=C5*RMR23*EX2*TAH3*DIF2+C5*RMR3*EX3*RKHC2
      DC41R1=0.0D0
      DC51R1=0.0D0
      DRR1=RR1/R1
      DRR2=RR2/R2
      DRR3=RR3/R3
      DDISP1=-4.0D0*C4OO/R14*D4R1**3*(D4R1/R1+(1.0D0-D4R1)*(-ADAMP(4)- &
             2.0D0*BDAMP(4)*RR1)*DRR1) &
            -5.0D0*C5OO/R15*D5R1**4*(D5R1/R1+(1.0D0-D5R1)*(-ADAMP(5)- &
             2.0D0*BDAMP(5)*RR1)*DRR1)
      DDISP2=-4.D0*C4OHR2/R24*D4R2**3*(D4R2/R2+(1.0D0-D4R2)*(-ADAMP(4)- &
             2.0D0*BDAMP(4)*RR2)*DRR2) &
            -5.0D0*C5OHR2/R25*D5R2**4*(D5R2/R2+(1.0D0-D5R2)*(-ADAMP(5)- &
             2.0D0*BDAMP(5)*RR2)*DRR2)
      DDISP3=-4.D0*C4OHR3/R34*D4R3**3*(D4R3/R3+(1.0D0-D4R3)*(-ADAMP(4)- &
             2.0D0*BDAMP(4)*RR3)*DRR3) &
            -5.0D0*C5OHR3/R35*D5R3**4*(D5R3/R3+(1.0D0-D5R3)*(-ADAMP(5)- &
             2.0D0*BDAMP(5)*RR3)*DRR3)
      DER1=DDISP1+DC42R1*T4R2+DC43R1*T4R3+DC52R1*T5R2+DC53R1*T5R3 &
           +DC41R1*T4R1+DC51R1*T5R1
      DER2=DDISP2+DC41R2*T4R1+DC43R2*T4R3+DC51R2*T5R1+DC53R2*T5R3 &
           +DC42R2*T4R2+DC52R2*T5R2
      DER3=DDISP3+DC42R3*T4R2+DC41R3*T4R1+DC52R3*T5R2+DC51R3*T5R1 &
           +DC43R3*T4R3+DC53R3*T5R3
      RETURN
      END


      FUNCTION DVOH_13(R)
!     ****************************************************************
!
!     TO COMPUTE THE DER. OF THE HFACE FOR  O...H
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DISPC_13/COO(10),COH(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DATA D,ASV/0.13825385D0,2.6564788D0,1.7450528D0,0.71014391D0, &
        2.5453276D0/
!     ****************************************************************
      X=R-RMOH
      R2=X*X
      R3=R2*X
      POL=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)
      DPOL=-D*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*R2)
      POT=DEXP(-ASV(4)*X)
      DVOH_13=-ASV(4)*POT*POL+DPOL*POT+DDISP_13(R,COH(6), &
         COH(8),COH(10),R0OH,RMOH)
      RETURN
      END

      FUNCTION DVOO_13(R)
!     ****************************************************************
!
!     TO COMPUTE THE DER. OF THE HFACE FOR  O...O
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DISPC_13/COO(10),COH(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
!     ****************************************************************
      DATA D,ASV/0.14291202D0,3.6445906D0,3.9281238D0,2.0986689D0, &
                 3.3522498D0/
!     ****************************************************************
      X=R-RMOO
      R2=X*X
      R3=R2*X
      POL=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)
      DPOL=-D*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*R2)
      POT=DEXP(-ASV(4)*X)
      DVOO_13=-ASV(4)*POT*POL+DPOL*POT+DDISP_13(R,COO(6), &
         COO(8),COO(10),R0OO,RMOO)
      RETURN
      END

      FUNCTION DDISP_13(R,C6,C8,C10,R0,RM)
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DAMPC_13/ADAMP(10),BDAMP(10)
!     ****************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      DRR=RR/R
      T6=1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR)
      T8=1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR)
      T10=1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR)
      DDISP_13=6.D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-ADAMP(6)- &
            2.0D0*BDAMP(6)*RR)*DRR)+ &
            8.D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-ADAMP(8)-2.0D0*BDAMP(8) &
            *RR)*DRR)+ &
            10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-ADAMP(10)-2.0D0* &
           BDAMP(10)*RR)*DRR)
      RETURN
      END





      BLOCK DATA HO2DAT_13
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/COEFF_13/C(52)
      COMMON/DISPC_13/COO(10),COH(10)
      COMMON/DIATDI_13/R0OO,RMOO,R0OH,RMOH
      COMMON/RKVAL_13/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/POLAR_13/C4,C5
      COMMON/DAMPC_13/ADAMP(10),BDAMP(10)
      COMMON/REFGEO_13/R10,R20,R30
!     ***************************************************************


      DATA C/ &
        .49040645D+01, -.86748216D+01,  .50555792D+01,  .42941301D+01, &
       -.41874792D+01,  .13461379D+00, -.99064922D+00,  .13358488D+01, &
        .13495231D+01, -.18529696D+00, -.23534213D+02,  .24289930D+02, &
       -.50209026D+01, -.10365484D+02,  .46692224D+01, -.14747138D+01, &
        .23119718D+01, -.18247842D+01, -.28472166D+00,  .51036509D+00, &
        .19124083D+00,  .45405729D+01,  .11087611D+00, -.19990481D+00, &
       -.37356178D+00,  .46142042D-01, -.20565580D+00, -.27015963D+00, &
        .34085281D+00,  .28321162D+00, -.11558481D+00, -.29448886D+00, &
       -.52932488D+00,  .58159523D-01, -.48649560D-02,  .11949167D-01, &
        .21409804D-01, -.20620608D-02,  .30177088D-01,  .27880291D-01, &
        .88458711D-02,  .13137410D-01, -.24705619D-01, -.31085889D-01, &
        .34317857D-02,  .52593878D-01,  .79500714D-01, -.79782216D-02, &
        .31164575D-01, -.28737598D-01,  .98201698D+00,  .62000000D+00/



      DATA R0OO,RMOO,R0OH,RMOH/5.661693D0,2.2818D0,6.294894D0,1.8344D0/
      DATA COO/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,15.40D0,0.0D0,235.219943D0, &
               0.0D0,4066.23929D0/
      DATA COH/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,10.00D0,0.0D0,180.447673D0, &
               0.0D0,3685.25842D0/
      DATA C4,C5/-0.92921D0,-1.79000D0/
      DATA RK0OO/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,-.27847758D0,0.0D0, &
                 -.46815641D0,0.0D0,-1.20506384D0/
      DATA RK1OO/0.0D0,0.0D0,0.0D0,3.35224980D0,3.35224980D0, &
                 0.95273753D0,0.0D0,0.94148408D0,0.0D0,0.72379129D0/
      DATA RK0OH/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.02465005D0,0.0D0, &
                 0.05036950D0,0.0D0,0.06294371D0/
      DATA RK1OH/0.0D0,0.0D0,0.0D0,2.54532760D0,2.54532760D0, &
                 0.68758097D0,0.0D0,0.82542359D0,0.0D0,0.94034225D0/
      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0951333D0, &
                 0.0D0,2.1999000D0,0.0D0,1.6880714D0/
      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7787895D0, &
                 0.0D0,7.2265123D0,0.0D0,5.9487108D0/
      DATA R10,R20,R30/2.5143000D0,2.6469057D0,2.6469057D0/
      END

!      *****************************************************************
!      *                  POTENCIAL O4                                 *
!      *****************************************************************

      FUNCTION VO4_13(R1,R2,R3,R4,R5,R6)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DM1_13/D1(4),D2(4)
      COMMON/TA_13/FPE,GAMATA,SUMTA
      COMMON/TA2_13/FPE2,GAMAT2,SUMT2
      COMMON/PC_13/CC(8)
      DATA RM,CRED,C5/2.2818,142.105,5462.358/
      SUM=R1+R2+R3+R4+R5+R6
      redK=sum-SUMTA
      REDV=SUM-SUMT2
      TK=FPE*EXP(-GAMATA*redK**2)
      TV=FPE2*EXP(-GAMAT2*REDV**2)
      SUM1=R2+R4+R5+R6
      SUM2=R1+R3+R5+R6
      SUM3=R1+R2+R3+R4
      x1=sum1*2./CRED
      x2=sum2*2./CRED
      x3=sum3*2./CRED
      E1=-((R1-RM)**2+(R3-RM)**2)
      E2=-((R2-RM)**2+(R4-RM)**2)
      E3=-((R5-RM)**2+(R6-RM)**2)
      DMP1=(1.-exp(-D1(4)*X1*(1.d0+D2(4)*X1)))**5
      DMP2=(1.-exp(-D1(4)*X2*(1.d0+D2(4)*X2)))**5
      DMP3=(1.-exp(-D1(4)*X3*(1.d0+D2(4)*X3)))**5
      SUMT=EXP(E1)*DMP1/SUM1**5+EXP(E2)*DMP2/SUM2**5 &
           +EXP(E3)*DMP3/SUM3**5
      VO4_13=TV+TK &
         +VO2_13(R2)+VO2_13(R4)+VO2_13(R5)+VO2_13(R6) &
         +VO3_13(R1,R2,R5)+VO3_13(R1,R4,R6)+VO3_13(R3,R4,R5)+ &
         VO3_13(R6,R2,R3) +C5*SUMT &
         +VO2_13(R1)+VO2_13(R3)
      RETURN
      END

      FUNCTION DVO4_13(R1,R2,R3,R4,R5,R6,I)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/TA_13/FPE,GAMATA,SUMTA
      COMMON/TA2_13/FPE2,GAMAT2,SUMT2
      SUM=R1+R2+R3+R4+R5+R6
      REDK=SUMTA-SUM
      REDV=SUMT2-SUM
      DTK=FPE*EXP(-GAMATA*REDK**2)*GAMATA*2.*REDK
      DTV=FPE2*EXP(-GAMAT2*REDV**2)*GAMAT2*2.*REDV
      DTA=DTK+DTV
      GOTO(1,2,3,4,5,6)I
    1 DVO4_13=DVO2_13(R1)+DVO3_13(R1,R2,R5,1)+DVO3_13(R1,R4,R6,1)+DTA+ &
       DERED_13(R1,R2,R3,R4,R5,R6)
      GOTO 7
    2 DVO4_13=DVO2_13(R2)+DVO3_13(R1,R2,R5,2)+DVO3_13(R6,R2,R3,2)+DTA+ &
       DERED_13(R2,R1,R4,R3,R5,R6)
      GOTO 7
    3 DVO4_13=DVO2_13(R3)+DVO3_13(R3,R4,R5,1)+DVO3_13(R6,R2,R3,3)+DTA+ &
       DERED_13(R3,R2,R1,R4,R5,R6)
      GOTO 7
    4 DVO4_13=DVO2_13(R4)+DVO3_13(R1,R4,R6,2)+DVO3_13(R3,R4,R5,2)+DTA+ &
       DERED_13(R4,R1,R2,R3,R5,R6)
      GOTO 7
    5 DVO4_13=DVO2_13(R5)+DVO3_13(R1,R2,R5,3)+DVO3_13(R3,R4,R5,3)+DTA+ &
       DERED_13(R5,R2,R6,R4,R1,R3)
      GOTO 7
    6 DVO4_13=DVO2_13(R6)+DVO3_13(R1,R4,R6,3)+DVO3_13(R6,R2,R3,1)+DTA+ &
       DERED_13(R6,R2,R5,R4,R1,R3)
    7 RETURN
      END

      FUNCTION VO2_13(R)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MAIN_13/CC(3),A(3),B(3),DD
      COMMON /V2MAS_13/VTOT,DVTOT,D2VTOT,RCSN
      RCSN=R
      CALL V2_13
      ACE=0.0
      DO 1 I=1,3
      ACE=ACE+CC(I)*DMR_13(I,R)
    1 CONTINUE
      VO2_13=ACE+VTOT
      RETURN
      END

      FUNCTION DVO2_13(R)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MAIN_13/CC(3),A(3),B(3),DD
      COMMON /V2MAS_13/VTOT,DVTOT,D2VTOT,RCSN
      RCSN=R
      CALL DV2_13
      DACE=0.0
      DO 1 I=1,3
      DACE=DACE+CC(I)*DDMR_13(I,R)
    1 CONTINUE
      DVO2_13=DACE+DVTOT
      RETURN
      END

      FUNCTION VO3_13(R1,R2,R3)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(3)
      EPSD=1.0D-20
      R(1)=R1
      R(2)=R2
      R(3)=R3
      DO 1 I=1,3
 1      IF (R(I).LT.EPSD) R(I)=EPSD
      CALL PRINC_13(R,YR)
      aaa=VAD_13(r1,r2,r3)
      bbb=PTH_13(r1,r2,r3)
      VO3_13=YR+aaa+bbb
      RETURN
      END

      FUNCTION DVO3_13(R1,R2,R3,I)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(3)
      R(1)=R1
      R(2)=R2
      R(3)=R3
      CALL DPRINC_13(R,DY,I)
      aaa=DVAD_13(r1,r2,r3,i)
      bbb=DPTH_13(r1,r2,r3,i)
      DVO3_13=DY+aaa+bbb
      RETURN
      END


      FUNCTION DVO4N_13(R1,R2,R3,R4,R5,R6,I)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(6)
      R(1)=R1
      R(2)=R2
      R(3)=R3
      R(4)=R4
      R(5)=R5
      R(6)=R6
      R(I)=R(I)+1.E-7
      VSUP=VO4_13(R(1),R(2),R(3),R(4),R(5),R(6))
      R(I)=R(I)-2.E-7
      VINF=VO4_13(R(1),R(2),R(3),R(4),R(5),R(6))
      DVO4N_13=(VSUP-VINF)/2.E-7
      RETURN
      END


      FUNCTION DERED_13(R1,R2,R3,R4,R5,R6)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DM1_13/D1(4),D2(4)
      DATA RM,CRED,C5/2.2818,142.105,5462.358/
      SUM1=R2+R4+R5+R6
      SUM2=R1+R3+R5+R6
      SUM3=R1+R2+R3+R4
      x1=sum1*2./CRED
      x2=sum2*2./CRED
      x3=sum3*2./CRED
      E1=-((R1-RM)**2+(R3-RM)**2)
      E2=-((R2-RM)**2+(R4-RM)**2)
      E3=-((R5-RM)**2+(R6-RM)**2)
      DMP1=(1.-exp(-D1(4)*X1*(1.d0+D2(4)*X1)))**5
      DX2=2.D0/CRED
      DX3=DX2
      POL2=-D1(4)*X2*(1.D0+D2(4)*X2)
      POL3=-D1(4)*X3*(1.D0+D2(4)*X3)
      DPOL2=-D1(4)*(1.D0+2.D0*D2(4)*X2)
      DPOL3=-D1(4)*(1.D0+2.D0*D2(4)*X3)
      EXPO2=EXP(POL2)
      EXPO3=EXP(POL3)
      RL2=SUM2**5
      RL3=SUM3**5
      C2=(1.-EXPO2)**5
      C3=(1.-EXPO3)**5
      EXPOEL2=C2/(1.-EXPO2)
      EXPOEL3=C3/(1.-EXPO3)
      DC2=5.D0*EXPOEL2*(-EXPO2)*DPOL2*DX2
      DC3=5.D0*EXPOEL3*(-EXPO3)*DPOL3*DX3
      DDMR2=DC2/RL2-5.D0*C2/(RL2*SUM2)
      DDMR3=DC3/RL3-5.D0*C3/(RL3*SUM3)
      DERED_13=C5*(2./SUM1**5*EXP(E1)*DMP1*(RM-R1)+ &
                DDMR2*EXP(E2)+DDMR3*EXP(E3))
      RETURN
      END



      FUNCTION DDVO4_13(R,ID,JD)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(6)
      R(JD)=R(JD)+1.E-5
      VSUP=DVO4_13(R(1),R(2),R(3),R(4),R(5),R(6),ID)
      R(JD)=R(JD)-2.E-5
      VINF=DVO4_13(R(1),R(2),R(3),R(4),R(5),R(6),ID)
      R(JD)=R(JD)+1.E-5
      DDVO4_13=(VSUP-VINF)/2.E-5
      RETURN
      END

      SUBROUTINE PRINC_13(RR,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BFF_13/VHF,ACE
      COMMON/COMN_13/VHF2,ACE2
      COMMON/PART2_13/VEQ,DVEQ,D2VEQ
      COMMON/MAIN_13/C2C(3),A0(3),A1(3),DD
      COMMON/COORD_13/R(3)
      COMMON/PART_13/VDP,DVDP,D2VDP
      COMMON/V2MAS_13/VTOT,DVTOT,D2VTOT,RCSN
      COMMON/MASMAIN_13/VR2,DVR2,D2VR2
      COMMON/CONTROL_13/K1
      COMMON/VVRR_13/ZP(3)
      DIMENSION RR(3)
      DO 5 I=1,3
      R(I)=RR(I)
    5 CONTINUE
      CALL CE_13
      CALL CEQ_13
      Y=VDP+VEQ
      RETURN
      END

      SUBROUTINE DPRINC_13(RR,DY,K1D)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BFF_13/VHF,ACE
      COMMON/COMN_13/VHF2,ACE2
      COMMON/PART2_13/VEQ,DVEQ,D2VEQ
      COMMON/MAIN_13/C2C(3),A0(3),A1(3),DD
      COMMON/COORD_13/R(3)
      COMMON/PART_13/VDP,DVDP,D2VDP
      COMMON/V2MAS_13/VTOT,DVTOT,D2VTOT,RCSN
      COMMON/MASMAIN_13/VR2,DVR2,D2VR2
      COMMON/CONTROL_13/K1
      COMMON/VVRR_13/ZP(3)
      DIMENSION RR(3)
      K1=K1D
      DO 5 I=1,3
      R(I)=RR(I)
    5 CONTINUE
      CALL DCE_13
      RCSN=R(K1D)
      CALL DCEQ_13
      DY=DVDP+DVEQ
      RETURN
      END

      SUBROUTINE V2_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/BFF_13/VHF,ACE
      COMMON/START_13/A,B,AK0,AK1,RZ,CXX
      COMMON/V2MAS_13/VTOT,DVTOT,D2VTOT,R
      COMMON/MAIN_13/CC(3),A00(3),A11(3),DD
      COMMON/DM1_13/D1(4),D2(4)
      COMMON/DM2_13/R0,RM
      COMMON/DM3_13/GAMA
      COMMON/DM4_13/A1,A2,A3,GM,RE,DE
      RR=R-RE
      EGM=EXP(-GM*RR)
      POL=DE*(1.+A1*RR+A2*RR**2+A3*RR**3)
      VTOT=POL*EGM
      RETURN
      END

      SUBROUTINE DV2_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/START_13/A,B,AK0,AK1,RZ,CXX
      COMMON/V2MAS_13/VTOT,DVTOT,D2VTOT,R
      COMMON/MAIN_13/CC(3),A00(3),A11(3),DD
      COMMON/DM1_13/D1(4),D2(4)
      COMMON/DM2_13/R0,RM
      COMMON/DM3_13/GAMA
      COMMON/DM4_13/A1,A2,A3,GM,RE,DE
      RR=R-RE
      EGM=EXP(-GM*RR)
      POL=DE*(1.+A1*RR+A2*RR**2+A3*RR**3)
      DPOL=DE*(A1+2.*A2*RR+3.*A3*RR**2)
      DVTOT=(DPOL*EGM-GM*POL*EGM)
      RETURN
      END

      SUBROUTINE CE_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CECF_13/CC,DCC,D2CC
      COMMON/MAIN_13/C2C(3),A0(3),A1(3),DD
      COMMON/START_13/A,B,AK0,AK1,R0,CXX
      COMMON/CONTROL_13/K1D
      COMMON/PART_13/VDP,DVDP,D2VDP
      COMMON/DUMMY_13/INDEX(3,2)
      COMMON/COORD_13/R(3)
      VDP=0.0
      DO 10 J=1,3
      I1=INDEX(J,1)
      I2=INDEX(J,2)
      A=R(I1)
      B=R(I2)
      DO 10 I=1,3
      CXX=C2C(I)
      AK0=A0(I)
      AK1=A1(I)
      CALL CF_13
      VDP=VDP+(CC-CXX)*DMR_13(I,R(J))
   10 CONTINUE
      RETURN
      END

      SUBROUTINE DCE_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CECF_13/CC,DCC,D2CC
      COMMON/MAIN_13/C2C(3),A0(3),A1(3),DD
      COMMON/START_13/A,B,AK0,AK1,R0,CXX
      COMMON/CONTROL_13/K1D
      COMMON/PART_13/VDP,DVDP,D2VDP
      COMMON/DUMMY_13/INDEX(3,2)
      COMMON/DMY_13/IX(3,3)
      COMMON/COORD_13/R(3)
      DVDP=0.0
      DO 10 I=1,3
      J=K1D
      I1=INDEX(J,1)
      I2=INDEX(J,2)
      A=R(I1)
      B=R(I2)
      CXX=C2C(I)
      AK0=A0(I)
      AK1=A1(I)
      CALL CF_13
      DVDP=DVDP+(CC-CXX)*DDMR_13(I,R(J))
      DO 10 K=1,2
      J=INDEX(K1D,K)
      IO=IX(K1D,J)
      A=R(K1D)
      B=R(IO)
      CALL DCF_13
      DVDP=DVDP+DCC*DMR_13(I,R(J))
   10 CONTINUE
      RETURN
      END

      SUBROUTINE CF_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/START_13/A,B,AK0,AK1,R0,CXX
      COMMON/CECF_13/CC,DCC,D2CC
      TANHA=TANH(AK1*A)
      TANHB=TANH(AK1*B)
      EA=AK0*EXP(-AK1*(A-R0))
      EB=AK0*EXP(-AK1*(B-R0))
      CC=0.5*CXX*((1.+EA)*TANHB+(1.+EB)*TANHA)
      RETURN
      END

      SUBROUTINE DCF_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/START_13/A,B,AK0,AK1,R0,CXX
      COMMON/CECF_13/CC,DCC,D2CC
      EA=AK0*EXP(-AK1*(A-R0))
      EB=AK0*EXP(-AK1*(B-R0))
      AKA=AK1*A
      TANHA=TANH(AKA)
      TANHB=TANH(AK1*B)
      SINHA=SINH(AKA)
      COSHA=SINHA/TANHA
      CSHA2=1./COSHA*1./COSHA
      DCC=0.5*CXX*AK1*(-EA*TANHB+(1.+EB)*CSHA2)
      RETURN
      END

      FUNCTION DMR_13(I,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PN_13/N(4)
      COMMON/DM1_13/D1(4),D2(4)
      COMMON/DM2_13/R0,RM
      COMMON/DM3_13/GAMA
        EPSD=1.0D-30
      IF (R.LT.EPSD) R=EPSD
       RO=(RM+GAMA*R0)/2.
      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      RL=R**IDXL
      DMR_13=-C/RL
      RETURN
      END

      FUNCTION DDMR_13(I,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/PN_13/N(4)
      COMMON/DM1_13/D1(4),D2(4)
      COMMON/DM2_13/R0,RM
      COMMON/DM3_13/GAMA
      RO=(RM+GAMA*R0)/2.
      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      DPOL=-D1(I)*(1.+2.*D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      EXPOEL=C/(1.-EXPO)
      DC=DXL*EXPOEL*(-EXPO)*DPOL*DX
      RL=R**IDXL
      DDMR_13=-(DC/RL-DXL*C/(RL*R))
      RETURN
      END

      SUBROUTINE CEQ_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MAIN2_13/C5Q,RK
      COMMON/START2_13/C
      COMMON/PART2_13/VEQ,DVEQ,D2VEQ
      COMMON/CECF_13/CC,DCC,D2CC
      COMMON/CECF2_13/DCCM
      COMMON/MAIN_13/C2C(3),A0(3),A1(3),DD
      COMMON/START_13/A,B,AK0,AK1,R0,CXX
      COMMON/CONTROL_13/K1D
      COMMON/PART_13/VDP,DVDP,D2VDP
      COMMON/DUMMY_13/INDEX(3,2)
      COMMON/COORD_13/R(3)
      VEQ=0.0
      DO 10 J=1,3
      I1=INDEX(J,1)
      I2=INDEX(J,2)
      A=R(I1)
      B=R(I2)
      C=R(J)
      CALL CFQ_13
      VEQ=VEQ+CC*DMR_13(4,R(J))
   10 CONTINUE
      RETURN
      END

      SUBROUTINE DCEQ_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MAIN2_13/C5Q,RK
      COMMON/START2_13/C
      COMMON/PART2_13/VEQ,DVEQ,D2VEQ
      COMMON/CECF_13/CC,DCC,D2CC
      COMMON/CECF2_13/DCCM
      COMMON/MAIN_13/C2C(3),A0(3),A1(3),DD
      COMMON/START_13/A,B,AK0,AK1,R0,CXX
      COMMON/CONTROL_13/K1D
      COMMON/PART_13/VDP,DVDP,D2VDP
      COMMON/DUMMY_13/INDEX(3,2)
      COMMON/DMY_13/IX(3,3)
      COMMON/COORD_13/R(3)
      DVEQ=0.0
      J=K1D
      I1=INDEX(J,1)
      I2=INDEX(J,2)
      A=R(I1)
      B=R(I2)
      C=R(J)
      CALL CFQ_13
      DVEQ=DVEQ+CC*DDMR_13(4,R(J))
      DO 10 K=1,2
      J=INDEX(K1D,K)
      IO=IX(K1D,J)
      A=R(K1D)
      B=R(IO)
      C=R(J)
      CALL DCFQ_13
      DVEQ=DVEQ+DCC*DMR_13(4,R(J))
   10 CONTINUE
      RETURN
      END

      SUBROUTINE CFQ_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MAIN2_13/C5Q,RK
      COMMON/START2_13/C
      COMMON/PART2_13/VEQ,DVEQ,D2VEQ
      COMMON/START_13/A,B,AK0,AK1,R0,CXX
      COMMON/CECF_13/CC,DCC,D2CC
      COMMON/CECF2_13/DCCM
      CCA=A**4/R0**4
      CCB=B**4/R0**4
      TANHA=TANH(RK*A)
      TANHB=TANH(RK*B)
      EA=EXP(-RK*(A-R0))
      EB=EXP(-RK*(B-R0))
      CC=0.5*C5Q*(CCA*EA*TANHB+CCB*EB*TANHA)
      RETURN
      END

      SUBROUTINE DCFQ_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/START2_13/C
      COMMON/MAIN2_13/C5Q,RK
      COMMON/START_13/A,B,AK0,AK1,R0,CXX
      COMMON/CECF_13/CC,DCC,D2CC
      COMMON/CECF2_13/DCCM
      RKA=RK*A
      CCA=A**4/R0**4
      DCCA=4.*A**3/R0**4
      CCB=B**4/R0**4
      EA=EXP(-RK*(A-R0))
      EB=EXP(-RK*(B-R0))
      TANHA=TANH(RKA)
      TANHB=TANH(RK*B)
      if(rka.lt.20.0)then
      SINHA=SINH(RKA)
      COSHA=SINHA/TANHA
      CSHA2=1./COSHA*1./COSHA
      else
      csha2=0.0
      endif
      DCC=0.5*C5Q*(-CCA*EA*TANHB*RK+DCCA*EA*TANHB+RK*CCB*EB*CSHA2)
      RETURN
      END

!      SUBROUTINE DCFQM_13
!     IMPLICIT REAL*8 (A-H,O-Z)
!      COMMON/MAIN2_13/C5Q,RK
!      COMMON/START2_13/C
!      COMMON/PART2_13/VEQ,DVEQ,D2VEQ
!      COMMON/START_13/A,B,AK0,AK1,R0,CXX
!      COMMON/CECF_13/CC,DCC,D2CC
!      COMMON/CECF2_13/DCCM
!      RRA=(B-C)/A
!      RRB=(A-C)/B
!      DCCA=-6.0*RRA*A**3/R0**4
!      DCCB=-6.0*RRB*B**3/R0**4
!      TANHA=TANH(RK*A)
!      TANHB=TANH(RK*B)
!      EA=EXP(-RK*(A-R0))
!      EB=EXP(-RK*(B-R0))
!      DCCM=0.5*C5Q*(DCCA*EA*TANHB+DCCB*EB*TANHA)
!      RETURN
!      END

      FUNCTION VAD_13(R1,R2,R3)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/DM8_13/DRDQ(3,3),DQDR(3,3)
      COMMON/VAR_13/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
      COMMON/MAIN_13/V(3),AK(3),BK(3),DD
      COMMON/PASS_13/Q1,Q2,Q3,TH,DTH
      S1=R1-DD
      S2=R2-DD
      S3=R3-DD
      Q1=S1*DQDR(1,1)+S2*DQDR(1,2)+S3*DQDR(1,3)
      Q2=S1*DQDR(2,1)+S2*DQDR(2,2)+S3*DQDR(2,3)
      Q3=S1*DQDR(3,1)+S2*DQDR(3,2)+S3*DQDR(3,3)
      PART=SUMAT_13(R1,R2,R3)
      TH=1.-TANH(GAMA*Q1/4.0)
      VAD_13=(B1+B2*Q1+(AL(1)+AL(2)*(Q2**2+Q3**2)+AL(3)* &
        (Q2**2+Q3**2)**2+AL4*(Q2**2+Q3**2))*(1.5+PART)**5)* &
       (1.2527+PART)**5*EXP(-BETAD*(Q2**2+Q3**2))*TH
      RETURN
      END

      FUNCTION SUMAT_13(R1,R2,R3)
      IMPLICIT REAL*8 (A-H,O-Z)
      SUMAT_13=(R1**2-R2**2-R3**2)/(2.*R2*R3)+ &
            (R2**2-R3**2-R1**2)/(2.*R1*R3)+ &
            (R3**2-R1**2-R2**2)/(2.*R1*R2)
      RETURN
      END

      FUNCTION PTH_13(R1,R2,R3)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/VAR_13/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
      COMMON/CONST2_13/A,C3,C6,C10
      COMMON/PASS_13/Q1,Q2,Q3,TH,DTH

      TH2=1.-TANH(GAMA*Q1/2.)
      P=A+C1*Q1+C2*Q1**2+C3*(Q2**2+Q3**2)+C5*Q1*(Q2**2+Q3**2) &
        +C6*(Q3**3-3.*Q2**2*Q3)+C10*(Q2**2+Q3**2)**2
      PTH_13=P*TH2
      RETURN
      END

      FUNCTION DVAD_13(R1,R2,R3,I)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DM8_13/DRDQ(3,3),DQDR(3,3)
      COMMON/VAR_13/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
      COMMON/MAIN_13/V(3),AK(3),BK(3),DD
      COMMON/PASS_13/Q1,Q2,Q3,TH,DTH
      S1=R1-DD
      S2=R2-DD
      S3=R3-DD
      Q1=S1*DQDR(1,1)+S2*DQDR(1,2)+S3*DQDR(1,3)
      Q2=S1*DQDR(2,1)+S2*DQDR(2,2)+S3*DQDR(2,3)
      Q3=S1*DQDR(3,1)+S2*DQDR(3,2)+S3*DQDR(3,3)
      TH=1.-TANH(GAMA*Q1/4.0)
      DTH=-GAMA/(4.0*COSH(GAMA*Q1/4.0)**2)*DQDR(1,I)
      PART=SUMAT_13(R1,R2,R3)
      DPART=DSUMAT_13(R1,R2,R3,I)
      T1=B1+B2*Q1
      T2=(AL(1)+AL(2)*(Q2**2+Q3**2)+AL(3)*(Q2**2+Q3**2)**2)
      T3=(1.5+PART)**5
      T4=(1.2527+PART)**5
      T5=EXP(-BETAD*(Q2**2+Q3**2))
      T6=TH
      DT1=B2*DQDR(1,I)
      DT2=2.*AL(2)*(Q2*DQDR(2,I)+Q3*DQDR(3,I)) &
               +4.*AL(3)*((Q2**3+Q2*Q3**2)*DQDR(2,I)+(Q3**3+Q2**2*Q3) &
               *DQDR(3,I))
      DT3=5.*(1.5000+PART)**4*DPART
      DT4=5.*(1.2527+PART)**4*DPART
      DT5=-2.*BETAD*(Q2*DQDR(2,I)+Q3*DQDR(3,I))*T5
      DT6=DTH
      DVAD_13=(DT1+DT2*T3+T2*DT3)*T4*T5*T6 &
          +(T1+T2*T3)*DT4*T5*T6 &
          +(T1+T2*T3)*T4*DT5*T6 &
          +(T1+T2*T3)*T4*T5*DT6
      RETURN
      END

      FUNCTION DSUMAT_13( R1,R2,R3,I)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DUMMY_13/INDEX(3,2)
      DIMENSION R(3)
      R(1)=R1
      R(2)=R2
      R(3)=R3
      I1=INDEX(I,1)
      I2=INDEX(I,2)
      DSUMAT_13=R(I)/(R(I1)*R(I2)) &
          +(R(I2)**3-R(I1)**2*R(I2)-R(I)**2*R(I2))/(2.*R(I)**2*R(I2)**2) &
          +(R(I1)**3-R(I2)**2*R(I1)-R(I)**2*R(I1))/(2.*R(I)**2*R(I1)**2)
      RETURN
      END

      FUNCTION DPTH_13(R1,R2,R3,I)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/DM8_13/DRDQ(3,3),DQDR(3,3)
      COMMON/VAR_13/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
      COMMON/CONST2_13/A,C3,C6,C10
      COMMON/PASS_13/Q1,Q2,Q3,TH,DTH

      DTH2=-GAMA/(2.*COSH(GAMA*Q1/2.)**2)*DQDR(1,I)
      TH2=1.-TANH(GAMA*Q1/2.)
      P=A+C1*Q1+C2*Q1**2+C3*(Q2**2+Q3**2)+C5*Q1*(Q2**2+Q3**2)+C6*(Q3**3- &
        3.*Q2**2*Q3)+C10*(Q2**2+Q3**2)**2
      DPVTH=DTH2*P
      DP=C1*DQDR(1,I)+2.*C2*Q1*DQDR(1,I)+2.*C3*(Q2*DQDR(2,I) &
         +Q3*DQDR(3,I))+C5*DQDR(1,I)*(Q2**2+Q3**2)+C5*Q1*(2.*Q2 &
         *DQDR(2,I)+2.*Q3*DQDR(3,I))+C6*(3.*Q3**2*DQDR(3,I)-6.*Q2*Q3 &
         *DQDR(2,I)-3.*Q2**2*DQDR(3,I))+C10*((4.*Q2**3+4.*Q2*Q3**2) &
         *DQDR(2,I)+(4.*Q3**3+4.*Q2**2*Q3)*DQDR(3,I))
      DPTH_13=DPVTH+DP*TH2
      RETURN
      END

      BLOCK DATA AAA1_13
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/MAIN_13/VW(3),AK0(3),AK1(3),DD
      COMMON/MAIN2_13/C5Q,RK
      COMMON/PN_13/N(4)
      COMMON/DM1_13/D1(4),D2(4)
      COMMON/DUMMY_13/INDEX(3,2)
      COMMON/DMY_13/IX(3,3)
      COMMON/DM2_13/R0,RM
      COMMON/DM3_13/GAMA
      COMMON/DM4_13/A1,A2,A3,GM,RE,DE
      COMMON/DM8_13/DRDQ(3,3),DQDR(3,3)
      COMMON/VAR_13/C1,C2,C5,GAMAP,B1,B2,AL(3),AL4,BETAD
      COMMON/CONST2_13/A,C3,C6,C10
      COMMON/START_13/AS,BS,ASK0,ASK1,R0PASS,CSXX
      COMMON/TA2_13/FPE2,GAMAT2,SUMT2
      COMMON/TA_13/FPE,GAMATA,SUMTA
      DATA FPE,GAMATA,SUMTA/2.9E-2,3.1410047d-2,22.323365/
      DATA FPE2,GAMAT2,SUMT2/1.4301464d-3,3.1410047d-2,35.549516/
      DATA R0PASS/2.2818/
      DATA C1,C2,C5,GAMAP,B1,B2,(AL(I),I=1,3),AL4,BETAD/.1157783, &
        .1682194, &
       -.6563872e-02,2.54,.1561575E+03,.4951110E+02,.11574778E7, &
       -.35624264E6,+.27725091E5,0.0,.224/
      DATA A,C3,C6,C10/.3058714,-.1502505,.1805864E-1,.3434551E-1/
      DATA(VW(I),I=1,3)/15.4,235.21994,4066.2393/
      DATA(AK0(I),I=1,3)/.2772806,.6233148,1.237123/
      DATA(AK1(I),I=1,3)/.9545174,.8116684,.7043832/
      DATA DD/2.967200/
      DATA RK,C5Q/3.3522,+1.3144/
      DATA(N(I),I=1,4)/6,8,10,5/
      DATA ((IX(I,J),I=1,3),J=1,3)/0,3,2,3,0,1,2,1,0/
      DATA ((INDEX(I,J),I=1,3),J=1,2)/2,1,1,3,3,2/
      DATA (D1(I),I=1,4)/3.0951333,2.1999,1.6880714,3.8428294/
      DATA (D2(I),I=1,4)/2.8363203,3.2849276,3.5239687,2.5178884/
      DATA GAMA/2.5/
      DATA R0,RM/5.661693,2.2818/
      DATA A1,A2,A3,GM,RE,DE/3.6445902,3.9281208,2.0986665,3.3522487, &
                             2.2818,-0.142912/
      DATA ((DRDQ(I,J),I=1,3),J=1,3)/.5773502,.5773502,.5773502, &
                                      0.0,.7071067,-.7071067, &
                                      .8164965,-.4082482,-.4082482/
      DATA ((DQDR(I,J),I=1,3),J=1,3)/.5773502,0.0,.8164965, &
                                    .5773502,.7071067,-.4082482, &
                                    .5773502,-.7071067,-.4082482/
      END




      FUNCTION ELECHO3_13(RAB,RAC,RAD,RBC,RBD,RCD)
      IMPLICIT REAL*8(A-H,O-Z)
      C4HOO=0.0D0
      C4OHO=-0.92921D0
      C5HOO=0.0D0
      C5OHO=-1.790D0
      C4OOO=0.0D0
      C5OOO=-1.3144D0
!        C4HO=-0.5904D0
!      C5HO=-0.30825D0
      C4HO=0.0D0
      C5HO=0.0D0
      CALL DAGSHOH_13(RAB,GAB,SHAB)
      CALL DAGSHOH_13(RAC,GAC,SHAC)
      CALL DAGSHOH_13(RAD,GAD,SHAD)
      CALL DAGSHOO_13(RBC,GBC,SHBC)
      CALL DAGSHOO_13(RBD,GBD,SHBD)
      CALL DAGSHOO_13(RCD,GCD,SHCD)
      e=1.0d0
        HAB=DAHOH_13(RAB)
        HAC=DAHOH_13(RAC)
        HAD=DAHOH_13(RAD)
        HBC=DAHOO_13(RBC)
        HBD=DAHOO_13(RBD)
        HCD=DAHOO_13(RCD)
!  THE DIPOLE-QUADRAPOLE INTERACTION
          XRAB=DMROH_13(4,RAB)
      XRAC=DMROH_13(4,RAC)
      XRAD=DMROH_13(4,RAD)
      XRBC=DMROO_13(4,RBC)
      XRBD=DMROO_13(4,RBD)
      XRCD=DMROO_13(4,RCD)
          SGAB=XIAOGOH_13(4,RAB)
      SGAC=XIAOGOH_13(4,RAC)
      SGAD=XIAOGOH_13(4,RAD)
      SGBC=XIAOGOO_13(4,RBC)
      SGBD=XIAOGOO_13(4,RBD)
      SGCD=XIAOGOO_13(4,RCD)
          CHOO=C4HOO
      COHO=C4OHO
      COOO=C4OOO
        CABCD=C4HO*BOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
        CACBD=C4HO*BOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD)
        CADBC=C4HO*BOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD)
      XRAB=0.5D0*XRAB*(COHO*(GAC*SHBC*(HBD*SGAD*SGCD-e)+ &
                GAD*SHBD*(HBC*SGAC*SGCD-e))+ &
                0.5D0*(CADBC*GAD*GBC*SHAC*SHBD*SHCD+ &
                CACBD*GAC*GBD*SHAD*SHBC*SHCD))
           XRAC=0.5D0*XRAC*(COHO*(GAB*SHBC*(HCD*SGAD*SGBD-e)+ &
                GAD*SHCD*(HBC*SGAB*SGBD-e))+ &
                0.5D0*(CABCD*GAB*GCD*SHAD*SHBC*SHBD+ &
                CADBC*GAD*GBC*SHAB*SHCD*SHBD))
           XRAD=0.5D0*XRAD*(COHO*(GAB*SHBD*(HCD*SGAC*SGBC-e)+ &
                GAC*SHCD*(HBD*SGAB*SGBC-e))+ &
                0.5D0*(CABCD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CACBD*GAC*GBD*SHAB*SHCD*SHBC))
        XRBC=0.5D0*XRBC*(COHO*(GAC*SHAB*(HBD*SGAD*SGCD-e)+ &
                GAB*SHAC*(HCD*SGAD*SGBD-e))+ &
                COOO*(GCD*SHBD*(HAB*SGAC*SGAD-e)+ &
                GBD*SHCD*(HAC*SGAB*SGAD-e))+0.5D0*( &
                CABCD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CACBD*GAC*GBD*SHAB*SHCD*SHAD))
        XRBD=0.5D0*XRBD*(COHO*(GAD*SHAB*(HBC*SGAC*SGCD-e)+ &
                GAB*SHAD*(HCD*SGAC*SGBC-e))+ &
                COOO*(GCD*SHBC*(HAB*SGAC*SGAD-e)+ &
                GBC*SHCD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CABCD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CADBC*GAD*GBC*SHAB*SHCD*SHAC))
        XRCD=0.5D0*XRCD*(COHO*(GAD*SHAC*(HBC*SGAB*SGBD-e)+ &
                GAC*SHAD*(HBD*SGAB*SGBC-e))+ &
                COOO*(GBD*SHBC*(HAC*SGAB*SGAD-e)+ &
                GBC*SHBD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CACBD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CADBC*GAD*GBC*SHAC*SHBD*SHAB))
            TERM4=XRAB+XRAC+XRAD+XRBC+XRBD+XRCD
!  THE QUADRAPOLE-QUADRAPOLE INTERACTION
            XRAB=DMROH_13(5,RAB)
      XRAC=DMROH_13(5,RAC)
      XRAD=DMROH_13(5,RAD)
      XRBC=DMROO_13(5,RBC)
      XRBD=DMROO_13(5,RBD)
      XRCD=DMROO_13(5,RCD)
          SGAB=XIAOGOH_13(5,RAB)
      SGAC=XIAOGOH_13(5,RAC)
      SGAD=XIAOGOH_13(5,RAD)
      SGBC=XIAOGOO_13(5,RBC)
      SGBD=XIAOGOO_13(5,RBD)
      SGCD=XIAOGOO_13(5,RCD)
          CHOO=C5HOO
      COHO=C5OHO
      COOO=C5OOO
        CABCD=C5HO*AOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
        CACBD=C5HO*AOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD)
        CADBC=C5HO*AOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD)
      XRAB=0.5D0*XRAB*(COHO*(GAC*SHBC*(HBD*SGAD*SGCD-e)+ &
                GAD*SHBD*(HBC*SGAC*SGCD-e))+ &
                0.5D0*(CADBC*GAD*GBC*SHAC*SHBD*SHCD+ &
                CACBD*GAC*GBD*SHAD*SHBC*SHCD))
           XRAC=0.5D0*XRAC*(COHO*(GAB*SHBC*(HCD*SGAD*SGBD-e)+ &
                GAD*SHCD*(HBC*SGAB*SGBD-e))+ &
                0.5D0*(CABCD*GAB*GCD*SHAD*SHBC*SHBD+ &
                CADBC*GAD*GBC*SHAB*SHCD*SHBD))
           XRAD=0.5D0*XRAD*(COHO*(GAB*SHBD*(HCD*SGAC*SGBC-e)+ &
                GAC*SHCD*(HBD*SGAB*SGBC-e))+ &
                0.5D0*(CABCD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CACBD*GAC*GBD*SHAB*SHCD*SHBC))
        XRBC=0.5D0*XRBC*(COHO*(GAC*SHAB*(HBD*SGAD*SGCD-e)+ &
                GAB*SHAC*(HCD*SGAD*SGBD-e))+ &
                COOO*(GCD*SHBD*(HAB*SGAC*SGAD-e)+ &
                GBD*SHCD*(HAC*SGAB*SGAD-e))+0.5D0*( &
                CABCD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CACBD*GAC*GBD*SHAB*SHCD*SHAD))
        XRBD=0.5D0*XRBD*(COHO*(GAD*SHAB*(HBC*SGAC*SGCD-e)+ &
                GAB*SHAD*(HCD*SGAC*SGBC-e))+ &
                COOO*(GCD*SHBC*(HAB*SGAC*SGAD-e)+ &
                GBC*SHCD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CABCD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CADBC*GAD*GBC*SHAB*SHCD*SHAC))
        XRCD=0.5D0*XRCD*(COHO*(GAD*SHAC*(HBC*SGAB*SGBD-e)+ &
                GAC*SHAD*(HBD*SGAB*SGBC-e))+ &
                COOO*(GBD*SHBC*(HAC*SGAB*SGAD-e)+ &
                GBC*SHBD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CACBD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CADBC*GAD*GBC*SHAC*SHBD*SHAB))
            TERM5=XRAB+XRAC+XRAD+XRBC+XRBD+XRCD
      TERMA=ELECHO3P_13(RAB,RAC,RAD,RBC,RBD,RCD)
      ELECHO3_13=TERM4+TERM5+TERMA
      RETURN
      END

!************************************************
        FUNCTION DMRDD_13(I,R)
        IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION D1(5),D2(5),N(5)
      D1(4)=5.0079875D0
      D1(5)=3.8428294D0
      D2(4)=2.129498D0
      D2(5)=2.5178884D0
      N(4)=4
      N(5)=5
      RO=15.15188

      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      RL=R**IDXL
      DMRDD_13=C/RL
      RETURN
      END

      FUNCTION DDMRDD_13(I,R)
      IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION D1(5),D2(5),N(5)
      D1(4)=5.0079875D0
      D1(5)=3.8428294D0
      D2(4)=2.129498D0
      D2(5)=2.5178884D0
      N(4)=4
      N(5)=5
      RO=15.15188

      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      DPOL=-D1(I)*(1.+2.*D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      EXPOEL=C/(1.-EXPO)
      DC=DXL*EXPOEL*(-EXPO)*DPOL*DX
      RL=R**IDXL
      DDMRDD_13=(DC/RL-DXL*C/(RL*R))
      RETURN
      END

      FUNCTION ELECHO3P_13(RAB,RAC,RAD,RBC,RBD,RCD)
      IMPLICIT REAL*8(A-H,O-Z)
      C4HOO=0.0D0
      C4OHO=-0.92921D0
      C4OOO=0.0D0
      C5HOO=0.0D0
      C5OHO=-1.790D0
      C5OOO=-1.3144D0
        C4HO=-0.5904D0
      C5HO=-0.30825D0
      CALL DAGSHOH_13(RAB,GAB,SHAB)
      CALL DAGSHOH_13(RAC,GAC,SHAC)
      CALL DAGSHOH_13(RAD,GAD,SHAD)
      CALL DAGSHOO_13(RBC,GBC,SHBC)
      CALL DAGSHOO_13(RBD,GBD,SHBD)
      CALL DAGSHOO_13(RCD,GCD,SHCD)
        HAB=DAHOH_13(RAB)
        HAC=DAHOH_13(RAC)
        HAD=DAHOH_13(RAD)
        HBC=DAHOO_13(RBC)
        HBD=DAHOO_13(RBD)
        HCD=DAHOO_13(RCD)
!---------------------------
        CABCD=C4HO*BOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
      CABCD=CABCD-(C4HOO+C4OOO)*HAB-(C4OHO+C4OHO)*HCD
      RABCD=(RAC+RAD+RBC+RBD)/4.0
        EABCD=CABCD*GAB*GCD*DMRDD_13(4,RABCD)

        CACBD=C4HO*BOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD)
      CACBD=CACBD-(C4HOO+C4OOO)*HAC-(C4OHO+C4OHO)*HBD
      RACBD=(RAB+RAD+RBC+RCD)/4.0
        EACBD=CACBD*GAC*GBD*DMRDD_13(4,RACBD)

        CADBC=C4HO*BOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD)
      CADBC=CADBC-(C4HOO+C4OOO)*HAD-(C4OHO+C4OHO)*HBC
      RADBC=(RAB+RAC+RBD+RCD)/4.0
        EADBC=CADBC*GAD*GBC*DMRDD_13(4,RADBC)
      TERM4=EABCD+EACBD+EADBC
!---------------------------
        CABCD=C5HO*AOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
      CABCD=CABCD-(C5HOO+C5OOO)*HAB-(C5OHO+C5OHO)*HCD
      RABCD=(RAC+RAD+RBC+RBD)/4.0
        EABCD=CABCD*GAB*GCD*DMRDD_13(5,RABCD)

        CACBD=C5HO*AOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD)
      CACBD=CACBD-(C5HOO+C5OOO)*HAC-(C5OHO+C5OHO)*HBD
      RACBD=(RAB+RAD+RBC+RCD)/4.0
        EACBD=CACBD*GAC*GBD*DMRDD_13(5,RACBD)

        CADBC=C5HO*AOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD)
      CADBC=CADBC-(C5HOO+C5OOO)*HAD-(C5OHO+C5OHO)*HBC
      RADBC=(RAB+RAC+RBD+RCD)/4.0
        EADBC=CADBC*GAD*GBC*DMRDD_13(5,RADBC)
      TERM5=EABCD+EACBD+EADBC

!--------------------------------------------------
      ELECHO3P_13=TERM4+TERM5
      return
      END

      SUBROUTINE DELECHO3P_13(RAB,RAC,RAD,RBC,RBD,RCD,DELP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DELP(6),DVB(6),DVA(6),TERM4(6),TERM5(6)
      C4HOO=0.0D0
      C4OHO=-0.92921D0
      C4OOO=0.0D0
      C5HOO=0.0D0
      C5OHO=-1.790D0
      C5OOO=-1.3144D0
        C4HO=-0.5904D0
      C5HO=-0.30825D0
      CALL DAGSHOH_13(RAB,GAB,SHAB)
      CALL DAGSHOH_13(RAC,GAC,SHAC)
      CALL DAGSHOH_13(RAD,GAD,SHAD)
      CALL DAGSHOO_13(RBC,GBC,SHBC)
      CALL DAGSHOO_13(RBD,GBD,SHBD)
      CALL DAGSHOO_13(RCD,GCD,SHCD)
        HAB=DAHOH_13(RAB)
        HAC=DAHOH_13(RAC)
        HAD=DAHOH_13(RAD)
        HBC=DAHOO_13(RBC)
        HBD=DAHOO_13(RBD)
        HCD=DAHOO_13(RCD)
      CALL DDAGSHOH_13(RAB,DGAB,DSHAB)
      CALL DDAGSHOH_13(RAC,DGAC,DSHAC)
      CALL DDAGSHOH_13(RAD,DGAD,DSHAD)
      CALL DDAGSHOO_13(RBC,DGBC,DSHBC)
      CALL DDAGSHOO_13(RBD,DGBD,DSHBD)
      CALL DDAGSHOO_13(RCD,DGCD,DSHCD)
        DHAB=DDAHOH_13(RAB)
        DHAC=DDAHOH_13(RAC)
        DHAD=DDAHOH_13(RAD)
        DHBC=DDAHOO_13(RBC)
        DHBD=DDAHOO_13(RBD)
        DHCD=DDAHOO_13(RCD)
!---------------------------
        CABCD=C4HO*BOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
      CALL DBOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD,DVB)
      CABCD=CABCD-(C4HOO+C4OOO)*HAB-(C4OHO+C4OHO)*HCD
       DCRAB=C4HO*DVB(1)-(C4HOO+C4OOO)*DHAB
      DCRCD=C4HO*DVB(2)-(C4OHO+C4OHO)*DHCD
      DCRAC=C4HO*DVB(3)
      DCRAD=C4HO*DVB(4)
      DCRBC=C4HO*DVB(5)
      DCRBD=C4HO*DVB(6)
        RABCD=(RAC+RAD+RBC+RBD)/4.0
      DMR_13=DMRDD_13(4,RABCD)
      DDMR_13=DDMRDD_13(4,RABCD)/4.0
        TERM4(1)=DCRAB*GAB*GCD*DMR_13+CABCD*DGAB*GCD*DMR_13
        TERM4(2)=DCRAC*GAB*GCD*DMR_13+CABCD*GAB*GCD*DDMR_13
        TERM4(3)=DCRAD*GAB*GCD*DMR_13+CABCD*GAB*GCD*DDMR_13
        TERM4(4)=DCRBC*GAB*GCD*DMR_13+CABCD*GAB*GCD*DDMR_13
        TERM4(5)=DCRBD*GAB*GCD*DMR_13+CABCD*GAB*GCD*DDMR_13
        TERM4(6)=DCRCD*GAB*GCD*DMR_13+CABCD*GAB*DGCD*DMR_13

        CACBD=C4HO*BOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD)
        CALL DBOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD,DVB)
      CACBD=CACBD-(C4HOO+C4OOO)*HAC-(C4OHO+C4OHO)*HBD
      DCRAC=C4HO*DVB(1)-(C4HOO+C4OOO)*DHAC
      DCRBD=C4HO*DVB(2)-(C4OHO+C4OHO)*DHBD
      DCRAB=C4HO*DVB(3)
      DCRAD=C4HO*DVB(4)
      DCRBC=C4HO*DVB(5)
      DCRCD=C4HO*DVB(6)
      RACBD=(RAB+RAD+RBC+RCD)/4.0
        DMR_13=DMRDD_13(4,RACBD)
      DDMR_13=DDMRDD_13(4,RACBD)/4.0
        TERM4(1)=TERM4(1)+DCRAB*GAC*GBD*DMR_13+CACBD*GAC*GBD*DDMR_13
        TERM4(2)=TERM4(2)+DCRAC*GAC*GBD*DMR_13+CACBD*DGAC*GBD*DMR_13
        TERM4(3)=TERM4(3)+DCRAD*GAC*GBD*DMR_13+CACBD*GAC*GBD*DDMR_13
        TERM4(4)=TERM4(4)+DCRBC*GAC*GBD*DMR_13+CACBD*GAC*GBD*DDMR_13
        TERM4(5)=TERM4(5)+DCRBD*GAC*GBD*DMR_13+CACBD*GAC*DGBD*DMR_13
        TERM4(6)=TERM4(6)+DCRCD*GAC*GBD*DMR_13+CACBD*GAC*GBD*DDMR_13

        CADBC=C4HO*BOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD)
      CALL DBOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD,DVB)
      CADBC=CADBC-(C4HOO+C4OOO)*HAD-(C4OHO+C4OHO)*HBC
      DCRAD=C4HO*DVB(1)-(C4HOO+C4OOO)*DHAD
      DCRBC=C4HO*DVB(2)-(C4OHO+C4OHO)*DHBC
      DCRAB=C4HO*DVB(3)
      DCRAC=C4HO*DVB(4)
      DCRBD=C4HO*DVB(5)
      DCRCD=C4HO*DVB(6)
      RADBC=(RAB+RAC+RBD+RCD)/4.0
        DMR_13=DMRDD_13(4,RADBC)
      DDMR_13=DDMRDD_13(4,RADBC)/4.0
        TERM4(1)=TERM4(1)+DCRAB*GAD*GBC*DMR_13+CADBC*GAD*GBC*DDMR_13
        TERM4(2)=TERM4(2)+DCRAC*GAD*GBC*DMR_13+CADBC*GAD*GBC*DDMR_13
        TERM4(3)=TERM4(3)+DCRAD*GAD*GBC*DMR_13+CADBC*DGAD*GBC*DMR_13
        TERM4(4)=TERM4(4)+DCRBC*GAD*GBC*DMR_13+CADBC*GAD*DGBC*DMR_13
        TERM4(5)=TERM4(5)+DCRBD*GAD*GBC*DMR_13+CADBC*GAD*GBC*DDMR_13
        TERM4(6)=TERM4(6)+DCRCD*GAD*GBC*DMR_13+CADBC*GAD*GBC*DDMR_13
!---------------------------
        CABCD=C5HO*AOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
      CALL DAOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD,DVA)
      CABCD=CABCD-(C5HOO+C5OOO)*HAB-(C5OHO+C5OHO)*HCD
      DCRAB=C5HO*DVA(1)-(C5HOO+C5OOO)*DHAB
      DCRCD=C5HO*DVA(2)-(C5OHO+C5OHO)*DHCD
      DCRAC=C5HO*DVA(3)
      DCRAD=C5HO*DVA(4)
      DCRBC=C5HO*DVA(5)
      DCRBD=C5HO*DVA(6)
      RABCD=(RAC+RAD+RBC+RBD)/4.0
        DMR_13=DMRDD_13(5,RABCD)
      DDMR_13=DDMRDD_13(5,RABCD)/4.0
        TERM5(1)=DCRAB*GAB*GCD*DMR_13+CABCD*DGAB*GCD*DMR_13
        TERM5(2)=DCRAC*GAB*GCD*DMR_13+CABCD*GAB*GCD*DDMR_13
        TERM5(3)=DCRAD*GAB*GCD*DMR_13+CABCD*GAB*GCD*DDMR_13
        TERM5(4)=DCRBC*GAB*GCD*DMR_13+CABCD*GAB*GCD*DDMR_13
        TERM5(5)=DCRBD*GAB*GCD*DMR_13+CABCD*GAB*GCD*DDMR_13
        TERM5(6)=DCRCD*GAB*GCD*DMR_13+CABCD*GAB*DGCD*DMR_13

        CACBD=C5HO*AOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD)
      CALL DAOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD,DVA)
      CACBD=CACBD-(C5HOO+C5OOO)*HAC-(C5OHO+C5OHO)*HBD
      DCRAC=C5HO*DVA(1)-(C5HOO+C5OOO)*DHAC
      DCRBD=C5HO*DVA(2)-(C5OHO+C5OHO)*DHBD
      DCRAB=C5HO*DVA(3)
      DCRAD=C5HO*DVA(4)
      DCRBC=C5HO*DVA(5)
      DCRCD=C5HO*DVA(6)
      RACBD=(RAB+RAD+RBC+RCD)/4.0
      DMR_13=DMRDD_13(5,RACBD)
      DDMR_13=DDMRDD_13(5,RACBD)/4.0
        TERM5(1)=TERM5(1)+DCRAB*GAC*GBD*DMR_13+CACBD*GAC*GBD*DDMR_13
        TERM5(2)=TERM5(2)+DCRAC*GAC*GBD*DMR_13+CACBD*DGAC*GBD*DMR_13
        TERM5(3)=TERM5(3)+DCRAD*GAC*GBD*DMR_13+CACBD*GAC*GBD*DDMR_13
        TERM5(4)=TERM5(4)+DCRBC*GAC*GBD*DMR_13+CACBD*GAC*GBD*DDMR_13
        TERM5(5)=TERM5(5)+DCRBD*GAC*GBD*DMR_13+CACBD*GAC*DGBD*DMR_13
        TERM5(6)=TERM5(6)+DCRCD*GAC*GBD*DMR_13+CACBD*GAC*GBD*DDMR_13

        CADBC=C5HO*AOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD)
      CALL DAOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD,DVA)
      CADBC=CADBC-(C5HOO+C5OOO)*HAD-(C5OHO+C5OHO)*HBC
      DCRAD=C5HO*DVA(1)-(C5HOO+C5OOO)*DHAD
      DCRBC=C5HO*DVA(2)-(C5OHO+C5OHO)*DHBC
      DCRAB=C5HO*DVA(3)
      DCRAC=C5HO*DVA(4)
      DCRBD=C5HO*DVA(5)
      DCRCD=C5HO*DVA(6)
      RADBC=(RAB+RAC+RBD+RCD)/4.0
        DMR_13=DMRDD_13(5,RADBC)
      DDMR_13=DDMRDD_13(5,RADBC)/4.0
        TERM5(1)=TERM5(1)+DCRAB*GAD*GBC*DMR_13+CADBC*GAD*GBC*DDMR_13
        TERM5(2)=TERM5(2)+DCRAC*GAD*GBC*DMR_13+CADBC*GAD*GBC*DDMR_13
        TERM5(3)=TERM5(3)+DCRAD*GAD*GBC*DMR_13+CADBC*DGAD*GBC*DMR_13
        TERM5(4)=TERM5(4)+DCRBC*GAD*GBC*DMR_13+CADBC*GAD*DGBC*DMR_13
        TERM5(5)=TERM5(5)+DCRBD*GAD*GBC*DMR_13+CADBC*GAD*GBC*DDMR_13
        TERM5(6)=TERM5(6)+DCRCD*GAD*GBC*DMR_13+CADBC*GAD*GBC*DDMR_13
!--------------------------------------------------
      DO 10 I=1,6
      DELP(I)=TERM4(I)+TERM5(I)
 10      CONTINUE
      return
      END

!********************************************



      SUBROUTINE DELECHO3_13(RAB,RAC,RAD,RBC,RBD,RCD,DVEL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DVEL(6),DEL(6),TERM4(6),TERM5(6),PEL(6), &
                PXRAB(6),PXRAC(6),PXRAD(6),PXRBC(6),PXRBD(6),PXRCD(6)
      C4HOO=0.0D0
      C4OHO=-0.92921D0
      C5HOO=0.0D0
      C5OHO=-1.790D0
      C4OOO=0.0D0
      C5OOO=-1.3144D0
!        C4HO=-0.5904D0
!      C5HO=-0.30825D0
      C4HO=0.0D0
      C5HO=0.0D0
      CALL DAGSHOH_13(RAB,GAB,SHAB)
      CALL DAGSHOH_13(RAC,GAC,SHAC)
      CALL DAGSHOH_13(RAD,GAD,SHAD)
      CALL DAGSHOO_13(RBC,GBC,SHBC)
      CALL DAGSHOO_13(RBD,GBD,SHBD)
      CALL DAGSHOO_13(RCD,GCD,SHCD)

      CALL DDAGSHOH_13(RAB,DGAB,DSHAB)
      CALL DDAGSHOH_13(RAC,DGAC,DSHAC)
      CALL DDAGSHOH_13(RAD,DGAD,DSHAD)
      CALL DDAGSHOO_13(RBC,DGBC,DSHBC)
      CALL DDAGSHOO_13(RBD,DGBD,DSHBD)
      CALL DDAGSHOO_13(RCD,DGCD,DSHCD)

        e=1.0d0
        HAB=DAHOH_13(RAB)
        HAC=DAHOH_13(RAC)
        HAD=DAHOH_13(RAD)
        HBC=DAHOO_13(RBC)
        HBD=DAHOO_13(RBD)
        HCD=DAHOO_13(RCD)

      DHAB=DDAHOH_13(RAB)
        DHAC=DDAHOH_13(RAC)
        DHAD=DDAHOH_13(RAD)
        DHBC=DDAHOO_13(RBC)
        DHBD=DDAHOO_13(RBD)
        DHCD=DDAHOO_13(RCD)

!  THE DIPOLE-QUADRAPOLE INTERACTION
          XRAB=DMROH_13(4,RAB)
      XRAC=DMROH_13(4,RAC)
      XRAD=DMROH_13(4,RAD)
      XRBC=DMROO_13(4,RBC)
      XRBD=DMROO_13(4,RBD)
      XRCD=DMROO_13(4,RCD)

          DXRAB=DDMROH_13(4,RAB)
      DXRAC=DDMROH_13(4,RAC)
      DXRAD=DDMROH_13(4,RAD)
      DXRBC=DDMROO_13(4,RBC)
      DXRBD=DDMROO_13(4,RBD)
      DXRCD=DDMROO_13(4,RCD)

          SGAB=XIAOGOH_13(4,RAB)
      SGAC=XIAOGOH_13(4,RAC)
      SGAD=XIAOGOH_13(4,RAD)
      SGBC=XIAOGOO_13(4,RBC)
      SGBD=XIAOGOO_13(4,RBD)
      SGCD=XIAOGOO_13(4,RCD)

          DSGAB=DXIAOGOH_13(4,RAB)
      DSGAC=DXIAOGOH_13(4,RAC)
      DSGAD=DXIAOGOH_13(4,RAD)
      DSGBC=DXIAOGOO_13(4,RBC)
      DSGBD=DXIAOGOO_13(4,RBD)
      DSGCD=DXIAOGOO_13(4,RCD)

          CHOO=C4HOO
      COHO=C4OHO
      COOO=C4OOO
      CFHO=C4HO
        CABCD=C4HO*BOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
        CALL DBOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD,DEL)
      CABRAB=DEL(1)*CFHO
      CABRCD=DEL(2)*CFHO
      CABRAC=DEL(3)*CFHO
      CABRAD=DEL(4)*CFHO
      CABRBC=DEL(5)*CFHO
      CABRBD=DEL(6)*CFHO

      CACBD=C4HO*BOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD)
        CALL DBOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD,DEL)
      CACRAC=DEL(1)*CFHO
      CACRBD=DEL(2)*CFHO
      CACRAB=DEL(3)*CFHO
      CACRAD=DEL(4)*CFHO
      CACRBC=DEL(5)*CFHO
      CACRCD=DEL(6)*CFHO

        CADBC=C4HO*BOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD)
        CALL DBOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD,DEL)
      CADRAD=DEL(1)*CFHO
      CADRBC=DEL(2)*CFHO
      CADRAB=DEL(3)*CFHO
      CADRAC=DEL(4)*CFHO
      CADRBD=DEL(5)*CFHO
      CADRCD=DEL(6)*CFHO

         PRAB=COHO*(GAC*SHBC*(HBD*SGAD*SGCD-e)+ &
                GAD*SHBD*(HBC*SGAC*SGCD-e))+ &
                0.5D0*(CADBC*GAD*GBC*SHAC*SHBD*SHCD+ &
                CACBD*GAC*GBD*SHAD*SHBC*SHCD)
      PXRAB(1)=0.5D0*(DXRAB*PRAB+XRAB* &
                0.5D0*(CADRAB*GAD*GBC*SHAC*SHBD*SHCD+ &
                CACRAB*GAC*GBD*SHAD*SHBC*SHCD))
      PXRAB(2)=0.5D0*XRAB*(COHO*(DGAC*SHBC*(HBD*SGAD*SGCD-e)+ &
                GAD*SHBD*HBC*DSGAC*SGCD)+ &
                0.5D0*(CADRAC*GAD*GBC*SHAC*SHBD*SHCD &
                +CADBC*GAD*GBC*DSHAC*SHBD*SHCD+ &
                CACRAC*GAC*GBD*SHAD*SHBC*SHCD+ &
                CACBD*DGAC*GBD*SHAD*SHBC*SHCD))
       PXRAB(3)=0.5D0*XRAB*(COHO*(GAC*SHBC*HBD*DSGAD*SGCD+ &
                DGAD*SHBD*(HBC*SGAC*SGCD-e))+ &
                0.5D0*(CADRAD*GAD*GBC*SHAC*SHBD*SHCD+ &
                CADBC*DGAD*GBC*SHAC*SHBD*SHCD+ &
                CACRAD*GAC*GBD*SHAD*SHBC*SHCD+ &
            CACBD*GAC*GBD*DSHAD*SHBC*SHCD))
      PXRAB(4)=0.5D0*XRAB*(COHO*(GAC*DSHBC*(HBD*SGAD*SGCD-e)+ &
                GAD*SHBD*DHBC*SGAC*SGCD)+ &
                0.5D0*(CADRBC*GAD*GBC*SHAC*SHBD*SHCD+ &
                CADBC*GAD*DGBC*SHAC*SHBD*SHCD+ &
                CACRBC*GAC*GBD*SHAD*SHBC*SHCD+ &
            CACBD*GAC*GBD*SHAD*DSHBC*SHCD))
          PXRAB(5)=0.5D0*XRAB*(COHO*(GAC*SHBC*DHBD*SGAD*SGCD+ &
                GAD*DSHBD*(HBC*SGAC*SGCD-e))+ &
                0.5D0*(CADRBD*GAD*GBC*SHAC*SHBD*SHCD+ &
                CADBC*GAD*GBC*SHAC*DSHBD*SHCD+ &
                CACRBD*GAC*GBD*SHAD*SHBC*SHCD+ &
                CACBD*GAC*DGBD*SHAD*SHBC*SHCD))
             PXRAB(6)=0.5D0*XRAB*(COHO*(GAC*SHBC*HBD*SGAD*DSGCD+ &
                GAD*SHBD*HBC*SGAC*DSGCD)+ &
                0.5D0*(CADRCD*GAD*GBC*SHAC*SHBD*SHCD+ &
                CADBC*GAD*GBC*SHAC*SHBD*DSHCD+ &
                CACRCD*GAC*GBD*SHAD*SHBC*SHCD+ &
            CACBD*GAC*GBD*SHAD*SHBC*DSHCD))


           PRAC=COHO*(GAB*SHBC*(HCD*SGAD*SGBD-e)+ &
                GAD*SHCD*(HBC*SGAB*SGBD-e))+ &
                0.5D0*(CABCD*GAB*GCD*SHAD*SHBC*SHBD+ &
                CADBC*GAD*GBC*SHAB*SHCD*SHBD)
           PXRAC(1)=0.5D0*XRAC*(COHO*(DGAB*SHBC*(HCD*SGAD*SGBD-e)+ &
                GAD*SHCD*HBC*DSGAB*SGBD)+ &
                0.5D0*(CABRAB*GAB*GCD*SHAD*SHBC*SHBD+ &
                CABCD*DGAB*GCD*SHAD*SHBC*SHBD+ &
                CADRAB*GAD*GBC*SHAB*SHCD*SHBD+ &
                CADBC*GAD*GBC*DSHAB*SHCD*SHBD))
           PXRAC(2)=0.5D0*(DXRAC*PRAC+XRAC* &
                0.5D0*(CABRAC*GAB*GCD*SHAD*SHBC*SHBD+ &
                CADRAC*GAD*GBC*SHAB*SHCD*SHBD))
           PXRAC(3)=0.5D0*XRAC*(COHO*(GAB*SHBC*HCD*DSGAD*SGBD+ &
                DGAD*SHCD*(HBC*SGAB*SGBD-e))+ &
                0.5D0*(CABRAD*GAB*GCD*SHAD*SHBC*SHBD+ &
                CABCD*GAB*GCD*DSHAD*SHBC*SHBD+ &
                CADRAD*GAD*GBC*SHAB*SHCD*SHBD+ &
                CADBC*DGAD*GBC*SHAB*SHCD*SHBD ))
           PXRAC(4)=0.5D0*XRAC*(COHO*(GAB*DSHBC*(HCD*SGAD*SGBD-e)+ &
                GAD*SHCD*DHBC*SGAB*SGBD)+ &
                0.5D0*(CABRBC*GAB*GCD*SHAD*SHBC*SHBD+ &
            CABCD*GAB*GCD*SHAD*DSHBC*SHBD+ &
                CADRBC*GAD*GBC*SHAB*SHCD*SHBD+ &
            CADBC*GAD*DGBC*SHAB*SHCD*SHBD))
           PXRAC(5)=0.5D0*XRAC*(COHO*(GAB*SHBC*HCD*SGAD*DSGBD+ &
                GAD*SHCD*HBC*SGAB*DSGBD)+ &
                0.5D0*(CABRBD*GAB*GCD*SHAD*SHBC*SHBD+ &
            CABCD*GAB*GCD*SHAD*SHBC*DSHBD+ &
                CADRBD*GAD*GBC*SHAB*SHCD*SHBD+ &
                CADBC*GAD*GBC*SHAB*SHCD*DSHBD))
           PXRAC(6)=0.5D0*XRAC*(COHO*(GAB*SHBC*DHCD*SGAD*SGBD+ &
                GAD*DSHCD*(HBC*SGAB*SGBD-e))+ &
                0.5D0*(CABRCD*GAB*GCD*SHAD*SHBC*SHBD+ &
                CABCD*GAB*DGCD*SHAD*SHBC*SHBD+ &
                CADRCD*GAD*GBC*SHAB*SHCD*SHBD+ &
                CADBC*GAD*GBC*SHAB*DSHCD*SHBD))


           PRAD=COHO*(GAB*SHBD*(HCD*SGAC*SGBC-e)+ &
                GAC*SHCD*(HBD*SGAB*SGBC-e))+ &
                0.5D0*(CABCD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CACBD*GAC*GBD*SHAB*SHCD*SHBC)
           PXRAD(1)=0.5D0*XRAD*(COHO*(DGAB*SHBD*(HCD*SGAC*SGBC-e)+ &
                GAC*SHCD*HBD*DSGAB*SGBC)+ &
                0.5D0*(CABRAB*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*DGAB*GCD*SHAC*SHBD*SHBC+ &
                CACRAB*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*GAC*GBD*DSHAB*SHCD*SHBC))
           PXRAD(2)=0.5D0*XRAD*(COHO*(GAB*SHBD*HCD*DSGAC*SGBC+ &
                DGAC*SHCD*(HBD*SGAB*SGBC-e))+ &
                0.5D0*(CABRAC*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*GAB*GCD*DSHAC*SHBD*SHBC+ &
                CACRAC*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*DGAC*GBD*SHAB*SHCD*SHBC))
             PXRAD(3)=0.5D0*(DXRAD*PRAD+XRAD* &
                0.5D0*(CABRAD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CACRAD*GAC*GBD*SHAB*SHCD*SHBC))
           PXRAD(4)=0.5D0*XRAD*(COHO*(GAB*SHBD*HCD*SGAC*DSGBC+ &
                GAC*SHCD*HBD*SGAB*DSGBC)+ &
                0.5D0*(CABRBC*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*GAB*GCD*SHAC*SHBD*DSHBC+ &
                CACRBC*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*GAC*GBD*SHAB*SHCD*DSHBC))
           PXRAD(5)=0.5D0*XRAD*(COHO*(GAB*DSHBD*(HCD*SGAC*SGBC-e)+ &
                GAC*SHCD*DHBD*SGAB*SGBC)+ &
                0.5D0*(CABRBD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*GAB*GCD*SHAC*DSHBD*SHBC+ &
                CACRBD*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*GAC*DGBD*SHAB*SHCD*SHBC))
           PXRAD(6)=0.5D0*XRAD*(COHO*(GAB*SHBD*DHCD*SGAC*SGBC+ &
                GAC*DSHCD*(HBD*SGAB*SGBC-e))+ &
                0.5D0*(CABRCD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*GAB*DGCD*SHAC*SHBD*SHBC+ &
                CACRCD*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*GAC*GBD*SHAB*DSHCD*SHBC))


        PRBC=COHO*(GAC*SHAB*(HBD*SGAD*SGCD-e)+ &
                GAB*SHAC*(HCD*SGAD*SGBD-e))+ &
                COOO*(GCD*SHBD*(HAB*SGAC*SGAD-e)+ &
                GBD*SHCD*(HAC*SGAB*SGAD-e))+0.5D0*( &
                CABCD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CACBD*GAC*GBD*SHAB*SHCD*SHAD)
        PXRBC(1)=0.5D0*XRBC*(COHO*(GAC*DSHAB*(HBD*SGAD*SGCD-e)+ &
                DGAB*SHAC*(HCD*SGAD*SGBD-e))+ &
                COOO*(GCD*SHBD*DHAB*SGAC*SGAD+ &
                GBD*SHCD*HAC*DSGAB*SGAD)+0.5D0*( &
                CABRAB*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*DGAB*GCD*SHAC*SHBD*SHAD+ &
                CACRAB*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*GAC*GBD*DSHAB*SHCD*SHAD))
        PXRBC(2)=0.5D0*XRBC*(COHO*(DGAC*SHAB*(HBD*SGAD*SGCD-e)+ &
                GAB*DSHAC*(HCD*SGAD*SGBD-e))+ &
                COOO*(GCD*SHBD*HAB*DSGAC*SGAD+ &
                GBD*SHCD*DHAC*SGAB*SGAD)+0.5D0*( &
                CABRAC*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*GAB*GCD*DSHAC*SHBD*SHAD+ &
                CACRAC*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*DGAC*GBD*SHAB*SHCD*SHAD))
        PXRBC(3)=0.5D0*XRBC*(COHO*(GAC*SHAB*HBD*DSGAD*SGCD+ &
                GAB*SHAC*HCD*DSGAD*SGBD)+ &
                COOO*(GCD*SHBD*HAB*SGAC*DSGAD+ &
                GBD*SHCD*HAC*SGAB*DSGAD)+0.5D0*( &
                CABRAD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*GAB*GCD*SHAC*SHBD*DSHAD+ &
                CACRAD*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*GAC*GBD*SHAB*SHCD*DSHAD ))
        PXRBC(4)=0.5D0*(DXRBC*PRBC+XRBC*0.5D0*( &
                CABRBC*GAB*GCD*SHAC*SHBD*SHAD+ &
                CACRBC*GAC*GBD*SHAB*SHCD*SHAD))
        PXRBC(5)=0.5D0*XRBC*(COHO*(GAC*SHAB*DHBD*SGAD*SGCD+ &
                GAB*SHAC*HCD*SGAD*DSGBD)+ &
                COOO*(GCD*DSHBD*(HAB*SGAC*SGAD-e)+ &
                DGBD*SHCD*(HAC*SGAB*SGAD-e))+0.5D0*( &
                CABRBD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*GAB*GCD*SHAC*DSHBD*SHAD+ &
                CACRBD*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*GAC*DGBD*SHAB*SHCD*SHAD))
        PXRBC(6)=0.5D0*XRBC*(COHO*(GAC*SHAB*HBD*SGAD*DSGCD+ &
                GAB*SHAC*DHCD*SGAD*SGBD)+ &
                COOO*(DGCD*SHBD*(HAB*SGAC*SGAD-e)+ &
                GBD*DSHCD*(HAC*SGAB*SGAD-e))+0.5D0*( &
                CABRCD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*GAB*DGCD*SHAC*SHBD*SHAD+ &
                CACRCD*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*GAC*GBD*SHAB*DSHCD*SHAD))


        PRBD=COHO*(GAD*SHAB*(HBC*SGAC*SGCD-e)+ &
                GAB*SHAD*(HCD*SGAC*SGBC-e))+ &
                COOO*(GCD*SHBC*(HAB*SGAC*SGAD-e)+ &
                GBC*SHCD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CABCD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CADBC*GAD*GBC*SHAB*SHCD*SHAC)
        PXRBD(1)=0.5D0*XRBD*(COHO*(GAD*DSHAB*(HBC*SGAC*SGCD-e)+ &
                DGAB*SHAD*(HCD*SGAC*SGBC-e))+ &
                COOO*(GCD*SHBC*DHAB*SGAC*SGAD+ &
                GBC*SHCD*HAD*DSGAB*SGAC)+0.5D0*( &
                CABRAB*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*DGAB*GCD*SHAD*SHBC*SHAC+ &
                CADRAB*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*GAD*GBC*DSHAB*SHCD*SHAC))
        PXRBD(2)=0.5D0*XRBD*(COHO*(GAD*SHAB*HBC*DSGAC*SGCD+ &
                GAB*SHAD*HCD*DSGAC*SGBC)+ &
                COOO*(GCD*SHBC*HAB*DSGAC*SGAD+ &
                GBC*SHCD*HAD*SGAB*DSGAC)+0.5D0*( &
                CABRAC*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*GAB*GCD*SHAD*SHBC*DSHAC+ &
                CADRAC*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*GAD*GBC*SHAB*SHCD*DSHAC))
        PXRBD(3)=0.5D0*XRBD*(COHO*(DGAD*SHAB*(HBC*SGAC*SGCD-e)+ &
                GAB*DSHAD*(HCD*SGAC*SGBC-e))+ &
                COOO*(GCD*SHBC*HAB*SGAC*DSGAD+ &
                GBC*SHCD*DHAD*SGAB*SGAC)+0.5D0*( &
                CABRAD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*GAB*GCD*DSHAD*SHBC*SHAC+ &
                CADRAD*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*DGAD*GBC*SHAB*SHCD*SHAC))
        PXRBD(4)=0.5D0*XRBD*(COHO*(GAD*SHAB*DHBC*SGAC*SGCD+ &
                GAB*SHAD*HCD*SGAC*DSGBC)+ &
                COOO*(GCD*DSHBC*(HAB*SGAC*SGAD-e)+ &
                DGBC*SHCD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CABRBC*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*GAB*GCD*SHAD*DSHBC*SHAC+ &
                CADRBC*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*GAD*DGBC*SHAB*SHCD*SHAC))
        PXRBD(5)=0.5D0*(DXRBD*PRBD+XRBD*0.5D0*( &
                CABRBD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CADRBD*GAD*GBC*SHAB*SHCD*SHAC))
        PXRBD(6)=0.5D0*XRBD*(COHO*(GAD*SHAB*HBC*SGAC*DSGCD+ &
                GAB*SHAD*DHCD*SGAC*SGBC)+ &
                COOO*(DGCD*SHBC*(HAB*SGAC*SGAD-e)+ &
                GBC*DSHCD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CABRCD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*GAB*DGCD*SHAD*SHBC*SHAC+ &
                CADRCD*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*GAD*GBC*SHAB*DSHCD*SHAC))

        PRCD=COHO*(GAD*SHAC*(HBC*SGAB*SGBD-e)+ &
                GAC*SHAD*(HBD*SGAB*SGBC-e))+ &
                COOO*(GBD*SHBC*(HAC*SGAB*SGAD-e)+ &
                GBC*SHBD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CACBD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CADBC*GAD*GBC*SHAC*SHBD*SHAB)
        PXRCD(1)=0.5D0*XRCD*(COHO*(GAD*SHAC*HBC*DSGAB*SGBD+ &
                GAC*SHAD*HBD*DSGAB*SGBC)+ &
                COOO*(GBD*SHBC*HAC*DSGAB*SGAD+ &
                GBC*SHBD*HAD*DSGAB*SGAC)+0.5D0*( &
                CACRAB*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*GAC*GBD*SHAD*SHBC*DSHAB+ &
                CADRAB*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*GAD*GBC*SHAC*SHBD*DSHAB))
        PXRCD(2)=0.5D0*XRCD*(COHO*(GAD*DSHAC*(HBC*SGAB*SGBD-e)+ &
                DGAC*SHAD*(HBD*SGAB*SGBC-e))+ &
                COOO*(GBD*SHBC*DHAC*SGAB*SGAD+ &
                GBC*SHBD*HAD*SGAB*DSGAC)+0.5D0*( &
                CACRAC*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*DGAC*GBD*SHAD*SHBC*SHAB+ &
                CADRAC*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*GAD*GBC*DSHAC*SHBD*SHAB))
        PXRCD(3)=0.5D0*XRCD*(COHO*(DGAD*SHAC*(HBC*SGAB*SGBD-e)+ &
                GAC*DSHAD*(HBD*SGAB*SGBC-e))+ &
                COOO*(GBD*SHBC*HAC*SGAB*DSGAD+ &
                GBC*SHBD*DHAD*SGAB*SGAC)+0.5D0*( &
                CACRAD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*GAC*GBD*DSHAD*SHBC*SHAB+ &
                CADRAD*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*DGAD*GBC*SHAC*SHBD*SHAB))
        PXRCD(4)=0.5D0*XRCD*(COHO*(GAD*SHAC*DHBC*SGAB*SGBD+ &
                GAC*SHAD*HBD*SGAB*DSGBC)+ &
                COOO*(GBD*DSHBC*(HAC*SGAB*SGAD-e)+ &
                DGBC*SHBD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CACRBC*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*GAC*GBD*SHAD*DSHBC*SHAB+ &
                CADRBC*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*GAD*DGBC*SHAC*SHBD*SHAB))
        PXRCD(5)=0.5D0*XRCD*(COHO*(GAD*SHAC*HBC*SGAB*DSGBD+ &
                GAC*SHAD*DHBD*SGAB*SGBC)+ &
                COOO*(DGBD*SHBC*(HAC*SGAB*SGAD-e)+ &
                GBC*DSHBD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CACRBD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*GAC*DGBD*SHAD*SHBC*SHAB+ &
                CADRBD*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*GAD*GBC*SHAC*DSHBD*SHAB))
        PXRCD(6)=0.5D0*(DXRCD*PRCD+XRCD*0.5D0*( &
                CACRCD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CADRCD*GAD*GBC*SHAC*SHBD*SHAB))

      DO 10 I=1,6
            TERM4(I)=PXRAB(I)+PXRAC(I)+PXRAD(I)+ &
                  PXRBC(I)+PXRBD(I)+PXRCD(I)
  10      CONTINUE
!  THE QUADRAPOLE-QUADRAPOLE INTERACTION
          XRAB=DMROH_13(5,RAB)
      XRAC=DMROH_13(5,RAC)
      XRAD=DMROH_13(5,RAD)
      XRBC=DMROO_13(5,RBC)
      XRBD=DMROO_13(5,RBD)
      XRCD=DMROO_13(5,RCD)

          DXRAB=DDMROH_13(5,RAB)
      DXRAC=DDMROH_13(5,RAC)
      DXRAD=DDMROH_13(5,RAD)
      DXRBC=DDMROO_13(5,RBC)
      DXRBD=DDMROO_13(5,RBD)
      DXRCD=DDMROO_13(5,RCD)

          SGAB=XIAOGOH_13(5,RAB)
      SGAC=XIAOGOH_13(5,RAC)
      SGAD=XIAOGOH_13(5,RAD)
      SGBC=XIAOGOO_13(5,RBC)
      SGBD=XIAOGOO_13(5,RBD)
      SGCD=XIAOGOO_13(5,RCD)

          DSGAB=DXIAOGOH_13(5,RAB)
      DSGAC=DXIAOGOH_13(5,RAC)
      DSGAD=DXIAOGOH_13(5,RAD)
      DSGBC=DXIAOGOO_13(5,RBC)
      DSGBD=DXIAOGOO_13(5,RBD)
      DSGCD=DXIAOGOO_13(5,RCD)

          CHOO=C5HOO
      COHO=C5OHO
      COOO=C5OOO
      CFHO=C5HO
        CABCD=C5HO*AOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
        CALL DAOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD,DEL)
      CABRAB=DEL(1)*CFHO
      CABRCD=DEL(2)*CFHO
      CABRAC=DEL(3)*CFHO
      CABRAD=DEL(4)*CFHO
      CABRBC=DEL(5)*CFHO
      CABRBD=DEL(6)*CFHO

      CACBD=C5HO*AOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD)
        CALL DAOMEGA_13(RAC,RBD,RAB,RAD,RBC,RCD,DEL)
      CACRAC=DEL(1)*CFHO
      CACRBD=DEL(2)*CFHO
      CACRAB=DEL(3)*CFHO
      CACRAD=DEL(4)*CFHO
      CACRBC=DEL(5)*CFHO
      CACRCD=DEL(6)*CFHO

        CADBC=C5HO*AOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD)
        CALL DAOMEGA_13(RAD,RBC,RAB,RAC,RBD,RCD,DEL)
      CADRAD=DEL(1)*CFHO
      CADRBC=DEL(2)*CFHO
      CADRAB=DEL(3)*CFHO
      CADRAC=DEL(4)*CFHO
      CADRBD=DEL(5)*CFHO
      CADRCD=DEL(6)*CFHO


         PRAB=COHO*(GAC*SHBC*(HBD*SGAD*SGCD-e)+ &
                GAD*SHBD*(HBC*SGAC*SGCD-e))+ &
                0.5D0*(CADBC*GAD*GBC*SHAC*SHBD*SHCD+ &
                CACBD*GAC*GBD*SHAD*SHBC*SHCD)
      PXRAB(1)=0.5D0*(DXRAB*PRAB+XRAB* &
                0.5D0*(CADRAB*GAD*GBC*SHAC*SHBD*SHCD+ &
                CACRAB*GAC*GBD*SHAD*SHBC*SHCD))
      PXRAB(2)=0.5D0*XRAB*(COHO*(DGAC*SHBC*(HBD*SGAD*SGCD-e)+ &
                GAD*SHBD*HBC*DSGAC*SGCD)+ &
                0.5D0*(CADRAC*GAD*GBC*SHAC*SHBD*SHCD &
                +CADBC*GAD*GBC*DSHAC*SHBD*SHCD+ &
                CACRAC*GAC*GBD*SHAD*SHBC*SHCD+ &
                CACBD*DGAC*GBD*SHAD*SHBC*SHCD))
       PXRAB(3)=0.5D0*XRAB*(COHO*(GAC*SHBC*HBD*DSGAD*SGCD+ &
                DGAD*SHBD*(HBC*SGAC*SGCD-e))+ &
                0.5D0*(CADRAD*GAD*GBC*SHAC*SHBD*SHCD+ &
                CADBC*DGAD*GBC*SHAC*SHBD*SHCD+ &
                CACRAD*GAC*GBD*SHAD*SHBC*SHCD+ &
            CACBD*GAC*GBD*DSHAD*SHBC*SHCD))
      PXRAB(4)=0.5D0*XRAB*(COHO*(GAC*DSHBC*(HBD*SGAD*SGCD-e)+ &
                GAD*SHBD*DHBC*SGAC*SGCD)+ &
                0.5D0*(CADRBC*GAD*GBC*SHAC*SHBD*SHCD+ &
                CADBC*GAD*DGBC*SHAC*SHBD*SHCD+ &
                CACRBC*GAC*GBD*SHAD*SHBC*SHCD+ &
            CACBD*GAC*GBD*SHAD*DSHBC*SHCD))
          PXRAB(5)=0.5D0*XRAB*(COHO*(GAC*SHBC*DHBD*SGAD*SGCD+ &
                GAD*DSHBD*(HBC*SGAC*SGCD-e))+ &
                0.5D0*(CADRBD*GAD*GBC*SHAC*SHBD*SHCD+ &
                CADBC*GAD*GBC*SHAC*DSHBD*SHCD+ &
                CACRBD*GAC*GBD*SHAD*SHBC*SHCD+ &
                CACBD*GAC*DGBD*SHAD*SHBC*SHCD))
             PXRAB(6)=0.5D0*XRAB*(COHO*(GAC*SHBC*HBD*SGAD*DSGCD+ &
                GAD*SHBD*HBC*SGAC*DSGCD)+ &
                0.5D0*(CADRCD*GAD*GBC*SHAC*SHBD*SHCD+ &
                CADBC*GAD*GBC*SHAC*SHBD*DSHCD+ &
                CACRCD*GAC*GBD*SHAD*SHBC*SHCD+ &
            CACBD*GAC*GBD*SHAD*SHBC*DSHCD))


           PRAC=COHO*(GAB*SHBC*(HCD*SGAD*SGBD-e)+ &
                GAD*SHCD*(HBC*SGAB*SGBD-e))+ &
                0.5D0*(CABCD*GAB*GCD*SHAD*SHBC*SHBD+ &
                CADBC*GAD*GBC*SHAB*SHCD*SHBD)
           PXRAC(1)=0.5D0*XRAC*(COHO*(DGAB*SHBC*(HCD*SGAD*SGBD-e)+ &
                GAD*SHCD*HBC*DSGAB*SGBD)+ &
                0.5D0*(CABRAB*GAB*GCD*SHAD*SHBC*SHBD+ &
                CABCD*DGAB*GCD*SHAD*SHBC*SHBD+ &
                CADRAB*GAD*GBC*SHAB*SHCD*SHBD+ &
                CADBC*GAD*GBC*DSHAB*SHCD*SHBD))
           PXRAC(2)=0.5D0*(DXRAC*PRAC+XRAC* &
                0.5D0*(CABRAC*GAB*GCD*SHAD*SHBC*SHBD+ &
                CADRAC*GAD*GBC*SHAB*SHCD*SHBD))
           PXRAC(3)=0.5D0*XRAC*(COHO*(GAB*SHBC*HCD*DSGAD*SGBD+ &
                DGAD*SHCD*(HBC*SGAB*SGBD-e))+ &
                0.5D0*(CABRAD*GAB*GCD*SHAD*SHBC*SHBD+ &
                CABCD*GAB*GCD*DSHAD*SHBC*SHBD+ &
                CADRAD*GAD*GBC*SHAB*SHCD*SHBD+ &
                CADBC*DGAD*GBC*SHAB*SHCD*SHBD ))
           PXRAC(4)=0.5D0*XRAC*(COHO*(GAB*DSHBC*(HCD*SGAD*SGBD-e)+ &
                GAD*SHCD*DHBC*SGAB*SGBD)+ &
                0.5D0*(CABRBC*GAB*GCD*SHAD*SHBC*SHBD+ &
            CABCD*GAB*GCD*SHAD*DSHBC*SHBD+ &
                CADRBC*GAD*GBC*SHAB*SHCD*SHBD+ &
            CADBC*GAD*DGBC*SHAB*SHCD*SHBD))
           PXRAC(5)=0.5D0*XRAC*(COHO*(GAB*SHBC*HCD*SGAD*DSGBD+ &
                GAD*SHCD*HBC*SGAB*DSGBD)+ &
                0.5D0*(CABRBD*GAB*GCD*SHAD*SHBC*SHBD+ &
            CABCD*GAB*GCD*SHAD*SHBC*DSHBD+ &
                CADRBD*GAD*GBC*SHAB*SHCD*SHBD+ &
                CADBC*GAD*GBC*SHAB*SHCD*DSHBD))
           PXRAC(6)=0.5D0*XRAC*(COHO*(GAB*SHBC*DHCD*SGAD*SGBD+ &
                GAD*DSHCD*(HBC*SGAB*SGBD-e))+ &
                0.5D0*(CABRCD*GAB*GCD*SHAD*SHBC*SHBD+ &
                CABCD*GAB*DGCD*SHAD*SHBC*SHBD+ &
                CADRCD*GAD*GBC*SHAB*SHCD*SHBD+ &
                CADBC*GAD*GBC*SHAB*DSHCD*SHBD))


           PRAD=COHO*(GAB*SHBD*(HCD*SGAC*SGBC-e)+ &
                GAC*SHCD*(HBD*SGAB*SGBC-e))+ &
                0.5D0*(CABCD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CACBD*GAC*GBD*SHAB*SHCD*SHBC)
           PXRAD(1)=0.5D0*XRAD*(COHO*(DGAB*SHBD*(HCD*SGAC*SGBC-e)+ &
                GAC*SHCD*HBD*DSGAB*SGBC)+ &
                0.5D0*(CABRAB*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*DGAB*GCD*SHAC*SHBD*SHBC+ &
                CACRAB*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*GAC*GBD*DSHAB*SHCD*SHBC))
           PXRAD(2)=0.5D0*XRAD*(COHO*(GAB*SHBD*HCD*DSGAC*SGBC+ &
                DGAC*SHCD*(HBD*SGAB*SGBC-e))+ &
                0.5D0*(CABRAC*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*GAB*GCD*DSHAC*SHBD*SHBC+ &
                CACRAC*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*DGAC*GBD*SHAB*SHCD*SHBC))
             PXRAD(3)=0.5D0*(DXRAD*PRAD+XRAD* &
                0.5D0*(CABRAD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CACRAD*GAC*GBD*SHAB*SHCD*SHBC))
           PXRAD(4)=0.5D0*XRAD*(COHO*(GAB*SHBD*HCD*SGAC*DSGBC+ &
                GAC*SHCD*HBD*SGAB*DSGBC)+ &
                0.5D0*(CABRBC*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*GAB*GCD*SHAC*SHBD*DSHBC+ &
                CACRBC*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*GAC*GBD*SHAB*SHCD*DSHBC))
           PXRAD(5)=0.5D0*XRAD*(COHO*(GAB*DSHBD*(HCD*SGAC*SGBC-e)+ &
                GAC*SHCD*DHBD*SGAB*SGBC)+ &
                0.5D0*(CABRBD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*GAB*GCD*SHAC*DSHBD*SHBC+ &
                CACRBD*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*GAC*DGBD*SHAB*SHCD*SHBC))
           PXRAD(6)=0.5D0*XRAD*(COHO*(GAB*SHBD*DHCD*SGAC*SGBC+ &
                GAC*DSHCD*(HBD*SGAB*SGBC-e))+ &
                0.5D0*(CABRCD*GAB*GCD*SHAC*SHBD*SHBC+ &
                CABCD*GAB*DGCD*SHAC*SHBD*SHBC+ &
                CACRCD*GAC*GBD*SHAB*SHCD*SHBC+ &
                CACBD*GAC*GBD*SHAB*DSHCD*SHBC))

        PRBC=COHO*(GAC*SHAB*(HBD*SGAD*SGCD-e)+ &
                GAB*SHAC*(HCD*SGAD*SGBD-e))+ &
                COOO*(GCD*SHBD*(HAB*SGAC*SGAD-e)+ &
                GBD*SHCD*(HAC*SGAB*SGAD-e))+0.5D0*( &
                CABCD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CACBD*GAC*GBD*SHAB*SHCD*SHAD)
        PXRBC(1)=0.5D0*XRBC*(COHO*(GAC*DSHAB*(HBD*SGAD*SGCD-e)+ &
                DGAB*SHAC*(HCD*SGAD*SGBD-e))+ &
                COOO*(GCD*SHBD*DHAB*SGAC*SGAD+ &
                GBD*SHCD*HAC*DSGAB*SGAD)+0.5D0*( &
                CABRAB*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*DGAB*GCD*SHAC*SHBD*SHAD+ &
                CACRAB*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*GAC*GBD*DSHAB*SHCD*SHAD))
        PXRBC(2)=0.5D0*XRBC*(COHO*(DGAC*SHAB*(HBD*SGAD*SGCD-e)+ &
                GAB*DSHAC*(HCD*SGAD*SGBD-e))+ &
                COOO*(GCD*SHBD*HAB*DSGAC*SGAD+ &
                GBD*SHCD*DHAC*SGAB*SGAD)+0.5D0*( &
                CABRAC*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*GAB*GCD*DSHAC*SHBD*SHAD+ &
                CACRAC*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*DGAC*GBD*SHAB*SHCD*SHAD))
        PXRBC(3)=0.5D0*XRBC*(COHO*(GAC*SHAB*HBD*DSGAD*SGCD+ &
                GAB*SHAC*HCD*DSGAD*SGBD)+ &
                COOO*(GCD*SHBD*HAB*SGAC*DSGAD+ &
                GBD*SHCD*HAC*SGAB*DSGAD)+0.5D0*( &
                CABRAD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*GAB*GCD*SHAC*SHBD*DSHAD+ &
                CACRAD*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*GAC*GBD*SHAB*SHCD*DSHAD ))
        PXRBC(4)=0.5D0*(DXRBC*PRBC+XRBC*0.5D0*( &
                CABRBC*GAB*GCD*SHAC*SHBD*SHAD+ &
                CACRBC*GAC*GBD*SHAB*SHCD*SHAD))
        PXRBC(5)=0.5D0*XRBC*(COHO*(GAC*SHAB*DHBD*SGAD*SGCD+ &
                GAB*SHAC*HCD*SGAD*DSGBD)+ &
                COOO*(GCD*DSHBD*(HAB*SGAC*SGAD-e)+ &
                DGBD*SHCD*(HAC*SGAB*SGAD-e))+0.5D0*( &
                CABRBD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*GAB*GCD*SHAC*DSHBD*SHAD+ &
                CACRBD*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*GAC*DGBD*SHAB*SHCD*SHAD))
        PXRBC(6)=0.5D0*XRBC*(COHO*(GAC*SHAB*HBD*SGAD*DSGCD+ &
                GAB*SHAC*DHCD*SGAD*SGBD)+ &
                COOO*(DGCD*SHBD*(HAB*SGAC*SGAD-e)+ &
                GBD*DSHCD*(HAC*SGAB*SGAD-e))+0.5D0*( &
                CABRCD*GAB*GCD*SHAC*SHBD*SHAD+ &
                CABCD*GAB*DGCD*SHAC*SHBD*SHAD+ &
                CACRCD*GAC*GBD*SHAB*SHCD*SHAD+ &
                CACBD*GAC*GBD*SHAB*DSHCD*SHAD))


        PRBD=COHO*(GAD*SHAB*(HBC*SGAC*SGCD-e)+ &
                GAB*SHAD*(HCD*SGAC*SGBC-e))+ &
                COOO*(GCD*SHBC*(HAB*SGAC*SGAD-e)+ &
                GBC*SHCD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CABCD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CADBC*GAD*GBC*SHAB*SHCD*SHAC)
        PXRBD(1)=0.5D0*XRBD*(COHO*(GAD*DSHAB*(HBC*SGAC*SGCD-e)+ &
                DGAB*SHAD*(HCD*SGAC*SGBC-e))+ &
                COOO*(GCD*SHBC*DHAB*SGAC*SGAD+ &
                GBC*SHCD*HAD*DSGAB*SGAC)+0.5D0*( &
                CABRAB*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*DGAB*GCD*SHAD*SHBC*SHAC+ &
                CADRAB*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*GAD*GBC*DSHAB*SHCD*SHAC))
        PXRBD(2)=0.5D0*XRBD*(COHO*(GAD*SHAB*HBC*DSGAC*SGCD+ &
                GAB*SHAD*HCD*DSGAC*SGBC)+ &
                COOO*(GCD*SHBC*HAB*DSGAC*SGAD+ &
                GBC*SHCD*HAD*SGAB*DSGAC)+0.5D0*( &
                CABRAC*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*GAB*GCD*SHAD*SHBC*DSHAC+ &
                CADRAC*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*GAD*GBC*SHAB*SHCD*DSHAC))
        PXRBD(3)=0.5D0*XRBD*(COHO*(DGAD*SHAB*(HBC*SGAC*SGCD-e)+ &
                GAB*DSHAD*(HCD*SGAC*SGBC-e))+ &
                COOO*(GCD*SHBC*HAB*SGAC*DSGAD+ &
                GBC*SHCD*DHAD*SGAB*SGAC)+0.5D0*( &
                CABRAD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*GAB*GCD*DSHAD*SHBC*SHAC+ &
                CADRAD*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*DGAD*GBC*SHAB*SHCD*SHAC))
        PXRBD(4)=0.5D0*XRBD*(COHO*(GAD*SHAB*DHBC*SGAC*SGCD+ &
                GAB*SHAD*HCD*SGAC*DSGBC)+ &
                COOO*(GCD*DSHBC*(HAB*SGAC*SGAD-e)+ &
                DGBC*SHCD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CABRBC*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*GAB*GCD*SHAD*DSHBC*SHAC+ &
                CADRBC*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*GAD*DGBC*SHAB*SHCD*SHAC))
        PXRBD(5)=0.5D0*(DXRBD*PRBD+XRBD*0.5D0*( &
                CABRBD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CADRBD*GAD*GBC*SHAB*SHCD*SHAC))
        PXRBD(6)=0.5D0*XRBD*(COHO*(GAD*SHAB*HBC*SGAC*DSGCD+ &
                GAB*SHAD*DHCD*SGAC*SGBC)+ &
                COOO*(DGCD*SHBC*(HAB*SGAC*SGAD-e)+ &
                GBC*DSHCD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CABRCD*GAB*GCD*SHAD*SHBC*SHAC+ &
                CABCD*GAB*DGCD*SHAD*SHBC*SHAC+ &
                CADRCD*GAD*GBC*SHAB*SHCD*SHAC+ &
                CADBC*GAD*GBC*SHAB*DSHCD*SHAC))

        PRCD=COHO*(GAD*SHAC*(HBC*SGAB*SGBD-e)+ &
                GAC*SHAD*(HBD*SGAB*SGBC-e))+ &
                COOO*(GBD*SHBC*(HAC*SGAB*SGAD-e)+ &
                GBC*SHBD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CACBD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CADBC*GAD*GBC*SHAC*SHBD*SHAB)
        PXRCD(1)=0.5D0*XRCD*(COHO*(GAD*SHAC*HBC*DSGAB*SGBD+ &
                GAC*SHAD*HBD*DSGAB*SGBC)+ &
                COOO*(GBD*SHBC*HAC*DSGAB*SGAD+ &
                GBC*SHBD*HAD*DSGAB*SGAC)+0.5D0*( &
                CACRAB*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*GAC*GBD*SHAD*SHBC*DSHAB+ &
                CADRAB*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*GAD*GBC*SHAC*SHBD*DSHAB))
        PXRCD(2)=0.5D0*XRCD*(COHO*(GAD*DSHAC*(HBC*SGAB*SGBD-e)+ &
                DGAC*SHAD*(HBD*SGAB*SGBC-e))+ &
                COOO*(GBD*SHBC*DHAC*SGAB*SGAD+ &
                GBC*SHBD*HAD*SGAB*DSGAC)+0.5D0*( &
                CACRAC*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*DGAC*GBD*SHAD*SHBC*SHAB+ &
                CADRAC*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*GAD*GBC*DSHAC*SHBD*SHAB))
        PXRCD(3)=0.5D0*XRCD*(COHO*(DGAD*SHAC*(HBC*SGAB*SGBD-e)+ &
                GAC*DSHAD*(HBD*SGAB*SGBC-e))+ &
                COOO*(GBD*SHBC*HAC*SGAB*DSGAD+ &
                GBC*SHBD*DHAD*SGAB*SGAC)+0.5D0*( &
                CACRAD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*GAC*GBD*DSHAD*SHBC*SHAB+ &
                CADRAD*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*DGAD*GBC*SHAC*SHBD*SHAB))
        PXRCD(4)=0.5D0*XRCD*(COHO*(GAD*SHAC*DHBC*SGAB*SGBD+ &
                GAC*SHAD*HBD*SGAB*DSGBC)+ &
                COOO*(GBD*DSHBC*(HAC*SGAB*SGAD-e)+ &
                DGBC*SHBD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CACRBC*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*GAC*GBD*SHAD*DSHBC*SHAB+ &
                CADRBC*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*GAD*DGBC*SHAC*SHBD*SHAB))
        PXRCD(5)=0.5D0*XRCD*(COHO*(GAD*SHAC*HBC*SGAB*DSGBD+ &
                GAC*SHAD*DHBD*SGAB*SGBC)+ &
                COOO*(DGBD*SHBC*(HAC*SGAB*SGAD-e)+ &
                GBC*DSHBD*(HAD*SGAB*SGAC-e))+0.5D0*( &
                CACRBD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CACBD*GAC*DGBD*SHAD*SHBC*SHAB+ &
                CADRBD*GAD*GBC*SHAC*SHBD*SHAB+ &
                CADBC*GAD*GBC*SHAC*DSHBD*SHAB))
        PXRCD(6)=0.5D0*(DXRCD*PRCD+XRCD*0.5D0*( &
                CACRCD*GAC*GBD*SHAD*SHBC*SHAB+ &
                CADRCD*GAD*GBC*SHAC*SHBD*SHAB))

      DO 20 I=1,6
            TERM5(I)=PXRAB(I)+PXRAC(I)+PXRAD(I)+ &
                  PXRBC(I)+PXRBD(I)+PXRCD(I)
  20      CONTINUE

      CALL DELECHO3P_13(RAB,RAC,RAD,RBC,RBD,RCD,PEL)
       DO 100 I=1,6
 100      DVEL(I)=TERM4(I)+TERM5(I)+PEL(I)
      RETURN
      END

      FUNCTION XIAOGOH_13(N,R)
      IMPLICIT REAL*8(A-H,O-Z)
      R0=2.6469057D0
      X=R-R0
      ALPH=2.5453276D0
      IF (N.EQ.4) THEN
      BETA=-1.1828544D-3
      ELSE
      BETA=0.0D0
      ENDIF
      XIAOGOH_13=1.0D0+BETA*DEXP(-ALPH*X)
      RETURN
      END

      FUNCTION DXIAOGOH_13(N,R)
      IMPLICIT REAL*8(A-H,O-Z)
      R0=2.6469057D0
      X=R-R0
      ALPH=2.5453276D0
      IF (N.EQ.4) THEN
      BETA=-1.1828544D-3
      ELSE
      BETA=0.0D0
      ENDIF
      DXIAOGOH_13=-ALPH*BETA*DEXP(-ALPH*X)
      RETURN
      END

      FUNCTION XIAOGOO_13(N,R)
      IMPLICIT REAL*8(A-H,O-Z)
      R0=2.5143D0
      X=R-R0
      ALPH=3.3522498D0
      IF (N.EQ.4) THEN
      BETA=-1.1918913D-4
      ELSE
      BETA=0.0D0
      ENDIF
      XIAOGOO_13=1.0D0+BETA*DEXP(-ALPH*X)
      RETURN
      END

      FUNCTION DXIAOGOO_13(N,R)
      IMPLICIT REAL*8(A-H,O-Z)
      R0=2.5143D0
      X=R-R0
      ALPH=3.3522498D0
      IF (N.EQ.4) THEN
      BETA=-1.1918913D-4
      ELSE
      BETA=0.0D0
      ENDIF
      DXIAOGOO_13=-ALPH*BETA*DEXP(-ALPH*X)
      RETURN
      END

             SUBROUTINE DAGSHOH_13(R,DAG,XIAOH)
      IMPLICIT REAL*8(A-H,O-Z)
      RM=1.8344D0
      ALPH=2.5453276D0
      DAG=((R/RM)**4)*DEXP(-ALPH*(R-RM))
      XIAOH=DTANH(ALPH*R)
      RETURN
      END

           SUBROUTINE DDAGSHOH_13(R,DAG,XIAOH)
      IMPLICIT REAL*8(A-H,O-Z)
      RM=1.8344D0
      ALPH=2.5453276D0
      DAG=((R/RM)**3)*DEXP(-ALPH*(R-RM))*(4.D0-ALPH*R)/RM
      XIAOH=ALPH*(1.0D0-(TANH(ALPH*R))**2)
      RETURN
      END

             SUBROUTINE DAGSHOO_13(R,DAG,XIAOH)
      IMPLICIT REAL*8(A-H,O-Z)
      RM=2.2818D0
      ALPH=3.3522498D0
      DAG=((R/RM)**4)*DEXP(-ALPH*(R-RM))
             XIAOH=DTANH(ALPH*R)
      RETURN
      END

           SUBROUTINE DDAGSHOO_13(R,DAG,XIAOH)
      IMPLICIT REAL*8(A-H,O-Z)
      RM=2.2818D0
      ALPH=3.3522498D0
      DAG=((R/RM)**3)*DEXP(-ALPH*(R-RM))*(4.D0-ALPH*R)/RM
            XIAOH=ALPH*(1.0D0-(TANH(ALPH*R))**2)
      RETURN
      END

      FUNCTION DAHOH_13(R)
      IMPLICIT REAL*8(A-H,O-Z)
           DAHOH_13=1.0-0.1688*EXP(-2.5453276*R)
      RETURN
      END

      FUNCTION DDAHOH_13(R)
      IMPLICIT REAL*8(A-H,O-Z)
      DDAHOH_13=0.1688*2.5453276*EXP(-2.5453276*R)
      RETURN
      END

      FUNCTION DAHOO_13(R)
      IMPLICIT REAL*8(A-H,O-Z)
      DAHOO_13=1.0+0.2875*EXP(-3.3522498*R)
      RETURN
      END

      FUNCTION DDAHOO_13(R)
      IMPLICIT REAL*8(A-H,O-Z)
      DDAHOO_13=-0.2875*3.3522498*EXP(-3.3522498*R)
      RETURN
      END


      FUNCTION DMROH_13(I,R)
      IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION D1(5),D2(5),N(5)
      D1(4)=5.0079875D0
      D1(5)=3.8428294D0
      D2(4)=2.129498D0
      D2(5)=2.5178884D0
      N(4)=4
      N(5)=5
      RM=1.8344D0
      R0=6.294894D0
      GAMA=2.5D0
      RO=(RM+GAMA*R0)/2.0D0
      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      RL=R**IDXL
      DMROH_13=C/RL
      RETURN
      END

      FUNCTION DDMROH_13(I,R)
      IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION D1(5),D2(5),N(5)
      D1(4)=5.0079875D0
      D1(5)=3.8428294D0
      D2(4)=2.129498D0
      D2(5)=2.5178884D0
      N(4)=4
      N(5)=5
      RM=1.8344D0
      R0=6.294894D0
      GAMA=2.5D0
      RO=(RM+GAMA*R0)/2.
      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      DPOL=-D1(I)*(1.+2.*D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      EXPOEL=C/(1.-EXPO)
      DC=DXL*EXPOEL*(-EXPO)*DPOL*DX
      RL=R**IDXL
      DDMROH_13=(DC/RL-DXL*C/(RL*R))
      RETURN
      END

      FUNCTION DMROO_13(I,R)
      IMPLICIT REAL*8 (A-H,O-Z)
           DIMENSION D1(5),D2(5),N(5)
      D1(4)=5.0079875D0
      D1(5)=3.8428294D0
      D2(4)=2.129498D0
      D2(5)=2.5178884D0
      N(4)=4
      N(5)=5
      RM=2.2818D0
      R0=5.661693D0
      GAMA=2.5D0
      RO=(RM+GAMA*R0)/2.
      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      RL=R**IDXL
      DMROO_13=C/RL
      RETURN
      END

      FUNCTION DDMROO_13(I,R)
      IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION D1(5),D2(5),N(5)
      D1(4)=5.0079875D0
      D1(5)=3.8428294D0
      D2(4)=2.129498D0
      D2(5)=2.5178884D0
      N(4)=4
      N(5)=5
      RM=2.2818D0
      R0=5.661693D0
      GAMA=2.5D0
      RO=(RM+GAMA*R0)/2.
      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      DPOL=-D1(I)*(1.+2.*D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      EXPOEL=C/(1.-EXPO)
      DC=DXL*EXPOEL*(-EXPO)*DPOL*DX
      RL=R**IDXL
      DDMROO_13=(DC/RL-DXL*C/(RL*R))
      RETURN
      END

      FUNCTION BOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
      IMPLICIT REAL*8(A-H,O-Z)
      CA=(RBD+RBC-RAC-RAD)/(2.*RAB)
      CB=(RBD+RAD-RAC-RBC)/(2.*RCD)
      SABCP=(RBD-RBC+RAC-RAD)*(RAD+RAC)/(2.*RAB*RCD)
      BOMEGA_13=CA*(3.0*CB*CB-1.0)-2.0*CB*SABCP
      RETURN
      END

      SUBROUTINE DBOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD,DB)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DB(6),DCA(6),DCB(6),DS(6)
      CA=(RBD+RBC-RAC-RAD)/(2.*RAB)
      DCA(1)=-CA/RAB
      DCA(2)=0.0
      DCA(3)=-1./(2.*RAB)
      DCA(4)=DCA(3)
      DCA(5)=-DCA(3)
      DCA(6)=DCA(5)
      CB=(RBD+RAD-RAC-RBC)/(2.*RCD)
      DCB(1)=0.0
      DCB(2)=-CB/RCD
      DCB(3)=-1./(2.*RCD)
      DCB(4)=-DCB(3)
      DCB(5)=DCB(3)
      DCB(6)=-DCB(3)
      CCC=2.*RAB*RCD
      CC1=RBD-RBC+RAC-RAD
      CC2=RAD+RAC
      SABCP=CC1*CC2/CCC
      DS(1)=-SABCP/RAB
      DS(2)=-SABCP/RCD
      DS(3)=(CC2+CC1)/CCC
      DS(4)=(CC1-CC2)/CCC
      DS(5)=-CC2/CCC
      DS(6)=CC2/CCC
      DO 10 I=1,6
!      BOMEGA=CA*(3.0*CB*CB-1.0)-2.0*CB*SABCP
      DB(I)=DCA(I)*(3.0*CB*CB-1.0)+CA*6.0*CB*DCB(I)-2.0* &
                 (DCB(I)*SABCP+CB*DS(I))
  10      CONTINUE
       RETURN
      END


      FUNCTION AOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD)
      IMPLICIT REAL*8(A-H,O-Z)
      CA=(RBD+RBC-RAC-RAD)/(2.*RAB)
      CB=(RBD+RAD-RAC-RBC)/(2.*RCD)
      SABCP=(RBD-RBC+RAC-RAD)*(RAD+RAC)/(2.*RAB*RCD)
      AOMEGA_13=1.-5.*CA*CA-5.*CB*CB-15.*(CA*CB)**2+ &
                  2.*(4.*CA*CB-SABCP)**2
      RETURN
      END

      SUBROUTINE DAOMEGA_13(RAB,RCD,RAC,RAD,RBC,RBD,DB)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DB(6),DCA(6),DCB(6),DS(6)
      CA=(RBD+RBC-RAC-RAD)/(2.*RAB)
      DCA(1)=-CA/RAB
      DCA(2)=0.0
      DCA(3)=-1./(2.*RAB)
      DCA(4)=DCA(3)
      DCA(5)=-DCA(3)
      DCA(6)=DCA(5)
      CB=(RBD+RAD-RAC-RBC)/(2.*RCD)
      DCB(1)=0.0
      DCB(2)=-CB/RCD
      DCB(3)=-1./(2.*RCD)
      DCB(4)=-DCB(3)
      DCB(5)=DCB(3)
      DCB(6)=-DCB(3)
      CCC=2.*RAB*RCD
      CC1=RBD-RBC+RAC-RAD
      CC2=RAD+RAC
      SABCP=CC1*CC2/CCC
      DS(1)=-SABCP/RAB
      DS(2)=-SABCP/RCD
      DS(3)=(CC2+CC1)/CCC
      DS(4)=(CC1-CC2)/CCC
      DS(5)=-CC2/CCC
      DS(6)=CC2/CCC
      DO 10 I=1,6
!      AOMEGA_13=1.-5.*CA*CA-5.*CB*CB-15.*(CA*CB)**2+
!     &            2.*(4.*CA*CB-SABCP)**2
      DB(I)=-10.0*(CA*DCA(I)+CB*DCB(I))-30.0*CA*CB*(CA*DCB(I)+ &
                    DCA(I)*CB)+4.0*(4.*CA*CB-SABCP)*(4.0*DCA(I)*CB+ &
                    4.0*CA*DCB(I)-DS(I))
 10      CONTINUE
      RETURN
      END


           FUNCTION FVHO3_13(RAB,RAC,RAD,RBC,RBD,RCD)
      IMPLICIT REAL*8(A-H,O-Z)
      EFM= DVTFM_13(RAB,RAC,RAD,RBC,RBD,RCD)
      EHF=FVHO3HF_13(RAB,RAC,RAD,RBC,RBD,RCD)
      Efit=vho3fit_13(RAB,RAC,RAD,RBC,RBD,RCD)
      fVHO3_13=efit+ehf+EFM
!w
!        PRINT*,'VS=',EFM
!        PRINT*,'VM=',EFIT
!        PRINT*,'PT=',EHF
!w
      RETURN
      END

      SUBROUTINE DFVHO3_13(RAB,RAC,RAD,RBC,RBD,RCD,DVR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DVR(6),ES(6),EB(6)
      CALL DFVHO3HF_13(RAB,RAC,RAD,RBC,RBD,RCD,DVR)
      CALL dvho3fit_10(RAB,RAC,RAD,RBC,RBD,RCD,EB)
      CALL DDVTFM_13(RAB,RAC,RAD,RBC,RBD,RCD,ES)
      DO 10 I=1,6
 10      DVR(I)=dvr(i)+EB(I)+es(i)
      RETURN
      END
!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
      FUNCTION vho3fit_13(rab,rac,rad,rbc,rbd,rcd)
        IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EXPO_13/GM1,gm2,AL,al2,BT,ROH,roh2,RO1,RO2,q30
!      dimension dvr(6)
         gz=    0.1250000000E+02
         AL=    0.7540057893E+00
      P=vfit_13(rab,rac,rad,rbc,rbd,rcd)+ &
       vfit_13(RAB,RAD,rac,RBD,RBC,RCD)+vfit_13(RAC,RAD,rab,RCD,RBC,RBD)
        vrep=vdvtrep_13(rab,rac,rad,rbc,rbd,rcd)
        sr=FSR_13(rbc,rbd,rcd)+FSR_13(rbd,rcd,rbc)+FSR_13(rcd,rbc,rbd)
      sr=(1.252705+sr)**2
      sr=exp(-gz*sr)
      q=rbc+rbd+rcd-8.90236
      g=al/4.0
      t=exp(-g*q**2)
      vho3fit_13=(p+vrep)*t*sr
!w
!        PRINT*,'VA=',p
!        PRINT*,'EXP=',t*sr
!        PRINT*,'VA*EXP=',p*t*sr
!        PRINT*,'VR*EXP=',vrep*t*sr
!w
      return
      end

      subroutine dvho3fit_10(rab,rac,rad,rbc,rbd,rcd,dvr)
        IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EXPO_13/GM1,gm2,AL,al2,BT,ROH,roh2,RO1,RO2,q30
      dimension dvr(6),vf1(6),vf2(6),vf3(6),vf4(6),dpv(6),dtsr(6)
         gz=    0.1250000000E+02
         AL=    0.7540057893E+00
      P=vfit_13(rab,rac,rad,rbc,rbd,rcd)+ &
        vfit_13(RAB,RAD,rac,RBD,RBC,RCD)+ &
        vfit_13(RAC,RAD,rab,RCD,RBC,RBD)
        vrep=vdvtrep_13(rab,rac,rad,rbc,rbd,rcd)
        pvrep=p+vrep
        sr0=FSR_13(rbc,rbd,rcd)+FSR_13(rbd,rcd,rbc)+FSR_13(rcd,rbc,rbd)
      sr1=(1.252705+sr0)**2
      sr=exp(-gz*sr1)
      dsr=-gz*2.*(1.252705+sr0)*sr
        CALL DFSR_13(RBC,RBD,RCD,DRBC1,DRBD1,DRCD1)
      CALL DFSR_13(RBD,RCD,RBC,DRBD2,DRCD2,DRBC2)
      CALL DFSR_13(RCD,RBC,RBD,DRCD3,DRBC3,DRBD3)
        SRBC=(DRBC1+DRBC2+DRBC3)*dsr
        SRBD=(DRBD1+DRBD2+DRBD3)*dsr
        SRCD=(DRCD1+DRCD2+DRCD3)*dsr
      q=rbc+rbd+rcd-8.90236
      g=al/4.0
      t=exp(-g*q**2)
      tq=-2.*g*q*t
      tsr=t*sr
      dtsr(1)=0.0
      dtsr(2)=0.0
      dtsr(3)=0.0
      dtsr(4)=sr*tq+t*srbc
      dtsr(5)=sr*tq+t*srbd
      dtsr(6)=sr*tq+t*srcd
        call dvdvtrep_13(rab,rac,rad,rbc,rbd,rcd,vf4)
      call dvfit_13(rab,rac,rad,rbc,rbd,rcd,vf1)
      call dvfit_13(RAB,RAD,rac,RBD,RBC,RCD,vf2)
           call dvfit_13(RAC,RAD,rab,RCD,RBC,RBD,vf3)
      dpv(1)=vf4(1)+vf1(1)+vf2(1)+vf3(3)
      dpv(2)=vf4(2)+vf1(2)+vf2(3)+vf3(1)
      dpv(3)=vf4(3)+vf1(3)+vf2(2)+vf3(2)
      dpv(4)=vf4(4)+vf1(4)+vf2(5)+vf3(5)
      dpv(5)=vf4(5)+vf1(5)+vf2(4)+vf3(6)
      dpv(6)=vf4(6)+vf1(6)+vf2(6)+vf3(4)
      dvr(1)=dpv(1)*tsr+pvrep*dtsr(1)
      dvr(2)=dpv(2)*tsr+pvrep*dtsr(2)
      dvr(3)=dpv(3)*tsr+pvrep*dtsr(3)
      dvr(4)=dpv(4)*tsr+pvrep*dtsr(4)
      dvr(5)=dpv(5)*tsr+pvrep*dtsr(5)
      dvr(6)=dpv(6)*tsr+pvrep*dtsr(6)
      return
      end


      function vfit_13(rab,rac,rad,rbc,rbd,rcd)
        IMPLICIT REAL*8(A-H,O-Z)
      parameter(np=19,nda=200)
      DIMENSION B(np),X(np)
      COMMON/EXPO_13/GM1,gm2,AL,al2,BT,ROH,roh2,RO1,RO2,q30
!         B0= -6.962208485341752E-002
         B0= -6.962208485341752E-002
         B( 1)=    0.1000000000E+01
         B( 2)=   -0.5489879243E+02
         B( 3)=    0.1120796600E+02
         B( 4)=   -0.1484596458E+03
         B( 5)=   -0.1338767542E+03
         B( 6)=    0.3335266723E+02
         B( 7)=    0.3253953639E+02
         B( 8)=   -0.3942563522E+02
         B( 9)=    0.2265036086E+01
         B(10)=    0.4642174751E+01
         B(11)=    0.1010778548E+02
         B(12)=    0.2957131995E+01
         B(13)=    0.5988692128E+00
         B(14)=    0.2334087256E+02
         B(15)=    0.3641005850E+02
         B(16)=   -0.1097782286E+00
         B(17)=    0.1772927025E+01
         B(18)=    0.1251542679E+01
         B(19)=    0.1414286015E+02
         GM1=   -0.7144442303E+00
         GM2=    0.2563402228E+00
         ROH=    0.5459891503E+01
         BT=    0.0000000000E+00
         ROH2=    0.5460209434E+01
         AL=    0.7540057893E+00
         gz=    0.1250000000E+02
         RO1=    0.2445710415E+01
         AL2=    0.0000000000E+00
         RO2=    0.4620820877E+01
         q30=    0.0000000000E+00
! SUM=  7.615499770707021E-003
! SUMreal=  8.387903806722165E-003

      CALL DVTf_13(RAB,RAC,rad,RBC,RBD,RCD,Q1,Q2,Q3,T)
      X(1)=1.0
      X(2)=Q1
      X(3)=Q1**2
            X(4)=Q1*Q3
            X(5)=Q1*Q3**2
            X(6)=Q1**2*Q3
            X(7)=(Q1*Q3)**2
            X(8)=Q1*Q3**3
       X(9)=Q1**3
            X(10)=Q1**3*Q3
            X(11)=Q1**2*Q3**3
            X(12)=Q1**3*Q3**2
            X(13)=Q1**3*Q3**3
      X(14)=Q2
        X(15)=Q2*Q3
        X(16)=Q2*Q1
      X(17)=Q2**2
      X(18)=Q2**2*Q3
        X(19)=Q2*Q3**2

      P=0.0
        DO 150 K=1,NP
      P=P+X(K)*B(K)
 150      CONTINUE
      vfit_13=b0*p*t
      return
      END

      subroutine dvfit_13(rab,rac,rad,rbc,rbd,rcd,dvr)
        IMPLICIT REAL*8(A-H,O-Z)
      parameter(np=19,nda=200)
      DIMENSION B(np),X(np),X1(np),X2(np),X3(np)
      dimension DQ1(6),DQ2(6),DQ3(6),DT(6),dvr(6)
      COMMON/EXPO_13/GM1,gm2,AL,al2,BT,ROH,roh2,RO1,RO2,q30
         B0= -6.962208485341752E-002
         B( 1)=    0.1000000000E+01
         B( 2)=   -0.5489879243E+02
         B( 3)=    0.1120796600E+02
         B( 4)=   -0.1484596458E+03
         B( 5)=   -0.1338767542E+03
         B( 6)=    0.3335266723E+02
         B( 7)=    0.3253953639E+02
         B( 8)=   -0.3942563522E+02
         B( 9)=    0.2265036086E+01
         B(10)=    0.4642174751E+01
         B(11)=    0.1010778548E+02
         B(12)=    0.2957131995E+01
         B(13)=    0.5988692128E+00
         B(14)=    0.2334087256E+02
         B(15)=    0.3641005850E+02
         B(16)=   -0.1097782286E+00
         B(17)=    0.1772927025E+01
         B(18)=    0.1251542679E+01
         B(19)=    0.1414286015E+02
         GM1=   -0.7144442303E+00
         GM2=    0.2563402228E+00
         ROH=    0.5459891503E+01
         BT=    0.0000000000E+00
         ROH2=    0.5460209434E+01
         AL=    0.7540057893E+00
         gz=    0.1250000000E+02
         RO1=    0.2445710415E+01
         AL2=    0.0000000000E+00
         RO2=    0.4620820877E+01
         q30=    0.0000000000E+00

      CALL DVTf_13(RAB,RAC,rad,RBC,RBD,RCD,Q1,Q2,Q3,T)
      X(1)=1.0
      X(2)=Q1
      X(3)=Q1**2
            X(4)=Q1*Q3
            X(5)=Q1*Q3**2
            X(6)=Q1**2*Q3
            X(7)=(Q1*Q3)**2
            X(8)=Q1*Q3**3
       X(9)=Q1**3
            X(10)=Q1**3*Q3
            X(11)=Q1**2*Q3**3
            X(12)=Q1**3*Q3**2
            X(13)=Q1**3*Q3**3
      X(14)=Q2
        X(15)=Q2*Q3
        X(16)=Q2*Q1
      X(17)=Q2**2
      X(18)=Q2**2*Q3
        X(19)=Q2*Q3**2

      X1(1)=0.0
      X1(2)=1.0
      X1(3)=2.*Q1
            X1(4)=Q3
            X1(5)=Q3**2
            X1(6)=2.*Q1*Q3
            X1(7)=2.*Q1*Q3**2
            X1(8)=Q3**3
       X1(9)=3.*Q1**2
            X1(10)=3.*Q1**2*Q3
            X1(11)=2.*Q1*Q3**3
            X1(12)=3.*Q1**2*Q3**2
            X1(13)=3.*Q1**2*Q3**3
      X1(14)=0.0
        X1(15)=0.0
        X1(16)=Q2
      X1(17)=0.0
      X1(18)=0.0
        X1(19)=0.0

        do 20 i=1,13
      X2(i)=0.0
 20      continue
      X2(14)=1.0
        X2(15)=Q3
        X2(16)=Q1
      X2(17)=2.*Q2
      X2(18)=2.*Q2*Q3
        X2(19)=Q3**2

      X3(1)=0.0
      X3(2)=0.0
      X3(3)=0.0
            X3(4)=Q1
            X3(5)=Q1*Q3*2.
            X3(6)=Q1**2
            X3(7)=2.*Q3*Q1**2
            X3(8)=3.*Q1*Q3**2
       X3(9)=0.0
            X3(10)=Q1**3
            X3(11)=3.*Q1**2*Q3**2
            X3(12)=Q1**3*Q3*2.
            X3(13)=3.*Q1**3*Q3**2
      X3(14)=0.0
        X3(15)=Q2
        X3(16)=0.0
      X3(17)=0.0
      X3(18)=Q2**2
        X3(19)=2.*Q2*Q3

      P=0.0
      pq1=0.0
      pq2=0.0
      pq3=0.0
        DO 150 K=1,NP
      P=P+X(K)*B(K)
      Pq1=Pq1+X1(K)*B(K)
      Pq2=Pq2+X2(K)*B(K)
      Pq3=Pq3+X3(K)*B(K)
 150      CONTINUE
      call DDVTf_13(RAB,RAC,rad,RBC,RBD,RCD,DQ1,DQ2,DQ3,DT)
      vfit_13=p*t
      do 99 i=1,6
      dvr(i)=(pq1*dq1(i)+pq2*dq2(i)+pq3*dq3(i))*t+p*dt(i)
      dvr(i)=dvr(i)*b0
 99      continue
      return
      END

      SUBROUTINE DVTf_13(RAB,RAC,rad,RBC,RBD,RCD,Q1,Q2,Q3,T)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EXPO_13/GM1,gm2,AL,al2,BT,ROH,roh2,RO1,RO2,q30
      Q1=RAB+RAC-ROH
      Q2=RBC-RO1
        Q3=FSR_13(RBC,RAB,RAC)+FSR_13(RAB,RAC,RBC)+FSR_13(RAC,RBC,RAB)
      q3=q3-q30
      q4=rbd+rcd-ro2
      q6=rad-roh2
      T=EXP(-(gm1*q1+GM2*Q1**2))
      T=T*EXP(-AL*Q2**2-al2*q4**2-bt*q6**2)
      RETURN
      END

      SUBROUTINE DDVTf_13(RAB,RAC,rad,RBC,RBD,RCD,DQ1,DQ2,DQ3,DT)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/EXPO_13/GM1,gm2,AL,al2,BT,ROH,roh2,RO1,RO2,q30
      dimension DQ1(6),DQ2(6),DQ3(6),DT(6)
      do 6 i=1,6
      dq1(i)=0.0
      dq2(i)=0.0
            dq3(i)=0.0
 6      dt(i)=0.0
      Q1=RAB+RAC-ROH
      Q2=RBC-RO1
        Q3=FSR_13(RBC,RAB,RAC)+FSR_13(RAB,RAC,RBC)+FSR_13(RAC,RBC,RAB)
      q3=q3-q30
      dq1(1)=1.0
      dq1(2)=1.0
        dq2(4)=1.0
      CALL DFSR_13(RBC,RAB,RAC,DRBC1,DRAB1,DRAC1)
      CALL DFSR_13(RAB,RAC,RBC,DRAB2,DRAC2,DRBC2)
      CALL DFSR_13(RAC,RBC,RAB,DRAC3,DRBC3,DRAB3)
      dQ3(1)=DRAB1+DRAB2+DRAB3
      dQ3(2)=DRAC1+DRAC2+DRAC3
      dQ3(4)=DRBC1+DRBC2+DRBC3
      q4=rbd+rcd-ro2
      q6=rad-roh2
      T=EXP(-(gm1*q1+GM2*Q1**2))
      T=T*EXP(-AL*Q2**2-al2*q4**2-bt*q6**2)
      tq1=-(gm1+2.*gm2*q1)*t
      tq2=-2.*al*q2*t
       tq4=-2.*al2*q4*t
       tq6=-2.*bt*q6*t
            dt(1)=tq1*dq1(1)+tq2*dq2(1)
            dt(2)=tq1*dq1(2)+tq2*dq2(2)
            dt(3)=tq1*dq1(3)+tq2*dq2(3)+tq6
            dt(4)=tq1*dq1(4)+tq2*dq2(4)
            dt(5)=tq1*dq1(5)+tq2*dq2(5)+tq4
            dt(6)=tq1*dq1(6)+tq2*dq2(6)+tq4
      RETURN
      END

       function vdvtrep_13(rab,rac,rad,rbc,rbd,rcd)
       implicit real*8(a-h,o-z)
       vdvtrep_13=DVTr_13(RAB,RAC,rad,RBC,RBD,RCD)+ &
           DVTr_13(RAB,RAD,rac,RBD,RBC,RCD) &
           +DVTr_13(RAC,RAD,rab,RCD,RBC,RBD)
       return
       end

       subroutine dvdvtrep_13(rab,rac,rad,rbc,rbd,rcd,dvr)
       implicit real*8(a-h,o-z)
       dimension dvr(6),dvr1(6),dvr2(6)
       call DDVTr_13(RAB,RAC,rad,RBC,RBD,RCD,dvr)
         call DDVTr_13(RAB,RAD,rac,RBD,RBC,RCD,dvr1)
       call DDVTr_13(RAC,RAD,rab,RCD,RBC,RBD,dvr2)
       dvr(1)=dvr(1)+dvr1(1)+dvr2(3)
       dvr(2)=dvr(2)+dvr1(3)+dvr2(1)
       dvr(3)=dvr(3)+dvr1(2)+dvr2(2)
       dvr(4)=dvr(4)+dvr1(5)+dvr2(5)
       dvr(5)=dvr(5)+dvr1(4)+dvr2(6)
       dvr(6)=dvr(6)+dvr1(6)+dvr2(4)
       do 6 i=1,6
 6       dvr(i)=dvr(i)
       return
       end

      function DVTr_13(RAB,RAC,rad,RBC,RBD,RCD)
      IMPLICIT REAL*8(A-H,O-Z)
      B1=    0.1717103975E-02
      B2=   -0.3744575433E-01
      B3=    0.2301744939E-01
      GM0=   -0.5330657145E-01
      GM1=    0.1582941367E+00
      GM2=    0.1598776430E+01
      Q20=    0.4262872754E+01
      Q30=    0.1252248542E+01
      Q40=    0.6845708355E+01
      ! SUM=  4.753287014195392E-003
      ! SUMreal=  3.665349272643972E-003
      
      Q4=RAB+RAc-rbc-q40
      Q2=RBC-q20
      Q3=FSR_13(RBC,RAB,RAC)+FSR_13(RAB,RAC,RBC)+FSR_13(RAC,RBC,RAB)
      q3=q3+q30
      t=exp(-gm2*q2**2-gm0*q4-gm1*q4**2)
      x=b1+b2*q3+b3*q3**2
      AL=0.7540057893E+00
      a0=-0.0014
      b0=0.0019
      Q2=RBC-2.4038
      Q3=RAB+RAc-9.763678
      q4=RAB+RAc-8.465767
      vgauss=(b0*exp(-2.*al*q4**2)+a0*exp(-4.*al*q3**2))*exp(-al*q2**2)
      dvtr_13=x*t+vgauss
      RETURN
      END

      subroutine DDVTr_13(RAB,RAC,rad,RBC,RBD,RCD,dvr)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension dvr(6)
      B1=    0.1717103975E-02
      B2=   -0.3744575433E-01
      B3=    0.2301744939E-01
      GM0=   -0.5330657145E-01
      GM1=    0.1582941367E+00
      GM2=    0.1598776430E+01
      Q20=    0.4262872754E+01
      Q30=    0.1252248542E+01
      Q40=    0.6845708355E+01
      Q4=RAB+RAc-rbc-q40
      Q2=RBC-q20
      Q3=FSR_13(RBC,RAB,RAC)+FSR_13(RAB,RAC,RBC)+FSR_13(RAC,RBC,RAB)
      q3=q3+q30
      CALL DFSR_13(RBC,RAB,RAC,DRBC1,DRAB1,DRAC1)
      CALL DFSR_13(RAB,RAC,RBC,DRAB2,DRAC2,DRBC2)
      CALL DFSR_13(RAC,RBC,RAB,DRAC3,DRBC3,DRAB3)
      Q3AB=DRAB1+DRAB2+DRAB3
      Q3AC=DRAC1+DRAC2+DRAC3
      Q3BC=DRBC1+DRBC2+DRBC3
      t=exp(-gm2*q2**2-gm0*q4-gm1*q4**2)
      tq2=-2.*gm2*q2*t
      tq4=-(gm0+2.*gm1*q4)*t
      x=b1+b2*q3+b3*q3**2
      xq3=b2+2.*b3*q3
      !      dvtr=x*t
      dvtrq2=x*tq2
      dvtrq3=xq3*t
      dvtrq4=x*tq4
      dvr(1)=dvtrq3*q3ab+dvtrq4
      dvr(2)=dvtrq3*q3ac+dvtrq4
      dvr(3)=0.0
      dvr(4)=dvtrq2+dvtrq3*q3bc-dvtrq4
      dvr(5)=0.0
      dvr(6)=0.0
      AL=0.7540057893E+00
      a0=-0.0014
      b0=0.0019
      Q2=RBC-2.4038
      Q3=RAB+RAc-9.763678
      q4=RAB+RAc-8.465767
      vex=exp(-al*q2**2)
      vexa0=a0*exp(-4.*al*q3**2)
      vexb0=b0*exp(-2.*al*q4**2)
      vgauss=(vexb0+vexa0)*vex
      VOO_13=-2.*al*q2*vgauss
      VOH_13=-al*(8.*q3*vexa0+4.*q4*vexb0)*vex
      dvr(1)=dvr(1)+VOH_13
      dvr(2)=dvr(2)+VOH_13
      dvr(4)=dvr(4)+VOO_13
      RETURN
      END

!============================================================
!-------------- The begin of 4-body Hartree-Fock Term ----------
           FUNCTION FVHO3HF_13(RAB,RAC,RAD,RBC,RBD,RCD)
      IMPLICIT REAL*8(A-H,O-Z)
!----- E=-0.355A.U.
      roh=   3.09282839000000
      roo= 3.07007816666667
      gam1=  3.50000000000
      gam2=  2.45000000000000
      v00=  1.005196250191485E-002
      v0=   1.00000000000000
      c1=   0.199385391766434
      c2=   3.67478464280366
      c3=   0.640803221096484
      c4=  -0.320433840343514
      c5=  -0.892125761340085
      c6=  -4.815695817096211E-002
      PAB=RAB-ROH
      PAC=RAC-ROH
      PAD=RAD-ROH
      PBC=RBC-ROO
      PBD=RBD-ROO
      PCD=RCD-ROO
      CC3=1.0D0/DSQRT(3.0D0)
      CC2=1.0D0/DSQRT(2.0D0)
      CC23=DSQRT(2.0D0/3.0D0)
      CC6=1.0D0/DSQRT(6.0D0)
      Q1=CC3*(PAB+PAC+PAD)
      Q2=CC2*(PAD-PAC)
      Q3=CC23*PAB-CC6*(PAD+PAC)
      Q4=CC3*(PCD+PBC+PBD)
      Q5=CC2*(PBC-PBD)
      Q6=CC23*PCD-CC6*(PBC+PBD)
      Q25=Q2+Q5
      Q36=Q3+Q6
      W1=Q1
      P0=V0+C1*W1+C2*Q4+C3*(Q2*Q2+Q3*Q3)+C4*(Q5*Q5+Q6*Q6)+ &
               C5*(Q2*Q5+Q3*Q6)+C6*Q36*(Q36*Q36-3.0D0*Q25*Q25)
      EQ=(1.0D0-TANH(GAM1*Q1/2.0D0))*(1.0D0-TANH(GAM2*Q4/2.0D0))
      FVHO3HF_13=v00*P0*EQ
      RETURN
      END

      SUBROUTINE DFVHO3HF_13(RAB,RAC,RAD,RBC,RBD,RCD,DVR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DVR(6)
!----- E=-0.355A.U.
      roh=   3.09282839000000
      roo= 3.07007816666667
      gam1=  3.50000000000
      gam2=  2.45000000000000
      v00=  1.005196250191485E-002
      v0=   1.00000000000000
      c1=   0.199385391766434
      c2=   3.67478464280366
      c3=   0.640803221096484
      c4=  -0.320433840343514
      c5=  -0.892125761340085
      c6=  -4.815695817096211E-002
      PAB=RAB-ROH
      PAC=RAC-ROH
      PAD=RAD-ROH
      PBC=RBC-ROO
      PBD=RBD-ROO
      PCD=RCD-ROO
      CC3=1.0D0/DSQRT(3.0D0)
      CC2=1.0D0/DSQRT(2.0D0)
      CC23=DSQRT(2.0D0/3.0D0)
      CC6=1.0D0/DSQRT(6.0D0)
      Q1=CC3*(PAB+PAC+PAD)
      Q2=CC2*(PAD-PAC)
      Q3=CC23*PAB-CC6*(PAD+PAC)
      Q4=CC3*(PCD+PBC+PBD)
      Q5=CC2*(PBC-PBD)
      Q6=CC23*PCD-CC6*(PBC+PBD)
      Q36=Q3+Q6
      Q25=Q2+Q5
      P0=V0+C1*Q1+C2*Q4+C3*(Q2*Q2+Q3*Q3)+C4*(Q5*Q5+Q6*Q6)+ &
               C5*(Q2*Q5+Q3*Q6)+C6*Q36*(Q36*Q36-3.0D0*Q25*Q25)
      EQ1=1.0D0-TANH(GAM1*Q1/2.0D0)
             DEQ1=-0.5D0*GAM1*(1.0D0-(TANH(GAM1*Q1/2.0D0))**2)
      EQ2=1.0D0-TANH(GAM2*Q4/2.0D0)
      DEQ2=-0.5D0*GAM2*(1.0D0-(TANH(GAM2*Q4/2.0D0))**2)
      EQ=EQ1*EQ2
      PQ1=C1*EQ+P0*EQ2*DEQ1
      PQ2=(2.0D0*C3*Q2+C5*Q5-6.*c6*q36*q25)*EQ
      PQ3=(2.D0*C3*Q3+C5*q6+3.*c6*(q36*q36-q25*q25))*EQ
      PQ4=C2*EQ+P0*EQ1*DEQ2
      PQ5=(2.0D0*C4*Q5+c5*q2-6.D0*C6*Q36*Q25)*EQ
      PQ6=(2.0D0*C4*Q6+c5*q3+3.D0*C6*(Q36*Q36-Q25*Q25))*EQ
      DVR(1)=CC3*PQ1+CC23*PQ3
      DVR(2)=CC3*PQ1-CC2*PQ2-CC6*PQ3
      DVR(3)=CC3*PQ1+CC2*PQ2-CC6*PQ3
      DVR(6)=CC3*PQ4+CC23*PQ6
      DVR(5)=CC3*PQ4-CC2*PQ5-CC6*PQ6
      DVR(4)=CC3*PQ4+CC2*PQ5-CC6*PQ6
      do 6 i=1,6
 6      dvr(i)=dvr(i)*v00
      RETURN
      END
!-------------- The end of 4-body Hartree-Fock Term ----------
      FUNCTION DVTFM_13(RAB,RAC,RAD,RBC,RBD,RCD)
      IMPLICIT REAL*8(A-H,O-Z)
      M=2
      A0=1.718D0
            AOH=0.035D0
      AOO=0.120D0
      AOO2=0.650D0
             R0OH=10.0882205963135
      R0OO=10.0937557220459
      RA=RAB+RAC+RAD
      RBCD=RBC+RBD+RCD
      Q1=(RA-R0OH)**2
      Q4=(RBCD-R0OO)**2
      Q5=(RBC-RBD)/DSQRT(2.0D0)
      Q6=DSQRT(2.0D0/3.0D0)*RCD-(RBC+RBD)/DSQRT(6.0D0)
      Q56=Q5*Q5+Q6*Q6
      SR=FSR_13(RBC,RBD,RCD)+FSR_13(RBD,RCD,RBC)+FSR_13(RCD,RBC,RBD)
      DVTFM_13=A0*((1.252705D0+SR)**M) &
                  *DEXP(-AOO2*Q56-AOH*Q1-AOO*Q4)
      RETURN
      END

      SUBROUTINE DDVTFM_13(RAB,RAC,RAD,RBC,RBD,RCD,EVR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION EVR(6)
      M=2
      A0=1.718D0
            AOH=0.035D0
      AOO=0.120D0
      AOO2=0.650D0
             R0OH=10.0882205963135
      R0OO=10.0937557220459
      RA=RAB+RAC+RAD
      RBCD=RBC+RBD+RCD
      Q1=RA-R0OH
      Q4=RBCD-R0OO
      S2=DSQRT(2.0D0)
      S23=DSQRT(2.0D0/3.0D0)
      S6=DSQRT(6.0D0)
      Q5=(RBC-RBD)/S2
      Q6=S23*RCD-(RBC+RBD)/S6
      Q56=Q5*Q5+Q6*Q6
      SR=1.252705+FSR_13(RBC,RBD,RCD)+FSR_13(RBD,RCD,RBC)+ &
       FSR_13(RCD,RBC,RBD)
        CALL DFSR_13(RBC,RBD,RCD,DRBC1,DRBD1,DRCD1)
      CALL DFSR_13(RBD,RCD,RBC,DRBD2,DRCD2,DRBC2)
      CALL DFSR_13(RCD,RBC,RBD,DRCD3,DRBC3,DRBD3)
        SRBC=DRBC1+DRBC2+DRBC3
        SRBD=DRBD1+DRBD2+DRBD3
        SRCD=DRCD1+DRCD2+DRCD3
!      DVTFM_13=A0*(SR**M)*DEXP(-AOO2*Q56-AOH*Q1*Q1-AOO*Q4*Q4)
      AE=A0*DEXP(-AOO2*Q56-AOH*Q1*Q1-AOO*Q4*Q4)
      FMQ5=-AOO2*2.*Q5
      FMQ6=-AOO2*2.*Q6
      FMQ1=-2.*AOH*Q1
      FMQ4=-2.*AOO*Q4
      FMRAB=FMQ1
      FMRAC=FMQ1
      FMRAD=FMQ1
      FMRBC=FMQ4+FMQ5/S2-FMQ6/S6
      FMRBD=FMQ4-FMQ5/S2-FMQ6/S6
      FMRCD=FMQ4+S23*FMQ6
      SRM=SR**M
      SRM1=FLOAT(M)*SR**(M-1)
      EVR(1)=AE*SRM*FMRAB
      EVR(2)=AE*SRM*FMRAC
      EVR(3)=AE*SRM*FMRAD
      EVR(4)=AE*(SRBC*SRM1+SRM*FMRBC)
      EVR(5)=AE*(SRBD*SRM1+SRM*FMRBD)
      EVR(6)=AE*(SRCD*SRM1+SRM*FMRCD)
      RETURN
      END

! -----------------------------------------------
             FUNCTION FSR_13(RI,RJ,RK)
      IMPLICIT REAL*8(A-H,O-Z)
      FSR_13=(RI*RI-RJ*RJ-RK*RK)/(2.0D0*RJ*RK)
      RETURN
      END
             SUBROUTINE DFSR_13(RI,RJ,RK,DRI,DRJ,DRK)
      IMPLICIT REAL*8(A-H,O-Z)
      RJK=RJ*RK
      RIJK=RI*RI-RJ*RJ-RK*RK
      DRI=RI/RJK
      DRJ=-(2.0D0*RJ*RJK+RK*RIJK)/(2.0D0*RJK**2)
        DRK=-(2.0D0*RK*RJK+RJ*RIJK)/(2.0D0*RJK**2)
      RETURN
      END

