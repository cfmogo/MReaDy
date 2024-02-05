      subroutine  vh2o2_switch(i,V,dv)
      use variables
      double precision, dimension(12) :: dv,dV1
      integer      :: gi3,gi4,gi5,gi6
      
      double precision :: rxij,ryij,rzij,rij
      double precision :: r1xij,r1yij,r1zij
      double precision :: r2xij,r2yij,r2zij
      double precision :: r3xij,r3yij,r3zij
      double precision :: r4xij,r4yij,r4zij
      double precision :: r5xij,r5yij,r5zij
      double precision :: r6xij,r6yij,r6zij
      double precision :: dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6
      double precision :: r1,r2,r3,r4,r5,r6
      
      double precision, dimension(6) :: rin,r_4atoms,dvdr_4atoms
      double precision ::fs
      double precision :: V2,V1,V,dV2,dV3,dvdiat,V4 ,vdiat ,v3
      double precision :: f3,f4,vh2o_f5,vh2ot_f6,vho2_f4
      double precision :: dfsdr1,dfsdr2,dfsdr3,dfsdr4,dfsdr5,dfsdr6
      
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
      
      
      
      
      
      !Geting the coordenates
      
      !The positions are predefined in the matrix "group" with H atoms first and Oxi later
      gi3=group(i,3)      ! H atom
      gi4=group(i,4)      ! H atom
      gi5=group(i,5)      ! H atom
      gi6=group(i,6)      ! O atom
      
      
      
      
      
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
      
      
      
      
      
      
      !       The quadriatomic surface
      
      call h2o2sur_10(6,r_4atoms,V1,dvdr_4atoms)
      !        !print*,'r_4atoms'
      !        !print*,r_4atoms
      !        !print*,'V1'
      !        !print*,V1
      !        !print*,'dvdr_4atoms'
      !        !print*,dvdr_4atoms
      
      
      
      
      
      
      
      dv1dr1=dvdr_4atoms(1)
      dv1dr2=dvdr_4atoms(2)
      dv1dr3=dvdr_4atoms(3)
      dv1dr4=dvdr_4atoms(4)
      dv1dr5=dvdr_4atoms(5)
      dv1dr6=dvdr_4atoms(6)
      
      
      
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
      
      
      
      
      
      
      
      
      !            !print*,'entrada',rin
      
      
      !'H2O2 -> O(4) + H2O'
      
      if ((r1.ge.in_bound).and.(r4.ge.in_bound).and. &
      (r5.ge.in_bound)) then
      !           !print*,'exit0'
      
      if(group(i,9).eq.0) group(i,9)=302
      
      
      
      if((group(i,9).ne.301).or.(group(i,9).ne.302))then
      !print*,'exchange in molecular complex: group',group(i,9),tempo
      group(i,9)=302
      end if
      
      if (group(i,9).eq.302) then
      ! H2O(1A') function
      V2=vh2o_f5(rin(2),rin(3),rin(6))+7.1955d-2*4.3597482D-18
      call dervh2o_f5(rin(2),rin(3),rin(6),dv2dr2,dv2dr3,dv2dr6)
      
      else if (group(i,9).eq.301) then
      ! H2O(3A") function
      V2=vh2ot_f6(rin(2),rin(3),rin(6))
      call dervh2ot_f6(rin(2),rin(3),rin(6),dv2dr2,dv2dr3,dv2dr6)
      
      end if
      
      ! Switch function
      
      
      fs=f3(rin(1),rin(4),rin(5))
      call df3(rin(1),rin(4),rin(5),dfsdr1,dfsdr4,dfsdr5)
      
      
      ! for the diatomic potencials
      
      call pothoq_f17(rin(4),Vdiat,dVdiat)
      
      
      dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
      dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
      dv(6) = dv(6)  + fs * dVdiat*r4zij/r4
      
      dv(10)= dv(10) - fs * dVdiat*r4xij/r4
      dv(11)= dv(11) - fs * dVdiat*r4yij/r4
      dv(12)= dv(12) - fs * dVdiat*r4zij/r4
      
      
      V4 = V4 + Vdiat
      
      call pothoq_f17(rin(5),Vdiat,dVdiat)
      
      
      dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
      dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
      dv(3) = dv(3)  + fs * dVdiat*r5zij/r5
      
      dv(10)= dv(10) - fs * dVdiat*r5xij/r5
      dv(11)= dv(11) - fs * dVdiat*r5yij/r5
      dv(12)= dv(12) - fs * dVdiat*r5zij/r5
      
      
      V4 = V4 + Vdiat
      
      
      
      call poto2q_f18(rin(1),Vdiat,dVdiat)
      
      !            !print*,'rin(3),Vdiat,dvdiat'
      !            !print*,rin(3),Vdiat,dvdiat
      
      
      
      dv(7) = dv(7) + fs * dVdiat*r1xij/r1
      dv(8) = dv(8) + fs * dVdiat*r1yij/r1
      dv(9) = dv(9) + fs * dVdiat*r1zij/r1
      
      dv(10)= dv(10)- fs * dVdiat*r1xij/r1
      dv(11)= dv(11)- fs * dVdiat*r1yij/r1
      dv(12)= dv(12)- fs * dVdiat*r1zij/r1
      
      
      V4 = V4 + Vdiat
      
      V=V1*(1-fs)+(V2+V4)*fs
      !        !print*,rin
      !        !print*,'V=',V,' V1=',V1,' V2=',V2
      
      dv=dv + dv1*(1-fs)
      
      
      dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr5*r5xij/r5)
      dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr5*r5yij/r5)
      dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr5*r5zij/r5)
      
      dv(4) = dv(4) + (v4 + v2 - v1) * ( dfsdr4*r4xij/r4)
      dv(5) = dv(5) + (v4 + v2 - v1) * ( dfsdr4*r4yij/r4)
      dv(6) = dv(6) + (v4 + v2 - v1) * ( dfsdr4*r4zij/r4)
      
      dv(7) = dv(7) + (v4 + v2 - v1) * ( dfsdr1*r1xij/r1)
      dv(8) = dv(8) + (v4 + v2 - v1) * ( dfsdr1*r1yij/r1)
      dv(9) = dv(9) + (v4 + v2 - v1) * ( dfsdr1*r1zij/r1)
      
      dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr1*r1xij/r1 - dfsdr4*r4xij/r4 - dfsdr5*r5xij/r5)
      dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr1*r1yij/r1 - dfsdr4*r4yij/r4 - dfsdr5*r5yij/r5)
      dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr1*r1zij/r1 - dfsdr4*r4zij/r4 - dfsdr5*r5zij/r5)
      
      
      dv(1) = dv(1) + fs*( dv2dr2*r2xij/r2 + dv2dr3*r3xij/r3)
      dv(2) = dv(2) + fs*( dv2dr2*r2yij/r2 + dv2dr3*r3yij/r3)
      dv(3) = dv(3) + fs*( dv2dr2*r2zij/r2 + dv2dr3*r3zij/r3)
      
      dv(4) = dv(4) + fs*(-dv2dr2*r2xij/r2 + dv2dr6*r6xij/r6)
      dv(5) = dv(5) + fs*(-dv2dr2*r2yij/r2 + dv2dr6*r6yij/r6)
      dv(6) = dv(6) + fs*(-dv2dr2*r2zij/r2 + dv2dr6*r6zij/r6)
      
      dv(7) = dv(7) + fs*(-dv2dr3*r3xij/r3 - dv2dr6*r6xij/r6)
      dv(8) = dv(8) + fs*(-dv2dr3*r3yij/r3 - dv2dr6*r6yij/r6)
      dv(9) = dv(9) + fs*(-dv2dr3*r3zij/r3 - dv2dr6*r6zij/r6)
      
      
      
      
      
      
      !'H2O2 -> O(3) + H2O'
      
      else if ((r1.ge.in_bound).and.(r3.ge.in_bound).and. &
      (r6.ge.in_bound)) then
      
      if(group(i,9).eq.0) group(i,9)=302
      
      
      
      if((group(i,9).ne.301).or.(group(i,9).ne.302))then
      print*,'exchange in molecular complex: group',group(i,9),tempo
      group(i,9)=302
      end if
      
      if (group(i,9).eq.302) then
      ! H2O(1A') function
      V2=vh2o_f5(rin(2),rin(4),rin(5))+7.1955d-2*4.3597482D-18
      call dervh2o_f5(rin(2),rin(4),rin(5),dv2dr2,dv2dr4,dv2dr5)
      
      else if (group(i,9).eq.301) then
      ! H2O(3A") function
      V2=vh2ot_f6(rin(2),rin(4),rin(5))
      call dervh2ot_f6(rin(2),rin(4),rin(5),dv2dr2,dv2dr4,dv2dr5)
      
      end if
      
      ! Switch function
      
      
      fs=f3(rin(1),rin(3),rin(6))
      call df3(rin(1),rin(3),rin(6),dfsdr1,dfsdr3,dfsdr6)
      
      
      ! for the diatomic potencials
      
      call pothoq_f17(rin(3),Vdiat,dVdiat)
      
      
      dv(1) = dv(1) + fs * dVdiat*r3xij/r3
      dv(2) = dv(2) + fs * dVdiat*r3yij/r3
      dv(3) = dv(3) + fs * dVdiat*r3zij/r3
      
      dv(7) = dv(7) - fs * dVdiat*r3xij/r3
      dv(8) = dv(8) - fs * dVdiat*r3yij/r3
      dv(9) = dv(9) - fs * dVdiat*r3zij/r3
      
      
      V4 = V4 + Vdiat
      
      call pothoq_f17(rin(6),Vdiat,dVdiat)
      
      
      dv(4) = dv(4) + fs * dVdiat*r6xij/r6
      dv(5) = dv(5) + fs * dVdiat*r6yij/r6
      dv(6) = dv(6) + fs * dVdiat*r6zij/r6
      
      dv(7) = dv(7) - fs * dVdiat*r6xij/r6
      dv(8) = dv(8) - fs * dVdiat*r6yij/r6
      dv(9) = dv(9) - fs * dVdiat*r6zij/r6
      
      
      V4 = V4 + Vdiat
      
      
      
      call poto2q_f18(rin(1),Vdiat,dVdiat)
      
      !            !print*,'rin(3),Vdiat,dvdiat'
      !            !print*,rin(3),Vdiat,dvdiat
      
      
      
      dv(7) = dv(7) + fs * dVdiat*r1xij/r1
      dv(8) = dv(8) + fs * dVdiat*r1yij/r1
      dv(9) = dv(9) + fs * dVdiat*r1zij/r1
      
      dv(10)= dv(10)- fs * dVdiat*r1xij/r1
      dv(11)= dv(11)- fs * dVdiat*r1yij/r1
      dv(12)= dv(12)- fs * dVdiat*r1zij/r1
      
      
      V4 = V4 + Vdiat
      
      V=V1*(1-fs)+(V2+V4)*fs
      !        !print*,rin
      !        !print*,'V=',V,' V1=',V1,' V2=',V2
      
      dv=dv + dv1*(1-fs)
      
      
      dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr3*r3xij/r3)
      dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr3*r3yij/r3)
      dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr3*r3zij/r3)
      
      dv(4) = dv(4) + (v4 + v2 - v1) * ( dfsdr6*r6xij/r6)
      dv(5) = dv(5) + (v4 + v2 - v1) * ( dfsdr6*r6yij/r6)
      dv(6) = dv(6) + (v4 + v2 - v1) * ( dfsdr6*r6zij/r6)
      
      dv(7) = dv(7) + (v4 + v2 - v1) * ( dfsdr1*r1xij/r1 - dfsdr3*r3xij/r3 - dfsdr6*r6xij/r6)
      dv(8) = dv(8) + (v4 + v2 - v1) * ( dfsdr1*r1yij/r1 - dfsdr3*r3yij/r3 - dfsdr6*r6yij/r6)
      dv(9) = dv(9) + (v4 + v2 - v1) * ( dfsdr1*r1zij/r1 - dfsdr3*r3zij/r3 - dfsdr6*r6zij/r6)
      
      dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr1*r1xij/r1)
      dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr1*r1yij/r1)
      dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr1*r1zij/r1)
      
      
      dv(1) = dv(1) + fs*( dv2dr2*r2xij/r2 + dv2dr5*r5xij/r5)
      dv(2) = dv(2) + fs*( dv2dr2*r2yij/r2 + dv2dr5*r5yij/r5)
      dv(3) = dv(3) + fs*( dv2dr2*r2zij/r2 + dv2dr5*r5zij/r5)
      
      dv(4) = dv(4) + fs*(-dv2dr2*r2xij/r2 + dv2dr4*r4xij/r4)
      dv(5) = dv(5) + fs*(-dv2dr2*r2yij/r2 + dv2dr4*r4yij/r4)
      dv(6) = dv(6) + fs*(-dv2dr2*r2zij/r2 + dv2dr4*r4zij/r4)
      
      dv(7) = dv(7) + fs*(-dv2dr4*r4xij/r4 - dv2dr5*r5xij/r5)
      dv(8) = dv(8) + fs*(-dv2dr4*r4yij/r4 - dv2dr5*r5yij/r5)
      dv(9) = dv(9) + fs*(-dv2dr4*r4zij/r4 - dv2dr5*r5zij/r5)
      
      
      
      
      
      
      !'H2O2 -> H(2) + HO2'
      
      else if ((r2.ge.in_bound).and.(r4.ge.in_bound).and. &
      (r6.ge.in_bound)) then
      !            !print*,'exit5'
      
      if(group(i,9).eq.0) group(i,9)=310
      
      
      
      
      !print*,'exchange in complex into 310 from',group(i,9),tempo
      group(i,9)=310
      
      
      
      ! HO2 function
      V2=vho2_f4(r1,r3,r5)
      call dervho2_f4(r1,r3,r5,dv2dr1,dv2dr3,dv2dr5)
      
      
      
      
      ! Switch function
      
      
      fs=f3(rin(2),rin(4),rin(6))
      call df3(rin(2),rin(4),rin(6),dfsdr2,dfsdr4,dfsdr6)
      
      
      ! for the diatomic potencials
      
      call pothoq_f17(rin(4),Vdiat,dVdiat)
      
      
      dv(4) = dv(4) + fs * dVdiat*r4xij/r4
      dv(5) = dv(5) + fs * dVdiat*r4yij/r4
      dv(6) = dv(6) + fs * dVdiat*r4zij/r4
      
      dv(10) = dv(10) - fs * dVdiat*r4xij/r4
      dv(11) = dv(11) - fs * dVdiat*r4yij/r4
      dv(12) = dv(12) - fs * dVdiat*r4zij/r4
      
      
      V4 = V4 + Vdiat
      
      call pothoq_f17(rin(6),Vdiat,dVdiat)
      
      
      dv(4) = dv(4) + fs * dVdiat*r6xij/r6
      dv(5) = dv(5) + fs * dVdiat*r6yij/r6
      dv(6) = dv(6) + fs * dVdiat*r6zij/r6
      
      dv(7) = dv(7) - fs * dVdiat*r6xij/r6
      dv(8) = dv(8) - fs * dVdiat*r6yij/r6
      dv(9) = dv(9) - fs * dVdiat*r6zij/r6
      
      
      V4 = V4 + Vdiat
      
      
      call pothht_f14(r2,Vdiat,dVdiat)
      
      
      
      dv(1) = dv(1) + fs * dVdiat*r2xij/r2
      dv(2) = dv(2) + fs * dVdiat*r2yij/r2
      dv(3) = dv(3) + fs * dVdiat*r2zij/r2
      
      dv(4) = dv(4) - fs * dVdiat*r2xij/r2
      dv(5) = dv(5) - fs * dVdiat*r2yij/r2
      dv(6) = dv(6) - fs * dVdiat*r2zij/r2
      
      
      V4 = V4 + Vdiat
      
      V=V1*(1-fs)+(V2+V4)*fs
      !        !print*,rin
      !        !print*,'V=',V,' V1=',V1,' V2=',V2
      
      dv=dv + dv1*(1-fs)
      
      
      dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr2*r2xij/r2)
      dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr2*r2yij/r2)
      dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr2*r2zij/r2)
      
      dv(4) = dv(4) + (v4 + v2 - v1) * (-dfsdr2*r2xij/r2 + dfsdr4*r4xij/r4 + dfsdr6*r6xij/r6)
      dv(5) = dv(5) + (v4 + v2 - v1) * (-dfsdr2*r2yij/r2 + dfsdr4*r4yij/r4 + dfsdr6*r6yij/r6)
      dv(6) = dv(6) + (v4 + v2 - v1) * (-dfsdr2*r2zij/r2 + dfsdr4*r4zij/r4 + dfsdr6*r6zij/r6)
      
      dv(7) = dv(7) + (v4 + v2 - v1) * (-dfsdr6*r6xij/r6)
      dv(8) = dv(8) + (v4 + v2 - v1) * (-dfsdr6*r6yij/r6)
      dv(9) = dv(9) + (v4 + v2 - v1) * (-dfsdr6*r6zij/r6)
      
      dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr4*r4xij/r4)
      dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr4*r4yij/r4)
      dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr4*r4zij/r4)
      
      
      dv(1) = dv(1) + fs*( dv2dr3*r3xij/r3 + dv2dr5*r5xij/r5)
      dv(2) = dv(2) + fs*( dv2dr3*r3yij/r3 + dv2dr5*r5yij/r5)
      dv(3) = dv(3) + fs*( dv2dr3*r3zij/r3 + dv2dr5*r5zij/r5)
      
      dv(7) = dv(7) + fs*( dv2dr1*r1xij/r1 - dv2dr3*r3xij/r3)
      dv(8) = dv(8) + fs*( dv2dr1*r1yij/r1 - dv2dr3*r3yij/r3)
      dv(9) = dv(9) + fs*( dv2dr1*r1zij/r1 - dv2dr3*r3zij/r3)
      
      dv(10)= dv(10)+ fs*(-dv2dr1*r1xij/r1 - dv2dr5*r5xij/r5)
      dv(11)= dv(11)+ fs*(-dv2dr1*r1yij/r1 - dv2dr5*r5yij/r5)
      dv(12)= dv(12)+ fs*(-dv2dr1*r1zij/r1 - dv2dr5*r5zij/r5)
      
      
      
      
      !'H2O2 -> H(1) + HO2'
      
      else if ((r2.ge.in_bound).and.(r3.ge.in_bound).and. &
      (r5.ge.in_bound)) then
      !            !print*,'exit9'
      
      if(group(i,9).eq.0) group(i,9)=310
      
      
      
      if (group(i,9).ne.310) then
      
      !print*,'exchange in complex into 310 from',group(i,9),tempo
      group(i,9)=310
      
      end if
      
      ! HO2 function
      V2=vho2_f4(r1,r4,r6)
      call dervho2_f4(r1,r4,r6,dv2dr1,dv2dr4,dv2dr6)
      
      
      
      
      ! Switch function
      
      
      fs=f3(rin(2),rin(3),rin(5))
      call df3(rin(2),rin(3),rin(5),dfsdr2,dfsdr3,dfsdr5)
      
      
      ! for the diatomic potencials
      
      call pothoq_f17(rin(3),Vdiat,dVdiat)
      
      
      dv(1) = dv(1) + fs * dVdiat*r3xij/r3
      dv(2) = dv(2) + fs * dVdiat*r3yij/r3
      dv(3) = dv(3) + fs * dVdiat*r3zij/r3
      
      dv(7) = dv(7) - fs * dVdiat*r3xij/r3
      dv(8) = dv(8) - fs * dVdiat*r3yij/r3
      dv(9) = dv(9) - fs * dVdiat*r3zij/r3
      
      
      V4 = V4 + Vdiat
      
      call pothoq_f17(rin(5),Vdiat,dVdiat)
      
      
      dv(1) = dv(1) + fs * dVdiat*r5xij/r5
      dv(2) = dv(2) + fs * dVdiat*r5yij/r5
      dv(3) = dv(3) + fs * dVdiat*r5zij/r5
      
      dv(10) = dv(10) - fs * dVdiat*r5xij/r5
      dv(11) = dv(11) - fs * dVdiat*r5yij/r5
      dv(12) = dv(12) - fs * dVdiat*r5zij/r5
      
      
      V4 = V4 + Vdiat
      
      
      call pothht_f14(r2,Vdiat,dVdiat)
      
      
      
      dv(1) = dv(1) + fs * dVdiat*r2xij/r2
      dv(2) = dv(2) + fs * dVdiat*r2yij/r2
      dv(3) = dv(3) + fs * dVdiat*r2zij/r2
      
      dv(4) = dv(4) - fs * dVdiat*r2xij/r2
      dv(5) = dv(5) - fs * dVdiat*r2yij/r2
      dv(6) = dv(6) - fs * dVdiat*r2zij/r2
      
      
      V4 = V4 + Vdiat
      
      V=V1*(1-fs)+(V2+V4)*fs
      !        !print*,rin
      !        !print*,'V=',V,' V1=',V1,' V2=',V2
      
      dv=dv + dv1*(1-fs)
      
      
      dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr2*r2xij/r2 + dfsdr3*r3xij/r3 + dfsdr5*r5xij/r5)
      dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr2*r2yij/r2 + dfsdr3*r3yij/r3 + dfsdr5*r5yij/r5)
      dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr2*r2zij/r2 + dfsdr3*r3zij/r3 + dfsdr5*r5zij/r5)
      
      dv(4) = dv(4) + (v4 + v2 - v1) * (-dfsdr2*r2xij/r2)
      dv(5) = dv(5) + (v4 + v2 - v1) * (-dfsdr2*r2yij/r2)
      dv(6) = dv(6) + (v4 + v2 - v1) * (-dfsdr2*r2zij/r2)
      
      dv(7) = dv(7) + (v4 + v2 - v1) * (-dfsdr3*r3xij/r3)
      dv(8) = dv(8) + (v4 + v2 - v1) * (-dfsdr3*r3yij/r3)
      dv(9) = dv(9) + (v4 + v2 - v1) * (-dfsdr3*r3zij/r3)
      
      dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr5*r5xij/r5)
      dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr5*r5yij/r5)
      dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr5*r5zij/r5)
      
      
      dv(4) = dv(4) + fs*( dv2dr4*r4xij/r4 + dv2dr6*r6xij/r6)
      dv(5) = dv(5) + fs*( dv2dr4*r4yij/r4 + dv2dr6*r6yij/r6)
      dv(6) = dv(6) + fs*( dv2dr4*r4zij/r4 + dv2dr6*r6zij/r6)
      
      dv(7) = dv(7) + fs*( dv2dr1*r1xij/r1 - dv2dr6*r6xij/r6)
      dv(8) = dv(8) + fs*( dv2dr1*r1yij/r1 - dv2dr6*r6yij/r6)
      dv(9) = dv(9) + fs*( dv2dr1*r1zij/r1 - dv2dr6*r6zij/r6)
      
      dv(10)= dv(10)+ fs*(-dv2dr1*r1xij/r1 - dv2dr4*r4xij/r4)
      dv(11)= dv(11)+ fs*(-dv2dr1*r1yij/r1 - dv2dr4*r4yij/r4)
      dv(12)= dv(12)+ fs*(-dv2dr1*r1zij/r1 - dv2dr4*r4zij/r4)
      
      
      
      else if ((r1.ge.in_bound).and.(r2.ge.in_bound).and. &
      (r5.ge.in_bound).and.(r6.ge.in_bound)) then
      !              !print*,'exit1',tempo
      if(group(i,9).eq.0) group(i,9)=320
      
      ! 'H2O2 -> HO(1-3) + HO(2-4)'
      if(group(i,9).ne.320) then
      !print*,'exchange in molecular complex 2'
      group(i,9)=320
      end if
      
      ! OH function
      call potohp_f3(r3,V2,dV2)
      ! OH function
      call potohp_f3(r4,V3,dV3)
      !            !print*,'r3,V2,dv2'
      !             !print*,r3,V2,dv2
      !              !print*,'r4,V3,dv3'
      !               !print*,r4,V3,dv3
      
      
      ! Switch function
      fs=f4(rin(1),rin(2),rin(5),rin(6))
      call df4(rin(1),rin(2),rin(5),rin(6),dfsdr1,dfsdr2,dfsdr5,dfsdr6)
      !               !print*,'rin(1),rin(2),rin(5),rin(6)'
      !               !print*,rin(1),rin(2),rin(5),rin(6)
      !               !print*,'dfsdr1,dfsdr2,dfsdr5,dfsdr6'
      !               !print*,dfsdr1,dfsdr2,dfsdr5,dfsdr6
      
      
      ! for the diatomic potencials
      
      call pothoq_f17(rin(5),Vdiat,dVdiat)
      
      !            !print*,'rin(5),Vdiat,dvdiat'
      !            !print*,rin(5),Vdiat,dvdiat
      !
      
      
      dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
      dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
      dv(3) = dv(3)  + fs * dVdiat*r5zij/r5
      
      dv(10)= dv(10) - fs * dVdiat*r5xij/r5
      dv(11)= dv(11) - fs * dVdiat*r5yij/r5
      dv(12)= dv(12) - fs * dVdiat*r5zij/r5
      
      
      V4 = V4 + Vdiat
      
      
      
      call pothht_f14(rin(2),Vdiat,dVdiat)
      
      
      !            !print*,'rin(2),Vdiat,dvdiat'
      !            !print*,rin(2),Vdiat,dvdiat
      
      dv(1) = dv(1) + fs * dVdiat*r2xij/r2
      dv(2) = dv(2) + fs * dVdiat*r2yij/r2
      dv(3) = dv(3) + fs * dVdiat*r2zij/r2
      
      dv(4) = dv(4) - fs * dVdiat*r2xij/r2
      dv(5) = dv(5) - fs * dVdiat*r2yij/r2
      dv(6) = dv(6) - fs * dVdiat*r2zij/r2
      
      
      V4 = V4 + Vdiat
      
      
      
      call pothoq_f17(rin(1),Vdiat,dVdiat)
      
      !           !print*,'rin(1),Vdiat,dvdiat'
      !           !print*,rin(1),Vdiat,dvdiat
      
      
      
      dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
      dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
      dv(9) = dv(9)  + fs * dVdiat*r1zij/r1
      
      dv(10)= dv(10) - fs * dVdiat*r1xij/r1
      dv(11)= dv(11) - fs * dVdiat*r1yij/r1
      dv(12)= dv(12) - fs * dVdiat*r1zij/r1
      
      
      V4 = V4 + Vdiat
      
      
      
      call pothoq_f17(rin(6),Vdiat,dVdiat)
      
      
      !print*,'rin(6),Vdiat,dvdiat'
      !print*,rin(6),Vdiat,dvdiat
      
      dv(4) = dv(4) + fs * dVdiat*r6xij/r6
      dv(5) = dv(5) + fs * dVdiat*r6yij/r6
      dv(6) = dv(6) + fs * dVdiat*r6zij/r6
      
      dv(7) = dv(7) - fs * dVdiat*r6xij/r6
      dv(8) = dv(8) - fs * dVdiat*r6yij/r6
      dv(9) = dv(9) - fs * dVdiat*r6zij/r6
      
      
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
      !print*,'exit5'
      
      
      
      ! 'H2O2 -> HO(1-4) + HO(2-3)'
      if(group(i,9).eq.0) group(i,9)=330
      
      if(group(i,9).ne.330) then
      print*,'exchange in molecular complex: group',group(i,9),tempo
      group(i,9)=330
      end if
      
      
      ! OH function
      call potohp_f3(r5,V3,dV3)
      ! OH function
      call potohp_f3(r6,V2,dV2)
      
      
      
      ! Switch function
      fs=f4(rin(1),rin(2),rin(3),rin(4))
      call df4(rin(1),rin(2),rin(3),rin(4),dfsdr1,dfsdr2,dfsdr3,dfsdr4)
      
      
      
      ! for the diatomic potencials
      
      call pothht_f14(rin(2),Vdiat,dVdiat)
      
      
      !            !print*,'rin(2),Vdiat,dvdiat'
      !            !print*,rin(2),Vdiat,dvdiat
      
      dv(1) = dv(1) + fs * dVdiat*r2xij/r2
      dv(2) = dv(2) + fs * dVdiat*r2yij/r2
      dv(3) = dv(3) + fs * dVdiat*r2zij/r2
      
      dv(4) = dv(4) - fs * dVdiat*r2xij/r2
      dv(5) = dv(5) - fs * dVdiat*r2yij/r2
      dv(6) = dv(6) - fs * dVdiat*r2zij/r2
      
      
      V4 = V4 + Vdiat
      
      
      call pothoq_f17(rin(3),Vdiat,dVdiat)
      
      
      !            !print*,'rin(3),Vdiat,dvdiat'
      !            !print*,rin(3),Vdiat,dvdiat
      
      dv(1) = dv(1) + fs * dVdiat*r3xij/r3
      dv(2) = dv(2) + fs * dVdiat*r3yij/r3
      dv(3) = dv(3) + fs * dVdiat*r3zij/r3
      
      dv(7) = dv(7) - fs * dVdiat*r3xij/r3
      dv(8) = dv(8) - fs * dVdiat*r3yij/r3
      dv(9) = dv(9) - fs * dVdiat*r3zij/r3
      
      
      V4 = V4 + Vdiat
      
      
      
      call pothoq_f17(rin(4),Vdiat,dVdiat)
      
      !            !print*,'rin(4),Vdiat,dvdiat'
      !            !print*,rin(4),Vdiat,dvdiat
      
      
      
      dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
      dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
      dv(6) = dv(6)  + fs * dVdiat*r4zij/r4
      
      dv(10)= dv(10) - fs * dVdiat*r4xij/r4
      dv(11)= dv(11) - fs * dVdiat*r4yij/r4
      dv(12)= dv(12) - fs * dVdiat*r4zij/r4
      
      
      V4 = V4 + Vdiat
      
      
      
      call poto2q_f18(rin(1),Vdiat,dVdiat)
      !           !print*,'rin(1),Vdiat,dvdiat'
      !           !print*,rin(1),Vdiat,dvdiat
      
      
      
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
      !              !print*,'exit4'
      
      
      
      !'H2O2 -> H2(1-2) + O2(3-4)'
      !!print*,'tempo',tempo
      if(group(i,9).eq.0) group(i,9)=340
      if(group(i,9).ne.340) then
      print*,'exchange in molecular complex: group',group(i,9),tempo
      group(i,9)=340
      end if
      
      
      ! H2 function
      call pothhs_f1(r2,V2,dV2)
      ! O2 function
      call potood_f22(r1,V3,dV3)
      
      !             !print*,'hop', r1,V3,dv3
      
      
      ! Switch function
      fs=f4(rin(3),rin(4),rin(5),rin(6))
      call df4(rin(3),rin(4),rin(5),rin(6),dfsdr3,dfsdr4,dfsdr5,dfsdr6)
      !        !print*,'dfsdr',dfsdr3,dfsdr4,dfsdr5,dfsdr6,fs
      
      
      ! for the diatomic potencials
      
      call pothoq_f17(rin(3),Vdiat,dVdiat)
      
      
      !            !print*,'rin(3),Vdiat,dvdiat'
      !            !print*,rin(3),Vdiat,dvdiat
      
      dv(1) = dv(1) + fs * dVdiat*r3xij/r3
      dv(2) = dv(2) + fs * dVdiat*r3yij/r3
      dv(3) = dv(3) + fs * dVdiat*r3zij/r3
      
      dv(7) = dv(7) - fs * dVdiat*r3xij/r3
      dv(8) = dv(8) - fs * dVdiat*r3yij/r3
      dv(9) = dv(9) - fs * dVdiat*r3zij/r3
      V4 = V4 + Vdiat
      
      
      call pothoq_f17(rin(5),Vdiat,dVdiat)
      
      !            !print*,'rin(5),Vdiat,dvdiat'
      !            !print*,rin(5),Vdiat,dvdiat
      
      
      
      dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
      dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
      dv(3) = dv(3)  + fs * dVdiat*r5zij/r5
      
      dv(10)= dv(10) - fs * dVdiat*r5xij/r5
      dv(11)= dv(11) - fs * dVdiat*r5yij/r5
      dv(12)= dv(12) - fs * dVdiat*r5zij/r5
      
      
      V4 = V4 + Vdiat
      
      
      
      
      
      
      call pothoq_f17(rin(6),Vdiat,dVdiat)
      
      
      !            !print*,'rin(6),Vdiat,dvdiat'
      !            !print*,rin(6),Vdiat,dvdiat
      
      dv(4) = dv(4) + fs * dVdiat*r6xij/r6
      dv(5) = dv(5) + fs * dVdiat*r6yij/r6
      dv(6) = dv(6) + fs * dVdiat*r6zij/r6
      
      dv(7) = dv(7) - fs * dVdiat*r6xij/r6
      dv(8) = dv(8) - fs * dVdiat*r6yij/r6
      dv(9) = dv(9) - fs * dVdiat*r6zij/r6
      
      
      V4 = V4 + Vdiat
      
      
      
      call pothoq_f17(rin(4),Vdiat,dVdiat)
      
      !            !print*,'rin(4),Vdiat,dvdiat'
      !            !print*,rin(4),Vdiat,dvdiat
      
      
      
      dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
      dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
      dv(6) = dv(6)  + fs * dVdiat*r4zij/r4
      
      dv(10)= dv(10) - fs * dVdiat*r4xij/r4
      dv(11)= dv(11) - fs * dVdiat*r4yij/r4
      dv(12)= dv(12) - fs * dVdiat*r4zij/r4
      
      
      V4 = V4 + Vdiat
      
      !            !print*,'dv depois  diat',dv
      
      
      
      
      V=V1*(1-fs)+(V2+V3+V4)*fs
      
      dv=dv + dv1*(1-fs)
      !       !print*,'dv depois 0',dv
      dv(1)  =   dv(1)  + (v2+v3+v4-v1) * ( dfsdr3 *  r3xij / r3 + dfsdr5 * r5xij / r5)
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
      !            !print*,'dv depois  1',dv
      
      
      
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
      
      !            !print*,'dv depois  2',dv
      
      
      
      else!  end if
      !print*,'inside',tempo
      
      ! H2O2 potential
      
      
      
      v=v1
      
      dv  = dv1
      
      group(i,9)=0
      
      end if
      
      !print*,'v',V
      !print*,dv(1:3)
      !print*,dv(4:6)
      !print*,dv(7:9)
      !print*,dv(10:12)
      
      
      
      
      end subroutine vh2o2_switch
      
      
      !*************************************************
      ! PROGRAM TO COMPUTE H2O2 POTENCIAL ENERGY SURFACE
      !*************************************************
      
      subroutine h2o2sur_10(n,rin,f,df)
      implicit none
      
      real*8 rin(6),r(6),f,df(6),vh2o2_10
      integer dummy,n,icalc 
      
      !COMMON/SURFAC_10/ICALC
      dummy=n 
      
      icalc=2
      
      r=rin/5.2917721092d-11
      !!print*,'from h2o2'
      !!print*,r
      
      f= vh2o2_10(r(1),r(2),r(3),r(4),r(5),r(6))
      
      
      if (icalc==1) return
      call dvh2o2_10(r(1),r(2),r(3),r(4),r(5),r(6),df(1),df(2),df(3),df(4),df(5),df(6))
      
      f=f*4.3597482D-18
      df=df*4.3597482D-18/5.2917721092d-11
      
      
      
      
      
      return
      end
      
      
      FUNCTION Vh2o2_10(R1,R2,R3,R4,R5,R6)
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
      real*8 v1,VOO_10,VHHS_10,VOHP_10,VHO23c_10,v2,VOOD_10,VOHS_10,vh2o3c2_10
      real*8 v12c,v,vh2o2_10,VH2O24c_10
      real*8 vo1d,Vinterm,V3,vh2o3c1_10,vh2o3c12_10,vh2ot3c_10,VHHT_10
      real*8 r1,r2,r3,r4,r5,r6
      real*8 Vdip_10
      vo1d=7.1955d-2
      V1=VOO_10(R1)+VHHS_10(R2)+VOHS_10(R3)+VOHS_10(R4)+VOHS_10(R5)+VOHS_10(R6)+VHO23c_10(R1,R3,R5)+VHO23c_10(R1,R4,R6) &
      +vh2o3c1_10(R2,R3,R6)+vh2o3c1_10(R2,R4,R5)+vo1d+vo1d
      
      V2=VOO_10(R1)+VHHT_10(R2)+VOHP_10(R3)+VOHP_10(R4)+VOHP_10(R5)+VOHP_10(R6)+VHO23c_10(R1,R3,R5)+ &
      VHO23c_10(R1,R4,R6)+vh2o3c2_10(R2,R3,R6)+vh2o3c2_10(R2,R4,R5)
      
      v12c=vh2o3c12_10(R2,R3,R6)+vh2o3c12_10(R2,R4,R5)
      
      Vinterm=((V1+V2)-SQRT((V1-V2)**2+4.0D0*v12c**2))/2.0D0
      
      V3=VOOD_10(R1)+VHHS_10(R2)+VOHP_10(R3)+VOHP_10(R4)+VOHP_10(R5)+VOHP_10(R6)+vh2ot3c_10(R2,R3,R6) &
      +vh2ot3c_10(R2,R4,R5)+VHO23c_10(R1,R3,R5)+VHO23c_10(R1,R4,R6)
      
      If (Vinterm < V3) then
      V = Vinterm
      else
      V = V3
      end if
      !     V=min(Vinterm,V3)
      
      
      Vh2o2_10=V+VH2O24c_10(r1,r2,r3,r4,r5,r6)+Vdip_10(r1,r2,r3,r4,r5,r6)
      
      
      
      !write(30,*) 'v1,v2,v12c,vinterm,v3,v'
      !write(30,*)  v1,v2,v12c,vinterm,v3,v
      
      
      END
      
      
      
      !***************************************************
      ! PROGRAM TO COMPUTE THE DERIVATIVES OF THE H2O2 PES
      !***************************************************
      
      subroutine dvh2o2_10 (R1,R2,R3,R4,R5,R6,dg1,dg2,dg3,dg4,dg5,dg6)
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      !vo1d=7.1955d-2
      implicit none
      real*8 v1,VOO_10,VHHS_10,VOHP_10,VHO23c_10,vh2o3c1_10,vh2o3c2_10,v2,VOOD_10,VOHS_10,vh2ot3c_10,v12c,v,vh2o2_10,VH2O24c_10
      real*8 Vinterm,V3,vh2o3c12_10,VHHT_10,dho2r3a,dho2r4b
      real*8 r1,r2,r3,r4,r5,r6,dr1,dr2,dr3,dr4,dr5,dr6,dv1dr1,dv1dr2,dv1dr3,dv1dr4,dv1dr5,dv1dr6
      real*8 dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6,dv12dr1,dv12dr2,dv12dr3,dv12dr4,dv12dr5,dv12dr6
      real*8 dvdr1,dvdr2,dvdr3,dvdr4,dvdr5,dvdr6,dho2r1a,dho2r1b
      real*8 dho2r5a,dho2r6b,dr2ta,dr3ta,dr6ta,dr2tb,dr4tb,dr5tb,dVOHP_10,dVOO_10,dVHHS_10,dVOHS_10,dVOOD_10
      real*8 dh2o2dr1,dh2o2dr2,dh2o2dr3,dh2o2dr4,dh2o2dr5,dh2o2dr6,dv4cdr1,dv4cdr2,dv4cdr3
      real*8 dv4cdr4,dv4cdr5,dv4cdr6,dg1,dg2,dg3,dg4,dg5,dg6,F1,F2,F3,F4,F5,F6
      real*8 dv3dr1,dv3dr2,dv3dr3,dv3dr4,dv3dr5,dv3dr6,dVHHT_10,dvintermdv1,dvintermdv2,dvintermdv12
      real*8 dr2tb1a,dr3tb1a,dr6tb1a,dr2tb1b,dr4tb1b,dr5tb1b,dr2rl1a,dr3rl1a,dr6rl1a,dr2rl1b,dr4rl1b,dr5rl1b
      real*8 dr2tb2a,dr3tb2a,dr6tb2a,dr2tb2b,dr4tb2b,dr5tb2b,dr2rl2a,dr3rl2a,dr6rl2a,dr2rl2b,dr4rl2b,dr5rl2b
      real*8 dr2tb12a,dr3tb12a,dr6tb12a,dr2tb12b,dr4tb12b,dr5tb12b
      real*8 dvdipdr1,dvdipdr2,dvdipdr3,dvdipdr4,dvdipdr5,dvdipdr6
      ! JB acrescentou
      real*8 R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      
      !$OMP THREADPRIVATE(/DIATDI_10/)
      
      V1=VOO_10(R1)+VHHS_10(R2)+VOHS_10(R3)+VOHS_10(R4)+VOHS_10(R5)+VOHS_10(R6)+VHO23c_10(R1,R3,R5)+VHO23c_10(R1,R4,R6) &
      +vh2o3c1_10(R2,R3,R6)+vh2o3c1_10(R2,R4,R5)+vo1d+vo1d
      
      V2=VOO_10(R1)+VHHT_10(R2)+VOHP_10(R3)+VOHP_10(R4)+VOHP_10(R5)+VOHP_10(R6)+VHO23c_10(R1,R3,R5)+VHO23c_10(R1,R4,R6) &
      +vh2o3c2_10(R2,R3,R6)+vh2o3c2_10(R2,R4,R5)
      
      v12c=vh2o3c12_10(R2,R3,R6)+vh2o3c12_10(R2,R4,R5)
      
      Vinterm=((V1+V2)-SQRT((V1-V2)**2+4.0D0*v12c**2))/2.0D0
      
      V3=VOOD_10(R1)+VHHS_10(R2)+VOHP_10(R3)+VOHP_10(R4)+VOHP_10(R5)+VOHP_10(R6)+vh2ot3c_10(R2,R3,R6)+vh2ot3c_10(R2,R4,R5) &
      + VHO23c_10(R1,R3,R5)+VHO23c_10(R1,R4,R6)
      
      
      If (Vinterm < V3) then
      V = Vinterm
      else
      V = V3
      end if
      
      !     V=min(Vinterm,V3)
      
      Vh2o2_10=V+VH2O24c_10(r1,r2,r3,r4,r5,r6)
      
      
      !derivadas de vinterm em ordem a v1,v2,v12
      
      dvintermdv1= 0.5d0-0.5d0*(v1-v2)/sqrt((v1-v2)**2 + 4.0d0*v12c**2)
      
      dvintermdv2= 0.5d0-0.5d0*(v2-v1)/sqrt((v1-v2)**2 + 4.0d0*v12c**2)
      
      dvintermdv12= -0.5d0*4.0d0*v12c/sqrt((v1-v2)**2 + 4.0d0*v12c**2)
      
      
      call derVHO23c_10(r1,r3,r5,dr1,dr3,dr5)
      dho2r1a=dr1
      dho2r3a=dr3
      dho2r5a=dr5
      call derVHO23c_10(r1,r4,r6,dr1,dr4,dr6)
      dho2r1b=dr1
      dho2r4b=dr4
      dho2r6b=dr6
      
      call DERVH2OT3c_10(r2,r3,r6,dr2,dr3,dr6)
      dr2ta=dr2
      dr3ta=dr3
      dr6ta=dr6
      call DERVH2OT3c_10(r2,r4,r5,dr2,dr4,dr5)
      dr2tb=dr2
      dr4tb=dr4
      dr5tb=dr5
      
      !derivadas de v1 em ordem aos rs
      
      !V1=VOO_10(R1)+VHHS_10(R2)+VOHS_10(R3)+VOHS_10(R4)+VOHS_10(R5)+VOHS_10(R6)+VHO23c_10(R1,R3,R5)+VHO23c_10(R1,R4,R6) &
      !      +vh2o3c1_10(R2,R3,R6)+vh2o3c1_10(R2,R4,R5)+vo1d+vo1d
      
      
      call dTHREBR1_10(r2,r3,r6,dr2,dr3,dr6) !superficieh2o.f
      dr2tb1a=dr2
      dr3tb1a=dr3
      dr6tb1a=dr6
      
      call dRL1_10(r2,r3,r6,dr2,dr3,dr6) !superficieh2o.f
      dr2rl1a=dr2
      dr3rl1a=dr3
      dr6rl1a=dr6
      
      call dTHREBR1_10(r2,r4,r5,dr2,dr4,dr5) !superficieh2o.f
      dr2tb1b=dr2
      dr4tb1b=dr4
      dr5tb1b=dr5
      
      call dRL1_10(r2,r4,r5,dr2,dr4,dr5) !superficieh2o.f
      dr2rl1b=dr2
      dr4rl1b=dr4
      dr5rl1b=dr5
      
      
      
      dv1dr1=dVOO_10(r1) + dho2r1a + dho2r1b
      dv1dr2=dVHHS_10(r2) + dr2tb1a + dr2rl1a + dr2tb1b + dr2rl1b
      dv1dr3=dVOHS_10(r3)+ dho2r3a + dr3tb1a + dr3rl1a
      dv1dr4=dVOHS_10(r4) + dho2r4b + dr4tb1b + dr4rl1b
      dv1dr5=dVOHS_10(r5) + dho2r5a + dr5tb1b + dr5rl1b
      dv1dr6=dVOHS_10(r6) + dho2r6b + dr6tb1a + dr6rl1a
      
      
      
      !derivadas de v2 em ordem aos rs
      
      !V2=VOO_10(R1)+VHHT_10(R2)+VOHP_10(R3)+VOHP_10(R4)+VOHP_10(R5)+VOHP_10(R6)+VHO23c_10(R1,R3,R5)+VHO23c_10(R1,R4,R6) &
      !      +vh2o3c2_10(R2,R3,R6)+vh2o3c2_10(R2,R4,R5)
      
      call dTHREBR2_10(r2,r3,r6,dr2,dr3,dr6)
      dr2tb2a=dr2
      dr3tb2a=dr3
      dr6tb2a=dr6
      
      call drl2_10(r2,r3,r6,dr2,dr3,dr6)
      dr2rl2a=dr2
      dr3rl2a=dr3
      dr6rl2a=dr6
      
      call dTHREBR2_10(r2,r4,r5,dr2,dr4,dr5)
      dr2tb2b=dr2
      dr4tb2b=dr4
      dr5tb2b=dr5
      
      call drl2_10(r2,r4,r5,dr2,dr4,dr5)
      dr2rl2b=dr2
      dr4rl2b=dr4
      dr5rl2b=dr5
      
      
      dv2dr1=dVOO_10(r1) + dho2r1a + dho2r1b
      dv2dr2=dVHHT_10(r2) + dr2tb2a + dr2rl2a + dr2tb2b + dr2rl2b
      dv2dr3=dVOHP_10(r3) + dho2r3a + dr3tb2a + dr3rl2a
      dv2dr4=dVOHP_10(r4) + dho2r4b + dr4tb2b + dr4rl2b
      dv2dr5=dVOHP_10(r5) + dho2r5a + dr5tb2b + dr5rl2b
      dv2dr6=dVOHP_10(r6) + dho2r6b + dr6tb2a + dr6rl2a
      
      
      !derivadas de v12c em ordem aos rs
      
      !v12c=vh2o3c12_10(R2,R3,R6)+vh2o3c12_10(R2,R4,R5)
      
      call DTHREB12_10(r2,r3,r6,dr2,dr3,dr6)
      dr2tb12a=dr2
      dr3tb12a=dr3
      dr6tb12a=dr6
      
      call DTHREB12_10(r2,r4,r5,dr2,dr4,dr5)
      dr2tb12b=dr2
      dr4tb12b=dr4
      dr5tb12b=dr5
      
      
      dv12dr1=0.0D0
      dv12dr2=dr2tb12a + dr2tb12b
      dv12dr3=dr3tb12a
      dv12dr4=dr4tb12b
      dv12dr5=dr5tb12b
      dv12dr6=dr6tb12a
      
      !derivadas de v3 em ordem aos rs
      
      !V3=VOOD_10(R1)+VHHS_10(R2)+VOHP_10(R3)+VOHP_10(R4)+VOHP_10(R5)+VOHP_10(R6)+vh2ot3c_10(R2,R3,R6)+vh2ot3c_10(R2,R4,R5) &
      !   + VHO23c_10(R1,R3,R5)+VHO23c_10(R1,R4,R6)
      
      
      dv3dr1=dVOOD_10(r1)+ dho2r1a + dho2r1b
      dv3dr2=dVHHS_10(r2) + dr2ta + dr2tb
      dv3dr3=dVOHP_10(r3) + dr3ta + dho2r3a
      dv3dr4=dVOHP_10(r4) + dr4tb + dho2r4b
      dv3dr5=dVOHP_10(r5) + dr5tb + dho2r5a
      dv3dr6=dVOHP_10(r6) + dr6ta + dho2r6b
      
      
      !derivadas de v em ordem aos rs quando v=vinterm
      
      if (vinterm.lt.v3) then
      
      dvdr1= dvintermdv1*dv1dr1 + dvintermdv2*dv2dr1 + dvintermdv12*dv12dr1
      
      dvdr2= dvintermdv1*dv1dr2 + dvintermdv2*dv2dr2 + dvintermdv12*dv12dr2
      
      dvdr3= dvintermdv1*dv1dr3+ dvintermdv2*dv2dr3 + dvintermdv12*dv12dr3
      
      dvdr4= dvintermdv1*dv1dr4 + dvintermdv2*dv2dr4 + dvintermdv12*dv12dr4
      
      dvdr5= dvintermdv1*dv1dr5 + dvintermdv2*dv2dr5 + dvintermdv12*dv12dr5
      
      dvdr6= dvintermdv1*dv1dr6 + dvintermdv2*dv2dr6 + dvintermdv12*dv12dr6
      
      else
      
      dvdr1= dv3dr1
      
      dvdr2= dv3dr2
      
      dvdr3= dv3dr3
      
      dvdr4= dv3dr4
      
      dvdr5= dv3dr5
      
      dvdr6= dv3dr6
      
      end if
      
      
      !derivadas de vh2o2
      
      call dVH2O24c_10(R1,R2,R3,R4,R5,R6,F1,F2,F3,F4,F5,F6)
      
      dv4cdr1= F1
      dv4cdr2= F2
      dv4cdr3= F3
      dv4cdr4= F4
      dv4cdr5= F5
      dv4cdr6= F6
      
      !derivadas de Vdip
      
      call dVdip_10(r1,r2,r3,r4,r5,r6,dr1,dr2,dr3,dr4,dr5,dr6)
      
      dvdipdr1=dr1
      dvdipdr2=dr2
      dvdipdr3=dr3
      dvdipdr4=dr4
      dvdipdr5=dr5
      dvdipdr6=dr6
      
      !derivada de h2o2
      
      dh2o2dr1= dvdr1 + dv4cdr1 + dvdipdr1
      dh2o2dr2= dvdr2 + dv4cdr2 + dvdipdr2
      dh2o2dr3= dvdr3 + dv4cdr3 + dvdipdr3
      dh2o2dr4= dvdr4 + dv4cdr4 + dvdipdr4
      dh2o2dr5= dvdr5 + dv4cdr5 + dvdipdr5
      dh2o2dr6= dvdr6 + dv4cdr6 + dvdipdr6
      
      
      
      dg1= dh2o2dr1
      dg2= dh2o2dr2
      dg3= dh2o2dr3
      dg4= dh2o2dr4
      dg5= dh2o2dr5
      dg6= dh2o2dr6
      
      
      END
      
      
      
      
      !**********************************
      FUNCTION VH2O24c_10(R1,R2,R3,R4,R5,R6)
      !**********************************
      use param4c_10
      
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
      real*8 rho1eq,rho2eq2,rho3eq2,rho4eq2,rho234eq,s1,s2,rho1,rho22,rho32,rho42,rho234,dec
      real*8 pol,VH2O24c_10,r1,r2,r3,r4,r5,r6
      real*8 rho1min1,rho2min12,rho3min12,rho4min12,rho234min1
      real*8 rho1ts2,rho2ts22,rho3ts22,rho4ts22,rho234ts2
      real*8 rho1ts7m,rho2ts7m2,rho3ts7m2,rho4ts7m2,rho234ts7m
      real*8 rho1ts7,rho2ts72,rho3ts72,rho4ts72,rho234ts7
      real*8 rho1ts9,rho2ts92,rho3ts92,rho4ts92,rho234ts9
      real*8 rho1ts13,rho2ts132,rho3ts132,rho4ts132,rho234ts13
      real*8 v4c1,v4c2,v4c3,v4c4,v4c5,v4c6,v4c7
      real*8 xc(6),CV(6),expoente
      integer ncoord,j,k
      
      parameter (ncoord=6)
      
      
      !geometria de referencia h2o2
      
      rho1eq=r3eq+r4eq+r5eq+r6eq
      rho2eq2=(r3eq+r4eq-r5eq-r6eq)**2
      rho3eq2=(r3eq-r4eq+r5eq-r6eq)**2
      rho4eq2=(r3eq-r4eq-r5eq+r6eq)**2
      rho234eq=(r3eq+r4eq-r5eq-r6eq)*(r3eq-r4eq+r5eq-r6eq)*(r3eq-r4eq-r5eq+r6eq)
      
      s1=r1-r1eq
      s2=r2-r2eq
      rho1=r3+r4+r5+r6-rho1eq
      rho22=(r3+r4-r5-r6)**2-rho2eq2
      rho32=(r3-r4+r5-r6)**2-rho3eq2
      rho42=(r3-r4-r5+r6)**2-rho4eq2
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234eq
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp1(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp1(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      deceq=dec
      
      
      
      pol=c(1) + c(2)*s1 + c(3)*s2 + c(4)*rho1 + c(5)*s1**2 + c(6)*s1*s2 + c(7)*s1*rho1 + c(8)*s2**2 + c(9)*s2*rho1 +&
      c(10)*rho1**2 + c(11)*rho22 + c(12)*rho32 + c(13)*rho42 + c(14)*s1**3 + c(15)*s1**2*s2 + c(16)*s1**2*rho1 +&
      c(17)*s1*s2**2 + c(18)*s1*s2*rho1 + c(19)*s1*rho1**2 + c(20)*s1*rho22 + c(21)*s1*rho32 + c(22)*s1*rho42 +&
      c(23)*s2**3 + c(24)*s2**2*rho1 + c(25)*s2*rho1**2 + c(26)*s2*rho22 + c(27)*s2*rho32 + c(28)*s2*rho42 +&
      c(29)*rho1**3 + c(30)*rho1*rho22 + c(31)*rho1*rho32 + c(32)*rho1*rho42 + c(33)*rho234 +&
      c(34)*s1**4 + c(35)*s1**3*s2 + c(36)*s1**3*rho1 + c(37)*s1**2*s2**2 + c(38)*s1**2*s2*rho1 +&
      c(39)*s1**2*rho1**2 + c(40)*s1**2*rho22 + c(41)*s1**2*rho32 + c(42)*s1**2*rho42 + c(43)*s1*s2**3 +&
      c(44)*s1*s2**2*rho1 + c(45)*s1*s2*rho1**2 + c(46)*s1*s2*rho22 + c(47)*s1*s2*rho32 + c(48)*s1*s2*rho42 +&
      c(49)*s1*rho1**3 + c(50)*s1*rho1*rho22 + c(51)*s1*rho1*rho32 + c(52)*s1*rho1*rho42 + c(53)*s1*rho234 +&
      c(54)*s2**4 + c(55)*s2**3*rho1 + c(56)*s2**2*rho1**2 + c(57)*s2**2*rho22 + c(58)*s2**2*rho32 +&
      c(59)*s2**2*rho42 + c(60)*s2*rho1**3 + c(61)*s2*rho1*rho22 + c(62)*s2*rho1*rho32 + c(63)*s2*rho1*rho42 +&
      c(64)*s2*rho234 + c(65)*rho1**4 + c(66)*rho1**2*rho22 + c(67)*rho1**2*rho32 + c(68)*rho1**2*rho42 +&
      c(69)*rho1*rho234 + c(70)*rho22*rho32 + c(71)*rho22*rho42 + c(72)*rho32*rho42 + c(73)*rho22*rho22 +&
      c(74)*rho32*rho32 + c(75)*rho42*rho42
      
      V4c1=pol*dec
      
      
      
      !   com outra geometria de referencia (min1)
      
      rho1min1=r3min1+r4min1+r5min1+r6min1
      rho2min12=(r3min1+r4min1-r5min1-r6min1)**2
      rho3min12=(r3min1-r4min1+r5min1-r6min1)**2
      rho4min12=(r3min1-r4min1-r5min1+r6min1)**2
      rho234min1=(r3min1+r4min1-r5min1-r6min1)*(r3min1-r4min1+r5min1-r6min1)*(r3min1-r4min1-r5min1+r6min1)
      
      s1=r1-r1min1
      s2=r2-r2min1
      rho1=r3+r4+r5+r6-rho1min1
      rho22=(r3+r4-r5-r6)**2-rho2min12
      rho32=(r3-r4+r5-r6)**2-rho3min12
      rho42=(r3-r4-r5+r6)**2-rho4min12
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234min1
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp2(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp2(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      decmin1=dec
      
      
      
      pol=c(76) + c(77)*s1 + c(78)*s2 + c(79)*rho1 + c(80)*s1**2 + c(81)*s1*s2 + c(82)*s1*rho1 +&
      c(83)*s2**2+c(84)*s2*rho1 + c(85)*rho1**2 + c(86)*rho22 + c(87)*rho32 + c(88)*rho42 + c(89)*s1**3 +&
      c(90)*s1**2*s2 + c(91)*s1**2*rho1 + c(92)*s1*s2**2 + c(93)*s1*s2*rho1 + c(94)*s1*rho1**2 +&
      c(95)*s1*rho22 + c(96)*s1*rho32 + c(97)*s1*rho42 + c(98)*s2**3 + c(99)*s2**2*rho1 + c(100)*s2*rho1**2 +&
      c(101)*s2*rho22 + c(102)*s2*rho32 + c(103)*s2*rho42 + c(104)*rho1**3 + c(105)*rho1*rho22 +&
      c(106)*rho1*rho32 + c(107)*rho1*rho42 + c(108)*rho234 +&
      c(109)*s1**4 + c(110)*s1**3*s2 + c(111)*s1**3*rho1 + c(112)*s1**2*s2**2 + c(113)*s1**2*s2*rho1 +&
      c(114)*s1**2*rho1**2 + c(115)*s1**2*rho22 + c(116)*s1**2*rho32 + c(117)*s1**2*rho42 + c(118)*s1*s2**3 +&
      c(119)*s1*s2**2*rho1 + c(120)*s1*s2*rho1**2 + c(121)*s1*s2*rho22 + c(122)*s1*s2*rho32 + c(123)*s1*s2*rho42 +&
      c(124)*s1*rho1**3 + c(125)*s1*rho1*rho22 + c(126)*s1*rho1*rho32 + c(127)*s1*rho1*rho42 + c(128)*s1*rho234 +&
      c(129)*s2**4 + c(130)*s2**3*rho1 + c(131)*s2**2*rho1**2 + c(132)*s2**2*rho22 + c(133)*s2**2*rho32 +&
      c(134)*s2**2*rho42 + c(135)*s2*rho1**3 + c(136)*s2*rho1*rho22 + c(137)*s2*rho1*rho32 + c(138)*s2*rho1*rho42 +&
      c(139)*s2*rho234 + c(140)*rho1**4 + c(141)*rho1**2*rho22 + c(142)*rho1**2*rho32 + c(143)*rho1**2*rho42 +&
      c(144)*rho1*rho234 + c(145)*rho22*rho32 + c(146)*rho22*rho42 + c(147)*rho32*rho42 + c(148)*rho22*rho22 +&
      c(149)*rho32*rho32 + c(150)*rho42*rho42
      
      V4c2=pol*dec
      
      
      
      !   com outra geometria de referencia (ts7)
      
      rho1ts7=r3ts7+r4ts7+r5ts7+r6ts7
      rho2ts72=(r3ts7+r4ts7-r5ts7-r6ts7)**2
      rho3ts72=(r3ts7-r4ts7+r5ts7-r6ts7)**2
      rho4ts72=(r3ts7-r4ts7-r5ts7+r6ts7)**2
      rho234ts7=(r3ts7+r4ts7-r5ts7-r6ts7)*(r3ts7-r4ts7+r5ts7-r6ts7)*(r3ts7-r4ts7-r5ts7+r6ts7)
      
      s1=r1-r1ts7
      s2=r2-r2ts7
      rho1=r3+r4+r5+r6-rho1ts7
      rho22=(r3+r4-r5-r6)**2-rho2ts72
      rho32=(r3-r4+r5-r6)**2-rho3ts72
      rho42=(r3-r4-r5+r6)**2-rho4ts72
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts7
      
      !write(30,*)
      !write(30,'(a9,7f15.7)') 'rts2 coor', s1,s2,rho1,rho22,rho32,rho42,rho234
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp3(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp3(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      dects7=dec
      
      
      
      
      pol=c(151) + c(152)*s1 + c(153)*s2 + c(154)*rho1 + c(155)*s1**2 + c(156)*s1*s2 + c(157)*s1*rho1 +&
      c(158)*s2**2 +c(159)*s2*rho1 + c(160)*rho1**2 + c(161)*rho22 + c(162)*rho32 + c(163)*rho42 + c(164)*s1**3 +&
      c(165)*s1**2*s2 + c(166)*s1**2*rho1 + c(167)*s1*s2**2 + c(168)*s1*s2*rho1 + c(169)*s1*rho1**2 +&
      c(170)*s1*rho22 + c(171)*s1*rho32 + c(172)*s1*rho42 + c(173)*s2**3 + c(174)*s2**2*rho1 + c(175)*s2*rho1**2 +&
      c(176)*s2*rho22 + c(177)*s2*rho32 + c(178)*s2*rho42 + c(179)*rho1**3 + c(180)*rho1*rho22 +&
      c(181)*rho1*rho32 + c(182)*rho1*rho42 + c(183)*rho234 +&
      c(184)*s1**4 + c(185)*s1**3*s2 + c(186)*s1**3*rho1 + c(187)*s1**2*s2**2 + c(188)*s1**2*s2*rho1 +&
      c(189)*s1**2*rho1**2 + c(190)*s1**2*rho22 + c(191)*s1**2*rho32 + c(192)*s1**2*rho42 + c(193)*s1*s2**3 +&
      c(194)*s1*s2**2*rho1 + c(195)*s1*s2*rho1**2 + c(196)*s1*s2*rho22 + c(197)*s1*s2*rho32 + c(198)*s1*s2*rho42 +&
      c(199)*s1*rho1**3 + c(200)*s1*rho1*rho22 + c(201)*s1*rho1*rho32 + c(202)*s1*rho1*rho42 + c(203)*s1*rho234 +&
      c(204)*s2**4 + c(205)*s2**3*rho1 + c(206)*s2**2*rho1**2 + c(207)*s2**2*rho22 + c(208)*s2**2*rho32 +&
      c(209)*s2**2*rho42 + c(210)*s2*rho1**3 + c(211)*s2*rho1*rho22 + c(212)*s2*rho1*rho32 + c(213)*s2*rho1*rho42 +&
      c(214)*s2*rho234 + c(215)*rho1**4 + c(216)*rho1**2*rho22 + c(217)*rho1**2*rho32 + c(218)*rho1**2*rho42 +&
      c(219)*rho1*rho234 + c(220)*rho22*rho32 + c(221)*rho22*rho42 + c(222)*rho32*rho42 + c(223)*rho22*rho22 +&
      c(224)*rho32*rho32 + c(225)*rho42*rho42
      
      
      
      V4c3=pol*dec
      
      
      
      
      !   com outra geometria de referencia (ts2)
      
      rho1ts2=r3ts2+r4ts2+r5ts2+r6ts2
      rho2ts22=(r3ts2+r4ts2-r5ts2-r6ts2)**2
      rho3ts22=(r3ts2-r4ts2+r5ts2-r6ts2)**2
      rho4ts22=(r3ts2-r4ts2-r5ts2+r6ts2)**2
      rho234ts2=(r3ts2+r4ts2-r5ts2-r6ts2)*(r3ts2-r4ts2+r5ts2-r6ts2)*(r3ts2-r4ts2-r5ts2+r6ts2)
      
      s1=r1-r1ts2
      s2=r2-r2ts2
      rho1=r3+r4+r5+r6-rho1ts2
      rho22=(r3+r4-r5-r6)**2-rho2ts22
      rho32=(r3-r4+r5-r6)**2-rho3ts22
      rho42=(r3-r4-r5+r6)**2-rho4ts22
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts2
      
      !write(30,*)
      !write(30,'(a7,6f15.7)') 'rts7m ', s1,s2,rho1,rho22,rho32,rho42
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp4(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp4(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      dects2=dec
      
      
      
      
      pol=c(226) + c(227)*s1 + c(228)*s2 + c(229)*rho1 + c(230)*s1**2 + c(231)*s1*s2 + c(232)*s1*rho1 +&
      c(233)*s2**2 +c(234)*s2*rho1 + c(235)*rho1**2 + c(236)*rho22 + c(237)*rho32 + c(238)*rho42 + c(239)*s1**3 +&
      c(240)*s1**2*s2 + c(241)*s1**2*rho1 + c(242)*s1*s2**2 + c(243)*s1*s2*rho1 + c(244)*s1*rho1**2 +&
      c(245)*s1*rho22 + c(246)*s1*rho32 + c(247)*s1*rho42 + c(248)*s2**3 + c(249)*s2**2*rho1 + c(250)*s2*rho1**2 +&
      c(251)*s2*rho22 + c(252)*s2*rho32 + c(253)*s2*rho42 + c(254)*rho1**3 + c(255)*rho1*rho22 +&
      c(256)*rho1*rho32 + c(257)*rho1*rho42 + c(258)*rho234+&
      c(259)*s1**4 + c(260)*s1**3*s2 + c(261)*s1**3*rho1 + c(262)*s1**2*s2**2 + c(263)*s1**2*s2*rho1 +&
      c(264)*s1**2*rho1**2 + c(265)*s1**2*rho22 + c(266)*s1**2*rho32 + c(267)*s1**2*rho42 + c(268)*s1*s2**3 +&
      c(269)*s1*s2**2*rho1 + c(270)*s1*s2*rho1**2 + c(271)*s1*s2*rho22 + c(272)*s1*s2*rho32 + c(273)*s1*s2*rho42 +&
      c(274)*s1*rho1**3 + c(275)*s1*rho1*rho22 + c(276)*s1*rho1*rho32 + c(277)*s1*rho1*rho42 + c(278)*s1*rho234 +&
      c(279)*s2**4 + c(280)*s2**3*rho1 + c(281)*s2**2*rho1**2 + c(282)*s2**2*rho22 + c(283)*s2**2*rho32 +&
      c(284)*s2**2*rho42 + c(285)*s2*rho1**3 + c(286)*s2*rho1*rho22 + c(287)*s2*rho1*rho32 + c(288)*s2*rho1*rho42 +&
      c(289)*s2*rho234 + c(290)*rho1**4 + c(291)*rho1**2*rho22 + c(292)*rho1**2*rho32 + c(293)*rho1**2*rho42 +&
      c(294)*rho1*rho234 + c(295)*rho22*rho32 + c(296)*rho22*rho42 + c(297)*rho32*rho42 + c(298)*rho22*rho22 +&
      c(299)*rho32*rho32 + c(300)*rho42*rho42
      
      
      V4c4=pol*dec
      
      !write(30,*) 'v4c4=', v4c4, 'pol=', pol, 'dec=', dec
      
      !   com outra geometria de referencia (ts9)
      
      rho1ts9=r3ts9+r4ts9+r5ts9+r6ts9
      rho2ts92=(r3ts9+r4ts9-r5ts9-r6ts9)**2
      rho3ts92=(r3ts9-r4ts9+r5ts9-r6ts9)**2
      rho4ts92=(r3ts9-r4ts9-r5ts9+r6ts9)**2
      rho234ts9=(r3ts9+r4ts9-r5ts9-r6ts9)*(r3ts9-r4ts9+r5ts9-r6ts9)*(r3ts9-r4ts9-r5ts9+r6ts9)
      
      s1=r1-r1ts9
      s2=r2-r2ts9
      rho1=r3+r4+r5+r6-rho1ts9
      rho22=(r3+r4-r5-r6)**2-rho2ts92
      rho32=(r3-r4+r5-r6)**2-rho3ts92
      rho42=(r3-r4-r5+r6)**2-rho4ts92
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts9
      
      !write(30,*)
      !write(30,'(a7,6f15.7)') 'rts7ts13 ', s1,s2,rho1,rho22,rho32,rho42
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp5(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp5(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      dects9=dec
      
      
      !dec1=exp(-gamma1ts7ts13*s1**2)
      !dec2=exp(-gamma2ts7ts13*s2**2)
      !dec3=exp(-gamma3ts7ts13*rho1**2)
      !dec4=exp(-gamma4ts7ts13*rho22**2)
      !dec5=exp(-gamma5ts7ts13*rho32**2)
      !dec6=exp(-gamma6ts7ts13*rho42**2)
      !dec=dec1*dec2*dec3*dec4*dec5*dec6
      !dects7ts13=dec
      
      !write(30,'(a7,6e15.7)') 'rts7ts13 ', dec1,dec2,dec3,dec4,dec5,dec6
      
      
      pol=c(301) + c(302)*s1 + c(303)*s2 + c(304)*rho1 + c(305)*s1**2 + c(306)*s1*s2 + c(307)*s1*rho1 +&
      c(308)*s2**2 +c(309)*s2*rho1 + c(310)*rho1**2 + c(311)*rho22 + c(312)*rho32 + c(313)*rho42 + c(314)*s1**3 +&
      c(315)*s1**2*s2 + c(316)*s1**2*rho1 + c(317)*s1*s2**2 + c(318)*s1*s2*rho1 + c(319)*s1*rho1**2 +&
      c(320)*s1*rho22 + c(321)*s1*rho32 + c(322)*s1*rho42 + c(323)*s2**3 + c(324)*s2**2*rho1 + c(325)*s2*rho1**2 +&
      c(326)*s2*rho22 + c(327)*s2*rho32 + c(328)*s2*rho42 + c(329)*rho1**3 + c(330)*rho1*rho22 +&
      c(331)*rho1*rho32 + c(332)*rho1*rho42 + c(333)*rho234+&
      c(334)*s1**4 + c(335)*s1**3*s2 + c(336)*s1**3*rho1 + c(337)*s1**2*s2**2 + c(338)*s1**2*s2*rho1 +&
      c(339)*s1**2*rho1**2 + c(340)*s1**2*rho22 + c(341)*s1**2*rho32 + c(342)*s1**2*rho42 + c(343)*s1*s2**3 +&
      c(344)*s1*s2**2*rho1 + c(345)*s1*s2*rho1**2 + c(346)*s1*s2*rho22 + c(347)*s1*s2*rho32 + c(348)*s1*s2*rho42 +&
      c(349)*s1*rho1**3 + c(350)*s1*rho1*rho22 + c(351)*s1*rho1*rho32 + c(352)*s1*rho1*rho42 + c(353)*s1*rho234 +&
      c(354)*s2**4 + c(355)*s2**3*rho1 + c(356)*s2**2*rho1**2 + c(357)*s2**2*rho22 + c(358)*s2**2*rho32 +&
      c(359)*s2**2*rho42 + c(360)*s2*rho1**3 + c(361)*s2*rho1*rho22 + c(362)*s2*rho1*rho32 + c(363)*s2*rho1*rho42 +&
      c(364)*s2*rho234 + c(365)*rho1**4 + c(366)*rho1**2*rho22 + c(367)*rho1**2*rho32 + c(368)*rho1**2*rho42 +&
      c(369)*rho1*rho234 + c(370)*rho22*rho32 + c(371)*rho22*rho42 + c(372)*rho32*rho42 + c(373)*rho22*rho22 +&
      c(374)*rho32*rho32 + c(375)*rho42*rho42
      
      
      V4c5=pol*dec
      !write(30,*) 'v4c5=', v4c5, 'pol=', pol, 'dec=', dec
      
      
      !   com outra geometria de referencia (ts7 media)
      
      rho1ts7m=r3ts7m+r4ts7m+r5ts7m+r6ts7m
      rho2ts7m2=(r3ts7m+r4ts7m-r5ts7m-r6ts7m)**2
      rho3ts7m2=(r3ts7m-r4ts7m+r5ts7m-r6ts7m)**2
      rho4ts7m2=(r3ts7m-r4ts7m-r5ts7m+r6ts7m)**2
      rho234ts7m=(r3ts7m+r4ts7m-r5ts7m-r6ts7m)*(r3ts7m-r4ts7m+r5ts7m-r6ts7m)*(r3ts7m-r4ts7m-r5ts7m+r6ts7m)
      
      s1=r1-r1ts7m
      s2=r2-r2ts7m
      rho1=r3+r4+r5+r6-rho1ts7m
      rho22=(r3+r4-r5-r6)**2-rho2ts7m2
      rho32=(r3-r4+r5-r6)**2-rho3ts7m2
      rho42=(r3-r4-r5+r6)**2-rho4ts7m2
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts7m
      
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp6(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp6(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      dects7m=dec
      
      
      pol=c(376) + c(377)*s1 + c(378)*s2 + c(379)*rho1 + c(380)*s1**2 + c(381)*s1*s2 + c(382)*s1*rho1 +&
      c(383)*s2**2 +c(384)*s2*rho1 + c(385)*rho1**2 + c(386)*rho22 + c(387)*rho32 + c(388)*rho42 + c(389)*s1**3 +&
      c(390)*s1**2*s2 + c(391)*s1**2*rho1 + c(392)*s1*s2**2 + c(393)*s1*s2*rho1 + c(394)*s1*rho1**2 +&
      c(395)*s1*rho22 + c(396)*s1*rho32 + c(397)*s1*rho42 + c(398)*s2**3 + c(399)*s2**2*rho1 + c(400)*s2*rho1**2 +&
      c(401)*s2*rho22 + c(402)*s2*rho32 + c(403)*s2*rho42 + c(404)*rho1**3 + c(405)*rho1*rho22 +&
      c(406)*rho1*rho32 + c(407)*rho1*rho42 + c(408)*rho234+&
      c(409)*s1**4 + c(410)*s1**3*s2 + c(411)*s1**3*rho1 + c(412)*s1**2*s2**2 + c(413)*s1**2*s2*rho1 +&
      c(414)*s1**2*rho1**2 + c(415)*s1**2*rho22 + c(416)*s1**2*rho32 + c(417)*s1**2*rho42 + c(418)*s1*s2**3 +&
      c(419)*s1*s2**2*rho1 + c(420)*s1*s2*rho1**2 + c(421)*s1*s2*rho22 + c(422)*s1*s2*rho32 + c(423)*s1*s2*rho42 +&
      c(424)*s1*rho1**3 + c(425)*s1*rho1*rho22 + c(426)*s1*rho1*rho32 + c(427)*s1*rho1*rho42 + c(428)*s1*rho234 +&
      c(429)*s2**4 + c(430)*s2**3*rho1 + c(431)*s2**2*rho1**2 + c(432)*s2**2*rho22 + c(433)*s2**2*rho32 +&
      c(434)*s2**2*rho42 + c(435)*s2*rho1**3 + c(436)*s2*rho1*rho22 + c(437)*s2*rho1*rho32 + c(438)*s2*rho1*rho42 +&
      c(439)*s2*rho234 + c(440)*rho1**4 + c(441)*rho1**2*rho22 + c(442)*rho1**2*rho32 + c(443)*rho1**2*rho42 +&
      c(444)*rho1*rho234 + c(445)*rho22*rho32 + c(446)*rho22*rho42 + c(447)*rho32*rho42 + c(448)*rho22*rho22 +&
      c(449)*rho32*rho32 + c(450)*rho42*rho42
      
      V4c6=pol*dec
      
      
      !com outra geometria de referencia ts13
      
      rho1ts13=r3ts13+r4ts13+r5ts13+r6ts13
      rho2ts132=(r3ts13+r4ts13-r5ts13-r6ts13)**2
      rho3ts132=(r3ts13-r4ts13+r5ts13-r6ts13)**2
      rho4ts132=(r3ts13-r4ts13-r5ts13+r6ts13)**2
      rho234ts13=(r3ts13+r4ts13-r5ts13-r6ts13)*(r3ts13-r4ts13+r5ts13-r6ts13)*(r3ts13-r4ts13-r5ts13+r6ts13)
      
      s1=r1-r1ts13
      s2=r2-r2ts13
      rho1=r3+r4+r5+r6-rho1ts13
      rho22=(r3+r4-r5-r6)**2-rho2ts132
      rho32=(r3-r4+r5-r6)**2-rho3ts132
      rho42=(r3-r4-r5+r6)**2-rho4ts132
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts13
      
      
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp7(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp7(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      dects13=dec
      
      
      pol=c(451) + c(452)*s1 + c(453)*s2 + c(454)*rho1 + c(455)*s1**2 + c(456)*s1*s2 + c(457)*s1*rho1 +&
      c(458)*s2**2 +c(459)*s2*rho1 + c(460)*rho1**2 + c(461)*rho22 + c(462)*rho32 + c(463)*rho42 + c(464)*s1**3 +&
      c(465)*s1**2*s2 + c(466)*s1**2*rho1 + c(467)*s1*s2**2 + c(468)*s1*s2*rho1 + c(469)*s1*rho1**2 +&
      c(470)*s1*rho22 + c(471)*s1*rho32 + c(472)*s1*rho42 + c(473)*s2**3 + c(474)*s2**2*rho1 + c(475)*s2*rho1**2 +&
      c(476)*s2*rho22 + c(477)*s2*rho32 + c(478)*s2*rho42 + c(479)*rho1**3 + c(480)*rho1*rho22 +&
      c(481)*rho1*rho32 + c(482)*rho1*rho42 + c(483)*rho234+&
      c(484)*s1**4 + c(485)*s1**3*s2 + c(486)*s1**3*rho1 + c(487)*s1**2*s2**2 + c(488)*s1**2*s2*rho1 +&
      c(489)*s1**2*rho1**2 + c(490)*s1**2*rho22 + c(491)*s1**2*rho32 + c(492)*s1**2*rho42 + c(493)*s1*s2**3 +&
      c(494)*s1*s2**2*rho1 + c(495)*s1*s2*rho1**2 + c(496)*s1*s2*rho22 + c(497)*s1*s2*rho32 + c(498)*s1*s2*rho42 +&
      c(499)*s1*rho1**3 + c(500)*s1*rho1*rho22 + c(501)*s1*rho1*rho32 + c(502)*s1*rho1*rho42 + c(503)*s1*rho234 +&
      c(504)*s2**4 + c(505)*s2**3*rho1 + c(506)*s2**2*rho1**2 + c(507)*s2**2*rho22 + c(508)*s2**2*rho32 +&
      c(509)*s2**2*rho42 + c(510)*s2*rho1**3 + c(511)*s2*rho1*rho22 + c(512)*s2*rho1*rho32 + c(513)*s2*rho1*rho42 +&
      c(514)*s2*rho234 + c(515)*rho1**4 + c(516)*rho1**2*rho22 + c(517)*rho1**2*rho32 + c(518)*rho1**2*rho42 +&
      c(519)*rho1*rho234 + c(520)*rho22*rho32 + c(521)*rho22*rho42 + c(522)*rho32*rho42 + c(523)*rho22*rho22 +&
      c(524)*rho32*rho32 + c(525)*rho42*rho42
      
      V4c7=pol*dec
      
      
      
      VH2O24c_10=V4c1 + V4c2 + V4c3 + V4c4 + V4c5 + V4c6 + V4c7
      
      
      
      END
      
      
      !***********************************************************
      subroutine DVH2O24c_10(R1,R2,R3,R4,R5,R6,F1,F2,F3,F4,F5,F6)
      !***********************************************************
      
      use param4c_10
      
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
      real*8 rho1eq,rho2eq2,rho3eq2,rho4eq2,rho234eq,s1,s2,rho1,rho22,rho32,rho42,rho234,dec
      real*8 pol,r1,r2,r3,r4,r5,r6,rho1min1,rho2min12,rho3min12,rho4min12,rho234min1
      real*8 F1,F2,F3,F4,F5,F6,dpolds1,dpolds2,dpoldrho1,dpoldrho22,dpoldrho32,dpoldrho42,dpoldrho234
      real*8 ddecds1,ddecds2,ddecdrho1,ddecdrho22,ddecdrho32,ddecdrho42
      real*8 dv4c1ds1,dv4c1ds2,dv4c1drho1,dv4c1drho22,dv4c1drho32,dv4c1drho42
      real*8 dv4c1drho234,ds1dr1,ds2dr2,drho1dr3,drho22dr3,drho32dr3,drho42dr3,drho234dr3,drho1dr4,drho22dr4,drho32dr4
      real*8 drho42dr4,drho234dr4,drho1dr5,drho22dr5,drho32dr5,drho42dr5,drho234dr5,drho1dr6,drho22dr6,drho32dr6
      real*8 drho42dr6,drho234dr6,dv4c1dr1,dv4c1dr2,dv4c1dr3,dv4c1dr4,dv4c1dr5,dv4c1dr6
      real*8 dv4c2ds1,dv4c2ds2,dv4c2drho1,dv4c2drho22,dv4c2drho32,dv4c2drho42,dv4c2drho234,dv4c2dr1,dv4c2dr2
      real*8 dv4c2dr3,dv4c2dr4,dv4c2dr5,dv4c2dr6
      real*8 dv4c3ds1,dv4c3ds2,dv4c3drho1,dv4c3drho22,dv4c3drho32,dv4c3drho42,dv4c3drho234,dv4c3dr1,dv4c3dr2
      real*8 dv4c3dr3,dv4c3dr4,dv4c3dr5,dv4c3dr6,rho1ts7,rho2ts72,rho3ts72,rho4ts72,rho234ts7
      real*8 dv4c4ds1,dv4c4ds2,dv4c4drho1,dv4c4drho22,dv4c4drho32,dv4c4drho42,dv4c4drho234,dv4c4dr1,dv4c4dr2
      real*8 dv4c4dr3,dv4c4dr4,dv4c4dr5,dv4c4dr6,rho1ts2,rho2ts22,rho3ts22,rho4ts22,rho234ts2
      real*8 dv4c5ds1,dv4c5ds2,dv4c5drho1,dv4c5drho22,dv4c5drho32,dv4c5drho42,dv4c5drho234,dv4c5dr1,dv4c5dr2
      real*8 dv4c5dr3,dv4c5dr4,dv4c5dr5,dv4c5dr6,rho1ts9,rho2ts92,rho3ts92,rho4ts92,rho234ts9
      real*8 dv4c6ds1,dv4c6ds2,dv4c6drho1,dv4c6drho22,dv4c6drho32,dv4c6drho42,dv4c6drho234,dv4c6dr1,dv4c6dr2
      real*8 dv4c6dr3,dv4c6dr4,dv4c6dr5,dv4c6dr6,rho1ts7m,rho2ts7m2,rho3ts7m2,rho4ts7m2,rho234ts7m
      real*8 dv4c7ds1,dv4c7ds2,dv4c7drho1,dv4c7drho22,dv4c7drho32,dv4c7drho42,dv4c7drho234,dv4c7dr1,dv4c7dr2
      real*8 dv4c7dr3,dv4c7dr4,dv4c7dr5,dv4c7dr6,rho1ts13,rho2ts132,rho3ts132,rho4ts132,rho234ts13
      real*8 xc(6),cv(6),expoente,escal
      integer j,k,ncoord
      parameter (ncoord=6)
      
      escal=4.0d0
      
      
      !   geometria de referencia (h2o2)
      
      rho1eq=r3eq+r4eq+r5eq+r6eq
      rho2eq2=(r3eq+r4eq-r5eq-r6eq)**2
      rho3eq2=(r3eq-r4eq+r5eq-r6eq)**2
      rho4eq2=(r3eq-r4eq-r5eq+r6eq)**2
      rho234eq=(r3eq+r4eq-r5eq-r6eq)*(r3eq-r4eq+r5eq-r6eq)*(r3eq-r4eq-r5eq+r6eq)
      
      s1=r1-r1eq
      s2=r2-r2eq
      rho1=r3+r4+r5+r6-rho1eq
      rho22=(r3+r4-r5-r6)**2-rho2eq2
      rho32=(r3-r4+r5-r6)**2-rho3eq2
      rho42=(r3-r4-r5+r6)**2-rho4eq2
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234eq
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp1(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp1(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      
      
      pol=c(1) + c(2)*s1 + c(3)*s2 + c(4)*rho1 + c(5)*s1**2 + c(6)*s1*s2 + c(7)*s1*rho1 + c(8)*s2**2 + c(9)*s2*rho1 +&
      c(10)*rho1**2 + c(11)*rho22 + c(12)*rho32 + c(13)*rho42 + c(14)*s1**3 + c(15)*s1**2*s2 + c(16)*s1**2*rho1 +&
      c(17)*s1*s2**2 + c(18)*s1*s2*rho1 + c(19)*s1*rho1**2 + c(20)*s1*rho22 + c(21)*s1*rho32 + c(22)*s1*rho42 +&
      c(23)*s2**3 + c(24)*s2**2*rho1 + c(25)*s2*rho1**2 + c(26)*s2*rho22 + c(27)*s2*rho32 + c(28)*s2*rho42 +&
      c(29)*rho1**3 + c(30)*rho1*rho22 + c(31)*rho1*rho32 + c(32)*rho1*rho42 + c(33)*rho234 +&
      c(34)*s1**4 + c(35)*s1**3*s2 + c(36)*s1**3*rho1 + c(37)*s1**2*s2**2 + c(38)*s1**2*s2*rho1 +&
      c(39)*s1**2*rho1**2 + c(40)*s1**2*rho22 + c(41)*s1**2*rho32 + c(42)*s1**2*rho42 + c(43)*s1*s2**3 +&
      c(44)*s1*s2**2*rho1 + c(45)*s1*s2*rho1**2 + c(46)*s1*s2*rho22 + c(47)*s1*s2*rho32 + c(48)*s1*s2*rho42 +&
      c(49)*s1*rho1**3 + c(50)*s1*rho1*rho22 + c(51)*s1*rho1*rho32 + c(52)*s1*rho1*rho42 + c(53)*s1*rho234 +&
      c(54)*s2**4 + c(55)*s2**3*rho1 + c(56)*s2**2*rho1**2 + c(57)*s2**2*rho22 + c(58)*s2**2*rho32 +&
      c(59)*s2**2*rho42 + c(60)*s2*rho1**3 + c(61)*s2*rho1*rho22 + c(62)*s2*rho1*rho32 + c(63)*s2*rho1*rho42 +&
      c(64)*s2*rho234 + c(65)*rho1**4 + c(66)*rho1**2*rho22 + c(67)*rho1**2*rho32 + c(68)*rho1**2*rho42 +&
      c(69)*rho1*rho234 + c(70)*rho22*rho32 + c(71)*rho22*rho42 + c(72)*rho32*rho42 + c(73)*rho22*rho22 +&
      c(74)*rho32*rho32 + c(75)*rho42*rho42
      
      
      !derivada em ordem a s1
      
      dpolds1=c(2) + 2.0d0*c(5)*s1 + c(6)*s2 + c(7)*rho1 + 3.0d0*c(14)*s1**2 + 2.0d0*c(15)*s1*s2 + 2.0d0*c(16)*s1*rho1 +&
      c(17)*s2**2 + c(18)*s2*rho1 + c(19)*rho1**2 + c(20)*rho22 + c(21)*rho32 + c(22)*rho42 +&
      4.0d0*c(34)*s1**3 + 3.0d0*c(35)*s1**2*s2 + 3.0d0*c(36)*s1**2*rho1 + 2.0d0*c(37)*s1*s2**2 + &
      2.0d0*c(38)*s1*s2*rho1 + 2.0d0*c(39)*s1*rho1**2 + 2.0d0*c(40)*s1*rho22 + 2.0d0*c(41)*s1*rho32 +&
      2.0d0*c(42)*s1*rho42 + c(43)*s2**3 + c(44)*s2**2*rho1 + c(45)*s2*rho1**2 + c(46)*s2*rho22 + &
      c(47)*s2*rho32 + c(48)*s2*rho42 + c(49)*rho1**3 + c(50)*rho1*rho22 + c(51)*rho1*rho32 +&
      c(52)*rho1*rho42 + c(53)*rho234
      
      
      
      !derivada em ordem a s2
      
      dpolds2=c(3) + c(6)*s1 + 2.0d0*c(8)*s2 + c(9)*rho1 + c(15)*s1**2 + 2.0d0*c(17)*s1*s2 + c(18)*s1*rho1 +&
      3.0d0*c(23)*s2**2 + 2.0d0*c(24)*s2*rho1 + c(25)*rho1**2 + c(26)*rho22 + c(27)*rho32 + c(28)*rho42+&
      c(35)*s1**3 + 2.0d0*c(37)*s1**2*s2 + c(38)*s1**2*rho1 + 3.0d0*c(43)*s1*s2**2 +&
      2.0d0*c(44)*s1*s2*rho1 + c(45)*s1*rho1**2 + c(46)*s1*rho22 + c(47)*s1*rho32 + c(48)*s1*rho42 +&
      4.0d0*c(54)*s2**3 + 3.0d0*c(55)*s2**2*rho1 + 2.0d0*c(56)*s2*rho1**2 + 2.0d0*c(57)*s2*rho22 +&
      2.0d0*c(58)*s2*rho32 + 2.0d0*c(59)*s2*rho42 + c(60)*rho1**3 + c(61)*rho1*rho22 +&
      c(62)*rho1*rho32 + c(63)*rho1*rho42 + c(64)*rho234
      
      
      
      
      !derivada em ordem a rho1
      
      dpoldrho1=c(4) + c(7)*s1 + c(9)*s2 + 2.0d0*c(10)*rho1 + c(16)*s1**2 + c(18)*s1*s2 + 2.0d0*c(19)*s1*rho1 +&
      c(24)*s2**2 + 2.0d0*c(25)*s2*rho1 + 3.0d0*c(29)*rho1**2 + c(30)*rho22 + c(31)*rho32 + c(32)*rho42+&
      c(36)*s1**3 + c(38)*s1**2*s2 + 2.0d0*c(39)*s1**2*rho1 + c(44)*s1*s2**2 +&
      2.0d0*c(45)*s1*s2*rho1 + 3.0d0*c(49)*s1*rho1**2 + c(50)*s1*rho22 + c(51)*s1*rho32 + c(52)*s1*rho42 +&
      c(55)*s2**3 + 2.0d0*c(56)*s2**2*rho1 + 3.0d0*c(60)*s2*rho1**2 + c(61)*s2*rho22 + c(62)*s2*rho32 +&
      c(63)*s2*rho42 + 4.0d0*c(65)*rho1**3 + 2.0d0*c(66)*rho1*rho22 + 2.0d0*c(67)*rho1*rho32 + &
      2.0d0*c(68)*rho1*rho42 + c(69)*rho234
      
      
      
      !derivada  em ordem a rho22
      
      
      dpoldrho22=c(11) + c(20)*s1 + c(26)*s2 + c(30)*rho1+ c(40)*s1**2 + c(46)*s1*s2 + c(50)*s1*rho1 +&
      c(57)*s2**2 + c(61)*s2*rho1 + c(66)*rho1**2 + c(70)*rho32 + c(71)*rho42 + 2.0d0*c(73)*rho22
      
      
      
      !derivada em ordem a rho32
      
      dpoldrho32=c(12) + c(21)*s1 + c(27)*s2 + c(31)*rho1 + c(41)*s1**2 + c(47)*s1*s2 + c(51)*s1*rho1 +&
      c(58)*s2**2 +c(62)*s2*rho1 + c(67)*rho1**2 + c(70)*rho22 + c(72)*rho42 + 2.0d0*c(74)*rho32
      
      
      !derivada em ordem a rho42
      
      dpoldrho42=c(13) + c(22)*s1 + c(28)*s2 + c(32)*rho1 + c(42)*s1**2 + c(48)*s1*s2 + c(52)*s1*rho1 +&
      c(59)*s2**2 + c(63)*s2*rho1 + c(68)*rho1**2 + c(71)*rho22 + c(72)*rho32 + 2.0d0*c(75)*rho42
      
      !derivada em ordem a rho234
      
      dpoldrho234=c(33) + c(53)*s1 + c(64)*s2 + c(69)*rho1
      
      ddecds1=0.0d0
      ddecds2=0.0d0
      ddecdrho1=0.0d0
      ddecdrho22=0.0d0
      ddecdrho32=0.0d0
      ddecdrho42=0.0d0
      
      
      do j=1,6
      
      ddecds1=ddecds1-(2.0d0/escal)*CV(j)*vecp1(1,j)/valp1(j)*exp(-expoente*0.25d0)
      ddecds2=ddecds2-(2.0d0/escal)*CV(j)*vecp1(2,j)/valp1(j)*exp(-expoente*0.25d0)
      ddecdrho1=ddecdrho1-(2.0d0/escal)*CV(j)*vecp1(3,j)/valp1(j)*exp(-expoente*0.25d0)
      ddecdrho22=ddecdrho22-(2.0d0/escal)*CV(j)*vecp1(4,j)/valp1(j)*exp(-expoente*0.25d0)
      ddecdrho32=ddecdrho32-(2.0d0/escal)*CV(j)*vecp1(5,j)/valp1(j)*exp(-expoente*0.25d0)
      ddecdrho42=ddecdrho42-(2.0d0/escal)*CV(j)*vecp1(6,j)/valp1(j)*exp(-expoente*0.25d0)
      
      end do
      
      
      dv4c1ds1= dpolds1*dec + pol*ddecds1
      dv4c1ds2= dpolds2*dec + pol*ddecds2
      dv4c1drho1= dpoldrho1*dec + pol*ddecdrho1
      dv4c1drho22= dpoldrho22*dec + pol*ddecdrho22
      dv4c1drho32= dpoldrho32*dec + pol*ddecdrho32
      dv4c1drho42= dpoldrho42*dec + pol*ddecdrho42
      dv4c1drho234= dpoldrho234*dec
      
      
      ds1dr1=1.0d0
      ds2dr2=1.0d0
      
      drho1dr3=1.0d0
      drho22dr3=2.0d0*(r3+r4-r5-r6)
      drho32dr3=2.0d0*(r3-r4+r5-r6)
      drho42dr3=2.0d0*(r3-r4-r5+r6)
      drho234dr3=3.0d0*r3**2 -2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r3*r5 + 2.0d0*r4*r6 + 2.0d0*r5*r6 -&
      r4**2 + 2.0d0*r4*r5 - r5**2 - r6**2
      
      drho1dr4=1.0d0
      drho22dr4=2.0d0*(r3+r4-r5-r6)
      drho32dr4=-2.0d0*(r3-r4+r5-r6)
      drho42dr4=-2.0d0*(r3-r4-r5+r6)
      drho234dr4=3.0d0*r4**2 - r3**2 + 2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r4*r5 - 2.0d0*r4*r6 +&
      2.0d0*r3*r5 + 2.0d0*r5*r6 - r5**2 - r6**2
      
      drho1dr5=1.0d0
      drho22dr5=-2.0d0*(r3+r4-r5-r6)
      drho32dr5=2.0d0*(r3-r4+r5-r6)
      drho42dr5=-2.0d0*(r3-r4-r5+r6)
      drho234dr5=3.0d0*r5**2 - r3**2 + 2.0d0*r3*r6 - r4**2 + 2.0d0*r3*r4 -2.0d0*r4*r5 + 2.0d0*r4*r6 -&
      2.0d0*r5*r3 - 2.0d0*r5*r6 - r6**2
      
      
      drho1dr6=1.0d0
      drho22dr6=-2.0d0*(r3+r4-r5-r6)
      drho32dr6=-2.0d0*(r3-r4+r5-r6)
      drho42dr6=2.0d0*(r3-r4-r5+r6)
      drho234dr6=3.0d0*r6**2 - r3**2 + 2.0d0*r3*r5 + 2.0d0*r3*r4 - r4**2 + 2.0d0*r4*r5 - r5**2 -&
      2.0d0*r6*r3 -2.0d0*r6*r4 - 2.0d0*r6*r5
      
      
      
      dv4c1dr1=dv4c1ds1*ds1dr1
      
      dv4c1dr2=dv4c1ds2*ds2dr2
      
      dv4c1dr3=dv4c1drho1*drho1dr3 + dv4c1drho22*drho22dr3 + dv4c1drho32*drho32dr3 + dv4c1drho42*drho42dr3 +&
      dv4c1drho234*drho234dr3
      
      dv4c1dr4=dv4c1drho1*drho1dr4 + dv4c1drho22*drho22dr4 + dv4c1drho32*drho32dr4 + dv4c1drho42*drho42dr4 +&
      dv4c1drho234*drho234dr4
      
      dv4c1dr5=dv4c1drho1*drho1dr5 + dv4c1drho22*drho22dr5 + dv4c1drho32*drho32dr5 + dv4c1drho42*drho42dr5 +&
      dv4c1drho234*drho234dr5
      
      dv4c1dr6=dv4c1drho1*drho1dr6 + dv4c1drho22*drho22dr6 + dv4c1drho32*drho32dr6 + dv4c1drho42*drho42dr6 +&
      dv4c1drho234*drho234dr6
      
      
      
      !   com outra geometria de referencia (min1)
      
      rho1min1=r3min1+r4min1+r5min1+r6min1
      rho2min12=(r3min1+r4min1-r5min1-r6min1)**2
      rho3min12=(r3min1-r4min1+r5min1-r6min1)**2
      rho4min12=(r3min1-r4min1-r5min1+r6min1)**2
      rho234min1=(r3min1+r4min1-r5min1-r6min1)*(r3min1-r4min1+r5min1-r6min1)*(r3min1-r4min1-r5min1+r6min1)
      
      s1=r1-r1min1
      s2=r2-r2min1
      rho1=r3+r4+r5+r6-rho1min1
      rho22=(r3+r4-r5-r6)**2-rho2min12
      rho32=(r3-r4+r5-r6)**2-rho3min12
      rho42=(r3-r4-r5+r6)**2-rho4min12
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234min1
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp2(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp2(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      
      
      pol=c(76) + c(77)*s1 + c(78)*s2 + c(79)*rho1 + c(80)*s1**2 + c(81)*s1*s2 + c(82)*s1*rho1 +&
      c(83)*s2**2 + c(84)*s2*rho1 + c(85)*rho1**2 + c(86)*rho22 + c(87)*rho32 + c(88)*rho42 + c(89)*s1**3 +&
      c(90)*s1**2*s2 + c(91)*s1**2*rho1 + c(92)*s1*s2**2 + c(93)*s1*s2*rho1 + c(94)*s1*rho1**2 +&
      c(95)*s1*rho22 + c(96)*s1*rho32 + c(97)*s1*rho42 + c(98)*s2**3 + c(99)*s2**2*rho1 + c(100)*s2*rho1**2 +&
      c(101)*s2*rho22 + c(102)*s2*rho32 + c(103)*s2*rho42 + c(104)*rho1**3 + c(105)*rho1*rho22 +&
      c(106)*rho1*rho32 + c(107)*rho1*rho42 + c(108)*rho234 +&
      c(109)*s1**4 + c(110)*s1**3*s2 + c(111)*s1**3*rho1 + c(112)*s1**2*s2**2 + c(113)*s1**2*s2*rho1 +&
      c(114)*s1**2*rho1**2 + c(115)*s1**2*rho22 + c(116)*s1**2*rho32 + c(117)*s1**2*rho42 + c(118)*s1*s2**3 +&
      c(119)*s1*s2**2*rho1 + c(120)*s1*s2*rho1**2 + c(121)*s1*s2*rho22 + c(122)*s1*s2*rho32 + c(123)*s1*s2*rho42 +&
      c(124)*s1*rho1**3 + c(125)*s1*rho1*rho22 + c(126)*s1*rho1*rho32 + c(127)*s1*rho1*rho42 + c(128)*s1*rho234 +&
      c(129)*s2**4 + c(130)*s2**3*rho1 + c(131)*s2**2*rho1**2 + c(132)*s2**2*rho22 + c(133)*s2**2*rho32 +&
      c(134)*s2**2*rho42 + c(135)*s2*rho1**3 + c(136)*s2*rho1*rho22 + c(137)*s2*rho1*rho32 + c(138)*s2*rho1*rho42 +&
      c(139)*s2*rho234 + c(140)*rho1**4 + c(141)*rho1**2*rho22 + c(142)*rho1**2*rho32 + c(143)*rho1**2*rho42 +&
      c(144)*rho1*rho234 + c(145)*rho22*rho32 + c(146)*rho22*rho42 + c(147)*rho32*rho42 + c(148)*rho22*rho22 +&
      c(149)*rho32*rho32 + c(150)*rho42*rho42
      
      !derivada em ordem a s1
      dpolds1=c(77) + 2.0d0*c(80)*s1 + c(81)*s2 + c(82)*rho1 + 3.0d0*c(89)*s1**2 +&
      2.0d0*c(90)*s1*s2 + 2.0d0*c(91)*s1*rho1 + c(92)*s2**2 + c(93)*s2*rho1 + c(94)*rho1**2 +&
      c(95)*rho22 + c(96)*rho32 + c(97)*rho42 +&
      4.0d0*c(109)*s1**3 + 3.0d0*c(110)*s1**2*s2 + 3.0d0*c(111)*s1**2*rho1 + 2.0d0*c(112)*s1*s2**2 + &
      2.0d0*c(113)*s1*s2*rho1 + 2.0d0*c(114)*s1*rho1**2 + 2.0d0*c(115)*s1*rho22 + 2.0d0*c(116)*s1*rho32 +&
      2.0d0*c(117)*s1*rho42 + c(118)*s2**3 + c(119)*s2**2*rho1 + c(120)*s2*rho1**2 + c(121)*s2*rho22 + &
      c(122)*s2*rho32 + c(123)*s2*rho42 + c(124)*rho1**3 + c(125)*rho1*rho22 + c(126)*rho1*rho32 +&
      c(127)*rho1*rho42 + c(128)*rho234
      
      !derivada em ordem a s2
      dpolds2=c(78) + c(81)*s1 + 2.0d0*c(83)*s2 +c(84)*rho1 + c(90)*s1**2 +  2.0d0*c(92)*s1*s2 + c(93)*s1*rho1 +&
      3.0d0*c(98)*s2**2 + 2.0d0*c(99)*s2*rho1 + c(100)*rho1**2 + c(101)*rho22 + c(102)*rho32 + c(103)*rho42+&
      c(110)*s1**3 + 2.0d0*c(112)*s1**2*s2 + c(113)*s1**2*rho1 + 3.0d0*c(118)*s1*s2**2 +&
      2.0d0*c(119)*s1*s2*rho1 + c(120)*s1*rho1**2 + c(121)*s1*rho22 + c(122)*s1*rho32 + c(123)*s1*rho42 +&
      4.0d0*c(129)*s2**3 + 3.0d0*c(130)*s2**2*rho1 + 2.0d0*c(131)*s2*rho1**2 + 2.0d0*c(132)*s2*rho22 +&
      2.0d0*c(133)*s2*rho32 + 2.0d0*c(134)*s2*rho42 + c(135)*rho1**3 + c(136)*rho1*rho22 +&
      c(137)*rho1*rho32 + c(138)*rho1*rho42 + c(139)*rho234
      
      
      
      !derivada em ordem a rho1
      dpoldrho1=c(79) + c(82)*s1 + c(84)*s2 + 2.0d0*c(85)*rho1 + c(91)*s1**2 + c(93)*s1*s2 + 2.0d0*c(94)*s1*rho1 +&
      c(99)*s2**2 + 2.0d0*c(100)*s2*rho1 + 3.0d0*c(104)*rho1**2 + c(105)*rho22 +&
      c(106)*rho32 + c(107)*rho42 +&
      c(111)*s1**3 + c(113)*s1**2*s2 + 2.0d0*c(114)*s1**2*rho1 + c(119)*s1*s2**2 +&
      2.0d0*c(120)*s1*s2*rho1 + 3.0d0*c(124)*s1*rho1**2 + c(125)*s1*rho22 + c(126)*s1*rho32 + c(127)*s1*rho42 +&
      c(130)*s2**3 + 2.0d0*c(131)*s2**2*rho1 + 3.0d0*c(135)*s2*rho1**2 + c(136)*s2*rho22 + c(137)*s2*rho32 +&
      c(138)*s2*rho42 + 4.0d0*c(140)*rho1**3 + 2.0d0*c(141)*rho1*rho22 + 2.0d0*c(142)*rho1*rho32 + &
      2.0d0*c(143)*rho1*rho42 + c(144)*rho234
      
      
      !derivada em ordem a rho22
      dpoldrho22=c(86) + c(95)*s1 + c(101)*s2 + c(105)*rho1 + c(115)*s1**2 + c(121)*s1*s2 + c(125)*s1*rho1 +&
      c(132)*s2**2 + c(136)*s2*rho1 + c(141)*rho1**2 + c(145)*rho32 + c(146)*rho42 + 2.0d0*c(148)*rho22
      
      
      !derivada em ordem a rho32
      dpoldrho32=c(87) + c(96)*s1 + c(102)*s2 + c(106)*rho1  + c(116)*s1**2 + c(122)*s1*s2 + c(126)*s1*rho1 +&
      c(133)*s2**2 +c(137)*s2*rho1 + c(142)*rho1**2 + c(145)*rho22 + c(147)*rho42 + 2.0d0*c(149)*rho32
      
      !derivada em ordem a rho42
      dpoldrho42=c(88) + c(97)*s1 + c(103)*s2 + c(107)*rho1  + c(117)*s1**2 + c(123)*s1*s2 + c(127)*s1*rho1 +&
      c(134)*s2**2 + c(138)*s2*rho1 + c(143)*rho1**2 + c(146)*rho22 + c(147)*rho32 + 2.0d0*c(150)*rho42
      
      
      !derivada em ordem a rho234
      
      dpoldrho234=c(108) + c(128)*s1 + c(139)*s2 + c(144)*rho1
      
      
      ddecds1=0.0d0
      ddecds2=0.0d0
      ddecdrho1=0.0d0
      ddecdrho22=0.0d0
      ddecdrho32=0.0d0
      ddecdrho42=0.0d0
      
      
      do j=1,6
      
      ddecds1=ddecds1-(2.0d0/escal)*(CV(j)*vecp2(1,j))/valp2(j)*exp(-expoente*0.25d0)
      ddecds2=ddecds2-(2.0d0/escal)*(CV(j)*vecp2(2,j))/valp2(j)*exp(-expoente*0.25d0)
      ddecdrho1=ddecdrho1-(2.0d0/escal)*(CV(j)*vecp2(3,j))/valp2(j)*exp(-expoente*0.25d0)
      ddecdrho22=ddecdrho22-(2.0d0/escal)*(CV(j)*vecp2(4,j))/valp2(j)*exp(-expoente*0.25d0)
      ddecdrho32=ddecdrho32-(2.0d0/escal)*(CV(j)*vecp2(5,j))/valp2(j)*exp(-expoente*0.25d0)
      ddecdrho42=ddecdrho42-(2.0d0/escal)*(CV(j)*vecp2(6,j))/valp2(j)*exp(-expoente*0.25d0)
      
      end do
      
      
      dv4c2ds1= dpolds1*dec + pol*ddecds1
      dv4c2ds2= dpolds2*dec + pol*ddecds2
      dv4c2drho1= dpoldrho1*dec + pol*ddecdrho1
      dv4c2drho22= dpoldrho22*dec + pol*ddecdrho22
      dv4c2drho32= dpoldrho32*dec + pol*ddecdrho32
      dv4c2drho42= dpoldrho42*dec + pol*ddecdrho42
      dv4c2drho234= dpoldrho234*dec
      
      
      ds1dr1=1.0d0
      ds2dr2=1.0d0
      
      drho1dr3=1.0d0
      drho22dr3=2.0d0*(r3+r4-r5-r6)
      drho32dr3=2.0d0*(r3-r4+r5-r6)
      drho42dr3=2.0d0*(r3-r4-r5+r6)
      drho234dr3=3.0d0*r3**2 -2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r3*r5 + 2.0d0*r4*r6 + 2.0d0*r5*r6 -&
      r4**2 + 2.0d0*r4*r5 - r5**2 - r6**2
      
      drho1dr4=1.0d0
      drho22dr4=2.0d0*(r3+r4-r5-r6)
      drho32dr4=-2.0d0*(r3-r4+r5-r6)
      drho42dr4=-2.0d0*(r3-r4-r5+r6)
      drho234dr4=3.0d0*r4**2 - r3**2 + 2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r4*r5 - 2.0d0*r4*r6 +&
      2.0d0*r3*r5 + 2.0d0*r5*r6 - r5**2 - r6**2
      
      drho1dr5=1.0d0
      drho22dr5=-2.0d0*(r3+r4-r5-r6)
      drho32dr5=2.0d0*(r3-r4+r5-r6)
      drho42dr5=-2.0d0*(r3-r4-r5+r6)
      drho234dr5=3.0d0*r5**2 - r3**2 + 2.0d0*r3*r6 - r4**2 + 2.0d0*r3*r4 -2.0d0*r4*r5 + 2.0d0*r4*r6 -&
      2.0d0*r5*r3 - 2.0d0*r5*r6 - r6**2
      
      
      drho1dr6=1.0d0
      drho22dr6=-2.0d0*(r3+r4-r5-r6)
      drho32dr6=-2.0d0*(r3-r4+r5-r6)
      drho42dr6=2.0d0*(r3-r4-r5+r6)
      drho234dr6=3.0d0*r6**2 - r3**2 + 2.0d0*r3*r5 + 2.0d0*r3*r4 - r4**2 + 2.0d0*r4*r5 - r5**2 -&
      2.0d0*r6*r3 -2.0d0*r6*r4 - 2.0d0*r6*r5
      
      
      
      dv4c2dr1=dv4c2ds1*ds1dr1
      
      dv4c2dr2=dv4c2ds2*ds2dr2
      
      dv4c2dr3=dv4c2drho1*drho1dr3 + dv4c2drho22*drho22dr3 + dv4c2drho32*drho32dr3 + dv4c2drho42*drho42dr3 +&
      dv4c2drho234*drho234dr3
      
      dv4c2dr4=dv4c2drho1*drho1dr4 + dv4c2drho22*drho22dr4 + dv4c2drho32*drho32dr4 + dv4c2drho42*drho42dr4 +&
      dv4c2drho234*drho234dr4
      
      dv4c2dr5=dv4c2drho1*drho1dr5 + dv4c2drho22*drho22dr5 + dv4c2drho32*drho32dr5 + dv4c2drho42*drho42dr5 +&
      dv4c2drho234*drho234dr5
      
      dv4c2dr6=dv4c2drho1*drho1dr6 + dv4c2drho22*drho22dr6 + dv4c2drho32*drho32dr6 + dv4c2drho42*drho42dr6 +&
      dv4c2drho234*drho234dr6
      
      
      !   com outra geometria de referencia (ts7)
      
      rho1ts7=r3ts7+r4ts7+r5ts7+r6ts7
      rho2ts72=(r3ts7+r4ts7-r5ts7-r6ts7)**2
      rho3ts72=(r3ts7-r4ts7+r5ts7-r6ts7)**2
      rho4ts72=(r3ts7-r4ts7-r5ts7+r6ts7)**2
      rho234ts7=(r3ts7+r4ts7-r5ts7-r6ts7)*(r3ts7-r4ts7+r5ts7-r6ts7)*(r3ts7-r4ts7-r5ts7+r6ts7)
      
      s1=r1-r1ts7
      s2=r2-r2ts7
      rho1=r3+r4+r5+r6-rho1ts7
      rho22=(r3+r4-r5-r6)**2-rho2ts72
      rho32=(r3-r4+r5-r6)**2-rho3ts72
      rho42=(r3-r4-r5+r6)**2-rho4ts72
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts7
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp3(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp3(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      
      pol=c(151) + c(152)*s1 + c(153)*s2 + c(154)*rho1 + c(155)*s1**2 + c(156)*s1*s2 + c(157)*s1*rho1 +&
      c(158)*s2**2 +c(159)*s2*rho1 + c(160)*rho1**2 + c(161)*rho22 + c(162)*rho32 + c(163)*rho42 + c(164)*s1**3 +&
      c(165)*s1**2*s2 + c(166)*s1**2*rho1 + c(167)*s1*s2**2 + c(168)*s1*s2*rho1 + c(169)*s1*rho1**2 +&
      c(170)*s1*rho22 + c(171)*s1*rho32 + c(172)*s1*rho42 + c(173)*s2**3 + c(174)*s2**2*rho1 + c(175)*s2*rho1**2 +&
      c(176)*s2*rho22 + c(177)*s2*rho32 + c(178)*s2*rho42 + c(179)*rho1**3 + c(180)*rho1*rho22 +&
      c(181)*rho1*rho32 + c(182)*rho1*rho42 + c(183)*rho234 +&
      c(184)*s1**4 + c(185)*s1**3*s2 + c(186)*s1**3*rho1 + c(187)*s1**2*s2**2 + c(188)*s1**2*s2*rho1 +&
      c(189)*s1**2*rho1**2 + c(190)*s1**2*rho22 + c(191)*s1**2*rho32 + c(192)*s1**2*rho42 + c(193)*s1*s2**3 +&
      c(194)*s1*s2**2*rho1 + c(195)*s1*s2*rho1**2 + c(196)*s1*s2*rho22 + c(197)*s1*s2*rho32 + c(198)*s1*s2*rho42 +&
      c(199)*s1*rho1**3 + c(200)*s1*rho1*rho22 + c(201)*s1*rho1*rho32 + c(202)*s1*rho1*rho42 + c(203)*s1*rho234 +&
      c(204)*s2**4 + c(205)*s2**3*rho1 + c(206)*s2**2*rho1**2 + c(207)*s2**2*rho22 + c(208)*s2**2*rho32 +&
      c(209)*s2**2*rho42 + c(210)*s2*rho1**3 + c(211)*s2*rho1*rho22 + c(212)*s2*rho1*rho32 + c(213)*s2*rho1*rho42 +&
      c(214)*s2*rho234 + c(215)*rho1**4 + c(216)*rho1**2*rho22 + c(217)*rho1**2*rho32 + c(218)*rho1**2*rho42 +&
      c(219)*rho1*rho234 + c(220)*rho22*rho32 + c(221)*rho22*rho42 + c(222)*rho32*rho42 + c(223)*rho22*rho22 +&
      c(224)*rho32*rho32 + c(225)*rho42*rho42
      
      
      
      !derivada em ordem a s1
      dpolds1=c(152) + 2.0d0*c(155)*s1 + c(156)*s2 + c(157)*rho1 + 3.0d0*c(164)*s1**2 +&
      2.0d0*c(165)*s1*s2 + 2.0d0*c(166)*s1*rho1 + c(167)*s2**2 + c(168)*s2*rho1 + c(169)*rho1**2 +&
      c(170)*rho22 + c(171)*rho32 + c(172)*rho42+&
      4.0d0*c(184)*s1**3 + 3.0d0*c(185)*s1**2*s2 + 3.0d0*c(186)*s1**2*rho1 + 2.0d0*c(187)*s1*s2**2 + &
      2.0d0*c(188)*s1*s2*rho1 + 2.0d0*c(189)*s1*rho1**2 + 2.0d0*c(190)*s1*rho22 + 2.0d0*c(191)*s1*rho32 +&
      2.0d0*c(192)*s1*rho42 + c(193)*s2**3 + c(194)*s2**2*rho1 + c(195)*s2*rho1**2 + c(196)*s2*rho22 + &
      c(197)*s2*rho32 + c(198)*s2*rho42 + c(199)*rho1**3 + c(200)*rho1*rho22 + c(201)*rho1*rho32 +&
      c(202)*rho1*rho42 + c(203)*rho234
      
      
      
      
      !derivada em ordem a s2
      dpolds2=c(153) + c(156)*s1 + 2.0d0*c(158)*s2 +c(159)*rho1 + c(165)*s1**2 +  2.0d0*c(167)*s1*s2 + c(168)*s1*rho1 +&
      3.0d0*c(173)*s2**2 + 2.0d0*c(174)*s2*rho1 + c(175)*rho1**2 + c(176)*rho22 + c(177)*rho32 + c(178)*rho42 +&
      c(185)*s1**3 + 2.0d0*c(187)*s1**2*s2 + c(188)*s1**2*rho1 + 3.0d0*c(193)*s1*s2**2 +&
      2.0d0*c(194)*s1*s2*rho1 + c(195)*s1*rho1**2 + c(196)*s1*rho22 + c(197)*s1*rho32 + c(198)*s1*rho42 +&
      4.0d0*c(204)*s2**3 + 3.0d0*c(205)*s2**2*rho1 + 2.0d0*c(206)*s2*rho1**2 + 2.0d0*c(207)*s2*rho22 +&
      2.0d0*c(208)*s2*rho32 + 2.0d0*c(209)*s2*rho42 + c(210)*rho1**3 + c(211)*rho1*rho22 +&
      c(212)*rho1*rho32 + c(213)*rho1*rho42 + c(214)*rho234
      
      
      !derivada em ordem a rho1
      dpoldrho1=c(154) + c(157)*s1 + c(159)*s2 + 2.0d0*c(160)*rho1 + c(166)*s1**2 + c(168)*s1*s2 + 2.0d0*c(169)*s1*rho1 +&
      c(174)*s2**2 + 2.0d0*c(175)*s2*rho1 + 3.0d0*c(179)*rho1**2 + c(180)*rho22 +&
      c(181)*rho32 + c(182)*rho42 +&
      c(186)*s1**3 + c(188)*s1**2*s2 + 2.0d0*c(189)*s1**2*rho1 + c(194)*s1*s2**2 +&
      2.0d0*c(195)*s1*s2*rho1 + 3.0d0*c(199)*s1*rho1**2 + c(200)*s1*rho22 + c(201)*s1*rho32 + c(202)*s1*rho42 +&
      c(205)*s2**3 + 2.0d0*c(206)*s2**2*rho1 + 3.0d0*c(210)*s2*rho1**2 + c(211)*s2*rho22 + c(212)*s2*rho32 +&
      c(213)*s2*rho42 + 4.0d0*c(215)*rho1**3 + 2.0d0*c(216)*rho1*rho22 + 2.0d0*c(217)*rho1*rho32 + &
      2.0d0*c(218)*rho1*rho42 + c(219)*rho234
      
      
      !derivada em ordem a rho22
      dpoldrho22=c(161) + c(170)*s1 + c(176)*s2 + c(180)*rho1 + c(190)*s1**2 + c(196)*s1*s2 + c(200)*s1*rho1 +&
      c(207)*s2**2 + c(211)*s2*rho1 + c(216)*rho1**2 + c(220)*rho32 + c(221)*rho42 + 2.0d0*c(223)*rho22
      
      !derivada em ordem a rho32
      dpoldrho32=c(162) + c(171)*s1 + c(177)*s2 + c(181)*rho1  + c(191)*s1**2 + c(197)*s1*s2 + c(201)*s1*rho1 +&
      c(208)*s2**2 +c(212)*s2*rho1 + c(217)*rho1**2 + c(220)*rho22 + c(222)*rho42 + 2.0d0*c(224)*rho32
      
      !derivada em ordem a rho42
      dpoldrho42=c(163) + c(172)*s1 + c(178)*s2 + c(182)*rho1 + c(192)*s1**2 + c(198)*s1*s2 + c(202)*s1*rho1 +&
      c(209)*s2**2 + c(213)*s2*rho1 + c(218)*rho1**2 + c(221)*rho22 + c(222)*rho32 + 2.0d0*c(225)*rho42
      
      
      !derivada em ordem a rho234
      dpoldrho234=c(183)+ c(203)*s1 + c(214)*s2 + c(219)*rho1
      
      
      ddecds1=0.0d0
      ddecds2=0.0d0
      ddecdrho1=0.0d0
      ddecdrho22=0.0d0
      ddecdrho32=0.0d0
      ddecdrho42=0.0d0
      
      
      do j=1,6
      
      ddecds1=ddecds1-(2.0d0/escal)*(CV(j)*vecp3(1,j))/valp3(j)*exp(-expoente*0.25d0)
      ddecds2=ddecds2-(2.0d0/escal)*(CV(j)*vecp3(2,j))/valp3(j)*exp(-expoente*0.25d0)
      ddecdrho1=ddecdrho1-(2.0d0/escal)*(CV(j)*vecp3(3,j))/valp3(j)*exp(-expoente*0.25d0)
      ddecdrho22=ddecdrho22-(2.0d0/escal)*(CV(j)*vecp3(4,j))/valp3(j)*exp(-expoente*0.25d0)
      ddecdrho32=ddecdrho32-(2.0d0/escal)*(CV(j)*vecp3(5,j))/valp3(j)*exp(-expoente*0.25d0)
      ddecdrho42=ddecdrho42-(2.0d0/escal)*(CV(j)*vecp3(6,j))/valp3(j)*exp(-expoente*0.25d0)
      
      end do
      
      
      dv4c3ds1= dpolds1*dec + pol*ddecds1
      dv4c3ds2= dpolds2*dec + pol*ddecds2
      dv4c3drho1= dpoldrho1*dec + pol*ddecdrho1
      dv4c3drho22= dpoldrho22*dec + pol*ddecdrho22
      dv4c3drho32= dpoldrho32*dec + pol*ddecdrho32
      dv4c3drho42= dpoldrho42*dec + pol*ddecdrho42
      dv4c3drho234= dpoldrho234*dec
      
      
      ds1dr1=1.0d0
      ds2dr2=1.0d0
      
      drho1dr3=1.0d0
      drho22dr3=2.0d0*(r3+r4-r5-r6)
      drho32dr3=2.0d0*(r3-r4+r5-r6)
      drho42dr3=2.0d0*(r3-r4-r5+r6)
      drho234dr3=3.0d0*r3**2 -2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r3*r5 + 2.0d0*r4*r6 + 2.0d0*r5*r6 -&
      r4**2 + 2.0d0*r4*r5 - r5**2 - r6**2
      
      drho1dr4=1.0d0
      drho22dr4=2.0d0*(r3+r4-r5-r6)
      drho32dr4=-2.0d0*(r3-r4+r5-r6)
      drho42dr4=-2.0d0*(r3-r4-r5+r6)
      drho234dr4=3.0d0*r4**2 - r3**2 + 2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r4*r5 - 2.0d0*r4*r6 +&
      2.0d0*r3*r5 + 2.0d0*r5*r6 - r5**2 - r6**2
      
      drho1dr5=1.0d0
      drho22dr5=-2.0d0*(r3+r4-r5-r6)
      drho32dr5=2.0d0*(r3-r4+r5-r6)
      drho42dr5=-2.0d0*(r3-r4-r5+r6)
      drho234dr5=3.0d0*r5**2 - r3**2 + 2.0d0*r3*r6 - r4**2 + 2.0d0*r3*r4 -2.0d0*r4*r5 + 2.0d0*r4*r6 -&
      2.0d0*r5*r3 - 2.0d0*r5*r6 - r6**2
      
      
      drho1dr6=1.0d0
      drho22dr6=-2.0d0*(r3+r4-r5-r6)
      drho32dr6=-2.0d0*(r3-r4+r5-r6)
      drho42dr6=2.0d0*(r3-r4-r5+r6)
      drho234dr6=3.0d0*r6**2 - r3**2 + 2.0d0*r3*r5 + 2.0d0*r3*r4 - r4**2 + 2.0d0*r4*r5 - r5**2 -&
      2.0d0*r6*r3 -2.0d0*r6*r4 - 2.0d0*r6*r5
      
      
      
      dv4c3dr1=dv4c3ds1*ds1dr1
      
      dv4c3dr2=dv4c3ds2*ds2dr2
      
      dv4c3dr3=dv4c3drho1*drho1dr3 + dv4c3drho22*drho22dr3 + dv4c3drho32*drho32dr3 + dv4c3drho42*drho42dr3 +&
      dv4c3drho234*drho234dr3
      
      dv4c3dr4=dv4c3drho1*drho1dr4 + dv4c3drho22*drho22dr4 + dv4c3drho32*drho32dr4 + dv4c3drho42*drho42dr4 +&
      dv4c3drho234*drho234dr4
      
      dv4c3dr5=dv4c3drho1*drho1dr5 + dv4c3drho22*drho22dr5 + dv4c3drho32*drho32dr5 + dv4c3drho42*drho42dr5 +&
      dv4c3drho234*drho234dr5
      
      dv4c3dr6=dv4c3drho1*drho1dr6 + dv4c3drho22*drho22dr6 + dv4c3drho32*drho32dr6 + dv4c3drho42*drho42dr6 +&
      dv4c3drho234*drho234dr6
      
      
      
      
      !   com outra geometria de referencia (ts2)
      
      rho1ts2=r3ts2+r4ts2+r5ts2+r6ts2
      rho2ts22=(r3ts2+r4ts2-r5ts2-r6ts2)**2
      rho3ts22=(r3ts2-r4ts2+r5ts2-r6ts2)**2
      rho4ts22=(r3ts2-r4ts2-r5ts2+r6ts2)**2
      rho234ts2=(r3ts2+r4ts2-r5ts2-r6ts2)*(r3ts2-r4ts2+r5ts2-r6ts2)*(r3ts2-r4ts2-r5ts2+r6ts2)
      
      s1=r1-r1ts2
      s2=r2-r2ts2
      rho1=r3+r4+r5+r6-rho1ts2
      rho22=(r3+r4-r5-r6)**2-rho2ts22
      rho32=(r3-r4+r5-r6)**2-rho3ts22
      rho42=(r3-r4-r5+r6)**2-rho4ts22
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts2
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp4(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp4(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      
      
      pol=c(226) + c(227)*s1 + c(228)*s2 + c(229)*rho1 + c(230)*s1**2 + c(231)*s1*s2 + c(232)*s1*rho1 +&
      c(233)*s2**2 +c(234)*s2*rho1 + c(235)*rho1**2 + c(236)*rho22 + c(237)*rho32 + c(238)*rho42 + c(239)*s1**3 +&
      c(240)*s1**2*s2 + c(241)*s1**2*rho1 + c(242)*s1*s2**2 + c(243)*s1*s2*rho1 + c(244)*s1*rho1**2 +&
      c(245)*s1*rho22 + c(246)*s1*rho32 + c(247)*s1*rho42 + c(248)*s2**3 + c(249)*s2**2*rho1 + c(250)*s2*rho1**2 +&
      c(251)*s2*rho22 + c(252)*s2*rho32 + c(253)*s2*rho42 + c(254)*rho1**3 + c(255)*rho1*rho22 +&
      c(256)*rho1*rho32 + c(257)*rho1*rho42 + c(258)*rho234+&
      c(259)*s1**4 + c(260)*s1**3*s2 + c(261)*s1**3*rho1 + c(262)*s1**2*s2**2 + c(263)*s1**2*s2*rho1 +&
      c(264)*s1**2*rho1**2 + c(265)*s1**2*rho22 + c(266)*s1**2*rho32 + c(267)*s1**2*rho42 + c(268)*s1*s2**3 +&
      c(269)*s1*s2**2*rho1 + c(270)*s1*s2*rho1**2 + c(271)*s1*s2*rho22 + c(272)*s1*s2*rho32 + c(273)*s1*s2*rho42 +&
      c(274)*s1*rho1**3 + c(275)*s1*rho1*rho22 + c(276)*s1*rho1*rho32 + c(277)*s1*rho1*rho42 + c(278)*s1*rho234 +&
      c(279)*s2**4 + c(280)*s2**3*rho1 + c(281)*s2**2*rho1**2 + c(282)*s2**2*rho22 + c(283)*s2**2*rho32 +&
      c(284)*s2**2*rho42 + c(285)*s2*rho1**3 + c(286)*s2*rho1*rho22 + c(287)*s2*rho1*rho32 + c(288)*s2*rho1*rho42 +&
      c(289)*s2*rho234 + c(290)*rho1**4 + c(291)*rho1**2*rho22 + c(292)*rho1**2*rho32 + c(293)*rho1**2*rho42 +&
      c(294)*rho1*rho234 + c(295)*rho22*rho32 + c(296)*rho22*rho42 + c(297)*rho32*rho42 + c(298)*rho22*rho22 +&
      c(299)*rho32*rho32 + c(300)*rho42*rho42
      
      !derivada em ordem a s1
      dpolds1=c(227) + 2.0d0*c(230)*s1 + c(231)*s2 + c(232)*rho1 + 3.0d0*c(239)*s1**2 +&
      2.0d0*c(240)*s1*s2 + 2.0d0*c(241)*s1*rho1 + c(242)*s2**2 + c(243)*s2*rho1 + c(244)*rho1**2 +&
      c(245)*rho22 + c(246)*rho32 + c(247)*rho42+&
      4.0d0*c(259)*s1**3 + 3.0d0*c(260)*s1**2*s2 + 3.0d0*c(261)*s1**2*rho1 + 2.0d0*c(262)*s1*s2**2 + &
      2.0d0*c(263)*s1*s2*rho1 + 2.0d0*c(264)*s1*rho1**2 + 2.0d0*c(265)*s1*rho22 + 2.0d0*c(266)*s1*rho32 +&
      2.0d0*c(267)*s1*rho42 + c(268)*s2**3 + c(269)*s2**2*rho1 + c(270)*s2*rho1**2 + c(271)*s2*rho22 + &
      c(272)*s2*rho32 + c(273)*s2*rho42 + c(274)*rho1**3 + c(275)*rho1*rho22 + c(276)*rho1*rho32 +&
      c(277)*rho1*rho42 + c(278)*rho234
      
      !derivada em ordem a s2
      dpolds2=c(228) + c(231)*s1 + 2.0d0*c(233)*s2 +c(234)*rho1 + c(240)*s1**2 +  2.0d0*c(242)*s1*s2 + c(243)*s1*rho1 +&
      3.0d0*c(248)*s2**2 + 2.0d0*c(249)*s2*rho1 + c(250)*rho1**2 + c(251)*rho22 + c(252)*rho32 + c(253)*rho42 +&
      c(260)*s1**3 + 2.0d0*c(262)*s1**2*s2 + c(263)*s1**2*rho1 + 3.0d0*c(268)*s1*s2**2 +&
      2.0d0*c(269)*s1*s2*rho1 + c(270)*s1*rho1**2 + c(271)*s1*rho22 + c(272)*s1*rho32 + c(273)*s1*rho42 +&
      4.0d0*c(279)*s2**3 + 3.0d0*c(280)*s2**2*rho1 + 2.0d0*c(281)*s2*rho1**2 + 2.0d0*c(282)*s2*rho22 +&
      2.0d0*c(283)*s2*rho32 + 2.0d0*c(284)*s2*rho42 + c(285)*rho1**3 + c(286)*rho1*rho22 +&
      c(287)*rho1*rho32 + c(288)*rho1*rho42 + c(289)*rho234
      
      
      !derivada em ordem a rho1
      
      dpoldrho1=c(229) + c(232)*s1 + c(234)*s2 + 2.0d0*c(235)*rho1 + c(241)*s1**2 + c(243)*s1*s2 + 2.0d0*c(244)*s1*rho1 +&
      c(249)*s2**2 + 2.0d0*c(250)*s2*rho1 + 3.0d0*c(254)*rho1**2 + c(255)*rho22 +&
      c(256)*rho32 + c(257)*rho42+&
      c(261)*s1**3 + c(263)*s1**2*s2 + 2.0d0*c(264)*s1**2*rho1 + c(269)*s1*s2**2 +&
      2.0d0*c(270)*s1*s2*rho1 + 3.0d0*c(274)*s1*rho1**2 + c(275)*s1*rho22 + c(276)*s1*rho32 + c(277)*s1*rho42 +&
      c(280)*s2**3 + 2.0d0*c(281)*s2**2*rho1 + 3.0d0*c(285)*s2*rho1**2 + c(286)*s2*rho22 + c(287)*s2*rho32 +&
      c(288)*s2*rho42 + 4.0d0*c(290)*rho1**3 + 2.0d0*c(291)*rho1*rho22 + 2.0d0*c(292)*rho1*rho32 + &
      2.0d0*c(293)*rho1*rho42 + c(294)*rho234
      
      
      !derivada em ordem a rho22
      dpoldrho22=c(236) + c(245)*s1 + c(251)*s2 + c(255)*rho1 + c(265)*s1**2 + c(271)*s1*s2 + c(275)*s1*rho1 +&
      c(282)*s2**2 +c(286)*s2*rho1 + c(291)*rho1**2 + c(295)*rho32 + c(296)*rho42 + 2.0d0*c(298)*rho22
      
      
      
      !derivada em ordem a rho32
      dpoldrho32=c(237) + c(246)*s1 + c(252)*s2 + c(256)*rho1 + c(266)*s1**2 + c(272)*s1*s2 + c(276)*s1*rho1 +&
      c(283)*s2**2 +c(287)*s2*rho1 + c(292)*rho1**2 + c(295)*rho22 + c(297)*rho42 + 2.0d0*c(299)*rho32
      
      
      !derivada em ordem a rho42
      dpoldrho42=c(238) + c(247)*s1 + c(253)*s2 + c(257)*rho1 + c(267)*s1**2 + c(273)*s1*s2 + c(277)*s1*rho1 +&
      c(284)*s2**2 + c(288)*s2*rho1 + c(293)*rho1**2 + c(296)*rho22 + c(297)*rho32 + 2.0d0*c(300)*rho42
      
      
      !derivada em ordem a rho234
      dpoldrho234=c(258) + c(278)*s1 + c(289)*s2 + c(294)*rho1
      
      ddecds1=0.0d0
      ddecds2=0.0d0
      ddecdrho1=0.0d0
      ddecdrho22=0.0d0
      ddecdrho32=0.0d0
      ddecdrho42=0.0d0
      
      
      do j=1,6
      
      ddecds1=ddecds1-(2.0d0/escal)*(CV(j)*vecp4(1,j))/valp4(j)*exp(-expoente*0.25d0)
      ddecds2=ddecds2-(2.0d0/escal)*(CV(j)*vecp4(2,j))/valp4(j)*exp(-expoente*0.250)
      ddecdrho1=ddecdrho1-(2.0d0/escal)*(CV(j)*vecp4(3,j))/valp4(j)*exp(-expoente*0.25d0)
      ddecdrho22=ddecdrho22-(2.0d0/escal)*(CV(j)*vecp4(4,j))/valp4(j)*exp(-expoente*0.25d0)
      ddecdrho32=ddecdrho32-(2.0d0/escal)*(CV(j)*vecp4(5,j))/valp4(j)*exp(-expoente*0.25d0)
      ddecdrho42=ddecdrho42-(2.0d0/escal)*(CV(j)*vecp4(6,j))/valp4(j)*exp(-expoente*0.25d0)
      
      end do
      
      
      
      dv4c4ds1= dpolds1*dec + pol*ddecds1
      dv4c4ds2= dpolds2*dec + pol*ddecds2
      dv4c4drho1= dpoldrho1*dec + pol*ddecdrho1
      dv4c4drho22= dpoldrho22*dec + pol*ddecdrho22
      dv4c4drho32= dpoldrho32*dec + pol*ddecdrho32
      dv4c4drho42= dpoldrho42*dec + pol*ddecdrho42
      dv4c4drho234= dpoldrho234*dec
      
      
      ds1dr1=1.0d0
      ds2dr2=1.0d0
      
      drho1dr3=1.0d0
      drho22dr3=2.0d0*(r3+r4-r5-r6)
      drho32dr3=2.0d0*(r3-r4+r5-r6)
      drho42dr3=2.0d0*(r3-r4-r5+r6)
      drho234dr3=3.0d0*r3**2 -2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r3*r5 + 2.0d0*r4*r6 + 2.0d0*r5*r6 -&
      r4**2 + 2.0d0*r4*r5 - r5**2 - r6**2
      
      drho1dr4=1.0d0
      drho22dr4=2.0d0*(r3+r4-r5-r6)
      drho32dr4=-2.0d0*(r3-r4+r5-r6)
      drho42dr4=-2.0d0*(r3-r4-r5+r6)
      drho234dr4=3.0d0*r4**2 - r3**2 + 2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r4*r5 - 2.0d0*r4*r6 +&
      2.0d0*r3*r5 + 2.0d0*r5*r6 - r5**2 - r6**2
      
      drho1dr5=1.0d0
      drho22dr5=-2.0d0*(r3+r4-r5-r6)
      drho32dr5=2.0d0*(r3-r4+r5-r6)
      drho42dr5=-2.0d0*(r3-r4-r5+r6)
      drho234dr5=3.0d0*r5**2 - r3**2 + 2.0d0*r3*r6 - r4**2 + 2.0d0*r3*r4 -2.0d0*r4*r5 + 2.0d0*r4*r6 -&
      2.0d0*r5*r3 - 2.0d0*r5*r6 - r6**2
      
      
      drho1dr6=1.0d0
      drho22dr6=-2.0d0*(r3+r4-r5-r6)
      drho32dr6=-2.0d0*(r3-r4+r5-r6)
      drho42dr6=2.0d0*(r3-r4-r5+r6)
      drho234dr6=3.0d0*r6**2 - r3**2 + 2.0d0*r3*r5 + 2.0d0*r3*r4 - r4**2 + 2.0d0*r4*r5 - r5**2 -&
      2.0d0*r6*r3 -2.0d0*r6*r4 - 2.0d0*r6*r5
      
      
      
      dv4c4dr1=dv4c4ds1*ds1dr1
      
      dv4c4dr2=dv4c4ds2*ds2dr2
      
      dv4c4dr3=dv4c4drho1*drho1dr3 + dv4c4drho22*drho22dr3 + dv4c4drho32*drho32dr3 + dv4c4drho42*drho42dr3 +&
      dv4c4drho234*drho234dr3
      
      dv4c4dr4=dv4c4drho1*drho1dr4 + dv4c4drho22*drho22dr4 + dv4c4drho32*drho32dr4 + dv4c4drho42*drho42dr4 +&
      dv4c4drho234*drho234dr4
      
      dv4c4dr5=dv4c4drho1*drho1dr5 + dv4c4drho22*drho22dr5 + dv4c4drho32*drho32dr5 + dv4c4drho42*drho42dr5 +&
      dv4c4drho234*drho234dr5
      
      dv4c4dr6=dv4c4drho1*drho1dr6 + dv4c4drho22*drho22dr6 + dv4c4drho32*drho32dr6 + dv4c4drho42*drho42dr6 +&
      dv4c4drho234*drho234dr6
      
      
      
      !   com outra geometria de referencia (ts9)
      
      rho1ts9=r3ts9+r4ts9+r5ts9+r6ts9
      rho2ts92=(r3ts9+r4ts9-r5ts9-r6ts9)**2
      rho3ts92=(r3ts9-r4ts9+r5ts9-r6ts9)**2
      rho4ts92=(r3ts9-r4ts9-r5ts9+r6ts9)**2
      rho234ts9=(r3ts9+r4ts9-r5ts9-r6ts9)*(r3ts9-r4ts9+r5ts9-r6ts9)*(r3ts9-r4ts9-r5ts9+r6ts9)
      
      s1=r1-r1ts9
      s2=r2-r2ts9
      rho1=r3+r4+r5+r6-rho1ts9
      rho22=(r3+r4-r5-r6)**2-rho2ts92
      rho32=(r3-r4+r5-r6)**2-rho3ts92
      rho42=(r3-r4-r5+r6)**2-rho4ts92
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts9
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp5(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp5(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      
      
      pol=c(301) + c(302)*s1 + c(303)*s2 + c(304)*rho1 + c(305)*s1**2 + c(306)*s1*s2 + c(307)*s1*rho1 +&
      c(308)*s2**2 +c(309)*s2*rho1 + c(310)*rho1**2 + c(311)*rho22 + c(312)*rho32 + c(313)*rho42 + c(314)*s1**3 +&
      c(315)*s1**2*s2 + c(316)*s1**2*rho1 + c(317)*s1*s2**2 + c(318)*s1*s2*rho1 + c(319)*s1*rho1**2 +&
      c(320)*s1*rho22 + c(321)*s1*rho32 + c(322)*s1*rho42 + c(323)*s2**3 + c(324)*s2**2*rho1 + c(325)*s2*rho1**2 +&
      c(326)*s2*rho22 + c(327)*s2*rho32 + c(328)*s2*rho42 + c(329)*rho1**3 + c(330)*rho1*rho22 +&
      c(331)*rho1*rho32 + c(332)*rho1*rho42 + c(333)*rho234+&
      c(334)*s1**4 + c(335)*s1**3*s2 + c(336)*s1**3*rho1 + c(337)*s1**2*s2**2 + c(338)*s1**2*s2*rho1 +&
      c(339)*s1**2*rho1**2 + c(340)*s1**2*rho22 + c(341)*s1**2*rho32 + c(342)*s1**2*rho42 + c(343)*s1*s2**3 +&
      c(344)*s1*s2**2*rho1 + c(345)*s1*s2*rho1**2 + c(346)*s1*s2*rho22 + c(347)*s1*s2*rho32 + c(348)*s1*s2*rho42 +&
      c(349)*s1*rho1**3 + c(350)*s1*rho1*rho22 + c(351)*s1*rho1*rho32 + c(352)*s1*rho1*rho42 + c(353)*s1*rho234 +&
      c(354)*s2**4 + c(355)*s2**3*rho1 + c(356)*s2**2*rho1**2 + c(357)*s2**2*rho22 + c(358)*s2**2*rho32 +&
      c(359)*s2**2*rho42 + c(360)*s2*rho1**3 + c(361)*s2*rho1*rho22 + c(362)*s2*rho1*rho32 + c(363)*s2*rho1*rho42 +&
      c(364)*s2*rho234 + c(365)*rho1**4 + c(366)*rho1**2*rho22 + c(367)*rho1**2*rho32 + c(368)*rho1**2*rho42 +&
      c(369)*rho1*rho234 + c(370)*rho22*rho32 + c(371)*rho22*rho42 + c(372)*rho32*rho42 + c(373)*rho22*rho22 +&
      c(374)*rho32*rho32 + c(375)*rho42*rho42
      
      
      !derivada em ordem a s1
      dpolds1=c(302) + 2.0d0*c(305)*s1 + c(306)*s2 + c(307)*rho1 + 3.0d0*c(314)*s1**2 +&
      2.0d0*c(315)*s1*s2 + 2.0d0*c(316)*s1*rho1 + c(317)*s2**2 + c(318)*s2*rho1 + c(319)*rho1**2 +&
      c(320)*rho22 + c(321)*rho32 + c(322)*rho42+&
      4.0d0*c(334)*s1**3 + 3.0d0*c(335)*s1**2*s2 + 3.0d0*c(336)*s1**2*rho1 + 2.0d0*c(337)*s1*s2**2 + &
      2.0d0*c(338)*s1*s2*rho1 + 2.0d0*c(339)*s1*rho1**2 + 2.0d0*c(340)*s1*rho22 + 2.0d0*c(341)*s1*rho32 +&
      2.0d0*c(342)*s1*rho42 + c(343)*s2**3 + c(344)*s2**2*rho1 + c(345)*s2*rho1**2 + c(346)*s2*rho22 + &
      c(347)*s2*rho32 + c(348)*s2*rho42 + c(349)*rho1**3 + c(350)*rho1*rho22 + c(351)*rho1*rho32 +&
      c(352)*rho1*rho42 + c(353)*rho234
      
      !derivada em ordem a s2
      dpolds2=c(303) + c(306)*s1 + 2.0d0*c(308)*s2 +c(309)*rho1 + c(315)*s1**2 +  2.0d0*c(317)*s1*s2 + c(318)*s1*rho1 +&
      3.0d0*c(323)*s2**2 + 2.0d0*c(324)*s2*rho1 + c(325)*rho1**2 + c(326)*rho22 + c(327)*rho32 + c(328)*rho42 +&
      c(335)*s1**3 + 2.0d0*c(337)*s1**2*s2 + c(338)*s1**2*rho1 + 3.0d0*c(343)*s1*s2**2 +&
      2.0d0*c(344)*s1*s2*rho1 + c(345)*s1*rho1**2 + c(346)*s1*rho22 + c(347)*s1*rho32 + c(348)*s1*rho42 +&
      4.0d0*c(354)*s2**3 + 3.0d0*c(355)*s2**2*rho1 + 2.0d0*c(356)*s2*rho1**2 + 2.0d0*c(357)*s2*rho22 +&
      2.0d0*c(358)*s2*rho32 + 2.0d0*c(359)*s2*rho42 + c(360)*rho1**3 + c(361)*rho1*rho22 +&
      c(362)*rho1*rho32 + c(363)*rho1*rho42 + c(364)*rho234
      
      
      !derivada em ordem a rho1
      
      dpoldrho1=c(304) + c(307)*s1 + c(309)*s2 + 2.0d0*c(310)*rho1 + c(316)*s1**2 + c(318)*s1*s2 + 2.0d0*c(319)*s1*rho1 +&
      c(324)*s2**2 + 2.0d0*c(325)*s2*rho1 + 3.0d0*c(329)*rho1**2 + c(330)*rho22 +&
      c(331)*rho32 + c(332)*rho42+&
      c(336)*s1**3 + c(338)*s1**2*s2 + 2.0d0*c(339)*s1**2*rho1 + c(344)*s1*s2**2 +&
      2.0d0*c(345)*s1*s2*rho1 + 3.0d0*c(349)*s1*rho1**2 + c(350)*s1*rho22 + c(351)*s1*rho32 + c(352)*s1*rho42 +&
      c(355)*s2**3 + 2.0d0*c(356)*s2**2*rho1 + 3.0d0*c(360)*s2*rho1**2 + c(361)*s2*rho22 + c(362)*s2*rho32 +&
      c(363)*s2*rho42 + 4.0d0*c(365)*rho1**3 + 2.0d0*c(366)*rho1*rho22 + 2.0d0*c(367)*rho1*rho32 + &
      2.0d0*c(368)*rho1*rho42 + c(369)*rho234
      
      
      !derivada em ordem a rho22
      dpoldrho22=c(311) + c(320)*s1 + c(326)*s2 + c(330)*rho1 + c(340)*s1**2 + c(346)*s1*s2 + c(350)*s1*rho1 +&
      c(357)*s2**2 +c(361)*s2*rho1 + c(366)*rho1**2 + c(370)*rho32 + c(371)*rho42 + 2.0d0*c(373)*rho22
      
      
      !derivada em ordem a rho32
      dpoldrho32=c(312) + c(321)*s1 + c(327)*s2 + c(331)*rho1 + c(341)*s1**2 + c(347)*s1*s2 + c(351)*s1*rho1 +&
      c(358)*s2**2 +c(362)*s2*rho1 + c(367)*rho1**2 + c(370)*rho22 + c(372)*rho42 + 2.0d0*c(374)*rho32
      
      
      !derivada em ordem a rho42
      dpoldrho42=c(313) + c(322)*s1 + c(328)*s2 + c(332)*rho1 + c(342)*s1**2 + c(348)*s1*s2 + c(352)*s1*rho1 +&
      c(359)*s2**2 + c(363)*s2*rho1 + c(368)*rho1**2 + c(371)*rho22 + c(372)*rho32 + 2.0d0*c(375)*rho42
      
      
      !derivada em ordem a rho234
      dpoldrho234=c(333) + c(353)*s1 + c(364)*s2 + c(369)*rho1
      
      ddecds1=0.0d0
      ddecds2=0.0d0
      ddecdrho1=0.0d0
      ddecdrho22=0.0d0
      ddecdrho32=0.0d0
      ddecdrho42=0.0d0
      
      
      do j=1,6
      
      ddecds1=ddecds1-(2.0d0/escal)*(CV(j)*vecp5(1,j))/valp5(j)*exp(-expoente*0.25d0)
      ddecds2=ddecds2-(2.0d0/escal)*(CV(j)*vecp5(2,j))/valp5(j)*exp(-expoente*0.25d0)
      ddecdrho1=ddecdrho1-(2.0d0/escal)*(CV(j)*vecp5(3,j))/valp5(j)*exp(-expoente*0.25d0)
      ddecdrho22=ddecdrho22-(2.0d0/escal)*(CV(j)*vecp5(4,j))/valp5(j)*exp(-expoente*0.25d0)
      ddecdrho32=ddecdrho32-(2.0d0/escal)*(CV(j)*vecp5(5,j))/valp5(j)*exp(-expoente*0.25d0)
      ddecdrho42=ddecdrho42-(2.0d0/escal)*(CV(j)*vecp5(6,j))/valp5(j)*exp(-expoente*0.25d0)
      
      end do
      
      
      dv4c5ds1= dpolds1*dec + pol*ddecds1
      dv4c5ds2= dpolds2*dec + pol*ddecds2
      dv4c5drho1= dpoldrho1*dec + pol*ddecdrho1
      dv4c5drho22= dpoldrho22*dec + pol*ddecdrho22
      dv4c5drho32= dpoldrho32*dec + pol*ddecdrho32
      dv4c5drho42= dpoldrho42*dec + pol*ddecdrho42
      dv4c5drho234= dpoldrho234*dec
      
      
      ds1dr1=1.0d0
      ds2dr2=1.0d0
      
      drho1dr3=1.0d0
      drho22dr3=2.0d0*(r3+r4-r5-r6)
      drho32dr3=2.0d0*(r3-r4+r5-r6)
      drho42dr3=2.0d0*(r3-r4-r5+r6)
      drho234dr3=3.0d0*r3**2 -2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r3*r5 + 2.0d0*r4*r6 + 2.0d0*r5*r6 -&
      r4**2 + 2.0d0*r4*r5 - r5**2 - r6**2
      
      drho1dr4=1.0d0
      drho22dr4=2.0d0*(r3+r4-r5-r6)
      drho32dr4=-2.0d0*(r3-r4+r5-r6)
      drho42dr4=-2.0d0*(r3-r4-r5+r6)
      drho234dr4=3.0d0*r4**2 - r3**2 + 2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r4*r5 - 2.0d0*r4*r6 +&
      2.0d0*r3*r5 + 2.0d0*r5*r6 - r5**2 - r6**2
      
      drho1dr5=1.0d0
      drho22dr5=-2.0d0*(r3+r4-r5-r6)
      drho32dr5=2.0d0*(r3-r4+r5-r6)
      drho42dr5=-2.0d0*(r3-r4-r5+r6)
      drho234dr5=3.0d0*r5**2 - r3**2 + 2.0d0*r3*r6 - r4**2 + 2.0d0*r3*r4 -2.0d0*r4*r5 + 2.0d0*r4*r6 -&
      2.0d0*r5*r3 - 2.0d0*r5*r6 - r6**2
      
      
      drho1dr6=1.0d0
      drho22dr6=-2.0d0*(r3+r4-r5-r6)
      drho32dr6=-2.0d0*(r3-r4+r5-r6)
      drho42dr6=2.0d0*(r3-r4-r5+r6)
      drho234dr6=3.0d0*r6**2 - r3**2 + 2.0d0*r3*r5 + 2.0d0*r3*r4 - r4**2 + 2.0d0*r4*r5 - r5**2 -&
      2.0d0*r6*r3 -2.0d0*r6*r4 - 2.0d0*r6*r5
      
      
      
      dv4c5dr1=dv4c5ds1*ds1dr1
      
      dv4c5dr2=dv4c5ds2*ds2dr2
      
      dv4c5dr3=dv4c5drho1*drho1dr3 + dv4c5drho22*drho22dr3 + dv4c5drho32*drho32dr3 + dv4c5drho42*drho42dr3 +&
      dv4c5drho234*drho234dr3
      
      dv4c5dr4=dv4c5drho1*drho1dr4 + dv4c5drho22*drho22dr4 + dv4c5drho32*drho32dr4 + dv4c5drho42*drho42dr4 +&
      dv4c5drho234*drho234dr4
      
      dv4c5dr5=dv4c5drho1*drho1dr5 + dv4c5drho22*drho22dr5 + dv4c5drho32*drho32dr5 + dv4c5drho42*drho42dr5 +&
      dv4c5drho234*drho234dr5
      
      dv4c5dr6=dv4c5drho1*drho1dr6 + dv4c5drho22*drho22dr6 + dv4c5drho32*drho32dr6 + dv4c5drho42*drho42dr6 +&
      dv4c5drho234*drho234dr6
      
      !   com outra geometria de referencia (ts7 media)
      
      rho1ts7m=r3ts7m+r4ts7m+r5ts7m+r6ts7m
      rho2ts7m2=(r3ts7m+r4ts7m-r5ts7m-r6ts7m)**2
      rho3ts7m2=(r3ts7m-r4ts7m+r5ts7m-r6ts7m)**2
      rho4ts7m2=(r3ts7m-r4ts7m-r5ts7m+r6ts7m)**2
      rho234ts7m=(r3ts7m+r4ts7m-r5ts7m-r6ts7m)*(r3ts7m-r4ts7m+r5ts7m-r6ts7m)*(r3ts7m-r4ts7m-r5ts7m+r6ts7m)
      
      s1=r1-r1ts7m
      s2=r2-r2ts7m
      rho1=r3+r4+r5+r6-rho1ts7m
      rho22=(r3+r4-r5-r6)**2-rho2ts7m2
      rho32=(r3-r4+r5-r6)**2-rho3ts7m2
      rho42=(r3-r4-r5+r6)**2-rho4ts7m2
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts7m
      
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp6(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp6(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      
      
      pol=c(376) + c(377)*s1 + c(378)*s2 + c(379)*rho1 + c(380)*s1**2 + c(381)*s1*s2 + c(382)*s1*rho1 +&
      c(383)*s2**2 +c(384)*s2*rho1 + c(385)*rho1**2 + c(386)*rho22 + c(387)*rho32 + c(388)*rho42 + c(389)*s1**3 +&
      c(390)*s1**2*s2 + c(391)*s1**2*rho1 + c(392)*s1*s2**2 + c(393)*s1*s2*rho1 + c(394)*s1*rho1**2 +&
      c(395)*s1*rho22 + c(396)*s1*rho32 + c(397)*s1*rho42 + c(398)*s2**3 + c(399)*s2**2*rho1 + c(400)*s2*rho1**2 +&
      c(401)*s2*rho22 + c(402)*s2*rho32 + c(403)*s2*rho42 + c(404)*rho1**3 + c(405)*rho1*rho22 +&
      c(406)*rho1*rho32 + c(407)*rho1*rho42 + c(408)*rho234+&
      c(409)*s1**4 + c(410)*s1**3*s2 + c(411)*s1**3*rho1 + c(412)*s1**2*s2**2 + c(413)*s1**2*s2*rho1 +&
      c(414)*s1**2*rho1**2 + c(415)*s1**2*rho22 + c(416)*s1**2*rho32 + c(417)*s1**2*rho42 + c(418)*s1*s2**3 +&
      c(419)*s1*s2**2*rho1 + c(420)*s1*s2*rho1**2 + c(421)*s1*s2*rho22 + c(422)*s1*s2*rho32 + c(423)*s1*s2*rho42 +&
      c(424)*s1*rho1**3 + c(425)*s1*rho1*rho22 + c(426)*s1*rho1*rho32 + c(427)*s1*rho1*rho42 + c(428)*s1*rho234 +&
      c(429)*s2**4 + c(430)*s2**3*rho1 + c(431)*s2**2*rho1**2 + c(432)*s2**2*rho22 + c(433)*s2**2*rho32 +&
      c(434)*s2**2*rho42 + c(435)*s2*rho1**3 + c(436)*s2*rho1*rho22 + c(437)*s2*rho1*rho32 + c(438)*s2*rho1*rho42 +&
      c(439)*s2*rho234 + c(440)*rho1**4 + c(441)*rho1**2*rho22 + c(442)*rho1**2*rho32 + c(443)*rho1**2*rho42 +&
      c(444)*rho1*rho234 + c(445)*rho22*rho32 + c(446)*rho22*rho42 + c(447)*rho32*rho42 + c(448)*rho22*rho22 +&
      c(449)*rho32*rho32 + c(450)*rho42*rho42
      
      
      !derivada em ordem a s1
      dpolds1=c(377) + 2.0d0*c(380)*s1 + c(381)*s2 + c(382)*rho1 + 3.0d0*c(389)*s1**2 +&
      2.0d0*c(390)*s1*s2 + 2.0d0*c(391)*s1*rho1 + c(392)*s2**2 + c(393)*s2*rho1 + c(394)*rho1**2 +&
      c(395)*rho22 + c(396)*rho32 + c(397)*rho42+&
      4.0d0*c(409)*s1**3 + 3.0d0*c(410)*s1**2*s2 + 3.0d0*c(411)*s1**2*rho1 + 2.0d0*c(412)*s1*s2**2 + &
      2.0d0*c(413)*s1*s2*rho1 + 2.0d0*c(414)*s1*rho1**2 + 2.0d0*c(415)*s1*rho22 + 2.0d0*c(416)*s1*rho32 +&
      2.0d0*c(417)*s1*rho42 + c(418)*s2**3 + c(419)*s2**2*rho1 + c(420)*s2*rho1**2 + c(421)*s2*rho22 + &
      c(422)*s2*rho32 + c(423)*s2*rho42 + c(424)*rho1**3 + c(425)*rho1*rho22 + c(426)*rho1*rho32 +&
      c(427)*rho1*rho42 + c(428)*rho234
      
      !derivada em ordem a s2
      dpolds2=c(378) + c(381)*s1 + 2.0d0*c(383)*s2 +c(384)*rho1 + c(390)*s1**2 +  2.0d0*c(392)*s1*s2 + c(393)*s1*rho1 +&
      3.0d0*c(398)*s2**2 + 2.0d0*c(399)*s2*rho1 + c(400)*rho1**2 + c(401)*rho22 + c(402)*rho32 + c(403)*rho42 +&
      c(410)*s1**3 + 2.0d0*c(412)*s1**2*s2 + c(413)*s1**2*rho1 + 3.0d0*c(418)*s1*s2**2 +&
      2.0d0*c(419)*s1*s2*rho1 + c(420)*s1*rho1**2 + c(421)*s1*rho22 + c(422)*s1*rho32 + c(423)*s1*rho42 +&
      4.0d0*c(429)*s2**3 + 3.0d0*c(430)*s2**2*rho1 + 2.0d0*c(431)*s2*rho1**2 + 2.0d0*c(432)*s2*rho22 +&
      2.0d0*c(433)*s2*rho32 + 2.0d0*c(434)*s2*rho42 + c(435)*rho1**3 + c(436)*rho1*rho22 +&
      c(437)*rho1*rho32 + c(438)*rho1*rho42 + c(439)*rho234
      
      
      
      !derivada em ordem a rho1
      
      dpoldrho1=c(379) + c(382)*s1 + c(384)*s2 + 2.0d0*c(385)*rho1 + c(391)*s1**2 + c(393)*s1*s2 + 2.0d0*c(394)*s1*rho1 +&
      c(399)*s2**2 + 2.0d0*c(400)*s2*rho1 + 3.0d0*c(404)*rho1**2 + c(405)*rho22 +&
      c(406)*rho32 + c(407)*rho42+&
      c(411)*s1**3 + c(413)*s1**2*s2 + 2.0d0*c(414)*s1**2*rho1 + c(419)*s1*s2**2 +&
      2.0d0*c(420)*s1*s2*rho1 + 3.0d0*c(424)*s1*rho1**2 + c(425)*s1*rho22 + c(426)*s1*rho32 + c(427)*s1*rho42 +&
      c(430)*s2**3 + 2.0d0*c(431)*s2**2*rho1 + 3.0d0*c(435)*s2*rho1**2 + c(436)*s2*rho22 + c(437)*s2*rho32 +&
      c(438)*s2*rho42 + 4.0d0*c(440)*rho1**3 + 2.0d0*c(441)*rho1*rho22 + 2.0d0*c(442)*rho1*rho32 + &
      2.0d0*c(443)*rho1*rho42 + c(444)*rho234
      
      
      !derivada em ordem a rho22
      dpoldrho22=c(386) + c(395)*s1 + c(401)*s2 + c(405)*rho1 + c(415)*s1**2 + c(421)*s1*s2 + c(425)*s1*rho1 +&
      c(432)*s2**2 +c(436)*s2*rho1 + c(441)*rho1**2 + c(445)*rho32 + c(446)*rho42 + 2.0d0*c(448)*rho22
      
      
      !derivada em ordem a rho32
      dpoldrho32=c(387) + c(396)*s1 + c(402)*s2 + c(406)*rho1 + c(416)*s1**2 + c(422)*s1*s2 + c(426)*s1*rho1 +&
      c(433)*s2**2 +c(437)*s2*rho1 + c(442)*rho1**2 + c(445)*rho22 + c(447)*rho42 + 2.0d0*c(449)*rho32
      
      
      !derivada em ordem a rho42
      dpoldrho42=c(388) + c(397)*s1 + c(403)*s2 + c(407)*rho1 + c(417)*s1**2 + c(423)*s1*s2 + c(427)*s1*rho1 +&
      c(434)*s2**2 + c(438)*s2*rho1 + c(443)*rho1**2 + c(446)*rho22 + c(447)*rho32 + 2.0d0*c(450)*rho42
      
      
      !derivada em ordem a rho234
      dpoldrho234=c(408) + c(428)*s1 + c(439)*s2 + c(444)*rho1
      
      ddecds1=0.0d0
      ddecds2=0.0d0
      ddecdrho1=0.0d0
      ddecdrho22=0.0d0
      ddecdrho32=0.0d0
      ddecdrho42=0.0d0
      
      
      do j=1,6
      
      ddecds1=ddecds1-(2.0d0/escal)*(CV(j)*vecp6(1,j))/valp6(j)*exp(-expoente*0.25d0)
      ddecds2=ddecds2-(2.0d0/escal)*(CV(j)*vecp6(2,j))/valp6(j)*exp(-expoente*0.25d0)
      ddecdrho1=ddecdrho1-(2.0d0/escal)*(CV(j)*vecp6(3,j))/valp6(j)*exp(-expoente*0.25d0)
      ddecdrho22=ddecdrho22-(2.0d0/escal)*(CV(j)*vecp6(4,j))/valp6(j)*exp(-expoente*0.25d0)
      ddecdrho32=ddecdrho32-(2.0d0/escal)*(CV(j)*vecp6(5,j))/valp6(j)*exp(-expoente*0.25d0)
      ddecdrho42=ddecdrho42-(2.0d0/escal)*(CV(j)*vecp6(6,j))/valp6(j)*exp(-expoente*0.25d0)
      
      end do
      
      
      dv4c6ds1= dpolds1*dec + pol*ddecds1
      dv4c6ds2= dpolds2*dec + pol*ddecds2
      dv4c6drho1= dpoldrho1*dec + pol*ddecdrho1
      dv4c6drho22= dpoldrho22*dec + pol*ddecdrho22
      dv4c6drho32= dpoldrho32*dec + pol*ddecdrho32
      dv4c6drho42= dpoldrho42*dec + pol*ddecdrho42
      dv4c6drho234= dpoldrho234*dec
      
      
      ds1dr1=1.0d0
      ds2dr2=1.0d0
      
      drho1dr3=1.0d0
      drho22dr3=2.0d0*(r3+r4-r5-r6)
      drho32dr3=2.0d0*(r3-r4+r5-r6)
      drho42dr3=2.0d0*(r3-r4-r5+r6)
      drho234dr3=3.0d0*r3**2 -2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r3*r5 + 2.0d0*r4*r6 + 2.0d0*r5*r6 -&
      r4**2 + 2.0d0*r4*r5 - r5**2 - r6**2
      
      drho1dr4=1.0d0
      drho22dr4=2.0d0*(r3+r4-r5-r6)
      drho32dr4=-2.0d0*(r3-r4+r5-r6)
      drho42dr4=-2.0d0*(r3-r4-r5+r6)
      drho234dr4=3.0d0*r4**2 - r3**2 + 2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r4*r5 - 2.0d0*r4*r6 +&
      2.0d0*r3*r5 + 2.0d0*r5*r6 - r5**2 - r6**2
      
      drho1dr5=1.0d0
      drho22dr5=-2.0d0*(r3+r4-r5-r6)
      drho32dr5=2.0d0*(r3-r4+r5-r6)
      drho42dr5=-2.0d0*(r3-r4-r5+r6)
      drho234dr5=3.0d0*r5**2 - r3**2 + 2.0d0*r3*r6 - r4**2 + 2.0d0*r3*r4 -2.0d0*r4*r5 + 2.0d0*r4*r6 -&
      2.0d0*r5*r3 - 2.0d0*r5*r6 - r6**2
      
      
      drho1dr6=1.0d0
      drho22dr6=-2.0d0*(r3+r4-r5-r6)
      drho32dr6=-2.0d0*(r3-r4+r5-r6)
      drho42dr6=2.0d0*(r3-r4-r5+r6)
      drho234dr6=3.0d0*r6**2 - r3**2 + 2.0d0*r3*r5 + 2.0d0*r3*r4 - r4**2 + 2.0d0*r4*r5 - r5**2 -&
      2.0d0*r6*r3 -2.0d0*r6*r4 - 2.0d0*r6*r5
      
      
      
      dv4c6dr1=dv4c6ds1*ds1dr1
      
      dv4c6dr2=dv4c6ds2*ds2dr2
      
      dv4c6dr3=dv4c6drho1*drho1dr3 + dv4c6drho22*drho22dr3 + dv4c6drho32*drho32dr3 + dv4c6drho42*drho42dr3 +&
      dv4c6drho234*drho234dr3
      
      dv4c6dr4=dv4c6drho1*drho1dr4 + dv4c6drho22*drho22dr4 + dv4c6drho32*drho32dr4 + dv4c6drho42*drho42dr4 +&
      dv4c6drho234*drho234dr4
      
      dv4c6dr5=dv4c6drho1*drho1dr5 + dv4c6drho22*drho22dr5 + dv4c6drho32*drho32dr5 + dv4c6drho42*drho42dr5 +&
      dv4c6drho234*drho234dr5
      
      dv4c6dr6=dv4c6drho1*drho1dr6 + dv4c6drho22*drho22dr6 + dv4c6drho32*drho32dr6 + dv4c6drho42*drho42dr6 +&
      dv4c6drho234*drho234dr6
      
      
      !derivada em ordem a geometria do estado de transicao ts13
      
      rho1ts13=r3ts13+r4ts13+r5ts13+r6ts13
      rho2ts132=(r3ts13+r4ts13-r5ts13-r6ts13)**2
      rho3ts132=(r3ts13-r4ts13+r5ts13-r6ts13)**2
      rho4ts132=(r3ts13-r4ts13-r5ts13+r6ts13)**2
      rho234ts13=(r3ts13+r4ts13-r5ts13-r6ts13)*(r3ts13-r4ts13+r5ts13-r6ts13)*(r3ts13-r4ts13-r5ts13+r6ts13)
      
      s1=r1-r1ts13
      s2=r2-r2ts13
      rho1=r3+r4+r5+r6-rho1ts13
      rho22=(r3+r4-r5-r6)**2-rho2ts132
      rho32=(r3-r4+r5-r6)**2-rho3ts132
      rho42=(r3-r4-r5+r6)**2-rho4ts132
      rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234ts13
      
      
      
      
      xc(1)=s1
      xc(2)=s2
      xc(3)=rho1
      xc(4)=rho22
      xc(5)=rho32
      xc(6)=rho42
      
      DO j=1,6
      CV(j)=0.0d0
      DO k=1,ncoord
      CV(j)=CV(j)+XC(k)*vecp7(k,j)
      END DO
      END DO
      
      expoente=0.0d0
      do J=1,6
      expoente=expoente + CV(j)**2/valp7(j)
      end do
      
      dec=exp(-expoente*0.25d0)
      
      
      pol=c(451) + c(452)*s1 + c(453)*s2 + c(454)*rho1 + c(455)*s1**2 + c(456)*s1*s2 + c(457)*s1*rho1 +&
      c(458)*s2**2 +c(459)*s2*rho1 + c(460)*rho1**2 + c(461)*rho22 + c(462)*rho32 + c(463)*rho42 + c(464)*s1**3 +&
      c(465)*s1**2*s2 + c(466)*s1**2*rho1 + c(467)*s1*s2**2 + c(468)*s1*s2*rho1 + c(469)*s1*rho1**2 +&
      c(470)*s1*rho22 + c(471)*s1*rho32 + c(472)*s1*rho42 + c(473)*s2**3 + c(474)*s2**2*rho1 + c(475)*s2*rho1**2 +&
      c(476)*s2*rho22 + c(477)*s2*rho32 + c(478)*s2*rho42 + c(479)*rho1**3 + c(480)*rho1*rho22 +&
      c(481)*rho1*rho32 + c(482)*rho1*rho42 + c(483)*rho234+&
      c(484)*s1**4 + c(485)*s1**3*s2 + c(486)*s1**3*rho1 + c(487)*s1**2*s2**2 + c(488)*s1**2*s2*rho1 +&
      c(489)*s1**2*rho1**2 + c(490)*s1**2*rho22 + c(491)*s1**2*rho32 + c(492)*s1**2*rho42 + c(493)*s1*s2**3 +&
      c(494)*s1*s2**2*rho1 + c(495)*s1*s2*rho1**2 + c(496)*s1*s2*rho22 + c(497)*s1*s2*rho32 + c(498)*s1*s2*rho42 +&
      c(499)*s1*rho1**3 + c(500)*s1*rho1*rho22 + c(501)*s1*rho1*rho32 + c(502)*s1*rho1*rho42 + c(503)*s1*rho234 +&
      c(504)*s2**4 + c(505)*s2**3*rho1 + c(506)*s2**2*rho1**2 + c(507)*s2**2*rho22 + c(508)*s2**2*rho32 +&
      c(509)*s2**2*rho42 + c(510)*s2*rho1**3 + c(511)*s2*rho1*rho22 + c(512)*s2*rho1*rho32 + c(513)*s2*rho1*rho42 +&
      c(514)*s2*rho234 + c(515)*rho1**4 + c(516)*rho1**2*rho22 + c(517)*rho1**2*rho32 + c(518)*rho1**2*rho42 +&
      c(519)*rho1*rho234 + c(520)*rho22*rho32 + c(521)*rho22*rho42 + c(522)*rho32*rho42 + c(523)*rho22*rho22 +&
      c(524)*rho32*rho32 + c(525)*rho42*rho42
      
      !derivada em ordem a s1
      dpolds1=c(452) + 2.0d0*c(455)*s1 + c(456)*s2 + c(457)*rho1 + 3.0d0*c(464)*s1**2 +&
      2.0d0*c(465)*s1*s2 + 2.0d0*c(466)*s1*rho1 + c(467)*s2**2 + c(468)*s2*rho1 + c(469)*rho1**2 +&
      c(470)*rho22 + c(471)*rho32 + c(472)*rho42+&
      4.0d0*c(484)*s1**3 + 3.0d0*c(485)*s1**2*s2 + 3.0d0*c(486)*s1**2*rho1 + 2.0d0*c(487)*s1*s2**2 + &
      2.0d0*c(488)*s1*s2*rho1 + 2.0d0*c(489)*s1*rho1**2 + 2.0d0*c(490)*s1*rho22 + 2.0d0*c(491)*s1*rho32 +&
      2.0d0*c(492)*s1*rho42 + c(493)*s2**3 + c(494)*s2**2*rho1 + c(495)*s2*rho1**2 + c(496)*s2*rho22 + &
      c(497)*s2*rho32 + c(498)*s2*rho42 + c(499)*rho1**3 + c(500)*rho1*rho22 + c(501)*rho1*rho32 +&
      c(502)*rho1*rho42 + c(503)*rho234
      
      !derivada em ordem a s2
      dpolds2=c(453) + c(456)*s1 + 2.0d0*c(458)*s2 +c(459)*rho1 + c(465)*s1**2 +  2.0d0*c(467)*s1*s2 + c(468)*s1*rho1 +&
      3.0d0*c(473)*s2**2 + 2.0d0*c(474)*s2*rho1 + c(475)*rho1**2 + c(476)*rho22 + c(477)*rho32 + c(478)*rho42 +&
      c(485)*s1**3 + 2.0d0*c(487)*s1**2*s2 + c(488)*s1**2*rho1 + 3.0d0*c(493)*s1*s2**2 +&
      2.0d0*c(494)*s1*s2*rho1 + c(495)*s1*rho1**2 + c(496)*s1*rho22 + c(497)*s1*rho32 + c(498)*s1*rho42 +&
      4.0d0*c(504)*s2**3 + 3.0d0*c(505)*s2**2*rho1 + 2.0d0*c(506)*s2*rho1**2 + 2.0d0*c(507)*s2*rho22 +&
      2.0d0*c(508)*s2*rho32 + 2.0d0*c(509)*s2*rho42 + c(510)*rho1**3 + c(511)*rho1*rho22 +&
      c(512)*rho1*rho32 + c(513)*rho1*rho42 + c(514)*rho234
      
      
      !derivada em ordem a rho1
      
      dpoldrho1=c(454) + c(457)*s1 + c(459)*s2 + 2.0d0*c(460)*rho1 + c(466)*s1**2 + c(468)*s1*s2 + 2.0d0*c(469)*s1*rho1 +&
      c(474)*s2**2 + 2.0d0*c(475)*s2*rho1 + 3.0d0*c(479)*rho1**2 + c(480)*rho22 +&
      c(481)*rho32 + c(482)*rho42+&
      c(486)*s1**3 + c(488)*s1**2*s2 + 2.0d0*c(489)*s1**2*rho1 + c(494)*s1*s2**2 +&
      2.0d0*c(495)*s1*s2*rho1 + 3.0d0*c(499)*s1*rho1**2 + c(500)*s1*rho22 + c(501)*s1*rho32 + c(502)*s1*rho42 +&
      c(505)*s2**3 + 2.0d0*c(506)*s2**2*rho1 + 3.0d0*c(510)*s2*rho1**2 + c(511)*s2*rho22 + c(512)*s2*rho32 +&
      c(513)*s2*rho42 + 4.0d0*c(515)*rho1**3 + 2.0d0*c(516)*rho1*rho22 + 2.0d0*c(517)*rho1*rho32 + &
      2.0d0*c(518)*rho1*rho42 + c(519)*rho234
      
      
      !derivada em ordem a rho22
      dpoldrho22=c(461) + c(470)*s1 + c(476)*s2 + c(480)*rho1 + c(490)*s1**2 + c(496)*s1*s2 + c(500)*s1*rho1 +&
      c(507)*s2**2 +c(511)*s2*rho1 + c(516)*rho1**2 + c(520)*rho32 + c(521)*rho42 + 2.0d0*c(523)*rho22
      
      
      
      !derivada em ordem a rho32
      dpoldrho32=c(462) + c(471)*s1 + c(477)*s2 + c(481)*rho1 + c(491)*s1**2 + c(497)*s1*s2 + c(501)*s1*rho1 +&
      c(508)*s2**2 +c(512)*s2*rho1 + c(517)*rho1**2 + c(520)*rho22 + c(522)*rho42 + 2.0d0*c(524)*rho32
      
      
      !derivada em ordem a rho42
      dpoldrho42=c(463) + c(472)*s1 + c(478)*s2 + c(482)*rho1 + c(492)*s1**2 + c(498)*s1*s2 + c(502)*s1*rho1 +&
      c(509)*s2**2 + c(513)*s2*rho1 + c(518)*rho1**2 + c(521)*rho22 + c(522)*rho32 + 2.0d0*c(525)*rho42
      
      
      !derivada em ordem a rho234
      dpoldrho234=c(483) + c(503)*s1 + c(514)*s2 + c(519)*rho1
      
      
      ddecds1=0.0d0
      ddecds2=0.0d0
      ddecdrho1=0.0d0
      ddecdrho22=0.0d0
      ddecdrho32=0.0d0
      ddecdrho42=0.0d0
      
      
      do j=1,6
      
      ddecds1=ddecds1-(2.0d0/escal)*(CV(j)*vecp7(1,j))/valp7(j)*exp(-expoente*0.25d0)
      ddecds2=ddecds2-(2.0d0/escal)*(CV(j)*vecp7(2,j))/valp7(j)*exp(-expoente*0.25d0)
      ddecdrho1=ddecdrho1-(2.0d0/escal)*(CV(j)*vecp7(3,j))/valp7(j)*exp(-expoente*0.25d0)
      ddecdrho22=ddecdrho22-(2.0d0/escal)*(CV(j)*vecp7(4,j))/valp7(j)*exp(-expoente*0.25d0)
      ddecdrho32=ddecdrho32-(2.0d0/escal)*(CV(j)*vecp7(5,j))/valp7(j)*exp(-expoente*0.25d0)
      ddecdrho42=ddecdrho42-(2.0d0/escal)*(CV(j)*vecp7(6,j))/valp7(j)*exp(-expoente*0.25d0)
      
      end do
      
      
      dv4c7ds1= dpolds1*dec + pol*ddecds1
      dv4c7ds2= dpolds2*dec + pol*ddecds2
      dv4c7drho1= dpoldrho1*dec + pol*ddecdrho1
      dv4c7drho22= dpoldrho22*dec + pol*ddecdrho22
      dv4c7drho32= dpoldrho32*dec + pol*ddecdrho32
      dv4c7drho42= dpoldrho42*dec + pol*ddecdrho42
      dv4c7drho234= dpoldrho234*dec
      
      
      ds1dr1=1.0d0
      ds2dr2=1.0d0
      
      drho1dr3=1.0d0
      drho22dr3=2.0d0*(r3+r4-r5-r6)
      drho32dr3=2.0d0*(r3-r4+r5-r6)
      drho42dr3=2.0d0*(r3-r4-r5+r6)
      drho234dr3=3.0d0*r3**2 -2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r3*r5 + 2.0d0*r4*r6 + 2.0d0*r5*r6 -&
      r4**2 + 2.0d0*r4*r5 - r5**2 - r6**2
      
      drho1dr4=1.0d0
      drho22dr4=2.0d0*(r3+r4-r5-r6)
      drho32dr4=-2.0d0*(r3-r4+r5-r6)
      drho42dr4=-2.0d0*(r3-r4-r5+r6)
      drho234dr4=3.0d0*r4**2 - r3**2 + 2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r4*r5 - 2.0d0*r4*r6 +&
      2.0d0*r3*r5 + 2.0d0*r5*r6 - r5**2 - r6**2
      
      drho1dr5=1.0d0
      drho22dr5=-2.0d0*(r3+r4-r5-r6)
      drho32dr5=2.0d0*(r3-r4+r5-r6)
      drho42dr5=-2.0d0*(r3-r4-r5+r6)
      drho234dr5=3.0d0*r5**2 - r3**2 + 2.0d0*r3*r6 - r4**2 + 2.0d0*r3*r4 -2.0d0*r4*r5 + 2.0d0*r4*r6 -&
      2.0d0*r5*r3 - 2.0d0*r5*r6 - r6**2
      
      
      drho1dr6=1.0d0
      drho22dr6=-2.0d0*(r3+r4-r5-r6)
      drho32dr6=-2.0d0*(r3-r4+r5-r6)
      drho42dr6=2.0d0*(r3-r4-r5+r6)
      drho234dr6=3.0d0*r6**2 - r3**2 + 2.0d0*r3*r5 + 2.0d0*r3*r4 - r4**2 + 2.0d0*r4*r5 - r5**2 -&
      2.0d0*r6*r3 -2.0d0*r6*r4 - 2.0d0*r6*r5
      
      
      
      dv4c7dr1=dv4c7ds1*ds1dr1
      
      dv4c7dr2=dv4c7ds2*ds2dr2
      
      dv4c7dr3=dv4c7drho1*drho1dr3 + dv4c7drho22*drho22dr3 + dv4c7drho32*drho32dr3 + dv4c7drho42*drho42dr3 +&
      dv4c7drho234*drho234dr3
      
      dv4c7dr4=dv4c7drho1*drho1dr4 + dv4c7drho22*drho22dr4 + dv4c7drho32*drho32dr4 + dv4c7drho42*drho42dr4 +&
      dv4c7drho234*drho234dr4
      
      dv4c7dr5=dv4c7drho1*drho1dr5 + dv4c7drho22*drho22dr5 + dv4c7drho32*drho32dr5 + dv4c7drho42*drho42dr5 +&
      dv4c7drho234*drho234dr5
      
      dv4c7dr6=dv4c7drho1*drho1dr6 + dv4c7drho22*drho22dr6 + dv4c7drho32*drho32dr6 + dv4c7drho42*drho42dr6 +&
      dv4c7drho234*drho234dr6
      
      
      
      
      F1= dv4c1dr1 + dv4c2dr1 + dv4c3dr1 + dv4c4dr1 + dv4c5dr1 + dv4c6dr1 + dv4c7dr1
      F2= dv4c1dr2 + dv4c2dr2 + dv4c3dr2 + dv4c4dr2 + dv4c5dr2 + dv4c6dr2 + dv4c7dr2
      F3= dv4c1dr3 + dv4c2dr3 + dv4c3dr3 + dv4c4dr3 + dv4c5dr3 + dv4c6dr3 + dv4c7dr3
      F4= dv4c1dr4 + dv4c2dr4 + dv4c3dr4 + dv4c4dr4 + dv4c5dr4 + dv4c6dr4 + dv4c7dr4
      F5= dv4c1dr5 + dv4c2dr5 + dv4c3dr5 + dv4c4dr5 + dv4c5dr5 + dv4c6dr5 + dv4c7dr5
      F6= dv4c1dr6 + dv4c2dr6 + dv4c3dr6 + dv4c4dr6 + dv4c5dr6 + dv4c6dr6 + dv4c7dr6
      
      END
      
      
      FUNCTION VOOD_10(R)
      !IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      implicit none
      real*8 de,a1,a2,a3,re,ra,x,VOOD_10,r
      
      DATA DE,A1,A2,A3,RE/0.155494053057397D0,6.000D0,10.075D0,8.714D0,1.21563D0/
      
      !como a1,a2 e a3 estao em angstrons mas o r que recebo da funcao esta em bohr
      !tenho de o converter para angstrons
      !o programa depois devolve o valor de vood em hartree(unidade do DE)
      !posso calcular esta funcao em angstrons porque
      !devolve vood em hartree(unidade do DE)
      
      ra=r*0.529177d0
      X=Ra-RE
      
      VOOD_10=(-DE*(1.0d0+A1*X+A2*X**2+A3*X**3)*DEXP(-A1*X))
      
      END
      
      !****************
      FUNCTION DVOOD_10(R)
      !****************
      
      implicit none
      real*8 de,a1,a2,a3,re,ra,x,VOOD_10,r,dvooddx,dxdra,DVOODDRa,dvood_10
      
      DATA DE,A1,A2,A3,RE/0.155494053057397D0,6.000D0,10.075D0,8.714D0,1.21563D0/
      
      ra=r*0.529177d0
      X=Ra-RE
      
      VOOD_10=(-DE*(1.0d0+A1*X+A2*X**2+A3*X**3)*DEXP(-A1*X))
      
      DVOODDX= -DE*(-A1**2*X*DEXP(-A1*X) + 2.0d0*A2*X*DEXP(-A1*X) - A2*A1*X**2*DEXP(-A1*X) +&
      3.0d0*A3*X**2*DEXP(-A1*X) - A3*A1*X**3*DEXP(-A1*X))
      
      DXDRa=1.0D0
      
      DVOODDRa= DVOODDX*DXDRa
      
      DVOOD_10=DVOODDRa*0.529177d0
      
      END
      
      
      
      !*******************************************************************
      !  To compute the electrostatic interaction between two OH diatomics
      !*******************************************************************
      
      FUNCTION Vdip_10(r1,r2,r3,r4,r5,r6)
      !     **************************************************************
      !  To compute the electrostatic interaction between two OH diatomics
      !
      IMPLICIT none
      real*8 vdipOHOH2,r1,r2,r3,r4,r5,r6
      real*8 dip1,dip2,FOHPdip_10,fdampdip_10
      real*8 qoa,qob,qha,qhb,vdipOHOH3,vdip_10
      
      !     **************************************************************
      !     ************V EM Eh***************
      dip1=FOHPdip_10(r3)
      dip2=FOHPdip_10(r4)
      qha=dip1/r3
      qoa=-qha
      qhb=dip2/r4
      qob=-qhb
      VdipOHOH2=qoa*qob/r1*fdampdip_10(r1)+qha*qhb/r2*fdampdip_10(r2)+  &
      qob*qha/r5*fdampdip_10(r5)+qoa*qhb/r6*fdampdip_10(r6)
      dip1=FOHPdip_10(r5)
      dip2=FOHPdip_10(r6)
      qha=dip1/r5
      qob=-qha
      qhb=dip2/r6
      qoa=-qhb
      VdipOHOH3=qoa*qob/r1*fdampdip_10(r1)+qha*qhb/r2*fdampdip_10(r2)+  &
      qoa*qha/r3*fdampdip_10(r3)+qob*qhb/r4*fdampdip_10(r4)
      Vdip_10=VdipOHOH2+VdipOHOH3
      
      RETURN
      END
      
      
      
      subroutine dVdip_10(r1,r2,r3,r4,r5,r6,dr1,dr2,dr3,dr4,dr5,dr6)
      !     **************************************************************
      !  To compute the electrostatic interaction between two OH diatomics
      !
      IMPLICIT none
      real*8 r1,r2,r3,r4,r5,r6
      real*8 dip1,dip2,FOHPdip_10,fdampdip_10
      real*8 dfdampdip_10,dFOHPdip_10
      real*8 qoa,qob,qha,qhb
      real*8 dr1,dr2,dr3,dr4,dr5,dr6
      real*8 dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6
      real*8 dv3dr1,dv3dr2,dv3dr3,dv3dr4,dv3dr5,dv3dr6
      real*8 dqhadr3,dqoadr3,dqhbdr4,dqobdr4
      real*8 dqhadr5, dqobdr5,dqhbdr6,dqoadr6
      
      !     **************************************************************
      !     ************V EM Eh***************
      dip1=FOHPdip_10(r3)
      dip2=FOHPdip_10(r4)
      qha=dip1/r3
      qoa=-qha
      qhb=dip2/r4
      qob=-qhb
      !      VdipOHOH2=qoa*qob/r1*fdampdip_10(r1)+qha*qhb/r2*fdampdip_10(r2)+  &
      !                qob*qha/r5*fdampdip_10(r5)+qoa*qhb/r6*fdampdip_10(r6)
      
      dqhadr3=dFOHPdip_10(r3)/r3-dip1/r3**2
      dqoadr3=-dqhadr3
      dqhbdr4=dFOHPdip_10(r4)/r4-dip2/r4**2
      dqobdr4=-dqhbdr4
      
      dv2dr1=-qoa*qob/r1**2*fdampdip_10(r1)+qoa*qob/r1*dfdampdip_10(r1)
      dv2dr2=-qha*qhb/r2**2*fdampdip_10(r2)+qha*qhb/r2*dfdampdip_10(r2)
      dv2dr3=dqoadr3*qob/r1*fdampdip_10(r1)+dqhadr3*qhb/r2*fdampdip_10(r2)+&
      dqhadr3*qob/r5*fdampdip_10(r5)+dqoadr3*qhb/r6*fdampdip_10(r6)
      dv2dr4=dqobdr4*qoa/r1*fdampdip_10(r1)+dqhbdr4*qha/r2*fdampdip_10(r2)+&
      dqobdr4*qha/r5*fdampdip_10(r5)+dqhbdr4*qoa/r6*fdampdip_10(r6)
      dv2dr5=-qob*qha/r5**2*fdampdip_10(r5)+qob*qha/r5*dfdampdip_10(r5)
      dv2dr6=-qoa*qhb/r6**2*fdampdip_10(r6)+qoa*qhb/r6*dfdampdip_10(r6)
      
      dip1=FOHPdip_10(r5)
      dip2=FOHPdip_10(r6)
      qha=dip1/r5
      qob=-qha
      qhb=dip2/r6
      qoa=-qhb
      !      VdipOHOH3=qoa*qob/r1*fdampdip_10(r1)+qha*qhb/r2*fdampdip_10(r2)+  &
      !                qoa*qha/r3*fdampdip_10(r3)+qob*qhb/r4*fdampdip_10(r4)
      
      
      dqhadr5=-dip1/r5**2+dFOHPdip_10(r5)/r5
      dqobdr5=-dqhadr5
      dqhbdr6=-dip2/r6**2+dFOHPdip_10(r6)/r6
      dqoadr6=-dqhbdr6
      
      
      dv3dr1=-qoa*qob/r1**2*fdampdip_10(r1)+qoa*qob/r1*dfdampdip_10(r1)
      dv3dr2=-qha*qhb/r2**2*fdampdip_10(r2)+qha*qhb/r2*dfdampdip_10(r2)
      dv3dr3=-qoa*qha/r3**2*fdampdip_10(r3)+qoa*qha/r3*dfdampdip_10(r3)
      dv3dr4=-qob*qhb/r4**2*fdampdip_10(r4)+qob*qhb/r4*dfdampdip_10(r4)
      dv3dr5=dqobdr5*qoa/r1*fdampdip_10(r1)+dqhadr5*qhb/r2*fdampdip_10(r2)+&
      dqhadr5*qoa/r3*fdampdip_10(r3)+dqobdr5*qhb/r4*fdampdip_10(r4)
      dv3dr6=dqoadr6*qob/r1*fdampdip_10(r1)+dqhbdr6*qha/r2*fdampdip_10(r2)+&
      dqoadr6*qha/r3*fdampdip_10(r3)+dqhbdr6*qob/r4*fdampdip_10(r4)
      
      
      !      Vdip_10=VdipOHOH2+VdipOHOH3
      
      dr1=dv2dr1+dv3dr1
      dr2=dv2dr2+dv3dr2
      dr3=dv2dr3+dv3dr3
      dr4=dv2dr4+dv3dr4
      dr5=dv2dr5+dv3dr5
      dr6=dv2dr6+dv3dr6
      
      RETURN
      END
      
      
      FUNCTION FOHPdip_10(R)
      !     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA C1,C2,C3,C4,C5,C6,C7/0.122307769D+01,-0.975516629D+00,   &
      0.447656609D+00,-0.928112515D-01,0.350125118D+00,  &
      0.371793971D-05,0.331141524D+03/
      !     ******************************************************************
      rr2=r*r
      rr3=rr2*r
      rr6=rr3*rr3
      rr8=rr6*rr2
      FOHPdip_10=(C1*R+C2*RR2+C3*RR3)*EXP(-c4*R-C5*RR2)+(1-EXP(-C6*   &
      RR8))*C7/RR6
      return
      end
      
      
      FUNCTION DFOHPdip_10(R)
      !     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA C1,C2,C3,C4,C5,C6,C7/0.122307769D+01,-0.975516629D+00,   &
      0.447656609D+00,-0.928112515D-01,0.350125118D+00,  &
      0.371793971D-05,0.331141524D+03/
      !     ******************************************************************
      rr2=r*r
      rr3=rr2*r
      rr6=rr3*rr3
      rr8=rr6*rr2
      
      DFOHPdip_10=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*R-C5*RR2)-(C4+   &
      2.0D0*C5*R)*(C1*R+C2*RR2+C3*RR3)*EXP(-C4*R-C5*RR2)  &
      -(6.0D0*C7/(RR6*R))*(1.0d0-EXP(-C6*RR8))+ &
      8.0D0*C6*C7*R*EXP(-C6*RR8)
      
      return
      end
      
      
      Function fdampdip_10(r)
      !     **************************************************************
      !  To damp the electrostatic interaction due to orbital overlap
      !
      IMPLICIT none
      real*8 fdampdip_10,r
      real*8 a(20),b(20),R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS
      real*8 rmohohp,r0ohohp,ro
      !
      DATA R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS/6.928203D0,   &
      1.40100D0,7.82D0,6.294894D0,6.334299622D0,1.8344D0,1.9086D0/
      !    **************************************************************
      !     AN E BN OBTIDOS A PARTIR DOS COEF.(A0,A1,B0,B1) DO ARTG. VAR82:857
      
      DATA A/ 25.9527760D+00,11.4004902D+00,7.0459516D+00   &
      ,0.500798749D+01,0.384282945D+01,0.309513333D+01   &
      ,0.257767037D+01,0.219990002D+01,0.191291419D+01   &
      ,0.168807142D+01,0.150753106D+01,0.135962478D+01   &
      ,0.123641324D+01,0.113231455D+01,0.104329456D+01   &
      ,0.966368248D+00,0.899281484D+00,0.840301611D+00   &
      ,0.788075808D+00,0.741533075D+00/
      DATA B/14.2790514D+00,12.9552683D+00,11.7542106D+00   &
      ,0.106645006D+02,0.967581549D+01,0.877878947D+01   &
      ,0.796492498D+01,0.722651228D+01,0.655655639D+01   &
      ,0.594871080D+01,0.539721740D+01,0.489685187D+01   &
      ,0.444287427D+01,0.403098404D+01,0.365727935D+01   &
      ,0.331822010D+01,0.301059437D+01,0.273148802D+01   &
      ,0.247825708D+01,0.224850269D+01/
      !    **************************************************************
      RMOHOHP=RMOHP*1.5d0
      R0OHOHP=r0OHP*1.5d0
      RO=0.5D0*(RMOHOHP+2.5D0*R0OHOHP)
      
      fdampdip_10=(1-exp(-A(3)*(r/ro)-B(3)*(r/ro)**2))
      
      return
      end
      
      
      Function Dfdampdip_10(r)
      !     **************************************************************
      !  To damp the electrostatic interaction due to orbital overlap
      !
      IMPLICIT none
      real*8 Dfdampdip_10,r
      real*8 a(20),b(20),R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS
      real*8 rmohohp,r0ohohp,ro
      !
      DATA R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS/6.928203D0,   &
      1.40100D0,7.82D0,6.294894D0,6.334299622D0,1.8344D0,1.9086D0/
      !    **************************************************************
      !     AN E BN OBTIDOS A PARTIR DOS COEF.(A0,A1,B0,B1) DO ARTG. VAR82:857
      
      DATA A/ 25.9527760D+00,11.4004902D+00,7.0459516D+00   &
      ,0.500798749D+01,0.384282945D+01,0.309513333D+01   &
      ,0.257767037D+01,0.219990002D+01,0.191291419D+01   &
      ,0.168807142D+01,0.150753106D+01,0.135962478D+01   &
      ,0.123641324D+01,0.113231455D+01,0.104329456D+01   &
      ,0.966368248D+00,0.899281484D+00,0.840301611D+00   &
      ,0.788075808D+00,0.741533075D+00/
      DATA B/14.2790514D+00,12.9552683D+00,11.7542106D+00   &
      ,0.106645006D+02,0.967581549D+01,0.877878947D+01   &
      ,0.796492498D+01,0.722651228D+01,0.655655639D+01   &
      ,0.594871080D+01,0.539721740D+01,0.489685187D+01   &
      ,0.444287427D+01,0.403098404D+01,0.365727935D+01   &
      ,0.331822010D+01,0.301059437D+01,0.273148802D+01   &
      ,0.247825708D+01,0.224850269D+01/
      !    **************************************************************
      RMOHOHP=RMOHP*1.5d0
      R0OHOHP=r0OHP*1.5d0
      RO=0.5D0*(RMOHOHP+2.5D0*R0OHOHP)
      
      Dfdampdip_10=(A(3)/ro+2.0d0*B(3)*(r/ro**2))*exp(-A(3)*(r/ro)-B(3)*  &
      (r/ro)**2)
      
      return
      end
      







