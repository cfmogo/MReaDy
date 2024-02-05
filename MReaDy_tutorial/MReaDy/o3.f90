        subroutine VO3_f7_switch(i,V,dV)
        use variables
        use constants

        double precision, dimension(12) :: dV1,dV
        integer      :: gi3,gi4,gi5
        integer      :: i
        double precision :: rxij,ryij,rzij,rij
        double precision :: r1xij,r1yij,r1zij
        double precision :: r2xij,r2yij,r2zij
        double precision :: r3xij,r3yij,r3zij
        double precision :: r1,r2,r3
        double precision :: V, V1,V2,dV2,V4,vdiat,dvdiat
        double precision :: f2,fs,dfsdr1,dfsdr2,dfsdr3
        double precision :: VO3_f7
        double precision :: dv1dr1,dv1dr2,dv1dr3


      dV1=0.0d0
      dV=0.0d0
      rxij=0.0d0
      ryij=0.0d0
      rzij=0.0d0
      rij=0.0d0
      r1xij=0.0d0
      r1yij=0.0d0
      r1zij=0.0d0
      r2xij=0.0d0
      r2yij=0.0d0
      r2zij=0.0d0
      r3xij=0.0d0
      r3yij=0.0d0
      r3zij=0.0d0
      r1 =0.0d0
      r2=0.0d0
      r3=0.0d0
      V =0.0d0
      V1=0.0d0
      V2=0.0d0
      dV2=0.0d0
      V4=0.0d0
      vdiat=0.0d0
      dvdiat=0.0d0
      dv1dr1=0.0d0
      dv1dr2=0.0d0
      dv1dr3=0.0d0



            gi3=group(i,3)      ! O atom
            gi4=group(i,4)      ! O atom
            gi5=group(i,5)      ! O atom

call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi4),ry(gi4),rz(gi4),rxij,ryij,rzij,rij)

            r1=rij
            r1xij=rxij
            r1yij=ryij
            r1zij=rzij

call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi5),ry(gi5),rz(gi5),rxij,ryij,rzij,rij)


            r2=rij
            r2xij=rxij
            r2yij=ryij
            r2zij=rzij

call mod_dist(rx(gi4),ry(gi4),rz(gi4),rx(gi5),ry(gi5),rz(gi5),rxij,ryij,rzij,rij)

            r3=rij
            r3xij=rxij
            r3yij=ryij
            r3zij=rzij

            V1=vo3_f7(r1,r2,r3)


           call dero3_f7(r1,r2,r3,dv1dr1,dv1dr2,dv1dr3)
           dv1(1)  =   + dv1dr1 * r1xij / r1 + dv1dr2 * r2xij / r2
           dv1(2)  =   + dv1dr1 * r1yij / r1 + dv1dr2 * r2yij / r2
           dv1(3)  =   + dv1dr1 * r1zij / r1 + dv1dr2 * r2zij / r2

           dv1(4)  =   - dv1dr1 * r1xij / r1 + dv1dr3 * r3xij / r3
           dv1(5)  =   - dv1dr1 * r1yij / r1 + dv1dr3 * r3yij / r3
           dv1(6)  =   - dv1dr1 * r1zij / r1 + dv1dr3 * r3zij / r3

           dv1(7)  =   - dv1dr2 * r2xij / r2 - dv1dr3 * r3xij / r3
           dv1(8)  =   - dv1dr2 * r2yij / r2 - dv1dr3 * r3yij / r3
           dv1(9)  =   - dv1dr2 * r2zij / r2 - dv1dr3 * r3zij / r3


            if ((r1.ge.in_bound).and.(r2.ge.in_bound).and. &
        (r3.le.in_bound)) then
!            print*,'118'


        ! O2 potential

        call potoop_f2(r3,V2,dV2)
       ! Switch function
        fs=f2(r1,r2)
        call df2(r1,r2,dfsdr1,dfsdr2)
       ! for the diatomic potencials

            call poto2q_f18(r1,Vdiat,dVdiat)



            dv(1) = dv(1) + fs * dVdiat*r1xij/r1
            dv(2) = dv(2) + fs * dVdiat*r1yij/r1
            dv(3) = dv(3) + fs * dVdiat*r1zij/r1

            dv(4) = dv(4) - fs * dVdiat*r1xij/r1
            dv(5) = dv(5) - fs * dVdiat*r1yij/r1
            dv(6) = dv(6) - fs * dVdiat*r1zij/r1


            V4 = V4 + Vdiat



            call poto2q_f18(r2,Vdiat,dVdiat)


            dv(1) = dv(1) + fs * dVdiat*r2xij/r2
            dv(2) = dv(2) + fs * dVdiat*r2yij/r2
            dv(3) = dv(3) + fs * dVdiat*r2zij/r2

            dv(7) = dv(7) - fs * dVdiat*r2xij/r2
            dv(8) = dv(8) - fs * dVdiat*r2yij/r2
            dv(9) = dv(9) - fs * dVdiat*r2zij/r2


            V4 = V4 + Vdiat



        V=V1*(1-fs)+(V2+V4)*fs

       dv=dv + dv1*(1-fs)

       dv(1)  =   dv(1)  + (v2+v4-v1) * ( dfsdr1 *  r1xij / r1 + dfsdr2 * r2xij / r2)
       dv(2)  =   dv(2)  + (v2+v4-v1) * ( dfsdr1 *  r1yij / r1 + dfsdr2 * r2yij / r2)
       dv(3)  =   dv(3)  + (v2+v4-v1) * ( dfsdr1 *  r1zij / r1 + dfsdr2 * r2zij / r2)

       dv(4)  =   dv(4)  + (v2+v4-v1) * (-dfsdr1 *  r1xij / r1)
       dv(5)  =   dv(5)  + (v2+v4-v1) * (-dfsdr1 *  r1yij / r1)
       dv(6)  =   dv(6)  + (v2+v4-v1) * (-dfsdr1 *  r1zij / r1)

       dv(7)  =   dv(7)  + (v2+v4-v1) * (-dfsdr2 *  r2xij / r2)
       dv(8)  =   dv(8)  + (v2+v4-v1) * (-dfsdr2 *  r2yij / r2)
       dv(9)  =   dv(9)  + (v2+v4-v1) * (-dfsdr2 *  r2zij / r2)


       dv(4)  =   dv(4)  + fs * ( dv2 * r3xij / r3)
       dv(5)  =   dv(5)  + fs * ( dv2 * r3yij / r3)
       dv(6)  =   dv(6)  + fs * ( dv2 * r3zij / r3)

       dv(7)  =   dv(7)  + fs * (-dv2 * r3xij / r3)
       dv(8)  =   dv(8)  + fs * (-dv2 * r3yij / r3)
       dv(9)  =   dv(9)  + fs * (-dv2 * r3zij / r3)





       else if ((r1.ge.in_bound).and.(r3.ge.in_bound).and. &
        (r2.le.in_bound)) then
!        print*,'208'


        ! O2 potential

        call potoop_f2(r2,V2,dV2)


       ! Switch function
        fs=f2(r1,r3)
        call df2(r1,r3,dfsdr1,dfsdr3)

        ! for the diatomic potencials

            call poto2q_f18(r1,Vdiat,dVdiat)



            dv(1) = dv(1) + fs * dVdiat*r1xij/r1
            dv(2) = dv(2) + fs * dVdiat*r1yij/r1
            dv(3) = dv(3) + fs * dVdiat*r1zij/r1

            dv(4) = dv(4) - fs * dVdiat*r1xij/r1
            dv(5) = dv(5) - fs * dVdiat*r1yij/r1
            dv(6) = dv(6) - fs * dVdiat*r1zij/r1


            V4 = V4 + Vdiat



            call poto2q_f18(r3,Vdiat,dVdiat)


            dv(4) = dv(4) + fs * dVdiat*r2xij/r2
            dv(5) = dv(5) + fs * dVdiat*r2yij/r2
            dv(6) = dv(6) + fs * dVdiat*r2zij/r2

            dv(7) = dv(7) - fs * dVdiat*r2xij/r2
            dv(8) = dv(8) - fs * dVdiat*r2yij/r2
            dv(9) = dv(9) - fs * dVdiat*r2zij/r2


            V4 = V4 + Vdiat



        V=V1*(1-fs)+(V2+V4)*fs

       dv=dv + dv1*(1-fs)

       dv(1)  =   dv(1)  + (v2+v4-v1) * ( dfsdr1 *  r1xij / r1)
       dv(2)  =   dv(2)  + (v2+v4-v1) * ( dfsdr1 *  r1yij / r1)
       dv(3)  =   dv(3)  + (v2+v4-v1) * ( dfsdr1 *  r1zij / r1)

       dv(4)  =   dv(4)  + (v2+v4-v1) * (-dfsdr1 *  r1xij / r1 + dfsdr3 *  r3xij / r3)
       dv(5)  =   dv(5)  + (v2+v4-v1) * (-dfsdr1 *  r1yij / r1 + dfsdr3 *  r3yij / r3)
       dv(6)  =   dv(6)  + (v2+v4-v1) * (-dfsdr1 *  r1zij / r1 + dfsdr3 *  r3zij / r3)

       dv(7)  =   dv(7)  + (v2+v4-v1) * (-dfsdr3 *  r3xij / r3)
       dv(8)  =   dv(8)  + (v2+v4-v1) * (-dfsdr3 *  r3yij / r3)
       dv(9)  =   dv(9)  + (v2+v4-v1) * (-dfsdr3 *  r3zij / r3)


       dv(1)  =   dv(1)  + fs * ( dv2 * r2xij / r2)
       dv(2)  =   dv(2)  + fs * ( dv2 * r2yij / r2)
       dv(3)  =   dv(3)  + fs * ( dv2 * r2zij / r2)

       dv(7)  =   dv(7)  + fs * (-dv2 * r2xij / r2)
       dv(8)  =   dv(8)  + fs * (-dv2 * r2yij / r2)
       dv(9)  =   dv(9)  + fs * (-dv2 * r2zij / r2)




       else if ((r2.ge.in_bound).and.(r3.ge.in_bound).and. &
        (r1.le.in_bound)) then
!        print*,'in  297'


        ! O2 potential

        call potoop_f2(r1,V2,dV2)
       ! Switch function
        fs=f2(r2,r3)
        call df2(r2,r3,dfsdr2,dfsdr3)



        ! for the diatomic potencials

            call poto2q_f18(r2,Vdiat,dVdiat)



            dv(1) = dv(1) + fs * dVdiat*r2xij/r2
            dv(2) = dv(2) + fs * dVdiat*r2yij/r2
            dv(3) = dv(3) + fs * dVdiat*r2zij/r2

            dv(7) = dv(7) - fs * dVdiat*r2xij/r2
            dv(8) = dv(8) - fs * dVdiat*r2yij/r2
            dv(9) = dv(9) - fs * dVdiat*r2zij/r2


            V4 = V4 + Vdiat



            call poto2q_f18(r3,Vdiat,dVdiat)


            dv(4) = dv(4) + fs * dVdiat*r3xij/r3
            dv(5) = dv(5) + fs * dVdiat*r3yij/r3
            dv(6) = dv(6) + fs * dVdiat*r3zij/r3

            dv(7) = dv(7) - fs * dVdiat*r3xij/r3
            dv(8) = dv(8) - fs * dVdiat*r3yij/r3
            dv(9) = dv(9) - fs * dVdiat*r3zij/r3


            V4 = V4 + Vdiat



        V=V1*(1-fs)+(V2+V4)*fs


       dv=dv + dv1*(1-fs)

       dv(1)  =   dv(1)  + (v2+v4-v1) * ( dfsdr2 *  r2xij / r2)
       dv(2)  =   dv(2)  + (v2+v4-v1) * ( dfsdr2 *  r2yij / r2)
       dv(3)  =   dv(3)  + (v2+v4-v1) * ( dfsdr2 *  r2zij / r2)

       dv(4)  =   dv(4)  + (v2+v4-v1) * ( dfsdr3 *  r3xij / r3)
       dv(5)  =   dv(5)  + (v2+v4-v1) * ( dfsdr3 *  r3yij / r3)
       dv(6)  =   dv(6)  + (v2+v4-v1) * ( dfsdr3 *  r3zij / r3)

       dv(7)  =   dv(7)  + (v2+v4-v1) * (-dfsdr2 *  r2xij / r2- dfsdr3 *  r3xij / r3)
       dv(8)  =   dv(8)  + (v2+v4-v1) * (-dfsdr2 *  r2yij / r2- dfsdr3 *  r3yij / r3)
       dv(9)  =   dv(9)  + (v2+v4-v1) * (-dfsdr2 *  r2zij / r2- dfsdr3 *  r3zij / r3)


       dv(1)  =   dv(1)  + fs * ( dv2 * r1xij / r1)
       dv(2)  =   dv(2)  + fs * ( dv2 * r1yij / r1)
       dv(3)  =   dv(3)  + fs * ( dv2 * r1zij / r1)

       dv(4)  =   dv(4)  + fs * (-dv2 * r1xij / r1)
       dv(5)  =   dv(5)  + fs * (-dv2 * r1yij / r1)
       dv(6)  =   dv(6)  + fs * (-dv2 * r1zij / r1)






         else!  end if

         v=v1

          dv  = dv1

          group(i,9)=0

      end if



      end subroutine









      FUNCTION VO3_f7(R1in,R2in,R3in)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION R(3)

      r1=r1in/5.2917721092d-11
      r2=r2in/5.2917721092d-11
      r3=r3in/5.2917721092d-11

      R(1)=R1
      R(2)=R2
      R(3)=R3
      CALL PRINC_f7(R,YR)
      VO3_f7=(YR+VAD_f7(R1,R2,R3,Q1,Q2,Q3)+ &
       PTH_f7(R1,R2,R3,Q1,Q2,Q3))*4.3597482D-18
      RETURN
      END





      SUBROUTINE PRINC_f7(RR,Y)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      double precision  :: Rout(3)
      DIMENSION RR(3)


      VR2=0.0
      DO 5 I=1,3
      Rout(I)=RR(I)
    5 CONTINUE
      CALL CE_f7(Rout,VDP)
      DO 10 I=1,3
      RCSN=Rout(I)
      CALL V2_f7(RCSN,VTOT)
      VR2=VR2+VTOT
   10 CONTINUE
      CALL CEQ_f7(Rout,VEQ)
      Y=VDP+VR2+VEQ
      RETURN
      END


      SUBROUTINE DPRINC_f7(RR,DY,K1D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      double precision :: Ro(3)
      DIMENSION RR(3)
      K1=K1D
      DO 5 I=1,3
      Ro(I)=RR(I)
    5 CONTINUE
      CALL DCE_f7(Ro,DVDP,K1)
      RCSN=Ro(K1D)
      CALL DV2_f7(RCSN,DVTOT)
      CALL DCEQ_f7(Ro,DVEQ,K1)
      DY=DVDP+DVTOT+DVEQ
      RETURN
      END


      SUBROUTINE V2_f7(Rin,VTOTout)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)


      COMMON/DM4_f7/A1,A2,A3,GM,RE,DE


      double precision, intent(out) :: VTOTout
      double precision, intent(in) :: Rin


      RR=Rin-RE
      EGM=EXP(-GM*RR)
      POL=DE*(1.+A1*RR+A2*RR**2+A3*RR**3)
      VTOTout=POL*EGM
      RETURN
      END


      SUBROUTINE DV2_f7(Ri,DVTOTout)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/DM4_f7/A1,A2,A3,GM,RE,DE


      double precision, intent(out) :: DVTOTout
      double precision, intent(in) :: Ri

      RR=Ri-RE
      EGM=EXP(-GM*RR)
      POL=DE*(1.+A1*RR+A2*RR**2+A3*RR**3)
      DPOL=DE*(A1+2.*A2*RR+3.*A3*RR**2)
      DVTOTout=(DPOL*EGM-GM*POL*EGM)
      RETURN
      END




      SUBROUTINE CE_f7(Rin,VDPout)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/MAIN_f7/C2C(3),A0(3),A1(3),DD
      COMMON/DUMMY_f7/INDEX_f7(3,2)

      double precision,intent(in) :: Rin(3)
      double precision,intent(out) :: VDPout

      VDPout=0.0d0
      DO 10 J=1,3
      I1=INDEX_f7(J,1)
      I2=INDEX_f7(J,2)
      Aout=Rin(I1)
      Bout=Rin(I2)
      DO 10 I=1,3
      CXX=C2C(I)
      AK0=A0(I)
      AK1=A1(I)
      CALL CF_f7(CCo,Aout,Bout,CXX,AK0,AK1)
      VDPout=VDPout+CCo*DMR_f7(I,Rin(J))
   10 CONTINUE
      RETURN
      END


      SUBROUTINE DCE_f7(Ri,DVDPout,K1Din)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)

      COMMON/MAIN_f7/C2C(3),A0(3),A1(3),DD
      COMMON/DUMMY_f7/INDEX_f7(3,2)
      COMMON/DMY_f7/IX(3,3)

      double precision, intent(in) :: Ri(3)
      double precision, intent(out) :: DVDPout

      DVDPout=0.0d0
      DO 10 I=1,3
      J=K1Din
!     INDEX_f7:DADA UMA INTERACAO DEFINE AS OUTRAS DUAS
      I1=INDEX_f7(J,1)
      I2=INDEX_f7(J,2)
      A=Ri(I1)
      B=Ri(I2)
      CXX=C2C(I)
      AK0=A0(I)
      AK1=A1(I)
      CALL CF_f7(CCo,A,B,CXX,AK0,AK1)
      DVDPout=DVDPout+CCo*DDMR_f7(I,Ri(J))
      DO 10 K=1,2
      J=INDEX_f7(K1Din,K)
      IO=IX(K1Din,J)
      A=Ri(K1Din)
      B=Ri(IO)
      CALL DCF_f7(DCC,A,B,CXX,AK0,AK1)
      DVDPout=DVDPout+DCC*DMR_f7(I,Ri(J))
   10 CONTINUE
      RETURN
      END





      SUBROUTINE CF_f7(CCout,Ain,Bin,CXX,AK0,AK1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      common/rzero/R0

      TANHA=TANH(AK1*Ain)
      TANHB=TANH(AK1*Bin)
      EA=AK0*EXP(-AK1*(Ain-R0))
      EB=AK0*EXP(-AK1*(Bin-R0))
      CCout=0.5*CXX*((1.+EA)*TANHB+(1.+EB)*TANHA)
      RETURN
      END


      SUBROUTINE DCF_f7(DCCout,Ai,Bi,CXX,AK0,AK1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      common/rzero/R0

      EA=AK0*EXP(-AK1*(Ai-R0))
      EB=AK0*EXP(-AK1*(Bi-R0))
      AKA=AK1*Ai
      TANHA=TANH(AKA)
      TANHB=TANH(AK1*Bi)
      SINHA=SINH(AKA)
      COSHA=SINHA/TANHA
      CSHA2=1./COSHA*1./COSHA
      DCCout=0.5*CXX*AK1*(-EA*TANHB+(1.+EB)*CSHA2)
      RETURN
      END




      FUNCTION DMR_f7(I,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
!    I=1:C6;I=2:C8;I=3:C8;I=4:C5

      COMMON/PN_f7/N(4)
      COMMON/DM1_f7/D1(4),D2(4)
      COMMON/DM2_f7/R0,RM
      COMMON/DM3_f7/GAMA

      RO=(RM+GAMA*R0)/2.
      X=R/RO
      DX=1./RO
      POL=-D1(I)*X*(1.+D2(I)*X)
      EXPO=EXP(POL)
      IDXL=N(I)
      DXL=FLOAT(IDXL)
      C=(1.-EXPO)**IDXL
      RL=R**IDXL
      DMR_f7=-C/RL
      RETURN
      END

      FUNCTION DDMR_f7(I,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
!    I=1:C6;I=2:C8;I=3:C8;I=4:C5

      COMMON/PN_f7/N(4)
      COMMON/DM1_f7/D1(4),D2(4)
      COMMON/DM2_f7/R0,RM
      COMMON/DM3_f7/GAMA

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
      DDMR_f7=-(DC/RL-DXL*C/(RL*R))
      RETURN
      END

      SUBROUTINE CEQ_f7(Rimp,VEQout)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      common/rzero/R0
      COMMON/DUMMY_f7/INDEX_f7(3,2)

      double precision, intent(in) :: Rimp(3)
      double precision, intent(out) :: VEQout

      VEQout=0.0
      DO 10 J=1,3
      I1=INDEX_f7(J,1)
      I2=INDEX_f7(J,2)
      Ao=Rimp(I1)
      Bo=Rimp(I2)
      C=Rimp(J)
      CALL CFQ_f7(CCi,Ao,Bo)
      VEQout=VEQout+CCi*DMR_f7(4,Rimp(J))
   10 CONTINUE
      RETURN
      END


      SUBROUTINE DCEQ_f7(Rim,DVEQout,K1Di)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/DUMMY_f7/INDEX_f7(3,2)
      COMMON/DMY_f7/IX(3,3)

      double precision, intent(in) :: Rim(3)
      double precision, intent(out) :: DVEQout

      DVEQout=0.0
      J=K1Di
!     INDEX_f7:DADA UMA INTERACAO DEFINE AS OUTRAS DUAS
      I1=INDEX_f7(J,1)
      I2=INDEX_f7(J,2)
      Ao=Rim(I1)
      Bo=Rim(I2)
      C=Rim(J)
      CALL CFQ_f7(CCin,Ao,Bo)
      DVEQout=DVEQout+CCin*DDMR_f7(4,Rim(J))
      DO 10 K=1,2
      J=INDEX_f7(K1Di,K)
      IO=IX(K1Di,J)
      A=Rim(K1Di)
      B=Rim(IO)
      C=Rim(J)
      CALL DCFQ_f7(DCC,A,B)
      DVEQout=DVEQout+DCC*DMR_f7(4,Rim(J))
   10 CONTINUE
      RETURN
      END



      SUBROUTINE CFQ_f7(CCout,Ai,Bi)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/MAIN2_f7/C5Q,RK
      common/rzero/R0

      CCA=Ai**4/R0**4
      CCB=Bi**4/R0**4
      TANHA=TANH(RK*Ai)
      TANHB=TANH(RK*Bi)
      EA=EXP(-RK*(Ai-R0))
      EB=EXP(-RK*(Bi-R0))
      CCout=0.5*C5Q*(CCA*EA*TANHB+CCB*EB*TANHA)
      RETURN
      END


      SUBROUTINE DCFQ_f7(DCCo,Ain,Bin)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/MAIN2_f7/C5Q,RK
      common/rzero/R0

      RKA=RK*Ain
      CCA=Ain**4/R0**4
      DCCA=4.*Ain**3/R0**4
      CCB=Bin**4/R0**4
      EA=EXP(-RK*(Ain-R0))
      EB=EXP(-RK*(Bin-R0))
      TANHA=TANH(RKA)
      TANHB=TANH(RK*Bin)
      SINHA=SINH(RKA)
      COSHA=SINHA/TANHA
      CSHA2=1.0D0/COSHA*1.0D0/COSHA
      DCCo=0.5D0*C5Q*(-CCA*EA*TANHB*RK+DCCA*EA*TANHB+RK*CCB*EB*CSHA2)
      RETURN
      END





      FUNCTION VAD_f7(R1,R2,R3,Q1o,Q2o,Q3o)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/DM8_f7/DRDQ(3,3),DQDR(3,3)
      COMMON/VAR_f7/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
      COMMON/MAIN_f7/V(3),AK(3),BK(3),DD

      S1=R1-DD
      S2=R2-DD
      S3=R3-DD
      Q1o=S1*DQDR(1,1)+S2*DQDR(1,2)+S3*DQDR(1,3)
      Q2o=S1*DQDR(2,1)+S2*DQDR(2,2)+S3*DQDR(2,3)
      Q3o=S1*DQDR(3,1)+S2*DQDR(3,2)+S3*DQDR(3,3)
      PART=SUMAT_f7(R1,R2,R3)
      TH=1.-TANH(GAMA*Q1o/4.0)
      VAD_f7=(B1+B2*Q1o+(AL(1)+AL(2)*(Q2o**2+Q3o**2)+AL(3)* &
      (Q2o**2+Q3o**2)**2+ &
      AL4*(Q2o**2+Q3o**2))*(1.5+PART)**5)*(1.2527+PART)**5 &
      *EXP(-BETAD*(Q2o**2+Q3o**2))*TH
      RETURN
      END

      FUNCTION SUMAT_f7(R1,R2,R3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      SUMAT_f7=(R1**2-R2**2-R3**2)/(2.*R2*R3)+ &
            (R2**2-R3**2-R1**2)/(2.*R1*R3)+ &
            (R3**2-R1**2-R2**2)/(2.*R1*R2)
      RETURN
      END


      FUNCTION PTH_f7(R1,R2,R3,Q1i,Q2i,Q3i)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/VAR_f7/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
      COMMON/CONST2_f7/A,C3,C6,C10
      TH2=1.-TANH(GAMA*Q1i/2.)
      P=A+C1*Q1i+C2*Q1i**2+C3*(Q2i**2+Q3i**2)+C5*Q1i* &
       (Q2i**2+Q3i**2)+C6*(Q3i**3-3.*Q2i**2*Q3i)+C10* &
       (Q2i**2+Q3i**2)**2

      PTH_f7=P*TH2
      RETURN
      END




      FUNCTION DVAD_f7(R1,R2,R3,I,Q1o,Q2o,Q3o)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/DM8_f7/DRDQ(3,3),DQDR(3,3)
      COMMON/VAR_f7/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
      COMMON/MAIN_f7/V(3),AK(3),BK(3),DD

      S1=R1-DD
      S2=R2-DD
      S3=R3-DD
      Q1o=S1*DQDR(1,1)+S2*DQDR(1,2)+S3*DQDR(1,3)
      Q2o=S1*DQDR(2,1)+S2*DQDR(2,2)+S3*DQDR(2,3)
      Q3o=S1*DQDR(3,1)+S2*DQDR(3,2)+S3*DQDR(3,3)
      TH=1.-TANH(GAMA*Q1o/4.0)
      DTH=-GAMA/(4.0*COSH(GAMA*Q1o/4.0)**2)*DQDR(1,I)
      PART=SUMAT_f7(R1,R2,R3)
      DPART=DSUMAT_f7(R1,R2,R3,I)
      T1=B1+B2*Q1o
      T2=(AL(1)+AL(2)*(Q2o**2+Q3o**2)+AL(3)*(Q2o**2+Q3o**2)**2)
      T3=(1.5+PART)**5
      T4=(1.2527+PART)**5
      T5=EXP(-BETAD*(Q2o**2+Q3o**2))
      T6=TH
      DT1=B2*DQDR(1,I)
      DT2=2.*AL(2)*(Q2o*DQDR(2,I)+Q3o*DQDR(3,I)) &
         +4.*AL(3)*((Q2o**3+Q2o*Q3o**2)*DQDR(2,I) &
         +(Q3o**3+Q2o**2*Q3o)*DQDR(3,I))
      DT3=5.*(1.5000+PART)**4*DPART
      DT4=5.*(1.2527+PART)**4*DPART
      DT5=-2.*BETAD*(Q2o*DQDR(2,I)+Q3o*DQDR(3,I))*T5
      DT6=DTH
      DVAD_f7=(DT1+DT2*T3+T2*DT3)*T4*T5*T6 &
          +(T1+T2*T3)*DT4*T5*T6 &
          +(T1+T2*T3)*T4*DT5*T6 &
          +(T1+T2*T3)*T4*T5*DT6
      RETURN
      END


      FUNCTION DSUMAT_f7( R1,R2,R3,I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/DUMMY_f7/INDEX_f7(3,2)

      DIMENSION R(3)
      R(1)=R1
      R(2)=R2
      R(3)=R3
      I1=INDEX_f7(I,1)
      I2=INDEX_f7(I,2)
      DSUMAT_f7=R(I)/(R(I1)*R(I2)) &
          +(R(I2)**3-R(I1)**2*R(I2)-R(I)**2*R(I2))/(2.*R(I)**2*R(I2)**2) &
          +(R(I1)**3-R(I2)**2*R(I1)-R(I)**2*R(I1))/(2.*R(I)**2*R(I1)**2)
      RETURN
      END

      FUNCTION DPTH_f7(R1,R2,R3,I,Q1i,Q2i,Q3i)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/DM8_f7/DRDQ(3,3),DQDR(3,3)
      COMMON/VAR_f7/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
      COMMON/CONST2_f7/A,C3,C6,C10
      DTH2=-GAMA/(2.*COSH(GAMA*Q1i/2.)**2)*DQDR(1,I)
      TH2=1.-TANH(GAMA*Q1i/2.)
      P=A+C1*Q1i+C2*Q1i**2+C3*(Q2i**2+Q3i**2)+C5*Q1i* &
        (Q2i**2+Q3i**2)+C6*(Q3i**3-3.*Q2i**2*Q3i)+C10* &
        (Q2i**2+Q3i**2)**2
      DPVTH=DTH2*P
      DP=C1*DQDR(1,I)+2.*C2*Q1i*DQDR(1,I)+2.*C3*(Q2i*DQDR(2,I) &
         +Q3i*DQDR(3,I))+C5*DQDR(1,I)* &
         (Q2i**2+Q3i**2)+C5*Q1i*(2.*Q2i &
         *DQDR(2,I)+2.*Q3i*DQDR(3,I))+C6* &
         (3.*Q3i**2*DQDR(3,I)-6.*Q2i*Q3i &
         *DQDR(2,I)-3.*Q2i**2*DQDR(3,I))+C10* &
         ((4.*Q2i**3+4.*Q2i*Q3i**2) &
         *DQDR(2,I)+(4.*Q3i**3+4.*Q2i**2*Q3i)*DQDR(3,I))
      DPTH_f7=DPVTH+DP*TH2
      RETURN
      END








!  *********************************************************************
!  DERIVATIVES OF THE POTENTIAL ENERGY FUNCTION FOR O3
!  *********************************************************************
      subroutine DERO3_f7(R1in,R2in,R3in,dvr1,dvr2,dvr3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION R(3)
        con=4.3597482D-18/5.2917721092d-11
       r1=r1in/5.2917721092d-11
       r2=r2in/5.2917721092d-11
       r3=r3in/5.2917721092d-11

      R(1)=R1
      R(2)=R2
      R(3)=R3
        call DPRINC_f7(R,DY,1)
        DVR1=(DY+DVAD_f7(R1,R2,R3,1,Q1,Q2,Q3)+ &
        DPTH_f7(R1,R2,R3,1,Q1,Q2,Q3))*con
        call DPRINC_f7(R,DY,2)
        DVR2=(DY+DVAD_f7(R1,R2,R3,2,Q1,Q2,Q3)+ &
        DPTH_f7(R1,R2,R3,2,Q1,Q2,Q3))*con
        call DPRINC_f7(R,DY,3)
        DVR3=(DY+DVAD_f7(R1,R2,R3,3,Q1,Q2,Q3)+ &
        DPTH_f7(R1,R2,R3,3,Q1,Q2,Q3))*con
      RETURN
      END






!     *****************************************************************
!     DATA FOR O3 SURFACE
!     ***************************************************************
      BLOCK DATA A1aa
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER(I-N)
      COMMON/MAIN_f7/VW(3),AK0(3),AK1(3),DD
      COMMON/MAIN2_f7/C5Q,RK
      COMMON/PN_f7/N(4)
      COMMON/DM1_f7/D1(4),D2(4)
      COMMON/DUMMY_f7/INDEX_f7(3,2)
      COMMON/DMY_f7/IX(3,3)
      COMMON/DM2_f7/R0,RM
      COMMON/DM3_f7/GAMA
      COMMON/DM4_f7/A1,A2,A3,GM,RE,DE
      COMMON/DM8_f7/DRDQ(3,3),DQDR(3,3)
      COMMON/VAR_f7/C1,C2,C5,GAMAP,B1,B2,AL(3),AL4,BETAD
      COMMON/CONST2_f7/A,C3,C6,C10

!  ----------------------------------------------------------------------
      common/damps_f7/dd1(3),dd2(3)
      common/rzero/R0PASS




      DATA (DD1(I),I=1,3)/.46547873E1,.38038890E1,.32525512E1/
      DATA (DD2(I),I=1,3)/.96802180E1,.79933059E1,.66003615E1/
!  ----------------------------------------------------------------------
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
      DATA ((INDEX_f7(I,J),I=1,3),J=1,2)/2,1,1,3,3,2/
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

