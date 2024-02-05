      subroutine pothht_f14(rin,VSout,dVSout)
        use phys_parameters

      implicit none


      double precision, intent(in)  :: rin
      double precision, intent(out) :: VSout,dVSout
      double precision :: VHHT_f14,DVHHT_f14,R
      double precision :: fs,dfs,f1_lr,df1_lr


      R=Rin/5.2917721092d-11
      VSout=VHHT_f14(R)*4.3597482D-18
      dVSout=DVHHT_f14(R)*4.3597482D-18/5.2917721092d-11

      if (Rin.gt.(rl_cut_off-0.2d-9)) then
        fs=f1_lr(Rin)
        dfs=df1_lr(Rin)

        VSout=VSout*(1.0d0-fs)
        dVSout=dVSout*(1.0d0-fs)-VSout*dfs
      end if

      end subroutine pothht_f14










      FUNCTION VHHT_f14(R)
!    **************************************************************
!     TO COMPUTE THE HFACE FOR triplet H...H
!    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!    **************************************************************
!     O AJUSTE - VHH(3SI) obteve-se PARTIR DOS PONTOS DE
!                       dois pontos dos nossos para 0.22 e 0.3
!                                    - KOL65:2429 (R=1.0<>5.9)
!                                    - KOL74:457  (R=6.0<>12.0)
!    **************************************************************
      VHHT_f14=EHFHHT_f14(R)+DISHHT_f14(R)

      RETURN
      END

      FUNCTION EHFHHT_f14(R)
!    **************************************************************
!     TO COMPUTE THE EHF FOR  triplet H...H
!    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION E(6)
!    **************************************************************
!     coef retirados de var82:857 pag 862
      DATA D,E/2.64103167D-05,-1.4993795d+00,2.9009680d+00, &
             4.7898562d-01,2.7407202d-02,-9.0948301d-01, &
             -2.0122436d-01/
      DATA AGEX,ALPHEX,GAMEX/-0.8205D0,2.5D0,2.0D0/
      COMMON/AN_f14/A(20)
      COMMON/BN_f14/B(20)
      COMMON/DIATDI_f14/R0HH,RMHHT
!      COMMON/CEHFHT_f14/R1,R2,R3,RHO,RR,EX,D6,TEXP,POL,EHF


!    **************************************************************
      R1=R-RMHHT
      R2=R1*R1
      R3=R2*R1

      RHO=(rmHHT+2.5D0*R0HH)/2.0d0
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6
      TEXP=exp(-E(2)*R1-E(3)*R2-E(4)*R3)
      POL=1.0d0-E(1)*R1-E(5)*R2-E(6)*R3
      EHF=-D/R*POL*TEXP
      EX=-agex*R**ALPHEX*EXP(-GAMEX*R)*D6
      EHFHHT_f14=EHF+EX


      RETURN
      END

      FUNCTION DISHHT_f14(R)
!    **************************************************************
!     TO COMPUTE THE DISPERSION FOR triplet H...H
!    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/AN_f14/A(20)
      COMMON/BN_f14/B(20)
      COMMON/DISPC1_f14/CHH(16)
      COMMON/DIATDI_f14/R0HH,RMHHT
!    **************************************************************
      R6=R**6
      R8=R6*R*R
      ADISP=-CHH(6)*DAMPHHT_f14(6,R)/R6-CHH(8)*DAMPHHT_f14(8,R)/R8
      DO N=10,16
           ADISP=ADISP-CHH(N)*DAMPHHT_f14(N,R)/(R8*R**(N-8.0D0))
      ENDDO
      DISHHT_f14=ADISP
      RETURN
      END

      FUNCTION DAMPHHT_f14(N,R)
!    ***************************************************************
!     CALCULATES VARANDAS-BRANDAO DAMPING FUNCTIONS
!    ***************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DIATDI_f14/R0HH,RMHHT
      COMMON/AN_f14/A(20)
      COMMON/BN_f14/B(20)
!    ***************************************************************
      RR=2.0D0*R/(RMHHT+2.5D0*R0HH)
      DAMPHHT_f14=(1.0D0-EXP(-A(N)*RR-B(N)*RR*RR))**N
      RETURN
      END

      FUNCTION DVHHT_f14(R)
!    **************************************************************
!     TO COMPUTE THE DER. OF THE HFACE FOR  triplet H...H
!    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!    **************************************************************
      DVHHT_f14=DEHFHHT_f14(R)+DDISHHT_f14(R)


      RETURN
      END

      FUNCTION DEHFHHT_f14(R)
!    **************************************************************
!     TO COMPUTE THE EHF FOR  triplet H...H
!    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION E(6)
!    **************************************************************
      DATA D,E/2.64103167D-05,-1.4993795d+00,2.9009680d+00, &
             4.7898562d-01,2.7407202d-02,-9.0948301d-01, &
             -2.0122436d-01/
      DATA AGEX,ALPHEX,GAMEX/-0.8205D0,2.5D0,2.0D0/
      COMMON/DIATDI_f14/R0HH,RMHHT
      COMMON/AN_f14/A(20)
      COMMON/BN_f14/B(20)
!      COMMON/CEHFHT_f14/R1,R2,R3,RHO,RR,EX,D6,TEXP,POL,EHF

!    **************************************************************

      R1=R-RMHHT
      R2=R1*R1
      R3=R2*R1

      RHO=(rmHHT+2.5D0*R0HH)/2.0d0
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6
      TEXP=exp(-E(2)*R1-E(3)*R2-E(4)*R3)
      POL=1.0d0-E(1)*R1-E(5)*R2-E(6)*R3
      EHF=-D/R*POL*TEXP
      EX=-agex*R**ALPHEX*EXP(-GAMEX*R)*D6






      rr2=rr*rr
      DBASE=(A(6)+2.0d0*B(6)*RR)/RHO*EXP(-A(6)*RR-B(6)*RR2)
      DD6=6.0d0*DBASE*(1.0d0-EXP(-A(6)*RR-B(6)*RR2))**5
      FEXC=EX/D6
      DFEXC=-agex*alphex*r**(ALPHEX-1.0d0)*EXP(-GAMEX*R)-GAMEX* &
             FEXC
      DEXCHHT=DFEXC*D6+DD6*FEXC
      DEHFHHT_f14=D/R**2*POL*Texp+E(1)*D/R*Texp+2.0d0*E(5)*r1*D/r*Texp+ &
             3.0d0*E(6)*r2*D/r*Texp-D/R*POL*(-E(2)-2.0D0* &
              E(3)*R1-3.0d0*E(4)*R2)*Texp
      DEHFHHT_f14=DEHFHHT_f14+DEXCHHT
      RETURN
      END

      FUNCTION DDISHHT_f14(R)
!    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC1_f14/CHH(16)
      COMMON/DIATDI_f14/R0HH,RMHHT
      COMMON/AN_f14/A(20)
      COMMON/BN_f14/B(20)
!    **************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      R11=R10*R
      R12=R11*R
      R13=R12*R
      R14=R13*R
      R15=R14*R
      R16=R15*R
      RR=2.0D0*R/(RMHHT+2.5D0*R0HH)
      rr2=rr*RR
      DRR=RR/R
      T6=1.0D0-EXP(-A(6)*RR-B(6)*RR2)
      T8=1.0D0-EXP(-A(8)*RR-B(8)*RR2)
      T10=1.0D0-EXP(-A(10)*RR-B(10)*RR2)
      T11=1.0D0-EXP(-A(11)*RR-B(11)*RR2)
      T12=1.0D0-EXP(-A(12)*RR-B(12)*RR2)
      T13=1.0D0-EXP(-A(13)*RR-B(13)*RR2)
      T14=1.0D0-EXP(-A(14)*RR-B(14)*RR2)
      T15=1.0D0-EXP(-A(15)*RR-B(15)*RR2)
      T16=1.0D0-EXP(-A(16)*RR-B(16)*RR2)

      DDISP=6.0D0*CHH(6)/R6*T6**5*(T6/R+(1.0D0-T6)*(-A(6)-2.0D0 &
           *B(6)*RR)*DRR)+ &
       8.0D0*CHH(8)/R8*T8**7*(T8/R+(1.0D0-T8)*(-A(8)-2.0D0* &
            B(8)*RR)*DRR)+ &
       10.0D0*CHH(10)/R10*T10**9*(T10/R+(1.0D0-T10)*(-A(10)-2.0D0* &
            B(10)*RR)*DRR)+ &
      11.0D0*CHH(11)/R11*T11**10*(T11/R+(1.0D0-T11)*(-A(11)-2.0D0* &
        B(11)*RR)*DRR)+ &
      12.0D0*CHH(12)/R12*T12**11*(T12/R+(1.0D0-T12)*(-A(12)-2.0D0* &
            B(12)*RR)*DRR)+ &
      13.0D0*CHH(13)/R13*T13**12*(T13/R+(1.0D0-T13)*(-A(13)-2.0D0* &
            B(13)*RR)*DRR)+ &
      14.0D0*CHH(14)/R14*T14**13*(T14/R+(1.0D0-T14)*(-A(14)-2.0D0* &
            B(14)*RR)*DRR)+ &
      15.0D0*CHH(15)/R15*T15**14*(T15/R+(1.0D0-T15)*(-A(15)-2.0D0* &
            B(15)*RR)*DRR)+ &
      16.0D0*CHH(16)/R16*T16**15*(T16/R+(1.0D0-T16)*(-A(16)-2.0D0* &
            B(16)*RR)*DRR)
      DDISHHT_f14=DDISP
      RETURN
      END




      BLOCK DATA HHDAT_f14
!     *************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/AN_f14/A(20)
      COMMON/BN_f14/B(20)

      COMMON/DISPC1_f14/CHH(16)


      COMMON/DIATDI_f14/R0HH,RMHHT
!    **************************************************************



      DATA R0HH,RMHHT/6.928203D0,7.82D0/
!    **************************************************************
!     AN E BN OBTIDOS A PARTIR DOS COEF.(A0,A1,B0,B1) DO ARTG. VAR82:857
      DATA A/0.000000000D+00,0.000000000D+00,0.000000000D+00 &
       ,0.500798749D+01,0.384282945D+01,0.309513333D+01 &
       ,0.257767037D+01,0.219990002D+01,0.191291419D+01 &
       ,0.168807142D+01,0.150753106D+01,0.135962478D+01 &
       ,0.123641324D+01,0.113231455D+01,0.104329456D+01 &
       ,0.966368248D+00,0.899281484D+00,0.840301611D+00 &
       ,0.788075808D+00,0.741533075D+00/
      DATA B/0.000000000D+00,0.000000000D+00,0.000000000D+00 &
       ,0.106645006D+02,0.967581549D+01,0.877878947D+01 &
       ,0.796492498D+01,0.722651228D+01,0.655655639D+01 &
       ,0.594871080D+01,0.539721740D+01,0.489685187D+01 &
       ,0.444287427D+01,0.403098404D+01,0.365727935D+01 &
       ,0.331822010D+01,0.301059437D+01,0.273148802D+01 &
       ,0.247825708D+01,0.224850269D+01/
!    **************************************************************
      DATA CHH/5*0.0D0,6.499027D0,0.0d0,124.3991D0,0.0d0,3285.828D0, &
            -3475.0D0,1.215D5,-2.914D5,6.061D6,-2.305D7,3.938D8/


      END





          real*8 function f1_hh(r)
          implicit none
          real*8 r,gamma,rmed
          gamma=5.0d0
          rmed=14.0d0/0.52917721092d0
          f1_hh=0.5d0*(1.0d0+tanh(gamma*(r-rmed)))
          return
          end

          real*8 function df1_hh(r)
          implicit none
          real*8 r,gamma,rmed
          gamma=5.0d0
          rmed=14.0d0/0.52917721092d0
          df1_hh=0.5d0*gamma/cosh(gamma*(r-rmed))**2
          return
          end



