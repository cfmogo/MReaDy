!      program teste
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      do i=30,100,1
!      r=i/10.0d0
!      print*, r,VHAr_f42(r),DVHAr_f42(r)
!      end do
!      do i=6330,6360,1
!      r=i/1000.0d0
!      print*, r,VHAr_f42(r),DVHAr_f42(r)
!      end do
!      do i=7050,7150,1
!      r=i/1000.0d0
!      print*, r,VHAr_f42(r),DVHAr_f42(r)
!      end do
!      end


      subroutine potArH_f42(rin,Vout,dVout)
      use phys_parameters

      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: Vout,dVout
      double precision :: VHAr_f42,DVHAr_f42
      double precision :: fs,dfs,f1_lr,df1_lr                                                      

      Vout=VHAr_f42(rin/5.2917721092d-11)*4.3597482D-18
      dVout=DVHAr_f42(rin/5.2917721092d-11)*4.3597482D-18/ &
       5.2917721092d-11

      if (Rin.gt.(rl_cut_off-0.2d-9)) then              
        fs=f1_lr(Rin)                                                   
        dfs=df1_lr(Rin)                                             
 
        Vout=Vout*(1.0d0-fs)                            
        dVout=dVout*(1.0d0-fs)-Vout*dfs
      end if
      end subroutine  potArH_f42

      function VHAr_f42(R)
!     ****************************************************************
!
!     TO COMPUTE THE HFACE FOR  H...Ar
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ****************************************************************
      VHAr_f42=EHFHAr_f42(R)+DISHAr_f42(R)
      RETURN
      END


      function EHFHAr_f42(R)
!     ****************************************************************
!
!     TO COMPUTE THE EHF FOR  H...Ar
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     a AND b HAVE BEEN CALCULATED IN ORDER TO GET THE VAN DER WAALS
!     MINIMUM AS Vmin=-1.5273E-4 Eh AND Rmin=6.84 a0. DAMP FUNCTIONS
!     FROM VARANDAS (1987) HAVE BEEN USED. THE DISPERSION COEFFICIENTS
!     HAVE BEEN OBTAINED FROM STANDARD & CERTAIN (1985)
      a=61.93993d0         !61.94d0
      b=1.907065d0       !1.907d0
      EHFHAr_f42=a*exp(-b*R)
      RETURN
      END


      function DISHAr_f42(R)
!     ****************************************************************
!
!     TO COMPUTE THE DISPERSION FOR  H...Ar
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPAr_f42/COAr(10),CHAr(10)
      COMMON/DIATAr_f42/R0OAr,R0HAr
!     ****************************************************************
      DISHAr_f42=DISPN_f42(R,CHAr(6),CHAr(8),CHAr(10),R0HAr)
      RETURN
      END


      FUNCTION DISPN_f42(R,C6,C8,C10,R0)
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension ADAMP(10),BDAMP(10)
      common/damps_f42/d1(3),d2(3)
!     ****************************************************************
      do i=1,3
      adamp(2*i+4)=d1(i)
      bdamp(2*i+4)=d2(i)
      enddo
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      rho=5.5+1.25*r0
      rr=r/rho
      D6=(1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6
      D8=(1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8
      D10=(1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10
      DISPN_f42=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

!  *********************************************************************
!  DERIVATIVES OF THE POTENTIAL ENERGY FUNCTION FOR ArHO2
!  *********************************************************************


      function DVHAr_f42(R)
!     ****************************************************************
!
!     TO COMPUTE THE DER. OF THE HFACE FOR  H...Ar
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPAr_f42/COAr(10),CHAr(10)
      COMMON/DIATAr_f42/R0OAr,R0HAr
!     ****************************************************************
!     a AND b HAVE BEEN CALCULATED IN ORDER TO GET THE VAN DER WAALS
!     MINIMUM AS Vmin=-1.5273E-4 Eh AND Rmin=6.84 a0. DAMP FUNCTIONS
!     FROM VARANDAS (1987) HAVE BEEN USED. THE DISPERSION COEFFICIENTS
!     HAVE BEEN OBTAINED FROM STANDARD & CERTAIN (1985)
      a=61.93993d0         !61.94d0
      b=1.907065d0       !1.907d0
      DVHAr_f42=-a*b*EXP(-b*R)+ &
      DDISPN_f42(R,CHAr(6),CHAr(8),CHAr(10),R0HAr)
      RETURN
      END


      FUNCTION DDISPN_f42(R,C6,C8,C10,R0)
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension ADAMP(10),BDAMP(10)
      common/damps_f42/d1(3),d2(3)
!     ****************************************************************
      do i=1,3
      adamp(2*i+4)=d1(i)
      bdamp(2*i+4)=d2(i)
      enddo
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      rho=5.5+1.25*r0
      rr=r/rho
      DRR=RR/R
      T6=1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR)
      T8=1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR)
      T10=1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR)
      DDISPN_f42=6.D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)  &
      *(-ADAMP(6)-2.0D0*BDAMP(6)*RR)*DRR)+      &
       8.D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)        &
       *(-ADAMP(8)-2.0D0*BDAMP(8)*RR)*DRR)+     &
       10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10) &
      *(-ADAMP(10)-2.0D0*                       &
      BDAMP(10)*RR)*DRR)
      RETURN
      END


      BLOCK DATA F32

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPAr_f42/COAr(10),CHAr(10)
      COMMON/DIATAr_f42/R0OAr,R0HAr
      common/damps_f42/d1(3),d2(3)
!     ***************************************************************

      DATA COAr/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,32.17D0,0.0D0,589.99661D0,&
      0.0D0,13176.4259D0/
!   valores obtidos do trabalho de Standard e Certain (1985)!
      DATA CHAr/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,20.00D0,0.0D0,426.00D0,&
      0.0D0,12300.0D0/
!   valores obtidos pela regra de combinaçao!

      DATA R0OAr,R0HAr/6.4651D0,7.09835D0/
      DATA (D1(I),I=1,3)/.46547873E1,.38038890E1,.32525512E1/
      DATA (D2(I),I=1,3)/.96802180E1,.79933059E1,.66003615E1/
      END

