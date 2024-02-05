!      program teste
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      do i=30,100,1
!        r=i/10.0d0
!        print*, r,VAr2_f40(r),DVAr2_f40(r)
!      end do
!      do i=6330,6360,1
!        r=i/1000.0d0
!        print*, r,VAr2_f40(r),DVAr2_f40(r)
!      end do
!      do i=7050,7150,1
!        r=i/1000.0d0
!        print*, r,VAr2_f40(r),DVAr2_f40(r)
!      end do
!      end


      subroutine potAr2_f40(rin,Vout,dVout)
      use phys_parameters

      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: Vout,dVout
      double precision :: VAr2_f40,DVAr2_f40
      double precision :: fs,dfs,f1_lr,df1_lr                                                      

      Vout=VAr2_f40(rin/5.2917721092d-11)*4.3597482D-18
      dVout=DVAr2_f40(rin/5.2917721092d-11)*4.3597482D-18/ &
       5.2917721092d-11

      if (Rin.gt.(rl_cut_off-0.2d-9)) then              
        fs=f1_lr(Rin)                                                   
        dfs=df1_lr(Rin)                                             
 
        Vout=Vout*(1.0d0-fs)                            
        dVout=dVout*(1.0d0-fs)-Vout*dfs
      end if

      end subroutine  potAr2_f40

      FUNCTION VAr2_f40(R)
!     ****************************************************************

!     TO COMPUTE THE HFACE FOR  Ar...Ar

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ****************************************************************
      VAr2_f40=EHFAr2_f40(R)+DISAr2_f40(R)
      RETURN
      END

      FUNCTION EHFAr2_f40(R)
!     ****************************************************************

!     TO COMPUTE THE EHF FOR  Ar...Ar

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(2)
!     ****************************************************************
      DATA D,ASV/41.1228D0,0.8005D0,0.082870D0/
!     ****************************************************************
      X=R
      X2=X*X
      EHFAr2_f40=D/X*DEXP(-ASV(1)*X-ASV(2)*X2)
      RETURN
      END


      FUNCTION DISAr2_f40(R)
!     ****************************************************************
!
!     TO COMPUTE THE DISPERSION FOR  Ar...Ar
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension CAr2(10)     
      DATA CAr2/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,66.0791D0,0.0D0,1506.144D0, &
              0.0D0,53733.63D0/
      DATA RMAr2,R0Ar2/7.10919D0,7.26845D0/
!     ****************************************************************
      DISAr2_f40=DISP_f40(R,CAr2(6),CAr2(8),CAr2(10),R0Ar2,RMAr2)
      RETURN
      END


      FUNCTION DVAr2_f40(R)
!     ****************************************************************
!
!     TO COMPUTE THE DER. OF THE HFACE FOR  Ar...Ar
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(2)
      Dimension CAr2(10)     
      DATA CAr2/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,66.0791D0,0.0D0,1506.144D0, &
              0.0D0,53733.63D0/
      DATA RMAr2,R0Ar2/7.10919D0,7.26845D0/
!     ****************************************************************
      DATA D,ASV/41.1228D0,0.8005D0,0.082870D0/
!     ****************************************************************
      X=R
      X2=X*X
      POT=DEXP(-ASV(1)*X-ASV(2)*X2)
!      EHFAr2=D/X*DEXP(-ASV(1)*X-ASV(2)*X2)
!      DEHFAr2=-D/X**2*POT+D/X*(-ASV(1)-2.0d0*ASV(2)*X)*POT
      DEHFAr2=D/X*POT*(-1.0d0/X+(-ASV(1)-2.0d0*ASV(2)*X))

      DVAr2_f40=DEHFAr2+DDISP_f40(R,CAr2(6),CAr2(8),CAr2(10),  &
        R0Ar2,RMAr2)
      RETURN
      END

      FUNCTION DISP_f40(R,C6,C8,C10,R0,RM)
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension ADAMP(10),BDAMP(10)
      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0951333D0,  &
                0.0D0,2.1999000D0,0.0D0,1.6880714D0/
      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7787895D0,  &
                0.0D0,7.2265123D0,0.0D0,5.9487108D0/
!      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0950974D0,  &
!                0.0D0,2.1998700D0,0.0D0,1.6880458D0/
!      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7788422D0,  &
!                0.0D0,7.2265701D0,0.0D0,5.9487703D0/

!     ****************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      D6=(1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6
      D8=(1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8
      D10=(1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10
      DISP_f40=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

      FUNCTION DDISP_f40(R,C6,C8,C10,R0,RM)
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension ADAMP(10),BDAMP(10)
      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0951333D0,&
                0.0D0,2.1999000D0,0.0D0,1.6880714D0/
      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7787895D0,&
                0.0D0,7.2265123D0,0.0D0,5.9487108D0/
!      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0950974D0,  &
!                0.0D0,2.1998700D0,0.0D0,1.6880458D0/
!      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7788422D0,  &
!                0.0D0,7.2265701D0,0.0D0,5.9487703D0/

!     ****************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      DRR=RR/R
      T6=1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR)
      T8=1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR)
      T10=1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR)
      DDISP_f40=6.D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)                  &
      *(-ADAMP(6)-2.0D0*BDAMP(6)*RR)*DRR)+                         &
      8.D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-ADAMP(8)-2.0D0*BDAMP(8)  &
      *RR)*DRR)+10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)           &
      *(-ADAMP(10)-2.0D0*BDAMP(10)*RR)*DRR) 
      RETURN
      END


