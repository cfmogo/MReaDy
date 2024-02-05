!       program teste
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!      do i=30,100,1
!      r=i/10.0d0
!      print*, r,VArC_f43(r),DVArC_f43(r)
!      end do
!      do i=6330,6360,1
!      r=i/1000.0d0
!      print*, r,VArC_f43(r),DVArC_f43(r)
!      end do
!      do i=7050,7150,1
!      r=i/1000.0d0
!      print*, r,VArC_f43(r),DVArC_f43(r)
!      end do
!      end

      subroutine potArC_f43(rin,Vout,dVout)
      use phys_parameters

      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: Vout,dVout
      double precision :: VArC_f43,DVArC_f43
      double precision :: fs,dfs,f1_lr,df1_lr                                                      

      Vout=VArC_f43(rin/5.2917721092d-11)*4.3597482D-18
      dVout=DVArC_f43(rin/5.2917721092d-11)*4.3597482D-18/ &
       5.2917721092d-11

      if (Rin.gt.(rl_cut_off-0.2d-9)) then              
        fs=f1_lr(Rin)                                                   
        dfs=df1_lr(Rin)                                             
 
        Vout=Vout*(1.0d0-fs)                            
        dVout=dVout*(1.0d0-fs)-Vout*dfs
      end if

      end subroutine  potArC_f43

      FUNCTION VArC_f43(R)
!     ****************************************************************

!     TO COMPUTE THE HFACE FOR  Ar...C   3Pi

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ****************************************************************
      VArC_f43=EHFArC_f43(R)+DISArC_f43(R)
      RETURN
      END

      FUNCTION EHFArC_f43(R)
!     ****************************************************************

!     TO COMPUTE THE EHF FOR  Ar...C   3Pi

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(2)
!     ****************************************************************
      DATA D,ASV/0.7102838D+01,-0.1001057D+01,-0.4574497D-01/
!     ****************************************************************
      X=R
      R2=X*X
      R3=R2*X
      R4=R3*X
      EHFArC_f43=D* Exp(ASV(1)*X+ASV(2)*R2) !+ASV(3)*R3+ASV(4)*R4)
!  
      RETURN
      END


      FUNCTION DISArC_f43(R)
!     ****************************************************************
!
!     TO COMPUTE THE DISPERSION FOR  Ar...C   3Pi
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension CArC(10)     
      DATA CArC/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,80.00D0,0.0D0,1937.5827D0, &
              0.0D0,53112.39667D0/
      DATA R0ArC,RMArC/7.590594D0,7.4D0/
!     ****************************************************************
      DISArC_f43=DISP_f43(R,CArC(6),CArC(8),CArC(10),R0ArC,RMArC)
      RETURN
      END


      FUNCTION DVArC_f43(R)
!     ****************************************************************
!
!     TO COMPUTE THE DER. OF THE HFACE FOR  Ar...C   3Pi
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(2)
      Dimension CArC(10)     
      DATA CArC/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,80.00D0,0.0D0,1937.5827D0, &
              0.0D0,53112.39667D0/
      DATA R0ArC,RMArC/7.590594D0,7.4D0/
!     ****************************************************************
      DATA D,ASV/0.7102838D+01,-0.1001057D+01,-0.4574497D-01/
!     ****************************************************************
      X=R
      R2=X*X
      R3=R2*X
      R4=R3*X
!      EHFArC_f43=D* Exp(ASV(1)*X+ASV(2)*R2+ASV(3)*R3+ASV(4)*R4)
      POL=ASV(1)*X+ASV(2)*R2  !+ASV(3)*R3+ASV(4)*R4
      DPOL=ASV(1)+2.0D0*ASV(2)*X   ! +3.0D0*ASV(3)*R2+4.0D0*ASV(4)*R3 
      POT=D*DEXP(pol)
      DVArC_f43=DPOL*POT+DDISP_f43(R,CArC(6),CArC(8),CArC(10),  &
        R0ArC,RMArC)
      RETURN
      END

      FUNCTION DISP_f43(R,C6,C8,C10,R0,RM)
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension ADAMP(10),BDAMP(10)
      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0951333D0,  &
                0.0D0,2.1999000D0,0.0D0,1.6880714D0/
      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7787895D0,  &
                0.0D0,7.2265123D0,0.0D0,5.9487108D0/

!     ****************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      D6=(1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6
      D8=(1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8
      D10=(1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10
      DISP_f43=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

      FUNCTION DDISP_f43(R,C6,C8,C10,R0,RM)
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension ADAMP(10),BDAMP(10)
      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0951333D0,  &
                0.0D0,2.1999000D0,0.0D0,1.6880714D0/
      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7787895D0,  &
                0.0D0,7.2265123D0,0.0D0,5.9487108D0/

!     ****************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      DRR=RR/R
      T6=1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR)
      T8=1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR)
      T10=1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR)
      DDISP_f43=6.D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-ADAMP(6)-2.0D0*BDAMP(6)  &
           *RR)*DRR)+     &
           8.D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-ADAMP(8)-2.0D0*BDAMP(8)  &
           *RR)*DRR)+     &
           10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-ADAMP(10)-2.0D0*  &
           BDAMP(10)*RR)*DRR) 
      RETURN
      END


