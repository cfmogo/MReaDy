      subroutine poto2q_f18(rin,VSout,dVSout)
		
      use phys_parameters
      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: VSout,dVSout
      double precision :: VO2Q_f18,DVO2Q_f18,R
      double precision :: fs,dfs,f1_lr,df1_lr                                                      
 

      R=Rin/5.2917721092d-11
      VSout=VO2Q_f18(R)*4.3597482D-18
      dVSout=DVO2Q_f18(R)*4.3597482D-18/5.2917721092d-11


      if (Rin.gt.(rl_cut_off-0.2d-9)) then              
        fs=f1_lr(Rin)                                                   
        dfs=df1_lr(Rin)                                             
 
        VSout=VSout*(1.0d0-fs)                            
        dVSout=dVSout*(1.0d0-fs)-VSout*dfs
      end if                                                      
 

      end subroutine poto2q_f18


      FUNCTION VO2Q_f18(R)
!     ****************************************************************

!     TO COMPUTE THE HFACE FOR  O...O Quintet

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ****************************************************************
      VO2Q_f18=EHFOOQ_f18(R)+DISOOQ_f18(R)


      RETURN
      END

      FUNCTION EHFOOQ_f18(R)
!     ****************************************************************

!     TO COMPUTE THE EHF FOR  O...O

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(2)
!     ****************************************************************
      DATA R0OO,RMOO/5.661693D0,6.85D0/
      DATA D,ASV/9.2794667846D-05,1.8202063065D+00,1.2816216351D-02/
!     ****************************************************************
      X=R-RMOO
      R2=X*X
      R3=R2*X
      EHFOOQ_f18=D*DEXP(-ASV(1)*X-ASV(2)*R2)
      RETURN
      END


      FUNCTION DISOOQ_f18(R)
!     ****************************************************************
!
!     TO COMPUTE THE DISP_f18ERSION FOR  O...O  Quintet
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension COO(10)     
      DATA R0OO,RMOO/5.661693D0,6.85D0/
      DATA COO/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,15.40D0,0.0D0,235.219943D0, &
               0.0D0,4066.23929D0/
!     ****************************************************************
      DISOOQ_f18=DISP_f18(R,COO(6),COO(8),COO(10),R0OO,RMOO)
      RETURN
      END


      FUNCTION DVO2Q_f18(R)
!     ****************************************************************
!
!     TO COMPUTE THE DER. OF THE HFACE FOR  O...O  Quintet
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(2)
      Dimension COO(10)     
      DATA R0OO,RMOO/5.661693D0,6.85D0/
      DATA COO/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,15.40D0,0.0D0,235.219943D0, &
               0.0D0,4066.23929D0/
!     ****************************************************************
      DATA D,ASV/9.2794667846D-05,1.8202063065D+00,1.2816216351D-02/
!     ****************************************************************
      X=R-RMOO
      R2=X*X
      R3=R2*X
      POL=D
      DPOL=1.0d0
      POT=DEXP(-ASV(1)*X-ASV(2)*R2)
      DVO2Q_f18=(-ASV(1)-2.0d0*ASV(2)*X)*POT*POL+  &
      DDISP_f18(R,COO(6),COO(8),COO(10),R0OO,RMOO)

      RETURN
      END

      FUNCTION DISP_f18(R,C6,C8,C10,R0,RM)
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
      DISP_f18=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

      FUNCTION DDISP_f18(R,C6,C8,C10,R0,RM)
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
      DDISP_f18=6.D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*  &
          (-ADAMP(6)-2.0D0*BDAMP(6) *RR)*DRR)+     &
           8.D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-ADAMP(8)-2.0D0*BDAMP(8)  &
           *RR)*DRR)+     &
           10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-ADAMP(10)-2.0D0*  &
           BDAMP(10)*RR)*DRR) 
      RETURN
      END

!          real*8 function f1_oo(r)                                   
!          implicit none                                          
!          real*8 r,gamma,rmed                                       
!          gamma=5.0d0                                            
!          rmed=14.0d0/0.52917721092d0                               
!          f1_oo=0.5d0*(1.0d0+tanh(gamma*(r-rmed)))          
!          return                                            
!          end                                               
!                                                            
!          real*8 function df1_oo(r)                         
!          implicit none                                     
!          real*8 r,gamma,rmed                               
!          gamma=5.0d0                                       
!          rmed=14.0d0/0.52917721092d0                       
!          df1_oo=0.5d0*gamma/cosh(gamma*(r-rmed))**2        
!          return                                            
!          end                                               
                                                            
                                                            
 
