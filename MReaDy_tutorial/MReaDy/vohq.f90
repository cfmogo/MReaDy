      subroutine pothoq_f17(rin,VSout,dVSout)

      use phys_parameters
           
      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: VSout,dVSout
      double precision :: VOH_f17,DVOH_f17,R
      double precision :: fs,dfs,f1_lr,df1_lr                                                      

      R=Rin/5.2917721092d-11
      VSout=VOH_f17(R)*4.3597482D-18
      dVSout=DVOH_f17(R)*4.3597482D-18/5.2917721092d-11

      if (Rin.gt.(rl_cut_off-0.2d-9)) then              
        fs=f1_lr(Rin)                                                   
        dfs=df1_lr(Rin)                                             
 
        VSout=VSout*(1.0d0-fs)                            
        dVSout=dVSout*(1.0d0-fs)-VSout*dfs
      end if                                                      

      end subroutine pothoq_f17


      FUNCTION VOH_f17(R)
!     ****************************************************************

!     TO COMPUTE THE HFACE FOR  O...H for the quartet state

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ****************************************************************
 
      VOH_f17=EHFOH_f17(R)+DISOH_f17(R)



      RETURN
      END

      FUNCTION EHFOH_f17(R)
!     ****************************************************************

!     TO COMPUTE THE EHF FOR  O...H for the quartet state

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(2)
!     ****************************************************************
      DATA R0OH,RMOH/6.294894D0,6.50D0/
      DATA D,ASV/0.1007194693D-03,0.1734748439D+01,0.4714414205D-01/
!     ****************************************************************
      X=R-RMOH
      R2=X*X
      R3=R2*X
      EHFOH_f17=D*DEXP(-ASV(1)*X-ASV(2)*R2)
      RETURN
      END


      FUNCTION DISOH_f17(R)
!     ****************************************************************
!
!     TO COMPUTE THE DISPERSION FOR  O...H for the quartet state
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension COH(10)     
      DATA COH/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,10.00D0,0.0D0,180.447673D0, &
              0.0D0,3685.25842D0/
      DATA R0OH,RMOH/6.294894D0,6.50D0/
!     ****************************************************************
      DISOH_f17=DISP_f17(R,COH(6),COH(8),COH(10),R0OH,RMOH)
      RETURN
      END


      FUNCTION DVOH_f17(R)
!     ****************************************************************
!
!     TO COMPUTE THE DER. OF THE HFACE FOR  O...H for the quartet state
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(2)
      Dimension COH(10)     
      DATA COH/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,10.00D0,0.0D0,180.447673D0, &
              0.0D0,3685.25842D0/
      DATA R0OH,RMOH/6.294894D0,6.50D0/
!     ****************************************************************
      DATA D,ASV/0.1007194693D-03,0.1734748439D+01,0.4714414205D-01/
!     ****************************************************************
       
 
      X=R-RMOH
      R2=X*X
      R3=R2*X
      POL=D
      DPOL=1.0d0
      POT=DEXP(-ASV(1)*X-ASV(2)*R2)
      DVOH_f17=(-ASV(1)-2.0d0*ASV(2)*X)*POT*POL+DDISP_f17(R,COH(6),COH(8),COH(10),  &
        R0OH,RMOH)


      RETURN
      END

      FUNCTION DISP_f17(R,C6,C8,C10,R0,RM)
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
      DISP_f17=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

      FUNCTION DDISP_f17(R,C6,C8,C10,R0,RM)
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
      DDISP_f17=6.D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-ADAMP(6)-2.0D0*BDAMP(6)  &
           *RR)*DRR)+     &
           8.D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-ADAMP(8)-2.0D0*BDAMP(8)  &
           *RR)*DRR)+     &
           10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-ADAMP(10)-2.0D0*  &
           BDAMP(10)*RR)*DRR) 
      RETURN
      END

          real*8 function f1_oh(r)                                   
          implicit none                                          
          real*8 r,gamma,rmed                                       
          gamma=5.0d0                                            
          rmed=14.0d0/0.52917721092d0                               
          f1_oh=0.5d0*(1.0d0+tanh(gamma*(r-rmed)))          
          return                                            
          end                                               
                                                            
          real*8 function df1_oh(r)                         
          implicit none                                     
          real*8 r,gamma,rmed                               
          gamma=5.0d0                                       
          rmed=14.0d0/0.52917721092d0                       
          df1_oh=0.5d0*gamma/cosh(gamma*(r-rmed))**2        
          return                                            
          end                                               
                                                            
                                                            
 
