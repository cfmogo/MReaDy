 
      subroutine potOH_pH_f19(rin,Vout,dVout)
		
      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: Vout,dVout
      double precision :: VOHpH,DVOHpH


      Vout=VOHpH(rin/5.2917721092d-11)*4.3597482D-18
      dVout=DVOHpH(rin/5.2917721092d-11)*4.3597482D-18/5.2917721092d-11
      end subroutine potOH_pH_f19

      FUNCTION VOHpH(R)
!     ****************************************************************

!     TO COMPUTE THE HFACE FOR  O...H

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!     ****************************************************************
      

      VOHpH=EHFOHpH(R)+DISOHpH(R)+VELEpH(R)
      
      RETURN
      END

      FUNCTION EHFOHpH(R)
!     ****************************************************************

!     TO COMPUTE THE EHF FOR  O...H

!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3)
!     ****************************************************************
      DATA R0OH,RMOH/6.294894D0,3.6D0/
      data D,ASV/1.221360716D-3,4.77340309D0,2.220897756D0,0.613384505D0/  
!     ****************************************************************
      X=R-RMOH
      R2=X*X
      R3=R2*X
      EHFOHpH=D*DEXP(-ASV(1)*X-ASV(2)*R2-ASV(3)*R3)

      RETURN
      END


      FUNCTION DISOHpH(R)
!     ****************************************************************
!
!     TO COMPUTE THE DISPERSION FOR  O...H
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Dimension COH(10)     
      DATA COH/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,10.00D0,0.0D0,180.447673D0, &
              0.0D0,3685.25842D0/
      DATA R0OH,RMOH/6.294894D0,3.6D0/
!     ****************************************************************
      DISOHpH=DISPpH(R,COH(6),COH(8),COH(10),R0OH,RMOH)
      RETURN
      END

  
      FUNCTION VELEpH(R)
!     ****************************************************************
!
!     TO COMPUTE THE electrostatic interaction FOR  O...H
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA R0OH,RMOH/6.294894D0,3.6D0/
!     ****************************************************************
      qO=-0.44d0/3.0d0
      qH=+0.25273d0/3.0d0
      VELEpH=qO*qH/R*fdampele(R)
      RETURN
      END

      FUNCTION DVELEpH(R)
!     ****************************************************************
!
!     TO COMPUTE THE electrostatic interaction FOR  O...H
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA R0OH,RMOH/6.294894D0,3.6D0/
!     ****************************************************************
      qO=-0.44d0/3.0d0
      qH=+0.25273d0/3.0d0
      DVELEpH=(-qO*qH/R**2*fdampele(R)+qO*qH/R*dfdampele(R))
      RETURN
      END

  


      FUNCTION DVOHpH(R)
!     ****************************************************************
!
!     TO COMPUTE THE DER. OF THE HFACE FOR  O...H
!
!     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3)
      Dimension COH(10)     
      DATA COH/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,10.00D0,0.0D0,180.447673D0, &
              0.0D0,3685.25842D0/
      DATA R0OH,RMOH/6.294894D0,3.6D0/
!     ****************************************************************
      data D,ASV/1.221360716D-3,4.77340309D0,2.220897756D0,0.613384505D0/  
!     ****************************************************************
      X=R-RMOH
      R2=X*X
      R3=R2*X
      POL=D
      DPOL=1.0d0
      POT=DEXP(-ASV(1)*X-ASV(2)*R2-ASV(3)*R3)
      DVOHpH=(-ASV(1)-2.0d0*ASV(2)*X-3.0d0*ASV(3)*R2)*POT*POL+   &
           DDISPpH(R,COH(6),COH(8),COH(10),R0OH,RMOH)+DVELEpH(R)
      RETURN
      END

      FUNCTION DISPpH(R,C6,C8,C10,R0,RM)
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
      DISPpH=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

      FUNCTION DDISPpH(R,C6,C8,C10,R0,RM)
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
      DDISPpH=6.D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-ADAMP(6)-2.0D0*BDAMP(6)  &
           *RR)*DRR)+     &
           8.D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-ADAMP(8)-2.0D0*BDAMP(8)  &
           *RR)*DRR)+     &
           10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-ADAMP(10)-2.0D0*  &
           BDAMP(10)*RR)*DRR) 
      RETURN
      END


Function fdampele(r)
!     **************************************************************
!  To damp the electrostatic interaction due to orbital overlap
!
IMPLICIT none
real*8 fdampele,r
real*8 a(20),b(20),R0OH,RMOH,r0
!
      DATA R0OH,RMOH/6.294894D0,3.6D0/   
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

R0=0.5D0*(RMOH+2.5D0*R0OH)

fdampele=(1-exp(-A(3)*(r/r0)-B(3)*(r/r0)**2))

return
end  


Function Dfdampele(r)
!     **************************************************************
!  To damp the electrostatic interaction due to orbital overlap
!
IMPLICIT none
real*8 Dfdampele,r
real*8 a(20),b(20),R0OH,RMOH,r0

!
      DATA R0OH,RMOH/6.294894D0,3.6D0/  
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
R0=0.5D0*(RMOH+2.5D0*R0OH)

Dfdampele=(A(3)/r0+2.0d0*B(3)*(r/r0**2))*exp(-A(3)*(r/r0)-B(3)*  &
          (r/r0)**2)   

return
end 


