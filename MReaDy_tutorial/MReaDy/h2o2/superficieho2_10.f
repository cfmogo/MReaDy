

C     ******************************************************************
C     TO COMPUTE THE HO2 SURFACE in atomic units
C     ******************************************************************
      FUNCTION VHO23c_10(R1J,R2J,R3J)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/THRBOD_10/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3
!$OMP THREADPRIVATE(/THRBOD_10/)



      R1=R1J
      R2=R2J
      R3=R3J
      Q1=1.0D0/DSQRT(3.0D0)*(R1+R2+R3)
      Q2=1.0D0/DSQRT(2.0D0)*(R2-R3)
      Q3=1.0D0/DSQRT(6.0D0)*(2.0D0*R1-R2-R3)

      gammar=3.0d0
      corr2=VOH_10(r2)-VOHP_10(r2)
      corr3=VOH_10(r3)-VOHP_10(r3)
      damp1=0.5d0*(1.0d0-tanh(gammar*(r1-5.0d0)))
      damp2=0.5d0*(1.0d0-tanh(gammar*(r2-5.0d0)))
      damp3=0.5d0*(1.0d0-tanh(gammar*(r3-5.0d0)))

      V=THREBQ_10(Q1,Q2,Q3)+
     1  EXDIS_10(R1,R2,R3)+ELECT_10(R1,R2,R3)+VSPECT_10(R1,R2,R3)
     2  +corr2*damp1*damp3+corr3*damp1*damp2

c      V=THREBQ_10(Q1,Q2,Q3)
     
      VHO23c_10=V
      END


      FUNCTION VSPECT_10(r1,r2,r3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Y)

      R1eq=2.5143d0
      R2eq=1.8346d0
      R3eq=3.4592d0

      r1d=r1-r1eq
      s1d=((r2+r3)-(r2eq+r3eq))/sqrt(2.0d0)
      s2d2=((r2-r3)**2-(r2eq-r3eq)**2)/2.0d0

      c1fix=2.0d0
      c2fix=2.0d0
      c3fix=2.0d0
      c4fix=2.0d0
     
      C1=1.0944219d-03 
      C2=2.9389603d-04
      C3=-9.4723705d-03
      C4=7.1459133d-02
      C5=2.2984321d-02
      C6=-4.8325504d-02
      C7=4.7234450d-04
      C8=-5.9777409d-02
      C9=1.8128239d-02
      C10=6.4184643d-02
      C11=-1.1384161d-02
      C12=1.2619628d-03
      C13=-2.9029563d-03

      t=c1+c2*r1d+c3*s1d+c4*r1d**2+c5*s1d**2+c6*r1d*s1d+
     1  c7*s2d2+c8*r1d**3+c9*r1d*s1d**2+c10*r1d**2*s1d+
     2  c11*r1d*s2d2+c12*s1d**3+c13*s1d*s2d2
      dec=c1fix*r1d**2+c2fix*s1d**2+c3fix*r1d*s1d+c4fix*s2d2**2

      VSPECT_10=t*exp(-dec)
      
      END

      

C     ****************************************************************
C
C     TO COMPUTE THE THREE BODY TERM IN SIMETRIC COORDINATES Q1,Q2,Q3
C
C     ****************************************************************
      FUNCTION THREBQ_10(Q1,Q2,Q3) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/COEFFHO2_10/C(52)
      COMMON/THRBOD_10/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3
      COMMON/REFGEOHO2_10/R10,R20,R30
!$OMP THREADPRIVATE(/THRBOD_10/)
     

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
      POLQ=C(1)*Q1+C(2)*Q12+C(3)*TQ1+C(4)*Q13+C(5)*Q1*TQ1+
     1C(6)*Q3*TQ2+C(7)*Q14+C(8)*Q12*TQ1+C(9)*TQ1**2+C(10)*Q1*Q3*TQ2+
     2C(11)*Q3+C(12)*Q1*Q3+C(13)*TQ3+C(14)*Q12*Q3+C(15)*Q1*TQ3+
     3C(16)*Q3*TQ1+C(17)*Q13*Q3+C(18)*Q12*TQ3+C(19)*Q1*Q3*TQ1+
     4C(20)*Q32*TQ2+C(21)*TQ1*TQ3+C(22)+C(23)*Q15+C(24)*Q13*TQ1+
     5C(25)*Q1*TQ12+C(26)*Q12*Q3*TQ2+C(27)*Q3*TQ1*TQ2+C(28)*Q14*Q3+
     6C(29)*Q13*TQ3+C(30)*Q12*Q3*TQ1+C(31)*Q1*Q32*TQ2+C(32)*Q1*TQ1*TQ3+
     7C(33)*Q3*TQ12+C(34)*Q3*TQ2*TQ3+C(35)*Q16+C(36)*Q14*TQ1+
     8C(37)*Q12*TQ12+C(38)*Q13*Q3*TQ2+C(39)*Q1*Q3*TQ1*TQ2+C(40)*TQ13+
     9C(41)*Q32*TQ22+C(42)*Q15*Q3+C(43)*Q14*TQ3+C(44)*Q13*Q3*TQ1+
     AC(45)*Q12*Q32*TQ2+C(46)*Q12*TQ1*TQ3+C(47)*Q1*Q3*TQ12+
     BC(48)*Q1*Q3*TQ2*TQ3+C(49)*Q32*TQ1*TQ2+C(50)*TQ12*TQ3
      DECAY1=1.0D0-DTANH(C(51)*S1)
      DECAY2=1.0D0-DTANH(C(52)*S2)
      DECAY3=1.0D0-DTANH(C(52)*S3)
      THREBQ_10=POLQ*DECAY1*DECAY2*DECAY3

      END


C     ****************************************************************
C
C     TO COMPUTE THE HFACE FOR  O...H
C
C     ****************************************************************
      FUNCTION VOH_10(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      VOH_10=EHFOH_10(R)+DISOH_10(R)
      
      END


C     ****************************************************************
C
C     TO COMPUTE THE EHF FOR  O...H
C
C     ****************************************************************
      FUNCTION EHFOH_10(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      SAVE /DIATDIHO2_10/,D,ASV
      DATA D,ASV/0.13825385D0,2.6564788D0,1.7450528D0,0.71014391D0,
     1           2.5453276D0/

      X=R-RMOH
      R2=X*X
      R3=R2*X
      EHFOH_10=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*DEXP(-ASV(4)*X)
      
      END


C     ****************************************************************
C
C     TO COMPUTE THE DISPERSION FOR  O...H
C
C     ****************************************************************
      FUNCTION DISOH_10(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_10/COO(10),COH(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      SAVE /DISPC_10/,/DIATDIHO2_10/
      
      DISOH_10=DISP_10(R,COH(6),COH(8),COH(10),R0OH,RMOH)

      END



C     ****************************************************************
C
C     TO COMPUTE THE HFACE FOR  O...O
C
C     ****************************************************************
      FUNCTION VOO_10(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      VOO_10=EHFOO_10(R)+DISOO_10(R)

      END


      
C     ****************************************************************
C
C     TO COMPUTE THE EHF FOR  O...O
C
C     ****************************************************************
      FUNCTION EHFOO_10(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      SAVE /DIATDIHO2_10/,D,ASV      
      DATA D,ASV/0.14291202D0,3.6445906D0,3.9281238D0,2.0986689D0,
     1           3.3522498D0/

      X=R-RMOO
      R2=X*X
      R3=R2*X
      EHFOO_10=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*DEXP(-ASV(4)*X)

      END


      
C     ****************************************************************
C
C     TO COMPUTE THE DISPERSION FOR  O...O
C
C     ****************************************************************
      FUNCTION DISOO_10(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_10/COO(10),COH(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      SAVE /DISPC_10/,/DIATDIHO2_10/
      
      DISOO_10=DISP_10(R,COO(6),COO(8),COO(10),R0OO,RMOO)

      END

      
C     ****************************************************************
C     TO COMPUTE THE EXCHANGE - DISPERSION TERM
C     ****************************************************************
      FUNCTION EXDIS_10 (R1,R2,R3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_10/COO(10),COH(10)
      COMMON/RKVAL_10/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DISCO_10/CEFOO(10),CEFOH2(10),CEFOH3(10),CEDOO(10),
     1 CEDOH2(10)
      COMMON/DISCO2_10/CEDOH3(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
!$OMP THREADPRIVATE(/DISCO_10/,/DISCO2_10/)


      SAVE /DISPC_10/,/DIATDIHO2_10/,/RKVAL_10/
      
      DO 10 IN=6,10,2
        CEFOO(IN)=CEF_10(COO(IN),RK0OH(IN),RK1OH(IN),RK0OH(IN)
     1          ,RK1OH(IN),RMOH,RMOH,R2,R3)
        CEDOO(IN)=CEFOO(IN)-COO(IN)
        CEFOH2(IN)=CEF_10(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),
     1        RK1OH(IN),RMOO,RMOH,R1,R3)
        CEDOH2(IN)=CEFOH2(IN)-COH(IN)
        CEFOH3(IN)=CEF_10(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),
     1           RK1OH(IN),RMOO,RMOH,R1,R2)
        CEDOH3(IN)=CEFOH3(IN)-COH(IN)
   10 CONTINUE
      EXDIS_10=DISP_10(R1,CEDOO(6),CEDOO(8),CEDOO(10),R0OO,RMOO)
     1     +DISP_10(R2,CEDOH2(6),CEDOH2(8),CEDOH2(10),R0OH,RMOH)
     2     +DISP_10(R3,CEDOH3(6),CEDOH3(8),CEDOH3(10),R0OH,RMOH)

      END


      
C     ****************************************************************
C     TO COMPUTE THE EFECTIVE Cn
C     ****************************************************************
      FUNCTION CEF_10(CAS,RK01,RK11,RK02,RK12,RE1,RE2,R1,R2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      CEF_10=0.5D0*CAS*((1.0D0-RK01*DEXP(-RK11*
     1  (R1-RE1)))*DTANH(RK12*R2)+
     1  (1.0D0-RK02*DEXP(-RK12*(R2-RE2)))*DTANH(RK11*R1))

      END


C     ****************************************************************
C     TO COMPUTE THE ELECTROSTATIC TERM
C     ****************************************************************
      FUNCTION ELECT_10(R1,R2,R3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/POLAR_10/C4,C5
      COMMON/RKVAL_10/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      COMMON/DAMPC_10/ADAMP(10),BDAMP(10)
!      COMMON/WELECT_10/C4OHR2,C5OHR2,C4OHR3,C5OHR3,C4OO,C5OO,TERM4,
!     1  TERM5
!!$OMP THREADPRIVATE(/WELECT_10/)

      SAVE /DIATDIHO2_10/,/RKVAL_10/,/POLAR_10/,/DAMPC_10/
      
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
      TERM4=C4OO/R14*(1.0D0-DEXP(-ADAMP(4)*RROO-BDAMP(4)*RROO**2))**4+
     1 C4OHR2/R24*(1.0D0-DEXP(-ADAMP(4)*RROH2-BDAMP(4)*RROH2**2))**4+
     2 C4OHR3/R34*(1.0D0-DEXP(-ADAMP(4)*RROH3-BDAMP(4)*RROH3**2))**4
      TERM5=C5OO/R15*(1.0D0-DEXP(-ADAMP(5)*RROO-BDAMP(5)*RROO**2))**5+
     1 C5OHR2/R25*(1.0D0-DEXP(-ADAMP(5)*RROH2-BDAMP(5)*RROH2**2))**5+
     2 C5OHR3/R35*(1.0D0-DEXP(-ADAMP(5)*RROH3-BDAMP(5)*RROH3**2))**5
      ELECT_10=TERM4+TERM5

      END

      
C     ***************************************************************      
      FUNCTION DISP_10(R,C6,C8,C10,R0,RM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DAMPC_10/ADAMP(10),BDAMP(10)

      SAVE /DAMPC_10/
      
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      D6=(1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6
      D8=(1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8
      D10=(1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10
      DISP_10=-C6/R6*D6-C8/R8*D8-C10/R10*D10

      END
      


      
C     ***************************************************************
C     DATA FOR HO2 SURFACE
C     ***************************************************************
      BLOCK DATA HO2DAT_10
      
C     ESTE BLOCK DATA CONTIENE LOS DATOS CORRESPONDIENTES A
C     LA SUPERFICIE CALCULADA CON BOO=BOH=1.55
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/COEFFHO2_10/C(52)
      COMMON/DISPC_10/COO(10),COH(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      COMMON/RKVAL_10/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/POLAR_10/C4,C5
      COMMON/DAMPC_10/ADAMP(10),BDAMP(10)
      COMMON/REFGEOHO2_10/R10,R20,R30
      SAVE /COEFFHO2_10/,/DISPC_10/,/DIATDIHO2_10/,/RKVAL_10/,
     1 /POLAR_10/,/DAMPC_10/,/REFGEOHO2_10/
      
      DATA C/
     1  .49040645D+01, -.86748216D+01,  .50555792D+01,  .42941301D+01,
     1 -.41874792D+01,  .13461379D+00, -.99064922D+00,  .13358488D+01,
     1  .13495231D+01, -.18529696D+00, -.23534213D+02,  .24289930D+02,
     1 -.50209026D+01, -.10365484D+02,  .46692224D+01, -.14747138D+01,
     1  .23119718D+01, -.18247842D+01, -.28472166D+00,  .51036509D+00,
     1  .19124083D+00,  .45405729D+01,  .11087611D+00, -.19990481D+00,
     1 -.37356178D+00,  .46142042D-01, -.20565580D+00, -.27015963D+00,
     1  .34085281D+00,  .28321162D+00, -.11558481D+00, -.29448886D+00,
     1 -.52932488D+00,  .58159523D-01, -.48649560D-02,  .11949167D-01,
     1  .21409804D-01, -.20620608D-02,  .30177088D-01,  .27880291D-01,
     1  .88458711D-02,  .13137410D-01, -.24705619D-01, -.31085889D-01,
     1  .34317857D-02,  .52593878D-01,  .79500714D-01, -.79782216D-02,
     2  .31164575D-01, -.28737598D-01,  .98201698D+00,  .62000000D+00/

      DATA R0OO,RMOO,R0OH,RMOH/5.661693D0,2.2818D0,6.294894D0,1.8344D0/
      DATA COO/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,15.40D0,0.0D0,235.219943D0,
     1         0.0D0,4066.23929D0/
      DATA COH/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,10.00D0,0.0D0,180.447673D0,
     1         0.0D0,3685.25842D0/
      DATA C4,C5/-0.92921D0,-1.79000D0/
      DATA RK0OO/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,-.27847758D0,0.0D0,
     1           -.46815641D0,0.0D0,-1.20506384D0/
      DATA RK1OO/0.0D0,0.0D0,0.0D0,3.35224980D0,3.35224980D0,
     1           0.95273753D0,0.0D0,0.94148408D0,0.0D0,0.72379129D0/
      DATA RK0OH/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.02465005D0,0.0D0,
     1           0.05036950D0,0.0D0,0.06294371D0/
      DATA RK1OH/0.0D0,0.0D0,0.0D0,2.54532760D0,2.54532760D0,
     1           0.68758097D0,0.0D0,0.82542359D0,0.0D0,0.94034225D0/
      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0951333D0,
     1           0.0D0,2.1999000D0,0.0D0,1.6880714D0/
      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7787895D0,
     1           0.0D0,7.2265123D0,0.0D0,5.9487108D0/
      DATA R10,R20,R30/2.5143000D0,2.6469057D0,2.6469057D0/
      END




c    ********************************************************************
c     Here starts the derivative part of the program
c     It needs the function THREBQ_10(Q1,Q2,Q3) and the BLOCK DATA HO2DAT_10
c    ********************************************************************




      SUBROUTINE DERVHO23c_10(r1j,r2j,r3j,dg1,dg2,dg3)
c    **************************************************************
C     TO COMPUTE THE DERIVATIVE H2O SURFACE in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(3),G(3)
      COMMON/COEFFHO2_10/C(52)
      COMMON/THRBOD_10/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3
      COMMON/REFGEOHO2_10/R10,R20,R30
!$OMP THREADPRIVATE(/THRBOD_10/)

C    **************************************************************

      R1=R1J
      R2=R2J
      R3=R3J
      Q1=1.0D0/DSQRT(3.0D0)*(R1+R2+R3)
      Q2=1.0D0/DSQRT(2.0D0)*(R2-R3)
      Q3=1.0D0/DSQRT(6.0D0)*(2.0D0*R1-R2-R3)
      X(1)=r1
      X(2)=r2
      X(3)=r3
      term3Q=THREBQ_10(Q1,Q2,Q3)
      S1=X(1)-R10
      S2=X(2)-R20
      S3=X(3)-R30
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
      DQ1R1=1.0D0/SQRT(3.0D0)
      DQ1R2=DQ1R1
      DQ1R3=DQ1R1
      DQ2R1=0.0D0
      DQ2R2=1.0D0/SQRT(2.0D0)
      DQ2R3=-DQ2R2
      DQ3R1=2.0D0/SQRT(6.0D0)
      DQ3R2=-0.5D0*DQ3R1
      DQ3R3=DQ3R2
      DPOQ1=C(1)+2.0D0*C(2)*Q1+3.0D0*C(4)*Q12+C(5)*TQ1+4.0D0*C(7)*Q13+
     1   2.0D0*C(8)*Q1*TQ1+C(10)*Q3*TQ2+C(12)*Q3+2.0D0*C(14)*Q1*Q3+
     2   C(15)*TQ3+3.0D0*C(17)*Q12*Q3+2.0D0*C(18)*Q1*TQ3+C(19)*Q3*TQ1+
     3   5.0D0*C(23)*Q14+3.0D0*C(24)*Q12*TQ1+C(25)*TQ12+2.0D0*C(26)*
     4   Q1*Q3*TQ2+4.0D0*C(28)*Q13*Q3+3.0D0*C(29)*Q12*TQ3+2.0D0*C(30)*
     5   Q1*Q3*TQ1+C(31)*Q32*TQ2+C(32)*TQ1*TQ3+6.0D0*C(35)*Q15+4.0D0*
     6   C(36)*Q13*TQ1+2.0D0*C(37)*Q1*TQ12+3.0D0*C(38)*Q12*Q3*TQ2+
     7   C(39)*Q3*TQ1*TQ2+5.0D0*C(42)*Q14*Q3+4.0D0*C(43)*Q13*TQ3+
     8   3.0D0*C(44)*Q12*Q3*TQ1+2.0D0*C(45)*Q1*Q32*TQ2+
     9   2.0D0*C(46)*Q1*TQ1*TQ3+C(47)*Q3*TQ12+C(48)*Q3*TQ2*TQ3
      DPOQ2=C(3)*2.0D0*Q2+2.0D0*C(5)*Q1*Q2+C(6)*Q3*(-6.0D0*Q2)+
     1   C(8)*Q12*2.0D0*Q2+C(9)*2.0D0*TQ1*2.0D0*Q2+C(10)*Q1*Q3*
     2   (-6.0D0*Q2)+C(13)*2.0D0*Q2+C(15)*Q1*2.0D0*Q2+C(16)*Q3*2.0D0*Q2+
     3   C(18)*Q12*2.0D0*Q2+C(19)*Q1*Q3*2.0D0*Q2+C(20)*Q32*(-6.0D0*Q2)+
     4   C(21)*2.0D0*Q2*TQ3+C(21)*TQ1*2.0D0*Q2+C(24)*Q13*2.0D0*Q2+
     5   C(25)*Q1*2.0D0*TQ1*2.0D0*Q2+C(26)*Q12*Q3*(-6.0D0*Q2)+
     6   C(27)*Q3*2.0D0*Q2*TQ2+C(27)*Q3*TQ1*(-6.0D0*Q2)+
     7   C(29)*Q13*2.0D0*Q2+C(30)*Q12*Q3*2.0D0*Q2+C(31)*Q1*Q32*(-6.0D0
     8   *Q2)+C(32)*Q1*2.0D0*Q2*(TQ1+TQ3)+C(33)*Q3*4.0D0*TQ1*Q2+C(34)*
     9   Q3*Q2*(2.0D0*TQ2-6.0D0*TQ3)+2.0D0*C(36)*Q14*Q2+4.0D0*C(37)*Q12*
     1   Q2*TQ1+C(38)*Q13*Q3*(-6.0D0*Q2)+C(39)*Q1*Q3*Q2*(2.0D0*TQ2-6.0D0
     2   *TQ1)+C(40)*3.0D0*TQ12*2.0D0*Q2+C(41)*Q32*2.0D0*TQ2*(-6.0D0*Q2)
     3   +C(43)*Q14*2.0D0*Q2+C(44)*Q13*Q3*2.0D0*Q2+C(45)*Q12*Q32*(-6.0D0
     4   *Q2)+C(46)*Q12*Q2*2.0D0*(TQ3+TQ1)+C(47)*Q1*Q3*2.0D0*TQ1*2.0D0*
     D   Q2+C(48)*Q1*Q3*Q2*(2.0D0*TQ2-6.0D0*TQ3)+C(49)*Q32*Q2*(2.0D0*TQ2
     E   -6.0D0*TQ1)+C(50)*2.0D0*Q2*TQ1*(2.0D0*TQ3+TQ1)
      DTQ=2.0D0*Q3
      DPOQ3=C(3)*DTQ+C(5)*Q1*DTQ+C(6)*TQ2+C(6)*Q3*DTQ+C(8)*Q12*DTQ+
     1   C(9)*TQ1*2.0D0*DTQ+C(10)*Q1*TQ2+C(10)*Q1*Q3*DTQ+C(11)+
     2   C(12)*Q1-C(13)*DTQ+C(14)*Q12+C(15)*Q1*(-DTQ)+C(16)*TQ1+
     3   C(16)*Q3*DTQ+C(17)*Q13+C(18)*Q12*(-DTQ)+C(19)*Q1*(TQ1+Q3*
     4   DTQ)+C(20)*2.0D0*Q3*TQ2+C(20)*Q32*DTQ+C(21)*(DTQ*TQ3-TQ1*DTQ)
     5   +C(24)*Q13*DTQ+C(25)*Q1*2.0D0*TQ1*DTQ+C(26)*Q12*(TQ2+Q3*DTQ)+
     6  C(27)*(TQ1*TQ2+Q3*DTQ*TQ2+Q3*TQ1*DTQ)+C(28)*Q14+C(29)*Q13*(-DTQ)
     7 +C(30)*Q12*(TQ1+Q3*DTQ)+C(31)*Q1*(2.0D0*Q3*TQ2+Q32*DTQ)+C(32)*Q1*
     8   (DTQ*TQ3-DTQ*TQ1)+C(33)*(TQ12+2.0D0*Q3*TQ1*DTQ)+C(34)*TQ2*TQ3+
     9   C(34)*Q3*DTQ*(TQ3-TQ2)+C(36)*Q14*DTQ+C(37)*Q12*2.0D0*TQ1*DTQ+
     A   C(38)*Q13*(TQ2+Q3*DTQ)+C(39)*Q1*TQ1*TQ2+C(39)*Q1*Q3*DTQ*
     B   (TQ1+TQ2)+C(40)*3.0D0*TQ12*DTQ+C(41)*DTQ*TQ2*(TQ2+2.0D0*Q32)+
     C   C(42)*Q15-C(43)*Q14*DTQ+C(44)*Q13*(TQ1+Q3*DTQ)+C(45)*Q12*DTQ*
     D  (TQ2+Q32)+C(46)*Q12*DTQ*(TQ3-TQ1)+C(47)*Q1*(TQ12+TQ1*4.0D0*Q32)+
     E   C(48)*Q1*TQ2*TQ3+C(48)*Q1*Q3*DTQ*(TQ3-TQ2)+C(49)*DTQ*TQ1*TQ2+
     F   C(49)*Q32*DTQ*(TQ1+TQ2)+C(50)*DTQ*(2.0D0*TQ1*TQ3-TQ12)
      CALL  DEREXDIS_10(X(1),X(2),X(3),DER1,DER2,DER3)
      CALL  DERELECT_10(X(1),X(2),X(3),DERI1,DERI2,DERI3)
      CALL  DERVSPECT_10(X(1),X(2),X(3),DERV1,DERV2,DERV3)
      gammar=3.0d0
      corr2=VOH_10(r2)-VOHP_10(r2)
      corr3=VOH_10(r3)-VOHP_10(r3)
      damp1=0.5d0*(1.0d0-tanh(gammar*(r1-5.0d0)))
      damp2=0.5d0*(1.0d0-tanh(gammar*(r2-5.0d0)))
      damp3=0.5d0*(1.0d0-tanh(gammar*(r3-5.0d0)))
      ddamp1=-0.5d0*gammar/cosh(gammar*(r1-5.0d0))**2
      ddamp2=-0.5d0*gammar/cosh(gammar*(r2-5.0d0))**2
      ddamp3=-0.5d0*gammar/cosh(gammar*(r3-5.0d0))**2


      G(1)=((DPOQ1*DQ1R1+DPOQ2*DQ2R1+DPOQ3*DQ3R1)*DECAY1-POLQ*C(51)/
     1   COSH(C(51)*S1)**2)*DECAY2*DECAY3+DER1+DERI1+DERV1+
     2   (corr2*damp3+corr3*damp2)*ddamp1
      G(2)=((DPOQ1*DQ1R2+DPOQ2*DQ2R2+DPOQ3*DQ3R2)*DECAY2-POLQ*C(52)/
     1   COSH(C(52)*S2)**2)*DECAY1*DECAY3+DER2+DERI2+DERV2+
     2   (dVOH_10(r2)-dVOHP_10(r2))*damp1*damp3+corr3*damp1*ddamp2
      G(3)=((DPOQ1*DQ1R3+DPOQ2*DQ2R3+DPOQ3*DQ3R3)*DECAY3-POLQ*C(52)/
     1   COSH(C(52)*S3)**2)*DECAY1*DECAY2+DER3+DERI3+DERV3+
     2   corr2*damp1*ddamp3+(dVOH_10(r3)-dVOHP_10(r3))*damp1*damp2


c      G(1)=((DPOQ1*DQ1R1+DPOQ2*DQ2R1+DPOQ3*DQ3R1)*DECAY1-POLQ*C(51)/
c     1   COSH(C(51)*S1)**2)*DECAY2*DECAY3+DER1+DERI1+DERV1
c      G(2)=((DPOQ1*DQ1R2+DPOQ2*DQ2R2+DPOQ3*DQ3R2)*DECAY2-POLQ*C(52)/
c     1   COSH(C(52)*S2)**2)*DECAY1*DECAY3+DER2+DERI2+DERV2
c      G(3)=((DPOQ1*DQ1R3+DPOQ2*DQ2R3+DPOQ3*DQ3R3)*DECAY3-POLQ*C(52)/
c     1   COSH(C(52)*S3)**2)*DECAY1*DECAY2+DER3+DERI3+DERV3

      dg1=G(1)
      dg2=G(2)
      dg3=G(3)

      END


      SUBROUTINE DERVSPECT_10(R1,R2,R3,DERV1,DERV2,DERV3)
C     ****************************************************************
C
C     TO COMPUTE THE DERIVATIVES OF A GAUSSIAN TERM TO FIX 
C     VIBRATIONAL SPECTRA
C
C     ****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Y)
     
      data r1eq,r2eq,r3eq/2.5143d0,1.8346d0,3.4592d0/
      
C     ****************************************************************
c

      r1d=r1-r1eq
      s1d=((r2+r3)-(r2eq+r3eq))/sqrt(2.0d0)
      s2d2=((r2-r3)**2-(r2eq-r3eq)**2)/2.0d0

      c1fix=2.0d0
      c2fix=2.0d0
      c3fix=2.0d0
      c4fix=2.0d0
     
      C1=1.0944219d-03 
      C2=2.9389603d-04
      C3=-9.4723705d-03
      C4=7.1459133d-02
      C5=2.2984321d-02
      C6=-4.8325504d-02
      C7=4.7234450d-04
      C8=-5.9777409d-02
      C9=1.8128239d-02
      C10=6.4184643d-02
      C11=-1.1384161d-02
      C12=1.2619628d-03
      C13=-2.9029563d-03


      t=c1+c2*r1d+c3*s1d+c4*r1d**2+c5*s1d**2+c6*r1d*s1d+
     1  c7*s2d2+c8*r1d**3+c9*r1d*s1d**2+c10*r1d**2*s1d+
     2  c11*r1d*s2d2+c12*s1d**3+c13*s1d*s2d2
      dec=c1fix*r1d**2+c2fix*s1d**2+c3fix*r1d*s1d+c4fix*s2d2**2

      VSPECT_10=t*exp(-dec)

  
C      VSPECT_10=t*exp(-dec)
 
       DERV1=(c2+2.0d0*c4*r1d+c6*s1d+
     2      3.0d0*c8*r1d**2+c9*s1d**2+2.0d0*c10*r1d*s1d+
     3      c11*s2d2
     1     -t*(2.0d0*c1fix*r1d+c3fix*s1d))*exp(-dec)

       DERVS1=(c3+2.0d0*c5*s1d+c6*r1d+
     2   2.0d0*c9*r1d*s1d+c10*r1d**2+
     3   3.0d0*c12*s1d**2+c13*s2d2
     1       -t*(2.0d0*c2fix*s1d+c3fix*r1d))*exp(-dec)

       DERVS2=(c7+c11*r1d+c13*s1d
     1   -t*c4fix*2.0D0*s2d2)*exp(-dec)

       DS1DR2=1.0D0/SQRT(2.0D0)
       DS1DR3=1.0D0/SQRT(2.0D0)
       DS2DR2=r2-r3
       DS2DR3=r3-r2
    
       DERV2=DERVS1*DS1DR2+DERVS2*DS2DR2
       DERV3=DERVS1*DS1DR3+DERVS2*DS2DR3

 
      RETURN
      END
      
      
C     ****************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  O...H
C     ****************************************************************
      FUNCTION DVOH_10(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DISPC_10/COO(10),COH(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      DATA D,ASV/0.13825385D0,2.6564788D0,1.7450528D0,0.71014391D0,
     1           2.5453276D0/
      SAVE /DISPC_10/,/DIATDIHO2_10/,D,ASV
      
      X=R-RMOH
      R2=X*X
      R3=R2*X
      POL=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)
      DPOL=-D*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*R2)
      POT=EXP(-ASV(4)*X)
      DVOH_10=-ASV(4)*POT*POL+DPOL*POT+DDISP_10(R,COH(6),COH(8),COH(10),
     1   R0OH,RMOH)
     
      END
      

C     ****************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  O...O
C     ****************************************************************
      FUNCTION DVOO_10(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DISPC_10/COO(10),COH(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      DATA D,ASV/0.14291202D0,3.6445906D0,3.9281238D0,2.0986689D0,
     1           3.3522498D0/
      SAVE /DISPC_10/,/DIATDIHO2_10/,D,ASV
      
      X=R-RMOO
      R2=X*X
      R3=R2*X
      POL=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)
      DPOL=-D*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*R2)
      POT=EXP(-ASV(4)*X)
      DVOO_10=-ASV(4)*POT*POL+DPOL*POT+DDISP_10(R,COO(6),COO(8),COO(10),
     1   R0OO,RMOO)

      END

      
C     ****************************************************************
C     TO COMPUTE THE DERIVATIVES OF THE EXCHANGE - DISPERSION TERM
C     ****************************************************************
      SUBROUTINE DEREXDIS_10(R1,R2,R3,DER1,DER2,DER3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC_10/COO(10),COH(10)
      COMMON/RKVAL_10/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DISCO_10/CEFOO(10),CEFOH2(10),CEFOH3(10),CEDOO(10),
     1 CEDOH2(10)
      COMMON/DISCO2_10/CEDOH3(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      COMMON/DAMPC_10/ADAMP(10),BDAMP(10)
!$OMP THREADPRIVATE(/DISCO_10/,/DISCO2_10/)

      SAVE /DISPC_10/,/DIATDIHO2_10/,/RKVAL_10/
      
      DO 10 IN=6,10,2
      CEFOO(IN)=CEF_10(COO(IN),RK0OH(IN),RK1OH(IN),RK0OH(IN),RK1OH(IN),
     1   RMOH,RMOH,R2,R3)
      CEDOO(IN)=CEFOO(IN)-COO(IN)
      CEFOH2(IN)=CEF_10(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN),
     1   RMOO,RMOH,R1,R3)
      CEDOH2(IN)=CEFOH2(IN)-COH(IN)
      CEFOH3(IN)=CEF_10(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN),
     1   RMOO,RMOH,R1,R2)
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
      CALL DCEF_10(COO(6),RK0OH(6),RK1OH(6),RK0OH(6),RK1OH(6),RMOH,
     1     RMOH,R2,R3,DC61R2,DC61R3)
      CALL DCEF_10(COO(8),RK0OH(8),RK1OH(8),RK0OH(8),RK1OH(8),RMOH,
     1     RMOH,R2,R3,DC81R2,DC81R3)
      CALL DCEF_10(COO(10),RK0OH(10),RK1OH(10),RK0OH(10),RK1OH(10),RMOH,
     1     RMOH,R2,R3,D101R2,D101R3)
      CALL DCEF_10(COH(6),RK0OO(6),RK1OO(6),RK0OH(6),RK1OH(6),RMOO,
     1     RMOH,R1,R3,DC62R1,DC62R3)
      CALL DCEF_10(COH(8),RK0OO(8),RK1OO(8),RK0OH(8),RK1OH(8),RMOO,
     1     RMOH,R1,R3,DC82R1,DC82R3)
      CALL DCEF_10(COH(10),RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10),RMOO,
     1     RMOH,R1,R3,D102R1,D102R3)
      CALL DCEF_10(COH(6),RK0OO(6),RK1OO(6),RK0OH(6),RK1OH(6),RMOO,
     1     RMOH,R1,R2,DC63R1,DC63R2)
      CALL DCEF_10(COH(8),RK0OO(8),RK1OO(8),RK0OH(8),RK1OH(8),RMOO,
     1     RMOH,R1,R2,DC83R1,DC83R2)
      CALL DCEF_10(COH(10),RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10),RMOO,
     1     RMOH,R1,R2,D103R1,D103R2)
      DER1=DDISP_10(R1,CEDOO(6),CEDOO(8),CEDOO(10),R0OO,RMOO)-DC62R1*
     1     T6R2-DC63R1*T6R3-DC82R1*T8R2-DC83R1*T8R3-D102R1*T10R2-
     2     D103R1*T10R3
      DER2=DDISP_10(R2,CEDOH2(6),CEDOH2(8),CEDOH2(10),R0OH,RMOH)-DC61R2*
     1     T6R1-DC63R2*T6R3-DC81R2*T8R1-DC83R2*T8R3-D101R2*T10R1-
     2     D103R2*T10R3
      DER3=DDISP_10(R3,CEDOH3(6),CEDOH3(8),CEDOH3(10),R0OH,RMOH)-DC61R3*
     1     T6R1-DC62R3*T6R2-DC81R3*T8R1-DC82R3*T8R2-D101R3*T10R1-
     2     D102R3*T10R2

      END
      

C     ****************************************************************
C     TO COMPUTE THE DERIVATIVES OF THE EFECTIVE Cn
C     ****************************************************************
      SUBROUTINE DCEF_10(CAS,RK01,RK11,RK02,RK12,RE1,RE2
     1 ,R1,R2,DCR1,DCR2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      DCR1=0.5D0*CAS*(RK11*RK01*EXP(-RK11*(R1-RE1))*TANH(RK12*R2)+
     1  (1.0D0-RK02*EXP(-RK12*(R2-RE2)))*RK11/COSH(RK11*R1)**2)
      DCR2=0.5D0*CAS*(RK12*RK02*EXP(-RK12*(R2-RE2))*TANH(RK11*R1)+
     1  (1.0D0-RK01*EXP(-RK11*(R1-RE1)))*RK12/COSH(RK12*R2)**2)
     
      END


C     ****************************************************************
C     TO COMPUTE THE DERIVATIVES OF THE ELECTROSTATIC TERM
C     ****************************************************************
      SUBROUTINE DERELECT_10(R1,R2,R3,DER1,DER2,DER3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/POLAR_10/C4,C5
      COMMON/RKVAL_10/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)
      COMMON/DIATDIHO2_10/R0OO,RMOO,R0OH,RMOH
      COMMON/DAMPC_10/ADAMP(10),BDAMP(10)
      SAVE /DIATDIHO2_10/,/RKVAL_10/,/POLAR_10/,/DAMPC_10/
      
C     RCR2=(R1-R3)/R2
C     RCR3=(R1-R2)/R3
C     C42=C4*RCR2
C     C43=C4*RCR3
      C42=C4
      C43=C4
C     PR2=3.0D0*RCR2**2-1.0D0
C     PR3=3.0D0*RCR3**2-1.0D0
C     C52=C5*PR2
C     C53=C5*PR3
      C52=C5
      C53=C5
C     RR2D1=1.0D0/R2
C     RR2D2=(R3-R1)/R2**2
C     RR2D3=-1.0D0/R2
C     RR3D1=1.0D0/R3
C     RR3D2=-1.0D0/R3
C     RR3D3=(R2-R1)/R3**2
      RR1=2.0D0*R1/(RMOO+2.5D0*R0OO)
      D4R1=1.0D0-EXP(-ADAMP(4)*RR1-BDAMP(4)*RR1**2)
      D5R1=1.0D0-EXP(-ADAMP(5)*RR1-BDAMP(5)*RR1**2)
      T4R1=(D4R1/R1)**4
      T5R1=(D5R1/R1)**5
      RR2=2.0D0*R2/(RMOH+2.5D0*R0OH)
      D4R2=1.0D0-EXP(-ADAMP(4)*RR2-BDAMP(4)*RR2**2)
      D5R2=1.0D0-EXP(-ADAMP(5)*RR2-BDAMP(5)*RR2**2)
      T4R2=(D4R2/R2)**4
      T5R2=(D5R2/R2)**5
      RR3=2.0D0*R3/(RMOH+2.5D0*R0OH)
      D4R3=1.0D0-EXP(-ADAMP(4)*RR3-BDAMP(4)*RR3**2)
      D5R3=1.0D0-EXP(-ADAMP(5)*RR3-BDAMP(5)*RR3**2)
      T4R3=(D4R3/R3)**4
      T5R3=(D5R3/R3)**5
      RMQ=RMOH**4
      C4OHR2=.5D0*C43/RMQ*R3**4*EXP(-RK1OH(4)*(R3-RMOH))*TANH(RK1OO(4)*
     1       R1)
      C5OHR2=.5D0*C53/RMQ*R3**4*EXP(-RK1OH(5)*(R3-RMOH))*TANH(RK1OO(5)*
     1       R1)
      C4OHR3=.5D0*C42/RMQ*R2**4*EXP(-RK1OH(4)*(R2-RMOH))*TANH(RK1OO(4)*
     1       R1)
      C5OHR3=.5D0*C52/RMQ*R2**4*EXP(-RK1OH(5)*(R2-RMOH))*TANH(RK1OO(5)*
     1       R1)
      C4OO=.5D0*C43/RMQ*R3**4*EXP(-RK1OH(4)*(R3-RMOH))*TANH(RK1OH(4)*R2)
     1   +.5D0*C42/RMQ*R2**4*EXP(-RK1OH(4)*(R2-RMOH))*TANH(RK1OH(4)*R3)
      C5OO=.5D0*C53/RMQ*R3**4*EXP(-RK1OH(5)*(R3-RMOH))*TANH(RK1OH(5)*R2)
     1   +.5D0*C52/RMQ*R2**4*EXP(-RK1OH(5)*(R2-RMOH))*TANH(RK1OH(5)*R3)
      DC42R3=C4OHR2*(4.0D0/R3-RK1OH(4))
      DC42R1=0.5D0*C4/RMQ*R3**4*EXP(-RK1OH(4)*(R3-RMOH))*(RK1OO(4)/
     1      COSH(RK1OO(4)*R1)**2)
      DC52R3=C5OHR2*(4.0D0/R3-RK1OH(5))
      DC52R1=0.5D0*C5/RMQ*R3**4*EXP(-RK1OH(5)*(R3-RMOH))*(RK1OO(5)/
     1      COSH(RK1OO(5)*R1)**2)
      DC43R2=C4OHR3*(4.0D0/R2-RK1OH(4))
      DC43R1=0.5D0*C4/RMQ*R2**4*EXP(-RK1OH(4)*(R2-RMOH))*(RK1OO(4)/
     1      COSH(RK1OO(4)*R1)**2)
      DC53R2=C5OHR3*(4.0D0/R2-RK1OH(5))
      DC53R1=0.5D0*C5/RMQ*R2**4*EXP(-RK1OH(5)*(R2-RMOH))*(RK1OO(5)/
     1      COSH(RK1OO(5)*R1)**2)
      DC42R2=0.0D0
      DC52R2=0.0D0
      DC43R3=0.0D0
      DC53R3=0.0D0
      DC41R3=0.5D0*C4/RMQ*R3**3*EXP(-RK1OH(4)*(R3-RMOH))*TANH(RK1OH(4)*
     1       R2)*(4.0D0-R3*RK1OH(4))+
     2.5D0*C4/RMQ*R2**4*EXP(-RK1OH(4)*(R2-RMOH))*(RK1OH(4)/COSH(RK1OH(4)
     3  *R3)**2)
      DC41R2=0.5D0*C4/RMQ*R2**3*EXP(-RK1OH(4)*(R2-RMOH))*TANH(RK1OH(4)*
     1       R3)*(4.0D0-R2*RK1OH(4))+
     2.5D0*C4/RMQ*R3**4*EXP(-RK1OH(4)*(R3-RMOH))*(RK1OH(4)/COSH(RK1OH(4)
     3  *R2)**2)
      DC51R3=0.5D0*C5/RMQ*R3**3*EXP(-RK1OH(5)*(R3-RMOH))*TANH(RK1OH(5)*
     1       R2)*(4.0D0-R3*RK1OH(5))+
     2.5D0*C5/RMQ*R2**4*EXP(-RK1OH(5)*(R2-RMOH))*(RK1OH(5)/COSH(RK1OH(5)
     3  *R3)**2)
      DC51R2=0.5D0*C5/RMQ*R2**3*EXP(-RK1OH(5)*(R2-RMOH))*TANH(RK1OH(5)*
     1       R3)*(4.0D0-R2*RK1OH(5))+
     2.5D0*C5/RMQ*R3**4*EXP(-RK1OH(5)*(R3-RMOH))*(RK1OH(5)/COSH(RK1OH(5)
     3  *R2)**2)
      DC41R1=0.0D0
      DC51R1=0.0D0
      DRR1=RR1/R1
      DRR2=RR2/R2
      DRR3=RR3/R3
      DDISP1=-4.0D0*C4OO/R1**4*D4R1**3*(D4R1/R1+(1.0D0-D4R1)*(-ADAMP(4)-
     1       2.0D0*BDAMP(4)*RR1)*DRR1)
     2      -5.0D0*C5OO/R1**5*D5R1**4*(D5R1/R1+(1.0D0-D5R1)*(-ADAMP(5)-
     3       2.0D0*BDAMP(5)*RR1)*DRR1)
      DDISP2=-4.0D0*C4OHR2/R2**4*D4R2**3*(D4R2/R2+(1.0D0-D4R2)*(-ADAMP
     1       (4)-2.0D0*BDAMP(4)*RR2)*DRR2)
     2    -5.0D0*C5OHR2/R2**5*D5R2**4*(D5R2/R2+(1.0D0-D5R2)*(-ADAMP(5)-
     3       2.0D0*BDAMP(5)*RR2)*DRR2)
      DDISP3=-4.0D0*C4OHR3/R3**4*D4R3**3*(D4R3/R3+(1.0D0-D4R3)*(-ADAMP
     1       (4)-2.0D0*BDAMP(4)*RR3)*DRR3)
     2     -5.0D0*C5OHR3/R3**5*D5R3**4*(D5R3/R3+(1.0D0-D5R3)*(-ADAMP(5)-
     3       2.0D0*BDAMP(5)*RR3)*DRR3)
      DER1=DDISP1+DC42R1*T4R2+DC43R1*T4R3+DC52R1*T5R2+DC53R1*T5R3
     1     +DC41R1*T4R1+DC51R1*T5R1
      DER2=DDISP2+DC41R2*T4R1+DC43R2*T4R3+DC51R2*T5R1+DC53R2*T5R3
     1     +DC42R2*T4R2+DC52R2*T5R2
      DER3=DDISP3+DC42R3*T4R2+DC41R3*T4R1+DC52R3*T5R2+DC51R3*T5R1
     1     +DC43R3*T4R3+DC53R3*T5R3

      END


      FUNCTION DDISP_10(R,C6,C8,C10,R0,RM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DAMPC_10/ADAMP(10),BDAMP(10)
      SAVE /DAMPC_10/
      
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      DRR=RR/R
      T6=1.0D0-EXP(-ADAMP(6)*RR-BDAMP(6)*RR**2)
      T8=1.0D0-EXP(-ADAMP(8)*RR-BDAMP(8)*RR**2)
      T10=1.0D0-EXP(-ADAMP(10)*RR-BDAMP(10)*RR**2)
      DDISP_10=6.0D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*
     1      (-ADAMP(6)-2.0D0*BDAMP(6)
     1      *RR)*DRR)+
     2      8.0D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-ADAMP(8)-2.0D0*BDAMP(8)
     3      *RR)*DRR)+
     4      10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-ADAMP(10)-2.0D0*
     5      BDAMP(10)*RR)*DRR)

      END



