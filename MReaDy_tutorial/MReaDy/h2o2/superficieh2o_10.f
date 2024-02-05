      Function vh2o_10(r1,r2,r3)
c    **************************************************************
C     TO COMPUTE THE H2O SURFACE - Sx in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
C    **************************************************************
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/CTHRB12_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q22,Q32,
     1 TQ1,TQ2,TQ3,POL,DECAY1,DECAY2,DECAY3,DECAY,DV1DR1,
     2 DV1DR2,DV1DR3,DV2DR1,DV2DR2,DV2DR3,F1,F2
!$OMP THREADPRIVATE(/DIATDI_10/,/CTHRB12_10/)

C    **************************************************************
      F1=VHHS_10(r1)+VOHS_10(r2)+VOHS_10(r3)+THREBR1_10(r1,r2,r3)+
     1  RL1_10(r1,r2,r3)+VO1D
      F2=VHHT_10(r1)+VOHP_10(r2)+VOHP_10(r3)+THREBR2_10(r1,r2,r3)+
     1  rl2_10(r1,r2,r3)
      F12=THREB12_10(r1,r2,r3)
      vh2o_10=((F1+F2)-SQRT((F1-F2)**2+4.0D0*F12**2))/2.0D0
!      print*,'vh2o_10'
      RETURN
      END
       
      Function vh2o3c1_10(r1,r2,r3)
c    **************************************************************
C     TO COMPUTE THE H2O SURFACE - Sx in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
C    **************************************************************
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/CTHRB12_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q22,Q32,
     1 TQ1,TQ2,TQ3,POL,DECAY1,DECAY2,DECAY3,DECAY,DV1DR1,
     2 DV1DR2,DV1DR3,DV2DR1,DV2DR2,DV2DR3,F1,F2
!$OMP THREADPRIVATE(/DIATDI_10/,/CTHRB12_10/)

c    ***************************************************************            
      vh2o3c1_10=THREBR1_10(r1,r2,r3)+RL1_10(r1,r2,r3)
      RETURN
      END
 
      Function vh2o3c2_10(r1,r2,r3)
c    **************************************************************
C     TO COMPUTE THE H2O SURFACE - Sx in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
C    **************************************************************
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/CTHRB12_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q22,Q32,
     1 TQ1,TQ2,TQ3,POL,DECAY1,DECAY2,DECAY3,DECAY,DV1DR1,
     2 DV1DR2,DV1DR3,DV2DR1,DV2DR2,DV2DR3,F1,F2
!$OMP THREADPRIVATE(/DIATDI_10/,/CTHRB12_10/)

c    ***************************************************************            
      vh2o3c2_10=THREBR2_10(r1,r2,r3)+rl2_10(r1,r2,r3)
      RETURN
      END

      Function vh2o3c12_10(r1,r2,r3)
c    **************************************************************
C     TO COMPUTE THE H2O SURFACE - Sx in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
C    **************************************************************
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/CTHRB12_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q22,Q32,
     1 TQ1,TQ2,TQ3,POL,DECAY1,DECAY2,DECAY3,DECAY,DV1DR1,
     2 DV1DR2,DV1DR3,DV2DR1,DV2DR2,DV2DR3,F1,F2
!$OMP THREADPRIVATE(/DIATDI_10/,/CTHRB12_10/)

c    ***************************************************************            
      F1=VHHS_10(r1)+VOHS_10(r2)+VOHS_10(r3)+THREBR1_10(r1,r2,r3)+
     1  RL1_10(r1,r2,r3)+VO1D
      F2=VHHT_10(r1)+VOHP_10(r2)+VOHP_10(r3)+THREBR2_10(r1,r2,r3)+
     1  rl2_10(r1,r2,r3)
c     os valores de F1 e F2 são necessárs para o cálculo do termo V12   
      vh2o3c12_10=THREB12_10(r1,r2,r3)

      RETURN
      END



    
      SUBROUTINE DERvh2o_10(r1,r2,r3,g1,g2,g3)
c    **************************************************************
C     TO COMPUTE THE DERIVATIVE H2O SURFACE - Sx in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
      DIMENSION GF1(3),GF2(3)      
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
C    **************************************************************
      COMMON/CTHRB12_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q22,Q32,
     1 TQ1,TQ2,TQ3,POL,DECAY1,DECAY2,DECAY3,DECAY,DV1DR1,
     2 DV1DR2,DV1DR3,DV2DR1,DV2DR2,DV2DR3,F1,F2
!$OMP THREADPRIVATE(/DIATDI_10/,/CTHRB12_10/)

c    ***************************************************************            
      F1=VHHS_10(r1)+VOHS_10(r2)+VOHS_10(r3)+THREBR1_10(r1,r2,r3)+
     1  RL1_10(r1,r2,r3)+VO1D
      F2=VHHT_10(r1)+VOHP_10(r2)+VOHP_10(r3)+THREBR2_10(r1,r2,r3)+
     1  rl2_10(r1,r2,r3)
      F12=THREB12_10(r1,r2,r3)
      
      CALL DTHREBR1_10(r1,r2,r3,DTDR1,DTDR2,DTDR3)
      CALL DRL1_10(r1,r2,r3,DER1,DER2,DER3)
      GF1(1)=DTDR1+DER1+DVHHS_10(r1)
      GF1(2)=DTDR2+DER2+DVOHS_10(r2)
      GF1(3)=DTDR3+DER3+DVOHS_10(r3)
      DV1DR1=GF1(1)
      DV1DR2=GF1(2)
      DV1DR3=GF1(3)
      CALL DTHREBR2_10(r1,r2,r3,DT2DR1,DT2DR2,DT2DR3)
      CALL Drl2_10(r1,r2,r3,DER1,DER2,DER3)
      GF2(1)=DT2DR1+DER1+DVHHT_10(r1)
      GF2(2)=DT2DR2+DER2+DVOHP_10(r2)
      GF2(3)=DT2DR3+DER3+DVOHP_10(r3)
      DV2DR1=GF2(1)
      DV2DR2=GF2(2)
      DV2DR3=GF2(3)
      CALL DTHREB12_10(r1,r2,r3,DT12R1,DT12R2,DT12R3)
      DFDF1=0.5D0*(1.0D0-(F1-F2)/SQRT((F1-F2)**2+4.0D0*F12**2))
      DFDF2=0.5D0*(1.0D0+(F1-F2)/SQRT((F1-F2)**2+4.0D0*F12**2))
      DFDF12=-2.0D0*F12/SQRT((F1-F2)**2+4.0D0*F12**2)
      G1=DFDF1*GF1(1)+DFDF2*GF2(1)+DFDF12*DT12R1
      G2=DFDF1*GF1(2)+DFDF2*GF2(2)+DFDF12*DT12R2
      G3=DFDF1*GF1(3)+DFDF2*GF2(3)+DFDF12*DT12R3

      RETURN
      END

!      SUBROUTINE DERVH2O3c(r1,r2,r3,g1,g2,g3)
c     **************************************************************
C      TO COMPUTE THE DERIVATIVE H2O SURFACE - Sx in atomic units
C     **************************************************************
!      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
!      DIMENSION GF1(3),GF2(3)      
!      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
C     **************************************************************
!      COMMON/CTHRB12_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q22,Q32,
!     1 TQ1,TQ2,TQ3,POL,DECAY1,DECAY2,DECAY3,DECAY,DV1DR1,
!     2 DV1DR2,DV1DR3,DV2DR1,DV2DR2,DV2DR3,F1,F2
c     ***************************************************************            
!      CALL DERvh2o_10(r1,r2,r3,dvh2o1,dvh2o2_10,dvh2o3)
!      g1=dvh2o1-dVHHS_10(r1)
!      g2=dvh2o2-dVOHP_10(r2)
!      g3=dvh2o3-dVOHP_10(r3)

c     vh2o3c=vh2o_10-VHHS_10(r1)-VOHP_10(r2)-VOHP_10(r3)

!      RETURN
!      END
      
      SUBROUTINE DTHREBR1_10(R1,R2,R3,DTDR1,DTDR2,DTDR3)
c    **************************************************************
C     TO COMPUTE THE DERIVATIVES OF THE THREE BODY TERM IN 
c     COORDINATES R1,R2,R3
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(nc1=95,nc2=34,nc3=15,nc=144)
      COMMON/COEFF_10/C(148)
      COMMON/REFGEO_10/R10,R20,R30
      COMMON/CTHRB1_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q14,Q15,
     1 Q16,Q17,Q18,Q22,Q23,Q32,Q33,TQ1,TQ12,TQ13,TQ14,TQ2,TQ3,TQ32,
     2 POLQ,DECAY1,DECAY2,DECAY3
!$OMP THREADPRIVATE(/CTHRB1_10/)

c    **************************************************************
      prov=THREBR1_10(r1,r2,r3)


      DECAY=DECAY1*DECAY2*DECAY3
      DQ1R1=RO3
      DQ1R2=RO3
      DQ1R3=RO3
      DQ2R1=0.0D0
      DQ2R2=RO2
      DQ2R3=-RO2
      DQ3R1=2.0D0*RO6
      DQ3R2=-RO6
      DQ3R3=-RO6
      
      DDECAYR1=DECAY2*DECAY3*C(145)/COSH(S1*C(145))**2
      DDECAYR2=DECAY1*DECAY3*C(146)/COSH(S2*C(146))**2  
      DDECAYR3=DECAY1*DECAY2*C(146)/COSH(S3*C(146))**2

      DTQ1Q2=2.0D0*Q2
      DTQ1Q3=2.0D0*Q3
      DTQ2Q2=2.0D0*Q2
      DTQ2Q3=-2.0D0*Q3
      DTQ3Q2=-6.0D0*Q2*Q3
      DTQ3Q3=3.0D0*Q32-3.0d0*q22
      
      DPOQ1=C(2)+2.0D0*C(4)*Q1+C(6)*Q3+3.0D0*C(8)*Q12+C(9)*TQ1+
     1      2.0D0*C(11)*Q1*Q3+C(12)*TQ2+4.0D0*C(14)*Q13+
     2      2.0D0*C(15)*Q1*TQ1+C(17)*TQ3+3.0D0*C(18)*Q12*Q3+
     3      2.0D0*C(19)*Q1*TQ2+C(20)*Q3*TQ1+5.0D0*C(23)*Q14+
     4      3.0D0*C(24)*Q12*TQ1+C(25)*TQ1**2+2.0D0*C(26)*Q1*TQ3+
     5      4.0D0*C(28)*Q13*Q3+3.0D0*C(29)*Q12*TQ2+2.0D0*C(30)*
     6      Q1*Q3*TQ1+C(31)*Q3*TQ3+C(32)*TQ1*TQ2+
     9 6.0D0*C(35)*Q15+4.0D0*C(36)*Q13*TQ1+2.0D0*C(37)*Q1*TQ12+
     1 3.0D0*C(38)*Q12*TQ3+C(39)*TQ3*TQ1+5.0D0*C(42)*Q14*Q3+4.0D0*
     2 C(43)*Q13*TQ2+3.0D0*C(44)*Q12*Q3*TQ1+2.0D0*C(45)*Q1*Q3*TQ3+
     3 2.0D0*C(46)*Q1*TQ1*TQ2+C(47)*Q3*TQ12+C(48)*TQ3*TQ2     
     4 +7.0d0*C(51)*Q16+5.0d0*C(52)*Q14*TQ1+3.0d0*C(53)*Q12*TQ12+
     5 4.0d0*C(54)*Q13*TQ3+2.0d0*C(55)*Q1*TQ1*TQ3+C(56)*TQ13+C(57)
     6 *TQ32+6.0d0*C(59)*Q3*Q15+4.0d0*C(60)*Q3*Q13*TQ1+2.0d0*C(61)
     7 *Q3*Q1*TQ12+3.0d0*C(62)*Q12*Q3*TQ3+C(63)*Q3*TQ3*TQ1+5.0d0*
     8 C(66)*Q14*TQ2+3.0D0*C(67)*Q12*TQ1*TQ2+2.0D0*C(68)*Q1*TQ3*
     9 TQ2+C(69)*TQ12*TQ2+
     9 8.0D0*c(71)*Q17+6.0D0*C(72)*Q15*TQ1+4.0D0*C(73)*Q13*TQ12+
     1 5.0D0*C(74)*Q14*TQ3+3.0D0*C(75)*Q12*TQ3*TQ1+2.0D0*C(76)*Q1*
     2 TQ13+2.0D0*C(77)*Q1*TQ32+C(78)*TQ3*TQ12+7.0D0*C(81)*Q16*Q3
     3 +5.0D0*C(82)*Q14*Q3*TQ1+3.0D0*C(83)*Q12*Q3*TQ12+4.0D0*C(84)
     4 *Q13*TQ32+2.0D0*C(85)*Q1*TQ32*TQ1+C(86)*Q3*TQ13+C(87)*
     5 Q3*TQ32+6.0D0*C(88)*Q15*TQ2+4.0D0*C(89)*Q13*TQ2*TQ1+3.0D0*
     6 C(90)*Q12*TQ3*TQ2+2.0D0*C(91)*Q1*TQ12*TQ2+C(92)*TQ3*TQ2*TQ1
     
      DPOQ2=2.0D0*C(5)*Q2+2.0D0*C(7)*Q2+2.0D0*C(9)*Q2*Q1+
     1      2.0D0*C(12)*Q1*Q2-6.0D0*C(10)*Q2*Q3+2.0D0*C(13)*
     2      Q2*Q3+2.0d0*C(15)*Q12*Q2+4.0D0*C(16)*Q23+4.0D0*C(16)*
     3      Q2*Q32-6.0D0*C(17)*Q1*Q3*Q2+2.0D0*C(19)*Q12*Q2+
     4      2.0D0*C(20)*Q1*Q2*Q3-6.0D0*C(21)*Q32*Q2+4.0D0*C(22)*Q23+
     5      2.0D0*C(24)*Q2*Q13+2.0D0*C(25)*Q1*DTQ1Q2*TQ1-6.0D0*
     6      C(26)*Q12*Q2*Q3+C(27)*(DTQ1Q2*TQ3+DTQ3Q2*TQ1)+
     7      C(29)*Q13*DTQ2Q2+C(30)*Q12*Q3*DTQ1Q2+C(31)*Q1*Q3*DTQ3Q2+
     8      C(32)*Q1*(DTQ1Q2*TQ2+DTQ2Q2*TQ1)+2.0D0*C(33)*Q3*DTQ1Q2*TQ1+
     9      C(34)*(DTQ2Q2*TQ3+DTQ3Q2*TQ2)+          
     9 C(36)*Q14*DTQ1Q2+2.0D0*C(37)*Q12*DTQ1Q2*TQ1+C(38)*Q13*DTQ3Q2+
     1 C(39)*Q1*(DTQ3Q2*TQ1+DTQ1Q2*TQ3)+3.0D0*C(40)*TQ12*DTQ1Q2+2.0D0*
     2 C(41)*DTQ3Q2*TQ3+C(43)*Q14*DTQ2Q2+C(44)*Q13*Q3*DTQ1Q2+C(45)*
     3 Q12*Q3*DTQ3Q2+C(46)*Q12*(DTQ1Q2*TQ2+DTQ2Q2*TQ1)+2.0D0*C(47)*Q1
     4 *Q3*DTQ1Q2*TQ1+C(48)*Q1*(DTQ3Q2*TQ2+DTQ2Q2*TQ3)+C(49)*Q3*
     5 (DTQ1Q2*TQ3+DTQ3Q2*TQ1)+C(50)*(DTQ2Q2*TQ12+2.0D0*TQ1*DTQ1Q2*TQ2)
     4 +C(52)*Q15*DTQ1Q2+2.0D0*C(53)*Q13*TQ1*DTQ1Q2+C(54)*Q14*DTQ3Q2+
     5 C(55)*Q12*(DTQ1Q2*TQ3+DTQ3Q2*TQ1)+3.0D0*C(56)*Q1*TQ12*DTQ1Q2+
     6 2.0D0*C(57)*Q1*TQ3*DTQ3Q2+C(58)*(DTQ3Q2*TQ12+2.0D0*TQ1*DTQ1Q2*
     7 TQ3)+C(60)*Q3*Q14*DTQ1Q2+2.0D0*C(61)*Q3*Q12*TQ1*DTQ1Q2+C(62)*
     8 Q13*Q3*DTQ3Q2+C(63)*Q1*Q3*(DTQ3Q2*TQ1+DTQ1Q2*TQ3)+3.0D0*C(64)*
     9 Q3*TQ12*DTQ1Q2+2.0D0*C(65)*Q3*TQ3*DTQ3Q2+C(66)*Q15*DTQ2Q2+C(67)
     1 *Q13*(DTQ1Q2*TQ2+DTQ2Q2*TQ1)+C(68)*Q12*(DTQ3Q2*TQ2+DTQ2Q2*TQ3)+
     2 C(69)*Q1*(2.0D0*TQ1*DTQ1Q2*TQ2+DTQ2Q2*TQ12)+C(70)*((DTQ3Q2*TQ1+
     3 DTQ1Q2*TQ3)*TQ2+DTQ2Q2*TQ1*TQ3)     
     9 +C(72)*Q16*DTQ1Q2+C(73)*Q14*2.0D0*TQ1*DTQ1Q2+C(74)*Q15*DTQ3Q2+
     1 C(75)*Q13*(DTQ3Q2*TQ1+DTQ1Q2*TQ3)+C(76)*Q12*3.0D0*TQ12*DTQ1Q2+
     2 C(77)*Q12*2.0D0*TQ3*DTQ3Q2+C(78)*Q1*(DTQ3Q2*TQ12+2.0D0*TQ1*
     3 DTQ1Q2*TQ3)+C(79)*4.0D0*TQ13*DTQ1Q2+C(80)*(2.0D0*TQ3*DTQ3Q2*
     4 TQ1+DTQ1Q2*TQ32)+C(82)*Q15*Q3*DTQ1Q2+C(83)*Q13*Q3*2.0D0*
     3 TQ1*DTQ1Q2+C(84)*Q14*2.0D0*TQ3*DTQ3Q2+C(85)*Q12*(2.0D0*DTQ3Q2*
     4 TQ3*TQ1+DTQ1Q2*TQ32)+C(86)*Q1*Q3*3.0D0*TQ12*DTQ1Q2+C(87)*
     5 Q1*Q3*2.0D0*TQ3*DTQ3Q2+C(88)*Q16*DTQ2Q2+C(89)*Q14*(DTQ2Q2*TQ1
     6 +DTQ1Q2*TQ2)+C(90)*Q13*(DTQ3Q2*TQ2+DTQ2Q2*TQ3)+
     6 C(91)*Q12*(2.0D0*TQ1*DTQ1Q2*TQ2+DTQ2Q2*TQ12)+C(92)*Q1*((DTQ3Q2
     7 *TQ2+DTQ2Q2*TQ3)*TQ1+DTQ1Q2*TQ3*TQ2)+C(93)*Q3*(DTQ3Q2*TQ12+
     8 2.0D0*TQ1*DTQ1Q2*TQ3)+C(94)*(2.0D0*TQ3*DTQ3Q2*TQ2+DTQ2Q2*TQ32)
     9 +C(95)*(3.0D0*TQ12*DTQ1Q2*TQ2+DTQ2Q2*TQ13)

     
      DPOQ3=C(3)+2.0D0*C(5)*Q3+C(6)*Q1-2.0D0*C(7)*Q3+2.0D0*C(9)*Q1*Q3+
     1      3.0D0*C(10)*Q32-3.0D0*C(10)*Q22+C(11)*Q12-2.0D0*C(12)*
     2      Q1*Q3+C(13)*Q22+3.0D0*C(13)*Q32+2.0D0*C(15)*Q12*Q3+
     3      4.0D0*C(16)*Q33+4.0D0*C(16)*Q22*Q3+3.0D0*C(17)*Q1*Q32-
     4      3.0D0*C(17)*Q1*Q22+
     5      C(18)*Q13-2.0D0*C(19)*Q12*Q3+C(20)*Q1*Q22+3.0D0*C(20)*
     6      Q1*Q32+4.0D0*C(21)*Q33-6.0D0*C(21)*Q3*Q22-4.0D0*C(22)*Q33+
     7      C(24)*Q13*DTQ1Q3+2.0D0*C(25)*Q1*DTQ1Q3*TQ1+
     8      C(26)*Q12*DTQ3Q3+C(27)*(DTQ1Q3*TQ3+DTQ3Q3*TQ1)+C(28)*
     9      Q14+C(29)*Q13*DTQ2Q3+C(30)*Q12*(TQ1+DTQ1Q3*Q3)+
     1      C(31)*Q1*(TQ3+DTQ3Q3*Q3)+C(32)*Q1*(DTQ1Q3*TQ2+DTQ2Q3*TQ1)+
     2      C(33)*(TQ1**2+2.0D0*TQ1*DTQ1Q3*Q3)+C(34)*(DTQ2Q3*TQ3+
     3      DTQ3Q3*TQ2)+
     9 C(36)*Q14*DTQ1Q3+2.0D0*C(37)*Q12*TQ1*DTQ1Q3+C(38)*Q13*DTQ3Q3+
     1 C(39)*Q1*(DTQ3Q3*TQ1+DTQ1Q3*TQ3)+3.0D0*C(40)*TQ12*DTQ1Q3+2.0D0*
     2 C(41)*TQ3*DTQ3Q3+C(42)*Q15+C(43)*Q14*DTQ2Q3+C(44)*
     3 Q13*(TQ1+Q3*DTQ1Q3)+C(45)*Q12*(TQ3+Q3*DTQ3Q3)+C(46)*Q12*
     4 (DTQ1Q3*TQ2+DTQ2Q3*TQ1)+C(47)*Q1*(TQ12+2.0D0*Q3*TQ1*DTQ1Q3)+
     5 C(48)*Q1*(DTQ3Q3*TQ2+DTQ2Q3*TQ3)+C(49)*(TQ3*TQ1+(DTQ1Q3*TQ3+
     6 DTQ3Q3*TQ1)*Q3)+C(50)*(DTQ2Q3*TQ12+2.0D0*TQ1*DTQ1Q3*TQ2)     
     4 +C(52)*Q15*DTQ1Q3+2.0D0*C(53)*Q13*TQ1*DTQ1Q3+C(54)*Q14*DTQ3Q3+
     5 C(55)*Q12*(DTQ1Q3*TQ3+DTQ3Q3*TQ1)+3.0D0*C(56)*Q1*TQ12*DTQ1Q3+
     6 2.0D0*C(57)*Q1*TQ3*DTQ3Q3+C(58)*(DTQ3Q3*TQ12+2.0D0*TQ1*DTQ1Q3*
     7 TQ3)+C(59)*Q16+C(60)*Q14*(TQ1+Q3*DTQ1Q3)+C(61)*Q12*(TQ12+Q3*
     8 2.0D0*TQ1*DTQ1Q3)+C(62)*Q13*(TQ3+Q3*DTQ3Q3)+C(63)*Q1*((TQ1+
     9 DTQ1Q3*Q3)*TQ3+DTQ3Q3*Q3*TQ1)+C(64)*(TQ13+Q3*3.0D0*TQ12*
     1 DTQ1Q3)+C(65)*(TQ32+2.0D0*Q3*TQ3*DTQ3Q3)+C(66)*Q15*DTQ2Q3+
     2 C(67)*Q13*(DTQ1Q3*TQ2+DTQ2Q3*TQ1)+C(68)*Q12*(DTQ3Q3*TQ2+
     3 DTQ2Q3*TQ3)+C(69)*Q1*(2.0D0*DTQ1Q3*TQ1*TQ2+DTQ2Q3*TQ12)+C(70)*
     4 ((DTQ3Q3*TQ1+DTQ1Q3*TQ3)*TQ2+DTQ2Q3*TQ1*TQ3)
     5 +C(72)*Q16*DTQ1Q3+C(73)*Q14*2.0D0*TQ1*DTQ1Q3+C(74)*Q15*DTQ3Q3
     6 +C(75)*Q13*(DTQ3Q3*TQ1+DTQ1Q3*TQ3)+C(76)*Q12*3.0D0*TQ12*
     7 DTQ1Q3+C(77)*Q12*2.0D0*TQ3*DTQ3Q3+C(78)*Q1*(DTQ3Q3*TQ12+2.0D0*
     8 TQ1*DTQ1Q3*tq3)+C(79)*4.0D0*TQ13*DTQ1Q3+C(80)*(DTQ3Q3*TQ3*2.0D0*
     9 TQ1+DTQ1Q3*TQ32)+C(81)*Q17+C(82)*Q15*(TQ1+Q3*DTQ1Q3)+C(83)*
     1 Q13*(2.0D0*TQ1*DTQ1Q3*Q3+TQ12)+C(84)*Q14*2.0D0*TQ3*DTQ3Q3+
     2 C(85)*Q12*(2.0D0*TQ3*DTQ3Q3*TQ1+DTQ1Q3*TQ32)+C(86)*Q1*(Q3*
     3 3.0D0*TQ12*DTQ1Q3+TQ13)+C(87)*Q1*(Q3*2.0D0*TQ3*DTQ3Q3+TQ32)+
     4 C(88)*Q16*DTQ2Q3+C(89)*Q14*(DTQ2Q3*TQ1+DTQ1Q3*TQ2)+C(90)*Q13*
     5 (DTQ3Q3*TQ2+DTQ2Q3*TQ3)+C(91)*Q12*(2.0D0*TQ1*DTQ1Q3*TQ2+
     6 DTQ2Q3*TQ12)+C(92)*Q1*((DTQ3Q3*TQ2+DTQ2Q3*TQ3)*TQ1+DTQ1Q3*
     7 TQ2*TQ3)+C(93)*((Q3*DTQ3Q3+TQ3)*TQ12+2.0D0*TQ1*DTQ1Q3*Q3*TQ3)
     8 +C(94)*(2.0D0*TQ3*DTQ3Q3*TQ2+DTQ2Q3*TQ32)+C(95)*(3.0D0*TQ12*
     9 DTQ1Q3*TQ2+DTQ2Q3*TQ13)
    
      DTDR1=(DPOQ1*DQ1R1+DPOQ3*DQ3R1)*DECAY-POLQ*DDECAYR1
      DTDR2=(DPOQ1*DQ1R2+DPOQ2*DQ2R2+DPOQ3*DQ3R2)*DECAY-POLQ*DDECAYR2
      DTDR3=(DPOQ1*DQ1R3+DPOQ2*DQ2R3+DPOQ3*DQ3R3)*DECAY-POLQ*DDECAYR3
      RETURN
      END

      SUBROUTINE DTHREBR2_10(R1,R2,R3,DT2DR1,DT2DR2,DT2DR3)
c    **************************************************************
C     TO COMPUTE THE DERIVATIVES OF THE THREE BODY TERM IN 
c     COORDINATES R1,R2,R3
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(nc1=95,nc2=34,nc3=15,nc=144)
      COMMON/COEFF_10/C(148)
      DIMENSION D(NC2+19)
      EQUIVALENCE (C(NC1+1),D(1))
      COMMON/REFGEO2_10/R102,R202,R302
      COMMON/CTHRB2_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q14,Q15,
     1 Q16,Q17,Q18,Q22,Q23,Q32,Q33,TQ1,TQ12,TQ13,TQ14,TQ2,TQ3,TQ32,
     2 POLQ,DECAY1,DECAY2,DECAY3
!$OMP THREADPRIVATE(/CTHRB2_10/)

c    **************************************************************
      prov=THREBR2_10(r1,r2,r3)




      DECAY=DECAY1*DECAY2*DECAY3
      DQ1R1=RO3
      DQ1R2=RO3
      DQ1R3=RO3
      DQ2R1=0.0D0
      DQ2R2=RO2
      DQ2R3=-RO2
      DQ3R1=2.0D0*RO6
      DQ3R2=-RO6
      DQ3R3=-RO6
      
      DDECAYR1=DECAY2*DECAY3*D(52)/COSH(S1*D(52))**2
      DDECAYR2=DECAY1*DECAY3*D(53)/COSH(S2*D(53))**2  
      DDECAYR3=DECAY1*DECAY2*D(53)/COSH(S3*D(53))**2

      DTQ1Q2=2.0D0*Q2
      DTQ1Q3=2.0D0*Q3
      DTQ2Q2=2.0D0*Q2
      DTQ2Q3=-2.0D0*Q3
      DTQ3Q2=-6.0D0*Q2*Q3
      DTQ3Q3=3.0D0*Q32-3.0d0*q22

      DPOQ1=D(2)+2.0D0*D(4)*Q1+D(6)*Q3+3.0D0*D(8)*Q12+D(9)*TQ1+
     1      2.0D0*D(11)*Q1*Q3+D(12)*TQ2+4.0D0*D(14)*Q13+
     2      2.0D0*D(15)*Q1*TQ1+D(17)*TQ3+3.0D0*D(18)*Q12*Q3+
     3      2.0D0*D(19)*Q1*TQ2+D(20)*Q3*TQ1+5.0D0*D(23)*Q14+
     4      3.0D0*D(24)*Q12*TQ1+D(25)*TQ1**2+2.0D0*D(26)*Q1*TQ3+
     5      4.0D0*D(28)*Q13*Q3+3.0D0*D(29)*Q12*TQ2+2.0D0*D(30)*
     6      Q1*Q3*TQ1+D(31)*Q3*TQ3+D(32)*TQ1*TQ2
     
      DPOQ2=2.0D0*D(5)*Q2+2.0D0*D(7)*Q2+2.0D0*D(9)*Q2*Q1+
     1      2.0D0*D(12)*Q1*Q2-6.0D0*D(10)*Q2*Q3+2.0D0*D(13)*
     2      Q2*Q3+2.0d0*D(15)*Q12*Q2+4.0D0*D(16)*Q23+4.0D0*D(16)*
     3      Q2*Q32-6.0D0*D(17)*Q1*Q3*Q2+2.0D0*D(19)*Q12*Q2+
     4      2.0D0*D(20)*Q1*Q2*Q3-6.0D0*D(21)*Q32*Q2+4.0D0*D(22)*Q23+
     5      2.0D0*D(24)*Q2*Q13+2.0D0*D(25)*Q1*DTQ1Q2*TQ1-6.0D0*
     6      D(26)*Q12*Q2*Q3+D(27)*(DTQ1Q2*TQ3+DTQ3Q2*TQ1)+
     7      D(29)*Q13*DTQ2Q2+D(30)*Q12*Q3*DTQ1Q2+D(31)*Q1*Q3*DTQ3Q2+
     8      D(32)*Q1*(DTQ1Q2*TQ2+DTQ2Q2*TQ1)+2.0D0*D(33)*Q3*DTQ1Q2*TQ1+
     9      D(34)*(DTQ2Q2*TQ3+DTQ3Q2*TQ2)
     
      DPOQ3=D(3)+2.0D0*D(5)*Q3+D(6)*Q1-2.0D0*D(7)*Q3+2.0D0*D(9)*Q1*Q3+
     1      3.0D0*D(10)*Q32-3.0D0*D(10)*Q22+D(11)*Q12-2.0D0*D(12)*
     2      Q1*Q3+D(13)*Q22+3.0D0*D(13)*Q32+2.0D0*D(15)*Q12*Q3+
     3      4.0D0*D(16)*Q33+4.0D0*D(16)*Q22*Q3+3.0D0*D(17)*Q1*Q32-
     4      3.0D0*D(17)*Q1*Q22+
     5      D(18)*Q13-2.0D0*D(19)*Q12*Q3+D(20)*Q1*Q22+3.0D0*D(20)*
     6      Q1*Q32+4.0D0*D(21)*Q33-6.0D0*D(21)*Q3*Q22-4.0D0*D(22)*Q33+
     7      D(24)*Q13*DTQ1Q3+2.0D0*D(25)*Q1*DTQ1Q3*TQ1+
     8      D(26)*Q12*DTQ3Q3+D(27)*(DTQ1Q3*TQ3+DTQ3Q3*TQ1)+D(28)*
     9      Q14+D(29)*Q13*DTQ2Q3+D(30)*Q12*(TQ1+DTQ1Q3*Q3)+
     1      D(31)*Q1*(TQ3+DTQ3Q3*Q3)+D(32)*Q1*(DTQ1Q3*TQ2+DTQ2Q3*TQ1)+
     2      D(33)*(TQ1**2+2.0D0*TQ1*DTQ1Q3*Q3)+D(34)*(DTQ2Q3*TQ3+
     3      DTQ3Q3*TQ2) 
    
      DT2DR1=(DPOQ1*DQ1R1+DPOQ3*DQ3R1)*DECAY-POLQ*DDECAYR1
      DT2DR2=(DPOQ1*DQ1R2+DPOQ2*DQ2R2+DPOQ3*DQ3R2)*DECAY-POLQ*DDECAYR2
      DT2DR3=(DPOQ1*DQ1R3+DPOQ2*DQ2R3+DPOQ3*DQ3R3)*DECAY-POLQ*DDECAYR3
      RETURN
      END      
      
      SUBROUTINE DTHREB12_10(R1,R2,R3,DT12R1,DT12R2,DT12R3)
c    **************************************************************
C     TO COMPUTE THE THREE BODY TERM IN COORDINATES R1,R2,R3
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(nc1=95,nc2=34,nc3=15,nc=144)
      DIMENSION d(NC3)
c      DIMENSION X(3),G(3)
c      DATA VO1D,EMIN/7.1955D-2,-0.3704003D0/
      COMMON/REFGEO12_10/R1012,R2012,R3012
      COMMON/COEFF_10/C(148)
      COMMON/CTHRB12_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q22,Q32,
     1 TQ1,TQ2,TQ3,POL,DECAY1,DECAY2,DECAY3,DECAY,DV1DR1,
     2 DV1DR2,DV1DR3,DV2DR1,DV2DR2,DV2DR3,V1,V2
!$OMP THREADPRIVATE(/CTHRB12_10/)




      EQUIVALENCE (C(NC1+NC2+1),d(1)) 
c    **************************************************************


c      provb=THREB12_10(r1,r2,r3)

      v1=VHHS_10(r1)+VOHS_10(r2)+VOHS_10(r3)+THREBR1_10(r1,r2,r3)+
     1  RL1_10(r1,r2,r3)+VO1D
      v2=VHHT_10(r1)+VOHP_10(r2)+VOHP_10(r3)+THREBR2_10(r1,r2,r3)+
     1  rl2_10(r1,r2,r3)



c      do i=1,nc3
c          d(i)=c(nc1+nc2+i)
c      end do

      term=(r2**2+r3**2-r1**2)/(2.0D0*r2*r3)
      
      if (abs(term).ge.1.0d0) then
         sina=0.0d0
       else
         sina=(1.0d0-term**2)
      end if
      

      RO2=1.0D0/SQRT(2.0D0)
      RO3=1.0D0/SQRT(3.0D0)
      RO6=1.0D0/SQRT(6.0D0)

      S1=R1-r1012
      S2=R2-r2012
      S3=R3-r3012
      Q1=RO3*(s1+s2+s3)
      Q2=RO2*(s2-s3)
      Q3=RO6*(2.0D0*s1-s2-s3)
      Q12=Q1*Q1
      Q13=Q12*Q1
      Q22=Q2*Q2
      Q32=Q3*Q3
      TQ1=Q22+Q32
      TQ2=Q22-Q32
      TQ3=(Q32-3.0D0*Q22)*Q3
      
      pol=D(1)+D(2)*Q1+D(3)*Q3+D(4)*Q12+
     1   D(5)*TQ1+D(6)*Q1*Q3+D(7)*TQ2+D(8)*Q13+D(9)*Q1*TQ1+
     2   D(10)*TQ3+D(11)*Q12*Q3+D(12)*Q1*TQ2+D(13)*Q3*TQ1

C      
C     ESTES DVS SAO POSTOS NO COMMON PARA SEREM USADOS NA DTHRB12
C
      GAMA1=100.0D0      
      decay1=EXP(-GAMA1*(v1-v2)**2)
      decay2=exp(-D(14)*((r1-3.0d0)**2))
      decay3=exp(-D(15)*((r2+r3-5.2d0)**2))
      decay=decay1*(decay2*decay3)


      DQ1R1=RO3
      DQ1R2=RO3
      DQ1R3=RO3
      DQ2R1=0.0D0
      DQ2R2=RO2
      DQ2R3=-RO2
      DQ3R1=2.0D0*RO6
      DQ3R2=-RO6
      DQ3R3=-RO6
      
      r12=r1*R1
      r22=r2*R2
      r32=r3*R3
      
      COSa=(r22+r32-r12)/(2.0d0*r2*r3)
      if (abs(COSa).ge.1.0d0) then
         sina=0.0d0
       else
         sina=(1.0d0-COSa**2)
      end if
            
      GAMA1=100.0D0      
      THREB12_10=sina*pol*decay
c********************************************************

      call dTHREBR1_10(r1,r2,r3,dr1,dr2,dr3) !superficieh2o.f
      dr1tb1=dr1
      dr2tb1=dr2
      dr3tb1=dr3

      call dRL1_10(r1,r2,r3,dr1,dr2,dr3) !superficieh2o.f
      dr1RL1_10=dr1
      dr2RL1_10=dr2
      dr3RL1_10=dr3


c      v1=VHHS_10(r1)+VOHS_10(r2)+VOHS_10(r3)+THREBR1_10(r1,r2,r3)+
c     1  RL1_10(r1,r2,r3)+VO1D


c     derivadas de v1 em ordem aos rs
      dv1dr1=dVHHS_10(r1) + dr1tb1 + dr1RL1_10
      dv1dr2=dVOHS_10(r2) + dr2tb1 + dr2RL1_10
      dv1dr3=dVOHS_10(r3) + dr3tb1 + dr3RL1_10

c     v2 +++++++

      call dTHREBR2_10(r1,r2,r3,dr1,dr2,dr3) !superficieh2o.f
      dr1tb2=dr1
      dr2tb2=dr2
      dr3tb2=dr3

      call drl2_10(r1,r2,r3,dr1,dr2,dr3) !superficieh2o.f
      dr1rl2=dr1
      dr2rl2=dr2
      dr3rl2=dr3


c      v2=VHHT_10(r1)+VOHP_10(r2)+VOHP_10(r3)+THREBR2_10(r1,r2,r3)+
c     1  rl2_10(r1,r2,r3)

c     derivadas de v2 em ordem aos rs

      dv2dr1=dVHHT_10(r1) + dr1tb2 + dr1rl2
      dv2dr2=dVOHP_10(r2) + dr2tb2 + dr2rl2
      dv2dr3=dVOHP_10(r3) + dr3tb2 + dr3rl2

c   ************************************************************            
      DSINR1=R1*(R22+R32-R12)/(R22*R32)
      DSINR2=-(R22+R32-R12)*(R22-R32+R12)/(2.0D0
     1                                          *R32*R2**3)
      DSINR3=-(R22+R32-R12)*(R32-R22+R12)/(2.0D0
     1                                          *R3**3*R22)
           
      DDECAY1DV1=-2.0D0*GAMA1*(V1-V2)*DECAY1
      DDECAY1DV2=-DDECAY1DV1
      
      ddecay2r1=-D(14)*2.0d0*(r1-3.0d0)*decay2
      ddecay2r2=0.0d0
      ddecay2r3=0.0d0
      
      ddecay3r1=0.0d0
      ddecay3r2=-D(15)*2.0d0*(r2+r3-5.2d0)*decay3
      ddecay3r3=-D(15)*2.0d0*(r2+r3-5.2d0)*decay3
      
      decay1p2=decay1*decay2     
      decay1p3=decay1*decay3
      decay2p3=decay2*decay3
            
      DEXPR1=(DDECAY1DV1*DV1DR1+DDECAY1DV2*DV2DR1)*(decay2p3)+
     1       ddecay2r1*decay1p3+ddecay3r1*decay1p2
      DEXPR2=(DDECAY1DV1*DV1DR2+DDECAY1DV2*DV2DR2)*(decay2p3)+
     1       ddecay2r2*decay1p3+ddecay3r2*decay1p2
      DEXPR3=(DDECAY1DV1*DV1DR3+DDECAY1DV2*DV2DR3)*(decay2p3)+
     1       ddecay2r3*decay1p3+ddecay3r3*decay1p2
      
      dpoq1=D(2)+2.0D0*D(4)*Q1+D(6)*Q3+3.0D0*D(8)*Q12+D(9)*TQ1+
     1      2.0D0*D(11)*Q1*Q3+D(12)*TQ2
      dpoq2=2.0D0*D(5)*Q2+2.0D0*D(7)*Q2+2.0D0*D(9)*Q2*Q1+
     1      2.0D0*D(12)*Q1*Q2-6.0D0*D(10)*Q2*Q3+2.0D0*D(13)*
     2      Q2*Q3
      dpoq3=D(3)+2.0D0*D(5)*Q3+D(6)*Q1-2.0D0*D(7)*Q3+2.0D0*D(9)*Q1*Q3+
     1      3.0D0*D(10)*Q32-3.0D0*D(10)*Q22+D(11)*Q12-2.0D0*D(12)*
     2      Q1*Q3+D(13)*Q22+3.0D0*D(13)*Q32
      
      dpolr1=DPOQ1*DQ1R1+DPOQ3*DQ3R1  
      dpolr2=DPOQ1*DQ1R2+DPOQ2*DQ2R2+DPOQ3*DQ3R2
      dpolr3=DPOQ1*DQ1R3+DPOQ2*DQ2R3+DPOQ3*DQ3R3
      
      POLDECAY=pol*DECAY
      SINaDECAY=SINa*DECAY
      SINapol=SINa*pol
      
      DT12R1=DSINR1*polDECAY+SINaDECAY*dpolr1+DEXPR1*SINapol
      DT12R2=DSINR2*polDECAY+SINaDECAY*dpolr2+DEXPR2*SINapol
      DT12R3=DSINR3*polDECAY+SINaDECAY*dpolr3+DEXPR3*SINapol
 
     
      RETURN
      END

      FUNCTION VOHP_10(R)
c    **************************************************************
C     TO COMPUTE THE HFACE FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c    **************************************************************
      VOHP_10=EHFOHP_10(R)+DISOHP_10(R)
      RETURN
      END

      FUNCTION EHFOHP_10(R)
c    **************************************************************
C     TO COMPUTE THE EHF FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
!$OMP THREADPRIVATE(/DIATDI_10/)


c    **************************************************************
      DATA D,ASV/0.240599442D0,2.4093609D0,1.1587677D0,5.3697044D-1/
      DATA (GAM(I),i=0,2)/1.8450419D0,2.7809326D3,5.3432177D-5/
      DATA AGEX,ALPHEX,A1EX,GAMEX/0.307D0,1.5D0,2.257329D0,2.0D0/
c    **************************************************************
      X=R-RMOHP
      X2=X*X
      X3=X2*X
      RHO=(RMOHP+2.5D0*R0OHP)/2.0D0      
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6
      GAMA=gam(0)*(1+gam(1)*TANH(gam(2)*X))
      EHF=-D/R*(1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3)*EXP(-GAMA*X)
      EX=-AGEX*R**ALPHEX*(1.0D0+A1EX*R)*EXP(-GAMEX*R)*D6
      EHFOHP_10=EHF+EX
      RETURN
      END


      FUNCTION DISOHP_10(R)
c    **************************************************************
C     TO COMPUTE THE DISPERSION FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************
      DISOHP_10=DISPOH_10(R,COHP(6),COHP(8),COHP(10),R0OHP,RMOHP)
      RETURN
      END

      FUNCTION DVOHP_10(R)
c    **************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************      DATA VO1D,EMIN/7.1955D-2,-0.3704003D0/

      DATA D,ASV/0.240599442D0,2.4093609D0,1.1587677D0,5.3697044D-1/
      DATA (GAM(I),i=0,2)/1.8450419D0,2.7809326D3,5.3432177D-5/
      DATA AGEX,ALPHEX,A1EX,GAMEX/0.307D0,1.5D0,2.257329D0,2.0D0/
c    **************************************************************
      X=R-RMOHP
      X2=X*X
      X3=X2*X
      RHO=(RMOHP+2.5D0*R0OHP)/2.0D0      
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6
      GAMA=gam(0)*(1+gam(1)*TANH(gam(2)*X))
      POL=-D/R*(1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3)
      DPOL=-D/R*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*X2)-POL/R
      DGAMA=GAM(0)*GAM(1)*GAM(2)/COSH(GAM(2)*X)**2
      POT=EXP(-GAMA*X)
      DEHF=-(DGAMA*X+GAMA)*POT*POL+DPOL*POT
      RR2=RR*RR
      DBASE=(+A(6)+2.0D0*B(6)*RR)/RHO*
     1        EXP(-A(6)*RR-B(6)*RR2)
      DD6=6.0D0*(1.0D0-EXP(-A(6)*RR-B(6)*RR2))**5*DBASE
      EX=-AGEX*R**ALPHEX*(1.0D0+A1EX*R)*EXP(-GAMEX*R)*D6
      DEXND=(-AGEX*ALPHEX*R**(ALPHEX-1.0D0)-AGEX*A1EX*(ALPHEX+1.0D0)
     1  *R**ALPHEX)*EXP(-GAMEX*R)+GAMEX*AGEX*R**ALPHEX*(1.0D0+A1EX*R)
     2    *EXP(-GAMEX*R)
      DEX=EX/D6*DD6+DEXND*D6
c     Dispersao
      DERDISP=DDISPOH_10(R,COHP(6),COHP(8),COHP(10),R0OHP,RMOHP)
c     Total
      DVOHP_10=DEHF+DEX+DERDISP
      RETURN
      END


      FUNCTION VOHS_10(R)
c    **************************************************************
C     TO COMPUTE THE HFACE FOR  O...H   2SIGMA
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c    **************************************************************
      VOHS_10=EHFOHS_10(R)+DISOHS_10(R)
      RETURN
      END

c     novo VOHS_10 novos coeficientes resultantes dos pontos de Fal61:167
c     e do De de HUB79:0
      
      FUNCTION EHFOHS_10(R)
c    **************************************************************
c     TO COMPUTE THE EHF FOR  O...H   2SIGMA
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************
      DATA D,ASV/0.103283115D0,3.0513954D0,1.8576639D0,1.4517489D0/
      DATA (GAM(I),i=0,2)/2.1847897D0,4.5080859D3,5.3464530D-5/
      DATA AGEX,ALPHEX,A1EX,GAMEX/1.022D0,1.5D0,-0.452055D0,2.0D0/
c    **************************************************************
      X=R-RMOHS
      X2=X*X
      X3=X2*X
      RHO=(rmohs+2.5D0*R0OHS)/2.0d0
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6
      GAMA=gam(0)*(1+gam(1)*TANH(gam(2)*X))
      EHF=-D/R*(1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3)*EXP(-GAMA*X)
      EX=-AGEX*R**ALPHEX*(1.0D0+A1EX*R)*EXP(-GAMEX*R)*D6
      EHFOHS_10=EHF+EX
      RETURN
      END

      FUNCTION DISOHS_10(R)
c    **************************************************************
C     TO COMPUTE THE DISPERSION FOR  O...H   2SIGMA
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************
      DISOHS_10=DISPOH_10(R,COHS(6),COHS(8),COHS(10),R0OHS,RMOHS)
      RETURN
      END

      FUNCTION DVOHS_10(R)
c    *****EXDIS_10*********************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  O...H   2SIGMA
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************
      DATA D,ASV/0.103283115D0,3.0513954D0,1.8576639D0,1.4517489D0/
      DATA (GAM(I),i=0,2)/2.1847897D0,4.5080859D3,5.3464530D-5/
      DATA AGEX,ALPHEX,A1EX,GAMEX/1.022D0,1.5D0,-0.452055D0,2.0D0/
c    **************************************************************
c     Termo sem exchange assimp.   
      X=R-RMOHS
      X2=X*X
      X3=X2*X
      
      POL=-D/R*(1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3)
      DPOL=-D/R*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*X2)-POL/R
      
      GAMA=gam(0)*(1+gam(1)*TANH(gam(2)*X))
      DGAMA=GAM(0)*GAM(1)*GAM(2)/COSH(GAM(2)*X)**2
      POT=EXP(-GAMA*X)
      
      DEHF=-(DGAMA*X+GAMA)*POT*POL+DPOL*POT
      RHO=(rmohs+2.5D0*R0OHS)/2.0d0
      RR=R/RHO
      RR2=RR*RR      
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR2))**6
      DBASE=(+A(6)+2.0D0*B(6)*RR)/RHO*
     1        EXP(-A(6)*RR-B(6)*RR2)
      DD6=6.0D0*(1.0D0-EXP(-A(6)*RR-B(6)*RR2))**5*DBASE
      EX=-AGEX*R**ALPHEX*(1.0D0+A1EX*R)*EXP(-GAMEX*R)*D6
      DEXND=(-AGEX*ALPHEX*R**(ALPHEX-1.0D0)-AGEX*A1EX*(ALPHEX+1.0D0)
     1  *R**ALPHEX)*EXP(-GAMEX*R)+GAMEX*AGEX*R**ALPHEX*(1.0D0+A1EX*R)
     2    *EXP(-GAMEX*R)
      DEX=EX/D6*DD6+DEXND*D6
c     Dispersao
      DERDISP=DDISPOH_10(R,COHS(6),COHS(8),COHS(10),R0OHS,RMOHS)
c     Total
      DVOHS_10=DEHF+DEX+DERDISP
      RETURN
      END

      FUNCTION VHHS_10(R)
c    **************************************************************
C     TO COMPUTE THE HFACE FOR singlet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c    **************************************************************
      VHHS_10=EHFHHS_10(R)+DISHHS_10(R)
      RETURN
      END

      FUNCTION EHFHHS_10(R)
c    **************************************************************
C     TO COMPUTE THE EHF FOR  singlet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CEHFHS_10/X,X2,X3,RHO,RR,GAMA,POL,D6
!$OMP THREADPRIVATE(/DIATDI_10/,/CEHFHS_10/)

c    **************************************************************
      DATA D,ASV/0.218973D0,1.91479d0,0.646041D0,0.346414D0/
      DATA (GAM(I),i=0,2)/1.22349D0,1.04334D0,0.208477D0/
      DATA AGEX,ALPHEX,GAMEX/0.8205D0,2.5D0,2.0D0/
c    **************************************************************
      X=R-RMHHS
      X2=X*X
      X3=X2*X
      POL=(1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3)
      Gama=Gam(0)*(1.0D0+Gam(1)*TANH(Gam(2)*X))      
      EHFHHS_10=-D/R*POL*EXP(-Gama*X)
      RHO=(RMHHS+2.5D0*R0HH)/2.0d0
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6  
      EXCHHS=-AGEX*R**ALPHEX*EXP(-GAMEX*R)*D6
      EHFHHS_10=EHFHHS_10+EXCHHS
      RETURN
      END

      FUNCTION DISHHS_10(R)
c    **************************************************************
C     TO COMPUTE THE DISPERSION FOR singlet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
c    **************************************************************
      R2=R*R
      R6=R2*R2*R2
      R8=R6*R2
      ADISP=-CHH(6)*DAMPHH_10(6,R)/R6-CHH(8)*DAMPHH_10(8,R)/R8
      DO N=10,16
        ADISP=ADISP-CHH(N)*DAMPHH_10(N,R)/(R8*R**(N-8.0d0))
      ENDDO
      DISHHS_10=ADISP 
      RETURN
      END

      FUNCTION DAMPHH_10(N,R)
c    *************************************************************** 
c     CALCULATES VARANDAS-BRANDAO DAMPING FUNCTIONS
c    ***************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
!$OMP THREADPRIVATE(/DIATDI_10/)

c    ***************************************************************
      RR=2.0D0*R/(RMHHS+2.5D0*R0HH)
      DAMPHH_10=(1.0D0-EXP(-A(N)*RR-B(N)*RR*RR))**N
      RETURN
      END

      FUNCTION DVHHS_10(R)
C     *************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  singlet H...H
C     *************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *************************************************************
      DVHHS_10=DEHFHHS_10(R)+DDISHHS_10(R)
      RETURN
      END

      FUNCTION DEHFHHS_10(R)
C     *************************************************************
C     TO COMPUTE THE EHF FOR  singlet H...H
C     *************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CEHFHS_10/X,X2,X3,RHO,RR,GAMA,POL,D6
!$OMP THREADPRIVATE(/DIATDI_10/,/CEHFHS_10/)

c     **************************************************************
      DATA D,ASV/0.218973D0,1.91479d0,0.646041D0,0.346414D0/
      DATA (GAM(I),i=0,2)/1.22349D0,1.04334D0,0.208477D0/
      DATA AGEX,ALPHEX,GAMEX/0.8205D0,2.5D0,2.0D0/
c     **************************************************************
      DGama=Gam(0)*Gam(1)*Gam(2)/COSH(Gam(2)*X)**2
      DPOLC=D/R**2*POL-D/R*(ASV(1)+2*ASV(2)*X+3*ASV(3)*X2)
      DEXP=-(Gama+DGama*X)*EXP(-Gama*X)
      DEHFHHS_10=DPOLC*EXP(-GAMA*X)-D/R*POL*DEXP
      RHH=RHO
      RR2=RR*RR
      DBASE=(A(6)+2*B(6)*RR)/RHH*EXP(-A(6)*RR-    
     1    B(6)*RR2)
      DD6=6*(1-EXP(-A(6)*RR-B(6)*RR2))**5*DBASE
      FEXC=-AGEX*R**ALPHEX*EXP(-GAMEX*R)
      DFEXC=-FEXC*GAMEX+ALPHEX*FEXC/R
      DEXCHHS=DFEXC*D6+DD6*FEXC
      DEHFHHS_10=DEHFHHS_10+DEXCHHS
      RETURN
      END

      FUNCTION DDISHHS_10(R)
c    ************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      R11=R10*R
      R12=R11*R
      R13=R12*R
      R14=R13*R
      R15=R14*R
      R16=R15*R
      RR=2.0D0*R/(RMHHS+2.5D0*R0HH)
      
      rr2=rr*rr
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
         
      DDISP=6.0D0*CHH(6)/R6*T6**5*(T6/R+(1.0D0-T6)*(-A(6)-2.0D0
     1     *B(6)*RR)*DRR)+
     2 8.0D0*CHH(8)/R8*T8**7*(T8/R+(1.0D0-T8)*(-A(8)-2.0D0*
     3      B(8)*RR)*DRR)+
     4 10.0D0*CHH(10)/R10*T10**9*(T10/R+(1.0D0-T10)*(-A(10)-2.0D0*
     5      B(10)*RR)*DRR)+
     611.0D0*CHH(11)/R11*T11**10*(T11/R+(1.0D0-T11)*(-A(11)-2.0D0*
     7  B(11)*RR)*DRR)+
     812.0D0*CHH(12)/R12*T12**11*(T12/R+(1.0D0-T12)*(-A(12)-2.0D0*
     9      B(12)*RR)*DRR)+
     113.0D0*CHH(13)/R13*T13**12*(T13/R+(1.0D0-T13)*(-A(13)-2.0D0*
     2      B(13)*RR)*DRR)+
     314.0D0*CHH(14)/R14*T14**13*(T14/R+(1.0D0-T14)*(-A(14)-2.0D0*
     4      B(14)*RR)*DRR)+
     515.0D0*CHH(15)/R15*T15**14*(T15/R+(1.0D0-T15)*(-A(15)-2.0D0*
     6      B(15)*RR)*DRR)+
     716.0D0*CHH(16)/R16*T16**15*(T16/R+(1.0D0-T16)*(-A(16)-2.0D0*
     8      B(16)*RR)*DRR)
      DDISHHS_10=DDISP
      RETURN
      END

      FUNCTION VHHT_10(R)
c    **************************************************************
C     TO COMPUTE THE HFACE FOR triplet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C    **************************************************************
C     O AJUSTE - VHH(3SI) obteve-se PARTIR DOS PONTOS DE
c                       dois pontos dos nossos para 0.22 e 0.3
C                                    - KOL65:2429 (R=1.0<>5.9)
C                                    - KOL74:457  (R=6.0<>12.0)
c    **************************************************************
      VHHT_10=EHFHHT_10(R)+DISHHT_10(R)
      RETURN
      END

      FUNCTION EHFHHT_10(R)
c    **************************************************************
C     TO COMPUTE THE EHF FOR  triplet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION E(6)
c    **************************************************************
c     coef retirados de var82:857 pag 862
      DATA D,E/2.64103167D-05,-1.4993795d+00,2.9009680d+00,
     1       4.7898562d-01,2.7407202d-02,-9.0948301d-01,
     2       -2.0122436d-01/
      DATA AGEX,ALPHEX,GAMEX/-0.8205D0,2.5D0,2.0D0/
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/CEHFHT_10/R1,R2,R3,RHO,RR,EX,D6,TEXP,POL,EHF
!$OMP THREADPRIVATE(/DIATDI_10/,/CEHFHT_10/)

c    **************************************************************
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
      EHFHHT_10=EHF+EX
      RETURN
      END

      FUNCTION DISHHT_10(R)
c    **************************************************************
C     TO COMPUTE THE DISPERSION FOR triplet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************
      R6=R**6
      R8=R6*R*R
      ADISP=-CHH(6)*DAMPHHT_10(6,R)/R6-CHH(8)*DAMPHHT_10(8,R)/R8
      DO N=10,16
           ADISP=ADISP-CHH(N)*DAMPHHT_10(N,R)/(R8*R**(N-8.0D0))
      ENDDO
      DISHHT_10=ADISP
      RETURN
      END
      
      FUNCTION DAMPHHT_10(N,R)
c    *************************************************************** 
c     CALCULATES VARANDAS-BRANDAO DAMPING FUNCTIONS
c    ***************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
!$OMP THREADPRIVATE(/DIATDI_10/)

c    ***************************************************************
      RR=2.0D0*R/(RMHHT+2.5D0*R0HH)
      DAMPHHT_10=(1.0D0-EXP(-A(N)*RR-B(N)*RR*RR))**N
      RETURN
      END

      FUNCTION DVHHT_10(R)
c    **************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  triplet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c    **************************************************************
      DVHHT_10=DEHFHHT_10(R)+DDISHHT_10(R)
      RETURN
      END

      FUNCTION DEHFHHT_10(R)
c    **************************************************************
C     TO COMPUTE THE EHF FOR  triplet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION E(6)
c    **************************************************************
      DATA D,E/2.64103167D-05,-1.4993795d+00,2.9009680d+00,
     1       4.7898562d-01,2.7407202d-02,-9.0948301d-01,
     2       -2.0122436d-01/
      DATA AGEX,ALPHEX,GAMEX/-0.8205D0,2.5D0,2.0D0/
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CEHFHT_10/R1,R2,R3,RHO,RR,EX,D6,TEXP,POL,EHF
!$OMP THREADPRIVATE(/DIATDI_10/,/CEHFHT_10/)

c    **************************************************************
      rr2=rr*rr
      DBASE=(A(6)+2.0d0*B(6)*RR)/RHO*EXP(-A(6)*RR-B(6)*RR2)
      DD6=6.0d0*DBASE*(1.0d0-EXP(-A(6)*RR-B(6)*RR2))**5
      FEXC=EX/D6
      DFEXC=-agex*alphex*r**(ALPHEX-1.0d0)*EXP(-GAMEX*R)-GAMEX*
     1       FEXC
      DEXCHHT=DFEXC*D6+DD6*FEXC
      DEHFHHT_10=D/R**2*POL*Texp+E(1)*D/R*Texp+2.0d0*E(5)*r1*D/r*Texp+
     1       3.0d0*E(6)*r2*D/r*Texp-D/R*POL*(-E(2)-2.0D0*
     1        E(3)*R1-3.0d0*E(4)*R2)*Texp      
      DEHFHHT_10=DEHFHHT_10+DEXCHHT
      RETURN
      END      

      FUNCTION DDISHHT_10(R)
c    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************
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
         
      DDISP=6.0D0*CHH(6)/R6*T6**5*(T6/R+(1.0D0-T6)*(-A(6)-2.0D0
     1     *B(6)*RR)*DRR)+
     2 8.0D0*CHH(8)/R8*T8**7*(T8/R+(1.0D0-T8)*(-A(8)-2.0D0*
     3      B(8)*RR)*DRR)+
     4 10.0D0*CHH(10)/R10*T10**9*(T10/R+(1.0D0-T10)*(-A(10)-2.0D0*
     5      B(10)*RR)*DRR)+
     611.0D0*CHH(11)/R11*T11**10*(T11/R+(1.0D0-T11)*(-A(11)-2.0D0*
     7  B(11)*RR)*DRR)+
     812.0D0*CHH(12)/R12*T12**11*(T12/R+(1.0D0-T12)*(-A(12)-2.0D0*
     9      B(12)*RR)*DRR)+
     113.0D0*CHH(13)/R13*T13**12*(T13/R+(1.0D0-T13)*(-A(13)-2.0D0*
     2      B(13)*RR)*DRR)+
     314.0D0*CHH(14)/R14*T14**13*(T14/R+(1.0D0-T14)*(-A(14)-2.0D0*
     4      B(14)*RR)*DRR)+
     515.0D0*CHH(15)/R15*T15**14*(T15/R+(1.0D0-T15)*(-A(15)-2.0D0*
     6      B(15)*RR)*DRR)+
     716.0D0*CHH(16)/R16*T16**15*(T16/R+(1.0D0-T16)*(-A(16)-2.0D0*
     8      B(16)*RR)*DRR)
      DDISHHT_10=DDISP
      RETURN
      END

      FUNCTION THREBR1_10(R1,R2,R3)
c    **************************************************************
c
C     TO COMPUTE THE THREE BODY TERM IN COORDINATES R1,R2,R3
c
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(nc1=95,nc2=34,nc3=15,nc=144)
      COMMON/COEFF_10/C(148)
      COMMON/REFGEO_10/R10,R20,R30
      COMMON/REFGEO2_10/R102,R202,R302
      COMMON/REFGEO12_10/R1012,R2012,R3012
      COMMON/CTHRB1_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q14,Q15,
     1 Q16,Q17,Q18,Q22,Q23,Q32,Q33,TQ1,TQ12,TQ13,TQ14,TQ2,TQ3,TQ32,
     2 POLQ,DECAY1,DECAY2,DECAY3
!$OMP THREADPRIVATE(/CTHRB1_10/)

c    **************************************************************
      RO2=1.0D0/SQRT(2.0D0)
      RO3=1.0D0/SQRT(3.0D0)
      RO6=1.0D0/SQRT(6.0D0)

      S1=R1-R10
      S2=R2-R20
      S3=R3-R30
      Q1=RO3*(s1+s2+s3)
      Q2=RO2*(s2-s3)
      Q3=RO6*(2.0D0*s1-s2-s3)
      Q12=Q1*Q1
      Q13=Q12*Q1
      Q14=Q13*Q1
      Q15=Q14*Q1
      Q16=Q15*Q1
      Q17=Q16*Q1
      Q18=Q17*Q1
      Q22=Q2*Q2
      Q23=Q22*Q2
      Q32=Q3*Q3
      Q33=Q32*Q3
      TQ1=Q22+Q32
      TQ12=TQ1*TQ1     
      TQ13=TQ12*TQ1     
      TQ14=TQ13*TQ1
      TQ2=Q22-Q32
      TQ3=(Q32-3.0D0*Q22)*Q3
      TQ32=TQ3*TQ3     
      
      POLQ=C(1)+C(2)*Q1+C(3)*Q3+C(4)*Q12+
     1   C(5)*TQ1+C(6)*Q1*Q3+C(7)*TQ2+C(8)*Q13+C(9)*Q1*TQ1+
     2   C(10)*TQ3+C(11)*Q12*Q3+C(12)*Q1*TQ2+C(13)*Q3*TQ1+
     3   C(14)*Q14+C(15)*Q12*TQ1+C(16)*TQ1**2+C(17)*Q1*TQ3+
     4   C(18)*Q13*Q3+C(19)*Q12*TQ2+C(20)*Q1*Q3*TQ1+        
     5   C(21)*Q3*TQ3+C(22)*TQ1*TQ2+C(23)*Q15+C(24)*Q13*TQ1+
     6 C(25)*Q1*TQ1**2+C(26)*Q12*TQ3+C(27)*TQ1*TQ3+C(28)*Q14*Q3+
     7 C(29)*Q13*TQ2+C(30)*Q12*Q3*TQ1+C(31)*Q1*Q3*TQ3+C(32)*Q1*TQ1*TQ2
     8 +C(33)*Q3*TQ1**2+C(34)*TQ2*TQ3+
     9 C(35)*Q16+C(36)*Q14*TQ1+C(37)*Q12*TQ12+C(38)*Q13*TQ3+C(39)*Q1*
     1 TQ3*TQ1+C(40)*TQ13+C(41)*TQ32+C(42)*Q15*Q3+C(43)*Q14*TQ2+C(44)*
     2 Q13*Q3*TQ1+C(45)*Q12*Q3*TQ3+C(46)*Q12*TQ1*TQ2+C(47)*Q1*Q3*TQ12+
     3 C(48)*Q1*TQ3*TQ2+C(49)*Q3*TQ1*TQ3+C(50)*TQ2*TQ12
     4 +C(51)*Q17+C(52)*Q15*TQ1+C(53)*Q13*TQ12+C(54)*Q14*TQ3+C(55)*Q12*
     5 TQ1*TQ3+C(56)*Q1*TQ13+C(57)*Q1*TQ32+C(58)*TQ3*TQ12+C(59)*Q3*Q16+
     6 C(60)*Q3*Q14*TQ1+C(61)*Q3*Q12*TQ12+C(62)*Q13*Q3*TQ3+C(63)*Q1*Q3*
     7 TQ3*TQ1+C(64)*Q3*TQ13+C(65)*Q3*TQ32+C(66)*Q15*TQ2+C(67)*Q13*TQ1*
     8 TQ2+C(68)*Q12*TQ3*TQ2+C(69)*Q1*TQ12*TQ2+C(70)*TQ3*TQ1*TQ2+
     9 c(71)*Q18+C(72)*Q16*TQ1+C(73)*Q14*TQ12+C(74)*Q15*TQ3+C(75)*Q13*
     1 TQ3*TQ1+C(76)*Q12*TQ13+C(77)*Q12*TQ32+C(78)*Q1*TQ3*TQ12+C(79)*
     2 TQ14+C(80)*TQ32*TQ1+C(81)*Q17*Q3+C(82)*Q15*Q3*TQ1+C(83)*Q13*Q3*
     3 TQ12+C(84)*Q14*TQ32+C(85)*Q12*TQ32*TQ1+C(86)*Q1*Q3*TQ13+C(87)*
     4 Q1*Q3*TQ32+C(88)*Q16*TQ2+C(89)*Q14*TQ2*TQ1+C(90)*Q13*TQ3*TQ2+
     5 C(91)*Q12*TQ12*TQ2+C(92)*Q1*TQ3*TQ2*TQ1+C(93)*Q3*TQ3*TQ12+C(94)
     6 *TQ32*TQ2+C(95)*TQ13*TQ2
     
      DECAY1=1.0D0-TANH(S1*C(145))
      DECAY2=1.0D0-TANH(S2*C(146))
      DECAY3=1.0D0-TANH(S3*C(146))
      THREBR1_10=POLQ*DECAY1*DECAY2*DECAY3
      RETURN
      END

      FUNCTION THREBR2_10(R1,R2,R3)
c    **************************************************************
c
C     TO COMPUTE THE THREE BODY TERM IN COORDINATES R1,R2,R3
c
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(nc1=95,nc2=34,nc3=15,nc=144)
      DIMENSION D(nc2+19)
      COMMON/COEFF_10/C(148)
      COMMON/REFGEO_10/R10,R20,R30
      COMMON/REFGEO2_10/R102,R202,R302
      COMMON/REFGEO12_10/R1012,R2012,R3012
      COMMON/CTHRB2_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q14,Q15,
     1 Q16,Q17,Q18,Q22,Q23,Q32,Q33,TQ1,TQ12,TQ13,TQ14,TQ2,TQ3,TQ32,
     2 POLQ,DECAY1,DECAY2,DECAY3
!$OMP THREADPRIVATE(/CTHRB2_10/)


      EQUIVALENCE (C(NC1+1),D(1))
c    *************************************************************
      RO2=1.0D0/SQRT(2.0D0)
      RO3=1.0D0/SQRT(3.0D0)
      RO6=1.0D0/SQRT(6.0D0)
      S1=R1-R102
      S2=R2-R202
      S3=R3-R302
      Q1=RO3*(s1+s2+s3)
      Q2=RO2*(s2-s3)
      Q3=RO6*(2.0D0*s1-s2-s3)
      Q12=Q1*Q1
      Q13=Q12*Q1
      Q14=Q13*Q1
      Q15=Q14*Q1
      Q16=Q15*Q1
      Q17=Q16*Q1
      Q18=Q17*Q1
      Q22=Q2*Q2
      Q23=Q22*Q2
      Q32=Q3*Q3
      Q33=Q32*Q3
      TQ1=Q22+Q32
      TQ12=TQ1*TQ1     
      TQ13=TQ12*TQ1     
      TQ14=TQ13*TQ1
      TQ2=Q22-Q32
      TQ3=(Q32-3.0D0*Q22)*Q3
      TQ32=TQ3*TQ3     
      
      POLQ=D(1)+D(2)*Q1+D(3)*Q3+D(4)*Q12+
     1   D(5)*TQ1+D(6)*Q1*Q3+D(7)*TQ2+D(8)*Q13+D(9)*Q1*TQ1+
     2   D(10)*TQ3+D(11)*Q12*Q3+D(12)*Q1*TQ2+D(13)*Q3*TQ1+
     3   D(14)*Q14+D(15)*Q12*TQ1+D(16)*TQ1**2+D(17)*Q1*TQ3+
     4   D(18)*Q13*Q3+D(19)*Q12*TQ2+D(20)*Q1*Q3*TQ1+        
     5   D(21)*Q3*TQ3+D(22)*TQ1*TQ2+D(23)*Q15+D(24)*Q13*TQ1+
     6 D(25)*Q1*TQ1**2+D(26)*Q12*TQ3+D(27)*TQ1*TQ3+D(28)*Q14*Q3+
     7 D(29)*Q13*TQ2+D(30)*Q12*Q3*TQ1+D(31)*Q1*Q3*TQ3+D(32)*Q1*TQ1*TQ2
     8 +D(33)*Q3*TQ1**2+D(34)*TQ2*TQ3
     
      DECAY1=1.0D0-TANH(S1*D(52))
      DECAY2=1.0D0-TANH(S2*D(53))
      DECAY3=1.0D0-TANH(S3*D(53))
      THREBR2_10=POLQ*DECAY1*DECAY2*DECAY3
      RETURN
      END

      FUNCTION THREB12_10(R1,R2,R3)
c    **************************************************************
C     TO COMPUTE THE THREE BODY TERM IN COORDINATES R1,R2,R3
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter(nc1=95,nc2=34,nc3=15,nc=144)
      DIMENSION d(NC3)
C      DIMENSION X(3),G(3)
C      DATA VO1D,EMIN/7.1955D-2,-0.3704003D0/
      COMMON/COEFF_10/C(148)
c      EQUIVALENCE (C(NC1+NC2+1),d(1))
      COMMON/REFGEO12_10/R1012,R2012,R3012
!      common/escrev_10/iwrit
      COMMON/CTHRB12_10/RO2,RO3,RO6,S1,S2,S3,Q1,Q2,Q3,Q12,Q13,Q22,Q32,
     1 TQ1,TQ2,TQ3,POL,DECAY1,DECAY2,DECAY3,DECAY,DV1DR1,
     2 DV1DR2,DV1DR3,DV2DR1,DV2DR2,DV2DR3,V1,V2
!$OMP THREADPRIVATE(/CTHRB12_10/)

c    **************************************************************

      v1=VHHS_10(r1)+VOHS_10(r2)+VOHS_10(r3)+THREBR1_10(r1,r2,r3)+
     1  RL1_10(r1,r2,r3)+VO1D
      v2=VHHT_10(r1)+VOHP_10(r2)+VOHP_10(r3)+THREBR2_10(r1,r2,r3)+
     1  rl2_10(r1,r2,r3)


      do i=1,nc3
          d(i)=c(nc1+nc2+i)
      end do

      term=(r2**2+r3**2-r1**2)/(2.0D0*r2*r3)
      
      if (abs(term).ge.1.0d0) then
         sina=0.0d0
       else
         sina=(1.0d0-term**2)
      end if
      

      RO2=1.0D0/SQRT(2.0D0)
      RO3=1.0D0/SQRT(3.0D0)
      RO6=1.0D0/SQRT(6.0D0)

      S1=R1-r1012
      S2=R2-r2012
      S3=R3-r3012
      Q1=RO3*(s1+s2+s3)
      Q2=RO2*(s2-s3)
      Q3=RO6*(2.0D0*s1-s2-s3)
      Q12=Q1*Q1
      Q13=Q12*Q1
      Q22=Q2*Q2
      Q32=Q3*Q3
      TQ1=Q22+Q32
      TQ2=Q22-Q32
      TQ3=(Q32-3.0D0*Q22)*Q3
      
      pol=D(1)+D(2)*Q1+D(3)*Q3+D(4)*Q12+
     1   D(5)*TQ1+D(6)*Q1*Q3+D(7)*TQ2+D(8)*Q13+D(9)*Q1*TQ1+
     2   D(10)*TQ3+D(11)*Q12*Q3+D(12)*Q1*TQ2+D(13)*Q3*TQ1
C      
C     ESTES DVS SAO POSTOS NO COMMON PARA SEREM USADOS NA DTHRB12
C
      GAMA1=100.0D0      
      decay1=EXP(-GAMA1*(v1-v2)**2)
      decay2=exp(-D(14)*((r1-3.0d0)**2))
      decay3=exp(-D(15)*((r2+r3-5.2d0)**2))
      decay=decay1*(decay2*decay3)
      THREB12_10=sina*pol*decay
C      PRINT*,'ESTOU NA TRH12'
C      PRINT*,'decay1,decay2,decay3',decay1,decay2,decay3
C      PRINT*,'decay,pol,sina,threb12',decay,pol,sina,threb12
C      PRINT*,'GAMA1,v1,v2,r1,r2,r3,D(14),D(15)'
C      PRINT*,GAMA1,v1,v2,r1,r2,r3,D(14),D(15)


      RETURN
      END

      SUBROUTINE GEOMRL_10(R1,R2,R3)  
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONST_10/PI
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
!$OMP THREADPRIVATE(/DIATDI_10/,/CVED_10/,/CVED31_10/,/CVED32_10/)

c     ******************************************************************     
      CONST=4.5d0
      
      r12=r1**2
      r22=r2**2
      r32=r3**2    
      
      COSB=(r12+r22-r32)/(2.0D0*r2*r1)
      if (cosb.gt.1.0d0) cosb=1.0d0
      if (cosb.lt.-1.0d0) cosb=-1.0d0
      
      COSC=(r12+r32-r22)/(2.0D0*r3*r1)
      if (cosc.gt.1.0d0) cosc=1.0d0
      if (cosc.lt.-1.0d0) cosc=-1.0d0
      
      rg1=SQRT((0.5D0*r1)**2+r22-r1*r2*COSB) 
      rg2=SQRT((0.5D0*r2)**2+r12-r2*r1*COSB) 
      rg3=SQRT((0.5D0*r3)**2+r12-r3*r1*COSC) 
      
      RG14=RG1**4      
      rg15=rg1*RG14
      
      DRG1R1=0.5D0/RG1*(2.0D0*0.5D0**2*R1-R2*COSB-
     1                            (R12-R22+R32)/(2.0D0*R1))
      DRG2R1=0.5D0/RG2*(2.0D0*R1-R2*COSB-
     1                            (R12-R22+R32)/(2.0D0*R1))
      DRG3R1=0.5D0/RG3*(2.0D0*R1-R3*COSC-
     1                            (R12-R32+R22)/(2.0D0*R1))
     
      DRG1R2=0.5D0/RG1*(2.0D0*R2-R1*COSB-(R22-R12+R32)
     1                                              /(2.0D0*R2))
      DRG2R2=0.5D0/RG2*(0.5D0**2*2.0D0*R2-R1*COSB-(R22-R12+R32)
     1                                              /(2.0D0*R2))
      DRG3R2=0.5D0/RG3*R2
      
      DRG1R3=0.5D0*r3/RG1
      DRG2R3=0.5D0*r3/RG2
      DRG3R3=0.5D0/RG3*(2.0D0*0.5D0**2*R3-R1*COSC-(R32-R12+
     1                                       R22)/(2.0D0*R3))
     
      cost1=(R3-R2)/R1
      if (cost1.gt.1.0d0) cost1=1.0d0
      if (cost1.lt.-1.0d0) cost1=-1.0d0
      DCOST1R1=-(R3-R2)/R12
      DCOST1R2=-1.0D0/R1
      DCOST1R3=1.0D0/R1
      
      cost2=(R3-R1)/R2
      if (cost2.gt.1.0d0) cost2=1.0d0
      if (cost2.lt.-1.0d0) cost2=-1.0d0
      DCOST2R1=-1.0D0/R2
      DCOST2R2=-(R3-R1)/R22
      DCOST2R3=1.0D0/R2
      
      cost3=(R2-R1)/R3
      if (cost3.gt.1.0d0) cost3=1.0d0
      if (cost3.lt.-1.0d0) cost3=-1.0d0
      DCOST3R1=-1.0D0/R3     
      DCOST3R2=1.0D0/R3
      DCOST1R3=1.0D0/R1
      DCOST3R3=-(R2-R1)/R32
      
C     **********
C     CONSTANTES ADICIONAIS USADAS NA VEOHHS, VEOHHT, VINDHOHS E VINDHOHP
      cost12=cost1*cost1
      cost13=cost12*cost1
      cost14=cost12*cost12
      cost15=cost14*cost1
      cost16=cost14*cost12
      cos2t1=2.0d0*cost12-1.0d0
      cos4t1=8.0d0*cost14-8.0d0*cost12+1.0d0
      cos6t1=32.0d0*cost16-48.0d0*cost14+18.0d0*cost12-1.0d0
      dcos2t=4.0d0*cost1
      dcos4t=32.0d0*cost13-16.0d0*cost1
      dcos6t=192.0d0*cost15-192.0d0*cost13+36.0d0*cost1
      DAMPOHHS5=DAMPOHHS_10(RG1,5,R1)/RG15
      DAMPOHHT5=DAMPOHHT_10(Rg1,5,r1)/rg15
C     *******************
 
      umsabc=(0.5d0+0.5d0*tanh(CONST*((rg1/r1)-2.0d0)))
      sabc=1.0d0-umsabc**2

      umscab=(0.5d0+0.5d0*tanh(CONST*((rg2/r2)-2.0d0)))
      scab=1.0d0-umscab**2

      umsbac=(0.5d0+0.5d0*tanh(CONST*((rg3/r3)-2.0d0)))
      sbac=1.0d0-umsbac**2

      FAB=sabc*sbac
      FAC=sabc*scab
      FbC=sbac*scab
      FABC=umsabc
      FCAB=umscab
      FBAC=umsbac
      
      
      DFABCR1=0.5D0*(CONST*DRG1R1*R1-CONST*RG1)/R12*(1.0D0/COSH
     1        (CONST*(RG1/R1-2.0D0)))**2
      DFCABR1=0.5D0*CONST*DRG2R1/R2*(1.0D0/COSH(CONST*(
     1        RG2/R2-2.0D0)))**2
      DFBACR1=0.5D0*CONST*DRG3R1/R3*(1.0D0/COSH(CONST*(
     1        RG3/R3-2.0D0)))**2

      DFABCR2=0.5d0*(1.0d0/cosH(CONST*(RG1/R1-2.0D0)))**2*
     1        CONST*DRG1R2/R1
      DFCABR2=0.5d0*(1.0d0/cosH(CONST*(RG2/R2-2.0D0)))**2*
     1        (CONST*DRG2R2*R2-CONST*RG2)/R22
      DFBACR2=0.5d0*(1.0d0/cosH(CONST*(RG3/R3-2.0D0)))**2*
     1        CONST*DRG3R2/R3

      DFABCR3=0.5D0*(1.0D0/COSH(CONST*(RG1/R1-2.0D0)))**2*
     1        CONST*DRG1R3/R1
      DFCABR3=0.5D0*(1.0D0/COSH(CONST*(RG2/R2-2.0D0)))**2*
     1        CONST*DRG2R3/R2
      DFBACR3=0.5D0*(1.0D0/COSH(CONST*(RG3/R3-2.0D0)))**2*
     1        (CONST*DRG3R3*R3-CONST*RG3)/R32
      
      DSABCR1=-2.0D0*UMSABC*DFABCR1
      DSBACR1=-2.0D0*UMSBAC*DFBACR1
      DSCABR1=-2.0D0*UMSCAB*DFCABR1
      
      DSABCR2=-2.0D0*UMSABC*DFABCR2
      DSBACR2=-2.0D0*UMSBAC*DFBACR2
      DSCABR2=-2.0D0*UMSCAB*DFCABR2
      
      DSABCR3=-2.0D0*UMSABC*DFABCR3
      DSBACR3=-2.0D0*UMSBAC*DFBACR3
      DSCABR3=-2.0D0*UMSCAB*DFCABR3
      
      
      DFABR1=DSABCR1*SBAC+SABC*DSBACR1
      DFACR1=DSABCR1*SCAB+SABC*DSCABR1
      DFBCR1=DSBACR1*SCAB+SBAC*DSCABR1
      
      DFABR2=DSABCR2*SBAC+SABC*DSBACR2
      DFACR2=DSABCR2*SCAB+SABC*DSCABR2
      DFBCR2=DSBACR2*SCAB+SBAC*DSCABR2
      
      DFABR3=DSABCR3*SBAC+SABC*DSBACR3
      DFACR3=DSABCR3*SCAB+SABC*DSCABR3
      DFBCR3=DSBACR3*SCAB+SBAC*DSCABR3


      C6BC=CHH(6)
      C8BC=CHH(8)
      C10BC=CHH(10)
       
      C6AC=COHS(6)
      C8AC=COHS(8)
      C10AC=COHS(10)
      C6AB=COHS(6)
      C8AB=COHS(8)
      C10AB=COHS(10)
      
      C6AC2=COHP(6)
      C8AC2=COHP(8)
      C10AC2=COHP(10)
      C6AB2=COHP(6)
      C8AB2=COHP(8)
      C10AB2=COHP(10)
      
c calculado a partir das dist. de eq.(POL94:7651) RMHHS=2.8619327675d0
c calculado a partir das dist. de eq.(POL94:7651) RMOHS=1.80965D0
      
      C6OHHT11=C6OHHT_10(r1,cost1)
      C8OHHT11=C8OHHT_10(r1,cost1)
      C10OHHT11=C10OHHT_10(r1,cost1)
      
      C6HOHP12=C6HOHP_10(r2,cost2)
      C8HOHP12=C8HOHP_10(r2,cost2)
      C10HOHP12=C10HOHP_10(r2,cost2)
      
      C6HOHP13=C6HOHP_10(r3,cost3)
      C8HOHP13=C8HOHP_10(r3,cost3)
      C10HOHP13=C10HOHP_10(r3,cost3)
      
      DAMPOHHT16=DAMPOHHT_10(RG1,6,r1)
      DAMPOHHT18=DAMPOHHT_10(RG1,8,r1)
      DAMPOHHT110=DAMPOHHT_10(RG1,10,r1)
      
      DAMPHOHP26=DAMPHOHP_10(RG2,6,r2)
      DAMPHOHP28=DAMPHOHP_10(RG2,8,r2)
      DAMPHOHP210=DAMPHOHP_10(RG2,10,r2)
      
      DAMPHOHP36=DAMPHOHP_10(RG3,6,r3)
      DAMPHOHP38=DAMPHOHP_10(RG3,8,r3)
      DAMPHOHP310=DAMPHOHP_10(RG3,10,r3)
      
      C6OHHS11=C6OHHS_10(r1,cost1)
      C8OHHS11=C8OHHS_10(r1,cost1)
      C10OHHS11=C10OHHS_10(r1,cost1)
      
      C6HOHS12=C6HOHS_10(r2,cost2)
      C8HOHS12=C8HOHS_10(r2,cost2)
      C10HOHS12=C10HOHS_10(r2,cost2)
      
      C6HOHS13=C6HOHS_10(r3,cost3)
      C8HOHS13=C8HOHS_10(r3,cost3)
      C10HOHS13=C10HOHS_10(r3,cost3)
      
      RG16=RG15*RG1
      RG17=RG16*RG1
      RG18=RG17*RG1
      RG19=RG18*RG1
      RG110=RG19*RG1
      RG112=RG16*RG16
      RG116=RG18*RG18
      RG120=RG110*RG110
      
      RG25=RG2**5      
      RG26=RG25*RG2
      RG27=RG26*RG2
      RG28=RG27*RG2
      RG29=RG28*RG2
      RG210=RG29*RG2
      RG212=RG26*RG26
      RG216=RG28*RG28
      RG220=RG210*RG210
      
      RG35=RG3**5      
      RG36=RG35*RG3
      RG37=RG36*RG3
      RG38=RG37*RG3
      RG39=RG38*RG3
      RG310=RG39*RG3
      RG312=RG36*RG36
      RG316=RG38*RG38
      RG320=RG310*RG310
            
      DAMPOHHS16=DAMPOHHS_10(RG1,6,r1)
      DAMPOHHS18=DAMPOHHS_10(RG1,8,r1)
      DAMPOHHS110=DAMPOHHS_10(RG1,10,r1)
      
      DAMPHOHS26=DAMPHOHS_10(RG2,6,r2)
      DAMPHOHS28=DAMPHOHS_10(RG2,8,r2)
      DAMPHOHS210=DAMPHOHS_10(RG2,10,r2)
      
      DAMPHOHS36=DAMPHOHS_10(RG3,6,r3)
      DAMPHOHS38=DAMPHOHS_10(RG3,8,r3)
      DAMPHOHS310=DAMPHOHS_10(RG3,10,r3)
      
      RETURN
      END

      

      FUNCTION VED31_10(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB

!$OMP THREADPRIVATE(/CVED_10/,/CVED31_10/)
c     ******************************************************************
      VED3=
     1     -C6OHHS11*FABC/rg16*DAMPOHHS16-
     1      C8OHHS11*FABC/rg18*DAMPOHHS18-
     2      C10OHHS11*FABC/rg110*DAMPOHHS110
     3     -C6HOHS12*FCAB/rg26*DAMPHOHS26-
     4      C8HOHS12*FCAB/rg28*DAMPHOHS28-
     5      C10HOHS12*FCAB/rg210*DAMPHOHS210
     6     -C6HOHS13*FBAC/rg36*DAMPHOHS36-
     7      C8HOHS13*FBAC/rg38*DAMPHOHS38-
     8      C10HOHS13*FBAC/rg310*DAMPHOHS310
     
      VED2=(FBC-1.0D0)*DISHHS_10(R1)+(FAB-1.0D0)*DISOHS_10(R2)+
     1                            (FAC-1.0D0)*DISOHS_10(R3)      
      
      VED31_10=VED3+VED2
      
      RETURN
      END

      FUNCTION DVED31R1(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/,/CVED_10/,/CVED31_10/)

c     ******************************************************************
      DVED31R1=-((DC6OHHSR1_10(r1,cost1,DCOST1R1)*FABC+DFABCR1*
     1      C6OHHS11)*
     1      DAMPOHHS16/rg16+(DDAMPOHHSR1_10(RG1,6,R1,DRG1R1)*
     1      RG16-6.0D0*
     1      DRG1R1*RG15*DAMPOHHS16)/RG112*C6OHHS11*FABC)-
     1      ((DC8OHHSR1_10(r1,cost1,DCOST1R1)*FABC+DFABCR1*C8OHHS11)*
     1      DAMPOHHS18/rg18+(DDAMPOHHSR1_10(RG1,8,r1,DRG1R1)*
     1      RG18-8.0D0*
     1      DRG1R1*RG17*DAMPOHHS18)/RG116*C8OHHS11*FABC)-
     1      ((DC10OHHSR1_10(r1,cost1,DCOST1R1)*FABC+DFABCR1*C10OHHS11)*
     1      DAMPOHHS110/rg110+(DDAMPOHHSR1_10(RG1,10,R1,DRG1R1)*RG110-
     1      10.0D0*DRG1R1*RG19*DAMPOHHS110)/RG120*C10OHHS11*FABC)-
     1      ((DC6HOHSR1_10(r2,cost2,DCOST2R1)*FCAB+DFCABR1*C6HOHS12)*
     1      DAMPHOHS26/rg26+(DDAMPHOHSR1_10(RG2,6,R2,DRG2R1)*RG26-6.0D0*
     1      DRG2R1*RG25*DAMPHOHS26)/RG212*C6HOHS12*FCAB)-
     1      ((DC8HOHSR1_10(r2,cost2,DCOST2R1)*FCAB+DFCABR1*C8HOHS12)*
     1      DAMPHOHS28/rg28+(DDAMPHOHSR1_10(RG2,8,R2,DRG2R1)*RG28-8.0D0*
     1      DRG2R1*RG27*DAMPHOHS28)/RG216*C8HOHS12*FCAB)-
     1      ((DC10HOHSR1_10(r2,cost2,DCOST2R1)*FCAB+DFCABR1*C10HOHS12)*
     1      DAMPHOHS210/rg210+(DDAMPHOHSR1_10(RG2,10,R2,DRG2R1)*RG210-
     1      10.0D0*DRG2R1*RG29*DAMPHOHS210)/RG220*C10HOHS12*FCAB)-
     1      ((DC6HOHSR1_10(r3,cost3,DCOST3R1)*FBAC+DFBACR1*C6HOHS13)*
     1      DAMPHOHS36/rg36+(DDAMPHOHSR1_10(RG3,6,R3,DRG3R1)*RG36-
     1      6.0D0*DRG3R1*RG35*DAMPHOHS36)/RG312*C6HOHS13*FBAC)-
     1      ((DC8HOHSR1_10(r3,cost3,DCOST3R1)*FBAC+DFBACR1*C8HOHS13)*
     1      DAMPHOHS38/rg38+(DDAMPHOHSR1_10(RG3,8,R3,DRG3R1)*RG38-
     1      8.0D0*DRG3R1*RG37*DAMPHOHS38)/RG316*C8HOHS13*FBAC)-
     1      ((DC10HOHSR1_10(r3,cost3,DCOST3R1)*FBAC+DFBACR1*C10HOHS13)*
     1      DAMPHOHS310/rg310+(DDAMPHOHSR1_10(RG3,10,R3,DRG3R1)*RG310-
     1    10.0D0*DRG3R1*RG39*DAMPHOHS310)/RG320*C10HOHS13*FBAC)
     
CC      VED2=(FBC-1.0D0)*DISHHS_10(R1)+(FAB-1.0D0)*DISOHS_10(R2)+
CC     1                            (FAC-1.0D0)*DISOHS_10(R3)

      DVED21R1=DFBCR1*DISHHS_10(R1)+(FBC-1.0D0)*DDISHHS_10(R1)+
     1        DFABR1*DISOHS_10(R2)+
     2        DFACR1*DISOHS_10(R3)
     
      DVED31R1=DVED31R1+DVED21R1
      
      RETURN
      END

      FUNCTION DVED31R2_10(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/,/CVED_10/,/CVED31_10/)

C     ******************************************************************
      DVED31R2_10=-((DC6OHHSR2_10(r1,cost1,DCOST1R2)*FABC+DFABCR2*
     1      C6OHHS11)*
     1      DAMPOHHS16/rg16+(DDAMPOHHSR2_10(RG1,6,R1,DRG1R2)*RG16-6.0D0*
     1      DRG1R2*RG15*DAMPOHHS16)/RG112*C6OHHS11*FABC)-
     1      ((DC8OHHSR2_10(r1,cost1,DCOST1R2)*FABC+DFABCR2*C8OHHS11)*
     1      DAMPOHHS18/rg18+(DDAMPOHHSR2_10(RG1,8,R1,DRG1R2)*RG18-8.0D0*
     1      DRG1R2*RG17*DAMPOHHS18)/RG116*C8OHHS11*FABC)-
     1      ((DC10OHHSR2_10(r1,cost1,DCOST1R2)*FABC+DFABCR2*C10OHHS11)*
     1      DAMPOHHS110/rg110+(DDAMPOHHSR2_10(RG1,10,R1,DRG1R2)*RG110-
     1      10.0D0*DRG1R2*RG19*DAMPOHHS110)/RG120*C10OHHS11*FABC)-
     1      ((DC6HOHSRAA_10(r2,cost2,DCOST2R2)*FCAB+DFCABR2*C6HOHS12)*
     1      DAMPHOHS26/rg26+(DDAMPHOHSRAA_10(RG2,6,R2,DRG2R2)*
     1      RG26-6.0D0*
     1      DRG2R2*RG25*DAMPHOHS26)/RG212*C6HOHS12*FCAB)-
     1      ((DC8HOHSRAA_10(r2,cost2,DCOST2R2)*FCAB+DFCABR2*C8HOHS12)*
     1      DAMPHOHS28/rg28+(DDAMPHOHSRAA_10(RG2,8,R2,DRG2R2)*
     1      RG28-8.0D0*
     1      DRG2R2*RG27*DAMPHOHS28)/RG216*C8HOHS12*FCAB)-
     1      ((DC10HOHSRAA_10(r2,cost2,DCOST2R2)*FCAB+DFCABR2*C10HOHS12)*
     1      DAMPHOHS210/rg210+(DDAMPHOHSRAA_10(RG2,10,R2,DRG2R2)*RG210-
     1      10.0D0*DRG2R2*RG29*DAMPHOHS210)/RG220*C10HOHS12*FCAB)-
     1      ((DC6HOHSRAB_10(r3,cost3,DCOST3R2)*FBAC+DFBACR2*C6HOHS13)*
     1      DAMPHOHS36/rg36+(DDAMPHOHSRAB_10(RG3,6,R3,DRG3R2)*
     1      RG36-6.0D0*
     1      DRG3R2*RG35*DAMPHOHS36)/RG312*C6HOHS13*FBAC)-
     1      ((DC8HOHSRAB_10(r3,cost3,DCOST3R2)*FBAC+DFBACR2*C8HOHS13)*
     1      DAMPHOHS38/rg38+(DDAMPHOHSRAB_10(RG3,8,R3,DRG3R2)*
     1      RG38-8.0D0*
     1      DRG3R2*RG37*DAMPHOHS38)/RG316*C8HOHS13*FBAC)-
     1      ((DC10HOHSRAB_10(r3,cost3,DCOST3R2)*FBAC+DFBACR2*C10HOHS13)*
     1      DAMPHOHS310/rg310+(DDAMPHOHSRAB_10(RG3,10,R3,DRG3R2)*RG310-
     1      10.0D0*DRG3R2*RG39*DAMPHOHS310)/RG320*C10HOHS13*FBAC)

CC      VED2=(FBC-1.0D0)*DISHHS_10(R1)+(FAB-1.0D0)*DISOHS_10(R2)+
CC     1                            (FAC-1.0D0)*DISOHS_10(R3)

      DVED21R2=DFBCR2*DISHHS_10(R1)+
     1        DFABR2*DISOHS_10(R2)+(FAB-1.0D0)*
     2        DDISPOH_10(R2,COHS(6),COHS(8),COHS(10),R0OHS,RMOHS)+
     2        DFACR2*DISOHS_10(R3)
     
      DVED31R2_10=DVED31R2_10+DVED21R2
      RETURN
      END

      FUNCTION DVED31R3_10(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/,/CVED_10/,/CVED31_10/)

C     ******************************************************************
      DVED31R3_10=-((DC6OHHSR3_10(r1,cost1,DCOST1R3)*FABC+DFABCR3*
     1      C6OHHS11)*
     1      DAMPOHHS16/rg16+(DDAMPOHHSR3_10(RG1,6,R1,DRG1R3)*RG16-6.0D0*
     1      DRG1R3*RG15*DAMPOHHS16)/RG112*C6OHHS11*FABC)-
     1      ((DC8OHHSR3_10(r1,cost1,DCOST1R3)*FABC+DFABCR3*C8OHHS11)*
     1      DAMPOHHS18/rg18+(DDAMPOHHSR3_10(RG1,8,R1,DRG1R3)*RG18-8.0D0*
     1      DRG1R3*RG17*DAMPOHHS18)/RG116*C8OHHS11*FABC)-
     1      ((DC10OHHSR3_10(r1,cost1,DCOST1R3)*FABC+DFABCR3*C10OHHS11)*
     1      DAMPOHHS110/rg110+(DDAMPOHHSR3_10(RG1,10,R1,DRG1R3)*RG110-
     1      10.0D0*DRG1R3*RG19*DAMPOHHS110)/RG120*C10OHHS11*FABC)-
     1      ((DC6HOHSRAB_10(r2,cost2,DCOST2R3)*FCAB+DFCABR3*C6HOHS12)*
     1      DAMPHOHS26/rg26+(DDAMPHOHSRAB_10(RG2,6,R2,DRG2R3)*
     1      RG26-6.0D0*
     1      DRG2R3*RG25*DAMPHOHS26)/RG212*C6HOHS12*FCAB)-
     1      ((DC8HOHSRAB_10(r2,cost2,DCOST2R3)*FCAB+DFCABR3*C8HOHS12)*
     1      DAMPHOHS28/rg28+(DDAMPHOHSRAB_10(RG2,8,R2,DRG2R3)*
     1      RG28-8.0D0*
     1      DRG2R3*RG27*DAMPHOHS28)/RG216*C8HOHS12*FCAB)-
     1      ((DC10HOHSRAB_10(r2,cost2,DCOST2R3)*FCAB+DFCABR3*C10HOHS12)*
     1      DAMPHOHS210/rg210+(DDAMPHOHSRAB_10(RG2,10,R2,DRG2R3)*RG210-
     1      10.0D0*DRG2R3*RG29*DAMPHOHS210)/RG220*C10HOHS12*FCAB)-
     1      ((DC6HOHSRAA_10(r3,cost3,DCOST3R3)*FBAC+DFBACR3*C6HOHS13)*
     1      DAMPHOHS36/rg36+(DDAMPHOHSRAA_10(RG3,6,R3,DRG3R3)*
     1      RG36-6.0D0*
     1      DRG3R3*RG35*DAMPHOHS36)/RG312*C6HOHS13*FBAC)-
     1      ((DC8HOHSRAA_10(r3,cost3,DCOST3R3)*FBAC+DFBACR3*C8HOHS13)*
     1      DAMPHOHS38/rg38+(DDAMPHOHSRAA_10(RG3,8,R3,DRG3R3)*
     1      RG38-8.0D0*
     1      DRG3R3*RG37*DAMPHOHS38)/RG316*C8HOHS13*FBAC)-
     1      ((DC10HOHSRAA_10(r3,cost3,DCOST3R3)*FBAC+DFBACR3*C10HOHS13)*
     1      DAMPHOHS310/rg310+(DDAMPHOHSRAA_10(RG3,10,R3,DRG3R3)*RG310-
     1      10.0D0*DRG3R3*RG39*DAMPHOHS310)/RG320*C10HOHS13*FBAC)
     
CC      VED2=(FBC-1.0D0)*DISHHS_10(R1)+(FAB-1.0D0)*DISOHS_10(R2)+
CC     1                            (FAC-1.0D0)*DISOHS_10(R3)

      DVED21R3=DFBCR3*DISHHS_10(R1)+
     1        DFABR3*DISOHS_10(R2)+
     2        DFACR3*DISOHS_10(R3)+(FAC-1.0D0)*
     2        DDISPOH_10(R3,COHS(6),COHS(8),COHS(10),R0OHS,RMOHS)
     
      DVED31R3_10=DVED31R3_10+DVED21R3
      RETURN
      END
      
      FUNCTION DAMPHOHS_10(R,N,X)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMHOHS=2.216766092d0
      RO=0.5D0*(RMHOHS+2.5D0*R0HOHSC_10(X))
      DAMPHOHS_10=(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO**2))**N
      RETURN
      END 
      
      FUNCTION DDAMPHOHSR1_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMHOHS=2.216766092d0
      RO=0.5D0*(RMHOHS+2.5D0*R0HOHSC_10(X))
      RO2=RO*RO
      DROR1=0.0D0
      DDAMPHOHSR1_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR1*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR1*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPHOHSRAA_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMHOHS=2.216766092d0
      RO=0.5D0*(RMHOHS+2.5D0*R0HOHSC_10(X))
      RO2=RO*RO
      DROR23=0.5D0*2.5D0*DR0HOHSCR23_10(X)
      DDAMPHOHSRAA_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR23*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR23*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPHOHSRAB_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMHOHS=2.216766092d0
      RO=0.5D0*(RMHOHS+2.5D0*R0HOHSC_10(X))
      RO2=RO*RO
      DROR23=0.0D0
      DDAMPHOHSRAB_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR23*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR23*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
                 
      FUNCTION DAMPOHHS_10(R,N,X)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMOHHS=1.1077763d0
      RO=0.5D0*(RMOHHS+2.5D0*R0OHHSC_10(X))
      DAMPOHHS_10=(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO**2))**N
      RETURN
      END 
      
      FUNCTION DDAMPOHHSR1_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMOHHS=1.1077763d0
      RO=0.5D0*(RMOHHS+2.5D0*R0OHHSC_10(X))
      RO2=RO*RO
      DROR1=0.5D0*2.5D0*DR0OHHSCR1_10(X)
      DDAMPOHHSR1_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR1*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR1*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPOHHSR2_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMOHHS=1.1077763d0
      RO=0.5D0*(RMOHHS+2.5D0*R0OHHSC_10(X))
      RO2=RO*RO
      DROR2=0.0D0
      DDAMPOHHSR2_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR2*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR2*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPOHHSR3_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMOHHS=1.1077763d0
      RO=0.5D0*(RMOHHS+2.5D0*R0OHHSC_10(X))
      RO2=RO*RO
      DROR3=0.0D0
      DDAMPOHHSR3_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR3*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR3*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
                  
      FUNCTION C6OHHS_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c6ohhsae_10/a6s,b6s
!$OMP THREADPRIVATE(/c6ohhsae_10/)

c     ******************************************************************
      A6s=C6OHHSPA_10(R)
      B6s=C6OHHSPE_10(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C6OHHS_10=(1.0D0/3.0D0)*(2.0D0*B6s+A6s)+(2.0D0/3.0D0)*(A6s-B6s)*PL
      RETURN
      END
      
      FUNCTION DC6OHHSR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c6ohhsae_10/a6s,b6s
!$OMP THREADPRIVATE(/c6ohhsae_10/)

c     ******************************************************************
      DC6PE=DC6OHHSPER1_10(R)
      DC6PA=DC6OHHSPAR1_10(R)
      DC6OHHSR1_10=2.0D0/3.0D0*DC6PE+1.0D0/3.0D0*DC6PA
     1        +1.0D0/3.0D0*(DC6PA-DC6PE)*(3.0D0*
     1         COST**2-1.0D0)+6.0D0*D*COST*1.0D0/3.0D0*(A6s-B6s)
      RETURN
      END
      
      FUNCTION DC6OHHSR2_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c6ohhsae_10/a6s,b6s
!$OMP THREADPRIVATE(/c6ohhsae_10/)
      dummy=r
c     ******************************************************************
      DC6OHHSR2_10=2.0D0*D*COST*(A6s-B6s)
      RETURN
      END      
      
      FUNCTION DC6OHHSR3_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c6ohhsae_10/a6s,b6s
!$OMP THREADPRIVATE(/c6ohhsae_10/)
      dummy=r
c     ******************************************************************
      DC6OHHSR3_10=2.0D0*D*COST*(A6s-B6s)
      RETURN
      END      
      
      FUNCTION C8OHHS_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c8ohhsae_10/a8s,b8s
!$OMP THREADPRIVATE(/c8ohhsae_10/)

c     ******************************************************************
      A8s=C8OHHSPA_10(R)
      B8s=C8OHHSPE_10(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C8OHHS_10=(1.0D0/3.0D0)*(2.0D0*B8s+A8s)+(2.0D0/3.0D0)*(A8s-B8s)*PL
      RETURN
      END
      
      FUNCTION DC8OHHSR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c8ohhsae_10/a8s,b8s
!$OMP THREADPRIVATE(/c8ohhsae_10/)
c     ******************************************************************
      DC8PE=DC8OHHSPE_10R1(R)
      DC8PA=DC8OHHSPAR1_10(R)
      DC8OHHSR1_10=2.0D0/3.0D0*DC8PE+1.0D0/3.0D0*DC8PA
     1        +1.0D0/3.0D0*(DC8PA-DC8PE)*(3.0D0*
     1         COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A8s-B8s)
      RETURN
      END
      
      FUNCTION DC8OHHSR2_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c8ohhsae_10/a8s,b8s
!$OMP THREADPRIVATE(/c8ohhsae_10/)
      dummy=r
c     ******************************************************************
      DC8OHHSR2_10=2.0D0*D*COST*(A8s-B8s)
      RETURN
      END
      
      FUNCTION DC8OHHSR3_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c8ohhsae_10/a8s,b8s
!$OMP THREADPRIVATE(/c8ohhsae_10/)
      dummy=r
c     ******************************************************************
      DC8OHHSR3_10=2.0D0*D*COST*(A8s-B8s)
      RETURN
      END
      
      FUNCTION C10OHHS_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhsae_10/a1s,b1s
!$OMP THREADPRIVATE(/c1ohhsae_10/)
c     ******************************************************************
      A1s=C10OHHSPA_10(R)
      B1s=C10OHHSPE_10(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C10OHHS_10=(1.0D0/3.0D0)*(2.0D0*B1s+A1s)+(2.0D0/3.0D0)*
     1  (A1s-B1s)*PL
      RETURN
      END
      
      FUNCTION DC10OHHSR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhsae_10/a1s,b1s
!$OMP THREADPRIVATE(/c1ohhsae_10/)

c     ******************************************************************
      DC10PE=DC10OHHSPE_10R1(R)
      DC10PA=DC10OHHSPAR1_10(R)
      DC10OHHSR1_10=(2.0D0/3.0D0)*DC10PE+(1.0D0/3.0D0)*
     1  DC10PA+(1.0D0/3.0D0)*(DC10PA-DC10PE)
     1   *(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A1s-B1s)
      RETURN
      END    
      
      FUNCTION DC10OHHSR2_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhsae_10/a1s,b1s
!$OMP THREADPRIVATE(/c1ohhsae_10/)
      dummy=r
c     ******************************************************************
      DC10OHHSR2_10=2.0D0*D*COST*(A1s-B1s)
      RETURN
      END
      
      FUNCTION DC10OHHSR3_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhsae_10/a1s,b1s
!$OMP THREADPRIVATE(/c1ohhsae_10/)
      dummy=r
c     ******************************************************************
      DC10OHHSR3_10=2.0D0*D*COST*(A1s-B1s)
      RETURN
      END
                  
      FUNCTION C6HOHS_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      A1=C6HOHSPA_10(R)
      B1=C6HOHSPE_10(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C6HOHS_10=(1.0D0/3.0D0)*(2.0D0*B1+A1)+(2.0D0/3.0D0)*(A1-B1)*PL1
      RETURN
      END
      
       FUNCTION DC6HOHSR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      A1=C6HOHSPA_10(R)
      B1=C6HOHSPE_10(R)
      DC6HOHSR1_10=2.0D0*D*COST*(A1-B1)     
      RETURN
      END
      
      FUNCTION DC6HOHSRAA_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      A1=C6HOHSPA_10(R)
      B1=C6HOHSPE_10(R)
      DC6PEA=DC6HOHSPERAA_10(R)
      DC6PAA=DC6HOHSPARAA_10(R)
      DC6HOHSRAA_10=(2.0D0/3.0D0)*DC6PEA+(1.0D0/3.0D0)*
     1   DC6PAA+(1.0D0/3.0D0)*(DC6PAA-DC6PEA
     1   )*(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A1-B1)     
      RETURN
      END
      
      FUNCTION DC6HOHSRAB_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      A1=C6HOHSPA_10(R)
      B1=C6HOHSPE_10(R)
      DC6HOHSRAB_10=2.0D0*D*COST*(A1-B1)     
      RETURN
      END
    
      FUNCTION C8HOHS_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C8HOHSPA_10(R)
      B1=C8HOHSPE_10(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C8HOHS_10=(1.0D0/3.0D0)*(2.0D0*B1+A1)+(2.0D0/3.0D0)*(A1-B1)*PL1
      RETURN
      END
      
      FUNCTION DC8HOHSR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C8HOHSPA_10(R)
      B1=C8HOHSPE_10(R)
      DC8HOHSR1_10=2.0D0*D*COST*(A1-B1)
      RETURN
      END
      
      FUNCTION DC8HOHSRAA_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C8HOHSPA_10(R)
      B1=C8HOHSPE_10(R)
      DC8PEA=DC8HOHSPERAA_10(R)
      DC8PAA=DC8HOHSPARAA_10(R)
      DC8HOHSRAA_10=(2.0D0/3.0D0)*DC8PEA+(1.0D0/3.0D0)*
     1   DC8PAA+(1.0D0/3.0D0)*(DC8PAA-DC8PEA)
     1   *(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A1-B1)     
      RETURN
      END
      
      FUNCTION DC8HOHSRAB_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C8HOHSPA_10(R)
      B1=C8HOHSPE_10(R)
      DC8HOHSRAB_10=2.0D0*D*COST*(A1-B1)     
      RETURN
      END
     
      FUNCTION C10HOHS_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      A1=C10HOHSPA_10(R)
      B1=C10HOHSPE_10(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C10HOHS_10=(1.0D0/3.0D0)*(2.0D0*B1+A1)+(2.0D0/3.0D0)*(A1-B1)*PL1
      RETURN
      END
      
      FUNCTION DC10HOHSR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      A1=C10HOHSPA_10(R)
      B1=C10HOHSPE_10(R)
      DC10HOHSR1_10=2.0D0*D*COST*(A1-B1)
      RETURN
      END
      
      FUNCTION DC10HOHSRAA_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C10HOHSPA_10(R)
      B1=C10HOHSPE_10(R)
      DC10PEA=DC10HOHSPERAA_10(R)
      DC10PAA=DC10HOHSPARAA_10(R)
      DC10HOHSRAA_10=(2.0D0/3.0D0)*DC10PEA+(1.0D0/3.0D0)*
     1   DC10PAA+(1.0D0/3.0D0)*(DC10PAA-DC10PEA
     1   )*(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A1-B1)     
      RETURN
      END
      
      FUNCTION DC10HOHSRAB_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C10HOHSPA_10(R)
      B1=C10HOHSPE_10(R)
      DC10HOHSRAB_10=2.0D0*D*COST*(A1-B1)     
      RETURN
      END
            
      FUNCTION C6OHHSPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_10/C6HSPA,R0O,ALP,RNH
!$OMP THREADPRIVATE(/C6PA_10/)

c     ******************************************************************
      alp=ALPHHSPA_10(R)
      R0O=R0OHHSC_10(R)
      RNH=RNHHS_10(R)
      C6HSPA=1.5d0*5.62136d0*alp/(sqrt(alp/RNH)+
     1         sqrt(5.62136d0/2.79d0))
      C6OHHSPA_10=C6HSPA   
      RETURN
      END
      
      FUNCTION DC6OHHSPAR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/DC6PA1_10/DC6HSPa1
      COMMON/C6PA_10/C6HSPA,R0O,ALP,RNH
!$OMP THREADPRIVATE(/C6PA_10/,/DC6PA1_10/)
c     ******************************************************************
      const1=1.5d0*5.62136d0 
      DALP=DALPHHSPAR1_10(r)      
      e=0.5d0/SQRT(ALP/RNH)*(dalp*rnh-DRNHHSR1_10(r)*alp)/rnh**2
      d=SQRT(ALP/RNH)+SQRT(5.62136d0/2.79d0)
      DC6HSPa1=(const1*dalp*d-const1*alp*e)/d**2
      DC6OHHSPAR1_10=DC6HSPa1     
      RETURN
      END
      
      FUNCTION C8OHHSPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_10/C6HSPA,R0O,ALP,RNH
!$OMP THREADPRIVATE(/C6PA_10/)
      dummy=r
c     ******************************************************************
      C8OHHSPA_10=1.0d0*C6HSPA*R0O**(1.57243d0)
      RETURN
      END
      
      FUNCTION DC8OHHSPAR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_10/C6HSPA,R0O,ALP,RNH
      COMMON/DC6PA1_10/DC6HSPa1
      COMMON/DRO1_10/DR0HSCR1
!$OMP THREADPRIVATE(/C6PA_10/,/DC6PA1_10/,/DRO1_10/)
c     ****************************************************************** 
      DR0HSCR1=DR0OHHSCR1_10(R)
      DC8OHHSPAR1_10=DC6HSPA1*R0O**(1.57243D0)+1.57243D0*
     1      R0O**(0.57243D0)*DR0HSCR1*C6HSPA
      RETURN
      END

      FUNCTION C10OHHSPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_10/C6HSPA,R0O,ALP,RNH
!$OMP THREADPRIVATE(/C6PA_10/)
      dummy=r
c     ******************************************************************
      C10OHHSPA_10=1.13178d0*C6HSPA*R0O**(1.57243d0*2.0d0)
      RETURN
      END
      
      FUNCTION DC10OHHSPAR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_10/C6HSPA,R0O,ALP,RNH
      COMMON/DC6PA1_10/DC6HSPa1
      COMMON/DRO1_10/DR0HSCR1
!$OMP THREADPRIVATE(/C6PA_10/,/DC6PA1_10/,/DRO1_10/)
c     ******************************************************************
c     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      dummy=r
      DC10OHHSPAR1_10=1.13178D0*DC6HSPA1*R0O**(1.57243D0*
     1     2.0D0)+(1.57243D0*2.0D0)*R0O**(1.57243D0*2.0D0-1.0D0)
     1     *DR0HSCR1*1.13178D0*C6HSPA
      RETURN
      END

      FUNCTION C6OHHSPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_10/C6HSPE,R0O,RNH,alp,DALP
!$OMP THREADPRIVATE(/C6PE_10/)
c     ******************************************************************
      alp=ALPHHSPE_10(R)
      R0O=R0OHHSC_10(R)
      RNH=RNHHS_10(R)      
      dalp=DALPHHSPER1_10(r)      
      C6HSPE=1.5d0*5.62136d0*alp/(sqrt(alp/RNH)+
     1         sqrt(5.62136d0/2.79d0))
      C6OHHSPE_10=C6HSPE    
      RETURN
      END
      
      FUNCTION DC6OHHSPER1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_10/C6HSPE,R0O,RNH,alp,DALP
      COMMON/DC6PE_10/DC6HSPE1
!$OMP THREADPRIVATE(/C6PE_10/,/DC6PE_10/)
c     ******************************************************************
      const1=1.5D0*5.62136d0  
      e=0.5d0/SQRT(ALP/RNH)*(dalp*rnh-DRNHHSR1_10(r)*alp)/rnh**2
      d=SQRT(ALP/RNH)+SQRT(5.62136d0/2.79d0)
      DC6HSPE1=(const1*dalp*d-const1*alp*e)/d**2
      DC6OHHSPER1_10=DC6HSPE1   
      RETURN
      END
      
      FUNCTION C8OHHSPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_10/C6HSPE,R0O,RNH,alp,DALP
!$OMP THREADPRIVATE(/C6PE_10/)
      dummy=r
c     ******************************************************************
c     substitui 1.54 por 1.57243
      C8OHHSPE_10=1.0d0*C6HSPE*R0O**1.57243d0
      RETURN
      END
      
      FUNCTION DC8OHHSPE_10R1(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_10/C6HSPE,R0O,RNH,alp,DALP
      COMMON/DC6PE_10/DC6HSPE1
      COMMON/DROE1_10/DR0HSCE1
!$OMP THREADPRIVATE(/C6PE_10/,/DC6PE_10/,/DROE1_10/)
c     ******************************************************************
c     substitui 1.54 por 1.57243
      DR0HSCE1=DR0OHHSCR1_10(R)
      DC8OHHSPE_10R1=DC6HSPE1*R0O**(1.57243D0)+1.57243D0*
     1      R0O**(1.57243D0-1.0D0)*DR0HSCE1*C6HSPE
      RETURN
      END
       
      FUNCTION C10OHHSPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_10/C6HSPE,R0O,RNH,alp,DALP
!$OMP THREADPRIVATE(/C6PE_10/)
      dummy=r
c     ******************************************************************        
c     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      C10OHHSPE_10=1.13178d0*C6HSPE*R0O**(1.57243d0*2.0d0)
      RETURN
      END
      
      FUNCTION DC10OHHSPE_10R1(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_10/C6HSPE,R0O,RNH,alp,DALP
      COMMON/DC6PE_10/DC6HSPE1
      COMMON/DROE1_10/DR0HSCE1
!$OMP THREADPRIVATE(/C6PE_10/,/DC6PE_10/,/DROE1_10/)
c     ******************************************************************
c     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      dummy=r
      DC10OHHSPE_10R1=1.13178D0*DC6HSPE1*R0O**(1.57243D0*
     1    2.0D0)+(1.57243D0*2.0D0)*R0O**(1.57243D0*2.0D0-1.0D0)*
     1     DR0HSCE1*1.13178D0*C6HSPE
      RETURN
      END

      FUNCTION C6HOHSPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      ALP=ALPOHSPA_10(R)
      C6HOHSPA_10=1.5d0*4.500d0*ALP/(SQRT(ALP/RNOHS_10(R))+2.25d0)
      RETURN
      END
      
      FUNCTION DC6HOHSPARAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      const1=1.5D0*4.500D0
      ALP=ALPOHSPA_10(R)
      RNO=RNOHS_10(R)
      dalp=DALPOHSPAR23_10(r)      
      d=SQRT(ALP/RNO)+2.25D0
      e=0.5d0/SQRT(ALP/RNO)*(dalp*rno-DRNOHSR23_10(r)*alp)/rno**2
      DC6HOHSPARAA_10=(const1*dalp*d-e*const1*alp)/d**2
      RETURN
      END
      
      FUNCTION C8HOHSPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      C8HOHSPA_10=1.0d0*C6HOHSPA_10(R)*R0HOHSC_10(R)**(1.57243d0)      
      RETURN
      END
      
      FUNCTION DC8HOHSPARAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      R0H=R0HOHSC_10(R)     
      DC8HOHSPARAA_10=DC6HOHSPARAA_10(R)*R0H**1.57243D0+1.57243D0*
     1      R0H**(1.57243D0-1.0D0)*DR0HOHSCR23_10(R)*C6HOHSPA_10(R)
      RETURN
      END
      
      FUNCTION C10HOHSPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
c     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      C10HOHSPA_10=1.13178d0*C6HOHSPA_10(r)*
     1 R0HOHSC10_10(r)**(2.0d0*1.57243d0)      
      RETURN
      END
      
      FUNCTION DC10HOHSPARAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      R0H=R0HOHSC10_10(R)      
      DC10HOHSPARAA_10=1.13178D0*DC6HOHSPARAA_10(R)*R0H**(2.0D0*
     1  1.57243D0)+2.0D0*1.57243D0*R0H**(2.0D0*1.57243D0-
     1    1.0D0)*DR0HOHSC10R23_10(R)*1.13178D0*C6HOHSPA_10(R)
      RETURN
      END

      FUNCTION C6HOHSPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      ALP=ALPOHSPE_10(r)      
      C6HOHSPE_10=1.5d0*4.500d0*ALP/(SQRT(ALP/RNOHS_10(R))+2.25d0)   
      RETURN
      END
      
      FUNCTION DC6HOHSPERAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      const1=1.5D0*4.500D0
      RNO=RNOHS_10(R)
      ALP=ALPOHSPE_10(R)
      dalp=DALPOHSPER23_10(r)      
      d=SQRT(ALP/RNO)+2.25D0
      e=0.5d0/SQRT(ALP/RNO)*(dalp*rno-DRNOHSR23_10(r)*alp)/rno**2
      DC6HOHSPERAA_10=(const1*dalp*d-e*const1*alp)/d**2
      RETURN
      END
      
      FUNCTION C8HOHSPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
c     substitui 1.54 por 1.57243
      C8HOHSPE_10=1.0d0*C6HOHSPE_10(r)*R0HOHSC_10(r)**(1.57243d0)      
      RETURN
      END

      FUNCTION DC8HOHSPERAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
c     substitui 1.54 por 1.57243
      R0H=R0HOHSC_10(R)      
      DC8HOHSPERAA_10=DC6HOHSPERAA_10(R)*R0H**1.57243D0+1.57243D0*
     1      R0H**(1.57243D0-1.0D0)*DR0HOHSCR23_10(R)*C6HOHSPE_10(R)
      RETURN
      END
      
      FUNCTION C10HOHSPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
c     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      C10HOHSPE_10=1.13178d0*C6HOHSPE_10(R)*
     1 R0HOHSC10_10(R)**(2.0d0*1.57243d0)      
      RETURN
      END
      
      FUNCTION DC10HOHSPERAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
C     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      R0H=R0HOHSC10_10(R)            
      DC10HOHSPERAA_10=1.13178D0*DC6HOHSPERAA_10(R)*R0H**(2.0D0*
     1  1.57243D0)+2.0D0*1.57243D0*R0H**(2.0D0*1.57243D0-
     1    1.0D0)*DR0HOHSC10R23_10(R)*1.13178D0*C6HOHSPE_10(R)
      RETURN
      END
      
      FUNCTION RNOHS_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      Ainf=4.2594d0
      re=1.8344d0
      a=-0.2963609d0
      c=0.584820271D0/1.04D0 
      RNOHS_10=ainf+a*exp(-c*(R-re))
      RETURN
      END
            
      FUNCTION DRNOHSR23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************      
      Ainf=4.2594d0
      re=1.8344d0
      a=-0.2963609d0
      c=0.584820271D0/1.04D0 
      DRNOHSR23_10=-A*C*exp(-c*(R-re))
      RETURN
      END
      
      FUNCTION RNHHS_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      Ainf=1.78d0
      re=1.449d0
      a=-0.1478d0
      c=0.60634691d0          
      RNHHS_10=ainf+a*exp(-c*(R-re))
      RETURN
      END
      
      FUNCTION DRNHHSR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      Ainf=1.78d0
      re=1.449d0
      a=-0.1478d0
      c=0.60634691d0          
      DRNHHSR1_10=-a*C*exp(-c*(R-re))
      RETURN
      END 

      FUNCTION ALPHHSPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7/9.0D0,6.3134223d0,-2.8677438d0,
     1 7.7770980d-1,-7.6966611d-1,2.8116781d-1,4.1567195d-3,1.5074557d2/
c     ******************************************************************
      azero=1.383d0-ainf1
      rr2=r*r
      rr3=rr2*r
      rr5=rr2*rr3
      ALPHHSPA_10=Ainf1+(azero+C1*r+C2*rr2+c3*rr3)*
     1   EXP(-C4*r-C5*rr2)+(1-EXP(-C6*rr5))*C7/rr3
      RETURN
      END
      
      FUNCTION DALPHHSPAR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7/9.0D0,6.3134223d0,-2.8677438d0,
     1 7.7770980d-1,-7.6966611d-1,2.8116781d-1,4.1567195d-3,1.5074557d2/
c     ******************************************************************
      azero=1.383d0-ainf1
      rr2=r*r
      rr3=rr2*r
      rr4=rr3*r
      rr5=rr4*r
      DALPHHSPAR1_10=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*r-C5*RR2)-
     1          (C4+2.0D0*C5*R)*(azero+C1*r+C2*RR2+c3*RR3)*EXP(-C4*r
     1          -C5*RR2)-3.0D0*C7/RR4+3.0D0*C7/RR4*EXP(-C6*RR5)+
     1          5.0D0*C6*C7*R*EXP(-C6*RR5)
      RETURN
      END
      
      FUNCTION ALPHHSPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7/9.0D0,4.8750306d0,-1.8684334d0,3.8
     1 222971d-1,-5.6354499d-1,2.8106756d-1,1.7861826d-3,-2.8274363d1/    
c     ******************************************************************
      azero=1.383d0-ainf1
      rr2=r*r
      rr3=rr2*r
      rr5=rr2*rr3
      ALPHHSPE_10=Ainf1+(azero+C1*r+C2*RR2+c3*RR3)*
     1   EXP(-C4*r-C5*RR2)+(1-EXP(-C6*RR5))*C7/RR3
      RETURN
      END
      
      FUNCTION DALPHHSPER1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7/9.0D0,4.8750306d0,-1.8684334d0,3.8
     1 222971d-1,-5.6354499d-1,2.8106756d-1,1.7861826d-3,-2.8274363d1/    
c     ******************************************************************
      azero=1.383d0-ainf1
      rr2=r*r
      rr3=rr2*r
      rr4=rr3*r
      rr5=rr4*r
      DALPHHSPER1_10=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*r-C5*RR2)-
     1          (C4+2.0D0*C5*R)*(azero+C1*r+C2*RR2+c3*RR3)*EXP(-C4*r
     1          -C5*RR2)-3.0D0*C7/RR4+3.0D0*C7/RR4*EXP(-C6*RR5)+
     1          5.0D0*C6*C7*R*EXP(-C6*RR5)
      RETURN
      END
      
      FUNCTION ALPOHSPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7,c8/8.7752955d0,7.7784665d0,-3.77
     1   08106d0,6.3231864d-1,-1.1225788d0,-7.2086523d-1,2.7382895d-1,
     2   1.7301967d-5,2.9244902d4/
c     ******************************************************************
      Azero=3.13939d0-Ainf1
      rr2=r*r
      rr3=rr2*r
      ALPOHSPA_10=Ainf1+(Azero+C1*r+C2*RR2+C3*RR3)*EXP(
     1 -C4*r-C5*RR2-C6*RR3)+(1-EXP(-C7*r**8))*C8/r**6
      RETURN
      END
      
      FUNCTION DALPOHSPAR23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7,c8/8.7752955d0,7.7784665d0,-3.77
     1   08106d0,6.3231864d-1,-1.1225788d0,-7.2086523d-1,2.7382895d-1,
     2   1.7301967d-5,2.9244902d4/
c     ******************************************************************
      Azero=3.13939d0-Ainf1
      rr2=r*r
      rr3=rr2*r
      rr7=rr2*RR2*RR3
      rr8=rr2*RR3*RR3
      DALPOHSPAR23_10=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*r-C5*RR2-C6
     1          *RR3)-(Azero+C1*r+C2*RR2+C3*RR3)*(C4+2.0D0*C5*R+
     1          3.0D0*C6*RR2)*EXP(-C4*r-C5*RR2-C6*RR3)-6.0D0*C8/RR7
     1          +6.0D0*C8/RR7*EXP(-C7*RR8)+8.0D0*R*C7*C8*EXP(-C7*RR8)
      RETURN
      END
      
      FUNCTION ALPOHSPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5/9.6378446d0,-1.5996623d-3,3.2213045d-1,
     1   1.4377753d-2,2.3615487d-3,-4.6276713d1/
c     ******************************************************************
      Azero=ainf1-3.5396940d0
      ALPOHSPE_10=Ainf1-Azero*EXP(-C1*r-C2*r**2-C3*
     1       r**3)+(1-EXP(-C4*r**5))*C5/r**3
      RETURN
      END
      
      FUNCTION DALPOHSPER23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5/9.6378446d0,-1.5996623d-3,3.2213045d-1,
     1   1.4377753d-2,2.3615487d-3,-4.6276713d1/
c     ******************************************************************
      Azero=ainf1-3.5396940d0
      rr2=r*r
      rr3=rr2*r
      RR4=RR2*RR2
      RR5=RR4*R
      DALPOHSPER23_10=Azero*(C1+2.0D0*C2*r+3.0D0*
     1         C3*RR2)*EXP(-C1*r-C2*RR2
     1         -C3*RR3)-3.0D0*c5/RR4+3.0D0*C5/RR4*EXP(-C4*RR5)
     1         +5.0D0*C4*C5*R*EXP(-C4*RR5)
      RETURN
      END 
      
      FUNCTION R0OHHSC_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************       
      finf=6.334299622d0/6.9916231d0
      fzero=5.047147656d0/5.0591314d0
      d=fzero-finf
      c=0.60634691d0      
      FCORR=finf+d*exp(-c*r)
      ALPHHSM=(1.0d0/3.0d0)*ALPHHSPA_10(r)+(2.0d0/3.0d0)*ALPHHSPE_10(r)       
      R0OHHSC_10=2.0d0*((alphhsm**(1.0d0/3.0d0))+sqrt(2.003423d0))*FCORR
      RETURN
      END
      
      FUNCTION DR0OHHSCR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************       
      finf=6.334299622d0/6.9916231d0
      fzero=5.047147656d0/5.0591314d0
      d=fzero-finf
      c=0.60634691d0      
      FCORR=finf+d*exp(-c*r)
      DFCORR=-C*D*exp(-c*r)
      ALPHHSM=(1.0d0/3.0d0)*ALPHHSPA_10(r)+(2.0d0/3.0d0)*ALPHHSPE_10(r)       
      DALPHHSMR1=(1.0d0/3.0d0)*DALPHHSPAR1_10(r)+(2.0d0/3.0d0)*
     1          DALPHHSPER1_10(r)
      DR0OHHSCR1_10=2.0D0*DFCORR*((alphhsm**(1.0d0/3.0d0))+sqrt(2.003423
     1       d0))+(1.0d0/3.0d0)*alphhsm**(-2.0d0/3.0d0)*DALPHHSMR1*
     1       2.0D0*FCORR
      RETURN
      END
      
      FUNCTION R0HOHSC_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      finf=7.86489501d0/7.6781675d0
      fzero=5.950512078d0/6.473336d0      
      d=fzero-finf
      c=0.584820271D0/1.04D0    
      FCORR=finf+d*exp(-c*r)
      ALPOHSM=(1.0d0/3.0d0)*ALPOHSPA_10(R)+(2.0d0/3.0d0)*ALPOHSPE_10(R)    
      R0HOHSC_10=2.0d0*((ALPOHSM**(1.0d0/3.0d0))+(6.9282d0/4.0d0))*FCORR
      RETURN
      END
      
      FUNCTION DR0HOHSCR23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      finf=7.86489501d0/7.6781675d0
      fzero=5.950512078d0/6.473336d0      
      d=fzero-finf
      c=0.584820271D0/1.04D0    
      FCORR=finf+d*exp(-c*r)
      ALPOHSM=(1.0d0/3.0d0)*ALPOHSPA_10(R)+(2.0d0/3.0d0)*ALPOHSPE_10(R)    
      DALPHAM=1.0d0/3.0d0*DALPOHSPAR23_10(R)+2.0d0/3.0d0*
     1 DALPOHSPER23_10(R)
      DR0HOHSCR23_10=2.0d0/3.0D0*alpOHSM**(-2.0d0/3.0d0)*DALPHAM*fcorr-
     1    2.0D0*c*d*exp(-c*r)*(alpOHSm**(1.0d0/3.0d0)+6.9282d0/4.0d0) 
      RETURN
      END 
     
      FUNCTION R0HOHSC10_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      finf=7.235196035d0/7.6781675d0
      fzero=5.950512078d0/6.473336d0      
      d=fzero-finf
      c=0.584820271D0/1.04D0    
      FCORR=finf+d*exp(-c*r)
      ALPOHSM=(1.0d0/3.0d0)*ALPOHSPA_10(R)+(2.0d0/3.0d0)*ALPOHSPE_10(R)    
      R0HOHSC10_10=2.0d0*((ALPOHSM**(1.0d0/3.0d0))+(6.9282d0/4.0d0))*
     1 FCORR      
      RETURN
      END
      
      FUNCTION DR0HOHSC10R23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************     
      finf=7.235196035d0/7.6781675d0
      fzero=5.950512078d0/6.473336d0      
      d=fzero-finf
      c=0.584820271D0/1.04D0    
      FCORR=finf+d*exp(-c*r)
      ALPOHSM=(1.0d0/3.0d0)*ALPOHSPA_10(R)+(2.0d0/3.0d0)*ALPOHSPE_10(R)    
      DALPHAM=1.0d0/3.0d0*DALPOHSPAR23_10(R)+2.0d0/3.0d0*
     1 DALPOHSPER23_10(R)
      DR0HOHSC10R23_10=2.0d0/(3.0D0*alpOHSM**(2.0d0/3.0d0))*
     1 DALPHAM*fcorr
     1  -2.0D0*c*d*(alpOHSM**(1.0d0/3.0d0)+(6.9282d0/4.0d0))*exp(-c*r) 
      RETURN
      END 
      
      FUNCTION FDOHS_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7,C8/-5.3537295E-01,6.4491094E-01,8.53312
     1     26E-03,1.2881749E+00,-7.4776246E-01,1.8391757E-01,1.3506306
     2     E-02,5.5330102E+00/
c     ******************************************************************
      RR2=R*R
      RR3=RR2*R
      RR5=RR2*RR3
      FDOHS_10=(C1*R+C2*RR2+C3*RR3)*EXP(-C4*R-C5*RR2-C6*RR3)+(1-EXP(-C7
     1       *RR5))*C8/RR3
      return
      end
      
      FUNCTION DFDOHSR23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7,C8/-5.3537295E-01,6.4491094E-01,8.53312
     1     26E-03,1.2881749E+00,-7.4776246E-01,1.8391757E-01,1.3506306
     2     E-02,5.5330102E+00/
c     ******************************************************************
c     DIPOLO PARA O OH(2SI)      
      RR2=R*R
      RR3=RR2*R
      RR4=RR3*R
      RR5=RR2*RR3
      DFDOHSR23_10=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*R-C5*RR2-C6*RR3)
     1        -(C4+2.0D0*C5*R+3.0D0*C6*RR2)*(C1*R+C2*RR2+C3*RR3)*EXP
     1        (-C4*R-C5*RR2-C6*RR3)-3.0D0*C8/RR4+3.0D0*C8/RR4*EXP(
     1        -C7*RR5)+5.0D0*C7*C8*R*EXP(-C7*RR5)
      return
      end
      
      FUNCTION VINDHOHS_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB

!$OMP THREADPRIVATE(/CVED_10/,/CVED31_10/)
      dummy=r1
c     ******************************************************************
      VINDHOHS_10=-FDOHS_10(R2)**2*4.5D0*(3.0D0*COST2**2+1.0D0)*
     1          DAMPHOHS26/(2.0D0*rG26)
     1         -FDOHS_10(R3)**2*4.5D0*(3.0D0*COST3**2+1.0D0)*
     1          DAMPHOHS36/(2.0D0*rG36)

      return
      end
      
      FUNCTION DVINDHOHSR1_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/DVINDR_10/TA,TB,XB,XB1,FDOHSR2,FDOHSR3,TA1,TA2,XA,XA1

!$OMP THREADPRIVATE(/CVED_10/,/CVED31_10/,/DVINDR_10/)

c     ******************************************************************
      TA=3.0D0*COST2**2+1.0D0    
      DTAR1=-6.0D0*(R3-R1)/R2**2
      TB=3.0D0*COST3**2+1.0D0    
      DTBR1=-6.0D0*(R2-R1)/R3**2
      DFDOHSR1=0.0D0
      xb=DAMPHOHS26/(2.0D0*rG26) 
      xb1=DAMPHOHS36/(2.0D0*rG36) 
      FDOHSR2=FDOHS_10(R2)
      FDOHSR3=FDOHS_10(R3)
      TA1=-4.5D0*FDOHSR2**2
      TA2=-4.5D0*FDOHSR3**2
      xa=TA1*TA  
      xa1=TA2*TB  
      DVINDHOHSR1_10=xb*(-4.5d0*2.0d0*FDOHSR2*DFDOHSR1*Ta+dTar1*TA1)
     1 +xa*(DDAMPHOHSR1_10(Rg2,6,r2,DRG2R1)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R1*DAMPHOHS26)/(4.0D0*RG212)
     1         +xb1*(-4.5d0*2.0d0*FDOHSR3*DFDOHSR1*TB+dTBr1*TA2)
     1 +Xa1*(DDAMPHOHSR1_10(Rg3,6,r3,DRG3R1)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R1*DAMPHOHS36)/(4.0D0*RG312) 
      return
      end
            
      FUNCTION DVINDHOHSR2_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/DVINDR_10/TA,TB,XB,XB1,FDOHSR2,FDOHSR3,TA1,TA2,XA,XA1
!$OMP THREADPRIVATE(/CVED_10/,/CVED31_10/,/DVINDR_10/)

c     ******************************************************************
      DTAR2=-6.0D0*(R3-R1)**2/R2**3
      DTBR2=6.0D0*(R2-R1)/R3**2
      DFDOHSR23R3=0.0D0
      DVINDHOHSR2_10=xb*(-4.5d0*2.0d0*FDOHSR2*DFDOHSR23_10(R2)*Ta+dTar2
     1 *TA1)+xa*(DDAMPHOHSRAA_10(Rg2,6,r2,DRG2R2)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R2*DAMPHOHS26)/(4.0D0*RG212)
     1      +xb1*(-4.5d0*2.0d0*FDOHSR3*DFDOHSR23R3*TB+dTBr2*TA2)
     1 +Xa1*(DDAMPHOHSRAB_10(Rg3,6,r3,DRG3R2)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R2*DAMPHOHS36)/(4.0D0*RG312)
      return
      end      
      
      FUNCTION DVINDHOHSR3_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31_10/C6OHHS11,C8OHHS11,C10OHHS11,C6HOHS12,C8HOHS12,
     1    C10HOHS12,C6HOHS13,C8HOHS13,C10HOHS13,DAMPOHHS16,
     1    DAMPOHHS18,DAMPOHHS110,DAMPHOHS26,DAMPHOHS28,DAMPHOHS210,
     1    DAMPHOHS36,DAMPHOHS38,DAMPHOHS310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/DVINDR_10/TA,TB,XB,XB1,FDOHSR2,FDOHSR3,TA1,TA2,XA,XA1
!$OMP THREADPRIVATE(/CVED_10/,/CVED31_10/,/DVINDR_10/)

c     ******************************************************************
      DTAR3=6.0D0*(R3-R1)/R2**2
      DTBR3=-6.0D0*(R2-R1)**2/R3**3
      DFDOHSR23R2=0.0D0
      DVINDHOHSR3_10=xb*(-4.5d0*2.0d0*FDOHSR2*DFDOHSR23R2*Ta+dTar3
     1 *TA1)+xa*(DDAMPHOHSRAB_10(Rg2,6,r2,DRG2R3)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R3*DAMPHOHS26)/(4.0D0*RG212)
     1    +xb1*(-4.5d0*2.0d0*FDOHSR3*DFDOHSR23_10(R3)*TB+dTBr3*TA2)
     1 +Xa1*(DDAMPHOHSRAA_10(Rg3,6,r3,DRG3R3)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R3*DAMPHOHS36)/(4.0D0*RG312)
      return
      end
    
      FUNCTION FQHHS_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7,C8/9.2292863d-02,-1.0187873d-02,1.03075
     1     21d-02,-1.0635640d0,2.8476217d-01,1.1961650d-08,1.3473360
     2     d-02,3.9037435d+00/
c     ******************************************************************
      RR2=R*R
      RR3=RR2*R
      RR5=RR3*RR2
      FQHHS_10=(C1*R+C2*RR2+C3*RR3)*EXP(-C4*R
     1      -C5*RR2-C6*RR3)+(1-EXP(-C7*RR5))*C8/RR3      
      return
      end
      
      FUNCTION DFQHHSR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7,C8/9.2292863d-02,-1.0187873d-02,1.03075
     1     21d-02,-1.0635640d0,2.8476217d-01,1.1961650d-08,1.3473360
     2     d-02,3.9037435d+00/
c     ******************************************************************
      RR2=R*R
      RR3=RR2*R
      RR4=RR3*R
      RR5=RR3*RR2
      DFQHHSR1_10=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*R-C5*RR2-C6*RR3)
     1        -(C4+2.0D0*C5*R+3.0D0*C6*RR2)*(C1*R+C2*RR2+C3*RR3)*EXP
     1        (-C4*R-C5*RR2-C6*RR3)-3.0D0*C8/RR4+3.0D0*C8/RR4*EXP(
     1        -C7*RR5)+5.0D0*C7*C8*R*EXP(-C7*RR5)
      return
      end
      
      FUNCTION VEOHHS_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      DATA b0,b2,b4,b6/-3.68666806d0,-0.0391352778d0,-0.277479722d0,
     1                 0.0402558333d0/
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3

!$OMP THREADPRIVATE(/CVED_10/)
      dummy=r2
      dummy=r3
c     ******************************************************************
c     o valor do quadrupolo do O(1D)=1.233793Buckingham=0.917288208au
      at=b0+b2*cos2t1+b4*cos4t1+b6*cos6t1
      VEOHHS_10=0.75d0*FQHHS_10(R1)*0.917288208D0*at*DAMPOHHS5
      return
      end
      
      FUNCTION DVEOHHSR1_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      DATA b0,b2,b4,b6/-3.68666806d0,-0.0391352778d0,-0.277479722d0,
     1                 0.0402558333d0/
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/DVEOHS_10/AT,DAT,FQHHSR1
!$OMP THREADPRIVATE(/CVED_10/,/DVEOHS_10/)
      dummy=r2
      dummy=r3
c     ******************************************************************
      at=b0+b2*cos2t1+b4*cos4t1+b6*cos6t1
      DAT=B2*dcos2t+b4*dcos4t+b6*dcos6t
      DATR1=DAT*dcost1r1
      CONSTV=0.75d0*0.917288208d0
      FQHHSR1=FQHHS_10(R1)
      DVEOHHSR1_10=(DFQHHSR1_10(R1)*AT+DATR1*FQHHSR1)*CONSTV*DAMPOHHS5+
     1     (DDAMPOHHSR1_10(Rg1,5,r1,drg1r1)*rg15-5.0D0*DRG1R1*
     1     RG14*DAMPOHHS5*RG15)/(rg110)*at*CONSTV*FQHHSR1
      return
      end
      
      FUNCTION DVEOHHSR2_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/DVEOHS_10/AT,DAT,FQHHSR1

!$OMP THREADPRIVATE(/CVED_10/,/DVEOHS_10/)
      dummy=r2
      dummy=r3
c     ******************************************************************
c     o valor do quadrupolo do O(1D)=1.233793Buckingham=0.917288208au
      DATR2=DAT*dcost1r2
      CONSTV=0.75d0*0.917288208d0
      DVEOHHSR2_10=(DATR2*FQHHSR1)*CONSTV*DAMPOHHS5+
     1        (DDAMPOHHSR2_10(Rg1,5,r1,DRG1R2)*rg15-5.0D0*DRG1R2*
     1        RG14*DAMPOHHS5*Rg15)/(rg110)*at*CONSTV*FQHHSR1
      return
      end
      
      FUNCTION DVEOHHSR3_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c      DATA b0,b2,b4,b6/-3.68666806d0,-0.0391352778d0,-0.277479722d0,
c     1                 0.0402558333d0/
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/DVEOHS_10/AT,DAT,FQHHSR1

!$OMP THREADPRIVATE(/CVED_10/,/DVEOHS_10/)
      dummy=r2      
      dummy=r3

c     ******************************************************************
c     o valor do quadrupolo do O(1D)=1.233793Buckingham=0.917288208au
      DATR3=DAT*dcost1r3
      CONSTV=0.75d0*0.917288208d0
      DVEOHHSR3_10=(DATR3*FQHHSR1)*CONSTV*DAMPOHHS5+
     1    (DDAMPOHHSR3_10(Rg1,5,r1,DRG1R3)*rg15-5.0D0*DRG1R3*
     1    RG14*DAMPOHHS5*RG15)/(rg110)*at*CONSTV*FQHHSR1
      return
      end
      
      FUNCTION RL1_10(r1E,r2E,r3E)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      R1=R1E+1.0D-8
      R2=R2E+1.0D-8
      R3=R3E+1.0D-8
      CALL GEOMRL_10(R1,R2,R3)
c    ***************************************************************            

      RL1_10=VED31_10(r1,r2,r3)+VEOHHS_10(r1,r2,r3)+
     1 VINDHOHS_10(r1,r2,r3)
      return
      end
      
      SUBROUTINE DRL1_10(R1E,R2E,R3E,DER1,DER2,DER3)
c     ******************************************************************
c     TO COMPUTE THE DERIVATIVES OF THE LONG RANGE TERM
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      R1=R1E+1.0D-8
      R2=R2E+1.0D-8
      R3=R3E+1.0D-8
      CALL GEOMRL_10(R1,R2,R3)
      DER1=DVED31R1(r1,r2,r3)+DVEOHHSR1_10(r1,r2,r3)+
     1 DVINDHOHSR1_10(r1,r2,r3)

      DER2=DVED31R2_10(r1,r2,r3)+DVEOHHSR2_10(r1,r2,r3)+
     1 DVINDHOHSR2_10(r1,r2,r3)

      DER3=DVED31R3_10(r1,r2,r3)+DVEOHHSR3_10(r1,r2,r3)+
     1 DVINDHOHSR3_10(r1,r2,r3)
      RETURN      
      END
      
      FUNCTION VED32_10(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
!$OMP THREADPRIVATE(/DIATDI_10/,/CVED_10/,/CVED32_10/)

c     ******************************************************************
      VED3=
     1     -C6OHHT11*FABC/rg16*DAMPOHHT16-
     1      C8OHHT11*FABC/rg18*DAMPOHHT18-
     2      C10OHHT11*FABC/rg110*DAMPOHHT110
     3     -C6HOHP12*FCAB/rg26*DAMPHOHP26-
     4      C8HOHP12*FCAB/rg28*DAMPHOHP28-
     5      C10HOHP12*FCAB/rg210*DAMPHOHP210
     6     -C6HOHP13*FBAC/rg36*DAMPHOHP36-
     7      C8HOHP13*FBAC/rg38*DAMPHOHP38-
     8      C10HOHP13*FBAC/rg310*DAMPHOHP310
     
      VED2=(FBC-1.0D0)*DISHHT_10(R1)+(FAB-1.0D0)*DISOHP_10(R2)+
     1                            (FAC-1.0D0)*DISOHP_10(R3)
      
      
      VED32_10=VED3+VED2
     
C      PRINT*, 'R1, R2, R3', R1,R2,R3
C      PRINT*, 'RG1, RG2,RG3', RG1,RG2,RG3
C      PRINT*, 'FAB,FAC,FBC',FAB,FAC,FBC
C      PRINT*,'FABC,FCAB,FBAC',FABC,FCAB,FBAC
C      PRINT*,'C6OHHT,C8OHHT,C10OHHT',C6OHHT11,C8OHHT11,C10OHHT11
C      PRINT*,'C6HOHP1,C8HOHP,C10HOHP',C6HOHP12,C8HOHP12,C10HOHP12
C      PRINT*,'C6HOHP2,C8HOHP,C10HOHP',C6HOHP13,C8HOHP13,C10HOHP13
C      PRINT*,'VED32_10',VED32_10      
C      PRINT*,'VED3,VED2',VED3,VED2
      RETURN
      END

      FUNCTION DVED32R1_10(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/,/CVED_10/,/CVED32_10/)

c     ******************************************************************
      DVED32R1_10=-((DC6OHHTR1_10(r1,cost1,DCOST1R1)*FABC+DFABCR1*
     1      C6OHHT11)*
     1      DAMPOHHT16/rg16+(DDAMPOHHTR1_10(RG1,6,R1,DRG1R1)*RG16-6.0D0*
     1      DRG1R1*RG15*DAMPOHHT16)/RG112*C6OHHT11*FABC)-
     1      ((DC8OHHTR1_10(r1,cost1,DCOST1R1)*FABC+DFABCR1*C8OHHT11)*
     1      DAMPOHHT18/rg18+(DDAMPOHHTR1_10(RG1,8,r1,DRG1R1)*RG18-8.0D0*
     1      DRG1R1*RG17*DAMPOHHT18)/RG116*C8OHHT11*FABC)-
     1      ((DC10OHHTR1_10(r1,cost1,DCOST1R1)*FABC+DFABCR1*C10OHHT11)*
     1      DAMPOHHT110/rg110+(DDAMPOHHTR1_10(RG1,10,R1,DRG1R1)*RG110-
     1      10.0D0*DRG1R1*RG19*DAMPOHHT110)/RG120*C10OHHT11*FABC)-
     1      ((DC6HOHPR1_10(r2,cost2,DCOST2R1)*FCAB+DFCABR1*C6HOHP12)*
     1      DAMPHOHP26/rg26+(DDAMPHOHPR1_10(RG2,6,R2,DRG2R1)*RG26-6.0D0*
     1      DRG2R1*RG25*DAMPHOHP26)/RG212*C6HOHP12*FCAB)-
     1      ((DC8HOHPR1_10(r2,cost2,DCOST2R1)*FCAB+DFCABR1*C8HOHP12)*
     1      DAMPHOHP28/rg28+(DDAMPHOHPR1_10(RG2,8,R2,DRG2R1)*RG28-8.0D0*
     1      DRG2R1*RG27*DAMPHOHP28)/RG216*C8HOHP12*FCAB)-
     1      ((DC10HOHPR1_10(r2,cost2,DCOST2R1)*FCAB+DFCABR1*C10HOHP12)*
     1      DAMPHOHP210/rg210+(DDAMPHOHPR1_10(RG2,10,R2,DRG2R1)*RG210-
     1      10.0D0*DRG2R1*RG29*DAMPHOHP210)/RG220*C10HOHP12*FCAB)-
     1      ((DC6HOHPR1_10(r3,cost3,DCOST3R1)*FBAC+DFBACR1*C6HOHP13)*
     1      DAMPHOHP36/rg36+(DDAMPHOHPR1_10(RG3,6,R3,DRG3R1)*RG36-
     1      6.0D0*DRG3R1*RG35*DAMPHOHP36)/RG312*C6HOHP13*FBAC)-
     1      ((DC8HOHPR1_10(r3,cost3,DCOST3R1)*FBAC+DFBACR1*C8HOHP13)*
     1      DAMPHOHP38/rg38+(DDAMPHOHPR1_10(RG3,8,R3,DRG3R1)*RG38-
     1      8.0D0*DRG3R1*RG37*DAMPHOHP38)/RG316*C8HOHP13*FBAC)-
     1      ((DC10HOHPR1_10(r3,cost3,DCOST3R1)*FBAC+DFBACR1*C10HOHP13)*
     1      DAMPHOHP310/rg310+(DDAMPHOHPR1_10(RG3,10,R3,DRG3R1)*RG310-
     1    10.0D0*DRG3R1*RG39*DAMPHOHP310)/RG320*C10HOHP13*FBAC)
     
CC      VED2=(FBC-1.0D0)*DISHHT_10(R1)+(FAB-1.0D0)*DISOHP_10(R2)+
CC     1                            (FAC-1.0D0)*DISOHP_10(R3)

      DVED22R1=DFBCR1*DISHHT_10(R1)+(FBC-1.0D0)*DDISHHT_10(R1)+
     1        DFABR1*DISOHP_10(R2)+
     2        DFACR1*DISOHP_10(R3)
     
      DVED32R1_10=DVED32R1_10+DVED22R1

      RETURN
      END

      FUNCTION DVED32R2_10(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/,/CVED_10/,/CVED32_10/)

c     ******************************************************************
      DVED32R2_10=-((DC6OHHTR2_10(r1,cost1,DCOST1R2)*
     1      FABC+DFABCR2*C6OHHT11)*
     1      DAMPOHHT16/rg16+(DDAMPOHHTR2_10(RG1,6,R1,DRG1R2)*RG16-6.0D0*
     1      DRG1R2*RG15*DAMPOHHT16)/RG112*C6OHHT11*FABC)-
     1      ((DC8OHHTR2_10(r1,cost1,DCOST1R2)*FABC+DFABCR2*C8OHHT11)*
     1      DAMPOHHT18/rg18+(DDAMPOHHTR2_10(RG1,8,R1,DRG1R2)*RG18-8.0D0*
     1      DRG1R2*RG17*DAMPOHHT18)/RG116*C8OHHT11*FABC)-
     1      ((DC10OHHTR2_10(r1,cost1,DCOST1R2)*FABC+DFABCR2*C10OHHT11)*
     1      DAMPOHHT110/rg1**10+(DDAMPOHHTR2_10(RG1,10,R1,DRG1R2)*RG110-
     1      10.0D0*DRG1R2*RG19*DAMPOHHT110)/RG120*C10OHHT11*FABC)-
     1      ((DC6HOHPRAA_10(r2,cost2,DCOST2R2)*FCAB+DFCABR2*C6HOHP12)*
     1      DAMPHOHP26/rg26+(DDAMPHOHPRAA_10(RG2,6,R2,DRG2R2)
     1      *RG26-6.0D0*
     1      DRG2R2*RG25*DAMPHOHP26)/RG212*C6HOHP12*FCAB)-
     1      ((DC8HOHPRAA_10(r2,cost2,DCOST2R2)*FCAB+DFCABR2*C8HOHP12)*
     1      DAMPHOHP28/rg28+(DDAMPHOHPRAA_10(RG2,8,R2,DRG2R2)*
     1      RG28-8.0D0*
     1      DRG2R2*RG27*DAMPHOHP28)/RG216*C8HOHP12*FCAB)-
     1      ((DC10HOHPRAA_10(r2,cost2,DCOST2R2)*FCAB+DFCABR2*C10HOHP12)*
     1      DAMPHOHP210/rg210+(DDAMPHOHPRAA_10(RG2,10,R2,DRG2R2)*RG210-
     1      10.0D0*DRG2R2*RG29*DAMPHOHP210)/RG220*C10HOHP12*FCAB)-
     1      ((DC6HOHPRAB_10(r3,cost3,DCOST3R2)*FBAC+DFBACR2*C6HOHP13)*
     1      DAMPHOHP36/rg36+(DDAMPHOHPRAB_10(RG3,6,R3,DRG3R2)*
     1      RG36-6.0D0*
     1      DRG3R2*RG35*DAMPHOHP36)/RG312*C6HOHP13*FBAC)-
     1      ((DC8HOHPRAB_10(r3,cost3,DCOST3R2)*FBAC+DFBACR2*C8HOHP13)*
     1      DAMPHOHP38/rg38+(DDAMPHOHPRAB_10(RG3,8,R3,DRG3R2)*
     1      RG38-8.0D0*
     1      DRG3R2*RG37*DAMPHOHP38)/RG316*C8HOHP13*FBAC)-
     1      ((DC10HOHPRAB_10(r3,cost3,DCOST3R2)*FBAC+DFBACR2*C10HOHP13)*
     1      DAMPHOHP310/rg310+(DDAMPHOHPRAB_10(RG3,10,R3,DRG3R2)*RG310-
     1      10.0D0*DRG3R2*RG39*DAMPHOHP310)/RG320*C10HOHP13*FBAC)
     
CC      VED2=(FBC-1.0D0)*DISHHT_10(R1)+(FAB-1.0D0)*DISOHP_10(R2)+
CC     1                            (FAC-1.0D0)*DISOHP_10(R3)

      DVED22R2=DFBCR2*DISHHT_10(R1)+
     1        DFABR2*DISOHP_10(R2)+(FAB-1.0D0)*
     2        DDISPOH_10(R2,COHP(6),COHP(8), COHP(10),R0OHP,RMOHP)+
     2        DFACR2*DISOHP_10(R3)
     
      DVED32R2_10=DVED32R2_10+DVED22R2

      RETURN
      END

      FUNCTION DVED32R3_10(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/,/CVED_10/,/CVED32_10/)

c     ******************************************************************
      DVED32R3_10=-((DC6OHHTR3_10(r1,cost1,DCOST1R3)*FABC+DFABCR3
     1      *C6OHHT11)*
     1      DAMPOHHT16/rg16+(DDAMPOHHTR3_10(RG1,6,R1,DRG1R3)*RG16-6.0D0*
     1      DRG1R3*RG15*DAMPOHHT16)/RG112*C6OHHT11*FABC)-
     1      ((DC8OHHTR3_10(r1,cost1,DCOST1R3)*FABC+DFABCR3*C8OHHT11)*
     1      DAMPOHHT18/rg18+(DDAMPOHHTR3_10(RG1,8,R1,DRG1R3)*RG18-8.0D0*
     1      DRG1R3*RG17*DAMPOHHT18)/RG116*C8OHHT11*FABC)-
     1      ((DC10OHHTR3_10(r1,cost1,DCOST1R3)*FABC+DFABCR3*C10OHHT11)*
     1      DAMPOHHT110/rg110+(DDAMPOHHTR3_10(RG1,10,R1,DRG1R3)*RG110-
     1      10.0D0*DRG1R3*RG19*DAMPOHHT110)/RG120*C10OHHT11*FABC)-
     1      ((DC6HOHPRAB_10(r2,cost2,DCOST2R3)*FCAB+DFCABR3*C6HOHP12)*
     1      DAMPHOHP26/rg26+(DDAMPHOHPRAB_10(RG2,6,R2,DRG2R3)*
     1      RG26-6.0D0*
     1      DRG2R3*RG25*DAMPHOHP26)/RG212*C6HOHP12*FCAB)-
     1      ((DC8HOHPRAB_10(r2,cost2,DCOST2R3)*FCAB+DFCABR3*C8HOHP12)*
     1      DAMPHOHP28/rg28+(DDAMPHOHPRAB_10(RG2,8,R2,DRG2R3)*
     1      RG28-8.0D0*
     1      DRG2R3*RG27*DAMPHOHP28)/RG216*C8HOHP12*FCAB)-
     1      ((DC10HOHPRAB_10(r2,cost2,DCOST2R3)*FCAB+DFCABR3*C10HOHP12)*
     1      DAMPHOHP210/rg210+(DDAMPHOHPRAB_10(RG2,10,R2,DRG2R3)*RG210-
     1      10.0D0*DRG2R3*RG29*DAMPHOHP210)/RG220*C10HOHP12*FCAB)-
     1      ((DC6HOHPRAA_10(r3,cost3,DCOST3R3)*FBAC+DFBACR3*C6HOHP13)*
     1      DAMPHOHP36/rg36+(DDAMPHOHPRAA_10(RG3,6,R3,DRG3R3)*
     1      RG36-6.0D0*
     1      DRG3R3*RG35*DAMPHOHP36)/RG312*C6HOHP13*FBAC)-
     1      ((DC8HOHPRAA_10(r3,cost3,DCOST3R3)*FBAC+DFBACR3*C8HOHP13)*
     1      DAMPHOHP38/rg38+(DDAMPHOHPRAA_10(RG3,8,R3,DRG3R3)*
     1      RG38-8.0D0*
     1      DRG3R3*RG37*DAMPHOHP38)/RG316*C8HOHP13*FBAC)-
     1      ((DC10HOHPRAA_10(r3,cost3,DCOST3R3)*FBAC+DFBACR3*C10HOHP13)*
     1      DAMPHOHP310/rg310+(DDAMPHOHPRAA_10(RG3,10,R3,DRG3R3)*RG310-
     1      10.0D0*DRG3R3*RG39*DAMPHOHP310)/RG320*C10HOHP13*FBAC)

     
CC      VED2=(FBC-1.0D0)*DISHHT_10(R1)+(FAB-1.0D0)*DISOHP_10(R2)+
CC     1                            (FAC-1.0D0)*DISOHP_10(R3)

      DVED22R3=DFBCR3*DISHHT_10(R1)+
     1        DFABR3*DISOHP_10(R2)+
     2        DFACR3*DISOHP_10(R3)+(FAC-1.0D0)*
     2        DDISPOH_10(R3,COHP(6),COHP(8),COHP(10),R0OHP,RMOHP)
     
      DVED32R3_10=DVED32R3_10+DVED22R3
      RETURN
      END
      
            
      FUNCTION DAMPHOHP_10(R,N,X)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
c     calculado a partir das dist. de eq.(POL94:7651) 
      RMHOHP=2.216766092d0
      RO=0.5D0*(RMHOHP+2.5D0*R0HOHPC_10(x))
      DAMPHOHP_10=(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO**2))**N      
      RETURN
      END 
      
      FUNCTION DDAMPHOHPR1_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMHOHP=2.216766092d0
      RO=0.5D0*(RMHOHP+2.5D0*R0HOHPC_10(X))
      RO2=RO*RO
      DROR1=0.0D0
      DDAMPHOHPR1_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR1*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR1*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPHOHPRAA_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMHOHP=2.216766092d0
      RO=0.5D0*(RMHOHP+2.5D0*R0HOHPC_10(x))
      RO2=RO*RO
      DROR23=0.5D0*2.5D0*DR0HOHPCR23_10(x)
      DDAMPHOHPRAA_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR23*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR23*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPHOHPRAB_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMHOHP=2.216766092d0
      RO=0.5D0*(RMHOHP+2.5D0*R0HOHPC_10(x))
      RO2=RO*RO
      DROR23=0.0D0
      DDAMPHOHPRAB_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR23*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR23*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
             
      FUNCTION DAMPOHHT_10(R,N,X)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMOHHT=1.1077763d0
      RO=0.5D0*(RMOHHT+2.5D0*R0OHHTC_10(X))
      DAMPOHHT_10=(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO**2))**N
      RETURN
      END 
      
      FUNCTION DDAMPOHHTR1_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMOHHT=1.1077763d0
      RO=0.5D0*(RMOHHT+2.5D0*R0OHHTC_10(X))
      RO2=RO*RO
      DROR1=0.5D0*2.5D0*DR0OHHTCR1_10(X)
      DDAMPOHHTR1_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR1*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR1*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPOHHTR2_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMOHHT=1.1077763d0
      RO=0.5D0*(RMOHHT+2.5D0*R0OHHTC_10(X))
      RO2=RO*RO
      DROR2=0.0D0
      DDAMPOHHTR2_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR2*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR2*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPOHHTR3_10(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c     ******************************************************************
      RMOHHT=1.1077763d0
      RO=0.5D0*(RMOHHT+2.5D0*R0OHHTC_10(X))
      RO2=RO*RO
      DROR3=0.0D0
      DDAMPOHHTR3_10=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1.0d0)
     1 *(A(N)*(D*RO-DROR3*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR3*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
                       
      FUNCTION C6OHHT_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c6ohhtae_10/a6t,b6t
!$OMP THREADPRIVATE(/c6ohhtae_10/)
c     ******************************************************************
      A6t=C6OHHTPA_10(R)
      B6t=C6OHHTPE_10(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C6OHHT_10=(1.0D0/3.0D0)*(2.0D0*B6t+A6t)+(2.0D0/3.0D0)*(A6t-B6t)*PL
      RETURN
      END
      
      FUNCTION DC6OHHTR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c6ohhtae_10/a6t,b6t
!$OMP THREADPRIVATE(/c6ohhtae_10/)
c     ******************************************************************
      DC6PE=DC6OHHTPER1_10(R)
      DC6PA=DC6OHHTPAR1_10(R)
      DC6OHHTR1_10=2.0D0/3.0D0*DC6PE+1.0D0/3.0D0*DC6PA
     1        +1.0D0/3.0D0*(DC6PA-DC6PE)*(3.0D0*
     1         COST**2-1.0D0)+6.0D0*D*COST*1.0D0/3.0D0*(A6t-B6t)
      RETURN
      END
      
      FUNCTION DC6OHHTR2_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c6ohhtae_10/a6t,b6t
!$OMP THREADPRIVATE(/c6ohhtae_10/)
c     ******************************************************************
      dummy=r
      DC6OHHTR2_10=6.0D0*D*COST*1.0D0/3.0D0*(A6t-B6t)
      RETURN
      END      
            
      FUNCTION DC6OHHTR3_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c6ohhtae_10/a6t,b6t
!$OMP THREADPRIVATE(/c6ohhtae_10/)
c     ******************************************************************
      dummy=r 
      DC6OHHTR3_10=6.0D0*D*COST*1.0D0/3.0D0*(A6t-B6t)
      RETURN
      END      
      
      FUNCTION C8OHHT_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c8ohhtae_10/a8t,b8t
!$OMP THREADPRIVATE(/c8ohhtae_10/)
c     ******************************************************************
      A8t=C8OHHTPA_10(R)
      B8t=C8OHHTPE_10(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C8OHHT_10=(1.0D0/3.0D0)*(2.0D0*B8t+A8t)+(2.0D0/3.0D0)*(A8t-B8t)*PL
      RETURN
      END
      
      FUNCTION DC8OHHTR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c8ohhtae_10/a8t,b8t
!$OMP THREADPRIVATE(/c8ohhtae_10/)
c     ******************************************************************
      DC8PE=DC8OHHTPER1_10(R)
      DC8PA=DC8OHHTPAR1_10(R)      
      DC8OHHTR1_10=2.0D0/3.0D0*DC8PE+1.0D0/3.0D0*DC8PA
     1        +1.0D0/3.0D0*(DC8PA-DC8PE)*(3.0D0*
     1         COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A8t-B8t)
      RETURN
      END
      
      FUNCTION DC8OHHTR2_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c8ohhtae_10/a8t,b8t
!$OMP THREADPRIVATE(/c8ohhtae_10/)
c     ******************************************************************
      dummy=r
      DC8OHHTR2_10=6.0D0*D*COST*(1.0D0/3.0D0)*(A8t-B8t)
      RETURN
      END
      
      FUNCTION DC8OHHTR3_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c8ohhtae_10/a8t,b8t
!$OMP THREADPRIVATE(/c8ohhtae_10/)
      dummy=r
c     ******************************************************************
      DC8OHHTR3_10=6.0D0*D*COST*(1.0D0/3.0D0)*(A8t-B8t)
      RETURN
      END
      
      FUNCTION C10OHHT_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhtae_10/a1t,b1t
!$OMP THREADPRIVATE(/c1ohhtae_10/)
c     ******************************************************************
      A1t=C10OHHTPA_10(R)
      B1t=C10OHHTPE_10(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C10OHHT_10=(1.0D0/3.0D0)*(2.0D0*B1t+A1t)+(2.0D0/3.0D0)*
     1 (A1t-B1t)*PL
      RETURN
      END
      
      FUNCTION DC10OHHTR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhtae_10/a1t,b1t
!$OMP THREADPRIVATE(/c1ohhtae_10/)

c     ******************************************************************
      DC10PE=DC10OHHTPER1_10(R)
      DC10PA=DC10OHHTPAR1_10(R)      
      DC10OHHTR1_10=(2.0D0/3.0D0)*DC10PE+(1.0D0/3.0D0)*
     1  DC10PA+(1.0D0/3.0D0)*(DC10PA-DC10PE)
     1   *(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A1t-B1t)
      RETURN
      END
      
      FUNCTION DC10OHHTR2_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhtae_10/a1t,b1t
!$OMP THREADPRIVATE(/c1ohhtae_10/)
      dummy=r
c     ******************************************************************
      DC10OHHTR2_10=2.0D0*D*COST*(A1t-B1t)
      RETURN
      END
      
      FUNCTION DC10OHHTR3_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhtae_10/a1t,b1t
!$OMP THREADPRIVATE(/c1ohhtae_10/)
      dummy=r
c     ******************************************************************
      DC10OHHTR3_10=2.0D0*D*COST*(A1t-B1t)
      RETURN
      END
      
      FUNCTION C6HOHP_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap=C6HOHPPA_10(R)
      Bp=C6HOHPPE_10(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C6HOHP_10=(1.0D0/3.0D0)*(2.0D0*Bp+Ap)+(2.0D0/3.0D0)*(Ap-Bp)*PL1
      RETURN
      END
      
      FUNCTION DC6HOHPR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap=C6HOHPPA_10(R)
      Bp=C6HOHPPE_10(R)
      DC6HOHPR1_10=2.0D0*D*COST*(Ap-Bp)    
      RETURN
      END
      
      FUNCTION DC6HOHPRAA_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap=C6HOHPPA_10(R)
      Bp=C6HOHPPE_10(R)
      DC6PEA=DC6HOHPPERAA_10(R)
      DC6PAA=DC6HOHPPARAA_10(R)      
      DC6HOHPRAA_10=(2.0D0/3.0D0)*DC6PEA+(1.0D0/3.0D0)*
     1   DC6PAA+(1.0D0/3.0D0)*(DC6PAA-DC6PEA
     1   )*(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(Ap-Bp)     
      RETURN
      END
      
      FUNCTION DC6HOHPRAB_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap=C6HOHPPA_10(R)
      Bp=C6HOHPPE_10(R)
      DC6HOHPRAB_10=2.0D0*D*COST*(Ap-Bp)     
      RETURN
      END
      
      FUNCTION C8HOHP_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      Ap8=C8HOHPPA_10(R)
      Bp8=C8HOHPPE_10(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C8HOHP_10=(1.0D0/3.0D0)*(2.0D0*Bp8+Ap8)+(2.0D0/3.0D0)*
     1 (Ap8-Bp8)*PL1
      RETURN
      END
      
      FUNCTION DC8HOHPR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      Ap8=C8HOHPPA_10(R)
      Bp8=C8HOHPPE_10(R)
      DC8HOHPR1_10=2.0D0*D*COST*(Ap8-Bp8)
      RETURN
      END
          
      FUNCTION DC8HOHPRAA_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      Ap8=C8HOHPPA_10(R)
      Bp8=C8HOHPPE_10(R)
      DC8PEA=DC8HOHPPERAA_10(R)
      DC8PAA=DC8HOHPPARAA_10(R)      
      DC8HOHPRAA_10=(2.0D0/3.0D0)*DC8PEA+(1.0D0/3.0D0)*
     1   DC8PAA+(1.0D0/3.0D0)*(DC8PAA-DC8PEA)
     1   *(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(Ap8-Bp8)  
   
      RETURN
      END
      
      FUNCTION DC8HOHPRAB_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      Ap8=C8HOHPPA_10(R)
      Bp8=C8HOHPPE_10(R)
      DC8HOHPRAB_10=2.0D0*D*COST*(Ap8-Bp8)     
      RETURN
      END
           
      FUNCTION C10HOHP_10(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap10=C10HOHPPA_10(R)
      Bp10=C10HOHPPE_10(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C10HOHP_10=(1.0D0/3.0D0)*(2.0D0*Bp10+Ap10)+
     1 (2.0D0/3.0D0)*(Ap10-Bp10)
     1         *PL1
      RETURN
      END
      
      FUNCTION DC10HOHPR1_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      A1=C10HOHPPA_10(R)
      B1=C10HOHPPE_10(R)
      DC10HOHPR1_10=2.0D0*D*COST*(A1-B1)
      RETURN
      END
      
      FUNCTION DC10HOHPRAA_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C10HOHPPA_10(R)
      B1=C10HOHPPE_10(R)
      DC10PEA=DC10HOHPPERAA_10(R)
      DC10PAA=DC10HOHPPARaa_10(R)      
      DC10HOHPRAA_10=(2.0D0/3.0D0)*DC10PEA+(1.0D0/3.0D0)*
     1   DC10PAA+(1.0D0/3.0D0)*(DC10PAA-DC10PEA
     1   )*(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*
     1   (A1-B1)     
      RETURN
      END
      
      FUNCTION DC10HOHPRAB_10(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C10HOHPPA_10(R)
      B1=C10HOHPPE_10(R)
      DC10HOHPRAB_10=2.0D0*D*COST*(A1-B1)     
      RETURN
      END
      
      FUNCTION C6OHHTPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6HTA_10/C6HTPA,ALP,RN,R0HTC
!$OMP THREADPRIVATE(/C6HTA_10/)

c     ******************************************************************
      ALP=ALPHHTPA_10(R)
      RN=RNHHT_10(R)
      R0HTC=R0OHHTC_10(R)
      C6HTPA=1.5D0*5.393D0*ALP/(SQRT(ALP/RN)+
     1         SQRT(5.393D0/2.79D0))
      C6OHHTPA_10=C6HTPA    
      RETURN
      END
      
      FUNCTION DC6OHHTPAR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6HTA_10/C6HTPA,ALP,RN,R0HTC
      COMMON/DC6TA_10/DC6HTA,DR0HT1
!$OMP THREADPRIVATE(/C6HTA_10/,/DC6TA_10/)
c     ******************************************************************
      DR0HT1=DR0OHHTCR1_10(R)
      const1=1.5D0*5.393D0
      dalp=DALPHHTPAR1_10(r)       
      e=0.5d0/SQRT(ALP/RN)*(dalp*rn-DRNHHTR1_10(r)*alp)/rn**2
      d=SQRT(ALP/RN)+SQRT(5.393D0/2.79D0)
      DC6HTa=(const1*dalp*d-const1*alp*e)/d**2
      DC6OHHTPaR1_10=DC6HTa   
      RETURN
      END
      
      FUNCTION C8OHHTPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6HTA_10/C6HTPA,ALP,RN,R0HTC
!$OMP THREADPRIVATE(/C6HTA_10/)
      dummy=r
c     ******************************************************************      
c     substitui 1.54 por 1.57243
      C8OHHTPA_10=1.0D0*C6HTPA*R0HTC**(1.57243D0)      
      RETURN
      END
      
      FUNCTION DC8OHHTPAR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6HTA_10/C6HTPA,ALP,RN,R0HTC
      COMMON/DC6TA_10/DC6HTA,DR0HT1
!$OMP THREADPRIVATE(/C6HTA_10/,/DC6TA_10/)
      dummy=r
c     ******************************************************************      
      DC8OHHTPAR1_10=DC6HTA*R0HTC**(1.57243D0)+1.57243D0*
     1      R0HTC**(0.57243D0)*DR0HT1*C6HTPA
      RETURN
      END
      
       FUNCTION C10OHHTPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6HTA_10/C6HTPA,ALP,RN,R0HTC
!$OMP THREADPRIVATE(/C6HTA_10/)
      dummy=r
c     ******************************************************************
      C10OHHTPA_10=1.13178D0*C6HTPA*R0HTC**(1.57243D0*2.0D0)      
      RETURN
      END
      
      FUNCTION DC10OHHTPAR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6HTA_10/C6HTPA,ALP,RN,R0HTC
      COMMON/DC6TA_10/DC6HTA,DR0HT1
!$OMP THREADPRIVATE(/C6HTA_10/,/DC6TA_10/)
      dummy=r
c     ******************************************************************
      DC10OHHTPAR1_10=1.13178D0*DC6HTA*R0HTC**(1.57243D0*
     1     2.0D0)+(1.57243D0*2.0D0)*R0HTC**(1.57243D0*2.0D0-1.0D0)
     1     *DR0HT1*1.13178D0*C6HTPA
      RETURN
      END

      FUNCTION C6OHHTPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6TPE_10/C6HTPE,ALP,RN,R0HTC
!$OMP THREADPRIVATE(/C6TPE_10/)

c     ******************************************************************
      ALP=ALPHHTPE_10(R)
      RN=RNHHT_10(R)
      R0HTC=R0OHHTC_10(R)
      C6HTPE=1.5D0*5.393D0*ALP/(SQRT(ALP/RN)+
     1         SQRT(5.393D0/2.79D0))
      C6OHHTPE_10=C6HTPE    
      RETURN
      END
      
      FUNCTION DC6OHHTPER1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6TPE_10/C6HTPE,ALP,RN,R0HTC
      COMMON/DC6TPE_10/DC6HTE1,DR0HTC1
!$OMP THREADPRIVATE(/C6TPE_10/,/DC6TPE_10/)

c     ******************************************************************
      const1=1.5D0*5.393D0
      dalp=DALPHHTPER1_10(r) 
      DR0HTC1=DR0OHHTCR1_10(R)     
      e=0.5d0/SQRT(ALP/RN)*(dalp*rn-DRNHHTR1_10(r)*alp)/rn**2
      d=SQRT(ALP/RN)+SQRT(5.393D0/2.79D0)
      DC6HTE1=(const1*dalp*d-const1*alp*e)/d**2      
      DC6OHHTPER1_10=DC6HTE1    
      RETURN
      END
       
      FUNCTION C8OHHTPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6TPE_10/C6HTPE,ALP,RN,R0HTC
!$OMP THREADPRIVATE(/C6TPE_10/)
      dummy=r
c     ******************************************************************
      C8OHHTPE_10=1.0D0*C6HTPE*R0HTC**(1.57243D0)       
      RETURN
      END
      
      FUNCTION DC8OHHTPER1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6TPE_10/C6HTPE,ALP,RN,R0HTC
      COMMON/DC6TPE_10/DC6HTE1,DR0HTC1
!$OMP THREADPRIVATE(/C6TPE_10/,/DC6TPE_10/)
      dummy=r
c     ******************************************************************
      DC8OHHTPER1_10=DC6HTE1*R0HTC**(1.57243D0)+1.57243D0*
     1      R0HTC**(1.57243D0-1.0D0)*DR0HTC1*C6HTPE
      RETURN
      END
                 
      FUNCTION C10OHHTPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6TPE_10/C6HTPE,ALP,RN,R0HTC
!$OMP THREADPRIVATE(/C6TPE_10/)
      dummy=r
c     ******************************************************************
      C10OHHTPE_10=1.13178D0*C6HTPE*R0HTC**(1.57243D0*2.0D0)      
      RETURN
      END
      
      FUNCTION DC10OHHTPER1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6TPE_10/C6HTPE,ALP,RN,R0HTC
      COMMON/DC6TPE_10/DC6HTE1,DR0HTC1
!$OMP THREADPRIVATE(/C6TPE_10/,/DC6TPE_10/)
      dummy=r
c     ******************************************************************
      DC10OHHTPER1_10=1.13178D0*DC6HTE1*R0HTC**(1.57243D0*
     1     2.0D0)+(1.57243D0*2.0D0)*R0HTC**(1.57243D0*2.0D0-1.0D0)*
     1     DR0HTC1*1.13178D0*C6HTPE
      RETURN
      END
      
      FUNCTION C6HOHPPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      ALP=ALPOHPPA_10(R)
      C6HOHPPA_10=1.5D0*4.500D0*ALP/(SQRT(ALP/RNOHP_10(R))+2.25D0)
      RETURN
      END
      
      FUNCTION DC6HOHPPARAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      const1=1.5D0*4.500D0
      ALP=ALPOHPPA_10(R)
      RN=RNOHP_10(R)
      dalp=DALPOHPPAR23_10(r)      
      d=SQRT(ALP/RN)+2.25D0
      e=0.5d0/SQRT(ALP/RN)*(dalp*rn-DRNOHPR23_10(r)*alp)/rn**2
      DC6HOHPPaRAA_10=(const1*dalp*d-e*const1*alp)/d**2
      RETURN
      END

      FUNCTION C8HOHPPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      C8HOHPPA_10=1.0D0*C6HOHPPA_10(R)*R0HOHPC_10(R)**(1.57243D0)      
      RETURN
      END
      
      FUNCTION DC8HOHPPARAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      R0H=R0HOHPC_10(R)      
      DC8HOHPPARAA_10=DC6HOHPPARAA_10(R)*R0H**1.57243D0+1.57243D0*
     1      R0H**(1.57243D0-1.0D0)*DR0HOHPCR23_10(R)*C6HOHPPA_10(R)
      RETURN
      END

      FUNCTION C10HOHPPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      C10HOHPPA_10=1.13178D0*C6HOHPPA_10(R)*R0HOHPC10_10(R)**(2.0D0*
     1 1.57243D0)      
      RETURN
      END
      
      FUNCTION DC10HOHPPARaa_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      R0H=R0HOHPC10_10(R)      
      DC10HOHPPARAA_10=1.13178D0*DC6HOHPPARAA_10(R)*R0H**(2.0D0*
     1  1.57243D0)+2.0D0*1.57243D0*R0H**(2.0D0*1.57243D0-
     1    1.0D0)*DR0HOHPC10R23_10(R)*1.13178D0*C6HOHPPA_10(R)
      RETURN
      END

      FUNCTION C6HOHPPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      ALP=ALPOHPPE_10(R)
      C6HOHPPE_10=1.5D0*4.500D0*ALP/(SQRT(ALP/RNOHP_10(R))
     1           +2.25D0)
      RETURN
      END
      
      FUNCTION DC6HOHPPERAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      const1=1.5D0*4.500D0
      ALP=ALPOHPPE_10(R)
      RNO=RNOHP_10(R)
      dalp=DALPOHPPER23_10(r)      
      d=SQRT(ALP/RNO)+2.25D0
      e=0.5d0/SQRT(ALP/RNO)*(dalp*rno-DRNOHPR23_10(r)*alp)/rno**2
      DC6HOHPPERAA_10=(const1*dalp*d-e*const1*alp)/d**2
      RETURN
      END
      
      FUNCTION C8HOHPPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      C8HOHPPE_10=1.0D0*C6HOHPPE_10(R)*R0HOHPC_10(R)**(1.57243D0)      
      RETURN
      END
      
      FUNCTION DC8HOHPPERAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      R0H=R0HOHPC_10(R)
      DC8HOHPPERAA_10=DC6HOHPPERAA_10(R)*R0H**1.57243D0+1.57243D0*
     1      R0H**(1.57243D0-1.0D0)*DR0HOHPCR23_10(R)*C6HOHPPE_10(R)
      RETURN
      END
      

      FUNCTION C10HOHPPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************  
      C10HOHPPE_10=1.13178d0*C6HOHPPE_10(R)*
     1 R0HOHPC10_10(R)**(1.57243D0*2.0D0)      
      RETURN
      END
      
      FUNCTION DC10HOHPPERAA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      R0H=R0HOHPC10_10(R)      
      DC10HOHPPERAA_10=1.13178D0*DC6HOHPPERAA_10(R)*R0H**(2.0D0*
     1  1.57243D0)+2.0D0*1.57243D0*R0H**(2.0D0*1.57243D0-
     1    1.0D0)*DR0HOHPC10R23_10(R)*1.13178D0*C6HOHPPE_10(R)
      RETURN
      END
 
      FUNCTION RNOHP_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************      
      Ainf=3.52d0
      re=1.8344d0
      a=1.04d0
      b=0.584820271D0
      c=0.584820271D0/1.04D0   
      RNOHP_10=ainf+(a+b*(R-re))*exp(-c*(R-re))
      RETURN
      END 
      
      FUNCTION DRNOHPR23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************      
      Ainf=3.52d0
      re=1.8344d0
      a=1.04d0
      b=0.584820271D0
      c=0.584820271D0/1.04D0   
      DRNOHPR23_10=b*exp(-c*(R-re))-c*(a+b*(R-re))*exp(-c*(R-re))
      RETURN
      END 

      FUNCTION RNHHT_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      c=0.60634691d0
      Ainf=1.78d0
      re=1.449d0
      a=-0.1478d0
      RNHHT_10=ainf+a*exp(-c*(R-re))
      RETURN
      END 
      
      FUNCTION DRNHHTR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      c=0.60634691d0
      Ainf=1.78d0
      re=1.449d0
      a=-0.1478d0
      DRNHHTR1_10=-a*C*exp(-c*(R-re))
      RETURN
      END 

      FUNCTION ALPHHTPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7,C8,C9/9.0D0,9.2091036d0,
     1     3.7378246d3,4.93109d0,-7.2538212d0,4.267343d0,3.8970194d-2,
     2     1.3447009d1,2.4769383d-2,8.9922858d2/     
c     ******************************************************************   
      ALPHHTPA_10=Ainf1+(C1+C2*r)*EXP(-C3*r-C4*r**2
     1    -C5*r**3)+(1-EXP(-C6*r**5))*C7/r**3
     2     +(1-EXP(-C8*r**8))*C9/r**6
      RETURN
      END 
      
      FUNCTION DALPHHTPAR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA C1,C2,C3,C4,C5,C6,C7,C8,C9/9.2091036d0,
     1     3.7378246d3,4.93109d0,-7.2538212d0,4.267343d0,3.8970194d-2,
     2     1.3447009d1,2.4769383d-2,8.9922858d2/     
      dummy=r
c     ******************************************************************
      RR2=R*R
      RR3=RR2*R
      RR4=RR3*R
      RR5=RR4*R
      RR6=RR5*R
      RR7=RR6*R
      RR8=RR7*R
      DALPHHTPAR1_10=C2*EXP(-C3*r-C4*RR2-C5*RR3)+(C1+C2*R)*
     1    (-C3-2.0D0*C4
     1    *R-3.0D0*C5*RR2)*EXP(-C3*r-C4*RR2-C5*RR3)+C7/RR3*5.0D0*C6
     1      *RR4*EXP(-C6*RR5)-3.0D0*C7/RR4*(1-EXP(-C6*RR5))+C9/RR6
     1     *8.0D0*C8*RR7*EXP(-C8*RR8)-6.0D0*C9/RR7*(1-EXP(-C8*RR8))
      RETURN
      END 

      FUNCTION ALPHHTPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7/9.2090967d0,7.7638102d1,
     1 1.6964103d0,4.3441020d-9,2.7534699d-1,4.7725385d-3,-5.8164722d1/
c     ******************************************************************
      ALPHHTPE_10=Ainf1+(C1+C2*r)*EXP(-C3*r-C4*r**2
     1    -C5*r**3)+(1-EXP(-C6*r**5))*C7/r**3
      RETURN
      END 
      
      FUNCTION DALPHHTPER1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7/9.2090967d0,7.7638102d1,
     1 1.6964103d0,4.3441020d-9,2.7534699d-1,4.7725385d-3,-5.8164722d1/
      dummy=r
c     ******************************************************************
      RR2=R*R
      RR3=RR2*R
      RR4=RR3*R
      RR5=RR4*R
      DALPHHTPER1_10=C2*EXP(-C3*r-C4*RR2-C5*RR3)+(C1+C2*r)*(-C3-2.0D0*C4
     1     *r-3.0D0*C5*RR2)*EXP(-C3*r-C4*RR2-C5*RR3)-3.0D0*C7/RR4+
     1         3.0D0*C7/RR4*EXP(-C6*RR5)+5.0D0*C6*C7*R*EXP(-C6*RR5)
      RETURN
      END 

      FUNCTION ALPOHPPA_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7,C8/9.3381D0,5.2957323d0,-1.36948
     1  90d0,1.1344069d-1,-7.3074372d-1,-2.8938141d-1,7.298064d-2,4.27
     2  10027d-3,1.2016556d2/  
c     ******************************************************************
      Azero=3.13939d0-Ainf1
      rr2=r*r
      rr3=rr2*r
      rr5=rr2*rr3
      ALPOHPPA_10=Ainf1+(Azero+C1*r+C2*rr2+C3*rr3)*
     1  EXP(-C4*r-C5*rr2-C6*rr3)+(1-EXP(-C7*rr5))*C8/rr3
      RETURN
      END 
      
      FUNCTION DALPOHPPAR23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7,C8/9.3381D0,5.2957323d0,-1.36948
     1  90d0,1.1344069d-1,-7.3074372d-1,-2.8938141d-1,7.298064d-2,4.27
     2  10027d-3,1.2016556d2/  
c     ******************************************************************
      Azero=3.13939d0-Ainf1
      rr2=r*r
      rr3=rr2*r
      rr4=rr3*r
      rr5=rr4*r      
      DALPOHPPAR23_10=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*R-C5*
     1         RR2-C6*RR3
     1         )-(C4+2.0D0*C5*R+3.0D0*C6*RR2)*(AZERO+C1*R+C2*RR2+C3*
     1         RR3)*EXP(-C4*r-C5*RR2-C6*RR3)-3.0D0*C8/RR4
     1         +3.0D0*C8/RR4*EXP(-C7*RR5)+5.0D0*C7*C8*R*EXP(-C7*RR5)
      RETURN
      END 
      
      FUNCTION ALPOHPPE_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5/9.6061296d0,4.8648692d-1,-1.3199234d-1,
     1  4.2427009d-2,2.1725130d-4,-2.5535932d1/
c     ******************************************************************
      Azero=ainf1-3.5396940d0
      ALPOHPPE_10=Ainf1-Azero*EXP(-C1*r-C2*r**2-C3*r**3)+
     1                 (1-EXP(-C4*r**5))  *C5/r**3
      RETURN
      END 

      FUNCTION DALPOHPPER23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5/9.6061296d0,4.8648692d-1,-1.3199234d-1,
     1  4.2427009d-2,2.1725130d-4,-2.5535932d1/
c     ******************************************************************
      Azero=ainf1-3.5396940d0
      rr2=r*r
      rr3=rr2*r
      rr4=rr3*r
      rr5=rr4*r      
      DALPOHPPER23_10=Azero*(C1+2.0D0*C2*r+3.0D0*C3*RR2)*EXP(-C1*
     1     r-C2*RR2 -C3*RR3)-3.0D0*c5/RR4+3.0D0*C5/RR4*EXP(-C4*RR5)
     1         +5.0D0*C4*C5*R*EXP(-C4*RR5)
      RETURN
      END 
      
      FUNCTION R0OHHTC_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************       
      finf=6.334299622d0/6.9898852d0
      fzero=7.971192026d0/8.0925471d0
      d=fzero-finf
      c=0.60634691d0
      fcorr=finf+d*exp(-c*r)      
      R0OHHTC_10=2.0d0*((ALPHHTPE_10(R)**(1.0d0/3.0d0))+
     1       sqrt(2.003423d0)) *fcorr
      RETURN
      END 
      
      FUNCTION DR0OHHTCR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************       
      finf=6.334299622d0/6.9898852d0
      fzero=7.971192026d0/8.0925471d0
      d=fzero-finf
      c=0.60634691d0
      fcorr=finf+d*exp(-c*r)      
      ALP=ALPHHTPE_10(R)
      DR0OHHTCR1_10=2.0D0/3.0D0*DALPHHTPER1_10(R)*ALP**(-2.0d0/3.0d0)*
     1       FCORR-2.0D0*D*C*exp(-c*r)*(ALP**(1.0d0/3.0d0)+
     1       sqrt(2.003423d0))
      RETURN
      END 
          
      FUNCTION R0HOHPC_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      alpham=(1.0d0/3.0d0)*ALPOHPPA_10(R)+(2.0d0/3.0d0)*ALPOHPPE_10(R)    
      finf=6.390259718d0/7.7028477d0
      fzero=5.950512078d0/6.4733336d0
      d=fzero-finf
      c=0.584820271D0/1.04D0         
      fcorr=finf+d*exp(-c*r)      
      R0HOHPC_10=2.0d0*((alpham**(1.0d0/3.0d0))+(6.9282d0/4.0d0))*fcorr  
      RETURN
      END 
      
      FUNCTION DR0HOHPCR23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      alpham=1.0d0/3.0d0*ALPOHPPA_10(R)+2.0d0/3.0d0*ALPOHPPE_10(R)    
      DALPHAM=1.0d0/3.0d0*DALPOHPPAR23_10(R)+2.0d0/3.0d0*
     1  DALPOHPPER23_10(R)
      finf=6.390259718d0/7.7028477d0
      fzero=5.950512078d0/6.4733336d0
      d=fzero-finf
      c=0.584820271D0/1.04D0         
      fcorr=finf+d*exp(-c*r)      
      DR0HOHPCR23_10=2.0d0/3.0D0*alpham**(-2.0d0/3.0d0)*DALPHAM*fcorr-
     1    2.0D0*c*d*exp(-c*r)*(alpham**(1.0d0/3.0d0)+6.9282d0/4.0d0) 
      RETURN
      END 
    
      FUNCTION R0HOHPC10_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************     
      ALPHAM=(1.0D0/3.0D0)*ALPOHPPA_10(R)+(2.0D0/3.0D0)*ALPOHPPE_10(R)  
      finf=6.574599892d0/7.7028477d0
      fzero=5.950512078d0/6.4733336d0
      d=fzero-finf
      c=0.584820271D0/1.04D0   
      FCORR=finf+d*exp(-c*r)      
      R0HOHPC10_10=2.0d0*((ALPHAM**(1.0D0/3.0D0))+(6.9282D0/4.0D0))*
     1 FCORR          
      RETURN
      END 
      
      FUNCTION DR0HOHPC10R23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************     
      ALPHAM=(1.0D0/3.0D0)*ALPOHPPA_10(R)+(2.0D0/3.0D0)*ALPOHPPE_10(R)  
      DALPHAM=1.0d0/3.0d0*DALPOHPPAR23_10(R)+2.0d0/3.0d0*
     1 DALPOHPPER23_10(R)
      finf=6.574599892d0/7.7028477d0
      fzero=5.950512078d0/6.4733336d0
      d=fzero-finf
      c=0.584820271D0/1.04D0   
      FCORR=finf+d*exp(-c*r)      
      DR0HOHPC10R23_10=2.0d0/(3.0D0*alpham**(2.0d0/3.0d0))*
     1    DALPHAM*fcorr
     1    -2.0D0*c*d*(alpham**(1.0d0/3.0d0)+(6.9282d0/4.0d0))*exp(-c*r) 
      RETURN
      END 
       
      FUNCTION FDOHP_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7/1.1014506E+00,-9.2434425E-01,3.471174
     1     2E-01,-3.6259720E-01,3.6732452E-01,6.8313966E-03,2.392610
     2     1E+00/
c     ******************************************************************
      an=3.34634328d-6
      rr2=r*r
      rr3=rr2*r
      rr5=rr3*rr2      
      FDOHP_10=(an+C1*R+C2*RR2+C3*RR3)*EXP(-c4*R-C5*RR2)+(1-EXP(-C6*
     2   RR5))*C7/RR3
      return
      end
      
      FUNCTION DFDOHPR23_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7/1.1014506E+00,-9.2434425E-01,3.471174
     1     2E-01,-3.6259720E-01,3.6732452E-01,6.8313966E-03,2.392610
     2     1E+00/
c     ******************************************************************
      an=3.34634328d-6
      rr2=r*r
      rr3=rr2*r
      rr4=rr3*r
      rr5=rr4*r      
      DFDOHPR23_10=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*R-C5*RR2)-(C4+
     1         2.0D0*C5*R)*(an+C1*R+C2*RR2+C3*RR3)*EXP(-c4*R-C5*RR2)
     1         -3.0D0*C7/RR4+5.0D0*C6*C7*R*EXP(-C6*RR5)+3.0D0*
     1         C7/RR4*EXP(-C6*RR5)
      return
      end
     
      FUNCTION VINDHOHP_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
!$OMP THREADPRIVATE(/CVED_10/,/CVED32_10/)
      dummy=r1
      dummy=r2      
      dummy=r3

c     ******************************************************************               
      VINDHOHP_10=-FDOHP_10(R2)**2*4.5D0*(3.0D0*COST2**2+1.0D0)*
     1          DAMPHOHP26/(2.0D0*rG26)
     1         -FDOHP_10(R3)**2*4.5D0*(3.0D0*COST3**2+1.0D0)*
     1          DAMPHOHP36/(2.0D0*rG36)
      return
      end
      
      FUNCTION DVINDHOHPR1_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
!$OMP THREADPRIVATE(/CVED_10/,/CVED32_10/)

c     ******************************************************************
      TA=3.0D0*COST2**2+1.0D0    
      DTAR1=-6.0D0*(R3-R1)/R2**2
      TB=3.0D0*COST3**2+1.0D0    
      DTBR1=-6.0D0*(R2-R1)/R3**2
      DFDOHPR1=0.0D0
      xb=DAMPHOHP26/(2.0D0*rG26) 
      xb1=DAMPHOHP36/(2.0D0*rG36) 
      FDOHPR2=FDOHP_10(R2)
      FDOHPR3=FDOHP_10(R3)
      TA1=-4.5D0*FDOHPR2**2
      TA2=-4.5D0*FDOHPR3**2
      xa=TA1*TA  
      xa1=TA2*TB  
      DVINDHOHPR1_10=xb*(-4.5d0*2.0d0*FDOHPR2*DFDOHPR1*Ta+dTar1*TA1)
     1 +xa*(DDAMPHOHPR1_10(Rg2,6,r2,DRG2R1)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R1*DAMPHOHP26)/(4.0D0*RG212)
     1         +xb1*(-4.5d0*2.0d0*FDOHPR3*DFDOHPR1*TB+dTBr1*TA2)
     1 +Xa1*(DDAMPHOHPR1_10(Rg3,6,r3,DRG3R1)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R1*DAMPHOHP36)/(4.0D0*RG312) 
      return
      end
      
      FUNCTION DVINDHOHPR2_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
!$OMP THREADPRIVATE(/CVED_10/,/CVED32_10/)

c     ******************************************************************
      TA=3.0D0*COST2**2+1.0D0    
      DTAR2=-6.0D0*(R3-R1)**2/R2**3
      TB=3.0D0*COST3**2+1.0D0    
      DTBR2=6.0D0*(R2-R1)/R3**2
      FDOHPR2=FDOHP_10(R2)
      FDOHPR3=FDOHP_10(R3)
      xb=DAMPHOHP26/(2.0D0*rG26) 
      xb1=DAMPHOHP36/(2.0D0*rG36) 
      TA1=-4.5D0*FDOHPR2**2
      TA2=-4.5D0*FDOHPR3**2
      xa=TA1*TA  
      xa1=TA2*TB  
      DFDOHPR23R3=0.0D0
      DVINDHOHPR2_10=xb*(-4.5d0*2.0d0*FDOHPR2*DFDOHPR23_10(R2)*Ta+dTar2
     1 *TA1)+xa*(DDAMPHOHPRAA_10(Rg2,6,r2,DRG2R2)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R2*DAMPHOHP26)/(4.0D0*RG212)
     1      +xb1*(-4.5d0*2.0d0*FDOHPR3*DFDOHPR23R3*TB+dTBr2*TA2)
     1 +Xa1*(DDAMPHOHPRAB_10(Rg3,6,r3,DRG3R2)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R2*DAMPHOHP36)/(4.0D0*RG312)
      return
      end      
      
      FUNCTION DVINDHOHPR3_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED32_10/C6OHHT11,C8OHHT11,C10OHHT11,C6HOHP12,C8HOHP12,
     1    C10HOHP12,C6HOHP13,C8HOHP13,C10HOHP13,DAMPOHHT16,
     1    DAMPOHHT18,DAMPOHHT110,DAMPHOHP26,DAMPHOHP28,DAMPHOHP210,
     1    DAMPHOHP36,DAMPHOHP38,DAMPHOHP310,C6AC2,C8AC2,C10AC2,C6AB2,
     1    C8AB2,C10AB2
!$OMP THREADPRIVATE(/CVED_10/,/CVED32_10/)

c     ******************************************************************
      TA=3.0D0*COST2**2+1.0D0    
      DTAR3=6.0D0*(R3-R1)/R2**2
      TB=3.0D0*COST3**2+1.0D0    
      DTBR3=-6.0D0*(R2-R1)**2/R3**3
      FDOHPR2=FDOHP_10(R2)
      FDOHPR3=FDOHP_10(R3)
      xb=DAMPHOHP26/(2.0D0*rG26) 
      xb1=DAMPHOHP36/(2.0D0*rG36) 
      TA1=-4.5D0*FDOHPR2**2
      TA2=-4.5D0*FDOHPR3**2
      xa=TA1*TA  
      xa1=TA2*TB  
      DFDOHPR23R2=0.0D0      
      DVINDHOHPR3_10=xb*(-4.5d0*2.0d0*FDOHPR2*DFDOHPR23R2*Ta+dTar3
     1 *TA1)+xa*(DDAMPHOHPRAB_10(Rg2,6,r2,DRG2R3)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R3*DAMPHOHP26)/(4.0D0*RG212)
     1    +xb1*(-4.5d0*2.0d0*FDOHPR3*DFDOHPR23_10(R3)*TB+dTBr3*TA2)
     1 +Xa1*(DDAMPHOHPRAA_10(Rg3,6,r3,DRG3R3)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R3*DAMPHOHP36)/(4.0D0*RG312)
      return
      end      
      
      FUNCTION FQHHT_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4/1.9930192E-02,6.0772551E-01,1.2408772E-02,-2.553
     1                 2637E+01/
c     ******************************************************************
      FQHHT_10=-12.0d0*EXP(-C1*r-C2*r**2)+(1-EXP(-C3*r**5))*C4/r**3 
      return
      end
      
      FUNCTION DFQHHTR1_10(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4/1.9930192E-02,6.0772551E-01,1.2408772E-02,-2.553
     1                 2637E+01/
c     ******************************************************************
      RR2=R*R
      RR3=RR2*R
      RR4=RR2*RR2
      RR5=RR3*RR2
      DFQHHTR1_10=12.0D0*(C1+2.0D0*C2*R)*EXP(-C1*r-C2*RR2)-3.0D0*C4/RR4
     1       +3.0D0*C4/RR4*EXP(-C3*RR5)+5.0D0*C4*C3*R*EXP(-C3*RR5)
      return
      end
     
      FUNCTION VEOHHT_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      DATA b0,b2,b4,b6/-3.68666806d0,-0.0391352778d0,-0.277479722d0,
     1                 0.0402558333d0/
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
!$OMP THREADPRIVATE(/CVED_10/)
      dummy=r2      
      dummy=r3

c     ******************************************************************
      at=b0+b2*cos2t1+b4*cos4t1+b6*cos6t1
      VEOHHT_10=0.75d0*FQHHT_10(R1)*(-0.992148212d0)*at*DAMPOHHT5
      return
      end
      
      FUNCTION DVEOHHTR1_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      DATA b0,b2,b4,b6/-3.68666806d0,-0.0391352778d0,-0.277479722d0,
     1                 0.0402558333d0/
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3

!$OMP THREADPRIVATE(/CVED_10/)
      dummy=r2
      dummy=r3
c     ******************************************************************
      at=b0+b2*cos2t1+b4*cos4t1+b6*cos6t1
      DATR1=(B2*dcos2t+b4*dcos4t+b6*dcos6t)*dcost1r1
      CONST2=0.75d0*(-0.992148212d0)
      FQHHTR1=FQHHT_10(R1)
      DVEOHHTR1_10=(DFQHHTR1_10(R1)*AT+DATR1*FQHHTR1)*CONST2*DAMPOHHT5+
     1     (DDAMPOHHTR1_10(Rg1,5,r1,drg1r1)*rg15-5.0D0*DRG1R1*
     1     RG14*DAMPOHHT5*RG15)/(rg110)*at*CONST2*FQHHTR1
      return
      end
      
      FUNCTION DVEOHHTR2_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      DATA b0,b2,b4,b6/-3.68666806d0,-0.0391352778d0,-0.277479722d0,
     1                 0.0402558333d0/
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3

!$OMP THREADPRIVATE(/CVED_10/)

      dummy=r2
      dummy=r3

c     ******************************************************************
      at=b0+b2*cos2t1+b4*cos4t1+b6*cos6t1
      DATR2=(B2*dcos2t+b4*dcos4t+b6*dcos6t)*dcost1r2
      CONST2=0.75d0*(-0.992148212d0)
      DFQHHTR2=0.0D0
      FQHHTR1=FQHHT_10(R1)
      DVEOHHTR2_10=(DFQHHTR2*AT+DATR2*FQHHTR1)*CONST2*DAMPOHHT5+
     1       (DDAMPOHHTR2_10(Rg1,5,r1,DRG1R2)*rg15-5.0D0*DRG1R2*
     1        RG14*DAMPOHHT5*Rg15)/(rg110)*at*CONST2*FQHHTR1
      return
      end
      
      FUNCTION DVEOHHTR3_10(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      DATA b0,b2,b4,b6/-3.68666806d0,-0.0391352778d0,-0.277479722d0,
     1                 0.0402558333d0/
      COMMON/CVED_10/COSB,COSC,RG1,RG2,RG3,COST1,COST2,COST3,FAB,FAC,
     1  FBC,FABC,FCAB,FBAC,CONST,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHS5,DAMPOHHT5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3

!$OMP THREADPRIVATE(/CVED_10/)
      dummy=r2
      dummy=r3


c     ******************************************************************
      at=b0+b2*cos2t1+b4*cos4t1+b6*cos6t1
      DATR3=(B2*dcos2t+b4*dcos4t+b6*dcos6t)*dcost1r3
      CONST2=0.75d0*(-0.992148212d0)
      DFQHHTR3=0.0D0
      FQHHTR1=FQHHT_10(R1)
      DVEOHHTR3_10=(DFQHHTR3*AT+DATR3*FQHHTR1)*CONST2*DAMPOHHT5+
     1    (DDAMPOHHTR3_10(Rg1,5,r1,DRG1R3)*rg15-5.0D0*DRG1R3*
     1    RG14*DAMPOHHT5*RG15)/(rg110)*at*CONST2*FQHHTR1
      return
      end
     
      FUNCTION rl2_10(r1E,r2E,r3E)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      R1=R1E+1.0D-8
      R2=R2E+1.0D-8
      R3=R3E+1.0D-8
      CALL GEOMRL_10(R1,R2,R3)
c    ***************************************************************  

      rl2_10=VED32_10(r1,r2,r3)+VEOHHT_10(r1,r2,r3)+
     1 VINDHOHP_10(r1,r2,r3)
      return
      end
      
      SUBROUTINE Drl2_10(R1E,R2E,R3E,DER1,DER2,DER3)
c     ******************************************************************
c     TO COMPUTE THE DERIVATIVES OF THE LONG RANGE TERM
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      R1=R1E+1.0D-8
      R2=R2E+1.0D-8
      R3=R3E+1.0D-8
      CALL GEOMRL_10(R1,R2,R3)
c    ***************************************************************  

      DER1=DVED32R1_10(r1,r2,r3)+DVEOHHTR1_10(r1,r2,r3)+
     1 DVINDHOHPR1_10(r1,r2,r3)

      DER2=DVED32R2_10(r1,r2,r3)+DVEOHHTR2_10(r1,r2,r3)+
     1 DVINDHOHPR2_10(r1,r2,r3)

      DER3=DVED32R3_10(r1,r2,r3)+DVEOHHTR3_10(r1,r2,r3)+
     1 DVINDHOHPR3_10(r1,r2,r3)
      RETURN
      END

     
      FUNCTION DISPOH_10(R,C6,C8,C10,R0,RM)
c    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c    **************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      RR2=RR*RR
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR2))**6
      D8=(1.0D0-EXP(-A(8)*RR-B(8)*RR2))**8
      D10=(1.0D0-EXP(-A(10)*RR-B(10)*RR2))**10
      DISPOH_10=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

      FUNCTION DDISPOH_10(R,C6,C8,C10,R0,RM)
c    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
c    **************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      RR2=RR*RR
      DRR=RR/R
      T6=1.0D0-EXP(-A(6)*RR-B(6)*RR2)
      T8=1.0D0-EXP(-A(8)*RR-B(8)*RR2)
      T10=1.0D0-EXP(-A(10)*RR-B(10)*RR2)
      DDISPOH_10=6.0D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-A(6)-
     1  2.0D0*B(6)*RR)*DRR)+
     2 8.0D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-A(8)-2.0D0*B(8)
     3      *RR)*DRR)+
     4 10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-A(10)-2.0D0*
     5      B(10)*RR)*DRR)
      RETURN
      END
      
      BLOCK DATA H2OSDAT_10
C     *************************************************************      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/CONST_10/PI
      COMMON/AN_10/A(20)
      COMMON/BN_10/B(20)
      COMMON/EMI_10/EMIN
      COMMON/DISPC1_10/CHH(16),COHP(10),COHS(10)
      COMMON/REFGEO_10/R10,R20,R30
      COMMON/REFGEO2_10/R102,R202,R302
      COMMON/REFGEO12_10/R1012,R2012,R3012
      COMMON/COEFF_10/C(148)
      COMMON/DIATDI_10/R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS,VO1D
!$OMP THREADPRIVATE(/DIATDI_10/)

c    **************************************************************
      DATA PI/3.141592653589793238462643D0/
      DATA R10,R20,R30/2.86194D0,1.80965D0,1.80965D0/
      DATA R102,R202,R302/4.04d0,2.717D0,2.717D0/
      DATA r1012,r2012,r3012/3.0d0,2.61660919D0,2.61660919D0/
      DATA VO1D,EMIN/7.1955D-2,-0.3704003D0/
      DATA R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS/6.928203D0,
     1   1.40100D0,7.82D0,6.294894D0,6.334299622D0,1.8344D0,1.9086D0/   
c    **************************************************************      
c     AN E BN OBTIDOS A PARTIR DOS COEF.(A0,A1,B0,B1) DO ARTG. VAR82:857 
      DATA A/0.000000000D+00,0.000000000D+00,0.000000000D+00
     1 ,0.500798749D+01,0.384282945D+01,0.309513333D+01
     1 ,0.257767037D+01,0.219990002D+01,0.191291419D+01
     1 ,0.168807142D+01,0.150753106D+01,0.135962478D+01
     1 ,0.123641324D+01,0.113231455D+01,0.104329456D+01
     1 ,0.966368248D+00,0.899281484D+00,0.840301611D+00
     1 ,0.788075808D+00,0.741533075D+00/
      DATA B/0.000000000D+00,0.000000000D+00,0.000000000D+00
     1 ,0.106645006D+02,0.967581549D+01,0.877878947D+01
     1 ,0.796492498D+01,0.722651228D+01,0.655655639D+01
     1 ,0.594871080D+01,0.539721740D+01,0.489685187D+01
     1 ,0.444287427D+01,0.403098404D+01,0.365727935D+01
     1 ,0.331822010D+01,0.301059437D+01,0.273148802D+01
     1 ,0.247825708D+01,0.224850269D+01/
c    **************************************************************
      DATA CHH/5*0.0D0,6.499027D0,0.0d0,124.3991D0,0.0d0,3285.828D0,
     1      -3475.0D0,1.215D5,-2.914D5,6.061D6,-2.305D7,3.938D8/ 
      DATA COHS/5*0.0D0,10.34D0,0.0D0,306.8585254D0,0.0D0,
     1 6328.709221D0/
      DATA COHP/5*0.0D0,10.00D0,0.0D0,180.45D0,0.0D0,3685.26D0/
C    ***************************************************************      
      DATA C/
     1   -0.1565690230D+00,   -0.2544931663D+00,   -0.4876804350D-01,
     1   -0.2690037685D+00,   -0.1535865259D-02,   -0.4890066230D-01,
     1   -0.7329374732D-01,   -0.2109344267D+00,   -0.6445301047D-01,
     1    0.2411710888D-01,   -0.1555868493D-01,   -0.5607487072D-01,
     1    0.6202541496D-01,   -0.1429634757D+00,   -0.7444150867D-01,
     1    0.6156263566D-01,   -0.2245935938D-02,    0.1330353546D-01,
     1   -0.6788598952D-01,    0.2271341007D-02,    0.8980070809D-02,
     1   -0.4986636030D-01,   -0.8429336120D-01,    0.8978555840D-02,
     1    0.3568207396D-01,    0.7898108250D-01,   -0.3299218814D-02,
     1    0.2266503020D-01,    0.4924099395D-01,   -0.6192243909D-01,
     1    0.1006236871D-01,   -0.1067837639D+00,    0.1965447873D-01,
     1    0.1755636267D-01,   -0.4571252263D-01,    0.1272353899D-01,
     1   -0.1199479525D+00,    0.9468891942D-01,   -0.2616965983D-01,
     1   -0.6339820029D-01,   -0.4934658276D-01,    0.3225361439D-01,
     1    0.7834891306D-01,   -0.4235299197D-01,   -0.8355138053D-01,
     1   -0.2899296365D+00,    0.1331932799D+00,   -0.5167976333D-01,
     1   -0.2347859542D-01,   -0.1777711260D-01,   -0.1299842803D-01,
     1    0.2970783506D-01,   -0.7347156570D-01,    0.7465634301D-02,
     1    0.9699947247D-02,    0.1649570767D+00,    0.3481731649D-01,
     1    0.6307806127D-01,    0.5148261697D-02,   -0.6304245065D-01,
     1    0.6882370932D-01,    0.4033730863D-01,    0.2565541551D-01,
     1   -0.6078121832D-01,   -0.2962569550D-01,    0.3084390610D-01,
     1   -0.2015458742D+00,    0.6721735197D-01,    0.1712406839D+00,
     1    0.6228997124D-01,   -0.3420843299D-03,    0.1340173021D-01,
     1   -0.4185837518D-01,   -0.1143824325D-01,    0.3392045669D-01,
     1    0.8324394198D-01,   -0.8282447683D-02,   -0.7865970804D-01,
     1   -0.2481027858D-01,   -0.1354741296D-01,   -0.5307323810D-02,
     1   -0.9247430789D-02,    0.2436516685D-01,   -0.3674205193D-02,
     1    0.3402165882D-02,    0.1899518378D-01,   -0.1741461802D-01,
     1    0.1174832037D-02,   -0.5667240702D-01,   -0.7274439089D-02,
     1    0.9909591032D-01,   -0.7952480341D-01,   -0.3072103204D-02,
     1    0.1970459378D-02,   -0.2362856551D-01,    0.1208067587D+00,
     1   -0.3461573091D-01,   -0.3750302864D-01,   -0.1024072632D+00,
     1    0.3886559362D-02,   -0.5565028247D-01,   -0.1387475546D+00,
     1   -0.5669833170D-01,    0.8381138633D-01,    0.2316598891D-02,
     1   -0.4894908326D-01,   -0.8779017665D-01,    0.5333761656D-01,
     1   -0.1313796927D-01,    0.5581602062D-01,    0.8072039816D-02,
     1    0.3448517557D-01,   -0.1509775116D-01,    0.6871987343D-02,
     1    0.4548979350D-01,    0.8392592250D-03,    0.2332029125D-01,
     1   -0.1285789631D-02,   -0.1553370286D-01,   -0.1519678996D-01,
     1   -0.2608303613D-01,   -0.2439369041D-02,   -0.3230616849D-03,
     1   -0.1622479086D-01,    0.4314598580D-01,    0.3349953640D-01,
     1    0.1267395642D-01,    0.3163071652D-02,    0.5184081419D-02,
     1   -0.1058457775D+00,   -0.2059102604D-01,    0.6661037391D-01,
     1   -0.6610026832D-01,    0.1843067372D-01,    0.1103968854D+00,
     1    0.5152243701D-01,    0.7878724664D-01,    0.8313253687D-01,
     1   -0.2827314866D-01,   -0.1888672360D+00,   -0.1251485382D+00,
     1   -0.7811582863D-01,    0.7702019726D-08,    0.3652832367D+00,
     1    0.8000000000D+00,    0.1000000000D+01,    0.8000000000D+00,
     1    0.1000000000D+01/
      END
