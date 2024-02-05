       
      Function VH2OT3C_f31(r1,r2,r3)
c    **************************************************************
C     TO COMPUTE THE H2O SURFACE - Sx in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
C    **************************************************************
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
C    **************************************************************
  
      
c      DIAT=VHHST_f31(r1)+VOHPT_f31(r2)+VOHPT_f31(r3)
c      F=THREBR_f31(r1,r2,r3)+RL_f31(r1,r2,r3)+DIAT
     
      F=THREBR_f31(r1,r2,r3)+RL_f31(r1,r2,r3)
      
      VH2OT3C_f31=F
      RETURN
      END
      
      SUBROUTINE DERVH2OT_f31(r1,r2,r3,g1,g2,g3)
c    **************************************************************
C     TO COMPUTE THE DERIVATIVE H2O SURFACE - Sx in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
c      DIMENSION GF1(3),GF2(3)      
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP

c    ***************************************************************            

c     if (r3.GT.r2) then  
c       TEMP=r2
c       r2=r3
c       r3=TEMP
c      endif

      F=VHHST_f31(r1)+VOHPT_f31(r2)+VOHPT_f31(r3)+THREBR_f31(r1,r2,r3)+
     1  RL_f31(r1,r2,r3)
                 
      CALL DTHREBR_f31(r1,r2,r3,DTDR1,DTDR2,DTDR3)

      CALL DRL_f31(r1,r2,r3,DER1,DER2,DER3)
            
      G1=DTDR1+DER1+DVHHST_f31(r1)	  
      G2=DTDR2+DER2+DVOHPT_f31(r2)	  
      G3=DTDR3+DER3+DVOHPT_f31(r3)
     
      RETURN
      END
      
       SUBROUTINE DERVH2OT3C_f31(r1,r2,r3,g1,g2,g3)
c    **************************************************************
C     TO COMPUTE THE DERIVATIVE H2O SURFACE - Sx in atomic units
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)     
c      DIMENSION GF1(3),GF2(3)      
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP

c    ***************************************************************            

c     if (r3.GT.r2) then  
c       TEMP=r2
c       r2=r3
c       r3=TEMP
c      endif

                 
      CALL DTHREBR_f31(r1,r2,r3,DTDR1,DTDR2,DTDR3)

      CALL DRL_f31(r1,r2,r3,DER1,DER2,DER3)
            
      G1=DTDR1+DER1	  
      G2=DTDR2+DER2	  
      G3=DTDR3+DER3
      RETURN
      END     
     

      
            FUNCTION VOHPT_f31(R)
c    **************************************************************
C     TO COMPUTE THE HFACE FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c    **************************************************************
      VOHPT_f31=EHFOHPT_f31(R)+DISOHPT_f31(R)
      RETURN
      END

      FUNCTION EHFOHPT_f31(R)
c    **************************************************************
C     TO COMPUTE THE EHF FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
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
      EHFOHPT_f31=EHF+EX
      RETURN
      END

      FUNCTION DISOHPT_f31(R)
c    **************************************************************
C     TO COMPUTE THE DISPERSION FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/dispc1t_f31/CHH(16),COHP(10)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
c    **************************************************************
      DISOHPT_f31=DISPOHT_f31(R,COHP(6),COHP(8),COHP(10),R0OHP,RMOHP)
      RETURN
      END

      FUNCTION DVOHPT_f31(R)
c    **************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/dispc1t_f31/CHH(16),COHP(10)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c    **************************************************************
      DATA D,ASV/0.240599442D0,2.4093609D0,1.1587677D0,5.3697044D-1/
      DATA (GAM(I),i=0,2)/1.8450419D0,2.7809326D3,5.3432177D-5/
      DATA AGEX,ALPHEX,A1EX,GAMEX/0.307D0,1.5D0,2.257329D0,2.0D0/
c    **************************************************************
      X=R-RMOHP
      X2=X*X
      X3=X2*X
            
      POL=-D/R*(1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3)
      
      DPOL=-D/R*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*X2)-POL/R
      GAMA=gam(0)*(1+gam(1)*TANH(gam(2)*X))
      DGAMA=GAM(0)*GAM(1)*GAM(2)/COSH(GAM(2)*X)**2
      POT=EXP(-GAMA*X)
      DEHF=-(DGAMA*X+GAMA)*POT*POL+DPOL*POT
      RHO=(RMOHP+2.5D0*R0OHP)/2.0D0
    
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6
      DBASE=(+A(6)+2.0D0*B(6)*RR)/RHO*
     1        EXP(-A(6)*RR-B(6)*RR**2)
      DD6=6.0D0*(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**5*DBASE
      EX=-AGEX*R**ALPHEX*(1.0D0+A1EX*R)*EXP(-GAMEX*R)*D6
      DEXND=(-AGEX*ALPHEX*R**(ALPHEX-1.0D0)-AGEX*A1EX*(ALPHEX+1.0D0)
     1  *R**ALPHEX)*EXP(-GAMEX*R)+GAMEX*AGEX*R**ALPHEX*(1.0D0+A1EX*R)
     2    *EXP(-GAMEX*R)
      
      DEX=EX/D6*DD6+DEXND*D6
c     Dispersao
      
      DERDISP=DDISPOHT_f31(R,COHP(6),COHP(8),COHP(10),R0OHP,RMOHP)
      
c     Total
         
      DVOHPT_f31=DEHF+DEX+DERDISP
     
      RETURN
      END 
      
            FUNCTION VHHST_f31(R)
c    **************************************************************
C     TO COMPUTE THE HFACE FOR singlet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c    **************************************************************
      VHHST_f31=EHFHHST_f31(R)+DISHHST_f31(R)
      RETURN
      END

      FUNCTION EHFHHST_f31(R)
c    **************************************************************
C     TO COMPUTE THE EHF FOR  singlet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c    **************************************************************
      DATA D,ASV/0.218973D0,1.91479d0,0.646041D0,0.346414D0/
      DATA (GAM(I),i=0,2)/1.22349D0,1.04334D0,0.208477D0/
      DATA AGEX,ALPHEX,GAMEX/0.8205D0,2.5D0,2.0D0/
c    **************************************************************
      X=R-RMHHS
      X2=X*X
      X3=X2*X
      POL=-D/R*(1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3)
      Gama=Gam(0)*(1.0D0+Gam(1)*TANH(Gam(2)*X))      
      EHFHHST_f31=POL*EXP(-Gama*X)
      RHO=(RMHHS+2.5D0*R0HH)/2.0d0
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6  
      EXCHHS=-AGEX*R**ALPHEX*EXP(-GAMEX*R)*D6
      EHFHHST_f31=EHFHHST_f31+EXCHHS
      RETURN
      END

      FUNCTION DISHHST_f31(R)
c    **************************************************************
C     TO COMPUTE THE DISPERSION FOR singlet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/dispc1t_f31/CHH(16),COHP(10)
c    **************************************************************
      R2=R*R
      R6=R2*R2*R2
      R8=R6*R2
      ADISP=-CHH(6)*DAMPHHTTRIP_f31(6,R)/R6-CHH(8)*
     1 DAMPHHTTRIP_f31(8,R)/R8
      DO N=10,16
        ADISP=ADISP-CHH(N)*DAMPHHTTRIP_f31(N,R)/(R8*R**(N-8.0d0))
      ENDDO
      DISHHST_f31=ADISP 
      RETURN
      END
      
            FUNCTION DAMPHHTTRIP_f31(N,R)
c    *************************************************************** 
c     CALCULATES VARANDAS-BRANDAO DAMPING FUNCTIONS
c    ***************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c    ***************************************************************
      RR=2.0D0*R/(RMHHS+2.5D0*R0HH)
      DAMPHHTTRIP_f31=(1.0D0-EXP(-A(N)*RR-B(N)*RR*RR))**N
      RETURN
      END

      FUNCTION DVHHST_f31(R)
C     *************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  singlet H...H
C     *************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *************************************************************
      DVHHST_f31=DEHFHHST_f31(R)+DDISHHSTT_f31(R)
      RETURN
      END

      FUNCTION DEHFHHST_f31(R)
C     *************************************************************
C     TO COMPUTE THE EHF FOR  singlet H...H
C     *************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     **************************************************************
      DATA D,ASV/0.218973D0,1.91479d0,0.646041D0,0.346414D0/
      DATA (GAM(I),i=0,2)/1.22349D0,1.04334D0,0.208477D0/
      DATA AGEX,ALPHEX,GAMEX/0.8205D0,2.5D0,2.0D0/
c     **************************************************************
c     Derivada de EHFHHST_f31=POL*EXP(-Gama*X)
      X=R-RMHHS
      X2=X*X
      X3=X2*X
      POL=1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3
      Gama=Gam(0)*(1.0D0+Gam(1)*TANH(Gam(2)*X))      
      DGama=Gam(0)*Gam(1)*Gam(2)/COSH(Gam(2)*X)**2
      DPOLC=D/R**2*POL-D/R*(ASV(1)+2*ASV(2)*X+3*ASV(3)*X2)
      DEXP=-(Gama+DGama*X)*EXP(-Gama*X)
      DEHFHHST_f31=DPOLC*EXP(-GAMA*X)-D/R*POL*DEXP

c     Derivada de   EXCHHS=-AGEX*R**ALPHEX*EXP(-GAMEX*R)*D6
      RHH=(RMHHS+2.5D0*R0HH)/2.0d0
      RR=R/RHH
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6  
      DBASE=(A(6)+2*B(6)*RR)/RHH*EXP(-A(6)*RR-    
     1    B(6)*RR**2)
      DD6=6*(1-EXP(-A(6)*RR-B(6)*RR**2))**5*DBASE
      FEXC=-AGEX*R**ALPHEX*EXP(-GAMEX*R)
      DFEXC=-FEXC*GAMEX+ALPHEX*FEXC/R
      DEXCHHS=DFEXC*D6+DD6*FEXC
c     Total
      DEHFHHST_f31=DEHFHHST_f31+DEXCHHS
      RETURN
      END

      FUNCTION DDISHHSTT_f31(R)
c    ************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/dispc1t_f31/CHH(16),COHP(10)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
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
      DRR=RR/R
      T6=1.0D0-EXP(-A(6)*RR-B(6)*RR**2)
      T8=1.0D0-EXP(-A(8)*RR-B(8)*RR**2)
      T10=1.0D0-EXP(-A(10)*RR-B(10)*RR**2)
      T11=1.0D0-EXP(-A(11)*RR-B(11)*RR**2)
      T12=1.0D0-EXP(-A(12)*RR-B(12)*RR**2)
      T13=1.0D0-EXP(-A(13)*RR-B(13)*RR**2)
      T14=1.0D0-EXP(-A(14)*RR-B(14)*RR**2)
      T15=1.0D0-EXP(-A(15)*RR-B(15)*RR**2)
      T16=1.0D0-EXP(-A(16)*RR-B(16)*RR**2)
         
      DDISPLOCAL=6.0D0*CHH(6)/R6*T6**5*(T6/R+(1.0D0-T6)*(-A(6)-2.0D0
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
      DDISHHSTT_f31=DDISPLOCAL
      RETURN
      END
       

      FUNCTION THREBR_f31(R1,R2,R3)
c **************************************************************
c
C TO COMPUTE THE THREE BODY TERM IN COORDINATES R1,R2,R3
c
C **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON/coefft_f31/C(127)
      COMMON/refgeot_f31/R10,R20,R30

c **************************************************************
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
      Q19=Q18*Q1
      
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
      TQ33=TQ32*TQ3
      
      POLQ=
     1 C(1)+
     2 C(2)*Q1+
     3 C(3)*Q3+
     4 C(4)*Q12+
     5 C(5)*TQ1+
     6 C(6)*Q1*Q3+
     7 C(7)*TQ2+
     8 C(8)*Q13+
     9 C(9)*Q1*TQ1+
     1 C(10)*TQ3+
     2 C(11)*Q12*Q3+
     3 C(12)*Q1*TQ2+
     4 C(13)*Q3*TQ1+
     5 C(14)*Q14+
     6 C(15)*Q12*TQ1+
     7 C(16)*TQ12+
     8 C(17)*Q1*TQ3+
     9 C(18)*Q13*Q3+
     1 C(19)*Q12*TQ2+
     2 C(20)*Q1*Q3*TQ1+
     3 C(21)*Q3*TQ3+
     4 C(22)*TQ1*TQ2+
     5 C(23)*Q15+
     6 C(24)*Q13*TQ1+
     6 C(25)*Q1*TQ12+
     8 C(26)*Q12*TQ3+
     9 C(27)*TQ1*TQ3+
     1 C(28)*Q14*Q3+
     2 C(29)*Q13*TQ2+
     3 C(30)*Q12*Q3*TQ1+
     4 C(31)*Q1*Q3*TQ3+
     5 C(32)*Q1*TQ1*TQ2+
     6 C(33)*Q3*TQ12+
     7 C(34)*TQ2*TQ3+
     8 C(35)*Q16+
     9 C(36)*Q14*TQ1+
     1 C(37)*Q12*TQ12+
     2 C(38)*Q13*TQ3+
     3 C(39)*Q1*TQ3*TQ1+
     4 C(40)*TQ13+
     5 C(41)*TQ32+
     6 C(42)*Q15*Q3+
     7 C(43)*Q14*TQ2+
     8 C(44)*Q13*Q3*TQ1+
     9 C(45)*Q12*Q3*TQ3+
     1 C(46)*Q12*TQ1*TQ2+
     2 C(47)*Q1*Q3*TQ12+
     3 C(48)*Q1*TQ3*TQ2+
     4 C(49)*Q3*TQ1*TQ3+
     5 C(50)*TQ2*TQ12+
     6 C(51)*Q17+
     7 C(52)*Q15*TQ1+
     8 C(53)*Q13*TQ12+
     9 C(54)*Q14*TQ3+
     1 C(55)*Q12*TQ1*TQ3+
     2 C(56)*Q1*TQ13+
     3 C(57)*Q1*TQ32+
     4 C(58)*TQ3*TQ12+
     5 C(59)*Q3*Q16+
     6 C(60)*Q3*Q14*TQ1+
     7 C(61)*Q3*Q12*TQ12+
     8 C(62)*Q13*Q3*TQ3+
     9 C(63)*Q1*Q3*TQ3*TQ1+
     1 C(64)*Q3*TQ13+
     2 C(65)*Q3*TQ32+
     3 C(66)*Q15*TQ2+
     4 C(67)*Q13*TQ1*TQ2+
     5 C(68)*Q12*TQ3*TQ2+
     6 C(69)*Q1*TQ12*TQ2+
     7 C(70)*TQ3*TQ1*TQ2+
c  8grau     
     1 C(71)*Q18+
     2 C(72)*Q16*TQ1+
     3 C(73)*Q14*TQ12+
     4 C(74)*Q15*TQ3+
     5 C(75)*Q13*TQ3*TQ1+
     6 C(76)*Q12*TQ13+
     7 C(77)*Q12*TQ32+
     8 C(78)*Q1*TQ3*TQ12+
     9 C(79)*TQ14+
     1 C(80)*TQ32*TQ1+
     2 C(81)*Q17*Q3+
     3 C(82)*Q15*Q3*TQ1+
     4 C(83)*Q13*Q3*TQ12+
     5 C(84)*Q14*TQ32+
     6 C(85)*Q12*TQ32*TQ1+
     7 C(86)*Q1*Q3*TQ13+
     8 C(87)*Q1*Q3*TQ32+
     9 C(88)*Q16*TQ2+
     1 C(89)*Q14*TQ2*TQ1+
     2 C(90)*Q13*TQ3*TQ2+
     3 C(91)*Q12*TQ12*TQ2+
     4 C(92)*Q1*TQ3*TQ2*TQ1+
     5 C(93)*Q3*TQ3*TQ12+
     6 C(94)*TQ32*TQ2+
     7 C(95)*TQ13*TQ2+
c   9Grau
     1 Q1*(
     2 C(96)*Q18+
     3 C(97)*Q16*TQ1+
     4 C(98)*Q14*TQ12+
     5 C(99)*Q15*TQ3+
     6 C(100)*Q13*TQ3*TQ1+
     7 C(101)*Q12*TQ13+
     8 C(102)*Q12*TQ32+
     9 C(103)*Q1*TQ3*TQ12+
     1 C(104)*TQ14+
     2 C(105)*TQ32*TQ1+
     3 C(106)*Q17*Q3+
     4 C(107)*Q15*Q3*TQ1+
     5 C(108)*Q13*Q3*TQ12+
     6 C(109)*Q14*TQ32+
     7 C(110)*Q12*TQ32*TQ1+
     8 C(111)*Q1*Q3*TQ13+
     9 C(112)*Q1*Q3*TQ32+
     1 C(113)*Q16*TQ2+
     2 C(114)*Q14*TQ2*TQ1+
     3 C(115)*Q13*TQ3*TQ2+
     4 C(116)*Q12*TQ12*TQ2+
     5 C(117)*Q1*TQ3*TQ2*TQ1+
     6 C(118)*Q3*TQ3*TQ12+
     7 C(119)*TQ32*TQ2+
     8 C(120)*TQ13*TQ2)+
c  9grau termos extra
     1 C(121)*Q3*TQ1*TQ32+
     2 C(122)*Q3*TQ14+
     3 C(123)*TQ33+
     4 C(124)*TQ13*TQ3+
     5 C(125)*TQ2*TQ12*TQ3
     

      arg1=S1*C(126)
      arg2=S2*C(127)
      arg3=S3*C(127)      
      DECAY1=2.0d0*exp(-arg1)/(exp(arg1)+exp(-arg1))
      DECAY2=2.0d0*exp(-arg2)/(exp(arg2)+exp(-arg2))
      DECAY3=2.0d0*exp(-arg3)/(exp(arg3)+exp(-arg3))      

c      DECAY1=1.0D0-TANH(S1*C(126))
c      DECAY2=1.0D0-TANH(S2*C(127))
c      DECAY3=1.0D0-TANH(S3*C(127))



      THREBR_f31=POLQ*DECAY1*DECAY2*DECAY3*0.5d0**3.0d0

      RETURN
      END

      SUBROUTINE DTHREBR_f31(R1,R2,R3,DTDR1,DTDR2,DTDR3)
c **************************************************************
C TO COMPUTE THE DERIVATIVES OF THE THREE BODY TERM IN
c COORDINATES R1,R2,R3
C **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)


      COMMON/coefft_f31/C(127)
      COMMON/refgeot_f31/R10,R20,R30

c **************************************************************
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
      Q19=Q18*Q1
      
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
      TQ33=TQ32*TQ3
      
      POLQ=
     1 C(1)+
     2 C(2)*Q1+
     3 C(3)*Q3+
     4 C(4)*Q12+
     5 C(5)*TQ1+
     6 C(6)*Q1*Q3+
     7 C(7)*TQ2+
     8 C(8)*Q13+
     9 C(9)*Q1*TQ1+
     1 C(10)*TQ3+
     2 C(11)*Q12*Q3+
     3 C(12)*Q1*TQ2+
     4 C(13)*Q3*TQ1+
     5 C(14)*Q14+
     6 C(15)*Q12*TQ1+
     7 C(16)*TQ12+
     8 C(17)*Q1*TQ3+
     9 C(18)*Q13*Q3+
     1 C(19)*Q12*TQ2+
     2 C(20)*Q1*Q3*TQ1+
     3 C(21)*Q3*TQ3+
     4 C(22)*TQ1*TQ2+
     5 C(23)*Q15+
     6 C(24)*Q13*TQ1+
     6 C(25)*Q1*TQ12+
     8 C(26)*Q12*TQ3+
     9 C(27)*TQ1*TQ3+
     1 C(28)*Q14*Q3+
     2 C(29)*Q13*TQ2+
     3 C(30)*Q12*Q3*TQ1+
     4 C(31)*Q1*Q3*TQ3+
     5 C(32)*Q1*TQ1*TQ2+
     6 C(33)*Q3*TQ12+
     7 C(34)*TQ2*TQ3+
     8 C(35)*Q16+
     9 C(36)*Q14*TQ1+
     1 C(37)*Q12*TQ12+
     2 C(38)*Q13*TQ3+
     3 C(39)*Q1*TQ3*TQ1+
     4 C(40)*TQ13+
     5 C(41)*TQ32+
     6 C(42)*Q15*Q3+
     7 C(43)*Q14*TQ2+
     8 C(44)*Q13*Q3*TQ1+
     9 C(45)*Q12*Q3*TQ3+
     1 C(46)*Q12*TQ1*TQ2+
     2 C(47)*Q1*Q3*TQ12+
     3 C(48)*Q1*TQ3*TQ2+
     4 C(49)*Q3*TQ1*TQ3+
     5 C(50)*TQ2*TQ12+
     6 C(51)*Q17+
     7 C(52)*Q15*TQ1+
     8 C(53)*Q13*TQ12+
     9 C(54)*Q14*TQ3+
     1 C(55)*Q12*TQ1*TQ3+
     2 C(56)*Q1*TQ13+
     3 C(57)*Q1*TQ32+
     4 C(58)*TQ3*TQ12+
     5 C(59)*Q3*Q16+
     6 C(60)*Q3*Q14*TQ1+
     7 C(61)*Q3*Q12*TQ12+
     8 C(62)*Q13*Q3*TQ3+
     9 C(63)*Q1*Q3*TQ3*TQ1+
     1 C(64)*Q3*TQ13+
     2 C(65)*Q3*TQ32+
     3 C(66)*Q15*TQ2+
     4 C(67)*Q13*TQ1*TQ2+
     5 C(68)*Q12*TQ3*TQ2+
     6 C(69)*Q1*TQ12*TQ2+
     7 C(70)*TQ3*TQ1*TQ2+
c  8grau     
     1 C(71)*Q18+
     2 C(72)*Q16*TQ1+
     3 C(73)*Q14*TQ12+
     4 C(74)*Q15*TQ3+
     5 C(75)*Q13*TQ3*TQ1+
     6 C(76)*Q12*TQ13+
     7 C(77)*Q12*TQ32+
     8 C(78)*Q1*TQ3*TQ12+
     9 C(79)*TQ14+
     1 C(80)*TQ32*TQ1+
     2 C(81)*Q17*Q3+
     3 C(82)*Q15*Q3*TQ1+
     4 C(83)*Q13*Q3*TQ12+
     5 C(84)*Q14*TQ32+
     6 C(85)*Q12*TQ32*TQ1+
     7 C(86)*Q1*Q3*TQ13+
     8 C(87)*Q1*Q3*TQ32+
     9 C(88)*Q16*TQ2+
     1 C(89)*Q14*TQ2*TQ1+
     2 C(90)*Q13*TQ3*TQ2+
     3 C(91)*Q12*TQ12*TQ2+
     4 C(92)*Q1*TQ3*TQ2*TQ1+
     5 C(93)*Q3*TQ3*TQ12+
     6 C(94)*TQ32*TQ2+
     7 C(95)*TQ13*TQ2+
c   9Grau
     1 Q1*(
     2 C(96)*Q18+
     3 C(97)*Q16*TQ1+
     4 C(98)*Q14*TQ12+
     5 C(99)*Q15*TQ3+
     6 C(100)*Q13*TQ3*TQ1+
     7 C(101)*Q12*TQ13+
     8 C(102)*Q12*TQ32+
     9 C(103)*Q1*TQ3*TQ12+
     1 C(104)*TQ14+
     2 C(105)*TQ32*TQ1+
     3 C(106)*Q17*Q3+
     4 C(107)*Q15*Q3*TQ1+
     5 C(108)*Q13*Q3*TQ12+
     6 C(109)*Q14*TQ32+
     7 C(110)*Q12*TQ32*TQ1+
     8 C(111)*Q1*Q3*TQ13+
     9 C(112)*Q1*Q3*TQ32+
     1 C(113)*Q16*TQ2+
     2 C(114)*Q14*TQ2*TQ1+
     3 C(115)*Q13*TQ3*TQ2+
     4 C(116)*Q12*TQ12*TQ2+
     5 C(117)*Q1*TQ3*TQ2*TQ1+
     6 C(118)*Q3*TQ3*TQ12+
     7 C(119)*TQ32*TQ2+
     8 C(120)*TQ13*TQ2)+
c  9grau termos extra
     1 C(121)*Q3*TQ1*TQ32+
     2 C(122)*Q3*TQ14+
     3 C(123)*TQ33+
     4 C(124)*TQ13*TQ3+
     5 C(125)*TQ2*TQ12*TQ3

      arg1=S1*C(126)
      arg2=S2*C(127)
      arg3=S3*C(127)      
      DECAY1=2.0d0*exp(-arg1)/(exp(arg1)+exp(-arg1))
      DECAY2=2.0d0*exp(-arg2)/(exp(arg2)+exp(-arg2))
      DECAY3=2.0d0*exp(-arg3)/(exp(arg3)+exp(-arg3))

c      DECAY1=1.0D0-TANH(S1*C(126))
c      DECAY2=1.0D0-TANH(S2*C(127))
c      DECAY3=1.0D0-TANH(S3*C(127))



      DECAY=DECAY1*DECAY2*DECAY3*0.5d0**3.0d0

      DQ1R1=RO3
      DQ1R2=RO3
      DQ1R3=RO3

      DQ2R1=0.0D0
      DQ2R2=RO2
      DQ2R3=-RO2

      DQ3R1=2.0D0*RO6
      DQ3R2=-RO6
      DQ3R3=-RO6


c DDECAYR1=DECAY2*DECAY3*2*C(126)*S1*EXP(-C(126)*S1**2)
c DDECAYR2=DECAY1*DECAY3*2*C(127)*S2*EXP(-C(127)*S2**2)
c DDECAYR3=DECAY1*DECAY2*2*C(127)*S3*EXP(-C(127)*S3**2)

      DDECAYR1=DECAY2*DECAY3*C(126)/
     1 COSH(S1*C(126))**2*0.5d0**3.0d0
      DDECAYR2=DECAY1*DECAY3*C(127)/
     1 COSH(S2*C(127))**2*0.5d0**3.0d0
      DDECAYR3=DECAY1*DECAY2*C(127)/
     1 COSH(S3*C(127))**2*0.5d0**3.0d0

      DTQ1Q2=2.0D0*Q2
      DTQ1Q3=2.0D0*Q3

      DTQ2Q2=2.0D0*Q2
      DTQ2Q3=-2.0D0*Q3

      DTQ3Q2=-6.0D0*Q2*Q3
      DTQ3Q3=3.0D0*Q32-3.0d0*q22

      DPOQ1=
     1 C(2)+
     2 C(4)*2.0D0*Q1+
     3 C(6)*Q3+3.0D0*
     4 C(8)*Q12+
     5 C(9)*TQ1+C(11)*2.0D0*Q1*Q3+
     6 C(12)*TQ2+
     7 C(14)*4.0D0*Q13+
     8 C(15)*2.0D0*Q1*TQ1+
     9 C(17)*TQ3+
     1 C(18)*3.0D0*Q12*Q3+C(19)*2.0D0*Q1*TQ2+
     2 C(20)*Q3*TQ1+
     3 C(23)*5.0D0*Q14+
     4 C(24)*3.0D0*Q12*TQ1+
     5 C(25)*TQ12+
     6 C(26)*2.0D0*Q1*TQ3+
     7 C(28)*4.0D0*Q13*Q3+
     8 C(29)*3.0D0*Q12*TQ2+
     9 C(30)*2.0D0*Q1*Q3*TQ1+
     1 C(31)*Q3*TQ3+
     2 C(32)*TQ1*TQ2+
     3 C(35)*6.0D0*Q15+
     4 C(36)*4.0D0*Q13*TQ1+
     5 C(37)*2.0D0*Q1*TQ12+
     6 C(38)*3.0D0*Q12*TQ3+
     7 C(39)*TQ3*TQ1+
     8 C(42)*5.0D0*Q14*Q3+
     9 C(43)*4.0D0*Q13*TQ2+
     1 C(44)*3.0D0*Q12*Q3*TQ1+
     2 C(45)*2.0D0*Q1*Q3*TQ3+
     3 C(46)*2.0D0*Q1*TQ1*TQ2+
     4 C(47)*Q3*TQ12+
     5 C(48)*TQ3*TQ2+
     6 C(51)*7.0d0*Q16+
     7 C(52)*5.0d0*Q14*TQ1+
     8 C(53)*3.0d0*Q12*TQ12+
     5 C(54)*4.0d0*Q13*TQ3+
     6 C(55)*2.0d0*Q1*TQ1*TQ3+
     7 C(56)*TQ13+
     8 C(57)*TQ32+
     9 C(59)*6.0d0*Q3*Q15+
     1 C(60)*4.0d0*Q3*Q13*TQ1+
     2 C(61)*2.0d0*Q3*Q1*TQ12+
     3 C(62)*3.0d0*Q12*Q3*TQ3+
     4 C(63)*Q3*TQ3*TQ1+5.0d0*
     5 C(66)*Q14*TQ2+
     6 C(67)*3.0D0*Q12*TQ1*TQ2+
     7 C(68)*2.0D0*Q1*TQ3*TQ2+
     8 C(69)*TQ12*TQ2+
     9 C(71)*8.0D0*Q17+
     1 C(72)*6.0D0*Q15*TQ1+
     2 C(73)*4.0D0*Q13*TQ12+
     3 C(74)*5.0D0*Q14*TQ3+
     4 C(75)*3.0D0*Q12*TQ3*TQ1+
     5 C(76)*2.0D0*Q1*TQ13+
     6 C(77)*2.0D0*Q1*TQ32+
     7 C(78)*TQ3*TQ12+
     8 C(81)*7.0D0*Q16*Q3+
     9 C(82)*5.0D0*Q14*Q3*TQ1+
     1 C(83)*3.0D0*Q12*Q3*TQ12+
     2 C(84)*4.0D0*Q13*TQ32+
     3 C(85)*2.0D0*Q1*TQ32*TQ1+
     4 C(86)*Q3*TQ13+
     5 C(87)*Q3*TQ32+
     6 C(88)*6.0D0*Q15*TQ2+
     7 C(89)*4.0D0*Q13*TQ2*TQ1+
     8 C(90)*3.0D0*Q12*TQ3*TQ2+
     9 C(91)*2.0D0*Q1*TQ12*TQ2+
     1 C(92)*TQ3*TQ2*TQ1+     
c   9grau
     1 C(96)*9.0D0*Q18+
     2 C(97)*7.0D0*Q16*TQ1+
     3 C(98)*5.0D0*Q14*TQ12+
     4 C(99)*6.0D0*Q15*TQ3+
     5 C(100)*4.0D0*Q13*TQ3*TQ1+
     6 C(101)*3.0D0*Q12*TQ13+
     7 C(102)*3.0D0*Q12*TQ32+
     8 C(103)*2.0D0*Q1*TQ3*TQ12+
     9 C(104)*TQ14+
     1 C(105)*TQ32*TQ1+
     2 C(106)*8.0D0*Q17*Q3+
     3 C(107)*6.0D0*Q15*Q3*TQ1+
     4 C(108)*4.0D0*Q13*Q3*TQ12+
     5 C(109)*5.0D0*Q14*TQ32+
     6 C(110)*3.0D0*Q12*TQ32*TQ1+
     7 C(111)*2.0D0*Q1*Q3*TQ13+
     8 C(112)*2.0D0*Q1*Q3*TQ32+
     9 C(113)*7.0D0*Q16*TQ2+
     1 C(114)*5.0D0*Q14*TQ2*TQ1+
     2 C(115)*4.0D0*Q13*TQ3*TQ2+
     3 C(116)*3.0D0*Q12*TQ12*TQ2+
     4 C(117)*2.0D0*Q1*TQ3*TQ2*TQ1+
     5 C(118)*Q3*TQ3*TQ12+
     6 C(119)*TQ32*TQ2+
     7 C(120)*TQ13*TQ2
c 9grau termos extra
C NAO EXISTEM     
        

      DPOQ2=
     1 C(5)*2.0D0*Q2+
     2 C(7)*2.0D0*Q2+
     3 C(9)*2.0D0*Q2*Q1+
     4 C(12)*2.0D0*Q1*Q2-
     5 C(10)*6.0D0*Q2*Q3+
     6 C(13)*2.0D0*Q2*Q3+
     7 C(15)*2.0d0*Q12*Q2+
     8 C(16)*4.0D0*Q23+
     9 C(16)*4.0D0*Q2*Q32-
     1 C(17)*6.0D0*Q1*Q3*Q2+
     2 C(19)*2.0D0*Q12*Q2+
     3 C(20)*2.0D0*Q1*Q2*Q3-
     4 C(21)*6.0D0*Q32*Q2+
     5 C(22)*4.0D0*Q23+
     6 C(24)*2.0D0*Q2*Q13+
     7 C(25)*2.0D0*Q1*DTQ1Q2*TQ1-
     8 C(26)*6.0D0*Q12*Q2*Q3+
     9 C(27)*(DTQ1Q2*TQ3+DTQ3Q2*TQ1)+
     7 C(29)*Q13*DTQ2Q2+
     1 C(30)*Q12*Q3*DTQ1Q2+
     2 C(31)*Q1*Q3*DTQ3Q2+
     3 C(32)*Q1*(DTQ1Q2*TQ2+DTQ2Q2*TQ1)+
     4 C(33)*2.0D0*Q3*DTQ1Q2*TQ1+
     5 C(34)*(DTQ2Q2*TQ3+DTQ3Q2*TQ2)+
     6 C(36)*Q14*DTQ1Q2+
     7 C(37)*2.0D0*Q12*DTQ1Q2*TQ1+
     8 C(38)*Q13*DTQ3Q2+
     9 C(39)*Q1*(DTQ3Q2*TQ1+DTQ1Q2*TQ3)+
     1 C(40)*3.0D0*TQ12*DTQ1Q2+
     2 C(41)*2.0D0*DTQ3Q2*TQ3+
     3 C(43)*Q14*DTQ2Q2+
     4 C(44)*Q13*Q3*DTQ1Q2+
     5 C(45)*Q12*Q3*DTQ3Q2+
     6 C(46)*Q12*(DTQ1Q2*TQ2+DTQ2Q2*TQ1)+
     7 C(47)*2.0D0*Q1*Q3*DTQ1Q2*TQ1+
     8 C(48)*Q1*(DTQ3Q2*TQ2+DTQ2Q2*TQ3)+
     9 C(49)*Q3*(DTQ1Q2*TQ3+DTQ3Q2*TQ1)+
     1 C(50)*(DTQ2Q2*TQ12+2.0D0*TQ1*DTQ1Q2*TQ2)+
     2 C(52)*Q15*DTQ1Q2+
     3 C(53)*2.0D0*Q13*TQ1*DTQ1Q2+
     4 C(54)*Q14*DTQ3Q2+
     5 C(55)*Q12*(DTQ1Q2*TQ3+DTQ3Q2*TQ1)+
     6 C(56)*3.0D0*Q1*TQ12*DTQ1Q2+
     7 C(57)*2.0D0*Q1*TQ3*DTQ3Q2+
     8 C(58)*(DTQ3Q2*TQ12+2.0D0*TQ1*DTQ1Q2*TQ3)+
     9 C(60)*Q3*Q14*DTQ1Q2+
     1 C(61)*2.0D0*Q3*Q12*TQ1*DTQ1Q2+
     2 C(62)*Q13*Q3*DTQ3Q2+
     3 C(63)*Q1*Q3*(DTQ3Q2*TQ1+DTQ1Q2*TQ3)+
     4 C(64)*3.0D0*Q3*TQ12*DTQ1Q2+
     5 C(65)*2.0D0*Q3*TQ3*DTQ3Q2+
     6 C(66)*Q15*DTQ2Q2+
     7 C(67)*Q13*(DTQ1Q2*TQ2+DTQ2Q2*TQ1)+
     8 C(68)*Q12*(DTQ3Q2*TQ2+DTQ2Q2*TQ3)+
     9 C(69)*Q1*(2.0D0*TQ1*DTQ1Q2*TQ2+DTQ2Q2*TQ12)+
     1 C(70)*((DTQ3Q2*TQ1+DTQ1Q2*TQ3)*TQ2+DTQ2Q2*TQ1*TQ3)+
     2 C(72)*Q16*DTQ1Q2+
     3 C(73)*Q14*2.0D0*TQ1*DTQ1Q2+
     4 C(74)*Q15*DTQ3Q2+
     5 C(75)*Q13*(DTQ3Q2*TQ1+DTQ1Q2*TQ3)+
     6 C(76)*Q12*3.0D0*TQ12*DTQ1Q2+
     7 C(77)*Q12*2.0D0*TQ3*DTQ3Q2+
     8 C(78)*Q1*(DTQ3Q2*TQ12+2.0D0*TQ1*DTQ1Q2*TQ3)+
     9 C(79)*4.0D0*TQ13*DTQ1Q2+
     1 C(80)*(2.0D0*TQ3*DTQ3Q2*TQ1+DTQ1Q2*TQ32)+
     2 C(82)*Q15*Q3*DTQ1Q2+
     3 C(83)*Q13*Q3*2.0D0*TQ1*DTQ1Q2+
     4 C(84)*Q14*2.0D0*TQ3*DTQ3Q2+
     5 C(85)*Q12*(2.0D0*DTQ3Q2*TQ3*TQ1+DTQ1Q2*TQ32)+
     6 C(86)*Q1*Q3*3.0D0*TQ12*DTQ1Q2+
     7 C(87)*Q1*Q3*2.0D0*TQ3*DTQ3Q2+
     8 C(88)*Q16*DTQ2Q2+
     9 C(89)*Q14*(DTQ2Q2*TQ1+DTQ1Q2*TQ2)+
     1 C(90)*Q13*(DTQ3Q2*TQ2+DTQ2Q2*TQ3)+
     2 C(91)*Q12*(2.0D0*TQ1*DTQ1Q2*TQ2+DTQ2Q2*TQ12)+
     3 C(92)*Q1*((DTQ3Q2*TQ2+DTQ2Q2*TQ3)*TQ1+DTQ1Q2*TQ3*TQ2)+
     4 C(93)*Q3*(DTQ3Q2*TQ12+2.0D0*TQ1*DTQ1Q2*TQ3)+
     5 C(94)*(2.0D0*TQ3*DTQ3Q2*TQ2+DTQ2Q2*TQ32)+
     6 C(95)*(3.0D0*TQ12*DTQ1Q2*TQ2+DTQ2Q2*TQ13)+
c   9grau
     1 Q1*(
     2 C(97)*Q16*DTQ1Q2+
     3 C(98)*Q14*2.0D0*TQ1*DTQ1Q2+
     4 C(99)*Q15*DTQ3Q2+
     5 C(100)*Q13*(DTQ3Q2*TQ1+TQ3*DTQ1Q2)+
     6 C(101)*Q12*3.0D0*TQ12*DTQ1Q2+
     7 C(102)*Q12*2.0D0*TQ3*DTQ3Q2+
     8 C(103)*Q1*(DTQ3Q2*TQ12+TQ3*2.0D0*TQ1*DTQ1Q2)+
     9 C(104)*4.0D0*TQ13*DTQ1Q2+
     1 C(105)*(2.0D0*TQ3*DTQ3Q2*TQ1+TQ32*DTQ1Q2)+
     2 C(107)*Q15*Q3*DTQ1Q2+
     3 C(108)*Q13*Q3*2.0D0*TQ1*DTQ1Q2+
     4 C(109)*Q14*2.0D0*TQ3*DTQ3Q2+
     5 C(110)*Q12*(2.0D0*TQ3*DTQ3Q2*TQ1+TQ32*DTQ1Q2)+
     6 C(111)*Q1*Q3*3.0D0*TQ12*DTQ1Q2+
     7 C(112)*Q1*Q3*2.0D0*TQ3*DTQ3Q2+
     8 C(113)*Q16*DTQ2Q2+
     9 C(114)*Q14*(DTQ2Q2*TQ1+TQ2*DTQ1Q2)+
     1 C(115)*Q13*(DTQ3Q2*TQ2+TQ3*DTQ2Q2)+
     2 C(116)*Q12*(2.0D0*TQ1*DTQ1Q2*TQ2+TQ12*DTQ2Q2)+
     3 C(117)*Q1*(DTQ3Q2*TQ2*TQ1+TQ3*DTQ2Q2*TQ1+TQ3*TQ2*DTQ1Q2)+
     4 C(118)*Q3*(DTQ3Q2*TQ12+TQ3*2.0D0*TQ1*DTQ1Q2)+
     5 C(119)*(2.0D0*TQ3*DTQ3Q2*TQ2+TQ32*DTQ2Q2)+
     6 C(120)*(3.0D0*TQ12*DTQ1Q2*TQ2+TQ13*DTQ2Q2))+
c   9grau termos extra
     1 C(121)*Q3*(DTQ1Q2*TQ32+TQ1*2.0D0*TQ3*DTQ3Q2)+
     2 C(122)*Q3*4.0D0*TQ13*DTQ1Q2+
     3 C(123)*3.0D0*TQ32*DTQ3Q2+
     4 C(124)*(3.0D0*TQ12*DTQ1Q2*TQ3+TQ13*DTQ3Q2)+
     5 C(125)*(DTQ2Q2*TQ12*TQ3+TQ2*2.0D0*TQ1*DTQ1Q2*TQ3+
     6 TQ2*TQ12*DTQ3Q2)




      DPOQ3=C(3)+
     1 C(5)*2.0D0*Q3+
     2 C(6)*Q1-
     3 C(7)*2.0D0*Q3+
     4 C(9)*2.0D0*Q1*Q3+
     5 C(10)*3.0D0*Q32-
     6 C(10)*3.0D0*Q22+
     7 C(11)*Q12-
     8 C(12)*2.0D0*Q1*Q3+
     9 C(13)*Q22+
     1 C(13)*3.0D0*Q32+
     2 C(15)*2.0D0*Q12*Q3+
     3 C(16)*4.0D0*Q33+
     4 C(16)*4.0D0*Q22*Q3+
     5 C(17)*3.0D0*Q1*Q32-
     6 C(17)*3.0D0*Q1*Q22+
     7 C(18)*Q13-
     8 C(19)*2.0D0*Q12*Q3+
     9 C(20)*Q1*Q22+
     1 C(20)*3.0D0*Q1*Q32+
     2 C(21)*4.0D0*Q33-
     3 C(21)*6.0D0*Q3*Q22-
     4 C(22)*4.0D0*Q33+
     5 C(24)*Q13*DTQ1Q3+
     6 C(25)*2.0D0*Q1*DTQ1Q3*TQ1+
     7 C(26)*Q12*DTQ3Q3+
     8 C(27)*(DTQ1Q3*TQ3+DTQ3Q3*TQ1)+
     9 C(28)*Q14+
     1 C(29)*Q13*DTQ2Q3+
     2 C(30)*Q12*(TQ1+DTQ1Q3*Q3)+
     3 C(31)*Q1*(TQ3+DTQ3Q3*Q3)+
     4 C(32)*Q1*(DTQ1Q3*TQ2+DTQ2Q3*TQ1)+
     5 C(33)*(TQ12+2.0D0*TQ1*DTQ1Q3*Q3)+
     6 C(34)*(DTQ2Q3*TQ3+DTQ3Q3*TQ2)+
     7 C(36)*Q14*DTQ1Q3+
     8 C(37)*2.0D0*Q12*TQ1*DTQ1Q3+
     9 C(38)*Q13*DTQ3Q3+
     1 C(39)*Q1*(DTQ3Q3*TQ1+DTQ1Q3*TQ3)+
     2 C(40)*3.0D0*TQ12*DTQ1Q3+
     3 C(41)*2.0D0*TQ3*DTQ3Q3+
     4 C(42)*Q15+
     5 C(43)*Q14*DTQ2Q3+
     6 C(44)*Q13*(TQ1+Q3*DTQ1Q3)+
     7 C(45)*Q12*(TQ3+Q3*DTQ3Q3)+
     8 C(46)*Q12*(DTQ1Q3*TQ2+DTQ2Q3*TQ1)+
     9 C(47)*Q1*(TQ12+2.0D0*Q3*TQ1*DTQ1Q3)+
     1 C(48)*Q1*(DTQ3Q3*TQ2+DTQ2Q3*TQ3)+
     2 C(49)*(TQ3*TQ1+(DTQ1Q3*TQ3+DTQ3Q3*TQ1)*Q3)+
     3 C(50)*(DTQ2Q3*TQ12+2.0D0*TQ1*DTQ1Q3*TQ2)+
     4 C(52)*Q15*DTQ1Q3+
     5 C(53)*2.0D0*Q13*TQ1*DTQ1Q3+
     6 C(54)*Q14*DTQ3Q3+
     7 C(55)*Q12*(DTQ1Q3*TQ3+DTQ3Q3*TQ1)+
     8 C(56)*3.0D0*Q1*TQ12*DTQ1Q3+
     9 C(57)*2.0D0*Q1*TQ3*DTQ3Q3+
     1 C(58)*(DTQ3Q3*TQ12+2.0D0*TQ1*DTQ1Q3*TQ3)+
     2 C(59)*Q16+
     3 C(60)*Q14*(TQ1+Q3*DTQ1Q3)+
     4 C(61)*Q12*(TQ12+Q3*2.0D0*TQ1*DTQ1Q3)+
     5 C(62)*Q13*(TQ3+Q3*DTQ3Q3)+
     6 C(63)*Q1*((TQ1+DTQ1Q3*Q3)*TQ3+DTQ3Q3*Q3*TQ1)+
     7 C(64)*(TQ13+Q3*3.0D0*TQ12*DTQ1Q3)+
     8 C(65)*(TQ32+2.0D0*Q3*TQ3*DTQ3Q3)+
     9 C(66)*Q15*DTQ2Q3+
     1 C(67)*Q13*(DTQ1Q3*TQ2+DTQ2Q3*TQ1)+
     2 C(68)*Q12*(DTQ3Q3*TQ2+DTQ2Q3*TQ3)+
     3 C(69)*Q1*(2.0D0*DTQ1Q3*TQ1*TQ2+DTQ2Q3*TQ12)+
     4 C(70)*((DTQ3Q3*TQ1+DTQ1Q3*TQ3)*TQ2+DTQ2Q3*TQ1*TQ3)+
     5 C(72)*Q16*DTQ1Q3+
     6 C(73)*Q14*2.0D0*TQ1*DTQ1Q3+
     7 C(74)*Q15*DTQ3Q3+
     8 C(75)*Q13*(DTQ3Q3*TQ1+DTQ1Q3*TQ3)+
     9 C(76)*Q12*3.0D0*TQ12*DTQ1Q3+
     1 C(77)*Q12*2.0D0*TQ3*DTQ3Q3+
     2 C(78)*Q1*(DTQ3Q3*TQ12+2.0D0*TQ1*DTQ1Q3*tq3)+
     3 C(79)*4.0D0*TQ13*DTQ1Q3+
     4 C(80)*(DTQ3Q3*TQ3*2.0D0*TQ1+DTQ1Q3*TQ32)+
     5 C(81)*Q17+
     6 C(82)*Q15*(TQ1+Q3*DTQ1Q3)+
     7 C(83)*Q13*(2.0D0*TQ1*DTQ1Q3*Q3+TQ12)+
     8 C(84)*Q14*2.0D0*TQ3*DTQ3Q3+
     9 C(85)*Q12*(2.0D0*TQ3*DTQ3Q3*TQ1+DTQ1Q3*TQ32)+
     1 C(86)*Q1*(Q3*3.0D0*TQ12*DTQ1Q3+TQ13)+
     2 C(87)*Q1*(Q3*2.0D0*TQ3*DTQ3Q3+TQ32)+
     3 C(88)*Q16*DTQ2Q3+
     4 C(89)*Q14*(DTQ2Q3*TQ1+DTQ1Q3*TQ2)+
     5 C(90)*Q13*(DTQ3Q3*TQ2+DTQ2Q3*TQ3)+
     6 C(91)*Q12*(2.0D0*TQ1*DTQ1Q3*TQ2+DTQ2Q3*TQ12)+
     7 C(92)*Q1*((DTQ3Q3*TQ2+DTQ2Q3*TQ3)*TQ1+DTQ1Q3*TQ2*TQ3)+
     8 C(93)*((Q3*DTQ3Q3+TQ3)*TQ12+2.0D0*TQ1*DTQ1Q3*Q3*TQ3)+
     9 C(94)*(2.0D0*TQ3*DTQ3Q3*TQ2+DTQ2Q3*TQ32)+
     1 C(95)*(3.0D0*TQ12*DTQ1Q3*TQ2+DTQ2Q3*TQ13)+
c   9grau
     1 Q1*(
     2 C(97)*Q16*DTQ1Q3+
     3 C(98)*Q14*2.0D0*TQ1*DTQ1Q3+
     4 C(99)*Q15*DTQ3Q3+
     5 C(100)*Q13*(DTQ3Q3*TQ1+TQ3*DTQ1Q3)+
     6 C(101)*Q12*3.0D0*TQ12*DTQ1Q3+
     7 C(102)*Q12*2.0D0*TQ3*DTQ3Q3+
     8 C(103)*Q1*(DTQ3Q3*TQ12+TQ3*2.0D0*TQ1*DTQ1Q3)+
     9 C(104)*4.0D0*TQ13*DTQ1Q3+
     1 C(105)*(2.0D0*TQ3*DTQ3Q3*TQ1+TQ32*DTQ1Q3)+
     2 C(106)*Q17+C(107)*Q15*(TQ1+Q3*DTQ1Q3)+
     3 C(108)*Q13*(TQ12+Q3*2.0D0*TQ1*DTQ1Q3)+
     4 C(109)*Q14*2.0D0*TQ3*DTQ3Q3+
     5 C(110)*Q12*(2.0D0*TQ3*DTQ3Q3*TQ1+TQ32*DTQ1Q3)+
     6 C(111)*Q1*(TQ13+Q3*3.0D0*TQ12*DTQ1Q3)+
     7 C(112)*Q1*(TQ32+Q3*2.0D0*TQ3*DTQ3Q3)+
     8 C(113)*Q16*DTQ2Q3+
     9 C(114)*Q14*(DTQ2Q3*TQ1+TQ2*DTQ1Q3)+
     1 C(115)*Q13*(DTQ3Q3*TQ2+TQ3*DTQ2Q3)+
     2 C(116)*Q12*(2.0D0*TQ1*DTQ1Q3*TQ2+TQ12*DTQ2Q3)+
     3 C(117)*Q1*(DTQ3Q3*TQ2*TQ1+TQ3*DTQ2Q3*TQ1+TQ3*TQ2*DTQ1Q3)+
     4 C(118)*(TQ3*TQ12+Q3*DTQ3Q3*TQ12+Q3*TQ3*2.0D0*TQ1*DTQ1Q3)+
     5 C(119)*(2.0D0*TQ3*DTQ3Q3*TQ2+TQ32*DTQ2Q3)+
     6 C(120)*(3.0D0*TQ12*DTQ1Q3*TQ2+TQ13*DTQ2Q3))+
c   9grau termos extra
     1 C(121)*(TQ1*TQ32+Q3*DTQ1Q3*TQ32+Q3*TQ1*2.0D0*TQ3*DTQ3Q3)+
     2 C(122)*(TQ14+Q3*4.0D0*TQ13*DTQ1Q3)+
     3 C(123)*3.0D0*TQ32*DTQ3Q3+
     4 C(124)*(3.0D0*TQ12*DTQ1Q3*TQ3+TQ13*DTQ3Q3)+
     5 C(125)*(DTQ2Q3*TQ12*TQ3+TQ2*2.0D0*TQ1*DTQ1Q3*TQ3+TQ2*TQ12*DTQ3Q3)





      DTDR1=(DPOQ1*DQ1R1+DPOQ3*DQ3R1)*DECAY-POLQ*DDECAYR1
      DTDR2=(DPOQ1*DQ1R2+DPOQ2*DQ2R2+DPOQ3*DQ3R2)*DECAY-POLQ*DDECAYR2
      DTDR3=(DPOQ1*DQ1R3+DPOQ2*DQ2R3+DPOQ3*DQ3R3)*DECAY-POLQ*DDECAYR3
      RETURN
      END



      
      SUBROUTINE GEOMRLT_f31(R1,R2,R3)  
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONSTT_f31/PI
      COMMON/dispc1t_f31/CHH(16),COHP(10)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,COST2,
     1  COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,C8OHHST11,C10OHHST11,
     1    C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
c     ******************************************************************     
c      CONSTT=4.5d0
c      CONSTT=3.0D0
c      CONSTT=2.5D0
c       CONSTT=2.8D0
         CONSTT=2.0D0
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
C     CONSTTANTES ADICIONAIS USADAS NA VEOHHST_f31, E VINDHOHPT
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
      DAMPOHHST5=DAMPOHHST_f31(RG1,5,R1)/RG15
C     *******************
 
      umsabc=(0.5d0+0.5d0*tanh(CONSTT*((rg1/r1)-2.0d0)))
      sabc=1.0d0-umsabc**2

      umscab=(0.5d0+0.5d0*tanh(CONSTT*((rg2/r2)-2.0d0)))
      scab=1.0d0-umscab**2

      umsbac=(0.5d0+0.5d0*tanh(CONSTT*((rg3/r3)-2.0d0)))
      sbac=1.0d0-umsbac**2

      FAB=sabc*sbac
      FAC=sabc*scab
      FbC=sbac*scab
      FABC=umsabc
      FCAB=umscab
      FBAC=umsbac
      
      
      DFABCR1=0.5D0*(CONSTT*DRG1R1*R1-CONSTT*RG1)/R12*(1.0D0/COSH
     1        (CONSTT*(RG1/R1-2.0D0)))**2
      DFCABR1=0.5D0*CONSTT*DRG2R1/R2*(1.0D0/COSH(CONSTT*(
     1        RG2/R2-2.0D0)))**2
      DFBACR1=0.5D0*CONSTT*DRG3R1/R3*(1.0D0/COSH(CONSTT*(
     1        RG3/R3-2.0D0)))**2

      DFABCR2=0.5d0*(1.0d0/cosH(CONSTT*(RG1/R1-2.0D0)))**2*
     1        CONSTT*DRG1R2/R1
      DFCABR2=0.5d0*(1.0d0/cosH(CONSTT*(RG2/R2-2.0D0)))**2*
     1        (CONSTT*DRG2R2*R2-CONSTT*RG2)/R22
      DFBACR2=0.5d0*(1.0d0/cosH(CONSTT*(RG3/R3-2.0D0)))**2*
     1        CONSTT*DRG3R2/R3

      DFABCR3=0.5D0*(1.0D0/COSH(CONSTT*(RG1/R1-2.0D0)))**2*
     1        CONSTT*DRG1R3/R1
      DFCABR3=0.5D0*(1.0D0/COSH(CONSTT*(RG2/R2-2.0D0)))**2*
     1        CONSTT*DRG2R3/R2
      DFBACR3=0.5D0*(1.0D0/COSH(CONSTT*(RG3/R3-2.0D0)))**2*
     1        (CONSTT*DRG3R3*R3-CONSTT*RG3)/R32
      
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
       
      
      C6AC2=COHP(6)
      C8AC2=COHP(8)
      C10AC2=COHP(10)
      C6AB2=COHP(6)
      C8AB2=COHP(8)
      C10AB2=COHP(10)
      
c calculado a partir das dist. de eq.(POL94:7651) RMHHS=2.8619327675d0

      
      
      C6HOHPT12=C6HOHPT_f31(r2,cost2)
      C8HOHPT12=C8HOHPT_f31(r2,cost2)
      C10HOHPT12=C10HOHPT_f31(r2,cost2)
      
      C6HOHPT13=C6HOHPT_f31(r3,cost3)
      C8HOHPT13=C8HOHPT_f31(r3,cost3)
      C10HOHPT13=C10HOHPT_f31(r3,cost3)
      
      
      DAMPHOHPT26=DAMPHOHPT_f31(RG2,6,r2)
      DAMPHOHPT28=DAMPHOHPT_f31(RG2,8,r2)
      DAMPHOHPT210=DAMPHOHPT_f31(RG2,10,r2)
      
      DAMPHOHPT36=DAMPHOHPT_f31(RG3,6,r3)
      DAMPHOHPT38=DAMPHOHPT_f31(RG3,8,r3)
      DAMPHOHPT310=DAMPHOHPT_f31(RG3,10,r3)
      
      C6OHHST11=C6OHHST_f31(r1,cost1)
      C8OHHST11=C8OHHST_f31(r1,cost1)
      C10OHHST11=C10OHHST_f31(r1,cost1)
      
      
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
            
      DAMPOHHST16=DAMPOHHST_f31(RG1,6,r1)
      DAMPOHHST18=DAMPOHHST_f31(RG1,8,r1)
      DAMPOHHST110=DAMPOHHST_f31(RG1,10,r1)
      
      
      
      RETURN
      END

      

      FUNCTION VED31T_f31(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,C8OHHST11,C10OHHST11,
     1    C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
c     ******************************************************************
      VED3=
     1     -C6OHHST11*FABC/rg16*DAMPOHHST16-
     1      C8OHHST11*FABC/rg18*DAMPOHHST18-
     2      C10OHHST11*FABC/rg110*DAMPOHHST110
     3     -C6HOHPT12*FCAB/rg26*DAMPHOHPT26-
     4      C8HOHPT12*FCAB/rg28*DAMPHOHPT28-
     5      C10HOHPT12*FCAB/rg210*DAMPHOHPT210
     6     -C6HOHPT13*FBAC/rg36*DAMPHOHPT36-
     7      C8HOHPT13*FBAC/rg38*DAMPHOHPT38-
     8      C10HOHPT13*FBAC/rg310*DAMPHOHPT310
     
      VED2=(FBC-1.0D0)*DISHHST_f31(R1)+(FAB-1.0D0)*DISOHPT_f31(R2)+
     1                            (FAC-1.0D0)*DISOHPT_f31(R3)      
      
      VED31T_f31=VED3+VED2
      
      RETURN
      END

      FUNCTION DVED31TR1_f31(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,C8OHHST11,
     1    C10OHHST11,C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/dispc1t_f31/CHH(16),COHP(10)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
c     ******************************************************************
      DVED31TR1_f31=-((DC6OHHSTR1_f31(r1,cost1,DCOST1R1)*FABC+DFABCR1*
     1      C6OHHST11)*DAMPOHHST16/rg16+
     1      (DDAMPOHHSTR1_f31(RG1,6,R1,DRG1R1)*
     1      RG16-6.0D0*DRG1R1*RG15*DAMPOHHST16)/RG112*C6OHHST11*FABC)-
     1      ((DC8OHHSTR1_f31(r1,cost1,DCOST1R1)*FABC+DFABCR1*C8OHHST11)*
     1      DAMPOHHST18/rg18+
     1      (DDAMPOHHSTR1_f31(RG1,8,r1,DRG1R1)*RG18-8.0D0*
     1      DRG1R1*RG17*DAMPOHHST18)/RG116*C8OHHST11*FABC)-
     1      ((DC10OHHSTR1_f31(r1,cost1,DCOST1R1)*
     1      FABC+DFABCR1*C10OHHST11)*
     1      DAMPOHHST110/rg110+
     1      (DDAMPOHHSTR1_f31(RG1,10,R1,DRG1R1)*RG110-
     1      10.0D0*DRG1R1*RG19*DAMPOHHST110)/RG120*C10OHHST11*FABC)-
     1      ((DC6HOHPTR1_f31(r2,cost2,DCOST2R1)*FCAB+DFCABR1*C6HOHPT12)*
     1      DAMPHOHPT26/rg26+(DDAMPHOHPTR1_f31(RG2,6,R2,DRG2R1)*
     1      RG26-6.0D0*
     1      DRG2R1*RG25*DAMPHOHPT26)/RG212*C6HOHPT12*FCAB)-
     1      ((DC8HOHPTR1_f31(r2,cost2,DCOST2R1)*FCAB+DFCABR1*C8HOHPT12)*
     1      DAMPHOHPT28/rg28+(DDAMPHOHPTR1_f31(RG2,8,R2,DRG2R1)*
     1      RG28-8.0D0*
     1      DRG2R1*RG27*DAMPHOHPT28)/RG216*C8HOHPT12*FCAB)-
     1      ((DC10HOHPTR1_f31(r2,cost2,DCOST2R1)*
     1      FCAB+DFCABR1*C10HOHPT12)*
     1      DAMPHOHPT210/rg210+(DDAMPHOHPTR1_f31(RG2,10,R2,DRG2R1)
     1      *RG210-
     1      10.0D0*DRG2R1*RG29*DAMPHOHPT210)/RG220*C10HOHPT12*FCAB)-
     1      ((DC6HOHPTR1_f31(r3,cost3,DCOST3R1)*FBAC+DFBACR1*C6HOHPT13)*
     1      DAMPHOHPT36/rg36+(DDAMPHOHPTR1_f31(RG3,6,R3,DRG3R1)*RG36-
     1      6.0D0*DRG3R1*RG35*DAMPHOHPT36)/RG312*C6HOHPT13*FBAC)-
     1      ((DC8HOHPTR1_f31(r3,cost3,DCOST3R1)*FBAC+DFBACR1*C8HOHPT13)*
     1      DAMPHOHPT38/rg38+(DDAMPHOHPTR1_f31(RG3,8,R3,DRG3R1)*RG38-
     1      8.0D0*DRG3R1*RG37*DAMPHOHPT38)/RG316*C8HOHPT13*FBAC)-
     1      ((DC10HOHPTR1_f31(r3,cost3,DCOST3R1)*
     1      FBAC+DFBACR1*C10HOHPT13)*
     1      DAMPHOHPT310/rg310+(DDAMPHOHPTR1_f31(RG3,10,R3,DRG3R1)*
     1      RG310-
     1    10.0D0*DRG3R1*RG39*DAMPHOHPT310)/RG320*C10HOHPT13*FBAC)
     
CC      VED2=(FBC-1.0D0)*DISHHST_f31(R1)+(FAB-1.0D0)*DISOHPT_f31(R2)+
CC     1                            (FAC-1.0D0)*DISOHPT_f31(R3)

      DVED21R1=DFBCR1*DISHHST_f31(R1)+(FBC-1.0D0)*DDISHHSTT_f31(R1)+
     1        DFABR1*DISOHPT_f31(R2)+
     2        DFACR1*DISOHPT_f31(R3)

      DVED31TR1_f31=DVED31TR1_f31+DVED21R1
    
      RETURN
      END

      FUNCTION DVED31TR2_f31(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,
     1    C8OHHST11,C10OHHST11,C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/dispc1t_f31/CHH(16),COHP(10)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
C     ******************************************************************
      DVED31TR2_f31=-((DC6OHHSTR2_f31(r1,cost1,DCOST1R2)*FABC+DFABCR2*
     1      C6OHHST11)*DAMPOHHST16/rg16+
     1      (DDAMPOHHSTR2_f31(RG1,6,R1,DRG1R2)*
     1      RG16-6.0D0*DRG1R2*RG15*DAMPOHHST16)/RG112*C6OHHST11*FABC)-
     1      ((DC8OHHSTR2_f31(r1,cost1,DCOST1R2)*FABC+DFABCR2*C8OHHST11)*
     1      DAMPOHHST18/rg18+
     1      (DDAMPOHHSTR2_f31(RG1,8,R1,DRG1R2)*RG18-8.0D0*
     1      DRG1R2*RG17*DAMPOHHST18)/RG116*C8OHHST11*FABC)-
     1      ((DC10OHHSTR2_f31(r1,cost1,DCOST1R2)*
     1      FABC+DFABCR2*C10OHHST11)*
     1      DAMPOHHST110/rg110+
     1      (DDAMPOHHSTR2_f31(RG1,10,R1,DRG1R2)*RG110-
     1      10.0D0*DRG1R2*RG19*DAMPOHHST110)/RG120*C10OHHST11*FABC)-
     1      ((DC6HOHPTRAA_f31(r2,cost2,DCOST2R2)*FCAB+
     1      DFCABR2*C6HOHPT12)*
     1      DAMPHOHPT26/rg26+(DDAMPHOHPTRAA_f31(RG2,6,R2,DRG2R2)*
     1      RG26-6.0D0*
     1      DRG2R2*RG25*DAMPHOHPT26)/RG212*C6HOHPT12*FCAB)-
     1      ((DC8HOHPTRAA_f31(r2,cost2,DCOST2R2)*
     1      FCAB+DFCABR2*C8HOHPT12)*
     1      DAMPHOHPT28/rg28+
     1      (DDAMPHOHPTRAA_f31(RG2,8,R2,DRG2R2)*RG28-8.0D0*
     1      DRG2R2*RG27*DAMPHOHPT28)/RG216*C8HOHPT12*FCAB)-
     1      ((DC10HOHPTRAA_f31(r2,cost2,DCOST2R2)*FCAB+
     1      DFCABR2*C10HOHPT12)*
     1      DAMPHOHPT210/rg210+(DDAMPHOHPTRAA_f31(RG2,10,R2,DRG2R2)*
     1      RG210-
     1      10.0D0*DRG2R2*RG29*DAMPHOHPT210)/RG220*C10HOHPT12*FCAB)-
     1      ((DC6HOHPTRAB_f31(r3,cost3,DCOST3R2)*
     1      FBAC+DFBACR2*C6HOHPT13)*
     1      DAMPHOHPT36/rg36+(DDAMPHOHPTRAB_f31(RG3,6,R3,DRG3R2)*
     1      RG36-6.0D0*
     1      DRG3R2*RG35*DAMPHOHPT36)/RG312*C6HOHPT13*FBAC)-
     1      ((DC8HOHPTRAB_f31(r3,cost3,DCOST3R2)*
     1      FBAC+DFBACR2*C8HOHPT13)*
     1      DAMPHOHPT38/rg38+(DDAMPHOHPTRAB_f31(RG3,8,R3,DRG3R2)*RG38-
     1      8.0D0*
     1      DRG3R2*RG37*DAMPHOHPT38)/RG316*C8HOHPT13*FBAC)-
     1      ((DC10HOHPTRAB_f31(r3,cost3,DCOST3R2)*
     1      FBAC+DFBACR2*C10HOHPT13)*
     1      DAMPHOHPT310/rg310+(DDAMPHOHPTRAB_f31(RG3,10,R3,DRG3R2)*
     1      RG310-
     1      10.0D0*DRG3R2*RG39*DAMPHOHPT310)/RG320*C10HOHPT13*FBAC)

CC      VED2=(FBC-1.0D0)*DISHHST_f31(R1)+(FAB-1.0D0)*DISOHPT_f31(R2)+
CC     1                            (FAC-1.0D0)*DISOHPT_f31(R3)

      DVED21R2=DFBCR2*DISHHST_f31(R1)+
     1        DFABR2*DISOHPT_f31(R2)+(FAB-1.0D0)*
     1        DDISPOHT_f31(R2,COHP(6),
     2        COHP(8),COHP(10),R0OHP,RMOHP)+
     2        DFACR2*DISOHPT_f31(R3)
     
      DVED31TR2_f31=DVED31TR2_f31+DVED21R2
      RETURN
      END

      FUNCTION DVED31TR3_f31(R1,R2,R3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,
     1    C8OHHST11,C10OHHST11,C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
      COMMON/dispc1t_f31/CHH(16),COHP(10)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
C     ******************************************************************
      DVED31TR3_f31=-((DC6OHHSTR3_f31(r1,cost1,DCOST1R3)*FABC+DFABCR3*
     1      C6OHHST11)*DAMPOHHST16/rg16+
     1      (DDAMPOHHSTR3_f31(RG1,6,R1,DRG1R3)*
     1      RG16-6.0D0*DRG1R3*RG15*DAMPOHHST16)/RG112*C6OHHST11*FABC)-
     1      ((DC8OHHSTR3_f31(r1,cost1,DCOST1R3)*FABC+DFABCR3*C8OHHST11)*
     1      DAMPOHHST18/rg18+(DDAMPOHHSTR3_f31(RG1,8,R1,DRG1R3)*
     1      RG18-8.0D0*
     1      DRG1R3*RG17*DAMPOHHST18)/RG116*C8OHHST11*FABC)-
     1      ((DC10OHHSTR3_f31(r1,cost1,DCOST1R3)*
     1      FABC+DFABCR3*C10OHHST11)*
     1      DAMPOHHST110/rg110+
     1      (DDAMPOHHSTR3_f31(RG1,10,R1,DRG1R3)*RG110-
     1      10.0D0*DRG1R3*RG19*DAMPOHHST110)/RG120*C10OHHST11*FABC)-
     1      ((DC6HOHPTRAB_f31(r2,cost2,DCOST2R3)*
     1      FCAB+DFCABR3*C6HOHPT12)*
     1      DAMPHOHPT26/rg26+
     1      (DDAMPHOHPTRAB_f31(RG2,6,R2,DRG2R3)*RG26-6.0D0*
     1      DRG2R3*RG25*DAMPHOHPT26)/RG212*C6HOHPT12*FCAB)-
     1      ((DC8HOHPTRAB_f31(r2,cost2,DCOST2R3)*FCAB+
     1      DFCABR3*C8HOHPT12)*
     1      DAMPHOHPT28/rg28+
     1      (DDAMPHOHPTRAB_f31(RG2,8,R2,DRG2R3)*RG28-8.0D0*
     1      DRG2R3*RG27*DAMPHOHPT28)/RG216*C8HOHPT12*FCAB)-
     1      ((DC10HOHPTRAB_f31(r2,cost2,DCOST2R3)*
     1      FCAB+DFCABR3*C10HOHPT12)*
     1      DAMPHOHPT210/rg210+
     1      (DDAMPHOHPTRAB_f31(RG2,10,R2,DRG2R3)*RG210-
     1      10.0D0*DRG2R3*RG29*DAMPHOHPT210)/RG220*C10HOHPT12*FCAB)-
     1      ((DC6HOHPTRAA_f31(r3,cost3,DCOST3R3)*
     1      FBAC+DFBACR3*C6HOHPT13)*
     1      DAMPHOHPT36/rg36+
     1      (DDAMPHOHPTRAA_f31(RG3,6,R3,DRG3R3)*RG36-6.0D0*
     1      DRG3R3*RG35*DAMPHOHPT36)/RG312*C6HOHPT13*FBAC)-
     1      ((DC8HOHPTRAA_f31(r3,cost3,DCOST3R3)*FBAC+
     1      DFBACR3*C8HOHPT13)*
     1      DAMPHOHPT38/rg38+
     1      (DDAMPHOHPTRAA_f31(RG3,8,R3,DRG3R3)*RG38-8.0D0*
     1      DRG3R3*RG37*DAMPHOHPT38)/RG316*C8HOHPT13*FBAC)-
     1      ((DC10HOHPTRAA_f31(r3,cost3,DCOST3R3)*
     1      FBAC+DFBACR3*C10HOHPT13)*
     1      DAMPHOHPT310/rg310+
     1      (DDAMPHOHPTRAA_f31(RG3,10,R3,DRG3R3)*RG310-
     1      10.0D0*DRG3R3*RG39*DAMPHOHPT310)/RG320*C10HOHPT13*FBAC)
     
CC      VED2=(FBC-1.0D0)*DISHHST_f31(R1)+(FAB-1.0D0)*DISOHPT_f31(R2)+
CC     1                            (FAC-1.0D0)*DISOHPT_f31(R3)

      DVED21R3=DFBCR3*DISHHST_f31(R1)+
     1        DFABR3*DISOHPT_f31(R2)+
     2        DFACR3*DISOHPT_f31(R3)+
     1        (FAC-1.0D0)*DDISPOHT_f31(R3,COHP(6),
     2          COHP(8),COHP(10),R0OHP,RMOHP)
     
      DVED31TR3_f31=DVED31TR3_f31+DVED21R3

      RETURN
      END
      
                  
      FUNCTION DAMPHOHPT_f31(R,N,X)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     ******************************************************************
c     calculado a partir das dist. de eq.(POL94:7651) 
c      RMHOHP=2.216766092d0
      RMHOHP=7.70d0
c      RMHOHP=7.40d0   
      RO=0.5D0*(RMHOHP+2.5D0*R0HOHPCT_f31(x))
      DAMPHOHPT_f31=(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO**2))**N      
      RETURN
      END 
      
      FUNCTION DDAMPHOHPTR1_f31(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     ******************************************************************
c      RMHOHP=2.216766092d0
      RMHOHP=7.70d0
c      RMHOHP=7.40d0
      RO=0.5D0*(RMHOHP+2.5D0*R0HOHPCT_f31(X))
      RO2=RO*RO
      DROR1=0.0D0
      DDAMPHOHPTR1_f31=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR1*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR1*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPHOHPTRAA_f31(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     ******************************************************************
c      RMHOHP=2.216766092d0
      RMHOHP=7.70d0
c      RMHOHP=7.40d0
      RO=0.5D0*(RMHOHP+2.5D0*R0HOHPCT_f31(x))
      RO2=RO*RO
      DROR23=0.5D0*2.5D0*DR0HOHPCTR23_f31(x)
      DDAMPHOHPTRAA_f31=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR23*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR23*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPHOHPTRAB_f31(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     ******************************************************************
c      RMHOHP=2.216766092d0
      RMHOHP=7.70d0
c      RMHOHP=7.40d0
      RO=0.5D0*(RMHOHP+2.5D0*R0HOHPCT_f31(x))
      RO2=RO*RO
      DROR23=0.0D0
      DDAMPHOHPTRAB_f31=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR23*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR23*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END  
                 
      FUNCTION DAMPOHHST_f31(R,N,X)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     ******************************************************************
c      RMOHHS=1.1077763d0
      RMOHHS=6.80d0
c      RMOHHS=7.40d0 
      RO=0.5D0*(RMOHHS+2.5D0*R0OHHSCT_f31(X))
      DAMPOHHST_f31=(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO**2))**N
      RETURN
      END 
      
      FUNCTION DDAMPOHHSTR1_f31(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     ******************************************************************
c      RMOHHS=1.1077763d0
      RMOHHS=6.80d0
c      RMOHHS=7.40d0 
      RO=0.5D0*(RMOHHS+2.5D0*R0OHHSCT_f31(X))
      RO2=RO*RO
      DROR1=0.5D0*2.5D0*DR0OHHSCTR1_f31(X)
      DDAMPOHHSTR1_f31=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR1*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR1*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPOHHSTR2_f31(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     ******************************************************************
c      RMOHHS=1.1077763d0
      RMOHHS=6.80d0
c      RMOHHS=7.40d0 

      RO=0.5D0*(RMOHHS+2.5D0*R0OHHSCT_f31(X))
      RO2=RO*RO
      DROR2=0.0D0
      DDAMPOHHSTR2_f31=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR2*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR2*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
      
      FUNCTION DDAMPOHHSTR3_f31(R,N,X,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     ******************************************************************
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c     ******************************************************************
c      RMOHHS=1.1077763d0
      RMOHHS=6.80d0
c      RMOHHS=7.40d0 

      RO=0.5D0*(RMOHHS+2.5D0*R0OHHSCT_f31(X))
      RO2=RO*RO
      DROR3=0.0D0
      DDAMPOHHSTR3_f31=N*(1.0D0-EXP(-A(N)*R/RO-B(N)*R**2/RO2))**(N-1)
     1 *(A(N)*(D*RO-DROR3*R)/RO2+2.0D0*B(N)*(R/RO)*(D*RO-DROR3*R)
     1          /RO2)*EXP(-A(N)*R/RO-B(N)*R**2/RO2)    
      RETURN
      END 
                  
      FUNCTION C6OHHST_f31(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/C6OHHSTae_f31/a6s,b6s
c     ******************************************************************
      A6s=C6OHHSTPA_f31(R)
      B6s=C6OHHSTPE_f31(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C6OHHST_f31=(1.0D0/3.0D0)*
     1 (2.0D0*B6s+A6s)+(2.0D0/3.0D0)*(A6s-B6s)*PL
      RETURN
      END
      
      FUNCTION DC6OHHSTR1_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/C6OHHSTae_f31/a6s,b6s
c     ******************************************************************
      DC6PE=DC6OHHSTPER1_f31(R)
      DC6PA=DC6OHHSTPAR1_f31(R)
      DC6OHHSTR1_f31=2.0D0/3.0D0*DC6PE+1.0D0/3.0D0*DC6PA
     1        +1.0D0/3.0D0*(DC6PA-DC6PE)*(3.0D0*
     1         COST**2-1.0D0)+6.0D0*D*COST*1.0D0/3.0D0*(A6s-B6s)
      RETURN
      END
      
      FUNCTION DC6OHHSTR2_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/C6OHHSTae_f31/a6s,b6s
c     ******************************************************************
      DC6OHHSTR2_f31=2.0D0*D*COST*(A6s-B6s)
      RETURN
      END      
      
      FUNCTION DC6OHHSTR3_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/C6OHHSTae_f31/a6s,b6s
c     ******************************************************************
      DC6OHHSTR3_f31=2.0D0*D*COST*(A6s-B6s)
      RETURN
      END      
      
      FUNCTION C8OHHST_f31(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/C8OHHSTae_f31/a8s,b8s
c     ******************************************************************
      A8s=C8OHHSTPA_f31(R)
      B8s=C8OHHSTPE_f31(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C8OHHST_f31=(1.0D0/3.0D0)*(2.0D0*B8s+A8s)+
     1 (2.0D0/3.0D0)*(A8s-B8s)*PL
      RETURN
      END
      
      FUNCTION DC8OHHSTR1_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/C8OHHSTae_f31/a8s,b8s
c     ******************************************************************
      DC8PE=DC8OHHSTPER1_f31(R)
      DC8PA=DC8OHHSTPAR1_f31(R)
      DC8OHHSTR1_f31=2.0D0/3.0D0*DC8PE+1.0D0/3.0D0*DC8PA
     1        +1.0D0/3.0D0*(DC8PA-DC8PE)*(3.0D0*
     1         COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A8s-B8s)
      RETURN
      END
      
      FUNCTION DC8OHHSTR2_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/C8OHHSTae_f31/a8s,b8s
c     ******************************************************************
      DC8OHHSTR2_f31=2.0D0*D*COST*(A8s-B8s)
      RETURN
      END
      
      FUNCTION DC8OHHSTR3_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/C8OHHSTae_f31/a8s,b8s
c     ******************************************************************
      DC8OHHSTR3_f31=2.0D0*D*COST*(A8s-B8s)
      RETURN
      END
      
      FUNCTION C10OHHST_f31(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhsae_f31/a1s,b1s
c     ******************************************************************
      A1s=C10OHHSTPA_f31(R)
      B1s=C10OHHSTPE_f31(R)
      PL=0.5D0*(3.0D0*COST**2-1.0D0)
      C10OHHST_f31=(1.0D0/3.0D0)*(2.0D0*B1s+A1s)+
     1  (2.0D0/3.0D0)*(A1s-B1s)*PL
      RETURN
      END
      
      FUNCTION DC10OHHSTR1_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhsae_f31/a1s,b1s
c     ******************************************************************
      DC10PE=DC10OHHSTPER1_f31(R)
      DC10PA=DC10OHHSTPAR1_f31(R)
      DC10OHHSTR1_f31=(2.0D0/3.0D0)*DC10PE+(1.0D0/3.0D0)*
     1  DC10PA+(1.0D0/3.0D0)*(DC10PA-DC10PE)
     1   *(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(A1s-B1s)
      RETURN
      END    
      
      FUNCTION DC10OHHSTR2_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhsae_f31/a1s,b1s
c     ******************************************************************
      DC10OHHSTR2_f31=2.0D0*D*COST*(A1s-B1s)
      RETURN
      END
      
      FUNCTION DC10OHHSTR3_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common/c1ohhsae_f31/a1s,b1s
c     ******************************************************************
      DC10OHHSTR3_f31=2.0D0*D*COST*(A1s-B1s)
      RETURN
      END
                  
      FUNCTION C6HOHPT_f31(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap=C6HOHPTPA_f31(R)
      Bp=C6HOHPTPE_f31(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C6HOHPT_f31=(1.0D0/3.0D0)*(2.0D0*Bp+Ap)+(2.0D0/3.0D0)*(Ap-Bp)*PL1
      RETURN
      END
      
      FUNCTION DC6HOHPTR1_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap=C6HOHPTPA_f31(R)
      Bp=C6HOHPTPE_f31(R)
      DC6HOHPTR1_f31=2.0D0*D*COST*(Ap-Bp)    
      RETURN
      END
      
      FUNCTION DC6HOHPTRAA_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap=C6HOHPTPA_f31(R)
      Bp=C6HOHPTPE_f31(R)
      DC6PEA=DC6HOHPTPERAA_f31(R)
      DC6PAA=DC6HOHPTPARAA_f31(R)
      DC6HOHPTRAA_f31=(2.0D0/3.0D0)*DC6PEA+(1.0D0/3.0D0)*
     1   DC6PAA+(1.0D0/3.0D0)*(DC6PAA-DC6PEA
     1   )*(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(Ap-Bp)     
      RETURN
      END
      
      FUNCTION DC6HOHPTRAB_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap=C6HOHPTPA_f31(R)
      Bp=C6HOHPTPE_f31(R)
      DC6HOHPTRAB_f31=2.0D0*D*COST*(Ap-Bp)     
      RETURN
      END
      
      FUNCTION C8HOHPT_f31(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      Ap8=C8HOHPTPA_f31(R)
      Bp8=C8HOHPTPE_f31(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C8HOHPT_f31=(1.0D0/3.0D0)*(2.0D0*Bp8+
     1  Ap8)+(2.0D0/3.0D0)*(Ap8-Bp8)*PL1
      RETURN
      END
      
      FUNCTION DC8HOHPTR1_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      Ap8=C8HOHPTPA_f31(R)
      Bp8=C8HOHPTPE_f31(R)
      DC8HOHPTR1_f31=2.0D0*D*COST*(Ap8-Bp8)
      RETURN
      END
          
      FUNCTION DC8HOHPTRAA_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      Ap8=C8HOHPTPA_f31(R)
      Bp8=C8HOHPTPE_f31(R)
      DC8PEA=DC8HOHPTPERAA_f31(R)
      DC8PAA=DC8HOHPTPARAA_f31(R)      
      DC8HOHPTRAA_f31=(2.0D0/3.0D0)*DC8PEA+(1.0D0/3.0D0)*
     1   DC8PAA+(1.0D0/3.0D0)*(DC8PAA-DC8PEA)
     1   *(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*(Ap8-Bp8)  
   
      RETURN
      END
      
      FUNCTION DC8HOHPTRAB_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      Ap8=C8HOHPTPA_f31(R)
      Bp8=C8HOHPTPE_f31(R)
      DC8HOHPTRAB_f31=2.0D0*D*COST*(Ap8-Bp8)     
      RETURN
      END
           
      FUNCTION C10HOHPT_f31(R,COST)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      Ap10=C10HOHPTPA_f31(R)
      Bp10=C10HOHPTPE_f31(R)
      PL1=0.5D0*(3.0D0*COST**2-1.0D0)
      C10HOHPT_f31=(1.0D0/3.0D0)*(2.0D0*Bp10+Ap10)+
     1 (2.0D0/3.0D0)*(Ap10-Bp10)
     1         *PL1
      RETURN
      END
      
      FUNCTION DC10HOHPTR1_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      A1=C10HOHPTPA_f31(R)
      B1=C10HOHPTPE_f31(R)
      DC10HOHPTR1_f31=2.0D0*D*COST*(A1-B1)
      RETURN
      END
      
      FUNCTION DC10HOHPTRAA_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C10HOHPTPA_f31(R)
      B1=C10HOHPTPE_f31(R)
      DC10PEA=DC10HOHPTPERAA_f31(R)
      DC10PAA=DC10HOHPTPARaa_f31(R)      
      DC10HOHPTRAA_f31=(2.0D0/3.0D0)*DC10PEA+(1.0D0/3.0D0)*
     1   DC10PAA+(1.0D0/3.0D0)*(DC10PAA-DC10PEA
     1   )*(3.0D0*COST**2-1.0D0)+6.0D0*D*COST*(1.0D0/3.0D0)*
     1   (A1-B1)     
      RETURN
      END
      
      FUNCTION DC10HOHPTRAB_f31(R,COST,D)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      A1=C10HOHPTPA_f31(R)
      B1=C10HOHPTPE_f31(R)
      DC10HOHPTRAB_f31=2.0D0*D*COST*(A1-B1)     
      RETURN
      END
            
      FUNCTION C6OHHSTPA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_f31/C6HSPA,R0O,ALP,RNH
c     ******************************************************************
      alp=ALPHHSPAT_f31(R)
      R0O=R0OHHSCT_f31(R)
      RNH=RNHHST_f31(R)
      C6HSPA=1.5d0*5.393d0*alp/(sqrt(alp/RNH)+
     1         sqrt(5.393d0/2.79d0))
      C6OHHSTPA_f31=C6HSPA   
      RETURN
      END
      
      FUNCTION DC6OHHSTPAR1_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/DC6PA1_f31/DC6HSPa1
      COMMON/C6PA_f31/C6HSPA,R0O,ALP,RNH
c     ******************************************************************
      CONSTT1=1.5d0*5.393d0 
      DALP=DALPHHSPATR1_f31(r)      
      e=0.5d0/SQRT(ALP/RNH)*(dalp*rnh-DRNHHSTR1_f31(r)*alp)/rnh**2
      d=SQRT(ALP/RNH)+SQRT(5.393d0/2.79d0)
      DC6HSPa1=(CONSTT1*dalp*d-CONSTT1*alp*e)/d**2
      DC6OHHSTPAR1_f31=DC6HSPa1     
      RETURN
      END
      
      FUNCTION C8OHHSTPA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_f31/C6HSPA,R0O,ALP,RNH
c     ******************************************************************
      C8OHHSTPA_f31=1.0d0*C6HSPA*R0O**(1.57243d0)
      RETURN
      END
      
      FUNCTION DC8OHHSTPAR1_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_f31/C6HSPA,R0O,ALP,RNH
      COMMON/DC6PA1_f31/DC6HSPa1
      COMMON/DRO1_f31/DR0HSCR1
c     ****************************************************************** 
      DR0HSCR1=DR0OHHSCTR1_f31(R)
      DC8OHHSTPAR1_f31=DC6HSPA1*R0O**(1.57243D0)+1.57243D0*
     1      R0O**(0.57243D0)*DR0HSCR1*C6HSPA
      RETURN
      END

      FUNCTION C10OHHSTPA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_f31/C6HSPA,R0O,ALP,RNH
c     ******************************************************************
      C10OHHSTPA_f31=1.13178d0*C6HSPA*R0O**(1.57243d0*2.0d0)
      RETURN
      END
      
      FUNCTION DC10OHHSTPAR1_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PA_f31/C6HSPA,R0O,ALP,RNH
      COMMON/DC6PA1_f31/DC6HSPa1
      COMMON/DRO1_f31/DR0HSCR1
c     ******************************************************************
c     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      DC10OHHSTPAR1_f31=1.13178D0*DC6HSPA1*R0O**(1.57243D0*
     1     2.0D0)+(1.57243D0*2.0D0)*R0O**(1.57243D0*2.0D0-1.0D0)
     1     *DR0HSCR1*1.13178D0*C6HSPA
      RETURN
      END

      FUNCTION C6OHHSTPE_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_f31/C6HSPE,R0O,RNH,alp,DALP
c     ******************************************************************
      alp=ALPHHSPET_f31(R)
      R0O=R0OHHSCT_f31(R)
      RNH=RNHHST_f31(R)      
      dalp=DALPHHSPETR1_f31(r)      
      C6HSPE=1.5d0*5.393d0*alp/(sqrt(alp/RNH)+
     1         sqrt(5.393d0/2.79d0))
      C6OHHSTPE_f31=C6HSPE    
      RETURN
      END
      
      FUNCTION DC6OHHSTPER1_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_f31/C6HSPE,R0O,RNH,alp,DALP
      COMMON/DC6PE_f31/DC6HSPE1
c     ******************************************************************
      CONSTT1=1.5D0*5.393d0  
      e=0.5d0/SQRT(ALP/RNH)*(dalp*rnh-DRNHHSTR1_f31(r)*alp)/rnh**2
      d=SQRT(ALP/RNH)+SQRT(5.393d0/2.79d0)
      DC6HSPE1=(CONSTT1*dalp*d-CONSTT1*alp*e)/d**2
      DC6OHHSTPER1_f31=DC6HSPE1   
      RETURN
      END
      
      FUNCTION C8OHHSTPE_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_f31/C6HSPE,R0O,RNH,alp,DALP
c     ******************************************************************
c     substitui 1.54 por 1.57243
      C8OHHSTPE_f31=1.0d0*C6HSPE*R0O**1.57243d0
      RETURN
      END
      
      FUNCTION DC8OHHSTPER1_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_f31/C6HSPE,R0O,RNH,alp,DALP
      COMMON/DC6PE_f31/DC6HSPE1
      COMMON/DROE1_f31/DR0HSCE1
c     ******************************************************************
c     substitui 1.54 por 1.57243
      DR0HSCE1=DR0OHHSCTR1_f31(R)
      DC8OHHSTPER1_f31=DC6HSPE1*R0O**(1.57243D0)+1.57243D0*
     1      R0O**(1.57243D0-1.0D0)*DR0HSCE1*C6HSPE
      RETURN
      END
       
      FUNCTION C10OHHSTPE_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_f31/C6HSPE,R0O,RNH,alp,DALP
c     ******************************************************************        
c     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      C10OHHSTPE_f31=1.13178d0*C6HSPE*R0O**(1.57243d0*2.0d0)
      RETURN
      END
      
      FUNCTION DC10OHHSTPER1_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/C6PE_f31/C6HSPE,R0O,RNH,alp,DALP
      COMMON/DC6PE_f31/DC6HSPE1
      COMMON/DROE1_f31/DR0HSCE1
c     ******************************************************************
c     substitui 1.31 por 1.13178 ; 1.54 por 1.57243
      DC10OHHSTPER1_f31=1.13178D0*DC6HSPE1*R0O**(1.57243D0*
     1    2.0D0)+(1.57243D0*2.0D0)*R0O**(1.57243D0*2.0D0-1.0D0)*
     1     DR0HSCE1*1.13178D0*C6HSPE
      RETURN
      END

      FUNCTION C6HOHPTPA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      ALP=ALPOHPPAT_f31(R)
      C6HOHPTPA_f31=1.5D0*4.500D0*ALP/(SQRT(ALP/RNOHPT_f31(R))+2.25D0)
      RETURN
      END
      
      FUNCTION DC6HOHPTPARAA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      CONSTT1=1.5D0*4.500D0
      ALP=ALPOHPPAT_f31(R)
      RN=RNOHPT_f31(R)
      dalp=DALPOHPPATR23_f31(r)      
      d=SQRT(ALP/RN)+2.25D0
      e=0.5d0/SQRT(ALP/RN)*(dalp*rn-DRNOHPTR23_f31(r)*alp)/rn**2
      DC6HOHPTPARAA_f31=(CONSTT1*dalp*d-e*CONSTT1*alp)/d**2
      RETURN
      END

      FUNCTION C8HOHPTPA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      C8HOHPTPA_f31=1.0D0*C6HOHPTPA_f31(R)*R0HOHPCT_f31(R)**(1.57243D0)      
      RETURN
      END
      
      FUNCTION DC8HOHPTPARAA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      R0H=R0HOHPCT_f31(R)      
      DC8HOHPTPARAA_f31=DC6HOHPTPARAA_f31(R)*R0H**1.57243D0+1.57243D0*
     1      R0H**(1.57243D0-1.0D0)*DR0HOHPCTR23_f31(R)*C6HOHPTPA_f31(R)
      RETURN
      END

      FUNCTION C10HOHPTPA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      C10HOHPTPA_f31=1.13178D0*C6HOHPTPA_f31(R)*
     1  R0HOHPCT10_f31(R)**(2.0D0*1.57243D0)      
      RETURN
      END
      
      FUNCTION DC10HOHPTPARaa_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      R0H=R0HOHPCT10_f31(R)      
      DC10HOHPTPARaa_f31=1.13178D0*DC6HOHPTPARAA_f31(R)*R0H**(2.0D0*
     1  1.57243D0)+2.0D0*1.57243D0*R0H**(2.0D0*1.57243D0-
     1    1.0D0)*DR0HOHPCT10R23_f31(R)*1.13178D0*C6HOHPTPA_f31(R)
      RETURN
      END

      FUNCTION C6HOHPTPE_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      ALP=ALPOHPPET_f31(R)
      C6HOHPTPE_f31=1.5D0*4.500D0*ALP/(SQRT(ALP/RNOHPT_f31(R))
     1           +2.25D0)
      RETURN
      END
      
      FUNCTION DC6HOHPTPERAA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      CONSTT1=1.5D0*4.500D0
      ALP=ALPOHPPET_f31(R)
      RNO=RNOHPT_f31(R)
      dalp=DALPOHPPETR23_f31(r)      
      d=SQRT(ALP/RNO)+2.25D0
      e=0.5d0/SQRT(ALP/RNO)*(dalp*rno-DRNOHPTR23_f31(r)*alp)/rno**2
      DC6HOHPTPERAA_f31=(CONSTT1*dalp*d-e*CONSTT1*alp)/d**2
      RETURN
      END
      
      FUNCTION C8HOHPTPE_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      C8HOHPTPE_f31=1.0D0*C6HOHPTPE_f31(R)*R0HOHPCT_f31(R)**(1.57243D0)      
      RETURN
      END
      
      FUNCTION DC8HOHPTPERAA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************
      R0H=R0HOHPCT_f31(R)
      DC8HOHPTPERAA_f31=DC6HOHPTPERAA_f31(R)*R0H**1.57243D0+1.57243D0*
     1      R0H**(1.57243D0-1.0D0)*DR0HOHPCTR23_f31(R)*C6HOHPTPE_f31(R)
      RETURN
      END
      

      FUNCTION C10HOHPTPE_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************            C10HOHPTPE_f31=1.13178d0*C6HOHPTPE_f31(R)*R0HOHPCT10_f31(R)**(1.57243D0*2.0D0)      
      C10HOHPTPE_f31=1.13178d0*C6HOHPTPE_f31(R)*
     1   R0HOHPCT10_f31(R)**(1.57243D0*2.0D0)      
      RETURN
      END
      
      FUNCTION DC10HOHPTPERAA_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     ******************************************************************      
      R0H=R0HOHPCT10_f31(R)      
      DC10HOHPTPERAA_f31=1.13178D0*DC6HOHPTPERAA_f31(R)*R0H**(2.0D0*
     1  1.57243D0)+2.0D0*1.57243D0*R0H**(2.0D0*1.57243D0-
     1    1.0D0)*DR0HOHPCT10R23_f31(R)*1.13178D0*C6HOHPTPE_f31(R)
      RETURN
      END
 
      FUNCTION RNOHPT_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************      
      Ainf=3.52d0
      re=1.8344d0
      a=1.04d0
      b=0.584820271D0
      c=0.584820271D0/1.04D0   
      RNOHPT_f31=ainf+(a+b*(R-re))*exp(-c*(R-re))
      RETURN
      END 
      
      FUNCTION DRNOHPTR23_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************      
      Ainf=3.52d0
      re=1.8344d0
      a=1.04d0
      b=0.584820271D0
      c=0.584820271D0/1.04D0   
      DRNOHPTR23_f31=b*exp(-c*(R-re))-c*(a+b*(R-re))*exp(-c*(R-re))
      RETURN
      END 
      
      FUNCTION RNHHST_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      Ainf=1.78d0
      re=1.449d0
      a=-0.1478d0
      c=0.60634691d0          
      RNHHST_f31=ainf+a*exp(-c*(R-re))
      RETURN
      END
      
      FUNCTION DRNHHSTR1_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      Ainf=1.78d0
      re=1.449d0
      a=-0.1478d0
      c=0.60634691d0          
      DRNHHSTR1_f31=-a*C*exp(-c*(R-re))
      RETURN
      END 

      FUNCTION ALPHHSPAT_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7/9.0D0,6.3134223d0,-2.8677438d0,
     1 7.7770980d-1,-7.6966611d-1,2.8116781d-1,4.1567195d-3,1.5074557d2/
c     ******************************************************************
      azero=1.383d0-ainf1
      rr2=r*r
      rr3=rr2*r
      rr5=rr2*rr3
      ALPHHSPAT_f31=Ainf1+(azero+C1*r+C2*rr2+c3*rr3)*
     1   EXP(-C4*r-C5*rr2)+(1.0D0-EXP(-C6*rr5))*C7/rr3
      RETURN
      END
      
      FUNCTION DALPHHSPATR1_f31(R)
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
      DALPHHSPATR1_f31=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*r-C5*RR2)-
     1          (C4+2.0D0*C5*R)*(azero+C1*r+C2*RR2+c3*RR3)*EXP(-C4*r
     1          -C5*RR2)-3.0D0*C7/RR4+3.0D0*C7/RR4*EXP(-C6*RR5)+
     1          5.0D0*C6*C7*R*EXP(-C6*RR5)
      RETURN
      END
      
      FUNCTION ALPHHSPET_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5,C6,C7/9.0D0,4.8750306d0,-1.8684334d0,3.8
     1 222971d-1,-5.6354499d-1,2.8106756d-1,1.7861826d-3,-2.8274363d1/    
c     ******************************************************************
      azero=1.383d0-ainf1
      rr2=r*r
      rr3=rr2*r
      rr5=rr2*rr3
      ALPHHSPET_f31=Ainf1+(azero+C1*r+C2*RR2+c3*RR3)*
     1   EXP(-C4*r-C5*RR2)+(1.0d0-EXP(-C6*RR5))*C7/RR3
      RETURN
      END
      
      FUNCTION DALPHHSPETR1_f31(R)
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
      DALPHHSPETR1_f31=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*r-C5*RR2)-
     1          (C4+2.0D0*C5*R)*(azero+C1*r+C2*RR2+c3*RR3)*EXP(-C4*r
     1          -C5*RR2)-3.0D0*C7/RR4+3.0D0*C7/RR4*EXP(-C6*RR5)+
     1          5.0D0*C6*C7*R*EXP(-C6*RR5)
      RETURN
      END
      
      FUNCTION ALPOHPPAT_f31(R)
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
      ALPOHPPAT_f31=Ainf1+(Azero+C1*r+C2*rr2+C3*rr3)*
     1  EXP(-C4*r-C5*rr2-C6*rr3)+(1.0d0-EXP(-C7*rr5))*C8/rr3
      RETURN
      END 
      
      FUNCTION DALPOHPPATR23_f31(R)
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
      DALPOHPPATR23_f31=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*R-C5*RR2-
     1      C6*RR3)-(C4+2.0D0*C5*R+3.0D0*C6*RR2)*(AZERO+C1*R+C2*RR2+C3*
     1         RR3)*EXP(-C4*r-C5*RR2-C6*RR3)-3.0D0*C8/RR4
     1         +3.0D0*C8/RR4*EXP(-C7*RR5)+5.0D0*C7*C8*R*EXP(-C7*RR5)
      RETURN
      END 
      
      FUNCTION ALPOHPPET_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA AINF1,C1,C2,C3,C4,C5/9.6061296d0,4.8648692d-1,-1.3199234d-1,
     1  4.2427009d-2,2.1725130d-4,-2.5535932d1/
c     ******************************************************************
      Azero=ainf1-3.5396940d0
      ALPOHPPET_f31=Ainf1-Azero*EXP(-C1*r-C2*r**2-C3*r**3)+
     1 (1.0d0-EXP(-C4*r**5))*C5/r**3

      RETURN
      END 

      FUNCTION DALPOHPPETR23_f31(R)
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
      DALPOHPPETR23_f31=Azero*(C1+2.0D0*C2*r+
     1         3.0D0*C3*RR2)*EXP(-C1*r-C2*RR2
     1         -C3*RR3)-3.0D0*c5/RR4+3.0D0*C5/RR4*EXP(-C4*RR5)
     1         +5.0D0*C4*C5*R*EXP(-C4*RR5)
      RETURN
      END 
      
      FUNCTION R0OHHSCT_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************       
      finf=6.294894d0/6.991014d0
      fzero=5.007742471d0/5.0591314d0
      d=fzero-finf
      c=0.60634691d0      
      FCORR=finf+d*exp(-c*r)
      ALPHHSM=(1.0d0/3.0d0)*ALPHHSPAT_f31(r)+(2.0d0/3.0d0)*
     1  ALPHHSPET_f31(r)       
      R0OHHSCT_f31=2.0d0*((alphhsm**(1.0d0/3.0d0))+
     1   sqrt(2.003423d0))*FCORR
      RETURN
      END
      
      FUNCTION DR0OHHSCTR1_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************       
      finf=6.294894438d0/6.9910141d0
      fzero=5.007742471d0/5.0591314d0
      d=fzero-finf
      c=0.60634691d0      
      FCORR=finf+d*exp(-c*r)
      DFCORR=-C*D*exp(-c*r)
      ALPHHSM=(1.0d0/3.0d0)*ALPHHSPAT_f31(r)+(2.0d0/3.0d0)*
     1  ALPHHSPET_f31(r)       
      DALPHHSMR1=(1.0d0/3.0d0)*DALPHHSPATR1_f31(r)+(2.0d0/3.0d0)*
     1          DALPHHSPETR1_f31(r)
      DR0OHHSCTR1_f31=2.0D0*DFCORR*((alphhsm**(1.0d0/3.0d0))+sqrt(
     1   2.003423d0))+(1.0d0/3.0d0)*alphhsm**(-2.0d0/3.0d0)*DALPHHSMR1*
     1       2.0D0*FCORR
      RETURN
      END
      
      FUNCTION R0HOHPCT_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      alpham=(1.0d0/3.0d0)*ALPOHPPAT_f31(R)+
     1  (2.0d0/3.0d0)*ALPOHPPET_f31(R)    
      finf=6.390259718d0/7.7028477d0
      fzero=5.950512078d0/6.4733336d0
      d=fzero-finf
      c=0.584820271D0/1.04D0         
      fcorr=finf+d*exp(-c*r)      
      R0HOHPCT_f31=2.0d0*((alpham**(1.0d0/3.0d0))+
     1  (6.9282d0/4.0d0))*fcorr  
      RETURN
      END 
      
      FUNCTION DR0HOHPCTR23_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      alpham=1.0d0/3.0d0*ALPOHPPAT_f31(R)+2.0d0/3.0d0*ALPOHPPET_f31(R)    
      DALPHAM=1.0d0/3.0d0*DALPOHPPATR23_f31(R)+
     1    2.0d0/3.0d0*DALPOHPPETR23_f31(R)
      finf=6.390259718d0/7.7028477d0
      fzero=5.950512078d0/6.4733336d0
      d=fzero-finf
      c=0.584820271D0/1.04D0         
      fcorr=finf+d*exp(-c*r)      
      DR0HOHPCTR23_f31=2.0d0/3.0D0*alpham**(-2.0d0/3.0d0)*DALPHAM*fcorr-
     1    2.0D0*c*d*exp(-c*r)*(alpham**(1.0d0/3.0d0)+6.9282d0/4.0d0) 
      RETURN
      END 
    
      FUNCTION R0HOHPCT10_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************     
      ALPHAM=(1.0D0/3.0D0)*ALPOHPPAT_f31(R)+
     1   (2.0D0/3.0D0)*ALPOHPPET_f31(R)  
      finf=6.574599892d0/7.7028477d0
      fzero=5.950512078d0/6.4733336d0
      d=fzero-finf
      c=0.584820271D0/1.04D0   
      FCORR=finf+d*exp(-c*r)      
      R0HOHPCT10_f31=2.0d0*((ALPHAM**(1.0D0/3.0D0))+
     1 (6.9282D0/4.0D0))*FCORR          
      RETURN
      END 
      
      FUNCTION DR0HOHPCT10R23_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************     
      ALPHAM=(1.0D0/3.0D0)*ALPOHPPAT_f31(R)+
     1   (2.0D0/3.0D0)*ALPOHPPET_f31(R)  
      DALPHAM=1.0d0/3.0d0*DALPOHPPATR23_f31(R)+ 
     1   2.0d0/3.0d0*DALPOHPPETR23_f31(R)
      finf=6.574599892d0/7.7028477d0
      fzero=5.950512078d0/6.4733336d0
      d=fzero-finf
      c=0.584820271D0/1.04D0   
      FCORR=finf+d*exp(-c*r)      
      DR0HOHPCT10R23_f31=2.0d0/(3.0D0*alpham**(2.0d0/3.0d0))*
     1     DALPHAM*fcorr
     1    -2.0D0*c*d*(alpham**(1.0d0/3.0d0)+(6.9282d0/4.0d0))*exp(-c*r) 
      RETURN
      END 
       
      FUNCTION FDOHPT_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7/1.1014506D+00,-9.2434425D-01,3.471174
     1     2D-01,-3.6259720D-01,3.6732452D-01,6.8313966D-03,2.392610
     2     1D+00/
c     ******************************************************************
      ANT_f31=3.34634328d-6
      rr2=r*r
      rr3=rr2*r
      rr5=rr3*rr2      
      FDOHPT_f31=(ANT_f31+C1*R+C2*RR2+C3*RR3)*
     1   EXP(-c4*R-C5*RR2)+(1.0d0-EXP(-C6*
     2   RR5))*C7/RR3
      return
      end
      
      FUNCTION DFDOHPTR23_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7/1.1014506D+00,-9.2434425D-01,3.471174
     1     2D-01,-3.6259720D-01,3.6732452D-01,6.8313966D-03,2.392610
     2     1D+00/
c     ******************************************************************
      ANT_f31=3.34634328d-6
      rr2=r*r
      rr3=rr2*r
      rr4=rr3*r
      rr5=rr4*r      
      DFDOHPTR23_f31=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*
     1    EXP(-C4*R-C5*RR2)-(C4+
     1    2.0D0*C5*R)*(ANT_f31+C1*R+C2*RR2+C3*RR3)*EXP(-c4*R-C5*RR2)
     1         -3.0D0*C7/RR4+5.0D0*C6*C7*R*EXP(-C6*RR5)+3.0D0*
     1         C7/RR4*EXP(-C6*RR5)
      return
      end
     
      FUNCTION VINDHOHPT_f31(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,
     1    C8OHHST11,C10OHHST11,C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
c     ******************************************************************               
      VINDHOHPT_f31=-FDOHPT_f31(R2)**2*4.5D0*(3.0D0*COST2**2+1.0D0)*
     1          DAMPHOHPT26/(2.0D0*rG26)
     1         -FDOHPT_f31(R3)**2*4.5D0*(3.0D0*COST3**2+1.0D0)*
     1          DAMPHOHPT36/(2.0D0*rG36)
      return
      end
      
      FUNCTION DVINDHOHPTR1_f31(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,
     1    C8OHHST11,C10OHHST11,C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
c     ******************************************************************
      TA=3.0D0*COST2**2+1.0D0    
      DTAR1=-6.0D0*(R3-R1)/R2**2
      TB=3.0D0*COST3**2+1.0D0    
      DTBR1=-6.0D0*(R2-R1)/R3**2
      DFDOHPTR1=0.0D0
      xb=DAMPHOHPT26/(2.0D0*rG26) 
      xb1=DAMPHOHPT36/(2.0D0*rG36) 
      FDOHPTR2=FDOHPT_f31(R2)
      FDOHPTR3=FDOHPT_f31(R3)
      TA1=-4.5D0*FDOHPTR2**2
      TA2=-4.5D0*FDOHPTR3**2
      xa=TA1*TA  
      xa1=TA2*TB  
      DVINDHOHPTR1_f31=xb*(-4.5d0*2.0d0*FDOHPTR2*DFDOHPTR1*Ta+dTar1*TA1)
     1 +xa*(DDAMPHOHPTR1_f31(Rg2,6,r2,DRG2R1)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R1*DAMPHOHPT26)/(4.0D0*RG212)
     1         +xb1*(-4.5d0*2.0d0*FDOHPTR3*DFDOHPTR1*TB+dTBr1*TA2)
     1 +Xa1*(DDAMPHOHPTR1_f31(Rg3,6,r3,DRG3R1)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R1*DAMPHOHPT36)/(4.0D0*RG312) 
      return
      end
      
      FUNCTION DVINDHOHPTR2_f31(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,
     2    C8OHHST11,C10OHHST11,C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
c     ******************************************************************
      TA=3.0D0*COST2**2+1.0D0    
      DTAR2=-6.0D0*(R3-R1)**2/R2**3
      TB=3.0D0*COST3**2+1.0D0    
      DTBR2=6.0D0*(R2-R1)/R3**2
      FDOHPTR2=FDOHPT_f31(R2)
      FDOHPTR3=FDOHPT_f31(R3)
      xb=DAMPHOHPT26/(2.0D0*rG26) 
      xb1=DAMPHOHPT36/(2.0D0*rG36) 
      TA1=-4.5D0*FDOHPTR2**2
      TA2=-4.5D0*FDOHPTR3**2
      xa=TA1*TA  
      xa1=TA2*TB  
      DFDOHPTR23R3=0.0D0
      DVINDHOHPTR2_f31=xb*(-4.5d0*2.0d0*
     1 FDOHPTR2*DFDOHPTR23_f31(R2)*Ta+dTar2
     1 *TA1)+xa*(DDAMPHOHPTRAA_f31(Rg2,6,r2,DRG2R2)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R2*DAMPHOHPT26)/(4.0D0*RG212)
     1      +xb1*(-4.5d0*2.0d0*FDOHPTR3*DFDOHPTR23R3*TB+dTBr2*TA2)
     1 +Xa1*(DDAMPHOHPTRAB_f31(Rg3,6,r3,DRG3R2)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R2*DAMPHOHPT36)/(4.0D0*RG312)
      return
      end      
      
      FUNCTION DVINDHOHPTR3_f31(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/CVED31T_f31/C6OHHST11,
     1    C8OHHST11,C10OHHST11,C6HOHPT12,C8HOHPT12,
     1    C10HOHPT12,C6HOHPT13,C8HOHPT13,C10HOHPT13,DAMPOHHST16,
     1    DAMPOHHST18,DAMPOHHST110,DAMPHOHPT26,DAMPHOHPT28,DAMPHOHPT210,
     1    DAMPHOHPT36,DAMPHOHPT38,DAMPHOHPT310,C6AC,C8AC,C10AC,C6AB,
     1    C8AB,C10AB
     
c     ******************************************************************
      TA=3.0D0*COST2**2+1.0D0    
      DTAR3=6.0D0*(R3-R1)/R2**2
      TB=3.0D0*COST3**2+1.0D0    
      DTBR3=-6.0D0*(R2-R1)**2/R3**3
      FDOHPTR2=FDOHPT_f31(R2)
      FDOHPTR3=FDOHPT_f31(R3)
      xb=DAMPHOHPT26/(2.0D0*rG26) 
      xb1=DAMPHOHPT36/(2.0D0*rG36) 
      TA1=-4.5D0*FDOHPTR2**2
      TA2=-4.5D0*FDOHPTR3**2
      xa=TA1*TA  
      xa1=TA2*TB  
      DFDOHPTR23R2=0.0D0      
      DVINDHOHPTR3_f31=xb*(-4.5d0*2.0d0*FDOHPTR2*DFDOHPTR23R2*Ta+dTar3
     1 *TA1)+xa*(DDAMPHOHPTRAB_f31(Rg2,6,r2,DRG2R3)*2.0D0*RG26
     1 -12.0D0*RG25*DRG2R3*DAMPHOHPT26)/(4.0D0*RG212)
     1    +xb1*(-4.5d0*2.0d0*FDOHPTR3*DFDOHPTR23_f31(R3)*TB+dTBr3*TA2)
     1 +Xa1*(DDAMPHOHPTRAA_f31(Rg3,6,r3,DRG3R3)*2.0D0*RG36
     1 -12.0D0*RG35*DRG3R3*DAMPHOHPT36)/(4.0D0*RG312)
      return
      end
    
      FUNCTION FQHHST_f31(R)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7,C8/9.2292863d-02,-1.0187873d-02,1.03075
     1     21d-02,-1.0635640d0,2.8476217d-01,1.1961650d-08,1.3473360
     2     d-02,3.9037435d+00/
c     ******************************************************************
      RR2=R*R
      RR3=RR2*R
      RR5=RR3*RR2
      FQHHST_f31=(C1*R+C2*RR2+C3*RR3)*EXP(-C4*R
     1      -C5*RR2-C6*RR3)+(1.0d0-EXP(-C7*RR5))*C8/RR3      
      return
      end
      
      FUNCTION DFQHHSTR1_f31(R)
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
      DFQHHSTR1_f31=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*
     1        EXP(-C4*R-C5*RR2-C6*RR3)
     1        -(C4+2.0D0*C5*R+3.0D0*C6*RR2)*(C1*R+C2*RR2+C3*RR3)*EXP
     1        (-C4*R-C5*RR2-C6*RR3)-3.0D0*C8/RR4+3.0D0*C8/RR4*EXP(
     1        -C7*RR5)+5.0D0*C7*C8*R*EXP(-C7*RR5)
      return
      end
      
      FUNCTION VEOHHST_f31(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
      DATA b0,b2,b4,b6/5.219002433D0,2.539416058d0,0.279854015d0,
     1                 -0.0406812652D0/
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
c     ******************************************************************
c     o valor do quadrupolo do O(1D)=1.233793Buckingham=0.917288208au
c     o valor do quadrupolo foi alterado para -0.992148212d0 
c     correspondendo ao do O(3p)

      at=b0+b2*cos2t1+b4*cos4t1+b6*cos6t1
      VEOHHST_f31=0.75d0*FQHHST_f31(R1)*(-0.992148212d0)*at*DAMPOHHST5
      return
      end
      
      FUNCTION DVEOHHSTR1_f31(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
      DATA b0,b2,b4,b6/5.219002433D0,2.539416058d0,0.279854015d0,
     1                 -0.0406812652D0/
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/DVEOHP_f31/AT,DAT,FQHHSTR1
c     ******************************************************************
      at=b0+b2*cos2t1+b4*cos4t1+b6*cos6t1
      DAT=B2*dcos2t+b4*dcos4t+b6*dcos6t
      DATR1=DAT*dcost1r1
      CONSTTV=0.75d0*(-0.992148212d0)
      FQHHSTR1=FQHHST_f31(R1)
      DVEOHHSTR1_f31=(DFQHHSTR1_f31(R1)*AT+DATR1*FQHHSTR1)*
     1     CONSTTV*DAMPOHHST5+
     1     (DDAMPOHHSTR1_f31(Rg1,5,r1,drg1r1)*rg15-5.0D0*DRG1R1*
     1     RG14*DAMPOHHST5*RG15)/(rg110)*at*CONSTTV*FQHHSTR1
      return
      end
      
      FUNCTION DVEOHHSTR2_f31(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c      DATA b0,b2,b4,b6/5.219002433D0,2.539416058d0,0.279854015d0,
C     1                 -0.0406812652D0/
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/DVEOHP_f31/AT,DAT,FQHHSTR1
c     ******************************************************************
c     o valor do quadrupolo do O(1D)=1.233793Buckingham=0.917288208au
      DATR2=DAT*dcost1r2
      CONSTTV=0.75d0*(-0.992148212d0)
      DVEOHHSTR2_f31=(DATR2*FQHHSTR1)*CONSTTV*DAMPOHHST5+
     1        (DDAMPOHHSTR2_f31(Rg1,5,r1,DRG1R2)*rg15-5.0D0*DRG1R2*
     1        RG14*DAMPOHHST5*Rg15)/(rg110)*at*CONSTTV*FQHHSTR1
      return
      end
      
      FUNCTION DVEOHHSTR3_f31(r1,r2,r3)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c      DATA b0,b2,b4,b6/5.219002433D0,2.539416058d0,0.279854015d0,
c     1                 -0.0406812652D0/
      COMMON/CVED_f31/COSB,COSC,RG1,RG2,RG3,COST1,
     1  COST2,COST3,FAB,FAC,FBC,
     1  FABC,FCAB,FBAC,CONSTT,RG14,RG15,RG16,RG17,RG18,RG19,RG110,
     1  RG112,RG116,RG120,RG25,RG26,RG27,RG28,RG29,RG210,RG212,
     1  RG216,RG220,RG35,RG36,RG37,RG38,RG39,RG310,RG312,RG316,
     1  RG320,DRG1R1,DRG2R1,DRG3R1,DRG1R2,DRG2R2,DRG3R2,DRG1R3,
     1  DRG2R3,DRG3R3,DCOST1R1,DCOST2R1,DCOST3R1,DCOST1R2,
     1  DCOST2R2,DCOST3R2,DCOST1R3,DCOST2R3,DCOST3R3,DFABCR1,
     1  DFCABR1,DFBACR1,DFABCR2,DFCABR2,DFBACR2,DFABCR3,DFCABR3,
     1  DFBACR3,C6BC,C8BC,C10BC,cost12,cost13,cost14,cost15,cost16,
     1  cos2t1,cos4t1,cos6t1,dcos2t,dcos4t,dcos6t,DAMPOHHST5,
     1  DFABR1,DFACR1,DFBCR1,DFABR2,DFACR2,DFBCR2,DFABR3,DFACR3,DFBCR3
      COMMON/DVEOHP_f31/AT,DAT,FQHHSTR1
c     ******************************************************************
c     o valor do quadrupolo do O(1D)=1.233793Buckingham=0.917288208au
      DATR3=DAT*dcost1r3
      CONSTTV=0.75d0*(-0.992148212d0)
      DVEOHHSTR3_f31=(DATR3*FQHHSTR1)*CONSTTV*DAMPOHHST5+
     1    (DDAMPOHHSTR3_f31(Rg1,5,r1,DRG1R3)*rg15-5.0D0*DRG1R3*
     1    RG14*DAMPOHHST5*RG15)/(rg110)*at*CONSTTV*FQHHSTR1
      return
      end
      
      FUNCTION RL_f31(r1E,r2E,r3E)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      R1=R1E+1.0D-8
      R2=R2E+1.0D-8
      R3=R3E+1.0D-8
      call GEOMRLT_f31(r1,r2,r3)
      RL_f31=VED31T_f31(r1,r2,r3)+VEOHHST_f31(r1,r2,r3)+
     1   VINDHOHPT_f31(r1,r2,r3)
      return
      end

      FUNCTION RLCAS_f31(r1E,r2E,r3E)
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      R1=R1E+1.0D-8
      R2=R2E+1.0D-8
      R3=R3E+1.0D-8
      call GEOMRLT_f31(r1,r2,r3)
      RLCAS_f31=VEOHHST_f31(r1,r2,r3)
      return
      end
      
      SUBROUTINE DRL_f31(R1E,R2E,R3E,DER1,DER2,DER3)
c     ******************************************************************
c     TO COMPUTE THE DERIVATIVES OF THE LONG RANGE TERM
c     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
c     ******************************************************************
      R1=R1E+1.0D-8
      R2=R2E+1.0D-8
      R3=R3E+1.0D-8
      call GEOMRLT_f31(r1,r2,r3)

      DER1=DVED31TR1_f31(r1,r2,r3)+DVEOHHSTR1_f31(r1,r2,r3)+
     1      DVINDHOHPTR1_f31(r1,r2,r3)

      DER2=DVED31TR2_f31(r1,r2,r3)+DVEOHHSTR2_f31(r1,r2,r3)+
     1      DVINDHOHPTR2_f31(r1,r2,r3)

      DER3=DVED31TR3_f31(r1,r2,r3)+DVEOHHSTR3_f31(r1,r2,r3)+
     1      DVINDHOHPTR3_f31(r1,r2,r3)
      RETURN
      END
     

      SUBROUTINE CODISPT_f31(R0,C6,C8,C10)
c    **************************************************************
c     TO PREDICT C8 AND C10 FROM R0 AND C6
c    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      C8=R0**1.5724274D0*C6
      C10=EXP(0.12379935D0)*R0**3.1448548D0*C6
      RETURN
      END
      
      FUNCTION DISPOHT_f31(R,C6,C8,C10,R0,RM)
c    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
c    **************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      RR2=RR*RR
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR2))**6
      D8=(1.0D0-EXP(-A(8)*RR-B(8)*RR2))**8
      D10=(1.0D0-EXP(-A(10)*RR-B(10)*RR2))**10
      DISPOHT_f31=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

      FUNCTION DDISPOHT_f31(R,C6,C8,C10,R0,RM)
c    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
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
      DDISPOHT_f31=6.0D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-A(6)-
     1  2.0D0*B(6)*RR)*DRR)+
     2 8.0D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-A(8)-2.0D0*B(8)
     3      *RR)*DRR)+
     4 10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-A(10)-2.0D0*
     5      B(10)*RR)*DRR)
      RETURN
      END
      
      BLOCK DATA H2OStDAT
C     *************************************************************      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/CONSTT_f31/PI
      COMMON/ANT_f31/A(20)
      COMMON/BNT_f31/B(20)
      COMMON/EMIT_f31/EMITN
      COMMON/dispc1t_f31/CHH(16),COHP(10)
      COMMON/refgeot_f31/R10,R20,R30
c apenas para ser utilizado na funcao
      COMMON/coefft_f31/C(127)
      COMMON/DIATDIT_f31/R0HH,RMHHS,R0OHP,RMOHP
c    **************************************************************
      DATA PI/3.141592653589793238462643D0/
c      DATA R10,R20,R30/3.0D0,3.0D0,3.0D0/

c      Geometria obtida pelo ponto cela      
c      DATA R10,R20,R30/4.5277815D0,2.2638907D0,2.2638907D0/ 
      DATA R10,R20,R30/4.5D0,3.0D0,3.0D0/
      DATA EMITN/-0.3704003D0/
      DATA R0HH,RMHHS,R0OHP,RMOHP/6.928203D0,
     1   1.40100D0,6.294894D0,1.8344D0/   
c    **************************************************************      
c     ANT E BNT OBTIDOS A PARTIR DOS COEF.(A0,A1,B0,B1) DO ARTG. VAR82:857 
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
     
      DATA COHP/5*0.0D0,10.00D0,0.0D0,180.45D0,0.0D0,3685.26D0/
C    ***************************************************************      

      DATA C/
     1    0.6525952591D+00,    0.6059213801D+00,    0.6286902475D-01,
     1    0.5382013948D+00,   -0.3596542428D-01,   -0.2855656568D-01,
     1   -0.3135159953D+00,    0.3354917216D+00,   -0.4575071108D-01,
     1   -0.8178111696D-02,   -0.6723402438D-01,   -0.3707215878D+00,
     1   -0.5635087779D-01,    0.2183034962D+00,   -0.2656352287D-01,
     1    0.1497850059D+00,   -0.6499599003D-01,   -0.1167749569D+00,
     1   -0.3941619360D+00,   -0.1388503877D+00,    0.6654681874D-01,
     1    0.1899065127D+00,    0.1301533766D+00,   -0.9702614739D-02,
     1    0.8716655111D-01,   -0.5746247403D-01,    0.1633507741D-01,
     1   -0.8809282609D-01,   -0.2633321889D+00,   -0.8468112681D-01,
     1    0.4449952690D-01,    0.1865964439D+00,    0.3596828642D-01,
     1    0.3965350663D-02,    0.5692200930D-01,   -0.7079258591D-02,
     1    0.3132496109D-01,   -0.1762698067D-01,   -0.2106706104D-01,
     1   -0.1718918022D-01,    0.8309588125D-03,   -0.3633214348D-01,
     1   -0.1152405251D+00,   -0.1143744408D-01,    0.3044705789D-01,
     1    0.1121034542D+00,    0.1474073807D-01,    0.1095018141D-01,
     1   -0.1110958764D-01,   -0.7666399845D-01,    0.1618057581D-01,
     1   -0.1274008893D-02,    0.5803431360D-02,    0.3293688787D-02,
     1    0.9745874276D-03,    0.5471927251D-02,   -0.1526060801D-01,
     1    0.1207500523D-01,   -0.1255867163D-01,    0.8122120012D-02,
     1   -0.4206922766D-02,   -0.3539766233D-02,   -0.2452018044D-01,
     1   -0.7813442331D-02,    0.1123666809D-01,   -0.3388875370D-01,
     1    0.3404878521D-01,   -0.1572427289D-01,   -0.3499345260D-01,
     1   -0.6012860371D-02,    0.2655206532D-02,    0.3253199009D-02,
     1    0.8528718643D-02,   -0.6138026337D-03,   -0.3379726639D-02,
     1    0.9481045184D-02,   -0.4437091686D-02,    0.5992085143D-02,
     1    0.3776988155D-03,    0.4563875092D-02,   -0.4412550077D-02,
     1   -0.6679641503D-02,   -0.1631485747D-01,   -0.3729271713D-03,
     1    0.3948351511D-03,   -0.2726824983D-02,    0.1339644845D-02,
     1   -0.9240515739D-02,   -0.7307921432D-03,   -0.7083880205D-02,
     1   -0.1249675257D-01,    0.3341157715D-02,    0.6183759881D-02,
     1   -0.1553550410D-02,    0.8393866566D-02,    0.1881743783D-03,
     1    0.7901600291D-03,    0.3016551340D-02,   -0.5642886937D-03,
     1    0.4185438087D-03,    0.3112883486D-02,   -0.2109989995D-02,
     1    0.1610284485D-02,   -0.1159061075D-02,    0.9290594429D-03,
     1   -0.6906487927D-03,   -0.2447727869D-02,   -0.5173593663D-02,
     1   -0.5484383992D-04,   -0.5184261919D-04,    0.5304455883D-03,
     1    0.1944119173D-02,   -0.1498803079D-02,   -0.9771064689D-03,
     1   -0.1459160168D-02,    0.2119279693D-03,   -0.1129835250D-03,
     1   -0.2429045396D-03,    0.6631102801D-04,   -0.8230573661D-03,
     1    0.8767205738D-04,   -0.6138736435D-03,    0.8261582000D-04,
     1   -0.6212609060D-03,   -0.1374132106D-02,    0.8000000000D+00,
     1    0.1100000000D+01/


           END
