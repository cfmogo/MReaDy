      subroutine pothhs_f1(rin,Vout,dVout)
      use phys_parameters
      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: Vout,dVout
      double precision :: VHHS_f1,DVHHS_f1

      double precision ::Vout_lr,dVout_lr,f1,df1,fs,dfs
       


      Vout=VHHS_f1(rin/5.2917721092d-11)*4.3597482D-18
      dVout=DVHHS_f1(rin/5.2917721092d-11)*4.3597482D-18/
     1 5.2917721092d-11

      if (rin.gt.in_bound) then      
         call pothht_f14(rin,Vout_lr,dVout_lr)      
         fs=f1(rin) 
         dfs=df1(rin)
   
         Vout=Vout*(1-fs)+Vout_lr*fs      
         dVout=dVout*(1-fs)-Vout_lr*dfs      
      end if

      end subroutine pothhs_f1




c    **************************************************************
c
C     TO COMPUTE THE HFACE FOR singlet H...H
c
C    **************************************************************
      FUNCTION VHHS_f1(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c    **************************************************************
      VHHS_f1=EHFHHS_f1(R)+DISHHS_f1(R)
      RETURN
      END

      FUNCTION EHFHHS_f1(R)
c    **************************************************************
C     TO COMPUTE THE EHF FOR  singlet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      implicit integer(i-n)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATHHS_f1/R0HH,RMHHS
      COMMON/ANHHS_f1/A(20)
      COMMON/BNHHS_f1/B(20)
!      COMMON/CEHFHS_f1/X,X2,X3,RHO,RR,GAMA,POL,D6

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
      EHFHHS_f1=-D/R*POL*EXP(-Gama*X)
      RHO=(RMHHS+2.5D0*R0HH)/2.0d0
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6  
      EXCHHS=-AGEX*R**ALPHEX*EXP(-GAMEX*R)*D6
      EHFHHS_f1=EHFHHS_f1+EXCHHS
      RETURN
      END

      FUNCTION DISHHS_f1(R)
c    **************************************************************
C     TO COMPUTE THE DISPERSION FOR singlet H...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      implicit integer(i-n)
      COMMON/DISPCHHS_f1/CHH(16)
c    **************************************************************
      R2=R*R
      R6=R2*R2*R2
      R8=R6*R2
      ADISP=-CHH(6)*DAMPHH_f1(6,R)/R6-CHH(8)*DAMPHH_f1(8,R)/R8
      DO N=10,16
        ADISP=ADISP-CHH(N)*DAMPHH_f1(N,R)/(R8*R**(N-8.0d0))
      ENDDO
      DISHHS_f1=ADISP 
      RETURN
      END

      FUNCTION DAMPHH_f1(N,R)
c    *************************************************************** 
c     CALCULATES VARANDAS-BRANDAO DAMPING FUNCTIONS
c    ***************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      implicit integer(i-n)
      COMMON/DIATHHS_f1/R0HH,RMHHS
      COMMON/ANHHS_f1/A(20)
      COMMON/BNHHS_f1/B(20)
c    ***************************************************************
      RR=2.0D0*R/(RMHHS+2.5D0*R0HH)
      DAMPHH_f1=(1.0D0-EXP(-A(N)*RR-B(N)*RR*RR))**N
      RETURN
      END

      FUNCTION DVHHS_f1(R)
C     *************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  singlet H...H
C     *************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C     *************************************************************
      DVHHS_f1=DEHFHHS_f1(R)+DDISHHS_f1(R)
      RETURN
      END

      FUNCTION DEHFHHS_f1(R)
C     *************************************************************
C     TO COMPUTE THE EHF FOR  singlet H...H
C     *************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      implicit integer(i-n)
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATHHS_f1/R0HH,RMHHS
      COMMON/ANHHS_f1/A(20)
      COMMON/BNHHS_f1/B(20)
!      COMMON/CEHFHS_f1/X,X2,X3,RHO,RR,GAMA,POL,D6

c     **************************************************************
      DATA D,ASV/0.218973D0,1.91479d0,0.646041D0,0.346414D0/
      DATA (GAM(I),i=0,2)/1.22349D0,1.04334D0,0.208477D0/
      DATA AGEX,ALPHEX,GAMEX/0.8205D0,2.5D0,2.0D0/
c     **************************************************************
      X=R-RMHHS
      X2=X*X
      X3=X2*X
      POL=(1.0D0+ASV(1)*X+ASV(2)*X2+ASV(3)*X3)
      Gama=Gam(0)*(1.0D0+Gam(1)*TANH(Gam(2)*X))      

      RHO=(RMHHS+2.5D0*R0HH)/2.0d0
      RR=R/RHO
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR**2))**6  
 

      DGama=Gam(0)*Gam(1)*Gam(2)/COSH(Gam(2)*X)**2
      DPOLC=D/R**2*POL-D/R*(ASV(1)+2*ASV(2)*X+3*ASV(3)*X2)
      DEXP=-(Gama+DGama*X)*EXP(-Gama*X)
      DEHFHHS_f1=DPOLC*EXP(-GAMA*X)-D/R*POL*DEXP
      RHH=RHO
      RR2=RR*RR
      DBASE=(A(6)+2*B(6)*RR)/RHH*EXP(-A(6)*RR-    
     1    B(6)*RR2)
      DD6=6*(1-EXP(-A(6)*RR-B(6)*RR2))**5*DBASE
      FEXC=-AGEX*R**ALPHEX*EXP(-GAMEX*R)
      DFEXC=-FEXC*GAMEX+ALPHEX*FEXC/R
      DEXCHHS=DFEXC*D6+DD6*FEXC
      DEHFHHS_f1=DEHFHHS_f1+DEXCHHS
      RETURN
      END

      FUNCTION DDISHHS_f1(R)
c    ************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPCHHS_f1/CHH(16)
      COMMON/DIATHHS_f1/R0HH,RMHHS
      COMMON/ANHHS_f1/A(20)
      COMMON/BNHHS_f1/B(20)
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
      DDISHHS_f1=DDISP
      RETURN
      END

       BLOCK DATA H2SING_f1
C     *************************************************************      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ANHHS_f1/A(20)
      COMMON/BNHHS_f1/B(20)
      COMMON/DISPCHHS_f1/CHH(16)
      COMMON/DIATHHS_f1/R0HH,RMHHS
c    **************************************************************
      DATA R0HH,RMHHS/6.928203D0,1.40100D0/   
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
      END
