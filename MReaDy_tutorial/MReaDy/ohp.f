      subroutine potohp_f3(rin,Vout,dVout)
      use phys_parameters
      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: Vout,dVout
      double precision :: VOHP_f3,DVOHP_f3
      
      double precision ::Vout_lr,dVout_lr,f1,df1,fs,dfs

      Vout=VOHP_f3(rin/5.2917721092d-11)*4.3597482D-18
      dVout=DVOHP_f3(rin/5.2917721092d-11)*4.3597482D-18/
     1 5.2917721092d-11

      if (rin.gt.in_bound) then      
         call pothoq_f17(rin,Vout_lr,dVout_lr)      
         fs=f1(rin) 
         dfs=df1(rin)
   
         Vout=Vout*(1-fs)+Vout_lr*fs      
         dVout=dVout*(1-fs)-Vout_lr*dfs      
      end if

      end subroutine potohp_f3



      FUNCTION VOHP_f3(R)
c    **************************************************************
C     TO COMPUTE THE HFACE FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c    **************************************************************
      VOHP_f3=EHFOHP_f3(R)+DISOHP_f3(R)
      RETURN
      END

      FUNCTION EHFOHP_f3(R)
c    **************************************************************
C     TO COMPUTE THE EHF FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      implicit integer(i-n)      
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DIATOHP_f3/R0OHP,RMOHP
      COMMON/ANOHP_f3/A(20)
      COMMON/BNOHP_f3/B(20)
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
      EHFOHP_f3=EHF+EX
      RETURN
      END


      FUNCTION DISOHP_f3(R)
c    **************************************************************
C     TO COMPUTE THE DISPERSION FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPCOHP_f3/COHP(10)
      COMMON/DIATOHP_f3/R0OHP,RMOHP
c    **************************************************************
      DISOHP_f3=DISPOH_f3(R,COHP(6),COHP(8),COHP(10),R0OHP,RMOHP)
      RETURN
      END

      FUNCTION DVOHP_f3(R)
c    **************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  O...H
C    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      implicit integer(i-n) 
      DIMENSION ASV(3),GAM(0:2)
      COMMON/DISPCOHP_f3/COHP(10)
      COMMON/DIATOHP_f3/R0OHP,RMOHP
      COMMON/ANOHP_f3/A(20)
      COMMON/BNOHP_f3/B(20)
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
     1  *R**ALPHEX)*EXP(-GAMEX*R)+GAMEX*AGEX*R**ALPHEX*
     2    (1.0D0+A1EX*R)*EXP(-GAMEX*R)
      DEX=EX/D6*DD6+DEXND*D6
c     Dispersao
      DERDISP=DDISPOH_f3(R,COHP(6),COHP(8),COHP(10),R0OHP,RMOHP)
c     Total
      DVOHP_f3=DEHF+DEX+DERDISP
      RETURN
      END



      FUNCTION DISPOH_f3(R,C6,C8,C10,R0,RM)
c    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ANOHP_f3/A(20)
      COMMON/BNOHP_f3/B(20)
c    **************************************************************
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      RR2=RR*RR
      D6=(1.0D0-EXP(-A(6)*RR-B(6)*RR2))**6
      D8=(1.0D0-EXP(-A(8)*RR-B(8)*RR2))**8
      D10=(1.0D0-EXP(-A(10)*RR-B(10)*RR2))**10
      DISPOH_f3=-C6/R6*D6-C8/R8*D8-C10/R10*D10
      RETURN
      END

      FUNCTION DDISPOH_f3(R,C6,C8,C10,R0,RM)
c    **************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ANOHP_f3/A(20)
      COMMON/BNOHP_f3/B(20)
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
      DDISPOH_f3=6.0D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-A(6)-
     1  2.0D0*B(6)*RR)*DRR)+
     2 8.0D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-A(8)-2.0D0*B(8)
     3      *RR)*DRR)+
     4 10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-A(10)-2.0D0*
     5      B(10)*RR)*DRR)
      RETURN
      END
      
      BLOCK DATA OHP_f3
C     *************************************************************      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ANOHP_f3/A(20)
      COMMON/BNOHP_f3/B(20)
      COMMON/DISPCOHP_f3/COHP(10)
      COMMON/DIATOHP_f3/R0OHP,RMOHP
c    **************************************************************
      DATA R0OHP,RMOHP/6.294894D0,1.8344D0/   
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
      DATA COHP/5*0.0D0,10.00D0,0.0D0,180.45D0,0.0D0,3685.26D0/
      
      END
