      subroutine potoop_f2(rin,Vout,dVout)
      use phys_parameters
      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: Vout,dVout
      double precision :: VOO_f2,DVOO_f2
  
      double precision ::Vout_lr,dVout_lr,f1,df1,fs,dfs

      Vout=VOO_f2(rin/5.2917721092d-11)*4.3597482D-18
      dVout=DVOO_f2(rin/5.2917721092d-11)*4.3597482D-18/
     1 5.2917721092d-11

      if (rin.gt.in_bound) then      
         call poto2q_f18(rin,Vout_lr,dVout_lr)      
         fs=f1(rin) 
         dfs=df1(rin)
   
         Vout=Vout*(1-fs)+Vout_lr*fs      
         dVout=dVout*(1-fs)-Vout_lr*dfs      
      end if

      end subroutine potoop_f2



C     ****************************************************************
C
C     TO COMPUTE THE HFACE FOR  O...O
C
C     ****************************************************************
      FUNCTION VOO_f2(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      VOO_f2=EHFOO_f2(R)+DISOO_f2(R)

      END


      
C     ****************************************************************
C
C     TO COMPUTE THE EHF FOR  O...O
C
C     ****************************************************************
      FUNCTION EHFOO_f2(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DIATDIOO_f2/R0OO,RMOO
      DATA D,ASV/0.14291202D0,3.6445906D0,3.9281238D0,2.0986689D0,
     1           3.3522498D0/

      X=R-RMOO
      R2=X*X
      R3=R2*X
      EHFOO_f2=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*DEXP(-ASV(4)*X)

      END





      
C     ****************************************************************
C
C     TO COMPUTE THE DISPERSION FOR  O...O
C
C     ****************************************************************
      FUNCTION DISOO_f2(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPCOO_f2/COO(10)
      COMMON/DIATDIOO_f2/R0OO,RMOO
      
      
      DISOO_f2=DISP_f2(R,COO(6),COO(8),COO(10),R0OO,RMOO)

      END

            
C     ***************************************************************      
      FUNCTION DISP_f2(R,C6,C8,C10,R0,RM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DAMPC_f2/ADAMP(10),BDAMP(10)
      
      
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      D6=(1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6
      D8=(1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8
      D10=(1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10
      DISP_f2=-C6/R6*D6-C8/R8*D8-C10/R10*D10

      END
      
C     ****************************************************************
C     TO COMPUTE THE DER. OF THE HFACE FOR  O...O
C     ****************************************************************
      FUNCTION DVOO_f2(R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ASV(4)
      COMMON/DISPCOO_f2/COO(10)
      COMMON/DIATDIOO_f2/R0OO,RMOO
      DATA D,ASV/0.14291202D0,3.6445906D0,3.9281238D0,2.0986689D0,
     1           3.3522498D0/
      
      X=R-RMOO
      R2=X*X
      R3=R2*X
      POL=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)
      DPOL=-D*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*R2)
      POT=EXP(-ASV(4)*X)
      DVOO_f2=-ASV(4)*POT*POL+DPOL*POT+DDISP_f2(R,COO(6),
     1   COO(8),COO(10),R0OO,RMOO)

      END




      FUNCTION DDISP_f2(R,C6,C8,C10,R0,RM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DAMPC_f2/ADAMP(10),BDAMP(10)
      
      R6=R**6
      R8=R6*R*R
      R10=R8*R*R
      RR=2.0D0*R/(RM+2.5D0*R0)
      DRR=RR/R
      T6=1.0D0-EXP(-ADAMP(6)*RR-BDAMP(6)*RR**2)
      T8=1.0D0-EXP(-ADAMP(8)*RR-BDAMP(8)*RR**2)
      T10=1.0D0-EXP(-ADAMP(10)*RR-BDAMP(10)*RR**2)
      DDISP_f2=6.0D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-ADAMP(6)-2.0D0*
     1      BDAMP(6)*RR)*DRR)+
     2      8.0D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*
     3      (-ADAMP(8)-2.0D0*BDAMP(8)*RR)*DRR)+
     4      10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*
     5      (-ADAMP(10)-2.0D0*BDAMP(10)*RR)*DRR)

      END


C     ***************************************************************
C     DATA FOR OOP SURFACE
C     ***************************************************************
      BLOCK DATA OOP_f2
      
C     ESTE BLOCK DATA CONTIENE LOS DATOS CORRESPONDIENTES A
C     LA SUPERFICIE CALCULADA CON BOO=BOH=1.55
      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DISPCOO_f2/COO(10)
      COMMON/DIATDIOO_f2/R0OO,RMOO
      COMMON/DAMPC_f2/ADAMP(10),BDAMP(10)


      DATA R0OO,RMOO /5.661693D0,2.2818D0/
      DATA COO/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,15.40D0,
     1        0.0D0,235.219943D0, 0.0D0,4066.23929D0/
   
      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,     
     1    3.0951333D0,0.0D0,2.1999000D0,0.0D0,1.6880714D0/
      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,
     1    8.7787895D0,0.0D0,7.2265123D0,0.0D0,5.9487108D0/
      END




