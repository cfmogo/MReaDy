      subroutine potood_f22(rin,Vout,dVout)
      use phys_parameters
      implicit none
      double precision, intent(in)  :: rin
      double precision, intent(out) :: Vout,dVout
      double precision :: VOOD_f22,DVOOD_f22
  
      double precision ::Vout_lr,dVout_lr,f1,df1,fs,dfs
    
      Vout=VOOD_f22(rin)*4.3597482D-18
      dVout=DVOOD_f22(rin)*4.3597482D-18/
     1 5.2917721092d-11

      if (rin.gt.in_bound) then      
         call poto2q_f18(rin,Vout_lr,dVout_lr)      
         fs=f1(rin) 
         dfs=df1(rin)
   
         Vout=Vout*(1-fs)+Vout_lr*fs      
         dVout=dVout*(1-fs)-Vout_lr*dfs      
      end if

      end subroutine potood_f22



        FUNCTION VOOD_f22(R)
        use omp_lib

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!        WRITE(7,*) ''
       ! DATA DE,A1,A2,A3,RE/0.155494053057397D0,6.000D0,10.075D0,
! 1 8.714D0,1.21563D0/

        DE=0.155494053057397D0
        A1=6.000D0
        A2=10.075D0
        A3=8.714D0
        RE=1.21563D0




 !como a1,a2 e a3 estao em angstrons mas o r que recebo da funcao esta em bohr
 !tenho de o converter para angstrons
 !o programa depois devolve o valor de vood em hartree(unidade do DE)
 !posso calcular esta funcao em angstrons porque 
 !devolve vood em hartree(unidade do DE)
        ! print*,'VOOD_f22'
 	ra=r*1.0D10
 	X=Ra-RE

 	VOOD_f22=(-DE*(1.0d0+A1*X+A2*X**2+A3*X**3)*EXP(-A1*X))
      !  print*,'OMP=',OMP_GET_THREAD_NUM(),R,VOOD_f22


 	END

	!****************
	FUNCTION DVOOD_f22(R)
	!****************
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!         WRITE(7,*) ''
	!implicit none
	!real*8 de,a1,a2,a3,re,ra,x,vood_f22,r,dvooddx,dxdra,DVOODDRa
        ! 1 ,dvood_f22

       ! DATA DE,A1,A2,A3,RE/0.155494053057397D0,6.000D0,10.075D0,
    ! 1 8.714D0,1.21563D0/

           DE=0.155494053057397D0
           A1=6.000D0
           A2=10.075D0
           A3=8.714D0
           RE=1.21563D0

!!          write(7,*)'from DVOOD_f22, fake write avoids parallel error'
           
	ra=r*1.0D10
	X=Ra-RE

	!vood_f22=(-DE*(1.0d0+A1*X+A2*X**2+A3*X**3)*EXP(-A1*X))

	!DVOODDX= -DE*(-A1**2*X*EXP(-A1*X) + 2.0d0*A2*X*EXP(-A1*X)- 
    ! 1 A2*A1*X**2*EXP(-A1*X)+
    ! 1     3.0d0*A3*X**2*EXP(-A1*X) - A3*A1*X**3*EXP(-A1*X))
       Dvood_f22= -DE*(-A1**2*X*EXP(-A1*X) + 2.0d0*A2*X*EXP(-A1*X)- 
     1 A2*A1*X**2*EXP(-A1*X)+3.0d0*A3*X**2*
     1 EXP(-A1*X) - A3*A1*X**3*EXP(-A1*X))


	!DXDRa=1.0D0

	!DVOODDRa= DVOODDX*DXDRa

	!Dvood_f22=DVOODDRa
        ! Dvood_f22=DVOODDX

	END
