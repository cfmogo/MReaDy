

!**********************************************************
! PROGRAM TO COMPUTE triplet H2O2 POTENCIAL ENERGY SURFACE
!**********************************************************

subroutine h2o2surt_f31(rin,f,df)
implicit none
real*8 r(6),f,df(6),vh2o2t_f31,rin(6)
integer icalc
!COMMON/SURFAC_f31/ICALC

r=rin/5.2917721092d-11
f= vh2o2t_f31(r(1),r(2),r(3),r(4),r(5),r(6))
icalc=2
!print*,'f=', f

if (icalc==1) return
call dvh2o2t_f31(r(1),r(2),r(3),r(4),r(5),r(6),df(1),df(2),df(3),df(4),df(5),df(6))

f=f*4.3597482D-18
df=df*4.3597482D-18/5.2917721092d-11

return
end


FUNCTION vh2o2t_f31(R1,R2,R3,R4,R5,R6)
!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
implicit none
real*8 v1,VOO_f31,vhhs_f31,vohp_f31,vho23c_f31,v2,vohs_f31,VH2O3c2_f31,v12c,v,vh2o2t_f31,VH2O24c_f31
real*8 vo1d,Vinterm,V3,VH2O3c1_f31,VH2O3c12_f31,VH2OT3C_f31,vhht_f31
real*8 r1,r2,r3,r4,r5,r6
real*8 Vdip_f31
vo1d=7.1955d-2

!CR
!FUNCTION VOO (^3Sigma^-_g) está na superficieho2.f
!FUNCTION VHHs (^1Sigma^+_g) está na superficieh2o.f
!FUNCTION VHHt (^3Sigma) está na superficieh2o.f
!FUNCTION VOHP (^2Pi) está na superficieh2o.f
!FUNCTION VOHs (^2Sigma) está na superficieh2o.f
!FUNCTION VHO23c_f31 (^2A'') está na superficieho2.f
!FUNCTION VH2O3c1 (^1A_1) está na superficieho2.f
!FUNCTION VH2O3c2 (^1A_1) está na superficieho2.f

V1=VOO_f31(R1)+vhhs_f31(R2)+vohs_f31(R3)+vohs_f31(R4)+vohs_f31(R5)+vohs_f31(R6)+VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6) &
      +VH2O3c1_f31(R2,R3,R6)+VH2O3c1_f31(R2,R4,R5)+vo1d

V2=VOO_f31(R1)+vhht_f31(R2)+vohp_f31(R3)+vohp_f31(R4)+vohp_f31(R5)+vohp_f31(R6)+VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6) &
      +VH2O3c2_f31(R2,R3,R6)+VH2O3c2_f31(R2,R4,R5)
 
v12c=VH2O3c12_f31(R2,R3,R6)+VH2O3c12_f31(R2,R4,R5)

Vinterm=((V1+V2)-SQRT((V1-V2)**2+4.0D0*v12c**2))/2.0D0

V3=VOO_f31(R1)+vhhs_f31(R2)+vohp_f31(R3)+vohp_f31(R4)+vohp_f31(R5)+vohp_f31(R6)+VH2OT3C_f31(R2,R3,R6)+VH2OT3C_f31(R2,R4,R5) &
   + VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6)
If (Vinterm < V3) then
   V = Vinterm
else
   V = V3
end if 




vh2o2t_f31=V +VH2O24c_f31(r1,r2,r3,r4,r5,r6)+Vdip_f31(r1,r2,r3,r4,r5,r6)
!	print*,'v',VH2O24c_f31(r1,r2,r3,r4,r5,r6)
!write(30,*) 'v1,v2,v12c,vinterm,v3,v'
!write(30,*)  v1,v2,v12c,vinterm,v3,v


END



!***************************************************
! PROGRAM TO COMPUTE THE DERIVATIVES OF THE H2O2 PES
!***************************************************

subroutine Dvh2o2t_f31(R1,R2,R3,R4,R5,R6,dg1,dg2,dg3,dg4,dg5,dg6)
!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!vo1d=7.1955d-2
implicit none
real*8 v1,VOO_f31,vhhs_f31,vohp_f31,vho23c_f31,VH2O3c1_f31,VH2O3c2_f31,v2,vohs_f31,VH2OT3C_f31,v12c,v,vh2o2,VH2O24c_f31
real*8 Vinterm,V3,VH2O3c12_f31,vhht_f31
real*8 r1,r2,r3,r4,r5,r6,dr1,dr2,dr3,dr4,dr5,dr6,dv1dr1,dv1dr2,dv1dr3,dv1dr4,dv1dr5,dv1dr6
real*8 dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6,dv12dr1,dv12dr2,dv12dr3,dv12dr4,dv12dr5,dv12dr6
real*8 dvdr1,dvdr2,dvdr3,dvdr4,dvdr5,dvdr6,dho2r1a,dho2r1b,dho2r3a,dho2r4b
real*8 dho2r5a,dho2r6b,dr2ta,dr3ta,dr6ta,dr2tb,dr4tb,dr5tb,dvohp_f31,dVOO_f31,dvhhs_f31,dvohs_f31
real*8 dh2o2dr1,dh2o2dr2,dh2o2dr3,dh2o2dr4,dh2o2dr5,dh2o2dr6,dv4cdr1,dv4cdr2,dv4cdr3
real*8 dv4cdr4,dv4cdr5,dv4cdr6,dg1,dg2,dg3,dg4,dg5,dg6,F1,F2,F3,F4,F5,F6
real*8 dv3dr1,dv3dr2,dv3dr3,dv3dr4,dv3dr5,dv3dr6,dvhht_f31,dvintermdv1,dvintermdv2,dvintermdv12
real*8 dr2tb1a,dr3tb1a,dr6tb1a,dr2tb1b,dr4tb1b,dr5tb1b,dr2rl1a,dr3rl1a,dr6rl1a,dr2rl1b,dr4rl1b,dr5rl1b
real*8 dr2tb2a,dr3tb2a,dr6tb2a,dr2tb2b,dr4tb2b,dr5tb2b,dr2rl2a,dr3rl2a,dr6rl2a,dr2rl2b,dr4rl2b,dr5rl2b
real*8 dr2tb12a,dr3tb12a,dr6tb12a,dr2tb12b,dr4tb12b,dr5tb12b
real*8 dvdipdr1,dvdipdr2,dvdipdr3,dvdipdr4,dvdipdr5,dvdipdr6
! JB acrescentou
real*8 VO1D

vo1d=7.1955d-2


!cr
V1=VOO_f31(R1)+vhhs_f31(R2)+vohs_f31(R3)+vohs_f31(R4)+vohs_f31(R5)+vohs_f31(R6)+VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6) &
      +VH2O3c1_f31(R2,R3,R6)+VH2O3c1_f31(R2,R4,R5)+vo1d

V2=VOO_f31(R1)+vhht_f31(R2)+vohp_f31(R3)+vohp_f31(R4)+vohp_f31(R5)+vohp_f31(R6)+VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6) &
      +VH2O3c2_f31(R2,R3,R6)+VH2O3c2_f31(R2,R4,R5)
 
v12c=VH2O3c12_f31(R2,R3,R6)+VH2O3c12_f31(R2,R4,R5)

Vinterm=((V1+V2)-SQRT((V1-V2)**2+4.0D0*v12c**2))/2.0D0

V3=VOO_f31(R1)+vhhs_f31(R2)+vohp_f31(R3)+vohp_f31(R4)+vohp_f31(R5)+vohp_f31(R6)+VH2OT3C_f31(R2,R3,R6)+VH2OT3C_f31(R2,R4,R5) &
   + VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6)

If (Vinterm < V3) then
   V = Vinterm
else
   V = V3
end if

VH2O2=V +VH2O24c_f31(r1,r2,r3,r4,r5,r6)


!derivadas de vinterm em ordem a v1,v2,v12

dvintermdv1= 0.5d0-0.5d0*(v1-v2)/sqrt((v1-v2)**2 + 4.0d0*v12c**2)

dvintermdv2= 0.5d0-0.5d0*(v2-v1)/sqrt((v1-v2)**2 + 4.0d0*v12c**2)

dvintermdv12= -0.5d0*4.0d0*v12c/sqrt((v1-v2)**2 + 4.0d0*v12c**2)


call dervho23c_f31(r1,r3,r5,dr1,dr3,dr5)
dho2r1a=dr1
dho2r3a=dr3
dho2r5a=dr5
call dervho23c_f31(r1,r4,r6,dr1,dr4,dr6)
dho2r1b=dr1
dho2r4b=dr4
dho2r6b=dr6
call derVH2OT3C_f31(r2,r3,r6,dr2,dr3,dr6)
dr2ta=dr2
dr3ta=dr3
dr6ta=dr6
call derVH2OT3C_f31(r2,r4,r5,dr2,dr4,dr5)
dr2tb=dr2
dr4tb=dr4
dr5tb=dr5

!derivadas de v1 em ordem aos rs

!V1=VOO_f31(R1)+vhhs_f31(R2)+vohs_f31(R3)+vohs_f31(R4)+vohs_f31(R5)+vohs_f31(R6)+VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6) &
!      +VH2O3c1_f31(R2,R3,R6)+VH2O3c1_f31(R2,R4,R5)+vo1d+vo1d
!cr
!V1=VOO_f31(R1)+vhhs_f31(R2)+vohs_f31(R3)+vohs_f31(R4)+vohs_f31(R5)+vohs_f31(R6)+VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6) &
!      +VH2O3c1_f31(R2,R3,R6)+VH2O3c1_f31(R2,R4,R5)+vo1d

call dthrebr1_f31(r2,r3,r6,dr2,dr3,dr6) !superficieh2o.f
dr2tb1a=dr2
dr3tb1a=dr3
dr6tb1a=dr6

call drl1_f31(r2,r3,r6,dr2,dr3,dr6) !superficieh2o.f
dr2rl1a=dr2
dr3rl1a=dr3
dr6rl1a=dr6

call dthrebr1_f31(r2,r4,r5,dr2,dr4,dr5) !superficieh2o.f
dr2tb1b=dr2
dr4tb1b=dr4
dr5tb1b=dr5

call drl1_f31(r2,r4,r5,dr2,dr4,dr5) !superficieh2o.f
dr2rl1b=dr2
dr4rl1b=dr4
dr5rl1b=dr5



dv1dr1=dVOO_f31(r1) + dho2r1a + dho2r1b 
dv1dr2=dvhhs_f31(r2) + dr2tb1a + dr2rl1a + dr2tb1b + dr2rl1b
dv1dr3=dvohs_f31(r3)+ dho2r3a + dr3tb1a + dr3rl1a
dv1dr4=dvohs_f31(r4) + dho2r4b + dr4tb1b + dr4rl1b
dv1dr5=dvohs_f31(r5) + dho2r5a + dr5tb1b + dr5rl1b
dv1dr6=dvohs_f31(r6) + dho2r6b + dr6tb1a + dr6rl1a


!print *, 'aqui dv1',dv1dr1,dv1dr2,dv1dr3,dv1dr4,dv1dr5,dv1dr6
!print *, dvhhs_f31(r2) , dr2tb1a ,dr2rl1a , dr2tb1b ,dr2rl1b


!derivadas de v2 em ordem aos rs

!V2=VOO_f31(R1)+vhht_f31(R2)+vohp_f31(R3)+vohp_f31(R4)+vohp_f31(R5)+vohp_f31(R6)+VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6) &
!      +VH2O3c2_f31(R2,R3,R6)+VH2O3c2_f31(R2,R4,R5)
!cr
!V2=VOO_f31(R1)+vhht_f31(R2)+vohp_f31(R3)+vohp_f31(R4)+vohp_f31(R5)+vohp_f31(R6)+VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6) &
!      +VH2O3c2_f31(R2,R3,R6)+VH2O3c2_f31(R2,R4,R5)

call dthrebr2_f31(r2,r3,r6,dr2,dr3,dr6) 
dr2tb2a=dr2
dr3tb2a=dr3
dr6tb2a=dr6

call drl2_f31(r2,r3,r6,dr2,dr3,dr6) 
dr2rl2a=dr2
dr3rl2a=dr3
dr6rl2a=dr6

call dthrebr2_f31(r2,r4,r5,dr2,dr4,dr5) 
dr2tb2b=dr2
dr4tb2b=dr4
dr5tb2b=dr5

call drl2_f31(r2,r4,r5,dr2,dr4,dr5)
dr2rl2b=dr2
dr4rl2b=dr4
dr5rl2b=dr5

dv2dr1=dVOO_f31(r1) + dho2r1a + dho2r1b
dv2dr2=dvhht_f31(r2) + dr2tb2a + dr2rl2a + dr2tb2b + dr2rl2b
dv2dr3=dvohp_f31(r3) + dho2r3a + dr3tb2a + dr3rl2a
dv2dr4=dvohp_f31(r4) + dho2r4b + dr4tb2b + dr4rl2b
dv2dr5=dvohp_f31(r5) + dho2r5a + dr5tb2b + dr5rl2b
dv2dr6=dvohp_f31(r6) + dho2r6b + dr6tb2a + dr6rl2a

  !print *, 'aqui dv2',dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6



!derivadas de v12c em ordem aos rs

!v12c=VH2O3c12_f31(R2,R3,R6)+VH2O3c12_f31(R2,R4,R5)
!cr
!v12c=VH2O3c12_f31(R2,R3,R6)+VH2O3c12_f31(R2,R4,R5)

call dthreb12_f31(r2,r3,r6,dr2,dr3,dr6)
dr2tb12a=dr2
dr3tb12a=dr3
dr6tb12a=dr6

call dthreb12_f31(r2,r4,r5,dr2,dr4,dr5)
dr2tb12b=dr2
dr4tb12b=dr4
dr5tb12b=dr5


dv12dr1=0.0D0
dv12dr2=dr2tb12a + dr2tb12b
dv12dr3=dr3tb12a
dv12dr4=dr4tb12b
dv12dr5=dr5tb12b
dv12dr6=dr6tb12a

!print *, 'aqui dv12',dv12dr1,dv12dr2

 !print *, 'aqui dv12',dv12dr1,dv12dr2,dv12dr3,dv12dr4,dv12dr5,dv12dr6
!derivadas de v3 em ordem aos rs

!V3=VOOD(R1)+vhhs_f31(R2)+vohp_f31(R3)+vohp_f31(R4)+vohp_f31(R5)+vohp_f31(R6)+VH2OT3C_f31(R2,R3,R6)+VH2OT3C_f31(R2,R4,R5) &
!   + VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6)
!cr
!V3=VOO_f31(R1)+vhhs_f31(R2)+vohp_f31(R3)+vohp_f31(R4)+vohp_f31(R5)+vohp_f31(R6)+VH2OT3C_f31(R2,R3,R6)+VH2OT3C_f31(R2,R4,R5) &
!   + VHO23c_f31(R1,R3,R5)+VHO23c_f31(R1,R4,R6)

!cr
!dv3dr1=dvood(r1)+ dho2r1a + dho2r1b

dv3dr1=dVOO_f31(r1)+ dho2r1a + dho2r1b
dv3dr2=dVhhs_f31(r2) + dr2ta + dr2tb
dv3dr3=dvohp_f31(r3) + dr3ta + dho2r3a
dv3dr4=dvohp_f31(r4) + dr4tb + dho2r4b
dv3dr5=dvohp_f31(r5) + dr5tb + dho2r5a
dv3dr6=dvohp_f31(r6) + dr6ta + dho2r6b

!derivadas de v em ordem aos rs quando v=vinterm

if (vinterm.lt.v3) then

dvdr1= dvintermdv1*dv1dr1 + dvintermdv2*dv2dr1 + dvintermdv12*dv12dr1

dvdr2= dvintermdv1*dv1dr2 + dvintermdv2*dv2dr2 + dvintermdv12*dv12dr2 

dvdr3= dvintermdv1*dv1dr3 + dvintermdv2*dv2dr3 + dvintermdv12*dv12dr3 

dvdr4= dvintermdv1*dv1dr4 + dvintermdv2*dv2dr4 + dvintermdv12*dv12dr4 

dvdr5= dvintermdv1*dv1dr5 + dvintermdv2*dv2dr5 + dvintermdv12*dv12dr5 

dvdr6= dvintermdv1*dv1dr6 + dvintermdv2*dv2dr6 + dvintermdv12*dv12dr6 

else

dvdr1= dv3dr1

dvdr2= dv3dr2

dvdr3= dv3dr3

dvdr4= dv3dr4

dvdr5= dv3dr5

dvdr6= dv3dr6

end if


!derivadas de vh2o2

 call dVH2O24c_f31(R1,R2,R3,R4,R5,R6,F1,F2,F3,F4,F5,F6)

dv4cdr1= F1
dv4cdr2= F2
dv4cdr3= F3
dv4cdr4= F4
dv4cdr5= F5
dv4cdr6= F6


!derivadas de Vdip_f31

call dVdip_f31(r1,r2,r3,r4,r5,r6,dr1,dr2,dr3,dr4,dr5,dr6)

dvdipdr1=dr1
dvdipdr2=dr2
dvdipdr3=dr3
dvdipdr4=dr4
dvdipdr5=dr5
dvdipdr6=dr6


!derivada de h2o2

dh2o2dr1= dvdr1 + dv4cdr1 + dvdipdr1
dh2o2dr2= dvdr2 + dv4cdr2 + dvdipdr2
dh2o2dr3= dvdr3 + dv4cdr3 + dvdipdr3
dh2o2dr4= dvdr4 + dv4cdr4 + dvdipdr4
dh2o2dr5= dvdr5 + dv4cdr5 + dvdipdr5
dh2o2dr6= dvdr6 + dv4cdr6 + dvdipdr6


dg1= dh2o2dr1
dg2= dh2o2dr2
dg3= dh2o2dr3
dg4= dh2o2dr4
dg5= dh2o2dr5
dg6= dh2o2dr6

END




!**********************************
FUNCTION VH2O24c_f31(R1,R2,R3,R4,R5,R6)
!**********************************
use param4ct_f31
!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
implicit none
integer i,j,k
integer iref
real*8 vh2o24c_f31,r1,r2,r3,r4,r5,r6
real*8 rho1r(nref),rho22r(nref),rho32r(nref),rho42r(nref),rho234r(nref)

real*8 s1,s2,rho1,rho22,rho32,rho42,rho234
real*8 xc(6),CV(6),expoente


references: Do iref=1,nref


rho1r(iref)=rref(iref,3)+rref(iref,4)+rref(iref,5)+rref(iref,6)
rho22r(iref)=(rref(iref,3)+rref(iref,4)-rref(iref,5)-rref(iref,6))**2
rho32r(iref)=(rref(iref,3)-rref(iref,4)+rref(iref,5)-rref(iref,6))**2
rho42r(iref)=(rref(iref,3)-rref(iref,4)-rref(iref,5)+rref(iref,6))**2
rho234r(iref)=(rref(iref,3)+rref(iref,4)-rref(iref,5)-rref(iref,6))* &
             (rref(iref,3)-rref(iref,4)+rref(iref,5)-rref(iref,6))* &
             (rref(iref,3)-rref(iref,4)-rref(iref,5)+rref(iref,6))


s1=r1-rref(iref,1)
s2=r2-rref(iref,2)
rho1=r3+r4+r5+r6-rho1r(iref)
rho22=(r3+r4-r5-r6)**2-rho22r(iref)
rho32=(r3-r4+r5-r6)**2-rho32r(iref)
rho42=(r3-r4-r5+r6)**2-rho42r(iref)
rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234r(iref)


xc(1)=s1
xc(2)=s2
xc(3)=rho1
xc(4)=rho22
xc(5)=rho32
xc(6)=rho42


DO j=1,ncoord
 CV(j)=0.0d0
     DO k=1,ncoord
        CV(j)=CV(j)+XC(k)*vecp(iref,k,j)
     END DO
END DO



expoente=0.0d0
do J=1,ncoord
    expoente=expoente + CV(j)**2/valp(iref,j)
end do

dec(iref)=exp(-expoente/escal)

pol(iref)=c(1,iref) + c(2,iref)*s1 + c(3,iref)*s2 + c(4,iref)*rho1 + &
    c(5,iref)*s1**2 + c(6,iref)*s1*s2 + c(7,iref)*s1*rho1 + c(8,iref)*s2**2 + &
    c(9,iref)*s2*rho1 + c(10,iref)*rho1**2 + c(11,iref)*rho22 + c(12,iref)*rho32 + &
    c(13,iref)*rho42 + c(14,iref)*s1**3 + c(15,iref)*s1**2*s2 + c(16,iref)*s1**2*rho1 + &
    c(17,iref)*s1*s2**2 + c(18,iref)*s1*s2*rho1 + c(19,iref)*s1*rho1**2 + &
    c(20,iref)*s1*rho22 + c(21,iref)*s1*rho32 + c(22,iref)*s1*rho42 + & 
    c(23,iref)*s2**3 + c(24,iref)*s2**2*rho1 + c(25,iref)*s2*rho1**2 + c(26,iref)*s2*rho22 + &
    c(27,iref)*s2*rho32 + c(28,iref)*s2*rho42 + c(29,iref)*rho1**3 + c(30,iref)*rho1*rho22 + &
    c(31,iref)*rho1*rho32 + c(32,iref)*rho1*rho42 + c(33,iref)*rho234 + &
    c(34,iref)*s1**4 + c(35,iref)*s1**3*s2 + c(36,iref)*s1**3*rho1 + & 
    c(37,iref)*s1**2*s2**2 + c(38,iref)*s1**2*s2*rho1 + &
    c(39,iref)*s1**2*rho1**2 + c(40,iref)*s1**2*rho22 + c(41,iref)*s1**2*rho32 + &
    c(42,iref)*s1**2*rho42 + c(43,iref)*s1*s2**3 + &
    c(44,iref)*s1*s2**2*rho1 + c(45,iref)*s1*s2*rho1**2 + c(46,iref)*s1*s2*rho22 + &
    c(47,iref)*s1*s2*rho32 + c(48,iref)*s1*s2*rho42 + &
    c(49,iref)*s1*rho1**3 + c(50,iref)*s1*rho1*rho22 + c(51,iref)*s1*rho1*rho32 + &
    c(52,iref)*s1*rho1*rho42 + c(53,iref)*s1*rho234 + &
    c(54,iref)*s2**4 + c(55,iref)*s2**3*rho1 + c(56,iref)*s2**2*rho1**2 + &
    c(57,iref)*s2**2*rho22 + c(58,iref)*s2**2*rho32 + &
    c(59,iref)*s2**2*rho42 + c(60,iref)*s2*rho1**3 + c(61,iref)*s2*rho1*rho22 + &
    c(62,iref)*s2*rho1*rho32 + c(63,iref)*s2*rho1*rho42 + &
    c(64,iref)*s2*rho234 + c(65,iref)*rho1**4 + c(66,iref)*rho1**2*rho22 + &
    c(67,iref)*rho1**2*rho32 + c(68,iref)*rho1**2*rho42 + &
    c(69,iref)*rho1*rho234 + c(70,iref)*rho22*rho32 + c(71,iref)*rho22*rho42 + &
    c(72,iref)*rho32*rho42 + c(73,iref)*rho22*rho22 + &
    c(74,iref)*rho32*rho32 + c(75,iref)*rho42*rho42 + &
    c(76,iref)*s1**5 + c(77,iref)*s1**4*s2 + c(78,iref)*s1**4*rho1 + &
    c(79,iref)*s1**3*s2**2 + c(80,iref)*s1**3*s2*rho1 + &
    c(81,iref)*s1**3*rho1**2 + c(82,iref)*s1**3*rho22 + c(83,iref)*s1**3*rho32 + &
    c(84,iref)*s1**3*rho42 + c(85,iref)*s1**2*s2**3 + &
    c(86,iref)*s1**2*s2**2*rho1 + c(87,iref)*s1**2*s2*rho1**2 + c(88,iref)*s1**2*s2*rho22 + &
    c(89,iref)*s1**2*s2*rho32 + c(90,iref)*s1**2*s2*rho42 + &
    c(91,iref)*s1**2*rho1**3 + c(92,iref)*s1**2*rho1*rho22 + c(93,iref)*s1**2*rho1*rho32 + &
    c(94,iref)*s1**2*rho1*rho42 + c(95,iref)*s1**2*rho234 + &
    c(96,iref)*s1*s2**4 + c(97,iref)*s1*s2**3*rho1 + c(98,iref)*s1*s2**2*rho1**2 + &
    c(99,iref)*s1*s2**2*rho22 + c(100,iref)*s1*s2**2*rho32 + &
    c(101,iref)*s1*s2**2*rho42 + c(102,iref)*s1*s2*rho1**3 + c(103,iref)*s1*s2*rho1*rho22 + &
    c(104,iref)*s1*s2*rho1*rho32 + c(105,iref)*s1*s2*rho1*rho42 + &
    c(106,iref)*s1*s2*rho234 + c(107,iref)*s1*rho1**4 + c(108,iref)*s1*rho1**2*rho22 + &
    c(109,iref)*s1*rho1**2*rho32 + c(110,iref)*s1*rho1**2*rho42 + &
    c(111,iref)*s1*rho1*rho234 + c(112,iref)*s1*rho22*rho32 + c(113,iref)*s1*rho22*rho42 + &
    c(114,iref)*s1*rho32*rho42 + c(115,iref)*s1*rho22*rho22 + &
    c(116,iref)*s1*rho32*rho32 + c(117,iref)*s1*rho42*rho42 + c(118,iref)*s2**5 + &
    c(119,iref)*s2**4*rho1 + c(120,iref)*s2**3*rho1**2 + c(121,iref)*s2**3*rho22 + &
    c(122,iref)*s2**3*rho32 + c(123,iref)*s2**3*rho42 + c(124,iref)*s2**2*rho1**3 + &
    c(125,iref)*s2**2*rho1*rho22 + c(126,iref)*s2**2*rho1*rho32 + &
    c(127,iref)*s2**2*rho1*rho42 + c(128,iref)*s2**2*rho234 + c(129,iref)*s2*rho1**4 + &
    c(130,iref)*s2*rho1**2*rho22 + c(131,iref)*s2*rho1**2*rho32 + &
    c(132,iref)*s2*rho1**2*rho42 + c(133,iref)*s2*rho1*rho234 + c(134,iref)*s2*rho22*rho32 + &
    c(135,iref)*s2*rho22*rho42 + c(136,iref)*s2*rho32*rho42 + &
    c(137,iref)*s2*rho22*rho22 + c(138,iref)*s2*rho32*rho32 + c(139,iref)*s2*rho42*rho42 + &
    c(140,iref)*rho1**5 + c(141,iref)*rho1**3*rho22 + c(142,iref)*rho1**3*rho32 + &
    c(143,iref)*rho1**3*rho42 + c(144,iref)*rho1**2*rho234 + &
    c(145,iref)*rho1*rho22**2 + c(146,iref)*rho1*rho32**2 + c(147,iref)*rho1*rho42**2 + &
    c(148,iref)*rho1*rho22*rho32 + c(149,iref)*rho1*rho22*rho42 + &
    c(150,iref)*rho1*rho32*rho42 + c(151,iref)*rho22*rho234 + &
    c(152,iref)*rho32*rho234 + c(153,iref)*rho42*rho234

V4c(iref)=pol(iref)*dec(iref)

 


end do references

VH2O24c_f31=0.0d0
do i=1,nref
      VH2O24c_f31=VH2O24c_f31+V4c(i)
end do 

END


!***********************************************************
subroutine DVH2O24c_f31(R1,R2,R3,R4,R5,R6,F1,F2,F3,F4,F5,F6)
!***********************************************************

use param4ct_f31

!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
implicit none
integer j,k
integer iref

real*8 r1,r2,r3,r4,r5,r6
real*8 F1,F2,F3,F4,F5,F6

real*8 s1,s2,rho1,rho22,rho32,rho42,rho234
real*8 rho1r(nref),rho22r(nref),rho32r(nref),rho42r(nref),rho234r(nref)

real*8 xc(6),CV(6),expoente

real*8 dpolds1,dpolds2,dpoldrho1,dpoldrho22
real*8 dpoldrho32,dpoldrho42,dpoldrho234
real*8 ddecds1,ddecds2,ddecdrho1,ddecdrho22
real*8 ddecdrho32,ddecdrho42

real*8 ds1dr1,ds2dr2,drho1dr3,drho22dr3,drho32dr3,drho42dr3,drho234dr3,drho1dr4,drho22dr4,drho32dr4
real*8 drho42dr4,drho234dr4,drho1dr5,drho22dr5,drho32dr5,drho42dr5,drho234dr5,drho1dr6,drho22dr6,drho32dr6
real*8 drho42dr6,drho234dr6

real*8 dv4ds1,dv4ds2,dv4drho1,dv4drho22,dv4drho32,dv4drho42,dv4drho234,dv4dr1(nref),dv4dr2(nref)
real*8 dv4dr3(nref),dv4dr4(nref),dv4dr5(nref),dv4dr6(nref)


! derivadas das coordenadas

ds1dr1=1.0d0
ds2dr2=1.0d0

drho1dr3=1.0d0
drho22dr3=2.0d0*(r3+r4-r5-r6)
drho32dr3=2.0d0*(r3-r4+r5-r6)
drho42dr3=2.0d0*(r3-r4-r5+r6)
drho234dr3=3.0d0*r3**2 -2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r3*r5 + 2.0d0*r4*r6 + 2.0d0*r5*r6 -&
           r4**2 + 2.0d0*r4*r5 - r5**2 - r6**2

drho1dr4=1.0d0
drho22dr4=2.0d0*(r3+r4-r5-r6)
drho32dr4=-2.0d0*(r3-r4+r5-r6)
drho42dr4=-2.0d0*(r3-r4-r5+r6)
drho234dr4=3.0d0*r4**2 - r3**2 + 2.0d0*r3*r6 - 2.0d0*r3*r4 - 2.0d0*r4*r5 - 2.0d0*r4*r6 + &
           2.0d0*r3*r5 + 2.0d0*r5*r6 - r5**2 - r6**2

drho1dr5=1.0d0
drho22dr5=-2.0d0*(r3+r4-r5-r6)
drho32dr5=2.0d0*(r3-r4+r5-r6)
drho42dr5=-2.0d0*(r3-r4-r5+r6)
drho234dr5=3.0d0*r5**2 - r3**2 + 2.0d0*r3*r6 - r4**2 + 2.0d0*r3*r4 -2.0d0*r4*r5 + 2.0d0*r4*r6 -&
           2.0d0*r5*r3 - 2.0d0*r5*r6 - r6**2


drho1dr6=1.0d0
drho22dr6=-2.0d0*(r3+r4-r5-r6)
drho32dr6=-2.0d0*(r3-r4+r5-r6)
drho42dr6=2.0d0*(r3-r4-r5+r6)
drho234dr6=3.0d0*r6**2 - r3**2 + 2.0d0*r3*r5 + 2.0d0*r3*r4 - r4**2 + 2.0d0*r4*r5 - r5**2 -&
           2.0d0*r6*r3 -2.0d0*r6*r4 - 2.0d0*r6*r5


! para as varias referencias


references: Do iref=1,nref


rho1r(iref)=rref(iref,3)+rref(iref,4)+rref(iref,5)+rref(iref,6)
rho22r(iref)=(rref(iref,3)+rref(iref,4)-rref(iref,5)-rref(iref,6))**2
rho32r(iref)=(rref(iref,3)-rref(iref,4)+rref(iref,5)-rref(iref,6))**2
rho42r(iref)=(rref(iref,3)-rref(iref,4)-rref(iref,5)+rref(iref,6))**2
rho234r(iref)=(rref(iref,3)+rref(iref,4)-rref(iref,5)-rref(iref,6))* &
             (rref(iref,3)-rref(iref,4)+rref(iref,5)-rref(iref,6))* &
             (rref(iref,3)-rref(iref,4)-rref(iref,5)+rref(iref,6))


s1=r1-rref(iref,1)
s2=r2-rref(iref,2)
rho1=r3+r4+r5+r6-rho1r(iref)
rho22=(r3+r4-r5-r6)**2-rho22r(iref)
rho32=(r3-r4+r5-r6)**2-rho32r(iref)
rho42=(r3-r4-r5+r6)**2-rho42r(iref)
rho234=(r3+r4-r5-r6)*(r3-r4+r5-r6)*(r3-r4-r5+r6)-rho234r(iref)


xc(1)=s1
xc(2)=s2
xc(3)=rho1
xc(4)=rho22
xc(5)=rho32
xc(6)=rho42

DO j=1,ncoord
 CV(j)=0.0d0
     DO k=1,ncoord
        CV(j)=CV(j)+XC(k)*vecp(iref,k,j)
     END DO
END DO

expoente=0.0d0
do J=1,ncoord
    expoente=expoente + CV(j)**2/valp(iref,j)
end do

dec(iref)=exp(-expoente/escal)


pol(iref)=c(1,iref) + c(2,iref)*s1 + c(3,iref)*s2 + c(4,iref)*rho1 + &
    c(5,iref)*s1**2 + c(6,iref)*s1*s2 + c(7,iref)*s1*rho1 + c(8,iref)*s2**2 + &
    c(9,iref)*s2*rho1 + c(10,iref)*rho1**2 + c(11,iref)*rho22 + c(12,iref)*rho32 + &
    c(13,iref)*rho42 + c(14,iref)*s1**3 + c(15,iref)*s1**2*s2 + c(16,iref)*s1**2*rho1 + &
    c(17,iref)*s1*s2**2 + c(18,iref)*s1*s2*rho1 + c(19,iref)*s1*rho1**2 + &
    c(20,iref)*s1*rho22 + c(21,iref)*s1*rho32 + c(22,iref)*s1*rho42 + & 
    c(23,iref)*s2**3 + c(24,iref)*s2**2*rho1 + c(25,iref)*s2*rho1**2 + c(26,iref)*s2*rho22 + &
    c(27,iref)*s2*rho32 + c(28,iref)*s2*rho42 + c(29,iref)*rho1**3 + c(30,iref)*rho1*rho22 + &
    c(31,iref)*rho1*rho32 + c(32,iref)*rho1*rho42 + c(33,iref)*rho234 + &
    c(34,iref)*s1**4 + c(35,iref)*s1**3*s2 + c(36,iref)*s1**3*rho1 + & 
    c(37,iref)*s1**2*s2**2 + c(38,iref)*s1**2*s2*rho1 + &
    c(39,iref)*s1**2*rho1**2 + c(40,iref)*s1**2*rho22 + c(41,iref)*s1**2*rho32 + &
    c(42,iref)*s1**2*rho42 + c(43,iref)*s1*s2**3 + &
    c(44,iref)*s1*s2**2*rho1 + c(45,iref)*s1*s2*rho1**2 + c(46,iref)*s1*s2*rho22 + &
    c(47,iref)*s1*s2*rho32 + c(48,iref)*s1*s2*rho42 + &
    c(49,iref)*s1*rho1**3 + c(50,iref)*s1*rho1*rho22 + c(51,iref)*s1*rho1*rho32 + &
    c(52,iref)*s1*rho1*rho42 + c(53,iref)*s1*rho234 + &
    c(54,iref)*s2**4 + c(55,iref)*s2**3*rho1 + c(56,iref)*s2**2*rho1**2 + &
    c(57,iref)*s2**2*rho22 + c(58,iref)*s2**2*rho32 + &
    c(59,iref)*s2**2*rho42 + c(60,iref)*s2*rho1**3 + c(61,iref)*s2*rho1*rho22 + &
    c(62,iref)*s2*rho1*rho32 + c(63,iref)*s2*rho1*rho42 + &
    c(64,iref)*s2*rho234 + c(65,iref)*rho1**4 + c(66,iref)*rho1**2*rho22 + &
    c(67,iref)*rho1**2*rho32 + c(68,iref)*rho1**2*rho42 + &
    c(69,iref)*rho1*rho234 + c(70,iref)*rho22*rho32 + c(71,iref)*rho22*rho42 + &
    c(72,iref)*rho32*rho42 + c(73,iref)*rho22*rho22 + &
    c(74,iref)*rho32*rho32 + c(75,iref)*rho42*rho42 + &
    c(76,iref)*s1**5 + c(77,iref)*s1**4*s2 + c(78,iref)*s1**4*rho1 + &
    c(79,iref)*s1**3*s2**2 + c(80,iref)*s1**3*s2*rho1 + &
    c(81,iref)*s1**3*rho1**2 + c(82,iref)*s1**3*rho22 + c(83,iref)*s1**3*rho32 + &
    c(84,iref)*s1**3*rho42 + c(85,iref)*s1**2*s2**3 + &
    c(86,iref)*s1**2*s2**2*rho1 + c(87,iref)*s1**2*s2*rho1**2 + c(88,iref)*s1**2*s2*rho22 + &
    c(89,iref)*s1**2*s2*rho32 + c(90,iref)*s1**2*s2*rho42 + &
    c(91,iref)*s1**2*rho1**3 + c(92,iref)*s1**2*rho1*rho22 + c(93,iref)*s1**2*rho1*rho32 + &
    c(94,iref)*s1**2*rho1*rho42 + c(95,iref)*s1**2*rho234 + &
    c(96,iref)*s1*s2**4 + c(97,iref)*s1*s2**3*rho1 + c(98,iref)*s1*s2**2*rho1**2 + &
    c(99,iref)*s1*s2**2*rho22 + c(100,iref)*s1*s2**2*rho32 + &
    c(101,iref)*s1*s2**2*rho42 + c(102,iref)*s1*s2*rho1**3 + c(103,iref)*s1*s2*rho1*rho22 + &
    c(104,iref)*s1*s2*rho1*rho32 + c(105,iref)*s1*s2*rho1*rho42 + &
    c(106,iref)*s1*s2*rho234 + c(107,iref)*s1*rho1**4 + c(108,iref)*s1*rho1**2*rho22 + &
    c(109,iref)*s1*rho1**2*rho32 + c(110,iref)*s1*rho1**2*rho42 + &
    c(111,iref)*s1*rho1*rho234 + c(112,iref)*s1*rho22*rho32 + c(113,iref)*s1*rho22*rho42 + &
    c(114,iref)*s1*rho32*rho42 + c(115,iref)*s1*rho22*rho22 + &
    c(116,iref)*s1*rho32*rho32 + c(117,iref)*s1*rho42*rho42 + c(118,iref)*s2**5 + &
    c(119,iref)*s2**4*rho1 + c(120,iref)*s2**3*rho1**2 + c(121,iref)*s2**3*rho22 + &
    c(122,iref)*s2**3*rho32 + c(123,iref)*s2**3*rho42 + c(124,iref)*s2**2*rho1**3 + &
    c(125,iref)*s2**2*rho1*rho22 + c(126,iref)*s2**2*rho1*rho32 + &
    c(127,iref)*s2**2*rho1*rho42 + c(128,iref)*s2**2*rho234 + c(129,iref)*s2*rho1**4 + &
    c(130,iref)*s2*rho1**2*rho22 + c(131,iref)*s2*rho1**2*rho32 + &
    c(132,iref)*s2*rho1**2*rho42 + c(133,iref)*s2*rho1*rho234 + c(134,iref)*s2*rho22*rho32 + &
    c(135,iref)*s2*rho22*rho42 + c(136,iref)*s2*rho32*rho42 + &
    c(137,iref)*s2*rho22*rho22 + c(138,iref)*s2*rho32*rho32 + c(139,iref)*s2*rho42*rho42 + &
    c(140,iref)*rho1**5 + c(141,iref)*rho1**3*rho22 + c(142,iref)*rho1**3*rho32 + &
    c(143,iref)*rho1**3*rho42 + c(144,iref)*rho1**2*rho234 + &
    c(145,iref)*rho1*rho22**2 + c(146,iref)*rho1*rho32**2 + c(147,iref)*rho1*rho42**2 + &
    c(148,iref)*rho1*rho22*rho32 + c(149,iref)*rho1*rho22*rho42 + &
    c(150,iref)*rho1*rho32*rho42 + c(151,iref)*rho22*rho234 + &
    c(152,iref)*rho32*rho234 + c(153,iref)*rho42*rho234

V4c(iref)=pol(iref)*dec(iref)

!derivada em ordem a s1

dpolds1=c(2,iref) + 2.0d0*c(5,iref)*s1 + c(6,iref)*s2 + c(7,iref)*rho1 + & 
         3.0d0*c(14,iref)*s1**2 + 2.0d0*c(15,iref)*s1*s2 + 2.0d0*c(16,iref)*s1*rho1 + &
         c(17,iref)*s2**2 + c(18,iref)*s2*rho1 + c(19,iref)*rho1**2 + &
         c(20,iref)*rho22 + c(21,iref)*rho32 + c(22,iref)*rho42 + &
         4.0d0*c(34,iref)*s1**3 + 3.0d0*c(35,iref)*s1**2*s2 + &
         3.0d0*c(36,iref)*s1**2*rho1 + 2.0d0*c(37,iref)*s1*s2**2 + &
         2.0d0*c(38,iref)*s1*s2*rho1 + 2.0d0*c(39,iref)*s1*rho1**2 + &
         2.0d0*c(40,iref)*s1*rho22 + 2.0d0*c(41,iref)*s1*rho32 + &
         2.0d0*c(42,iref)*s1*rho42 + c(43,iref)*s2**3 + c(44,iref)*s2**2*rho1 + &
         c(45,iref)*s2*rho1**2 + c(46,iref)*s2*rho22 + &
         c(47,iref)*s2*rho32 + c(48,iref)*s2*rho42 + c(49,iref)*rho1**3 + &
         c(50,iref)*rho1*rho22 + c(51,iref)*rho1*rho32 + &
         c(52,iref)*rho1*rho42 + c(53,iref)*rho234 + &
         5.0d0*c(76,iref)*s1**4 + 4.0d0*c(77,iref)*s1**3*s2 + &
         4.0d0*c(78,iref)*s1**3*rho1 + 3.0d0*c(79,iref)*s1**2*s2**2 + &
         3.0d0*c(80,iref)*s1**2*s2*rho1 + 3.0d0*c(81,iref)*s1**2*rho1**2 + &
         3.0d0*c(82,iref)*s1**2*rho22 + &
         3.0d0*c(83,iref)*s1**2*rho32 + 3.0d0*c(84,iref)*s1**2*rho42 + &
         2.0d0*c(85,iref)*s1*s2**3 + &
         2.0d0*c(86,iref)*s1*s2**2*rho1 + 2.0d0*c(87,iref)*s1*s2*rho1**2 + &
         2.0d0*c(88,iref)*s1*s2*rho22 + &
         2.0d0*c(89,iref)*s1*s2*rho32 + 2.0d0*c(90,iref)*s1*s2*rho42 + &
         2.0d0*c(91,iref)*s1*rho1**3 + 2.0d0*c(92,iref)*s1*rho1*rho22 + &
         2.0d0*c(93,iref)*s1*rho1*rho32 + &
         2.0d0*c(94,iref)*s1*rho1*rho42 + 2.0d0*c(95,iref)*s1*rho234 + &
         c(96,iref)*s2**4 + c(97,iref)*s2**3*rho1 + c(98,iref)*s2**2*rho1**2 + &
         c(99,iref)*s2**2*rho22 + c(100,iref)*s2**2*rho32 + &
         c(101,iref)*s2**2*rho42 + c(102,iref)*s2*rho1**3 + &
         c(103,iref)*s2*rho1*rho22 + c(104,iref)*s2*rho1*rho32 + &
         c(105,iref)*s2*rho1*rho42 + c(106,iref)*s2*rho234 + &
         c(107,iref)*rho1**4 + c(108,iref)*rho1**2*rho22 + &
         c(109,iref)*rho1**2*rho32 + c(110,iref)*rho1**2*rho42 + &
         c(111,iref)*rho1*rho234 + c(112,iref)*rho22*rho32 + c(113,iref)*rho22*rho42 + &
         c(114,iref)*rho32*rho42 + c(115,iref)*rho22*rho22 + &
         c(116,iref)*rho32*rho32 + c(117,iref)*rho42*rho42 
  
 
!derivada em ordem a s2    
  
dpolds2=c(3,iref) + c(6,iref)*s1 + 2.0d0*c(8,iref)*s2 + c(9,iref)*rho1 + &
         c(15,iref)*s1**2 + 2.0d0*c(17,iref)*s1*s2 + c(18,iref)*s1*rho1 + &
         3.0d0*c(23,iref)*s2**2 + 2.0d0*c(24,iref)*s2*rho1 + c(25,iref)*rho1**2 + &
         c(26,iref)*rho22 + c(27,iref)*rho32 + c(28,iref)*rho42+ &
         c(35,iref)*s1**3 + 2.0d0*c(37,iref)*s1**2*s2 + c(38,iref)*s1**2*rho1 + &
         3.0d0*c(43,iref)*s1*s2**2 + &
         2.0d0*c(44,iref)*s1*s2*rho1 + c(45,iref)*s1*rho1**2 + c(46,iref)*s1*rho22 + &
         c(47,iref)*s1*rho32 + c(48,iref)*s1*rho42 + &
         4.0d0*c(54,iref)*s2**3 + 3.0d0*c(55,iref)*s2**2*rho1 + &
         2.0d0*c(56,iref)*s2*rho1**2 + 2.0d0*c(57,iref)*s2*rho22 + &
         2.0d0*c(58,iref)*s2*rho32 + 2.0d0*c(59,iref)*s2*rho42 + &
         c(60,iref)*rho1**3 + c(61,iref)*rho1*rho22 + &
         c(62,iref)*rho1*rho32 + c(63,iref)*rho1*rho42 + c(64,iref)*rho234 + &
         c(77,iref)*s1**4 + 2.0d0*c(79,iref)*s1**3*s2 + c(80,iref)*s1**3*rho1 + &
         3.0d0*c(85,iref)*s1**2*s2**2 + &
         2.0d0*c(86,iref)*s1**2*s2*rho1 + c(87,iref)*s1**2*rho1**2 + &
         c(88,iref)*s1**2*rho22 + c(89,iref)*s1**2*rho32 + &
         c(90,iref)*s1**2*rho42 + 4.0d0*c(96,iref)*s1*s2**3 + &
         3.0d0*c(97,iref)*s1*s2**2*rho1 + 2.0d0*c(98,iref)*s1*s2*rho1**2 + &
         2.0d0*c(99,iref)*s1*s2*rho22 + 2.0d0*c(100,iref)*s1*s2*rho32 + &
         2.0d0*c(101,iref)*s1*s2*rho42 + c(102,iref)*s1*rho1**3 + &
         c(103,iref)*s1*rho1*rho22 + c(104,iref)*s1*rho1*rho32 + &
         c(105,iref)*s1*rho1*rho42 + c(106,iref)*s1*rho234 + 5.0d0*c(118,iref)*s2**4 + &
         4.0d0*c(119,iref)*s2**3*rho1 + &
         3.0d0*c(120,iref)*s2**2*rho1**2 + 3.0d0*c(121,iref)*s2**2*rho22 + &
         3.0d0*c(122,iref)*s2**2*rho32 + 3.0d0*c(123,iref)*s2**2*rho42 + &
         2.0d0*c(124,iref)*s2*rho1**3 + &
         2.0d0*c(125,iref)*s2*rho1*rho22 + 2.0d0*c(126,iref)*s2*rho1*rho32 + &
         2.0d0*c(127,iref)*s2*rho1*rho42 + 2.0d0*c(128,iref)*s2*rho234 + &
         c(129,iref)*rho1**4 + c(130,iref)*rho1**2*rho22 + c(131,iref)*rho1**2*rho32 + &
         c(132,iref)*rho1**2*rho42 + c(133,iref)*rho1*rho234 + c(134,iref)*rho22*rho32 + &
         c(135,iref)*rho22*rho42 + c(136,iref)*rho32*rho42 + &
         c(137,iref)*rho22*rho22 + c(138,iref)*rho32*rho32 + c(139,iref)*rho42*rho42
  
 
!derivada em ordem a rho1
   
dpoldrho1=c(4,iref) + c(7,iref)*s1 + c(9,iref)*s2 + 2.0d0*c(10,iref)*rho1 + &
         c(16,iref)*s1**2 + c(18,iref)*s1*s2 + 2.0d0*c(19,iref)*s1*rho1 + &
         c(24,iref)*s2**2 + 2.0d0*c(25,iref)*s2*rho1 + 3.0d0*c(29,iref)*rho1**2 + &
         c(30,iref)*rho22 + c(31,iref)*rho32 + c(32,iref)*rho42+ &
         c(36,iref)*s1**3 + c(38,iref)*s1**2*s2 + 2.0d0*c(39,iref)*s1**2*rho1 + &
         c(44,iref)*s1*s2**2 + &
         2.0d0*c(45,iref)*s1*s2*rho1 + 3.0d0*c(49,iref)*s1*rho1**2 + c(50,iref)*s1*rho22 + &
         c(51,iref)*s1*rho32 + c(52,iref)*s1*rho42 + &
         c(55,iref)*s2**3 + 2.0d0*c(56,iref)*s2**2*rho1 + 3.0d0*c(60,iref)*s2*rho1**2 + &
         c(61,iref)*s2*rho22 + c(62,iref)*s2*rho32 + &
         c(63,iref)*s2*rho42 + 4.0d0*c(65,iref)*rho1**3 + 2.0d0*c(66,iref)*rho1*rho22 + &
         2.0d0*c(67,iref)*rho1*rho32 + &
         2.0d0*c(68,iref)*rho1*rho42 + c(69,iref)*rho234 + &
         c(78,iref)*s1**4 + c(80,iref)*s1**3*s2 + 2.0d0*c(81,iref)*s1**3*rho1 + &
         c(86,iref)*s1**2*s2**2 + 2.0d0*c(87,iref)*s1**2*s2*rho1 + &
         3.0d0*c(91,iref)*s1**2*rho1**2 + c(92,iref)*s1**2*rho22 + c(93,iref)*s1**2*rho32 + &
         c(94,iref)*s1**2*rho42 + c(97,iref)*s1*s2**3 + &
         2.0d0*c(98,iref)*s1*s2**2*rho1 + 3.0d0*c(102,iref)*s1*s2*rho1**2 + &
         c(103,iref)*s1*s2*rho22 + c(104,iref)*s1*s2*rho32 + c(105,iref)*s1*s2*rho42 + &
         4.0d0*c(107,iref)*s1*rho1**3 + 2.0d0*c(108,iref)*s1*rho1*rho22 + &
         2.0d0*c(109,iref)*s1*rho1*rho32 + 2.0d0*c(110,iref)*s1*rho1*rho42 + &
         c(111,iref)*s1*rho234 + c(119,iref)*s2**4 + 2.0d0*c(120,iref)*s2**3*rho1 + &
         3.0d0*c(124,iref)*s2**2*rho1**2 + c(125,iref)*s2**2*rho22 + &
         c(126,iref)*s2**2*rho32 + c(127,iref)*s2**2*rho42 + 4.0d0*c(129,iref)*s2*rho1**3 + &
         2.0d0*c(130,iref)*s2*rho1*rho22 + &
         2.0d0*c(131,iref)*s2*rho1*rho32 + 2.0d0*c(132,iref)*s2*rho1*rho42 + &
         c(133,iref)*s2*rho234 + 5.0d0*c(140,iref)*rho1**4 + &
         3.0d0*c(141,iref)*rho1**2*rho22 + 3.0d0*c(142,iref)*rho1**2*rho32 + &
         3.0d0*c(143,iref)*rho1**2*rho42 + 2.0d0*c(144,iref)*rho1*rho234 + &
         c(145,iref)*rho22**2 + c(146,iref)*rho32**2 + &
         c(147,iref)*rho42**2 + c(148,iref)*rho22*rho32 + &
         c(149,iref)*rho22*rho42 + c(150,iref)*rho32*rho42 

!derivada  em ordem a rho22


dpoldrho22=c(11,iref) + c(20,iref)*s1 + c(26,iref)*s2 + c(30,iref)*rho1+ &
         c(40,iref)*s1**2 + c(46,iref)*s1*s2 + c(50,iref)*s1*rho1 + & 
         c(57,iref)*s2**2 + c(61,iref)*s2*rho1 + c(66,iref)*rho1**2 + &
         c(70,iref)*rho32 + c(71,iref)*rho42 + 2.0d0*c(73,iref)*rho22 + &
         c(82,iref)*s1**3 + c(88,iref)*s1**2*s2 + c(92,iref)*s1**2*rho1 + &
         c(99,iref)*s1*s2**2 + c(103,iref)*s1*s2*rho1 + & 
         c(108,iref)*s1*rho1**2 + c(112,iref)*s1*rho32 + c(113,iref)*s1*rho42 + &
         2.0d0*c(115,iref)*s1*rho22 + c(121,iref)*s2**3 + &
         c(125,iref)*s2**2*rho1 + c(130,iref)*s2*rho1**2 + c(134,iref)*s2*rho32 + &
         c(135,iref)*s2*rho42 + 2.0d0*c(137,iref)*s2*rho22 + & 
         c(141,iref)*rho1**3 + 2.0d0*c(145,iref)*rho1*rho22  + &
         c(148,iref)*rho1*rho32 + c(149,iref)*rho1*rho42 + &
         c(151,iref)*rho234 

  
!derivada em ordem a rho32

dpoldrho32=c(12,iref) + c(21,iref)*s1 + c(27,iref)*s2 + c(31,iref)*rho1 + &
         c(41,iref)*s1**2 + c(47,iref)*s1*s2 + c(51,iref)*s1*rho1 + &
         c(58,iref)*s2**2 +c(62,iref)*s2*rho1 + c(67,iref)*rho1**2 + &
         c(70,iref)*rho22 + c(72,iref)*rho42 + 2.0d0*c(74,iref)*rho32 + &
         c(83,iref)*s1**3 + c(89,iref)*s1**2*s2 + c(93,iref)*s1**2*rho1 + &
         c(100,iref)*s1*s2**2 + c(104,iref)*s1*s2*rho1 + &
         c(109,iref)*s1*rho1**2 + c(112,iref)*s1*rho22 + c(114,iref)*s1*rho42 + &
         2.0d0*c(116,iref)*s1*rho32 + c(122,iref)*s2**3 + &
         c(126,iref)*s2**2*rho1 + c(131,iref)*s2*rho1**2 + c(134,iref)*s2*rho22 + &
         c(136,iref)*s2*rho42 + 2.0d0*c(138,iref)*s2*rho32 + &
         c(142,iref)*rho1**3 + &
         2.0d0*c(146,iref)*rho1*rho32 + c(148,iref)*rho1*rho22 + &
         c(150,iref)*rho1*rho42 + c(152,iref)*rho234 

!derivada em ordem a rho42

dpoldrho42=c(13,iref) + c(22,iref)*s1 + c(28,iref)*s2 + c(32,iref)*rho1 + &
         c(42,iref)*s1**2 + c(48,iref)*s1*s2 + c(52,iref)*s1*rho1 + &
         c(59,iref)*s2**2 + c(63,iref)*s2*rho1 + c(68,iref)*rho1**2 + &
         c(71,iref)*rho22 + c(72,iref)*rho32 + 2.0d0*c(75,iref)*rho42+ &
         c(84,iref)*s1**3 + c(90,iref)*s1**2*s2 + c(94,iref)*s1**2*rho1 + &
         c(101,iref)*s1*s2**2 + c(105,iref)*s1*s2*rho1 + &
         c(110,iref)*s1*rho1**2 + c(113,iref)*s1*rho22+ c(114,iref)*s1*rho32 + &
         2.0d0*c(117,iref)*s1*rho42 + c(123,iref)*s2**3 + &
         c(127,iref)*s2**2*rho1 + c(132,iref)*s2*rho1**2 + &
         c(135,iref)*s2*rho22 + c(136,iref)*s2*rho32 + &
         2.0d0*c(139,iref)*s2*rho42 + c(143,iref)*rho1**3 + &
         2.0d0*c(147,iref)*rho1*rho42 + c(149,iref)*rho1*rho22 + &
         c(150,iref)*rho1*rho32 +  c(153,iref)*rho234


!derivada em ordem a rho234

dpoldrho234=c(33,iref) + c(53,iref)*s1 + c(64,iref)*s2 + c(69,iref)*rho1 + &
         c(95,iref)*s1**2 + &
         c(106,iref)*s1*s2 + c(111,iref)*s1*rho1 + c(128,iref)*s2**2 + &
         c(133,iref)*s2*rho1 + &
         c(144,iref)*rho1**2 + c(151,iref)*rho22 + c(152,iref)*rho32 + &
         c(153,iref)*rho42


ddecds1=0.0d0
ddecds2=0.0d0
ddecdrho1=0.0d0
ddecdrho22=0.0d0
ddecdrho32=0.0d0
ddecdrho42=0.0d0

do j=1,ncoord

   ddecds1=ddecds1-(2.0d0/escal)*CV(j)*vecp(iref,1,j)/valp(iref,j)*exp(-expoente/escal)
   ddecds2=ddecds2-(2.0d0/escal)*CV(j)*vecp(iref,2,j)/valp(iref,j)*exp(-expoente/escal)
   ddecdrho1=ddecdrho1-(2.0d0/escal)*CV(j)*vecp(iref,3,j)/valp(iref,j)*exp(-expoente/escal)
   ddecdrho22=ddecdrho22-(2.0d0/escal)*CV(j)*vecp(iref,4,j)/valp(iref,j)*exp(-expoente/escal)
   ddecdrho32=ddecdrho32-(2.0d0/escal)*CV(j)*vecp(iref,5,j)/valp(iref,j)*exp(-expoente/escal)
   ddecdrho42=ddecdrho42-(2.0d0/escal)*CV(j)*vecp(iref,6,j)/valp(iref,j)*exp(-expoente/escal)

end do

dv4ds1= dpolds1*dec(iref) + pol(iref)*ddecds1
dv4ds2= dpolds2*dec(iref) + pol(iref)*ddecds2
dv4drho1= dpoldrho1*dec(iref) + pol(iref)*ddecdrho1
dv4drho22= dpoldrho22*dec(iref) + pol(iref)*ddecdrho22
dv4drho32= dpoldrho32*dec(iref) + pol(iref)*ddecdrho32
dv4drho42= dpoldrho42*dec(iref) + pol(iref)*ddecdrho42
dv4drho234= dpoldrho234*dec(iref)


dv4dr1(iref)=dv4ds1*ds1dr1
dv4dr2(iref)=dv4ds2*ds2dr2
dv4dr3(iref)=dv4drho1*drho1dr3 + dv4drho22*drho22dr3 + dv4drho32*drho32dr3 + dv4drho42*drho42dr3 + &
         dv4drho234*drho234dr3
dv4dr4(iref)=dv4drho1*drho1dr4 + dv4drho22*drho22dr4 + dv4drho32*drho32dr4 + dv4drho42*drho42dr4 + &
         dv4drho234*drho234dr4
dv4dr5(iref)=dv4drho1*drho1dr5 + dv4drho22*drho22dr5 + dv4drho32*drho32dr5 + dv4drho42*drho42dr5 + &
         dv4drho234*drho234dr5
dv4dr6(iref)=dv4drho1*drho1dr6 + dv4drho22*drho22dr6 + dv4drho32*drho32dr6 + dv4drho42*drho42dr6 + &
         dv4drho234*drho234dr6

end do references

F1=0.0d0  
F2=0.0d0  
F3=0.0d0  
F4=0.0d0  
F5=0.0d0  
F6=0.0d0  
  
do iref=1,nref

   F1= F1+ dv4dr1(iref)
   F2= F2+ dv4dr2(iref)
   F3= F3+ dv4dr3(iref)
   F4= F4+ dv4dr4(iref)
   F5= F5+ dv4dr5(iref)
   F6= F6+ dv4dr6(iref)

end do

END


!*******************************************************************
!  To compute the electrostatic interaction between two OH diatomics
!*******************************************************************
  
      FUNCTION Vdip_f31(r1,r2,r3,r4,r5,r6)
!     **************************************************************
!  To compute the electrostatic interaction between two OH diatomics
!
      IMPLICIT none
      real*8 vdipOHOH2,r1,r2,r3,r4,r5,r6
      real*8 dip1,dip2,FOHPdip_f31,fdampdip_f31
      real*8 qoa,qob,qha,qhb,vdipOHOH3,Vdip_f31

!     **************************************************************
!     ************V EM Eh***************
      dip1=FOHPdip_f31(r3)
      dip2=FOHPdip_f31(r4)
      qha=dip1/r3
      qoa=-qha
      qhb=dip2/r4
      qob=-qhb
      VdipOHOH2=qoa*qob/r1*fdampdip_f31(r1)+qha*qhb/r2*fdampdip_f31(r2)+  &
                qob*qha/r5*fdampdip_f31(r5)+qoa*qhb/r6*fdampdip_f31(r6)
      dip1=FOHPdip_f31(r5)
      dip2=FOHPdip_f31(r6)
      qha=dip1/r5
      qob=-qha
      qhb=dip2/r6
      qoa=-qhb
      VdipOHOH3=qoa*qob/r1*fdampdip_f31(r1)+qha*qhb/r2*fdampdip_f31(r2)+  &
                qoa*qha/r3*fdampdip_f31(r3)+qob*qhb/r4*fdampdip_f31(r4)
      Vdip_f31=VdipOHOH2+VdipOHOH3

RETURN
END



  subroutine dVdip_f31(r1,r2,r3,r4,r5,r6,dr1,dr2,dr3,dr4,dr5,dr6)
!     **************************************************************
!  To compute the electrostatic interaction between two OH diatomics
!
      IMPLICIT none
      real*8 r1,r2,r3,r4,r5,r6
      real*8 dip1,dip2,FOHPdip_f31,fdampdip_f31
      real*8 dfdampdip_f31,dFOHPdip_f31
      real*8 qoa,qob,qha,qhb
      real*8 dr1,dr2,dr3,dr4,dr5,dr6
      real*8 dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6
      real*8 dv3dr1,dv3dr2,dv3dr3,dv3dr4,dv3dr5,dv3dr6
      real*8 dqhadr3,dqoadr3,dqhbdr4,dqobdr4
      real*8 dqhadr5, dqobdr5,dqhbdr6,dqoadr6

!     **************************************************************
!     ************V EM Eh***************
      dip1=FOHPdip_f31(r3)
      dip2=FOHPdip_f31(r4)
      qha=dip1/r3
      qoa=-qha
      qhb=dip2/r4
      qob=-qhb
!      VdipOHOH2=qoa*qob/r1*fdampdip_f31(r1)+qha*qhb/r2*fdampdip_f31(r2)+  &
!                qob*qha/r5*fdampdip_f31(r5)+qoa*qhb/r6*fdampdip_f31(r6)
      
      dqhadr3=dFOHPdip_f31(r3)/r3-dip1/r3**2
      dqoadr3=-dqhadr3
      dqhbdr4=dFOHPdip_f31(r4)/r4-dip2/r4**2
      dqobdr4=-dqhbdr4

      dv2dr1=-qoa*qob/r1**2*fdampdip_f31(r1)+qoa*qob/r1*dfdampdip_f31(r1)
      dv2dr2=-qha*qhb/r2**2*fdampdip_f31(r2)+qha*qhb/r2*dfdampdip_f31(r2)
      dv2dr3=dqoadr3*qob/r1*fdampdip_f31(r1)+dqhadr3*qhb/r2*fdampdip_f31(r2)+&
             dqhadr3*qob/r5*fdampdip_f31(r5)+dqoadr3*qhb/r6*fdampdip_f31(r6)
      dv2dr4=dqobdr4*qoa/r1*fdampdip_f31(r1)+dqhbdr4*qha/r2*fdampdip_f31(r2)+&
             dqobdr4*qha/r5*fdampdip_f31(r5)+dqhbdr4*qoa/r6*fdampdip_f31(r6)
      dv2dr5=-qob*qha/r5**2*fdampdip_f31(r5)+qob*qha/r5*dfdampdip_f31(r5)
      dv2dr6=-qoa*qhb/r6**2*fdampdip_f31(r6)+qoa*qhb/r6*dfdampdip_f31(r6)

      dip1=FOHPdip_f31(r5)
      dip2=FOHPdip_f31(r6)
      qha=dip1/r5
      qob=-qha
      qhb=dip2/r6
      qoa=-qhb
!      VdipOHOH3=qoa*qob/r1*fdampdip_f31(r1)+qha*qhb/r2*fdampdip_f31(r2)+  &
!                qoa*qha/r3*fdampdip_f31(r3)+qob*qhb/r4*fdampdip_f31(r4)

      
      dqhadr5=-dip1/r5**2+dFOHPdip_f31(r5)/r5
      dqobdr5=-dqhadr5
      dqhbdr6=-dip2/r6**2+dFOHPdip_f31(r6)/r6
      dqoadr6=-dqhbdr6
      

      dv3dr1=-qoa*qob/r1**2*fdampdip_f31(r1)+qoa*qob/r1*dfdampdip_f31(r1)
      dv3dr2=-qha*qhb/r2**2*fdampdip_f31(r2)+qha*qhb/r2*dfdampdip_f31(r2)
      dv3dr3=-qoa*qha/r3**2*fdampdip_f31(r3)+qoa*qha/r3*dfdampdip_f31(r3)
      dv3dr4=-qob*qhb/r4**2*fdampdip_f31(r4)+qob*qhb/r4*dfdampdip_f31(r4)
      dv3dr5=dqobdr5*qoa/r1*fdampdip_f31(r1)+dqhadr5*qhb/r2*fdampdip_f31(r2)+&
             dqhadr5*qoa/r3*fdampdip_f31(r3)+dqobdr5*qhb/r4*fdampdip_f31(r4)
      dv3dr6=dqoadr6*qob/r1*fdampdip_f31(r1)+dqhbdr6*qha/r2*fdampdip_f31(r2)+&
             dqoadr6*qha/r3*fdampdip_f31(r3)+dqhbdr6*qob/r4*fdampdip_f31(r4)


!      Vdip_f31=VdipOHOH2+VdipOHOH3

      dr1=dv2dr1+dv3dr1
      dr2=dv2dr2+dv3dr2
      dr3=dv2dr3+dv3dr3
      dr4=dv2dr4+dv3dr4
      dr5=dv2dr5+dv3dr5
      dr6=dv2dr6+dv3dr6
      
RETURN
END


      FUNCTION FOHPdip_f31(R)
!     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7/0.122307769D+01,-0.975516629D+00,   &
           0.447656609D+00,-0.928112515D-01,0.350125118D+00,  & 
           0.371793971D-05,0.331141524D+03/ 
!     ******************************************************************
      rr2=r*r
      rr3=rr2*r
      rr6=rr3*rr3
      rr8=rr6*rr2      
      FOHPdip_f31=(C1*R+C2*RR2+C3*RR3)*EXP(-c4*R-C5*RR2)+(1-EXP(-C6*   &
         RR8))*C7/RR6
      return
      end

      
      FUNCTION DFOHPdip_f31(R)
!     ******************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)      
      DATA C1,C2,C3,C4,C5,C6,C7/0.122307769D+01,-0.975516629D+00,   &
           0.447656609D+00,-0.928112515D-01,0.350125118D+00,  & 
           0.371793971D-05,0.331141524D+03/ 
!     ******************************************************************
      rr2=r*r
      rr3=rr2*r
      rr6=rr3*rr3
      rr8=rr6*rr2 

      DFOHPdip_f31=(C1+2.0D0*C2*R+3.0D0*C3*RR2)*EXP(-C4*R-C5*RR2)-(C4+   &
               2.0D0*C5*R)*(C1*R+C2*RR2+C3*RR3)*EXP(-C4*R-C5*RR2)  &
               -(6.0D0*C7/(RR6*R))*(1.0d0-EXP(-C6*RR8))+ &
               8.0D0*C6*C7*R*EXP(-C6*RR8)
     
      return
      end


Function fdampdip_f31(r)
!     **************************************************************
!  To damp the electrostatic interaction due to orbital overlap
!
IMPLICIT none
real*8 fdampdip_f31,r
real*8 a(20),b(20),R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS
real*8 rmohohp,r0ohohp,ro
!
      DATA R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS/6.928203D0,   &
         1.40100D0,7.82D0,6.294894D0,6.334299622D0,1.8344D0,1.9086D0/   
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
RMOHOHP=RMOHP*1.5d0
R0OHOHP=r0OHP*1.5d0
RO=0.5D0*(RMOHOHP+2.5D0*R0OHOHP)

fdampdip_f31=(1-exp(-A(3)*(r/ro)-B(3)*(r/ro)**2))

return
end  


Function Dfdampdip_f31(r)
!     **************************************************************
!  To damp the electrostatic interaction due to orbital overlap
!
IMPLICIT none
real*8 Dfdampdip_f31,r
real*8 a(20),b(20),R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS
real*8 rmohohp,r0ohohp,ro
!
      DATA R0HH,RMHHS,RMHHT,R0OHP,R0OHS,RMOHP,RMOHS/6.928203D0,   &
         1.40100D0,7.82D0,6.294894D0,6.334299622D0,1.8344D0,1.9086D0/   
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
RMOHOHP=RMOHP*1.5d0
R0OHOHP=r0OHP*1.5d0
RO=0.5D0*(RMOHOHP+2.5D0*R0OHOHP)

Dfdampdip_f31=(A(3)/ro+2.0d0*B(3)*(r/ro**2))*exp(-A(3)*(r/ro)-B(3)*  &
          (r/ro)**2)   

return
end 




























