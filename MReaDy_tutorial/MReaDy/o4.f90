subroutine VO4_f11_switch(i,V,dv)
use variables
use phys_parameters
use sim_variables
use constants
double precision, dimension(12) :: dV1,dV
integer      :: gi3,gi4,gi5,gi6
integer      :: i!,j,k,l
!double precision, dimension(3) :: acc,bcc,ccc,dcc
!double precision :: dist
double precision :: rxij,ryij,rzij,rij
double precision :: r1xij,r1yij,r1zij
double precision :: r2xij,r2yij,r2zij
double precision :: r3xij,r3yij,r3zij
double precision :: r4xij,r4yij,r4zij
double precision :: r5xij,r5yij,r5zij
double precision :: r6xij,r6yij,r6zij
double precision :: dv2dr1,dv2dr2,dv2dr3,dv2dr4,dv2dr5,dv2dr6
double precision :: r1,r2,r3,r4,r5,r6

double precision, dimension(6) :: rin!,fn
double precision :: fs!gama,rzero,fs
!double precision, dimension(3) :: rH3,dvdrV2,dvdrH3
double precision :: V2,V1,V,dV2,dV3,dvdiat,V4 ,vdiat ,v3
double precision :: f3,f4
!double precision :: dr1,dr2,dr3,dr4
double precision :: dfsdr1,dfsdr2,dfsdr3,dfsdr4,dfsdr5,dfsdr6
double precision, dimension(4,3) :: vect_4atoms ,dvdcc
!double precision, dimension(6) ::r_4atoms,dvdr_4atoms
double precision :: VO4_11,DVO4_11,vo3_f7
double precision :: dv1dr1,dv1dr2,dv1dr3,dv1dr4,dv1dr5,dv1dr6



fs=0.0d0
V=0.0d0
V1=0.0d0
V2=0.0d0
V3=0.0d0
dv1=0.0d0
dv2=0.0d0
dv3=0.0d0
dv4=0.0d0
dv=0.0d0
V4=0.0d0
vect_4atoms=0.0d0
dvdcc=0.0d0

dfsdr1=0.0d0
dfsdr2=0.0d0
dfsdr3=0.0d0
dfsdr4=0.0d0
dfsdr5=0.0d0
dfsdr6=0.0d0

r1xij=0.0d0
r1yij=0.0d0
r1zij=0.0d0

r2xij=0.0d0
r2yij=0.0d0
r2zij=0.0d0

r3xij=0.0d0
r3yij=0.0d0
r3zij=0.0d0

r4xij=0.0d0
r4yij=0.0d0
r4zij=0.0d0

r5xij=0.0d0
r5yij=0.0d0
r5zij=0.0d0

r6xij=0.0d0
r6yij=0.0d0
r6zij=0.0d0

dv2dr1=0.0d0
dv2dr2=0.0d0
dv2dr3=0.0d0
dv2dr4=0.0d0
dv2dr5=0.0d0
dv2dr6=0.0d0

dv1dr1=0.0d0
dv1dr2=0.0d0
dv1dr3=0.0d0
dv1dr4=0.0d0
dv1dr5=0.0d0
dv1dr6=0.0d0




!Geting the coordenates

!The positions are predefined in the matrix "group"
            !The positions are predefined in the matrix "group"
            gi3=group(i,3)      ! O atom
            gi4=group(i,4)      ! O atom
            gi5=group(i,5)      ! O atom
            gi6=group(i,6)      ! O atom



!      Print*,rx(gi3),ry(gi3),rz(gi3)
!      Print*,rx(gi4),ry(gi4),rz(gi4)
!      Print*,rx(gi5),ry(gi5),rz(gi5)
!      Print*,rx(gi6),ry(gi6),rz(gi6)



















! Checking inner distance



call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi4),ry(gi4), &
rz(gi4),rxij,ryij,rzij,rij)

rin(2)=rij
r2=rij
r2xij=rxij
r2yij=ryij
r2zij=rzij

call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi5),ry(gi5), &
rz(gi5),rxij,ryij,rzij,rij)

rin(3)=rij
r3=rij
r3xij=rxij
r3yij=ryij
r3zij=rzij


call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi6),ry(gi6), &
rz(gi6),rxij,ryij,rzij,rij)

rin(5)=rij
r5=rij
r5xij=rxij
r5yij=ryij
r5zij=rzij


call mod_dist(rx(gi4),ry(gi4),rz(gi4),rx(gi5),ry(gi5), &
rz(gi5),rxij,ryij,rzij,rij)

rin(6)=rij
r6=rij
r6xij=rxij
r6yij=ryij
r6zij=rzij



call mod_dist(rx(gi4),ry(gi4),rz(gi4),rx(gi6),ry(gi6), &
rz(gi6),rxij,ryij,rzij,rij)

rin(4)=rij
r4=rij
r4xij=rxij
r4yij=ryij
r4zij=rzij

call mod_dist(rx(gi5),ry(gi5),rz(gi5),rx(gi6),ry(gi6), &
rz(gi6),rxij,ryij,rzij,rij)

rin(1)=rij
r1=rij
r1xij=rxij
r1yij=ryij
r1zij=rzij




V1=VO4_11(R2,R4,R1,R3,R5,R6)

dv1dr2=DVO4_11(R2,R4,R1,R3,R5,R6,1)
dv1dr4=DVO4_11(R2,R4,R1,R3,R5,R6,2)
dv1dr1=DVO4_11(R2,R4,R1,R3,R5,R6,3)
dv1dr3=DVO4_11(R2,R4,R1,R3,R5,R6,4)
dv1dr5=DVO4_11(R2,R4,R1,R3,R5,R6,5)
dv1dr6=DVO4_11(R2,R4,R1,R3,R5,R6,6)



dv1(1)  =  dv1(1)  + dv1dr2 * r2xij / r2 + dv1dr3 * r3xij / r3 + dv1dr5 * r5xij / r5
dv1(2)  =  dv1(2)  + dv1dr2 * r2yij / r2 + dv1dr3 * r3yij / r3 + dv1dr5 * r5yij / r5
dv1(3)  =  dv1(3)  + dv1dr2 * r2zij / r2 + dv1dr3 * r3zij / r3 + dv1dr5 * r5zij / r5

dv1(4)  =  dv1(4)  - dv1dr2 * r2xij / r2 + dv1dr4 * r4xij / r4 + dv1dr6 * r6xij / r6
dv1(5)  =  dv1(5)  - dv1dr2 * r2yij / r2 + dv1dr4 * r4yij / r4 + dv1dr6 * r6yij / r6
dv1(6)  =  dv1(6)  - dv1dr2 * r2zij / r2 + dv1dr4 * r4zij / r4 + dv1dr6 * r6zij / r6

dv1(7)  =  dv1(7)  + dv1dr1 * r1xij / r1 - dv1dr3 * r3xij / r3 - dv1dr6 * r6xij / r6
dv1(8)  =  dv1(8)  + dv1dr1 * r1yij / r1 - dv1dr3 * r3yij / r3 - dv1dr6 * r6yij / r6
dv1(9)  =  dv1(9)  + dv1dr1 * r1zij / r1 - dv1dr3 * r3zij / r3 - dv1dr6 * r6zij / r6

dv1(10) =  dv1(10) - dv1dr1 * r1xij / r1 - dv1dr4 * r4xij / r4 - dv1dr5 * r5xij / r5
dv1(11) =  dv1(11) - dv1dr1 * r1yij / r1 - dv1dr4 * r4yij / r4 - dv1dr5 * r5yij / r5
dv1(12) =  dv1(12) - dv1dr1 * r1zij / r1 - dv1dr4 * r4zij / r4 - dv1dr5 * r5zij / r5


!          print*,'V1A',V1/con
!          stop


!'O4 -> O(4) + O3'

if ((r1.ge.in_bound).and.(r4.ge.in_bound).and. &
(r5.ge.in_bound)) then
if(group(i,9).eq.0) group(i,9)=11
if(group(i,9).ne.11) then

print*,'exchange in molecular complex 1 ',tempo
group(i,9)=11
end if

!            print*,'tempo',tempo

! O3 potential

V2=vo3_f7(r2,r3,r6)
call dero3_f7(r2,r3,r6,dv2dr2,dv2dr3,dv2dr6)


! Switch function


fs=f3(rin(1),rin(4),rin(5))
call df3(rin(1),rin(4),rin(5),dfsdr1,dfsdr4,dfsdr5)



! for the diatomic potencials

call poto2q_f18(rin(5),Vdiat,dVdiat)

!            print*,'rin(5),Vdiat'
!            print*,rin(5),Vdiat



dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
dv(3) = dv(3)  + fs * dVdiat*r5zij/r5

dv(10)= dv(10) - fs * dVdiat*r5xij/r5
dv(11)= dv(11) - fs * dVdiat*r5yij/r5
dv(12)= dv(12) - fs * dVdiat*r5zij/r5


V4 = V4 + Vdiat


call poto2q_f18(rin(4),Vdiat,dVdiat)

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat



dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
dv(6) = dv(6)  + fs * dVdiat*r4zij/r4

dv(10)= dv(10) - fs * dVdiat*r4xij/r4
dv(11)= dv(11) - fs * dVdiat*r4yij/r4
dv(12)= dv(12) - fs * dVdiat*r4zij/r4


V4 = V4 + Vdiat



call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(1),Vdiat'
!            print*,rin(1),Vdiat



dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

dv(10)= dv(10) - fs * dVdiat*r1xij/r1
dv(11)= dv(11) - fs * dVdiat*r1yij/r1
dv(12)= dv(12) - fs * dVdiat*r1zij/r1


V4 = V4 + Vdiat

V=V1*(1-fs)+(V2+V4)*fs
!
dv=dv + dv1*(1-fs)




dv(1)  = dv(1) + (v4 + v2 - v1) * ( dfsdr5 * r5xij / r5)
dv(2)  = dv(2) + (v4 + v2 - v1) * ( dfsdr5 * r5yij / r5)
dv(3)  = dv(3) + (v4 + v2 - v1) * ( dfsdr5 * r5zij / r5)

dv(4)  = dv(4) + (v4 + v2 - v1) * ( dfsdr4 * r4xij / r4)
dv(5)  = dv(5) + (v4 + v2 - v1) * ( dfsdr4 * r4yij / r4)
dv(6)  = dv(6) + (v4 + v2 - v1) * ( dfsdr4 * r4zij / r4)

dv(7)  = dv(7) + (v4 + v2 - v1) * ( dfsdr1 * r1xij / r1)
dv(8)  = dv(8) + (v4 + v2 - v1) * ( dfsdr1 * r1yij / r1)
dv(9)  = dv(9) + (v4 + v2 - v1) * ( dfsdr1 * r1zij / r1)

dv(10) = dv(10) + (v4 +v2 - v1) * (-dfsdr1 * r1xij / r1 - dfsdr4 * r4xij / r4 - dfsdr5 * r5xij / r5)
dv(11) = dv(11) + (v4 +v2 - v1) * (-dfsdr1 * r1yij / r1 - dfsdr4 * r4yij / r4 - dfsdr5 * r5yij / r5)
dv(12) = dv(12) + (v4 +v2 - v1) * (-dfsdr1 * r1zij / r1 - dfsdr4 * r4zij / r4 - dfsdr5 * r5zij / r5)


dv(1)  = dv(1)  + fs * ( dv2dr2 * r2xij / r2 + dv2dr3 * r3xij /r3)
dv(2)  = dv(2)  + fs * ( dv2dr2 * r2yij / r2 + dv2dr3 * r3yij /r3)
dv(3)  = dv(3)  + fs * ( dv2dr2 * r2zij / r2 + dv2dr3 * r3zij /r3)

dv(4)  = dv(4)  + fs * (-dv2dr2 * r2xij / r2 + dv2dr6 * r6xij /r6)
dv(5)  = dv(5)  + fs * (-dv2dr2 * r2yij / r2 + dv2dr6 * r6yij /r6)
dv(6)  = dv(6)  + fs * (-dv2dr2 * r2zij / r2 + dv2dr6 * r6zij /r6)

dv(7)  = dv(7)  + fs * (-dv2dr3 * r3xij / r3 - dv2dr6 * r6xij /r6)
dv(8)  = dv(8)  + fs * (-dv2dr3 * r3yij / r3 - dv2dr6 * r6yij /r6)
dv(9)  = dv(9)  + fs * (-dv2dr3 * r3zij / r3 - dv2dr6 * r6zij /r6)


!       print*,'dv',dv


!O4 -> O(1) + O3
else if ((r2.ge.in_bound).and.(r3.ge.in_bound).and. &
(r5.ge.in_bound)) then

if(group(i,9).eq.0)  group(i,9)=12

if(group(i,9).ne.12) then
print*,'Change in group'
group(i,9)=12
end if

!        print*,'here we are'
!O4 -> O(1) + O3

!        print*,"mid O(1) + O3",tempo
!        print*,'Exchange in channels. Time:',t



! O3 potential

V2=vo3_f7(r6,r4,r1)

call dero3_f7(r6,r4,r1,dv2dr6,dv2dr4,dv2dr1)




! Switch function


fs=f3(rin(2),rin(3),rin(5))
call df3(rin(2),rin(3),rin(5),dfsdr2,dfsdr3,dfsdr5)


!        print*,'fs=',fs,dfsdr2,dfsdr3,dfsdr5





! for the diatomic potencials

call poto2q_f18(rin(5),Vdiat,dVdiat)

!            print*,'rin(5),Vdiat'
!            print*,rin(5),Vdiat



dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
dv(3) = dv(3)  + fs * dVdiat*r5zij/r5

dv(10)= dv(10) - fs * dVdiat*r5xij/r5
dv(11)= dv(11) - fs * dVdiat*r5yij/r5
dv(12)= dv(12) - fs * dVdiat*r5zij/r5


V4 = V4 + Vdiat


call poto2q_f18(rin(2),Vdiat,dVdiat)

!            print*,'rin(2),Vdiat'
!            print*,rin(2),Vdiat

dv(1) = dv(1) + fs * dVdiat*r2xij/r2
dv(2) = dv(2) + fs * dVdiat*r2yij/r2
dv(3) = dv(3) + fs * dVdiat*r2zij/r2

dv(4) = dv(4) - fs * dVdiat*r2xij/r2
dv(5) = dv(5) - fs * dVdiat*r2yij/r2
dv(6) = dv(6) - fs * dVdiat*r2zij/r2


V4 = V4 + Vdiat


call poto2q_f18(rin(3),Vdiat,dVdiat)
!            print*,'rin(3),Vdiat'
!            print*,rin(3),Vdiat



dv(1) = dv(1) + fs * dVdiat*r3xij/r3
dv(2) = dv(2) + fs * dVdiat*r3yij/r3
dv(3) = dv(3) + fs * dVdiat*r3zij/r3

dv(7) = dv(7) - fs * dVdiat*r3xij/r3
dv(8) = dv(8) - fs * dVdiat*r3yij/r3
dv(9) = dv(9) - fs * dVdiat*r3zij/r3


V4 = V4 + Vdiat

V=V1*(1-fs)+(V2+V4)*fs
!        print*,rin
!        print*,'V=',V,' V1=',V1,' V2=',V2

dv=dv + dv1*(1-fs)


dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr2 *  r2xij / r2 + dfsdr3 * r3xij / r3 + dfsdr5 * r5xij / r5)
dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr2 *  r2yij / r2 + dfsdr3 * r3yij / r3 + dfsdr5 * r5yij / r5)
dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr2 *  r2zij / r2 + dfsdr3 * r3zij / r3 + dfsdr5 * r5zij / r5)

dv(4) = dv(4) + (v4 + v2 - v1) * (-dfsdr2 *  r2xij / r2)
dv(5) = dv(5) + (v4 + v2 - v1) * (-dfsdr2 *  r2yij / r2)
dv(6) = dv(6) + (v4 + v2 - v1) * (-dfsdr2 *  r2zij / r2)

dv(7) = dv(7) + (v4 + v2 - v1) * (-dfsdr3 *  r3xij / r3)
dv(8) = dv(8) + (v4 + v2 - v1) * (-dfsdr3 *  r3yij / r3)
dv(9) = dv(9) + (v4 + v2 - v1) * (-dfsdr3 *  r3zij / r3)

dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr5 *  r5xij / r5)
dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr5 *  r5yij / r5)
dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr5 *  r5zij / r5)


dv(4) = dv(4) + fs * ( dv2dr4 * r4xij / r4 + dv2dr6 * r6xij / r6)
dv(5) = dv(5) + fs * ( dv2dr4 * r4yij / r4 + dv2dr6 * r6yij / r6)
dv(6) = dv(6) + fs * ( dv2dr4 * r4zij / r4 + dv2dr6 * r6zij / r6)

dv(7) = dv(7)  + fs * ( dv2dr1 * r1xij / r1 - dv2dr6 * r6xij / r6)
dv(8) = dv(8)  + fs * ( dv2dr1 * r1yij / r1 - dv2dr6 * r6yij / r6)
dv(9) = dv(9)  + fs * ( dv2dr1 * r1zij / r1 - dv2dr6 * r6zij / r6)

dv(10)= dv(10)  + fs * (-dv2dr1 * r1xij / r1 - dv2dr4 * r4xij / r4)
dv(11)= dv(11)  + fs * (-dv2dr1 * r1yij / r1 - dv2dr4 * r4yij / r4)
dv(12)= dv(12)  + fs * (-dv2dr1 * r1zij / r1 - dv2dr4 * r4zij / r4)




!'O4 -> O(2) + O3
else if ((r2.ge.in_bound).and.(r4.ge.in_bound).and. &
(r6.ge.in_bound)) then

if(group(i,9).eq.0)  group(i,9)=13

if(group(i,9).ne.13) then
print*,'Change in group'
group(i,9)=13
end if


!        print*,'here we are'
!O4 -> O(2) + O3

!             print*,"mid O(2) + O3",tempo
!print*,'Exchange in channels. Time:',t

! O3 potential




V2=vo3_f7(r6,r1,r4)
call dero3_f7(r1,r3,r5,dv2dr1,dv2dr3,dv2dr5)












! Switch function


fs=f3(rin(2),rin(4),rin(6))
call df3(rin(2),rin(4),rin(6),dfsdr2,dfsdr4,dfsdr6)


!        print*,'fs=',fs,dfsdr2,dfsdr4,dfsdr6





! for the diatomic potencials

call poto2q_f18(rin(4),Vdiat,dVdiat)

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat



dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
dv(6) = dv(6)  + fs * dVdiat*r4zij/r4

dv(10)= dv(10) - fs * dVdiat*r4xij/r4
dv(11)= dv(11) - fs * dVdiat*r4yij/r4
dv(12)= dv(12) - fs * dVdiat*r4zij/r4


V4 = V4 + Vdiat


call poto2q_f18(rin(2),Vdiat,dVdiat)


!            print*,'rin(2),Vdiat'
!            print*,rin(2),Vdiat

dv(1) = dv(1) + fs * dVdiat*r2xij/r2
dv(2) = dv(2) + fs * dVdiat*r2yij/r2
dv(3) = dv(3) + fs * dVdiat*r2zij/r2

dv(4) = dv(4) - fs * dVdiat*r2xij/r2
dv(5) = dv(5) - fs * dVdiat*r2yij/r2
dv(6) = dv(6) - fs * dVdiat*r2zij/r2


V4 = V4 + Vdiat


call poto2q_f18(rin(6),Vdiat,dVdiat)

!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat



dv(4) = dv(4) + fs * dVdiat*r6xij/r6
dv(5) = dv(5) + fs * dVdiat*r6yij/r6
dv(6) = dv(6) + fs * dVdiat*r6zij/r6

dv(7) = dv(7) - fs * dVdiat*r6xij/r6
dv(8) = dv(8) - fs * dVdiat*r6yij/r6
dv(9) = dv(9) - fs * dVdiat*r6zij/r6


V4 = V4 + Vdiat

V=V1*(1-fs)+(V2+V4)*fs
!        print*,rin
!        print*,'V=',V,' V1=',V1,' V2=',V2

dv=dv + dv1*(1-fs)


dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr2 *  r2xij / r2)
dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr2 *  r2yij / r2)
dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr2 *  r2zij / r2)

dv(4) = dv(4) + (v4 + v2 - v1) * (-dfsdr2 *  r2xij / r2 + dfsdr4 *  r4xij / r4 + dfsdr6 *  r6xij / r6)
dv(5) = dv(5) + (v4 + v2 - v1) * (-dfsdr2 *  r2yij / r2 + dfsdr4 *  r4yij / r4 + dfsdr6 *  r6yij / r6)
dv(6) = dv(6) + (v4 + v2 - v1) * (-dfsdr2 *  r2zij / r2 + dfsdr4 *  r4zij / r4 + dfsdr6 *  r6zij / r6)

dv(7) = dv(7) + (v4 + v2 - v1) * (-dfsdr6 *  r6xij / r6)
dv(8) = dv(8) + (v4 + v2 - v1) * (-dfsdr6 *  r6yij / r6)
dv(9) = dv(9) + (v4 + v2 - v1) * (-dfsdr6 *  r6zij / r6)

dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr4 *  r4xij / r4)
dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr4 *  r4yij / r4)
dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr4 *  r4zij / r4)


dv(1) = dv(1) + fs * ( dv2dr3 * r3xij / r3 + dv2dr5 * r5xij / r5)
dv(2) = dv(2) + fs * ( dv2dr3 * r3yij / r3 + dv2dr5 * r5yij / r5)
dv(3) = dv(3) + fs * ( dv2dr3 * r3zij / r3 + dv2dr5 * r5zij / r5)

dv(7)  = dv(7)  + fs * ( dv2dr1 * r1xij / r1 - dv2dr3 * r3xij / r3)
dv(8)  = dv(8)  + fs * ( dv2dr1 * r1yij / r1 - dv2dr3 * r3yij / r3)
dv(9)  = dv(9)  + fs * ( dv2dr1 * r1zij / r1 - dv2dr3 * r3zij / r3)

dv(10)  = dv(10)  + fs * (-dv2dr1 * r1xij / r1 - dv2dr5 * r5xij / r5)
dv(11)  = dv(11)  + fs * (-dv2dr1 * r1yij / r1 - dv2dr5 * r5yij / r5)
dv(12)  = dv(12)  + fs * (-dv2dr1 * r1zij / r1 - dv2dr5 * r5zij / r5)

!'O -> O(3) + O3
else if ((r1.ge.in_bound).and.(r3.ge.in_bound).and. &
(r6.ge.in_bound)) then

if(group(i,9).eq.0)  group(i,9)=14

if(group(i,9).ne.14) then
print*,'Change in group'
group(i,9)=14
end if
!        print*,'here we are'
!O4 -> O(3) + O3




!             print*,"mid O(2) + O3",tempo
!print*,'Exchange in channels. Time:',t

! O3 potential

V2=vo3_f7(r2,r4,r5)
call dero3_f7(r2,r4,r5,dv2dr2,dv2dr4,dv2dr5)


! Switch function


fs=f3(rin(1),rin(3),rin(6))
call df3(rin(1),rin(3),rin(6),dfsdr1,dfsdr3,dfsdr6)


!        print*,'fs=',fs,dfsdr1,dfsdr3,dfsdr6





! for the diatomic potencials

call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(1),Vdiat'
!            print*,rin(1),Vdiat



dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

dv(10)= dv(10) - fs * dVdiat*r1xij/r1
dv(11)= dv(11) - fs * dVdiat*r1yij/r1
dv(12)= dv(12) - fs * dVdiat*r1zij/r1


V4 = V4 + Vdiat


call poto2q_f18(rin(3),Vdiat,dVdiat)


!            print*,'rin(3),Vdiat'
!            print*,rin(3),Vdiat

dv(1) = dv(1) + fs * dVdiat*r3xij/r3
dv(2) = dv(2) + fs * dVdiat*r3yij/r3
dv(3) = dv(3) + fs * dVdiat*r3zij/r3

dv(7) = dv(7) - fs * dVdiat*r3xij/r3
dv(8) = dv(8) - fs * dVdiat*r3yij/r3
dv(9) = dv(9) - fs * dVdiat*r3zij/r3


V4 = V4 + Vdiat


call poto2q_f18(rin(6),Vdiat,dVdiat)

!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat



dv(4) = dv(4) + fs * dVdiat*r6xij/r6
dv(5) = dv(5) + fs * dVdiat*r6yij/r6
dv(6) = dv(6) + fs * dVdiat*r6zij/r6

dv(7) = dv(7) - fs * dVdiat*r6xij/r6
dv(8) = dv(8) - fs * dVdiat*r6yij/r6
dv(9) = dv(9) - fs * dVdiat*r6zij/r6


V4 = V4 + Vdiat

V=V1*(1-fs)+(V2+V4)*fs
!        print*,rin
!        print*,'V=',V,' V1=',V1,' V2=',V2

dv=dv + dv1*(1-fs)


dv(1) = dv(1) + (v4 + v2 - v1) * ( dfsdr3 *  r3xij / r3)
dv(2) = dv(2) + (v4 + v2 - v1) * ( dfsdr3 *  r3yij / r3)
dv(3) = dv(3) + (v4 + v2 - v1) * ( dfsdr3 *  r3zij / r3)

dv(4) = dv(4) + (v4 + v2 - v1) * ( dfsdr6 *  r6xij / r6)
dv(5) = dv(5) + (v4 + v2 - v1) * ( dfsdr6 *  r6yij / r6)
dv(6) = dv(6) + (v4 + v2 - v1) * ( dfsdr6 *  r6zij / r6)

dv(7) = dv(7) + (v4 + v2 - v1) * ( dfsdr1 *  r1xij / r1 - dfsdr3 *  r3xij / r3 - dfsdr6 *  r6xij / r6)
dv(8) = dv(8) + (v4 + v2 - v1) * ( dfsdr1 *  r1yij / r1 - dfsdr3 *  r3yij / r3 - dfsdr6 *  r6yij / r6)
dv(9) = dv(9) + (v4 + v2 - v1) * ( dfsdr1 *  r1zij / r1 - dfsdr3 *  r3zij / r3 - dfsdr6 *  r6zij / r6)

dv(10)= dv(10)+ (v4 + v2 - v1) * (-dfsdr1 *  r1xij / r1)
dv(11)= dv(11)+ (v4 + v2 - v1) * (-dfsdr1 *  r1yij / r1)
dv(12)= dv(12)+ (v4 + v2 - v1) * (-dfsdr1 *  r1zij / r1)


dv(1) = dv(1) + fs * ( dv2dr2 * r2xij / r2 + dv2dr5 * r5xij / r5)
dv(2) = dv(2) + fs * ( dv2dr2 * r2yij / r2 + dv2dr5 * r5yij / r5)
dv(3) = dv(3) + fs * ( dv2dr2 * r2zij / r2 + dv2dr5 * r5zij / r5)

dv(4)  = dv(4)  + fs * (-dv2dr2 * r2xij / r2 + dv2dr4 * r4xij / r4)
dv(5)  = dv(5)  + fs * (-dv2dr2 * r2yij / r2 + dv2dr4 * r4yij / r4)
dv(6)  = dv(6)  + fs * (-dv2dr2 * r2zij / r2 + dv2dr4 * r4zij / r4)

dv(10)  = dv(10)  + fs * (-dv2dr4 * r4xij / r4 - dv2dr5 * r5xij / r5)
dv(11)  = dv(11)  + fs * (-dv2dr4 * r4yij / r4 - dv2dr5 * r5yij / r5)
dv(12)  = dv(12)  + fs * (-dv2dr4 * r4zij / r4 - dv2dr5 * r5zij / r5)


else if ((r1.ge.in_bound).and.(r2.ge.in_bound).and. &
(r5.ge.in_bound).and.(r6.ge.in_bound)) then
!              print*,'exit'



! 'O4 -> O2(1-3) + O2(2-4)'
if(group(i,9).eq.0)  group(i,9)=15

if(group(i,9).ne.15) then
print*,'Change in group'
!        stop 'here'
group(i,9)=15
end if

! O2 function
call potoop_f2(r3,V2,dV2)
! O2 function
call potoop_f2(r4,V3,dV3)



! Switch function
fs=f4(rin(1),rin(2),rin(5),rin(6))
call df4(rin(1),rin(2),rin(5),rin(6),dfsdr1,dfsdr2,dfsdr5,dfsdr6)



! for the diatomic potencials

call poto2q_f18(rin(5),Vdiat,dVdiat)

!            print*,'rin(5),Vdiat'
!            print*,rin(5),Vdiat



dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
dv(3) = dv(3)  + fs * dVdiat*r5zij/r5

dv(10)= dv(10) - fs * dVdiat*r5xij/r5
dv(11)= dv(11) - fs * dVdiat*r5yij/r5
dv(12)= dv(12) - fs * dVdiat*r5zij/r5


V4 = V4 + Vdiat



call poto2q_f18(rin(2),Vdiat,dVdiat)


!            print*,'rin(2),Vdiat'
!            print*,rin(2),Vdiat

dv(1) = dv(1) + fs * dVdiat*r2xij/r2
dv(2) = dv(2) + fs * dVdiat*r2yij/r2
dv(3) = dv(3) + fs * dVdiat*r2zij/r2

dv(4) = dv(4) - fs * dVdiat*r2xij/r2
dv(5) = dv(5) - fs * dVdiat*r2yij/r2
dv(6) = dv(6) - fs * dVdiat*r2zij/r2


V4 = V4 + Vdiat



!call poto2q_f18(rin(4),Vdiat,dVdiat)
!
!!            print*,'rin(4),Vdiat'
!!            print*,rin(4),Vdiat
!
!
!
!dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
!dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
!dv(6) = dv(6)  + fs * dVdiat*r4zij/r4
!
!dv(10)= dv(10) - fs * dVdiat*r4xij/r4
!dv(11)= dv(11) - fs * dVdiat*r4yij/r4
!dv(12)= dv(12) - fs * dVdiat*r4zij/r4
!
!
!V4 = V4 + Vdiat



call poto2q_f18(rin(6),Vdiat,dVdiat)


!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat

dv(4) = dv(4) + fs * dVdiat*r6xij/r6
dv(5) = dv(5) + fs * dVdiat*r6yij/r6
dv(6) = dv(6) + fs * dVdiat*r6zij/r6

dv(7) = dv(7) - fs * dVdiat*r6xij/r6
dv(8) = dv(8) - fs * dVdiat*r6yij/r6
dv(9) = dv(9) - fs * dVdiat*r6zij/r6


V4 = V4 + Vdiat



call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(1),Vdiat'
!            print*,rin(1),Vdiat



dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

dv(10)= dv(10) - fs * dVdiat*r1xij/r1
dv(11)= dv(11) - fs * dVdiat*r1yij/r1
dv(12)= dv(12) - fs * dVdiat*r1zij/r1


V4 = V4 + Vdiat





V=V1*(1-fs)+(V2+V3+V4)*fs

dv=dv + dv1*(1-fs)

dv(1)  =   dv(1)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2xij / r2 + dfsdr5 * r5xij / r5)
dv(2)  =   dv(2)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2yij / r2 + dfsdr5 * r5yij / r5)
dv(3)  =   dv(3)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2zij / r2 + dfsdr5 * r5zij / r5)

dv(4)  =   dv(4)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2xij / r2 + dfsdr6 * r6xij / r6)
dv(5)  =   dv(5)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2yij / r2 + dfsdr6 * r6yij / r6)
dv(6)  =   dv(6)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2zij / r2 + dfsdr6 * r6zij / r6)

dv(7)  =   dv(7)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1xij / r1 - dfsdr6 * r6xij / r6)
dv(8)  =   dv(8)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1yij / r1 - dfsdr6 * r6yij / r6)
dv(9)  =   dv(9)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1zij / r1 - dfsdr6 * r6zij / r6)

dv(10) =   dv(10) + (v2+v3+v4-v1) * (-dfsdr1 *  r1xij / r1 - dfsdr5 * r5xij / r5)
dv(11) =   dv(11) + (v2+v3+v4-v1) * (-dfsdr1 *  r1yij / r1 - dfsdr5 * r5yij / r5)
dv(12) =   dv(12) + (v2+v3+v4-v1) * (-dfsdr1 *  r1zij / r1 - dfsdr5 * r5zij / r5)


dv(1)  =   dv(1)  + fs * ( dv2 * r3xij / r3)
dv(2)  =   dv(2)  + fs * ( dv2 * r3yij / r3)
dv(3)  =   dv(3)  + fs * ( dv2 * r3zij / r3)

dv(4)  =   dv(4)  + fs * ( dv3 * r4xij / r4)
dv(5)  =   dv(5)  + fs * ( dv3 * r4yij / r4)
dv(6)  =   dv(6)  + fs * ( dv3 * r4zij / r4)

dv(7)  =   dv(7)  + fs * (-dv2 * r3xij / r3)
dv(8)  =   dv(8)  + fs * (-dv2 * r3yij / r3)
dv(9)  =   dv(9)  + fs * (-dv2 * r3zij / r3)

dv(10) =   dv(10) + fs * (-dv3 * r4xij / r4)
dv(11) =   dv(11) + fs * (-dv3 * r4yij / r4)
dv(12) =   dv(12) + fs * (-dv3 * r4zij / r4)









else if ((r1.ge.in_bound).and.(r2.ge.in_bound).and. &
(r3.ge.in_bound).and.(r4.ge.in_bound)) then
!              print*,'exit5'



! 'O4 -> O2(2-3) + O2(1-4)'
if(group(i,9).eq.0)  group(i,9)=16

if(group(i,9).ne.16) then
print*,'Change in group'
group(i,9)=16
end if

! O2 function
call potoop_f2(r6,V2,dV2)
! O2 function
call potoop_f2(r5,V3,dV3)



! Switch function
fs=f4(rin(1),rin(2),rin(3),rin(4))
call df4(rin(1),rin(2),rin(3),rin(4),dfsdr1,dfsdr2,dfsdr3,dfsdr4)



! for the diatomic potencials

call poto2q_f18(rin(2),Vdiat,dVdiat)


!            print*,'rin(2),Vdiat'
!            print*,rin(2),Vdiat

dv(1) = dv(1) + fs * dVdiat*r2xij/r2
dv(2) = dv(2) + fs * dVdiat*r2yij/r2
dv(3) = dv(3) + fs * dVdiat*r2zij/r2

dv(4) = dv(4) - fs * dVdiat*r2xij/r2
dv(5) = dv(5) - fs * dVdiat*r2yij/r2
dv(6) = dv(6) - fs * dVdiat*r2zij/r2


V4 = V4 + Vdiat


call poto2q_f18(rin(3),Vdiat,dVdiat)


!            print*,'rin(3),Vdiat'
!            print*,rin(3),Vdiat

dv(1) = dv(1) + fs * dVdiat*r3xij/r3
dv(2) = dv(2) + fs * dVdiat*r3yij/r3
dv(3) = dv(3) + fs * dVdiat*r3zij/r3

dv(7) = dv(7) - fs * dVdiat*r3xij/r3
dv(8) = dv(8) - fs * dVdiat*r3yij/r3
dv(9) = dv(9) - fs * dVdiat*r3zij/r3


V4 = V4 + Vdiat



call poto2q_f18(rin(4),Vdiat,dVdiat)

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat



dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
dv(6) = dv(6)  + fs * dVdiat*r4zij/r4

dv(10)= dv(10) - fs * dVdiat*r4xij/r4
dv(11)= dv(11) - fs * dVdiat*r4yij/r4
dv(12)= dv(12) - fs * dVdiat*r4zij/r4


V4 = V4 + Vdiat



call poto2q_f18(rin(6),Vdiat,dVdiat)


!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat

dv(4) = dv(4) + fs * dVdiat*r6xij/r6
dv(5) = dv(5) + fs * dVdiat*r6yij/r6
dv(6) = dv(6) + fs * dVdiat*r6zij/r6

dv(7) = dv(7) - fs * dVdiat*r6xij/r6
dv(8) = dv(8) - fs * dVdiat*r6yij/r6
dv(9) = dv(9) - fs * dVdiat*r6zij/r6


V4 = V4 + Vdiat






call poto2q_f18(rin(1),Vdiat,dVdiat)

!            print*,'rin(1),Vdiat'
!            print*,rin(1),Vdiat



dv(7) = dv(7)  + fs * dVdiat*r1xij/r1
dv(8) = dv(8)  + fs * dVdiat*r1yij/r1
dv(9) = dv(9)  + fs * dVdiat*r1zij/r1

dv(10)= dv(10) - fs * dVdiat*r1xij/r1
dv(11)= dv(11) - fs * dVdiat*r1yij/r1
dv(12)= dv(12) - fs * dVdiat*r1zij/r1


V4 = V4 + Vdiat





V=V1*(1-fs)+(V2+V3+V4)*fs

dv=dv + dv1*(1-fs)

dv(1)  =   dv(1)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2xij / r2 + dfsdr3 * r3xij / r3)
dv(2)  =   dv(2)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2yij / r2 + dfsdr3 * r3yij / r3)
dv(3)  =   dv(3)  + (v2+v3+v4-v1) * ( dfsdr2 *  r2zij / r2 + dfsdr3 * r3zij / r3)

dv(4)  =   dv(4)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2xij / r2 + dfsdr4 * r4xij / r4)
dv(5)  =   dv(5)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2yij / r2 + dfsdr4 * r4yij / r4)
dv(6)  =   dv(6)  + (v2+v3+v4-v1) * (-dfsdr2 *  r2zij / r2 + dfsdr4 * r4zij / r4)

dv(7)  =   dv(7)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1xij / r1 - dfsdr3 * r3xij / r3)
dv(8)  =   dv(8)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1yij / r1 - dfsdr3 * r3yij / r3)
dv(9)  =   dv(9)  + (v2+v3+v4-v1) * ( dfsdr1 *  r1zij / r1 - dfsdr3 * r3zij / r3)

dv(10) =   dv(10) + (v2+v3+v4-v1) * (-dfsdr1 *  r1xij / r1 - dfsdr4 * r4xij / r4)
dv(11) =   dv(11) + (v2+v3+v4-v1) * (-dfsdr1 *  r1yij / r1 - dfsdr4 * r4yij / r4)
dv(12) =   dv(12) + (v2+v3+v4-v1) * (-dfsdr1 *  r1zij / r1 - dfsdr4 * r4zij / r4)


dv(1)  =   dv(1)  + fs * ( dv3 * r5xij / r5)
dv(2)  =   dv(2)  + fs * ( dv3 * r5yij / r5)
dv(3)  =   dv(3)  + fs * ( dv3 * r5zij / r5)

dv(4)  =   dv(4)  + fs * ( dv2 * r6xij / r6)
dv(5)  =   dv(5)  + fs * ( dv2 * r6yij / r6)
dv(6)  =   dv(6)  + fs * ( dv2 * r6zij / r6)

dv(7)  =   dv(7)  + fs * (-dv2 * r6xij / r6)
dv(8)  =   dv(8)  + fs * (-dv2 * r6yij / r6)
dv(9)  =   dv(9)  + fs * (-dv2 * r6zij / r6)

dv(10) =   dv(10) + fs * (-dv3 * r5xij / r5)
dv(11) =   dv(11) + fs * (-dv3 * r5yij / r5)
dv(12) =   dv(12) + fs * (-dv3 * r5zij / r5)


else if ((r3.ge.in_bound).and.(r4.ge.in_bound).and. &
(r5.ge.in_bound).and.(r6.ge.in_bound)) then
!             print*,'exit4'



!             print*,'O4 -> O2(1-2) + O2(3-4)'
!             print*,'tempo',tempo



if(group(i,9).eq.0)  group(i,9)=17

if(group(i,9).ne.17) then
print*,'Change in group'
group(i,9)=17
end if
! O2 function
call potoop_f2(r2,V2,dV2)
! O2 function
call potoop_f2(r1,V3,dV3)

!             print*,'r2,V2,r1,V3', r2,V2/con,r1,V3/con


! Switch function
fs=f4(rin(3),rin(4),rin(5),rin(6))
call df4(rin(3),rin(4),rin(5),rin(6),dfsdr3,dfsdr4,dfsdr5,dfsdr6)
!        print*,'dfsdr1A',dfsdr3,dfsdr4,dfsdr5,dfsdr6,fs


! for the diatomic potencials

call poto2q_f18(rin(3),Vdiat,dVdiat)


!            print*,'rin(3),Vdiat'
!            print*,rin(3),Vdiat/con

dv(1) = dv(1) + fs * dVdiat*r3xij/r3
dv(2) = dv(2) + fs * dVdiat*r3yij/r3
dv(3) = dv(3) + fs * dVdiat*r3zij/r3

dv(7) = dv(7) - fs * dVdiat*r3xij/r3
dv(8) = dv(8) - fs * dVdiat*r3yij/r3
dv(9) = dv(9) - fs * dVdiat*r3zij/r3
V4 = V4 + Vdiat


call poto2q_f18(rin(5),Vdiat,dVdiat)

!            print*,'rin(5),Vdiat'
!            print*,rin(5),Vdiat/con



dv(1) = dv(1)  + fs * dVdiat*r5xij/r5
dv(2) = dv(2)  + fs * dVdiat*r5yij/r5
dv(3) = dv(3)  + fs * dVdiat*r5zij/r5

dv(10)= dv(10) - fs * dVdiat*r5xij/r5
dv(11)= dv(11) - fs * dVdiat*r5yij/r5
dv(12)= dv(12) - fs * dVdiat*r5zij/r5


V4 = V4 + Vdiat






call poto2q_f18(rin(6),Vdiat,dVdiat)


!            print*,'rin(6),Vdiat'
!            print*,rin(6),Vdiat/con

dv(4) = dv(4) + fs * dVdiat*r6xij/r6
dv(5) = dv(5) + fs * dVdiat*r6yij/r6
dv(6) = dv(6) + fs * dVdiat*r6zij/r6

dv(7) = dv(7) - fs * dVdiat*r6xij/r6
dv(8) = dv(8) - fs * dVdiat*r6yij/r6
dv(9) = dv(9) - fs * dVdiat*r6zij/r6


V4 = V4 + Vdiat



call poto2q_f18(rin(4),Vdiat,dVdiat)

!            print*,'rin(4),Vdiat'
!            print*,rin(4),Vdiat/con



dv(4) = dv(4)  + fs * dVdiat*r4xij/r4
dv(5) = dv(5)  + fs * dVdiat*r4yij/r4
dv(6) = dv(6)  + fs * dVdiat*r4zij/r4

dv(10)= dv(10) - fs * dVdiat*r4xij/r4
dv(11)= dv(11) - fs * dVdiat*r4yij/r4
dv(12)= dv(12) - fs * dVdiat*r4zij/r4


V4 = V4 + Vdiat


!            print*,'V2',V2/con
!            print*,'V3',V3/con
!            print*,'V4',V4/con



V=V1*(1-fs)+(V2+V3+V4)*fs
!       print*,'V1,(V2+V3+V4),fs'
!       print*,V1/con,(V2+V3+V4)/con,fs

dv=dv + dv1*(1-fs)
!       print*,'dv depois 0',dv
dv(1)  =   dv(1)  + (v2+v3+v4-v1) * ( dfsdr3 *  r3xij / r3 + dfsdr5 * r5xij / r5)
!        print*
!        print*, v2,v3,v4,v1,  dfsdr3 ,  r3xij , r3 , dfsdr5 , r5xij , r5

dv(2)  =   dv(2)  + (v2+v3+v4-v1) * ( dfsdr3 *  r3yij / r3 + dfsdr5 * r5yij / r5)
dv(3)  =   dv(3)  + (v2+v3+v4-v1) * ( dfsdr3 *  r3zij / r3 + dfsdr5 * r5zij / r5)

dv(4)  =   dv(4)  + (v2+v3+v4-v1) * ( dfsdr4 * r4xij / r4 + dfsdr6 * r6xij / r6)
dv(5)  =   dv(5)  + (v2+v3+v4-v1) * ( dfsdr4 * r4yij / r4 + dfsdr6 * r6yij / r6)
dv(6)  =   dv(6)  + (v2+v3+v4-v1) * ( dfsdr4 * r4zij / r4 + dfsdr6 * r6zij / r6)

dv(7)  =   dv(7)  + (v2+v3+v4-v1) * (-dfsdr3 * r3xij / r3 - dfsdr6 * r6xij / r6)
dv(8)  =   dv(8)  + (v2+v3+v4-v1) * (-dfsdr3 * r3yij / r3 - dfsdr6 * r6yij / r6)
dv(9)  =   dv(9)  + (v2+v3+v4-v1) * (-dfsdr3 * r3zij / r3 - dfsdr6 * r6zij / r6)

dv(10) =   dv(10) + (v2+v3+v4-v1) * (-dfsdr4 *  r4xij / r4 - dfsdr5 * r5xij / r5)
dv(11) =   dv(11) + (v2+v3+v4-v1) * (-dfsdr4 *  r4yij / r4 - dfsdr5 * r5yij / r5)
dv(12) =   dv(12) + (v2+v3+v4-v1) * (-dfsdr4 *  r4zij / r4 - dfsdr5 * r5zij / r5)
!            print*,'dv depois  1',dv



dv(1)  =   dv(1)  + fs * ( dv2 * r2xij / r2)
dv(2)  =   dv(2)  + fs * ( dv2 * r2yij / r2)
dv(3)  =   dv(3)  + fs * ( dv2 * r2zij / r2)

dv(4)  =   dv(4)  + fs * (-dv2 * r2xij / r2)
dv(5)  =   dv(5)  + fs * (-dv2 * r2yij / r2)
dv(6)  =   dv(6)  + fs * (-dv2 * r2zij / r2)

dv(7)  =   dv(7)  + fs * ( dv3 * r1xij / r1)
dv(8)  =   dv(8)  + fs * ( dv3 * r1yij / r1)
dv(9)  =   dv(9)  + fs * ( dv3 * r1zij / r1)

dv(10) =   dv(10) + fs * (-dv3 * r1xij / r1)
dv(11) =   dv(11) + fs * (-dv3 * r1yij / r1)
dv(12) =   dv(12) + fs * (-dv3 * r1zij / r1)

!            print*,'dv depois  2',dv



else!  end if

v=v1

dv  = dv1

group(i,9)=0

end if

!           print*,'dv(1) ',dv(1)
!           print*,'dv(2) ',dv(2)
!           print*,'dv(3) ',dv(3)
!
!           print*,'dv(4) ',dv(4)
!           print*,'dv(5) ',dv(5)
!           print*,'dv(6) ',dv(6)
!
!           print*,'dv(7) ',dv(7)
!           print*,'dv(8) ',dv(8)
!           print*,'dv(9) ',dv(9)
!
!           print*,'dv(10)',dv(10)
!           print*,'dv(11)',dv(11)
!           print*,'dv(12)',dv(12)
!

end subroutine



!     ******************************************************************
!
!     CAIXA POTENCIAL COM:
!           FUNÇÃO VO4_11(R1,R2,R3,R4,R5,R6) CALCULA A ENERGIA POTENCIAL
!           FUNÇÃO DVO4_11(R1,R2,R3,R4,R5,R6,I) CALCULA A DERIVADA
!
!     ******************************************************************



FUNCTION VO4_11(R1in,R2in,R3in,R4in,R5in,R6in)
IMPLICIT REAL*8(A-H,O-Z)

COMMON/DM1_11/D1(4),D2(4)
COMMON/TA_11/FPE,GAMATA,SUMTA
COMMON/TA2_11/FPE2,GAMAT2,SUMT2
COMMON/PC_11/CC(8)
DATA RM,CRED,C5/2.2818,142.105,5462.358/

r1=r1in/5.2917721092d-11
r2=r2in/5.2917721092d-11
r3=r3in/5.2917721092d-11
r4=r4in/5.2917721092d-11
r5=r5in/5.2917721092d-11
r6=r6in/5.2917721092d-11

SUM=R1+R2+R3+R4+R5+R6
redK=sum-SUMTA
REDV=SUM-SUMT2
TK=FPE*EXP(-GAMATA*redK**2)
TV=FPE2*EXP(-GAMAT2*REDV**2)
SUM1=R2+R4+R5+R6
SUM2=R1+R3+R5+R6
SUM3=R1+R2+R3+R4
x1=sum1*2./CRED
x2=sum2*2./CRED
x3=sum3*2./CRED
E1=-((R1-RM)**2+(R3-RM)**2)
E2=-((R2-RM)**2+(R4-RM)**2)
E3=-((R5-RM)**2+(R6-RM)**2)
DMP1=(1.-exp(-D1(4)*X1*(1.d0+D2(4)*X1)))**5
DMP2=(1.-exp(-D1(4)*X2*(1.d0+D2(4)*X2)))**5
DMP3=(1.-exp(-D1(4)*X3*(1.d0+D2(4)*X3)))**5
SUMT=EXP(E1)*DMP1/SUM1**5+EXP(E2)*DMP2/SUM2**5 &
+EXP(E3)*DMP3/SUM3**5
VO4_11=TV+TK &
+VO2_11(R2)+VO2_11(R4)+VO2_11(R5)+VO2_11(R6) &
+VO3_11(R1,R2,R5)+VO3_11(R1,R4,R6)+VO3_11(R3,R4,R5)+ &
VO3_11(R6,R2,R3)+C5*SUMT &
+VO2_11(R1)+VO2_11(R3)


VO4_11=VO4_11*4.3597482D-18

RETURN
END

FUNCTION DVO4_11(R1in,R2in,R3in,R4in,R5in,R6in,I)
IMPLICIT double precision (A-H,O-Z)
COMMON/TA_11/FPE,GAMATA,SUMTA
COMMON/TA2_11/FPE2,GAMAT2,SUMT2

r1=r1in/5.2917721092d-11
r2=r2in/5.2917721092d-11
r3=r3in/5.2917721092d-11
r4=r4in/5.2917721092d-11
r5=r5in/5.2917721092d-11
r6=r6in/5.2917721092d-11


SUM=R1+R2+R3+R4+R5+R6
REDK=SUMTA-SUM
REDV=SUMT2-SUM
DTK=FPE*EXP(-GAMATA*REDK**2)*GAMATA*2.*REDK
DTV=FPE2*EXP(-GAMAT2*REDV**2)*GAMAT2*2.*REDV
DTA=DTK+DTV
GOTO(1,2,3,4,5,6)I
1 DVO4_11=DVO2_11(R1)+DVO3_11(R1,R2,R5,1)+DVO3_11(R1,R4,R6,1)+DTA+ &
DERED_11(R1,R2,R3,R4,R5,R6)!*4.3597482D-18/5.2917721092d-11
GOTO 7
2 DVO4_11=DVO2_11(R2)+DVO3_11(R1,R2,R5,2)+DVO3_11(R6,R2,R3,2)+DTA+ &
DERED_11(R2,R1,R4,R3,R5,R6)!*4.3597482D-18/5.2917721092d-11
GOTO 7
3 DVO4_11=DVO2_11(R3)+DVO3_11(R3,R4,R5,1)+DVO3_11(R6,R2,R3,3)+DTA+ &
DERED_11(R3,R2,R1,R4,R5,R6)!*4.3597482D-18/5.2917721092d-11
GOTO 7
4 DVO4_11=DVO2_11(R4)+DVO3_11(R1,R4,R6,2)+DVO3_11(R3,R4,R5,2)+DTA+ &
DERED_11(R4,R1,R2,R3,R5,R6)!*4.3597482D-18/5.2917721092d-11
GOTO 7
5 DVO4_11=DVO2_11(R5)+DVO3_11(R1,R2,R5,3)+DVO3_11(R3,R4,R5,3)+DTA+ &
DERED_11(R5,R2,R6,R4,R1,R3)!*4.3597482D-18/5.2917721092d-11
GOTO 7
6 DVO4_11=DVO2_11(R6)+DVO3_11(R1,R4,R6,3)+DVO3_11(R6,R2,R3,1)+DTA+ &
DERED_11(R6,R2,R5,R4,R1,R3)!*4.3597482D-18/5.2917721092d-11
7 DVO4_11=DVO4_11*4.3597482D-18/5.2917721092d-11
RETURN


END

FUNCTION VO2_11(R)
IMPLICIT REAL*8(A-H,O-Z)
COMMON/MAIN_11/CC(3),A(3),B(3),DD
COMMON /V2MAS_11/VTOT,DVTOT,D2VTOT,RCSN
RCSN=R
CALL V2_11
ACE=0.0
DO 1 I=1,3
ACE=ACE+CC(I)*DMR_11(I,R)
1 CONTINUE
VO2_11=ACE+VTOT
RETURN
END

FUNCTION DVO2_11(R)
IMPLICIT REAL*8(A-H,O-Z)
COMMON/MAIN_11/CC(3),A(3),B(3),DD
COMMON /V2MAS_11/VTOT,DVTOT,D2VTOT,RCSN
RCSN=R
CALL DV2_11
DACE=0.0
DO 1 I=1,3
DACE=DACE+CC(I)*DDMR_11(I,R)
1 CONTINUE
DVO2_11=DACE+DVTOT
RETURN
END

FUNCTION VO3_11(R1,R2,R3)
IMPLICIT REAL*8(A-H,O-Z)
DIMENSION R(3)
R(1)=R1
R(2)=R2
R(3)=R3
CALL PRINC_11(R,YR)
aaa=VAD_11(r1,r2,r3)
bbb=PTH_11(r1,r2,r3)
VO3_11=YR+aaa+bbb
RETURN
END

FUNCTION DVO3_11(R1,R2,R3,I)
IMPLICIT REAL*8(A-H,O-Z)
DIMENSION R(3)
R(1)=R1
R(2)=R2
R(3)=R3
CALL DPRINC_11(R,DY,I)
aaa=dVAD_11(r1,r2,r3,i)
bbb=dPTH_11(r1,r2,r3,i)
DVO3_11=DY+aaa+bbb
RETURN
END


FUNCTION DVO4N_11(R1,R2,R3,R4,R5,R6,I)
IMPLICIT REAL*8(A-H,O-Z)
DIMENSION R(6)
R(1)=R1
R(2)=R2
R(3)=R3
R(4)=R4
R(5)=R5
R(6)=R6
R(I)=R(I)+1.E-7
VSUP=VO4_11(R(1),R(2),R(3),R(4),R(5),R(6))
R(I)=R(I)-2.E-7
VINF=VO4_11(R(1),R(2),R(3),R(4),R(5),R(6))
DVO4N_11=(VSUP-VINF)/2.E-7
RETURN
END


FUNCTION DERED_11(R1,R2,R3,R4,R5,R6)
IMPLICIT REAL*8(A-H,O-Z)
COMMON/DM1_11/D1(4),D2(4)
DATA RM,CRED,C5/2.2818,142.105,5462.358/
SUM1=R2+R4+R5+R6
SUM2=R1+R3+R5+R6
SUM3=R1+R2+R3+R4
x1=sum1*2./CRED
x2=sum2*2./CRED
x3=sum3*2./CRED
E1=-((R1-RM)**2+(R3-RM)**2)
E2=-((R2-RM)**2+(R4-RM)**2)
E3=-((R5-RM)**2+(R6-RM)**2)
DMP1=(1.-exp(-D1(4)*X1*(1.d0+D2(4)*X1)))**5
DX2=2.D0/CRED
DX3=DX2
POL2=-D1(4)*X2*(1.D0+D2(4)*X2)
POL3=-D1(4)*X3*(1.D0+D2(4)*X3)
DPOL2=-D1(4)*(1.D0+2.D0*D2(4)*X2)
DPOL3=-D1(4)*(1.D0+2.D0*D2(4)*X3)
EXPO2=EXP(POL2)
EXPO3=EXP(POL3)
RL2=SUM2**5
RL3=SUM3**5
C2=(1.-EXPO2)**5
C3=(1.-EXPO3)**5
EXPOEL2=C2/(1.-EXPO2)
EXPOEL3=C3/(1.-EXPO3)
DC2=5.D0*EXPOEL2*(-EXPO2)*DPOL2*DX2
DC3=5.D0*EXPOEL3*(-EXPO3)*DPOL3*DX3
DDMR2=DC2/RL2-5.D0*C2/(RL2*SUM2)
DDMR3=DC3/RL3-5.D0*C3/(RL3*SUM3)
DERED_11=C5*(2./SUM1**5*EXP(E1)*DMP1*(RM-R1)+ &
DDMR2*EXP(E2)+DDMR3*EXP(E3))
RETURN
END



FUNCTION DDVO4_11(R,ID,JD)
IMPLICIT REAL*8(A-H,O-Z)
DIMENSION R(6)
R(JD)=R(JD)+1.E-5
VSUP=DVO4_11(R(1),R(2),R(3),R(4),R(5),R(6),ID)
R(JD)=R(JD)-2.E-5
VINF=DVO4_11(R(1),R(2),R(3),R(4),R(5),R(6),ID)
R(JD)=R(JD)+1.E-5
DDVO4_11=(VSUP-VINF)/2.E-5
RETURN
END

SUBROUTINE PRINC_11(RR,Y)
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/BFF_11/VHF,ACE
COMMON/COMN_11/VHF2,ACE2
COMMON/PART2_11/VEQ,DVEQ,D2VEQ
COMMON/MAIN_11/C2C(3),A0(3),A1(3),DD
COMMON/COORD_11/R(3)
COMMON/PART_11/VDP,DVDP,D2VDP
COMMON/V2MAS_11/VTOT,DVTOT,D2VTOT,RCSN
COMMON/MASMAIN_11/VR2,DVR2,D2VR2
COMMON/CONTROL_11/K1
COMMON/VVRR_11/ZP(3)
DIMENSION RR(3)
DO 5 I=1,3
R(I)=RR(I)
5 CONTINUE
CALL CE_11
CALL CEQ_11
Y=VDP+VEQ
RETURN
END

SUBROUTINE DPRINC_11(RR,DY,K1D)
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/BFF_11/VHF,ACE
COMMON/COMN_11/VHF2,ACE2
COMMON/PART2_11/VEQ,DVEQ,D2VEQ
COMMON/MAIN_11/C2C(3),A0(3),A1(3),DD
COMMON/COORD_11/R(3)
COMMON/PART_11/VDP,DVDP,D2VDP
COMMON/V2MAS_11/VTOT,DVTOT,D2VTOT,RCSN
COMMON/MASMAIN_11/VR2,DVR2,D2VR2
COMMON/CONTROL_11/K1
COMMON/VVRR_11/ZP(3)
DIMENSION RR(3)
K1=K1D
DO 5 I=1,3
R(I)=RR(I)
5 CONTINUE
CALL DCE_11
RCSN=R(K1D)
CALL DCEQ_11
DY=DVDP+DVEQ
RETURN
END

SUBROUTINE V2_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/BFF_11/VHF,ACE
COMMON/START_11/A,B,AK0,AK1,RZ,CXX
COMMON/V2MAS_11/VTOT,DVTOT,D2VTOT,R
COMMON/MAIN_11/CC(3),A00(3),A11(3),DD
COMMON/DM1_11/D1(4),D2(4)
COMMON/DM2_11/R0,RM
COMMON/DM3_11/GAMA
COMMON/DM4_11/A1,A2,A3,GM,RE,DE
RR=R-RE
EGM=EXP(-GM*RR)
POL=DE*(1.+A1*RR+A2*RR**2+A3*RR**3)
VTOT=POL*EGM
RETURN
END

SUBROUTINE DV2_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/START_11/A,B,AK0,AK1,RZ,CXX
COMMON/V2MAS_11/VTOT,DVTOT,D2VTOT,R
COMMON/MAIN_11/CC(3),A00(3),A11(3),DD
COMMON/DM1_11/D1(4),D2(4)
COMMON/DM2_11/R0,RM
COMMON/DM3_11/GAMA
COMMON/DM4_11/A1,A2,A3,GM,RE,DE
RR=R-RE
EGM=EXP(-GM*RR)
POL=DE*(1.+A1*RR+A2*RR**2+A3*RR**3)
DPOL=DE*(A1+2.*A2*RR+3.*A3*RR**2)
DVTOT=(DPOL*EGM-GM*POL*EGM)
RETURN
END

SUBROUTINE CE_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/CECF_11/CC,DCC,D2CC
COMMON/MAIN_11/C2C(3),A0(3),A1(3),DD
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/PART_11/VDP,DVDP,D2VDP
COMMON/DUMMY_11/INDEX(3,2)
COMMON/COORD_11/R(3)
VDP=0.0
DO 10 J=1,3
I1=INDEX(J,1)
I2=INDEX(J,2)
A=R(I1)
B=R(I2)
DO 10 I=1,3
CXX=C2C(I)
AK0=A0(I)
AK1=A1(I)
CALL CF_11
VDP=VDP+(CC-CXX)*DMR_11(I,R(J))
10 CONTINUE
RETURN
END

SUBROUTINE DCE_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/CECF_11/CC,DCC,D2CC
COMMON/MAIN_11/C2C(3),A0(3),A1(3),DD
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/CONTROL_11/K1D
COMMON/PART_11/VDP,DVDP,D2VDP
COMMON/DUMMY_11/INDEX(3,2)
COMMON/DMY_11/IX(3,3)
COMMON/COORD_11/R(3)
DVDP=0.0
DO 10 I=1,3
J=K1D
I1=INDEX(J,1)
I2=INDEX(J,2)
A=R(I1)
B=R(I2)
CXX=C2C(I)
AK0=A0(I)
AK1=A1(I)
CALL CF_11
DVDP=DVDP+(CC-CXX)*DDMR_11(I,R(J))
DO 10 K=1,2
J=INDEX(K1D,K)
IO=IX(K1D,J)
A=R(K1D)
B=R(IO)
CALL DCF_11
DVDP=DVDP+DCC*DMR_11(I,R(J))
10 CONTINUE
RETURN
END

SUBROUTINE CF_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/CECF_11/CC,DCC,D2CC
TANHA=TANH(AK1*A)
TANHB=TANH(AK1*B)
EA=AK0*EXP(-AK1*(A-R0))
EB=AK0*EXP(-AK1*(B-R0))
CC=0.5*CXX*((1.+EA)*TANHB+(1.+EB)*TANHA)
RETURN
END

SUBROUTINE DCF_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/CECF_11/CC,DCC,D2CC
EA=AK0*EXP(-AK1*(A-R0))
EB=AK0*EXP(-AK1*(B-R0))
AKA=AK1*A
TANHA=TANH(AKA)
TANHB=TANH(AK1*B)
SINHA=SINH(AKA)
COSHA=SINHA/TANHA
CSHA2=1./COSHA*1./COSHA
DCC=0.5*CXX*AK1*(-EA*TANHB+(1.+EB)*CSHA2)
RETURN
END

FUNCTION DMR_11(I,R)
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/PN_11/N(4)
COMMON/DM1_11/D1(4),D2(4)
COMMON/DM2_11/R0,RM
COMMON/DM3_11/GAMA
RO=(RM+GAMA*R0)/2.
X=R/RO
DX=1./RO
POL=-D1(I)*X*(1.+D2(I)*X)
EXPO=EXP(POL)
IDXL=N(I)
DXL=FLOAT(IDXL)
C=(1.-EXPO)**IDXL
RL=R**IDXL
DMR_11=-C/RL
RETURN
END

FUNCTION DDMR_11(I,R)
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/PN_11/N(4)
COMMON/DM1_11/D1(4),D2(4)
COMMON/DM2_11/R0,RM
COMMON/DM3_11/GAMA
RO=(RM+GAMA*R0)/2.
X=R/RO
DX=1./RO
POL=-D1(I)*X*(1.+D2(I)*X)
DPOL=-D1(I)*(1.+2.*D2(I)*X)
EXPO=EXP(POL)
IDXL=N(I)
DXL=FLOAT(IDXL)
C=(1.-EXPO)**IDXL
EXPOEL=C/(1.-EXPO)
DC=DXL*EXPOEL*(-EXPO)*DPOL*DX
RL=R**IDXL
DDMR_11=-(DC/RL-DXL*C/(RL*R))
RETURN
END

SUBROUTINE CEQ_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/MAIN2_11/C5Q,RK
COMMON/START2_11/C
COMMON/PART2_11/VEQ,DVEQ,D2VEQ
COMMON/CECF_11/CC,DCC,D2CC
COMMON/CECF2_11/DCCM
COMMON/MAIN_11/C2C(3),A0(3),A1(3),DD
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/CONTROL_11/K1D
COMMON/PART_11/VDP,DVDP,D2VDP
COMMON/DUMMY_11/INDEX(3,2)
COMMON/COORD_11/R(3)
VEQ=0.0
DO 10 J=1,3
I1=INDEX(J,1)
I2=INDEX(J,2)
A=R(I1)
B=R(I2)
C=R(J)
CALL CFQ_11
VEQ=VEQ+CC*DMR_11(4,R(J))
10 CONTINUE
RETURN
END

SUBROUTINE DCEQ_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/MAIN2_11/C5Q,RK
COMMON/START2_11/C
COMMON/PART2_11/VEQ,DVEQ,D2VEQ
COMMON/CECF_11/CC,DCC,D2CC
COMMON/CECF2_11/DCCM
COMMON/MAIN_11/C2C(3),A0(3),A1(3),DD
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/CONTROL_11/K1D
COMMON/PART_11/VDP,DVDP,D2VDP
COMMON/DUMMY_11/INDEX(3,2)
COMMON/DMY_11/IX(3,3)
COMMON/COORD_11/R(3)
DVEQ=0.0
J=K1D
I1=INDEX(J,1)
I2=INDEX(J,2)
A=R(I1)
B=R(I2)
C=R(J)
CALL CFQ_11
DVEQ=DVEQ+CC*DDMR_11(4,R(J))
DO 10 K=1,2
J=INDEX(K1D,K)
IO=IX(K1D,J)
A=R(K1D)
B=R(IO)
C=R(J)
CALL DCFQ_11
DVEQ=DVEQ+DCC*DMR_11(4,R(J))
10 CONTINUE
RETURN
END

SUBROUTINE CFQ_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/MAIN2_11/C5Q,RK
COMMON/START2_11/C
COMMON/PART2_11/VEQ,DVEQ,D2VEQ
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/CECF_11/CC,DCC,D2CC
COMMON/CECF2_11/DCCM
CCA=A**4/R0**4
CCB=B**4/R0**4
TANHA=TANH(RK*A)
TANHB=TANH(RK*B)
EA=EXP(-RK*(A-R0))
EB=EXP(-RK*(B-R0))
CC=0.5*C5Q*(CCA*EA*TANHB+CCB*EB*TANHA)
RETURN
END

SUBROUTINE DCFQ_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/START2_11/C
COMMON/MAIN2_11/C5Q,RK
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/CECF_11/CC,DCC,D2CC
COMMON/CECF2_11/DCCM
RKA=RK*A
CCA=A**4/R0**4
DCCA=4.*A**3/R0**4
CCB=B**4/R0**4
EA=EXP(-RK*(A-R0))
EB=EXP(-RK*(B-R0))
TANHA=TANH(RKA)
TANHB=TANH(RK*B)
if(rka.lt.20.0)then
SINHA=SINH(RKA)
COSHA=SINHA/TANHA
CSHA2=1./COSHA*1./COSHA
else
csha2=0.0
endif
DCC=0.5*C5Q*(-CCA*EA*TANHB*RK+DCCA*EA*TANHB+RK*CCB*EB*CSHA2)
RETURN
END

SUBROUTINE DCFQM_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/MAIN2_11/C5Q,RK
COMMON/START2_11/C
COMMON/PART2_11/VEQ,DVEQ,D2VEQ
COMMON/START_11/A,B,AK0,AK1,R0,CXX
COMMON/CECF_11/CC,DCC,D2CC
COMMON/CECF2_11/DCCM
RRA=(B-C)/A
RRB=(A-C)/B
DCCA=-6.0*RRA*A**3/R0**4
DCCB=-6.0*RRB*B**3/R0**4
TANHA=TANH(RK*A)
TANHB=TANH(RK*B)
EA=EXP(-RK*(A-R0))
EB=EXP(-RK*(B-R0))
DCCM=0.5*C5Q*(DCCA*EA*TANHB+DCCB*EB*TANHA)
RETURN
END

FUNCTION VAD_11(R1,R2,R3)
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/DM8_11/DRDQ(3,3),DQDR(3,3)
COMMON/VAR_11/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
COMMON/MAIN_11/V(3),AK(3),BK(3),DD
COMMON/PASS_11/Q1,Q2,Q3,TH,DTH
S1=R1-DD
S2=R2-DD
S3=R3-DD
Q1=S1*DQDR(1,1)+S2*DQDR(1,2)+S3*DQDR(1,3)
Q2=S1*DQDR(2,1)+S2*DQDR(2,2)+S3*DQDR(2,3)
Q3=S1*DQDR(3,1)+S2*DQDR(3,2)+S3*DQDR(3,3)
PART=SUMAT_11(R1,R2,R3)
TH=1.-TANH(GAMA*Q1/4.0)
VAD_11=(B1+B2*Q1+(AL(1)+AL(2)*(Q2**2+Q3**2)+AL(3)* &
(Q2**2+Q3**2)**2+ AL4*(Q2**2+Q3**2))*(1.5+PART)**5)* &
(1.2527+PART)**5*EXP(-BETAD*(Q2**2+Q3**2))*TH
RETURN
END

FUNCTION SUMAT_11(R1,R2,R3)
IMPLICIT REAL*8 (A-H,O-Z)
SUMAT_11=(R1**2-R2**2-R3**2)/(2.*R2*R3)+ &
(R2**2-R3**2-R1**2)/(2.*R1*R3)+ &
(R3**2-R1**2-R2**2)/(2.*R1*R2)
RETURN
END

FUNCTION PTH_11(R1,R2,R3)
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/VAR_11/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
COMMON/CONST2_11/A,C3,C6,C10
COMMON/PASS_11/Q1,Q2,Q3,TH,DTH


TH2=1.-TANH(GAMA*Q1/2.)
P=A+C1*Q1+C2*Q1**2+C3*(Q2**2+Q3**2)+C5*Q1*(Q2**2+Q3**2) &
+C6*(Q3**3-3.*Q2**2*Q3)+C10*(Q2**2+Q3**2)**2
PTH_11=P*TH2
RETURN
END

FUNCTION DVAD_11(R1,R2,R3,I)
IMPLICIT REAL*8(A-H,O-Z)
COMMON/DM8_11/DRDQ(3,3),DQDR(3,3)
COMMON/VAR_11/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
COMMON/MAIN_11/V(3),AK(3),BK(3),DD
COMMON/PASS_11/Q1,Q2,Q3,TH,DTH
S1=R1-DD
S2=R2-DD
S3=R3-DD
Q1=S1*DQDR(1,1)+S2*DQDR(1,2)+S3*DQDR(1,3)
Q2=S1*DQDR(2,1)+S2*DQDR(2,2)+S3*DQDR(2,3)
Q3=S1*DQDR(3,1)+S2*DQDR(3,2)+S3*DQDR(3,3)
TH=1.-TANH(GAMA*Q1/4.0)
DTH=-GAMA/(4.0*COSH(GAMA*Q1/4.0)**2)*DQDR(1,I)
PART=SUMAT_11(R1,R2,R3)
DPART=DSUMAT_11(R1,R2,R3,I)
T1=B1+B2*Q1
T2=(AL(1)+AL(2)*(Q2**2+Q3**2)+AL(3)*(Q2**2+Q3**2)**2)
T3=(1.5+PART)**5
T4=(1.2527+PART)**5
T5=EXP(-BETAD*(Q2**2+Q3**2))
T6=TH
DT1=B2*DQDR(1,I)
DT2=2.*AL(2)*(Q2*DQDR(2,I)+Q3*DQDR(3,I)) &
+4.*AL(3)*((Q2**3+Q2*Q3**2)*DQDR(2,I)+(Q3**3+Q2**2*Q3) &
*DQDR(3,I))
DT3=5.*(1.5000+PART)**4*DPART
DT4=5.*(1.2527+PART)**4*DPART
DT5=-2.*BETAD*(Q2*DQDR(2,I)+Q3*DQDR(3,I))*T5
DT6=DTH
DVAD_11=(DT1+DT2*T3+T2*DT3)*T4*T5*T6 &
+(T1+T2*T3)*DT4*T5*T6 &
+(T1+T2*T3)*T4*DT5*T6 &
+(T1+T2*T3)*T4*T5*DT6
RETURN
END

FUNCTION DSUMAT_11( R1,R2,R3,I)
IMPLICIT REAL*8(A-H,O-Z)
COMMON/DUMMY_11/INDEX(3,2)
DIMENSION R(3)
R(1)=R1
R(2)=R2
R(3)=R3
I1=INDEX(I,1)
I2=INDEX(I,2)
DSUMAT_11=R(I)/(R(I1)*R(I2)) &
+(R(I2)**3-R(I1)**2*R(I2)-R(I)**2*R(I2))/(2.*R(I)**2*R(I2)**2) &
+(R(I1)**3-R(I2)**2*R(I1)-R(I)**2*R(I1))/(2.*R(I)**2*R(I1)**2)
RETURN
END

FUNCTION DPTH_11(R1,R2,R3,I)
IMPLICIT REAL*8(A-H,O-Z)
COMMON/DM8_11/DRDQ(3,3),DQDR(3,3)
COMMON/VAR_11/C1,C2,C5,GAMA,B1,B2,AL(3),AL4,BETAD
COMMON/CONST2_11/A,C3,C6,C10
COMMON/PASS_11/Q1,Q2,Q3,TH,DTH


DTH2=-GAMA/(2.*COSH(GAMA*Q1/2.)**2)*DQDR(1,I)
TH2=1.-TANH(GAMA*Q1/2.)
P=A+C1*Q1+C2*Q1**2+C3*(Q2**2+Q3**2)+C5*Q1*(Q2**2+Q3**2)+C6*(Q3**3- &
3.*Q2**2*Q3)+C10*(Q2**2+Q3**2)**2
DPVTH=DTH2*P
DP=C1*DQDR(1,I)+2.*C2*Q1*DQDR(1,I)+2.*C3*(Q2*DQDR(2,I) &
+Q3*DQDR(3,I))+C5*DQDR(1,I)*(Q2**2+Q3**2)+C5*Q1*(2.*Q2 &
*DQDR(2,I)+2.*Q3*DQDR(3,I))+C6*(3.*Q3**2*DQDR(3,I)-6.*Q2*Q3 &
*DQDR(2,I)-3.*Q2**2*DQDR(3,I))+C10*((4.*Q2**3+4.*Q2*Q3**2) &
*DQDR(2,I)+(4.*Q3**3+4.*Q2**2*Q3)*DQDR(3,I))
DPTH_11=DPVTH+DP*TH2
RETURN
END

BLOCK DATA AAA1_11
IMPLICIT REAL*8 (A-H,O-Z)
COMMON/MAIN_11/VW(3),AK0(3),AK1(3),DD
COMMON/MAIN2_11/C5Q,RK
COMMON/PN_11/N(4)
COMMON/DM1_11/D1(4),D2(4)
COMMON/DUMMY_11/INDEX(3,2)
COMMON/DMY_11/IX(3,3)
COMMON/DM2_11/R0,RM
COMMON/DM3_11/GAMA
COMMON/DM4_11/A1,A2,A3,GM,RE,DE
COMMON/DM8_11/DRDQ(3,3),DQDR(3,3)
COMMON/VAR_11/C1,C2,C5,GAMAP,B1,B2,AL(3),AL4,BETAD
COMMON/CONST2_11/A,C3,C6,C10
COMMON/START_11/AS,BS,ASK0,ASK1,R0PASS,CSXX
COMMON/TA2_11/FPE2,GAMAT2,SUMT2
COMMON/TA_11/FPE,GAMATA,SUMTA
DATA FPE,GAMATA,SUMTA/2.9E-2,3.1410047d-2,22.323365/
DATA FPE2,GAMAT2,SUMT2/1.4301464d-3,3.1410047d-2,35.549516/
DATA R0PASS/2.2818/
DATA C1,C2,C5,GAMAP,B1,B2,(AL(I),I=1,3),AL4,BETAD/.1157783, &
.1682194, &
-.6563872e-02,2.54,.1561575E+03,.4951110E+02,.11574778E7, &
-.35624264E6,+.27725091E5,0.0,.224/
DATA A,C3,C6,C10/.3058714,-.1502505,.1805864E-1,.3434551E-1/
DATA(VW(I),I=1,3)/15.4,235.21994,4066.2393/
DATA(AK0(I),I=1,3)/.2772806,.6233148,1.237123/
DATA(AK1(I),I=1,3)/.9545174,.8116684,.7043832/
DATA DD/2.967200/
DATA RK,C5Q/3.3522,+1.3144/
DATA(N(I),I=1,4)/6,8,10,5/
DATA ((IX(I,J),I=1,3),J=1,3)/0,3,2,3,0,1,2,1,0/
DATA ((INDEX(I,J),I=1,3),J=1,2)/2,1,1,3,3,2/
DATA (D1(I),I=1,4)/3.0951333,2.1999,1.6880714,3.8428294/
DATA (D2(I),I=1,4)/2.8363203,3.2849276,3.5239687,2.5178884/
DATA GAMA/2.5/
DATA R0,RM/5.661693,2.2818/
DATA A1,A2,A3,GM,RE,DE/3.6445902,3.9281208,2.0986665,3.3522487, &
2.2818,-0.142912/
DATA ((DRDQ(I,J),I=1,3),J=1,3)/.5773502,.5773502,.5773502, &
0.0,.7071067,-.7071067, &
.8164965,-.4082482,-.4082482/
DATA ((DQDR(I,J),I=1,3),J=1,3)/.5773502,0.0,.8164965, &
.5773502,.7071067,-.4082482, &
.5773502,-.7071067,-.4082482/
END


!     ******************************************************************
!
!     FIM DE CAIXA POTENCIAL COM
!
!     ******************************************************************



