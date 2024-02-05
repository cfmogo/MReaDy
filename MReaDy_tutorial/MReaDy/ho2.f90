        subroutine vho2_f4_switch(i,V,dV)
        use variables
        use constants
        use phys_parameters 
        implicit none
        double precision, dimension(12) :: dV1,dV
        integer      :: gi3,gi4,gi5
        integer      :: i
        double precision :: rxij,ryij,rzij,rij
        double precision :: r1xij,r1yij,r1zij
        double precision :: r2xij,r2yij,r2zij
        double precision :: r3xij,r3yij,r3zij
        double precision :: r1,r2,r3
        double precision :: V, V1,V2,dV2,V4,vdiat,dvdiat        
        double precision :: f2,fs,dfsdr1,dfsdr2,dfsdr3
        double precision :: vho2_f4
        double precision :: dv1dr1,dv1dr2,dv1dr3

!      print*,'whats up?1'
      dV1=0.0d0
      dV=0.0d0
      rxij=0.0d0
      ryij=0.0d0
      rzij=0.0d0
      rij=0.0d0
      r1xij=0.0d0
      r1yij=0.0d0
      r1zij=0.0d0
      r2xij=0.0d0
      r2yij=0.0d0
      r2zij=0.0d0
      r3xij=0.0d0
      r3yij=0.0d0
      r3zij=0.0d0
      r1 =0.0d0
      r2=0.0d0
      r3=0.0d0
      V =0.0d0
      V1=0.0d0
      V2=0.0d0
      dV2=0.0d0
      V4=0.0d0
      vdiat=0.0d0
      dvdiat=0.0d0
      dv1dr1=0.0d0
      dv1dr2=0.0d0
      dv1dr3=0.0d0













 

            gi3=group(i,3)      ! H atom
            gi4=group(i,4)      ! O atom
            gi5=group(i,5)      ! O atom
!            print*,gi3,gi4,gi5,'ho2 1'
call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi4),ry(gi4),rz(gi4),rxij,ryij,rzij,rij)

            r1=rij
            r1xij=rxij            
            r1yij=ryij
            r1zij=rzij

call mod_dist(rx(gi3),ry(gi3),rz(gi3),rx(gi5),ry(gi5),rz(gi5),rxij,ryij,rzij,rij)


            r2=rij
            r2xij=rxij            
            r2yij=ryij
            r2zij=rzij

call mod_dist(rx(gi4),ry(gi4),rz(gi4),rx(gi5),ry(gi5),rz(gi5),rxij,ryij,rzij,rij)

            r3=rij
            r3xij=rxij            
            r3yij=ryij
            r3zij=rzij
!            print*,r1,r2,r3,'ho2 2'
      



            V1=vho2_f4(r3,r2,r1)
         
            
            call dervho2_f4(r3,r2,r1,dv1dr3,dv1dr2,dv1dr1)

 

           
           dv1(1)  =   + dv1dr1 * r1xij / r1 + dv1dr2 * r2xij / r2
           dv1(2)  =   + dv1dr1 * r1yij / r1 + dv1dr2 * r2yij / r2
           dv1(3)  =   + dv1dr1 * r1zij / r1 + dv1dr2 * r2zij / r2
                      
           dv1(4)  =   - dv1dr1 * r1xij / r1 + dv1dr3 * r3xij / r3
           dv1(5)  =   - dv1dr1 * r1yij / r1 + dv1dr3 * r3yij / r3
           dv1(6)  =   - dv1dr1 * r1zij / r1 + dv1dr3 * r3zij / r3
                      
           dv1(7)  =   - dv1dr2 * r2xij / r2 - dv1dr3 * r3xij / r3
           dv1(8)  =   - dv1dr2 * r2yij / r2 - dv1dr3 * r3yij / r3
           dv1(9)  =   - dv1dr2 * r2zij / r2 - dv1dr3 * r3zij / r3
 
            if ((r1.ge.in_bound).and.(r2.ge.in_bound).and. &        
        (r3.le.in_bound)) then             
!            print*,'118'


        ! O2 potential                      
                    
        call potoop_f2(r3,V2,dV2)        
       ! Switch function                    
        fs=f2(r1,r2)  
        call df2(r1,r2,dfsdr1,dfsdr2)
       ! for the diatomic potencials

            call pothoq_f17(r1,Vdiat,dVdiat)
 


            dv(1) = dv(1) + fs * dVdiat*r1xij/r1 
            dv(2) = dv(2) + fs * dVdiat*r1yij/r1 
            dv(3) = dv(3) + fs * dVdiat*r1zij/r1 
                          
            dv(4) = dv(4) - fs * dVdiat*r1xij/r1 
            dv(5) = dv(5) - fs * dVdiat*r1yij/r1 
            dv(6) = dv(6) - fs * dVdiat*r1zij/r1 


            V4 = V4 + Vdiat


 
            call pothoq_f17(r2,Vdiat,dVdiat)


            dv(1) = dv(1) + fs * dVdiat*r2xij/r2 
            dv(2) = dv(2) + fs * dVdiat*r2yij/r2 
            dv(3) = dv(3) + fs * dVdiat*r2zij/r2 
                                    
            dv(7) = dv(7) - fs * dVdiat*r2xij/r2 
            dv(8) = dv(8) - fs * dVdiat*r2yij/r2 
            dv(9) = dv(9) - fs * dVdiat*r2zij/r2 


            V4 = V4 + Vdiat



        V=V1*(1-fs)+(V2+V4)*fs

       dv=dv + dv1*(1-fs) 
              
       dv(1)  =   dv(1)  + (v2+v4-v1) * ( dfsdr1 *  r1xij / r1 + dfsdr2 * r2xij / r2)
       dv(2)  =   dv(2)  + (v2+v4-v1) * ( dfsdr1 *  r1yij / r1 + dfsdr2 * r2yij / r2)
       dv(3)  =   dv(3)  + (v2+v4-v1) * ( dfsdr1 *  r1zij / r1 + dfsdr2 * r2zij / r2)
                              
       dv(4)  =   dv(4)  + (v2+v4-v1) * (-dfsdr1 *  r1xij / r1)
       dv(5)  =   dv(5)  + (v2+v4-v1) * (-dfsdr1 *  r1yij / r1)
       dv(6)  =   dv(6)  + (v2+v4-v1) * (-dfsdr1 *  r1zij / r1)
                              
       dv(7)  =   dv(7)  + (v2+v4-v1) * (-dfsdr2 *  r2xij / r2)
       dv(8)  =   dv(8)  + (v2+v4-v1) * (-dfsdr2 *  r2yij / r2)
       dv(9)  =   dv(9)  + (v2+v4-v1) * (-dfsdr2 *  r2zij / r2)
                              
                              
       dv(4)  =   dv(4)  + fs * ( dv2 * r3xij / r3)
       dv(5)  =   dv(5)  + fs * ( dv2 * r3yij / r3)
       dv(6)  =   dv(6)  + fs * ( dv2 * r3zij / r3)
                               
       dv(7)  =   dv(7)  + fs * (-dv2 * r3xij / r3)
       dv(8)  =   dv(8)  + fs * (-dv2 * r3yij / r3)
       dv(9)  =   dv(9)  + fs * (-dv2 * r3zij / r3)
                                



       else if ((r1.ge.in_bound).and.(r3.ge.in_bound).and. &        
        (r2.le.in_bound)) then             
        !print*,'208'


        ! O2 potential                      
                    
        call potohp_f3(r2,V2,dV2)        
 

       ! Switch function                    
        fs=f2(r1,r3)  
        call df2(r1,r3,dfsdr1,dfsdr3)
 
        ! for the diatomic potencials

            call pothoq_f17(r1,Vdiat,dVdiat)



            dv(1) = dv(1) + fs * dVdiat*r1xij/r1 
            dv(2) = dv(2) + fs * dVdiat*r1yij/r1 
            dv(3) = dv(3) + fs * dVdiat*r1zij/r1 
                          
            dv(4) = dv(4) - fs * dVdiat*r1xij/r1 
            dv(5) = dv(5) - fs * dVdiat*r1yij/r1 
            dv(6) = dv(6) - fs * dVdiat*r1zij/r1 


            V4 = V4 + Vdiat


 
            call poto2q_f18(r3,Vdiat,dVdiat)


            dv(4) = dv(4) + fs * dVdiat*r2xij/r2 
            dv(5) = dv(5) + fs * dVdiat*r2yij/r2 
            dv(6) = dv(6) + fs * dVdiat*r2zij/r2 
                                    
            dv(7) = dv(7) - fs * dVdiat*r2xij/r2 
            dv(8) = dv(8) - fs * dVdiat*r2yij/r2 
            dv(9) = dv(9) - fs * dVdiat*r2zij/r2 


            V4 = V4 + Vdiat



        V=V1*(1-fs)+(V2+V4)*fs

       dv=dv + dv1*(1-fs) 
              
       dv(1)  =   dv(1)  + (v2+v4-v1) * ( dfsdr1 *  r1xij / r1)
       dv(2)  =   dv(2)  + (v2+v4-v1) * ( dfsdr1 *  r1yij / r1)
       dv(3)  =   dv(3)  + (v2+v4-v1) * ( dfsdr1 *  r1zij / r1)
                              
       dv(4)  =   dv(4)  + (v2+v4-v1) * (-dfsdr1 *  r1xij / r1 + dfsdr3 *  r3xij / r3)
       dv(5)  =   dv(5)  + (v2+v4-v1) * (-dfsdr1 *  r1yij / r1 + dfsdr3 *  r3yij / r3)
       dv(6)  =   dv(6)  + (v2+v4-v1) * (-dfsdr1 *  r1zij / r1 + dfsdr3 *  r3zij / r3)
                              
       dv(7)  =   dv(7)  + (v2+v4-v1) * (-dfsdr3 *  r3xij / r3)
       dv(8)  =   dv(8)  + (v2+v4-v1) * (-dfsdr3 *  r3yij / r3)
       dv(9)  =   dv(9)  + (v2+v4-v1) * (-dfsdr3 *  r3zij / r3)
                              
                              
       dv(1)  =   dv(1)  + fs * ( dv2 * r2xij / r2)
       dv(2)  =   dv(2)  + fs * ( dv2 * r2yij / r2)
       dv(3)  =   dv(3)  + fs * ( dv2 * r2zij / r2)
                              
       dv(7)  =   dv(7)  + fs * (-dv2 * r2xij / r2)
       dv(8)  =   dv(8)  + fs * (-dv2 * r2yij / r2)
       dv(9)  =   dv(9)  + fs * (-dv2 * r2zij / r2)
 



       else if ((r2.ge.in_bound).and.(r3.ge.in_bound).and. &        
        (r1.le.in_bound)) then             
        !print*,'in  297'


        ! O2 potential                      
                    
        call potohp_f3(r1,V2,dV2)        
       ! Switch function                    
        fs=f2(r2,r3)  
        call df2(r2,r3,dfsdr2,dfsdr3)



        ! for the diatomic potencials

            call pothoq_f17(r2,Vdiat,dVdiat)



            dv(1) = dv(1) + fs * dVdiat*r2xij/r2 
            dv(2) = dv(2) + fs * dVdiat*r2yij/r2 
            dv(3) = dv(3) + fs * dVdiat*r2zij/r2 
                          
            dv(7) = dv(7) - fs * dVdiat*r2xij/r2 
            dv(8) = dv(8) - fs * dVdiat*r2yij/r2 
            dv(9) = dv(9) - fs * dVdiat*r2zij/r2 


            V4 = V4 + Vdiat


 
            call poto2q_f18(r3,Vdiat,dVdiat)


            dv(4) = dv(4) + fs * dVdiat*r3xij/r3 
            dv(5) = dv(5) + fs * dVdiat*r3yij/r3 
            dv(6) = dv(6) + fs * dVdiat*r3zij/r3 
                                    
            dv(7) = dv(7) - fs * dVdiat*r3xij/r3 
            dv(8) = dv(8) - fs * dVdiat*r3yij/r3 
            dv(9) = dv(9) - fs * dVdiat*r3zij/r3 


            V4 = V4 + Vdiat



        V=V1*(1-fs)+(V2+V4)*fs


       dv=dv + dv1*(1-fs) 
              
       dv(1)  =   dv(1)  + (v2+v4-v1) * ( dfsdr2 *  r2xij / r2)
       dv(2)  =   dv(2)  + (v2+v4-v1) * ( dfsdr2 *  r2yij / r2)
       dv(3)  =   dv(3)  + (v2+v4-v1) * ( dfsdr2 *  r2zij / r2)
                              
       dv(4)  =   dv(4)  + (v2+v4-v1) * ( dfsdr3 *  r3xij / r3)
       dv(5)  =   dv(5)  + (v2+v4-v1) * ( dfsdr3 *  r3yij / r3)
       dv(6)  =   dv(6)  + (v2+v4-v1) * ( dfsdr3 *  r3zij / r3)
                              
       dv(7)  =   dv(7)  + (v2+v4-v1) * (-dfsdr2 *  r2xij / r2- dfsdr3 *  r3xij / r3)
       dv(8)  =   dv(8)  + (v2+v4-v1) * (-dfsdr2 *  r2yij / r2- dfsdr3 *  r3yij / r3)
       dv(9)  =   dv(9)  + (v2+v4-v1) * (-dfsdr2 *  r2zij / r2- dfsdr3 *  r3zij / r3)
                              
                              
       dv(1)  =   dv(1)  + fs * ( dv2 * r1xij / r1)
       dv(2)  =   dv(2)  + fs * ( dv2 * r1yij / r1)
       dv(3)  =   dv(3)  + fs * ( dv2 * r1zij / r1)
                              
       dv(4)  =   dv(4)  + fs * (-dv2 * r1xij / r1)
       dv(5)  =   dv(5)  + fs * (-dv2 * r1yij / r1)
       dv(6)  =   dv(6)  + fs * (-dv2 * r1zij / r1)


   




 
         else!  end if
               
         v=v1
           
          dv  = dv1   
        
          group(i,9)=0               

      end if               



!print*,'whatsup2'


      return

      end subroutine


                                                                                                   
!    ******************************************************************                             
!    TO COMPUTE THE HO2 SURFACE in atomic units                                                     
!    ******************************************************************                             
      FUNCTION VHO2_f4(R1J,R2J,R3J)                                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                                           
      COMMON/THRBOD_f4/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3                                           
                                                                                                    
                                                                                                    
                                                                                                    
!      print*,'funcao HO2 antes',r1j,r2j,r3j                                                        
      R1=R1J/5.2917721092d-11                                                                       
      R2=R2J/5.2917721092d-11                                                                       
      R3=R3J/5.2917721092d-11                                                                       
!      print*,'funcao HO2 depois',r1,r2,r3                                                          
      Q1=1.0D0/DSQRT(3.0D0)*(R1+R2+R3)                                                              
      Q2=1.0D0/DSQRT(2.0D0)*(R2-R3)                                                                 
      Q3=1.0D0/DSQRT(6.0D0)*(2.0D0*R1-R2-R3)                                                        
      V=voo_f4(R1)+VOH_f4(R2)+VOH_f4(R3)+THREBQ_f4(Q1,Q2,Q3)+ &                                     
        EXDIS_f4(R1,R2,R3)+ELECT_f4(R1,R2,R3)+VSPECT_f4(R1,R2,R3)                                   
                                                                                                    
      VHO2_f4=V*4.3597482D-18                                                                       
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
      FUNCTION VSPECT_f4(r1,r2,r3)                                                                  
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
                                                                                                    
      t=c1+c2*r1d+c3*s1d+c4*r1d**2+c5*s1d**2+c6*r1d*s1d+ &                                          
        c7*s2d2+c8*r1d**3+c9*r1d*s1d**2+c10*r1d**2*s1d+ &                                           
        c11*r1d*s2d2+c12*s1d**3+c13*s1d*s2d2                                                        
      dec=c1fix*r1d**2+c2fix*s1d**2+c3fix*r1d*s1d+c4fix*s2d2**2                                     
                                                                                                    
      VSPECT_f4=t*exp(-dec)                                                                         
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!                                                                                                   
!    TO COMPUTE THE THREE BODY TERM IN SIMETRIC COORDINATES Q1,Q2,Q3                                
!                                                                                                   
!    ****************************************************************                               
      FUNCTION THREBQ_f4(Q1,Q2,Q3)                                                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      COMMON/COEFF_f4/C(52)                                                                         
      COMMON/THRBOD_f4/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3                                           
      COMMON/REFGEO_f4/R10,R20,R30                                                                  
                                                                                                    
                                                                                                    
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
      POLQ=C(1)*Q1+C(2)*Q12+C(3)*TQ1+C(4)*Q13+C(5)*Q1*TQ1+ &                                        
      C(6)*Q3*TQ2+C(7)*Q14+C(8)*Q12*TQ1+C(9)*TQ1**2+C(10)*Q1*Q3*TQ2+ &                              
      C(11)*Q3+C(12)*Q1*Q3+C(13)*TQ3+C(14)*Q12*Q3+C(15)*Q1*TQ3+ &                                   
      C(16)*Q3*TQ1+C(17)*Q13*Q3+C(18)*Q12*TQ3+C(19)*Q1*Q3*TQ1+ &                                    
      C(20)*Q32*TQ2+C(21)*TQ1*TQ3+C(22)+C(23)*Q15+C(24)*Q13*TQ1+ &                                  
      C(25)*Q1*TQ12+C(26)*Q12*Q3*TQ2+C(27)*Q3*TQ1*TQ2+C(28)*Q14*Q3+ &                               
      C(29)*Q13*TQ3+C(30)*Q12*Q3*TQ1+C(31)*Q1*Q32*TQ2+C(32)*Q1*TQ1*TQ3+ &                           
      C(33)*Q3*TQ12+C(34)*Q3*TQ2*TQ3+C(35)*Q16+C(36)*Q14*TQ1+ &                                     
      C(37)*Q12*TQ12+C(38)*Q13*Q3*TQ2+C(39)*Q1*Q3*TQ1*TQ2+C(40)*TQ13+ &                             
      C(41)*Q32*TQ22+C(42)*Q15*Q3+C(43)*Q14*TQ3+C(44)*Q13*Q3*TQ1+ &                                 
      C(45)*Q12*Q32*TQ2+C(46)*Q12*TQ1*TQ3+C(47)*Q1*Q3*TQ12+ &                                       
      C(48)*Q1*Q3*TQ2*TQ3+C(49)*Q32*TQ1*TQ2+C(50)*TQ12*TQ3                                          
      DECAY1=1.0D0-DTANH(C(51)*S1)                                                                  
      DECAY2=1.0D0-DTANH(C(52)*S2)                                                                  
      DECAY3=1.0D0-DTANH(C(52)*S3)                                                                  
      THREBQ_f4=POLQ*DECAY1*DECAY2*DECAY3                                                           
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!                                                                                                   
!    TO COMPUTE THE HFACE FOR  O...H                                                                
!                                                                                                   
!    ****************************************************************                               
      FUNCTION VOH_f4(R)                                                                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
                                                                                                    
      VOH_f4=EHFOH_f4(R)+DISOH_f4(R)                                                                
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!                                                                                                   
!    TO COMPUTE THE EHF FOR  O...H                                                                  
!                                                                                                   
!    ****************************************************************                               
      FUNCTION EHFOH_f4(R)                                                                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      DIMENSION ASV(4)                                                                              
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      SAVE /DIATDI_f4/,D,ASV                                                                        
      DATA D,ASV/0.13825385D0,2.6564788D0,1.7450528D0,0.71014391D0, &                               
                 2.5453276D0/                                                                       
                                                                                                    
      X=R-RMOH                                                                                      
      R2=X*X                                                                                        
      R3=R2*X                                                                                       
      EHFOH_f4=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*DEXP(-ASV(4)*X)                              
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!                                                                                                   
!    TO COMPUTE THE DISPERSION FOR  O...H                                                           
!                                                                                                   
!    ****************************************************************                               
      FUNCTION DISOH_f4(R)                                                                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      COMMON/DISPC_f4/COO(10),COH(10)                                                               
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      SAVE /DISPC_f4/,/DIATDI_f4/                                                                   
                                                                                                    
      DISOH_f4=disp_f4(R,COH(6),COH(8),COH(10),R0OH,RMOH)                                           
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!                                                                                                   
!    TO COMPUTE THE HFACE FOR  O...O                                                                
!                                                                                                   
!    ****************************************************************                               
      FUNCTION voo_f4(R)                                                                            
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
                                                                                                    
      voo_f4=ehfoo_f4(R)+disoo_f4(R)                                                                
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!                                                                                                   
!    TO COMPUTE THE EHF FOR  O...O                                                                  
!                                                                                                   
!    ****************************************************************                               
      FUNCTION ehfoo_f4(R)                                                                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      DIMENSION ASV(4)                                                                              
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      SAVE /DIATDI_f4/,D,ASV                                                                        
      DATA D,ASV/0.14291202D0,3.6445906D0,3.9281238D0,2.0986689D0, &                                
                 3.3522498D0/                                                                       
                                                                                                    
      X=R-RMOO                                                                                      
      R2=X*X                                                                                        
      R3=R2*X                                                                                       
      ehfoo_f4=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)*DEXP(-ASV(4)*X)                              
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!                                                                                                   
!    TO COMPUTE THE DISPERSION FOR  O...O                                                           
!                                                                                                   
!    ****************************************************************                               
      FUNCTION disoo_f4(R)                                                                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      COMMON/DISPC_f4/COO(10),COH(10)                                                               
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      SAVE /DISPC_f4/,/DIATDI_f4/                                                                   
                                                                                                    
      disoo_f4=disp_f4(R,COO(6),COO(8),COO(10),R0OO,RMOO)                                           
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!    TO COMPUTE THE EXCHANGE - DISPERSION TERM                                                      
!    ****************************************************************                               
      FUNCTION EXDIS_f4 (R1,R2,R3)                                                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      IMPLICIT INTEGER(I-N)                                                                         
      COMMON/DISPC_f4/COO(10),COH(10)                                                               
      COMMON/RKVAL_f4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)                                       
      COMMON/DISCO_f4/CEFOO(10),CEFOH2(10),CEFOH3(10), &                                            
       CEDOO(10),CEDOH2(10)                                                                         
      COMMON/DISCO2_f4/CEDOH3(10)                                                                   
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      SAVE /DISPC_f4/,/DIATDI_f4/,/RKVAL_f4/                                                        
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    
      DO 10 IN=6,10,2                                                                               
        CEFOO(IN)=CEF_f4(COO(IN),RK0OH(IN),RK1OH(IN),RK0OH(IN), &                                   
                  RK1OH(IN),RMOH,RMOH,R2,R3)                                                        
        CEDOO(IN)=CEFOO(IN)-COO(IN)                                                                 
        CEFOH2(IN)=CEF_f4(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN), &                                  
                   RK1OH(IN),RMOO,RMOH,R1,R3)                                                       
        CEDOH2(IN)=CEFOH2(IN)-COH(IN)                                                               
        CEFOH3(IN)=CEF_f4(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN), &                                  
                   RK1OH(IN),RMOO,RMOH,R1,R2)                                                       
        CEDOH3(IN)=CEFOH3(IN)-COH(IN)                                                               
   10 CONTINUE                                                                                      
      EXDIS_f4=disp_f4(R1,CEDOO(6),CEDOO(8),CEDOO(10),R0OO,RMOO) &                                  
           +disp_f4(R2,CEDOH2(6),CEDOH2(8),CEDOH2(10),R0OH,RMOH) &                                  
           +disp_f4(R3,CEDOH3(6),CEDOH3(8),CEDOH3(10),R0OH,RMOH)                                    
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!    TO COMPUTE THE EFECTIVE Cn                                                                     
!    ****************************************************************                               
      FUNCTION CEF_f4(CAS,RK01,RK11,RK02,RK12,RE1,RE2,R1,R2)                                        
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
                                                                                                    
      CEF_f4=0.5D0*CAS*((1.0D0-RK01*DEXP(-RK11*(R1-RE1)))* &                                        
        DTANH(RK12*R2)+(1.0D0-RK02*DEXP(-RK12*(R2-RE2)))*DTANH(RK11*R1))                            
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!    TO COMPUTE THE ELECT_f4ROSTATIC TERM                                                           
!    ****************************************************************                               
      FUNCTION ELECT_f4(R1,R2,R3)                                                                   
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      COMMON/POLAR_f4/C4,C5                                                                         
      COMMON/RKVAL_f4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)                                       
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      COMMON/DAMPC_f4/ADAMP(10),BDAMP(10)                                                           
      COMMON/WELECT_f4/C4OHR2,C5OHR2,C4OHR3,C5OHR3,C4OO,C5OO,TERM4,TERM5                            
      SAVE /DIATDI_f4/,/RKVAL_f4/,/POLAR_f4/,/DAMPC_f4/                                             
                                                                                                    
                                                                                                    
                                                                                                    
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
      TERM4=C4OO/R14*(1.0D0-DEXP(-ADAMP(4)*RROO-BDAMP(4)*RROO**2))**4+ &                            
       C4OHR2/R24*(1.0D0-DEXP(-ADAMP(4)*RROH2-BDAMP(4)*RROH2**2))**4+ &                             
       C4OHR3/R34*(1.0D0-DEXP(-ADAMP(4)*RROH3-BDAMP(4)*RROH3**2))**4                                
      TERM5=C5OO/R15*(1.0D0-DEXP(-ADAMP(5)*RROO-BDAMP(5)*RROO**2))**5+ &                            
       C5OHR2/R25*(1.0D0-DEXP(-ADAMP(5)*RROH2-BDAMP(5)*RROH2**2))**5+ &                             
       C5OHR3/R35*(1.0D0-DEXP(-ADAMP(5)*RROH3-BDAMP(5)*RROH3**2))**5                                
      ELECT_f4=TERM4+TERM5                                                                          
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ***************************************************************                                
      FUNCTION disp_f4(R,C6,C8,C10,R0,RM)                                                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      COMMON/DAMPC_f4/ADAMP(10),BDAMP(10)                                                           
      SAVE /DAMPC_f4/                                                                               
                                                                                                    
      R6=R**6                                                                                       
      R8=R6*R*R                                                                                     
      R10=R8*R*R                                                                                    
      RR=2.0D0*R/(RM+2.5D0*R0)                                                                      
      D6=(1.0D0-DEXP(-ADAMP(6)*RR-BDAMP(6)*RR*RR))**6                                               
      D8=(1.0D0-DEXP(-ADAMP(8)*RR-BDAMP(8)*RR*RR))**8                                               
      D10=(1.0D0-DEXP(-ADAMP(10)*RR-BDAMP(10)*RR*RR))**10                                           
      disp_f4=-C6/R6*D6-C8/R8*D8-C10/R10*D10                                                        
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    
!    ***************************************************************                                
!    DATA FOR HO2 SURFACE                                                                           
!    ***************************************************************                                
      BLOCK DATA HO2DAT_f4                                                                          
                                                                                                    
!    ESTE BLOCK DATA CONTIENE LOS DATOS CORRESPONDIENTES A                                          
!    LA SUPERFICIE CALCULADA CON BOO=BOH=1.55                                                       
                                                                                                    
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      COMMON/COEFF_f4/C(52)                                                                         
      COMMON/DISPC_f4/COO(10),COH(10)                                                               
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      COMMON/RKVAL_f4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)                                       
      COMMON/POLAR_f4/C4,C5                                                                         
      COMMON/DAMPC_f4/ADAMP(10),BDAMP(10)                                                           
      COMMON/REFGEO_f4/R10,R20,R30                                                                  
      SAVE /COEFF_f4/,/DISPC_f4/,/DIATDI_f4/,/RKVAL_f4/ &                                           
       ,/POLAR_f4/,/DAMPC_f4/,/REFGEO_f4/                                                           
                                                                                                    
      DATA C/ &                                                                                     
        .49040645D+01, -.86748216D+01,  .50555792D+01,  .42941301D+01, &                            
       -.41874792D+01,  .13461379D+00, -.99064922D+00,  .13358488D+01, &                            
        .13495231D+01, -.18529696D+00, -.23534213D+02,  .24289930D+02, &                            
       -.50209026D+01, -.10365484D+02,  .46692224D+01, -.14747138D+01, &                            
        .23119718D+01, -.18247842D+01, -.28472166D+00,  .51036509D+00, &                            
        .19124083D+00,  .45405729D+01,  .11087611D+00, -.19990481D+00, &                            
       -.37356178D+00,  .46142042D-01, -.20565580D+00, -.27015963D+00, &                            
        .34085281D+00,  .28321162D+00, -.11558481D+00, -.29448886D+00, &                            
       -.52932488D+00,  .58159523D-01, -.48649560D-02,  .11949167D-01, &                            
        .21409804D-01, -.20620608D-02,  .30177088D-01,  .27880291D-01, &                            
        .88458711D-02,  .13137410D-01, -.24705619D-01, -.31085889D-01, &                            
        .34317857D-02,  .52593878D-01,  .79500714D-01, -.79782216D-02, &                            
        .31164575D-01, -.28737598D-01,  .98201698D+00,  .62000000D+00/                              
                                                                                                    
      DATA R0OO,RMOO,R0OH,RMOH/5.661693D0,2.2818D0,6.294894D0,1.8344D0/                             
      DATA COO/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,15.40D0,0.0D0,235.219943D0, &                           
               0.0D0,4066.23929D0/                                                                  
      DATA COH/0.0D0,0.0D0,0.0D0,0.0D0,0.D0,10.00D0,0.0D0,180.447673D0, &                           
               0.0D0,3685.25842D0/                                                                  
      DATA C4,C5/-0.92921D0,-1.79000D0/                                                             
      DATA RK0OO/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,-.27847758D0,0.0D0, &                                
                 -.46815641D0,0.0D0,-1.20506384D0/                                                  
      DATA RK1OO/0.0D0,0.0D0,0.0D0,3.35224980D0,3.35224980D0, &                                     
                 0.95273753D0,0.0D0,0.94148408D0,0.0D0,0.72379129D0/                                
      DATA RK0OH/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.02465005D0,0.0D0, &                                
                 0.05036950D0,0.0D0,0.06294371D0/                                                   
      DATA RK1OH/0.0D0,0.0D0,0.0D0,2.54532760D0,2.54532760D0, &                                     
                 0.68758097D0,0.0D0,0.82542359D0,0.0D0,0.94034225D0/                                
      DATA ADAMP/0.0D0,0.0D0,0.0D0,5.0079875D0,3.8428294D0,3.0951333D0, &                           
                 0.0D0,2.1999000D0,0.0D0,1.6880714D0/                                               
      DATA BDAMP/0.0D0,0.0D0,0.D0,10.6645006D0,9.6758155D0,8.7787895D0, &                           
                 0.0D0,7.2265123D0,0.0D0,5.9487108D0/                                               
      DATA R10,R20,R30/2.5143000D0,2.6469057D0,2.6469057D0/                                         
      END                                                                                           
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    
!   ********************************************************************                            
!    Here starts the derivative part of the program                                                 
!    It needs the function THREBQ_f4(Q1,Q2,Q3) and the BLOCK DATA HO2DAT                            
!   ********************************************************************                            
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    
      SUBROUTINE DERVHO2_f4(r1j,r2j,r3j,dg1,dg2,dg3)                                                
!   **************************************************************                                  
!    TO COMPUTE THE DERIVATIVE H2O SURFACE in atomic units                                          
!   **************************************************************                                  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      DIMENSION X(3),G(3)                                                                           
      COMMON/COEFF_f4/C(52)                                                                         
      COMMON/THRBOD_f4/POLQ,DECAY1,DECAY2,DECAY3,R1,R2,R3                                           
      COMMON/REFGEO_f4/R10,R20,R30                                                                  
                                                                                                    
                                                                                                    
                                                                                                    
!   **************************************************************                                  
!      print*,'funcao derHO2 antes',r1J,r2J,r3J                                                     
      R1=R1J/5.2917721092d-11                                                                       
      R2=R2J/5.2917721092d-11                                                                       
      R3=R3J/5.2917721092d-11                                                                       
                                                                                                    
!      print*,'funcao derHO2 depois ',r1,r2,r3                                                      
      Q1=1.0D0/DSQRT(3.0D0)*(R1+R2+R3)                                                              
      Q2=1.0D0/DSQRT(2.0D0)*(R2-R3)                                                                 
      Q3=1.0D0/DSQRT(6.0D0)*(2.0D0*R1-R2-R3)                                                        
      X(1)=r1                                                                                       
      X(2)=r2                                                                                       
      X(3)=r3                                                                                       
      term3Q=THREBQ_f4(Q1,Q2,Q3)                                                                    
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
      DPOQ1=C(1)+2.0D0*C(2)*Q1+3.0D0*C(4)*Q12+C(5)*TQ1+4.0D0*C(7)*Q13+ &                            
         2.0D0*C(8)*Q1*TQ1+C(10)*Q3*TQ2+C(12)*Q3+2.0D0*C(14)*Q1*Q3+ &                               
         C(15)*TQ3+3.0D0*C(17)*Q12*Q3+2.0D0*C(18)*Q1*TQ3+C(19)*Q3*TQ1+ &                            
         5.0D0*C(23)*Q14+3.0D0*C(24)*Q12*TQ1+C(25)*TQ12+2.0D0*C(26)* &                              
         Q1*Q3*TQ2+4.0D0*C(28)*Q13*Q3+3.0D0*C(29)*Q12*TQ3+2.0D0*C(30)* &                            
         Q1*Q3*TQ1+C(31)*Q32*TQ2+C(32)*TQ1*TQ3+6.0D0*C(35)*Q15+4.0D0* &                             
         C(36)*Q13*TQ1+2.0D0*C(37)*Q1*TQ12+3.0D0*C(38)*Q12*Q3*TQ2+ &                                
         C(39)*Q3*TQ1*TQ2+5.0D0*C(42)*Q14*Q3+4.0D0*C(43)*Q13*TQ3+ &                                 
         3.0D0*C(44)*Q12*Q3*TQ1+2.0D0*C(45)*Q1*Q32*TQ2+ &                                           
         2.0D0*C(46)*Q1*TQ1*TQ3+C(47)*Q3*TQ12+C(48)*Q3*TQ2*TQ3                                      
      DPOQ2=C(3)*2.0D0*Q2+2.0D0*C(5)*Q1*Q2+C(6)*Q3*(-6.0D0*Q2)+ &                                   
         C(8)*Q12*2.0D0*Q2+C(9)*2.0D0*TQ1*2.0D0*Q2+C(10)*Q1*Q3* &                                   
         (-6.0D0*Q2)+C(13)*2.0D0*Q2+C(15)*Q1*2.0D0*Q2+C(16)*Q3*2.0D0*Q2+ &                          
         C(18)*Q12*2.0D0*Q2+C(19)*Q1*Q3*2.0D0*Q2+C(20)*Q32*(-6.0D0*Q2)+ &                           
         C(21)*2.0D0*Q2*TQ3+C(21)*TQ1*2.0D0*Q2+C(24)*Q13*2.0D0*Q2+ &                                
         C(25)*Q1*2.0D0*TQ1*2.0D0*Q2+C(26)*Q12*Q3*(-6.0D0*Q2)+ &                                    
         C(27)*Q3*2.0D0*Q2*TQ2+C(27)*Q3*TQ1*(-6.0D0*Q2)+ &                                          
         C(29)*Q13*2.0D0*Q2+C(30)*Q12*Q3*2.0D0*Q2+C(31)*Q1*Q32*(-6.0D0 &                            
         *Q2)+C(32)*Q1*2.0D0*Q2*(TQ1+TQ3)+C(33)*Q3*4.0D0*TQ1*Q2+C(34)* &                            
         Q3*Q2*(2.0D0*TQ2-6.0D0*TQ3)+2.0D0*C(36)*Q14*Q2+4.0D0*C(37)*Q12* &                          
         Q2*TQ1+C(38)*Q13*Q3*(-6.0D0*Q2)+C(39)*Q1*Q3*Q2*(2.0D0*TQ2-6.0D0 &                          
         *TQ1)+C(40)*3.0D0*TQ12*2.0D0*Q2+C(41)*Q32*2.0D0*TQ2*(-6.0D0*Q2) &                          
         +C(43)*Q14*2.0D0*Q2+C(44)*Q13*Q3*2.0D0*Q2+C(45)*Q12*Q32*(-6.0D0 &                          
         *Q2)+C(46)*Q12*Q2*2.0D0*(TQ3+TQ1)+C(47)*Q1*Q3*2.0D0*TQ1*2.0D0* &                           
         Q2+C(48)*Q1*Q3*Q2*(2.0D0*TQ2-6.0D0*TQ3)+C(49)*Q32*Q2*(2.0D0*TQ2 &                          
         -6.0D0*TQ1)+C(50)*2.0D0*Q2*TQ1*(2.0D0*TQ3+TQ1)                                             
      DTQ=2.0D0*Q3                                                                                  
      DPOQ3=C(3)*DTQ+C(5)*Q1*DTQ+C(6)*TQ2+C(6)*Q3*DTQ+C(8)*Q12*DTQ+ &                               
         C(9)*TQ1*2.0D0*DTQ+C(10)*Q1*TQ2+C(10)*Q1*Q3*DTQ+C(11)+ &                                   
         C(12)*Q1-C(13)*DTQ+C(14)*Q12+C(15)*Q1*(-DTQ)+C(16)*TQ1+ &                                  
         C(16)*Q3*DTQ+C(17)*Q13+C(18)*Q12*(-DTQ)+C(19)*Q1*(TQ1+Q3* &                                
         DTQ)+C(20)*2.0D0*Q3*TQ2+C(20)*Q32*DTQ+C(21)*(DTQ*TQ3-TQ1*DTQ) &                            
         +C(24)*Q13*DTQ+C(25)*Q1*2.0D0*TQ1*DTQ+C(26)*Q12*(TQ2+Q3*DTQ)+ &                            
        C(27)*(TQ1*TQ2+Q3*DTQ*TQ2+Q3*TQ1*DTQ)+C(28)*Q14+C(29)*Q13*(-DTQ) &                          
       +C(30)*Q12*(TQ1+Q3*DTQ)+C(31)*Q1*(2.0D0*Q3*TQ2+Q32*DTQ)+C(32)*Q1* &                          
         (DTQ*TQ3-DTQ*TQ1)+C(33)*(TQ12+2.0D0*Q3*TQ1*DTQ)+C(34)*TQ2*TQ3+ &                           
         C(34)*Q3*DTQ*(TQ3-TQ2)+C(36)*Q14*DTQ+C(37)*Q12*2.0D0*TQ1*DTQ+ &                            
         C(38)*Q13*(TQ2+Q3*DTQ)+C(39)*Q1*TQ1*TQ2+C(39)*Q1*Q3*DTQ* &                                 
         (TQ1+TQ2)+C(40)*3.0D0*TQ12*DTQ+C(41)*DTQ*TQ2*(TQ2+2.0D0*Q32)+ &                            
        C(42)*Q15-C(43)*Q14*DTQ+C(44)*Q13*(TQ1+Q3*DTQ)+C(45)*Q12*DTQ* &                             
        (TQ2+Q32)+C(46)*Q12*DTQ*(TQ3-TQ1)+C(47)*Q1*(TQ12+TQ1*4.0D0*Q32)+ &                          
         C(48)*Q1*TQ2*TQ3+C(48)*Q1*Q3*DTQ*(TQ3-TQ2)+C(49)*DTQ*TQ1*TQ2+ &                            
         C(49)*Q32*DTQ*(TQ1+TQ2)+C(50)*DTQ*(2.0D0*TQ1*TQ3-TQ12)                                     
      CALL  DEREXDIS_f4(X(1),X(2),X(3),DER1,DER2,DER3)                                              
      CALL  DERELECT_f4(X(1),X(2),X(3),DERI1,DERI2,DERI3)                                           
      CALL  DERVSPECT_f4(X(1),X(2),X(3),DERV1,DERV2,DERV3)                                          
      G(1)=((DPOQ1*DQ1R1+DPOQ2*DQ2R1+DPOQ3*DQ3R1)*DECAY1-POLQ*C(51)/ &                              
         COSH(C(51)*S1)**2)*DECAY2*DECAY3+Dvoo_f4(X(1))+DER1+DERI1+DERV1                            
      G(2)=((DPOQ1*DQ1R2+DPOQ2*DQ2R2+DPOQ3*DQ3R2)*DECAY2-POLQ*C(52)/ &                              
         COSH(C(52)*S2)**2)*DECAY1*DECAY3+DVOH_f4(X(2))+DER2+DERI2+DERV2                            
      G(3)=((DPOQ1*DQ1R3+DPOQ2*DQ2R3+DPOQ3*DQ3R3)*DECAY3-POLQ*C(52)/ &                              
         COSH(C(52)*S3)**2)*DECAY1*DECAY2+DVOH_f4(X(3))+DER3+DERI3+DERV3                            
                                                                                                    
                                                                                                    
      dg1=G(1)*4.3597482D-18/5.2917721092d-11                                                       
      dg2=G(2)*4.3597482D-18/5.2917721092d-11                                                       
      dg3=G(3)*4.3597482D-18/5.2917721092d-11                                                       
                                                                                                    
      END                                                                                           
                                                                                                    
      SUBROUTINE DERVSPECT_f4(R1,R2,R3,DERV1,DERV2,DERV3)                                           
!    ****************************************************************                               
!                                                                                                   
!    TO COMPUTE THE DERIVATIVES OF A GAUSSIAN TERM TO FIX                                           
!    VIBRATIONAL SPECTRA                                                                            
!                                                                                                   
!    ****************************************************************                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Y)                                                            
                                                                                                    
      data r1eq,r2eq,r3eq/2.5143d0,1.8346d0,3.4592d0/                                               
                                                                                                    
!    ****************************************************************                               
!                                                                                                   
                                                                                                    
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
                                                                                                    
                                                                                                    
      t=c1+c2*r1d+c3*s1d+c4*r1d**2+c5*s1d**2+c6*r1d*s1d+ &                                          
        c7*s2d2+c8*r1d**3+c9*r1d*s1d**2+c10*r1d**2*s1d+ &                                           
        c11*r1d*s2d2+c12*s1d**3+c13*s1d*s2d2                                                        
      dec=c1fix*r1d**2+c2fix*s1d**2+c3fix*r1d*s1d+c4fix*s2d2**2                                     
                                                                                                    
      VSPECT_f4=t*exp(-dec)                                                                         
                                                                                                    
                                                                                                    
!     VSPECT_f4=t*exp(-dec)                                                                         
                                                                                                    
       DERV1=(c2+2.0d0*c4*r1d+c6*s1d+ &                                                             
            3.0d0*c8*r1d**2+c9*s1d**2+2.0d0*c10*r1d*s1d+ &                                          
            c11*s2d2 &                                                                              
           -t*(2.0d0*c1fix*r1d+c3fix*s1d))*exp(-dec)                                                
                                                                                                    
       DERVS1=(c3+2.0d0*c5*s1d+c6*r1d+ &                                                            
         2.0d0*c9*r1d*s1d+c10*r1d**2+ &                                                             
         3.0d0*c12*s1d**2+c13*s2d2 &                                                                
             -t*(2.0d0*c2fix*s1d+c3fix*r1d))*exp(-dec)                                              
                                                                                                    
       DERVS2=(c7+c11*r1d+c13*s1d &                                                                 
         -t*c4fix*2.0D0*s2d2)*exp(-dec)                                                             
                                                                                                    
       DS1DR2=1.0D0/SQRT(2.0D0)                                                                     
       DS1DR3=1.0D0/SQRT(2.0D0)                                                                     
       DS2DR2=r2-r3                                                                                 
       DS2DR3=r3-r2                                                                                 
                                                                                                    
       DERV2=DERVS1*DS1DR2+DERVS2*DS2DR2                                                            
       DERV3=DERVS1*DS1DR3+DERVS2*DS2DR3                                                            
                                                                                                    
                                                                                                    
      RETURN                                                                                        
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!    TO COMPUTE THE DER. OF THE HFACE FOR  O...H                                                    
!    ****************************************************************                               
      FUNCTION DVOH_f4(R)                                                                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      DIMENSION ASV(4)                                                                              
      COMMON/DISPC_f4/COO(10),COH(10)                                                               
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      DATA D,ASV/0.13825385D0,2.6564788D0,1.7450528D0,0.71014391D0, &                               
                 2.5453276D0/                                                                       
      SAVE /DISPC_f4/,/DIATDI_f4/,D,ASV                                                             
                                                                                                    
      X=R-RMOH                                                                                      
      R2=X*X                                                                                        
      R3=R2*X                                                                                       
      POL=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)                                                   
      DPOL=-D*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*R2)                                               
      POT=EXP(-ASV(4)*X)                                                                            
      DVOH_f4=-ASV(4)*POT*POL+DPOL*POT+Ddisp_f4(R,COH(6),COH(8),COH(10), &                          
         R0OH,RMOH)                                                                                 
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!    TO COMPUTE THE DER. OF THE HFACE FOR  O...O                                                    
!    ****************************************************************                               
      FUNCTION Dvoo_f4(R)                                                                           
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      DIMENSION ASV(4)                                                                              
      COMMON/DISPC_f4/COO(10),COH(10)                                                               
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      DATA D,ASV/0.14291202D0,3.6445906D0,3.9281238D0,2.0986689D0, &                                
                 3.3522498D0/                                                                       
      SAVE /DISPC_f4/,/DIATDI_f4/,D,ASV                                                             
                                                                                                    
      X=R-RMOO                                                                                      
      R2=X*X                                                                                        
      R3=R2*X                                                                                       
      POL=-D*(1.0D0+ASV(1)*X+ASV(2)*R2+ASV(3)*R3)                                                   
      DPOL=-D*(ASV(1)+2.0D0*ASV(2)*X+3.0D0*ASV(3)*R2)                                               
      POT=EXP(-ASV(4)*X)                                                                            
      Dvoo_f4=-ASV(4)*POT*POL+DPOL*POT+Ddisp_f4(R,COO(6),COO(8),COO(10), &                          
         R0OO,RMOO)                                                                                 
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!    TO COMPUTE THE DERIVATIVES OF THE EXCHANGE - DISPERSION TERM                                   
!    ****************************************************************                               
      SUBROUTINE DEREXDIS_f4(R1,R2,R3,DER1,DER2,DER3)                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      IMPLICIT INTEGER(I-N)                                                                         
      COMMON/DISPC_f4/COO(10),COH(10)                                                               
      COMMON/RKVAL_f4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)                                       
      COMMON/DISCO_f4/CEFOO(10),CEFOH2(10),CEFOH3(10), &                                            
        CEDOO(10),CEDOH2(10)                                                                        
      COMMON/DISCO2_f4/CEDOH3(10)                                                                   
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      COMMON/DAMPC_f4/ADAMP(10),BDAMP(10)                                                           
      SAVE /DISPC_f4/,/DIATDI_f4/,/RKVAL_f4/                                                        
                                                                                                    
                                                                                                    
                                                                                                    
                                                                                                    
      DO 10 IN=6,10,2                                                                               
      CEFOO(IN)=CEF_f4(COO(IN),RK0OH(IN),RK1OH(IN),RK0OH(IN),RK1OH(IN), &                           
         RMOH,RMOH,R2,R3)                                                                           
      CEDOO(IN)=CEFOO(IN)-COO(IN)                                                                   
      CEFOH2(IN)=CEF_f4(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN), &                          
         RMOO,RMOH,R1,R3)                                                                           
      CEDOH2(IN)=CEFOH2(IN)-COH(IN)                                                                 
      CEFOH3(IN)=CEF_f4(COH(IN),RK0OO(IN),RK1OO(IN),RK0OH(IN),RK1OH(IN), &                          
         RMOO,RMOH,R1,R2)                                                                           
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
      CALL DCEF_f4(COO(6),RK0OH(6),RK1OH(6),RK0OH(6),RK1OH(6),RMOH, &                               
           RMOH,R2,R3,DC61R2,DC61R3)                                                                
      CALL DCEF_f4(COO(8),RK0OH(8),RK1OH(8),RK0OH(8),RK1OH(8),RMOH, &                               
           RMOH,R2,R3,DC81R2,DC81R3)                                                                
      CALL DCEF_f4(COO(10),RK0OH(10),RK1OH(10),RK0OH(10),RK1OH(10),RMOH, &                          
           RMOH,R2,R3,D101R2,D101R3)                                                                
      CALL DCEF_f4(COH(6),RK0OO(6),RK1OO(6),RK0OH(6),RK1OH(6),RMOO, &                               
           RMOH,R1,R3,DC62R1,DC62R3)                                                                
      CALL DCEF_f4(COH(8),RK0OO(8),RK1OO(8),RK0OH(8),RK1OH(8),RMOO, &                               
           RMOH,R1,R3,DC82R1,DC82R3)                                                                
      CALL DCEF_f4(COH(10),RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10),RMOO, &                          
           RMOH,R1,R3,D102R1,D102R3)                                                                
      CALL DCEF_f4(COH(6),RK0OO(6),RK1OO(6),RK0OH(6),RK1OH(6),RMOO, &                               
           RMOH,R1,R2,DC63R1,DC63R2)                                                                
      CALL DCEF_f4(COH(8),RK0OO(8),RK1OO(8),RK0OH(8),RK1OH(8),RMOO, &                               
           RMOH,R1,R2,DC83R1,DC83R2)                                                                
      CALL DCEF_f4(COH(10),RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10),RMOO, &                          
           RMOH,R1,R2,D103R1,D103R2)                                                                
      DER1=Ddisp_f4(R1,CEDOO(6),CEDOO(8),CEDOO(10),R0OO,RMOO)-DC62R1* &                             
           T6R2-DC63R1*T6R3-DC82R1*T8R2-DC83R1*T8R3-D102R1*T10R2- &                                 
           D103R1*T10R3                                                                             
      DER2=Ddisp_f4(R2,CEDOH2(6),CEDOH2(8),CEDOH2(10),R0OH,RMOH)-DC61R2* &                          
           T6R1-DC63R2*T6R3-DC81R2*T8R1-DC83R2*T8R3-D101R2*T10R1- &                                 
           D103R2*T10R3                                                                             
      DER3=Ddisp_f4(R3,CEDOH3(6),CEDOH3(8),CEDOH3(10),R0OH,RMOH)-DC61R3* &                          
           T6R1-DC62R3*T6R2-DC81R3*T8R1-DC82R3*T8R2-D101R3*T10R1- &                                 
           D102R3*T10R2                                                                             
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!    TO COMPUTE THE DERIVATIVES OF THE EFECTIVE Cn                                                  
!    ****************************************************************                               
      SUBROUTINE DCEF_f4(CAS,RK01,RK11,RK02,RK12,RE1,RE2,R1,R2, &                                   
         DCR1,DCR2)                                                                                 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
                                                                                                    
      DCR1=0.5D0*CAS*(RK11*RK01*EXP(-RK11*(R1-RE1))*TANH(RK12*R2)+ &                                
        (1.0D0-RK02*EXP(-RK12*(R2-RE2)))*RK11/COSH(RK11*R1)**2)                                     
      DCR2=0.5D0*CAS*(RK12*RK02*EXP(-RK12*(R2-RE2))*TANH(RK11*R1)+ &                                
        (1.0D0-RK01*EXP(-RK11*(R1-RE1)))*RK12/COSH(RK12*R2)**2)                                     
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
!    ****************************************************************                               
!    TO COMPUTE THE DERIVATIVES OF THE ELECT_f4ROSTATIC TERM                                        
!    ****************************************************************                               
      SUBROUTINE DERELECT_f4(R1,R2,R3,DER1,DER2,DER3)                                               
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      COMMON/POLAR_f4/C4,C5                                                                         
      COMMON/RKVAL_f4/RK0OO(10),RK1OO(10),RK0OH(10),RK1OH(10)                                       
      COMMON/DIATDI_f4/R0OO,RMOO,R0OH,RMOH                                                          
      COMMON/DAMPC_f4/ADAMP(10),BDAMP(10)                                                           
      SAVE /DIATDI_f4/,/RKVAL_f4/,/POLAR_f4/,/DAMPC_f4/                                             
                                                                                                    
!    RCR2=(R1-R3)/R2                                                                                
!    RCR3=(R1-R2)/R3                                                                                
!    C42=C4*RCR2                                                                                    
!    C43=C4*RCR3                                                                                    
      C42=C4                                                                                        
      C43=C4                                                                                        
!    PR2=3.0D0*RCR2**2-1.0D0                                                                        
!    PR3=3.0D0*RCR3**2-1.0D0                                                                        
!    C52=C5*PR2                                                                                     
!    C53=C5*PR3                                                                                     
      C52=C5                                                                                        
      C53=C5                                                                                        
!    RR2D1=1.0D0/R2                                                                                 
!    RR2D2=(R3-R1)/R2**2                                                                            
!    RR2D3=-1.0D0/R2                                                                                
!    RR3D1=1.0D0/R3                                                                                 
!    RR3D2=-1.0D0/R3                                                                                
!    RR3D3=(R2-R1)/R3**2                                                                            
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
      C4OHR2=.5D0*C43/RMQ*R3**4*EXP(-RK1OH(4)*(R3-RMOH))*TANH(RK1OO(4)* &                           
             R1)                                                                                    
      C5OHR2=.5D0*C53/RMQ*R3**4*EXP(-RK1OH(5)*(R3-RMOH))*TANH(RK1OO(5)* &                           
             R1)                                                                                    
      C4OHR3=.5D0*C42/RMQ*R2**4*EXP(-RK1OH(4)*(R2-RMOH))*TANH(RK1OO(4)* &                           
             R1)                                                                                    
      C5OHR3=.5D0*C52/RMQ*R2**4*EXP(-RK1OH(5)*(R2-RMOH))*TANH(RK1OO(5)* &                           
             R1)                                                                                    
      C4OO=.5D0*C43/RMQ*R3**4*EXP(-RK1OH(4)*(R3-RMOH))*TANH(RK1OH(4)*R2) &                          
         +.5D0*C42/RMQ*R2**4*EXP(-RK1OH(4)*(R2-RMOH))*TANH(RK1OH(4)*R3)                             
      C5OO=.5D0*C53/RMQ*R3**4*EXP(-RK1OH(5)*(R3-RMOH))*TANH(RK1OH(5)*R2) &                          
         +.5D0*C52/RMQ*R2**4*EXP(-RK1OH(5)*(R2-RMOH))*TANH(RK1OH(5)*R3)                             
      DC42R3=C4OHR2*(4.0D0/R3-RK1OH(4))                                                             
      DC42R1=0.5D0*C4/RMQ*R3**4*EXP(-RK1OH(4)*(R3-RMOH))*(RK1OO(4)/ &                               
            COSH(RK1OO(4)*R1)**2)                                                                   
      DC52R3=C5OHR2*(4.0D0/R3-RK1OH(5))                                                             
      DC52R1=0.5D0*C5/RMQ*R3**4*EXP(-RK1OH(5)*(R3-RMOH))*(RK1OO(5)/ &                               
            COSH(RK1OO(5)*R1)**2)                                                                   
      DC43R2=C4OHR3*(4.0D0/R2-RK1OH(4))                                                             
      DC43R1=0.5D0*C4/RMQ*R2**4*EXP(-RK1OH(4)*(R2-RMOH))*(RK1OO(4)/ &                               
            COSH(RK1OO(4)*R1)**2)                                                                   
      DC53R2=C5OHR3*(4.0D0/R2-RK1OH(5))                                                             
      DC53R1=0.5D0*C5/RMQ*R2**4*EXP(-RK1OH(5)*(R2-RMOH))*(RK1OO(5)/ &                               
            COSH(RK1OO(5)*R1)**2)                                                                   
      DC42R2=0.0D0                                                                                  
      DC52R2=0.0D0                                                                                  
      DC43R3=0.0D0                                                                                  
      DC53R3=0.0D0                                                                                  
      DC41R3=0.5D0*C4/RMQ*R3**3*EXP(-RK1OH(4)*(R3-RMOH))*TANH(RK1OH(4)* &                           
             R2)*(4.0D0-R3*RK1OH(4))+ &                                                             
      .5D0*C4/RMQ*R2**4*EXP(-RK1OH(4)*(R2-RMOH))*(RK1OH(4)/COSH(RK1OH(4) &                          
        *R3)**2)                                                                                    
      DC41R2=0.5D0*C4/RMQ*R2**3*EXP(-RK1OH(4)*(R2-RMOH))*TANH(RK1OH(4)* &                           
             R3)*(4.0D0-R2*RK1OH(4))+ &                                                             
      .5D0*C4/RMQ*R3**4*EXP(-RK1OH(4)*(R3-RMOH))*(RK1OH(4)/COSH(RK1OH(4) &                          
        *R2)**2)                                                                                    
      DC51R3=0.5D0*C5/RMQ*R3**3*EXP(-RK1OH(5)*(R3-RMOH))*TANH(RK1OH(5)* &                           
             R2)*(4.0D0-R3*RK1OH(5))+ &                                                             
      .5D0*C5/RMQ*R2**4*EXP(-RK1OH(5)*(R2-RMOH))*(RK1OH(5)/COSH(RK1OH(5) &                          
        *R3)**2)                                                                                    
      DC51R2=0.5D0*C5/RMQ*R2**3*EXP(-RK1OH(5)*(R2-RMOH))*TANH(RK1OH(5)* &                           
             R3)*(4.0D0-R2*RK1OH(5))+ &                                                             
      .5D0*C5/RMQ*R3**4*EXP(-RK1OH(5)*(R3-RMOH))*(RK1OH(5)/COSH(RK1OH(5) &                          
        *R2)**2)                                                                                    
      DC41R1=0.0D0                                                                                  
      DC51R1=0.0D0                                                                                  
      DRR1=RR1/R1                                                                                   
      DRR2=RR2/R2                                                                                   
      DRR3=RR3/R3                                                                                   
      Ddisp1_f4=-4.0D0*C4OO/R1**4*D4R1**3*(D4R1/R1+(1.0D0-D4R1)* &                                  
             (-ADAMP(4)-2.0D0*BDAMP(4)*RR1)*DRR1) &                                                 
            -5.0D0*C5OO/R1**5*D5R1**4*(D5R1/R1+(1.0D0-D5R1)*(-ADAMP(5)- &                           
             2.0D0*BDAMP(5)*RR1)*DRR1)                                                              
      Ddisp2_f4=-4.0D0*C4OHR2/R2**4*D4R2**3*(D4R2/R2+(1.0D0-D4R2) &                                 
             *(-ADAMP(4)-2.0D0*BDAMP(4)*RR2)*DRR2) &                                                
          -5.0D0*C5OHR2/R2**5*D5R2**4*(D5R2/R2+(1.0D0-D5R2)*(-ADAMP(5)- &                           
             2.0D0*BDAMP(5)*RR2)*DRR2)                                                              
      Ddisp3_f4=-4.0D0*C4OHR3/R3**4*D4R3**3*(D4R3/R3+(1.0D0-D4R3)* &                                
             (-ADAMP(4)-2.0D0*BDAMP(4)*RR3)*DRR3) &                                                 
           -5.0D0*C5OHR3/R3**5*D5R3**4*(D5R3/R3+(1.0D0-D5R3)*(-ADAMP(5)- &                          
             2.0D0*BDAMP(5)*RR3)*DRR3)                                                              
      DER1=Ddisp1_f4+DC42R1*T4R2+DC43R1*T4R3+DC52R1*T5R2+DC53R1*T5R3 &                              
           +DC41R1*T4R1+DC51R1*T5R1                                                                 
      DER2=Ddisp2_f4+DC41R2*T4R1+DC43R2*T4R3+DC51R2*T5R1+DC53R2*T5R3 &                              
           +DC42R2*T4R2+DC52R2*T5R2                                                                 
      DER3=Ddisp3_f4+DC42R3*T4R2+DC41R3*T4R1+DC52R3*T5R2+DC51R3*T5R1 &                              
           +DC43R3*T4R3+DC53R3*T5R3                                                                 
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
      FUNCTION Ddisp_f4(R,C6,C8,C10,R0,RM)                                                          
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                                            
      COMMON/DAMPC_f4/ADAMP(10),BDAMP(10)                                                           
      SAVE /DAMPC_f4/                                                                               
                                                                                                    
      R6=R**6                                                                                       
      R8=R6*R*R                                                                                     
      R10=R8*R*R                                                                                    
      RR=2.0D0*R/(RM+2.5D0*R0)                                                                      
      DRR=RR/R                                                                                      
      T6=1.0D0-EXP(-ADAMP(6)*RR-BDAMP(6)*RR**2)                                                     
      T8=1.0D0-EXP(-ADAMP(8)*RR-BDAMP(8)*RR**2)                                                     
      T10=1.0D0-EXP(-ADAMP(10)*RR-BDAMP(10)*RR**2)                                                  
      DDISP_f4=6.0D0*C6/R6*T6**5*(T6/R+(1.0D0-T6)*(-ADAMP(6)- &                                     
         2.0D0*BDAMP(6)*RR)*DRR)+ &                                                                 
            8.0D0*C8/R8*T8**7*(T8/R+(1.0D0-T8)*(-ADAMP(8)-2.0D0*BDAMP(8) &                          
            *RR)*DRR)+ &                                                                            
            10.0D0*C10/R10*T10**9*(T10/R+(1.0D0-T10)*(-ADAMP(10)-2.0D0* &                           
            BDAMP(10)*RR)*DRR)                                                                      
                                                                                                    
      END                                                                                           
                                                                                                    
                                                                                                    
                                                                                                    
