 real*8 function f1_lr(r)
 use variables
 use phys_parameters
 implicit none
 real*8 r,rin,rou,pi
 data pi/3.1415926535897931d0/
 rin=rl_cut_off-0.2d-9
 rou=rl_cut_off
 if (r.le.rin) then 
    f1_lr=0.0d0
 else if (r.ge.rou) then 
    f1_lr=1.0d0
 else 
    f1_lr=0.5d0*(1.0d0+sin(pi*((r-rin)/(rou-rin)-0.5d0)))
 end if
 return
 end

 real*8 function df1_lr(r)
 use variables
 use phys_parameters

 implicit none
 real*8 r,rin,rou,pi
 data pi/3.1415926535897931d0/

 rin=rl_cut_off-0.2d-9
 rou=rl_cut_off

 if (r.le.rin) then 
    df1_lr=0.0d0
 else if (r.ge.rou) then 
    df1_lr=0.0d0
 else 
    df1_lr=0.5d0*(cos(pi*((r-rin)/(rou-rin)-0.5d0)))*pi/(rou-rin)
 end if
 return
 end




real*8 function f1(r)
 use variables
 use phys_parameters
 implicit none
 real*8 r,rin,rou,pi
 data pi/3.1415926535897931d0/
 rin=in_bound
 rou=out_bound
 if (r.le.rin) then 
    f1=0.0d0
 else if (r.ge.rou) then 
    f1=1.0d0
 else 
    f1=0.5d0*(1.0d0+sin(pi*((r-rin)/(rou-rin)-0.5d0)))
 end if
 return
 end

 real*8 function df1(r)
 use variables
 use phys_parameters

 implicit none
 real*8 r,rin,rou,pi
 data pi/3.1415926535897931d0/

 rin=in_bound
 rou=out_bound

 if (r.le.rin) then 
    df1=0.0d0
 else if (r.ge.rou) then 
    df1=0.0d0
 else 
    df1=0.5d0*(cos(pi*((r-rin)/(rou-rin)-0.5d0)))*pi/(rou-rin)
 end if
 return
 end


          
              
      real*8 function f2(r1,r2)             
      implicit none                         
      real*8 r1,r2,f1,fc1,fc2
       fc1=f1(r1)
       fc2=f1(r2)
       if ((fc1.lt.1.0d-13).or.(fc2.lt.1.0d-13)) then
          f2=0.0d0
          return
      end if
                       
!      f2=2.0d0/(1.0d0/f1(r1)+1.0d0/f1(r2))  
      f2=2.0d0*fc1*fc2/(fc1+fc2)
      return    
      end       
                
      subroutine df2(r1,r2,dr1,dr2)         
      implicit none                         
      real*8 r1,r2,dr1,dr2,f1,df1,fc1,fc2           
      fc1=f1(r1)
      fc2=f1(r2) 
       if ((fc1.lt.1.0d-13).or.(fc2.lt.1.0d-13)) then
          dr1=0.0d0
          dr2=0.0d0
          return
      end if
      dr1=2.0d0*df1(r1)*(fc2/(fc1+fc2))**2                     
      dr2=2.0d0*df1(r2)*(fc1/(fc1+fc2))**2                     

!      dr1=2.0d0*df1(r1)*(f1(r2)/(f1(r1)+f1(r2)))**2 
!      dr2=2.0d0*df1(r2)*(f1(r1)/(f1(r1)+f1(r2)))**2 


      return    
      end       
                
      real*8 function f3(r1,r2,r3)          
      implicit none                         
      real*8 r1,r2,r3,f1,fc1,fc2,fc3                    
      fc1=f1(r1)
      fc2=f1(r2) 
      fc3=f1(r3) 
       if ((fc1.lt.1.0d-13).or.(fc2.lt.1.0d-13).or.(fc3.lt.1.0d-13)) then
          f3=0.0d0
          return
      end if

!      f3=3.0d0/(1.0d0/f1(r1)+1.0d0/f1(r2)+1.0d0/f1(r3))   
      f3=3.0d0*fc1*fc2*fc3/(fc1*fc2+fc1*fc3+fc2*fc3)             


      return    
      end       
                
      subroutine df3(r1,r2,r3,dr1,dr2,dr3)  
      implicit none                         
      real*8 r1,r2,r3,dr1,dr2,dr3,f1,df1,u,fc1,fc2,fc3      
      fc1=f1(r1)
      fc2=f1(r2) 
      fc3=f1(r3) 
       if ((fc1.lt.1.0d-13).or.(fc2.lt.1.0d-13).or.(fc3.lt.1.0d-13)) then
          dr1=0.0d0
          dr2=0.0d0
          dr3=0.0d0
          return
      end if

!      u2=(f1(r1)*f1(r2)+f1(r1)*f1(r3)+f1(r2)*f1(r3))**2 
!      dr1=3.0d0*df1(r1)*(f1(r2)*f1(r3))**2/u2                     
!      dr2=3.0d0*df1(r2)*(f1(r1)*f1(r3))**2/u2        
!      dr3=3.0d0*df1(r3)*(f1(r1)*f1(r2))**2/u2  
                  
      u=fc1*fc2+fc1*fc3+fc2*fc3
      dr1=3.0d0*df1(r1)*(fc2*fc3/u)**2                     
      dr2=3.0d0*df1(r2)*(fc1*fc3/u)**2        
      dr3=3.0d0*df1(r3)*(fc1*fc2/u)**2                    

      return    
      end       
               
      real*8 function f4(r1,r2,r3,r4)       
      implicit none                         
      real*8 r1,r2,r3,r4,f1,fc1,fc2,fc3,fc4                   
!      f4=4.0d0/(1.0d0/f1(r1)+1.0d0/f1(r2)+1.0d0/f1(r3)+1.0d0/f1(r4)) 
      fc1=f1(r1)
      fc2=f1(r2) 
      fc3=f1(r3) 
      fc4=f1(r4) 
       if ((fc1.lt.1.0d-13).or.(fc2.lt.1.0d-13).or.(fc3.lt.1.0d-13).or.(fc4.lt.1.0d-13)) then
          f4=0.0d0
          return
      end if
   
      f4=4.0d0*fc1*fc2*fc3*fc4/(fc2*fc3*fc4+fc1*fc3*fc4+fc1*fc2*fc4+fc1*fc2*fc3) 
      return    
      end       
                
      subroutine df4(r1,r2,r3,r4,dr1,dr2,dr3,dr4)                       
      implicit none                         
      real*8 r1,r2,r3,r4,dr1,dr2,dr3,dr4,f1,df1,u,fc1,fc2,fc3,fc4                         
  
      fc1=f1(r1)
      fc2=f1(r2) 
      fc3=f1(r3) 
      fc4=f1(r4) 
       if ((fc1.lt.1.0d-13).or.(fc2.lt.1.0d-13).or.(fc3.lt.1.0d-13).or.(fc3.lt.1.0d-13)) then
          dr1=0.0d0
          dr2=0.0d0
          dr3=0.0d0
          dr4=0.0d0
          return
      end if

      u=fc2*fc3*fc4+fc1*fc3*fc4+fc1*fc2*fc4+fc1*fc2*fc3

!      dr1=4.0d0*df1(r1)*(f1(r2)*f1(r3)*f1(r4)/(f1(r2)*f1(r3)*f1(r4)+f1(r1)*f1(r3)*f1(r4)+f1(r1)*f1(r2)*f1(r4)+f1(r1)*f1(r2)*f1(r3)))**2
!      dr2=4.0d0*df1(r2)*(f1(r1)*f1(r3)*f1(r4)/(f1(r2)*f1(r3)*f1(r4)+f1(r1)*f1(r3)*f1(r4)+f1(r1)*f1(r2)*f1(r4)+f1(r1)*f1(r2)*f1(r3)))**2
!      dr3=4.0d0*df1(r3)*(f1(r1)*f1(r2)*f1(r4)/(f1(r2)*f1(r3)*f1(r4)+f1(r1)*f1(r3)*f1(r4)+f1(r1)*f1(r2)*f1(r4)+f1(r1)*f1(r2)*f1(r3)))**2
!      dr4=4.0d0*df1(r4)*(f1(r1)*f1(r2)*f1(r3)/(f1(r2)*f1(r3)*f1(r4)+f1(r1)*f1(r3)*f1(r4)+f1(r1)*f1(r2)*f1(r4)+f1(r1)*f1(r2)*f1(r3)))**2
      dr1=4.0d0*df1(r1)*(fc2*fc3*fc4/u)**2
      dr2=4.0d0*df1(r2)*(fc1*fc3*fc4/u)**2
      dr3=4.0d0*df1(r3)*(fc1*fc2*fc4/u)**2
      dr4=4.0d0*df1(r4)*(fc1*fc2*fc3/u)**2



      return    
      end       
 
