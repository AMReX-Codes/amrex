#include <CONSTANTS.H>
#include <REAL.H>
#include <DIMS.H>

#include "main_F.H"


#define TEST_EXT 1
#define POWER 3
#define MUPOWER 3


      subroutine FORT_FILLCRSE (
     $     crse, DIMS(crse),
     $     h, nc
     $     )
      integer nc
      integer DIMDEC(crse)
      REAL_T crse(DIMV(crse),nc)
      REAL_T h(BL_SPACEDIM)
c
      integer  i, j, k, n
      REAL_T x, y, z

      REAL_T cux,cuy,cuz
      REAL_T cvx,cvy,cvz
      REAL_T cwx,cwy,cwz
      REAL_T alpha,beta,gamma
      namelist /fortin/ cux,cuy,cuz,
     $     cvx,cvy,cvz,cwx,cwy,cwz,alpha,beta,gamma
c      
      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(9)

      do k=ARG_L3(crse), ARG_H3(crse)
         z = (k+0.5)*h(3)
         do j=ARG_L2(crse), ARG_H2(crse)
            y = (j + 0.5)*h(2)
            do i=ARG_L1(crse), ARG_H1(crse)
               x = (i + 0.5)*h(1)
               crse(i,j,k,1)=cux*x**POWER+cuy*y**POWER+cuz*z**POWER
               crse(i,j,k,2)=cvx*x**POWER+cvy*y**POWER+cvz*z**POWER
               crse(i,j,k,3)=cwx*x**POWER+cwy*y**POWER+cwz*z**POWER
            enddo
         enddo
      enddo
      
      end
      

      
      subroutine FORT_FILLFINE (
     $     fine, DIMS(fine),
     $     h, nc
     $     )
      integer nc
      integer DIMDEC(fine)
      REAL_T fine(DIMV(fine),nc)
      REAL_T h(BL_SPACEDIM)
c
      integer  i, j, k, n
      REAL_T x, y,z
c
      REAL_T cux,cuy,cuz
      REAL_T cvx,cvy,cvz
      REAL_T cwx,cwy,cwz
      REAL_T alpha,beta,gamma
      namelist /fortin/ cux,cuy,cuz,
     $     cvx,cvy,cvz,cwx,cwy,cwz,alpha,beta,gamma
c      
      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(9)
c
      do k=ARG_L3(fine), ARG_H3(fine)
         z = (k + 0.5)*h(3)
#if TEST_EXT
         if( k.lt.0) z = 0
#endif
         do j=ARG_L2(fine), ARG_H2(fine)
            y = (j + 0.5)*h(2)
#if TEST_EXT
            if( j.lt.0) y = 0
#endif
            do i=ARG_L1(fine), ARG_H1(fine)
               x = (i + 0.5)*h(1)
#if TEST_EXT
               if( i.lt.0) x = 0
#endif
               fine(i,j,k,1)=cux*x**POWER+cuy*y**POWER+cuz*z**POWER
               fine(i,j,k,2)=cvx*x**POWER+cvy*y**POWER+cvz*z**POWER
               fine(i,j,k,3)=cwx*x**POWER+cwy*y**POWER+cwz*z**POWER
            enddo
         enddo
      enddo


      end
      
      subroutine FORT_FILLRHS (
     $     rhs, DIMS(rhs),
     $     h, nc
     $     )
      integer nc
      integer DIMDEC(rhs)
      REAL_T rhs(DIMV(rhs),nc)
      REAL_T h(BL_SPACEDIM)
c
      integer  i, j, k, n
      REAL_T x, y,z
c
      REAL_T cux,cuy,cuz
      REAL_T cvx,cvy,cvz
      REAL_T cwx,cwy,cwz
      REAL_T alpha,beta,gamma
      namelist /fortin/ cux,cuy,cuz,
     $     cvx,cvy,cvz,cwx,cwy,cwz,alpha,beta,gamma
c      
      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(9)
c
      do k=ARG_L3(rhs), ARG_H3(rhs)
         z = (k + 0.5)*h(3)
         do j=ARG_L2(rhs), ARG_H2(rhs)
            y = (j + 0.5)*h(2)
            do i=ARG_L1(rhs), ARG_H1(rhs)
               x = (i + 0.5)*h(1)
#if (POWER==3)&(MUPOWER==3)
               rhs(i,j,k,1)= -(
     -  18*alpha*cux*x**4 + 3*beta*y**2*(3*cvx*x**2 + 3*cuy*y**2) + 
     -  3*gamma*z**2*(3*cwx*x**2 + 3*cuz*z**2) + 
     -  12*cux*x*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  6*cuy*y*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  6*cuz*z*(1 + alpha*x**3 + beta*y**3 + gamma*z**3)
     &              )
               rhs(i,j,k,2)= -(
     -  18*beta*cvy*y**4 + 3*alpha*x**2*(3*cvx*x**2 + 3*cuy*y**2) + 
     -  3*gamma*z**2*(3*cwy*y**2 + 3*cvz*z**2) + 
     -  6*cvx*x*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  12*cvy*y*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  6*cvz*z*(1 + alpha*x**3 + beta*y**3 + gamma*z**3)
     &              )
               rhs(i,j,k,3)= -(
     -  18*cwz*gamma*z**4 + 3*alpha*x**2*(3*cwx*x**2 + 3*cuz*z**2) + 
     -  3*beta*y**2*(3*cwy*y**2 + 3*cvz*z**2) + 
     -  6*cwx*x*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  6*cwy*y*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  12*cwz*z*(1 + alpha*x**3 + beta*y**3 + gamma*z**3)
     &              )
#endif
#if (POWER==2)&(MUPOWER==3)
               rhs(i,j,k,1)= -(
     -   12*alpha*cux*x**3 + 3*beta*y**2*(2*cvx*x + 2*cuy*y) + 
     -  3*gamma*z**2*(2*cwx*x + 2*cuz*z) + 
     -  4*cux*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  2*cuy*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  2*cuz*(1 + alpha*x**3 + beta*y**3 + gamma*z**3)
     &              )
               rhs(i,j,k,2)= -(
     -  12*beta*cvy*y**3 + 3*alpha*x**2*(2*cvx*x + 2*cuy*y) + 
     -  3*gamma*z**2*(2*cwy*y + 2*cvz*z) + 
     -  2*cvx*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  4*cvy*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  2*cvz*(1 + alpha*x**3 + beta*y**3 + gamma*z**3)
     &              )
               rhs(i,j,k,3)= -(
     -  12*cwz*gamma*z**3 + 3*alpha*x**2*(2*cwx*x + 2*cuz*z) + 
     -  3*beta*y**2*(2*cwy*y + 2*cvz*z) + 
     -  2*cwx*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  2*cwy*(1 + alpha*x**3 + beta*y**3 + gamma*z**3) + 
     -  4*cwz*(1 + alpha*x**3 + beta*y**3 + gamma*z**3)
     &              )
#endif
            enddo
         enddo
      enddo


      end
      

      subroutine FORT_MAKEMU(
     $     mu,DIMS(mu),
     $     h, idir )

      integer idir
      integer DIMDEC(mu)
      REAL_T mu(DIMV(mu))
      REAL_T h(BL_SPACEDIM)

      integer i,j,k
      REAL_T x,y,z
c
      REAL_T cux,cuy,cuz
      REAL_T cvx,cvy,cvz
      REAL_T cwx,cwy,cwz
      REAL_T alpha,beta,gamma
      namelist /fortin/ cux,cuy,cuz,
     $     cvx,cvy,cvz,cwx,cwy,cwz,alpha,beta,gamma
c      
      open(9,file='probin',form='formatted',status='old')
      read(9,fortin)
      close(9)
c

      if( idir .eq. 0 ) then
         do k=ARG_L3(mu),ARG_H3(mu)
            z = (k+0.5)*h(3)
            do j=ARG_L2(mu),ARG_H2(mu)
               y = (j+0.5)*h(2)
               do i=ARG_L1(mu),ARG_H1(mu)
                  x = (i)*h(1)
                  mu(i,j,k) = 1+alpha*x**MUPOWER+beta*y**MUPOWER+gamma*z**MUPOWER
               enddo
            enddo
         enddo
      else if( idir .eq. 1 ) then
         do k=ARG_L3(mu),ARG_H3(mu)
            z = (k+0.5)*h(3)
            do j=ARG_L2(mu),ARG_H2(mu)
               y = (j)*h(2)
               do i=ARG_L1(mu),ARG_H1(mu)
                  x = (i+0.5)*h(1)
                  mu(i,j,k) = 1+alpha*x**MUPOWER+beta*y**MUPOWER+gamma*z**MUPOWER
               enddo
            enddo
         enddo
      else if( idir .eq. 2 ) then
         do k=ARG_L3(mu),ARG_H3(mu)
            z = (k)*h(3)
            do j=ARG_L2(mu),ARG_H2(mu)
               y = (j+0.5)*h(2)
               do i=ARG_L1(mu),ARG_H1(mu)
                  x = (i+0.5)*h(1)
                  mu(i,j,k) = 1+alpha*x**MUPOWER+beta*y**MUPOWER+gamma*z**MUPOWER
               enddo
            enddo
         enddo
      else
         write(6,*)'FORT_MAKEMU: bad idir'
         stop
      endif

      return
      end
