#include <CONSTANTS.H>
#include <REAL.H>
#include <DIMS.H>

#include "main_F.H"

#if 0
c .. 1 is 0
c .. 2 is x
c .. 3 is xy
#define DOGU 1

c .. 1 is y**2
c .. 2 is y
c .. 3 is (24-y)**2
c .. 4 is 0
c .. 5 is x**2
#define DOGV 1

      subroutine FORT_FILLCRSE (
     $     crse, DIMS(crse),
     $     h, nc
     $     )
      integer nc
      integer DIMDEC(crse)
      REAL_T crse(DIMV(crse),nc)
      REAL_T h(BL_SPACEDIM)
c
      integer  i, j, n
      REAL_T x, y
c      
      do n=1,nc
         if( n .eq. 1 ) then
            do j=ARG_L2(crse), ARG_H2(crse)
               y = (j + 0.5)*h(2)
               do i=ARG_L1(crse), ARG_H1(crse)
                  x = (i + 0.5)*h(1)
#if DOGU==1
                  crse(i,j,n)=0
#elif DOGU==2
                  crse(i,j,n)=x
#elif DOGU==3
                  crse(i,j,n)=x*y
#else
                  error is error
#endif
               enddo
            enddo
         else
            do j=ARG_L2(crse), ARG_H2(crse)
               y = (j + 0.5)*h(2)
               do i=ARG_L1(crse), ARG_H1(crse)
                  x = (i + 0.5)*h(1)
#if DOGV==1
                  crse(i,j,n)=y*y
#elif DOGV==2
                  crse(i,j,n)=y
#elif DOGV==3
                  crse(i,j,n)=(24-y)*(24-y)
#elif DOGV==4
                  crse(i,j,n)=0
#elif DOGV==5
                  crse(i,j,n)=x*x
#else
                  error is error
#endif
               enddo
            enddo
         endif
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
      integer  i, j, n
      REAL_T x, y
c
      do n=1,nc
         if( n .eq. 1 ) then
            do j=ARG_L2(fine), ARG_H2(fine)
               y = (j + 0.5)*h(2)
               do i=ARG_L1(fine), ARG_H1(fine)
                  x = (i + 0.5)*h(1)
                  if( x .gt. 0 ) then
#if DOGU==1
                  fine(i,j,n)=0
#elif DOGU==2
                  fine(i,j,n)=x
#elif DOGU==3
                  fine(i,j,n)=x*y
#else
                  error is error
#endif
                  else
c....................at physical boundary, use homog bc
                     fine(i,j,n) = 0
                  endif
               enddo
            enddo
         else
            do j=ARG_L2(fine), ARG_H2(fine)
               y = (j + 0.5)*h(2)
               do i=ARG_L1(fine), ARG_H1(fine)
                  x = (i + 0.5)*h(1)
                  if( x .gt. 0 ) then
#if DOGV==1
                  fine(i,j,n)=y*y
#elif DOGV==2
                  fine(i,j,n)=y
#elif DOGV==3
                  fine(i,j,n)=(24-y)*(24-y)
#elif DOGV==4
                  fine(i,j,n)=0
#elif DOGV==5
                  fine(i,j,n)=x*x
#else
                  error is error
#endif
                  else
                     fine(i,j,n) = 0
                  endif
               enddo
            enddo
         endif
      enddo
      end
#else
      
#define TEST_EXT 0
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
      integer  i, j, n
      REAL_T x, y

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

      do j=ARG_L2(crse), ARG_H2(crse)
         y = (j + 0.5)*h(2)
         do i=ARG_L1(crse), ARG_H1(crse)
            x = (i + 0.5)*h(1)
            crse(i,j,1)=cux*x**POWER+cuy*y**POWER
            crse(i,j,2)=cvx*x**POWER+cvy*y**POWER
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
      integer  i, j, n
      REAL_T x, y
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
            fine(i,j,1)=cux*x**POWER+cuy*y**POWER
            fine(i,j,2)=cvx*x**POWER+cvy*y**POWER
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
      integer  i, j, n
      REAL_T x, y
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
      do j=ARG_L2(rhs), ARG_H2(rhs)
         y = (j + 0.5)*h(2)
         do i=ARG_L1(rhs), ARG_H1(rhs)
            x = (i + 0.5)*h(1)
#if (POWER==3)&(MUPOWER==3)
            rhs(i,j,1)= -(
     -           18*alpha*cux*x**4 + 3*beta*y**2*(3*cvx*x**2 + 3*cuy*y**2) + 
     -           12*cux*x*(1 + alpha*x**3 + beta*y**3) + 
     -           6*cuy*y*(1 + alpha*x**3 + beta*y**3)
     &           )
            rhs(i,j,2)= -(
     -           18*beta*cvy*y**4 + 3*alpha*x**2*(3*cvx*x**2 + 3*cuy*y**2) + 
     -           6*cvx*x*(1 + alpha*x**3 + beta*y**3) + 
     -           12*cvy*y*(1 + alpha*x**3 + beta*y**3) 
     &           )
#endif
#if (POWER==2)&(MUPOWER==3)
            rhs(i,j,1)= -(
     -           12*alpha*cux*x**3 + 3*beta*y**2*(2*cvx*x + 2*cuy*y) + 
     -           4*cux*(1 + alpha*x**3 + beta*y**3) + 
     -           2*cuy*(1 + alpha*x**3 + beta*y**3)
     &           )
            rhs(i,j,2)= -(
     -           12*beta*cvy*y**3 + 3*alpha*x**2*(2*cvx*x + 2*cuy*y) + 
     -           2*cvx*(1 + alpha*x**3 + beta*y**3) + 
     -           4*cvy*(1 + alpha*x**3 + beta*y**3)
     &           )
#endif
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

      integer i,j
      REAL_T x,y
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
         do j=ARG_L2(mu),ARG_H2(mu)
            y = (j+0.5)*h(2)
            do i=ARG_L1(mu),ARG_H1(mu)
               x = (i)*h(1)
               mu(i,j) = 1+alpha*x**MUPOWER+beta*y**MUPOWER
            enddo
         enddo
      else if( idir .eq. 1 ) then
         do j=ARG_L2(mu),ARG_H2(mu)
            y = (j)*h(2)
            do i=ARG_L1(mu),ARG_H1(mu)
               x = (i+0.5)*h(1)
               mu(i,j) = 1+alpha*x**MUPOWER+beta*y**MUPOWER
            enddo
         enddo
      else
         write(6,*)'FORT_MAKEMU: bad idir'
         stop
      endif
      
      return
      end
      
#endif
