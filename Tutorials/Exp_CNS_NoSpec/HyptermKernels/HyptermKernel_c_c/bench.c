#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#include <omp.h>
#ifdef _MPI
#include <mpi.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#define irho 0
#define imx  1
#define imy  2
#define imz  3
#define iene 4

#define qu    1
#define qv    2
#define qw    3
#define qpres 4

#define ALP  0.8
#define BET -0.2
#define GAM  0.038095238095238 //  4.0/105.0
#define DEL -0.003571428571429 // -1.0/280.0

#define JBlockSize 11

//------------------------------------------------------------------------------------------------------------------------------
#include "timer.h"
//------------------------------------------------------------------------------------------------------------------------------
uint64_t _total_run_time;
uint64_t _total_time_hypterm;
uint64_t _total_time_hypterm_L1;
uint64_t _total_time_hypterm_L2;
uint64_t _total_time_hypterm_L3;

//------------------------------------------------------------------------------------------------------------------------------
void init_data(int lo[3], int hi[3], int ng, double dx[3], double ** __restrict__ _cons){

  int c;
  int lo0=lo[0];int hi0=hi[0];
  int lo1=lo[1];int hi1=hi[1];
  int lo2=lo[2];int hi2=hi[2];
  int pencilg =  (hi0-lo0)+2*ng;
  int  planeg = ((hi1-lo1)+2*ng)*pencilg;

  double * __restrict__ cons[5];for(c=0;c<5;c++){cons[c] = _cons[c] + (ng+ng*pencilg+ng*planeg);}

  int i, j, k;
  double scale, xloc, yloc, zloc, rholoc, uvel, vvel, wvel, eloc;

  printf("lo = %d %d %d\n", lo[0], lo[1], lo[2]);
  printf("hi = %d %d %d\n", hi[0], hi[1], hi[2]);

  scale = 3.14e0;

  #pragma omp parallel for private(i,j,k) reduction(+ : L1iters)
  for(k=lo2;k<=hi2;k++){
         zloc = k*dx[3]/scale;
    for(j=lo1;j<=hi1;j++){
            yloc = j*dx[2]/scale;
      for(i=lo0;i<=hi0;i++){
               xloc = i*dx[1]/scale;

               uvel   = 1.1e4*sin(1*xloc)*sin(2*yloc)*sin(3*zloc);
               vvel   = 1.0e4*sin(2*xloc)*sin(4*yloc)*sin(1*zloc);
               wvel   = 1.2e4*sin(3*xloc)*cos(2*yloc)*sin(2*zloc);
               rholoc = 1.0e-3 + 1.0e-5*sin(1*xloc)*cos(2*yloc)*cos(3*zloc);
               eloc   = 2.5e9  + 1.0e-3*sin(2*xloc)*cos(2*yloc)*sin(2*zloc);

        int ijkg = i + j*pencilg + k*planeg;
        cons[irho][ijkg] = rholoc;
        cons[irho][ijkg] = rholoc*uvel;
        cons[irho][ijkg] = rholoc*vvel;
        cons[irho][ijkg] = rholoc*wvel;
        cons[irho][ijkg] = rholoc*(eloc + (uvel*uvel+vvel*vvel+wvel*wvel)/2.0);
      }
    }
  }
}


//------------------------------------------------------------------------------------------------------------------------------
void hypterm_naive(int lo[3], int hi[3], int ng, double dx[3], double ** __restrict__ _cons, double ** __restrict__ _q, double ** __restrict__ _flux){ // (lo,hi,ng,dx,cons,q,flux)
  //    double precision, intent(in ) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,5)
  //    double precision, intent(in ) ::    q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,6)
  //    double precision, intent(out) :: flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)
  uint64_t _time_start      = CycleTime();
  int c;
  double dxinv0=1.0/dx[0];
  double dxinv1=1.0/dx[1];
  double dxinv2=1.0/dx[2];
  int lo0=lo[0];int hi0=hi[0];
  int lo1=lo[1];int hi1=hi[1];
  int lo2=lo[2];int hi2=hi[2];

  int pencil  =  (hi0-lo0);
  int  plane  = ((hi1-lo1))*pencil;
  int pencilg =  (hi0-lo0)+2*ng;
  int  planeg = ((hi1-lo1)+2*ng)*pencilg;

  double * __restrict__ cons[5];for(c=0;c<5;c++){cons[c] = _cons[c] + (ng+ng*pencilg+ng*planeg);}
  double * __restrict__    q[6];for(c=0;c<6;c++){   q[c] =    _q[c] + (ng+ng*pencilg+ng*planeg);}
  double * __restrict__ flux[5];for(c=0;c<5;c++){flux[c] = _flux[c];}

  int i,j,k;
  int jb,JBlocks;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  int L1iters = 0;
  printf("lo0 hi0 = %d %d\n", lo0, hi0);
  printf("lo1 hi1 = %d %d\n", lo1, hi1);
  printf("lo2 hi2 = %d %d\n", lo2, hi2);
  #pragma omp parallel for private(i,j,k) reduction(+ : L1iters)
  for(k=lo2;k<=hi2;k++){
    for(j=lo1;j<=hi1;j++){
      for(i=lo0;i<=hi0;i++){
	++L1iters;
        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;

        double qum4 = q[qu][ijkg-4];
        double qum3 = q[qu][ijkg-3];
        double qum2 = q[qu][ijkg-2];
        double qum1 = q[qu][ijkg-1];
        double qup1 = q[qu][ijkg+1];
        double qup2 = q[qu][ijkg+2];
        double qup3 = q[qu][ijkg+3];
        double qup4 = q[qu][ijkg+4];

        flux[irho][ijk] = 
             - dxinv0*( ALP*( cons[ imx][ijkg+1] - cons[ imx][ijkg-1] ) +
                        BET*( cons[ imx][ijkg+2] - cons[ imx][ijkg-2] ) +
                        GAM*( cons[ imx][ijkg+3] - cons[ imx][ijkg-3] ) +
                        DEL*( cons[ imx][ijkg+4] - cons[ imx][ijkg-4] ) );

        flux[imx][ijk] = 
             - dxinv0*( ALP*( cons[ imx][ijkg+1]*qup1 - cons[ imx][ijkg-1]*qum1 + q[qpres][ijkg+1] - q[qpres][ijkg-1] ) +
                        BET*( cons[ imx][ijkg+2]*qup2 - cons[ imx][ijkg-2]*qum2 + q[qpres][ijkg+2] - q[qpres][ijkg-2] ) +
                        GAM*( cons[ imx][ijkg+3]*qup3 - cons[ imx][ijkg-3]*qum3 + q[qpres][ijkg+3] - q[qpres][ijkg-3] ) +
                        DEL*( cons[ imx][ijkg+4]*qup4 - cons[ imx][ijkg-4]*qum4 + q[qpres][ijkg+4] - q[qpres][ijkg-4] ) );

        flux[imy][ijk] = 
             - dxinv0*( ALP*( cons[ imy][ijkg+1]*qup1 - cons[ imy][ijkg-1]*qum1 ) +
                        BET*( cons[ imy][ijkg+2]*qup2 - cons[ imy][ijkg-2]*qum2 ) +
                        GAM*( cons[ imy][ijkg+3]*qup3 - cons[ imy][ijkg-3]*qum3 ) +
                        DEL*( cons[ imy][ijkg+4]*qup4 - cons[ imy][ijkg-4]*qum4 ) );

        flux[imz][ijk] = 
             - dxinv0*( ALP*( cons[ imz][ijkg+1]*qup1 - cons[ imz][ijkg-1]*qum1 ) +
                        BET*( cons[ imz][ijkg+2]*qup2 - cons[ imz][ijkg-2]*qum2 ) +
                        GAM*( cons[ imz][ijkg+3]*qup3 - cons[ imz][ijkg-3]*qum3 ) +
                        DEL*( cons[ imz][ijkg+4]*qup4 - cons[ imz][ijkg-4]*qum4 ) );

        flux[iene][ijk] = 
             - dxinv0*( ALP*( cons[iene][ijkg+1]*qup1 - cons[iene][ijkg-1]*qum1 + q[qpres][ijkg+1]*qup1 - q[qpres][ijkg-1]*qum1 ) +
                        BET*( cons[iene][ijkg+2]*qup2 - cons[iene][ijkg-2]*qum2 + q[qpres][ijkg+2]*qup2 - q[qpres][ijkg-2]*qum2 ) +
                        GAM*( cons[iene][ijkg+3]*qup3 - cons[iene][ijkg-3]*qum3 + q[qpres][ijkg+3]*qup3 - q[qpres][ijkg-3]*qum3 ) +
                        DEL*( cons[iene][ijkg+4]*qup4 - cons[iene][ijkg-4]*qum4 + q[qpres][ijkg+4]*qup4 - q[qpres][ijkg-4]*qum4 ) );
  }}}
  uint64_t _time_L1 = CycleTime();
  _total_time_hypterm_L1 += (_time_L1 - _time_start);
  printf("L1 iters = %d\n",L1iters    );
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  int pencilg2 = pencilg*2;
  int pencilg3 = pencilg*3;
  int pencilg4 = pencilg*4;
  #pragma omp parallel for private(i,j,k)
  for(k=lo2;k<=hi2;k++){
    for(j=lo1;j<=hi1;j++){
      for(i=lo0;i<=hi0;i++){

        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;

        double qvp1 = q[qv][ijkg+pencilg  ];
        double qvp2 = q[qv][ijkg+pencilg2 ];
        double qvp3 = q[qv][ijkg+pencilg3 ];
        double qvp4 = q[qv][ijkg+pencilg4 ];
        double qvm1 = q[qv][ijkg-pencilg  ];
        double qvm2 = q[qv][ijkg-pencilg2 ];
        double qvm3 = q[qv][ijkg-pencilg3 ];
        double qvm4 = q[qv][ijkg-pencilg4 ];

        flux[irho][ijk] = flux[irho][ijk]
             - dxinv1*( ALP*( cons[ imy][ijkg+pencilg  ] - cons[ imy][ijkg-pencilg  ] ) +
                        BET*( cons[ imy][ijkg+pencilg2 ] - cons[ imy][ijkg-pencilg2 ] ) +
                        GAM*( cons[ imy][ijkg+pencilg3 ] - cons[ imy][ijkg-pencilg3 ] ) +
                        DEL*( cons[ imy][ijkg+pencilg4 ] - cons[ imy][ijkg-pencilg4 ] ) );

      }
      for(i=lo0;i<=hi0;i++){

        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;

        double qvp1 = q[qv][ijkg+pencilg  ];
        double qvp2 = q[qv][ijkg+pencilg2 ];
        double qvp3 = q[qv][ijkg+pencilg3 ];
        double qvp4 = q[qv][ijkg+pencilg4 ];
        double qvm1 = q[qv][ijkg-pencilg  ];
        double qvm2 = q[qv][ijkg-pencilg2 ];
        double qvm3 = q[qv][ijkg-pencilg3 ];
        double qvm4 = q[qv][ijkg-pencilg4 ];

        flux[imx][ijk] = flux[imx][ijk]
             - dxinv1*( ALP*( cons[ imx][ijkg+pencilg  ]*qvp1 - cons[ imx][ijkg-pencilg  ]*qvm1 ) +
                        BET*( cons[ imx][ijkg+pencilg2 ]*qvp2 - cons[ imx][ijkg-pencilg2 ]*qvm2 ) +
                        GAM*( cons[ imx][ijkg+pencilg3 ]*qvp3 - cons[ imx][ijkg-pencilg3 ]*qvm3 ) +
                        DEL*( cons[ imx][ijkg+pencilg4 ]*qvp4 - cons[ imx][ijkg-pencilg4 ]*qvm4 ) );

      }
      for(i=lo0;i<=hi0;i++){

        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;

        double qvp1 = q[qv][ijkg+pencilg  ];
        double qvp2 = q[qv][ijkg+pencilg2 ];
        double qvp3 = q[qv][ijkg+pencilg3 ];
        double qvp4 = q[qv][ijkg+pencilg4 ];
        double qvm1 = q[qv][ijkg-pencilg  ];
        double qvm2 = q[qv][ijkg-pencilg2 ];
        double qvm3 = q[qv][ijkg-pencilg3 ];
        double qvm4 = q[qv][ijkg-pencilg4 ];

        flux[imy][ijk] = flux[imy][ijk]
             - dxinv1*( ALP*( cons[ imy][ijkg+pencilg  ]*qvp1 - cons[ imy][ijkg-pencilg  ]*qvm1 + q[qpres][ijkg+pencilg  ] - q[qpres][ijkg-pencilg  ] ) +
                        BET*( cons[ imy][ijkg+pencilg2 ]*qvp2 - cons[ imy][ijkg-pencilg2 ]*qvm2 + q[qpres][ijkg+pencilg2 ] - q[qpres][ijkg-pencilg2 ] ) +
                        GAM*( cons[ imy][ijkg+pencilg3 ]*qvp3 - cons[ imy][ijkg-pencilg3 ]*qvm3 + q[qpres][ijkg+pencilg3 ] - q[qpres][ijkg-pencilg3 ] ) +
                        DEL*( cons[ imy][ijkg+pencilg4 ]*qvp4 - cons[ imy][ijkg-pencilg4 ]*qvm4 + q[qpres][ijkg+pencilg4 ] - q[qpres][ijkg-pencilg4 ] ) );

      }
      for(i=lo0;i<=hi0;i++){

        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;

        double qvp1 = q[qv][ijkg+pencilg  ];
        double qvp2 = q[qv][ijkg+pencilg2 ];
        double qvp3 = q[qv][ijkg+pencilg3 ];
        double qvp4 = q[qv][ijkg+pencilg4 ];
        double qvm1 = q[qv][ijkg-pencilg  ];
        double qvm2 = q[qv][ijkg-pencilg2 ];
        double qvm3 = q[qv][ijkg-pencilg3 ];
        double qvm4 = q[qv][ijkg-pencilg4 ];

        flux[imz][ijk] = flux[imz][ijk]
             - dxinv1*( ALP*( cons[ imz][ijkg+pencilg  ]*qvp1 - cons[ imz][ijkg-pencilg  ]*qvm1 ) +
                        BET*( cons[ imz][ijkg+pencilg2 ]*qvp2 - cons[ imz][ijkg-pencilg2 ]*qvm2 ) +
                        GAM*( cons[ imz][ijkg+pencilg3 ]*qvp3 - cons[ imz][ijkg-pencilg3 ]*qvm3 ) +
                        DEL*( cons[ imz][ijkg+pencilg4 ]*qvp4 - cons[ imz][ijkg-pencilg4 ]*qvm4 ) );

      }
      for(i=lo0;i<=hi0;i++){

        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;

        double qvp1 = q[qv][ijkg+pencilg  ];
        double qvp2 = q[qv][ijkg+pencilg2 ];
        double qvp3 = q[qv][ijkg+pencilg3 ];
        double qvp4 = q[qv][ijkg+pencilg4 ];
        double qvm1 = q[qv][ijkg-pencilg  ];
        double qvm2 = q[qv][ijkg-pencilg2 ];
        double qvm3 = q[qv][ijkg-pencilg3 ];
        double qvm4 = q[qv][ijkg-pencilg4 ];

        flux[iene][ijk] = flux[iene][ijk]
             - dxinv1*( ALP*( cons[iene][ijkg+pencilg  ]*qvp1 - cons[iene][ijkg-pencilg  ]*qvm1 + q[qpres][ijkg+pencilg  ]*qvp1 - q[qpres][ijkg-pencilg  ]*qvm1 ) +
                        BET*( cons[iene][ijkg+pencilg2 ]*qvp2 - cons[iene][ijkg-pencilg2 ]*qvm2 + q[qpres][ijkg+pencilg2 ]*qvp2 - q[qpres][ijkg-pencilg2 ]*qvm2 ) +
                        GAM*( cons[iene][ijkg+pencilg3 ]*qvp3 - cons[iene][ijkg-pencilg3 ]*qvm3 + q[qpres][ijkg+pencilg3 ]*qvp3 - q[qpres][ijkg-pencilg3 ]*qvm3 ) +
                        DEL*( cons[iene][ijkg+pencilg4 ]*qvp4 - cons[iene][ijkg-pencilg4 ]*qvm4 + q[qpres][ijkg+pencilg4 ]*qvp4 - q[qpres][ijkg-pencilg4 ]*qvm4 ) );
        }
  }}
  uint64_t _time_L2 = CycleTime();
  _total_time_hypterm_L2 += (_time_L2 - _time_L1);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  int planeg2 = planeg*2;
  int planeg3 = planeg*3;
  int planeg4 = planeg*4;
  JBlocks = (hi1-lo1 + JBlockSize-1 )/JBlockSize;
  #pragma omp parallel for private(i,j,k,jb) schedule(static,1)
  for(jb=0;jb<JBlocks;jb++){
    for(k=lo2;k<=hi2;k++){
      for(j=jb*JBlockSize;j<(jb+1)*JBlockSize;j++)if(j<=hi1){
  //#pragma omp parallel for private(i,j,k) schedule(static,1)
  //for(k=lo2;k<hi2;k++){
  //  for(j=lo1;j<hi1;j++){

      for(i=lo0;i<=hi0;i++){
  
        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;

        double qwp1 = q[qw][ijkg+planeg ];
        double qwp2 = q[qw][ijkg+planeg2];
        double qwp3 = q[qw][ijkg+planeg3];
        double qwp4 = q[qw][ijkg+planeg4];
        double qwm1 = q[qw][ijkg-planeg ];
        double qwm2 = q[qw][ijkg-planeg2];
        double qwm3 = q[qw][ijkg-planeg3];
        double qwm4 = q[qw][ijkg-planeg4];

        flux[irho][ijk] = flux[irho][ijk]
             - dxinv2*( ALP*( cons[ imz][ijkg+planeg ] - cons[ imz][ijkg-planeg ]) +
                        BET*( cons[ imz][ijkg+planeg2] - cons[ imz][ijkg-planeg2]) +
                        GAM*( cons[ imz][ijkg+planeg3] - cons[ imz][ijkg-planeg3]) +
                        DEL*( cons[ imz][ijkg+planeg4] - cons[ imz][ijkg-planeg4]) );

      }
      for(i=lo0;i<=hi0;i++){

        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;

        double qwp1 = q[qw][ijkg+planeg ];
        double qwp2 = q[qw][ijkg+planeg2];
        double qwp3 = q[qw][ijkg+planeg3];
        double qwp4 = q[qw][ijkg+planeg4];
        double qwm1 = q[qw][ijkg-planeg ];
        double qwm2 = q[qw][ijkg-planeg2];
        double qwm3 = q[qw][ijkg-planeg3];
        double qwm4 = q[qw][ijkg-planeg4];

        flux[imx][ijk] = flux[imx][ijk]
             - dxinv2*( ALP*( cons[ imx][ijkg+planeg ]*qwp1 - cons[ imx][ijkg-planeg ]*qwm1 ) +
                        BET*( cons[ imx][ijkg+planeg2]*qwp2 - cons[ imx][ijkg-planeg2]*qwm2 ) +
                        GAM*( cons[ imx][ijkg+planeg3]*qwp3 - cons[ imx][ijkg-planeg3]*qwm3 ) +
                        DEL*( cons[ imx][ijkg+planeg4]*qwp4 - cons[ imx][ijkg-planeg4]*qwm4 ) );

      }
      for(i=lo0;i<=hi0;i++){
        
        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;
                        
        double qwp1 = q[qw][ijkg+planeg ];
        double qwp2 = q[qw][ijkg+planeg2];
        double qwp3 = q[qw][ijkg+planeg3];
        double qwp4 = q[qw][ijkg+planeg4];
        double qwm1 = q[qw][ijkg-planeg ];
        double qwm2 = q[qw][ijkg-planeg2];
        double qwm3 = q[qw][ijkg-planeg3];
        double qwm4 = q[qw][ijkg-planeg4];

        flux[imy][ijk] = flux[imy][ijk]
             - dxinv2*( ALP*( cons[ imy][ijkg+planeg ]*qwp1 - cons[ imy][ijkg-planeg ]*qwm1 ) +
                        BET*( cons[ imy][ijkg+planeg2]*qwp2 - cons[ imy][ijkg-planeg2]*qwm2 ) +
                        GAM*( cons[ imy][ijkg+planeg3]*qwp3 - cons[ imy][ijkg-planeg3]*qwm3 ) +
                        DEL*( cons[ imy][ijkg+planeg4]*qwp4 - cons[ imy][ijkg-planeg4]*qwm4 ) );

      }
      for(i=lo0;i<=hi0;i++){
        
        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;
                        
        double qwp1 = q[qw][ijkg+planeg ];
        double qwp2 = q[qw][ijkg+planeg2];
        double qwp3 = q[qw][ijkg+planeg3];
        double qwp4 = q[qw][ijkg+planeg4];
        double qwm1 = q[qw][ijkg-planeg ];
        double qwm2 = q[qw][ijkg-planeg2];
        double qwm3 = q[qw][ijkg-planeg3];
        double qwm4 = q[qw][ijkg-planeg4];

        flux[imz][ijk] = flux[imz][ijk]
             - dxinv2*( ALP*( cons[ imz][ijkg+planeg ]*qwp1 - cons[ imz][ijkg-planeg ]*qwm1 + q[qpres][ijkg+planeg ] - q[qpres][ijkg-planeg ] ) +
                        BET*( cons[ imz][ijkg+planeg2]*qwp2 - cons[ imz][ijkg-planeg2]*qwm2 + q[qpres][ijkg+planeg2] - q[qpres][ijkg-planeg2] ) +
                        GAM*( cons[ imz][ijkg+planeg3]*qwp3 - cons[ imz][ijkg-planeg3]*qwm3 + q[qpres][ijkg+planeg3] - q[qpres][ijkg-planeg3] ) +
                        DEL*( cons[ imz][ijkg+planeg4]*qwp4 - cons[ imz][ijkg-planeg4]*qwm4 + q[qpres][ijkg+planeg4] - q[qpres][ijkg-planeg4] ) );

      }
      for(i=lo0;i<=hi0;i++){
        
        int ijk  = i + j*pencil  + k*plane;
        int ijkg = i + j*pencilg + k*planeg;
                        
        double qwp1 = q[qw][ijkg+planeg ];
        double qwp2 = q[qw][ijkg+planeg2];
        double qwp3 = q[qw][ijkg+planeg3];
        double qwp4 = q[qw][ijkg+planeg4];
        double qwm1 = q[qw][ijkg-planeg ];
        double qwm2 = q[qw][ijkg-planeg2];
        double qwm3 = q[qw][ijkg-planeg3];
        double qwm4 = q[qw][ijkg-planeg4];

        flux[iene][ijk] = flux[iene][ijk]
             - dxinv2*( ALP*( cons[iene][ijkg+planeg ]*qwp1 - cons[iene][ijkg-planeg ]*qwm1 + q[qpres][ijkg+planeg ]*qwp1 - q[qpres][ijkg-planeg ]*qwm1 ) +
                        BET*( cons[iene][ijkg+planeg2]*qwp2 - cons[iene][ijkg-planeg2]*qwm2 + q[qpres][ijkg+planeg2]*qwp2 - q[qpres][ijkg-planeg2]*qwm2 ) +
                        GAM*( cons[iene][ijkg+planeg3]*qwp3 - cons[iene][ijkg-planeg3]*qwm3 + q[qpres][ijkg+planeg3]*qwp3 - q[qpres][ijkg-planeg3]*qwm3 ) +
                        DEL*( cons[iene][ijkg+planeg4]*qwp4 - cons[iene][ijkg-planeg4]*qwm4 + q[qpres][ijkg+planeg4]*qwp4 - q[qpres][ijkg-planeg4]*qwm4 ) );
      }
  }}}
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  uint64_t _time_L3 = CycleTime();
  _total_time_hypterm_L3 += (_time_L3 - _time_L2);
  _total_time_hypterm    += (_time_L3 - _time_start);


}
//------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv){
  int MPI_Rank=0;
  int MPI_Tasks=1;
  int OMP_Threads = 1;

  #pragma omp parallel 
  {
    #pragma omp master
    {
//////////////////////      OMP_Threads = omp_get_num_threads();
    }
  }


  int max_grid_size = 64;
  int NG = 4;
  int NC = 5;

  int lo[3] = {0,0,0};
  int hi[3] = {max_grid_size-1,max_grid_size-1,max_grid_size-1};
  double dx[3] = {1.0,1.0,1.0};
  double * tmpbuf;
  int volume,c;

  double **U; // NC, NG
  posix_memalign((void**)&U,64,NC*sizeof(double*));
  volume = (max_grid_size+2*NG)*(max_grid_size+2*NG)*(max_grid_size+2*NG);
  posix_memalign((void**)&tmpbuf,64,volume*(NC)*sizeof(double));memset(tmpbuf,0,volume*(NC)*sizeof(double));
  for(c=0;c<NC;c++){
    U[c] = tmpbuf + c*volume;
  }
  double **D = (double**)malloc(NC*sizeof(double*)); // NC, 0
  posix_memalign((void**)&D,64,NC*sizeof(double*));
  volume = max_grid_size*max_grid_size*max_grid_size;
  posix_memalign((void**)&tmpbuf,64,volume*(NC)*sizeof(double));memset(tmpbuf,0,volume*(NC)*sizeof(double));
  for(c=0;c<NC;c++){
    D[c] = tmpbuf + c*volume;
  }
  double **F = (double**)malloc(NC*sizeof(double*)); // NC, 0
  posix_memalign((void**)&F,64,NC*sizeof(double*));
  volume = max_grid_size*max_grid_size*max_grid_size;
  posix_memalign((void**)&tmpbuf,64,volume*(NC)*sizeof(double));memset(tmpbuf,0,volume*(NC)*sizeof(double));
  for(c=0;c<NC;c++){
    F[c] = tmpbuf + c*volume;
  }
  double **Q = (double**)malloc((NC+1)*sizeof(double*)); // NC+1, NG
  volume = (max_grid_size+2*NG)*(max_grid_size+2*NG)*(max_grid_size+2*NG);
  posix_memalign((void**)&tmpbuf,64,volume*(NC+1)*sizeof(double));memset(tmpbuf,0,volume*(NC+1)*sizeof(double));
  for(c=0;c<NC+1;c++){
    Q[c] = tmpbuf + c*volume;
  }
  double **Unew = (double**)malloc(NC*sizeof(double*)); // NC, NG
  volume = (max_grid_size+2*NG)*(max_grid_size+2*NG)*(max_grid_size+2*NG);
  posix_memalign((void**)&tmpbuf,64,volume*(NC)*sizeof(double));memset(tmpbuf,0,volume*(NC)*sizeof(double));
  for(c=0;c<NC;c++){
    Unew[c] = tmpbuf + c*volume;
  }


  init_data(lo,hi,NG,dx,U);

  uint64_t t0 = CycleTime();
  sleep(1);
  uint64_t t1 = CycleTime();
  double frequency = (double)(t1-t0);

  _total_run_time        = 0;
  _total_time_hypterm    = 0;
  _total_time_hypterm_L1 = 0;
  _total_time_hypterm_L2 = 0;
  _total_time_hypterm_L3 = 0;
  int iteration,NIterations = 10;
  uint64_t _run_time_start      = CycleTime();
  for(iteration=0;iteration<NIterations;iteration++){
   // for all boxes...
   hypterm_naive(lo,hi,NG,dx,U,Q,F);
  }
  uint64_t _run_time_end = CycleTime();
  _total_run_time += (_run_time_end - _run_time_start);

  #ifdef JBlockSize
  printf("JBlockSize = %d\n",JBlockSize);
  printf("volume     = %d\n",volume    );
  #endif
  printf("frequency = %6.6e\n",frequency);
  printf("runtime = %9.6fs\n",(double)_total_run_time   /frequency);
  printf("hypterm = %9.6fs\n",(double)_total_time_hypterm   /(double)NIterations/frequency);
  printf("     L1 = %9.6fs\n",(double)_total_time_hypterm_L1/(double)NIterations/frequency);
  printf("     L2 = %9.6fs\n",(double)_total_time_hypterm_L2/(double)NIterations/frequency);
  printf("     L3 = %9.6fs\n",(double)_total_time_hypterm_L3/(double)NIterations/frequency);

}
