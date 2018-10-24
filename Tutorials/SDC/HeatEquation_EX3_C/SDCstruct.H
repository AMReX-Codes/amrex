#ifndef SDCSTRUCT_H_
#define SDCSTRUCT_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
// #include <AMReX_LO_BCTYPES.H>

using namespace amrex;

const int SDC_NNODES=5;
const int SDC_NPIECES=2;

extern "C"  {
void pf_quadrature(int* qtype_in, int* nnodes, int* nnodes0,
		   amrex_real* nodes, int* nflags,amrex_real* qmats);
}

class SDCstruct
 {
   

 public:
   Real qnodes [SDC_NNODES];
   int  nflags [SDC_NNODES];
   Real qmats [4][SDC_NNODES-1][SDC_NNODES];
   
   int qtype_in=1;
   int nnodes=SDC_NNODES;
   int nsweeps =8;

   SDCstruct()
      {
       int qtype_in=1;
       pf_quadrature(&qtype_in, &nnodes, &nnodes,qnodes,nflags, &qmats[0][0][0]);
     }

 };

#endif
