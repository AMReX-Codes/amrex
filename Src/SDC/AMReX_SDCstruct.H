#ifndef SDCSTRUCT_H_
#define SDCSTRUCT_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
// #include <AMReX_LO_BCTYPES.H>
   

using namespace amrex;

extern "C"  {
void pf_quadrature(int* qtype_in, int* nnodes, int* nnodes0,
		   amrex_real* nodes, int* nflags,amrex_real* qmats);
}
const int SDC_NNODES=5;
class SDCstruct
 {

 public:
   //   Real qnodes [];
   //   int  Nflags [];
   //   Real qmats [4][SDC_NNODES-1][SDC_NNODES];
   
   Real qnodes [SDC_NNODES];
   int  Nflags [SDC_NNODES];
   Real qmats [4][SDC_NNODES-1][SDC_NNODES];
   
   int qtype_in=1;
   int Nnodes;
   int Nsweeps=8;
   int Npieces;
   
   Vector<MultiFab> phi_sdc;
   Vector<MultiFab> res_sdc;
   Vector<MultiFab> f_sdc; // access by [npieces][node]
   //   Vector<Vector<MultiFab> > f_sdc; // access by [npieces][node]   

   SDCstruct(int Nnodes_in,int Npieces_in, MultiFab& phi_in);


 };

#endif
