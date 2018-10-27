#ifndef SDCSTRUCT_H_
#define SDCSTRUCT_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>

using namespace amrex;

//  External routine for computing quadrature tables
extern "C"  {
void SDC_quadrature(int* qtype_in, int* nnodes, int* nnodes0,
		   amrex_real* nodes, int* nflags,amrex_real* qmats);
}
const int SDC_NNODES=5;
class SDCstruct
 {

 public:
   //   Real qnodes [];
   //   int  Nflags [];
   Real qmats [4][SDC_NNODES-1][SDC_NNODES];
   
   //   Real qnodes [SDC_NNODES];
   Real* qnodes;
   int*  Nflags;   //  Flags for quadrature rules
   //   Real* qmats;    //  Quadrature tables [4][NNODES-1][NNODES]; 
   
   int Nnodes;      //  Number of quadrature nodes
   int qtype=1;     //  Type of quadrature nodes
   int Nsweeps=8;   //  Number of sweeps per time step
   int Npieces;     //  Number of terms in RHS
   
   Vector<MultiFab> sol;  //  Solution at the nodes
   Vector<Vector<MultiFab> > f; // Functions a nodes access by [npieces][node]   
   Vector<MultiFab> res;        // Temp storage 

   //  Constructor
   SDCstruct(int Nnodes_in,int Npieces_in, MultiFab& sol_in);

   //  Sweeper routines
   void SDC_rhs_integrals(Real dt);
   void SDC_rhs_k_plus_one(MultiFab& rhs, Real dt,int m);

 };

#endif
