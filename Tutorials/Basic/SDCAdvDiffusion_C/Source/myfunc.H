#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
// #include <AMReX_LO_BCTYPES.H>

using namespace amrex;

const int SDC_NNODES=5;
const int SDC_NPIECES=2;
 class SDCstuff
 {
   //   void pf_quadrature(int* qtype_in, int* nnodes, int* nnodes0,
   //		      amrex_real* nodes, int* nflags,amrex_real* qmats);


 public:
   Real qnodes [SDC_NNODES];
   int  nflags [SDC_NNODES];
   Real qmats [4][SDC_NNODES-1][SDC_NNODES];
   
   int qtype_in=1;
   int nnodes=SDC_NNODES;
   int nsweeps =8;

   //   SDCstuff()
   //   {
   //     int qtype_in=1;
   //     pf_quadrature(&qtype_in, &nnodes, &nnodes,qnodes,nflags, &qmats[0][0][0]);
   //   }

 };

void main_main ();

void sweep (MultiFab& phi_old,
	    MultiFab& phi_new,
	    std::array<MultiFab, AMREX_SPACEDIM>& flux,
	    Vector<MultiFab>& phi_sdc,
	    Vector<MultiFab>& res_sdc,
	    Vector<Vector<MultiFab> >& f_sdc,
	    Real dt,Real v,Real nu,
	    const Geometry& geom, 
	    const BoxArray& grids,
	    const DistributionMapping& dmap, 
	    const Vector<BCRec>& bc,
            SDCstuff sdcmats);

void compute_integral_term(Vector<MultiFab>& phi_sdc,
			     Vector<MultiFab>& res_sdc,
			     Vector<Vector<MultiFab> >& f_sdc,
			     Real dt, SDCstuff sdcmats);

void evaluate_f_m(Vector<MultiFab>& phi_sdc,
		    Vector<Vector<MultiFab> >& f_sdc,
		    std::array<MultiFab, AMREX_SPACEDIM>& flux,
		    const Geometry& geom, Real v,Real nu,
		    int m);

#endif
