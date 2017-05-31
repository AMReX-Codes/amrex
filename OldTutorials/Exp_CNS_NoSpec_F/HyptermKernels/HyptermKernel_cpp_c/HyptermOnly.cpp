// --------------------------------------------------------------------
// HyptermOnly.cpp
// --------------------------------------------------------------------

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using std::ios;
using std::cout;
using std::endl;

#include <unistd.h>

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#ifdef SHOWVAL
#undef SHOWVAL
#endif
#define SHOWVAL(val) { cout << #val << " = " << val << endl; }

extern "C" {
  void init_data(int lo[3], int hi[3], const int fablo[3], const int fabhi[3], int ng, double dx[3],
                 double ** __restrict__ _cons, double ** __restrict__ _q);

  void hypterm_naive(int lo[3], int hi[3], int ng, double dx[3], double ** __restrict__ _cons,
                     double ** __restrict__ _q, double ** __restrict__ _flux);
  void init_timer();
};


using namespace amrex;

// --------------------------------------------------------------------
int main(int argc, char *argv[]) {
    amrex::Initialize(argc,argv);    

    int maxgrid(64), nComps(5), nGhost(4);
    int ioProc(ParallelDescriptor::IOProcessorNumber());

    FArrayBox::setFormat(FABio::FAB_NATIVE);

    IntVect ivlo(0, 0, 0);
    IntVect ivhi(maxgrid-1, maxgrid-1, maxgrid-1);
    Box probDomain(ivlo, ivhi);
    BoxArray pdBoxArray(probDomain);
    MultiFab mfFlux(pdBoxArray, nComps, 0, Fab_allocate);
    MultiFab mfU(pdBoxArray, nComps, nGhost, Fab_allocate);
    MultiFab mfQ(pdBoxArray, nComps+1, nGhost, Fab_allocate);

    Array<Real> probLo(BL_SPACEDIM), probHi(BL_SPACEDIM);
    for(int i(0); i < BL_SPACEDIM; ++i) {
      probLo[i] = -2.3;
      probHi[i] =  2.3;
    }
    
    Real dx[BL_SPACEDIM];

    for(int i(0); i < BL_SPACEDIM; ++i) {
      dx[i] = (probHi[i] - probLo[i]) / (static_cast<Real> (probDomain.length(i)));
    }

    for(MFIter mfi(mfU); mfi.isValid(); ++mfi) {
        int idx = mfi.index();
        FArrayBox &myFabU = mfU[mfi];
        FArrayBox &myFabQ = mfQ[mfi];

        const int  *dlo     = myFabU.loVect();
        const int  *dhi     = myFabU.hiVect();
        const Real *dxptr   = dx;

        int lo[BL_SPACEDIM], hi[BL_SPACEDIM];
        for(int i(0); i < BL_SPACEDIM; ++i) {
          lo[i] = mfU.boxArray()[idx].loVect()[i];
          hi[i] = mfU.boxArray()[idx].hiVect()[i];
        }
        int NG(nGhost);
        double *U[nComps];
        double *Q[nComps+1];
        for(int i(0); i < nComps; ++i) {
          U[i] = myFabU.dataPtr(i);
          Q[i] = myFabQ.dataPtr(i);
        }
        Q[nComps] = myFabQ.dataPtr(nComps);

        init_data(lo,hi,dlo,dhi,NG,dx,U,Q);
    }

    VisMF::Write(mfU, "mfUInit");


    double tsleepstart = ParallelDescriptor::second();
    sleep(1);
    double tsleepend = ParallelDescriptor::second();
    if(ParallelDescriptor::IOProcessor()) {
      cout << "sleep(1) time = " << tsleepend - tsleepstart << endl;
    }

    for(int i(0); i < BL_SPACEDIM; ++i) {
      cout << "dx[" << i << "] = " << dx[i] << endl;
    }

    init_timer();

    double tstart = ParallelDescriptor::second();

    int nSteps(1);
    for(int iStep(0); iStep < nSteps; ++iStep) {
      for(MFIter mfi(mfU); mfi.isValid(); ++mfi) {
        int idx = mfi.index();
        FArrayBox &myFabU = mfU[mfi];
        FArrayBox &myFabQ = mfQ[mfi];
        FArrayBox &myFabFlux = mfFlux[mfi];

        const int  *dlo     = myFabU.loVect();
        const int  *dhi     = myFabU.hiVect();
        const Real *dxptr   = dx;

        int lo[BL_SPACEDIM], hi[BL_SPACEDIM];
        for(int i(0); i < BL_SPACEDIM; ++i) {
          lo[i] = mfU.boxArray()[idx].loVect()[i];
          hi[i] = mfU.boxArray()[idx].hiVect()[i];
        }
        int NG(nGhost);
        double *U[nComps];
        double *Q[nComps+1];
        double *F[nComps];
        for(int i(0); i < nComps; ++i) {
          U[i] = myFabU.dataPtr(i);
          Q[i] = myFabQ.dataPtr(i);
          F[i] = myFabFlux.dataPtr(i);
        }
        Q[nComps] = myFabQ.dataPtr(nComps);

        if(ParallelDescriptor::IOProcessor()) {
	  cout << "-----------------" << endl;
	}

        double tstart = ParallelDescriptor::second();

        hypterm_naive(lo,hi,NG,dx,U,Q,F);

        double tend = ParallelDescriptor::second();
        if(ParallelDescriptor::IOProcessor()) {
	  cout << "-----------------" << endl;
          cout << "hypterm =  " << tend - tstart << endl;
	}

      }
    }

    double tend = ParallelDescriptor::second();
    if(ParallelDescriptor::IOProcessor()) {
      cout << "-----------------" << endl;
      cout << "runtime tot =  " << tend - tstart << endl;
      cout << "runtime /it =  " << (tend - tstart) / (static_cast<double> (nSteps)) << endl;
      cout << "-----------------" << endl;
    }

    VisMF::Write(mfFlux, "mfFluxFinal");

    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------
