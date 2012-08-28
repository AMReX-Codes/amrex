// --------------------------------------------------------------------
// HyptermOnly.cpp
// --------------------------------------------------------------------
#include <winstd.H>

#include <new>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using std::ios;
using std::cout;
using std::endl;

#ifndef WIN32
#include <unistd.h>
#endif

#include <IntVect.H>
#include <Box.H>
#include <FArrayBox.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <Utility.H>
#include <VisMF.H>
#include <HyptermOnly_F.H>

#ifdef SHOWVAL
#undef SHOWVAL
#endif
#define SHOWVAL(val) { cout << #val << " = " << val << endl; }

#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

extern "C" {
  void hypterm_naive(int lo[3], int hi[3], int ng, double dx[3], double ** __restrict__ _cons, double ** __restrict__ _q, double ** __restrict__ _flux);
};

// --------------------------------------------------------------------
int main(int argc, char *argv[]) {
    BoxLib::Initialize(argc,argv);    

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
      FArrayBox &myFabU = mfU[mfi];
      FArrayBox &myFabQ = mfQ[mfi];
      const Real *dataPtrU = myFabU.dataPtr();
      const Real *dataPtrQ = myFabQ.dataPtr();
      const int  *dlo     = myFabU.loVect();
      const int  *dhi     = myFabU.hiVect();
      const Real *dxptr   = dx;
      FORT_INITDATA(dataPtrU, ARLIM(dlo), ARLIM(dhi), dataPtrQ, dxptr, &nComps);
    }

    VisMF::Write(mfU, "mfUInit");


    {
      double tstart = BoxLib::wsecond();
      sleep(1);
      double tend = BoxLib::wsecond();
      if(ParallelDescriptor::IOProcessor()) {
        cout << "sleep(1) time = " << tend - tstart << endl;
      }
    }

    double tstart = BoxLib::wsecond();

    int nSteps(1);
    for(int iStep(0); iStep < nSteps; ++iStep) {
    for(MFIter mfi(mfU); mfi.isValid(); ++mfi) {

      FArrayBox &myFabU = mfU[mfi];
      FArrayBox &myFabQ = mfQ[mfi];
      FArrayBox &myFabFlux = mfFlux[mfi];

      int idx = mfi.index();
      cout << "myproc idx = " << ParallelDescriptor::MyProc() << "  " << idx << endl;

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

      hypterm_naive(lo,hi,NG,dx,U,Q,F);
    }
    }

    double tend = BoxLib::wsecond();
    if(ParallelDescriptor::IOProcessor()) {
      tstart += (static_cast<double> (nSteps));  // adjust sleep time from bench.c
      cout << "Hypterm time    =  " << tend - tstart << endl;
      cout << "Hypterm time/it =  " << (tend - tstart) / (static_cast<double> (nSteps)) << endl;
    }

    VisMF::Write(mfFlux, "mfFluxFinal");

    BoxLib::Finalize();
    return 0;
}
// --------------------------------------------------------------------
