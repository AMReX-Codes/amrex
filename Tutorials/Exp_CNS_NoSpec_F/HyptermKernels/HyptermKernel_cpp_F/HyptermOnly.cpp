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

    mfU.setVal(0.0);
    mfQ.setVal(0.0);
    mfFlux.setVal(0.2);

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

      int idx = mfi.index();

      const Real *dataPtrU = myFabU.dataPtr();
      const Real *dataPtrQ = myFabQ.dataPtr();
      const int  *dlo      = myFabU.loVect();
      const int  *dhi      = myFabU.hiVect();
      const Real *dxptr    = dx;

      FORT_INITDATA(dataPtrU, ARLIM(dlo), ARLIM(dhi), dataPtrQ, dxptr, &nComps);
    }

    VisMF::Write(mfU, "mfUInit");


    double tsleepstart = BoxLib::wsecond();
    sleep(1);
    double tsleepend = BoxLib::wsecond();
    if(ParallelDescriptor::IOProcessor()) {
      cout << "sleep(1) time = " << tsleepend - tsleepstart << endl;
    }
    if(ParallelDescriptor::IOProcessor()) {
      for(int i(0); i < BL_SPACEDIM; ++i) {
        cout << "dx[" << i << "] = " << dx[i] << endl;
      }
    }

    double tstart = BoxLib::wsecond();

    int nSteps(1);
    for(int iStep(0); iStep < nSteps; ++iStep) {
      for(MFIter mfi(mfU); mfi.isValid(); ++mfi) {
        int idx = mfi.index();
        FArrayBox &myFabU = mfU[mfi];
        FArrayBox &myFabQ = mfQ[mfi];
        FArrayBox &myFabFlux = mfFlux[mfi];

        const Real *dataPtrU = myFabU.dataPtr();
        const Real *dataPtrQ = myFabQ.dataPtr();
        const Real *dataPtrFlux = myFabFlux.dataPtr();
        const int  *dlo     = myFabU.loVect();
        const int  *dhi     = myFabU.hiVect();
        const int  *lo      = mfU.boxArray()[idx].loVect();
        const int  *hi      = mfU.boxArray()[idx].hiVect();
        const Real *dxptr   = dx;

        FORT_HYPTERM(dataPtrU, ARLIM(dlo), ARLIM(dhi), ARLIM(lo), ARLIM(hi),
                     dataPtrQ, dataPtrFlux, dxptr, &nComps);
        //FORT_HYPTERM_UNOPT(dataPtrU, ARLIM(dlo), ARLIM(dhi), ARLIM(lo), ARLIM(hi),
                     //dataPtrQ, dataPtrFlux, dxptr, &nComps);

        if(ParallelDescriptor::IOProcessor()) {
	  for(int icomp(0); icomp < nComps; ++icomp) {
            cout << "minmax flux[" << icomp << "] = "
	         << myFabFlux.min(icomp) << "  " << myFabFlux.max(icomp) << endl;
	  }
        }
      }
    }

    double tend = BoxLib::wsecond();
    if(ParallelDescriptor::IOProcessor()) {
      cout << "-----------------" << endl;
      cout << "Hypterm time    =  " << tend - tstart << endl;
      cout << "Hypterm time/it =  " << (tend - tstart) / (static_cast<double> (nSteps)) << endl;
      cout << "-----------------" << endl;
    }

    VisMF::Write(mfFlux, "mfFluxFinal");

    BoxLib::Finalize();
    return 0;
}
// --------------------------------------------------------------------
