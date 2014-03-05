// --------------------------------------------------------------------------
// MultiFabFillBoundaryTest.cpp
// --------------------------------------------------------------------------
//   this file tests fillboundary.
// --------------------------------------------------------------------------
#include <winstd.H>

#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <cmath>

#ifndef WIN32
#include <unistd.h>
#endif

#include <IntVect.H>
#include <Box.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <Utility.H>

#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

const int maxGrid(64);

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    BoxLib::Initialize(argc,argv);    

    BL_PROFILE_VAR("main()", pmain);

    // ---- First use the number of processors to decide how many grids you have.
    // ---- We arbitrarily decide to have one grid per MPI process in a uniform
    // ---- cubic domain, so we require that the number of processors be N^3. 
    // ---- This requirement is somewhat arbitrary, but convenient for now.

    int nprocs = ParallelDescriptor::NProcs();

    // This is the cube root of the number of processors
    int N = exp(log(abs(nprocs))/3.0);

    if (N*N*N != nprocs) 
    {
        BoxLib::Error("We require that the number of processors be a perfect cube");
    }

    // std::cout << "Cube root of " << nprocs << " is " << N << std::endl;


    // Don't restrict ourselves to on-processor communication
    bool local(false);

    // ---- If cross == true then only the faces are exchanged
    // ---- If cross == false then the faces, edges and corners are all exchanged

    // ---- make a box, then a boxarray with maxSize
    int domain_hi = (N*maxGrid) - 1;
    Box Domain(IntVect(0,0,0), IntVect(domain_hi,domain_hi,domain_hi));
    BoxArray ba(Domain);
    ba.maxSize(maxGrid);

    // ---- Below we will make a MultiFab and set the values to 1.0
    // ----  (this is arbitrary, just makes the values not be undefined)
    // ---- We will do this for nGhost = 1,...,4 and nComp = 1, 4, 20
    // ----  (these are also arbitrary, just meant to be representative)
    // ---- and for "cross" = true or false.


    const std::vector<bool> cross = { true, false };
    const std::vector<int>  nComp = { 1, 4, 20 };
    const std::vector<int> nGhost = { 1, 2, 3, 4 };
    // do this for older compilers
    //static const int nGhostarr[] = {1, 2, 3, 4};
    //std::vector<int> nGhost(nGhostarr, nGhostarr + sizeof(nGhostarr) / sizeof(nGhostarr[0]) );

    for(int icross(0); icross < cross.size(); ++icross) {
      for(int icomp(0); icomp < nComp.size(); ++icomp) {
        for(int ighost(0); ighost < nGhost.size(); ++ighost) {

          std::ostringstream nametag;
          nametag << "FB_nGhost" << nGhost[ighost] << "_nComp" << nComp[icomp]
	          << "_cross" << (cross[icross] ? "True":"False");
	  if(ParallelDescriptor::IOProcessor()) {
	    std::cout << "Working on:"
	              << "  Ghost = " << nGhost[ighost]
	              << "  nComp = " << nComp[icomp]
		      << "  cross = " << (cross[icross] ? "true":"false")
		      << std::endl;
	  }

          MultiFab mf(ba, nComp[icomp], nGhost[ighost]);
          mf.setVal(1.0);

          BL_COMM_PROFILE_NAMETAG(nametag.str() + "_Start");

          mf.FillBoundary(local, cross[icross]);

          BL_COMM_PROFILE_NAMETAG(nametag.str() + "_End");

        }
      }
    }

    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "Finished." << std::endl;
    }

    BL_PROFILE_VAR_STOP(pmain);

    BoxLib::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
