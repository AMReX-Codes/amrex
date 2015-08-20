// An example demonstrating the usage of the abstract base class "Analysis" for
// doing data analysis during a simulation. It also makes use of the sidecar
// (in-transit) technology.

#include <new>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <Geometry.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>
#include <RealBox.H>
#include <Utility.H>

#include <InSituAnalysis.H>

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    // Unfortunately the # of sidecars is currently a compile-time constant
    // because ParmParse must come AFTER Initialize(), but setting the # of
    // sidecars must come BEFORE.

    // TODO: change Initialize() so that we can read # of sidecars from an
    // inputs file.

    const int nSidecarProcs(2);
    ParallelDescriptor::SetNProcsSidecar(nSidecarProcs);

    BoxLib::Initialize(argc,argv);

    int maxGrid;
    int nComp;
    int nGhost;
    int maxSize;

    ParmParse pp;
    pp.get("maxGrid", maxGrid);
    pp.get("nComp", nComp);
    pp.get("nGhost", nGhost);
    pp.get("maxSize", maxSize);

    MultiFab *mf;
    Geometry* geom;

    if (ParallelDescriptor::InCompGroup())
    {
      // Make a Box, then a BoxArray with maxSize.
      Box baseBox(IntVect(0,0,0), IntVect(maxGrid-1, maxGrid-1, maxGrid-1));
      BoxArray ba(baseBox);
      ba.maxSize(maxSize);

      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "ba.size() = " << ba.size() << std::endl;
        std::cout << "ba[0] = " << ba[0] << std::endl;
        std::cout << "ba[1] = " << ba[1] << std::endl;
        std::cout << "ba[last] = " << ba[ba.size() - 1] << std::endl;
      }

      // This is the DM for the compute processes.
      DistributionMapping comp_DM;
      comp_DM.define(ba, ParallelDescriptor::NProcsComp());

      // Make a MultiFab and populate it with a bunch of random numbers.
      mf = new MultiFab;
      mf->define(ba, nComp, nGhost, comp_DM, Fab_allocate);
      for(int i(0); i < mf->nComp(); ++i) {
        mf->setVal(rand()%100, i, 1);
      }

      // This defines the physical size of the box.
      RealBox real_box;
      for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 42.0);
      }

      // This says we are using Cartesian coordinates
      int coord(0);

      // This sets the boundary conditions to be non-periodic.
      int is_per[BL_SPACEDIM];
      for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 0;

      // This defines a Geometry object which is useful for writing the plotfiles
      geom = new Geometry(baseBox, &real_box, coord, is_per);
    }
    else // on the sidecars, just allocate the memory but don't do anything else
    {
      geom = new Geometry;
      mf = new MultiFab;
    }

    // whoooosh
    MultiFab::SendMultiFabToSidecars (mf);
    Geometry::SendGeometryToSidecars (geom);

    if (ParallelDescriptor::InSidecarGroup())
    {
      InSituAnalysis analysis(*mf, *geom);
      analysis.DoAnalysis();
      if (ParallelDescriptor::IOProcessor()) {
        std::cout << "From sidecar nodes, here is the analysis:" << std::endl;
        analysis.PrintResults();
      }
    }

    delete mf;
    delete geom;

    BoxLib::Finalize();

    return 0;
}
