// An example showing how to clone a MultiFab from the compute MPI group to the
// "sidecar" group. The tricky part about doing this is pointer management. The
// functions SendMultiFabToSidecars() and SendGeometryToSidecars() do no
// pointer allocation; they assume that the memory space has already been
// allocated beforehard. However, some (most) of these pointers will point to
// actual data only one MPI group, while on the other they will point to
// nothing. But they still need to exist on both MPI groups because everybody
// must call the Send...() functions. So be mindful of your pointers; if you do
// then this can be a very useful tool.

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

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {

    // Unfortunately the # of sidecars is currently a compile-time constant
    // because ParmParse must come AFTER Initialize(), but setting the # of
    // sidecars must come BEFORE.

    // TODO: change Initialize() so that we can read # of sidecars from an
    // inputs file.

    const int nSidecarProcs(6);
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

    MultiFab *mf = new MultiFab;
    MultiFab *MF_sidecar = new MultiFab;

    Box *baseBox;
    RealBox* real_box;
    int* coord;
    Geometry* geom;
    int is_per[BL_SPACEDIM];

    double t_start, t_stop;

    if (ParallelDescriptor::InCompGroup())
    {
      // Make a Box, then a BoxArray with maxSize.
      baseBox = new Box(IntVect(0,0,0), IntVect(maxGrid-1, maxGrid-1, maxGrid-1));
      BoxArray ba(*baseBox);
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
      mf->define(ba, nComp, nGhost, comp_DM, Fab_allocate);
      for(int i(0); i < mf->nComp(); ++i) {
        mf->setVal(rand()%100, i, 1);
      }

      // This defines the physical size of the box.
      real_box = new RealBox;
      for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box->setLo(n, 0.0);
        real_box->setHi(n, 42.0);
      }

      // This says we are using Cartesian coordinates
      coord = new int(0);

      // This sets the boundary conditions to be non-periodic.
      for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 0;

      // This defines a Geometry object which is useful for writing the plotfiles
      geom = new Geometry(*baseBox, real_box, *coord, is_per);
    }
    else // on the sidecars, just allocate the memory but don't do anything else
    {
      baseBox = new Box;
      real_box = new RealBox;
      coord = new int;
      geom = new Geometry;
    }

    // whoooosh
    t_start = MPI_Wtime();
    MultiFab::SendMultiFabToSidecars (mf, MF_sidecar);
    Geometry::SendGeometryToSidecars (baseBox, real_box, coord, is_per, geom);
    t_stop = MPI_Wtime();

    if (ParallelDescriptor::InSidecarGroup())
    {
      const double norm0 = MF_sidecar->norm0();
      const double norm1 = MF_sidecar->norm1();
      const double norm2 = MF_sidecar->norm2();
      const Real probsize = geom->ProbSize();
      if (ParallelDescriptor::IOProcessor()) {
        std::cout << "From sidecar nodes, norms (L0, L1, L2) of MF are (" << norm0 << ", " << norm1 << ", " << norm2 << ")" << std::endl;
        std::cout << "From sidecar nodes, probsize = " << probsize << std::endl;
      }
    }
    else
    {
      const double norm0 = mf->norm0();
      const double norm1 = mf->norm1();
      const double norm2 = mf->norm2();
      const Real probsize = geom->ProbSize();
      if (ParallelDescriptor::IOProcessor()) {
        std::cout << "From compute nodes, norms (L0, L1, L2) of MF are (" << norm0 << ", " << norm1 << ", " << norm2 << ")" << std::endl;
        std::cout << "From compute nodes, probsize = " << probsize << std::endl;

        std::cout << std::endl << std::endl << "total time for MultiFab and Geometry transfer: " << t_stop - t_start << " sec" << std::endl;
      }
    }

    delete mf;
    delete MF_sidecar;
    delete baseBox;
    delete real_box;
    delete coord;
    delete geom;

    BoxLib::Finalize();

    return 0;
}
