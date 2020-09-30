#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>

#include "myfunc_F.H"

// This tutorial demonstrates how to use CVODE to integrate a system of 3
// coupled ODEs per cell. It repeats the ODE integration twice, once with the
// analytic Jacobian matrix defined by the user, and again with no Jacobian
// matrix; in the latter case, CVODE constructs the Jacobian (assuming Newton
// iteration is used) by numerical differentiation. The latter case should be
// slower than the former, due to the additional computational expense of
// building and rebuilding the Jacobian during the iteration.
//
// The system of equations is taken from the fcvRoberts_dns.f example code in
// CVODE.

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    Real strt_time, stop_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    std::cout << std::setprecision(15);

    int n_cell, max_grid_size;
    int cvode_meth, cvode_itmeth, write_plotfile;
    bool do_tiling;

    // inputs parameters
    {
      // ParmParse is way of reading inputs from the inputs file
      ParmParse pp;

      // We need to get n_cell from the inputs file - this is the number of
      // cells on each side of a square (or cubic) domain.
      pp.get("n_cell",n_cell);

      // Default nsteps to 0, allow us to set it to something else in the
      // inputs file
      pp.get("max_grid_size",max_grid_size);

      // Select CVODE solve method.
      //   1 = Adams (for non-stiff problems)
      //   2 = BDF (for stiff problems)
      pp.get("cvode_meth",cvode_meth);
      // Select CVODE solver iteration method.
      //   1 = Functional iteration
      //   2 = Newton iteration
      pp.get("cvode_itmeth",cvode_itmeth);

      pp.get("write_plotfile",write_plotfile);
      pp.get("do_tiling",do_tiling);
    }

    if (cvode_meth < 1)
      amrex::Abort("Unknown cvode_meth");
    if (cvode_itmeth < 1)
      amrex::Abort("Unknown cvode_itmeth");

    amrex::Print() << "This is AMReX version " << amrex::Version() << std::endl;
    amrex::Print() << "Problem domain size: nx = ny = nz = " << n_cell << std::endl;
    amrex::Print() << "Max grid size: " << max_grid_size << std::endl;
    amrex::Print() << "CVODE method: ";
    if (cvode_meth == 1) {
      amrex::Print() << "Adams (non-stiff)";
    } else if (cvode_meth == 2) {
        amrex::Print() << "BDF (stiff)";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "CVODE iteration method: ";
    if (cvode_itmeth == 1) {
      amrex::Print() << "Functional";
    } else if (cvode_itmeth == 2) {
        amrex::Print() << "Newton";
    }
    amrex::Print() << std::endl;

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
      IntVect dom_lo(IntVect(D_DECL(0,0,0)));
      IntVect dom_hi(IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
      Box domain(dom_lo, dom_hi);

      // Initialize the boxarray "ba" from the single box "bx"
      ba.define(domain);

      // Break up boxarray "ba" into chunks no larger than "max_grid_size"
      // along a direction
      ba.maxSize(max_grid_size);

      // This defines the physical size of the box.  Right now the box is
      // [-1,1] in each direction.
      RealBox real_box;
      for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n,-1.0);
        real_box.setHi(n, 1.0);
      }

      // This sets the boundary conditions to be doubly or triply periodic
      int is_periodic[BL_SPACEDIM];
      for (int i = 0; i < BL_SPACEDIM; i++) {
        is_periodic[i] = 1;
      }

      // This defines a Geometry object
      geom.define(domain,&real_box,CoordSys::cartesian,is_periodic);
    }

    // Ncomp = number of components for each array
    int Ncomp  = 3;

    // time = starting time in the simulation
    Real time = 0.0;

    DistributionMapping dm(ba);

    // Create MultiFab with no ghost cells.
    MultiFab mf(ba, dm, Ncomp, 0);

    amrex::Print() << "Solving ODEs without a user-supplied Jacobian ..." << std::endl;
    strt_time = amrex::second();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(mf, do_tiling); mfi.isValid(); ++mfi )
    {
      const Box& tbx = mfi.tilebox();

      integrate_ode_no_jac(mf[mfi].dataPtr(),
        tbx.loVect(),
        tbx.hiVect(),
        &Ncomp,
        &cvode_meth,
        &cvode_itmeth);
    }
    stop_time = amrex::second();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);
    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Time = " << stop_time - strt_time << std::endl;

    if (write_plotfile)
    {
      amrex::WriteSingleLevelPlotfile("PLT_OUTPUT_NO_JAC",
                                      mf,
                                      {"y1", "y2", "y3"},
                                      geom,
                                      time,
                                      0);
    }

    amrex::Print() << "Solving ODEs with a user-supplied Jacobian ..." << std::endl;
    strt_time = amrex::second();
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(mf, do_tiling); mfi.isValid(); ++mfi )
    {
      const Box& tbx = mfi.tilebox();

      integrate_ode_with_jac(mf[mfi].dataPtr(),
        tbx.loVect(),
        tbx.hiVect(),
        &Ncomp,
        &cvode_meth,
        &cvode_itmeth);
    }
    stop_time = amrex::second();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);
    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Time = " << stop_time - strt_time << std::endl;

    if (write_plotfile)
    {
      amrex::WriteSingleLevelPlotfile("PLT_OUTPUT_WITH_JAC",
                                      mf,
                                      {"y1", "y2", "y3"},
                                      geom,
                                      time,
                                      0);
    }

    amrex::Finalize();
    return 0;
}
