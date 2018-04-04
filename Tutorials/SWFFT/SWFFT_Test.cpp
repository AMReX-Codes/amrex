
#include <SWFFT_Test.H>
#include <SWFFT_Test_F.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#define ALIGN 16

using namespace amrex;

void swfft_solver(MultiFab& rhs, MultiFab& soln, Geometry& geom, int verbose);

SWFFT_Test::SWFFT_Test ()
{
    static_assert(AMREX_SPACEDIM == 3, "3D only");

    // runtime parameters
    {
        ParmParse pp;

        // Read in n_cell.  Use defaults if not explicitly defined.
        int cnt = pp.countval("n_cell");
        if (cnt > 1) {
            Vector<int> ncs;
            pp.getarr("n_cell",ncs);
            n_cell = IntVect{ncs[0],ncs[1],ncs[2]};
        } else if (cnt > 0) {
            int ncs;
            pp.get("n_cell",ncs);
            n_cell = IntVect{ncs,ncs,ncs};
        } else {
           n_cell = IntVect{32,32,32};
        }

        // Read in max_grid_size.  Use defaults if not explicitly defined.
        cnt = pp.countval("max_grid_size");
        if (cnt > 1) {
            Vector<int> mgs;
            pp.getarr("max_grid_size",mgs);
            max_grid_size = IntVect{mgs[0],mgs[1],mgs[2]};
        } else if (cnt > 0) {
            int mgs;
            pp.get("max_grid_size",mgs);
            max_grid_size = IntVect{mgs,mgs,mgs};
        } else {
           max_grid_size = IntVect{32,32,32};
        }

        pp.query("verbose", verbose);
    }
    
    BoxArray ba;
    {
        // Make up a dx that is not 1
        Real dx = 1./double(n_cell[0]);

        IntVect dom_lo(0,0,0);
        IntVect dom_hi(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1);
        Box domain(dom_lo,dom_hi);
        ba.define(domain);
        ba.maxSize(max_grid_size);

        // We assume a unit box (just for convenience)
        Real x_hi = n_cell[0]*dx;
        Real y_hi = n_cell[1]*dx;
        Real z_hi = n_cell[2]*dx;
        RealBox real_box({0.0,0.0,0.0}, {x_hi,y_hi,z_hi});

        // The FFT assumes fully periodic boundaries
        std::array<int,3> is_periodic {1,1,1};

        geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());
    }

    // Make sure we define both the soln and the rhs with the same DistributionMapping 
    DistributionMapping dmap{ba};

    // Note that we are defining rhs with NO ghost cells
    rhs.define(ba, dmap, 1, 0);
    init_rhs();

    // Note that we are defining soln with NO ghost cells
    soln.define(ba, dmap, 1, 0);
    the_soln.define(ba, dmap, 1, 0);

    comp_the_solution();
}

void
SWFFT_Test::init_rhs ()
{
    const Real* dx = geom.CellSize();
    Box domain(geom.Domain());

    for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        fort_init_rhs(BL_TO_FORTRAN_BOX(tbx),
                      BL_TO_FORTRAN_ANYD(rhs[mfi]),
                      BL_TO_FORTRAN_BOX(domain),
                      geom.ProbLo(), geom.ProbHi(), dx);
    }

    Real sum_rhs = rhs.sum();
    amrex::Print() << "Sum of rhs over the domain was    " << sum_rhs << std::endl;

         sum_rhs = sum_rhs / domain.numPts();
    rhs.plus(-sum_rhs,0,1,0);

         sum_rhs = rhs.sum();
    amrex::Print() << "Sum of rhs over the domain is now " << sum_rhs << std::endl;
}

void
SWFFT_Test::comp_the_solution ()
{
    const Real* dx = geom.CellSize();
    Box domain(geom.Domain());

    for (MFIter mfi(the_soln); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        fort_comp_asol(BL_TO_FORTRAN_BOX(tbx),
                       BL_TO_FORTRAN_ANYD(the_soln[mfi]),
                       BL_TO_FORTRAN_BOX(domain),
                       geom.ProbLo(), geom.ProbHi(), dx);
    }
}

void
SWFFT_Test::solve ()
{
    Real start_time = amrex::second();
    swfft_solver(rhs, soln, geom, verbose);
    Real total_time = amrex::second() - start_time;

    bool write_data = false;

    if (write_data)
    {
       VisMF::Write(rhs,"RHS");
       VisMF::Write(soln,"SOL_COMP");
       VisMF::Write(the_soln,"SOL_EXACT");
    }

    if (verbose)
    {
       amrex::Print() << "MAX / MIN VALUE OF COMP  SOLN " <<     soln.max(0) << " " 
                      <<     soln.min(0) << std::endl;
       amrex::Print() << "MAX / MIN VALUE OF EXACT SOLN " << the_soln.max(0) << " " 
                      << the_soln.min(0) << std::endl;
    }

    BoxArray ba = rhs.boxArray();
    const DistributionMapping& dm = rhs.DistributionMap();
    MultiFab diff(ba, dm, 1, 0);
    MultiFab::Copy(diff, soln, 0, 0, 1, 0);
    MultiFab::Subtract(diff, the_soln, 0, 0, 1, 0);
    amrex::Print() << "\nMax-norm of the error is " << diff.norm0() << "\n";
    amrex::Print() << "Time spent in solve: " << total_time << std::endl;

    if (write_data)
       VisMF::Write(diff,"DIFF");
}
