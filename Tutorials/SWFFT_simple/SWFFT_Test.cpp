
#include <SWFFT_Test.H>
#include <SWFFT_Test_F.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#include <AMReX_PlotFileUtil.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#define ALIGN 16

using namespace amrex;

void swfft_compute(MultiFab& phi_spatial, MultiFab& phi, Geometry& geom, int verbose);

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
	pp.query("prob_type", prob_type);
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

    // Make sure we define both the phi and the phi_spatial with the same DistributionMapping 
    DistributionMapping dmap{ba};

    // Note that we are defining phi_spatial with NO ghost cells
    phi_spatial.define(ba, dmap, 1, 0);
    init_phi_spatial();

    // Note that we are defining phi with NO ghost cells
    phi.define(ba, dmap, 1, 0);

}

void
SWFFT_Test::init_phi_spatial ()
{
    const Real* dx = geom.CellSize();
    Box domain(geom.Domain());

    for (MFIter mfi(phi_spatial,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        fort_init_phi_spatial(BL_TO_FORTRAN_BOX(tbx),
                      BL_TO_FORTRAN_ANYD(phi_spatial[mfi]),
                      BL_TO_FORTRAN_BOX(domain),
                      geom.ProbLo(), geom.ProbHi(), dx, &prob_type);
    }
}

void
SWFFT_Test::computeFFT ()
{
    Real start_time = amrex::second();
    swfft_compute(phi_spatial, phi, geom, verbose);
    Real total_time = amrex::second() - start_time;

    WritePlotFile();
    bool write_data = false;

    if (write_data)
    {
       VisMF::Write(phi_spatial,"PHI_SPATIAL");
       VisMF::Write(phi,"PHI_DFT");
    }

    if (verbose)
    {
       amrex::Print() << "MAX / MIN VALUE OF DFT " <<  phi.max(0) << " " 
                      <<  phi.min(0) << std::endl;
    }

    amrex::Print() << "Time spent taking FFT: " << total_time << std::endl;
}

void
SWFFT_Test::WritePlotFile (const int step, const amrex::Real time)
{
    
    MultiFab plotfile;
    Vector<std::string> varNames;
    int nPlot;

    //////////////////////////////////////////////////////////////////////////////////
    // Write out FFT to plot file
    //////////////////////////////////////////////////////////////////////////////////
    const std::string plotfilename = "plt_fft";
    nPlot = 2;
    plotfile.define(phi.boxArray(), phi.DistributionMap(), nPlot, 0);
    varNames.resize(nPlot);

    // keep a counter for plotfile variables
    int cnt = 0;

    varNames[cnt++] = "FFT_of_phi";
    varNames[cnt++] = "phi";

    // reset plotfile variable counter
    cnt = 0;

    // copy into plotfile
    MultiFab::Copy(plotfile, phi, 0, cnt, 1, 0);
    cnt++;

    MultiFab::Copy(plotfile, phi_spatial,  0, cnt, 1, 0);
    cnt++;

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);
}
