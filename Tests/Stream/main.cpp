#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

void main_main ()
{
    // This test is not configured to use multiple ranks.
    if (ParallelDescriptor::NProcs() > 1) {
        amrex::Abort("This test is not configured to use multiple ranks.");
    }

    int n_cell, max_grid_size, nsteps, nvar;
    std::string test;

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // Number of cells on each side of a cubic domain.
        pp.get("n_cell", n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size", max_grid_size);

        // Number of "steps" to take by default.
        nsteps = 10;
        pp.query("nsteps", nsteps);

        // Number of variables to represent.
        nvar = 1;
        pp.query("nvar", nvar);

        // Stream test to run.
        test = "triad";
        pp.query("test", test);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "bx"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [-1,1] in each direction.
    RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                     {AMREX_D_DECL( 1.0,  1.0,  1.0)});

    // Specify periodic domain (does not matter for this test).
    Vector<int> is_periodic(AMREX_SPACEDIM, 1);

    // This defines a Geometry object
    geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // Allocate and initialize our data, stored as MultiFabs.
    MultiFab a(ba, dm, nvar, 0);
    MultiFab b(ba, dm, nvar, 0);
    MultiFab c(ba, dm, nvar, 0);

    a.setVal(0.0);
    b.setVal(1.0);
    c.setVal(2.0);

    Real scal = 3.0;

    amrex::Print() << std::endl;
    amrex::Print() << "Beginning stream benchmark with: " << std::endl;
    amrex::Print() << "n_cell = " << n_cell << std::endl;
    amrex::Print() << "nvar = " << nvar << std::endl;
    amrex::Print() << "nsteps = " << nsteps << std::endl;
    amrex::Print() << "max_grid_size = " << max_grid_size << std::endl;
    amrex::Print() << std::endl;

        // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    for (int s = 1; s <= nsteps; ++s)
    {
        for (MFIter mfi(a, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto a_arr = a.array(mfi);
            auto b_arr = b.array(mfi);
            auto c_arr = c.array(mfi);

            if (test == "copy") {

                AMREX_PARALLEL_FOR_4D(bx, nvar, i, j, k, n,
                {
                    a_arr(i,j,k,n) = b_arr(i,j,k,n);
                });

            }
            else if (test == "scale") {

                AMREX_PARALLEL_FOR_4D(bx, nvar, i, j, k, n,
                {
                    a_arr(i,j,k,n) = scal * b_arr(i,j,k,n);
                });

            }
            else if (test == "sum") {

                AMREX_PARALLEL_FOR_4D(bx, nvar, i, j, k, n,
                {
                    a_arr(i,j,k,n) = b_arr(i,j,k,n) + c_arr(i,j,k,n);
                });

            }
            else if (test == "triad") {

                AMREX_PARALLEL_FOR_4D(bx, nvar, i, j, k, n,
                {
                    a_arr(i,j,k,n) = b_arr(i,j,k,n) + scal * c_arr(i,j,k,n);
                });

            }
            else {

                amrex::Abort("Unknown stream test");

            }
        }
    }

    // Call the timer again, compute maximum run time across ranks.
    Real run_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(run_time);

    amrex::Print() << "Stream benchmark completed!" << std::endl;
    amrex::Print() << std::endl;

    // Write out runtime and effective average bandwidth.

    int nBytesScale = 1;
    int flopScale = 1;

    if (test == "copy") {
        nBytesScale = 2;
        flopScale = 0;
    }
    else if (test == "scale") {
        nBytesScale = 2;
        flopScale = 1;
    }
    else if (test == "sum") {
        nBytesScale = 3;
        flopScale = 1;
    }
    else if (test == "triad") {
        nBytesScale = 3;
        flopScale = 2;
    }

    size_t nBytes = nvar * domain.numPts() * nBytesScale * sizeof(Real);
    Real nBytesGB = nBytes / ((Real) 1024*1024*1024);

    Real bandwidthGBs = nBytesGB / (run_time / nsteps);

    Real flops = nvar * domain.numPts() * flopScale;
    Real Gflops = flops / ((Real) (1024*1024*1024));

    Real flopsGFs = Gflops / (run_time / nsteps);

    amrex::Print() << "Run time = " << run_time << " s" << std::endl;
    amrex::Print() << "Data size = " << nBytesGB << " GB" << std::endl;
    amrex::Print() << "Bandwidth = " << bandwidthGBs << " GB/s" << std::endl;
    amrex::Print() << "Compute = " << flopsGFs << " GF/s" << std::endl;
    amrex::Print() << std::endl;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        main_main();
    }
    amrex::Finalize();
    return 0;
}
