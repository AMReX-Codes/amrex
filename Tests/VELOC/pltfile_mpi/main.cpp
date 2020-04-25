#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_PlotFileUtil.H>

#include <thread>
#include <future>

using namespace amrex;

void main_main ();

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();

}

void main_main ()
{
    BL_PROFILE("main");

    int n_cell = 0, n_ghost = 0, n_comp = 1, n_boxes_per_rank = 0;
    amrex::Vector<int> n_cell_3d (AMREX_SPACEDIM, 512);
    
    int max_grid_size = 64;
    int n_files = 256;
    int nwork = 10;

    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.queryarr("n_cell_3d", n_cell_3d, 0, AMREX_SPACEDIM);
        pp.query("n_boxes_per_rank", n_boxes_per_rank);
        pp.query("max_grid_size", max_grid_size);
        pp.query("nghost", n_ghost);
        pp.query("ncomp", n_comp); 
        pp.query("noutfiles", n_files);
        pp.query("nwork", nwork);

        // inputs hierarchy: 
        // n_cell > n_boxes_per_rank > n_cell_3d

        if (n_cell != 0)
        {
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
            { n_cell_3d[i] = n_cell; }
        }
        else if (n_boxes_per_rank != 0)
        {
           n_cell_3d[0] = (max_grid_size) - 1;
           n_cell_3d[1] = (max_grid_size * n_boxes_per_rank) - 1;
           n_cell_3d[2] = (max_grid_size * ParallelDescriptor::NProcs()) - 1;
        }
    }

    VisMF::SetNOutFiles(n_files);
    Box domain(IntVect(0), IntVect(n_cell_3d));

    BoxArray ba;
    Geometry geom;
    {
        ba.define(domain);
        ba.maxSize(max_grid_size);

        RealBox real_box (0.0, 0.0, 0.0, 
                          1.0, Real(n_boxes_per_rank), Real(ParallelDescriptor::NProcs()));
        geom.define(domain, &real_box);
    }

    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, n_comp, n_ghost);

//    amrex::ResetRandomSeed(33344455666);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto& arr = mf.array(mfi);
        amrex::ParallelFor (bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k) = amrex::Random();
        });
        Gpu::streamSynchronize(); // because of random nubmer generator
    }

    { // touch pinned memory so that the one-time cost is removed from timers.
        MultiFab mfcpu(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect(),
                       MFInfo().SetArena(The_Pinned_Arena()));
        amrex::dtoh_memcpy(mfcpu, mf);
    }

    amrex::Print() << "I/O printing randomly filled plotfile with: "
                   << "\n  dimensions = "    << ba.minimalBox() 
                   << "\n  max_grid_size = " << max_grid_size
                   << "\n  nghost = "        << n_ghost
                   << "\n  noutfiles = "     << n_files
                   << "\n  boxes = "         << ba.size()
                   << "\n  and nwork = "     << nwork << std::endl;

    double mf_min = mf.min(0);
    double mf_max = mf.max(0);

// ***************************************************************

    amrex::Print() << " Time Write and Work separately. " << std::endl;
    {
        BL_PROFILE_REGION("pltfile-time");
        {
            BL_PROFILE_VAR("pltfile-time-write", t1);
            const std::string& pltfile = amrex::Concatenate("PLTtime",0,5);
            WriteSingleLevelPlotfile(pltfile, mf, {"random"}, geom, 0, 0);
            ParallelDescriptor::Barrier();
        }
        {
            BL_PROFILE_VAR("pltfile-time-work", t2);
            for (int i = 0; i < nwork*2; ++i) {
                double min = mf.min(0);
                double max = mf.max(0);
                if (mf_min != min)
                    { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min << std::endl; }
                if (mf_max != max)
                    { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max << std::endl; }
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Standard without Barriers. " << std::endl;
    {
        BL_PROFILE_REGION("pltfile-standard");
        const std::string& pltfile = amrex::Concatenate("PLTstandard",0,5);
        WriteSingleLevelPlotfile(pltfile, mf, {"random"}, geom, 0, 0);
        {
            BL_PROFILE_VAR("pltfile-standard-work", t3);
            for (int i = 0; i < nwork*2; ++i) {
                double min = mf.min(0);
                double max = mf.max(0);
                if (mf_min != min)
                    { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min << std::endl; }
                if (mf_max != max)
                    { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max << std::endl; }
            }
        }
    }
    ParallelDescriptor::Barrier();

#ifdef AMREX_MPI_MULTIPLE

// ***************************************************************

    // For timing of Async Methodology
    amrex::Print() << " Async Plotfile w/ Immediate wait. " << std::endl; 
    WriteAsyncStatus status_pltfile_async_timing;
    {
        BL_PROFILE_REGION("pltfile-async-timing");
        std::future<WriteAsyncStatus> wrt_future;

        {
            BL_PROFILE_VAR("pltfile-async-timing-launch", t4);
            const std::string& pltfile = amrex::Concatenate("PLTasyncTest",0,5);
            WriteAsyncSingleLevelPlotfile(pltfile, mf, {"random"}, geom, 0, 0);
        }
        {
//            BL_PROFILE_VAR("pltfile-async-timing-wait", t5);
//            status_pltfile_async_timing = VisMF::asyncWaitAll();
        }
    }


// ***************************************************************

    amrex::Print() << " Async Plotfile. " << std::endl; 
    WriteAsyncStatus status_pltfile_async;
    {
        BL_PROFILE_REGION("pltfile-async");
        std::future<WriteAsyncStatus> wrt_future;

        {
            BL_PROFILE_VAR("pltfile-async-launch", t4);
            const std::string& pltfile = amrex::Concatenate("PLTasync",0,5);
            WriteAsyncSingleLevelPlotfile(pltfile, mf, {"random"}, geom, 0, 0);
        }
        {
            BL_PROFILE_VAR("pltfile-async-work", t3);
            for (int i = 0; i < nwork*2; ++i) {
                double min = mf.min(0);
                double max = mf.max(0);
                if (mf_min != min)
                    { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min << std::endl; }
                if (mf_max != max)
                    { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max << std::endl; }
            }
        }
        {
//            BL_PROFILE_VAR("pltfile-async-wait", t5);
//            status_pltfile_async = VisMF::asyncWaitAll(); 
        }
    }

// ***************************************************************

#endif


#ifdef AMREX_MPI_MULTIPLE
    for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
        if (ip == ParallelDescriptor::MyProc()) {
            amrex::AllPrint() << "Proc. " << ip << std::endl;

            amrex::AllPrint() << "MPI-Async-Timing: " << status_pltfile_async_timing << std::endl;
            amrex::AllPrint() << "MPI-Async: "        << status_pltfile_async        << std::endl;
        }
        amrex::USleep(0.001);
        ParallelDescriptor::Barrier();
    }
#endif

}
