#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>

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

    int n_cell = 0;
    int n_boxes_per_rank = 0;
    amrex::Vector<int> n_cell_3d (AMREX_SPACEDIM, 512);
    int max_grid_size = 64;

    int n_files = 256;
    int nwork = 50;

    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.queryarr("n_cell_3d", n_cell_3d, 0, AMREX_SPACEDIM);
        pp.query("n_boxes_per_rank", n_boxes_per_rank);
        pp.query("max_grid_size", max_grid_size);
        pp.query("noutfiles", n_files);
        pp.query("nwork", nwork);

        // inputs hierarchy: 
        // n_cell > n_boxes_per_rank > n_cell_3d

        if (n_cell != 0)
        {
            for (int i; i < AMREX_SPACEDIM; ++i)
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
    BoxArray ba(Box(IntVect(0),IntVect(n_cell_3d)));
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, 1, 0);

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

    amrex::Print() << "I/O printing randomly filled multifab with: "
                   << "\n  dimensions = "    << ba.minimalBox() 
                   << "\n  max_grid_size = " << max_grid_size
                   << "\n  noutfiles = "     << n_files
                   << "\n  boxes = "         << ba.size()
                   << "\n  and nwork = "     << nwork << std::endl;

    double mf_min = mf.min(0);
    double mf_max = mf.max(0);

    for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
        if (ip == ParallelDescriptor::MyProc()) {
            amrex::AllPrint() << "Proc. " << ip << " number of boxes = " << mf.local_size() << std::endl;
        }
        amrex::USleep(0.001);
        ParallelDescriptor::Barrier();
    }

    amrex::UtilCreateDirectoryDestructive("vismfdata");

// ***************************************************************

    amrex::Print() << " Time Write and Work separately. " << std::endl;
    {
        BL_PROFILE_REGION("vismf-time");
        {
            BL_PROFILE_VAR("vismf-time-write", blp1);
            VisMF::Write(mf, "vismfdata/mf0");
            ParallelDescriptor::Barrier();
        }
        {
            BL_PROFILE_VAR("vismf-time-work", blp2);
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

    amrex::Print() << " No Async " << std::endl;
    {
        BL_PROFILE_REGION("vismf-orig");
        VisMF::Write(mf, "vismfdata/mf1");
        {
            BL_PROFILE_VAR("vismf-orig-work", blp2);
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

    amrex::Print() << " Async-file " << std::endl; 
    WriteAsyncStatus status_file;
    {
        BL_PROFILE_REGION("vismf-async-file-overlap");
        auto wrt_future = VisMF::WriteAsync(mf, "vismfdata/mf2");
        {
            BL_PROFILE_VAR("vismf-async-file-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-file-wait", blp3);
            wrt_future.wait();
            status_file = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

#ifdef AMREX_MPI_MULTIPLE

    amrex::Print() << " Async-MPI " << std::endl; 
    WriteAsyncStatus status_mpi_basic;
    {
        BL_PROFILE_REGION("vismf-async-mpi-basic-overlap");
        auto wrt_future = VisMF::WriteAsyncMPI(mf, "vismfdata/mf3");
        {
            BL_PROFILE_VAR("vismf-async-mpi-basic-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-mpi-basic-wait", blp3);
            wrt_future.wait();
            status_mpi_basic = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Comm " << std::endl; 
    WriteAsyncStatus status_mpi_comm;
    {
        BL_PROFILE_REGION("vismf-async-mpi-comm-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIComm(mf, "vismfdata/mf4");
        {
            BL_PROFILE_VAR("vismf-async-mpi-comm-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-mpi-comm-wait", blp3);
            wrt_future.wait();
            status_mpi_comm = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Wait " << std::endl; 
    WriteAsyncStatus status_mpi_wait;
    {
        BL_PROFILE_REGION("vismf-async-mpi-wait-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIWait(mf, "vismfdata/mf5");
        {
            BL_PROFILE_VAR("vismf-async-mpi-wait-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-mpi-wait-wait", blp3);
            wrt_future.wait();
            status_mpi_wait = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Barrier w/ Comm " << std::endl; 
    WriteAsyncStatus status_mpi_barrier;
    {
        BL_PROFILE_REGION("vismf-async-mpi-barrier-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIBarrier(mf, "vismfdata/mf6");
        {
            BL_PROFILE_VAR("vismf-async-mpi-barrier-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-mpi-barrier-wait", blp3);
            wrt_future.wait();
            status_mpi_barrier = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI IBarrier w/ Comm " << std::endl; 
    WriteAsyncStatus status_mpi_ibarrier;
    {
        BL_PROFILE_REGION("vismf-async-mpi-ibarrier-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIABarrier(mf, "vismfdata/mf7");
        {
            BL_PROFILE_VAR("vismf-async-mpi-ibarrier-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-mpi-ibarrier-wait", blp3);
            wrt_future.wait();
            status_mpi_ibarrier = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI IBarrier Waitall w/ Comm " << std::endl; 
    WriteAsyncStatus status_mpi_ibarrier_waitall;
    {
        BL_PROFILE_REGION("vismf-async-mpi-ibarrier-waitall-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIABarrier(mf, "vismfdata/mf8");
        {
            BL_PROFILE_VAR("vismf-async-mpi-ibarrier-waitall-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-mpi-ibarrier-waitall-wait", blp3);
            wrt_future.wait();
            status_mpi_ibarrier_waitall = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();


// ***************************************************************

    amrex::Print() << " Async-MPI Fence " << std::endl; 
    WriteAsyncStatus status_mpi_fence;
    {
        BL_PROFILE_REGION("vismf-async-mpi-fence-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIOneSidedFence(mf, "vismfdata/mf9");
        {
            BL_PROFILE_VAR("vismf-async-mpi-fence-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-mpi-fence-wait", blp3);
            wrt_future.wait();
            status_mpi_fence = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Post " << std::endl; 
    WriteAsyncStatus status_mpi_post;
    {
        BL_PROFILE_REGION("vismf-async-mpi-post-overlap");
        auto wrt_future = VisMF::WriteAsyncMPIOneSidedPost(mf, "vismfdata/mf10");
        {
            BL_PROFILE_VAR("vismf-async-mpi-post-work", blp2);
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
            BL_PROFILE_VAR("vismf-async-mpi-post-wait", blp3);
            wrt_future.wait();
            status_mpi_post = wrt_future.get();
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

#endif

    for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
        if (ip == ParallelDescriptor::MyProc()) {
            amrex::AllPrint() << "Proc. " << ip << std::endl;
            amrex::AllPrint() << "File: " << status_file << std::endl;
#ifdef AMREX_MPI_MULTIPLE
            amrex::AllPrint() << "MPI-Basic: "           << status_mpi_basic             << std::endl;
            amrex::AllPrint() << "MPI-Comm: "            << status_mpi_comm              << std::endl;
            amrex::AllPrint() << "MPI-Wait: "            << status_mpi_wait              << std::endl;
            amrex::AllPrint() << "MPI-Barrier: "         << status_mpi_barrier           << std::endl;
            amrex::AllPrint() << "MPI-IBarrier: "        << status_mpi_ibarrier          << std::endl;
            amrex::AllPrint() << "MPI-IBarrierWaitall: " << status_mpi_ibarrier_waitall  << std::endl;
            amrex::AllPrint() << "MPI-Fence: "           << status_mpi_fence             << std::endl;
            amrex::AllPrint() << "MPI-Post: "            << status_mpi_post              << std::endl;
#endif
        }
        amrex::USleep(0.001);
        ParallelDescriptor::Barrier();
    }
}
