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
    int nwrites = 1;

    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.queryarr("n_cell_3d", n_cell_3d, 0, AMREX_SPACEDIM);
        pp.query("n_boxes_per_rank", n_boxes_per_rank);
        pp.query("max_grid_size", max_grid_size);
        pp.query("noutfiles", n_files);
        pp.query("nwork", nwork);
        pp.query("nwrites", nwrites);

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

    Vector<MultiFab> mfs(nwrites);
    Vector< Array4<Real> > arrs(nwrites);

    for (int m = 0; m < nwrites; ++m) {
        mfs[m].define(ba, dm, 1, 0);   
    }

//    amrex::ResetRandomSeed(33344455666);
    for (MFIter mfi(mfs[0]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        for (int m = 0; m < nwrites; ++m) {
            arrs[m] = mfs[m].array(mfi);
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int m = 0; m < nwrites; ++m) {
                arrs[m](i,j,k) = amrex::Random();
            }
        });
        Gpu::streamSynchronize(); // because of random nubmer generator
    }

    { // touch pinned memory so that the one-time cost is removed from timers.
        Vector<MultiFab> mfcpus(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            mfcpus[m].define(mfs[m].boxArray(), mfs[m].DistributionMap(),
                             mfs[m].nComp(), mfs[m].nGrowVect(),
                             MFInfo().SetArena(The_Pinned_Arena()));
            amrex::dtoh_memcpy((mfcpus[m]), (mfs[m]));
        }
    }

    amrex::Print() << "I/O printing randomly filled multifab with: "
                   << "\n  dimensions = "    << ba.minimalBox() 
                   << "\n  max_grid_size = " << max_grid_size
                   << "\n  noutfiles = "     << n_files
                   << "\n  boxes = "         << ba.size()
                   << "\n  and nwork = "     << nwork << std::endl;

    Vector<Real> mf_min(nwrites), mf_max(nwrites);

    for (int m = 0; m < nwrites; ++m) {
        mf_min[m] = mfs[m].min(0);
        mf_max[m] = mfs[m].max(0);
    }

    for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
        if (ip == ParallelDescriptor::MyProc()) {
            amrex::AllPrint() << "Proc. " << ip << " number of boxes = " << mfs[0].local_size() << std::endl;
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
            for (int m = 0; m < nwrites; ++m) {
                VisMF::Write(mfs[m], std::string("vismfdata/time-" + std::to_string(m)));
            }
            ParallelDescriptor::Barrier();
        }
        {
            BL_PROFILE_VAR("vismf-time-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " No Async " << std::endl;
    {
        BL_PROFILE_REGION("vismf-orig");
        for (int m = 0; m < nwrites; ++m) {
            VisMF::Write(mfs[m], std::string("vismfdata/plain-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-orig-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-file " << std::endl; 
    Vector<WriteAsyncStatus> status_file(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-file-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsync(mfs[m], std::string("vismfdata/file-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-file-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-file-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_file[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

#ifdef AMREX_MPI_MULTIPLE

    amrex::Print() << " Async-MPI " << std::endl; 
    Vector<WriteAsyncStatus> status_mpi_basic(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-mpi-basic-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsyncMPI(mfs[m], std::string("vismfdata/mpi-basic-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-basic-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-basic-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_mpi_basic[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Comm" << std::endl; 
    Vector<WriteAsyncStatus> status_mpi_comm(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-mpi-comm-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsyncMPIComm(mfs[m], std::string("vismfdata/mpi-comm-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-comm-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-comm-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_mpi_comm[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Wait" << std::endl; 
    Vector<WriteAsyncStatus> status_mpi_wait(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-mpi-wait-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsyncMPIWait(mfs[m], std::string("vismfdata/mpi-wait-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-wait-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-wait-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_mpi_wait[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Barrier" << std::endl; 
    Vector<WriteAsyncStatus> status_mpi_barrier(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-mpi-barrier-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsyncMPIBarrier(mfs[m], std::string("vismfdata/mpi-barrier-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-barrier-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-barrier-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_mpi_barrier[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

   amrex::Print() << " Async-MPI IBarrier" << std::endl; 
    Vector<WriteAsyncStatus> status_mpi_ibarrier(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-mpi-ibarrier-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsyncMPIABarrier(mfs[m], std::string("vismfdata/ibarrier-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-ibarrier-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-ibarrier-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_mpi_ibarrier[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

   amrex::Print() << " Async-MPI IBarrier Waitall" << std::endl;
    Vector<WriteAsyncStatus> status_mpi_ibarrier_waitall(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-mpi-ibarrier-waitall-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsyncMPIABarrierWaitall(mfs[m], std::string("vismfdata/ibarrier-waitall-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-ibarrier-waitall-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-ibarrier-waitall-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_mpi_ibarrier_waitall[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Fence" << std::endl;
    Vector<WriteAsyncStatus> status_mpi_fence(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-mpi-fence-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsyncMPIOneSidedFence(mfs[m], std::string("vismfdata/fence-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-fence-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-fence-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_mpi_fence[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

    amrex::Print() << " Async-MPI Post" << std::endl;
    Vector<WriteAsyncStatus> status_mpi_post(nwrites);
    {
        BL_PROFILE_REGION("vismf-async-mpi-post-overlap");
        Vector<std::future<WriteAsyncStatus> > wrt_futures(nwrites);
        for (int m = 0; m < nwrites; ++m) {
            wrt_futures[m] = VisMF::WriteAsyncMPIOneSidedPost(mfs[m], std::string("vismfdata/post-" + std::to_string(m)));
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-post-work", blp2);
            for (int m = 0; m < nwrites; ++m) {
                for (int i = 0; i < nwork*2; ++i) {
                    Real min = mfs[m].min(0);
                    Real max = mfs[m].max(0);
                    if (mf_min[m] != min)
                        { amrex::AllPrint() << "Min failed: " << min << " != " << mf_min[m] << std::endl; }
                    if (mf_max[m] != max)
                        { amrex::AllPrint() << "Max failed: " << max << " != " << mf_max[m] << std::endl; }
                }
            }
        }
        {
            BL_PROFILE_VAR("vismf-async-mpi-post-wait", blp3);
            for (int m = 0; m < nwrites; ++m) {
                wrt_futures[m].wait();
                status_mpi_post[m] = wrt_futures[m].get();
            }
        }
    }
    ParallelDescriptor::Barrier();

// ***************************************************************

#endif

    for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
        if (ip == ParallelDescriptor::MyProc()) {
            amrex::AllPrint() << " ====== " << std::endl;
            amrex::AllPrint() << " Proc. " << ip << std::endl;

            for (int m = 0; m < nwrites; ++m)
            {
                amrex::AllPrint() << " MultiFab #" << m << std::endl;
                amrex::AllPrint() << " File: " << status_file[m] << std::endl;
#ifdef AMREX_MPI_MULTIPLE
                amrex::AllPrint() << " MPI-Basic: "           << status_mpi_basic[m]             << std::endl;
                amrex::AllPrint() << " MPI-Comm: "            << status_mpi_comm[m]              << std::endl;
                amrex::AllPrint() << " MPI-Wait: "            << status_mpi_wait[m]              << std::endl;
                amrex::AllPrint() << " MPI-Barrier: "         << status_mpi_barrier[m]           << std::endl;
                amrex::AllPrint() << " MPI-IBarrier: "        << status_mpi_ibarrier[m]          << std::endl;
                amrex::AllPrint() << " MPI-IBarrierWaitall: " << status_mpi_ibarrier_waitall[m]  << std::endl;
                amrex::AllPrint() << " MPI-Fence: "           << status_mpi_fence[m]             << std::endl;
                amrex::AllPrint() << " MPI-Post: "            << status_mpi_post[m]              << std::endl;
                amrex::AllPrint() << std::endl;
#endif
            }
        }
        amrex::USleep(0.001);
        ParallelDescriptor::Barrier();
    }



}
