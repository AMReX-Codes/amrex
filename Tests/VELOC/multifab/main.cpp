#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>

#include <thread>
#include <future>

#ifdef AMREX_USE_VELOC
#include <veloc.h>
#endif

using namespace amrex;

void main_main ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

#ifdef AMREX_USE_VELOC
    std::string veloc_cfg("veloc.cfg");
    {
        ParmParse pp("veloc");
        pp.get("config_file", veloc_cfg);
    }
    const auto veloc_status = VELOC_Init(ParallelDescriptor::Communicator(),
                                         veloc_cfg.c_str());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(veloc_status == VELOC_SUCCESS, "VELOC_Init failed");
#endif

    main_main();

#ifdef AMREX_USE_VELOC
    VELOC_Finalize(0); // no clean up
#endif

    amrex::Finalize();
}

#ifdef AMREX_USE_VELOC
template <class T>
struct VeloCDeleter {
    void operator() (Vector<std::unique_ptr<T,DataDeleter> > ptrs) {
        if (VELOC_SUCCESS == VELOC_Checkpoint_wait()) {
            amrex::Print() << "VELOC_Checkpoint_wait() finished" << std::endl; // Here for testing only
            ptrs.clear();
        } else {
            amrex::Abort("VeloCDeleter failed");
        }
    }
};

std::thread WriteVeloC (MultiFab& mf, std::string const& mf_name, int step)
{
    if (amrex::Verbose() > 0) {
        amrex::Print() << "Writing VeloC " << mf_name << " " << step << "\n";
    }

    BL_PROFILE("WriteVeloC()");
    if (ParallelDescriptor::IOProcessor())
    {
        std::string file_name("veloc");
        ParmParse pp("veloc");
        pp.query("directory", file_name);
        file_name += "/" + amrex::Concatenate(mf_name, step) + "-metadata";
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream ofs;
        ofs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        ofs.open(file_name.c_str(), std::ios::out | std::ios::trunc);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ofs.good(), "Failed to open file");
        ofs << ParallelDescriptor::NProcs() << '\n';
        ofs << mf.nComp() << '\n';
        ofs << mf.nGrowVect() << '\n';
        const auto& ba = mf.boxArray();
        ba.writeOn(ofs);
        ofs << '\n';
        const auto& dm = mf.DistributionMap();
        dm.writeOn(ofs);
        ofs << '\n';
    }

    MultiFab mfcpu(mf.boxArray(), mf.DistributionMap(), mf.nComp(), mf.nGrowVect(),
                   MFInfo().SetArena(The_Pinned_Arena()));
    {
        BL_PROFILE("dtoh_memcpy");
        amrex::dtoh_memcpy(mfcpu, mf);
    }

    Vector<std::unique_ptr<Real,DataDeleter> > ptrs;
    {
        BL_PROFILE("WriteVeloC()-memprotect");
        for (MFIter mfi(mfcpu); mfi.isValid(); ++mfi) {
            auto p = mfcpu[mfi].release();
            VELOC_Mem_protect(mfi.LocalIndex(), p.get(), mfcpu[mfi].size(), sizeof(Real));
            ptrs.push_back(std::move(p));
        }
    }

    BL_PROFILE_VAR("WriteVeloC()-checkpoint", blp);
    const auto veloc_status = VELOC_Checkpoint(mf_name.c_str(), step);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(veloc_status == VELOC_SUCCESS,
                                     "VELOC_Checkpoint failed");

    std::thread t(VeloCDeleter<Real>(), std::move(ptrs));

    if (amrex::Verbose() > 0) {
        amrex::Print() << "WriteVeloC finished." << std::endl;
    }

    return t;
}

MultiFab ReadVeloC (std::string const& mf_name, int step)
{
    if (amrex::Verbose() > 0) {
        amrex::Print() << "Reading VeloC " << mf_name << " " << step << "\n";
    }

    BL_PROFILE("ReadVeloC()");

    int ncomp;
    IntVect ngrow;
    BoxArray ba;
    DistributionMapping dm;
    {
        std::string file_name("veloc");
        ParmParse pp("veloc");
        pp.query("directory", file_name);
        file_name += "/" + amrex::Concatenate(mf_name, step) + "-metadata";
        Vector<char> medata_buffer;
        ParallelDescriptor::ReadAndBcastFile(file_name, medata_buffer);
        std::istringstream iss(medata_buffer.data(), std::istringstream::in);

        int nprocs;
        iss >> nprocs;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nprocs == ParallelDescriptor::NProcs(),
                                         "ReadVeloc: must use the same number of MPI processes");
        iss >> ncomp;
        iss >> ngrow;
        ba.readFrom(iss);
        dm.readFrom(iss);
    }

    MultiFab mf(ba, dm, ncomp, ngrow);

    BL_PROFILE_VAR("ReadVeloC()-data",blp_data);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto p = mf[mfi].dataPtr();
        VELOC_Mem_protect(mfi.LocalIndex(), p, mf[mfi].size(), sizeof(Real));
    }

    const auto veloc_status = VELOC_Restart(mf_name.c_str(), step);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(veloc_status == VELOC_SUCCESS,
                                     "VeloC_Restart failed");

    return mf;
}
#endif

void main_main ()
{
    BL_PROFILE("main");

    int n_cell = 512;
    int max_grid_size = 64;
    std::string check_file("chk");
    int restart_step = -1;
    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("max_grid_size", max_grid_size);
        pp.query("check_file", check_file);
        pp.query("restart_step", restart_step);
    }

    if (restart_step < 0)
    {
        BoxArray ba(Box(IntVect(0),IntVect(n_cell-1)));
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        MultiFab mf(ba, dm, 1, 0);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const auto& arr = mf.array(mfi);
            amrex::CheckSeedArraySizeAndResize(bx.numPts());
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

        amrex::UtilCreateDirectoryDestructive("vismfdata");

        const int nwork = 50;

        {
            BL_PROFILE_REGION("vismf-orig");
            VisMF::Write(mf, "vismfdata/mf1");
            {
                BL_PROFILE_VAR("vismf-orig-work", blp2);
                for (int i = 0; i < nwork; ++i) {
                    amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                    amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                }
            }
        }

        WriteAsyncStatus status;
        {
            BL_PROFILE_REGION("vismf-async-overlap");
            auto wrt_future = VisMF::WriteAsync(mf, "vismfdata/mf2");
            {
                BL_PROFILE_VAR("vismf-async-work", blp2);
                for (int i = 0; i < nwork; ++i) {
                    amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                    amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                }
            }
            {
                BL_PROFILE_VAR("vismf-async-wait", blp3);
                wrt_future.wait();
                status = wrt_future.get();
            }
        }
        for (int ip = 0; ip < ParallelDescriptor::NProcs(); ++ip) {
            if (ip == ParallelDescriptor::MyProc()) {
                amrex::AllPrint() << "Proc. " << ip << ": " << status << std::endl;
            }
            ParallelDescriptor::Barrier();
        }

#ifdef AMREX_USE_VELOC
        {
            BL_PROFILE_REGION("VeloCTotal");
            auto t = WriteVeloC(mf, check_file.c_str(), 0);

            {
                BL_PROFILE_VAR("VeloC-work", blp2);
                for (int i = 0; i < nwork; ++i) {
                    amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                    amrex::Print() << "mf min = " << mf.min(0) << ",  max = " << mf.max(0) << "\n";
                }
            }
            {
                BL_PROFILE_VAR("VeloCWait", blp3);
                t.join();
            }
        }
#endif
    }
    else
    {
#ifdef AMREX_USE_VELOC
        MultiFab mf_veloc = ReadVeloC(check_file, restart_step);
        amrex::prefetchToDevice(mf_veloc);

        MultiFab mf_vismf(mf_veloc.boxArray(), mf_veloc.DistributionMap(),
                          mf_veloc.nComp(), mf_veloc.nGrowVect());
        VisMF::Read(mf_vismf, "vismfdata/mf");

        MultiFab::Subtract(mf_vismf, mf_veloc, 0, 0, mf_veloc.nComp(), mf_veloc.nGrowVect());
        Real dmin = mf_vismf.min(0);
        Real dmax = mf_vismf.max(0);
        amrex::Print() << "diff min = " << dmin << ", max = " << dmax << "\n";
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dmin == 0.0 and dmax == 0.0,
                                         "Restart failed.");
#endif
    }
}

