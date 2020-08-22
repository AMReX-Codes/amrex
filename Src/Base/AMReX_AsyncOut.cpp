#include <AMReX_AsyncOut.H>
#include <AMReX_BackgroundThread.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX.H>

namespace amrex {
namespace AsyncOut {

namespace {

#if defined(AMREX_USE_DPCPP) || defined(AMREX_USE_HIP)
int s_asyncout = true; // Have this on by default for DPC++ for now so that
                       // I/O writing plotfile does not depend on unified
                       // memory.
#else
int s_asyncout = false;
#endif
int s_noutfiles = 64;
MPI_Comm s_comm = MPI_COMM_NULL;

std::unique_ptr<BackgroundThread> s_thread;

WriteInfo s_info;

}

void Initialize ()
{
    amrex::ignore_unused(s_comm,s_info);

    ParmParse pp("amrex");
    pp.query("async_out", s_asyncout);
    pp.query("async_out_nfiles", s_noutfiles);

    int nprocs = ParallelDescriptor::NProcs();
    s_noutfiles = std::min(s_noutfiles, nprocs);

    if (s_asyncout and s_noutfiles < nprocs)
    {
#ifdef AMREX_MPI_THREAD_MULTIPLE
        int myproc = ParallelDescriptor::MyProc();
        s_info = GetWriteInfo(myproc);
        MPI_Comm_split(ParallelDescriptor::Communicator(), s_info.ifile, myproc, &s_comm);
#else
        amrex::Abort("AsyncOut with " + std::to_string(s_noutfiles) + " and "
                     +std::to_string(nprocs) + " processes requires MPI_THREAD_MULTIPLE");
#endif
    }

    if (s_asyncout) s_thread.reset(new BackgroundThread());

    ExecOnFinalize(Finalize);
}

void Finalize ()
{
    if (s_thread) {
        s_thread.reset();
    }

#ifdef AMREX_USE_MPI
    if (s_comm != MPI_COMM_NULL) MPI_Comm_free(&s_comm);
    s_comm = MPI_COMM_NULL;
#endif
}

bool UseAsyncOut () { return s_asyncout; }

WriteInfo GetWriteInfo (int rank)
{
    const int nfiles = s_noutfiles;
    const int nprocs = ParallelDescriptor::NProcs();
    const int nmaxspots = (nprocs + (nfiles-1)) / nfiles;  // max spots per file
    const int nfull = nfiles + nprocs - nmaxspots*nfiles;  // the first nfull files are full

    int ifile, ispot, nspots;
    if (rank < nfull*nmaxspots) {
        ifile = rank / nmaxspots;
        ispot = rank - ifile*nmaxspots;
        nspots = nmaxspots;
    } else {
        int tmpproc = rank-nfull*nmaxspots;
        ifile = tmpproc/(nmaxspots-1);
        ispot = tmpproc - ifile*(nmaxspots-1);
        ifile += nfull;
        nspots = nmaxspots - 1;
    }

    return WriteInfo{ifile, ispot, nspots};
}

void Submit (std::function<void()>&& a_f)
{
    s_thread->Submit(std::move(a_f));
}

void Submit (std::function<void()> const& a_f)
{
    s_thread->Submit(a_f);
}

void Finish ()
{
    s_thread->Finish();
}

void Wait ()
{
#ifdef AMREX_USE_MPI
    const int N = s_info.ispot;
    if (N > 0) {
        Vector<MPI_Request> reqs(N);
        Vector<MPI_Status> stats(N);
        for (int i = 0; i < N; ++i) {
            reqs[i] = ParallelDescriptor::Abarrier(s_comm).req();
        }
        ParallelDescriptor::Waitall(reqs, stats);
    }
#endif
}

void Notify ()
{
#ifdef AMREX_USE_MPI
    const int N = s_info.nspots - 1 - s_info.ispot;
    if (N > 0) {
        Vector<MPI_Request> reqs(N);
        Vector<MPI_Status> stats(N);
        for (int i = 0; i < N; ++i) {
            reqs[i] = ParallelDescriptor::Abarrier(s_comm).req();
        }
        ParallelDescriptor::Waitall(reqs, stats);
    }
#endif
}

}}
