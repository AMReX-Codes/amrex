#include <AMReX_AsyncOut.H>
#include <AMReX_BackgroundThread.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX.H>

namespace amrex::AsyncOut {

namespace {

int s_asyncout = false;
int s_noutfiles = 64;
MPI_Comm s_comm = MPI_COMM_NULL;

std::unique_ptr<BackgroundThread> s_thread;

WriteInfo s_info;

}

void Initialize ()
{
    amrex::ignore_unused(s_comm,s_info);

    ParmParse pp("amrex");
    pp.queryAdd("async_out", s_asyncout);
    pp.queryAdd("async_out_nfiles", s_noutfiles);

    int nprocs = ParallelDescriptor::NProcs();
    s_noutfiles = std::min(s_noutfiles, nprocs);

#ifdef AMREX_USE_MPI
    if (s_asyncout && s_noutfiles < nprocs)
    {
        int provided = -1;
        MPI_Query_thread(&provided);
        if (provided < MPI_THREAD_MULTIPLE) {
            amrex::Abort("AsyncOut with " + std::to_string(s_noutfiles) + " and "
                         + std::to_string(nprocs) + " processes requires "
                         + "MPI_THREAD_MULTIPLE at runtime, but got "
                         + ParallelDescriptor::mpi_level_to_string(provided));
        }
        int myproc = ParallelDescriptor::MyProc();
        s_info = GetWriteInfo(myproc);
        MPI_Comm_split(ParallelDescriptor::Communicator(), s_info.ifile, myproc, &s_comm);
    }
#endif

    if (s_asyncout) {
        s_thread = std::make_unique<BackgroundThread>();
    }

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
    if (s_thread) {
        s_thread->Finish();
    }
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

}
