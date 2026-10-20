#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    // split_comm = 1: run the test on each half of the ranks (MPI builds
    // with at least two ranks), e.g. to check solvers on sub-communicators
    bool split_comm = false;
    {
        amrex::ParmParse pp;
        pp.query("split_comm", split_comm);
        // The halves would write the same plotfile
        if (split_comm) { pp.add("do_plots", false); }
    }
#ifdef AMREX_USE_MPI
    MPI_Comm subcomm = MPI_COMM_NULL;
    if (split_comm) {
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const int myproc = amrex::ParallelDescriptor::MyProc();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nprocs >= 2, "split_comm needs at least two ranks");
        const int color = (myproc < nprocs/2) ? 0 : 1;
        MPI_Comm_split(amrex::ParallelDescriptor::Communicator(), color, myproc, &subcomm);
        amrex::ParallelContext::push(subcomm);
    }
#else
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!split_comm, "split_comm needs MPI");
#endif

    {
        BL_PROFILE("main");
        MyTest mytest;
        bool first = true;
        for (auto const& mgt : mytest.getMultigridTypes()) {
            if (!mytest.setMultigridType(mgt)) { continue; }
            if (!first) { mytest.initData(); } // restore the initial guess
            first = false;
            mytest.solve();
        }
        mytest.writePlotfile();
    }

#ifdef AMREX_USE_MPI
    if (split_comm) {
        amrex::ParallelContext::pop();
        MPI_Comm_free(&subcomm);
    }
#endif

    amrex::Finalize();
}
