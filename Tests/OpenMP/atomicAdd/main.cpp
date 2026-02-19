#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_OpenMP.H>
#include <AMReX_Print.H>

using namespace amrex;

void test_atomicAdd (MultiFab& mf)
{

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        FArrayBox tmp;
        for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
            Box const& tbx = mfi.growntilebox();
            tmp.resize(tbx);
            tmp.template setVal<RunOn::Host>(0.2);
            mf[mfi].template atomicAdd<RunOn::Host>(tmp, tbx, tbx, 0, 0, 1);
        }
    }
}

void test_lockAdd (MultiFab& mf)
{

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        FArrayBox tmp;
        for (MFIter mfi(mf,true); mfi.isValid(); ++mfi) {
            Box const& tbx = mfi.growntilebox();
            tmp.resize(tbx);
            tmp.template setVal<RunOn::Host>(0.2);
            mf[mfi].template lockAdd<RunOn::Host>(tmp, tbx, tbx, 0, 0, 1);
        }
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BoxArray ba(Box(IntVect(0),IntVect(127)));
        ba.maxSize(32);
        DistributionMapping dm(ba);

        MultiFab mf(ba,dm,1,2);
        mf.setVal(0.0);

        test_atomicAdd(mf);
        double t = amrex::second();
        test_atomicAdd(mf);
        double t_aa = amrex::second() - t;

        test_lockAdd(mf);
        t = amrex::second();
        test_lockAdd(mf);
        double t_la = amrex::second() - t;

        amrex::Print() << "  atomicAdd time is " << t_aa << "\n"
                       << "    lockAdd time is " << t_la << "\n";
    }
    amrex::Finalize();
}
