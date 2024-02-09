#include "MyTest.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        MyTest mytest;
        mytest.solve();
    }

    amrex::Finalize();
}
