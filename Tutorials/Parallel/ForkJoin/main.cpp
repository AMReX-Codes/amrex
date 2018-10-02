
#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        MyTest mytest;
        mytest.runTest();
    }

    amrex::Finalize();
}
