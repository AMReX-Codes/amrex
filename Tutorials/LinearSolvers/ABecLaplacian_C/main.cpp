
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        amrex::Gpu::setGraphRegion(false);
        MyTest mytest;
        mytest.solve();
        mytest.writePlotfile();
    }

    amrex::Finalize();
}
