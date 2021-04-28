#include <AMReX.H>
#include <AMReX_Config.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        MyTest mytest;

        mytest.compute_gradient();
        mytest.writePlotfile();
    }

    amrex::Finalize();
}
