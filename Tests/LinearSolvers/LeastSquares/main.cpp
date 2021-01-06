#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        MyTest mytest;

        mytest.compute_gradient();

        for (int i = 0; i < 1; ++i) {
            mytest.solve();
            mytest.writePlotfile();
        }
    }

    amrex::Finalize();
}
