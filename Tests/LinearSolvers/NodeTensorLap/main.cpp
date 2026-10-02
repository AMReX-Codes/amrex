
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        MyTest mytest;
        for (auto const& mgt : mytest.getMultigridTypes()) {
            mytest.setMultigridType(mgt);
            mytest.solve();
            mytest.compute_norms();
        }
        mytest.writePlotfile();
    }

    amrex::Finalize();
}
