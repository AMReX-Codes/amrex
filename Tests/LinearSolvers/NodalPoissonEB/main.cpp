
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
            for (int i = 0; i < mytest.getNumTrials(); ++i) {
                mytest.solve();
                mytest.compute_norms();
            }
        }
        if (mytest.getDoPlots()) {
            mytest.writePlotfile();
        }
    }

    amrex::Finalize();
}
