
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

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

    amrex::Finalize();
}
