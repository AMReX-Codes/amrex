
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        MyTest mytest;
        bool first = true;
        for (auto const& mgt : mytest.getMultigridTypes()) {
            mytest.setMultigridType(mgt);
            if (!first) { mytest.initData(); } // restore the initial guess
            first = false;
            mytest.solve();
        }
        mytest.writePlotfile();
    }

    amrex::Finalize();
}
