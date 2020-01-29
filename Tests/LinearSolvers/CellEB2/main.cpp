
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        int scaling_test = 0;
        amrex::ParmParse pp;
        pp.query("scaling_test", scaling_test);

        BL_PROFILE("main");
        MyTest mytest;
        mytest.solve();
        if (scaling_test) {
            BL_PROFILE_REGION("SOLVE");
            mytest.solve();
        } else {
            mytest.writePlotfile();
        }
    }

    amrex::Finalize();
}
