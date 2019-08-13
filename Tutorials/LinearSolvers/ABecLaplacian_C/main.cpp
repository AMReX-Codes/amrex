
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        {
            MyTest mytest;
            mytest.solve();
        }
        {
            MyTest mytest;
            {
                BL_PROFILE_REGION("LinearSolver");
                mytest.solve();
            }
        }
//        mytest.writePlotfile();
    }

    amrex::Finalize();
}
