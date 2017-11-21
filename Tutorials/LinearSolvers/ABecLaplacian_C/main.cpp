
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        MyTest mytest;
        mytest.solve();
        mytest.writePlotfile();
    }

    amrex::Finalize();
}
