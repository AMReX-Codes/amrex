
#include <AMReX_EB2.H>
#include "MyTest.H"

using namespace amrex;

void
MyTest::initializeEB ()
{
    EB2::Build(geom, 0, 0);
}
