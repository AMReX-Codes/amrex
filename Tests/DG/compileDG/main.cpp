#include <iostream>

#include <AMReX.H>
#include <AMReX_DG.H>

using namespace amrex;

int main( int argc, char* argv[] )
{
    amrex::Initialize( argc, argv );
    amrex::Finalize();
}
