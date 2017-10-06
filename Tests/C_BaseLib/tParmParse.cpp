
#include <iostream>

#include <AMReX_Vector.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int
main (int argc, char** argv)
{
    amrex::Initialize(argc,argv);

    ParmParse pp;

    Vector<int> arr;

    pp.queryarr("arr",arr);

    for (int i = 0; i < arr.size(); i++)
    {
        std::cout << arr[i] << std::endl;
    }

    amrex::Finalize();
}
