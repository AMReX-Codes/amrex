
#include <iostream>

#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>

int
main (int argc, char** argv)
{
    BoxLib::Initialize(argc,argv);

    ParmParse pp;

    Array<int> arr;

    pp.queryarr("arr",arr);

    for (int i = 0; i < arr.size(); i++)
    {
        std::cout << arr[i] << std::endl;
    }

    BoxLib::Finalize();
}
