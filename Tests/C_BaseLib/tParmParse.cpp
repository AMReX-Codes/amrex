
//
// $Id: tParmParse.cpp,v 1.1 2004-01-21 21:04:13 lijewski Exp $
//

#include <iostream>

#include <Array.H>
#include <ParmParse.H>

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
