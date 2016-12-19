#include <iostream>
#include <iomanip>
#include <AMReX_Utility.H>

using namespace amrex;

// this mersenne_ran_main() outputs first 1000 generated numbers
// compare against the output of mt19937int.out
int
main(int argc, char** argv)
{
  amrex::Initialize(argc,argv);
  amrex::mt19937 rr(4357UL);
    std::ios::fmtflags ofmtflags = std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout << std::setprecision(8);
    for ( int j=0; j<1000; j++ )
    {
	std::cout << std::setw(10) << rr.u_value() << ' ';
	if ( j%5==4 ) std::cout << std::endl;
    }
    std::cout << std::endl;
    amrex::Finalize();
}
