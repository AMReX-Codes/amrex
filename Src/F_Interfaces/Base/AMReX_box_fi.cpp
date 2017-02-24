
#include <AMReX_Box.H>
#include <AMReX_Print.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_print_box (const int* lo, const int* hi, int all)
    {
	Box box {IntVect(lo), IntVect(hi)};
	if (all) {
	    AllPrint() << box << "\n";
	} else {
	    Print() << box << "\n";
	}
    }
}
