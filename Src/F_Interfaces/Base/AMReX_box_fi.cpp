
#include <AMReX_Box.H>
#include <AMReX_Print.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_print_box (const int* lo, const int* hi, const int* nodal)
    {
	Box box {IntVect(lo), IntVect(hi), IntVect(nodal)};
	AllPrint() << box << "\n";
    }
}
