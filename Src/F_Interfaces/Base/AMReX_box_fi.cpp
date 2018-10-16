
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
 
    void amrex_fi_get_box_intersection(const int lo1[3], const int hi1[3], const int nodal1[3], 
                                       const int lo2[3], const int hi2[3], const int nodal2[3], 
                                       int rlo[3], int rhi[3])
    {
	Box box1 {IntVect(lo1), IntVect(hi1), IntVect(nodal1)};
	Box box2 {IntVect(lo2), IntVect(hi2), IntVect(nodal2)};

    Box bx = box1 & box2;

	const int* lov = bx.loVect();
    const int* hiv = bx.hiVect();
    for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
        rlo[idim] = lov[idim];
        rhi[idim] = hiv[idim];
    }
    }
}
