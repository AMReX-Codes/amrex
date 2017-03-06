
#include <AMReX_BoxArray.H>
#include <AMReX_Print.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_boxarray (BoxArray*& ba, int lo[3], int hi[3])
    {
	IntVect small(lo), big(hi);
	ba = new BoxArray(Box(small,big));
    }

    void amrex_fi_delete_boxarray (BoxArray* ba)
    {
	delete ba;
    }

    void amrex_fi_clone_boxarray (BoxArray*& bao, const BoxArray* bai)
    {
	delete bao;
	bao = new BoxArray(*bai);
    }

    void amrex_fi_boxarray_maxsize (BoxArray* ba, int sz)
    {
	ba->maxSize(sz);
    }

    void amrex_fi_print_boxarray (const BoxArray* ba, int all)
    {
	if (all) {
	    AllPrint() << *ba;
	} else {
	    Print() << *ba;
	}
    }
}
