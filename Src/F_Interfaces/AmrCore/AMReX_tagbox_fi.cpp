
#include <AMReX_TagBox.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_tagboxarray_dataptr (TagBoxArray* tag, MFIter* mfi, char*& dp, int lo[3], int hi[3])
    {
	TagBox& fab = (*tag)[*mfi];
	dp = fab.dataPtr();
	const Box& bx = fab.box();
	const int* lov = bx.loVect();
	const int* hiv = bx.hiVect();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = lov[i];
	    hi[i] = hiv[i];
	}
    }
}
