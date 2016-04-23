
#include <MultiFab.H>

extern "C" {

    void fi_new_multifab (MultiFab*& mf, BoxArray*& bao, const BoxArray* bai, int nc, int ng)
    {
	mf = new MultiFab(*bai, nc, ng);
	bao = (BoxArray*)&(mf->boxArray());
    }

    void fi_delete_multifab (MultiFab* mf)
    {
	delete mf;
    }

    void fi_multifab_dataptr (MultiFab* mf, MFIter* mfi, double*& dp, int lo[3], int hi[3])
    {
	FArrayBox& fab = (*mf)[*mfi];
	dp = fab.dataPtr();
	const Box& bx = fab.box();
	const int* lov = bx.loVect();
	const int* hiv = bx.hiVect();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = lov[i];
	    hi[i] = hiv[i];
	}
    }

    // MFIter routines

    void fi_new_mfiter (MFIter*& mfi, MultiFab* mf, int tiling)
    {
	mfi = new MFIter(*mf, (bool)tiling);
    }

    void fi_delete_mfiter (MFIter* mfi)
    {
	delete mfi;
    }

    void fi_increment_mfiter (MFIter* mfi, int* isvalid)
    {
	++(*mfi);
	*isvalid = mfi->isValid();
    }

    void fi_mfiter_is_valid (MFIter* mfi, int* isvalid)
    {
	*isvalid = mfi->isValid();
    }

    void fi_mfiter_tilebox (MFIter* mfi, int lo[3], int hi[3])
    {
	const Box& bx = mfi->tilebox();
	const int* lov = bx.loVect();
	const int* hiv = bx.hiVect();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = lov[i];
	    hi[i] = hiv[i];
	}
    }

}

