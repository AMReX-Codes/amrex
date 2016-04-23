
#include <MultiFab.H>

extern "C" {

    void fi_new_multifab (void*& mf, void* ba, int nc, int ng)
    {
	mf = (void*) new MultiFab(*(BoxArray*)ba, nc, ng);
    }

    void fi_delete_multifab (void* mf)
    {
	delete (MultiFab*) mf;
    }

    void fi_multifab_dataptr (void* mf_, void* mfi_, void*& dp, int lo[3], int hi[3])
    {
	MultiFab& mf = *(MultiFab*)mf_;
	MFIter& mfi = *(MFIter*)mfi_;
	FArrayBox& fab = mf[mfi];
	dp = (void*) (fab.dataPtr());
	const Box& bx = fab.box();
	const int* lov = bx.loVect();
	const int* hiv = bx.hiVect();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = lov[i];
	    hi[i] = hiv[i];
	}
    }

    // MFIter routines

    void fi_new_mfiter (void*& mfi, void* mf, int tiling)
    {
	mfi = new MFIter(*(MultiFab*)mf, (bool)tiling);
    }

    void fi_delete_mfiter (void* mfi)
    {
	delete (MFIter*) mfi;
    }

    void fi_increment_mfiter (void* mfi_, int* isvalid)
    {
	MFIter& mfi = *(MFIter*)mfi_;
	++mfi;
	*isvalid = mfi.isValid();
    }

    void fi_mfiter_is_valid (void* mfi_, int* isvalid)
    {
	MFIter& mfi = *(MFIter*)mfi_;
	*isvalid = mfi.isValid();
    }

    void fi_mfiter_tilebox (void* mfi_, int lo[3], int hi[3])
    {
	MFIter& mfi = *(MFIter*)mfi_;
	const Box& bx = mfi.tilebox();
	const int* lov = bx.loVect();
	const int* hiv = bx.hiVect();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = lov[i];
	    hi[i] = hiv[i];
	}
    }

}

